"""
AI-powered result interpretation module for AllEnricher v2.0

Supports multiple AI backends:
- OpenAI (GPT-4, GPT-3.5)
- Anthropic (Claude)
- Local models (via Ollama)

中文模块说明：
    AI驱动的结果解读模块（AllEnricher v2.0）

    本模块提供基于人工智能的基因集富集分析结果解读功能，支持多种AI后端：
    - OpenAI（GPT-4、GPT-3.5）：云端大语言模型，提供高质量的生物学解读
    - Anthropic（Claude）：Anthropic公司的Claude系列模型，擅长长文本分析
    - Ollama：本地部署的开源模型，无需API密钥，适合离线或隐私敏感场景
    - Mock：测试用模拟后端，无需任何AI服务即可运行

    模块架构：
    - AIInterpreterBase：抽象基类，定义统一的解读接口
    - OpenAIInterpreter / ClaudeInterpreter / OllamaInterpreter：具体后端实现
    - MockInterpreter：用于测试的模拟实现
    - AIInterpreter：主入口类（门面模式），统一管理不同后端
    - create_interpreter()：工厂函数，便捷创建解释器实例
    - get_available_backends()：获取所有可用后端列表
"""

import os
import json
import logging
from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Any
import pandas as pd

logger = logging.getLogger(__name__)


class AIInterpreterBase(ABC):
    """
    AI解释器抽象基类

    定义了所有AI解释器后端必须实现的统一接口。
    所有具体的后端实现（OpenAI、Claude、Ollama、Mock）都需要继承此类，
    并实现 interpret() 和 summarize_term() 两个抽象方法。

    设计模式：策略模式（Strategy Pattern），通过统一的接口支持不同的AI后端切换。
    """

    @abstractmethod
    def interpret(self, results: Dict[str, pd.DataFrame], context: str = "") -> Dict[str, str]:
        """
        对富集分析结果生成AI解读

        参数:
            results (Dict[str, pd.DataFrame]): 富集分析结果字典，
                键为数据库名称（如 "GO_Biological_Process"），值为对应的DataFrame结果
            context (str): 额外的分析上下文信息，用于帮助AI更好地理解分析背景，默认为空字符串

        返回:
            Dict[str, str]: 解读结果字典，键为数据库名称，值为AI生成的解读文本
        """
        pass

    @abstractmethod
    def summarize_term(self, term_name: str, gene_list: List[str]) -> str:
        """
        对单个富集条目进行简要总结

        参数:
            term_name (str): 富集条目名称（如GO术语、KEGG通路名称等）
            gene_list (List[str]): 与该条目关联的基因列表

        返回:
            str: AI生成的该条目的生物学意义总结文本
        """
        pass


class OpenAIInterpreter(AIInterpreterBase):
    """
    OpenAI GPT后端解释器

    使用OpenAI的GPT系列模型（如GPT-4、GPT-3.5）对富集分析结果进行AI解读。
    需要提供有效的OpenAI API密钥，可通过构造函数参数或环境变量 OPENAI_API_KEY 配置。

    特点：
    - 支持GPT-4和GPT-3.5等多种模型
    - 可自定义生成参数（最大token数、温度等）
    - 使用系统提示词设定生物信息学专家角色
    """

    def __init__(
        self,
        api_key: str = None,
        model: str = "gpt-4",
        max_tokens: int = 2000,
        temperature: float = 0.7
    ):
        """
        初始化OpenAI解释器

        参数:
            api_key (str): OpenAI API密钥，若为None则从环境变量 OPENAI_API_KEY 获取
            model (str): 使用的模型名称，默认为 "gpt-4"，可选 "gpt-3.5-turbo" 等
            max_tokens (int): 生成文本的最大token数量，默认为2000
            temperature (float): 生成温度，控制输出随机性，范围0-1，默认为0.7
        """
        self.api_key = api_key or os.getenv("OPENAI_API_KEY")
        self.model = model
        self.max_tokens = max_tokens
        self.temperature = temperature

        # 未提供API密钥时发出警告，后续调用将跳过AI解读
        if not self.api_key:
            logger.warning("OpenAI API key not provided. AI interpretation will be disabled.")

    def _call_api(self, prompt: str) -> str:
        """
        调用OpenAI ChatCompletion API

        构建包含系统角色和用户提示的消息列表，发送至OpenAI API并返回生成的文本内容。
        系统提示将AI角色设定为"基因集富集分析领域的生物信息学专家"。

        参数:
            prompt (str): 用户提示词，包含需要解读的富集分析数据

        返回:
            str: AI生成的解读文本；若调用失败则返回错误信息字符串
        """
        try:
            import openai

            # 创建OpenAI客户端实例（openai>=1.0.0 新版API）
            client = openai.OpenAI(api_key=self.api_key)

            # 调用ChatCompletion接口，使用系统提示词设定专家角色
            response = client.chat.completions.create(
                model=self.model,
                messages=[
                    # 系统提示：设定AI为生物信息学专家角色
                    {"role": "system", "content": "You are a expert bioinformatician specializing in gene set enrichment analysis. Provide clear, accurate, and insightful interpretations of enrichment results."},
                    {"role": "user", "content": prompt}
                ],
                max_tokens=self.max_tokens,
                temperature=self.temperature
            )

            # 提取并返回生成的文本内容
            return response.choices[0].message.content

        except ImportError:
            logger.error("openai package not installed. Run: pip install openai")
            return "Error: openai package not installed"
        except Exception as e:
            logger.error(f"OpenAI API error: {e}")
            return f"Error: {str(e)}"

    def interpret(self, results: Dict[str, pd.DataFrame], context: str = "") -> Dict[str, str]:
        """
        对所有富集分析结果生成AI解读

        遍历每个数据库的富集结果，提取前10个最显著的富集条目，
        构建包含P值和基因数量的详细提示词，调用OpenAI API生成生物学解读。

        参数:
            results (Dict[str, pd.DataFrame]): 富集分析结果字典，键为数据库名称，值为结果DataFrame
            context (str): 额外的分析上下文信息，默认为空字符串

        返回:
            Dict[str, str]: 解读结果字典，键为数据库名称，值为AI生成的生物学解读文本
        """
        interpretations = {}

        # 若未配置API密钥，直接返回空字典
        if not self.api_key:
            return interpretations

        for db_name, df in results.items():
            # 跳过空结果
            if len(df) == 0:
                continue

            # 提取前10个最显著的富集结果，构建摘要信息
            top_results = df.head(10)
            summary_lines = []

            for _, row in top_results.iterrows():
                summary_lines.append(
                    f"- {row.get('Term_Name', 'N/A')}: "
                    f"P-value={row.get('P_Value', 1):.2e}, "
                    f"Genes={row.get('Gene_Count', 0)}"
                )

            # 构建详细的提示词，包含分析上下文和具体要求
            prompt = f"""
Please provide a biological interpretation of the following {db_name} enrichment analysis results.

Context: {context if context else "Gene set enrichment analysis"}

Top 10 enriched terms:
{chr(10).join(summary_lines)}

Please:
1. Summarize the main biological themes represented by these enriched terms
2. Identify any potential biological processes or pathways that are significantly overrepresented
3. Suggest potential biological implications of these findings
4. Note any interesting patterns or relationships between the enriched terms

Keep the interpretation concise (2-3 paragraphs) and focus on biological insights.
"""

            # 调用API并存储解读结果
            interpretation = self._call_api(prompt)
            interpretations[db_name] = interpretation

        return interpretations

    def summarize_term(self, term_name: str, gene_list: List[str]) -> str:
        """
        对单个富集条目生成简要总结

        构建包含条目名称和关联基因（最多显示前10个）的提示词，
        调用OpenAI API生成2-3句话的生物学意义总结。

        参数:
            term_name (str): 富集条目名称
            gene_list (List[str]): 与该条目关联的基因列表

        返回:
            str: AI生成的条目总结文本；若未配置API密钥则返回空字符串
        """
        if not self.api_key:
            return ""

        # 构建提示词，包含条目名称和关联基因（超过10个基因时截断并添加省略号）
        prompt = f"""
Please provide a brief description of the following biological term and its relevance:

Term: {term_name}
Associated genes: {', '.join(gene_list[:10])}{'...' if len(gene_list) > 10 else ''}

Provide a 2-3 sentence summary explaining what this term represents and its biological significance.
"""

        return self._call_api(prompt)


class ClaudeInterpreter(AIInterpreterBase):
    """
    Anthropic Claude后端解释器

    使用Anthropic公司的Claude系列模型（如claude-3-opus）对富集分析结果进行AI解读。
    需要提供有效的Anthropic API密钥，可通过构造函数参数或环境变量 ANTHROPIC_API_KEY 配置。

    特点：
    - 使用Anthropic Messages API
    - 擅长长文本分析和细致的推理
    - 支持claude-3-opus、claude-3-sonnet等多种模型
    """

    def __init__(
        self,
        api_key: str = None,
        model: str = "claude-3-opus-20240229",
        max_tokens: int = 2000
    ):
        """
        初始化Claude解释器

        参数:
            api_key (str): Anthropic API密钥，若为None则从环境变量 ANTHROPIC_API_KEY 获取
            model (str): 使用的模型名称，默认为 "claude-3-opus-20240229"
            max_tokens (int): 生成文本的最大token数量，默认为2000
        """
        self.api_key = api_key or os.getenv("ANTHROPIC_API_KEY")
        self.model = model
        self.max_tokens = max_tokens

        # 未提供API密钥时发出警告，后续调用将跳过AI解读
        if not self.api_key:
            logger.warning("Anthropic API key not provided. AI interpretation will be disabled.")

    def _call_api(self, prompt: str) -> str:
        """
        调用Anthropic Messages API

        使用Anthropic SDK创建客户端并发送消息请求，返回Claude生成的文本内容。

        参数:
            prompt (str): 用户提示词，包含需要解读的富集分析数据

        返回:
            str: Claude生成的解读文本；若调用失败则返回错误信息字符串
        """
        try:
            import anthropic

            # 创建Anthropic客户端实例
            client = anthropic.Anthropic(api_key=self.api_key)

            # 调用Messages API生成回复
            message = client.messages.create(
                model=self.model,
                max_tokens=self.max_tokens,
                messages=[
                    {"role": "user", "content": prompt}
                ]
            )

            # 提取并返回生成的文本内容（content数组的第一个元素）
            return message.content[0].text

        except ImportError:
            logger.error("anthropic package not installed. Run: pip install anthropic")
            return "Error: anthropic package not installed"
        except Exception as e:
            logger.error(f"Anthropic API error: {e}")
            return f"Error: {str(e)}"

    def interpret(self, results: Dict[str, pd.DataFrame], context: str = "") -> Dict[str, str]:
        """
        对所有富集分析结果生成AI解读

        遍历每个数据库的富集结果，提取前10个最显著的富集条目，
        构建提示词并调用Claude API生成简洁的生物学解读。

        参数:
            results (Dict[str, pd.DataFrame]): 富集分析结果字典，键为数据库名称，值为结果DataFrame
            context (str): 额外的分析上下文信息，默认为空字符串

        返回:
            Dict[str, str]: 解读结果字典，键为数据库名称，值为Claude生成的生物学解读文本
        """
        interpretations = {}

        # 若未配置API密钥，直接返回空字典
        if not self.api_key:
            return interpretations

        for db_name, df in results.items():
            # 跳过空结果
            if len(df) == 0:
                continue

            # 提取前10个最显著的富集结果，构建摘要信息（包含P值）
            top_results = df.head(10)
            summary_lines = []

            for _, row in top_results.iterrows():
                summary_lines.append(
                    f"- {row.get('Term_Name', 'N/A')}: "
                    f"P-value={row.get('P_Value', 1):.2e}"
                )

            # 构建提示词，将生物信息学专家角色嵌入用户提示中
            prompt = f"""
You are a bioinformatics expert. Please interpret the following {db_name} enrichment results:

{chr(10).join(summary_lines)}

Provide a concise biological interpretation (2-3 paragraphs).
"""

            # 调用API并存储解读结果
            interpretation = self._call_api(prompt)
            interpretations[db_name] = interpretation

        return interpretations

    def summarize_term(self, term_name: str, gene_list: List[str]) -> str:
        """
        对单个富集条目生成简要总结

        构建包含条目名称的简洁提示词，调用Claude API生成2-3句话的生物学意义总结。

        参数:
            term_name (str): 富集条目名称
            gene_list (List[str]): 与该条目关联的基因列表（本实现中未直接使用）

        返回:
            str: Claude生成的条目总结文本；若未配置API密钥则返回空字符串
        """
        if not self.api_key:
            return ""

        # 构建简洁的提示词，仅包含条目名称
        prompt = f"Describe the biological term '{term_name}' and its significance in 2-3 sentences."
        return self._call_api(prompt)


class OllamaInterpreter(AIInterpreterBase):
    """
    本地Ollama后端解释器

    使用本地部署的Ollama服务运行开源大语言模型（如llama2、mistral等）进行富集分析结果解读。
    无需API密钥，适合离线环境或对数据隐私有较高要求的场景。

    使用前提：
    - 需要本地安装并运行Ollama服务（默认地址 http://localhost:11434）
    - 需要预先拉取所需的模型（如 ollama pull llama2）

    特点：
    - 完全本地运行，数据不离开本机
    - 无需API密钥和付费订阅
    - 支持多种开源模型
    - 仅提取前5个富集条目（相比云端后端更少），以适配本地模型能力
    """

    def __init__(self, model: str = "llama2", base_url: str = "http://localhost:11434"):
        """
        初始化Ollama解释器

        参数:
            model (str): 使用的模型名称，默认为 "llama2"，需确保该模型已通过Ollama拉取
            base_url (str): Ollama服务的API地址，默认为 "http://localhost:11434"
        """
        self.model = model
        self.base_url = base_url

    def _call_api(self, prompt: str) -> str:
        """
        调用Ollama本地API

        通过HTTP POST请求调用Ollama的 /api/generate 接口，
        设置 stream=False 以一次性获取完整响应。

        参数:
            prompt (str): 用户提示词，包含需要解读的富集分析数据

        返回:
            str: 模型生成的解读文本；若调用失败则返回错误信息字符串
        """
        try:
            import requests

            # 向Ollama的生成接口发送POST请求，stream=False表示非流式返回完整结果
            response = requests.post(
                f"{self.base_url}/api/generate",
                json={
                    "model": self.model,
                    "prompt": prompt,
                    "stream": False  # 非流式模式，等待完整响应返回
                }
            )

            # 检查HTTP状态码，成功时提取生成的文本
            if response.status_code == 200:
                return response.json().get("response", "")
            else:
                return f"Error: {response.status_code}"

        except ImportError:
            logger.error("requests package not installed. Run: pip install requests")
            return "Error: requests package not installed"
        except Exception as e:
            logger.error(f"Ollama API error: {e}")
            return f"Error: {str(e)}"

    def interpret(self, results: Dict[str, pd.DataFrame], context: str = "") -> Dict[str, str]:
        """
        对所有富集分析结果生成AI解读

        遍历每个数据库的富集结果，提取前5个最显著的富集条目（适配本地模型能力），
        构建简洁的提示词并调用Ollama API生成生物学解读。

        参数:
            results (Dict[str, pd.DataFrame]): 富集分析结果字典，键为数据库名称，值为结果DataFrame
            context (str): 额外的分析上下文信息，默认为空字符串（本实现中未使用）

        返回:
            Dict[str, str]: 解读结果字典，键为数据库名称，值为模型生成的生物学解读文本
        """
        interpretations = {}

        for db_name, df in results.items():
            # 跳过空结果
            if len(df) == 0:
                continue

            # 提取前5个最显著的富集结果（比云端后端少，适配本地模型能力）
            top_results = df.head(5)
            summary_lines = []

            for _, row in top_results.iterrows():
                # 仅提取条目名称，不包含详细统计信息
                summary_lines.append(f"- {row.get('Term_Name', 'N/A')}")

            # 构建简洁的提示词
            prompt = f"""
Interpret these {db_name} enrichment results:
{chr(10).join(summary_lines)}

Provide a brief biological interpretation.
"""

            # 调用API并存储解读结果
            interpretation = self._call_api(prompt)
            interpretations[db_name] = interpretation

        return interpretations

    def summarize_term(self, term_name: str, gene_list: List[str]) -> str:
        """
        对单个富集条目生成简要总结

        构建极简提示词，调用Ollama API生成条目描述。

        参数:
            term_name (str): 富集条目名称
            gene_list (List[str]): 与该条目关联的基因列表（本实现中未直接使用）

        返回:
            str: 模型生成的条目描述文本
        """
        # 构建极简提示词，仅包含条目名称
        prompt = f"Describe: {term_name}"
        return self._call_api(prompt)


class MockInterpreter(AIInterpreterBase):
    """
    测试用模拟解释器

    无需任何AI服务或API密钥，生成固定格式的模拟解读结果。
    主要用于：
    - 单元测试和集成测试
    - 演示和开发调试
    - 在没有AI服务时的降级方案

    注意：生成的解读内容为模板化文本，不包含真实的生物学分析。
    """

    def interpret(self, results: Dict[str, pd.DataFrame], context: str = "") -> Dict[str, str]:
        """
        生成模拟的富集分析解读

        为每个数据库生成包含以下内容的模板化解读：
        - 数据库名称和富集条目总数
        - 前5个最显著的富集条目列表
        - 通用的生物学解读模板
        - 分析建议

        参数:
            results (Dict[str, pd.DataFrame]): 富集分析结果字典，键为数据库名称，值为结果DataFrame
            context (str): 额外的分析上下文信息，默认为空字符串（本实现中未使用）

        返回:
            Dict[str, str]: 模拟解读结果字典，键为数据库名称，值为模板化的解读文本
        """
        interpretations = {}

        for db_name, df in results.items():
            # 跳过空结果
            if len(df) == 0:
                continue

            # 提取前5个富集条目名称
            top_terms = df.head(5)['Term_Name'].tolist()

            # 生成模板化的模拟解读（Markdown格式）
            interpretation = f"""
**{db_name} Enrichment Analysis Summary**

The gene set shows significant enrichment in {len(df)} terms from the {db_name} database.

**Top enriched terms:**
{chr(10).join([f"- {term}" for term in top_terms])}

**Biological Interpretation:**
The enrichment results suggest that the input gene set is involved in several biological processes. 
The top enriched terms indicate potential functional themes that may be relevant to the biological 
context of your study.

**Recommendations:**
1. Validate key findings with literature review
2. Consider pathway analysis for deeper insights
3. Cross-reference with other omics data if available

*Note: This is a mock interpretation. Enable AI integration for detailed analysis.*
"""

            interpretations[db_name] = interpretation

        return interpretations

    def summarize_term(self, term_name: str, gene_list: List[str]) -> str:
        """
        生成模拟的单个条目总结

        返回包含条目名称和关联基因数量的简单描述。

        参数:
            term_name (str): 富集条目名称
            gene_list (List[str]): 与该条目关联的基因列表

        返回:
            str: 模拟的条目总结文本，格式为 "The term '{term_name}' is associated with {n} genes..."
        """
        return f"The term '{term_name}' is associated with {len(gene_list)} genes from your input set."


class AIInterpreter:
    """
    AI解释器主入口类（门面模式 / Facade Pattern）

    作为模块的统一对外接口，屏蔽不同AI后端的实现细节。
    用户只需通过此类即可使用所有支持的AI后端进行富集分析结果解读。

    使用方式：
        # 使用OpenAI后端
        interpreter = AIInterpreter(backend="openai", api_key="sk-xxx")
        result = interpreter.interpret_results(enrichment_results)

        # 使用本地Ollama后端
        interpreter = AIInterpreter(backend="ollama", model="mistral")
        result = interpreter.interpret_results(enrichment_results)

        # 使用测试模拟后端
        interpreter = AIInterpreter(backend="mock")
        result = interpreter.interpret_results(enrichment_results)

    支持的后端：
        - "openai": OpenAI GPT系列（需要API密钥）
        - "claude": Anthropic Claude系列（需要API密钥）
        - "ollama": 本地Ollama模型（无需API密钥）
        - "mock": 测试模拟后端（无需任何配置）
    """

    # 后端名称到实现类的映射表
    BACKENDS = {
        "openai": OpenAIInterpreter,
        "claude": ClaudeInterpreter,
        "ollama": OllamaInterpreter,
        "mock": MockInterpreter
    }

    def __init__(
        self,
        backend: str = "openai",
        api_key: str = None,
        model: str = None,
        **kwargs
    ):
        """
        初始化AI解释器

        根据指定的后端名称创建对应的解释器实例。
        不同后端使用不同的默认模型：OpenAI默认gpt-4，Claude默认claude-3-opus，Ollama默认llama2。

        参数:
            backend (str): AI后端名称，可选 "openai"、"claude"、"ollama"、"mock"，默认为 "openai"
            api_key (str): API密钥，用于openai和claude后端；若为None则从对应环境变量获取
            model (str): 模型名称，若为None则使用各后端的默认模型
            **kwargs: 传递给具体后端构造函数的额外参数（如 max_tokens、temperature 等）

        异常:
            ValueError: 当指定的后端名称不在支持列表中时抛出
        """
        if backend not in self.BACKENDS:
            raise ValueError(f"Unknown backend: {backend}. Available: {list(self.BACKENDS.keys())}")

        self.backend_name = backend
        interpreter_class = self.BACKENDS[backend]

        # 根据后端类型，使用适当的参数初始化对应的解释器实例
        if backend in ["openai", "claude"]:
            # 云端后端：需要API密钥，使用对应默认模型
            self.interpreter = interpreter_class(
                api_key=api_key,
                model=model or ("gpt-4" if backend == "openai" else "claude-3-opus-20240229"),
                **kwargs
            )
        elif backend == "ollama":
            # 本地Ollama后端：无需API密钥，默认使用llama2模型
            self.interpreter = interpreter_class(
                model=model or "llama2",
                **kwargs
            )
        else:
            # Mock后端：无需任何参数
            self.interpreter = interpreter_class(**kwargs)
    
    def interpret_results(
        self,
        results: Dict[str, pd.DataFrame],
        context: str = "",
        include_term_summaries: bool = False
    ) -> Dict[str, Any]:
        """
        生成富集分析结果的综合AI解读

        首先调用后端解释器生成各数据库的整体解读，可选地为每个数据库的前5个富集条目
        生成单独的条目总结。

        参数:
            results (Dict[str, pd.DataFrame]): 富集分析结果字典，键为数据库名称，值为结果DataFrame
            context (str): 额外的分析上下文信息（如实验设计、研究目的等），默认为空字符串
            include_term_summaries (bool): 是否为每个数据库的前5个富集条目生成单独总结，默认为False

        返回:
            Dict[str, Any]: 综合解读结果字典，包含：
                - 各数据库名称对应的整体解读文本
                - 若 include_term_summaries=True，还包含 "{数据库名}_term_summaries" 键，
                  值为该数据库前5个条目的总结字典
        """
        # 调用后端解释器生成各数据库的整体解读
        interpretations = self.interpreter.interpret(results, context)

        # 可选：为每个数据库的前5个富集条目生成单独的条目总结
        if include_term_summaries:
            for db_name, df in results.items():
                if len(df) > 0:
                    term_summaries = {}
                    # 遍历前5个富集条目
                    for _, row in df.head(5).iterrows():
                        term_name = row.get('Term_Name', '')
                        # 从 'Genes' 列提取基因列表（分号分隔）
                        genes = row.get('Genes', '').split(';')
                        if term_name:
                            term_summaries[term_name] = self.interpreter.summarize_term(term_name, genes)

                    # 将条目总结以 "{数据库名}_term_summaries" 为键存入结果
                    interpretations[f"{db_name}_term_summaries"] = term_summaries

        return interpretations
    
    def generate_report_section(
        self,
        results: Dict[str, pd.DataFrame],
        context: str = ""
    ) -> str:
        """
        生成AI解读的HTML报告段落

        将AI解读结果渲染为HTML格式的报告段落，包含：
        - 段落标题和AI免责声明（提示用户需由领域专家审核）
        - 各数据库的解读内容（跳过条目总结，仅展示整体解读）
        - 使用Font Awesome图标增强视觉效果

        参数:
            results (Dict[str, pd.DataFrame]): 富集分析结果字典，键为数据库名称，值为结果DataFrame
            context (str): 额外的分析上下文信息，默认为空字符串

        返回:
            str: 完整的HTML字符串，可直接嵌入报告页面中
        """
        # 获取AI解读结果
        interpretations = self.interpret_results(results, context)

        # 构建HTML段落的开头：包含标题和AI免责声明
        html_parts = ['''
        <div class="section" id="ai-interpretation">
            <h2><i class="fas fa-brain"></i> AI-Powered Interpretation</h2>
            <p class="ai-disclaimer">
                <i class="fas fa-info-circle"></i>
                The following interpretations are generated by AI ({}) and should be reviewed by domain experts.
            </p>
        '''.format(self.backend_name)]

        # 遍历解读结果，为每个数据库生成HTML展示块
        for db_name, interpretation in interpretations.items():
            # 跳过条目总结（仅展示整体解读）
            if db_name.endswith('_term_summaries'):
                continue

            # 将解读文本中的换行符替换为HTML换行标签
            html_parts.append(f'''
            <div class="ai-interpretation">
                <h3><i class="fas fa-robot"></i> {db_name}</h3>
                <div class="interpretation-content">
                    {interpretation.replace(chr(10), '<br>')}
                </div>
            </div>
            ''')

        # 闭合HTML段落
        html_parts.append('</div>')
        return ''.join(html_parts)


def get_available_backends() -> List[str]:
    """
    获取所有可用的AI后端名称列表

    返回:
        List[str]: 支持的后端名称列表，包含 ["openai", "claude", "ollama", "mock"]
    """
    return list(AIInterpreter.BACKENDS.keys())


def create_interpreter(
    backend: str = "mock",
    api_key: str = None,
    model: str = None,
    **kwargs
) -> AIInterpreter:
    """
    工厂函数：创建AI解释器实例

    提供便捷的方式来创建AIInterpreter实例，无需直接导入AIInterpreter类。
    默认使用 "mock" 后端，确保在没有AI服务配置的情况下也能正常使用。

    参数:
        backend (str): AI后端名称，可选 "openai"、"claude"、"ollama"、"mock"，默认为 "mock"
        api_key (str): API密钥，用于openai和claude后端
        model (str): 模型名称，若为None则使用各后端的默认模型
        **kwargs: 传递给具体后端构造函数的额外参数

    返回:
        AIInterpreter: 初始化完成的AI解释器实例

    使用示例:
        # 快速创建一个OpenAI解释器
        interpreter = create_interpreter(backend="openai", api_key="sk-xxx")

        # 创建一个本地Ollama解释器
        interpreter = create_interpreter(backend="ollama", model="mistral")

        # 创建一个测试用模拟解释器
        interpreter = create_interpreter()
    """
    return AIInterpreter(backend=backend, api_key=api_key, model=model, **kwargs)
