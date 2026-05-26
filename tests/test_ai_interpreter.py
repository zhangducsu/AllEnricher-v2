"""
AI 解读模块全面单元测试

覆盖所有后端解释器、门面类、工厂函数的测试。
使用 unittest.mock 模拟外部依赖（openai、anthropic、requests），
确保测试可独立运行，无需真实的 API 密钥或网络连接。
"""

import pytest
import sys
from pathlib import Path
from unittest.mock import patch, MagicMock, PropertyMock
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))


# ---------------------------------------------------------------------------
# 测试数据工厂
# ---------------------------------------------------------------------------

def _make_results(n_terms=5, include_genes_col=False):
    """创建模拟富集分析结果 DataFrame。

    参数:
        n_terms: 富集条目数量
        include_genes_col: 是否包含 'Genes' 列（分号分隔的基因字符串）
    """
    data = {
        "Term_ID": [f"GO:{1000 + i}" for i in range(n_terms)],
        "Term_Name": [f"biological_process_{i}" for i in range(n_terms)],
        "Gene_Count": [10 - i for i in range(n_terms)],
        "P_Value": [1e-5 * (i + 1) for i in range(n_terms)],
        "Adjusted_P_Value": [1e-3 * (i + 1) for i in range(n_terms)],
    }
    if include_genes_col:
        data["Genes"] = [";".join([f"GENE{j}" for j in range(5 - i)]) for i in range(n_terms)]
    df = pd.DataFrame(data)
    return {"GO_Biological_Process": df}


def _make_empty_results():
    """创建空结果的 DataFrame。"""
    return {"GO_Biological_Process": pd.DataFrame()}


# ===========================================================================
# 1. MockInterpreter 测试
# ===========================================================================

class TestMockInterpreter:
    """测试 MockInterpreter 模拟解释器。"""

    def test_interpret_normal(self):
        """interpret() 正常工作（传入 mock DataFrame）。"""
        from allenricher.ai.interpreter import MockInterpreter

        interpreter = MockInterpreter()
        results = _make_results(n_terms=5)
        interpretations = interpreter.interpret(results)

        assert "GO_Biological_Process" in interpretations
        assert "enrichment" in interpretations["GO_Biological_Process"].lower()
        assert "5 terms" in interpretations["GO_Biological_Process"]

    def test_interpret_empty_results(self):
        """interpret() 空结果时不返回该数据库的解读。"""
        from allenricher.ai.interpreter import MockInterpreter

        interpreter = MockInterpreter()
        results = _make_empty_results()
        interpretations = interpreter.interpret(results)

        # MockInterpreter 对空结果直接 continue，不放入字典
        assert "GO_Biological_Process" not in interpretations

    def test_interpret_fewer_than_20(self):
        """interpret() 不足 20 条时按实际数量展示。"""
        from allenricher.ai.interpreter import MockInterpreter

        interpreter = MockInterpreter()
        n = 8
        results = _make_results(n_terms=n)
        interpretations = interpreter.interpret(results)

        text = interpretations["GO_Biological_Process"]
        # 应包含全部 n 个条目名称
        for i in range(n):
            assert f"biological_process_{i}" in text

    def test_interpret_exactly_20(self):
        """interpret() 恰好 20 条时全部展示。"""
        from allenricher.ai.interpreter import MockInterpreter

        interpreter = MockInterpreter()
        results = _make_results(n_terms=20)
        interpretations = interpreter.interpret(results)

        text = interpretations["GO_Biological_Process"]
        for i in range(20):
            assert f"biological_process_{i}" in text

    def test_interpret_more_than_20(self):
        """interpret() 超过 20 条时只展示前 20 条。"""
        from allenricher.ai.interpreter import MockInterpreter

        interpreter = MockInterpreter()
        results = _make_results(n_terms=25)
        interpretations = interpreter.interpret(results)

        text = interpretations["GO_Biological_Process"]
        # 前 20 个应出现
        for i in range(20):
            assert f"biological_process_{i}" in text
        # 第 21 个不应出现（head(20) 截断）
        assert "biological_process_20" not in text

    def test_summarize_term(self):
        """summarize_term() 正常工作。"""
        from allenricher.ai.interpreter import MockInterpreter

        interpreter = MockInterpreter()
        summary = interpreter.summarize_term("cell division", ["GENE1", "GENE2", "GENE3"])

        assert "cell division" in summary
        assert "3 genes" in summary

    def test_summarize_term_empty_genes(self):
        """summarize_term() 空基因列表。"""
        from allenricher.ai.interpreter import MockInterpreter

        interpreter = MockInterpreter()
        summary = interpreter.summarize_term("apoptosis", [])

        assert "apoptosis" in summary
        assert "0 genes" in summary

    def test_create_interpreter_mock(self):
        """create_interpreter("mock") 返回 MockInterpreter 后端。"""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter("mock")
        assert interpreter.backend_name == "mock"
        assert isinstance(interpreter.interpreter, type(interpreter.interpreter))


# ===========================================================================
# 2. OpenAIInterpreter 测试
# ===========================================================================

class TestOpenAIInterpreter:
    """测试 OpenAIInterpreter（mock openai 包）。"""

    def test_call_api_success(self):
        """_call_api() 正常调用返回结果。"""
        from allenricher.ai.interpreter import OpenAIInterpreter

        interpreter = OpenAIInterpreter(api_key="test-key")

        mock_response = MagicMock()
        mock_response.choices = [MagicMock()]
        mock_response.choices[0].message.content = "This is a test interpretation."

        with patch.dict("sys.modules", {"openai": MagicMock()}):
            import openai as mock_openai
            mock_client = MagicMock()
            mock_client.chat.completions.create.return_value = mock_response
            mock_openai.OpenAI.return_value = mock_client

            result = interpreter._call_api("test prompt")
            assert result == "This is a test interpretation."
            mock_client.chat.completions.create.assert_called_once()

    def test_call_api_import_error(self):
        """_call_api() ImportError 降级处理。"""
        from allenricher.ai.interpreter import OpenAIInterpreter

        interpreter = OpenAIInterpreter(api_key="test-key")

        with patch.dict("sys.modules", {"openai": None}):
            # 让 import openai 抛出 ImportError
            with patch("builtins.__import__", side_effect=ImportError("No module named openai")):
                result = interpreter._call_api("test prompt")
                assert "Error" in result
                assert "not installed" in result

    def test_call_api_exception(self):
        """_call_api() 其他异常处理。"""
        from allenricher.ai.interpreter import OpenAIInterpreter

        interpreter = OpenAIInterpreter(api_key="test-key")

        with patch.dict("sys.modules", {"openai": MagicMock()}):
            import openai as mock_openai
            mock_openai.OpenAI.side_effect = Exception("API connection failed")

            result = interpreter._call_api("test prompt")
            assert "Error" in result
            assert "API connection failed" in result

    def test_interpret_extracts_top_20(self):
        """interpret() 提取 Top 20 条目。"""
        from allenricher.ai.interpreter import OpenAIInterpreter

        interpreter = OpenAIInterpreter(api_key="test-key")
        results = _make_results(n_terms=25)

        mock_response = MagicMock()
        mock_response.choices = [MagicMock()]
        mock_response.choices[0].message.content = "Interpretation of top terms."

        with patch.dict("sys.modules", {"openai": MagicMock()}):
            import openai as mock_openai
            mock_client = MagicMock()
            mock_client.chat.completions.create.return_value = mock_response
            mock_openai.OpenAI.return_value = mock_client

            interpretations = interpreter.interpret(results)

            assert "GO_Biological_Process" in interpretations
            # 验证 prompt 中只包含前 20 条
            call_args = mock_client.chat.completions.create.call_args
            prompt = call_args[1]["messages"][1]["content"]
            assert "Top 20" in prompt
            assert "biological_process_0" in prompt
            assert "biological_process_19" in prompt
            assert "biological_process_20" not in prompt

    def test_interpret_empty_results_no_api_call(self):
        """interpret() 空结果不调用 API。"""
        from allenricher.ai.interpreter import OpenAIInterpreter

        interpreter = OpenAIInterpreter(api_key="test-key")
        results = _make_empty_results()

        with patch.object(interpreter, "_call_api") as mock_call:
            interpretations = interpreter.interpret(results)

            assert "GO_Biological_Process" in interpretations
            assert "No significant enrichment" in interpretations["GO_Biological_Process"]
            mock_call.assert_not_called()

    def test_interpret_no_api_key(self):
        """interpret() 无 API key 时返回空字典。"""
        from allenricher.ai.interpreter import OpenAIInterpreter

        interpreter = OpenAIInterpreter(api_key=None)
        results = _make_results(n_terms=5)

        with patch.dict("os.environ", {"OPENAI_API_KEY": ""}):
            interpreter.api_key = None
            interpretations = interpreter.interpret(results)
            assert interpretations == {}

    def test_summarize_term_with_api_key(self):
        """summarize_term() 有 API key 时正常调用。"""
        from allenricher.ai.interpreter import OpenAIInterpreter

        interpreter = OpenAIInterpreter(api_key="test-key")

        mock_response = MagicMock()
        mock_response.choices = [MagicMock()]
        mock_response.choices[0].message.content = "Cell division summary."

        with patch.dict("sys.modules", {"openai": MagicMock()}):
            import openai as mock_openai
            mock_client = MagicMock()
            mock_client.chat.completions.create.return_value = mock_response
            mock_openai.OpenAI.return_value = mock_client

            result = interpreter.summarize_term("cell division", ["GENE1", "GENE2"])
            assert result == "Cell division summary."

    def test_summarize_term_no_api_key(self):
        """summarize_term() 无 API key 时返回空字符串。"""
        from allenricher.ai.interpreter import OpenAIInterpreter

        interpreter = OpenAIInterpreter(api_key=None)
        interpreter.api_key = None

        result = interpreter.summarize_term("cell division", ["GENE1"])
        assert result == ""

    def test_summarize_term_truncates_long_gene_list(self):
        """summarize_term() 超过 10 个基因时截断。"""
        from allenricher.ai.interpreter import OpenAIInterpreter

        interpreter = OpenAIInterpreter(api_key="test-key")
        genes = [f"GENE{i}" for i in range(15)]

        mock_response = MagicMock()
        mock_response.choices = [MagicMock()]
        mock_response.choices[0].message.content = "Summary."

        with patch.dict("sys.modules", {"openai": MagicMock()}):
            import openai as mock_openai
            mock_client = MagicMock()
            mock_client.chat.completions.create.return_value = mock_response
            mock_openai.OpenAI.return_value = mock_client

            interpreter.summarize_term("test_term", genes)

            call_args = mock_client.chat.completions.create.call_args
            prompt = call_args[1]["messages"][1]["content"]
            # 前 10 个基因应出现，第 11 个不应出现
            assert "GENE9" in prompt
            assert "GENE10" not in prompt
            assert "..." in prompt


# ===========================================================================
# 3. ClaudeInterpreter 测试
# ===========================================================================

class TestClaudeInterpreter:
    """测试 ClaudeInterpreter（mock anthropic 包）。"""

    def test_call_api_success(self):
        """_call_api() 正常调用返回结果。"""
        from allenricher.ai.interpreter import ClaudeInterpreter

        interpreter = ClaudeInterpreter(api_key="test-key")

        mock_content = MagicMock()
        mock_content.text = "Claude interpretation result."

        mock_message = MagicMock()
        mock_message.content = [mock_content]

        with patch.dict("sys.modules", {"anthropic": MagicMock()}):
            import anthropic as mock_anthropic
            mock_client = MagicMock()
            mock_client.messages.create.return_value = mock_message
            mock_anthropic.Anthropic.return_value = mock_client

            result = interpreter._call_api("test prompt")
            assert result == "Claude interpretation result."
            mock_client.messages.create.assert_called_once()

    def test_call_api_import_error(self):
        """_call_api() ImportError 降级处理。"""
        from allenricher.ai.interpreter import ClaudeInterpreter

        interpreter = ClaudeInterpreter(api_key="test-key")

        with patch("builtins.__import__", side_effect=ImportError("No module named anthropic")):
            result = interpreter._call_api("test prompt")
            assert "Error" in result
            assert "not installed" in result

    def test_call_api_exception(self):
        """_call_api() 其他异常处理。"""
        from allenricher.ai.interpreter import ClaudeInterpreter

        interpreter = ClaudeInterpreter(api_key="test-key")

        with patch.dict("sys.modules", {"anthropic": MagicMock()}):
            import anthropic as mock_anthropic
            mock_anthropic.Anthropic.side_effect = Exception("Claude API error")

            result = interpreter._call_api("test prompt")
            assert "Error" in result
            assert "Claude API error" in result

    def test_interpret_no_api_key(self):
        """interpret() 无 API key 时返回空字典。"""
        from allenricher.ai.interpreter import ClaudeInterpreter

        interpreter = ClaudeInterpreter(api_key=None)
        interpreter.api_key = None
        results = _make_results(n_terms=3)

        interpretations = interpreter.interpret(results)
        assert interpretations == {}

    def test_interpret_empty_results(self):
        """interpret() 空结果不调用 API。"""
        from allenricher.ai.interpreter import ClaudeInterpreter

        interpreter = ClaudeInterpreter(api_key="test-key")
        results = _make_empty_results()

        with patch.object(interpreter, "_call_api") as mock_call:
            interpretations = interpreter.interpret(results)

            assert "GO_Biological_Process" in interpretations
            assert "No significant enrichment" in interpretations["GO_Biological_Process"]
            mock_call.assert_not_called()

    def test_summarize_term_no_api_key(self):
        """summarize_term() 无 API key 时返回空字符串。"""
        from allenricher.ai.interpreter import ClaudeInterpreter

        interpreter = ClaudeInterpreter(api_key=None)
        interpreter.api_key = None

        result = interpreter.summarize_term("apoptosis", ["GENE1"])
        assert result == ""


# ===========================================================================
# 4. OllamaInterpreter 测试
# ===========================================================================

class TestOllamaInterpreter:
    """测试 OllamaInterpreter（mock requests.post）。"""

    def test_call_api_success(self):
        """_call_api() 正常调用返回结果。"""
        from allenricher.ai.interpreter import OllamaInterpreter

        interpreter = OllamaInterpreter(model="llama2")

        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {"response": "Ollama generated text."}

        with patch("requests.post", return_value=mock_response) as mock_post:
            result = interpreter._call_api("test prompt")
            assert result == "Ollama generated text."
            mock_post.assert_called_once()
            # 验证请求 URL 和参数
            call_args = mock_post.call_args
            assert "localhost:11434/api/generate" in call_args[0][0]
            assert call_args[1]["json"]["model"] == "llama2"
            assert call_args[1]["json"]["stream"] is False

    def test_call_api_connection_failure(self):
        """_call_api() 连接失败处理。"""
        from allenricher.ai.interpreter import OllamaInterpreter

        interpreter = OllamaInterpreter()

        with patch("requests.post", side_effect=Exception("Connection refused")):
            result = interpreter._call_api("test prompt")
            assert "Error" in result
            assert "Connection refused" in result

    def test_call_api_http_error(self):
        """_call_api() HTTP 错误状态码处理。"""
        from allenricher.ai.interpreter import OllamaInterpreter

        interpreter = OllamaInterpreter()

        mock_response = MagicMock()
        mock_response.status_code = 500

        with patch("requests.post", return_value=mock_response):
            result = interpreter._call_api("test prompt")
            assert "Error" in result
            assert "500" in result

    def test_call_api_import_error(self):
        """_call_api() requests 未安装时降级处理。"""
        from allenricher.ai.interpreter import OllamaInterpreter

        interpreter = OllamaInterpreter()

        with patch("builtins.__import__", side_effect=ImportError("No module named requests")):
            result = interpreter._call_api("test prompt")
            assert "Error" in result
            assert "not installed" in result

    def test_interpret_empty_results(self):
        """interpret() 空结果不调用 API。"""
        from allenricher.ai.interpreter import OllamaInterpreter

        interpreter = OllamaInterpreter()
        results = _make_empty_results()

        with patch.object(interpreter, "_call_api") as mock_call:
            interpretations = interpreter.interpret(results)

            assert "GO_Biological_Process" in interpretations
            assert "No significant enrichment" in interpretations["GO_Biological_Process"]
            mock_call.assert_not_called()

    def test_interpret_normal(self):
        """interpret() 正常工作。"""
        from allenricher.ai.interpreter import OllamaInterpreter

        interpreter = OllamaInterpreter()
        results = _make_results(n_terms=3)

        with patch.object(interpreter, "_call_api", return_value="Ollama interpretation."):
            interpretations = interpreter.interpret(results)

            assert "GO_Biological_Process" in interpretations
            assert interpretations["GO_Biological_Process"] == "Ollama interpretation."

    def test_summarize_term(self):
        """summarize_term() 正常工作。"""
        from allenricher.ai.interpreter import OllamaInterpreter

        interpreter = OllamaInterpreter()

        with patch.object(interpreter, "_call_api", return_value="Term description."):
            result = interpreter.summarize_term("apoptosis", ["GENE1"])
            assert result == "Term description."


# ===========================================================================
# 5. DeepSeek / GLM / MiniMaxInterpreter 测试
# ===========================================================================

class TestDeepSeekInterpreter:
    """测试 DeepSeekInterpreter。"""

    def test_init_no_api_key(self):
        """初始化时无 API key 的降级行为。"""
        from allenricher.ai.interpreter import DeepSeekInterpreter

        interpreter = DeepSeekInterpreter(api_key=None)
        assert interpreter.api_key is None
        assert interpreter.model == "deepseek-chat"

    def test_interpret_no_api_key(self):
        """无 API key 时 interpret() 返回空字典。"""
        from allenricher.ai.interpreter import DeepSeekInterpreter

        interpreter = DeepSeekInterpreter(api_key=None)
        interpreter.api_key = None
        results = _make_results(n_terms=3)

        interpretations = interpreter.interpret(results)
        assert interpretations == {}

    def test_summarize_term_no_api_key(self):
        """无 API key 时 summarize_term() 返回空字符串。"""
        from allenricher.ai.interpreter import DeepSeekInterpreter

        interpreter = DeepSeekInterpreter(api_key=None)
        interpreter.api_key = None

        result = interpreter.summarize_term("test", ["GENE1"])
        assert result == ""

    def test_call_api_success(self):
        """_call_api() 正常调用（使用 openai SDK 兼容格式）。"""
        from allenricher.ai.interpreter import DeepSeekInterpreter

        interpreter = DeepSeekInterpreter(api_key="test-key")

        mock_response = MagicMock()
        mock_response.choices = [MagicMock()]
        mock_response.choices[0].message.content = "DeepSeek result."

        with patch.dict("sys.modules", {"openai": MagicMock()}):
            import openai as mock_openai
            mock_client = MagicMock()
            mock_client.chat.completions.create.return_value = mock_response
            mock_openai.OpenAI.return_value = mock_client

            result = interpreter._call_api("test prompt")
            assert result == "DeepSeek result."
            # 验证使用了 DeepSeek 的 base_url
            mock_openai.OpenAI.assert_called_once_with(
                api_key="test-key",
                base_url="https://api.deepseek.com"
            )

    def test_call_api_import_error(self):
        """_call_api() ImportError 降级。"""
        from allenricher.ai.interpreter import DeepSeekInterpreter

        interpreter = DeepSeekInterpreter(api_key="test-key")

        with patch("builtins.__import__", side_effect=ImportError("No module named openai")):
            result = interpreter._call_api("test prompt")
            assert "Error" in result
            assert "not installed" in result


class TestGLMInterpreter:
    """测试 GLMInterpreter。"""

    def test_init_no_api_key(self):
        """初始化时无 API key 的降级行为。"""
        from allenricher.ai.interpreter import GLMInterpreter

        interpreter = GLMInterpreter(api_key=None)
        assert interpreter.api_key is None
        assert interpreter.model == "glm-4"

    def test_interpret_no_api_key(self):
        """无 API key 时 interpret() 返回空字典。"""
        from allenricher.ai.interpreter import GLMInterpreter

        interpreter = GLMInterpreter(api_key=None)
        interpreter.api_key = None
        results = _make_results(n_terms=3)

        interpretations = interpreter.interpret(results)
        assert interpretations == {}

    def test_summarize_term_no_api_key(self):
        """无 API key 时 summarize_term() 返回空字符串。"""
        from allenricher.ai.interpreter import GLMInterpreter

        interpreter = GLMInterpreter(api_key=None)
        interpreter.api_key = None

        result = interpreter.summarize_term("test", ["GENE1"])
        assert result == ""

    def test_call_api_success(self):
        """_call_api() 正常调用。"""
        from allenricher.ai.interpreter import GLMInterpreter

        interpreter = GLMInterpreter(api_key="test-key")

        mock_response = MagicMock()
        mock_response.choices = [MagicMock()]
        mock_response.choices[0].message.content = "GLM result."

        with patch.dict("sys.modules", {"openai": MagicMock()}):
            import openai as mock_openai
            mock_client = MagicMock()
            mock_client.chat.completions.create.return_value = mock_response
            mock_openai.OpenAI.return_value = mock_client

            result = interpreter._call_api("test prompt")
            assert result == "GLM result."
            # 验证使用了 GLM 的 base_url
            mock_openai.OpenAI.assert_called_once_with(
                api_key="test-key",
                base_url="https://open.bigmodel.cn/api/paas/v4"
            )


class TestMiniMaxInterpreter:
    """测试 MiniMaxInterpreter。"""

    def test_init_no_api_key(self):
        """初始化时无 API key 的降级行为。"""
        from allenricher.ai.interpreter import MiniMaxInterpreter

        interpreter = MiniMaxInterpreter(api_key=None, group_id=None)
        assert interpreter.api_key is None
        assert interpreter.group_id is None
        assert interpreter.model == "abab6.5s-chat"

    def test_interpret_no_api_key(self):
        """无 API key 或 group_id 时 interpret() 返回空字典。"""
        from allenricher.ai.interpreter import MiniMaxInterpreter

        interpreter = MiniMaxInterpreter(api_key=None, group_id=None)
        results = _make_results(n_terms=3)

        interpretations = interpreter.interpret(results)
        assert interpretations == {}

    def test_interpret_no_group_id(self):
        """有 API key 但无 group_id 时 interpret() 返回空字典。"""
        from allenricher.ai.interpreter import MiniMaxInterpreter

        interpreter = MiniMaxInterpreter(api_key="test-key", group_id=None)
        results = _make_results(n_terms=3)

        interpretations = interpreter.interpret(results)
        assert interpretations == {}

    def test_summarize_term_no_api_key(self):
        """无 API key 时 summarize_term() 返回空字符串。"""
        from allenricher.ai.interpreter import MiniMaxInterpreter

        interpreter = MiniMaxInterpreter(api_key=None, group_id=None)

        result = interpreter.summarize_term("test", ["GENE1"])
        assert result == ""

    def test_call_api_no_credentials(self):
        """_call_api() 无凭证时返回错误。"""
        from allenricher.ai.interpreter import MiniMaxInterpreter

        interpreter = MiniMaxInterpreter(api_key=None, group_id=None)

        # 需要 mock openai 包使 import 成功，但 api_key/group_id 检查在调用前
        with patch.dict("sys.modules", {"openai": MagicMock()}):
            result = interpreter._call_api("test prompt")
            assert "Error" in result
            assert "not configured" in result

    def test_call_api_success(self):
        """_call_api() 正常调用。"""
        from allenricher.ai.interpreter import MiniMaxInterpreter

        interpreter = MiniMaxInterpreter(api_key="test-key", group_id="test-group")

        mock_response = MagicMock()
        mock_response.choices = [MagicMock()]
        mock_response.choices[0].message.content = "MiniMax result."

        with patch.dict("sys.modules", {"openai": MagicMock()}):
            import openai as mock_openai
            mock_client = MagicMock()
            mock_client.chat.completions.create.return_value = mock_response
            mock_openai.OpenAI.return_value = mock_client

            result = interpreter._call_api("test prompt")
            assert result == "MiniMax result."


# ===========================================================================
# 6. AIInterpreter 门面类测试
# ===========================================================================

class TestAIInterpreterFacade:
    """测试 AIInterpreter 门面类。"""

    def test_interpret_results_calls_backend(self):
        """interpret_results() 调用后端并返回结果。"""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=3)

        interpretations = interpreter.interpret_results(results)

        assert "GO_Biological_Process" in interpretations
        assert "enrichment" in interpretations["GO_Biological_Process"].lower()

    def test_interpret_results_empty(self):
        """interpret_results() 空结果处理。"""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_empty_results()

        interpretations = interpreter.interpret_results(results)

        # MockInterpreter 对空结果不放入字典
        assert len(interpretations) == 0

    def test_interpret_results_with_context(self):
        """interpret_results() 传入 context 参数。"""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=2)

        # MockInterpreter 不使用 context，但不应抛出异常
        interpretations = interpreter.interpret_results(results, context="cancer research")
        assert "GO_Biological_Process" in interpretations

    def test_interpret_results_include_term_summaries(self):
        """include_term_summaries=True 时生成条目总结。"""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=3, include_genes_col=True)

        interpretations = interpreter.interpret_results(
            results, include_term_summaries=True
        )

        # 应包含整体解读
        assert "GO_Biological_Process" in interpretations
        # 应包含条目总结
        assert "GO_Biological_Process_term_summaries" in interpretations
        term_summaries = interpretations["GO_Biological_Process_term_summaries"]
        assert isinstance(term_summaries, dict)
        assert len(term_summaries) == 3
        # 验证每个条目总结包含基因数量信息
        for term_name, summary in term_summaries.items():
            assert "genes" in summary.lower()

    def test_interpret_results_include_term_summaries_empty(self):
        """include_term_summaries=True 但结果为空时不生成条目总结。"""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_empty_results()

        interpretations = interpreter.interpret_results(
            results, include_term_summaries=True
        )

        assert "GO_Biological_Process_term_summaries" not in interpretations

    def test_generate_report_section(self):
        """generate_report_section() 生成 HTML。"""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=2)

        html = interpreter.generate_report_section(results)

        assert "<div" in html
        assert "</div>" in html
        assert "AI-Powered Interpretation" in html
        assert "mock" in html
        assert "GO_Biological_Process" in html
        assert "ai-disclaimer" in html

    def test_generate_report_section_empty_results(self):
        """generate_report_section() 空结果时生成 HTML（无解读内容）。"""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_empty_results()

        html = interpreter.generate_report_section(results)

        assert "<div" in html
        assert "AI-Powered Interpretation" in html
        # 不应包含任何数据库解读块
        assert "GO_Biological_Process" not in html

    def test_generate_report_section_with_context(self):
        """generate_report_section() 传入 context 参数。"""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=2)

        html = interpreter.generate_report_section(results, context="diabetes study")
        assert "AI-Powered Interpretation" in html

    def test_generate_report_section_excludes_term_summaries(self):
        """generate_report_section() 不包含条目总结（仅整体解读）。"""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=2, include_genes_col=True)

        html = interpreter.generate_report_section(results)

        # HTML 中不应出现 _term_summaries 相关内容
        assert "_term_summaries" not in html

    def test_backend_name_stored(self):
        """backend_name 正确存储。"""
        from allenricher.ai.interpreter import AIInterpreter

        for backend in ["mock", "openai", "claude", "ollama", "deepseek", "glm", "minimax"]:
            kwargs = {}
            if backend == "minimax":
                kwargs["group_id"] = "test"
            interpreter = AIInterpreter(backend=backend, **kwargs)
            assert interpreter.backend_name == backend


# ===========================================================================
# 7. 工厂函数测试
# ===========================================================================

class TestFactoryFunctions:
    """测试 create_interpreter() 和 get_available_backends() 工厂函数。"""

    def test_create_interpreter_mock(self):
        """create_interpreter("mock") 返回正确实例。"""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter("mock")
        assert interpreter.backend_name == "mock"

    def test_create_interpreter_openai(self):
        """create_interpreter("openai") 返回正确实例。"""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter("openai", api_key="test-key")
        assert interpreter.backend_name == "openai"
        assert interpreter.interpreter.api_key == "test-key"

    def test_create_interpreter_claude(self):
        """create_interpreter("claude") 返回正确实例。"""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter("claude", api_key="test-key")
        assert interpreter.backend_name == "claude"
        assert interpreter.interpreter.api_key == "test-key"

    def test_create_interpreter_ollama(self):
        """create_interpreter("ollama") 返回正确实例。"""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter("ollama", model="mistral")
        assert interpreter.backend_name == "ollama"
        assert interpreter.interpreter.model == "mistral"

    def test_create_interpreter_deepseek(self):
        """create_interpreter("deepseek") 返回正确实例。"""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter("deepseek", api_key="test-key")
        assert interpreter.backend_name == "deepseek"
        assert interpreter.interpreter.api_key == "test-key"

    def test_create_interpreter_glm(self):
        """create_interpreter("glm") 返回正确实例。"""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter("glm", api_key="test-key")
        assert interpreter.backend_name == "glm"
        assert interpreter.interpreter.api_key == "test-key"

    def test_create_interpreter_minimax(self):
        """create_interpreter("minimax") 返回正确实例。"""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter(
            "minimax", api_key="test-key", group_id="test-group"
        )
        assert interpreter.backend_name == "minimax"
        assert interpreter.interpreter.api_key == "test-key"
        assert interpreter.interpreter.group_id == "test-group"

    def test_create_interpreter_invalid(self):
        """create_interpreter("invalid") 抛出 ValueError。"""
        from allenricher.ai.interpreter import create_interpreter

        with pytest.raises(ValueError, match="Unknown backend"):
            create_interpreter("invalid")

    def test_create_interpreter_default_backend(self):
        """create_interpreter() 默认使用 mock 后端。"""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter()
        assert interpreter.backend_name == "mock"

    def test_get_available_backends(self):
        """get_available_backends() 返回列表。"""
        from allenricher.ai.interpreter import get_available_backends

        backends = get_available_backends()
        assert isinstance(backends, list)
        assert len(backends) == 7
        expected = {"openai", "claude", "deepseek", "glm", "minimax", "ollama", "mock"}
        assert set(backends) == expected


# ===========================================================================
# 8. 抽象基类接口一致性测试
# ===========================================================================

class TestAbstractBaseInterface:
    """验证所有后端实现都遵循统一的抽象接口。"""

    @pytest.fixture(params=["mock", "ollama"])
    def concrete_interpreter(self, request):
        """创建不需要 API key 的具体解释器实例。"""
        from allenricher.ai.interpreter import AIInterpreter

        return AIInterpreter(backend=request.param)

    def test_interpret_signature(self, concrete_interpreter):
        """interpret() 方法接受 Dict[str, DataFrame] 参数。"""
        results = _make_results(n_terms=2)
        # 不应抛出异常
        result = concrete_interpreter.interpreter.interpret(results)
        assert isinstance(result, dict)

    def test_summarize_term_signature(self, concrete_interpreter):
        """summarize_term() 方法接受 str 和 List[str] 参数。"""
        result = concrete_interpreter.interpreter.summarize_term("test_term", ["GENE1"])
        assert isinstance(result, str)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
