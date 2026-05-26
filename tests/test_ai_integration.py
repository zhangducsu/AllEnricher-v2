"""
AI 集成测试

端到端测试 AI 解读模块与富集分析管线的集成，覆盖：
- MockInterpreter 完整流程：analyze → AI 解读 → JSON 文件 → HTML 报告
- 0 条显著结果：直接报告"没有显著富集结果"，不调用 AI
- 不足 20 条：按实际条目数发送给 AI
- 恰好 20 条：发送全部 20 条
- 超过 20 条：仅发送前 20 条
- AI JSON 输出格式：ai_interpretation.json 结构验证
- HTML 报告嵌入：AI 解读段落正确嵌入 HTML

使用 MockInterpreter 作为 AI 后端，无需真实 API 密钥或网络连接。
"""

import json
import pytest
import sys
from pathlib import Path
from unittest.mock import patch, MagicMock
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))


# ---------------------------------------------------------------------------
# 测试数据工厂
# ---------------------------------------------------------------------------

def _make_results(n_terms=5, include_genes_col=False):
    """创建模拟富集分析结果。

    参数:
        n_terms: 富集条目数量
        include_genes_col: 是否包含 'Genes' 列
    """
    data = {
        "Term_ID": [f"GO:{1000 + i}" for i in range(n_terms)],
        "Term_Name": [f"biological_process_{i}" for i in range(n_terms)],
        "Gene_Count": [10 - i for i in range(n_terms)],
        "P_Value": [1e-5 * (i + 1) for i in range(n_terms)],
        "Adjusted_P_Value": [1e-3 * (i + 1) for i in range(n_terms)],
    }
    if include_genes_col:
        data["Genes"] = [
            ";".join([f"GENE{j}" for j in range(5 - min(i, 4))])
            for i in range(n_terms)
        ]
    df = pd.DataFrame(data)
    return {"GO_Biological_Process": df}


def _make_empty_results():
    """创建空结果。"""
    return {"GO_Biological_Process": pd.DataFrame()}


def _make_multi_db_results(n_terms=5):
    """创建多数据库结果。"""
    data = {
        "Term_ID": [f"GO:{1000 + i}" for i in range(n_terms)],
        "Term_Name": [f"biological_process_{i}" for i in range(n_terms)],
        "Gene_Count": [10 - i for i in range(n_terms)],
        "P_Value": [1e-5 * (i + 1) for i in range(n_terms)],
        "Adjusted_P_Value": [1e-3 * (i + 1) for i in range(n_terms)],
    }
    return {
        "GO_Biological_Process": pd.DataFrame(data),
        "KEGG": pd.DataFrame({
            "Term_ID": [f"path:map{i:03d}" for i in range(n_terms)],
            "Term_Name": [f"kegg_pathway_{i}" for i in range(n_terms)],
            "Gene_Count": [8 - i for i in range(n_terms)],
            "P_Value": [1e-4 * (i + 1) for i in range(n_terms)],
            "Adjusted_P_Value": [1e-2 * (i + 1) for i in range(n_terms)],
        }),
    }


# ===========================================================================
# 1. MockInterpreter 端到端流程测试
# ===========================================================================

class TestMockInterpreterEndToEnd:
    """MockInterpreter 完整流程：analyze → AI 解读 → JSON 文件 → HTML 报告。"""

    def test_e2e_interpret_to_json(self, tmp_path):
        """端到端：interpret_results → 保存为 JSON → 验证文件结构。"""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter(backend="mock")
        results = _make_results(n_terms=10)

        # 1. 生成 AI 解读
        ai_interpretation = interpreter.interpret_results(results)

        # 2. 保存为 JSON 文件
        json_path = tmp_path / "ai_interpretation.json"
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(ai_interpretation, f, indent=2)

        # 3. 验证文件存在且可读
        assert json_path.exists()

        # 4. 验证 JSON 结构
        with open(json_path, 'r', encoding='utf-8') as f:
            loaded = json.load(f)

        assert isinstance(loaded, dict)
        assert "GO_Biological_Process" in loaded
        assert isinstance(loaded["GO_Biological_Process"], str)
        assert len(loaded["GO_Biological_Process"]) > 0

    def test_e2e_interpret_to_html(self, tmp_path):
        """端到端：interpret_results → HTML 报告嵌入 → 验证 HTML 结构。"""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter(backend="mock")
        results = _make_results(n_terms=5)

        # 1. 生成 AI 解读
        ai_interpretation = interpreter.interpret_results(results)

        # 2. 生成 HTML 报告段落
        html = interpreter.generate_report_section(results)

        # 3. 验证 HTML 结构
        assert "<div" in html
        assert "</div>" in html
        assert "AI-Powered Interpretation" in html
        assert "GO_Biological_Process" in html
        assert "ai-disclaimer" in html
        assert "mock" in html  # 后端名称

        # 4. 验证 AI 解读内容嵌入 HTML
        for i in range(5):
            assert f"biological_process_{i}" in html

    def test_e2e_multi_database(self, tmp_path):
        """端到端：多数据库同时解读 → JSON + HTML。"""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter(backend="mock")
        results = _make_multi_db_results(n_terms=5)

        # 生成 AI 解读
        ai_interpretation = interpreter.interpret_results(results)

        # 验证两个数据库都有解读
        assert "GO_Biological_Process" in ai_interpretation
        assert "KEGG" in ai_interpretation

        # 保存 JSON
        json_path = tmp_path / "ai_interpretation.json"
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(ai_interpretation, f, indent=2)

        with open(json_path, 'r', encoding='utf-8') as f:
            loaded = json.load(f)

        assert len(loaded) == 2

        # 验证 HTML 包含两个数据库
        html = interpreter.generate_report_section(results)
        assert "GO_Biological_Process" in html
        assert "KEGG" in html

    def test_e2e_with_term_summaries(self, tmp_path):
        """端到端：包含 term_summaries 的完整流程。"""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter(backend="mock")
        results = _make_results(n_terms=5, include_genes_col=True)

        # 生成包含 term_summaries 的 AI 解读
        ai_interpretation = interpreter.interpret_results(
            results, include_term_summaries=True
        )

        # 验证包含整体解读和条目总结
        assert "GO_Biological_Process" in ai_interpretation
        assert "GO_Biological_Process_term_summaries" in ai_interpretation

        term_summaries = ai_interpretation["GO_Biological_Process_term_summaries"]
        assert isinstance(term_summaries, dict)
        assert len(term_summaries) == 5

        # 保存 JSON 并验证结构
        json_path = tmp_path / "ai_interpretation.json"
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(ai_interpretation, f, indent=2)

        with open(json_path, 'r', encoding='utf-8') as f:
            loaded = json.load(f)

        # term_summaries 应为嵌套字典
        assert isinstance(loaded["GO_Biological_Process_term_summaries"], dict)
        for term_name, summary in loaded["GO_Biological_Process_term_summaries"].items():
            assert isinstance(term_name, str)
            assert isinstance(summary, str)
            assert len(summary) > 0


# ===========================================================================
# 2. 0 条显著结果测试
# ===========================================================================

class TestZeroResults:
    """0 条显著结果场景：直接报告"没有显著富集结果"，不调用 AI。"""

    def test_mock_zero_results_skips_interpretation(self):
        """MockInterpreter：0 条结果时不生成解读。"""
        from allenricher.ai.interpreter import MockInterpreter

        interpreter = MockInterpreter()
        results = _make_empty_results()

        interpretations = interpreter.interpret(results)

        # MockInterpreter 对空结果直接 continue，不放入字典
        assert len(interpretations) == 0

    def test_openai_zero_results_no_api_call(self):
        """OpenAI：0 条结果时不调用 API，返回固定消息。"""
        from allenricher.ai.interpreter import OpenAIInterpreter

        interpreter = OpenAIInterpreter(api_key="test-key")
        results = _make_empty_results()

        with patch.object(interpreter, "_call_api") as mock_call:
            interpretations = interpreter.interpret(results)

            # 应返回"没有显著富集结果"的消息
            assert "GO_Biological_Process" in interpretations
            assert "No significant enrichment" in interpretations["GO_Biological_Process"]
            # API 不应被调用
            mock_call.assert_not_called()

    def test_claude_zero_results_no_api_call(self):
        """Claude：0 条结果时不调用 API。"""
        from allenricher.ai.interpreter import ClaudeInterpreter

        interpreter = ClaudeInterpreter(api_key="test-key")
        results = _make_empty_results()

        with patch.object(interpreter, "_call_api") as mock_call:
            interpretations = interpreter.interpret(results)

            assert "GO_Biological_Process" in interpretations
            assert "No significant enrichment" in interpretations["GO_Biological_Process"]
            mock_call.assert_not_called()

    def test_ollama_zero_results_no_api_call(self):
        """Ollama：0 条结果时不调用 API。"""
        from allenricher.ai.interpreter import OllamaInterpreter

        interpreter = OllamaInterpreter()
        results = _make_empty_results()

        with patch.object(interpreter, "_call_api") as mock_call:
            interpretations = interpreter.interpret(results)

            assert "GO_Biological_Process" in interpretations
            assert "No significant enrichment" in interpretations["GO_Biological_Process"]
            mock_call.assert_not_called()

    def test_facade_zero_results(self):
        """门面类：0 条结果时返回空字典（MockInterpreter 行为）。"""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_empty_results()

        interpretations = interpreter.interpret_results(results)
        assert len(interpretations) == 0

    def test_facade_zero_results_no_term_summaries(self):
        """门面类：0 条结果时不生成 term_summaries。"""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_empty_results()

        interpretations = interpreter.interpret_results(
            results, include_term_summaries=True
        )

        assert len(interpretations) == 0
        assert "GO_Biological_Process_term_summaries" not in interpretations

    def test_zero_results_json_output(self, tmp_path):
        """0 条结果时 JSON 文件结构验证。"""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter(backend="mock")
        results = _make_empty_results()

        ai_interpretation = interpreter.interpret_results(results)

        json_path = tmp_path / "ai_interpretation.json"
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(ai_interpretation, f, indent=2)

        with open(json_path, 'r', encoding='utf-8') as f:
            loaded = json.load(f)

        # MockInterpreter 返回空字典
        assert loaded == {}

    def test_zero_results_html_empty_section(self):
        """0 条结果时 HTML 报告不包含解读内容。"""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_empty_results()

        html = interpreter.generate_report_section(results)

        # 应有 AI 标题但无数据库解读块
        assert "AI-Powered Interpretation" in html
        assert "GO_Biological_Process" not in html


# ===========================================================================
# 3. 不足 20 条结果测试
# ===========================================================================

class TestFewerThan20Results:
    """不足 20 条：按实际条目数发送给 AI。"""

    def test_5_results_sends_5(self):
        """5 条结果时 prompt 中包含 5 条。"""
        from allenricher.ai.interpreter import OpenAIInterpreter

        interpreter = OpenAIInterpreter(api_key="test-key")
        results = _make_results(n_terms=5)

        mock_response = MagicMock()
        mock_response.choices = [MagicMock()]
        mock_response.choices[0].message.content = "Interpretation."

        with patch.dict("sys.modules", {"openai": MagicMock()}):
            import openai as mock_openai
            mock_client = MagicMock()
            mock_client.chat.completions.create.return_value = mock_response
            mock_openai.OpenAI.return_value = mock_client

            interpreter.interpret(results)

            call_args = mock_client.chat.completions.create.call_args
            prompt = call_args[1]["messages"][1]["content"]

            # 验证 prompt 中说 "Top 5"（不是 Top 20）
            assert "Top 5" in prompt
            # 验证包含全部 5 条
            for i in range(5):
                assert f"biological_process_{i}" in prompt
            # 不应包含第 6 条
            assert "biological_process_5" not in prompt

    def test_15_results_sends_15(self):
        """15 条结果时 prompt 中包含 15 条。"""
        from allenricher.ai.interpreter import OpenAIInterpreter

        interpreter = OpenAIInterpreter(api_key="test-key")
        results = _make_results(n_terms=15)

        mock_response = MagicMock()
        mock_response.choices = [MagicMock()]
        mock_response.choices[0].message.content = "Interpretation."

        with patch.dict("sys.modules", {"openai": MagicMock()}):
            import openai as mock_openai
            mock_client = MagicMock()
            mock_client.chat.completions.create.return_value = mock_response
            mock_openai.OpenAI.return_value = mock_client

            interpreter.interpret(results)

            call_args = mock_client.chat.completions.create.call_args
            prompt = call_args[1]["messages"][1]["content"]

            assert "Top 15" in prompt
            for i in range(15):
                assert f"biological_process_{i}" in prompt
            assert "biological_process_15" not in prompt

    def test_1_result_sends_1(self):
        """1 条结果时 prompt 中包含 1 条。"""
        from allenricher.ai.interpreter import OpenAIInterpreter

        interpreter = OpenAIInterpreter(api_key="test-key")
        results = _make_results(n_terms=1)

        mock_response = MagicMock()
        mock_response.choices = [MagicMock()]
        mock_response.choices[0].message.content = "Interpretation."

        with patch.dict("sys.modules", {"openai": MagicMock()}):
            import openai as mock_openai
            mock_client = MagicMock()
            mock_client.chat.completions.create.return_value = mock_response
            mock_openai.OpenAI.return_value = mock_client

            interpreter.interpret(results)

            call_args = mock_client.chat.completions.create.call_args
            prompt = call_args[1]["messages"][1]["content"]

            assert "Top 1" in prompt
            assert "biological_process_0" in prompt

    def test_mock_5_results_shows_all(self):
        """MockInterpreter：5 条结果时展示全部 5 条。"""
        from allenricher.ai.interpreter import MockInterpreter

        interpreter = MockInterpreter()
        results = _make_results(n_terms=5)

        interpretations = interpreter.interpret(results)
        text = interpretations["GO_Biological_Process"]

        for i in range(5):
            assert f"biological_process_{i}" in text

    def test_facade_fewer_than_20_term_summaries(self):
        """门面类：不足 20 条时 term_summaries 按实际数量生成。"""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=8, include_genes_col=True)

        interpretations = interpreter.interpret_results(
            results, include_term_summaries=True
        )

        term_summaries = interpretations["GO_Biological_Process_term_summaries"]
        assert len(term_summaries) == 8


# ===========================================================================
# 4. 恰好 20 条结果测试
# ===========================================================================

class TestExactly20Results:
    """恰好 20 条：发送全部 20 条。"""

    def test_openai_20_results_sends_all(self):
        """OpenAI：恰好 20 条时 prompt 包含全部 20 条。"""
        from allenricher.ai.interpreter import OpenAIInterpreter

        interpreter = OpenAIInterpreter(api_key="test-key")
        results = _make_results(n_terms=20)

        mock_response = MagicMock()
        mock_response.choices = [MagicMock()]
        mock_response.choices[0].message.content = "Interpretation."

        with patch.dict("sys.modules", {"openai": MagicMock()}):
            import openai as mock_openai
            mock_client = MagicMock()
            mock_client.chat.completions.create.return_value = mock_response
            mock_openai.OpenAI.return_value = mock_client

            interpreter.interpret(results)

            call_args = mock_client.chat.completions.create.call_args
            prompt = call_args[1]["messages"][1]["content"]

            assert "Top 20" in prompt
            for i in range(20):
                assert f"biological_process_{i}" in prompt

    def test_mock_20_results_shows_all(self):
        """MockInterpreter：恰好 20 条时展示全部 20 条。"""
        from allenricher.ai.interpreter import MockInterpreter

        interpreter = MockInterpreter()
        results = _make_results(n_terms=20)

        interpretations = interpreter.interpret(results)
        text = interpretations["GO_Biological_Process"]

        for i in range(20):
            assert f"biological_process_{i}" in text

    def test_facade_20_term_summaries(self):
        """门面类：恰好 20 条时 term_summaries 生成 20 个。"""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=20, include_genes_col=True)

        interpretations = interpreter.interpret_results(
            results, include_term_summaries=True
        )

        term_summaries = interpretations["GO_Biological_Process_term_summaries"]
        assert len(term_summaries) == 20


# ===========================================================================
# 5. 超过 20 条结果测试
# ===========================================================================

class TestMoreThan20Results:
    """超过 20 条：仅发送前 20 条（Top 20 截断）。"""

    def test_openai_25_results_truncates_to_20(self):
        """OpenAI：25 条结果时 prompt 仅包含前 20 条。"""
        from allenricher.ai.interpreter import OpenAIInterpreter

        interpreter = OpenAIInterpreter(api_key="test-key")
        results = _make_results(n_terms=25)

        mock_response = MagicMock()
        mock_response.choices = [MagicMock()]
        mock_response.choices[0].message.content = "Interpretation."

        with patch.dict("sys.modules", {"openai": MagicMock()}):
            import openai as mock_openai
            mock_client = MagicMock()
            mock_client.chat.completions.create.return_value = mock_response
            mock_openai.OpenAI.return_value = mock_client

            interpreter.interpret(results)

            call_args = mock_client.chat.completions.create.call_args
            prompt = call_args[1]["messages"][1]["content"]

            # 应为 "Top 20"（不是 Top 25）
            assert "Top 20" in prompt
            # 前 20 条应出现
            for i in range(20):
                assert f"biological_process_{i}" in prompt
            # 第 21 条不应出现
            assert "biological_process_20" not in prompt

    def test_mock_25_results_truncates_to_20(self):
        """MockInterpreter：25 条结果时仅展示前 20 条。"""
        from allenricher.ai.interpreter import MockInterpreter

        interpreter = MockInterpreter()
        results = _make_results(n_terms=25)

        interpretations = interpreter.interpret(results)
        text = interpretations["GO_Biological_Process"]

        for i in range(20):
            assert f"biological_process_{i}" in text
        assert "biological_process_20" not in text

    def test_facade_25_term_summaries_truncates_to_20(self):
        """门面类：25 条时 term_summaries 仅生成前 20 个。"""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=25, include_genes_col=True)

        interpretations = interpreter.interpret_results(
            results, include_term_summaries=True
        )

        term_summaries = interpretations["GO_Biological_Process_term_summaries"]
        assert len(term_summaries) == 20
        # 前 20 个条目名称应在
        for i in range(20):
            assert f"biological_process_{i}" in term_summaries
        # 第 21 个不应在
        assert "biological_process_20" not in term_summaries


# ===========================================================================
# 6. AI JSON 输出格式验证
# ===========================================================================

class TestAIJsonOutputFormat:
    """ai_interpretation.json 结构验证。"""

    def test_json_structure_single_db(self, tmp_path):
        """单数据库 JSON 结构。"""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter(backend="mock")
        results = _make_results(n_terms=5)

        ai_interpretation = interpreter.interpret_results(results)

        json_path = tmp_path / "ai_interpretation.json"
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(ai_interpretation, f, indent=2)

        with open(json_path, 'r', encoding='utf-8') as f:
            loaded = json.load(f)

        # 顶层是 dict
        assert isinstance(loaded, dict)
        # 键是数据库名
        assert "GO_Biological_Process" in loaded
        # 值是字符串
        assert isinstance(loaded["GO_Biological_Process"], str)
        # 内容不为空
        assert len(loaded["GO_Biological_Process"]) > 50

    def test_json_structure_multi_db(self, tmp_path):
        """多数据库 JSON 结构。"""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter(backend="mock")
        results = _make_multi_db_results(n_terms=5)

        ai_interpretation = interpreter.interpret_results(results)

        json_path = tmp_path / "ai_interpretation.json"
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(ai_interpretation, f, indent=2)

        with open(json_path, 'r', encoding='utf-8') as f:
            loaded = json.load(f)

        assert len(loaded) == 2
        for db_name in ["GO_Biological_Process", "KEGG"]:
            assert db_name in loaded
            assert isinstance(loaded[db_name], str)
            assert len(loaded[db_name]) > 0

    def test_json_with_term_summaries_structure(self, tmp_path):
        """包含 term_summaries 的 JSON 结构。"""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter(backend="mock")
        results = _make_results(n_terms=3, include_genes_col=True)

        ai_interpretation = interpreter.interpret_results(
            results, include_term_summaries=True
        )

        json_path = tmp_path / "ai_interpretation.json"
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(ai_interpretation, f, indent=2)

        with open(json_path, 'r', encoding='utf-8') as f:
            loaded = json.load(f)

        # 整体解读
        assert "GO_Biological_Process" in loaded
        assert isinstance(loaded["GO_Biological_Process"], str)

        # term_summaries 嵌套字典
        assert "GO_Biological_Process_term_summaries" in loaded
        summaries = loaded["GO_Biological_Process_term_summaries"]
        assert isinstance(summaries, dict)
        assert len(summaries) == 3

        # 每个 summary 的值是字符串
        for key, value in summaries.items():
            assert isinstance(key, str)
            assert isinstance(value, str)

    def test_json_is_valid_utf8(self, tmp_path):
        """JSON 文件使用 UTF-8 编码，可正确读写。"""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter(backend="mock")
        results = _make_results(n_terms=5)

        ai_interpretation = interpreter.interpret_results(results)

        json_path = tmp_path / "ai_interpretation.json"
        # 写入 UTF-8
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(ai_interpretation, f, indent=2, ensure_ascii=False)

        # 读取 UTF-8
        with open(json_path, 'r', encoding='utf-8') as f:
            loaded = json.load(f)

        assert loaded == ai_interpretation

    def test_json_empty_results(self, tmp_path):
        """空结果时 JSON 为空字典。"""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter(backend="mock")
        results = _make_empty_results()

        ai_interpretation = interpreter.interpret_results(results)

        json_path = tmp_path / "ai_interpretation.json"
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(ai_interpretation, f, indent=2)

        with open(json_path, 'r', encoding='utf-8') as f:
            loaded = json.load(f)

        assert loaded == {}


# ===========================================================================
# 7. HTML 报告嵌入验证
# ===========================================================================

class TestHTMLReportEmbedding:
    """AI 解读段落正确嵌入 HTML 报告。"""

    def test_html_contains_ai_section(self):
        """HTML 包含 AI 解读段落。"""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=3)

        html = interpreter.generate_report_section(results)

        # 验证 AI 段落标记
        assert 'id="ai-interpretation"' in html
        assert "AI-Powered Interpretation" in html

    def test_html_contains_disclaimer(self):
        """HTML 包含 AI 免责声明。"""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=3)

        html = interpreter.generate_report_section(results)

        assert "ai-disclaimer" in html
        assert "reviewed by domain experts" in html

    def test_html_contains_backend_name(self):
        """HTML 包含后端名称。"""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=3)

        html = interpreter.generate_report_section(results)

        assert "mock" in html

    def test_html_contains_db_interpretation(self):
        """HTML 包含数据库解读内容。"""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=5)

        html = interpreter.generate_report_section(results)

        # 应包含数据库解读块
        assert "GO_Biological_Process" in html
        # 应包含富集条目
        for i in range(5):
            assert f"biological_process_{i}" in html

    def test_html_multi_db_sections(self):
        """HTML 包含多个数据库的解读段落。"""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_multi_db_results(n_terms=3)

        html = interpreter.generate_report_section(results)

        assert "GO_Biological_Process" in html
        assert "KEGG" in html

    def test_html_newlines_converted_to_br(self):
        """HTML 中换行符被转换为 <br>。"""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=3)

        html = interpreter.generate_report_section(results)

        # MockInterpreter 生成的文本包含换行符，应被转换为 <br>
        assert "<br>" in html

    def test_html_empty_results_no_db_block(self):
        """空结果时 HTML 不包含数据库解读块。"""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_empty_results()

        html = interpreter.generate_report_section(results)

        # 应有 AI 标题和免责声明
        assert "AI-Powered Interpretation" in html
        assert "ai-disclaimer" in html
        # 不应有数据库解读块
        assert "GO_Biological_Process" not in html

    def test_html_excludes_term_summaries(self):
        """HTML 不包含 term_summaries（仅展示整体解读）。"""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=3, include_genes_col=True)

        # 即使请求了 term_summaries，HTML 中也不应出现
        html = interpreter.generate_report_section(results)

        assert "_term_summaries" not in html

    def test_html_is_valid_fragment(self):
        """HTML 片段结构完整（有开闭标签）。"""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=3)

        html = interpreter.generate_report_section(results)

        # 验证 <div> 标签配对
        assert html.strip().startswith("<div")
        assert html.strip().endswith("</div>")
        # 验证 <h2> 标签
        assert "<h2>" in html
        assert "</h2>" in html


# ===========================================================================
# 8. 混合场景测试（部分数据库有结果，部分为空）
# ===========================================================================

class TestMixedResults:
    """混合场景：部分数据库有结果，部分为空。"""

    def test_mixed_empty_and_non_empty(self):
        """部分数据库有结果，部分为空。"""
        from allenricher.ai.interpreter import MockInterpreter

        interpreter = MockInterpreter()
        results = {
            "GO_Biological_Process": pd.DataFrame({
                "Term_Name": ["process_1", "process_2"],
                "P_Value": [1e-5, 1e-3],
                "Gene_Count": [10, 5],
            }),
            "KEGG": pd.DataFrame(),  # 空
        }

        interpretations = interpreter.interpret(results)

        # GO 有解读
        assert "GO_Biological_Process" in interpretations
        # KEGG 为空，MockInterpreter 不放入字典
        assert "KEGG" not in interpretations

    def test_mixed_openai_empty_skips_api(self):
        """OpenAI：混合场景中空数据库不调用 API。"""
        from allenricher.ai.interpreter import OpenAIInterpreter

        interpreter = OpenAIInterpreter(api_key="test-key")
        results = {
            "GO_Biological_Process": pd.DataFrame({
                "Term_Name": ["process_1"],
                "P_Value": [1e-5],
                "Gene_Count": [10],
            }),
            "KEGG": pd.DataFrame(),
        }

        mock_response = MagicMock()
        mock_response.choices = [MagicMock()]
        mock_response.choices[0].message.content = "GO interpretation."

        with patch.dict("sys.modules", {"openai": MagicMock()}):
            import openai as mock_openai
            mock_client = MagicMock()
            mock_client.chat.completions.create.return_value = mock_response
            mock_openai.OpenAI.return_value = mock_client

            interpretations = interpreter.interpret(results)

            # GO 应调用 API 一次
            assert mock_client.chat.completions.create.call_count == 1
            # KEGG 应返回固定消息
            assert "KEGG" in interpretations
            assert "No significant enrichment" in interpretations["KEGG"]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
