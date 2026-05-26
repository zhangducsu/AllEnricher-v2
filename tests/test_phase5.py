"""
阶段5模块单元测试：API、AI解读、HTML报告
"""

import pytest
import sys
import tempfile
import os
from pathlib import Path
from unittest.mock import patch, MagicMock
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))


class TestAPIEndpoints:
    """测试 API 端点"""

    def test_app_creation(self):
        """测试 FastAPI 应用创建"""
        from allenricher.api.server import app
        assert app is not None
        assert app.title == "AllEnricher API"

    def test_root_endpoint(self):
        """测试根端点"""
        from allenricher.api.server import app
        from fastapi.testclient import TestClient
        
        client = TestClient(app)
        response = client.get("/")
        assert response.status_code == 200
        # 根路由现在返回 HTML (Web 界面)
        content_type = response.headers.get("content-type", "")
        assert "text/html" in content_type or "application/json" in content_type

    def test_species_endpoint(self):
        """测试物种列表端点"""
        from allenricher.api.server import app
        from fastapi.testclient import TestClient
        
        client = TestClient(app)
        response = client.get("/api/species")
        assert response.status_code == 200
        data = response.json()
        assert isinstance(data, list)
        assert any(s["code"] == "hsa" for s in data)

    def test_databases_endpoint(self):
        """测试数据库列表端点"""
        from allenricher.api.server import app
        from fastapi.testclient import TestClient
        
        client = TestClient(app)
        response = client.get("/api/databases")
        assert response.status_code == 200
        data = response.json()
        assert "databases" in data
        db_names = [db["name"] for db in data["databases"]]
        assert "GO" in db_names
        assert "KEGG" in db_names

    def test_status_not_found(self):
        """测试查询不存在的任务"""
        from allenricher.api.server import app
        from fastapi.testclient import TestClient
        
        client = TestClient(app)
        response = client.get("/api/status/nonexistent-job-id")
        assert response.status_code == 404


class TestAIInterpreter:
    """测试 AI 解读模块"""

    def test_mock_interpreter(self):
        """测试 Mock 解释器"""
        from allenricher.ai.interpreter import MockInterpreter
        
        interpreter = MockInterpreter()
        
        # 创建模拟结果
        results = {
            "GO": pd.DataFrame({
                "Term_ID": ["GO:0005576", "GO:0051301"],
                "Term_Name": ["extracellular region", "cell division"],
                "Gene_Count": [10, 8],
                "P_Value": [1e-5, 1e-4],
                "Adjusted_P_Value": [1e-3, 1e-2]
            })
        }
        
        interpretations = interpreter.interpret(results)
        assert "GO" in interpretations
        assert "enrichment" in interpretations["GO"].lower()

    def test_mock_summarize_term(self):
        """测试 Mock 条目总结"""
        from allenricher.ai.interpreter import MockInterpreter
        
        interpreter = MockInterpreter()
        summary = interpreter.summarize_term("cell division", ["GENE1", "GENE2"])
        assert "cell division" in summary
        assert "2 genes" in summary

    def test_ai_interpreter_facade(self):
        """测试 AIInterpreter 门面类"""
        from allenricher.ai.interpreter import AIInterpreter
        
        interpreter = AIInterpreter(backend="mock")
        assert interpreter.backend_name == "mock"
        
        results = {
            "KEGG": pd.DataFrame({
                "Term_ID": ["hsa04110"],
                "Term_Name": ["Cell Cycle"],
                "Gene_Count": [15],
                "P_Value": [1e-6],
                "Adjusted_P_Value": [1e-4]
            })
        }
        
        interpretations = interpreter.interpret_results(results)
        assert "KEGG" in interpretations

    def test_create_interpreter_factory(self):
        """测试工厂函数"""
        from allenricher.ai.interpreter import create_interpreter, get_available_backends
        
        backends = get_available_backends()
        assert "mock" in backends
        assert "openai" in backends
        assert "claude" in backends
        assert "deepseek" in backends
        assert "glm" in backends
        assert "minimax" in backends
        
        interpreter = create_interpreter(backend="mock")
        assert interpreter.backend_name == "mock"

    def test_deepseek_interpreter(self):
        """测试 DeepSeek 解释器"""
        from allenricher.ai.interpreter import DeepSeekInterpreter
        
        # 无 API key 时初始化
        interpreter = DeepSeekInterpreter(api_key=None)
        assert interpreter.model == "deepseek-chat"
        
        # 无 API key 时 interpret 返回空字典
        results = {"GO": pd.DataFrame({"Term_ID": ["GO:0005576"], "Term_Name": ["test"], "Gene_Count": [1], "P_Value": [0.01], "Adjusted_P_Value": [0.01]})}
        interpretations = interpreter.interpret(results)
        assert interpretations == {}

    def test_glm_interpreter(self):
        """测试 GLM 解释器"""
        from allenricher.ai.interpreter import GLMInterpreter
        
        interpreter = GLMInterpreter(api_key=None)
        assert interpreter.model == "glm-4"
        
        # 无 API key 时 summarize_term 返回空字符串
        summary = interpreter.summarize_term("test term", ["GENE1"])
        assert summary == ""

    def test_minimax_interpreter(self):
        """测试 MiniMax 解释器"""
        from allenricher.ai.interpreter import MiniMaxInterpreter
        
        interpreter = MiniMaxInterpreter(api_key=None, group_id=None)
        assert interpreter.model == "abab6.5s-chat"
        
        # 无 API key 或 group_id 时 interpret 返回空字典
        results = {"GO": pd.DataFrame({"Term_ID": ["GO:0005576"], "Term_Name": ["test"], "Gene_Count": [1], "P_Value": [0.01], "Adjusted_P_Value": [0.01]})}
        interpretations = interpreter.interpret(results)
        assert interpretations == {}

    def test_ai_interpreter_new_backends(self):
        """测试 AIInterpreter 门面类支持新后端"""
        from allenricher.ai.interpreter import AIInterpreter
        
        # 测试 DeepSeek
        interpreter = AIInterpreter(backend="deepseek")
        assert interpreter.backend_name == "deepseek"
        
        # 测试 GLM
        interpreter = AIInterpreter(backend="glm")
        assert interpreter.backend_name == "glm"
        
        # 测试 MiniMax
        interpreter = AIInterpreter(backend="minimax", group_id="test-group")
        assert interpreter.backend_name == "minimax"

    def test_invalid_backend(self):
        """测试无效后端"""
        from allenricher.ai.interpreter import AIInterpreter
        
        with pytest.raises(ValueError):
            AIInterpreter(backend="invalid_backend")


class TestReportGenerator:
    """测试 HTML 报告生成器"""

    def test_generate_report(self, tmp_path):
        """测试生成报告"""
        from allenricher.report.generator import ReportGenerator
        
        generator = ReportGenerator(str(tmp_path))
        
        results = {
            "GO": pd.DataFrame({
                "Term_ID": ["GO:0005576"],
                "Term_Name": ["extracellular region"],
                "Gene_Count": [10],
                "Rich_Factor": [0.05],
                "P_Value": [1e-5],
                "Adjusted_P_Value": [1e-3],
                "Genes": ["GENE1;GENE2;GENE3"]
            })
        }
        
        output_file = str(tmp_path / "test_report.html")
        result_path = generator.generate(results, output_file, gene_list=["GENE1", "GENE2"])
        
        assert Path(result_path).exists()
        
        with open(result_path) as f:
            content = f.read()
            assert "AllEnricher Report" in content
            assert "GO" in content
            assert "extracellular region" in content

    def test_generate_empty_report(self, tmp_path):
        """测试生成空结果报告"""
        from allenricher.report.generator import ReportGenerator
        
        generator = ReportGenerator(str(tmp_path))
        
        results = {"GO": pd.DataFrame()}
        output_file = str(tmp_path / "empty_report.html")
        result_path = generator.generate(results, output_file)
        
        assert Path(result_path).exists()
        
        with open(result_path) as f:
            content = f.read()
            assert "未找到显著富集的结果" in content or "No enrichment" in content.lower()

    def test_generate_with_ai_interpretation(self, tmp_path):
        """测试带 AI 解读的报告"""
        from allenricher.report.generator import ReportGenerator
        
        generator = ReportGenerator(str(tmp_path))
        
        results = {
            "GO": pd.DataFrame({
                "Term_ID": ["GO:0005576"],
                "Term_Name": ["extracellular region"],
                "Gene_Count": [10],
                "Rich_Factor": [0.05],
                "P_Value": [1e-5],
                "Adjusted_P_Value": [1e-3],
                "Genes": ["GENE1;GENE2"]
            })
        }
        
        ai_interpretation = {"GO": "This is a mock AI interpretation."}
        output_file = str(tmp_path / "ai_report.html")
        result_path = generator.generate(results, output_file, ai_interpretation=ai_interpretation)
        
        with open(result_path) as f:
            content = f.read()
            assert "AI Interpretation" in content
            assert "mock AI interpretation" in content


class TestModels:
    """测试 Pydantic 模型"""

    def test_enrichment_request(self):
        """测试富集分析请求模型"""
        from allenricher.api.server import EnrichmentRequest
        
        request = EnrichmentRequest(
            genes=["TP53", "BRCA1"],
            species="hsa",
            databases=["GO", "KEGG"]
        )
        
        assert request.genes == ["TP53", "BRCA1"]
        assert request.species == "hsa"
        assert request.method == "fisher"
        assert request.pvalue_cutoff == 0.05

    def test_enrichment_response(self):
        """测试富集分析响应模型"""
        from allenricher.api.server import EnrichmentResponse
        
        response = EnrichmentResponse(
            job_id="test-job-id",
            status="pending",
            message="Analysis started"
        )
        
        assert response.job_id == "test-job-id"
        assert response.status == "pending"
        assert response.results is None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
