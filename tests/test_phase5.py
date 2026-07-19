"""
Phase 5 module testing: API, AI Interpretation, HTML Report
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
    """Test API Endpoint"""

    def test_app_creation(self):
        """Test FastAPI application creation"""
        from allenricher.api.server import app
        assert app is not None
        assert app.title == "AllEnricher API"

    def test_root_endpoint(self):
        """Test Roots"""
        from allenricher.api.server import app
        from fastapi.testclient import TestClient
        
        client = TestClient(app)
        response = client.get("/")
        assert response.status_code == 200
        # Root route now returns HTML (Web interface)
        content_type = response.headers.get("content-type", "")
        assert "text/html" in content_type or "application/json" in content_type

    def test_species_endpoint(self):
        """Test species list endpoint"""
        from allenricher.api.server import app
        from fastapi.testclient import TestClient
        
        client = TestClient(app)
        response = client.get("/api/species")
        assert response.status_code == 200
        data = response.json()
        assert isinstance(data, list)
        assert any(s["code"] == "hsa" for s in data)

    def test_databases_endpoint(self):
        """Test database list endpoints"""
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
        """Could not close temporary folder: %s"""
        from allenricher.api.server import app
        from fastapi.testclient import TestClient
        
        client = TestClient(app)
        response = client.get("/api/status/nonexistent-job-id")
        assert response.status_code == 404


class TestAIInterpreter:
    """Test AI Interpretation Module"""

    def test_mock_interpreter(self):
        """Tests Mock interpreter"""
        from allenricher.ai.interpreter import MockInterpreter
        
        interpreter = MockInterpreter()
        
        # Create simulation results
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
        """Test Lock Entry Summary"""
        from allenricher.ai.interpreter import MockInterpreter
        
        interpreter = MockInterpreter()
        summary = interpreter.summarize_term("cell division", ["GENE1", "GENE2"])
        assert "cell division" in summary
        assert "2 genes" in summary

    def test_ai_interpreter_facade(self):
        """Test AAIInterpreter"""
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
        """Test Factor function"""
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
        """Test DeepSeek Interpreter"""
        from allenricher.ai.interpreter import DeepSeekInterpreter
        
        # Initialize without API key
        interpreter = DeepSeekInterpreter(api_key=None)
        assert interpreter.model == "deepseek-chat"
        
        # returns empty dictionary when no API key
        results = {"GO": pd.DataFrame({"Term_ID": ["GO:0005576"], "Term_Name": ["test"], "Gene_Count": [1], "P_Value": [0.01], "Adjusted_P_Value": [0.01]})}
        interpretations = interpreter.interpret(results)
        assert interpretations == {}

    def test_glm_interpreter(self):
        """Test GLM Interpreters"""
        from allenricher.ai.interpreter import GLMInterpreter
        
        interpreter = GLMInterpreter(api_key=None)
        assert interpreter.model == "glm-4"
        
        # Sommarize_term returns empty string when no API key
        summary = interpreter.summarize_term("test term", ["GENE1"])
        assert summary == ""

    def test_minimax_interpreter(self):
        """Test MiniMax interpreter"""
        from allenricher.ai.interpreter import MiniMaxInterpreter
        
        interpreter = MiniMaxInterpreter(api_key=None, group_id=None)
        assert interpreter.model == "abab6.5s-chat"
        
        # returns empty dictionary without API key or group_id
        results = {"GO": pd.DataFrame({"Term_ID": ["GO:0005576"], "Term_Name": ["test"], "Gene_Count": [1], "P_Value": [0.01], "Adjusted_P_Value": [0.01]})}
        interpretations = interpreter.interpret(results)
        assert interpretations == {}

    def test_ai_interpreter_new_backends(self):
        """Test AIInterpreter Gate support new backend"""
        from allenricher.ai.interpreter import AIInterpreter
        
        # Test DeepSeek
        interpreter = AIInterpreter(backend="deepseek")
        assert interpreter.backend_name == "deepseek"
        
        # Test GLM
        interpreter = AIInterpreter(backend="glm")
        assert interpreter.backend_name == "glm"
        
        # Test MiniMax
        interpreter = AIInterpreter(backend="minimax", group_id="test-group")
        assert interpreter.backend_name == "minimax"

    def test_invalid_backend(self):
        """Test invalid backend"""
        from allenricher.ai.interpreter import AIInterpreter
        
        with pytest.raises(ValueError):
            AIInterpreter(backend="invalid_backend")


class TestReportGenerator:
    """Test HTML Report Generator"""

    def test_generate_report(self, tmp_path):
        """Test Generation Report"""
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
        
        with open(result_path, encoding="utf-8") as f:
            content = f.read()
            assert "AllEnricher Report" in content
            assert "GO" in content
            assert "extracellular region" in content

    def test_generate_empty_report(self, tmp_path):
        """Test produces empty results report"""
        from allenricher.report.generator import ReportGenerator
        
        generator = ReportGenerator(str(tmp_path))
        
        results = {"GO": pd.DataFrame()}
        output_file = str(tmp_path / "empty_report.html")
        result_path = generator.generate(results, output_file)
        
        assert Path(result_path).exists()
        
        with open(result_path, encoding="utf-8") as f:
            content = f.read()
            assert "No significant enrichment results found" in content or "No enrichment" in content.lower()

    def test_generate_with_ai_interpretation(self, tmp_path):
        """Test belt AI reading report"""
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
        
        with open(result_path, encoding="utf-8") as f:
            content = f.read()
            assert "AI Interpretation" in content
            assert "mock AI interpretation" in content


class TestModels:
    """Test Pydantic Model"""

    def test_enrichment_request(self):
        """Test the enrichment analysis request model"""
        from allenricher.api.server import EnrichmentRequest
        
        request = EnrichmentRequest(
            genes=["TP53", "BRCA1"],
            species="hsa",
            databases=["GO", "KEGG"]
        )
        
        assert request.genes == ["TP53", "BRCA1"]
        assert request.species == "hsa"
        assert request.method == "hypergeometric"
        assert request.pvalue_cutoff == 0.05

    def test_enrichment_response(self):
        """Test the Fuzzy Analytical Response Model"""
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
