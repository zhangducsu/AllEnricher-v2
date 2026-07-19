"""
AI Integrated Test

End-to-end testing AI Integrating module with enrichment-analysis pipes, Overwrite:
- MockiInterpreter Full Process: analyze  AI Interpretation  JSON File  HTML Report
- 0 A notable result.: Direct report"No significant enrichment.", Do Not Call AI
- Insufficient 20 Article: Send to Actual Entry AI
- Just in time. 20 Article: Send All 20 Article
- More than 20 Article: Only before sending 20 Article
- AI JSON Output Format: ai_interpretation.json Structure Validation
- HTML Report Embedded: AI Delineating the right paragraph embedded HTML

Use MockInterpreter As AI Backend, No need to be true. API Key or network connection.
"""

import json
import pytest
import sys
from pathlib import Path
from unittest.mock import patch, MagicMock
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))


# ---------------------------------------------------------------------------
# Test data plant
# ---------------------------------------------------------------------------

def _make_results(n_terms=5, include_genes_col=False):
    """Create simulation of enrichment analysis.

Parameters:
n_terms: Number of rich entries
include_genes_col: Does it contain the 'Genes' column?
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
    """Creates an empty result."""
    return {"GO_Biological_Process": pd.DataFrame()}


def _make_multi_db_results(n_terms=5):
    """Creates multi-database results."""
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
# 1. MockInterpreter End-to-end process testing
# ===========================================================================

class TestMockInterpreterEndToEnd:
    """MockiInterpreter Full Process: analyze  AI Interpretation  JSON Document  HTML Report."""

    def test_e2e_interpret_to_json(self, tmp_path):
        """End-to-end: interpret_results & Save as JSON & validation file structure."""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter(backend="mock")
        results = _make_results(n_terms=10)

        # 1. Generate AI interpretation
        ai_interpretation = interpreter.interpret_results(results)

        # 2. Save as JSON file
        json_path = tmp_path / "ai_interpretation.json"
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(ai_interpretation, f, indent=2)

        # 3. Validation of the existence and readability of documents
        assert json_path.exists()

        # 4. Validation of the JPON structure
        with open(json_path, 'r', encoding='utf-8') as f:
            loaded = json.load(f)

        assert isinstance(loaded, dict)
        assert "GO_Biological_Process" in loaded
        assert isinstance(loaded["GO_Biological_Process"], str)
        assert len(loaded["GO_Biological_Process"]) > 0

    def test_e2e_interpret_to_html(self, tmp_path):
        """End-to-end: interpret_results * HTML report embedded * Validation HTML structure."""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter(backend="mock")
        results = _make_results(n_terms=5)

        # 1. Generate AI interpretation
        ai_interpretation = interpreter.interpret_results(results)

        # 2. Generate HTML report paragraph
        html = interpreter.generate_report_section(results)

        # 3. Validation of HTML structures
        assert "<div" in html
        assert "</div>" in html
        assert "AI-Powered Interpretation" in html
        assert "GO_Biological_Process" in html
        assert "ai-disclaimer" in html
        assert "mock" in html  # Backend Name

        # 4. Validation AI interpret content embedded in HTML
        for i in range(5):
            assert f"biological_process_{i}" in html

    def test_e2e_multi_database(self, tmp_path):
        """End-to-end: Multi-databases read JSON+HTML at the same time."""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter(backend="mock")
        results = _make_multi_db_results(n_terms=5)

        # Generate AI reading
        ai_interpretation = interpreter.interpret_results(results)

        # Both databases are read.
        assert "GO_Biological_Process" in ai_interpretation
        assert "KEGG" in ai_interpretation

        # Save JSON
        json_path = tmp_path / "ai_interpretation.json"
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(ai_interpretation, f, indent=2)

        with open(json_path, 'r', encoding='utf-8') as f:
            loaded = json.load(f)

        assert len(loaded) == 2

        # Verify HTML contains two databases
        html = interpreter.generate_report_section(results)
        assert "GO_Biological_Process" in html
        assert "KEGG" in html

    def test_e2e_with_term_summaries(self, tmp_path):
        """End-to-end: full process with term_summaries."""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter(backend="mock")
        results = _make_results(n_terms=5, include_genes_col=True)

        # Generate AI interpretation with term_summaries
        ai_interpretation = interpreter.interpret_results(
            results, include_term_summaries=True
        )

        # Validation contains a holistic interpretation and a summary of entries
        assert "GO_Biological_Process" in ai_interpretation
        assert "GO_Biological_Process_term_summaries" in ai_interpretation

        term_summaries = ai_interpretation["GO_Biological_Process_term_summaries"]
        assert isinstance(term_summaries, dict)
        assert len(term_summaries) == 5

        # Save JSON and verify structure
        json_path = tmp_path / "ai_interpretation.json"
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(ai_interpretation, f, indent=2)

        with open(json_path, 'r', encoding='utf-8') as f:
            loaded = json.load(f)

        # Term_summaries should be an embedded dictionary
        assert isinstance(loaded["GO_Biological_Process_term_summaries"], dict)
        for term_name, summary in loaded["GO_Biological_Process_term_summaries"].items():
            assert isinstance(term_name, str)
            assert isinstance(summary, str)
            assert len(summary) > 0


# ===========================================================================
# 2.0 Tests for Visible Results
# ===========================================================================

class TestZeroResults:
    """0 scenes of significant results: report directly "no significant enrichment" and do not call AI."""

    def test_mock_zero_results_skips_interpretation(self):
        """MockInterpreter: 0 results do not generate an interpretation."""
        from allenricher.ai.interpreter import MockInterpreter

        interpreter = MockInterpreter()
        results = _make_empty_results()

        interpretations = interpreter.interpret(results)

        # MockInterpreter directs the result to empty
        assert len(interpretations) == 0

    def test_openai_zero_results_no_api_call(self):
        """OpenAI: 0 results do not call API when you want to return the fixed message."""
        from allenricher.ai.interpreter import OpenAIInterpreter

        interpreter = OpenAIInterpreter(api_key="test-key")
        results = _make_empty_results()

        with patch.object(interpreter, "_call_api") as mock_call:
            interpretations = interpreter.interpret(results)

            # I'm going to return to "no significant enrichment."
            assert "GO_Biological_Process" in interpretations
            assert "No enrichment terms were available" in interpretations["GO_Biological_Process"]
            # API should not be called
            mock_call.assert_not_called()

    def test_claude_zero_results_no_api_call(self):
        """Claude: 0 results do not call API."""
        from allenricher.ai.interpreter import ClaudeInterpreter

        interpreter = ClaudeInterpreter(api_key="test-key")
        results = _make_empty_results()

        with patch.object(interpreter, "_call_api") as mock_call:
            interpretations = interpreter.interpret(results)

            assert "GO_Biological_Process" in interpretations
            assert "No enrichment terms were available" in interpretations["GO_Biological_Process"]
            mock_call.assert_not_called()

    def test_ollama_zero_results_no_api_call(self):
        """Ollama: 0 results do not call API."""
        from allenricher.ai.interpreter import OllamaInterpreter

        interpreter = OllamaInterpreter()
        results = _make_empty_results()

        with patch.object(interpreter, "_call_api") as mock_call:
            interpretations = interpreter.interpret(results)

            assert "GO_Biological_Process" in interpretations
            assert "No enrichment terms were available" in interpretations["GO_Biological_Process"]
            mock_call.assert_not_called()

    def test_facade_zero_results(self):
        """Door class: 0 results return to empty dictionary (MockInterpreter behaviour)."""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_empty_results()

        interpretations = interpreter.interpret_results(results)
        assert len(interpretations) == 0

    def test_facade_zero_results_no_term_summaries(self):
        """Dock: 0 results do not generate term_summaries."""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_empty_results()

        interpretations = interpreter.interpret_results(
            results, include_term_summaries=True
        )

        assert len(interpretations) == 0
        assert "GO_Biological_Process_term_summaries" not in interpretations

    def test_zero_results_json_output(self, tmp_path):
        """0 When result is verified by JSON file structure."""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter(backend="mock")
        results = _make_empty_results()

        ai_interpretation = interpreter.interpret_results(results)

        json_path = tmp_path / "ai_interpretation.json"
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(ai_interpretation, f, indent=2)

        with open(json_path, 'r', encoding='utf-8') as f:
            loaded = json.load(f)

        # MockInterpreter returns empty dictionary
        assert loaded == {}

    def test_zero_results_html_empty_section(self):
        """0 The HTML report does not contain an interpretation when it comes to the outcome."""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_empty_results()

        html = interpreter.generate_report_section(results)

        # There should be AI headers but no database reading blocks
        assert "AI-Powered Interpretation" in html
        assert "GO_Biological_Process" not in html


# ===========================================================================
# 3. 20 results tests
# ===========================================================================

class TestFewerThan20Results:
    """Less than 20: Send to AI by actual entry."""

    def test_5_results_sends_5(self):
        """Include all five terms when a result table contains five terms."""
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

            # Validate "Top 5" in prompt
            assert "Top 5" in prompt
            # Every term is represented in the prompt.
            for i in range(5):
                assert f"biological_process_{i}" in prompt
            # Should not contain article 6
            assert "biological_process_5" not in prompt

    def test_15_results_sends_15(self):
        """Include all 15 terms when the result is below the prompt limit."""
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
        """When a result is reached, the prompt contains one."""
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
        """MockInterpreter: 5 results show all 5 outcomes."""
        from allenricher.ai.interpreter import MockInterpreter

        interpreter = MockInterpreter()
        results = _make_results(n_terms=5)

        interpretations = interpreter.interpret(results)
        text = interpretations["GO_Biological_Process"]

        for i in range(5):
            assert f"biological_process_{i}" in text

    def test_facade_fewer_than_20_term_summaries(self):
        """Generate one facade summary per term when fewer than 20 terms are available."""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=8, include_genes_col=True)

        interpretations = interpreter.interpret_results(
            results, include_term_summaries=True
        )

        term_summaries = interpretations["GO_Biological_Process_term_summaries"]
        assert len(term_summaries) == 8


# ===========================================================================
# 4. Just in time 20 results tests
# ===========================================================================

class TestExactly20Results:
    """Just 20: send all 20."""

    def test_openai_20_results_sends_all(self):
        """OpenAI: When 20 are right, the prompt contains all 20."""
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
        """MockInterpreter: Show all 20 at exactly 20 times."""
        from allenricher.ai.interpreter import MockInterpreter

        interpreter = MockInterpreter()
        results = _make_results(n_terms=20)

        interpretations = interpreter.interpret(results)
        text = interpretations["GO_Biological_Process"]

        for i in range(20):
            assert f"biological_process_{i}" in text

    def test_facade_20_term_summaries(self):
        """Door class: 20 when 20 times the term_summaries are created."""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=20, include_genes_col=True)

        interpretations = interpreter.interpret_results(
            results, include_term_summaries=True
        )

        term_summaries = interpretations["GO_Biological_Process_term_summaries"]
        assert len(term_summaries) == 20


# ===========================================================================
# 5. More than 20 results tested
# ===========================================================================

class TestMoreThan20Results:
    """More than 20: Only 20 before sending (Top 20 cut)."""

    def test_openai_25_results_truncates_to_20(self):
        """Limit the OpenAI prompt to the top 20 of 25 terms."""
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

            # For Top 20
            assert "Top 20" in prompt
            # Only the top 20 terms should appear.
            for i in range(20):
                assert f"biological_process_{i}" in prompt
            # Article 21
            assert "biological_process_20" not in prompt

    def test_mock_25_results_truncates_to_20(self):
        """MockiInterpreter: 25 results show only the first 20."""
        from allenricher.ai.interpreter import MockInterpreter

        interpreter = MockInterpreter()
        results = _make_results(n_terms=25)

        interpretations = interpreter.interpret(results)
        text = interpretations["GO_Biological_Process"]

        for i in range(20):
            assert f"biological_process_{i}" in text
        assert "biological_process_20" not in text

    def test_facade_25_term_summaries_truncates_to_20(self):
        """Domain: 25 times time timeterm_summaries only 20 before being created."""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=25, include_genes_col=True)

        interpretations = interpreter.interpret_results(
            results, include_term_summaries=True
        )

        term_summaries = interpretations["GO_Biological_Process_term_summaries"]
        assert len(term_summaries) == 20
        # The first 20 entries should be named
        for i in range(20):
            assert f"biological_process_{i}" in term_summaries
        # 21 should not be
        assert "biological_process_20" not in term_summaries


# ===========================================================================
# 6. AI JSON Output Format Validation
# ===========================================================================

class TestAIJsonOutputFormat:
    """i_interpretation.json structure authenticated."""

    def test_json_structure_single_db(self, tmp_path):
        """Single database JSON structure."""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter(backend="mock")
        results = _make_results(n_terms=5)

        ai_interpretation = interpreter.interpret_results(results)

        json_path = tmp_path / "ai_interpretation.json"
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(ai_interpretation, f, indent=2)

        with open(json_path, 'r', encoding='utf-8') as f:
            loaded = json.load(f)

        # Top floor is dict
        assert isinstance(loaded, dict)
        # Key is the database name.
        assert "GO_Biological_Process" in loaded
        # Value is string
        assert isinstance(loaded["GO_Biological_Process"], str)
        # Nothing is empty.
        assert len(loaded["GO_Biological_Process"]) > 50

    def test_json_structure_multi_db(self, tmp_path):
        """Multi-database JSON structure."""
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
        """The JSON structure with term_summaries."""
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

        # Overall reading
        assert "GO_Biological_Process" in loaded
        assert isinstance(loaded["GO_Biological_Process"], str)

        # Term_summaries Embedded Dictionary
        assert "GO_Biological_Process_term_summaries" in loaded
        summaries = loaded["GO_Biological_Process_term_summaries"]
        assert isinstance(summaries, dict)
        assert len(summaries) == 3

        # Each value for a sumary is a string
        for key, value in summaries.items():
            assert isinstance(key, str)
            assert isinstance(value, str)

    def test_json_is_valid_utf8(self, tmp_path):
        """JSON files are coded using UTF-8 and can be read and write correctly."""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter(backend="mock")
        results = _make_results(n_terms=5)

        ai_interpretation = interpreter.interpret_results(results)

        json_path = tmp_path / "ai_interpretation.json"
        # Writing UTF-8
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(ai_interpretation, f, indent=2, ensure_ascii=False)

        # Read UTF-8
        with open(json_path, 'r', encoding='utf-8') as f:
            loaded = json.load(f)

        assert loaded == ai_interpretation

    def test_json_empty_results(self, tmp_path):
        """The result when JSON is an empty dictionary."""
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
# 7. Embedding of HTML reports for validation
# ===========================================================================

class TestHTMLReportEmbedding:
    """AI interprets the paragraph correctly embedded in the HTML report."""

    def test_html_contains_ai_section(self):
        """HTML contains AI's reading paragraph."""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=3)

        html = interpreter.generate_report_section(results)

        # Validation AI paragraph tag
        assert 'id="ai-interpretation"' in html
        assert "AI-Powered Interpretation" in html

    def test_html_contains_disclaimer(self):
        """HTML contains AI disclaimer."""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=3)

        html = interpreter.generate_report_section(results)

        assert "ai-disclaimer" in html
        assert "reviewed by domain experts" in html

    def test_html_contains_backend_name(self):
        """HTML contains the backend name."""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=3)

        html = interpreter.generate_report_section(results)

        assert "mock" in html

    def test_html_contains_db_interpretation(self):
        """HTML contains a database interpretation."""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=5)

        html = interpreter.generate_report_section(results)

        # Should contain database reading blocks
        assert "GO_Biological_Process" in html
        # Should contain rich entries
        for i in range(5):
            assert f"biological_process_{i}" in html

    def test_html_multi_db_sections(self):
        """HTML contains reading paragraphs in multiple databases."""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_multi_db_results(n_terms=3)

        html = interpreter.generate_report_section(results)

        assert "GO_Biological_Process" in html
        assert "KEGG" in html

    def test_html_newlines_converted_to_br(self):
        """HTML Centre line breaks converted to <br>."""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=3)

        html = interpreter.generate_report_section(results)

        # The text generated by MockiInterpreter contains line breaks and should be converted to<br>
        assert "<br>" in html

    def test_html_empty_results_no_db_block(self):
        """When you have no results, HTML does not contain the database interpretation block."""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_empty_results()

        html = interpreter.generate_report_section(results)

        # There should be AI title and disclaimer
        assert "AI-Powered Interpretation" in html
        assert "ai-disclaimer" in html
        # There should be no database reading blocks
        assert "GO_Biological_Process" not in html

    def test_html_excludes_term_summaries(self):
        """HTML does not contain term_summaries (shows only overall reading)."""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=3, include_genes_col=True)

        # Even if I asked for it, it shouldn't appear in HTML.
        html = interpreter.generate_report_section(results)

        assert "_term_summaries" not in html

    def test_html_is_valid_fragment(self):
        """HTML snippet is complete (with a closed label)."""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=3)

        html = interpreter.generate_report_section(results)

        # Validation<div>Tab Match
        assert html.strip().startswith("<div")
        assert html.strip().endswith("</div>")
        # Validate <h2> label
        assert "<h2>" in html
        assert "</h2>" in html


# ===========================================================================
# 8. Mixed scenario test (some databases have results and some are empty)
# ===========================================================================

class TestMixedResults:
    """Mixed scene: some databases have results and some are empty."""

    def test_mixed_empty_and_non_empty(self):
        """Some databases have results and some are empty."""
        from allenricher.ai.interpreter import MockInterpreter

        interpreter = MockInterpreter()
        results = {
            "GO_Biological_Process": pd.DataFrame({
                "Term_Name": ["process_1", "process_2"],
                "P_Value": [1e-5, 1e-3],
                "Gene_Count": [10, 5],
            }),
            "KEGG": pd.DataFrame(),  # Empty
        }

        interpretations = interpreter.interpret(results)

        # Go has an explanation.
        assert "GO_Biological_Process" in interpretations
        # KEG is empty, MockInterpreter does not enter dictionary
        assert "KEGG" not in interpretations

    def test_mixed_openai_empty_skips_api(self):
        """OpenAI: API is not called in the Fusion Database for mixed scenarios."""
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

            # GO should call API once
            assert mock_client.chat.completions.create.call_count == 1
            # Kegg should return fixed message
            assert "KEGG" in interpretations
            assert "No enrichment terms were available" in interpretations["KEGG"]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
