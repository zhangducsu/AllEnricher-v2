"""Tests for AI interpretation backends and the interpreter facade.

External clients are mocked so the suite runs without network access or API keys.
"""

import pytest
import sys
from pathlib import Path
from unittest.mock import patch, MagicMock, PropertyMock
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))


# ---------------------------------------------------------------------------
# Shared result fixtures
# ---------------------------------------------------------------------------

def _make_results(n_terms=5, include_genes_col=False):
    """Create a synthetic enrichment result table.

    Args:
        n_terms: Number of enrichment terms.
        include_genes_col: Include a semicolon-delimited ``Genes`` column.
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
    """Create an empty enrichment result table."""
    return {"GO_Biological_Process": pd.DataFrame()}


# ===========================================================================
# 1. MockInterpreter Test
# ===========================================================================

class TestMockInterpreter:
    """Tests for deterministic mock interpretations."""

    def test_interpret_normal(self):
        """Verify interpretation of a non-empty result table."""
        from allenricher.ai.interpreter import MockInterpreter

        interpreter = MockInterpreter()
        results = _make_results(n_terms=5)
        interpretations = interpreter.interpret(results)

        assert "GO_Biological_Process" in interpretations
        assert "enrichment" in interpretations["GO_Biological_Process"].lower()
        assert "5 terms" in interpretations["GO_Biological_Process"]

    def test_interpret_empty_results(self):
        """Verify that empty database results do not produce an interpretation."""
        from allenricher.ai.interpreter import MockInterpreter

        interpreter = MockInterpreter()
        results = _make_empty_results()
        interpretations = interpreter.interpret(results)

        # Empty database results are omitted from the interpretation mapping.
        assert "GO_Biological_Process" not in interpretations

    def test_interpret_fewer_than_20(self):
        """When the number of interpret() is not enough, the number of items is displayed in real terms."""
        from allenricher.ai.interpreter import MockInterpreter

        interpreter = MockInterpreter()
        n = 8
        results = _make_results(n_terms=n)
        interpretations = interpreter.interpret(results)

        text = interpretations["GO_Biological_Process"]
        # should contain all n entry names
        for i in range(n):
            assert f"biological_process_{i}" in text

    def test_interpret_exactly_20(self):
        """Interpret() is displayed at exactly 20 times."""
        from allenricher.ai.interpreter import MockInterpreter

        interpreter = MockInterpreter()
        results = _make_results(n_terms=20)
        interpretations = interpreter.interpret(results)

        text = interpretations["GO_Biological_Process"]
        for i in range(20):
            assert f"biological_process_{i}" in text

    def test_interpret_more_than_20(self):
        """Only the previous 20 when the interpret() exceeds 20."""
        from allenricher.ai.interpreter import MockInterpreter

        interpreter = MockInterpreter()
        results = _make_results(n_terms=25)
        interpretations = interpreter.interpret(results)

        text = interpretations["GO_Biological_Process"]
        # Top 20 should appear
        for i in range(20):
            assert f"biological_process_{i}" in text
        # 21 should not appear (head (20) cut)
        assert "biological_process_20" not in text

    def test_summarize_term(self):
        """Summarize_term() is working normally."""
        from allenricher.ai.interpreter import MockInterpreter

        interpreter = MockInterpreter()
        summary = interpreter.summarize_term("cell division", ["GENE1", "GENE2", "GENE3"])

        assert "cell division" in summary
        assert "3 genes" in summary

    def test_summarize_term_empty_genes(self):
        """Qualified by the name of the name of the person."""
        from allenricher.ai.interpreter import MockInterpreter

        interpreter = MockInterpreter()
        summary = interpreter.summarize_term("apoptosis", [])

        assert "apoptosis" in summary
        assert "0 genes" in summary

    def test_create_interpreter_mock(self):
        """Create_interpreter("mock") returns the MockInterpreter backend."""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter("mock")
        assert interpreter.backend_name == "mock"
        assert isinstance(interpreter.interpreter, type(interpreter.interpreter))


# ===========================================================================
# OpenAIInterpreter test
# ===========================================================================

class TestOpenAIInterpreter:
    """Test OpenAIInterpreter (mock openaii package)."""

    def test_call_api_success(self):
        """_call_api() returns the normal call."""
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
        """_call_api() ImportError downgraded."""
        from allenricher.ai.interpreter import OpenAIInterpreter

        interpreter = OpenAIInterpreter(api_key="test-key")

        with patch.dict("sys.modules", {"openai": None}):
            # Let the import openaii throw ImportError
            with patch("builtins.__import__", side_effect=ImportError("No module named openai")):
                result = interpreter._call_api("test prompt")
                assert "Error" in result
                assert "not installed" in result

    def test_call_api_exception(self):
        """_call_api() other anomalies."""
        from allenricher.ai.interpreter import OpenAIInterpreter

        interpreter = OpenAIInterpreter(api_key="test-key")

        with patch.dict("sys.modules", {"openai": MagicMock()}):
            import openai as mock_openai
            mock_openai.OpenAI.side_effect = Exception("API connection failed")

            result = interpreter._call_api("test prompt")
            assert "Error" in result
            assert "API connection failed" in result

    def test_interpret_extracts_top_20(self):
        """Interpret() extracts Top 20 entries."""
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
            # The prompt includes at most the top 20 terms.
            call_args = mock_client.chat.completions.create.call_args
            prompt = call_args[1]["messages"][1]["content"]
            assert "Top 20" in prompt
            assert "biological_process_0" in prompt
            assert "biological_process_19" in prompt
            assert "biological_process_20" not in prompt

    def test_interpret_empty_results_no_api_call(self):
        """An empty result does not call API."""
        from allenricher.ai.interpreter import OpenAIInterpreter

        interpreter = OpenAIInterpreter(api_key="test-key")
        results = _make_empty_results()

        with patch.object(interpreter, "_call_api") as mock_call:
            interpretations = interpreter.interpret(results)

            assert "GO_Biological_Process" in interpretations
            assert "No enrichment terms were available" in interpretations["GO_Biological_Process"]
            mock_call.assert_not_called()

    def test_interpret_no_api_key(self):
        """returns an empty dictionary when you have no API key."""
        from allenricher.ai.interpreter import OpenAIInterpreter

        interpreter = OpenAIInterpreter(api_key=None)
        results = _make_results(n_terms=5)

        with patch.dict("os.environ", {"OPENAI_API_KEY": ""}):
            interpreter.api_key = None
            interpretations = interpreter.interpret(results)
            assert interpretations == {}

    def test_summarize_term_with_api_key(self):
        """The submarize_term() is normally called when API key is available."""
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
        """The empty string is returned when the sumarize_term() does not have API key."""
        from allenricher.ai.interpreter import OpenAIInterpreter

        interpreter = OpenAIInterpreter(api_key=None)
        interpreter.api_key = None

        result = interpreter.summarize_term("cell division", ["GENE1"])
        assert result == ""

    def test_summarize_term_truncates_long_gene_list(self):
        """Summarize_term() is cut when more than 10 genes are passed."""
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
            # The first 10 genes should appear, the eleventh should not appear.
            assert "GENE9" in prompt
            assert "GENE10" not in prompt
            assert "..." in prompt


# ===========================================================================
# 3. ClaudeInterpreter test
# ===========================================================================

class TestClaudeInterpreter:
    """Test Claude Interpreter (mock anthropological package)."""

    def test_call_api_success(self):
        """_call_api() returns the normal call."""
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
        """_call_api() ImportError downgraded."""
        from allenricher.ai.interpreter import ClaudeInterpreter

        interpreter = ClaudeInterpreter(api_key="test-key")

        with patch("builtins.__import__", side_effect=ImportError("No module named anthropic")):
            result = interpreter._call_api("test prompt")
            assert "Error" in result
            assert "not installed" in result

    def test_call_api_exception(self):
        """_call_api() other anomalies."""
        from allenricher.ai.interpreter import ClaudeInterpreter

        interpreter = ClaudeInterpreter(api_key="test-key")

        with patch.dict("sys.modules", {"anthropic": MagicMock()}):
            import anthropic as mock_anthropic
            mock_anthropic.Anthropic.side_effect = Exception("Claude API error")

            result = interpreter._call_api("test prompt")
            assert "Error" in result
            assert "Claude API error" in result

    def test_interpret_no_api_key(self):
        """returns an empty dictionary when you have no API key."""
        from allenricher.ai.interpreter import ClaudeInterpreter

        interpreter = ClaudeInterpreter(api_key=None)
        interpreter.api_key = None
        results = _make_results(n_terms=3)

        interpretations = interpreter.interpret(results)
        assert interpretations == {}

    def test_interpret_empty_results(self):
        """An empty result does not call API."""
        from allenricher.ai.interpreter import ClaudeInterpreter

        interpreter = ClaudeInterpreter(api_key="test-key")
        results = _make_empty_results()

        with patch.object(interpreter, "_call_api") as mock_call:
            interpretations = interpreter.interpret(results)

            assert "GO_Biological_Process" in interpretations
            assert "No enrichment terms were available" in interpretations["GO_Biological_Process"]
            mock_call.assert_not_called()

    def test_summarize_term_no_api_key(self):
        """The empty string is returned when the sumarize_term() does not have API key."""
        from allenricher.ai.interpreter import ClaudeInterpreter

        interpreter = ClaudeInterpreter(api_key=None)
        interpreter.api_key = None

        result = interpreter.summarize_term("apoptosis", ["GENE1"])
        assert result == ""


# ===========================================================================
# 4. OlamaInterpreter test
# ===========================================================================

class TestOllamaInterpreter:
    """Test OlamaInterpreter (mock recests.post)."""

    def test_call_api_success(self):
        """_call_api() returns the normal call."""
        from allenricher.ai.interpreter import OllamaInterpreter

        interpreter = OllamaInterpreter(model="llama2")

        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {"response": "Ollama generated text."}

        with patch("requests.post", return_value=mock_response) as mock_post:
            result = interpreter._call_api("test prompt")
            assert result == "Ollama generated text."
            mock_post.assert_called_once()
            # Validation request URL and parameters
            call_args = mock_post.call_args
            assert "localhost:11434/api/generate" in call_args[0][0]
            assert call_args[1]["json"]["model"] == "llama2"
            assert call_args[1]["json"]["stream"] is False

    def test_call_api_connection_failure(self):
        """_call_api() Connection failed to process."""
        from allenricher.ai.interpreter import OllamaInterpreter

        interpreter = OllamaInterpreter()

        with patch("requests.post", side_effect=Exception("Connection refused")):
            result = interpreter._call_api("test prompt")
            assert "Error" in result
            assert "Connection refused" in result

    def test_call_api_http_error(self):
        """_call_api() HTTP error state code processing."""
        from allenricher.ai.interpreter import OllamaInterpreter

        interpreter = OllamaInterpreter()

        mock_response = MagicMock()
        mock_response.status_code = 500

        with patch("requests.post", return_value=mock_response):
            result = interpreter._call_api("test prompt")
            assert "Error" in result
            assert "500" in result

    def test_call_api_import_error(self):
        """_call_api() downgrade when not installed."""
        from allenricher.ai.interpreter import OllamaInterpreter

        interpreter = OllamaInterpreter()

        with patch("builtins.__import__", side_effect=ImportError("No module named requests")):
            result = interpreter._call_api("test prompt")
            assert "Error" in result
            assert "not installed" in result

    def test_interpret_empty_results(self):
        """An empty result does not call API."""
        from allenricher.ai.interpreter import OllamaInterpreter

        interpreter = OllamaInterpreter()
        results = _make_empty_results()

        with patch.object(interpreter, "_call_api") as mock_call:
            interpretations = interpreter.interpret(results)

            assert "GO_Biological_Process" in interpretations
            assert "No enrichment terms were available" in interpretations["GO_Biological_Process"]
            mock_call.assert_not_called()

    def test_interpret_normal(self):
        """Interpret() is working."""
        from allenricher.ai.interpreter import OllamaInterpreter

        interpreter = OllamaInterpreter()
        results = _make_results(n_terms=3)

        with patch.object(interpreter, "_call_api", return_value="Ollama interpretation."):
            interpretations = interpreter.interpret(results)

            assert "GO_Biological_Process" in interpretations
            assert interpretations["GO_Biological_Process"] == "Ollama interpretation."

    def test_summarize_term(self):
        """Summarize_term() is working normally."""
        from allenricher.ai.interpreter import OllamaInterpreter

        interpreter = OllamaInterpreter()

        with patch.object(interpreter, "_call_api", return_value="Term description."):
            result = interpreter.summarize_term("apoptosis", ["GENE1"])
            assert result == "Term description."


# ===========================================================================
# DeepSeek / GLM / MiniMaxInterpreter
# ===========================================================================

class TestDeepSeekInterpreter:
    """Test DeepSeekInterpreter."""

    def test_init_no_api_key(self):
        """API key is not downgraded at initialization."""
        from allenricher.ai.interpreter import DeepSeekInterpreter

        interpreter = DeepSeekInterpreter(api_key=None)
        assert interpreter.api_key is None
        assert interpreter.model == "deepseek-chat"

    def test_interpret_no_api_key(self):
        """Returns empty dictionary when no API key() is available."""
        from allenricher.ai.interpreter import DeepSeekInterpreter

        interpreter = DeepSeekInterpreter(api_key=None)
        interpreter.api_key = None
        results = _make_results(n_terms=3)

        interpretations = interpreter.interpret(results)
        assert interpretations == {}

    def test_summarize_term_no_api_key(self):
        """The empty string returns the summarize_term() when there is no API key."""
        from allenricher.ai.interpreter import DeepSeekInterpreter

        interpreter = DeepSeekInterpreter(api_key=None)
        interpreter.api_key = None

        result = interpreter.summarize_term("test", ["GENE1"])
        assert result == ""

    def test_call_api_success(self):
        """_call_api() is normally called (in openaii SDK compatible format)."""
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
            # Validation uses base URL for DeepSeek
            mock_openai.OpenAI.assert_called_once_with(
                api_key="test-key",
                base_url="https://api.deepseek.com"
            )

    def test_call_api_import_error(self):
        """_call_api() ImportError downgrade."""
        from allenricher.ai.interpreter import DeepSeekInterpreter

        interpreter = DeepSeekInterpreter(api_key="test-key")

        with patch("builtins.__import__", side_effect=ImportError("No module named openai")):
            result = interpreter._call_api("test prompt")
            assert "Error" in result
            assert "not installed" in result


class TestGLMInterpreter:
    """Test GLMInterpreter."""

    def test_init_no_api_key(self):
        """API key is not downgraded at initialization."""
        from allenricher.ai.interpreter import GLMInterpreter

        interpreter = GLMInterpreter(api_key=None)
        assert interpreter.api_key is None
        assert interpreter.model == "glm-4"

    def test_interpret_no_api_key(self):
        """Returns empty dictionary when no API key() is available."""
        from allenricher.ai.interpreter import GLMInterpreter

        interpreter = GLMInterpreter(api_key=None)
        interpreter.api_key = None
        results = _make_results(n_terms=3)

        interpretations = interpreter.interpret(results)
        assert interpretations == {}

    def test_summarize_term_no_api_key(self):
        """The empty string returns the summarize_term() when there is no API key."""
        from allenricher.ai.interpreter import GLMInterpreter

        interpreter = GLMInterpreter(api_key=None)
        interpreter.api_key = None

        result = interpreter.summarize_term("test", ["GENE1"])
        assert result == ""

    def test_call_api_success(self):
        """_call_api() is on normal call."""
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
            # Validate base_url using GLM
            mock_openai.OpenAI.assert_called_once_with(
                api_key="test-key",
                base_url="https://open.bigmodel.cn/api/paas/v4"
            )


class TestMiniMaxInterpreter:
    """Test MiniMax Interpreter."""

    def test_init_no_api_key(self):
        """API key is not downgraded at initialization."""
        from allenricher.ai.interpreter import MiniMaxInterpreter

        interpreter = MiniMaxInterpreter(api_key=None, group_id=None)
        assert interpreter.api_key is None
        assert interpreter.group_id is None
        assert interpreter.model == "abab6.5s-chat"

    def test_interpret_no_api_key(self):
        """Returns empty dictionary when no API key or group_id is available."""
        from allenricher.ai.interpreter import MiniMaxInterpreter

        interpreter = MiniMaxInterpreter(api_key=None, group_id=None)
        results = _make_results(n_terms=3)

        interpretations = interpreter.interpret(results)
        assert interpretations == {}

    def test_interpret_no_group_id(self):
        """API key returns an empty dictionary when there is no group_id"""
        from allenricher.ai.interpreter import MiniMaxInterpreter

        interpreter = MiniMaxInterpreter(api_key="test-key", group_id=None)
        results = _make_results(n_terms=3)

        interpretations = interpreter.interpret(results)
        assert interpretations == {}

    def test_summarize_term_no_api_key(self):
        """The empty string returns the summarize_term() when there is no API key."""
        from allenricher.ai.interpreter import MiniMaxInterpreter

        interpreter = MiniMaxInterpreter(api_key=None, group_id=None)

        result = interpreter.summarize_term("test", ["GENE1"])
        assert result == ""

    def test_call_api_no_credentials(self):
        """_call_api() returns error when no documents are available."""
        from allenricher.ai.interpreter import MiniMaxInterpreter

        interpreter = MiniMaxInterpreter(api_key=None, group_id=None)

        # %1 requires a lock openanai package to make the import successful, but api_key/group_idCheck before calling
        with patch.dict("sys.modules", {"openai": MagicMock()}):
            result = interpreter._call_api("test prompt")
            assert "Error" in result
            assert "not configured" in result

    def test_call_api_success(self):
        """_call_api() is on normal call."""
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
# 6. AIInterpreter Dock Test
# ===========================================================================

class TestAIInterpreterFacade:
    """Test the AIInterpreter frontal class."""

    def test_interpret_results_calls_backend(self):
        """Interpret_revers() calls backend and returns result."""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=3)

        interpretations = interpreter.interpret_results(results)

        assert "GO_Biological_Process" in interpretations
        assert "enrichment" in interpretations["GO_Biological_Process"].lower()

    def test_interpret_results_empty(self):
        """An empty result process."""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_empty_results()

        interpretations = interpreter.interpret_results(results)

        # MockiInterpreter does not put results in dictionaries.
        assert len(interpretations) == 0

    def test_interpret_results_with_context(self):
        """Enter_ret_results() into the context parameter."""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=2)

        # MockiInterpreter does not use context, but should not throw out an anomaly.
        interpretations = interpreter.interpret_results(results, context="cancer research")
        assert "GO_Biological_Process" in interpretations

    def test_interpret_results_include_term_summaries(self):
        """Include_term_summaries=True generates a summary of entries."""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=3, include_genes_col=True)

        interpretations = interpreter.interpret_results(
            results, include_term_summaries=True
        )

        # It should be read in its entirety.
        assert "GO_Biological_Process" in interpretations
        # Summary of entries to be included
        assert "GO_Biological_Process_term_summaries" in interpretations
        term_summaries = interpretations["GO_Biological_Process_term_summaries"]
        assert isinstance(term_summaries, dict)
        assert len(term_summaries) == 3
        # Verify that each entry summary contains gene-set information
        for term_name, summary in term_summaries.items():
            assert "genes" in summary.lower()

    def test_interpret_results_include_term_summaries_empty(self):
        """Include_term_summaries=True but results are empty and do not result in an entry summary."""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_empty_results()

        interpretations = interpreter.interpret_results(
            results, include_term_summaries=True
        )

        assert "GO_Biological_Process_term_summaries" not in interpretations

    def test_generate_report_section(self):
        """Generates HTML."""
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
        """......................"""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_empty_results()

        html = interpreter.generate_report_section(results)

        assert "<div" in html
        assert "AI-Powered Interpretation" in html
        # Should not contain any database reading blocks
        assert "GO_Biological_Process" not in html

    def test_generate_report_section_with_context(self):
        """Generate_report_section() to enter parameters."""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=2)

        html = interpreter.generate_report_section(results, context="diabetes study")
        assert "AI-Powered Interpretation" in html

    def test_generate_report_section_excludes_term_summaries(self):
        """Generate_report_section() does not include the summary of entries (read as a whole only)."""
        from allenricher.ai.interpreter import AIInterpreter

        interpreter = AIInterpreter(backend="mock")
        results = _make_results(n_terms=2, include_genes_col=True)

        html = interpreter.generate_report_section(results)

        # _term_summaries
        assert "_term_summaries" not in html

    def test_backend_name_stored(self):
        """The right memory."""
        from allenricher.ai.interpreter import AIInterpreter

        for backend in ["mock", "openai", "claude", "ollama", "deepseek", "glm", "minimax"]:
            kwargs = {}
            if backend == "minimax":
                kwargs["group_id"] = "test"
            interpreter = AIInterpreter(backend=backend, **kwargs)
            assert interpreter.backend_name == backend


# ===========================================================================
# 7. Factorial function testing
# ===========================================================================

class TestFactoryFunctions:
    """Tests the Create_interpreter() and get_available_backends() factory functions."""

    def test_create_interpreter_mock(self):
        """Create_interpreter("mock") returns the correct example."""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter("mock")
        assert interpreter.backend_name == "mock"

    def test_create_interpreter_openai(self):
        """returns correct examples."""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter("openai", api_key="test-key")
        assert interpreter.backend_name == "openai"
        assert interpreter.interpreter.api_key == "test-key"

    def test_create_interpreter_claude(self):
        """returns correct examples."""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter("claude", api_key="test-key")
        assert interpreter.backend_name == "claude"
        assert interpreter.interpreter.api_key == "test-key"

    def test_create_interpreter_ollama(self):
        """Create_interpreter("olama") returns the correct example."""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter("ollama", model="mistral")
        assert interpreter.backend_name == "ollama"
        assert interpreter.interpreter.model == "mistral"

    def test_create_interpreter_deepseek(self):
        """Create_interpret ("deepseek") returns the correct example."""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter("deepseek", api_key="test-key")
        assert interpreter.backend_name == "deepseek"
        assert interpreter.interpreter.api_key == "test-key"

    def test_create_interpreter_glm(self):
        """returns correct examples."""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter("glm", api_key="test-key")
        assert interpreter.backend_name == "glm"
        assert interpreter.interpreter.api_key == "test-key"

    def test_create_interpreter_minimax(self):
        """returns correct examples."""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter(
            "minimax", api_key="test-key", group_id="test-group"
        )
        assert interpreter.backend_name == "minimax"
        assert interpreter.interpreter.api_key == "test-key"
        assert interpreter.interpreter.group_id == "test-group"

    def test_create_interpreter_invalid(self):
        """"Create_interpreter" throws down ValueError."""
        from allenricher.ai.interpreter import create_interpreter

        with pytest.raises(ValueError, match="Unknown backend"):
            create_interpreter("invalid")

    def test_create_interpreter_default_backend(self):
        """Create_interpreter() uses the mock backend by default."""
        from allenricher.ai.interpreter import create_interpreter

        interpreter = create_interpreter()
        assert interpreter.backend_name == "mock"

    def test_get_available_backends(self):
        """get_available_backends() returns the list."""
        from allenricher.ai.interpreter import get_available_backends

        backends = get_available_backends()
        assert isinstance(backends, list)
        assert len(backends) == 7
        expected = {"openai", "claude", "deepseek", "glm", "minimax", "ollama", "mock"}
        assert set(backends) == expected


# ===========================================================================
# 8. Consistency testing for the abstract base-based interface
# ===========================================================================

class TestAbstractBaseInterface:
    """Validates that all backends are performed following a unified abstract interface."""

    @pytest.fixture(params=["mock", "ollama"])
    def concrete_interpreter(self, request):
        """Creates specific examples of interpretors that do not require API key."""
        from allenricher.ai.interpreter import AIInterpreter

        return AIInterpreter(backend=request.param)

    def test_interpret_signature(self, concrete_interpreter):
        """Method accepts the Dit [str, DataFrame] parameter."""
        results = _make_results(n_terms=2)
        # There's no need to throw out the anomaly.
        result = concrete_interpreter.interpreter.interpret(results)
        assert isinstance(result, dict)

    def test_summarize_term_signature(self, concrete_interpreter):
        """The summarize_term() method accepts the st and list[str] parameters."""
        result = concrete_interpreter.interpreter.summarize_term("test_term", ["GENE1"])
        assert isinstance(result, str)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
