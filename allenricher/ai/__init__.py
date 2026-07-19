"""AI module - Provides an AI-driven interpretation of the enrichment analysis results in support of backends such as OpenAI, Claude, DeepSeek, GLM, MiniMax, Ollama"""
from allenricher.ai.interpreter import (
    AIInterpreter,
    create_interpreter,
    get_available_backends,
    OpenAIInterpreter,
    ClaudeInterpreter,
    DeepSeekInterpreter,
    GLMInterpreter,
    MiniMaxInterpreter,
    OllamaInterpreter,
    MockInterpreter,
    build_structured_evidence,
    build_interpretation_prompt,
    validate_interpretation,
)

__all__ = [
    "AIInterpreter",
    "create_interpreter",
    "get_available_backends",
    "OpenAIInterpreter",
    "ClaudeInterpreter",
    "DeepSeekInterpreter",
    "GLMInterpreter",
    "MiniMaxInterpreter",
    "OllamaInterpreter",
    "MockInterpreter",
    "build_structured_evidence",
    "build_interpretation_prompt",
    "validate_interpretation",
]
