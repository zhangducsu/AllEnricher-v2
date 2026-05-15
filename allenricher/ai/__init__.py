"""AI模块 - 提供AI驱动的富集分析结果解读，支持OpenAI、Claude、DeepSeek、GLM、MiniMax、Ollama等后端"""
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
    MockInterpreter
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
    "MockInterpreter"
]
