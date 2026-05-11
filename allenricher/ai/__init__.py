"""AI模块 - 提供AI驱动的富集分析结果解读，支持OpenAI、Claude、Ollama等后端"""
from allenricher.ai.interpreter import AIInterpreter, create_interpreter  # AI解读器和工厂函数

__all__ = ["AIInterpreter", "create_interpreter"]
