"""下载工具函数模块

提供文件完整性校验、进度显示、格式化等工具函数。
"""
import gzip
import hashlib
import time
from pathlib import Path


def verify_gzip_integrity(filepath: Path, sample_lines: int = 100) -> tuple:
    """验证 gzip 文件完整性

    采用采样策略：验证文件头、中间段、尾部，避免大文件全量解压。

    Args:
        filepath: gzip 文件路径
        sample_lines: 采样验证的行数（0 = 全量验证，较慢）

    Returns:
        (是否有效, 错误信息)
    """
    try:
        with gzip.open(filepath, 'rt', encoding='utf-8', errors='ignore') as f:
            if sample_lines > 0:
                # 采样验证：从头到尾完整遍历（不全量保存），检查能否读到文件末尾
                # 同时对前 N 行做格式检查
                checked = 0
                for i, line in enumerate(f):
                    checked += 1
                    # 只对前 sample_lines 行做格式检查
                    if i < sample_lines:
                        stripped = line.strip()
                        if stripped and '\t' not in stripped and not stripped.startswith('#'):
                            return False, f"Invalid line format at line {i}: {stripped[:80]}"
            else:
                # 全量验证（慢）
                for _ in f:
                    pass

        return True, "OK"

    except gzip.BadGzipFile as e:
        return False, f"Bad gzip file: {e}"
    except EOFError as e:
        return False, f"Unexpected EOF (truncated file): {e}"
    except Exception as e:
        return False, f"Verification error: {e}"


def calculate_file_hash(filepath: Path, algorithm: str = "md5") -> str:
    """计算文件哈希值

    Args:
        filepath: 文件路径
        algorithm: 哈希算法 ('md5' 或 'sha256')

    Returns:
        十六进制哈希字符串
    """
    h = hashlib.md5() if algorithm == "md5" else hashlib.sha256()
    with open(filepath, 'rb') as f:
        for chunk in iter(lambda: f.read(65536), b''):
            h.update(chunk)
    return h.hexdigest()


def format_size(size_bytes: int) -> str:
    """格式化文件大小显示

    Args:
        size_bytes: 字节数

    Returns:
        人类可读的大小字符串，如 "1.3 GB"
    """
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024
    return f"{size_bytes:.1f} PB"


def format_speed(bytes_per_sec: float) -> str:
    """格式化下载速度

    Args:
        bytes_per_sec: 每秒字节数

    Returns:
        人类可读的速度字符串，如 "5.2 MB/s"
    """
    if bytes_per_sec <= 0:
        return "---"
    return f"{format_size(bytes_per_sec)}/s"


def format_duration(seconds: float) -> str:
    """格式化时长

    Args:
        seconds: 秒数

    Returns:
        人类可读的时长字符串，如 "2m 30s"
    """
    if seconds < 60:
        return f"{seconds:.0f}s"
    elif seconds < 3600:
        m = int(seconds // 60)
        s = int(seconds % 60)
        return f"{m}m {s}s"
    else:
        h = int(seconds // 3600)
        m = int((seconds % 3600) // 60)
        return f"{h}h {m}m"


class SimpleProgressBar:
    """轻量级终端进度条（零外部依赖）

    显示格式: desc |████████░░░░░░░░░░░░| 45.2% 1.2 GB/s 3m 20s
    """

    def __init__(self, total: int, desc: str = "", width: int = 40):
        """
        Args:
            total: 总字节数
            desc: 描述文本
            width: 进度条字符宽度
        """
        self.total = total
        self.desc = desc
        self.width = width
        self.n = 0
        self.start_time = time.time()
        self._last_update = 0  # 上次更新的字节数

    def update(self, n: int = 1):
        """更新进度

        Args:
            n: 本次新增字节数
        """
        self.n += n
        # 每 0.5 秒或完成时刷新一次
        now = time.time()
        if now - self._last_update >= 0.5 or self.n >= self.total:
            self._last_update = now
            self._render()

    def _render(self):
        if self.total <= 0:
            return
        percent = min(self.n / self.total, 1.0)
        filled = int(self.width * percent)
        bar = '█' * filled + '░' * (self.width - filled)

        # 速度和剩余时间
        elapsed = time.time() - self.start_time
        speed = self.n / elapsed if elapsed > 0 else 0
        remaining = (self.total - self.n) / speed if speed > 0 else 0

        line = (
            f"\r{self.desc} |{bar}| {percent*100:5.1f}% "
            f"{format_speed(speed)} {format_duration(remaining)}"
        )
        print(line, end='', flush=True)

    def close(self):
        """关闭进度条（换行）"""
        if self.total > 0:
            # 确保最终状态渲染
            self.n = self.total
            self._render()
        print()
