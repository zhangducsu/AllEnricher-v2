"""Utility functions for verified downloads and terminal progress reporting."""
import gzip
import hashlib
import time
from pathlib import Path


def verify_gzip_integrity(filepath: Path, sample_lines: int = 100) -> tuple:
    """Return whether a gzip file can be read to completion."""
    try:
        with gzip.open(filepath, 'rt', encoding='utf-8', errors='ignore') as f:
            if sample_lines > 0:
                # Sample validation: Full from beginning to end (incomplete preservation), check if you can read to end of file
                # Check the front N lines with formatting
                checked = 0
                for i, line in enumerate(f):
                    checked += 1
                    # Format checks only for former sample_lines
                    if i < sample_lines:
                        stripped = line.strip()
                        if stripped and '\t' not in stripped and not stripped.startswith('#'):
                            return False, f"Invalid line format at line {i}: {stripped[:80]}"
            else:
                # Full-scale validation (slow)
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
    """Calculate a cryptographic digest for one file."""
    h = hashlib.md5() if algorithm == "md5" else hashlib.sha256()
    with open(filepath, 'rb') as f:
        for chunk in iter(lambda: f.read(65536), b''):
            h.update(chunk)
    return h.hexdigest()


def format_size(size_bytes: int) -> str:
    """Format a byte count for terminal display."""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024
    return f"{size_bytes:.1f} PB"


def format_speed(bytes_per_sec: float) -> str:
    """Format a transfer rate for terminal display."""
    if bytes_per_sec <= 0:
        return "---"
    return f"{format_size(bytes_per_sec)}/s"


def format_duration(seconds: float) -> str:
    """Format elapsed seconds for terminal display."""
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
    """Render dependency-free download progress in a terminal."""

    def __init__(self, total: int, desc: str = "", width: int = 40):
        """
        Args:
total: Total bytes
desc: Description of Text
width: Progress bar character width
        """
        self.total = total
        self.desc = desc
        self.width = width
        self.n = 0
        self.start_time = time.time()
        self._last_update = 0  # Number of bytes last updated

    def update(self, n: int = 1):
        """Update the displayed progress state."""
        self.n += n
        # Update every 0.5 seconds or when finished
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

        # Speed and time remaining
        elapsed = time.time() - self.start_time
        speed = self.n / elapsed if elapsed > 0 else 0
        remaining = (self.total - self.n) / speed if speed > 0 else 0

        line = (
            f"\r{self.desc} |{bar}| {percent*100:5.1f}% "
            f"{format_speed(speed)} {format_duration(remaining)}"
        )
        print(line, end='', flush=True)

    def close(self):
        """Finish the progress line."""
        if self.total > 0:
            # Ensure final state rendering
            self.n = self.total
            self._render()
        print()
