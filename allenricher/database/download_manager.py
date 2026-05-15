"""下载管理器模块

提供高级下载功能：多线程分片下载、镜像源自动切换、
gzip 完整性校验、断点续传、实时进度显示。

对应 v1: update_GOdb / update_ReactomeDB 的下载部分。
"""
import shutil
import threading
import urllib.request
import urllib.error
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Optional

from .download_utils import (
    SimpleProgressBar,
    format_size,
    verify_gzip_integrity,
)
from .mirrors import MirrorSource


class DownloadManager:
    """高级下载管理器

    特性：
    - 多线程分片下载大文件（>100MB）
    - 镜像源自动切换（主源失败 → 备用源）
    - gzip 完整性校验（采样验证，不全量解压）
    - 断点续传（单线程 Range + 多线程分片级）
    - 实时进度条（速度 + 剩余时间）

    Usage::

        dm = DownloadManager(max_workers=4)
        dm.download_file(url, dest_path)
        dm.download_with_mirror_fallback(mirrors, filename, dest_path)
    """

    # 启用多线程的文件大小阈值 (100 MB)
    MULTI_THREAD_THRESHOLD = 100 * 1024 * 1024
    # 多线程分片大小 (10 MB)
    CHUNK_SIZE = 10 * 1024 * 1024
    # 单线程读取块大小 (64 KB)
    READ_BLOCK = 64 * 1024
    # 连接超时 (秒)
    CONNECT_TIMEOUT = 30
    # 读取超时 (秒)
    READ_TIMEOUT = 600
    # User-Agent
    UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"

    def __init__(
        self,
        root_dir: str = "./database",
        overwrite: bool = False,
        max_workers: int = 4,
        use_multi_thread: bool = True,
        verify_integrity: bool = True,
        show_progress: bool = True,
    ):
        """
        Args:
            root_dir: 数据库根目录
            overwrite: 是否覆盖已存在文件
            max_workers: 多线程下载线程数
            use_multi_thread: 是否对大文件启用多线程
            verify_integrity: 下载后是否验证 gzip 完整性
            show_progress: 是否显示进度条
        """
        self.root_dir = Path(root_dir)
        self.overwrite = overwrite
        self.max_workers = max_workers
        self.use_multi_thread = use_multi_thread
        self.verify_integrity = verify_integrity
        self.show_progress = show_progress

    # ================================================================
    # 公共接口
    # ================================================================
    def download_file(
        self,
        url: str,
        dest_path: Path,
        expected_size: Optional[int] = None,
        desc: Optional[str] = None,
    ) -> Path:
        """下载单个文件（自动选择单线程/多线程）

        Args:
            url: 下载 URL
            dest_path: 目标路径
            expected_size: 预期文件大小（字节），用于进度条
            desc: 进度条描述文字

        Returns:
            下载完成的文件路径

        Raises:
            RuntimeError: 下载或校验失败
        """
        dest_path = Path(dest_path)

        # ---- 已存在文件处理 ----
        if dest_path.exists() and not self.overwrite:
            if self.verify_integrity and dest_path.suffix == '.gz':
                valid, msg = verify_gzip_integrity(dest_path)
                if valid:
                    print(f"|--- 已存在且有效，跳过: {dest_path.name}")
                    return dest_path
                else:
                    print(f"|--- 已存在但已损坏，将重新下载: {dest_path.name} ({msg})")
            else:
                print(f"|--- 已存在，跳过: {dest_path.name}")
                return dest_path

        # ---- 获取远程文件大小 ----
        if expected_size is None:
            expected_size = self._head_content_length(url)

        # ---- 选择策略 ----
        use_mt = (
            self.use_multi_thread
            and expected_size is not None
            and expected_size > self.MULTI_THREAD_THRESHOLD
        )
        if use_mt:
            return self._download_multi_thread(url, dest_path, expected_size, desc)
        else:
            return self._download_single_thread(url, dest_path, expected_size, desc)

    def download_with_mirror_fallback(
        self,
        mirrors: List[MirrorSource],
        filename: str,
        dest_path: Path,
        desc: Optional[str] = None,
    ) -> Path:
        """从多个镜像源依次尝试下载，自动切换

        Args:
            mirrors: 镜像源列表（已按优先级排序）
            filename: 相对于 base_url 的文件名
            dest_path: 目标路径
            desc: 进度条描述

        Returns:
            下载完成的文件路径

        Raises:
            RuntimeError: 所有镜像源均失败
        """
        last_error: Optional[Exception] = None

        for mirror in mirrors:
            url = f"{mirror.base_url.rstrip('/')}/{filename}"
            try:
                print(f"|--- 尝试镜像: {mirror.name} ({mirror.region})")
                return self.download_file(url, dest_path, desc=desc or filename)
            except Exception as e:
                last_error = e
                print(f"    ✗ 失败: {e}")
                # 失败后标记需要覆盖（避免残留文件被跳过）
                self.overwrite = True
                continue

        raise RuntimeError(
            f"所有镜像源均失败。最后错误: {last_error}"
        )

    # ================================================================
    # 内部方法
    # ================================================================
    def _head_content_length(self, url: str) -> Optional[int]:
        """HEAD 请求获取 Content-Length"""
        try:
            req = urllib.request.Request(url, method='HEAD')
            req.add_header("User-Agent", self.UA)
            with urllib.request.urlopen(req, timeout=self.CONNECT_TIMEOUT) as resp:
                length = resp.headers.get('Content-Length')
                return int(length) if length else None
        except Exception:
            return None

    # ---- 单线程下载（带断点续传 + 进度条） ----
    def _download_single_thread(
        self,
        url: str,
        dest_path: Path,
        expected_size: Optional[int],
        desc: Optional[str],
    ) -> Path:
        dest_path = Path(dest_path)
        dest_path.parent.mkdir(parents=True, exist_ok=True)

        # 断点续传
        downloaded = 0
        if dest_path.exists():
            downloaded = dest_path.stat().st_size
            if downloaded > 0:
                print(f"    断点续传: 已有 {format_size(downloaded)}")

        req = urllib.request.Request(url)
        req.add_header("User-Agent", self.UA)
        if downloaded > 0:
            req.add_header("Range", f"bytes={downloaded}-")

        mode = "ab" if downloaded > 0 else "wb"
        total = downloaded
        label = (desc or dest_path.name)[:30]

        progress = None
        if self.show_progress and expected_size:
            progress = SimpleProgressBar(expected_size, desc=label)
            progress.update(downloaded)

        try:
            with urllib.request.urlopen(
                req, timeout=self.READ_TIMEOUT
            ) as resp:
                # 服务器不支持断点续传 → 全量重下
                if downloaded > 0 and resp.status != 206:
                    print("    服务器不支持断点续传，从头下载")
                    downloaded = 0
                    total = 0
                    mode = "wb"
                    req2 = urllib.request.Request(url)
                    req2.add_header("User-Agent", self.UA)
                    resp = urllib.request.urlopen(
                        req2, timeout=self.READ_TIMEOUT
                    )

                with open(dest_path, mode) as f:
                    while True:
                        chunk = resp.read(self.READ_BLOCK)
                        if not chunk:
                            break
                        f.write(chunk)
                        total += len(chunk)
                        if progress:
                            progress.update(len(chunk))

            if progress:
                progress.close()

            # 完整性校验
            self._post_verify(dest_path)
            return dest_path

        except Exception as e:
            if progress:
                progress.close()
            raise RuntimeError(f"下载失败 {url}: {e}") from e

    # ---- 多线程分片下载 ----
    def _download_multi_thread(
        self,
        url: str,
        dest_path: Path,
        expected_size: int,
        desc: Optional[str],
    ) -> Path:
        dest_path = Path(dest_path)
        dest_path.parent.mkdir(parents=True, exist_ok=True)

        num_chunks = (expected_size + self.CHUNK_SIZE - 1) // self.CHUNK_SIZE
        label = (desc or dest_path.name)[:30]
        print(
            f"    多线程下载: {format_size(expected_size)}, "
            f"{num_chunks} 分片, {self.max_workers} 线程"
        )

        # 临时分片文件
        parts = [
            dest_path.parent / f"{dest_path.name}.part{i}"
            for i in range(num_chunks)
        ]

        # 共享进度
        chunk_bytes = [0] * num_chunks
        lock = threading.Lock()

        def _download_chunk(idx: int) -> bool:
            """下载第 idx 个分片"""
            start = idx * self.CHUNK_SIZE
            end = min(start + self.CHUNK_SIZE - 1, expected_size - 1)
            part = parts[idx]

            # 已完成？
            if part.exists() and part.stat().st_size == end - start + 1:
                with lock:
                    chunk_bytes[idx] = part.stat().st_size
                return True

            # 断点
            offset = part.stat().st_size if part.exists() else 0
            actual_start = start + offset

            req = urllib.request.Request(url)
            req.add_header("User-Agent", self.UA)
            req.add_header("Range", f"bytes={actual_start}-{end}")

            try:
                with urllib.request.urlopen(
                    req, timeout=self.READ_TIMEOUT
                ) as resp:
                    m = "ab" if part.exists() else "wb"
                    with open(part, m) as f:
                        while True:
                            block = resp.read(self.READ_BLOCK)
                            if not block:
                                break
                            f.write(block)
                            with lock:
                                chunk_bytes[idx] += len(block)
                return True
            except Exception as exc:
                print(f"\n    分片 {idx} 失败: {exc}")
                return False

        # 线程池
        progress = (
            SimpleProgressBar(expected_size, desc=label)
            if self.show_progress
            else None
        )

        with ThreadPoolExecutor(max_workers=self.max_workers) as pool:
            futures = {pool.submit(_download_chunk, i): i for i in range(num_chunks)}
            for future in as_completed(futures):
                idx = futures[future]
                if not future.result():
                    # 清理已下载分片
                    for p in parts:
                        if p.exists():
                            p.unlink()
                    raise RuntimeError(f"分片 {idx} 下载失败，已清理临时文件")

                # 更新进度
                if progress:
                    with lock:
                        progress.n = sum(chunk_bytes)
                    progress.update(0)

        if progress:
            progress.close()

        # 合并分片
        print("    合并分片...")
        with open(dest_path, 'wb') as out:
            for part in parts:
                if part.exists():
                    with open(part, 'rb') as inp:
                        shutil.copyfileobj(inp, out)
                    part.unlink()

        # 完整性校验
        self._post_verify(dest_path)
        return dest_path

    # ---- 下载后校验 ----
    def _post_verify(self, filepath: Path):
        """下载后 gzip 完整性校验"""
        if not self.verify_integrity:
            return
        if filepath.suffix != '.gz':
            return

        print("    验证完整性...", end='', flush=True)
        valid, msg = verify_gzip_integrity(filepath)
        if valid:
            print(" ✓ 通过")
        else:
            # 删除损坏文件，让上层可以重试
            filepath.unlink(missing_ok=True)
            raise RuntimeError(
                f"文件完整性验证失败: {msg}，已删除损坏文件"
            )
