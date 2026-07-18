"""Provide verified single-stream and parallel file downloads with mirror fallback."""
import shutil
import threading
import time
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
    """Download files with integrity checks, retries, and optional parallel ranges."""

    # Enable multilined file size threshold (100 MB)
    MULTI_THREAD_THRESHOLD = 100 * 1024 * 1024
    # Multiline Split Size (10 MB)
    CHUNK_SIZE = 10 * 1024 * 1024
    # Line Read Block Size (64 KB)
    READ_BLOCK = 64 * 1024
    # Connection timed out (sec)
    CONNECT_TIMEOUT = 30
    # Read Timeout (sec)
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
root_dir: database root
overwrite: Whether to overwrite existing files
Max_workers: number of threads downloaded over multiple threads
use_multi_thread: Whether to enable multi-lined processes for large files
vereify_integrity: Whether gzip integrity is verified after downloading
Show_progress: Whether progress bar is displayed
        """
        self.root_dir = Path(root_dir)
        self.overwrite = overwrite
        self.max_workers = max_workers
        self.use_multi_thread = use_multi_thread
        self.verify_integrity = verify_integrity
        self.show_progress = show_progress

    # ================================================================
    # Public interface
    # ================================================================
    def download_file(
        self,
        url: str,
        dest_path: Path,
        expected_size: Optional[int] = None,
        desc: Optional[str] = None,
    ) -> Path:
        """Download one file with integrity verification."""
        dest_path = Path(dest_path)

        # ----File Processing Exists----
        if dest_path.exists() and not self.overwrite:
            if self.verify_integrity and dest_path.suffix == '.gz':
                valid, msg = verify_gzip_integrity(dest_path)
                if valid:
                    print(f"|---Existence and validity; skipping: {dest_path.name}")
                    return dest_path
                else:
                    print(f"|---Existing but damaged, will be re-downloaded: {dest_path.name} ({msg})")
            else:
                print(f"|---Existence, Skipping: {dest_path.name}")
                return dest_path

        # ----Retrieving remote file size----
        if expected_size is None:
            expected_size = self._head_content_length(url)

        # ----Select Policy----
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
        """Try ordered mirrors until one verified download succeeds."""
        last_error: Optional[Exception] = None

        for mirror in mirrors:
            url = f"{mirror.base_url.rstrip('/')}/{filename}"
            try:
                print(f"|---Try mirror image: {mirror.name} ({mirror.region})")
                return self.download_file(url, dest_path, desc=desc or filename)
            except Exception as e:
                last_error = e
                print(f"    [FAILED] {e}")
                # The markup needs to be overwritten (avoid the residue file is skipped)
                self.overwrite = True
                continue

        raise RuntimeError(
            f"All mirror sources failed. Last error: {last_error}"
        )

    # ================================================================
    # Internal methodology
    # ================================================================
    def _head_content_length(self, url: str) -> Optional[int]:
        """Read Content-Length from an HTTP HEAD response when available."""
        try:
            req = urllib.request.Request(url, method='HEAD')
            req.add_header("User-Agent", self.UA)
            with urllib.request.urlopen(req, timeout=self.CONNECT_TIMEOUT) as resp:
                length = resp.headers.get('Content-Length')
                return int(length) if length else None
        except Exception:
            return None

    # ----Single-line download (with breakpoints + progress bar)----
    def _download_single_thread(
        self,
        url: str,
        dest_path: Path,
        expected_size: Optional[int],
        desc: Optional[str],
    ) -> Path:
        dest_path = Path(dest_path)
        dest_path.parent.mkdir(parents=True, exist_ok=True)

        # Pass it on.
        downloaded = 0
        if dest_path.exists():
            downloaded = dest_path.stat().st_size
            if downloaded > 0:
                print(f"Breakpoints: Already{format_size(downloaded)}")

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
                # The server does not support stop-point transmissions.
                if downloaded > 0 and resp.status != 206:
                    print("The server does not support intermittent transmissions, download from the header")
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

            # Complete Validation
            self._post_verify(dest_path)
            return dest_path

        except Exception as e:
            if progress:
                progress.close()
            raise RuntimeError(f"Download Failed{url}: {e}") from e

    # ----Download multilined fractions----
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
            f"Multiline download: {format_size(expected_size)},"
            f"{num_chunks}The film is a piece of the paper.{self.max_workers}Thread"
        )

        # Temporary Dissect Files
        parts = [
            dest_path.parent / f"{dest_path.name}.part{i}"
            for i in range(num_chunks)
        ]

        # Share Progress
        chunk_bytes = [0] * num_chunks
        lock = threading.Lock()

        def _download_chunk(idx: int) -> bool:
            """Download one byte range for a parallel transfer."""
            start = idx * self.CHUNK_SIZE
            end = min(start + self.CHUNK_SIZE - 1, expected_size - 1)
            part = parts[idx]
            expected_chunk_size = end - start + 1

            # Done?
            if part.exists() and part.stat().st_size == expected_chunk_size:
                with lock:
                    chunk_bytes[idx] = part.stat().st_size
                return True

            if part.exists() and part.stat().st_size > expected_chunk_size:
                part.unlink()

            last_error: Optional[Exception] = None
            for attempt in range(1, 4):
                offset = part.stat().st_size if part.exists() else 0
                actual_start = start + offset
                with lock:
                    chunk_bytes[idx] = offset
                req = urllib.request.Request(url)
                req.add_header("User-Agent", self.UA)
                req.add_header("Range", f"bytes={actual_start}-{end}")
                try:
                    with urllib.request.urlopen(
                        req, timeout=self.READ_TIMEOUT
                    ) as resp:
                        if getattr(resp, "status", 206) != 206:
                            raise RuntimeError("The server did not press the Range returns the partition")
                        with open(part, "ab" if offset else "wb") as f:
                            while True:
                                block = resp.read(self.READ_BLOCK)
                                if not block:
                                    break
                                f.write(block)
                                with lock:
                                    chunk_bytes[idx] += len(block)
                        if part.stat().st_size == expected_chunk_size:
                            return True
                        raise RuntimeError(
                            f"Split size does not match: {part.stat().st_size}!= {expected_chunk_size}"
                        )
                except Exception as exc:
                    last_error = exc
                if attempt < 3:
                    time.sleep(2 ** (attempt - 1))
            print(f"\nSplit{idx}Failed after retrying 3 times: {last_error}")
            return False

        # Thread Pool
        progress = (
            SimpleProgressBar(expected_size, desc=label)
            if self.show_progress
            else None
        )

        failed_chunks = []
        with ThreadPoolExecutor(max_workers=self.max_workers) as pool:
            futures = {pool.submit(_download_chunk, i): i for i in range(num_chunks)}
            for future in as_completed(futures):
                idx = futures[future]
                if not future.result():
                    failed_chunks.append(idx)

                # Update Progress
                if progress:
                    with lock:
                        progress.n = sum(chunk_bytes)
                    progress.update(0)

        if failed_chunks:
            if progress:
                progress.close()
            raise RuntimeError(
                f"Split{sorted(failed_chunks)}Failed to download; finished fragments are reserved for renewal"
            )

        if progress:
            progress.close()

        # Merge Sub-scopic
        print("Merge Disperses...")
        with open(dest_path, 'wb') as out:
            for part in parts:
                if part.exists():
                    with open(part, 'rb') as inp:
                        shutil.copyfileobj(inp, out)
                    part.unlink()

        # Complete Validation
        self._post_verify(dest_path)
        return dest_path

    # ----Check after download----
    def _post_verify(self, filepath: Path):
        """Verify downloaded size and optional gzip integrity."""
        if not self.verify_integrity:
            return
        if filepath.suffix != '.gz':
            return

        print("Validate Complete...", end='', flush=True)
        valid, msg = verify_gzip_integrity(filepath)
        if valid:
            print("[OK] Pass")
        else:
            # Remove corrupt file so the upper level can be retries
            filepath.unlink(missing_ok=True)
            raise RuntimeError(
                f"Could not close temporary folder: %s{msg}, Deleted corrupt file"
            )
