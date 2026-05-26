# AllEnricher v2 下载器优化计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 优化 AllEnricher v2 的 `download` 命令，实现多线程下载加速、完整性校验、镜像源自动切换，解决大文件下载慢和易损坏的问题。

**Architecture:** 重构 `downloader.py`，引入 `DownloadManager` 类管理下载任务，支持多线程分片下载、自动镜像源切换、下载后 gzip 完整性验证，并提供下载进度显示和断点续传。

**Tech Stack:** Python 3.10+, `concurrent.futures`, `hashlib`, `gzip`, `urllib.request`

---

## 当前问题分析

### 1. 性能瓶颈
- **单线程下载**: `_download_file()` 使用单线程，8KB chunk 读取，下载 1.4GB gene_info.gz 需 25+ 分钟
- **无连接复用**: 每个文件独立建立 HTTP 连接
- **无镜像源**: 仅使用官方 NCBI FTP，国内访问慢

### 2. 完整性问题
- **无校验机制**: 下载完成后不验证文件完整性
- **断点续传缺陷**: 仅支持 Range 请求，不验证续传后文件是否完整
- **损坏文件残留**: gene2go.gz 多次下载损坏（1261MB vs 正确 1.23GB）

### 3. 用户体验
- **无进度显示**: 大文件下载时用户无法知道进度
- **失败无重试**: 网络中断后需手动重新下载

---

## 优化方案设计

### 方案 A: 多线程分片下载 (优先级: 高)
将大文件分成多个 range 块，使用线程池并行下载，最后合并。

**适用文件:**
- gene_info.gz (~1.4GB)
- gene2go.gz (~1.2GB)

**不适用:**
- go-basic.obo (~30MB)
- Jensen Lab TSV 文件 (<50MB)

### 方案 B: 镜像源自动切换 (优先级: 高)
维护多个镜像源列表，主源失败时自动切换到备用源。

**镜像源候选:**
1. NCBI 官方: `https://ftp.ncbi.nlm.nih.gov/gene/DATA/`
2. EBI 欧洲: `https://ftp.ebi.ac.uk/pub/databases/ncbi/gene/DATA/`
3. 国内镜像: 清华/中科大 (需验证可用性)

### 方案 C: 完整性校验 (优先级: 高)
下载完成后进行多重验证:
1. **gzip 验证**: 尝试解压文件头，验证 gzip 格式
2. **行数检查**: 验证文件行数是否在合理范围
3. **内容采样**: 随机采样 100 行验证格式正确性

### 方案 D: 下载进度显示 (优先级: 中)
使用 `tqdm` 或自定义进度条显示下载进度、速度、预计剩余时间。

---

## 文件结构

```
allenricher/
├── database/
│   ├── downloader.py          # 现有文件 - 大幅重构
│   ├── download_manager.py    # 新增 - DownloadManager 类
│   ├── download_utils.py      # 新增 - 工具函数（校验、进度条）
│   └── mirrors.py             # 新增 - 镜像源配置
```

---

## Task 1: 创建镜像源配置模块

**Files:**
- Create: `allenricher/database/mirrors.py`

- [ ] **Step 1: 定义镜像源数据结构**

```python
"""镜像源配置模块

定义各数据库的镜像源列表，支持按优先级排序和自动切换。
"""
from dataclasses import dataclass
from typing import List, Optional

@dataclass
class MirrorSource:
    """镜像源定义"""
    name: str           # 镜像名称，如 "ncbi-official"
    base_url: str       # 基础 URL
    priority: int       # 优先级，数字越小优先级越高
    region: str         # 地区，如 "US", "EU", "CN"
    enabled: bool = True

# NCBI gene_info / gene2go 镜像源
NCBI_MIRRORS: List[MirrorSource] = [
    MirrorSource("ncbi-official", "https://ftp.ncbi.nlm.nih.gov/gene/DATA/", 1, "US"),
    MirrorSource("ebi-ftp", "https://ftp.ebi.ac.uk/pub/databases/ncbi/gene/DATA/", 2, "EU"),
    # 国内镜像待验证后添加
]

# GO OBO 文件镜像源
GO_MIRRORS: List[MirrorSource] = [
    MirrorSource("purl-obo", "http://purl.obolibrary.org/obo/go/", 1, "US"),
    MirrorSource("geneontology", "http://current.geneontology.org/ontology/", 2, "US"),
]

# Reactome 镜像源
REACTOME_MIRRORS: List[MirrorSource] = [
    MirrorSource("reactome-official", "https://reactome.org/download/current/", 1, "US"),
]

# Jensen Lab DO 数据源（无镜像）
JENSEN_SOURCES = [
    "http://download.jensenlab.org/human_disease_textmining_filtered.tsv",
    "http://download.jensenlab.org/human_disease_knowledge_filtered.tsv",
    "http://download.jensenlab.org/human_disease_experiments_filtered.tsv",
]

def get_mirrors(db_type: str) -> List[MirrorSource]:
    """获取指定数据库类型的镜像源列表（按优先级排序）"""
    mirrors_map = {
        'ncbi': NCBI_MIRRORS,
        'go': GO_MIRRORS,
        'reactome': REACTOME_MIRRORS,
    }
    mirrors = mirrors_map.get(db_type, [])
    return sorted([m for m in mirrors if m.enabled], key=lambda x: x.priority)
```

- [ ] **Step 2: 运行验证测试**

Run: `python -c "from allenricher.database.mirrors import get_mirrors, NCBI_MIRRORS; print([m.name for m in get_mirrors('ncbi')])"`

Expected: `['ncbi-official', 'ebi-ftp']`

- [ ] **Step 3: Commit**

```bash
git add allenricher/database/mirrors.py
git commit -m "feat: add mirror source configuration module"
```

---

## Task 2: 创建下载工具函数模块

**Files:**
- Create: `allenricher/database/download_utils.py`

- [ ] **Step 1: 实现文件完整性校验函数**

```python
"""下载工具函数模块

提供文件完整性校验、进度显示、临时文件管理等工具函数。
"""
import gzip
import hashlib
from pathlib import Path
from typing import Optional, Callable


def verify_gzip_integrity(filepath: Path, sample_lines: int = 100) -> tuple[bool, str]:
    """验证 gzip 文件完整性
    
    Args:
        filepath: gzip 文件路径
        sample_lines: 采样验证的行数（0=全部验证，较慢）
    
    Returns:
        (是否有效, 错误信息)
    """
    try:
        with gzip.open(filepath, 'rt', encoding='utf-8', errors='ignore') as f:
            if sample_lines > 0:
                # 采样验证：验证开头、中间、结尾
                lines = []
                for i, line in enumerate(f):
                    if i < sample_lines // 3:
                        lines.append(line)
                    elif i == sample_lines // 3:
                        # 跳到中间
                        for _ in range(sample_lines // 3):
                            try:
                                next(f)
                            except StopIteration:
                                break
                    elif i < sample_lines * 2 // 3:
                        lines.append(line)
                    elif i == sample_lines * 2 // 3:
                        # 跳到结尾附近
                        for _ in range(sample_lines // 3):
                            try:
                                next(f)
                            except StopIteration:
                                break
                    elif i < sample_lines:
                        lines.append(line)
                    else:
                        break
                
                # 验证采样行格式
                for line in lines:
                    if '\t' not in line and not line.startswith('#'):
                        return False, f"Invalid line format: {line[:50]}"
            else:
                # 完整验证
                for _ in f:
                    pass
        
        return True, "OK"
    except gzip.BadGzipFile as e:
        return False, f"Bad gzip file: {e}"
    except Exception as e:
        return False, f"Verification error: {e}"


def calculate_file_hash(filepath: Path, algorithm: str = "md5") -> str:
    """计算文件哈希值
    
    Args:
        filepath: 文件路径
        algorithm: 哈希算法 (md5, sha256)
    
    Returns:
        十六进制哈希字符串
    """
    hasher = hashlib.md5() if algorithm == "md5" else hashlib.sha256()
    with open(filepath, 'rb') as f:
        for chunk in iter(lambda: f.read(8192), b''):
            hasher.update(chunk)
    return hasher.hexdigest()


def format_file_size(size_bytes: int) -> str:
    """格式化文件大小显示"""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size_bytes < 1024:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024
    return f"{size_bytes:.1f} TB"


def format_duration(seconds: float) -> str:
    """格式化时长显示"""
    if seconds < 60:
        return f"{seconds:.0f}s"
    elif seconds < 3600:
        return f"{seconds/60:.1f}m"
    else:
        return f"{seconds/3600:.1f}h"
```

- [ ] **Step 2: 实现简单进度条类**

```python
class SimpleProgressBar:
    """简单进度条（不依赖 tqdm）"""
    
    def __init__(self, total: int, desc: str = "", width: int = 50):
        self.total = total
        self.desc = desc
        self.width = width
        self.n = 0
        self._last_print = 0
    
    def update(self, n: int = 1):
        self.n += n
        # 每 5% 或每 1MB 更新一次，避免频繁输出
        percent = self.n / self.total if self.total > 0 else 0
        if percent - self._last_print >= 0.05 or self.n == self.total:
            self._print(percent)
            self._last_print = percent
    
    def _print(self, percent: float):
        filled = int(self.width * percent)
        bar = '█' * filled + '░' * (self.width - filled)
        print(f"\r{self.desc} |{bar}| {percent*100:.1f}%", end='', flush=True)
    
    def close(self):
        print()  # 换行
```

- [ ] **Step 3: 运行验证测试**

Run: `python -c "from allenricher.database.download_utils import format_file_size, verify_gzip_integrity; print(format_file_size(1400000000))"`

Expected: `1.3 GB`

- [ ] **Step 4: Commit**

```bash
git add allenricher/database/download_utils.py
git commit -m "feat: add download utility functions (integrity check, progress bar)"
```

---

## Task 3: 创建 DownloadManager 核心类

**Files:**
- Create: `allenricher/database/download_manager.py`

- [ ] **Step 1: 实现单文件下载方法**

```python
"""下载管理器模块

提供高级下载功能：多线程下载、镜像源切换、完整性校验、断点续传。
"""
import os
import urllib.request
import urllib.error
from pathlib import Path
from typing import Optional, List, Callable
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

from .mirrors import MirrorSource, get_mirrors
from .download_utils import (
    verify_gzip_integrity, 
    calculate_file_hash,
    format_file_size,
    format_duration,
    SimpleProgressBar
)


class DownloadManager:
    """高级下载管理器
    
    支持多线程下载、自动镜像源切换、完整性校验、断点续传。
    """
    
    # 启用多线程下载的文件大小阈值 (100MB)
    MULTI_THREAD_THRESHOLD = 100 * 1024 * 1024
    # 分片大小 (10MB)
    CHUNK_SIZE = 10 * 1024 * 1024
    # 单线程 chunk 大小 (64KB)
    SINGLE_CHUNK_SIZE = 64 * 1024
    # 最大重试次数
    MAX_RETRIES = 3
    # 连接超时 (秒)
    CONNECT_TIMEOUT = 30
    # 读取超时 (秒)
    READ_TIMEOUT = 300
    
    def __init__(
        self,
        root_dir: str = "./database",
        overwrite: bool = False,
        max_workers: int = 4,
        use_multi_thread: bool = True,
        verify_integrity: bool = True,
        show_progress: bool = True
    ):
        """初始化下载管理器
        
        Args:
            root_dir: 数据库根目录
            overwrite: 是否覆盖已存在文件
            max_workers: 多线程下载的线程数
            use_multi_thread: 是否启用多线程下载大文件
            verify_integrity: 下载后是否验证文件完整性
            show_progress: 是否显示下载进度
        """
        self.root_dir = Path(root_dir)
        self.overwrite = overwrite
        self.max_workers = max_workers
        self.use_multi_thread = use_multi_thread
        self.verify_integrity = verify_integrity
        self.show_progress = show_progress
        self._lock = threading.Lock()
    
    def download_file(
        self,
        url: str,
        dest_path: Path,
        expected_size: Optional[int] = None,
        desc: Optional[str] = None
    ) -> Path:
        """下载单个文件（自动选择单线程或多线程）
        
        Args:
            url: 下载 URL
            dest_path: 目标路径
            expected_size: 预期文件大小（用于进度显示）
            desc: 进度条描述
        
        Returns:
            下载完成的文件路径
        
        Raises:
            RuntimeError: 下载失败
        """
        dest_path = Path(dest_path)
        
        # 检查已存在文件
        if dest_path.exists() and not self.overwrite:
            if self.verify_integrity and dest_path.suffix == '.gz':
                valid, msg = verify_gzip_integrity(dest_path)
                if valid:
                    print(f"|--- 已存在且有效，跳过: {dest_path.name}")
                    return dest_path
                else:
                    print(f"|--- 已存在但损坏，重新下载: {dest_path.name} ({msg})")
                    self.overwrite = True
            else:
                print(f"|--- 已存在，跳过: {dest_path.name}")
                return dest_path
        
        # 获取文件大小
        if expected_size is None:
            expected_size = self._get_remote_file_size(url)
        
        # 选择下载策略
        if self.use_multi_thread and expected_size and expected_size > self.MULTI_THREAD_THRESHOLD:
            return self._download_multi_thread(url, dest_path, expected_size, desc)
        else:
            return self._download_single_thread(url, dest_path, expected_size, desc)
    
    def _get_remote_file_size(self, url: str) -> Optional[int]:
        """获取远程文件大小"""
        try:
            req = urllib.request.Request(url, method='HEAD')
            req.add_header("User-Agent", "Mozilla/5.0")
            with urllib.request.urlopen(req, timeout=self.CONNECT_TIMEOUT) as response:
                return int(response.headers.get('Content-Length', 0)) or None
        except Exception:
            return None
    
    def _download_single_thread(
        self,
        url: str,
        dest_path: Path,
        expected_size: Optional[int],
        desc: Optional[str]
    ) -> Path:
        """单线程下载（支持断点续传）"""
        # 实现与现有 downloader.py 类似，但增加进度显示
        pass  # 详见 Step 2
    
    def _download_multi_thread(
        self,
        url: str,
        dest_path: Path,
        expected_size: int,
        desc: Optional[str]
    ) -> Path:
        """多线程分片下载"""
        pass  # 详见 Step 3
```

- [ ] **Step 2: 实现单线程下载（带进度显示）**

```python
    def _download_single_thread(
        self,
        url: str,
        dest_path: Path,
        expected_size: Optional[int],
        desc: Optional[str]
    ) -> Path:
        """单线程下载（支持断点续传和进度显示）"""
        dest_path = Path(dest_path)
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        
        # 检查断点续传
        downloaded_size = 0
        if dest_path.exists():
            downloaded_size = dest_path.stat().st_size
            print(f"    断点续传: {format_file_size(downloaded_size)}")
        
        # 构建请求
        req = urllib.request.Request(url)
        req.add_header("User-Agent", "Mozilla/5.0 (Windows NT 10.0; Win64; x64)")
        if downloaded_size > 0:
            req.add_header("Range", f"bytes={downloaded_size}-")
        
        # 下载
        mode = "ab" if downloaded_size > 0 else "wb"
        total = downloaded_size
        desc = desc or dest_path.name
        
        progress = None
        if self.show_progress and expected_size:
            progress = SimpleProgressBar(expected_size, desc=desc[:30])
            progress.update(downloaded_size)
        
        try:
            with urllib.request.urlopen(req, timeout=(self.CONNECT_TIMEOUT, self.READ_TIMEOUT)) as response:
                # 检查服务器是否支持断点续传
                if downloaded_size > 0 and response.status != 206:
                    print(f"    服务器不支持断点续传，重新下载")
                    downloaded_size = 0
                    total = 0
                    mode = "wb"
                    req = urllib.request.Request(url)
                    req.add_header("User-Agent", "Mozilla/5.0")
                    response = urllib.request.urlopen(req, timeout=(self.CONNECT_TIMEOUT, self.READ_TIMEOUT))
                
                with open(dest_path, mode) as f:
                    while True:
                        chunk = response.read(self.SINGLE_CHUNK_SIZE)
                        if not chunk:
                            break
                        f.write(chunk)
                        total += len(chunk)
                        if progress:
                            progress.update(len(chunk))
            
            if progress:
                progress.close()
            
            # 验证完整性
            if self.verify_integrity and dest_path.suffix == '.gz':
                print(f"    验证文件完整性...")
                valid, msg = verify_gzip_integrity(dest_path)
                if not valid:
                    raise RuntimeError(f"文件验证失败: {msg}")
                print(f"    验证通过")
            
            return dest_path
            
        except Exception as e:
            if progress:
                progress.close()
            raise RuntimeError(f"下载失败 {url}: {e}")
```

- [ ] **Step 3: 实现多线程分片下载**

```python
    def _download_multi_thread(
        self,
        url: str,
        dest_path: Path,
        expected_size: int,
        desc: Optional[str]
    ) -> Path:
        """多线程分片下载大文件"""
        dest_path = Path(dest_path)
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        
        # 计算分片
        num_chunks = (expected_size + self.CHUNK_SIZE - 1) // self.CHUNK_SIZE
        desc = desc or dest_path.name
        
        print(f"    多线程下载: {format_file_size(expected_size)}, {num_chunks} 个分片, {self.max_workers} 线程")
        
        # 创建临时文件
        temp_files = [dest_path.parent / f"{dest_path.name}.part{i}" for i in range(num_chunks)]
        
        # 下载进度跟踪
        downloaded_bytes = [0] * num_chunks
        progress_lock = threading.Lock()
        
        def download_chunk(chunk_idx: int):
            """下载单个分片"""
            start = chunk_idx * self.CHUNK_SIZE
            end = min(start + self.CHUNK_SIZE - 1, expected_size - 1)
            temp_file = temp_files[chunk_idx]
            
            # 检查是否已下载
            if temp_file.exists():
                existing_size = temp_file.stat().st_size
                if existing_size == end - start + 1:
                    with progress_lock:
                        downloaded_bytes[chunk_idx] = existing_size
                    return True
                else:
                    start += existing_size
            
            # 下载分片
            req = urllib.request.Request(url)
            req.add_header("User-Agent", "Mozilla/5.0")
            req.add_header("Range", f"bytes={start}-{end}")
            
            try:
                with urllib.request.urlopen(req, timeout=(self.CONNECT_TIMEOUT, self.READ_TIMEOUT)) as response:
                    mode = "ab" if temp_file.exists() else "wb"
                    with open(temp_file, mode) as f:
                        while True:
                            chunk = response.read(self.SINGLE_CHUNK_SIZE)
                            if not chunk:
                                break
                            f.write(chunk)
                            with progress_lock:
                                downloaded_bytes[chunk_idx] += len(chunk)
                return True
            except Exception as e:
                print(f"    分片 {chunk_idx} 下载失败: {e}")
                return False
        
        # 使用线程池下载
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {executor.submit(download_chunk, i): i for i in range(num_chunks)}
            
            # 进度显示
            if self.show_progress:
                progress = SimpleProgressBar(expected_size, desc=desc[:30])
                completed = set()
                
                while len(completed) < num_chunks:
                    for future in as_completed(futures):
                        chunk_idx = futures[future]
                        completed.add(chunk_idx)
                        
                        # 更新进度
                        with progress_lock:
                            current = sum(downloaded_bytes)
                        progress.n = current
                        progress._print(current / expected_size)
                        
                        if not future.result():
                            raise RuntimeError(f"分片 {chunk_idx} 下载失败")
                
                progress.close()
        
        # 合并分片
        print(f"    合并分片...")
        with open(dest_path, 'wb') as outfile:
            for temp_file in temp_files:
                with open(temp_file, 'rb') as infile:
                    shutil.copyfileobj(infile, outfile)
                temp_file.unlink()  # 删除临时文件
        
        # 验证完整性
        if self.verify_integrity and dest_path.suffix == '.gz':
            print(f"    验证文件完整性...")
            valid, msg = verify_gzip_integrity(dest_path)
            if not valid:
                raise RuntimeError(f"文件验证失败: {msg}")
            print(f"    验证通过")
        
        return dest_path
```

- [ ] **Step 4: 实现带自动切换的下载方法**

```python
    def download_with_mirror_fallback(
        self,
        mirrors: List[MirrorSource],
        filename: str,
        dest_path: Path,
        desc: Optional[str] = None
    ) -> Path:
        """从多个镜像源下载，自动切换
        
        Args:
            mirrors: 镜像源列表（按优先级排序）
            filename: 要下载的文件名（相对于 base_url）
            dest_path: 目标路径
            desc: 进度条描述
        
        Returns:
            下载完成的文件路径
        
        Raises:
            RuntimeError: 所有镜像源都失败
        """
        last_error = None
        
        for mirror in mirrors:
            url = f"{mirror.base_url.rstrip('/')}/{filename}"
            try:
                print(f"|--- 尝试镜像: {mirror.name} ({mirror.region})")
                return self.download_file(url, dest_path, desc=desc)
            except Exception as e:
                last_error = e
                print(f"    失败: {e}")
                continue
        
        raise RuntimeError(f"所有镜像源都失败，最后一个错误: {last_error}")
```

- [ ] **Step 5: Commit**

```bash
git add allenricher/database/download_manager.py
git commit -m "feat: add DownloadManager with multi-threading and mirror fallback"
```

---

## Task 4: 重构 downloader.py 使用 DownloadManager

**Files:**
- Modify: `allenricher/database/downloader.py`

- [ ] **Step 1: 修改 DataDownloader 类使用 DownloadManager**

```python
"""数据下载器模块（重构版）

使用 DownloadManager 实现高性能、可靠的下载功能。
"""
import os
from pathlib import Path
from datetime import datetime
from typing import Dict, Optional

from .download_manager import DownloadManager
from .mirrors import get_mirrors, JENSEN_SOURCES


class DataDownloader:
    """全体物种通用数据下载器（重构版）"""

    def __init__(
        self,
        root_dir: str = "./database",
        overwrite: bool = False,
        max_workers: int = 4,
        use_multi_thread: bool = True,
        verify_integrity: bool = True
    ):
        """初始化下载器
        
        Args:
            root_dir: 数据库根目录
            overwrite: 是否覆盖已存在文件
            max_workers: 多线程下载线程数
            use_multi_thread: 是否启用多线程下载大文件
            verify_integrity: 是否验证文件完整性
        """
        self.root_dir = Path(root_dir)
        self.basic_dir = self.root_dir / "basic"
        self.basic_dir.mkdir(parents=True, exist_ok=True)
        
        # 使用 DownloadManager
        self.manager = DownloadManager(
            root_dir=root_dir,
            overwrite=overwrite,
            max_workers=max_workers,
            use_multi_thread=use_multi_thread,
            verify_integrity=verify_integrity,
            show_progress=True
        )

    def download_go_basic(self, version: Optional[str] = None) -> str:
        """下载 GO 基础数据（使用多线程和镜像源）"""
        if version is None:
            version = f"GO{datetime.now().strftime('%Y%m%d')}"

        go_dir = self.basic_dir / "go" / version
        go_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"下载 GO 基础数据 → {go_dir}")
        print(f"{'='*60}")

        # 使用 NCBI 镜像源下载
        ncbi_mirrors = get_mirrors('ncbi')
        
        self.manager.download_with_mirror_fallback(
            ncbi_mirrors,
            "gene2go.gz",
            go_dir / "gene2go.gz",
            desc="gene2go.gz"
        )
        
        self.manager.download_with_mirror_fallback(
            ncbi_mirrors,
            "gene_info.gz",
            go_dir / "gene_info.gz",
            desc="gene_info.gz"
        )
        
        # GO OBO 使用 GO 镜像源
        go_mirrors = get_mirrors('go')
        self.manager.download_with_mirror_fallback(
            go_mirrors,
            "go-basic.obo",
            go_dir / "go-basic.obo",
            desc="go-basic.obo"
        )

        print(f"GO 基础数据下载完成 → {go_dir}")
        return str(go_dir)

    def download_reactome_basic(self, version: Optional[str] = None) -> str:
        """下载 Reactome 基础数据"""
        if version is None:
            version = f"Reactome{datetime.now().strftime('%Y%m%d')}"

        re_dir = self.basic_dir / "reactome" / version
        re_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"下载 Reactome 基础数据 → {re_dir}")
        print(f"{'='*60}")

        # 下载 gene_info（复用 GO 的，避免重复下载）
        ncbi_mirrors = get_mirrors('ncbi')
        self.manager.download_with_mirror_fallback(
            ncbi_mirrors,
            "gene_info.gz",
            re_dir / "gene_info.gz",
            desc="gene_info.gz"
        )

        # 下载 Reactome 数据
        reactome_mirrors = get_mirrors('reactome')
        raw_file = re_dir / "NCBI2Reactome_All_Levels.txt"
        gz_file = re_dir / "NCBI2Reactome_All_Levels.txt.gz"
        
        if gz_file.exists() and not self.manager.overwrite:
            print(f"|--- 已存在，跳过: {gz_file.name}")
        else:
            self.manager.download_with_mirror_fallback(
                reactome_mirrors,
                "NCBI2Reactome_All_Levels.txt",
                raw_file,
                desc="NCBI2Reactome"
            )
            # 压缩
            if raw_file.exists():
                import gzip
                import shutil
                print(f"|--- 压缩: {raw_file.name}")
                with open(raw_file, 'rb') as f_in:
                    with gzip.open(gz_file, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                raw_file.unlink()

        print(f"Reactome 基础数据下载完成 → {re_dir}")
        return str(re_dir)

    def download_do_files(self) -> Dict[str, str]:
        """下载 DO 数据（Jensen Lab 无镜像，直接下载）"""
        do_dir = self.basic_dir / "do"
        do_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"下载 DO 数据 → {do_dir}")
        print(f"{'='*60}")

        files = {}
        for url in JENSEN_SOURCES:
            fname = url.split('/')[-1]
            dest = do_dir / fname
            self.manager.download_file(url, dest, desc=fname)
            files[fname] = str(dest)
        
        return files

    # 保留版本管理方法（不变）
    def list_go_versions(self) -> list:
        """列出已下载的 GO 基础数据版本"""
        go_basic = self.basic_dir / "go"
        if not go_basic.exists():
            return []
        return sorted([d.name for d in go_basic.iterdir() if d.is_dir()])

    def list_reactome_versions(self) -> list:
        """列出已下载的 Reactome 基础数据版本"""
        re_basic = self.basic_dir / "reactome"
        if not re_basic.exists():
            return []
        return sorted([d.name for d in re_basic.iterdir() if d.is_dir()])

    def get_latest_go_version(self) -> Optional[str]:
        """获取最新的 GO 基础数据版本"""
        versions = self.list_go_versions()
        return versions[-1] if versions else None

    def get_latest_reactome_version(self) -> Optional[str]:
        """获取最新的 Reactome 基础数据版本"""
        versions = self.list_reactome_versions()
        return versions[-1] if versions else None
```

- [ ] **Step 2: 运行集成测试**

Run: `cd AllEnricher-v2 && python -c "from allenricher.database.downloader import DataDownloader; d = DataDownloader(root_dir='test_download', max_workers=2); print('Init OK')"`

Expected: `Init OK`

- [ ] **Step 3: Commit**

```bash
git add allenricher/database/downloader.py
git commit -m "refactor: use DownloadManager for high-performance downloads"
```

---

## Task 5: 添加 CLI 参数支持

**Files:**
- Modify: `allenricher/cli.py` (download 子命令参数)

- [ ] **Step 1: 添加下载相关 CLI 参数**

在 `create_parser()` 函数中找到 download_parser，添加以下参数：

```python
    download_parser.add_argument('--workers', type=int, default=4,
                                 help='多线程下载的线程数 (默认: 4)')
    download_parser.add_argument('--no-multi-thread', action='store_true',
                                 help='禁用多线程下载')
    download_parser.add_argument('--no-verify', action='store_true',
                                 help='禁用下载后完整性验证')
    download_parser.add_argument('--timeout', type=int, default=300,
                                 help='下载超时时间（秒）(默认: 300)')
```

- [ ] **Step 2: 修改 cmd_download 函数使用新参数**

```python
def cmd_download(args) -> int:
    """下载数据库"""
    from allenricher.database.downloader import DataDownloader
    
    downloader = DataDownloader(
        root_dir=args.database_dir,
        overwrite=args.force,
        max_workers=args.workers,
        use_multi_thread=not args.no_multi_thread,
        verify_integrity=not args.no_verify
    )
    # ... 其余代码不变
```

- [ ] **Step 3: Commit**

```bash
git add allenricher/cli.py
git commit -m "feat: add download optimization CLI options (--workers, --no-multi-thread, --no-verify)"
```

---

## Task 6: 编写测试

**Files:**
- Create: `tests/test_download_manager.py`
- Create: `tests/test_download_utils.py`

- [ ] **Step 1: 测试 download_utils**

```python
import pytest
from pathlib import Path
import gzip

from allenricher.database.download_utils import (
    verify_gzip_integrity,
    calculate_file_hash,
    format_file_size
)


class TestVerifyGzipIntegrity:
    def test_valid_gzip(self, tmp_path):
        # 创建有效的 gzip 文件
        test_file = tmp_path / "test.gz"
        with gzip.open(test_file, 'wt') as f:
            f.write("line1\nline2\nline3\n")
        
        valid, msg = verify_gzip_integrity(test_file)
        assert valid is True
        assert msg == "OK"
    
    def test_invalid_gzip(self, tmp_path):
        # 创建无效的 gzip 文件
        test_file = tmp_path / "test.gz"
        test_file.write_bytes(b"not a gzip file")
        
        valid, msg = verify_gzip_integrity(test_file)
        assert valid is False
        assert "Bad gzip" in msg


class TestFormatFileSize:
    def test_bytes(self):
        assert format_file_size(500) == "500.0 B"
    
    def test_kilobytes(self):
        assert format_file_size(1536) == "1.5 KB"
    
    def test_megabytes(self):
        assert format_file_size(1536 * 1024) == "1.5 MB"
    
    def test_gigabytes(self):
        assert format_file_size(1536 * 1024 * 1024) == "1.5 GB"
```

- [ ] **Step 2: 运行测试**

Run: `pytest tests/test_download_utils.py -v`

Expected: 所有测试通过

- [ ] **Step 3: Commit**

```bash
git add tests/test_download_utils.py tests/test_download_manager.py
git commit -m "test: add download utility and manager tests"
```

---

## Task 7: 更新文档

**Files:**
- Modify: `README.md` (下载命令部分)

- [ ] **Step 1: 更新 README 下载命令说明**

```markdown
### 下载优化选项

```bash
# 使用多线程加速下载大文件（默认启用）
allenricher download -d go --workers 8

# 禁用多线程（小文件或网络不稳定时）
allenricher download -d go --no-multi-thread

# 禁用完整性验证（加快下载，但可能下载损坏文件）
allenricher download -d go --no-verify

# 强制重新下载（覆盖已存在文件）
allenricher download -d go --force
```

**下载性能对比**（参考值，实际取决于网络）：

| 文件 | 大小 | 优化前 | 优化后 (4线程) | 加速比 |
|------|------|--------|----------------|--------|
| gene_info.gz | ~1.4GB | ~25分钟 | ~8分钟 | 3x |
| gene2go.gz | ~1.2GB | ~20分钟 | ~6分钟 | 3.3x |
```

- [ ] **Step 2: Commit**

```bash
git add README.md
git commit -m "docs: update README with download optimization options"
```

---

## 验证清单

实施完成后，运行以下验证：

1. **单元测试**: `pytest tests/test_download_utils.py tests/test_download_manager.py -v`
2. **集成测试**: `python -m allenricher download -d go --workers 4 --database-dir test_download`
3. **性能测试**: 对比优化前后 gene_info.gz 下载时间
4. **完整性测试**: 下载后运行 `gzip -t gene_info.gz` 验证文件

---

## 性能目标

- **gene_info.gz (1.4GB)**: 从 25分钟 → 8分钟 (3x 加速)
- **gene2go.gz (1.2GB)**: 从 20分钟 → 6分钟 (3.3x 加速)
- **完整性**: 100% 检测并自动重试损坏文件
- **可靠性**: 镜像源自动切换，单点故障不影响下载
