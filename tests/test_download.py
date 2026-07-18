"""Download Tool Function and Download Manager Module Test"""
import gzip
import pytest
import urllib.error
from pathlib import Path

from allenricher.database.download_utils import (
    verify_gzip_integrity,
    calculate_file_hash,
    format_size,
    format_speed,
    format_duration,
    SimpleProgressBar,
)
from allenricher.database.mirrors import (
    MirrorSource,
    get_mirrors,
    NCBI_MIRRORS,
    GO_MIRRORS,
    REACTOME_MIRRORS,
    JENSEN_SOURCES,
)
from allenricher.database.download_manager import DownloadManager


# ================================================================
# Download_utils test
# ================================================================

class TestVerifyGzipIntegrity:
    def test_valid_gzip(self, tmp_path):
        """Normal gzip file should be authenticated"""
        f = tmp_path / "valid.gz"
        with gzip.open(f, 'wt') as gf:
            gf.write("col1\tcol2\tcol3\n")
            gf.write("a\tb\tc\n")
            gf.write("d\te\tf\n")
        valid, msg = verify_gzip_integrity(f)
        assert valid is True
        assert msg == "OK"

    def test_invalid_gzip_not_gz(self, tmp_path):
        """Non gzip file should fail"""
        f = tmp_path / "fake.gz"
        f.write_bytes(b"this is not a gzip file at all")
        valid, msg = verify_gzip_integrity(f)
        assert valid is False
        assert "Bad gzip" in msg

    def test_truncated_gzip(self, tmp_path):
        """Disclosed gzip file should fail (Verification with full amount)"""
        f = tmp_path / "trunc.gz"
        with gzip.open(f, 'wt') as gf:
            for i in range(10000):
                gf.write(f"line{i}\tcol2\tcol3\n")
        # Cut File
        size = f.stat().st_size
        f.write_bytes(f.read_bytes()[:size // 2])
        # Sample verification may miss detection, full-scale verification
        valid, msg = verify_gzip_integrity(f, sample_lines=0)
        assert valid is False

    def test_gzip_with_header_only(self, tmp_path):
        """Gzip files only for comment headers should be adopted"""
        f = tmp_path / "header.gz"
        with gzip.open(f, 'wt') as gf:
            gf.write("#tax_id\tGeneID\tSymbol\n")
            gf.write("9606\t1\tA1BG\n")
        valid, msg = verify_gzip_integrity(f)
        assert valid is True


class TestFormatSize:
    def test_bytes(self):
        assert format_size(500) == "500.0 B"

    def test_kilobytes(self):
        assert format_size(1536) == "1.5 KB"

    def test_megabytes(self):
        assert format_size(1536 * 1024) == "1.5 MB"

    def test_gigabytes(self):
        assert format_size(1.4 * 1024 ** 3) == "1.4 GB"

    def test_zero(self):
        assert format_size(0) == "0.0 B"


class TestFormatSpeed:
    def test_normal(self):
        assert "MB/s" in format_speed(5 * 1024 * 1024)

    def test_zero(self):
        assert format_speed(0) == "---"


class TestFormatDuration:
    def test_seconds(self):
        assert format_duration(30) == "30s"

    def test_minutes(self):
        assert format_duration(150) == "2m 30s"

    def test_hours(self):
        assert format_duration(3661) == "1h 1m"


class TestCalculateFileHash:
    def test_md5(self, tmp_path):
        f = tmp_path / "test.txt"
        f.write_text("hello world")
        h = calculate_file_hash(f, "md5")
        assert len(h) == 32  # MD5 hex = 32 chars

    def test_sha256(self, tmp_path):
        f = tmp_path / "test.txt"
        f.write_text("hello world")
        h = calculate_file_hash(f, "sha256")
        assert len(h) == 64  # SHA256 hex = 64 chars

    def test_deterministic(self, tmp_path):
        f = tmp_path / "test.txt"
        f.write_text("same content")
        h1 = calculate_file_hash(f)
        h2 = calculate_file_hash(f)
        assert h1 == h2


class TestSimpleProgressBar:
    def test_init(self):
        p = SimpleProgressBar(total=1000, desc="test")
        assert p.total == 1000
        assert p.n == 0

    def test_update(self):
        p = SimpleProgressBar(total=1000, desc="test")
        p.update(500)
        assert p.n == 500

    def test_close(self, capsys):
        p = SimpleProgressBar(total=100, desc="test")
        p.close()
        captured = capsys.readouterr()
        assert '\n' in captured.out


class TestMultiThreadDownload:
    def test_chunk_retry_completes_without_deleting_other_parts(self, tmp_path, monkeypatch):
        source = b"abcdefghijkl"
        attempts = {}

        class Response:
            status = 206

            def __init__(self, payload):
                self.payload = payload

            def __enter__(self):
                return self

            def __exit__(self, *_args):
                return False

            def read(self, _size):
                payload, self.payload = self.payload, b""
                return payload

        def fake_urlopen(request, timeout):
            byte_range = request.get_header("Range")
            start, end = map(int, byte_range.removeprefix("bytes=").split("-"))
            attempts[start] = attempts.get(start, 0) + 1
            if start == 0 and attempts[start] == 1:
                raise urllib.error.URLError("temporary")
            return Response(source[start:end + 1])

        monkeypatch.setattr("urllib.request.urlopen", fake_urlopen)
        manager = DownloadManager(
            overwrite=True,
            max_workers=2,
            verify_integrity=False,
            show_progress=False,
        )
        manager.CHUNK_SIZE = 4
        output = tmp_path / "data.bin"

        manager._download_multi_thread("https://example.test/data", output, len(source), "data")

        assert output.read_bytes() == source
        assert attempts[0] == 2


# ================================================================
# mmirors test
# ================================================================

class TestMirrors:
    def test_ncbi_mirrors_exist(self):
        assert len(NCBI_MIRRORS) >= 1

    def test_ncbi_mirrors_priority(self):
        priorities = [m.priority for m in NCBI_MIRRORS]
        assert priorities == sorted(priorities)

    def test_get_mirrors_ncbi(self):
        mirrors = get_mirrors('ncbi')
        assert mirrors[0].name == "ncbi-official"

    def test_get_mirrors_go(self):
        mirrors = get_mirrors('go')
        assert len(mirrors) >= 1
        assert all('obo' in m.base_url.lower() or 'geneontology' in m.base_url.lower() for m in mirrors)

    def test_get_mirrors_reactome(self):
        mirrors = get_mirrors('reactome')
        assert len(mirrors) >= 1

    def test_get_mirrors_unknown(self):
        mirrors = get_mirrors('nonexistent')
        assert mirrors == []

    def test_jensen_sources(self):
        assert len(JENSEN_SOURCES) == 3
        assert all('jensenlab.org' in url for url in JENSEN_SOURCES)

    def test_mirror_source_dataclass(self):
        m = MirrorSource("test", "https://example.com/", 1, "US")
        assert m.name == "test"
        assert m.enabled is True
