"""Download the human hTFtarget transcription-factor target dataset."""

from pathlib import Path
from typing import Dict, List, Optional
import requests
import logging

logger = logging.getLogger(__name__)

# HTFtarget Download Link
HTFTARGET_DOWNLOAD_URL = (
    "https://guolab.wchscu.cn/static/hTFtarget/file_download/tf-target-infomation.txt"
)


class HTFtargetFetcher:
    """Retrieve the human hTFtarget regulatory interaction file."""

    REQUEST_TIMEOUT = 300  # 56MB file takes longer

    def __init__(self, basic_dir: str):
        """
        Args:
Basic_dir: Base Cache Directory
        """
        self.basic_dir = Path(basic_dir)
        self.basic_dir.mkdir(parents=True, exist_ok=True)

    def _get_cache_dir(self) -> Path:
        """Return the local cache directory for this data source."""
        cache_dir = self.basic_dir / "htftarget"
        cache_dir.mkdir(parents=True, exist_ok=True)
        return cache_dir

    @staticmethod
    def _is_valid_tabular_file(path: Path) -> bool:
        """Return whether a download is tabular data rather than an HTML error page."""
        if not path.is_file() or path.stat().st_size == 0:
            return False
        head = path.read_bytes()[:4096].lstrip().lower()
        return not head.startswith((b"<!doctype", b"<html")) and b"\t" in head

    def download(self, overwrite: bool = False) -> Path:
        """Download and validate the hTFtarget interaction file."""
        cache_dir = self._get_cache_dir()
        local_path = cache_dir / "tf-target-information.txt"

        if local_path.exists() and not overwrite and self._is_valid_tabular_file(local_path):
            logger.info(f"hTFtarget cached, skipping: {local_path}")
            return local_path
        if local_path.exists() and not overwrite:
            logger.warning("hTFtarget Cache is not a valid table and will be re-downloaded: %s", local_path)

        logger.info(f"Download hTFtarget: {HTFTARGET_DOWNLOAD_URL}")

        try:
            resp = requests.get(
                HTFTARGET_DOWNLOAD_URL,
                timeout=self.REQUEST_TIMEOUT,
                stream=True,
                headers={"User-Agent": "Mozilla/5.0"},
            )
            resp.raise_for_status()
        except requests.RequestException as e:
            raise RuntimeError(f"Could not close temporary folder: %s{e}") from e

        # Fluid Writing
        temporary_path = local_path.with_name(local_path.name + ".part")
        with open(temporary_path, "wb") as f:
            for chunk in resp.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)

        if not self._is_valid_tabular_file(temporary_path):
            temporary_path.unlink(missing_ok=True)
            raise RuntimeError("hTFtarget Downloading is not a valid tab data")
        temporary_path.replace(local_path)

        size_mb = local_path.stat().st_size / 1024 / 1024
        logger.info(f"hTFtarget saved: {local_path} ({size_mb: .1f} MB)")
        return local_path

    @staticmethod
    def get_info() -> Dict[str, str]:
        """Return source metadata suitable for provenance reporting."""
        return {
            "name": "hTFtarget",
            "version": "2020",
            "url": "https://guolab.wchscu.cn/hTFtarget/",
            "species": "Homo sapiens",
            "description": "Human TF-target based on ChIP-Seq",
        }

    @staticmethod
    def get_supported_species_records() -> List[tuple[int, str]]:
        """Return TaxID-keyed species coverage derived from downloaded source data."""
        return [(9606, "Homo sapiens")]
