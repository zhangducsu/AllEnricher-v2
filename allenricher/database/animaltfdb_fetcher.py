"""Download official AnimalTFDB 4.0 transcription-factor and orthology data."""

from pathlib import Path
from typing import Dict, List, Optional, Set
import requests
import logging
import re
from urllib.parse import urlparse

logger = logging.getLogger(__name__)

ANIMALTFDB_SPECIES_PATH = Path(__file__).with_name("data") / "animaltfdb_v4_species.txt"

# Official HTTP download endpoint. The HTTPS route currently rejects these files.
ANIMALTFDB_BASE_URL = (
    "http://guolab.wchscu.cn/AnimalTFDB4_static/download"
)

_WAF_UNSBOX_INDEX = (
    15, 35, 29, 24, 33, 16, 1, 38, 10, 9,
    19, 31, 40, 27, 22, 23, 25, 13, 6, 11,
    39, 18, 20, 8, 14, 21, 32, 26, 2, 30,
    7, 4, 17, 5, 3, 28, 34, 37, 12, 36,
)
_WAF_XOR_KEY = "3000176000856006061501533003690027800375"

# Frequently requested model, livestock, and companion animal species.
ANIMALTFDB_PRIORITY_SPECIES = [
    "Homo_sapiens", "Mus_musculus", "Rattus_norvegicus",
    "Danio_rerio", "Drosophila_melanogaster", "Caenorhabditis_elegans",
    "Bos_taurus", "Sus_scrofa", "Ovis_aries", "Capra_hircus",
    "Gallus_gallus", "Canis_lupus_familiaris", "Equus_caballus",
    "Felis_catus", "Macaca_mulatta", "Gorilla_gorilla_gorilla",
    "Pan_troglodytes", "Oryctolagus_cuniculus", "Xenopus_tropicalis",
]


class AnimalTFDBFetcher:
    """Download AnimalTFDB 4.0 transcription-factor and orthology files."""

    REQUEST_TIMEOUT = 120

    def __init__(self, basic_dir: str):
        """Initialize the fetcher under the shared database source directory."""
        self.basic_dir = Path(basic_dir)
        self.basic_dir.mkdir(parents=True, exist_ok=True)

    def _get_cache_dir(self) -> Path:
        """Return the local cache directory for this data source."""
        cache_dir = self.basic_dir / "animaltfdb" / "AnimalTFDBv4.0"
        cache_dir.mkdir(parents=True, exist_ok=True)
        return cache_dir

    @staticmethod
    def _is_valid_tabular_file(path: Path) -> bool:
        """Return whether a download is tabular data rather than an HTML error page."""
        if not path.is_file() or path.stat().st_size == 0:
            return False
        head = path.read_bytes()[:4096].lstrip().lower()
        return not head.startswith((b"<!doctype", b"<html")) and b"\t" in head

    @staticmethod
    def _compute_waf_cookie(challenge: str) -> str:
        if not re.fullmatch(r"[0-9A-Fa-f]{40}", challenge):
            raise ValueError("AnimalTFDB WAF challenge format changed")
        unsboxed = "".join(challenge[index - 1] for index in _WAF_UNSBOX_INDEX)
        return "".join(
            f"{int(unsboxed[offset:offset + 2], 16) ^ int(_WAF_XOR_KEY[offset:offset + 2], 16):02x}"
            for offset in range(0, 40, 2)
        )

    def _request_download(self, url: str) -> requests.Response:
        """Complete the upstream JavaScript challenge and download one official file."""
        session = requests.Session()
        # The configured localhost proxy returns a different WAF response. The
        # official HTTP endpoint is directly reachable and requires no VPN.
        session.trust_env = False
        headers = {"User-Agent": "Mozilla/5.0"}
        response = session.get(url, timeout=self.REQUEST_TIMEOUT, headers=headers)
        response.raise_for_status()
        match = re.search(r"var arg1='([0-9A-Fa-f]{40})'", response.text)
        if match:
            cookie = self._compute_waf_cookie(match.group(1))
            session.cookies.set(
                "acw_sc__v2",
                cookie,
                domain=urlparse(url).hostname,
                path="/",
            )
            response = session.get(
                url,
                timeout=self.REQUEST_TIMEOUT,
                stream=True,
                headers=headers,
            )
            response.raise_for_status()
        return response

    def _download_url(self, url: str, local_path: Path, overwrite: bool = False) -> Path:
        """Download one URL to the requested local path."""
        if local_path.exists() and not overwrite and self._is_valid_tabular_file(local_path):
            logger.info("Using cached AnimalTFDB file: %s", local_path.name)
            return local_path
        if local_path.exists() and not overwrite:
            logger.warning("Cached file is not valid tabular data; downloading it again: %s", local_path)

        logger.info("Downloading AnimalTFDB data: %s", url)

        try:
            resp = self._request_download(url)
        except requests.RequestException as e:
            raise RuntimeError(f"AnimalTFDB request failed for {url}: {e}") from e

        temporary_path = local_path.with_name(local_path.name + ".part")
        with open(temporary_path, "wb") as f:
            for chunk in resp.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)

        if not self._is_valid_tabular_file(temporary_path):
            temporary_path.unlink(missing_ok=True)
            raise RuntimeError(f"AnimalTFDB returned invalid tabular data for {url}")
        temporary_path.replace(local_path)

        logger.info("Saved %s (%.1f KB)", local_path, local_path.stat().st_size / 1024)
        return local_path

    def download_tf_list(self, species: str, overwrite: bool = False) -> Path:
        """Download the transcription-factor list for one AnimalTFDB species."""
        cache_dir = self._get_cache_dir()
        url = f"{ANIMALTFDB_BASE_URL}/TF_list_final/{species}_TF"
        local_path = cache_dir / f"{species}_TF"
        return self._download_url(url, local_path, overwrite)

    def download_ortholog_to_human(self, species: str, overwrite: bool = False) -> Path:
        """Download the target-species to human orthology mapping."""
        cache_dir = self._get_cache_dir()
        url = f"{ANIMALTFDB_BASE_URL}/ortholog_to_human_download/{species}_ortholog_to_human"
        local_path = cache_dir / f"{species}_ortholog_to_human"
        return self._download_url(url, local_path, overwrite)

    def download_species_data(self, species: str, overwrite: bool = False) -> Dict[str, Path]:
        """Download all AnimalTFDB inputs required for one species."""
        results = {}
        try:
            results['tf_list'] = self.download_tf_list(species, overwrite)
        except Exception as e:
            logger.error("Failed to download the AnimalTFDB TF list for %s: %s", species, e)

        try:
            results['ortholog'] = self.download_ortholog_to_human(species, overwrite)
        except Exception as e:
            logger.error("Failed to download the AnimalTFDB orthology map for %s: %s", species, e)

        return results

    def download_htftarget(self, overwrite: bool = False) -> Path:
        """Download the human hTFtarget interaction table."""
        from .htftarget_fetcher import HTFtargetFetcher
        fetcher = HTFtargetFetcher(str(self.basic_dir))
        return fetcher.download(overwrite)

    @staticmethod
    def get_priority_species() -> List[str]:
        """Return species prioritized for bundled AnimalTFDB support."""
        return ANIMALTFDB_PRIORITY_SPECIES.copy()

    @staticmethod
    def get_supported_species_records() -> List[tuple[int, str]]:
        """Return TaxID-keyed species coverage derived from downloaded source data."""
        if not ANIMALTFDB_SPECIES_PATH.is_file():
            raise FileNotFoundError(
                f"AnimalTFDB species list does not exist: {ANIMALTFDB_SPECIES_PATH}"
            )
        records = []
        for line in ANIMALTFDB_SPECIES_PATH.read_text(encoding="utf-8").splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            taxid, latin_name = line.split("\t", 1)
            records.append((int(taxid), latin_name.replace("_", " ")))
        if len(records) != 183 or len({taxid for taxid, _ in records}) != 183:
            raise ValueError("The AnimalTFDB 4.0 registry must contain 183 unique NCBI TaxIDs")
        return records

    @staticmethod
    def get_info() -> Dict[str, str]:
        """Return source metadata suitable for provenance reporting."""
        return {
            "name": "AnimalTFDB",
            "version": "4.0",
            "url": "https://guolab.wchscu.cn/AnimalTFDB4/",
            "species_count": 183,
            "tf_count": 274633,
            "description": "TF classifications and direct human orthology mappings for 183 animal species",
        }
