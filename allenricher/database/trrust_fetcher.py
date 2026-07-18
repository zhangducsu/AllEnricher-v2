"""Download the human and mouse TRRUST v2 regulatory interaction datasets."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional

import requests

logger = logging.getLogger(__name__)

# *TRUST species code
TRRUST_SPECIES_MAP: Dict[str, str] = {
    "Homo sapiens": "hsa",
    "Mus musculus": "mmu",
}

TRRUST_SPECIES_TAXIDS: Dict[str, int] = {
    "Homo sapiens": 9606,
    "Mus musculus": 10090,
}

# Species Code * Download URL
TRRUST_DOWNLOAD_URLS: Dict[str, str] = {
    "hsa": "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv",
    "mmu": "https://www.grnpedia.org/trrust/data/trrust_rawdata.mouse.tsv",
}

# Request timeout (sec)
_TIMEOUT = 60

# User-Agent
_UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"


class TRRUSTFetcher:
    """Retrieve TRRUST v2 regulatory interactions for human or mouse."""

    def __init__(self, basic_dir: str):
        """
        Args:
basic_dir: Basic Data Directory (Cache Root Directory)
        """
        self._basic_dir = Path(basic_dir)

    def _get_cache_dir(self) -> Path:
        """Return the local cache directory for this data source."""
        cache_dir = self._basic_dir / "trrust" / "TRRUSTv2"
        cache_dir.mkdir(parents=True, exist_ok=True)
        return cache_dir

    def download_species(self, species: str, overwrite: bool = False) -> Path:
        """Download TRRUST data for one supported species."""
        species_code = self.get_species_code(species)
        if species_code is None:
            raise ValueError(
                f"TTRUST unsupported species: '{species}',"
                f"Supported species: {self.get_supported_species()}"
            )

        url = TRRUST_DOWNLOAD_URLS[species_code]
        cache_dir = self._get_cache_dir()
        local_file = cache_dir / f"trrust_rawdata.{species_code}.tsv"

        if local_file.exists() and not overwrite:
            logger.info("TRUST data cached and overloaded: %s", local_file)
            return local_file

        logger.info("Download TRRUST Data: %s -> %s", url, local_file)

        resp = requests.get(
            url,
            headers={"User-Agent": _UA},
            timeout=_TIMEOUT,
        )
        resp.raise_for_status()

        local_file.write_bytes(resp.content)
        logger.info(
            "TRRUST data download completed: %s (%d bytes)",
            local_file,
            len(resp.content),
        )

        return local_file

    def download_all(self, overwrite: bool = False) -> Dict[str, Path]:
        """Download TRRUST data for all supported species."""
        results: Dict[str, Path] = {}
        for species in self.get_supported_species():
            results[species] = self.download_species(species, overwrite=overwrite)
        return results

    @staticmethod
    def get_species_code(latin_name: str) -> Optional[str]:
        """Resolve the TRRUST species code from a scientific name."""
        return TRRUST_SPECIES_MAP.get(latin_name)

    @staticmethod
    def get_latin_name(species_code: str) -> Optional[str]:
        """Resolve a scientific name from the project species code."""
        code_to_name = {v: k for k, v in TRRUST_SPECIES_MAP.items()}
        return code_to_name.get(species_code)

    @staticmethod
    def get_supported_species() -> List[str]:
        """Return the species supported by this data source."""
        return list(TRRUST_SPECIES_MAP.keys())

    @staticmethod
    def get_supported_species_records() -> List[tuple[int, str]]:
        """Return TaxID-keyed species coverage derived from downloaded source data."""
        return [(taxid, name) for name, taxid in TRRUST_SPECIES_TAXIDS.items()]
