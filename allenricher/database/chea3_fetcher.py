"""Download ChEA3 libraries and access the optional online enrichment service."""

from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Any, Dict, List

import requests

logger = logging.getLogger(__name__)


class ChEA3Fetcher:
    """Retrieve ChEA3 GMT libraries and submit optional API requests."""

    # Name of species  Code map
    CHEA3_SPECIES_MAP: Dict[str, str] = {
        "Homo sapiens": "hsa",
    }

    # Name of 6 GMT libraries  URL mapping
    CHEA3_GMT_LIBS: Dict[str, str] = {
        "ENCODE": "https://maayanlab.cloud/chea3/assets/tflibs/ENCODE_ChIP-seq.gmt",
        "ReMap": "https://maayanlab.cloud/chea3/assets/tflibs/ReMap_ChIP-seq.gmt",
        "LiteratureChIP": "https://maayanlab.cloud/chea3/assets/tflibs/Literature_ChIP-seq.gmt",
        "GTExCoexpression": "https://maayanlab.cloud/chea3/assets/tflibs/GTEx_Coexpression.gmt",
        "ARCHS4Coexpression": "https://maayanlab.cloud/chea3/assets/tflibs/ARCHS4_Coexpression.gmt",
        "EnrichrQueries": "https://maayanlab.cloud/chea3/assets/tflibs/Enrichr_Queries.gmt",
    }

    # API Endpoint
    API_URL = "https://maayanlab.cloud/chea3/api/enrich/"

    # Request timeout (sec)
    TIMEOUT = 120

    # User-Agent
    UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"

    def __init__(self, basic_dir: str):
        """
        Args:
basic_dir: Database Basic Directory (Like ./database/basic)
        """
        self.basic_dir = Path(basic_dir)

    def _get_cache_dir(self) -> Path:
        """Return the local cache directory for this data source."""
        cache_dir = self.basic_dir / "chea3" / "ChEA3v2024"
        cache_dir.mkdir(parents=True, exist_ok=True)
        return cache_dir

    def download_gmt_library(self, lib_name: str, overwrite: bool = False) -> Path:
        """Download one named ChEA3 GMT library."""
        if lib_name not in self.CHEA3_GMT_LIBS:
            raise ValueError(
                f"Unknown library '{lib_name}', supported library: {list(self.CHEA3_GMT_LIBS.keys())}"
            )

        url = self.CHEA3_GMT_LIBS[lib_name]
        cache_dir = self._get_cache_dir()
        local_file = cache_dir / f"{lib_name}_tf.gmt"

        # Cache Check
        if local_file.exists() and not overwrite:
            logger.info("ChEA3 GMT library cached, skipping download: %s", local_file)
            return local_file

        logger.info("Download ChEA3 GMT Library: %s -> %s", url, local_file)

        resp = requests.get(
            url,
            headers={"User-Agent": self.UA},
            timeout=self.TIMEOUT,
        )
        resp.raise_for_status()

        content = resp.content
        first_line = next((line for line in content.splitlines() if line.strip()), b"")
        if first_line.lstrip().startswith(b"<") or len(first_line.split(b"\t")) < 3:
            raise ValueError(f"ChEA3 Return content is not valid GMT: {url}")
        local_file.write_bytes(content)
        logger.info(
            "ChEA3 GMT library download completed: %s (%d bytes)",
            local_file, len(content),
        )

        return local_file

    def download_all_gmt_libraries(
        self, overwrite: bool = False
    ) -> Dict[str, Path]:
        """Download all supported ChEA3 GMT libraries."""
        results: Dict[str, Path] = {}
        for lib_name in self.CHEA3_GMT_LIBS:
            results[lib_name] = self.download_gmt_library(lib_name, overwrite=overwrite)
        return results

    def enrich_api(
        self,
        gene_set: List[str],
        query_name: str | None = None,
        max_retries: int = 3,
        backoff_factor: float = 1.0,
    ) -> Dict[str, Any]:
        """Submit a gene list to the ChEA3 online enrichment endpoint."""
        if not gene_set:
            raise ValueError("gene_set cannot be empty")

        if query_name is None:
            query_name = ", ".join(gene_set[:5])

        payload = {
            "query_name": query_name,
            "gene_set": gene_set,
        }

        logger.info("Call ChEA3 API: %dGenome.", len(gene_set))

        # Index de-fielding the test mechanism
        last_exception = None
        for attempt in range(max_retries):
            try:
                resp = requests.post(
                    self.API_URL,
                    json=payload,
                    headers={"User-Agent": self.UA},
                    timeout=self.TIMEOUT,
                )
                resp.raise_for_status()

                result = resp.json()

                # Validate the response schema before accepting the library.
                self._validate_api_response(result)

                logger.info("ChEA3 API returns results, containing %d library", len(result))
                return result

            except (requests.RequestException, ValueError) as e:
                last_exception = e
                logger.warning(
                    "ChEA3 API request failed (trying)%d/%d): %s",
                    attempt + 1, max_retries, str(e)
                )
                if attempt < max_retries - 1:
                    sleep_time = backoff_factor * (2 ** attempt)
                    logger.info("Wait%.1fRetry in seconds...", sleep_time)
                    time.sleep(sleep_time)

        # All retry failed.
        logger.error("ChEA3 API request failed after the %d", max_retries)
        raise last_exception

    def _validate_api_response(self, result: Any) -> None:
        """Validate the structure of a ChEA3 API response."""
        # Verify whether or not to be a dictionary
        if not isinstance(result, dict):
            raise ValueError(
                f"API response format error: expected dic, actual{type(result).__name__}"
            )

        # Verify whether the expected library key is included
        expected_libs = set(self.CHEA3_GMT_LIBS.keys())
        actual_libs = set(result.keys())

        # Allows some of the libraries to return, but at least one library is needed
        if not actual_libs:
            raise ValueError("API response is empty, no library data returned")

        # Verify data formats for each library
        for lib_name, lib_data in result.items():
            if not isinstance(lib_data, list):
                raise ValueError(
                    f"\"Cook.\"{lib_name}'Data format error: expected list,"
                    f"Actual{type(lib_data).__name__}"
                )
            # Verify whether each entry in the list is a dictionary
            for i, entry in enumerate(lib_data):
                if not isinstance(entry, dict):
                    raise ValueError(
                        f"Error in formatting the {i} entry from `{lib_name} ':"
                        f"Expect, actual{type(entry).__name__}"
                    )

    @staticmethod
    def get_supported_species() -> List[str]:
        """Return the species supported by this data source."""
        return list(ChEA3Fetcher.CHEA3_SPECIES_MAP.keys())

    @staticmethod
    def get_supported_species_records() -> List[tuple[int, str]]:
        """Return TaxID-keyed species coverage derived from downloaded source data."""
        return [(9606, "Homo sapiens")]

    @staticmethod
    def get_library_names() -> List[str]:
        """Return the available ChEA3 library names."""
        return list(ChEA3Fetcher.CHEA3_GMT_LIBS.keys())
