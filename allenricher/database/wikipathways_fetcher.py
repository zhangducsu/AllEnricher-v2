"""Download versioned, species-specific WikiPathways GMT files."""

import re
from pathlib import Path
from typing import Dict, List, Optional, Set
import requests


# The map of allEnricher species code for the Latin name of WikiPathways
# Comprising 38 species: 18 with GMT files + 20 only GPML
SPECIES_NAME_MAP: Dict[str, str] = {
    # = GMT species (18) = =
    "Anopheles gambiae": "aga",
    "Arabidopsis thaliana": "ath",
    "Bos taurus": "bta",
    "Caenorhabditis elegans": "cel",
    "Canis familiaris": "cfa",
    "Danio rerio": "dre",
    "Drosophila melanogaster": "dme",
    "Equus caballus": "eca",
    "Gallus gallus": "gga",
    "Homo sapiens": "hsa",
    "Mus musculus": "mmu",
    "Pan troglodytes": "ptr",
    "Populus trichocarpa": "ptc",
    "Rattus norvegicus": "rno",
    "Saccharomyces cerevisiae": "sce",
    "Solanum lycopersicum": "sly",
    "Sus scrofa": "ssc",
    "Zea mays": "zma",
    # == @elder_man
    "Acetobacterium woodii": "awo",
    "Bacillus subtilis": "bsu",
    "Beta vulgaris": "bvu",
    "Brassica napus": "bna",
    "Caulobacter vibrioides": "cvi",
    "Citrus sinensis": "csi",
    "Coffea arabica": "car",
    "Daphnia magna": "dma",
    "Escherichia coli": "eco",
    "Gibberella zeae": "gze",
    "Hordeum vulgare": "hvu",
    "Ilex paraguariensis": "ipa",
    "Mycobacterium tuberculosis": "mtu",
    "Oryza sativa": "osa",
    "Paullinia cupana": "pcu",
    "Perilla frutescens": "pfr",
    "Plasmodium falciparum": "pfa",
    "Theobroma cacao": "tcc",
    "Triticum aestivum": "tae",
    "Vitis vinifera": "vvi",
}

# 18 species with GMT files (Latin names in aggregate)
SPECIES_WITH_GMT: Set[str] = {
    "Anopheles gambiae",
    "Arabidopsis thaliana",
    "Bos taurus",
    "Caenorhabditis elegans",
    "Canis familiaris",
    "Danio rerio",
    "Drosophila melanogaster",
    "Equus caballus",
    "Gallus gallus",
    "Homo sapiens",
    "Mus musculus",
    "Pan troglodytes",
    "Populus trichocarpa",
    "Rattus norvegicus",
    "Saccharomyces cerevisiae",
    "Solanum lycopersicum",
    "Sus scrofa",
    "Zea mays",
}

# Only 20 species with GPML files (Latin names combined)
SPECIES_GPML_ONLY: Set[str] = {
    "Acetobacterium woodii",
    "Bacillus subtilis",
    "Beta vulgaris",
    "Brassica napus",
    "Caulobacter vibrioides",
    "Citrus sinensis",
    "Coffea arabica",
    "Daphnia magna",
    "Escherichia coli",
    "Gibberella zeae",
    "Hordeum vulgare",
    "Ilex paraguariensis",
    "Mycobacterium tuberculosis",
    "Oryza sativa",
    "Paullinia cupana",
    "Perilla frutescens",
    "Plasmodium falciparum",
    "Theobroma cacao",
    "Triticum aestivum",
    "Vitis vinifera",
}


class WikiPathwaysFetcher:
    """Retrieve versioned WikiPathways GMT data by species."""

    BASE_URL = "https://data.wikipathways.org"
    REQUEST_TIMEOUT = 60
    UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"

    def __init__(self, basic_dir: str):
        """
        Args:
basic_dir: Basic Cache Directory (To be created under this directory wikipathways Subdirectories)
        """
        self.basic_dir = Path(basic_dir)
        self.basic_dir.mkdir(parents=True, exist_ok=True)

    def _get_cache_dir(self, version: str) -> Path:
        """Return the local cache directory for this data source."""
        cache_dir = self.basic_dir / "wikipathways" / f"WP{version}"
        cache_dir.mkdir(parents=True, exist_ok=True)
        return cache_dir

    def _detect_latest_version(self) -> str:
        """Resolve the newest available WikiPathways release."""
        url = self.BASE_URL
        headers = {"User-Agent": self.UA}

        try:
            resp = requests.get(url, headers=headers, timeout=self.REQUEST_TIMEOUT)
            resp.raise_for_status()
        except requests.RequestException as e:
            raise RuntimeError(f"Could not close temporary folder: %s{e}") from e

        # Parsing the date folder in HTML (format: YYYYMMMDD or YYYYMMMDD/)
        # Matches similar: <a href="20240510">20240510</a>or<a href="20240510/">20240510/</a>
        pattern = r'href=["\'](\d{8})/?["\']'
        matches = re.findall(pattern, resp.text)

        if not matches:
            raise RuntimeError(f"Can't get from{url}Parsing version number")

        # Returns the latest version number (sortling by string is enough, as it is in YYYYMMDD format)
        latest = sorted(matches)[-1]
        return latest

    def get_available_species(self, version: Optional[str] = None) -> List[str]:
        """Return species available in one WikiPathways release."""
        if version is None:
            version = self._detect_latest_version()

        url = f"{self.BASE_URL}/{version}/gmt/"
        headers = {"User-Agent": self.UA}

        try:
            resp = requests.get(url, headers=headers, timeout=self.REQUEST_TIMEOUT)
            resp.raise_for_status()
        except requests.RequestException as e:
            raise RuntimeError(f"Can not get folder: %s: %s{e}") from e

        # Parsing GMT file links
        # Filename Format: wikipathways-{version}-gmt-{Species_Latin_Name}.gmt
        pattern = rf'wikipathways-{version}-gmt-([\w\s]+)\.gmt'
        matches = re.findall(pattern, resp.text)

        # Clear species names (remove URL code, etc.)
        species_list: List[str] = []
        for match in matches:
            # Replace underlined with spaces
            species_name = match.replace("_", " ")
            species_list.append(species_name)

        return sorted(species_list)

    def download_gmt(
        self,
        species_latin_name: str,
        version: Optional[str] = None,
        overwrite: bool = False
    ) -> Path:
        """Download the WikiPathways GMT file for one species."""
        if version is None:
            version = self._detect_latest_version()

        # Build filename (space replaced by underlined)
        species_filename = species_latin_name.replace(" ", "_")
        filename = f"wikipathways-{version}-gmt-{species_filename}.gmt"

        url = f"{self.BASE_URL}/{version}/gmt/{filename}"
        cache_dir = self._get_cache_dir(version)
        local_path = cache_dir / filename

        # Check if there is a presence
        if local_path.exists() and not overwrite:
            print(f"|---Cached, Skipped: {species_latin_name}")
            return local_path

        # Download the species-specific GMT archive.
        print(f"|--- Downloading WikiPathways for {species_latin_name} ({filename})")
        headers = {"User-Agent": self.UA}

        try:
            resp = requests.get(url, headers=headers, timeout=self.REQUEST_TIMEOUT)
            resp.raise_for_status()
        except requests.RequestException as e:
            raise RuntimeError(f"Download Failed{url}: {e}") from e

        # Save File
        with open(local_path, "wb") as f:
            f.write(resp.content)

        print(f"|--- Saved WikiPathways archive: {local_path}")
        return local_path

    def download_all_gmt(
        self,
        version: Optional[str] = None,
        overwrite: bool = False
    ) -> Dict[str, Path]:
        """Download WikiPathways GMT files for every available species."""
        if version is None:
            version = self._detect_latest_version()

        print(f"\n{'='*60}")
        print(f"WikiPathways GMT Download (version={version})")
        print(f"{'='*60}")

        # Get list of available species
        available_species = self.get_available_species(version)
        print(f"|---Available species: {len(available_species)}One.")

        # Download GMT files for each species
        results: Dict[str, Path] = {}
        for species in available_species:
            try:
                path = self.download_gmt(species, version, overwrite)
                results[species] = path
            except RuntimeError as e:
                print(f"--- Error: {e}")

        print(f"|--- Download complete: {len(results)}/{len(available_species)} species")
        return results

    @staticmethod
    def get_species_code(latin_name: str) -> Optional[str]:
        """Resolve the TRRUST species code from a scientific name."""
        return SPECIES_NAME_MAP.get(latin_name)

    @staticmethod
    def get_latin_name(species_code: str) -> Optional[str]:
        """Resolve a scientific name from the project species code."""
        for latin, code in SPECIES_NAME_MAP.items():
            if code == species_code:
                return latin
        return None

    @staticmethod
    def get_species_data_type(latin_name: str) -> str:
        """Return the WikiPathways data format available for a species."""
        if latin_name in SPECIES_WITH_GMT:
            return "gmt"
        elif latin_name in SPECIES_GPML_ONLY:
            return "gpml"
        else:
            return "none"

    @staticmethod
    def get_supported_species() -> Dict[str, List[str]]:
        """Return the species supported by this data source."""
        return {
            "gmt": sorted(list(SPECIES_WITH_GMT)),
            "gpml": sorted(list(SPECIES_GPML_ONLY)),
            "all": sorted(list(SPECIES_NAME_MAP.keys())),
        }

    def get_local_gmt_path(self, species_latin_name: str, version: str) -> Optional[Path]:
        """Return the cached GMT path without downloading it."""
        species_filename = species_latin_name.replace(" ", "_")
        filename = f"wikipathways-{version}-gmt-{species_filename}.gmt"
        cache_dir = self._get_cache_dir(version)
        local_path = cache_dir / filename

        if local_path.exists():
            return local_path
        return None

    def list_cached_versions(self) -> List[str]:
        """Return locally cached WikiPathways releases."""
        wp_dir = self.basic_dir / "wikipathways"
        if not wp_dir.exists():
            return []

        versions: List[str] = []
        for item in wp_dir.iterdir():
            if item.is_dir() and item.name.startswith("WP"):
                version = item.name[2:]  # Remove the prefix "WP".
                if re.match(r"^\d{8}$", version):
                    versions.append(version)

        return sorted(versions)
