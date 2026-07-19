"""Retrieve species-specific pathway annotations through the KEGG REST API."""

import gzip
import re
import time
import urllib.request
import urllib.error
from pathlib import Path
from typing import Dict, List, Optional, Tuple


class KEGGFetcher:
    """Fetch pathway membership, names, and categories for one KEGG organism."""

    BASE_URL = "https://rest.kegg.jp"
    # Kegg requires no more than 10 requests per second
    REQUEST_INTERVAL = 0.15  # sec
    MAX_RETRIES = 3
    UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"

    def __init__(self, cache_dir: str, overwrite: bool = False):
        """
        Args:
cache_dir: Cache Directory (Store downloaded raw data and generated files)
overwrite: Whether to overwrite the cached data
        """
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.overwrite = overwrite

    # ================================================================
    # Public interface
    # ================================================================

    def fetch_species_data(
        self,
        species: str,
        gene_info_path: str,
        taxid: Optional[int] = None,
    ) -> Tuple[str, str]:
        """Retrieve and convert source annotations for one species."""
        print(f"\n{'='*60}")
        print(f"KEGREST API Data Acquisition (SPECies={species})")
        print(f"{'='*60}")

        # Step 1: Retrieve the pathway list.
        pathways = self._list_pathways(species)
        print(f"|--- Retrieved {len(pathways)} KEGG pathways.")

        # Step 2: Retrieve pathway-to-gene associations.
        gene_pathway_links = self._get_gene_pathway_links(species)
        print(f"|---Gene-pathway associations: {len(gene_pathway_links)}")

        # Step 3: Get KEGID *NCBI Gene ID map
        kegg_to_ncbi = self._get_kegg_ncbi_mapping(species)
        print(f"|---ID map: {len(kegg_to_ncbi)}KEGG Gene")

        # Step 4: Build NCBI Gene ID & Symbol Map
        ncbi_to_symbol = self._ncbi_id_to_symbol(gene_info_path, taxid)
        print(f"|---Symbol: {len(ncbi_to_symbol)}One.")

        # Step 5: Generate Gene2pathway.txt
        gene2pathway_path = self._build_gene2pathway(
            species, gene_pathway_links, kegg_to_ncbi, ncbi_to_symbol, pathways
        )

        # Step 6: Generate path_summary.txt
        pathway_summary_path = self._build_pathway_summary(species, pathways)

        print(f"|---Kegg Data Acquisition Completed")
        return str(gene2pathway_path), str(pathway_summary_path)

    # ================================================================
    # API Call
    # ================================================================

    def _api_get(self, endpoint: str) -> str:
        """Call a KEGG REST endpoint with retry handling."""
        url = f"{self.BASE_URL}/{endpoint}"

        # Use priority (if any); the big response will occasionally be cut in advance by KEGG.
        requests_error = None
        try:
            import requests
        except ImportError:
            requests = None
        if requests is not None:
            for attempt in range(1, self.MAX_RETRIES + 1):
                try:
                    resp = requests.get(url, headers={"User-Agent": self.UA}, timeout=60)
                    resp.raise_for_status()
                    return resp.text
                except requests.RequestException as exc:
                    requests_error = exc
                    if attempt < self.MAX_RETRIES:
                        time.sleep(0.5 * (2 ** (attempt - 1)))

        # Back to urllib
        req = urllib.request.Request(url)
        req.add_header("User-Agent", self.UA)

        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                return resp.read().decode("utf-8")
        except urllib.error.HTTPError as e:
            raise RuntimeError(f"KEGG API HTTP {e.code}: {url}") from e
        except urllib.error.URLError as e:
            detail = f"error: {requests_error}" if requests_error else ""
            raise RuntimeError(f"KEGAPI error: {e.reason}: {url}{detail}") from e

    def _list_pathways(self, species: str) -> List[Tuple[str, str]]:
        """Return all pathways assigned to the selected KEGG organism."""
        cache_file = self.cache_dir / f"{species}_pathways.txt"

        if cache_file.exists() and not self.overwrite:
            print(f"|---Cached, Skipping: list/pathway/{species}")
            pathways = []
            with open(cache_file, "r") as f:
                for line in f:
                    parts = line.strip().split("\t", 1)
                    if len(parts) == 2:
                        pathway_id = parts[0]
                        if pathway_id.startswith(species):
                            pathway_id = pathway_id[len(species):]
                        pathways.append((pathway_id, parts[1]))
            return pathways

        print(f"|--- API: list/pathway/{species}")
        data = self._api_get(f"list/pathway/{species}")
        time.sleep(self.REQUEST_INTERVAL)

        # Solve path lists, extract classification levels
        pathways = []
        pathway_names: Dict[str, str] = {}  # pathway_id → full name

        for line in data.strip().split("\n"):
            parts = line.split("\t", 1)
            if len(parts) != 2:
                continue
            full_id = parts[0].strip()
            full_name = parts[1].strip()

            # Extracting a pure circuit number (specify prefixes)
            pathway_id = full_id.replace(species, "", 1)
            pathway_names[pathway_id] = full_name

        # Save Cache
        with open(cache_file, "w") as f:
            for pid, pname in pathway_names.items():
                f.write(f"{pid}\t{pname}\n")

        return [(pid, pname) for pid, pname in pathway_names.items()]

    def _get_gene_pathway_links(self, species: str) -> Dict[str, List[str]]:
        """Retrieve gene-to-pathway memberships from KEGG."""
        cache_file = self.cache_dir / f"{species}_gene_pathway_links.txt"

        if cache_file.exists() and not self.overwrite:
            print(f"|---Cached, Skipping: link/{species}/pathway")
            links = {}
            with open(cache_file, "r") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) == 2:
                        gene, pw = parts
                        links.setdefault(gene, []).append(pw)
            return links

        print(f"|--- API: link/{species}/pathway")
        data = self._api_get(f"link/{species}/pathway")
        time.sleep(self.REQUEST_INTERVAL)

        links: Dict[str, List[str]] = {}
        for line in data.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) != 2:
                continue
            # path:hsa00010 → 00010
            pw = parts[0].replace(f"path:{species}", "")
            # hsa:10327 → 10327
            gene = parts[1].replace(f"{species}:", "")
            links.setdefault(gene, []).append(pw)

        # Save Cache
        with open(cache_file, "w") as f:
            for gene, pws in links.items():
                for pw in pws:
                    f.write(f"{gene}\t{pw}\n")

        return links

    def _get_kegg_ncbi_mapping(self, species: str) -> Dict[str, str]:
        """Retrieve KEGG-to-NCBI gene identifier mappings."""
        cache_file = self.cache_dir / f"{species}_kegg_ncbi_map.txt"

        if cache_file.exists() and not self.overwrite:
            print(f"|---Cached, Skipped: conv/{species}/ncbi-geneid")
            mapping = {}
            with open(cache_file, "r") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) == 2:
                        ncbi, kegg = parts
                        mapping[kegg] = ncbi
            return mapping

        print(f"|--- API: conv/{species}/ncbi-geneid")
        data = self._api_get(f"conv/{species}/ncbi-geneid")
        time.sleep(self.REQUEST_INTERVAL)

        mapping: Dict[str, str] = {}
        for line in data.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) != 2:
                continue
            ncbi = parts[0].replace("ncbi-geneid:", "")
            kegg = parts[1].replace(f"{species}:", "")
            mapping[kegg] = ncbi

        # Save Cache
        with open(cache_file, "w") as f:
            for kegg, ncbi in mapping.items():
                f.write(f"{ncbi}\t{kegg}\n")

        return mapping

    def _ncbi_id_to_symbol(
        self,
        gene_info_path: str,
        taxid: Optional[int] = None,
    ) -> Dict[str, str]:
        """Build an NCBI Gene ID-to-symbol mapping for one TaxID."""
        mapping: Dict[str, str] = {}
        opener = gzip.open if gene_info_path.endswith(".gz") else open

        with opener(gene_info_path, "rt", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 3 and (taxid is None or parts[0] == str(taxid)):
                    mapping[parts[1]] = parts[2]  # gene_id → symbol

        return mapping

    # ================================================================
    # File Generation
    # ================================================================

    def _build_gene2pathway(
        self,
        species: str,
        gene_pathway_links: Dict[str, List[str]],
        kegg_to_ncbi: Dict[str, str],
        ncbi_to_symbol: Dict[str, str],
        pathways: List[Tuple[str, str]],
    ) -> Path:
        """Write normalized gene-to-pathway memberships."""
        pathway_names = {pid: pname for pid, pname in pathways}
        out_file = self.cache_dir / f"{species}_gene2pathway.txt"

        n = 0
        with open(out_file, "w", encoding="utf-8") as f:
            for kegg_gene_id, pw_ids in sorted(gene_pathway_links.items()):
                # KEGG ID → NCBI ID → Symbol
                ncbi_id = kegg_to_ncbi.get(kegg_gene_id)
                if not ncbi_id:
                    continue
                symbol = ncbi_to_symbol.get(ncbi_id)
                if not symbol:
                    continue

                for pw_id in pw_ids:
                    pw_name = pathway_names.get(pw_id, pw_id)
                    # Remove the suffix " -Homo sapiens (human)" from the name
                    pw_name = self._clean_pathway_name(pw_name)
                    f.write(f"{symbol}\t{ncbi_id}\t{pw_id}\t{pw_name}\n")
                    n += 1

        print(f"|--- gene2pathway.txt: {n}Article Link")
        return out_file

    def _build_pathway_summary(
        self,
        species: str,
        pathways: List[Tuple[str, str]],
    ) -> Path:
        """Write pathway identifiers, names, and category metadata."""
        out_file = self.cache_dir / f"{species}_pathway_summary.txt"

        # Solve path level (from list/pathwaythe indentation of the
        # In the list output of KEGG, the row starts with a path ID without a number
        # Actually, list/pathway/{org}Only return to the access entry, not the classification level
        # Classification information requires a general name to extrapolate or use a global classification

        # Fetch real classified information from KEGAPI
        brite_categories = self._get_brite_categories(species, pathways)

        n = 0
        with open(out_file, "w", encoding="utf-8") as f:
            for pw_id, pw_name in pathways:
                pw_name_clean = self._clean_pathway_name(pw_name)
                url = f"https://www.kegg.jp/entry/{species}{pw_id}"

                # Retrieve access classifications (from KEGG API CLASS field)
                category, subcategory = brite_categories.get(pw_id, ("Uncategorized", "Uncategorized"))

                f.write(f"{category}\t{subcategory}\t{pw_id}\t{pw_name_clean}\t{url}\n")
                n += 1

        print(f"|--- Wrote {n} pathways to pathway_summary.txt.")
        return out_file

    @staticmethod
    def _clean_pathway_name(name: str) -> str:
        """Remove organism suffixes from a KEGG pathway name."""
        # Match " -Species name (common name)" after
        # For example: "Glycolysis / Gluconeogenesis - Homo sapiens (human)"
        cleaned = re.sub(r"\s*-\s*[\w\s]+\([^)]*\)\s*$", "", name)
        return cleaned.strip()

    def fetch_organism_list(self) -> List[Tuple[str, str, int, int]]:
        """Retrieve the KEGG organism catalogue for registry generation."""
        import re

        data = self._api_get("list/organism")
        result: List[Tuple[str, str, int, int]] = []

        for line in data.strip().split("\n"):
            cols = line.split("\t")
            if len(cols) < 6:
                continue
            kegg_code = cols[0]
            name = cols[1]
            taxid = int(cols[2])
            gene_count = int(cols[5])
            # Remove the brackets from the brackets (e.g. "Homo sapiens" - "Homo sapiens")
            latin_name = re.sub(r"\s*\(.*?\)\s*$", "", name).strip()
            result.append((kegg_code, latin_name, taxid, gene_count))

        return result

    def _get_brite_categories(self, species: str, pathways: List[Tuple[str, str]]) -> Dict[str, Tuple[str, str]]:
        """Retrieve KEGG BRITE pathway category assignments."""
        cache_file = self.cache_dir / f"{species}_pathway_classes.txt"

        # Read it directly if there is an existing cache
        if cache_file.exists() and not self.overwrite:
            print(f"|---Loading access from cache")
            categories: Dict[str, Tuple[str, str]] = {}
            pathway_ids = {pathway_id for pathway_id, _ in pathways}
            with open(cache_file, 'r', encoding='utf-8') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) == 3:
                        pw_id, cat, subcat = parts
                        bare_id = re.sub(r"^[A-Za-z]{3}(?=\d{5}$)", "", pw_id)
                        if bare_id in pathway_ids:
                            current = categories.get(bare_id)
                            if current is None or (
                                current[0] == "Uncategorized" and cat != "Uncategorized"
                            ):
                                categories[bare_id] = (cat, subcat)
            return categories

        print(f"|---Fetch access classification information from KEGAPI")
        categories: Dict[str, Tuple[str, str]] = {}

        try:
            # The KEGG API accepts at most 10 pathway entries per batch request.
            # Format: get/hsa04110+hsa00010
            pw_ids = [pw_id for pw_id, _ in pathways]
            batch_size = 10
            fetched = 0

            for i in range(0, len(pw_ids), batch_size):
                batch = pw_ids[i:i + batch_size]
                ids_str = '+'.join(
                    pw_id if pw_id.startswith(species) else f"{species}{pw_id}"
                    for pw_id in batch
                )
                data = self._api_get(f"get/{ids_str}")

                # Parsing CLASS fields
                # Batch Back Format: Multiple circuit data integration, each starting with ENTRY, with///Over.
                # ENTRY -> NAME -> CLASS -> ... -> ///
                lines = data.split('\n')

                for pw_id in batch:
                    found = False
                    in_entry = False
                    class_line = None

                    for line in lines:
                        # Confirm that the returned ENTRY matches the requested pathway.
                        if line.startswith("ENTRY") and pw_id in line:
                            in_entry = True
                            continue

                        # Stop at the KEGG record terminator.
                        if line.startswith("///"):
                            in_entry = False
                            continue

                        # Find CLASS in current access
                        if in_entry and line.startswith("CLASS"):
                            class_line = line.replace("CLASS", "").strip()
                            break

                    if class_line:
                        # Format: "Cellular Services; Cell Growth and Death"
                        parts = class_line.split("; ")
                        category = parts[0].replace(" ", "_") if parts else "Uncategorized"
                        subcategory = parts[1].replace(" ", "_") if len(parts) > 1 else "Uncategorized"
                        categories[pw_id] = (category, subcategory)
                    else:
                        # CLASS not found, using default classification
                        categories[pw_id] = ("Uncategorized", "Uncategorized")

                    fetched += 1
                    if fetched % 50 == 0:
                        print(f"|--- Retrieved categories for {fetched}/{len(pw_ids)} pathways")

                time.sleep(self.REQUEST_INTERVAL)

        except Exception as e:
            raise RuntimeError(f"KEGG Access Failed: {e}") from e

        # Save Cache
        with open(cache_file, 'w', encoding='utf-8') as f:
            for pw_id, (cat, subcat) in categories.items():
                f.write(f"{pw_id}\t{cat}\t{subcat}\n")

        print(f"|---Got it.{len(categories)}A traffic class.")
        return categories

    @staticmethod
    def _get_hardcoded_categories() -> Dict[str, Tuple[str, str]]:
        """Return fallback KEGG categories for known pathways."""
        # Hard-coded classification map (using hsa prefix)
        categories = {
            # Metabolic pathways
            "hsa00010": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00020": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00030": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00040": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00051": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00052": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00053": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00061": ("Metabolism", "Lipid_Metabolism"),
            "hsa00062": ("Metabolism", "Lipid_Metabolism"),
            "hsa00071": ("Metabolism", "Lipid_Metabolism"),
            "hsa00100": ("Metabolism", "Lipid_Metabolism"),
            "hsa00120": ("Metabolism", "Lipid_Metabolism"),
            "hsa00130": ("Metabolism", "Lipid_Metabolism"),
            "hsa00140": ("Metabolism", "Lipid_Metabolism"),
            "hsa00190": ("Metabolism", "Energy_Metabolism"),
            "hsa00220": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00230": ("Metabolism", "Nucleotide_Metabolism"),
            "hsa00232": ("Metabolism", "Metabolism_of_Cofactors"),
            "hsa00240": ("Metabolism", "Nucleotide_Metabolism"),
            "hsa00250": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00260": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00270": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00280": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00290": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00300": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00310": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00330": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00340": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00350": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00360": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00380": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00400": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00410": ("Metabolism", "Metabolism_of_Other_Amino_Acids"),
            "hsa00430": ("Metabolism", "Metabolism_of_Other_Amino_Acids"),
            "hsa00440": ("Metabolism", "Metabolism_of_Other_Amino_Acids"),
            "hsa00450": ("Metabolism", "Metabolism_of_Other_Amino_Acids"),
            "hsa00460": ("Metabolism", "Metabolism_of_Other_Amino_Acids"),
            "hsa00470": ("Metabolism", "Metabolism_of_Other_Amino_Acids"),
            "hsa00471": ("Metabolism", "Metabolism_of_Other_Amino_Acids"),
            "hsa00472": ("Metabolism", "Metabolism_of_Other_Amino_Acids"),
            "hsa00480": ("Metabolism", "Metabolism_of_Other_Amino_Acids"),
            "hsa00500": ("Metabolism", "Glycan_Biosynthesis"),
            "hsa00510": ("Metabolism", "Glycan_Biosynthesis"),
            "hsa00511": ("Metabolism", "Glycan_Biosynthesis"),
            "hsa00512": ("Metabolism", "Glycan_Biosynthesis"),
            "hsa00513": ("Metabolism", "Glycan_Biosynthesis"),
            "hsa00514": ("Metabolism", "Glycan_Biosynthesis"),
            "hsa00520": ("Metabolism", "Glycan_Biosynthesis"),
            "hsa00531": ("Metabolism", "Glycan_Biosynthesis"),
            "hsa00532": ("Metabolism", "Glycan_Biosynthesis"),
            "hsa00533": ("Metabolism", "Glycan_Biosynthesis"),
            "hsa00534": ("Metabolism", "Glycan_Biosynthesis"),
            "hsa00540": ("Metabolism", "Glycan_Biosynthesis"),
            "hsa00550": ("Metabolism", "Glycan_Biosynthesis"),
            "hsa00561": ("Metabolism", "Lipid_Metabolism"),
            "hsa00562": ("Metabolism", "Lipid_Metabolism"),
            "hsa00563": ("Metabolism", "Lipid_Metabolism"),
            "hsa00564": ("Metabolism", "Lipid_Metabolism"),
            "hsa00565": ("Metabolism", "Lipid_Metabolism"),
            "hsa00590": ("Metabolism", "Lipid_Metabolism"),
            "hsa00591": ("Metabolism", "Lipid_Metabolism"),
            "hsa00592": ("Metabolism", "Lipid_Metabolism"),
            "hsa00600": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00601": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00620": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00625": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00626": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00630": ("Metabolism", "One_Carbon_Metabolism"),
            "hsa00640": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00650": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00670": ("Metabolism", "One_Carbon_Metabolism"),
            "hsa00790": ("Metabolism", "Metabolism_of_Cofactors"),
            "hsa00830": ("Metabolism", "Metabolism_of_Cofactors"),
            "hsa00860": ("Metabolism", "Metabolism_of_Cofactors"),
            "hsa00900": ("Metabolism", "Terpenoid_Backbone_Biosynthesis"),
            "hsa00910": ("Metabolism", "Nitrogen_Metabolism"),
            "hsa00920": ("Metabolism", "Sulfur_Metabolism"),
            "hsa00970": ("Metabolism", "Translation"),
            "hsa00980": ("Metabolism", "Metabolism_of_Cofactors"),
            "hsa00982": ("Metabolism", "Xenobiotics_Biodegradation"),
            "hsa00983": ("Metabolism", "Xenobiotics_Biodegradation"),
            # Cell process
            "hsa04110": ("Cellular_Processes", "Cell_Cycle"),
            "hsa04111": ("Cellular_Processes", "Cell_Cycle"),
            "hsa04112": ("Cellular_Processes", "Cell_Cycle"),
            "hsa04113": ("Cellular_Processes", "Cell_Cycle"),
            "hsa04114": ("Cellular_Processes", "Cell_Cycle"),
            "hsa04115": ("Cellular_Processes", "Cell_Cycle"),
            "hsa04120": ("Cellular_Processes", "Cellular_Senescence"),
            "hsa04122": ("Cellular_Processes", "Transport"),
            "hsa04130": ("Cellular_Processes", "Folding_and_Degradation"),
            "hsa04140": ("Cellular_Processes", "Transport"),
            "hsa04141": ("Cellular_Processes", "Transport"),
            "hsa04142": ("Cellular_Processes", "Transport"),
            "hsa04144": ("Cellular_Processes", "Endocytosis"),
            "hsa04145": ("Cellular_Processes", "Phagocytosis"),
            "hsa04146": ("Cellular_Processes", "Autophagy"),
            "hsa04150": ("Cellular_Processes", "Signal_Transduction"),
            "hsa04151": ("Cellular_Processes", "Signal_Transduction"),
            # Signal transduction
            "hsa04010": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04012": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04014": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04015": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04020": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04022": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04024": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04066": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04068": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04070": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04071": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04072": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04010": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04210": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04211": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04213": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04215": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04217": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04218": ("Environmental_Information_Processing", "Signal_Transduction"),
            # DNA, copy and repair.
            "hsa03030": ("Genetic_Information_Processing", "Replication_and_Repair"),
            "hsa03040": ("Genetic_Information_Processing", "Replication_and_Repair"),
            "hsa03410": ("Genetic_Information_Processing", "Replication_and_Repair"),
            "hsa03420": ("Genetic_Information_Processing", "Replication_and_Repair"),
            "hsa03430": ("Genetic_Information_Processing", "Replication_and_Repair"),
            "hsa03440": ("Genetic_Information_Processing", "Replication_and_Repair"),
            "hsa03450": ("Genetic_Information_Processing", "Replication_and_Repair"),
            "hsa03460": ("Genetic_Information_Processing", "Replication_and_Repair"),
            # Transcription
            "hsa03020": ("Genetic_Information_Processing", "Transcription"),
            # Translation
            "hsa03010": ("Genetic_Information_Processing", "Translation"),
            "hsa03013": ("Genetic_Information_Processing", "Translation"),
            "hsa03015": ("Genetic_Information_Processing", "Translation"),
            # Collapse and degradation
            "hsa04130": ("Genetic_Information_Processing", "Folding_and_Degradation"),
            # Immunization system
            "hsa04612": ("Organismal_Systems", "Immune_System"),
            "hsa04620": ("Organismal_Systems", "Immune_System"),
            "hsa04621": ("Organismal_Systems", "Immune_System"),
            "hsa04622": ("Organismal_Systems", "Immune_System"),
            "hsa04623": ("Organismal_Systems", "Immune_System"),
            "hsa04625": ("Organismal_Systems", "Immune_System"),
            "hsa04630": ("Organismal_Systems", "Immune_System"),
            "hsa04640": ("Organismal_Systems", "Immune_System"),
            "hsa04650": ("Organismal_Systems", "Immune_System"),
            "hsa04611": ("Organismal_Systems", "Immune_System"),
            "hsa04657": ("Human_Diseases", "Infectious_Disease"),
            # Endocrine system
            "hsa04910": ("Organismal_Systems", "Endocrine_System"),
            "hsa04911": ("Organismal_Systems", "Endocrine_System"),
            "hsa04912": ("Organismal_Systems", "Endocrine_System"),
            "hsa04913": ("Organismal_Systems", "Endocrine_System"),
            "hsa04914": ("Organismal_Systems", "Endocrine_System"),
            "hsa04915": ("Organismal_Systems", "Endocrine_System"),
            "hsa04916": ("Organismal_Systems", "Endocrine_System"),
            "hsa04917": ("Organismal_Systems", "Endocrine_System"),
            "hsa04918": ("Organismal_Systems", "Endocrine_System"),
            "hsa04919": ("Organismal_Systems", "Endocrine_System"),
            "hsa04920": ("Organismal_Systems", "Endocrine_System"),
            "hsa04921": ("Organismal_Systems", "Endocrine_System"),
            "hsa04922": ("Organismal_Systems", "Endocrine_System"),
            "hsa04923": ("Organismal_Systems", "Endocrine_System"),
            "hsa04924": ("Organismal_Systems", "Endocrine_System"),
            "hsa04925": ("Organismal_Systems", "Endocrine_System"),
            "hsa04926": ("Organismal_Systems", "Endocrine_System"),
            "hsa04927": ("Organismal_Systems", "Endocrine_System"),
            "hsa04928": ("Organismal_Systems", "Endocrine_System"),
            "hsa04929": ("Organismal_Systems", "Endocrine_System"),
            "hsa04930": ("Organismal_Systems", "Endocrine_System"),
            "hsa04931": ("Organismal_Systems", "Endocrine_System"),
            "hsa04932": ("Organismal_Systems", "Endocrine_System"),
            "hsa04933": ("Organismal_Systems", "Endocrine_System"),
            "hsa04934": ("Organismal_Systems", "Endocrine_System"),
            "hsa04935": ("Organismal_Systems", "Endocrine_System"),
            # Imposing system
            "hsa04970": ("Organismal_Systems", "Digestive_System"),
            "hsa04971": ("Organismal_Systems", "Digestive_System"),
            "hsa04972": ("Organismal_Systems", "Digestive_System"),
            "hsa04973": ("Organismal_Systems", "Digestive_System"),
            "hsa04974": ("Organismal_Systems", "Digestive_System"),
            "hsa04975": ("Organismal_Systems", "Digestive_System"),
            "hsa04976": ("Organismal_Systems", "Digestive_System"),
            "hsa04977": ("Organismal_Systems", "Digestive_System"),
            "hsa04978": ("Organismal_Systems", "Digestive_System"),
            "hsa04979": ("Organismal_Systems", "Digestive_System"),
            # The nervous system.
            "hsa04710": ("Organismal_Systems", "Nervous_System"),
            "hsa04711": ("Organismal_Systems", "Nervous_System"),
            "hsa04720": ("Organismal_Systems", "Nervous_System"),
            "hsa04721": ("Organismal_Systems", "Nervous_System"),
            "hsa04722": ("Organismal_Systems", "Nervous_System"),
            "hsa04723": ("Organismal_Systems", "Nervous_System"),
            "hsa04724": ("Organismal_Systems", "Nervous_System"),
            "hsa04725": ("Organismal_Systems", "Nervous_System"),
            "hsa04726": ("Organismal_Systems", "Nervous_System"),
            "hsa04727": ("Organismal_Systems", "Nervous_System"),
            "hsa04728": ("Organismal_Systems", "Nervous_System"),
            "hsa04730": ("Organismal_Systems", "Nervous_System"),
            "hsa04740": ("Organismal_Systems", "Sensory_System"),
            "hsa04742": ("Organismal_Systems", "Sensory_System"),
            "hsa04744": ("Organismal_Systems", "Sensory_System"),
            "hsa04750": ("Organismal_Systems", "Sensory_System"),
            "hsa04710": ("Organismal_Systems", "Nervous_System"),
            # Epidemic diseases
            "hsa05160": ("Human_Diseases", "Infectious_Disease"),
            "hsa05161": ("Human_Diseases", "Infectious_Disease"),
            "hsa05162": ("Human_Diseases", "Infectious_Disease"),
            "hsa05164": ("Human_Diseases", "Infectious_Disease"),
            "hsa05165": ("Human_Diseases", "Infectious_Disease"),
            "hsa05166": ("Human_Diseases", "Infectious_Disease"),
            "hsa05167": ("Human_Diseases", "Infectious_Disease"),
            "hsa05168": ("Human_Diseases", "Infectious_Disease"),
            "hsa05169": ("Human_Diseases", "Infectious_Disease"),
            "hsa05170": ("Human_Diseases", "Infectious_Disease"),
            "hsa05171": ("Human_Diseases", "Infectious_Disease"),
            # Cancer
            "hsa05200": ("Human_Diseases", "Cancer"),
            "hsa05210": ("Human_Diseases", "Cancer"),
            "hsa05211": ("Human_Diseases", "Cancer"),
            "hsa05212": ("Human_Diseases", "Cancer"),
            "hsa05213": ("Human_Diseases", "Cancer"),
            "hsa05214": ("Human_Diseases", "Cancer"),
            "hsa05215": ("Human_Diseases", "Cancer"),
            "hsa05216": ("Human_Diseases", "Cancer"),
            "hsa05217": ("Human_Diseases", "Cancer"),
            "hsa05218": ("Human_Diseases", "Cancer"),
            "hsa05219": ("Human_Diseases", "Cancer"),
            "hsa05220": ("Human_Diseases", "Cancer"),
            "hsa05221": ("Human_Diseases", "Cancer"),
            "hsa05222": ("Human_Diseases", "Cancer"),
            "hsa05223": ("Human_Diseases", "Cancer"),
            "hsa05224": ("Human_Diseases", "Cancer"),
            "hsa05225": ("Human_Diseases", "Cancer"),
            "hsa05226": ("Human_Diseases", "Cancer"),
            "hsa05230": ("Human_Diseases", "Cancer"),
            "hsa05231": ("Human_Diseases", "Cancer"),
            "hsa05232": ("Human_Diseases", "Cancer"),
            "hsa05235": ("Human_Diseases", "Cancer"),
        }
        return categories
