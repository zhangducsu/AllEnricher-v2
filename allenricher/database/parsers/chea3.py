"""Parse and normalize ChEA3 gene-set libraries and API responses."""

import gzip
import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)


class ChEA3Parser:
    """Normalize ChEA3 libraries without collapsing their provenance."""

    @staticmethod
    def parse_gmt(gmt_path: str, library_name: str = "unknown"
                  ) -> Tuple[Dict[str, Set[str]], Dict[str, str]]:
        """Parse gene sets from GMT format."""
        gmt_file = Path(gmt_path)
        if not gmt_file.exists():
            raise FileNotFoundError(f"[Error] The GMT file does not exist: {gmt_path}")

        tf_to_targets: Dict[str, Set[str]] = {}
        tf_descriptions: Dict[str, str] = {}

        open_func = gzip.open if gmt_path.endswith('.gz') else open
        mode = 'rt' if gmt_path.endswith('.gz') else 'r'

        with open_func(gmt_path, mode, encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) < 2:
                    continue

                tf_name = parts[0]
                # ChEA3 official document is "term"<TAB>gene1<TAB>Gene2...", no.
                # Standard GMT 's Description column. Old realization mischaracterizes the first target gene.
                description = tf_name
                targets = {target for target in parts[1:] if target}

                tf_to_targets[tf_name] = targets
                tf_descriptions[tf_name] = description

        logger.info("ChEA3Parser: parsing GMT [%s_ Completed, total%dTF",
                    library_name, len(tf_to_targets))

        return tf_to_targets, tf_descriptions

    @staticmethod
    def _library_name(gmt_path: str) -> str:
        name = Path(gmt_path).stem
        return name[:-3] if name.endswith("_tf") else name

    @staticmethod
    def _tf_name(term: str) -> str:
        return re.split(r"[_\s]", term, maxsplit=1)[0]

    @staticmethod
    def _display_name(tf: str, library: str, original_term: str) -> str:
        context = original_term[len(tf):].strip("_ ").replace("_", " ")
        if context.upper() in {"", "ARCHS4 PEARSON"}:
            return f"{tf} [{library}]"
        return f"{tf} [{library}; {context}]"

    @staticmethod
    def parse_api_result(api_result: Dict[str, List[Dict[str, str]]]
                        ) -> Dict[str, List[Dict[str, str]]]:
        """Normalize a ChEA3 API response."""
        standardized: Dict[str, List[Dict[str, str]]] = {}

        for lib_name, entries in api_result.items():
            standardized[lib_name] = []
            for entry in entries:
                standardized[lib_name].append({
                    'TF': str(entry.get('TF', '')),
                    'Rank': str(entry.get('Rank', '')),
                    'Pvalue': str(entry.get('Pvalue', '')),
                    'Overlap': str(entry.get('Overlap', '')),
                    'TargetCount': str(entry.get('TargetCount', '')),
                })

        return standardized

    @staticmethod
    def merge_libraries(libraries: Dict[str, Dict[str, Set[str]]],
                        method: str = "separate") -> Dict[str, Set[str]]:
        """Preserve each ChEA3 library as independent TF-target terms."""
        if method != "separate":
            raise ValueError(
                "ChEA3 libraries must remain separate; use method='separate'"
            )

        if not libraries:
            return {}

        merged = {
            f"{lib_name}|{term}": set(targets)
            for lib_name, tf_data in libraries.items()
            for term, targets in tf_data.items()
        }

        logger.info("ChEA3Parser: Merged complete (method=)%s(Same)%dTF",
                    method, len(merged))

        return merged

    @staticmethod
    def build_database(gmt_paths: List[str], output_dir: str, species: str,
                        merge_method: str = "separate",
                        valid_genes: Optional[Set[str]] = None) -> None:
        """Build normalized database artifacts from parsed source data."""
        outdir_path = Path(output_dir)
        outdir_path.mkdir(parents=True, exist_ok=True)

        logger.info("ChEA3Parser: Start building ChEA3 database (species=%s)", species)

        if merge_method != "separate":
            raise ValueError(
                "ChEA3 no longer supports cross-record integration of target genes; use merge_method='separate'"
            )

        # Each library and experiment term is independently maintained, avoiding the same name TF being combined with evidence sources.
        term_to_targets: Dict[str, Set[str]] = {}
        metadata: Dict[str, Dict[str, str]] = {}
        for gmt_path in gmt_paths:
            lib_name = ChEA3Parser._library_name(gmt_path)
            tf_to_targets, _ = ChEA3Parser.parse_gmt(gmt_path, lib_name)
            for original_term, targets in tf_to_targets.items():
                term_id = f"{lib_name}|{original_term}"
                tf = ChEA3Parser._tf_name(original_term)
                term_to_targets[term_id] = targets
                metadata[term_id] = {
                    "Term_Name": ChEA3Parser._display_name(tf, lib_name, original_term),
                    "TF": tf,
                    "Library": lib_name,
                    "Context": original_term,
                    "Evidence_Type": "coexpression" if "Coexpression" in lib_name
                    else "cooccurrence" if lib_name == "EnrichrQueries"
                    else "ChIP-seq",
                }

        if not term_to_targets:
            raise ValueError("[Error] No valid TF data is parsed!")

        # Filter for valid genes (if available)
        if valid_genes is not None:
            filtered_terms: Dict[str, Set[str]] = {}
            for term_id, targets in term_to_targets.items():
                filtered_targets = targets & valid_genes
                if filtered_targets:
                    filtered_terms[term_id] = filtered_targets
            term_to_targets = filtered_terms
            metadata = {term: metadata[term] for term in term_to_targets}
            logger.info("Gene filtering retained %d ChEA3 terms", len(term_to_targets))

        # Collect all genes and term
        all_genes: Set[str] = set()
        for targets in term_to_targets.values():
            all_genes.update(targets)

        sorted_genes = sorted(all_genes)
        sorted_terms = sorted(term_to_targets)

        # Step 5: Write in ChEA3_2gene.tab.gz
        # Format: Gene\tTF1\tTF2\t... (0/1(Atlas)
        tab_file = outdir_path / f"{species}.ChEA3_2gene.tab.gz"
        logger.info("Writing file: %s", tab_file)

        with gzip.open(tab_file, 'wt', encoding='utf-8') as f:
            header = ["Gene"] + sorted_terms
            f.write('\t'.join(header) + '\n')

            for gene in sorted_genes:
                row = [gene]
                for term_id in sorted_terms:
                    if gene in term_to_targets[term_id]:
                        row.append('1')
                    else:
                        row.append('0')
                f.write('\t'.join(row) + '\n')

        # Writes a description table with the source of the evidence.
        disc_file = outdir_path / f"{species}.ChEA3_2disc.gz"
        logger.info("Writing file: %s", disc_file)

        with gzip.open(disc_file, 'wt', encoding='utf-8') as f:
            f.write(
                "Term_ID\tTerm_Name\tTF\tLibrary\tContext\tEvidence_Type\t"
                "Target_Set_Size\tSource\n"
            )
            for term_id in sorted_terms:
                item = metadata[term_id]
                f.write(
                    f"{term_id}\t{item['Term_Name']}\t{item['TF']}\t"
                    f"{item['Library']}\t{item['Context']}\t{item['Evidence_Type']}\t"
                    f"{len(term_to_targets[term_id])}\tChEA3\n"
                )

        logger.info("Parsed %d ChEA3 TF gene sets covering %d genes", len(sorted_terms), len(sorted_genes))
        logger.info("ChEA3Parser: ChEA3 database built")
