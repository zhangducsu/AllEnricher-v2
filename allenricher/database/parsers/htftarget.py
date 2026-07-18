"""Parse hTFtarget regulatory contexts into independent TF-target gene sets."""

import gzip
import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from collections import defaultdict
import logging

import pandas as pd

logger = logging.getLogger(__name__)


class HTFtargetParser:
    """Build context-specific hTFtarget gene sets."""

    @staticmethod
    def _term_id(tf: str, tissue: str) -> str:
        context = re.sub(r"\s+", "_", tissue.strip()).replace("|", "/")
        return f"{tf}|{context or 'unspecified'}"

    @staticmethod
    def parse_context_terms(tsv_path: str) -> Tuple[
        Dict[str, Set[str]], Dict[str, Set[str]], Dict[str, Dict[str, str]]
    ]:
        """Build independent TF-target gene sets for each regulatory context."""
        term_to_targets: Dict[str, Set[str]] = defaultdict(set)
        gene_to_terms: Dict[str, Set[str]] = defaultdict(set)
        metadata: Dict[str, Dict[str, str]] = {}

        count = 0
        with open(tsv_path, 'r', encoding='utf-8') as handle:
            next(handle, None)
            for line in handle:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 3:
                    continue
                tf, target, tissue_field = (part.strip() for part in parts[:3])
                if not tf or not target:
                    continue
                tissues = [item.strip() for item in tissue_field.split(',') if item.strip()]
                for tissue in tissues or ["unspecified"]:
                    term_id = HTFtargetParser._term_id(tf, tissue)
                    term_to_targets[term_id].add(target)
                    gene_to_terms[target].add(term_id)
                    metadata[term_id] = {
                        "Term_Name": f"{tf} [hTFtarget; {tissue}]",
                        "TF": tf,
                        "Library": "hTFtarget",
                        "Context": tissue,
                        "Evidence_Type": "ChIP-seq",
                        "Inference_Type": "direct",
                    }
                count += 1

        if not term_to_targets:
            raise ValueError("No valid TF-target associations were found in the hTFtarget file")
        logger.info(
            "hTFtarget: parsed %d associations across %d TF-tissue terms",
            count, len(term_to_targets),
        )
        return dict(term_to_targets), dict(gene_to_terms), metadata

    @staticmethod
    def parse_tsv(tsv_path: str) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]], Dict[str, Set[str]]]:
        """Parse the source TSV into normalized regulatory gene sets."""
        term_to_targets, _, metadata = HTFtargetParser.parse_context_terms(tsv_path)
        tf_to_targets: Dict[str, Set[str]] = defaultdict(set)
        gene_to_tfs: Dict[str, Set[str]] = defaultdict(set)
        tf_to_tissues: Dict[str, Set[str]] = defaultdict(set)
        for term_id, targets in term_to_targets.items():
            item = metadata[term_id]
            tf = item["TF"]
            tf_to_targets[tf].update(targets)
            tf_to_tissues[tf].add(item["Context"])
            for target in targets:
                gene_to_tfs[target].add(tf)

        logger.info(f"hTFtarget: {len(tf_to_targets)}One TF, {len(gene_to_tfs)}A target gene.")
        return dict(tf_to_targets), dict(gene_to_tfs), dict(tf_to_tissues)

    @staticmethod
    def build_database(tsv_path: str, output_dir: str, species: str,
                       valid_genes: Optional[Set[str]] = None) -> None:
        """Build normalized database artifacts from parsed source data."""
        outdir = Path(output_dir)
        outdir.mkdir(parents=True, exist_ok=True)

        logger.info(f"HTFtargetParser: Start building the database (species={species})")

        term_to_targets, gene_to_terms, metadata = HTFtargetParser.parse_context_terms(tsv_path)

        # Restrict targets to genes present in the species database when available.
        if valid_genes:
            term_to_targets = {
                term_id: {g for g in targets if g in valid_genes}
                for term_id, targets in term_to_targets.items()
            }
            term_to_targets = {
                term_id: targets for term_id, targets in term_to_targets.items() if targets
            }
            metadata = {term_id: metadata[term_id] for term_id in term_to_targets}
            gene_to_terms = defaultdict(set)
            for term_id, targets in term_to_targets.items():
                for target in targets:
                    gene_to_terms[target].add(term_id)
            logger.info(f"Filters: {len(term_to_targets)}TF-tassue term")

        # Get all the genes.
        all_genes = set()
        for targets in term_to_targets.values():
            all_genes.update(targets)
        all_terms = set(term_to_targets)

        # Build Gene x TF Matrix
        gene_list = sorted(all_genes)
        term_list = sorted(all_terms)

        gene2tf_file = outdir / f"{species}.hTF_2gene.tab.gz"
        logger.info(f"Writing file: {gene2tf_file}")

        with gzip.open(gene2tf_file, 'wt') as f:
            f.write('Gene\t' + '\t'.join(term_list) + '\n')
            for gene in gene_list:
                regulating_terms = gene_to_terms.get(gene, set())
                values = ['1' if term in regulating_terms else '0' for term in term_list]
                f.write(gene + '\t' + '\t'.join(values) + '\n')

        # Build TF description file
        disc_file = outdir / f"{species}.hTF_2disc.gz"
        logger.info(f"Writing file: {disc_file}")

        with gzip.open(disc_file, 'wt') as f:
            f.write(
                "Term_ID\tTerm_Name\tTF\tLibrary\tContext\tEvidence_Type\t"
                "Inference_Type\tTarget_Set_Size\tSource\n"
            )
            for term_id in term_list:
                item = metadata[term_id]
                f.write(
                    f"{term_id}\t{item['Term_Name']}\t{item['TF']}\t"
                    f"{item['Library']}\t{item['Context']}\t{item['Evidence_Type']}\t"
                    f"{item['Inference_Type']}\t{len(term_to_targets[term_id])}\thTFtarget\n"
                )

        logger.info(f"HTFtargetParser: Database built")
