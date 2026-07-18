"""Parse TRRUST regulatory interactions into TF-target gene sets."""

import gzip
import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)


class TRRUSTParser:
    """Build TF-target gene sets from TRRUST regulatory edges."""

    @staticmethod
    def parse_edges(tsv_path: str) -> List[Tuple[str, str, str, str]]:
        """Read TRRUST regulatory edges with direction and PubMed identifiers."""
        open_func = gzip.open if tsv_path.endswith('.gz') else open
        edges: List[Tuple[str, str, str, str]] = []
        with open_func(tsv_path, 'rt', encoding='utf-8') as handle:
            for line in handle:
                if not line.strip() or line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 3:
                    continue
                tf, target, mode = (part.strip() for part in parts[:3])
                pmid = parts[3].strip() if len(parts) > 3 else ""
                if tf and target:
                    edges.append((tf, target, mode, pmid))
        if not edges:
            raise ValueError("No valid TF-target interactions were found in the TRRUST file")
        return edges

    @staticmethod
    def parse_tsv(tsv_path: str) -> Tuple[
        Dict[str, Set[str]],
        Dict[str, Set[str]],
        Dict[str, str]
    ]:
        """Parse the source TSV into normalized regulatory gene sets."""
        # AutoSelect Open
        if tsv_path.endswith('.gz'):
            f_open = gzip.open(tsv_path, 'rt', encoding='utf-8')
        else:
            f_open = open(tsv_path, 'r', encoding='utf-8')

        tf_to_targets: Dict[str, Set[str]] = {}
        gene_to_tfs: Dict[str, Set[str]] = {}
        tf_mode_counts: Dict[str, Dict[str, int]] = {}  # {TF: {'activation': n, 'repression': m, 'unknown': k}}
        n = 0

        with f_open:
            for line in f_open:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 3:
                    continue

                tf = parts[0].strip()
                target = parts[1].strip()
                mode = parts[2].strip()

                if not tf or not target:
                    continue

                # Build TF-> towers map
                if tf not in tf_to_targets:
                    tf_to_targets[tf] = set()
                tf_to_targets[tf].add(target)

                # Build Gene-> TFs Reverse Map
                if target not in gene_to_tfs:
                    gene_to_tfs[target] = set()
                gene_to_tfs[target].add(tf)

                # TF-stating mode count
                if tf not in tf_mode_counts:
                    tf_mode_counts[tf] = {'activation': 0, 'repression': 0, 'unknown': 0}
                mode_lower = mode.lower()
                if mode_lower == 'activation':
                    tf_mode_counts[tf]['activation'] += 1
                elif mode_lower == 'repression':
                    tf_mode_counts[tf]['repression'] += 1
                else:
                    tf_mode_counts[tf]['unknown'] += 1

                n += 1

        if n == 0:
            raise ValueError(
                "No valid TF-target interactions were found in the TRRUST file"
            )

        # Determine the main TF mode of regulation
        tf_modes: Dict[str, str] = {}
        for tf, counts in tf_mode_counts.items():
            act = counts['activation']
            rep = counts['repression']
            unk = counts['unknown']

            if act > 0 and rep > 0:
                tf_modes[tf] = 'mixed'
            elif act > 0:
                tf_modes[tf] = 'activator'
            elif rep > 0:
                tf_modes[tf] = 'repressor'
            else:
                tf_modes[tf] = 'unknown'

        logger.info("TRRUSTParser: solves %d Article TF-target association", n)
        logger.info("TRRUSTPARSER: %d TF, %d TREE",
                    len(tf_to_targets), len(gene_to_tfs))

        return tf_to_targets, gene_to_tfs, tf_modes

    @staticmethod
    def build_database(tsv_path: str, output_dir: str, species: str,
                       valid_genes: Optional[Set[str]] = None) -> None:
        """Build normalized database artifacts from parsed source data."""
        outdir_path = Path(output_dir)
        outdir_path.mkdir(parents=True, exist_ok=True)

        logger.info("TRRUSTParser: Start building TRRUST databases (species=%s)", species)

        # Step 1: parsing TSV files
        edges = TRRUSTParser.parse_edges(tsv_path)
        tf_to_targets, gene_to_tfs, tf_modes = TRRUSTParser.parse_tsv(tsv_path)

        # Step 2: If a valid_genes is provided, filtering genes not in the pool
        if valid_genes is not None:
            logger.info("Filtering TRRUST targets against %d valid genes", len(valid_genes))

            # Filter TF: TF itself must be a valid gene
            filtered_tf_to_targets: Dict[str, Set[str]] = {}
            for tf, targets in tf_to_targets.items():
                if tf not in valid_genes:
                    continue
                filtered_targets = {t for t in targets if t in valid_genes}
                if filtered_targets:
                    filtered_tf_to_targets[tf] = filtered_targets

            # Rebuild Gene_to_tfs
            filtered_gene_to_tfs: Dict[str, Set[str]] = {}
            for tf, targets in filtered_tf_to_targets.items():
                for target in targets:
                    if target not in filtered_gene_to_tfs:
                        filtered_gene_to_tfs[target] = set()
                    filtered_gene_to_tfs[target].add(tf)

            # Filter tf_modes
            filtered_tf_modes: Dict[str, str] = {
                tf: mode for tf, mode in tf_modes.items()
                if tf in filtered_tf_to_targets
            }

            tf_to_targets = filtered_tf_to_targets
            gene_to_tfs = filtered_gene_to_tfs
            tf_modes = filtered_tf_modes
            edges = [
                edge for edge in edges
                if edge[0] in tf_to_targets and edge[1] in tf_to_targets[edge[0]]
            ]

            logger.info("Filtering retained %d TFs and %d target genes",
                        len(tf_to_targets), len(gene_to_tfs))

        if not tf_to_targets:
            raise ValueError(
                "No valid TF-target associations remained after filtering"
            )

        # gene2TF rows define the statistical background and therefore include
        # only genes that occur as targets. Regulator-only TFs remain as columns
        # and in the TF2target file.
        all_genes = set(gene_to_tfs.keys())
        sorted_tfs = sorted(tf_to_targets.keys())
        sorted_genes = sorted(all_genes)

        # Step 4: Write TF2target.tab.gz
        # Format: TF\\ttarget1\\ttarget2\\t... (0/1(Atlas)
        sorted_targets = sorted(gene_to_tfs.keys())
        tab_file = outdir_path / f"{species}.TF2target.tab.gz"
        logger.info("Writing file: %s", tab_file)

        with gzip.open(tab_file, 'wt', encoding='utf-8') as f:
            header = ["TF"] + sorted_targets
            f.write('\t'.join(header) + '\n')

            for tf in sorted_tfs:
                row = [tf]
                for target in sorted_targets:
                    if target in tf_to_targets[tf]:
                        row.append('1')
                    else:
                        row.append('0')
                f.write('\t'.join(row) + '\n')

        # Step 5: write gene2TF.tab.gz
        # Format: Gene\\tTF1\\tTF2\\t... (0/1(Atlas)
        gene2tf_file = outdir_path / f"{species}.gene2TF.tab.gz"
        logger.info("Writing file: %s", gene2tf_file)

        with gzip.open(gene2tf_file, 'wt', encoding='utf-8') as f:
            header = ["Gene"] + sorted_tfs
            f.write('\t'.join(header) + '\n')

            for gene in sorted_genes:
                row = [gene]
                for tf in sorted_tfs:
                    if tf in gene_to_tfs.get(gene, set()):
                        row.append('1')
                    else:
                        row.append('0')
                f.write('\t'.join(row) + '\n')

        # Step 6: Write the TF metadata and the regularized border table with the header.
        disc_file = outdir_path / f"{species}.TF2disc.gz"
        logger.info("Writing file: %s", disc_file)

        with gzip.open(disc_file, 'wt', encoding='utf-8') as f:
            f.write("Term_ID\tTerm_Name\tTF\tMode\tTarget_Set_Size\tSource\n")
            for tf in sorted_tfs:
                mode = tf_modes.get(tf, 'unknown')
                target_count = len(tf_to_targets[tf])
                f.write(
                    f"{tf}\t{tf} targets [TRRUST]\t{tf}\t{mode}\t{target_count}\tTRRUST\n"
                )

        edge_file = outdir_path / f"{species}.TRRUST_edges.tsv.gz"
        with gzip.open(edge_file, 'wt', encoding='utf-8') as f:
            f.write("TF\tTarget\tMode\tPMID\tSource\n")
            for tf, target, mode, pmid in sorted(set(edges)):
                f.write(f"{tf}\t{target}\t{mode}\t{pmid}\tTRRUST\n")

        logger.info("TRRUSTParser: TRRUST database built up")
