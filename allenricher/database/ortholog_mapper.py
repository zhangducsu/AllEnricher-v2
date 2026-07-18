"""Map human TF-target relationships to a target species through orthology."""

import gzip
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple
from collections import defaultdict
import logging

import pandas as pd

logger = logging.getLogger(__name__)


class OrthologMapper:
    """Map human regulatory targets to orthologs in another species."""

    # Exclude ambiguous many-to-one mappings by default.
    DEFAULT_DEDUP_STRATEGY = "one_to_one"

    def __init__(
        self,
        human_tf_to_targets: Dict[str, Set[str]],
        species_to_human: Dict[str, str],
        species_tf_set: Optional[Set[str]] = None,
        dedup_strategy: str = "one_to_one",
        human_term_metadata: Optional[Dict[str, Dict[str, str]]] = None,
    ):
        """Initialize the ortholog-mapping engine.

        Args:
            human_tf_to_targets: Human TF gene sets keyed by term ID.
            species_to_human: Target-species gene to human ortholog mapping.
            species_tf_set: Optional set of TFs annotated in the target species.
            dedup_strategy: Policy for ambiguous human-to-species mappings:
                ``one_to_one`` excludes ambiguous mappings, ``none`` and ``all``
                retain every mapping, and ``first`` keeps the alphabetically first
                target-species gene.
            human_term_metadata: Optional provenance for each human TF gene set.
        """
        self.human_tf_to_targets = human_tf_to_targets
        self.species_to_human = species_to_human
        self.species_tf_set = species_tf_set
        self.human_term_metadata = human_term_metadata or {}
        self.mapped_term_metadata: Dict[str, Dict[str, str]] = {}

        # Validate the ambiguity policy.
        valid_strategies = {"one_to_one", "none", "first", "all"}
        if dedup_strategy not in valid_strategies:
            logger.warning(
                "Unknown dedup_strategy '%s'; using '%s'",
                dedup_strategy,
                self.DEFAULT_DEDUP_STRATEGY,
            )
            dedup_strategy = self.DEFAULT_DEDUP_STRATEGY
        self.dedup_strategy = dedup_strategy

        # Build the reverse human-to-species ortholog map.
        self.human_to_species: Dict[str, Set[str]] = defaultdict(set)
        for sp_gene, hu_gene in species_to_human.items():
            self.human_to_species[hu_gene].add(sp_gene)

        # Map each target-species TF to its human ortholog.
        self.species_tf_to_human_tf: Dict[str, str] = {}
        if species_tf_set:
            for sp_tf in species_tf_set:
                if sp_tf in species_to_human:
                    human_tf = species_to_human[sp_tf]
                    if self.dedup_strategy != "one_to_one" or len(self.human_to_species[human_tf]) == 1:
                        self.species_tf_to_human_tf[sp_tf] = human_tf

    def get_duplicate_stats(self) -> Dict[str, Any]:
        """Return ambiguity statistics for the orthology mapping."""
        # Count human genes with more than one target-species ortholog.
        multi_mappings = {
            hu: len(sp_genes)
            for hu, sp_genes in self.human_to_species.items()
            if len(sp_genes) > 1
        }
        return {
            "total_human_genes": len(self.human_to_species),
            "multi_mapping_count": len(multi_mappings),
            "multi_mapping_genes": multi_mappings,
            "total_species_genes": len(self.species_to_human),
        }

    def _get_deduplicated_targets(
        self, human_target: str, sp_targets: Set[str]
    ) -> Set[str]:
        """Resolve mapped target genes according to the selected ambiguity policy."""
        if self.dedup_strategy == "one_to_one":
            return sp_targets if len(sp_targets) == 1 else set()
        if self.dedup_strategy == "none" or len(sp_targets) <= 1:
            return sp_targets
        elif self.dedup_strategy == "first":
            # Keep one deterministic representative.
            return {sorted(sp_targets)[0]}
        else:  # ``all`` retains every mapping and records ambiguity statistics.
            return sp_targets

    def map_tf_targets(self) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]]]:
        """Map human TF-target gene sets to target-species orthologs."""
        tf_to_targets: Dict[str, Set[str]] = defaultdict(set)
        gene_to_tfs: Dict[str, Set[str]] = defaultdict(set)

        mapped_tf_count = 0
        unmapped_tf_count = 0
        dedup_removed_count = 0

        # Retain examples of ambiguous mappings for diagnostic logs.
        duplicate_mappings: Dict[str, List[str]] = {}

        terms_by_tf: Dict[str, List[str]] = defaultdict(list)
        if self.human_term_metadata:
            for term_id, item in self.human_term_metadata.items():
                terms_by_tf[item.get("TF", "")].append(term_id)

        for sp_tf, human_tf in self.species_tf_to_human_tf.items():
            human_terms = terms_by_tf.get(human_tf, [human_tf])
            if not any(self.human_tf_to_targets.get(term_id) for term_id in human_terms):
                unmapped_tf_count += 1
                continue

            for human_term in human_terms:
                human_targets = self.human_tf_to_targets.get(human_term, set())
                source_meta = self.human_term_metadata.get(human_term, {})
                context = source_meta.get("Context", "all")
                sp_term = f"{sp_tf}|{context.replace('|', '/').replace(' ', '_')}"
                if not self.human_term_metadata:
                    sp_term = sp_tf

                # Project each human target to orthologs in the target species.
                for human_target in human_targets:
                    sp_targets = self.human_to_species.get(human_target, set())
                    if not sp_targets:
                        continue

                    original_count = len(sp_targets)
                    deduped_targets = self._get_deduplicated_targets(human_target, sp_targets)
                    dedup_removed_count += original_count - len(deduped_targets)

                    if self.dedup_strategy == "all" and len(sp_targets) > 1:
                        duplicate_mappings[human_target] = sorted(sp_targets)

                    for sp_target in deduped_targets:
                        if sp_target != sp_tf:
                            tf_to_targets[sp_term].add(sp_target)
                            gene_to_tfs[sp_target].add(sp_term)

                if tf_to_targets.get(sp_term):
                    self.mapped_term_metadata[sp_term] = {
                        "Term_Name": f"{sp_tf} [hTFtarget inferred; {context}]",
                        "TF": sp_tf,
                        "Library": "AnimalTFDB_hTFtarget",
                        "Context": context,
                        "Evidence_Type": source_meta.get("Evidence_Type", "ChIP-seq"),
                        "Inference_Type": "ortholog-inferred",
                        "Human_TF": human_tf,
                    }

            mapped_tf_count += 1

        # Report ambiguity handling without changing the generated gene sets.
        if self.dedup_strategy != "none":
            dup_stats = self.get_duplicate_stats()
            logger.info(
                "Ortholog ambiguity policy '%s': %d human genes have multiple "
                "target-species mappings; %d mappings were excluded",
                self.dedup_strategy,
                dup_stats["multi_mapping_count"],
                dedup_removed_count,
            )
            if self.dedup_strategy == "all" and duplicate_mappings:
                logger.info(
                    "Examples of retained ambiguous mappings: %s",
                    dict(list(duplicate_mappings.items())[:5]),
                )

        logger.info(
            "Ortholog projection complete: %d TFs mapped, %d TFs unmapped, "
            "%d TF gene sets, %d target genes",
            mapped_tf_count,
            unmapped_tf_count,
            len(tf_to_targets),
            len(gene_to_tfs),
        )

        return dict(tf_to_targets), dict(gene_to_tfs)

    @staticmethod
    def build_mapped_database(
        tf_to_targets: Dict[str, Set[str]],
        gene_to_tfs: Dict[str, Set[str]],
        species_tf_df: pd.DataFrame,
        output_dir: str,
        species: str,
        term_metadata: Optional[Dict[str, Dict[str, str]]] = None,
    ) -> None:
        """Write an orthology-mapped TF-target database."""
        outdir = Path(output_dir)
        outdir.mkdir(parents=True, exist_ok=True)

        all_genes = set(gene_to_tfs.keys())
        all_tfs = set(tf_to_targets.keys())
        gene_list = sorted(all_genes)
        tf_list = sorted(all_tfs)

        # Write the gene-by-TF membership matrix.
        gene2tf_file = outdir / f"{species}.AnimalTFDB_2gene.tab.gz"
        logger.info("Writing TF membership matrix: %s", gene2tf_file)

        with gzip.open(gene2tf_file, 'wt') as f:
            f.write('Gene\t' + '\t'.join(tf_list) + '\n')
            for gene in gene_list:
                regulating_tfs = gene_to_tfs.get(gene, set())
                values = ['1' if tf in regulating_tfs else '0' for tf in tf_list]
                f.write(gene + '\t' + '\t'.join(values) + '\n')

        # Write term descriptions and source provenance.
        disc_file = outdir / f"{species}.AnimalTFDB_mapped_2disc.gz"
        logger.info("Writing TF term metadata: %s", disc_file)

        # Index AnimalTFDB family annotations by TF symbol.
        tf_family_map = {}
        if species_tf_df is not None and 'Symbol' in species_tf_df.columns:
            for _, row in species_tf_df.iterrows():
                tf_family_map[row['Symbol']] = row.get('Family', 'Unknown')

        term_metadata = term_metadata or {}
        with gzip.open(disc_file, 'wt', encoding='utf-8') as f:
            f.write(
                "Term_ID\tTerm_Name\tTF\tLibrary\tContext\tEvidence_Type\t"
                "Inference_Type\tHuman_TF\tFamily\tTarget_Set_Size\tSource\n"
            )
            for term_id in tf_list:
                item = term_metadata.get(term_id, {})
                tf = item.get('TF', term_id)
                family = tf_family_map.get(tf, 'Unknown')
                f.write(
                    f"{term_id}\t{item.get('Term_Name', tf)}\t{tf}\t"
                    f"{item.get('Library', 'AnimalTFDB_hTFtarget')}\t"
                    f"{item.get('Context', 'all')}\t{item.get('Evidence_Type', 'ChIP-seq')}\t"
                    f"{item.get('Inference_Type', 'ortholog-inferred')}\t"
                    f"{item.get('Human_TF', '')}\t{family}\t{len(tf_to_targets[term_id])}\t"
                    "AnimalTFDB_ortholog_inferred_hTFtarget\n"
                )

        logger.info("Ortholog-mapped TF database build complete: %s", outdir)
