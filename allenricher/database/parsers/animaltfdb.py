"""Parse AnimalTFDB transcription-factor lists and orthology mappings."""

import gzip
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from collections import defaultdict
import logging

import pandas as pd

logger = logging.getLogger(__name__)


class AnimalTFDBParser:
    """Parse AnimalTFDB TF lists and orthology mappings into gene sets."""

    @staticmethod
    def parse_tf_list(tf_list_path: str) -> pd.DataFrame:
        """Parse an AnimalTFDB transcription-factor list."""
        df = pd.read_csv(tf_list_path, sep='\t', low_memory=False)

        # Standardized listing (remove possible spaces)
        df.columns = df.columns.str.strip()

        logger.info(f"List of AnimamalTFDB TF: {len(df)}TF")
        return df

    _external_id_cache: Dict[str, Dict[int, Dict[str, str]]] = {}

    @classmethod
    def load_external_id_symbol_maps(
        cls, gene_info_path: str, taxids: Set[int]
    ) -> Dict[int, Dict[str, str]]:
        """Load external gene identifiers for several TaxIDs in one pass."""
        cache_key = str(Path(gene_info_path).resolve())
        cached = cls._external_id_cache.setdefault(cache_key, {})
        missing = {int(taxid) for taxid in taxids} - set(cached)
        if not missing:
            return {taxid: cached[taxid] for taxid in taxids}
        loaded: Dict[int, Dict[str, str]] = {taxid: {} for taxid in missing}
        opener = gzip.open if str(gene_info_path).endswith('.gz') else open
        with opener(gene_info_path, 'rt', encoding='utf-8') as handle:
            for line in handle:
                if not line or line.startswith('#'):
                    continue
                parts = line.rstrip('\r\n').split('\t')
                if len(parts) < 6:
                    continue
                try:
                    taxid = int(parts[0])
                except ValueError:
                    continue
                if taxid not in missing:
                    continue
                mapping = loaded[taxid]
                symbol = parts[2].strip()
                if not symbol or symbol == '-':
                    continue
                mapping[parts[1].strip()] = symbol
                mapping[symbol] = symbol
                for cross_reference in parts[5].split('|'):
                    external_id = cross_reference.partition(':')[2].strip()
                    if external_id and external_id != '-':
                        mapping[external_id] = symbol
        cached.update(loaded)
        for taxid, mapping in loaded.items():
            logger.info("gene_info external ID mapping: taxid=%s, IDs=%d", taxid, len(mapping))
        return {taxid: cached[taxid] for taxid in taxids}

    @classmethod
    def load_external_id_symbols(cls, gene_info_path: str, taxid: int) -> Dict[str, str]:
        """Load external gene identifiers for one TaxID."""
        return cls.load_external_id_symbol_maps(gene_info_path, {taxid})[taxid]

    @staticmethod
    def parse_ortholog_to_human(
        ortholog_path: str,
        species_id_to_symbol: Optional[Dict[str, str]] = None,
        human_id_to_symbol: Optional[Dict[str, str]] = None,
    ) -> Dict[str, str]:
        """Parse target-species to human orthology relationships."""
        species_id_to_symbol = species_id_to_symbol or {}
        human_id_to_symbol = human_id_to_symbol or {}
        frame = pd.read_csv(ortholog_path, sep='\t', dtype=str, low_memory=False)
        frame.columns = frame.columns.str.strip()
        official_columns = {'Ensembl ID', 'Ortholog ID'}
        ortholog_map: Dict[str, str] = {}

        if official_columns.issubset(frame.columns):
            if not species_id_to_symbol or not human_id_to_symbol:
                raise ValueError(
                    "AnimalTFDB v4 ortholog files contain external IDs; "
                    "species and human gene_info mappings are required"
                )
            score_columns = [
                column for column in frame.columns
                if column == 'Coverage' or column.startswith('Coverage.')
                or column == 'Identity' or column.startswith('Identity.')
            ]
            scores = frame[score_columns].apply(pd.to_numeric, errors='coerce').fillna(0).sum(axis=1)
            ranked = frame.assign(_score=scores).sort_values('_score', ascending=False, kind='mergesort')
            unmapped_species = 0
            unmapped_human = 0
            for _, row in ranked.iterrows():
                species_id = str(row['Ensembl ID']).strip()
                human_id = str(row['Ortholog ID']).strip()
                species_gene = species_id_to_symbol.get(species_id)
                human_gene = human_id_to_symbol.get(human_id)
                if not species_gene:
                    unmapped_species += 1
                    continue
                if not human_gene:
                    unmapped_human += 1
                    continue
                ortholog_map.setdefault(species_gene, human_gene)
            logger.info(
                "AnimalTFDB ortholog ID conversion: mapped=%d, unmapped species rows=%d, "
                "unmapped human rows=%d",
                len(ortholog_map), unmapped_species, unmapped_human,
            )
        else:
            # Backward compatibility for historical two-column symbol files.
            legacy = pd.read_csv(ortholog_path, sep='\t', dtype=str, header=None, usecols=[0, 1])
            for species_gene, human_gene in legacy.itertuples(index=False, name=None):
                species_gene = str(species_gene).strip()
                human_gene = str(human_gene).strip()
                if species_gene and human_gene:
                    ortholog_map.setdefault(species_gene, human_gene)

        logger.info("Loaded %d AnimalTFDB ortholog mappings", len(ortholog_map))
        return ortholog_map

    @staticmethod
    def build_database(
        tf_list_path: str,
        ortholog_path: str,
        output_dir: str,
        species: str,
        valid_genes: Optional[Set[str]] = None,
        gene_info_path: Optional[str] = None,
        species_taxid: Optional[int] = None,
    ) -> Tuple[pd.DataFrame, Dict[str, str]]:
        """Build normalized database artifacts from parsed source data."""
        outdir = Path(output_dir)
        outdir.mkdir(parents=True, exist_ok=True)

        logger.info(f"AnimalTFDFBParser: Start building database (species={species})")

        # Parsing TF List
        tf_df = AnimalTFDBParser.parse_tf_list(tf_list_path)

        # Restrict targets to genes present in the species database when available.
        if valid_genes:
            tf_df = tf_df[tf_df['Symbol'].isin(valid_genes)]
            logger.info(f"Filters: {len(tf_df)}TF")

        species_id_to_symbol: Dict[str, str] = {}
        human_id_to_symbol: Dict[str, str] = {}
        if gene_info_path and species_taxid is not None:
            mappings = AnimalTFDBParser.load_external_id_symbol_maps(
                gene_info_path, {species_taxid, 9606}
            )
            species_id_to_symbol = mappings[species_taxid]
            human_id_to_symbol = mappings[9606]
        if {'Ensembl', 'Symbol'}.issubset(tf_df.columns):
            species_id_to_symbol.update(
                {
                    str(external_id).strip(): str(symbol).strip()
                    for external_id, symbol in zip(tf_df['Ensembl'], tf_df['Symbol'])
                    if str(external_id).strip() and str(symbol).strip()
                }
            )

        # The official v4 file contains external IDs, not gene symbols.
        ortholog_map = AnimalTFDBParser.parse_ortholog_to_human(
            ortholog_path,
            species_id_to_symbol=species_id_to_symbol,
            human_id_to_symbol=human_id_to_symbol,
        )

        # Save TF Information Table
        tf_file = outdir / f"{species}.AnimalTFDB_2tf.tab.gz"
        logger.info(f"Writing file: {tf_file}")

        with gzip.open(tf_file, 'wt') as f:
            f.write('\t'.join(tf_df.columns) + '\n')
            for _, row in tf_df.iterrows():
                f.write('\t'.join(str(v) for v in row.values) + '\n')

        # Save TF description file
        disc_file = outdir / f"{species}.AnimalTFDB_2disc.gz"
        logger.info(f"Writing file: {disc_file}")

        with gzip.open(disc_file, 'wt') as f:
            f.write("TF\tFamily\tEntrez_ID\tEnsembl\tsource\n")
            for _, row in tf_df.iterrows():
                symbol = row.get('Symbol', '')
                family = row.get('Family', 'Unknown')
                entrez = row.get('Entrez_ID', 'NA')
                ensembl = row.get('Ensembl', 'NA')
                f.write(f"{symbol}\t{family}\t{entrez}\t{ensembl}\tAnimalTFDB\n")

        # Save Homogenesis Map
        ortholog_file = outdir / f"{species}.AnimalTFDB_ortholog.gz"
        logger.info(f"Writing file: {ortholog_file}")

        with gzip.open(ortholog_file, 'wt') as f:
            f.write("Species_Gene\tHuman_Gene\n")
            for sp_gene, hu_gene in ortholog_map.items():
                f.write(f"{sp_gene}\t{hu_gene}\n")

        logger.info(f"AnimalTFDFBParser: Database build complete")
        return tf_df, ortholog_map
