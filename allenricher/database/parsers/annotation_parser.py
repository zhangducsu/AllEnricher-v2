"""Parse user annotations with optional hierarchical term metadata."""

import gzip
from pathlib import Path
from typing import Dict, List, Optional, Set


class AnnotationRecord:
    """Store one parsed gene-to-term annotation."""

    def __init__(self, gene: str, term_id: str, term_name: str,
                 hierarchy: Optional[str] = None):
        self.gene = gene
        self.term_id = term_id
        self.term_name = term_name
        self.hierarchy = hierarchy

    @property
    def hierarchy_levels(self) -> List[str]:
        """Return non-empty hierarchy levels in order."""
        if not self.hierarchy:
            return []
        return self.hierarchy.split('|')

    def __repr__(self) -> str:
        return (
            f"AnnotationRecord(gene={self.gene!r}, term_id={self.term_id!r}, "
            f"term_name={self.term_name!r}, hierarchy={self.hierarchy!r})"
        )


class AnnotationParser:
    """Parse tabular annotations and optional hierarchy columns."""

    def __init__(self, file_path: Optional[str] = None, format_type: Optional[str] = None,
                 hierarchy_separator: str = '|', filepath: Optional[str] = None):
        # ``filepath`` is retained as a compatibility alias for ``file_path``.
        _path = filepath or file_path
        if _path is None:
            raise ValueError("Either file_path or filepath must be provided")
        self.file_path = Path(_path)
        self.format_type = format_type
        self.hierarchy_separator = hierarchy_separator
        self._records: Optional[List[AnnotationRecord]] = None

    def _detect_format(self) -> str:
        """Detect a supported annotation table layout."""
        if not self.file_path.exists():
            raise FileNotFoundError(f"Annotation file does not exist: {self.file_path}")

        opener = gzip.open if str(self.file_path).endswith('.gz') else open
        with opener(self.file_path, 'rt', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                col_count = len(line.split('\t'))
                if col_count >= 4:
                    return 'four_column'
                elif col_count == 3:
                    return 'three_column'
                elif col_count == 2:
                    return 'two_column'
                else:
                    raise ValueError(
                        f"Unsupported annotation format: the first data row has "
                        f"{col_count} columns; expected 2, 3, or at least 4"
                    )

        raise ValueError("The annotation file is empty or contains only headers and comments")

    def parse(self) -> List[AnnotationRecord]:
        """Parse the annotation table into normalized records."""
        if self._records is not None:
            return self._records

        if not self.file_path.exists():
            raise FileNotFoundError(f"Annotation file does not exist: {self.file_path}")

        fmt = self.format_type or self._detect_format()
        self._records = []

        opener = gzip.open if str(self.file_path).endswith('.gz') else open
        with opener(self.file_path, 'rt', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = [part.strip() for part in line.split('\t')]

                if parts[0].lower() in {'gene', 'gene_id', 'gene_symbol'}:
                    continue

                if fmt == 'four_column' and len(parts) >= 4:
                    gene = parts[0]
                    term_id = parts[1]
                    term_name = parts[2]
                    hierarchy = '|'.join(
                        level.strip()
                        for level in parts[3].split(self.hierarchy_separator)
                        if level.strip()
                    )
                    if not gene or not term_id or not term_name:
                        continue
                    self._records.append(
                        AnnotationRecord(gene, term_id, term_name, hierarchy)
                    )
                elif fmt == 'three_column' and len(parts) >= 3:
                    gene = parts[0]
                    term_id = parts[1]
                    term_name = parts[2]
                    if not gene or not term_id or not term_name:
                        continue
                    self._records.append(
                        AnnotationRecord(gene, term_id, term_name)
                    )
                elif fmt == 'two_column' and len(parts) >= 2:
                    gene = parts[0]
                    term = parts[1]
                    if not gene or not term:
                        continue
                    self._records.append(
                        AnnotationRecord(gene, term, term)
                    )

        return self._records

    def get_term_genes(self) -> Dict[str, Set[str]]:
        """Return term-to-gene memberships."""
        records = self.parse()
        term_genes: Dict[str, Set[str]] = {}
        for rec in records:
            if rec.term_id not in term_genes:
                term_genes[rec.term_id] = set()
            term_genes[rec.term_id].add(rec.gene)
        return term_genes

    def get_term_names(self) -> Dict[str, str]:
        """Return term identifiers mapped to descriptive names."""
        records = self.parse()
        term_names: Dict[str, str] = {}
        for rec in records:
            if rec.term_id not in term_names:
                term_names[rec.term_id] = rec.term_name
        return term_names

    def get_term_hierarchies(self) -> Dict[str, str]:
        """Return term identifiers mapped to hierarchy paths."""
        records = self.parse()
        term_hierarchies: Dict[str, str] = {}
        for rec in records:
            if rec.hierarchy and rec.term_id not in term_hierarchies:
                term_hierarchies[rec.term_id] = rec.hierarchy
        return term_hierarchies

    def get_hierarchy_tree(self) -> Dict:
        """Return the parsed annotation hierarchy as a nested tree."""
        records = self.parse()
        tree: Dict = {}

        # Press term_id for the polymer gene first.
        term_data: Dict[str, Dict] = {}
        for rec in records:
            if rec.term_id not in term_data:
                term_data[rec.term_id] = {
                    'genes': set(),
                    'term_name': rec.term_name,
                }
            term_data[rec.term_id]['genes'].add(rec.gene)

        # Build Tier Tree
        for rec in records:
            levels = rec.hierarchy_levels
            if not levels:
                # No hierarchical information. Put it on top.
                if rec.term_id not in tree:
                    tree[rec.term_id] = term_data[rec.term_id]
                continue

            # Build Trees From Layer to Layer
            current = tree
            for level in levels:
                if level not in current:
                    current[level] = {}
                current = current[level]

            # Place term data at the bottom
            if rec.term_id not in current:
                current[rec.term_id] = term_data[rec.term_id]

        return tree
