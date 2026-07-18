"""Build Gene Ontology gene sets and term descriptions from NCBI and OBO data."""

import gzip
import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


class GOParser:
    """Build Gene Ontology matrices and term descriptions."""

    @staticmethod
    def _open_gz_or_text(filepath: str):
        """Open plain text or gzip-compressed input transparently."""
        if filepath.endswith('.gz'):
            return gzip.open(filepath, 'rt', encoding='utf-8')
        else:
            return open(filepath, 'r', encoding='utf-8')

    @staticmethod
    def parse_gene2go(gene2go_path: str, gene_info_path: str,
                      taxid: int, species: str, outdir: str) -> None:
        """Build GO memberships from NCBI gene2go and gene_info."""
        outdir_path = Path(outdir)
        outdir_path.mkdir(parents=True, exist_ok=True)

        print(f"|---Govarser: Start parsing gene2go (taxid=){taxid}, species={species})")

        # Step 1: Read gene_info.gz, create geneid-> symbol map
        # Preserve the v1 mapping rule: filter gene_info by TaxID and extract GeneID-symbol pairs.
        gene_id_to_symbol: Dict[str, str] = {}
        print(f"|---Read file: {gene_info_path}")

        with GOParser._open_gz_or_text(gene_info_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 3:
                    continue
                file_taxid = parts[0]
                gene_id = parts[1]
                symbol = parts[2]
                if int(file_taxid) == taxid:
                    gene_id_to_symbol[gene_id] = symbol

        print(f"|---Found it.{len(gene_id_to_symbol)}Genome (taxid={taxid})")

        # Step 2: read gene2go.gz, filter specified tabid
        # v1 Logical: Column [0] = taxid, [1] =geneid, [2] =go_id, [5] =go_name, [7] =category
        # 2026 Add new Qualifier column in new format: [0] = taxid, [1] =geneid, [2] =go_id,
        #   [3]=Evidence, [4]=Qualifier, [5]=GO_term, [6]=PubMed, [7]=Category
        # Compatibility in two formats: Auto-judgment by column number
        # Output gene2go.txt: symbol\\tgeneid\\tgo_id\\tcategory\\tgo_name
        gene2go_file = outdir_path / f"{species}.gene2go.txt"
        tab_file = outdir_path / f"{species}.GO2gene.tab.gz"

        all_go: Set[str] = set()       # All GO IDs
        all_symbols: Set[str] = set()  # All genes symbol
        tab: Dict[str, Dict[str, int]] = {}  # {symbol: {go_id: 1}}

        print(f"|---Read file: {gene2go_path}")
        n = 0

        with GOParser._open_gz_or_text(gene2go_path) as f_in, \
             open(gene2go_file, 'w', encoding='utf-8') as f_txt:
            for line in f_in:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 7:
                    continue
                file_taxid = parts[0]
                gene_id = parts[1]
                go_id = parts[2]

                # The current NCBI gene2go is 8 columns and the old format is 7 columns.
                if len(parts) >= 8:
                    qualifier = parts[4]
                    if 'NOT' in qualifier.split('|'):
                        continue
                    go_name = parts[5]
                    category = parts[7]
                else:
                    go_name = parts[4] if len(parts) > 4 else ""
                    category = parts[6] if len(parts) > 6 else ""

                if int(file_taxid) != taxid:
                    continue

                # Get agene symbol or use agene_id if not
                symbol = gene_id_to_symbol.get(gene_id, gene_id)

                all_symbols.add(symbol)
                all_go.add(go_id)

                # Write gene2go.txt
                f_txt.write(f"{symbol}\t{gene_id}\t{go_id}\t{category}\t{go_name}\n")

                # Build tab matrix data
                if symbol not in tab:
                    tab[symbol] = {}
                tab[symbol][go_id] = 1
                n += 1

        if n == 0:
            raise ValueError(
                f"[Error] No GO comment information for taxid={taxid} found in NCBI gene2go.gz file!"
            )

        print(f"|--- Parsed {n} GO term-gene associations.")

        # Step 3: Write GO2gene.tab.gz
        # v1 Logic: Header Gene\\tGO_ID1\\tGO_ID2\\t...
        # Data Rows symbol\\t0/1\\t0/1\\t...
        sorted_go = sorted(all_go)
        print(f"|---Writing file: {tab_file}")

        with gzip.open(tab_file, 'wt', encoding='utf-8') as f:
            # Write to table headers
            header = ["Gene"] + sorted_go
            f.write('\t'.join(header) + '\n')

            # Writing Data Line
            for symbol in sorted(all_symbols):
                row = [symbol]
                for go_id in sorted_go:
                    val = tab.get(symbol, {}).get(go_id, 0)
                    row.append(str(val))
                f.write('\t'.join(row) + '\n')

        print(f"|---Govarser: gene2go parsing complete")

    @staticmethod
    def parse_obo(obo_path: str, outdir: str) -> None:
        """Write Gene Ontology term names from go-basic.obo."""
        outdir_path = Path(outdir)
        outdir_path.mkdir(parents=True, exist_ok=True)

        print(f"|---Govarser: Start deciphering go-basic.obo")

        # v1 Logical: Split, extract id, name, namespace, is_a
        # Output format: GO_ID\\tnamespace: name\\tfather1; father2; ...
        disc_file = outdir_path / "GO2disc.gz"

        go_id_pattern = re.compile(r'^id:\s(GO:\d+)')
        name_pattern = re.compile(r'^name:\s(.*)')
        namespace_pattern = re.compile(r'^namespace:\s(.*)')
        is_a_pattern = re.compile(r'^is_a:\s(GO:\d+)')

        term_count = 0

        with open(obo_path, 'r', encoding='utf-8') as f_in, \
             gzip.open(disc_file, 'wt', encoding='utf-8') as f_out:

            content = f_in.read()
            # Click [Term] to split and skip the first empty block
            # v1 Sets record separator with $/ = '[Term]'
            blocks = content.split('[Term]')

            for block in blocks[1:]:  # Skip the first empty block
                # Remove [Typedef] and subsequent content
                block = re.sub(r'\[Typedef\].*$', '', block, flags=re.DOTALL)

                lines = block.strip().split('\n')
                go_id = None
                name = None
                namespace = None
                fathers = []

                for line in lines:
                    line = line.strip()
                    if not line:
                        continue

                    m = go_id_pattern.match(line)
                    if m:
                        go_id = m.group(1)
                        continue

                    m = name_pattern.match(line)
                    if m:
                        name = m.group(1)
                        continue

                    m = namespace_pattern.match(line)
                    if m:
                        namespace = m.group(1)
                        continue

                    m = is_a_pattern.match(line)
                    if m:
                        fathers.append(m.group(1))

                if go_id and name and namespace:
                    # Format: GO_ID\\tnamespace: name\\tfather1; father2; ...
                    father_str = ";".join(fathers) if fathers else ""
                    f_out.write(f"{go_id}\t{namespace}:{name}\t{father_str}\n")
                    term_count += 1

        print(f"|---Parsing{term_count}Go Term")
        print(f"|---Writing file: {disc_file}")
        print(f"|---GoParser: obo parsing complete")
