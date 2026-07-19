"""End-to-end tests for GMT generation and gene-set extraction."""

import gzip
import json
import os
import unittest
from pathlib import Path

from allenricher.database.gmt_generator import GMTGenerator


class TestGMTGenerationE2E(unittest.TestCase):
    """Validate generated GMT files against the public format contract."""

    @classmethod
    def setUpClass(cls):
        """Test Class Initialization"""
        cls.base_dir = Path(__file__).parent.parent
        cls.organism_dir = cls.base_dir / "database" / "organism" / "v20260515" / "hsa"
        cls.test_data_dir = cls.base_dir / "test_data"
        cls.species = "hsa"

        # Ensure that test_data directory exists
        cls.test_data_dir.mkdir(exist_ok=True)

        # Initialize GMT Generator
        cls.generator = GMTGenerator(str(cls.organism_dir))

        # Expected GMT file generation
        cls.expected_gmt_files = {
            "GO": cls.organism_dir / f"{cls.species}.GO.gmt.gz",
            "KEGG": cls.organism_dir / f"{cls.species}.KEGG.gmt.gz",
            "Reactome": cls.organism_dir / f"{cls.species}.Reactome.gmt.gz",
            "DO": cls.organism_dir / f"{cls.species}.DO.gmt.gz",
        }

    def test_01_gmt_files_exist(self):
        """Test whether GMT files are correctly generated and exist"""
        for db_name, filepath in self.expected_gmt_files.items():
            with self.subTest(database=db_name):
                self.assertTrue(
                    filepath.exists(),
                    f"{db_name}The GMT file does not exist: {filepath}"
                )
                # Check file is not empty
                self.assertGreater(
                    filepath.stat().st_size,
                    0,
                    f"{db_name} GMT file empty"
                )

    def test_02_gmt_format_valid(self):
        """TestGMTDo you think the format of the document is in line with the norm?

        Format requirements: pathway_name<TAB>description<TAB>gene1<TAB>gene2...
        """
        for db_name, filepath in self.expected_gmt_files.items():
            with self.subTest(database=db_name):
                with gzip.open(filepath, 'rt', encoding='utf-8') as f:
                    line_count = 0
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue

                        parts = line.split('\t')
                        # At least include: path_name, decrition, a gene
                        self.assertGreaterEqual(
                            len(parts),
                            3,
                            f"{db_name}GMT file format error, line requires at least three columns: {line[: 100]}"
                        )

                        # Validate path_name is not empty
                        self.assertTrue(
                            parts[0],
                            f"{db_name}The GMT file is empty_pathway_name"
                        )

                        # Test at least one gene.
                        genes = parts[2:]
                        self.assertGreater(
                            len(genes),
                            0,
                            f"{db_name}Access to GMT files{parts[0]}No genes associated with it."
                        )

                        line_count += 1

                    # Ensure that the document contains content
                    self.assertGreater(
                        line_count,
                        0,
                        f"{db_name}GMT file does not have a valid data line"
                    )

    def test_03_gene_pool_extraction(self):
        """Test gene pool extraction"""
        all_genes = set()
        db_gene_stats = {}

        for db_name, filepath in self.expected_gmt_files.items():
            db_genes = set()
            pathway_count = 0
            total_gene_associations = 0

            with gzip.open(filepath, 'rt', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue

                    parts = line.split('\t')
                    if len(parts) < 3:
                        continue

                    pathway_count += 1
                    genes = parts[2:]
                    total_gene_associations += len(genes)
                    db_genes.update(genes)
                    all_genes.update(genes)

            db_gene_stats[db_name] = {
                "pathway_count": pathway_count,
                "unique_genes": len(db_genes),
                "total_associations": total_gene_associations
            }

        # Verify that gene pool is not empty
        self.assertGreater(
            len(all_genes),
            0,
            "The gene pool is empty."
        )

        # Print statistical information
        print("\n== sync, corrected by elderman ==")
        print(f"The only total number of genes for all databases combined: {len(all_genes)}")
        for db_name, stats in db_gene_stats.items():
            print(f"\n{db_name}:")
            print(f"Number of routes: {stats['pathway_count']}")
            print(f"Only number of genes: {stats['unique_genes']}")
            print(f"Total gene association: {stats['total_associations']}")

        # Save gene pool to file
        gene_pool_path = self.test_data_dir / "gmt_gene_pool.json"
        gene_pool_data = {
            "species": self.species,
            "database_version": "v20260515",
            "total_unique_genes": len(all_genes),
            "genes": sorted(list(all_genes))
        }
        with open(gene_pool_path, 'w', encoding='utf-8') as f:
            json.dump(gene_pool_data, f, indent=2, ensure_ascii=False)

        print(f"\nThe gene pool has been saved to: {gene_pool_path}")

    def test_04_generate_all_gmt_function(self):
        """Test the Generate_all_gmt function"""
        # Call Generation Functions
        results = self.generator.generate_all_gmt(self.species)

        # Verify return results include expected databases
        expected_dbs = ["GO", "KEGG", "Reactome", "DO"]
        for db in expected_dbs:
            self.assertIn(
                db,
                results,
                f"Missing in results from generate_all_gmt return{db}"
            )
            # Verify path exists
            self.assertTrue(
                Path(results[db]).exists(),
                f"{db}The GMT file path does not exist: {results[db]}"
            )

    def test_05_gmt_content_sample(self):
        """Test the GMT file contents"""
        for db_name, filepath in self.expected_gmt_files.items():
            with self.subTest(database=db_name):
                with gzip.open(filepath, 'rt', encoding='utf-8') as f:
                    lines = f.readlines()

                    # Get the first three lines as a sample
                    sample_lines = []
                    for line in lines[:3]:
                        line = line.strip()
                        if line:
                            sample_lines.append(line)

                    # Verify that the sample exists
                    self.assertGreaterEqual(
                        len(sample_lines),
                        1,
                        f"{db_name}GMT file does not have enough data lines"
                    )

                    # Print Sample
                    print(f"\n=== {db_name}GMT sample (first three lines) = = =")
                    for i, sample in enumerate(sample_lines, 1):
                        truncated = sample[:150] + "..." if len(sample) > 150 else sample
                        print(f"Okay.{i}: {truncated}")

    def test_06_specific_pathway_content(self):
        """Test the content integrity of a given route"""
        # Test specific access routes in Go files
        go_file = self.expected_gmt_files["GO"]
        with gzip.open(go_file, 'rt', encoding='utf-8') as f:
            found_go_terms = False
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3 and parts[0].startswith("GO:"):
                    found_go_terms = True
                    # Authenticate GO ID format
                    self.assertRegex(
                        parts[0],
                        r'^GO:\d{7}$',
                        f"GO ID format error: {parts[0]}"
                    )
            self.assertTrue(found_go_terms, "No valid GO term found in the GO file")

        # Test specific routes in the Kegg file
        kegg_file = self.expected_gmt_files["KEGG"]
        with gzip.open(kegg_file, 'rt', encoding='utf-8') as f:
            found_kegg_paths = False
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3 and parts[0].startswith("hsa"):
                    found_kegg_paths = True
                    # Validate KEG path ID format
                    self.assertRegex(
                        parts[0],
                        r'^hsa\d{5}$',
                        f"KEG path ID error: {parts[0]}"
                    )
            self.assertTrue(found_kegg_paths, "No valid pathway found in the Kegg file")


class TestGMTGeneratorIntegration(unittest.TestCase):
    """GMT Generator Integrated Test"""

    def test_gmt_generator_initialization(self):
        """Test GMT generator initialization"""
        base_dir = Path(__file__).parent.parent
        organism_dir = base_dir / "database" / "organism" / "v20260515" / "hsa"

        generator = GMTGenerator(str(organism_dir))
        self.assertEqual(
            str(generator.organism_dir),
            str(organism_dir)
        )

    def test_gmt_generator_with_nonexistent_dir(self):
        """Test GMT generator to handle non-existent directories"""
        nonexistent_dir = "/path/that/does/not/exist"
        generator = GMTGenerator(nonexistent_dir)

        # Try to generate a filenotfoundError
        with self.assertRaises(FileNotFoundError):
            generator.generate_go_gmt("hsa")


if __name__ == "__main__":
    # Set Test Output
    unittest.main(verbosity=2)
