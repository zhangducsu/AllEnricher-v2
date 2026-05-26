"""
GMT文件生成端到端测试

测试内容：
1. GMT文件是否正确生成
2. GMT格式是否符合规范
3. 基因池提取功能

Author: AllEnricher
Date: 2026-05-26
"""

import gzip
import json
import os
import unittest
from pathlib import Path

from allenricher.database.gmt_generator import GMTGenerator


class TestGMTGenerationE2E(unittest.TestCase):
    """GMT文件生成端到端测试类"""

    @classmethod
    def setUpClass(cls):
        """测试类初始化"""
        cls.base_dir = Path(__file__).parent.parent
        cls.organism_dir = cls.base_dir / "database" / "organism" / "v20260515" / "hsa"
        cls.test_data_dir = cls.base_dir / "test_data"
        cls.species = "hsa"

        # 确保test_data目录存在
        cls.test_data_dir.mkdir(exist_ok=True)

        # 初始化GMT生成器
        cls.generator = GMTGenerator(str(cls.organism_dir))

        # 预期生成的GMT文件
        cls.expected_gmt_files = {
            "GO": cls.organism_dir / f"{cls.species}.GO.gmt.gz",
            "KEGG": cls.organism_dir / f"{cls.species}.KEGG.gmt.gz",
            "Reactome": cls.organism_dir / f"{cls.species}.Reactome.gmt.gz",
            "DO": cls.organism_dir / f"{cls.species}.DO.gmt.gz",
        }

    def test_01_gmt_files_exist(self):
        """测试GMT文件是否正确生成且存在"""
        for db_name, filepath in self.expected_gmt_files.items():
            with self.subTest(database=db_name):
                self.assertTrue(
                    filepath.exists(),
                    f"{db_name} GMT文件不存在: {filepath}"
                )
                # 检查文件非空
                self.assertGreater(
                    filepath.stat().st_size,
                    0,
                    f"{db_name} GMT文件为空"
                )

    def test_02_gmt_format_valid(self):
        """测试GMT文件格式是否符合规范

        格式要求: pathway_name<TAB>description<TAB>gene1<TAB>gene2...
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
                        # 至少包含: pathway_name, description, 一个基因
                        self.assertGreaterEqual(
                            len(parts),
                            3,
                            f"{db_name} GMT文件格式错误，行至少需要3列: {line[:100]}"
                        )

                        # 验证pathway_name非空
                        self.assertTrue(
                            parts[0],
                            f"{db_name} GMT文件中pathway_name为空"
                        )

                        # 验证至少有一个基因
                        genes = parts[2:]
                        self.assertGreater(
                            len(genes),
                            0,
                            f"{db_name} GMT文件中通路 {parts[0]} 没有关联基因"
                        )

                        line_count += 1

                    # 确保文件有内容
                    self.assertGreater(
                        line_count,
                        0,
                        f"{db_name} GMT文件没有有效数据行"
                    )

    def test_03_gene_pool_extraction(self):
        """测试基因池提取功能"""
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

        # 验证基因池非空
        self.assertGreater(
            len(all_genes),
            0,
            "基因池为空"
        )

        # 打印统计信息
        print("\n=== 基因池统计 ===")
        print(f"所有数据库合并后的唯一基因总数: {len(all_genes)}")
        for db_name, stats in db_gene_stats.items():
            print(f"\n{db_name}:")
            print(f"  通路数量: {stats['pathway_count']}")
            print(f"  唯一基因数: {stats['unique_genes']}")
            print(f"  总基因关联数: {stats['total_associations']}")

        # 保存基因池到文件
        gene_pool_path = self.test_data_dir / "gmt_gene_pool.json"
        gene_pool_data = {
            "species": self.species,
            "database_version": "v20260515",
            "total_unique_genes": len(all_genes),
            "genes": sorted(list(all_genes))
        }
        with open(gene_pool_path, 'w', encoding='utf-8') as f:
            json.dump(gene_pool_data, f, indent=2, ensure_ascii=False)

        print(f"\n基因池已保存到: {gene_pool_path}")

    def test_04_generate_all_gmt_function(self):
        """测试generate_all_gmt函数"""
        # 调用生成函数
        results = self.generator.generate_all_gmt(self.species)

        # 验证返回结果包含预期的数据库
        expected_dbs = ["GO", "KEGG", "Reactome", "DO"]
        for db in expected_dbs:
            self.assertIn(
                db,
                results,
                f"generate_all_gmt返回结果中缺少 {db}"
            )
            # 验证路径存在
            self.assertTrue(
                Path(results[db]).exists(),
                f"{db} GMT文件路径不存在: {results[db]}"
            )

    def test_05_gmt_content_sample(self):
        """测试GMT文件内容样本"""
        for db_name, filepath in self.expected_gmt_files.items():
            with self.subTest(database=db_name):
                with gzip.open(filepath, 'rt', encoding='utf-8') as f:
                    lines = f.readlines()

                    # 获取前3行作为样本
                    sample_lines = []
                    for line in lines[:3]:
                        line = line.strip()
                        if line:
                            sample_lines.append(line)

                    # 验证样本存在
                    self.assertGreaterEqual(
                        len(sample_lines),
                        1,
                        f"{db_name} GMT文件没有足够的数据行"
                    )

                    # 打印样本
                    print(f"\n=== {db_name} GMT样本 (前3行) ===")
                    for i, sample in enumerate(sample_lines, 1):
                        truncated = sample[:150] + "..." if len(sample) > 150 else sample
                        print(f"  行{i}: {truncated}")

    def test_06_specific_pathway_content(self):
        """测试特定通路的内容完整性"""
        # 测试GO文件中的特定通路
        go_file = self.expected_gmt_files["GO"]
        with gzip.open(go_file, 'rt', encoding='utf-8') as f:
            found_go_terms = False
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3 and parts[0].startswith("GO:"):
                    found_go_terms = True
                    # 验证GO ID格式
                    self.assertRegex(
                        parts[0],
                        r'^GO:\d{7}$',
                        f"GO ID格式错误: {parts[0]}"
                    )
            self.assertTrue(found_go_terms, "GO文件中没有找到有效的GO term")

        # 测试KEGG文件中的特定通路
        kegg_file = self.expected_gmt_files["KEGG"]
        with gzip.open(kegg_file, 'rt', encoding='utf-8') as f:
            found_kegg_paths = False
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3 and parts[0].startswith("hsa"):
                    found_kegg_paths = True
                    # 验证KEGG pathway ID格式
                    self.assertRegex(
                        parts[0],
                        r'^hsa\d{5}$',
                        f"KEGG pathway ID格式错误: {parts[0]}"
                    )
            self.assertTrue(found_kegg_paths, "KEGG文件中没有找到有效的pathway")


class TestGMTGeneratorIntegration(unittest.TestCase):
    """GMT生成器集成测试"""

    def test_gmt_generator_initialization(self):
        """测试GMT生成器初始化"""
        base_dir = Path(__file__).parent.parent
        organism_dir = base_dir / "database" / "organism" / "v20260515" / "hsa"

        generator = GMTGenerator(str(organism_dir))
        self.assertEqual(
            str(generator.organism_dir),
            str(organism_dir)
        )

    def test_gmt_generator_with_nonexistent_dir(self):
        """测试GMT生成器处理不存在的目录"""
        nonexistent_dir = "/path/that/does/not/exist"
        generator = GMTGenerator(nonexistent_dir)

        # 尝试生成应该抛出FileNotFoundError
        with self.assertRaises(FileNotFoundError):
            generator.generate_go_gmt("hsa")


if __name__ == "__main__":
    # 设置测试输出
    unittest.main(verbosity=2)
