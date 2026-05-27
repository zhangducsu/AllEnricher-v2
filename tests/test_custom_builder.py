"""
CustomDatabaseBuilder 单元测试

测试从不同格式的注释文件构建自定义数据库的完整流程。
"""

import gzip
import os
import tempfile
from datetime import datetime
from pathlib import Path

import pandas as pd
import pytest

from allenricher.database.custom_builder import CustomDatabaseBuilder


@pytest.fixture
def tmp_root(tmp_path):
    """创建临时数据库根目录"""
    return str(tmp_path / "database")


@pytest.fixture
def four_col_annotation(tmp_path):
    """创建四列注释文件: gene<TAB>term_id<TAB>term_name<TAB>hierarchy"""
    fpath = tmp_path / "four_col.txt"
    fpath.write_text(
        "GENE1\tTERM001\tCell Cycle\tBiology|Cell Biology|Cell Cycle\n"
        "GENE2\tTERM001\tCell Cycle\tBiology|Cell Biology|Cell Cycle\n"
        "GENE3\tTERM001\tCell Cycle\tBiology|Cell Biology|Cell Cycle\n"
        "GENE2\tTERM002\tApoptosis\tBiology|Cell Biology|Apoptosis\n"
        "GENE4\tTERM002\tApoptosis\tBiology|Cell Biology|Apoptosis\n"
        "GENE5\tTERM003\tMetabolism\tBiology|Metabolism\n"
        "GENE1\tTERM003\tMetabolism\tBiology|Metabolism\n"
    )
    return str(fpath)


@pytest.fixture
def three_col_annotation(tmp_path):
    """创建三列注释文件: gene<TAB>term_id<TAB>term_name"""
    fpath = tmp_path / "three_col.txt"
    fpath.write_text(
        "GENEA\tPATH_A\tPathway A\n"
        "GENEB\tPATH_A\tPathway A\n"
        "GENEB\tPATH_B\tPathway B\n"
        "GENEC\tPATH_B\tPathway B\n"
    )
    return str(fpath)


@pytest.fixture
def two_col_annotation(tmp_path):
    """创建两列注释文件: gene<TAB>term_id"""
    fpath = tmp_path / "two_col.txt"
    fpath.write_text(
        "G1\tT1\n"
        "G2\tT1\n"
        "G2\tT2\n"
        "G3\tT2\n"
    )
    return str(fpath)


@pytest.fixture
def empty_annotation(tmp_path):
    """创建空注释文件"""
    fpath = tmp_path / "empty.txt"
    fpath.write_text("")
    return str(fpath)


class TestBuildFromAnnotation:
    """build_from_annotation 核心流程测试"""

    def test_build_from_four_column_annotation(
        self, tmp_root, four_col_annotation
    ):
        """四列注释构建完整数据库"""
        builder = CustomDatabaseBuilder(root_dir=tmp_root)
        outdir = builder.build_from_annotation(
            annotation_file=four_col_annotation,
            species="hsa",
            taxid=9606,
            db_name="CustomDB"
        )

        # 验证返回路径
        assert Path(outdir).exists()
        assert outdir.endswith("hsa")

        # 验证三个输出文件都存在
        assert os.path.exists(os.path.join(outdir, "hsa.CustomDB2gene.tab.gz"))
        assert os.path.exists(os.path.join(outdir, "CustomDB2disc.gz"))
        assert os.path.exists(os.path.join(outdir, "hsa.CustomDB.gmt.gz"))

    def test_build_from_three_column_annotation(
        self, tmp_root, three_col_annotation
    ):
        """三列注释构建"""
        builder = CustomDatabaseBuilder(root_dir=tmp_root)
        outdir = builder.build_from_annotation(
            annotation_file=three_col_annotation,
            species="mmu",
            taxid=10090,
            db_name="MyPathway"
        )

        assert Path(outdir).exists()
        assert os.path.exists(os.path.join(outdir, "mmu.MyPathway2gene.tab.gz"))
        assert os.path.exists(os.path.join(outdir, "MyPathway2disc.gz"))
        assert os.path.exists(os.path.join(outdir, "mmu.MyPathway.gmt.gz"))

    def test_build_from_two_column_annotation(
        self, tmp_root, two_col_annotation
    ):
        """两列注释构建"""
        builder = CustomDatabaseBuilder(root_dir=tmp_root)
        outdir = builder.build_from_annotation(
            annotation_file=two_col_annotation,
            species="hsa",
            taxid=9606,
            db_name="SimpleDB"
        )

        assert Path(outdir).exists()
        assert os.path.exists(os.path.join(outdir, "hsa.SimpleDB2gene.tab.gz"))
        assert os.path.exists(os.path.join(outdir, "SimpleDB2disc.gz"))
        assert os.path.exists(os.path.join(outdir, "hsa.SimpleDB.gmt.gz"))

    def test_empty_annotation_file(self, tmp_root, empty_annotation):
        """空注释文件处理"""
        builder = CustomDatabaseBuilder(root_dir=tmp_root)
        with pytest.raises(ValueError, match="没有有效的基因-条目映射"):
            builder.build_from_annotation(
                annotation_file=empty_annotation,
                species="hsa",
                taxid=9606,
                db_name="EmptyDB"
            )


class TestGeneMatrix:
    """基因矩阵格式验证"""

    def test_gene_matrix_format(self, tmp_root, four_col_annotation):
        """验证矩阵格式（Gene列+0/1值）"""
        builder = CustomDatabaseBuilder(root_dir=tmp_root)
        outdir = builder.build_from_annotation(
            annotation_file=four_col_annotation,
            species="hsa",
            taxid=9606,
            db_name="MatrixTest"
        )

        matrix_path = os.path.join(outdir, "hsa.MatrixTest2gene.tab.gz")
        df = pd.read_csv(matrix_path, sep='\t', compression='gzip')

        # 验证列: Gene + 各 term 列
        assert "Gene" in df.columns
        term_cols = [c for c in df.columns if c != "Gene"]
        assert sorted(term_cols) == ["TERM001", "TERM002", "TERM003"]

        # 验证基因列表（去重有序，按出现顺序）
        assert sorted(df["Gene"].tolist()) == ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"]

        # 验证 0/1 值
        assert set(df["TERM001"].unique()).issubset({0, 1})
        assert set(df["TERM002"].unique()).issubset({0, 1})
        assert set(df["TERM003"].unique()).issubset({0, 1})

        # 验证具体映射正确性
        # TERM001: GENE1, GENE2, GENE3
        assert df.loc[df["Gene"] == "GENE1", "TERM001"].values[0] == 1
        assert df.loc[df["Gene"] == "GENE2", "TERM001"].values[0] == 1
        assert df.loc[df["Gene"] == "GENE3", "TERM001"].values[0] == 1
        assert df.loc[df["Gene"] == "GENE4", "TERM001"].values[0] == 0
        assert df.loc[df["Gene"] == "GENE5", "TERM001"].values[0] == 0

        # TERM002: GENE2, GENE4
        assert df.loc[df["Gene"] == "GENE2", "TERM002"].values[0] == 1
        assert df.loc[df["Gene"] == "GENE4", "TERM002"].values[0] == 1
        assert df.loc[df["Gene"] == "GENE1", "TERM002"].values[0] == 0


class TestGMTFile:
    """GMT 文件内容验证"""

    def test_gmt_auto_generated_content(self, tmp_root, four_col_annotation):
        """验证GMT内容正确"""
        builder = CustomDatabaseBuilder(root_dir=tmp_root)
        outdir = builder.build_from_annotation(
            annotation_file=four_col_annotation,
            species="hsa",
            taxid=9606,
            db_name="GMTTest"
        )

        gmt_path = os.path.join(outdir, "hsa.GMTTest.gmt.gz")
        with gzip.open(gmt_path, 'rt', encoding='utf-8') as f:
            lines = [l.strip() for l in f if l.strip()]

        # 应该有 3 个基因集（TERM001, TERM002, TERM003）
        assert len(lines) == 3

        # 解析每行
        gmt_data = {}
        for line in lines:
            parts = line.split('\t')
            term_id = parts[0]
            term_name = parts[1]
            genes = parts[2:]
            gmt_data[term_id] = (term_name, genes)

        # TERM001: Cell Cycle, genes: GENE1, GENE2, GENE3
        assert gmt_data["TERM001"][0] == "Cell Cycle"
        assert sorted(gmt_data["TERM001"][1]) == ["GENE1", "GENE2", "GENE3"]

        # TERM002: Apoptosis, genes: GENE2, GENE4
        assert gmt_data["TERM002"][0] == "Apoptosis"
        assert sorted(gmt_data["TERM002"][1]) == ["GENE2", "GENE4"]

        # TERM003: Metabolism, genes: GENE1, GENE5
        assert gmt_data["TERM003"][0] == "Metabolism"
        assert sorted(gmt_data["TERM003"][1]) == ["GENE1", "GENE5"]


class TestDescriptionFile:
    """描述文件验证"""

    def test_description_with_hierarchy(self, tmp_root, four_col_annotation):
        """验证描述文件含层级"""
        builder = CustomDatabaseBuilder(root_dir=tmp_root)
        outdir = builder.build_from_annotation(
            annotation_file=four_col_annotation,
            species="hsa",
            taxid=9606,
            db_name="DescTest"
        )

        disc_path = os.path.join(outdir, "DescTest2disc.gz")
        with gzip.open(disc_path, 'rt', encoding='utf-8') as f:
            lines = [l.strip() for l in f if l.strip()]

        assert len(lines) == 3

        disc_data = {}
        for line in lines:
            parts = line.split('\t')
            disc_data[parts[0]] = (parts[1], parts[2] if len(parts) > 2 else "")

        # 验证层级信息
        assert disc_data["TERM001"][1] == "Biology|Cell Biology|Cell Cycle"
        assert disc_data["TERM002"][1] == "Biology|Cell Biology|Apoptosis"
        assert disc_data["TERM003"][1] == "Biology|Metabolism"

    def test_description_without_hierarchy(self, tmp_root, two_col_annotation):
        """无层级时使用term_name"""
        builder = CustomDatabaseBuilder(root_dir=tmp_root)
        outdir = builder.build_from_annotation(
            annotation_file=two_col_annotation,
            species="hsa",
            taxid=9606,
            db_name="NoHierDB"
        )

        disc_path = os.path.join(outdir, "NoHierDB2disc.gz")
        with gzip.open(disc_path, 'rt', encoding='utf-8') as f:
            lines = [l.strip() for l in f if l.strip()]

        assert len(lines) == 2

        for line in lines:
            parts = line.split('\t')
            term_id = parts[0]
            term_name = parts[1]
            hierarchy = parts[2] if len(parts) > 2 else ""
            # 两列格式没有 term_name，所以 term_name == term_id
            # 层级应该回退为 term_name（即 term_id）
            assert hierarchy == term_id


class TestOutputStructure:
    """输出目录结构验证"""

    def test_output_directory_structure(self, tmp_root, four_col_annotation):
        """验证输出目录结构"""
        builder = CustomDatabaseBuilder(root_dir=tmp_root)
        outdir = builder.build_from_annotation(
            annotation_file=four_col_annotation,
            species="hsa",
            taxid=9606,
            db_name="StructTest"
        )

        outdir_path = Path(outdir)

        # 验证目录层级: database/organism/v{YYYYMMDD}/hsa/
        assert outdir_path.parent.name.startswith("v")
        assert outdir_path.name == "hsa"
        assert outdir_path.parent.parent.name == "organism"

        # 验证日期格式
        date_str = outdir_path.parent.name[1:]  # 去掉 'v' 前缀
        datetime.strptime(date_str, "%Y%m%d")  # 不抛异常即通过
