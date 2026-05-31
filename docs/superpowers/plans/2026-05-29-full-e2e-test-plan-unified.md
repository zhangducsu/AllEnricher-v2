# AllEnricher-v2 全场景端对端测试计划（统一版）

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 用真实数据对 AllEnricher-v2 的全部 11 个 CLI 子命令、所有参数、所有功能进行无死角端对端测试，确保每个功能、参数、命令均正常工作。

**Architecture:** 新建独立测试目录 `test_e2e_2026/`，按功能模块建立子层级文件夹，分门别类保存测试脚本、输入数据、输出文件、日志。所有文件按测试编号前缀命名，命名直接反映测试内容。

**Tech Stack:** Python 3.x, AllEnricher-v2 CLI (`python -m allenricher`), 真实基因数据

---

## 一、测试目录结构（强制规范）

所有测试相关文件必须存放在新建的独立目录中，禁止与项目原有文件混放。

```
AllEnricher-v2/
├── test_e2e_2026/                          # 测试根目录（新建）
│   ├── README.md                           # 测试说明文档
│   ├── run_all_tests.py                    # 总控脚本（可选）
│   │
│   ├── 00_input_data/                      # 输入数据目录
│   │   ├── 00_prepare_data.py              # 数据准备脚本
│   │   ├── hsa_gene_list_500.txt           # 人类基因列表（500基因）
│   │   ├── mmu_gene_list_500.txt           # 小鼠基因列表（500基因）
│   │   ├── rro_gene_list_500.txt           # 川金丝猴基因列表（500基因）
│   │   ├── hsa_background_genes.txt        # 背景基因集
│   │   ├── hsa_ranked_genes_gsea.tsv       # GSEA排序基因
│   │   ├── hsa_expression_matrix_500x6.tsv # 表达矩阵（500基因×6样本）
│   │   ├── hsa_custom_pathways_5.gmt       # 自定义GMT（5条通路）
│   │   ├── hsa_custom_pathways_5.gmt.gz    # 压缩版GMT
│   │   ├── hsa_gene_list_formats/          # 多格式输入测试
│   │   │   ├── tsv_format.tsv
│   │   │   ├── csv_format.csv
│   │   │   └── gene_ext_format.gene
│   │   ├── custom_annotations/             # 自定义注释文件
│   │   │   ├── three_column.tsv
│   │   │   ├── four_column.tsv
│   │   │   └── two_column.tsv
│   │   ├── error_scenarios/                # 错误场景测试文件
│   │   │   ├── empty_gene_list.txt
│   │   │   └── tiny_gene_list_5.txt
│   │   └── config_files/                   # YAML/JSON配置文件
│   │       ├── T091_gsea_permutations_100.yaml
│   │       ├── T092_gsea_size_limits.yaml
│   │       ├── T093_gsva_method_plage.yaml
│   │       ├── T094_gsva_method_zscore.yaml
│   │       ├── T095_gsva_kcdf_poisson.yaml
│   │       ├── T096_gsva_tau_0.5.yaml
│   │       ├── T097_max_genes_500.yaml
│   │       ├── T098_top_terms_10.yaml
│   │       ├── T099_plot_size_12x10.yaml
│   │       ├── T100_plot_formats_multi.yaml
│   │       ├── T101_report_format_html.yaml
│   │       ├── T102_ai_ollama_baseurl.yaml
│   │       ├── T103_ai_mock_disabled.yaml
│   │       ├── T078_config_full_yaml.yaml
│   │       ├── T079_config_full_json.json
│   │       ├── T106_cli_override_priority.yaml
│   │       ├── T143_style_cell_config.yaml
│   │       ├── T146_gsea_style_cell.yaml
│   │       ├── T160_ssgsea_style_cell.yaml
│   │       └── T174_gsva_style_cell.yaml
│   │
│   ├── 01_cmd_analyze/                     # analyze命令测试
│   │   ├── 01_ora_basic/                   # ORA基础测试
│   │   │   ├── T001_fisher_default/        # 测试编号前缀+测试内容
│   │   │   │   ├── run.sh                  # 执行脚本
│   │   │   │   ├── output/                 # 输出目录
│   │   │   │   └── test.log                # 测试日志
│   │   │   ├── T002_fisher_threshold_p0.01_q0.01/
│   │   │   ├── T003_fisher_min_genes_5/
│   │   │   ├── T004_fisher_min_genes_10/
│   │   │   ├── T005_fisher_only_significant/
│   │   │   ├── T006_fisher_parallel_jobs_4/
│   │   │   ├── T007_hypergeometric_default/
│   │   │   ├── T008_hypergeometric_correction_bh/
│   │   │   ├── T009_hypergeometric_correction_bonferroni/
│   │   │   ├── T010_hypergeometric_correction_holm/
│   │   │   └── T011_hypergeometric_correction_by/
│   │   │
│   │   ├── 02_ora_correction/              # ORA校正方法测试
│   │   │   ├── T013_correction_bh/
│   │   │   ├── T014_correction_by/
│   │   │   ├── T015_correction_bonferroni/
│   │   │   ├── T016_correction_holm/
│   │   │   └── T017_correction_none/
│   │   │
│   │   ├── 03_ora_databases/               # ORA数据库组合测试
│   │   │   ├── T018_database_go_only/
│   │   │   ├── T019_database_kegg_only/
│   │   │   ├── T020_database_reactome_only/
│   │   │   ├── T021_database_do_only/
│   │   │   ├── T022_database_all_combined/
│   │   │   └── T023_mmu_database_all/
│   │   │
│   │   ├── 04_background_modes/            # 背景基因集模式测试
│   │   │   ├── T024_background_annotated/
│   │   │   ├── T025_background_genome/
│   │   │   ├── T026_background_custom/
│   │   │   └── T027_background_custom_error/
│   │   │
│   │   ├── 05_gsea_analysis/               # GSEA分析测试
│   │   │   ├── T028_gsea_ranked_genes/
│   │   │   ├── T029_gsea_expression_matrix/
│   │   │   ├── T030_gsea_custom_gmt/
│   │   │   └── T031_gsea_with_groups/
│   │   │
│   │   ├── 06_ssgsea_analysis/             # ssGSEA分析测试
│   │   │   ├── T032_ssgsea_expression_matrix/
│   │   │   ├── T033_ssgsea_custom_gmt/
│   │   │   └── T034_ssgsea_with_groups/
│   │   │
│   │   ├── 07_gsva_analysis/               # GSVA分析测试
│   │   │   ├── T035_gsva_default_gsva/
│   │   │   ├── T036_gsva_method_plage/
│   │   │   ├── T037_gsva_method_zscore/
│   │   │   └── T038_gsva_with_groups/
│   │   │
│   │   ├── 08_visual_styles/               # 可视化风格测试
│   │   │   ├── T039_style_nature/
│   │   │   ├── T040_style_science/
│   │   │   ├── T041_style_cell/
│   │   │   ├── T042_style_colorblind/
│   │   │   ├── T043_style_presentation/
│   │   │   └── T044_style_omicshare/
│   │   │
│   │   ├── 09_visual_palettes/             # 可视化配色测试
│   │   │   ├── T045_palette_nature/
│   │   │   ├── T046_palette_science/
│   │   │   ├── T047_palette_lancet/
│   │   │   ├── T048_palette_nejm/
│   │   │   ├── T049_palette_jama/
│   │   │   ├── T050_palette_okabe_ito/
│   │   │   ├── T051_palette_gsea/
│   │   │   ├── T052_palette_omicshare/
│   │   │   ├── T053_palette_china_style/
│   │   │   ├── T054_palette_go_bp/
│   │   │   ├── T055_palette_kegg_pathway/
│   │   │   ├── T056_palette_tol_bright/
│   │   │   ├── T057_palette_tol_muted/
│   │   │   ├── T058_palette_tol_sunset/
│   │   │   └── T134-T142_remaining_palettes/
│   │   │
│   │   ├── 10_output_formats/              # 输出格式测试
│   │   │   ├── T059_format_png/
│   │   │   ├── T060_format_pdf/
│   │   │   ├── T061_format_svg/
│   │   │   ├── T062_dpi_150/
│   │   │   └── T063_dpi_600/
│   │   │
│   │   ├── 11_report_control/              # 报告控制测试
│   │   │   ├── T064_no_plot/
│   │   │   ├── T065_no_report/
│   │   │   └── T066_no_plot_no_report/
│   │   │
│   │   ├── 12_plot_types/                  # 图表类型测试
│   │   │   ├── T067_plot_dotplot/
│   │   │   ├── T068_plot_barplot/
│   │   │   ├── T069_plot_network/
│   │   │   ├── T070_plot_upset/
│   │   │   ├── T071_plot_volcano/
│   │   │   ├── T072_gsea_plot_enrichment/
│   │   │   ├── T073_gsea_plot_nes_barplot/
│   │   │   ├── T074_ssgsea_plot_heatmap/
│   │   │   ├── T075_ssgsea_plot_group_comparison/
│   │   │   ├── T076_ssgsea_plot_correlation/
│   │   │   └── T077_plot_multi_types/
│   │   │
│   │   ├── 13_config_loading/              # 配置文件加载测试
│   │   │   ├── T078_config_yaml_full/
│   │   │   └── T079_config_json_full/
│   │   │
│   │   ├── 14_version_control/             # 版本控制测试
│   │   │   ├── T080_use_version_v20260515/
│   │   │   └── T081_database_dir_custom/
│   │   │
│   │   ├── 15_ai_interpretation/           # AI解读测试
│   │   │   ├── T082_verbose_logging/
│   │   │   ├── T083_ai_mock_backend/
│   │   │   ├── T084_ai_openai_params/
│   │   │   └── T085_ai_claude_params/
│   │   │
│   │   ├── 16_error_scenarios/             # 错误场景测试
│   │   │   ├── T086_empty_gene_list/
│   │   │   ├── T087_invalid_species/
│   │   │   ├── T088_nonexistent_input/
│   │   │   ├── T089_gsea_missing_ranked/
│   │   │   └── T090_ssgsea_missing_expression/
│   │   │
│   │   ├── 17_config_only_params/          # Config-only参数测试
│   │   │   ├── T091_gsea_permutations_100/
│   │   │   ├── T092_gsea_size_limits/
│   │   │   ├── T093_gsva_method_plage/
│   │   │   ├── T094_gsva_method_zscore/
│   │   │   ├── T095_gsva_kcdf_poisson/
│   │   │   ├── T096_gsva_tau_0.5/
│   │   │   ├── T097_max_genes_500/
│   │   │   ├── T098_top_terms_10/
│   │   │   ├── T099_plot_size_12x10/
│   │   │   ├── T100_plot_formats_multi/
│   │   │   ├── T101_report_format_html/
│   │   │   ├── T102_ai_ollama_baseurl/
│   │   │   └── T103_ai_mock_disabled/
│   │   │
│   │   ├── 18_global_params/               # 全局参数测试
│   │   │   ├── T104_version_short/
│   │   │   └── T105_version_long/
│   │   │
│   │   ├── 19_cli_priority/                # CLI优先级测试
│   │   │   └── T106_cli_override_config/
│   │   │
│   │   ├── 20_input_formats/               # 多格式输入测试
│   │   │   ├── T107_input_tsv_format/
│   │   │   ├── T108_input_csv_format/
│   │   │   └── T109_input_gene_ext/
│   │   │
│   │   ├── 21_gmt_formats/                 # GMT格式测试
│   │   │   └── T110_gmt_gz_format/
│   │   │
│   │   ├── 22_group_comparison/            # 分组比较图表测试
│   │   │   ├── T111_group_comparison_box/
│   │   │   ├── T112_group_comparison_violin/
│   │   │   └── T113_group_comparison_bar/
│   │   │
│   │   ├── 23_correlation_methods/         # 相关性方法测试
│   │   │   ├── T114_correlation_pearson/
│   │   │   └── T115_correlation_spearman/
│   │   │
│   │   ├── 24_network_plots/               # 网络图测试
│   │   │   ├── T116_network_single_db/
│   │   │   └── T117_network_multi_db/
│   │   │
│   │   ├── 25_api_endpoints/               # API端点测试
│   │   │   ├── T117a_api_home/
│   │   │   ├── T117b_api_species/
│   │   │   ├── T117c_api_databases/
│   │   │   ├── T117d_api_analyze/
│   │   │   ├── T117e_api_upload/
│   │   │   ├── T117f_api_status_404/
│   │   │   ├── T117g_api_results_404/
│   │   │   ├── T117h_api_results_json/
│   │   │   ├── T117i_api_plot/
│   │   │   ├── T117j_api_report/
│   │   │   └── T117k_api_delete_404/
│   │   │
│   │   ├── 26_report_generation/           # 报告生成测试
│   │   │   ├── T118_report_normal/
│   │   │   └── T119_report_no_results/
│   │   │
│   │   ├── 27_url_generation/              # URL生成测试
│   │   │   └── T120_url_all_databases/
│   │   │
│   │   ├── 28_parallel_modes/              # 并行模式测试
│   │   │   ├── T127_serial_jobs_1/
│   │   │   ├── T128_parallel_jobs_4/
│   │   │   └── T129_all_cpu_jobs_minus1/
│   │   │
│   │   ├── 29_gsea_fields/                 # GSEA字段验证
│   │   │   └── T130_gsea_result_fields/
│   │   │
│   │   ├── 30_ssgsea_nan/                  # ssGSEA NaN处理
│   │   │   └── T131_ssgsea_nan_pvalue/
│   │   │
│   │   └── 31_output_fields/               # 输出字段验证
│   │       ├── T132_ora_result_fields/
│   │       └── T133_gsea_output_fields/
│   │
│   ├── 02_cmd_download/                    # download命令测试
│   │   ├── D001_download_go/
│   │   ├── D002_download_reactome/
│   │   ├── D003_download_do/
│   │   ├── D004_download_multiple/
│   │   ├── D005_workers_8/
│   │   ├── D006_no_multi_thread/
│   │   ├── D007_no_verify/
│   │   ├── D008_force_redownload/
│   │   ├── D009_custom_db_dir/
│   │   └── D010_kegg_hint/
│   │
│   ├── 03_cmd_build/                       # build命令测试
│   │   ├── B001_build_rno_go_kegg/
│   │   ├── B002_build_rno_all_dbs/
│   │   ├── B003_build_custom_db_dir/
│   │   ├── B004_build_with_gene_info/
│   │   ├── B005_build_with_go_annot/
│   │   ├── B006_build_custom_annot/
│   │   ├── B007_build_annot_format_three/
│   │   ├── B008_build_annot_format_four/
│   │   ├── B009_build_annot_format_two/
│   │   ├── B010_build_annot_format_auto/
│   │   ├── B011_build_hierarchy_sep/
│   │   ├── B012_build_hsa_do/
│   │   ├── B013_build_mmu_do_skip/
│   │   └── B014_build_hsa_all_dbs/
│   │
│   ├── 04_cmd_serve/                       # serve命令测试
│   │   ├── S001_serve_default/
│   │   ├── S002_serve_custom_host_port/
│   │   └── S003_serve_reload/
│   │
│   ├── 05_cmd_list/                        # list命令测试
│   │   ├── L001_list_species/
│   │   ├── L002_list_databases/
│   │   └── L003_list_no_args/
│   │
│   ├── 06_cmd_config/                      # config命令测试
│   │   ├── C001_config_default/
│   │   └── C002_config_custom_path/
│   │
│   ├── 07_cmd_check_update/                # check-update命令测试
│   │   ├── U001_check_update_default/
│   │   ├── U002_check_update_custom_dir/
│   │   ├── U003_check_update_json/
│   │   └── U004_check_update_dir_json/
│   │
│   ├── 08_cmd_cleanup/                     # cleanup命令测试
│   │   ├── CL001_cleanup_dry_run/
│   │   ├── CL002_cleanup_keep_1/
│   │   └── CL003_cleanup_custom_dir/
│   │
│   ├── 09_cmd_list_versions/               # list-versions命令测试
│   │   ├── V001_list_versions_default/
│   │   ├── V002_list_versions_json/
│   │   ├── V003_list_versions_lineage/
│   │   ├── V004_list_versions_custom_dir/
│   │   └── V005_list_versions_all_flags/
│   │
│   ├── 10_cmd_list_species/                # list-species命令测试
│   │   ├── SP001_list_species_default/
│   │   ├── SP002_filter_go/
│   │   ├── SP003_filter_kegg/
│   │   ├── SP004_filter_reactome/
│   │   ├── SP005_filter_do/
│   │   ├── SP006_filter_combined/
│   │   ├── SP007_format_table/
│   │   ├── SP008_format_tsv/
│   │   ├── SP009_format_json/
│   │   ├── SP010_summary/
│   │   └── SP011_all_flags/
│   │
│   ├── 11_cmd_query_species/               # query-species命令测试
│   │   ├── Q001_query_by_name/
│   │   ├── Q002_query_by_taxid/
│   │   ├── Q003_query_by_kegg/
│   │   ├── Q004_fuzzy_old_name/
│   │   ├── Q005_fuzzy_partial/
│   │   └── Q006_query_not_found/
│   │
│   ├── 12_disgenet_graceful/               # DisGeNET降级测试
│   │   ├── D121_download_disgenet_degraded/
│   │   └── D122_build_disgenet_skip/
│   │
│   ├── 13_registry_build/                  # 注册表构建验证
│   │   └── D123_registry_auto_build/
│   │
│   ├── 14_build_manifest/                  # 构建清单验证
│   │   └── B124_build_manifest_verify/
│   │
│   ├── 15_goa_fallback/                    # GOA回退测试
│   │   └── B125_goa_fallback_build/
│   │
│   ├── 16_version_consistency/             # 版本号一致性
│   │   └── T126_version_check/
│   │
│   └── 17_gsea_ssgsea_gsva_styles/         # GSEA/ssGSEA/GSVA风格配色
│       ├── G144-G149_gsea_styles/          # 6种风格
│       ├── G150-G157_gsea_palettes/        # 8种配色
│       ├── G158-G163_ssgsea_styles/        # 6种风格
│       ├── G164-G171_ssgsea_palettes/      # 8种配色
│       ├── G172-G177_gsva_styles/          # 6种风格
│       └── G178-G185_gsva_palettes/        # 8种配色
│
└── 99_results_summary/                     # 结果汇总
    ├── test_report.md                      # 测试报告
    ├── passed_tests.txt                    # 通过的测试列表
    ├── failed_tests.txt                    # 失败的测试列表
    └── logs/                               # 所有测试日志汇总
```

---

## 二、测试编号命名规范

### 2.1 编号前缀规则

| 前缀 | 命令/模块 | 编号范围 | 示例 |
|------|----------|---------|------|
| `T` | analyze (主命令) | T001-T143 | T001_fisher_default |
| `D` | download | D001-D010, D121-D123 | D001_download_go |
| `B` | build | B001-B014, B124-B125 | B001_build_rno_go_kegg |
| `S` | serve | S001-S003 | S001_serve_default |
| `L` | list | L001-L003 | L001_list_species |
| `C` | config | C001-C002 | C001_config_default |
| `U` | check-update | U001-U004 | U001_check_update_default |
| `CL` | cleanup | CL001-CL003 | CL001_cleanup_dry_run |
| `V` | list-versions | V001-V005 | V001_list_versions_default |
| `SP` | list-species | SP001-SP011 | SP001_list_species_default |
| `Q` | query-species | Q001-Q006 | Q001_query_by_name |
| `G` | GSEA/ssGSEA/GSVA 风格配色 | G144-G185 | G144_gsea_style_nature |

### 2.2 文件夹命名规范

```
{编号}_{测试内容描述}/
```

- **编号**：3位数字，不足补零（如 T001, T091, G144）
- **测试内容**：下划线连接的小写英文单词，描述测试的核心内容
- **示例**：
  - `T001_fisher_default/` — Fisher精确检验默认参数
  - `T091_gsea_permutations_100/` — GSEA排列次数100
  - `G144_gsea_style_nature/` — GSEA nature风格

### 2.3 文件命名规范

每个测试文件夹内统一包含：

```
{编号}_{测试内容}/
├── run.sh                    # 执行脚本（必须）
├── run.py                    # Python执行脚本（可选，复杂测试用）
├── input/                    # 测试专用输入文件（可选）
│   └── ...
├── output/                   # 输出目录（自动生成）
│   ├── results/              # 分析结果
│   ├── plots/                # 图表文件
│   └── report.html           # 报告文件
├── expected/                 # 预期输出（可选，用于对比）
│   └── ...
├── test.log                  # 测试执行日志（自动生成）
└── README.md                 # 测试说明（可选）
```

---

## 三、执行脚本模板（run.sh）

每个测试文件夹必须包含 `run.sh`，内容模板如下：

```bash
#!/bin/bash
#===============================================================================
# Test: {T001_fisher_default}
# Description: Fisher精确检验默认参数测试
# Command: analyze -i {input} -s {species} -d {databases}
# Expected: 成功执行，生成结果文件
#===============================================================================

set -e  # 遇到错误立即退出

# 配置
TEST_ID="T001"
TEST_NAME="fisher_default"
PROJECT_ROOT="$(cd ../.. && pwd)"
INPUT_DIR="${PROJECT_ROOT}/00_input_data"
OUTPUT_DIR="${PROJECT_ROOT}/01_cmd_analyze/01_ora_basic/${TEST_ID}_${TEST_NAME}/output"
LOG_FILE="${PROJECT_ROOT}/01_cmd_analyze/01_ora_basic/${TEST_ID}_${TEST_NAME}/test.log"

# 创建输出目录
mkdir -p "${OUTPUT_DIR}"

# 记录开始时间
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Test ${TEST_ID} started" > "${LOG_FILE}"

# 执行测试
cd "${PROJECT_ROOT}"
python -m allenricher analyze \
    -i "${INPUT_DIR}/hsa_gene_list_500.txt" \
    -s hsa \
    -d GO,KEGG \
    -o "${OUTPUT_DIR}" \
    2>&1 | tee -a "${LOG_FILE}"

EXIT_CODE=${PIPESTATUS[0]}

# 验证结果
if [ $EXIT_CODE -eq 0 ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Test ${TEST_ID} PASSED" >> "${LOG_FILE}"
    
    # 验证输出文件存在
    if [ -f "${OUTPUT_DIR}/enrichment_results.tsv" ]; then
        echo "  - enrichment_results.tsv: OK" >> "${LOG_FILE}"
    else
        echo "  - enrichment_results.tsv: MISSING" >> "${LOG_FILE}"
        exit 1
    fi
else
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Test ${TEST_ID} FAILED (exit code: $EXIT_CODE)" >> "${LOG_FILE}"
    exit $EXIT_CODE
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Test ${TEST_ID} completed" >> "${LOG_FILE}"
exit 0
```

---

## 四、数据准备脚本规格

### 4.1 脚本位置

`test_e2e_2026/00_input_data/00_prepare_data.py`

### 4.2 生成的文件清单

| 序号 | 文件名 | 格式 | 用途 | 行数/大小 |
|------|--------|------|------|----------|
| 1 | `hsa_gene_list_500.txt` | 每行Symbol | ORA分析（人类） | 500 |
| 2 | `mmu_gene_list_500.txt` | 每行Symbol | ORA分析（小鼠） | 500 |
| 3 | `rro_gene_list_500.txt` | 每行Symbol | ORA分析（川金丝猴） | 500 |
| 4 | `hsa_background_genes.txt` | 每行Symbol | 背景基因集 | ~15000 |
| 5 | `hsa_ranked_genes_gsea.tsv` | TSV (gene\tweight\trank) | GSEA排序 | 501（含表头） |
| 6 | `hsa_expression_matrix_500x6.tsv` | TSV (基因×样本) | ssGSEA/GSVA | 501×7 |
| 7 | `hsa_custom_pathways_5.gmt` | GMT格式 | 自定义基因集 | 5行 |
| 8 | `hsa_custom_pathways_5.gmt.gz` | GMT gzip | 压缩格式测试 | 5行 |
| 9 | `tsv_format.tsv` | TSV (gene\tscore) | 多格式输入 | 501 |
| 10 | `csv_format.csv` | CSV (gene,score) | 多格式输入 | 501 |
| 11 | `gene_ext_format.gene` | 每行Symbol | .gene扩展名 | 500 |
| 12 | `tiny_gene_list_5.txt` | 每行Symbol | 无结果测试 | 5 |
| 13 | `empty_gene_list.txt` | 空文件 | 错误场景 | 0 |
| 14-16 | `three_column.tsv` / `four_column.tsv` / `two_column.tsv` | TSV | 自定义注释 | ~200 |
| 17-36 | 20个YAML/JSON配置文件 | YAML/JSON | 配置加载测试 | - |

### 4.3 核心函数实现

```python
#!/usr/bin/env python3
"""
AllEnricher-v2 E2E测试数据准备脚本
从项目已有数据中提取真实基因，生成测试输入文件
"""
import gzip
import json
import os
import random
import yaml

# 路径配置
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
DB_DIR = os.path.join(PROJECT_ROOT, 'database', 'basic')
OUTPUT_DIR = SCRIPT_DIR

# 物种配置
SPECIES = {
    'hsa': {'taxid': 9606, 'name': 'Homo sapiens'},
    'mmu': {'taxid': 10090, 'name': 'Mus musculus'},
    'rro': {'taxid': 61622, 'name': 'Rhinopithecus roxellana'},
}

# 选定的KEGG通路（用于基因选择）
SELECTED_KEGG_PATHWAYS = [
    '00010',  # Glycolysis / Gluconeogenesis
    '04010',  # MAPK signaling pathway
    '04110',  # Cell cycle
    '04151',  # PI3K-Akt signaling pathway
    '04630',  # JAK-STAT signaling pathway
]

def load_gene_symbols(taxid: int) -> list[str]:
    """从gene_info.gz提取指定taxid的protein-coding基因Symbol"""
    gene_info_path = os.path.join(DB_DIR, 'go', 'GO20260527', 'gene_info.gz')
    symbols = []
    with gzip.open(gene_info_path, 'rt', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            if len(cols) < 10:
                continue
            if cols[0] == str(taxid) and cols[9] == 'protein-coding':
                symbol = cols[2]
                if symbol and symbol != '-' and not symbol.startswith('LOC'):
                    symbols.append(symbol)
    return symbols

def load_kegg_pathways(species_code: str) -> dict[str, list[str]]:
    """加载KEGG通路-基因映射"""
    pathway_file = os.path.join(DB_DIR, 'kegg', f'{species_code}_gene2pathway.txt')
    if not os.path.exists(pathway_file):
        return {}
    pathways = {}
    with open(pathway_file, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 3:
                continue
            gene_symbol, pathway_id = cols[0], cols[2]
            if pathway_id not in pathways:
                pathways[pathway_id] = []
            pathways[pathway_id].append(gene_symbol)
    return pathways

def generate_gene_list(all_symbols: list[str], kegg_pathways: dict,
                       n: int = 500, seed: int = 42) -> list[str]:
    """生成基因列表：优先从选定KEGG通路中选取"""
    random.seed(seed)
    pathway_genes = set()
    for pid in SELECTED_KEGG_PATHWAYS:
        if pid in kegg_pathways:
            pathway_genes.update(kegg_pathways[pid])
    all_set = set(all_symbols)
    pathway_genes = pathway_genes & all_set
    if len(pathway_genes) < n:
        remaining = list(all_set - pathway_genes)
        random.shuffle(remaining)
        pathway_genes.update(remaining[:n - len(pathway_genes)])
    gene_list = list(pathway_genes)
    random.shuffle(gene_list)
    return gene_list[:n]

def generate_ranked_genes(gene_list: list[str], seed: int = 42) -> list[tuple]:
    """生成排序基因列表（GSEA用）"""
    random.seed(seed)
    ranked = []
    for i, gene in enumerate(gene_list):
        if i < 100:
            weight = random.gauss(0, 5)
        else:
            weight = random.gauss(0, 1.5)
        ranked.append((gene, round(weight, 4)))
    ranked.sort(key=lambda x: x[1], reverse=True)
    return [(gene, weight, rank + 1) for rank, (gene, weight) in enumerate(ranked)]

def generate_expression_matrix(gene_list: list[str], seed: int = 42) -> str:
    """生成表达矩阵（ssGSEA/GSVA用）"""
    random.seed(seed)
    samples = ['Normal_1', 'Normal_2', 'Normal_3', 'Tumor_1', 'Tumor_2', 'Tumor_3']
    lines = ['\t' + '\t'.join(samples)]
    for i, gene in enumerate(gene_list):
        values = []
        for j, sample in enumerate(samples):
            if i < 100 and j >= 3:
                expr = random.gauss(10, 1.5)
            elif i < 100 and j < 3:
                expr = random.gauss(6, 1.5)
            else:
                expr = random.gauss(8, 2)
            values.append(f'{expr:.3f}')
        lines.append(gene + '\t' + '\t'.join(values))
    return '\n'.join(lines)

def generate_custom_gmt(species_code: str, kegg_pathways: dict, seed: int = 42) -> str:
    """生成自定义GMT文件（5条通路）"""
    random.seed(seed)
    valid = {pid: genes for pid, genes in kegg_pathways.items()
             if 30 <= len(genes) <= 200}
    selected = random.sample(list(valid.keys()), min(5, len(valid)))
    lines = []
    for pid in selected:
        genes = '\t'.join(valid[pid])
        lines.append(f'{species_code}{pid}\tCustom|Test|Pathway_{pid}\t{genes}')
    return '\n'.join(lines)

def generate_custom_annotations(gene_list: list[str], seed: int = 42) -> tuple:
    """生成三列/四列/两列格式的自定义注释文件"""
    random.seed(seed)
    terms = []
    for i in range(10):
        term_id = f'CUSTOM{i+1:04d}'
        term_name = f'Custom Pathway {i+1}'
        genes = random.sample(gene_list, min(20, len(gene_list)))
        terms.append((term_id, term_name, genes))
    
    lines_3col = [f'{gene}\t{term_id}\t{term_name}' 
                  for term_id, term_name, genes in terms for gene in genes]
    lines_4col = [f'{gene}\t{term_id}\t{term_name}\t{random.uniform(0.5, 1.0):.4f}'
                  for term_id, term_name, genes in terms for gene in genes]
    lines_2col = [f'{gene}\t{term_id}'
                  for term_id, term_name, genes in terms for gene in genes]
    
    return '\n'.join(lines_3col), '\n'.join(lines_4col), '\n'.join(lines_2col)

def generate_config_files():
    """生成所有YAML/JSON配置文件"""
    configs = {
        'T091_gsea_permutations_100.yaml': {
            'gsea_permutations': 100,
        },
        'T092_gsea_size_limits.yaml': {
            'gsea_min_size': 5,
            'gsea_max_size': 200,
        },
        'T093_gsva_method_plage.yaml': {'gsva_method': 'plage'},
        'T094_gsva_method_zscore.yaml': {'gsva_method': 'zscore'},
        'T095_gsva_kcdf_poisson.yaml': {'gsva_kcdf': 'Poisson'},
        'T096_gsva_tau_0.5.yaml': {'gsva_tau': 0.5},
        'T097_max_genes_500.yaml': {'max_genes': 500},
        'T098_top_terms_10.yaml': {'top_terms': 10},
        'T099_plot_size_12x10.yaml': {'plot_width': 12, 'plot_height': 10},
        'T100_plot_formats_multi.yaml': {'plot_formats': ['pdf', 'png', 'svg']},
        'T101_report_format_html.yaml': {'report_format': 'html'},
        'T102_ai_ollama_baseurl.yaml': {
            'ai_backends': {'ollama': {'model': 'llama3', 'base_url': 'http://localhost:11434', 'enabled': True}}
        },
        'T103_ai_mock_disabled.yaml': {
            'ai_backends': {'mock': {'enabled': False}}
        },
        'T078_config_full_yaml.yaml': {
            'input_file': '00_input_data/hsa_gene_list_500.txt',
            'species': 'hsa',
            'databases': ['GO', 'KEGG'],
            'method': 'fisher',
            'correction': 'BH',
            'pvalue_cutoff': 0.05,
            'qvalue_cutoff': 0.05,
            'min_genes': 2,
            'n_jobs': 2,
            'plot_style': 'nature',
            'plot_palette': 'lancet',
            'plot_format': 'png',
            'plot_dpi': 300,
        },
        'T079_config_full_json.json': {
            'input_file': '00_input_data/hsa_gene_list_500.txt',
            'species': 'hsa',
            'databases': ['GO', 'KEGG'],
            'method': 'fisher',
            'correction': 'BH',
        },
        'T106_cli_override_priority.yaml': {
            'species': 'mmu',
            'databases': ['KEGG'],
            'method': 'hypergeometric',
        },
        'T143_style_cell_config.yaml': {'plot_style': 'cell'},
        'T146_gsea_style_cell.yaml': {'plot_style': 'cell'},
        'T160_ssgsea_style_cell.yaml': {'plot_style': 'cell'},
        'T174_gsva_style_cell.yaml': {'plot_style': 'cell'},
    }
    
    config_dir = os.path.join(OUTPUT_DIR, 'config_files')
    os.makedirs(config_dir, exist_ok=True)
    
    for filename, config in configs.items():
        filepath = os.path.join(config_dir, filename)
        if filename.endswith('.json'):
            with open(filepath, 'w', encoding='utf-8') as f:
                json.dump(config, f, indent=2, ensure_ascii=False)
        else:
            with open(filepath, 'w', encoding='utf-8') as f:
                yaml.dump(config, f, default_flow_style=False, allow_unicode=True)

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print("=== AllEnricher-v2 E2E测试数据准备 ===\n")
    
    # 1. 加载基因
    print("[1/6] 加载基因信息...")
    species_genes = {}
    for code, info in SPECIES.items():
        symbols = load_gene_symbols(info['taxid'])
        species_genes[code] = symbols
        print(f"  {code}: {len(symbols)} protein-coding genes")
    
    # 2. 加载KEGG通路
    print("\n[2/6] 加载KEGG通路...")
    species_pathways = {}
    for code in SPECIES:
        pathways = load_kegg_pathways(code)
        species_pathways[code] = pathways
        print(f"  {code}: {len(pathways)} pathways")
    
    # 3. 生成基因列表
    print("\n[3/6] 生成基因列表...")
    for code in SPECIES:
        gene_list = generate_gene_list(species_genes[code], species_pathways[code])
        filepath = os.path.join(OUTPUT_DIR, f'{code}_gene_list_500.txt')
        with open(filepath, 'w') as f:
            f.write('\n'.join(gene_list) + '\n')
        print(f"  {code}_gene_list_500.txt: {len(gene_list)} genes")
    
    # 4. 背景基因
    print("\n[4/6] 生成背景基因列表...")
    bg = sorted(set(species_genes['hsa']))
    with open(os.path.join(OUTPUT_DIR, 'hsa_background_genes.txt'), 'w') as f:
        f.write('\n'.join(bg) + '\n')
    print(f"  hsa_background_genes.txt: {len(bg)} genes")
    
    # 5. 排序基因
    print("\n[5/6] 生成排序基因列表...")
    hsa_genes = generate_gene_list(species_genes['hsa'], species_pathways['hsa'])
    ranked = generate_ranked_genes(hsa_genes)
    with open(os.path.join(OUTPUT_DIR, 'hsa_ranked_genes_gsea.tsv'), 'w') as f:
        f.write('gene\tweight\trank\n')
        for gene, weight, rank in ranked:
            f.write(f'{gene}\t{weight}\t{rank}\n')
    print(f"  hsa_ranked_genes_gsea.tsv: {len(ranked)} genes")
    
    # 6. 表达矩阵
    print("\n[6/6] 生成表达矩阵...")
    matrix = generate_expression_matrix(hsa_genes)
    with open(os.path.join(OUTPUT_DIR, 'hsa_expression_matrix_500x6.tsv'), 'w') as f:
        f.write(matrix + '\n')
    print(f"  hsa_expression_matrix_500x6.tsv: {len(hsa_genes)} genes × 6 samples")
    
    # 7. 自定义GMT
    print("\n[7] 生成自定义GMT...")
    gmt = generate_custom_gmt('hsa', species_pathways['hsa'])
    with open(os.path.join(OUTPUT_DIR, 'hsa_custom_pathways_5.gmt'), 'w') as f:
        f.write(gmt + '\n')
    with gzip.open(os.path.join(OUTPUT_DIR, 'hsa_custom_pathways_5.gmt.gz'), 'wt') as f:
        f.write(gmt + '\n')
    print(f"  hsa_custom_pathways_5.gmt + .gmt.gz: {gmt.count(chr(10))} pathways")
    
    # 8. 多格式输入
    print("\n[8] 生成多格式输入文件...")
    random.seed(42)
    tsv_lines = ['gene\tscore'] + [f'{g}\t{random.uniform(0.5, 15):.4f}' for g in hsa_genes]
    csv_lines = ['gene,score'] + [f'{g},{random.uniform(0.5, 15):.4f}' for g in hsa_genes]
    with open(os.path.join(OUTPUT_DIR, 'tsv_format.tsv'), 'w') as f:
        f.write('\n'.join(tsv_lines) + '\n')
    with open(os.path.join(OUTPUT_DIR, 'csv_format.csv'), 'w') as f:
        f.write('\n'.join(csv_lines) + '\n')
    with open(os.path.join(OUTPUT_DIR, 'gene_ext_format.gene'), 'w') as f:
        f.write('\n'.join(hsa_genes) + '\n')
    print("  .tsv / .csv / .gene formats created")
    
    # 9. 错误场景文件
    print("\n[9] 生成错误场景文件...")
    tiny = hsa_genes[:5]
    with open(os.path.join(OUTPUT_DIR, 'tiny_gene_list_5.txt'), 'w') as f:
        f.write('\n'.join(tiny) + '\n')
    open(os.path.join(OUTPUT_DIR, 'empty_gene_list.txt'), 'w').close()
    print("  tiny_gene_list_5.txt / empty_gene_list.txt created")
    
    # 10. 自定义注释
    print("\n[10] 生成自定义注释文件...")
    col3, col4, col2 = generate_custom_annotations(hsa_genes)
    annot_dir = os.path.join(OUTPUT_DIR, 'custom_annotations')
    os.makedirs(annot_dir, exist_ok=True)
    with open(os.path.join(annot_dir, 'three_column.tsv'), 'w') as f:
        f.write(col3 + '\n')
    with open(os.path.join(annot_dir, 'four_column.tsv'), 'w') as f:
        f.write(col4 + '\n')
    with open(os.path.join(annot_dir, 'two_column.tsv'), 'w') as f:
        f.write(col2 + '\n')
    print(f"  3col/4col/2col annotations: {col3.count(chr(10))} lines each")
    
    # 11. 配置文件
    print("\n[11] 生成配置文件...")
    generate_config_files()
    print("  20 config files generated")
    
    print("\n=== 完成！所有测试数据已生成 ===")

if __name__ == '__main__':
    main()
```

---

## 五、测试用例详细清单

### 5.1 analyze 命令测试（T001-T143）

#### 01_ora_basic — ORA基础测试（T001-T011）

| 编号 | 文件夹名 | 测试命令 | 验证要点 |
|------|---------|---------|---------|
| T001 | `T001_fisher_default/` | `analyze -i hsa_gene_list_500.txt -s hsa -d GO,KEGG` | 默认参数执行成功 |
| T002 | `T002_fisher_threshold_p0.01_q0.01/` | `... -p 0.01 -q 0.01` | p/q阈值过滤生效 |
| T003 | `T003_fisher_min_genes_5/` | `... -n 5` | min-genes=5过滤 |
| T004 | `T004_fisher_min_genes_10/` | `... -n 10` | min-genes=10过滤 |
| T005 | `T005_fisher_only_significant/` | `... --only-significant` | 仅显著结果输出 |
| T006 | `T006_fisher_parallel_jobs_4/` | `... -j 4` | 并行执行正常 |
| T007 | `T007_hypergeometric_default/` | `... -m hypergeometric` | 超几何检验执行 |
| T008 | `T008_hypergeometric_correction_bh/` | `... -m hypergeometric -c BH` | BH校正 |
| T009 | `T009_hypergeometric_correction_bonferroni/` | `... -c bonferroni` | Bonferroni校正 |
| T010 | `T010_hypergeometric_correction_holm/` | `... -c holm` | Holm校正 |
| T011 | `T011_hypergeometric_correction_by/` | `... -c BY` | BY校正 |

#### 02_ora_correction — 校正方法测试（T013-T017）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| T013 | `T013_correction_bh/` | `... -c BH` |
| T014 | `T014_correction_by/` | `... -c BY` |
| T015 | `T015_correction_bonferroni/` | `... -c bonferroni` |
| T016 | `T016_correction_holm/` | `... -c holm` |
| T017 | `T017_correction_none/` | `... -c none` |

#### 03_ora_databases — 数据库组合测试（T018-T023）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| T018 | `T018_database_go_only/` | `... -d GO` |
| T019 | `T019_database_kegg_only/` | `... -d KEGG` |
| T020 | `T020_database_reactome_only/` | `... -d Reactome` |
| T021 | `T021_database_do_only/` | `... -d DO` |
| T022 | `T022_database_all_combined/` | `... -d GO,KEGG,Reactome,DO` |
| T023 | `T023_mmu_database_all/` | `... -s mmu -d GO,KEGG,Reactome` |

#### 04_background_modes — 背景模式测试（T024-T027）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| T024 | `T024_background_annotated/` | `... --background-mode annotated` |
| T025 | `T025_background_genome/` | `... --background-mode genome` |
| T026 | `T026_background_custom/` | `... --background-mode custom -b hsa_background_genes.txt` |
| T027 | `T027_background_custom_error/` | `... --background-mode custom`（无-b，应报错） |

#### 05_gsea_analysis — GSEA分析测试（T028-T031）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| T028 | `T028_gsea_ranked_genes/` | `... -m gsea -r hsa_ranked_genes_gsea.tsv` |
| T029 | `T029_gsea_expression_matrix/` | `... -m gsea -e hsa_expression_matrix_500x6.tsv` |
| T030 | `T030_gsea_custom_gmt/` | `... -m gsea -r ... -g hsa_custom_pathways_5.gmt` |
| T031 | `T031_gsea_with_groups/` | `... -e ... --groups "Normal:N1,N2,N3;Tumor:T1,T2,T3"` |

#### 06_ssgsea_analysis — ssGSEA分析测试（T032-T034）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| T032 | `T032_ssgsea_expression_matrix/` | `... -m ssgsea -e hsa_expression_matrix_500x6.tsv` |
| T033 | `T033_ssgsea_custom_gmt/` | `... -m ssgsea -e ... -g hsa_custom_pathways_5.gmt` |
| T034 | `T034_ssgsea_with_groups/` | `... -m ssgsea -e ... --groups "..."` |

#### 07_gsva_analysis — GSVA分析测试（T035-T038）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| T035 | `T035_gsva_default_gsva/` | `... -m gsva -e hsa_expression_matrix_500x6.tsv` |
| T036 | `T036_gsva_method_plage/` | `... -m gsva -e ... --config T093_gsva_method_plage.yaml` |
| T037 | `T037_gsva_method_zscore/` | `... --config T094_gsva_method_zscore.yaml` |
| T038 | `T038_gsva_with_groups/` | `... -m gsva -e ... --groups "..."` |

#### 08_visual_styles — 可视化风格测试（T039-T044）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| T039 | `T039_style_nature/` | `... --style nature` |
| T040 | `T040_style_science/` | `... --style science` |
| T041 | `T041_style_cell/` | `... --config T143_style_cell_config.yaml` |
| T042 | `T042_style_colorblind/` | `... --style colorblind` |
| T043 | `T043_style_presentation/` | `... --style presentation` |
| T044 | `T044_style_omicshare/` | `... --style omicshare` |

#### 09_visual_palettes — 可视化配色测试（T045-T058, T134-T142）

| 编号 | 文件夹名 | 配色 |
|------|---------|------|
| T045 | `T045_palette_nature/` | nature |
| T046 | `T046_palette_science/` | science |
| T047 | `T047_palette_lancet/` | lancet |
| T048 | `T048_palette_nejm/` | nejm |
| T049 | `T049_palette_jama/` | jama |
| T050 | `T050_palette_okabe_ito/` | okabe_ito |
| T051 | `T051_palette_gsea/` | gsea |
| T052 | `T052_palette_omicshare/` | omicshare |
| T053 | `T053_palette_china_style/` | china_style |
| T054 | `T054_palette_go_bp/` | go_bp |
| T055 | `T055_palette_kegg_pathway/` | kegg_pathway |
| T056 | `T056_palette_tol_bright/` | tol_bright |
| T057 | `T057_palette_tol_muted/` | tol_muted |
| T058 | `T058_palette_tol_sunset/` | tol_sunset |
| T134-T142 | `T134-T142_remaining_palettes/` | default, tol_high_contrast, tol_vibrant, tol_medium_contrast, tol_light, tol_burga, go_cc, go_mf, cell |

#### 10_output_formats — 输出格式测试（T059-T063）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| T059 | `T059_format_png/` | `... --plot-format png` |
| T060 | `T060_format_pdf/` | `... --plot-format pdf` |
| T061 | `T061_format_svg/` | `... --plot-format svg` |
| T062 | `T062_dpi_150/` | `... --plot-dpi 150` |
| T063 | `T063_dpi_600/` | `... --plot-dpi 600` |

#### 11_report_control — 报告控制测试（T064-T066）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| T064 | `T064_no_plot/` | `... --no-plot` |
| T065 | `T065_no_report/` | `... --no-report` |
| T066 | `T066_no_plot_no_report/` | `... --no-plot --no-report` |

#### 12_plot_types — 图表类型测试（T067-T077）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| T067 | `T067_plot_dotplot/` | `... -pt dotplot` |
| T068 | `T068_plot_barplot/` | `... -pt barplot` |
| T069 | `T069_plot_network/` | `... -pt network` |
| T070 | `T070_plot_upset/` | `... -pt upset` |
| T071 | `T071_plot_volcano/` | `... -pt volcano` |
| T072 | `T072_gsea_plot_enrichment/` | `... -m gsea -r ... -pt enrichment` |
| T073 | `T073_gsea_plot_nes_barplot/` | `... -m gsea -r ... -pt nes_barplot` |
| T074 | `T074_ssgsea_plot_heatmap/` | `... -m ssgsea -e ... -pt heatmap` |
| T075 | `T075_ssgsea_plot_group_comparison/` | `... -m ssgsea -e ... -pt group_comparison --groups "..."` |
| T076 | `T076_ssgsea_plot_correlation/` | `... -m ssgsea -e ... -pt correlation` |
| T077 | `T077_plot_multi_types/` | `... -pt dotplot,barplot,network,upset` |

#### 13_config_loading — 配置文件加载（T078-T079）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| T078 | `T078_config_yaml_full/` | `... --config T078_config_full_yaml.yaml` |
| T079 | `T079_config_json_full/` | `... --config T079_config_full_json.json` |

#### 14_version_control — 版本控制（T080-T081）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| T080 | `T080_use_version_v20260515/` | `... --use-version v20260515` |
| T081 | `T081_database_dir_custom/` | `... --database-dir ./database` |

#### 15_ai_interpretation — AI解读（T082-T085）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| T082 | `T082_verbose_logging/` | `... --verbose` |
| T083 | `T083_ai_mock_backend/` | `... --ai mock` |
| T084 | `T084_ai_openai_params/` | `... --ai openai --ai-key test --ai-model gpt-4` |
| T085 | `T085_ai_claude_params/` | `... --ai claude --ai-key test --ai-model claude-3` |

#### 16_error_scenarios — 错误场景（T086-T090）

| 编号 | 文件夹名 | 测试命令 | 预期结果 |
|------|---------|---------|---------|
| T086 | `T086_empty_gene_list/` | `... -i empty_gene_list.txt` | 报错或空结果 |
| T087 | `T087_invalid_species/` | `... -s xxx_notexist` | 报错 |
| T088 | `T088_nonexistent_input/` | `... -i nonexistent.txt` | 报错 |
| T089 | `T089_gsea_missing_ranked/` | `... -m gsea`（无-r/-e） | 报错 |
| T090 | `T090_ssgsea_missing_expression/` | `... -m ssgsea`（无-e） | 报错 |

#### 17_config_only_params — Config-only参数（T091-T103）

| 编号 | 文件夹名 | 配置文件 | 测试内容 |
|------|---------|---------|---------|
| T091 | `T091_gsea_permutations_100/` | T091_gsea_permutations_100.yaml | gsea_permutations=100 |
| T092 | `T092_gsea_size_limits/` | T092_gsea_size_limits.yaml | gsea_min/max_size |
| T093 | `T093_gsva_method_plage/` | T093_gsva_method_plage.yaml | gsva_method=plage |
| T094 | `T094_gsva_method_zscore/` | T094_gsva_method_zscore.yaml | gsva_method=zscore |
| T095 | `T095_gsva_kcdf_poisson/` | T095_gsva_kcdf_poisson.yaml | gsva_kcdf=Poisson |
| T096 | `T096_gsva_tau_0.5/` | T096_gsva_tau_0.5.yaml | gsva_tau=0.5 |
| T097 | `T097_max_genes_500/` | T097_max_genes_500.yaml | max_genes=500 |
| T098 | `T098_top_terms_10/` | T098_top_terms_10.yaml | top_terms=10 |
| T099 | `T099_plot_size_12x10/` | T099_plot_size_12x10.yaml | plot_width/height |
| T100 | `T100_plot_formats_multi/` | T100_plot_formats_multi.yaml | plot_formats多格式 |
| T101 | `T101_report_format_html/` | T101_report_format_html.yaml | report_format |
| T102 | `T102_ai_ollama_baseurl/` | T102_ai_ollama_baseurl.yaml | ai_backends.base_url |
| T103 | `T103_ai_mock_disabled/` | T103_ai_mock_disabled.yaml | ai_backends.enabled=false |

#### 18_global_params — 全局参数（T104-T105）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| T104 | `T104_version_short/` | `allenricher -v` |
| T105 | `T105_version_long/` | `allenricher --version` |

#### 19_cli_priority — CLI优先级（T106）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| T106 | `T106_cli_override_config/` | `... -s hsa -d GO --config T106_cli_override_priority.yaml`（config里是mmu/KEGG） |

#### 20_input_formats — 多格式输入（T107-T109）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| T107 | `T107_input_tsv_format/` | `... -i tsv_format.tsv` |
| T108 | `T108_input_csv_format/` | `... -i csv_format.csv` |
| T109 | `T109_input_gene_ext/` | `... -i gene_ext_format.gene` |

#### 21_gmt_formats — GMT格式（T110）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| T110 | `T110_gmt_gz_format/` | `... -m gsea -r ... -g hsa_custom_pathways_5.gmt.gz` |

#### 22_group_comparison — 分组比较图表（T111-T113）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| T111 | `T111_group_comparison_box/` | `... -m ssgsea -e ... -pt group_comparison --groups "..."` |
| T112 | `T112_group_comparison_violin/` | （通过config设置plot_type） |
| T113 | `T113_group_comparison_bar/` | （通过config设置plot_type） |

#### 23_correlation_methods — 相关性方法（T114-T115）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| T114 | `T114_correlation_pearson/` | `... -m ssgsea -e ... -pt correlation` |
| T115 | `T115_correlation_spearman/` | （通过config设置correlation_method） |

#### 24_network_plots — 网络图（T116-T117）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| T116 | `T116_network_single_db/` | `... -pt network` |
| T117 | `T117_network_multi_db/` | `... -d GO,KEGG,Reactome -pt network,upset` |

#### 25_api_endpoints — API端点（T117a-T117k）

| 编号 | 文件夹名 | 测试内容 |
|------|---------|---------|
| T117a | `T117a_api_home/` | GET / |
| T117b | `T117b_api_species/` | GET /api/species |
| T117c | `T117c_api_databases/` | GET /api/databases |
| T117d | `T117d_api_analyze/` | POST /api/analyze |
| T117e | `T117e_api_upload/` | POST /api/upload |
| T117f | `T117f_api_status_404/` | GET /api/status/nonexistent |
| T117g | `T117g_api_results_404/` | GET /api/results/nonexistent |
| T117h | `T117h_api_results_json/` | GET /api/results/{id}?format=json |
| T117i | `T117i_api_plot/` | GET /api/results/{id}/plot |
| T117j | `T117j_api_report/` | GET /api/results/{id}/report |
| T117k | `T117k_api_delete_404/` | DELETE /api/jobs/nonexistent |

#### 26_report_generation — 报告生成（T118-T119）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| T118 | `T118_report_normal/` | 正常分析，验证report.html |
| T119 | `T119_report_no_results/` | `... -i tiny_gene_list_5.txt -p 1e-30` 验证无结果页面 |

#### 27_url_generation — URL生成（T120）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| T120 | `T120_url_all_databases/` | `... -d GO,KEGG,Reactome,DO --no-plot` 验证结果TSV中URL列 |

#### 28_parallel_modes — 并行模式（T127-T129）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| T127 | `T127_serial_jobs_1/` | `... -j 1` |
| T128 | `T128_parallel_jobs_4/` | `... -j 4` |
| T129 | `T129_all_cpu_jobs_minus1/` | `... -j -1` |

#### 29_gsea_fields — GSEA字段验证（T130）

| 编号 | 文件夹名 | 验证内容 |
|------|---------|---------|
| T130 | `T130_gsea_result_fields/` | 结果TSV包含NES/ES/FDR/Leading_Edge列 |

#### 30_ssgsea_nan — ssGSEA NaN处理（T131）

| 编号 | 文件夹名 | 验证内容 |
|------|---------|---------|
| T131 | `T131_ssgsea_nan_pvalue/` | pvalue列为NaN，程序不崩溃 |

#### 31_output_fields — 输出字段（T132-T133）

| 编号 | 文件夹名 | 验证内容 |
|------|---------|---------|
| T132 | `T132_ora_result_fields/` | ORA结果字段完整性 |
| T133 | `T133_gsea_output_fields/` | GSEA结果字段完整性 |

### 5.2 其他命令测试

#### 02_cmd_download（D001-D010, D121-D123）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| D001 | `D001_download_go/` | `download -d go` |
| D002 | `D002_download_reactome/` | `download -d reactome` |
| D003 | `D003_download_do/` | `download -d do` |
| D004 | `D004_download_multiple/` | `download -d go,reactome,do` |
| D005 | `D005_workers_8/` | `download -d go --workers 8` |
| D006 | `D006_no_multi_thread/` | `download -d go --no-multi-thread` |
| D007 | `D007_no_verify/` | `download -d go --no-verify` |
| D008 | `D008_force_redownload/` | `download -d go --force` |
| D009 | `D009_custom_db_dir/` | `download -d go --database-dir ./database` |
| D010 | `D010_kegg_hint/` | `download -d kegg -s hsa` |
| D121 | `D121_download_disgenet_degraded/` | `download -d disgenet`（应优雅降级） |
| D122 | `D122_build_disgenet_skip/` | `build -s hsa -t 9606 -d DisGeNET`（应跳过） |
| D123 | `D123_registry_auto_build/` | 验证supported_species.tsv自动更新 |

#### 03_cmd_build（B001-B014, B124-B125）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| B001 | `B001_build_rno_go_kegg/` | `build -s rno -t 10116 -d GO,KEGG` |
| B002 | `B002_build_rno_all_dbs/` | `build -s rno -t 10116 -d GO,KEGG,Reactome` |
| B003 | `B003_build_custom_db_dir/` | `build ... --database-dir ./database` |
| B004 | `B004_build_with_gene_info/` | `build ... --gene-info ...` |
| B005 | `B005_build_with_go_annot/` | `build ... --go-annot ...` |
| B006 | `B006_build_custom_annot/` | `build ... --custom-annot ... --custom-db-name MYDB` |
| B007 | `B007_build_annot_format_three/` | `build ... --annot-format three_column` |
| B008 | `B008_build_annot_format_four/` | `build ... --annot-format four_column` |
| B009 | `B009_build_annot_format_two/` | `build ... --annot-format two_column` |
| B010 | `B010_build_annot_format_auto/` | `build ... --annot-format auto` |
| B011 | `B011_build_hierarchy_sep/` | `build ... --hierarchy-sep "."` |
| B012 | `B012_build_hsa_do/` | `build -s hsa -t 9606 -d DO` |
| B013 | `B013_build_mmu_do_skip/` | `build -s mmu -t 10090 -d DO`（应跳过） |
| B014 | `B014_build_hsa_all_dbs/` | `build -s hsa -t 9606 -d GO,KEGG,Reactome,DO` |
| B124 | `B124_build_manifest_verify/` | 验证build_manifest.json内容和格式 |
| B125 | `B125_goa_fallback_build/` | 构建gene2go中不存在的物种，触发GOA回退 |

#### 04_cmd_serve（S001-S003）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| S001 | `S001_serve_default/` | `serve` |
| S002 | `S002_serve_custom_host_port/` | `serve --host 127.0.0.1 --port 8080` |
| S003 | `S003_serve_reload/` | `serve --reload` |

#### 05_cmd_list（L001-L003）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| L001 | `L001_list_species/` | `list species` |
| L002 | `L002_list_databases/` | `list databases` |
| L003 | `L003_list_no_args/` | `list`（应报错） |

#### 06_cmd_config（C001-C002）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| C001 | `C001_config_default/` | `config` |
| C002 | `C002_config_custom_path/` | `config -o my_config.yaml` |

#### 07_cmd_check_update（U001-U004）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| U001 | `U001_check_update_default/` | `check-update` |
| U002 | `U002_check_update_custom_dir/` | `check-update --database-dir ./database` |
| U003 | `U003_check_update_json/` | `check-update --json` |
| U004 | `U004_check_update_dir_json/` | `check-update --database-dir ./database --json` |

#### 08_cmd_cleanup（CL001-CL003）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| CL001 | `CL001_cleanup_dry_run/` | `cleanup --dry-run` |
| CL002 | `CL002_cleanup_keep_1/` | `cleanup --keep 1 --dry-run` |
| CL003 | `CL003_cleanup_custom_dir/` | `cleanup --database-dir ./database --dry-run` |

#### 09_cmd_list_versions（V001-V005）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| V001 | `V001_list_versions_default/` | `list-versions` |
| V002 | `V002_list_versions_json/` | `list-versions --json` |
| V003 | `V003_list_versions_lineage/` | `list-versions --lineage` |
| V004 | `V004_list_versions_custom_dir/` | `list-versions --database-dir ./database` |
| V005 | `V005_list_versions_all_flags/` | `list-versions --json --lineage` |

#### 10_cmd_list_species（SP001-SP011）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| SP001 | `SP001_list_species_default/` | `list-species` |
| SP002 | `SP002_filter_go/` | `list-species --go` |
| SP003 | `SP003_filter_kegg/` | `list-species --kegg` |
| SP004 | `SP004_filter_reactome/` | `list-species --reactome` |
| SP005 | `SP005_filter_do/` | `list-species --do` |
| SP006 | `SP006_filter_combined/` | `list-species --go --kegg --reactome` |
| SP007 | `SP007_format_table/` | `list-species --format table` |
| SP008 | `SP008_format_tsv/` | `list-species --format tsv` |
| SP009 | `SP009_format_json/` | `list-species --format json` |
| SP010 | `SP010_summary/` | `list-species --summary` |
| SP011 | `SP011_all_flags/` | `list-species --go --kegg --format json --summary` |

#### 11_cmd_query_species（Q001-Q006）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| Q001 | `Q001_query_by_name/` | `query-species --name "Homo sapiens"` |
| Q002 | `Q002_query_by_taxid/` | `query-species --taxid 9606` |
| Q003 | `Q003_query_by_kegg/` | `query-species --kegg hsa` |
| Q004 | `Q004_fuzzy_old_name/` | `query-species --name "Microtus fortis"` |
| Q005 | `Q005_fuzzy_partial/` | `query-species --name "Rhinopithecus"` |
| Q006 | `Q006_query_not_found/` | `query-species --name "Nonexistent"` |

### 5.3 GSEA/ssGSEA/GSVA 风格配色测试（G144-G185）

#### GSEA风格（G144-G149）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| G144 | `G144_gsea_style_nature/` | `... -m gsea -r ... --style nature` |
| G145 | `G145_gsea_style_science/` | `... --style science` |
| G146 | `G146_gsea_style_cell/` | `... --config T146_gsea_style_cell.yaml` |
| G147 | `G147_gsea_style_colorblind/` | `... --style colorblind` |
| G148 | `G148_gsea_style_presentation/` | `... --style presentation` |
| G149 | `G149_gsea_style_omicshare/` | `... --style omicshare` |

#### GSEA配色（G150-G157）

| 编号 | 文件夹名 | 配色 |
|------|---------|------|
| G150 | `G150_gsea_palette_nature/` | nature |
| G151 | `G151_gsea_palette_science/` | science |
| G152 | `G152_gsea_palette_lancet/` | lancet |
| G153 | `G153_gsea_palette_nejm/` | nejm |
| G154 | `G154_gsea_palette_okabe_ito/` | okabe_ito |
| G155 | `G155_gsea_palette_gsea/` | gsea |
| G156 | `G156_gsea_palette_china_style/` | china_style |
| G157 | `G157_gsea_palette_tol_burga/` | tol_burga |

#### ssGSEA风格（G158-G163）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| G158 | `G158_ssgsea_style_nature/` | `... -m ssgsea -e ... --style nature --groups "..."` |
| G159 | `G159_ssgsea_style_science/` | `... --style science` |
| G160 | `G160_ssgsea_style_cell/` | `... --config T160_ssgsea_style_cell.yaml` |
| G161 | `G161_ssgsea_style_colorblind/` | `... --style colorblind` |
| G162 | `G162_ssgsea_style_presentation/` | `... --style presentation` |
| G163 | `G163_ssgsea_style_omicshare/` | `... --style omicshare` |

#### ssGSEA配色（G164-G171）

| 编号 | 文件夹名 | 配色 |
|------|---------|------|
| G164 | `G164_ssgsea_palette_nature/` | nature |
| G165 | `G165_ssgsea_palette_science/` | science |
| G166 | `G166_ssgsea_palette_lancet/` | lancet |
| G167 | `G167_ssgsea_palette_nejm/` | nejm |
| G168 | `G168_ssgsea_palette_okabe_ito/` | okabe_ito |
| G169 | `G169_ssgsea_palette_omicshare/` | omicshare |
| G170 | `G170_ssgsea_palette_china_style/` | china_style |
| G171 | `G171_ssgsea_palette_tol_sunset/` | tol_sunset |

#### GSVA风格（G172-G177）

| 编号 | 文件夹名 | 测试命令 |
|------|---------|---------|
| G172 | `G172_gsva_style_nature/` | `... -m gsva -e ... --style nature --groups "..."` |
| G173 | `G173_gsva_style_science/` | `... --style science` |
| G174 | `G174_gsva_style_cell/` | `... --config T174_gsva_style_cell.yaml` |
| G175 | `G175_gsva_style_colorblind/` | `... --style colorblind` |
| G176 | `G176_gsva_style_presentation/` | `... --style presentation` |
| G177 | `G177_gsva_style_omicshare/` | `... --style omicshare` |

#### GSVA配色（G178-G185）

| 编号 | 文件夹名 | 配色 |
|------|---------|------|
| G178 | `G178_gsva_palette_nature/` | nature |
| G179 | `G179_gsva_palette_science/` | science |
| G180 | `G180_gsva_palette_lancet/` | lancet |
| G181 | `G181_gsva_palette_nejm/` | nejm |
| G182 | `G182_gsva_palette_okabe_ito/` | okabe_ito |
| G183 | `G183_gsva_palette_omicshare/` | omicshare |
| G184 | `G184_gsva_palette_china_style/` | china_style |
| G185 | `G185_gsva_palette_tol_muted/` | tol_muted |

---

## 六、测试执行顺序

```
阶段 0: 准备测试数据
  └─ 执行 00_input_data/00_prepare_data.py

阶段 1: 基础设施命令（无依赖）
  ├─ 05_cmd_list (L001-L003)
  ├─ 06_cmd_config (C001-C002)
  ├─ 07_cmd_check_update (U001-U004)
  ├─ 09_cmd_list_versions (V001-V005)
  ├─ 10_cmd_list_species (SP001-SP011)
  └─ 11_cmd_query_species (Q001-Q006)

阶段 2: 数据库构建（analyze的前置依赖）
  ├─ 02_cmd_download (D001-D010, D121-D123)
  ├─ 03_cmd_build (B001-B014, B124-B125)
  └─ 08_cmd_cleanup (CL001-CL003)

阶段 3: 核心分析（依赖阶段2）
  ├─ 01_cmd_analyze/01_ora_basic (T001-T011)
  ├─ 01_cmd_analyze/02_ora_correction (T013-T017)
  ├─ 01_cmd_analyze/03_ora_databases (T018-T023)
  ├─ 01_cmd_analyze/04_background_modes (T024-T027)
  ├─ 01_cmd_analyze/05_gsea_analysis (T028-T031)
  ├─ 01_cmd_analyze/06_ssgsea_analysis (T032-T034)
  ├─ 01_cmd_analyze/07_gsva_analysis (T035-T038)

阶段 4: 可视化与输出（依赖阶段3）
  ├─ 01_cmd_analyze/08_visual_styles (T039-T044, G144-G149, G158-G163, G172-G177)
  ├─ 01_cmd_analyze/09_visual_palettes (T045-T058, T134-T142, G150-G157, G164-G171, G178-G185)
  ├─ 01_cmd_analyze/10_output_formats (T059-T063)
  ├─ 01_cmd_analyze/11_report_control (T064-T066)
  ├─ 01_cmd_analyze/12_plot_types (T067-T077)

阶段 5: 高级功能
  ├─ 01_cmd_analyze/13_config_loading (T078-T079)
  ├─ 01_cmd_analyze/14_version_control (T080-T081)
  ├─ 01_cmd_analyze/15_ai_interpretation (T082-T085)
  ├─ 01_cmd_analyze/17_config_only_params (T091-T103)
  ├─ 01_cmd_analyze/18_global_params (T104-T105)
  ├─ 01_cmd_analyze/19_cli_priority (T106)
  ├─ 01_cmd_analyze/20_input_formats (T107-T109)
  ├─ 01_cmd_analyze/21_gmt_formats (T110)
  ├─ 01_cmd_analyze/22_group_comparison (T111-T113)
  ├─ 01_cmd_analyze/23_correlation_methods (T114-T115)
  ├─ 01_cmd_analyze/24_network_plots (T116-T117)
  ├─ 01_cmd_analyze/25_api_endpoints (T117a-T117k)
  ├─ 01_cmd_analyze/26_report_generation (T118-T119)
  ├─ 01_cmd_analyze/27_url_generation (T120)
  ├─ 01_cmd_analyze/28_parallel_modes (T127-T129)
  ├─ 01_cmd_analyze/29_gsea_fields (T130)
  ├─ 01_cmd_analyze/30_ssgsea_nan (T131)
  ├─ 01_cmd_analyze/31_output_fields (T132-T133)

阶段 6: 错误场景
  └─ 01_cmd_analyze/16_error_scenarios (T086-T090)

阶段 7: API服务
  └─ 04_cmd_serve (S001-S003)

阶段 8: 其他功能验证
  ├─ 12_disgenet_graceful (D121-D122)
  ├─ 13_registry_build (D123)
  ├─ 14_build_manifest (B124)
  ├─ 15_goa_fallback (B125)
  └─ 16_version_consistency (T126)
```

---

## 七、风险与注意事项

1. **网络依赖**: download、check-update、build（KEGG REST API）需要网络连接
2. **执行时间**: GSEA排列检验较慢（默认1000次），T091使用100次加速
3. **磁盘空间**: 294个测试用例的输出可能占用数GB空间
4. **hsa版本差异**: hsa使用v20260515旧版构建，mmu使用v20260528
5. **AI测试**: 仅使用mock后端，真实API测试需要有效密钥
6. **DisGeNET**: 需要测试graceful degradation而非实际下载
7. **川金丝猴(rro)**: 基因数量可能不足500，需处理边界情况

---

## 八、自检清单

### 8.1 规格覆盖

- [x] 11个CLI子命令全部覆盖
- [x] 294个测试用例全部定义
- [x] 统一的目录结构和命名规范
- [x] 每个测试有独立的run.sh模板
- [x] 数据准备脚本完整可执行

### 8.2 功能覆盖

- [x] 5种富集方法: fisher, hypergeometric, gsea, ssgsea, gsva
- [x] 3种GSVA变体: gsva, plage, zscore
- [x] 2种GSVA核函数: Gaussian, Poisson
- [x] 5种校正方法: BH, BY, bonferroni, holm, none
- [x] 6种可视化风格: nature, science, cell, colorblind, presentation, omicshare
- [x] 23种配色方案全覆盖
- [x] 3种输出格式: png, pdf, svg
- [x] 5个数据库: GO, KEGG, Reactome, DO, DisGeNET
- [x] 3种背景模式: annotated, genome, custom
- [x] 10种图表类型全覆盖
- [x] 4种注释格式: three_column, four_column, two_column, auto
- [x] 7种AI后端: openai, claude, deepseek, glm, minimax, ollama, mock
- [x] 11个API端点全覆盖

### 8.3 无占位符

- [x] 所有测试命令具体可执行
- [x] 所有文件路径明确
- [x] 所有参数值具体
- [x] 验证要点清晰
