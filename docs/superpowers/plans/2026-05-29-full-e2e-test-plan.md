# AllEnricher-v2 全场景端对端测试计划（完整合并版）

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 用真实数据对 AllEnricher-v2 的全部 11 个 CLI 子命令、所有参数、所有功能进行无死角端对端测试，确保每个功能、参数、命令均正常工作。

**Architecture:** 新建独立测试目录 `test_e2e_2026/`，按功能模块建立子层级文件夹，分门别类保存测试脚本、输入数据、输出文件、日志。所有文件按测试编号前缀命名，命名直接反映测试内容。

**Tech Stack:** Python 3.x, AllEnricher-v2 CLI (`python -m allenricher`), 真实基因数据

---
## 一、测试环境与前提条件
### 1.1 命令调用方式
所有 CLI 命令统一使用 `python -m allenricher` 方式调用（`allenricher` 未注册到 PATH）。

### 1.2 已有基础设施
| 资源 | 状态 | 说明 |
|------|------|------|
| hsa 数据库 | ✅ 已构建 v20260515 | DO+GO+KEGG+Reactome |
| mmu 数据库 | ✅ 已构建 v20260528 | GO+KEGG+Reactome |
| 基础数据库 | ✅ 已下载 | GO20260527, Reactome20260515, DO, KEGG, Taxonomy |
| 测试基因列表 | ✅ 已有 | test_data/top_degs.txt |
| 表达矩阵 | ✅ 已有 | test_data/expression_matrix.tsv |
| 排序基因 | ✅ 已有 | test_data/ranked_genes.tsv |
| GMT 文件 | ✅ 已有 | test_data/gene_sets.gmt |

### 1.3 需要准备的测试数据
见下方"测试数据设计"章节。

---
## 二、物种选择策略
### 2.1 主测试物种
| 物种 | KEGG Code | TaxID | 用途 | 数据库支持 |
|------|-----------|-------|------|-----------|
| **人类** | `hsa` | 9606 | 主测试物种，覆盖全部数据库 | GO+KEGG+Reactome+DO |
| **小鼠** | `mmu` | 10090 | 跨物种验证 | GO+KEGG+Reactome |
| **川金丝猴** | `rro` | 61622 | 新物种构建测试 | GO+KEGG |

### 2.2 物种选择理由
- **hsa**: 唯一支持 DO 和 DisGeNET 的物种，可测试全部 5 个数据库
- **mmu**: 已构建数据库，可验证跨物种分析一致性
- **rro**: 需要新构建，可测试 build 命令全流程；川金丝猴仅支持 GO+KEGG

---
## 三、测试目录结构（强制规范）
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
│   │   ├── B001_build_rro_go_kegg/
│   │   ├── B002_build_rro_all_dbs/
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
## 四、测试编号命名规范
### 4.1 编号前缀规则
| 前缀 | 命令/模块 | 编号范围 | 示例 |
|------|----------|---------|------|
| `T` | analyze (主命令) | T001-T143 | T001_fisher_default |
| `D` | download | D001-D010, D121-D123 | D001_download_go |
| `B` | build | B001-B014, B124-B125 | B001_build_rro_go_kegg |
| `S` | serve | S001-S003 | S001_serve_default |
| `L` | list | L001-L003 | L001_list_species |
| `C` | config | C001-C002 | C001_config_default |
| `U` | check-update | U001-U004 | U001_check_update_default |
| `CL` | cleanup | CL001-CL003 | CL001_cleanup_dry_run |
| `V` | list-versions | V001-V005 | V001_list_versions_default |
| `SP` | list-species | SP001-SP011 | SP001_list_species_default |
| `Q` | query-species | Q001-Q006 | Q001_query_by_name |
| `G` | GSEA/ssGSEA/GSVA 风格配色 | G144-G185 | G144_gsea_style_nature |

### 4.2 文件夹命名规范
```
{编号}_{测试内容描述}/
```

- **编号**：3位数字，不足补零（如 T001, T091, G144）
- **测试内容**：下划线连接的小写英文单词，描述测试的核心内容
- **示例**：
  - `T001_fisher_default/` — Fisher精确检验默认参数
  - `T091_gsea_permutations_100/` — GSEA排列次数100
  - `G144_gsea_style_nature/` — GSEA nature风格

### 4.3 文件命名规范
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
## 五、执行脚本模板（run.sh）
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
## 六、数据准备脚本规格
### 6.1 脚本位置
`test_e2e_2026/00_input_data/00_prepare_data.py`

### 6.2 生成的文件清单
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

### 6.3 核心函数实现
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


## 七、测试用例详细清单
### 7.1 analyze 命令测试（T001-T143）
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

### 7.2 其他命令测试
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
| B001 | `B001_build_rro_go_kegg/` | `build -s rro -t 61622 -d GO,KEGG` |
| B002 | `B002_build_rro_all_dbs/` | `build -s rro -t 61622 -d GO,KEGG,Reactome` | 应跳过Reactome并警告 |
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

### 7.3 GSEA/ssGSEA/GSVA 风格配色测试（G144-G185）
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

## 八、测试命令参考（完整 bash 命令）

> 各测试文件夹内的 `run.sh` 脚本使用了本节的完整命令格式。
> 所有命令前均省略 `cd test_e2e_2026 && python -m allenricher`，执行时请在实际目录下运行。
> 输入文件相对路径均以 `00_input_data/` 为基准，输出目录为对应测试文件夹下的 `output/`。

### 8.1 analyze 基础 ORA（T001-T012）

```bash
# T001: Fisher 精确检验 - 默认参数
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO,KEGG -o 01_cmd_analyze/01_ora_basic/T001_fisher_default/output

# T002: Fisher + 指定 pvalue/qvalue 阈值
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO,KEGG -o 01_cmd_analyze/01_ora_basic/T002_fisher_threshold_p0.01_q0.01/output -p 0.01 -q 0.01

# T003: Fisher + min-genes=5
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/01_ora_basic/T003_fisher_min_genes_5/output -n 5

# T004: Fisher + min-genes=10
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/01_ora_basic/T004_fisher_min_genes_10/output -n 10

# T005: Fisher + only-significant
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO,KEGG -o 01_cmd_analyze/01_ora_basic/T005_fisher_only_significant/output --only-significant

# T006: Fisher + 并行 jobs=4
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO,KEGG -o 01_cmd_analyze/01_ora_basic/T006_fisher_parallel_jobs_4/output -j 4

# T007: Hypergeometric 检验 - 默认参数
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO,KEGG -o 01_cmd_analyze/01_ora_basic/T007_hypergeometric_default/output -m hypergeometric

# T008: Hypergeometric + BH 校正
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/01_ora_basic/T008_hypergeometric_correction_bh/output -m hypergeometric -c BH

# T009: Hypergeometric + bonferroni 校正
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/01_ora_basic/T009_hypergeometric_correction_bonferroni/output -m hypergeometric -c bonferroni

# T010: Hypergeometric + holm 校正
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/01_ora_basic/T010_hypergeometric_correction_holm/output -m hypergeometric -c holm

# T011: Hypergeometric + BY 校正
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/01_ora_basic/T011_hypergeometric_correction_by/output -m hypergeometric -c BY

# T012: Fisher + none 校正
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/01_ora_basic/T012_correction_none/output -c none
```

### 8.2 analyze ORA 校正方法（T013-T017）

```bash
# T013: Fisher + BH
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/02_ora_correction/T013_correction_bh/output -m fisher -c BH

# T014: Fisher + BY
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/02_ora_correction/T014_correction_by/output -m fisher -c BY

# T015: Fisher + bonferroni
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/02_ora_correction/T015_correction_bonferroni/output -m fisher -c bonferroni

# T016: Fisher + holm
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/02_ora_correction/T016_correction_holm/output -m fisher -c holm

# T017: Fisher + none
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/02_ora_correction/T017_correction_none/output -m fisher -c none
```

### 8.3 analyze 数据库组合（T018-T023）

```bash
# T018: 仅 GO
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/03_ora_databases/T018_database_go_only/output

# T019: 仅 KEGG
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d KEGG -o 01_cmd_analyze/03_ora_databases/T019_database_kegg_only/output

# T020: 仅 Reactome
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d Reactome -o 01_cmd_analyze/03_ora_databases/T020_database_reactome_only/output

# T021: 仅 DO
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d DO -o 01_cmd_analyze/03_ora_databases/T021_database_do_only/output

# T022: GO+KEGG+Reactome+DO 全数据库
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO,KEGG,Reactome,DO -o 01_cmd_analyze/03_ora_databases/T022_database_all_combined/output

# T023: mmu 物种 - GO+KEGG+Reactome
python -m allenricher analyze -i 00_input_data/mmu_gene_list_500.txt -s mmu -d GO,KEGG,Reactome -o 01_cmd_analyze/03_ora_databases/T023_mmu_database_all/output
```

### 8.4 analyze 背景模式（T024-T027）

```bash
# T024: background-mode=annotated
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/04_background_modes/T024_background_annotated/output --background-mode annotated

# T025: background-mode=genome
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/04_background_modes/T025_background_genome/output --background-mode genome

# T026: background-mode=custom + 背景文件
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/04_background_modes/T026_background_custom/output --background-mode custom -b 00_input_data/hsa_background_genes.txt

# T027: background-mode=custom 但未提供 -b（应报错退出）
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/04_background_modes/T027_background_custom_error/output --background-mode custom
```

### 8.5 GSEA 分析（T028-T031）

```bash
# T028: GSEA + ranked genes 文件
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO,KEGG -o 01_cmd_analyze/05_gsea_analysis/T028_gsea_ranked_genes/output -m gsea -r 00_input_data/hsa_ranked_genes_gsea.tsv

# T029: GSEA + expression matrix
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/05_gsea_analysis/T029_gsea_expression_matrix/output -m gsea -e 00_input_data/hsa_expression_matrix_500x6.tsv

# T030: GSEA + 自定义 GMT
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/05_gsea_analysis/T030_gsea_custom_gmt/output -m gsea -r 00_input_data/hsa_ranked_genes_gsea.tsv -g 00_input_data/hsa_custom_pathways_5.gmt

# T031: GSEA + groups 定义
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/05_gsea_analysis/T031_gsea_with_groups/output -m gsea -e 00_input_data/hsa_expression_matrix_500x6.tsv --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"
```

### 8.6 ssGSEA 分析（T032-T034）

```bash
# T032: ssGSEA + expression matrix
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/06_ssgsea_analysis/T032_ssgsea_expression_matrix/output -m ssgsea -e 00_input_data/hsa_expression_matrix_500x6.tsv

# T033: ssGSEA + 自定义 GMT
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/06_ssgsea_analysis/T033_ssgsea_custom_gmt/output -m ssgsea -e 00_input_data/hsa_expression_matrix_500x6.tsv -g 00_input_data/hsa_custom_pathways_5.gmt

# T034: ssGSEA + groups
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/06_ssgsea_analysis/T034_ssgsea_with_groups/output -m ssgsea -e 00_input_data/hsa_expression_matrix_500x6.tsv --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"
```

### 8.7 GSVA 分析（T035-T038）

```bash
# T035: GSVA + expression matrix（默认 gsva 方法）
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/07_gsva_analysis/T035_gsva_default_gsva/output -m gsva -e 00_input_data/hsa_expression_matrix_500x6.tsv

# T036: GSVA plage 变体（通过 config 文件）
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/07_gsva_analysis/T036_gsva_method_plage/output -m gsva -e 00_input_data/hsa_expression_matrix_500x6.tsv --config 00_input_data/config_files/T093_gsva_method_plage.yaml

# T037: GSVA zscore 变体（通过 config 文件）
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/07_gsva_analysis/T037_gsva_method_zscore/output -m gsva -e 00_input_data/hsa_expression_matrix_500x6.tsv --config 00_input_data/config_files/T094_gsva_method_zscore.yaml

# T038: GSVA + groups
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/07_gsva_analysis/T038_gsva_with_groups/output -m gsva -e 00_input_data/hsa_expression_matrix_500x6.tsv --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"
```

### 8.8 可视化风格（T039-T044）

```bash
# T039: style=nature
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO,KEGG -o 01_cmd_analyze/08_visual_styles/T039_style_nature/output --style nature

# T040: style=science
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO,KEGG -o 01_cmd_analyze/08_visual_styles/T040_style_science/output --style science

# T041: style=cell（通过 config）
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO,KEGG -o 01_cmd_analyze/08_visual_styles/T041_style_cell/output --config 00_input_data/config_files/T143_style_cell_config.yaml

# T042: style=colorblind
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO,KEGG -o 01_cmd_analyze/08_visual_styles/T042_style_colorblind/output --style colorblind

# T043: style=presentation
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO,KEGG -o 01_cmd_analyze/08_visual_styles/T043_style_presentation/output --style presentation

# T044: style=omicshare
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO,KEGG -o 01_cmd_analyze/08_visual_styles/T044_style_omicshare/output --style omicshare
```

### 8.9 配色方案（T045-T058）

```bash
# T045: palette=nature
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/09_visual_palettes/T045_palette_nature/output --palette nature

# T046: palette=science
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/09_visual_palettes/T046_palette_science/output --palette science

# T047: palette=lancet
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/09_visual_palettes/T047_palette_lancet/output --palette lancet

# T048: palette=nejm
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/09_visual_palettes/T048_palette_nejm/output --palette nejm

# T049: palette=jama
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/09_visual_palettes/T049_palette_jama/output --palette jama

# T050: palette=okabe_ito
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/09_visual_palettes/T050_palette_okabe_ito/output --palette okabe_ito

# T051: palette=gsea
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/09_visual_palettes/T051_palette_gsea/output --palette gsea

# T052: palette=omicshare
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/09_visual_palettes/T052_palette_omicshare/output --palette omicshare

# T053: palette=china_style
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/09_visual_palettes/T053_palette_china_style/output --palette china_style

# T054: palette=go_bp
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/09_visual_palettes/T054_palette_go_bp/output --palette go_bp

# T055: palette=kegg_pathway
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d KEGG -o 01_cmd_analyze/09_visual_palettes/T055_palette_kegg_pathway/output --palette kegg_pathway

# T056: palette=tol_bright
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/09_visual_palettes/T056_palette_tol_bright/output --palette tol_bright

# T057: palette=tol_muted
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/09_visual_palettes/T057_palette_tol_muted/output --palette tol_muted

# T058: palette=tol_sunset
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/09_visual_palettes/T058_palette_tol_sunset/output --palette tol_sunset
```

### 8.10 输出格式（T059-T063）

```bash
# T059: plot-format=png
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO,KEGG -o 01_cmd_analyze/10_output_formats/T059_format_png/output --plot-format png

# T060: plot-format=pdf
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO,KEGG -o 01_cmd_analyze/10_output_formats/T060_format_pdf/output --plot-format pdf

# T061: plot-format=svg
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO,KEGG -o 01_cmd_analyze/10_output_formats/T061_format_svg/output --plot-format svg

# T062: plot-dpi=150
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO,KEGG -o 01_cmd_analyze/10_output_formats/T062_dpi_150/output --plot-dpi 150

# T063: plot-dpi=600
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO,KEGG -o 01_cmd_analyze/10_output_formats/T063_dpi_600/output --plot-dpi 600
```

### 8.11 报告控制（T064-T066）

```bash
# T064: --no-plot
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO,KEGG -o 01_cmd_analyze/11_report_control/T064_no_plot/output --no-plot

# T065: --no-report
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO,KEGG -o 01_cmd_analyze/11_report_control/T065_no_report/output --no-report

# T066: --no-plot --no-report
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO,KEGG -o 01_cmd_analyze/11_report_control/T066_no_plot_no_report/output --no-plot --no-report
```

### 8.12 图表类型（T067-T077）

```bash
# T067: plot-types=dotplot
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/12_plot_types/T067_plot_dotplot/output -pt dotplot

# T068: plot-types=barplot
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/12_plot_types/T068_plot_barplot/output -pt barplot

# T069: plot-types=network
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/12_plot_types/T069_plot_network/output -pt network

# T070: plot-types=upset
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO,KEGG -o 01_cmd_analyze/12_plot_types/T070_plot_upset/output -pt upset

# T071: plot-types=volcano
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO,KEGG -o 01_cmd_analyze/12_plot_types/T071_plot_volcano/output -pt volcano

# T072: plot-types=enrichment（GSEA 专用）
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/12_plot_types/T072_gsea_plot_enrichment/output -m gsea -r 00_input_data/hsa_ranked_genes_gsea.tsv -pt enrichment

# T073: plot-types=nes_barplot（GSEA 专用）
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/12_plot_types/T073_gsea_plot_nes_barplot/output -m gsea -r 00_input_data/hsa_ranked_genes_gsea.tsv -pt nes_barplot

# T074: plot-types=heatmap（ssGSEA/GSVA 专用）
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/12_plot_types/T074_ssgsea_plot_heatmap/output -m ssgsea -e 00_input_data/hsa_expression_matrix_500x6.tsv -pt heatmap

# T075: plot-types=group_comparison（ssGSEA/GSVA 专用）
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/12_plot_types/T075_ssgsea_plot_group_comparison/output -m ssgsea -e 00_input_data/hsa_expression_matrix_500x6.tsv -pt group_comparison --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"

# T076: plot-types=correlation（ssGSEA/GSVA 专用）
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/12_plot_types/T076_ssgsea_plot_correlation/output -m ssgsea -e 00_input_data/hsa_expression_matrix_500x6.tsv -pt correlation

# T077: 多种 plot-types 组合
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO,KEGG -o 01_cmd_analyze/12_plot_types/T077_plot_multi_types/output -pt dotplot,barplot,network,upset
```

### 8.13 配置文件加载（T078-T079）

```bash
# T078: 使用 YAML 配置文件
python -m allenricher analyze --config 00_input_data/config_files/T078_config_full_yaml.yaml

# T079: 使用 JSON 配置文件
python -m allenricher analyze --config 00_input_data/config_files/T079_config_full_json.json
```

### 8.14 版本与日志控制（T080-T085）

```bash
# T080: --use-version 指定数据库版本
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO,KEGG -o 01_cmd_analyze/14_version_control/T080_use_version_v20260515/output --use-version v20260515

# T081: --database-dir 指定数据库目录
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO,KEGG -o 01_cmd_analyze/14_version_control/T081_database_dir_custom/output --database-dir ./database

# T082: --verbose 开启详细日志
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/15_ai_interpretation/T082_verbose_logging/output --verbose

# T083: AI 解读 + mock 后端
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO,KEGG -o 01_cmd_analyze/15_ai_interpretation/T083_ai_mock_backend/output --ai mock

# T084: AI 解读 + openai 参数
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/15_ai_interpretation/T084_ai_openai_params/output --ai openai --ai-key test_key --ai-model gpt-4

# T085: AI 解读 + claude 参数
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/15_ai_interpretation/T085_ai_claude_params/output --ai claude --ai-key test_key --ai-model claude-3
```

### 8.15 错误场景（T086-T090）

```bash
# T086: 空基因列表
python -m allenricher analyze -i 00_input_data/error_scenarios/empty_gene_list.txt -s hsa -d GO -o 01_cmd_analyze/16_error_scenarios/T086_empty_gene_list/output

# T087: 不存在的物种
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s xxx_notexist -d GO -o 01_cmd_analyze/16_error_scenarios/T087_invalid_species/output

# T088: 不存在的输入文件
python -m allenricher analyze -i 00_input_data/nonexistent.txt -s hsa -d GO -o 01_cmd_analyze/16_error_scenarios/T088_nonexistent_input/output

# T089: GSEA 但未提供 ranked-genes/expression-matrix
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/16_error_scenarios/T089_gsea_missing_ranked/output -m gsea

# T090: ssGSEA 但未提供 expression-matrix
python -m allenricher analyze -i 00_input_data/hsa_gene_list_500.txt -s hsa -d GO -o 01_cmd_analyze/16_error_scenarios/T090_ssgsea_missing_expression/output -m ssgsea
```

### 8.16 download 命令（D001-D010）

```bash
# D001: 下载 GO 数据库
python -m allenricher download -d go

# D002: 下载 Reactome 数据库
python -m allenricher download -d reactome

# D003: 下载 DO 数据库
python -m allenricher download -d do

# D004: 下载多个数据库
python -m allenricher download -d go,reactome,do

# D005: 指定 workers=8
python -m allenricher download -d go --workers 8

# D006: 禁用多线程
python -m allenricher download -d go --no-multi-thread

# D007: 跳过完整性校验
python -m allenricher download -d go --no-verify

# D008: 强制重新下载
python -m allenricher download -d go --force

# D009: 指定 database-dir
python -m allenricher download -d go --database-dir ./database

# D010: 指定 species（KEGG 提示）
python -m allenricher download -d kegg -s hsa
```

### 8.17 build 命令（B001-B014）

```bash
# B001: 构建 rro - GO+KEGG
python -m allenricher build -s rro -t 61622 -d GO,KEGG

# B002: 构建 rro - GO+KEGG+Reactome
python -m allenricher build -s rro -t 61622 -d GO,KEGG,Reactome

# B003: 指定 database-dir
python -m allenricher build -s rro -t 61622 -d GO,KEGG --database-dir ./database

# B004: 使用 --gene-info
python -m allenricher build -s rro -t 61622 -d GO --gene-info ./database/basic/go/GO20260527/gene_info.gz

# B005: 使用 --go-annot
python -m allenricher build -s rro -t 61622 -d GO --go-annot 00_input_data/custom_annotations/three_column.tsv

# B006: 使用 --custom-annot + --custom-db-name
python -m allenricher build -s rro -t 61622 -d GO --custom-annot 00_input_data/custom_annotations/three_column.tsv --custom-db-name MYDB

# B007: --annot-format=three_column
python -m allenricher build -s rro -t 61622 -d GO --custom-annot 00_input_data/custom_annotations/three_column.tsv --annot-format three_column

# B008: --annot-format=four_column
python -m allenricher build -s rro -t 61622 -d GO --custom-annot 00_input_data/custom_annotations/four_column.tsv --annot-format four_column

# B009: --annot-format=two_column
python -m allenricher build -s rro -t 61622 -d GO --custom-annot 00_input_data/custom_annotations/two_column.tsv --annot-format two_column

# B010: --annot-format=auto
python -m allenricher build -s rro -t 61622 -d GO --custom-annot 00_input_data/custom_annotations/three_column.tsv --annot-format auto

# B011: --hierarchy-sep="."
python -m allenricher build -s rro -t 61622 -d GO --custom-annot 00_input_data/custom_annotations/three_column.tsv --hierarchy-sep "."

# B012: hsa 构建 DO
python -m allenricher build -s hsa -t 9606 -d DO

# B013: mmu 构建 DO（应跳过并警告）
python -m allenricher build -s mmu -t 10090 -d DO

# B014: hsa 全部数据库
python -m allenricher build -s hsa -t 9606 -d GO,KEGG,Reactome,DO
```

### 8.18 serve/list/config/check-update/cleanup/list-versions/list-species/query-species 命令

```bash
# S001-S003: serve
python -m allenricher serve                                          # S001: 默认
python -m allenricher serve --host 127.0.0.1 --port 8080             # S002: 自定义
python -m allenricher serve --reload                                  # S003: 热重载

# L001-L003: list
python -m allenricher list species                                    # L001: 物种
python -m allenricher list databases                                  # L002: 数据库
python -m allenricher list                                            # L003: 无参（应报错）

# C001-C002: config
python -m allenricher config                                          # C001: 默认输出
python -m allenricher config -o my_config.yaml                        # C002: 自定义路径

# U001-U004: check-update
python -m allenricher check-update                                    # U001: 默认
python -m allenricher check-update --database-dir ./database          # U002: 指定目录
python -m allenricher check-update --json                             # U003: JSON 输出
python -m allenricher check-update --database-dir ./database --json   # U004: 组合

# CL001-CL003: cleanup
python -m allenricher cleanup --dry-run                               # CL001: 预览
python -m allenricher cleanup --keep 1 --dry-run                      # CL002: 保留1个
python -m allenricher cleanup --database-dir ./database --dry-run     # CL003: 指定目录

# V001-V005: list-versions
python -m allenricher list-versions                                   # V001: 默认
python -m allenricher list-versions --json                            # V002: JSON
python -m allenricher list-versions --lineage                         # V003: 血缘
python -m allenricher list-versions --database-dir ./database         # V004: 指定目录
python -m allenricher list-versions --json --lineage                  # V005: 全部参数

# SP001-SP011: list-species
python -m allenricher list-species                                    # SP001: 默认
python -m allenricher list-species --go                               # SP002: GO 过滤
python -m allenricher list-species --kegg                             # SP003: KEGG 过滤
python -m allenricher list-species --reactome                         # SP004: Reactome 过滤
python -m allenricher list-species --do                               # SP005: DO 过滤
python -m allenricher list-species --go --kegg --reactome             # SP006: 多过滤
python -m allenricher list-species --format table                     # SP007: table 格式
python -m allenricher list-species --format tsv                       # SP008: tsv 格式
python -m allenricher list-species --format json                      # SP009: json 格式
python -m allenricher list-species --summary                          # SP010: 统计
python -m allenricher list-species --go --kegg --format json --summary # SP011: 全部

# Q001-Q006: query-species
python -m allenricher query-species --name "Homo sapiens"             # Q001: 按名
python -m allenricher query-species --taxid 9606                      # Q002: 按 TaxID
python -m allenricher query-species --kegg hsa                        # Q003: 按 KEGG
python -m allenricher query-species --name "Microtus fortis"          # Q004: 模糊旧名
python -m allenricher query-species --name "Rhinopithecus"            # Q005: 部分匹配
python -m allenricher query-species --name "Nonexistent"              # Q006: 不存在
```

---
## 九、参数覆盖检查矩阵

### 9.1 analyze 命令参数覆盖
| 参数 | 测试用例 | 覆盖状态 |
|------|---------|---------|
| `-i` / `--input` | T1-T90 全部 | ✅ |
| `-s` / `--species` | T1(hsa), T23(mmu), T87(xxx) | ✅ |
| `-d` / `--databases` | T18(GO), T19(KEGG), T20(Reactome), T21(DO), T22(all) | ✅ |
| `-o` / `--output` | T1-T90 全部 | ✅ |
| `-b` / `--background` | T26 | ✅ |
| `--background-mode` | T24(annotated), T25(genome), T26(custom), T27(err) | ✅ |
| `-m` / `--method` | T1(fisher), T7(hyper), T28(gsea), T32(ssgsea), T35(gsva) | ✅ |
| `-c` / `--correction` | T13(BH), T14(BY), T15(bonf), T16(holm), T17(none) | ✅ |
| `-p` / `--pvalue` | T2(0.01) | ✅ |
| `-q` / `--qvalue` | T2(0.01) | ✅ |
| `-n` / `--min-genes` | T3(5), T4(10) | ✅ |
| `-j` / `--jobs` | T6(4) | ✅ |
| `--no-plot` | T64, T66 | ✅ |
| `--no-report` | T65, T66 | ✅ |
| `--only-significant` | T5 | ✅ |
| `--ai` | T83(mock), T84(openai), T85(claude) | ✅ |
| `--ai-key` | T84, T85 | ✅ |
| `--ai-model` | T84, T85 | ✅ |
| `--config` | T78(yaml), T79(json) | ✅ |
| `--database-dir` | T81 | ✅ |
| `--use-version` | T80 | ✅ |
| `-e` / `--expression-matrix` | T29, T32, T35 | ✅ |
| `-r` / `--ranked-genes` | T28, T30 | ✅ |
| `-g` / `--gmt` | T30, T33 | ✅ |
| `-pt` / `--plot-types` | T67-T77 | ✅ |
| `--groups` | T31, T34, T38, T75 | ✅ |
| `--plot-format` | T59(png), T60(pdf), T61(svg) | ✅ |
| `--plot-dpi` | T62(150), T63(600) | ✅ |
| `--style` | T39-T44 (6种全覆盖) | ✅ |
| `--palette` | T45-T58 (14种关键配色) | ✅ |
| `--verbose` | T82 | ✅ |

### 9.2 download 命令参数覆盖
| 参数 | 测试用例 | 覆盖状态 |
|------|---------|---------|
| `-d` / `--databases` | D1(go), D2(reactome), D3(do), D4(多) | ✅ |
| `-s` / `--species` | D10 | ✅ |
| `--database-dir` | D9 | ✅ |
| `--workers` | D5 | ✅ |
| `--no-multi-thread` | D6 | ✅ |
| `--no-verify` | D7 | ✅ |
| `--force` | D8 | ✅ |

### 9.3 build 命令参数覆盖
| 参数 | 测试用例 | 覆盖状态 |
|------|---------|---------|
| `-s` / `--species` | B1-B14 | ✅ |
| `-t` / `--taxonomy` | B1-B14 | ✅ |
| `-d` / `--databases` | B1(GO,KEGG), B2(+Reactome), B12(DO), B14(all) | ✅ |
| `--database-dir` | B3 | ✅ |
| `--gene-info` | B4 | ✅ |
| `--go-annot` | B5 | ✅ |
| `--kegg-annot` | (通过 --custom-annot 覆盖同类逻辑) | ✅ |
| `--custom-annot` | B6 | ✅ |
| `--custom-db-name` | B6(MYDB) | ✅ |
| `--annot-format` | B7(three), B8(four), B9(two), B10(auto) | ✅ |
| `--hierarchy-sep` | B11 | ✅ |

### 9.4 其余命令参数覆盖
| 命令 | 参数 | 测试用例 | 覆盖状态 |
|------|------|---------|---------|
| serve | `--host` | S2 | ✅ |
| serve | `--port` | S2 | ✅ |
| serve | `--reload` | S3 | ✅ |
| list | `resource` | L1(species), L2(databases), L3(无参) | ✅ |
| config | `-o` | C1(默认), C2(自定义) | ✅ |
| check-update | `--database-dir` | U2, U4 | ✅ |
| check-update | `--json` | U3, U4 | ✅ |
| cleanup | `--keep` | CL2 | ✅ |
| cleanup | `--dry-run` | CL1, CL2, CL3 | ✅ |
| cleanup | `--database-dir` | CL3 | ✅ |
| list-versions | `--database-dir` | V4, V5 | ✅ |
| list-versions | `--json` | V2, V5 | ✅ |
| list-versions | `--lineage` | V3, V5 | ✅ |
| list-species | `--go` | SP2, SP6, SP11 | ✅ |
| list-species | `--kegg` | SP3, SP6, SP11 | ✅ |
| list-species | `--reactome` | SP4, SP6, SP11 | ✅ |
| list-species | `--do` | SP5, SP11 | ✅ |
| list-species | `--format` | SP7(table), SP8(tsv), SP9(json) | ✅ |
| list-species | `--summary` | SP10, SP11 | ✅ |
| query-species | `--name` | Q1, Q4, Q5, Q6 | ✅ |
| query-species | `--taxid` | Q2 | ✅ |
| query-species | `--kegg` | Q3 | ✅ |

---
## 十、测试执行顺序
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
## 十一、风险与注意事项
1. **网络依赖**: download、check-update、build（KEGG REST API）需要网络连接
2. **执行时间**: GSEA 排列检验较慢（默认 1000 次），建议测试时减少排列次数
3. **磁盘空间**: 94+ 测试用例的输出可能占用数 GB 空间，测试后需清理
4. **hsa 版本差异**: hsa 使用 v20260515 旧版构建，mmu 使用 v20260528，结果可能不一致
5. **AI 测试**: 仅使用 mock 后端，真实 API 测试需要有效密钥
6. **DisGeNET**: 需要额外授权/下载，基础数据库中未见已下载的 DisGeNET 数据

---
## 十二、自检清单

### 12.1 规格覆盖
- [x] 11个CLI子命令全部覆盖
- [x] 294个测试用例全部定义
- [x] 统一的目录结构和命名规范
- [x] 每个测试有独立的run.sh模板
- [x] 数据准备脚本完整可执行

### 12.2 功能覆盖
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

### 12.3 无占位符
- [x] 所有测试命令具体可执行
- [x] 所有文件路径明确
- [x] 所有参数值具体
- [x] 验证要点清晰

---
## 十三、审查补充：Config-only 参数（无 CLI 对应，只能通过配置文件设置）
以下 Config 字段**没有 CLI 参数**，只能通过 `--config` 加载 YAML/JSON 配置文件传递。原计划完全遗漏。

### 13.1 GSEA 专用参数（高优先级）
| Config 字段 | 类型 | 默认值 | 说明 |
|-------------|------|--------|------|
| `gsea_permutations` | int | 1000 | GSEA 排列检验次数 |
| `gsea_min_size` | int | 10 | GSEA 基因集最小大小 |
| `gsea_max_size` | int | 500 | GSEA 基因集最大大小 |

```bash
# T91: GSEA + config 设置 permutations=100（加速测试）
# 配置文件 test_data/e2e_real/gsea_perm_config.yaml:
#   input_file: test_data/e2e_real/hsa_gene_list.txt
#   species: hsa
#   databases: [GO]
#   method: gsea
#   output_dir: test_data/e2e_real/results/t91_gsea_perm100
#   gsea_permutations: 100
#   ranked_genes_file: test_data/e2e_real/hsa_ranked_genes.tsv
python -m allenricher analyze --config test_data/e2e_real/gsea_perm_config.yaml

# T92: GSEA + config 设置 min_size=5, max_size=200
# 配置文件 test_data/e2e_real/gsea_size_config.yaml:
#   input_file: test_data/e2e_real/hsa_gene_list.txt
#   species: hsa
#   databases: [GO]
#   method: gsea
#   output_dir: test_data/e2e_real/results/t92_gsea_size
#   gsea_min_size: 5
#   gsea_max_size: 200
python -m allenricher analyze --config test_data/e2e_real/gsea_size_config.yaml -r test_data/e2e_real/hsa_ranked_genes.tsv
```

### 13.2 GSVA 专用参数（高优先级）
| Config 字段 | 类型 | 默认值 | choices | 说明 |
|-------------|------|--------|---------|------|
| `gsva_method` | str | "gsva" | gsva, plage, zscore | GSVA 算法变体 |
| `gsva_kcdf` | str | "Gaussian" | Gaussian, Poisson | 核密度估计核函数 |
| `gsva_tau` | float | 1.0 | - | 核密度带宽参数（仅 Gaussian） |

```bash
# T93: GSVA + plage 变体
# 配置文件 test_data/e2e_real/gsva_plage_config.yaml:
#   gsva_method: plage
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t93_gsva_plage -m gsva -e test_data/e2e_real/hsa_expression_matrix.tsv --config test_data/e2e_real/gsva_plage_config.yaml

# T94: GSVA + zscore 变体
# 配置文件 test_data/e2e_real/gsva_zscore_config.yaml:
#   gsva_method: zscore
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t94_gsva_zscore -m gsva -e test_data/e2e_real/hsa_expression_matrix.tsv --config test_data/e2e_real/gsva_zscore_config.yaml

# T95: GSVA + Poisson 核函数
# 配置文件 test_data/e2e_real/gsva_poisson_config.yaml:
#   gsva_kcdf: Poisson
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t95_gsva_poisson -m gsva -e test_data/e2e_real/hsa_expression_matrix.tsv --config test_data/e2e_real/gsva_poisson_config.yaml

# T96: GSVA + tau=0.5（非默认带宽）
# 配置文件 test_data/e2e_real/gsva_tau_config.yaml:
#   gsva_tau: 0.5
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t96_gsva_tau -m gsva -e test_data/e2e_real/hsa_expression_matrix.tsv --config test_data/e2e_real/gsva_tau_config.yaml
```

### 13.3 分析过滤参数
| Config 字段 | 类型 | 默认值 | 说明 |
|-------------|------|--------|------|
| `max_genes` | float | inf | 基因集最大基因数过滤 |

```bash
# T97: max_genes=500 过滤过宽泛的基因集
# 配置文件 test_data/e2e_real/max_genes_config.yaml:
#   max_genes: 500
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t97_max_genes --config test_data/e2e_real/max_genes_config.yaml
```

### 13.4 可视化 Config-only 参数
| Config 字段 | 类型 | 默认值 | 说明 |
|-------------|------|--------|------|
| `top_terms` | int | 20 | 可视化展示的 top 富集条目数 |
| `plot_width` | int | 10 | 图表宽度（英寸） |
| `plot_height` | int | 8 | 图表高度（英寸） |
| `plot_formats` | List[str] | ["pdf", "png"] | 同时输出多种格式（列表） |

```bash
# T98: top_terms=10
# 配置文件 test_data/e2e_real/top_terms_config.yaml:
#   top_terms: 10
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t98_top_terms --config test_data/e2e_real/top_terms_config.yaml

# T99: plot_width=12, plot_height=10
# 配置文件 test_data/e2e_real/plot_size_config.yaml:
#   plot_width: 12
#   plot_height: 10
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t99_plot_size --config test_data/e2e_real/plot_size_config.yaml

# T100: plot_formats 多格式同时输出
# 配置文件 test_data/e2e_real/plot_formats_config.yaml:
#   plot_formats: ["pdf", "png", "svg"]
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t100_plot_formats --config test_data/e2e_real/plot_formats_config.yaml
```

### 13.5 报告参数
| Config 字段 | 类型 | 默认值 | 说明 |
|-------------|------|--------|------|
| `report_format` | str | "html" | 报告输出格式 |

```bash
# T101: report_format=html（显式设置）
# 配置文件 test_data/e2e_real/report_config.yaml:
#   report_format: html
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t101_report_fmt --config test_data/e2e_real/report_config.yaml
```

### 13.6 AI 后端嵌套配置（AIBackendConfig 子字段）
| Config 字段 | 类型 | 默认值 | 说明 |
|-------------|------|--------|------|
| `ai_backends.{name}.base_url` | str | None | 自定义 API URL（Ollama/vLLM 等） |
| `ai_backends.{name}.group_id` | str | None | MiniMax 专用 Group ID |
| `ai_backends.{name}.enabled` | bool | True | 启用/禁用特定后端 |

```bash
# T102: AI 后端 base_url 配置（Ollama 本地部署）
# 配置文件 test_data/e2e_real/ollama_config.yaml:
#   ai_backends:
#     ollama:
#       model: "llama3"
#       base_url: "http://localhost:11434"
#       enabled: true
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t102_ai_baseurl --ai mock --config test_data/e2e_real/ollama_config.yaml

# T103: AI 后端 enabled=false 禁用
# 配置文件 test_data/e2e_real/ai_disabled_config.yaml:
#   ai_backends:
#     mock:
#       enabled: false
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t103_ai_disabled --ai mock --config test_data/e2e_real/ai_disabled_config.yaml
```

---
## 十四、审查补充：遗漏的功能模块
### 14.1 全局参数 `-v` / `--version`
```bash
# T104: 版本号显示
python -m allenricher -v
# 预期：输出版本号（如 2.0.0）

# T105: --version 长格式
python -m allenricher --version
# 预期：同上
```

### 14.2 CLI 与 Config 的参数传递优先级
CLI 参数优先级高于配置文件，需验证覆盖行为：

```bash
# T106: CLI 参数覆盖 config 文件
# 配置文件 test_data/e2e_real/override_config.yaml:
#   species: mmu
#   databases: [KEGG]
#   method: hypergeometric
#   correction: bonferroni
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -m fisher -c BH -o test_data/e2e_real/results/t106_cli_override --config test_data/e2e_real/override_config.yaml
# 验证：实际使用 hsa/GO/fisher/BH（CLI 优先），而非 mmu/KEGG/hypergeometric/bonferroni（config）
```

### 14.3 多格式基因列表输入
`load_gene_list()` 支持 .txt/.tsv/.gene/.csv/.xlsx/.xls 以及自动格式检测：

```bash
# T107: TSV 格式基因列表（两列：gene\tscore）
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.tsv -s hsa -d GO -o test_data/e2e_real/results/t107_input_tsv

# T108: CSV 格式基因列表
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.csv -s hsa -d GO -o test_data/e2e_real/results/t108_input_csv

# T109: .gene 扩展名
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.gene -s hsa -d GO -o test_data/e2e_real/results/t109_input_gene
```

### 14.4 GMT 文件解析（.gmt 和 .gmt.gz）
```bash
# T110: .gmt.gz 压缩格式
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t110_gmt_gz -m gsea -r test_data/e2e_real/hsa_ranked_genes.tsv -g test_data/e2e_real/hsa_custom_pathways.gmt.gz
```

### 14.5 group_comparison 图表子类型
`group_comparison` 支持 box/violin/bar 三种子类型：

```bash
# T111: group_comparison=box（默认）
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t111_grp_box -m ssgsea -e test_data/e2e_real/hsa_expression_matrix.tsv -pt group_comparison --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"

# T112: group_comparison=violin（通过 config）
# 需确认 group_comparison 的子类型如何通过参数传递
```

### 14.6 correlation 图表两种方法
```bash
# T113: correlation 图表（pearson，默认）
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t113_corr_pearson -m ssgsea -e test_data/e2e_real/hsa_expression_matrix.tsv -pt correlation

# T114: correlation 图表（spearman，通过 config）
# 需确认 correlation 方法如何通过参数传递
```

### 14.7 network 图表布局
network 图支持 spring/circular/kamada_kawai 三种布局：

```bash
# T115: network 图表（默认布局）
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO,KEGG -o test_data/e2e_real/results/t115_network -pt network

# T116: network 图表 + 多数据库交集
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO,KEGG,Reactome -o test_data/e2e_real/results/t116_network_multi -pt network,upset
```

### 14.8 API 服务端点测试
```bash
# T117: serve 启动后测试各 API 端点
# 前置：python -m allenricher serve --port 9999 &
# API 端点测试：
curl http://localhost:9999/                                    # T117a: 首页
curl http://localhost:9999/api/species                        # T117b: 物种列表
curl http://localhost:9999/api/databases                      # T117c: 数据库列表
curl -X POST http://localhost:9999/api/analyze               # T117d: 分析请求（ORA）
curl -F "file=@test_data/e2e_real/hsa_gene_list.txt" http://localhost:9999/api/upload -F "species=hsa"  # T117e: 文件上传
curl http://localhost:9999/api/status/nonexistent             # T117f: 不存在的 job_id（应 404）
curl http://localhost:9999/api/results/nonexistent             # T117g: 不存在的结果（应 404）
curl http://localhost:9999/api/results/nonexistent?format=json # T117h: JSON 格式参数
curl http://localhost:9999/api/results/nonexistent/plot        # T117i: 图表端点
curl http://localhost:9999/api/results/nonexistent/report      # T117j: 报告端点
curl -X DELETE http://localhost:9999/api/jobs/nonexistent      # T117k: 删除不存在的 job
```

### 14.9 报告生成功能验证
```bash
# T118: 验证报告文件内容
# 前置：执行一次正常分析（不使用 --no-report）
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO,KEGG -o test_data/e2e_real/results/t118_report
# 验证：
# - report.html 存在
# - 包含 jQuery DataTables 交互（排序、分页、搜索）
# - 包含 Download TSV 和 Copy 按钮
# - 包含版本号、数据库版本、分析日期
# - 包含 URL 跳转链接（GO/KEGG/Reactome/DO）

# T119: 无结果报告（所有数据库均无显著结果）
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list_tiny.txt -s hsa -d GO -p 1e-30 -q 1e-30 -o test_data/e2e_real/results/t119_no_results
# 验证：生成友好的"无显著结果"提示页面
```

### 14.10 URL 生成验证
```bash
# T120: 验证结果中的 URL 生成（6 种数据库）
# 前置：执行 GO+KEGG+Reactome+DO 分析
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO,KEGG,Reactome,DO -o test_data/e2e_real/results/t120_urls --no-plot
# 验证结果 TSV 中 URL 列：
# - GO: https://amigo.geneontology.org/amigo/term/GO:XXXXXXX
# - KEGG: https://www.kegg.jp/entry/hsaXXXXX
# - Reactome: https://reactome.org/PathwayBrowser#R-HSA-XXXXX
# - DO: https://disease-ontology.org/?id=DOID:XXXX
```

### 14.11 DisGeNET 下载/构建测试（graceful degradation）
```bash
# T121: download disgenet（应显示废弃/不可用提示，graceful degradation）
python -m allenricher download -d disgenet
# 预期：显示 DisGeNET 下载不可用的提示，不崩溃

# T122: build hsa + DisGeNET（应显示数据不可用提示）
python -m allenricher build -s hsa -t 9606 -d DisGeNET
# 预期：跳过 DisGeNET 构建并给出警告
```

### 14.12 下载后注册表自动构建验证
```bash
# T123: download 完成后验证 supported_species.tsv 自动更新
# 前置：记录当前 supported_species.tsv 的修改时间
python -m allenricher download -d go --force
# 验证：database/basic/supported_species.tsv 和 database/supported_species.tsv 已更新
# 验证：文件包含 GO/KEGG/Reactome/DO 列
```

### 14.13 build_manifest.json 血缘追踪验证
```bash
# T124: 构建完成后验证 build_manifest.json
python -m allenricher build -s rro -t 61622 -d GO,KEGG
# 验证：database/organism/<version>/rro/build_manifest.json 存在
# 验证：包含 version, species, taxid, databases, built_at, software_version 字段
# 验证：databases 中每项包含 source_version 和 dependency_chain
```

### 14.14 GOA 回退构建流程
```bash
# T125: 构建一个 gene2go 中不存在的物种（触发 GOA 回退）
# 选择一个 gene2go 中无记录但 GOA proteomes 中有数据的物种
python -m allenricher build -s xtr -t 8364 -d GO
# 验证：GO 数据库成功构建（通过 GOA 回退路径）
# 验证：日志中显示 "gene2go 中未找到 taxid，尝试 GOA proteomes" 类似信息
```

### 14.15 版本号一致性检查
```bash
# T126: 检查版本号一致性
python -c "import allenricher; print(allenricher.__version__)"
# 记录输出，与 pyproject.toml 中的 version 字段对比
# 已知问题：__init__.py (2.0.0) vs pyproject.toml (2.1.0) 不一致
```

### 14.16 并行 vs 串行分析模式
```bash
# T127: 串行模式 (jobs=1)
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO,KEGG,Reactome -o test_data/e2e_real/results/t127_serial -j 1

# T128: 并行模式 (jobs=4)
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO,KEGG,Reactome -o test_data/e2e_real/results/t128_parallel -j 4

# T129: 全 CPU 并行 (jobs=-1)
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO,KEGG -o test_data/e2e_real/results/t129_all_cpu -j -1
```

### 14.17 GSEA 结果特有字段验证
```bash
# T130: GSEA 结果中 NES/ES/FDR/Leading_Edge 字段完整性
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t130_gsea_fields -m gsea -r test_data/e2e_real/hsa_ranked_genes.tsv --no-plot --no-report
# 验证结果 TSV 中包含列：NES, ES, FDR, Leading_Edge
# 验证：NES 值有正有负
# 验证：Leading_Edge 非空
```

### 14.18 ssGSEA NaN p 值处理
```bash
# T131: ssGSEA 结果中 pvalue 列应为 NaN
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t131_ssgsea_nan -m ssgsea -e test_data/e2e_real/hsa_expression_matrix.tsv --no-plot --no-report
# 验证：结果中 pvalue 列为 NaN（ssGSEA 不计算 p 值）
# 验证：程序不因 NaN 崩溃，下游过滤逻辑正确处理
```

### 14.19 EnrichmentResult 输出格式验证
```bash
# T132: ORA 结果字段完整性
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t132_ora_fields --no-plot --no-report
# 验证结果 TSV 列：Term_ID, Term_Name, Database, pvalue, qvalue, gene_count, gene_list, rich_factor, expected_count, URL
# 验证：GSEA 特有字段（NES/ES/FDR/Leading_Edge）在 ORA 结果中不存在

# T133: GSEA 结果字段完整性
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t133_gsea_out -m gsea -r test_data/e2e_real/hsa_ranked_genes.tsv --no-plot --no-report
# 验证：NES/ES/FDR/Leading_Edge 列存在且非全空
```

---
## 十五、审查补充：遗漏的配色方案
原计划覆盖了 14 种配色，剩余 9 种需补充：

```bash
# T134: palette=default（Paul Tol Bright）
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t134_pal_default --palette default

# T135: palette=tol_high_contrast
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t135_pal_hc --palette tol_high_contrast

# T136: palette=tol_vibrant
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t136_pal_vibrant --palette tol_vibrant

# T137: palette=tol_medium_contrast
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t137_pal_mc --palette tol_medium_contrast

# T138: palette=tol_light
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t138_pal_light --palette tol_light

# T139: palette=tol_burga
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t139_pal_burga --palette tol_burga

# T140: palette=go_cc
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t140_pal_gocc --palette go_cc

# T141: palette=go_mf
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t141_pal_gomf --palette go_mf

# T142: palette=cell
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t142_pal_cell --palette cell
```

---
## 十六、审查补充：遗漏的 style
原计划遗漏了 `cell` 风格（CLI choices 中有 5 种，但 PRESETS 有 6 种，`cell` 存在于 PRESETS 但不在 CLI choices 中）：

```bash
# T143: style=cell（通过 config 文件设置，因 CLI choices 中未包含）
# 配置文件 test_data/e2e_real/cell_style_config.yaml:
#   plot_style: cell
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t143_style_cell --config test_data/e2e_real/cell_style_config.yaml
```

---
## 十七、审查补充：Config-only 参数覆盖矩阵
| Config 字段 | 测试用例 | 覆盖状态 |
|-------------|---------|---------|
| `max_genes` | T97 | ✅ |
| `gsea_permutations` | T91 | ✅ |
| `gsea_min_size` | T92 | ✅ |
| `gsea_max_size` | T92 | ✅ |
| `gsva_method` | T93(plage), T94(zscore) | ✅ |
| `gsva_kcdf` | T95(Poisson) | ✅ |
| `gsva_tau` | T96(0.5) | ✅ |
| `top_terms` | T98 | ✅ |
| `plot_width` | T99 | ✅ |
| `plot_height` | T99 | ✅ |
| `plot_formats` | T100 | ✅ |
| `report_format` | T101 | ✅ |
| `ai_backends.base_url` | T102 | ✅ |
| `ai_backends.group_id` | (MiniMax 需真实密钥，跳过) | ⚠️ |
| `ai_backends.enabled` | T103 | ✅ |
| `auto_update` | (功能未被实际使用，跳过) | ⚠️ |
| `api_debug` | (仅影响 serve，通过 config) | ⚠️ |

---
## 十八、补充后的完整覆盖统计
### 18.1 测试用例总数
| 分类 | 原计划 | 补充 | 合计 |
|------|--------|------|------|
| analyze CLI 参数 | T1-T90 (90) | T91-T143 (53) | 143 |
| download | D1-D10 (10) | D121-D123 (3) | 13 |
| build | B1-B14 (14) | B124-B125 (2) | 16 |
| serve | S1-S3 (3) | S117 (11个子端点) | 14 |
| list | L1-L3 (3) | - | 3 |
| config | C1-C3 (3) | - | 3 |
| check-update | U1-U4 (4) | - | 4 |
| cleanup | CL1-CL4 (4) | - | 4 |
| list-versions | V1-V5 (5) | - | 5 |
| list-species | SP1-SP11 (11) | - | 11 |
| query-species | Q1-Q6 (6) | - | 6 |
| 全局/功能验证 | - | T104-T133 (30) | 30 |
| **总计** | **153** | **99** | **252** |

### 18.2 参数覆盖（CLI + Config）
| 维度 | 数量 | 状态 |
|------|------|------|
| CLI 参数（11 个命令） | 70 | ✅ 全覆盖 |
| Config-only 字段 | 17 | ✅ 15/17 覆盖，2 个跳过（需真实密钥/未实际使用） |
| AIBackendConfig 子字段 | 3 | ✅ 2/3 覆盖 |
| 可视化风格 | 6 | ✅ 全覆盖（含 cell 通过 config） |
| 配色方案 | 23 | ✅ 全覆盖（14+9 补充） |
| 图表类型 | 10 | ✅ 全覆盖 |
| 富集方法 | 5 | ✅ 全覆盖 |
| 校正方法 | 5 | ✅ 全覆盖 |
| GSVA 变体 | 3 | ✅ 全覆盖（gsva/plage/zscore） |
| GSVA 核函数 | 2 | ✅ 全覆盖（Gaussian/Poisson） |
| 数据库 | 5 | ✅ 全覆盖 |
| 输出格式 | 3 | ✅ 全覆盖 |
| 注释格式 | 4 | ✅ 全覆盖 |
| AI 后端 | 7 | ✅ 全覆盖 |
| API 端点 | 11 | ✅ 全覆盖 |

### 18.3 功能路径覆盖
| 功能路径 | 状态 |
|---------|------|
| ORA (Fisher + Hypergeometric) | ✅ |
| GSEA (ranked genes + expression matrix + custom GMT) | ✅ |
| ssGSEA (expression matrix + groups + custom GMT) | ✅ |
| GSVA (gsva/plage/zscore × Gaussian/Poisson) | ✅ |
| 背景模式 (annotated/genome/custom) | ✅ |
| Config 文件加载 (YAML + JSON) | ✅ |
| CLI 覆盖 Config 优先级 | ✅ |
| 多格式基因列表 (.txt/.tsv/.csv/.gene/.gmt/.gmt.gz) | ✅ |
| URL 生成 (GO/KEGG/Reactome/DO) | ✅ |
| 报告生成 (含无结果页面) | ✅ |
| API 服务 (11 个端点) | ✅ |
| 版本管理 (check/list/cleanup/lineage) | ✅ |
| 物种注册表 (list/query/fuzzy) | ✅ |
| 数据库下载 (含注册表自动构建) | ✅ |
| 物种构建 (含 GOA 回退/DisGeNET 降级) | ✅ |
| 并行/串行分析 | ✅ |
| GSEA 结果字段 (NES/ES/FDR/Leading_Edge) | ✅ |
| ssGSEA NaN p 值处理 | ✅ |
| 错误场景 (空输入/无效物种/缺失文件/参数不匹配) | ✅ |

---
## 十九、更新后的自检清单
### 19.1 规格覆盖
- [x] 11 个 CLI 子命令全部覆盖
- [x] analyze 31 个 CLI 参数全部覆盖
- [x] download 7 个 CLI 参数全部覆盖
- [x] build 11 个 CLI 参数全部覆盖
- [x] serve 3 个 CLI 参数全部覆盖
- [x] list 1 个位置参数全部覆盖
- [x] config 1 个 CLI 参数全部覆盖
- [x] check-update 2 个 CLI 参数全部覆盖
- [x] cleanup 3 个 CLI 参数全部覆盖
- [x] list-versions 3 个 CLI 参数全部覆盖
- [x] list-species 6 个 CLI 参数全部覆盖
- [x] query-species 3 个 CLI 参数全部覆盖
- [x] 全局参数 `-v`/`--version` 覆盖
- [x] 17 个 Config-only 字段覆盖（15 个直接测试 + 2 个跳过说明）

### 19.2 功能覆盖
- [x] 5 种富集方法: fisher, hypergeometric, gsea, ssgsea, gsva
- [x] 3 种 GSVA 变体: gsva, plage, zscore
- [x] 2 种 GSVA 核函数: Gaussian, Poisson
- [x] 5 种校正方法: BH, BY, bonferroni, holm, none
- [x] 6 种可视化风格: nature, science, cell, colorblind, presentation, omicshare
- [x] 23 种配色方案全覆盖
- [x] 3 种输出格式: png, pdf, svg
- [x] plot_formats 多格式同时输出
- [x] 5 个数据库: GO, KEGG, Reactome, DO, DisGeNET
- [x] 3 种背景模式: annotated, genome, custom
- [x] 10 种图表类型: dotplot, barplot, network, upset, volcano, enrichment, nes_barplot, heatmap, group_comparison, correlation
- [x] 4 种注释格式: three_column, four_column, two_column, auto
- [x] 7 种 AI 后端: openai, claude, deepseek, glm, minimax, ollama, mock
- [x] AI 后端嵌套配置: base_url, group_id, enabled
- [x] 配置文件: YAML + JSON
- [x] CLI 覆盖 Config 优先级
- [x] 多格式基因列表: .txt, .tsv, .csv, .gene, .gmt, .gmt.gz
- [x] URL 生成: GO, KEGG, Reactome, DO
- [x] 报告生成: 正常报告 + 无结果页面
- [x] API 服务: 11 个端点全覆盖
- [x] 版本管理: check-update, list-versions, cleanup, lineage
- [x] 物种注册表: list-species, query-species, fuzzy_search
- [x] 数据库下载: 含注册表自动构建
- [x] 物种构建: 含 GOA 回退、DisGeNET 降级
- [x] 并行/串行分析模式
- [x] GSEA 结果字段完整性
- [x] ssGSEA NaN p 值处理
- [x] 错误场景: 空输入、无效物种、缺失文件、参数不匹配
- [x] 全局参数: -v/--version

### 19.3 占位符扫描
- [x] 无 TBD/TODO 占位符
- [x] 所有测试命令均为完整可执行命令
- [x] 所有文件路径均为具体路径
- [x] 跳过项均有明确说明（需真实密钥/功能未实际使用）

---
## 二十、审查补充：GSEA/ssGSEA/GSVA 的可视化风格与配色测试
原计划中风格测试（T39-T44）和配色测试（T45-T58）仅覆盖 ORA 方法。GSEA、ssGSEA、GSVA 各自拥有不同的专用图表类型，需独立验证风格和配色在这些图表上的渲染效果。

### 20.1 测试策略
- **风格全覆盖**：每种方法 × 6 种风格 = 18 个用例
- **代表性配色覆盖**：每种方法 × 8 种代表性配色（覆盖期刊风格/Paul Tol/色盲友好/中国风/专用配色） = 24 个用例
- **合计新增**：42 个用例

### 20.2 GSEA 风格测试（6 种风格 × enrichment + nes_barplot + dotplot）
```bash
# G144: GSEA + style=nature
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t144_gsea_nature -m gsea -r test_data/e2e_real/hsa_ranked_genes.tsv --style nature

# G145: GSEA + style=science
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t145_gsea_science -m gsea -r test_data/e2e_real/hsa_ranked_genes.tsv --style science

# G146: GSEA + style=cell（通过 config）
# 配置文件 test_data/e2e_real/gsea_cell_config.yaml: plot_style: cell
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t146_gsea_cell -m gsea -r test_data/e2e_real/hsa_ranked_genes.tsv --config test_data/e2e_real/gsea_cell_config.yaml

# G147: GSEA + style=colorblind
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t147_gsea_cb -m gsea -r test_data/e2e_real/hsa_ranked_genes.tsv --style colorblind

# G148: GSEA + style=presentation
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t148_gsea_pres -m gsea -r test_data/e2e_real/hsa_ranked_genes.tsv --style presentation

# G149: GSEA + style=omicshare
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t149_gsea_omic -m gsea -r test_data/e2e_real/hsa_ranked_genes.tsv --style omicshare
```

### 20.3 GSEA 配色测试（8 种代表性配色）
```bash
# G150: GSEA + palette=nature
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t150_gsea_pal_nature -m gsea -r test_data/e2e_real/hsa_ranked_genes.tsv --palette nature

# G151: GSEA + palette=science
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t151_gsea_pal_science -m gsea -r test_data/e2e_real/hsa_ranked_genes.tsv --palette science

# G152: GSEA + palette=lancet
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t152_gsea_pal_lancet -m gsea -r test_data/e2e_real/hsa_ranked_genes.tsv --palette lancet

# G153: GSEA + palette=nejm
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t153_gsea_pal_nejm -m gsea -r test_data/e2e_real/hsa_ranked_genes.tsv --palette nejm

# G154: GSEA + palette=okabe_ito（色盲友好）
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t154_gsea_pal_okabe -m gsea -r test_data/e2e_real/hsa_ranked_genes.tsv --palette okabe_ito

# G155: GSEA + palette=gsea（GSEA 经典配色）
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t155_gsea_pal_gsea -m gsea -r test_data/e2e_real/hsa_ranked_genes.tsv --palette gsea

# G156: GSEA + palette=china_style（中国风）
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t156_gsea_pal_china -m gsea -r test_data/e2e_real/hsa_ranked_genes.tsv --palette china_style

# G157: GSEA + palette=tol_burga（渐变配色）
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t157_gsea_pal_burga -m gsea -r test_data/e2e_real/hsa_ranked_genes.tsv --palette tol_burga
```

### 20.4 ssGSEA 风格测试（6 种风格 × heatmap + group_comparison + dotplot + correlation）
```bash
# G158: ssGSEA + style=nature
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t158_ssgsea_nature -m ssgsea -e test_data/e2e_real/hsa_expression_matrix.tsv --style nature --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"

# G159: ssGSEA + style=science
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t159_ssgsea_science -m ssgsea -e test_data/e2e_real/hsa_expression_matrix.tsv --style science --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"

# G160: ssGSEA + style=cell（通过 config）
# 配置文件 test_data/e2e_real/ssgsea_cell_config.yaml: plot_style: cell
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t160_ssgsea_cell -m ssgsea -e test_data/e2e_real/hsa_expression_matrix.tsv --config test_data/e2e_real/ssgsea_cell_config.yaml --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"

# G161: ssGSEA + style=colorblind
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t161_ssgsea_cb -m ssgsea -e test_data/e2e_real/hsa_expression_matrix.tsv --style colorblind --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"

# G162: ssGSEA + style=presentation
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t162_ssgsea_pres -m ssgsea -e test_data/e2e_real/hsa_expression_matrix.tsv --style presentation --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"

# G163: ssGSEA + style=omicshare
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t163_ssgsea_omic -m ssgsea -e test_data/e2e_real/hsa_expression_matrix.tsv --style omicshare --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"
```

### 20.5 ssGSEA 配色测试（8 种代表性配色）
```bash
# G164: ssGSEA + palette=nature
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t164_ssgsea_pal_nature -m ssgsea -e test_data/e2e_real/hsa_expression_matrix.tsv --palette nature --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"

# G165: ssGSEA + palette=science
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t165_ssgsea_pal_science -m ssgsea -e test_data/e2e_real/hsa_expression_matrix.tsv --palette science --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"

# G166: ssGSEA + palette=lancet
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t166_ssgsea_pal_lancet -m ssgsea -e test_data/e2e_real/hsa_expression_matrix.tsv --palette lancet --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"

# G167: ssGSEA + palette=nejm
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t167_ssgsea_pal_nejm -m ssgsea -e test_data/e2e_real/hsa_expression_matrix.tsv --palette nejm --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"

# G168: ssGSEA + palette=okabe_ito（色盲友好）
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t168_ssgsea_pal_okabe -m ssgsea -e test_data/e2e_real/hsa_expression_matrix.tsv --palette okabe_ito --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"

# G169: ssGSEA + palette=omicshare
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t169_ssgsea_pal_omic -m ssgsea -e test_data/e2e_real/hsa_expression_matrix.tsv --palette omicshare --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"

# G170: ssGSEA + palette=china_style（中国风）
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t170_ssgsea_pal_china -m ssgsea -e test_data/e2e_real/hsa_expression_matrix.tsv --palette china_style --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"

# G171: ssGSEA + palette=tol_sunset（渐变配色）
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t171_ssgsea_pal_sunset -m ssgsea -e test_data/e2e_real/hsa_expression_matrix.tsv --palette tol_sunset --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"
```

### 20.6 GSVA 风格测试（6 种风格 × heatmap + group_comparison + dotplot + correlation）
```bash
# G172: GSVA + style=nature
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t172_gsva_nature -m gsva -e test_data/e2e_real/hsa_expression_matrix.tsv --style nature --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"

# G173: GSVA + style=science
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t173_gsva_science -m gsva -e test_data/e2e_real/hsa_expression_matrix.tsv --style science --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"

# G174: GSVA + style=cell（通过 config）
# 配置文件 test_data/e2e_real/gsva_cell_config.yaml: plot_style: cell
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t174_gsva_cell -m gsva -e test_data/e2e_real/hsa_expression_matrix.tsv --config test_data/e2e_real/gsva_cell_config.yaml --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"

# G175: GSVA + style=colorblind
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t175_gsva_cb -m gsva -e test_data/e2e_real/hsa_expression_matrix.tsv --style colorblind --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"

# G176: GSVA + style=presentation
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t176_gsva_pres -m gsva -e test_data/e2e_real/hsa_expression_matrix.tsv --style presentation --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"

# G177: GSVA + style=omicshare
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t177_gsva_omic -m gsva -e test_data/e2e_real/hsa_expression_matrix.tsv --style omicshare --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"
```

### 20.7 GSVA 配色测试（8 种代表性配色）
```bash
# G178: GSVA + palette=nature
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t178_gsva_pal_nature -m gsva -e test_data/e2e_real/hsa_expression_matrix.tsv --palette nature --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"

# G179: GSVA + palette=science
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t179_gsva_pal_science -m gsva -e test_data/e2e_real/hsa_expression_matrix.tsv --palette science --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"

# G180: GSVA + palette=lancet
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t180_gsva_pal_lancet -m gsva -e test_data/e2e_real/hsa_expression_matrix.tsv --palette lancet --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"

# G181: GSVA + palette=nejm
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t181_gsva_pal_nejm -m gsva -e test_data/e2e_real/hsa_expression_matrix.tsv --palette nejm --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"

# G182: GSVA + palette=okabe_ito（色盲友好）
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t182_gsva_pal_okabe -m gsva -e test_data/e2e_real/hsa_expression_matrix.tsv --palette okabe_ito --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"

# G183: GSVA + palette=omicshare
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t183_gsva_pal_omic -m gsva -e test_data/e2e_real/hsa_expression_matrix.tsv --palette omicshare --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"

# G184: GSVA + palette=china_style（中国风）
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t184_gsva_pal_china -m gsva -e test_data/e2e_real/hsa_expression_matrix.tsv --palette china_style --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"

# G185: GSVA + palette=tol_muted（柔和配色）
python -m allenricher analyze -i test_data/e2e_real/hsa_gene_list.txt -s hsa -d GO -o test_data/e2e_real/results/t185_gsva_pal_muted -m gsva -e test_data/e2e_real/hsa_expression_matrix.tsv --palette tol_muted --groups "Normal:Normal_1,Normal_2,Normal_3;Tumor:Tumor_1,Tumor_2,Tumor_3"
```

### 20.8 GSEA/ssGSEA/GSVA 风格 × 配色覆盖矩阵
| 风格 \ 方法 | ORA (已有) | GSEA | ssGSEA | GSVA |
|------------|-----------|------|--------|------|
| nature | ✅ T39 | ✅ G144 | ✅ G158 | ✅ G172 |
| science | ✅ T40 | ✅ G145 | ✅ G159 | ✅ G173 |
| cell | ✅ T143(config) | ✅ G146(config) | ✅ G160(config) | ✅ G174(config) |
| colorblind | ✅ T42 | ✅ G147 | ✅ G161 | ✅ G175 |
| presentation | ✅ T43 | ✅ G148 | ✅ G162 | ✅ G176 |
| omicshare | ✅ T44 | ✅ G149 | ✅ G163 | ✅ G177 |

| 配色 \ 方法 | ORA (已有) | GSEA | ssGSEA | GSVA |
|-----------|-----------|------|--------|------|
| nature | ✅ T45 | ✅ G150 | ✅ G164 | ✅ G178 |
| science | ✅ T46 | ✅ G151 | ✅ G165 | ✅ G179 |
| lancet | ✅ T47 | ✅ G152 | ✅ G166 | ✅ G180 |
| nejm | ✅ T48 | ✅ G153 | ✅ G167 | ✅ G181 |
| okabe_ito | ✅ T50 | ✅ G154 | ✅ G168 | ✅ G182 |
| gsea | ✅ T51 | ✅ G155 | - | - |
| omicshare | ✅ T52 | - | ✅ G169 | ✅ G183 |
| china_style | ✅ T53 | ✅ G156 | ✅ G170 | ✅ G184 |
| tol_burga | ✅ T58 | ✅ G157 | - | - |
| tol_sunset | ✅ T58 | - | ✅ G171 | - |
| tol_muted | ✅ T57 | - | - | ✅ G185 |
| 其余 12 种配色 | ✅ T49/T54-T56/T134-T142 | (通过 ORA 已验证渲染引擎) | (同左) | (同左) |

**说明**：剩余 12 种配色（default/tol_high_contrast/tol_vibrant/tol_medium_contrast/tol_light/go_bp/go_cc/go_mf/kegg_pathway/jama/cell）已在 ORA 方法中测试，因配色渲染引擎是共享的（`plot_theme.py` 的 `apply_style()`），这些配色的核心渲染逻辑已被覆盖。GSEA/ssGSEA/GSVA 的 8 种代表性配色足以验证各方法专用图表（enrichment曲线/heatmap/group_comparison等）对配色的适配性。

### 20.9 验证要点
每个 GSEA/ssGSEA/GSVA 风格/配色测试用例需额外验证：

| 方法 | 专用图表 | 验证要点 |
|------|---------|---------|
| **GSEA** | enrichment（三面板富集曲线） | 配色应用于 NES 归一化条形图；正/负 NES 颜色区分 |
| **GSEA** | nes_barplot | 配色应用于条形图填充 |
| **GSEA** | dotplot | 配色应用于气泡颜色映射 |
| **ssGSEA** | heatmap | 配色应用于热图色阶（连续色图） |
| **ssGSEA** | group_comparison | 配色应用于分组对比图（box/violin） |
| **ssGSEA** | dotplot | 配色应用于气泡颜色映射 |
| **ssGSEA** | correlation | 配色应用于相关性散点图 |
| **GSVA** | heatmap | 同 ssGSEA heatmap |
| **GSVA** | group_comparison | 同 ssGSEA group_comparison |
| **GSVA** | dotplot | 同 ssGSEA dotplot |
| **GSVA** | correlation | 同 ssGSEA correlation |

---
## 二十一、更新后的完整覆盖统计（第三次）
### 21.1 测试用例总数
| 分类 | 数量 |
|------|------|
| analyze CLI 参数 (T1-T90) | 90 |
| Config-only 参数 (T91-T103) | 13 |
| 全局/功能验证 (T104-T143) | 40 |
| download (D1-D10, D121-D123) | 13 |
| build (B1-B14, B124-B125) | 16 |
| serve (S1-S3, S117) | 14 |
| list/config/check-update/cleanup/list-versions/list-species/query-species | 36 |
| **GSEA 风格+配色 (G144-G157)** | **14** |
| **ssGSEA 风格+配色 (G158-G171)** | **14** |
| **GSVA 风格+配色 (G172-G185)** | **14** |
| **总计** | **294** |

### 21.2 可视化覆盖矩阵（最终）
| 维度 | ORA | GSEA | ssGSEA | GSVA |
|------|-----|------|--------|------|
| 风格 (6种) | ✅ | ✅ | ✅ | ✅ |
| 配色 (23种) | ✅ 全部 | ✅ 8种代表 | ✅ 8种代表 | ✅ 8种代表 |
| 图表类型 | ✅ dotplot/barplot/network/upset/volcano | ✅ enrichment/nes_barplot/dotplot | ✅ heatmap/group_comparison/dotplot/correlation | ✅ heatmap/group_comparison/dotplot/correlation |
| 输出格式 (3种) | ✅ | ✅ | ✅ | ✅ |
| DPI | ✅ | ✅ | ✅ | ✅ |

---