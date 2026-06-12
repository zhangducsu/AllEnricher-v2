# R 脚本出图系统 + Python 调用层 — 测试报告

**测试时间**: 2026-06-04 19:22:58

**测试统计**: 43/43 通过, 0 失败, 0 警告

## 通过项 (43)

- **R 脚本目录**: f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2\allenricher\visualization\r_scripts 存在
- **R 脚本文件**: common.R
- **R 脚本文件**: gsea_dotplot.R
- **R 脚本文件**: gsea_barplot.R
- **R 脚本文件**: gsea_nes_plot.R
- **R 脚本文件**: gsea_ridgeplot.R
- **R 脚本文件**: gsea_heatmap.R
- **R 脚本文件**: gsea_circos.R
- **R 脚本文件**: gsea_emapplot.R
- **R 脚本文件**: gsea_cnetplot.R
- **R 脚本文件**: gsea_enrichment_plot.R
- **R 脚本文件**: gsea_enrichment_plot2.R
- **Python 文件**: r_plotter.py
- **py_compile**: r_plotter.py 语法正确
- **模块导入**: allenricher.visualization.r_plotter 导入成功
- **R_PLOT_TYPES**: 包含 10 种图表类型
- **R_PLOT_FUNC_MAP**: 包含 10 个映射
- **便捷函数**: plot_gsea_dotplot_r
- **便捷函数**: plot_gsea_barplot_r
- **便捷函数**: plot_gsea_nes_plot_r
- **便捷函数**: plot_gsea_ridgeplot_r
- **便捷函数**: plot_gsea_heatmap_r
- **便捷函数**: plot_gsea_emapplot_r
- **便捷函数**: plot_gsea_cnetplot_r
- **便捷函数**: plot_gsea_circos_r
- **便捷函数**: plot_gsea_enrichment_r
- **便捷函数**: plot_gsea_enrichment2_r
- **核心函数**: check_r_environment
- **核心函数**: run_r_script
- **R 环境**: Rscript 可用
- **R 语法**: common.R
- **R 语法**: gsea_dotplot.R
- **R 语法**: gsea_barplot.R
- **R 语法**: gsea_nes_plot.R
- **R 语法**: gsea_ridgeplot.R
- **R 语法**: gsea_heatmap.R
- **R 语法**: gsea_circos.R
- **R 语法**: gsea_emapplot.R
- **R 语法**: gsea_cnetplot.R
- **R 语法**: gsea_enrichment_plot.R
- **R 语法**: gsea_enrichment_plot2.R
- **common.R 加载**: source 成功，函数可用
- **错误处理**: 不存在的脚本返回 False

## 失败项

无失败项。

## 创建的文件清单

### R 脚本 (`allenricher/visualization/r_scripts/`)

- `common.R`
- `gsea_dotplot.R`
- `gsea_barplot.R`
- `gsea_nes_plot.R`
- `gsea_ridgeplot.R`
- `gsea_heatmap.R`
- `gsea_circos.R`
- `gsea_emapplot.R`
- `gsea_cnetplot.R`
- `gsea_enrichment_plot.R`
- `gsea_enrichment_plot2.R`

### Python 文件 (`allenricher/visualization/`)

- `r_plotter.py`

## 结论

所有测试通过，R 脚本出图系统和 Python 调用层创建成功。
