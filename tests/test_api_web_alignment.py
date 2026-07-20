"""Contract tests for the API, CLI bridge, and local Web workbench."""

import subprocess
from pathlib import Path
from unittest.mock import patch

from fastapi.testclient import TestClient

from allenricher.api import server


def test_json_api_accepts_all_analysis_methods(monkeypatch, tmp_path):
    monkeypatch.setenv("ALLENRICHER_API_JOB_DIR", str(tmp_path / "jobs"))
    server.jobs.clear()
    client = TestClient(server.app)

    requests = [
        {"genes": ["TP53"], "databases": ["GO"], "method": "hypergeometric"},
        {
            "ranked_genes": [{"gene": "TP53", "weight": 2.0}, {"gene": "BRCA1", "weight": -1.0}],
            "databases": ["GO"],
            "method": "gsea",
        },
        {
            "expression_matrix": {"TP53": {"S1": 2.0, "S2": 3.0}, "BRCA1": {"S1": 1.0, "S2": 4.0}},
            "databases": ["GO"],
            "method": "ssgsea",
        },
        {
            "expression_matrix": {"TP53": {"S1": 2.0, "S2": 3.0}, "BRCA1": {"S1": 1.0, "S2": 4.0}},
            "databases": ["GO"],
            "method": "gsva",
        },
    ]
    with patch("allenricher.api.server.run_analysis"):
        for payload in requests:
            response = client.post("/api/analyze", json=payload)
            assert response.status_code == 200, response.text


def test_upload_api_accepts_method_specific_files(monkeypatch, tmp_path):
    monkeypatch.setenv("ALLENRICHER_API_JOB_DIR", str(tmp_path / "jobs"))
    server.jobs.clear()
    client = TestClient(server.app)
    cases = [
        (
            "gsea",
            {"ranked_file": ("ranked.tsv", b"gene\tweight\nTP53\t2\nBRCA1\t-1\n", "text/tab-separated-values")},
            {},
        ),
        (
            "ssgsea",
            {"expression_file": ("expression.tsv", b"gene\tS1\tS2\nTP53\t2\t3\nBRCA1\t1\t4\n", "text/tab-separated-values")},
            {"groups": "Control:S1;Treatment:S2"},
        ),
        (
            "gsva",
            {"expression_file": ("expression.tsv", b"gene\tS1\tS2\nTP53\t2\t3\nBRCA1\t1\t4\n", "text/tab-separated-values")},
            {"groups": "Control:S1;Treatment:S2"},
        ),
    ]
    with patch("allenricher.api.server.run_analysis"):
        for method, files, extra in cases:
            response = client.post("/api/upload", files=files, data={"method": method, "databases": "GO", **extra})
            assert response.status_code == 200, response.text
            request = server.jobs[response.json()["job_id"]]["request"]
            assert request["method"] == method
            assert request["databases"] == ["GO"]
            if extra:
                assert request["groups"] == extra["groups"]


def test_upload_api_accepts_ai_interpretation_options(monkeypatch, tmp_path):
    monkeypatch.setenv("ALLENRICHER_API_JOB_DIR", str(tmp_path / "jobs"))
    server.jobs.clear()
    client = TestClient(server.app)
    files = {"gene_file": ("genes.txt", b"TP53\nBRCA1\n", "text/plain")}
    data = {
        "method": "hypergeometric",
        "databases": "GO",
        "ai_backend": "mock",
        "ai_mode": "caption",
        "ai_top_n": "9",
    }
    with patch("allenricher.api.server.run_analysis"):
        response = client.post("/api/upload", files=files, data=data)
    assert response.status_code == 200, response.text
    request = server.jobs[response.json()["job_id"]]["request"]
    assert request["ai_backend"] == "mock"
    assert request["ai_mode"] == "caption"
    assert request["ai_top_n"] == 9


def test_method_specific_input_validation(monkeypatch, tmp_path):
    monkeypatch.setenv("ALLENRICHER_API_JOB_DIR", str(tmp_path / "jobs"))
    server.jobs.clear()
    client = TestClient(server.app)
    with patch("allenricher.api.server.run_analysis"):
        assert client.post("/api/analyze", json={"method": "hypergeometric"}).status_code == 422
        assert client.post("/api/analyze", json={"genes": ["TP53"], "method": "gsea"}).status_code == 422
        assert client.post("/api/analyze", json={"genes": ["TP53"], "method": "gsva"}).status_code == 422


def test_api_requests_r_plots_by_default(tmp_path):
    request = server.EnrichmentRequest(genes=["A"], ranked_genes=[server.RankedGene(gene="A", weight=1.0)], databases=["GO"], method="gsea")
    files = {"input": str(tmp_path / "genes.txt"), "ranked": str(tmp_path / "ranked.tsv")}
    command = server.build_cli_command(request, files, tmp_path / "output")
    assert "--use-r-plots" in command


def test_cli_command_contains_method_inputs_and_plot_options(tmp_path):
    request = server.EnrichmentRequest(
        genes=["A"],
        databases=["GO", "KEGG"],
        method="gsea",
        plot_types="enrichment,lollipop",
        use_r_plots=True,
        groups="Control:S1;Treatment:S2",
        categorical_palette="okabe_ito",
        ai_backend="mock",
        ai_mode="reviewer",
        ai_top_n=7,
    )
    files = {
        "input": str(tmp_path / "genes.txt"),
        "ranked": str(tmp_path / "ranked.tsv"),
        "expression": str(tmp_path / "expression.tsv"),
    }
    command = server.build_cli_command(request, files, tmp_path / "output")
    assert command[:4] == [server.sys.executable, "-m", "allenricher", "analyze"]
    assert command[command.index("-d") + 1] == "GO,KEGG"
    assert command[command.index("--ranked-genes") + 1] == files["ranked"]
    assert command[command.index("--expression-matrix") + 1] == files["expression"]
    assert command[command.index("--plot-types") + 1] == "enrichment,lollipop"
    assert command[command.index("--methods-language") + 1] == "en"
    assert "--use-r-plots" in command
    assert command[command.index("--ai") + 1] == "mock"
    assert command[command.index("--ai-mode") + 1] == "reviewer"
    assert command[command.index("--ai-top-n") + 1] == "7"


def test_run_analysis_collects_real_cli_artifacts(monkeypatch, tmp_path):
    monkeypatch.setenv("ALLENRICHER_API_JOB_DIR", str(tmp_path / "jobs"))
    server.jobs.clear()
    request = server.EnrichmentRequest(genes=["TP53", "BRCA1"], databases=["GO"])
    job_id, job = server._create_job(request)
    server._prepare_json_inputs(job, request)

    def fake_run(command, **kwargs):
        output_dir = Path(command[command.index("-o") + 1])
        (output_dir / "plots").mkdir(parents=True, exist_ok=True)
        (output_dir / "GO_enrichment.tsv").write_text(
            "Database\tTerm_ID\tTerm_Name\tGene_Count\tP_Value\tAdjusted_P_Value\n"
            "GO\tGO:0001\tCell cycle\t2\t0.001\t0.01\n",
            encoding="utf-8",
        )
        (output_dir / "plots" / "GO_barplot.png").write_bytes(b"PNG")
        (output_dir / "report.html").write_text("<html>AllEnricher</html>", encoding="utf-8")
        (output_dir / "ai_interpretation.json").write_text(
            '{"schema_version": 1, "profile": "summary", "databases": {}, "evidence": {}}',
            encoding="utf-8",
        )
        return subprocess.CompletedProcess(command, 0, stdout="completed\n", stderr="")

    monkeypatch.setattr(server.subprocess, "run", fake_run)
    server.run_analysis(job_id, request)
    completed = server.jobs[job_id]
    assert completed["status"] == "completed"
    assert completed["results"]["GO"][0]["Term_Name"] == "Cell cycle"
    assert Path(completed["results_file"]).is_file()
    assert Path(completed["report_file"]).is_file()
    assert completed["ai_interpretation"]["profile"] == "summary"
    assert {item["category"] for item in completed["artifacts"]} >= {"table", "figure", "report", "log"}

    server.jobs.clear()
    restored = server._get_job(job_id)
    assert restored is not None
    assert restored["status"] == "completed"

    client = TestClient(server.app)
    response = client.get(f"/api/results/{job_id}/ai-interpretation")
    assert response.status_code == 200
    assert response.json()["schema_version"] == 1


def test_ai_interpretation_failure_keeps_analysis_completed(monkeypatch, tmp_path):
    monkeypatch.setenv("ALLENRICHER_API_JOB_DIR", str(tmp_path / "jobs"))
    server.jobs.clear()
    request = server.EnrichmentRequest(genes=["TP53"], databases=["GO"], ai_backend="mock")
    job_id, job = server._create_job(request)
    server._prepare_json_inputs(job, request)

    def fake_run(command, **kwargs):
        output_dir = Path(command[command.index("-o") + 1])
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / "GO_enrichment.tsv").write_text(
            "Database\tTerm_ID\tTerm_Name\tGene_Count\tP_Value\tAdjusted_P_Value\n"
            "GO\tGO:0001\tCell cycle\t1\t0.001\t0.01\n",
            encoding="utf-8",
        )
        (output_dir / "ai_interpretation_error.json").write_text(
            '{"error_code":"AI_INTERPRETATION_FAILED","backend":"deepseek",'
            '"mode":"summary","message":"Error: Missing credentials"}\n',
            encoding="utf-8",
        )
        return subprocess.CompletedProcess(command, 0, stdout="completed\n", stderr="")

    monkeypatch.setattr(server.subprocess, "run", fake_run)
    server.run_analysis(job_id, request)
    completed = server.jobs[job_id]
    assert completed["status"] == "completed"
    assert completed["error"] is None
    assert completed["ai_interpretation_error"]["error_code"] == "AI_INTERPRETATION_FAILED"

    client = TestClient(server.app)
    status = client.get(f"/api/status/{job_id}").json()
    assert status["status"] == "completed"
    assert status["ai_interpretation_error"]["message"] == "Error: Missing credentials"


def test_persisted_ai_error_is_loaded_for_existing_job(monkeypatch, tmp_path):
    monkeypatch.setenv("ALLENRICHER_API_JOB_DIR", str(tmp_path / "jobs"))
    server.jobs.clear()
    request = server.EnrichmentRequest(genes=["TP53"], databases=["GO"])
    job_id, job = server._create_job(request)
    error_file = Path(job["output_dir"]) / "ai_interpretation_error.json"
    error_file.write_text('{"error_code":"AI_INTERPRETATION_FAILED"}\n', encoding="utf-8")
    server.jobs.clear()
    restored = server._get_job(job_id)
    assert restored["ai_interpretation_error"]["error_code"] == "AI_INTERPRETATION_FAILED"


def test_unconfigured_ai_backend_is_rejected_before_job_creation(monkeypatch, tmp_path):
    monkeypatch.setenv("ALLENRICHER_API_JOB_DIR", str(tmp_path / "jobs"))
    monkeypatch.delenv("ALLENRICHER_CONFIG", raising=False)
    monkeypatch.delenv("DEEPSEEK_API_KEY", raising=False)
    server.jobs.clear()
    client = TestClient(server.app)
    response = client.post(
        "/api/analyze",
        json={"genes": ["TP53"], "databases": ["GO"], "ai_backend": "deepseek"},
    )
    assert response.status_code == 422
    assert "not configured" in response.json()["detail"]
    assert not server.jobs


def test_artifact_endpoint_blocks_input_files(monkeypatch, tmp_path):
    monkeypatch.setenv("ALLENRICHER_API_JOB_DIR", str(tmp_path / "jobs"))
    server.jobs.clear()
    request = server.EnrichmentRequest(genes=["TP53"], databases=["GO"])
    job_id, job = server._create_job(request)
    server._prepare_json_inputs(job, request)
    client = TestClient(server.app)
    response = client.get(f"/api/results/{job_id}/files/input/genes.txt")
    assert response.status_code == 400


def test_default_job_root_is_persistent_user_storage():
    assert server.DEFAULT_JOB_ROOT == Path.home() / ".allenricher" / "api_jobs"


def test_web_omitted_runtime_fields_use_safe_backend_defaults():
    request = server.EnrichmentRequest(genes=["TP53"], databases=["GO"])
    assert request.jobs == 1
    assert request.database_dir is None
    assert request.use_version is None
    assert request.verbose is False


def test_species_database_support_reflects_installed_files(monkeypatch, tmp_path):
    database_dir = tmp_path / "database"
    species_dir = database_dir / "organism" / "v1" / "hsa"
    species_dir.mkdir(parents=True)
    for name in (
        "hsa.GO.gmt",
        "hsa.ChEA3_2gene.tab.gz",
        "hsa.AnimalTFDB_2gene.tab.gz",
        "hsa.CUSTOM.gmt",
    ):
        (species_dir / name).touch()
    monkeypatch.setenv("ALLENRICHER_DATABASE_DIR", str(database_dir))
    client = TestClient(server.app)
    expected = ["GO", "ChEA3", "CUSTOM"]
    assert client.get("/api/species/hsa/databases").json()["databases"] == expected
    species = {item["code"]: item for item in client.get("/api/species").json()}
    assert species["hsa"]["databases"] == expected
    assert client.get("/api/species/bad%20code/databases").status_code == 400


def test_web_workbench_uses_dynamic_safe_four_method_ui():
    html = (Path(server.__file__).parent / "static" / "index.html").read_text(encoding="utf-8")
    for method in ("hypergeometric", "gsea", "ssgsea", "gsva"):
        assert f'value="{method}"' in html
    assert "fisher" not in html.lower()
    assert 'fetch("/api/species")' in html
    assert 'fetch("/api/databases")' in html
    assert 'id="methodsLanguage"' in html
    assert 'form.append("methods_language", "en")' in html
    assert "/methods-reference" in html
    assert "innerHTML" not in html
    assert "--accent: #2c5282" in html
    assert 'font-family: Georgia, "Times New Roman", serif' in html
    assert html.count('value="0.05" min="0.000001" max="1" step="any"') == 3
    assert ".output { min-width: 0; width: 100%; }" in html
    assert ".workspace { grid-template-columns: minmax(0, 1fr); }" in html
    assert ".empty-state[hidden] { display: none; }" in html
    assert "[hidden] { display: none !important; }" in html
    assert 'id="tsvLink" class="button" href="#" hidden' in html
    assert 'id="reportLink" class="button" href="#" target="_blank" rel="noopener" hidden' in html
    assert 'artifact.category === "report"' in html


def test_web_workbench_footer_links_repository_and_shows_citation():
    html = (Path(server.__file__).parent / "static" / "index.html").read_text(encoding="utf-8")

    assert 'class="site-footer"' in html
    assert 'href="https://github.com/zhangducsu/AllEnricher-v2"' in html
    assert '>https://github.com/zhangducsu/AllEnricher-v2</a>' in html
    assert "Please cite:" in html
    assert "BMC Bioinformatics. 2020;21:106." in html


def test_web_workbench_uses_progressive_method_specific_workflow():
    html = (Path(server.__file__).parent / "static" / "index.html").read_text(encoding="utf-8")
    for heading in ("Choose an analysis method", "Provide input data", "Define the analysis scope"):
        assert heading in html
    assert 'name="oraMode" value="paste" checked' in html
    assert 'name="oraMode" value="file"' in html
    assert 'name="gseaMode" value="ranked" checked' in html
    assert 'name="gseaMode" value="differential"' in html
    assert 'id="rankedFile"' in html
    assert 'id="differentialFile"' in html
    assert 'id="gseaExpressionFile"' not in html
    assert 'id="expressionFile"' in html
    assert 'id="groups"' in html
    assert 'id="advancedDialog" class="settings-drawer"' in html
    assert 'byId("significanceSettings").hidden = method === "ssgsea" || method === "gsva"' in html
    assert 'if (method === "ssgsea" || method === "gsva") {' in html
    assert 'appendValue(form, "groups", "groups")' in html
    assert 'plotTypes.includes("enrichment")) {' in html
    assert 'appendFile(form, "gmt_file", "gmtFile")' in html


def test_web_explains_each_method_and_reports_registry_and_installed_support():
    html = (Path(server.__file__).parent / "static" / "index.html").read_text(encoding="utf-8")
    assert 'id="methodSummary" class="method-summary" aria-live="polite"' in html
    for method in ("ORA", "GSEA", "ssGSEA", "GSVA"):
        assert f'["{method}",' in html
    assert 'state.species.filter((item) => item.databases?.length)' in html
    assert 'fetch("/api/species/summary")' in html
    assert 'facts.className = "method-facts"' in html
    assert '["Registered species", registry ? registry.total_species.toLocaleString("en-US") : "Unavailable"]' in html
    assert '["Coverage by database", registryCoverage || "Unavailable"]' in html
    assert '["Installed locally", state.species.length ?' in html
    assert 'value.count.toLocaleString("zh-CN")' not in html
    assert "function databaseDisplayName(name)" in html
    assert '`${databaseDisplayName(key)} ${value.count.toLocaleString("en-US")} species`' in html
    assert 'supported.flatMap((item) => item.databases)' in html
    assert "If the required database is unavailable for your species, contact the administrator to build it." in html


def test_web_workbench_has_single_tf_database_entry_and_structured_plot_controls():
    html = (Path(server.__file__).parent / "static" / "index.html").read_text(encoding="utf-8")
    for database in ("TRRUST", "ChEA3", "AnimalTFDB", "hTFtarget"):
        assert database in html
    for legacy_control in ('id="tfDatabase"', 'id="tfOnly"', 'id="tfLibrary"', 'id="tfTissue"'):
        assert legacy_control not in html
    assert 'const TF_DATABASES = new Set(' in html
    assert 'databases.some((name) => TF_DATABASES.has(name))' in html
    assert 'name = "plotType"' in html
    assert 'form.append("plot_types", plotTypes.join(","))' in html
    assert 'id="plotTypes"' not in html
    assert 'form.append("no_plot", !generatePlots)' in html
    assert 'form.append("no_report", !byId("generateReport").checked)' in html
    assert 'return Boolean(selected?.databases?.includes(database.name))' in html
    assert 'byId("species").addEventListener("input", refreshSpeciesDatabases)' in html


def test_web_inputs_explain_the_actual_file_formats():
    html = (Path(server.__file__).parent / "static" / "index.html").read_text(encoding="utf-8")
    for hint in (
        "Enter one gene ID per line",
        "TSV/CSV with <code>gene</code> and numeric <code>weight</code> columns",
        "a gene column, and at least one numeric ranking statistic",
        "gene IDs in the first column, sample names in the header",
        "Samples are read from the matrix header",
        "GMT file with a gene-set ID, description, and tab-delimited gene list",
    ):
        assert hint in html
    for field in ("geneText", "groups"):
        assert f'id="{field}"' in html
        assert f'aria-describedby="{field}Hint"' in html
    for field in ("geneFile", "rankedFile", "differentialFile", "expressionFile", "backgroundFile", "gmtFile"):
        assert f'id="{field}" class="native-file-input" type="file"' in html
        assert f'aria-describedby="{field}Hint {field}Name"' in html
        assert f'data-file-input="{field}">Choose file</button>' in html
        assert f'id="{field}Name" class="file-name" aria-live="polite">No file selected</span>' in html
    assert "function initializeFilePickers()" in html
    assert 'id="expressionFileHint" class="field-hint">TSV/CSV with ' in html
    assert 'id="expressionFileHint" class="field-hint">UTF-8' not in html


def test_web_uses_clear_labels_and_species_display_without_changing_species_code():
    html = (Path(server.__file__).parent / "static" / "index.html").read_text(encoding="utf-8")

    assert "<title>AllEnricher Gene Function Enrichment Workbench</title>" in html
    assert "<h1>AllEnricher <small>v{{ ALLENRICHER_VERSION }}</small></h1>" in html
    assert "<span>Gene Function Enrichment Workbench</span>" in html
    assert '["Custom databases", ["CUSTOM"]' in html
    assert 'value="Human (Homo sapiens) - 9606"' in html
    assert "function speciesLabel(item)" in html
    assert "function selectedSpeciesCode()" in html
    assert 'form.append("species", species)' in html
    assert "label.textContent = databaseDisplayName(database.name)" in html


def test_web_converts_differential_results_to_a_standard_ranked_file():
    html = (Path(server.__file__).parent / "static" / "index.html").read_text(encoding="utf-8")
    for element in ("differentialGeneColumn", "differentialMetricColumn", "differentialDirection"):
        assert f'id="{element}"' in html
    assert "async function buildRankedFileFromDifferential()" in html
    assert 'form.append("ranked_file", rankedFile, "ranked_from_differential.tsv")' in html
    assert 'const form = await buildRequestForm();' in html
    assert 'byId("differentialFile").addEventListener("change", readDifferentialColumns)' in html
    assert 'const metricNames = ["stat", "waldstatistic", "t", "tvalue", "score", "signedscore", "log2foldchange", "logfc"]' in html


def test_web_builds_sample_groups_from_the_expression_header():
    html = (Path(server.__file__).parent / "static" / "index.html").read_text(encoding="utf-8")
    for element in ("groupEditor", "sampleGroupName", "assignSampleGroup", "clearSampleGroup", "sampleGroupList"):
        assert f'id="{element}"' in html
    assert 'byId("expressionFile").addEventListener("change", readExpressionSamples)' in html
    assert 'function setExpressionSamples(samples)' in html
    assert 'state.sampleGroups = new Map(samples.map((sample, index) => [sample, reliable ? inferred[index] : ""]))' in html
    assert 'throw new Error("Group-comparison figures require a group for every sample")' in html


def test_web_plot_settings_are_context_aware_and_hide_runtime_internals():
    html = (Path(server.__file__).parent / "static" / "index.html").read_text(encoding="utf-8")
    assert '<option value="nature">Nature</option>' in html
    assert '<option value="science">Science</option>' in html
    assert '<option value="presentation">Presentation</option>' in html
    assert '<option value="cell">' not in html
    assert '<option value="omicshare">' not in html
    for control in ("databaseDir", "useVersion", "jobs", "verbose", "useRPlots"):
        assert f'id="{control}"' not in html
    assert "Operating Environment" not in html
    assert 'form.append("use_r_plots", true)' not in html
    assert 'const PLOT_PALETTE_ROLES = {' in html
    assert 'class="palette-preview"' not in html
    for role in ("categorical", "sequential", "diverging"):
        assert f'class="palette-picker" data-palette-role="{role}"' in html
        assert f'id="{role}PaletteTrigger"' in html
        assert f'id="{role}PaletteMenu" class="palette-menu"' in html
    assert "function initializePalettePickers()" in html
    assert "palette-option-preview" in html
    assert 'colorbrewer_blues: ["#9ECAE1", "#08519C"]' in html
    assert 'element.style.background = `linear-gradient(to right, ${colors.join(", ")})`' in html
    assert 'byId(`${role}PaletteField`).hidden = !roles.has(role)' in html
    assert 'byId("plotOutputSettings").hidden = !enabled' in html
    assert '"Select at least one figure type"' in html


def test_web_species_selection_refreshes_installed_database_support():
    html = (Path(server.__file__).parent / "static" / "index.html").read_text(encoding="utf-8")
    assert "async function refreshSpeciesDatabases()" in html
    assert 'fetch(`/api/species/${encodeURIComponent(species)}/databases`)' in html
    assert "No installed databases are available for the selected species" in html
    assert "request !== state.speciesSupportRequest" in html
    assert 'if (database.name === "CUSTOM") return true' not in html
    assert "input.disabled = !supported" in html
    assert 'wrapper.classList.toggle("unavailable", !supported)' in html


def test_web_gsea_settings_follow_selected_plot_types():
    html = (Path(server.__file__).parent / "static" / "index.html").read_text(encoding="utf-8")
    for setting in ("gseaEnrichmentSettings", "gseaMultiSettings", "gseaEmapplotSettings"):
        assert f'id="{setting}"' in html
    assert '!selected.has("enrichment")' in html
    assert '!selected.has("enrichment2")' in html
    assert '!selected.has("emapplot")' in html
    assert 'if (plotTypes.includes("emapplot")) {' in html
    assert 'if (plotTypes.includes("enrichment2")) {' in html


def test_web_workbench_uses_auditable_result_views_and_figure_browser():
    html = (Path(server.__file__).parent / "static" / "index.html").read_text(encoding="utf-8")
    for view in ("overview", "table", "figures", "files"):
        assert f'data-view="{view}"' in html
        assert f'data-result-view="{view}"' in html
        assert f'aria-controls="{view}View"' in html
        assert f'aria-labelledby="{view}ViewTab"' in html
    assert 'id="resultSearch"' in html
    assert 'id="figureList"' in html
    assert 'id="figurePreview"' in html
    assert 'link.textContent = "Open PDF Figures"' in html
    assert 'const labels = {table: "Results", report: "Reports", figure: "Figures", log: "Logs", other: "Other"}' in html
    assert 'window.confirm("Delete this analysis and all of its output files?")' in html


def test_web_workbench_submits_and_renders_structured_ai_interpretation():
    html = (Path(server.__file__).parent / "static" / "index.html").read_text(encoding="utf-8")
    for element in ("generateAi", "aiBackend", "aiMode", "aiTopN", "aiInterpretation"):
        assert f'id="{element}"' in html
    assert 'form.append("ai_backend", byId("aiBackend").value)' in html
    assert 'appendValue(form, "ai_top_n", "aiTopN")' in html
    assert "/ai-interpretation" in html
    assert "function renderAiInterpretation()" in html
    assert 'Model: ${backend} | Profile: ${interpretation.profile || "summary"}' in html
    assert "function showEvidence(evidenceId)" in html
    assert "evidenceElementId(evidenceId)" in html
    assert "aiInterpretationError" in html
    assert 'fetch("/api/ai/backends")' in html
    assert "option.disabled = !backend.configured" in html


def test_static_web_assets_are_packaged():
    pyproject = Path(server.__file__).resolve().parents[2] / "pyproject.toml"
    assert '"api/static/*"' in pyproject.read_text(encoding="utf-8")
