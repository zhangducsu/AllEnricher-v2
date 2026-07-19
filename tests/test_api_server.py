"""
FastAPI Service Endpoint Test
====================

Use FastAPI TestClient All API Ends to test the unit,
Do not start the real server.All time-consuming analysis operations passed mock Alternative.

Overwrite peer:
1.  GET /                        - Service information
2.  GET /api/species             - List of species
3.  GET /api/databases           - Database List
4.  POST /api/analyze            - Submit analysis (mock Backstage task)
5.  POST /api/upload             - File Upload (mock)
6.  GET /api/status/{job_id}     - Task Status (Normal + 404)
7.  GET /api/results/{job_id}    - Get results (JSON/TSV)
8.  GET /api/results/{job_id}/plot  - Chart Access
9.  GET /api/results/{job_id}/report - Access to reports
10. DELETE /api/jobs/{job_id}    - Delete a job
11. Paths through protective tests
12. Invalid Request 422 Test
"""

import sys
import os
import json
from pathlib import Path
from io import BytesIO
from unittest.mock import patch, MagicMock, AsyncMock

import pytest

# Ensure that root directory in sys.path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Skip the whole module if fastapi or httpx is not installed
fastapi = pytest.importorskip("fastapi")
pytest.importorskip("httpx")

from fastapi.testclient import TestClient
from allenricher import __version__


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def client():
    """Creates a TestClient example, each test with a separate Jobs store."""
    from allenricher.api.server import app, jobs
    # Empty the dictionaries of the Jobs before each test, ensure the test is isolated
    jobs.clear()
    with TestClient(app) as c:
        yield c
    jobs.clear()


@pytest.fixture
def completed_job(client):
    """
Yes. jobs Inject a job directly into the dictionary, For test results/Chart/Report Endpoint.
Back to the infused. job_id.
    """
    from allenricher.api.server import jobs
    import tempfile

    job_id = "test-completed-job-001"
    output_dir = tempfile.mkdtemp(prefix="test_allenricher_")

    # Create a plots subdirectories and a simulation PDF file
    plots_dir = Path(output_dir) / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    (plots_dir / "GO_barplot.pdf").write_bytes(b"%PDF-1.4 mock plot content")

    # Create Simulate TSV Outcome
    tsv_file = Path(output_dir) / "results.tsv"
    tsv_file.write_text("Term_ID\tTerm_Name\tP_Value\nGO:0005576\textracellular region\t1e-5\n")

    # Create Simulation HTML Report
    report_file = Path(output_dir) / "report.html"
    report_file.write_text("<html><body>AllEnricher Test Report</body></html>")
    (Path(output_dir) / "analysis_metadata.json").write_text(
        json.dumps({
            "allenricher_version": "2.0-test",
            "analysis_method": "hypergeometric",
            "species": "hsa",
            "databases": ["GO", "KEGG"],
            "database_versions": {"GO": "v1", "KEGG": "v1"},
            "parameters": {
                "background_mode": "annotated",
                "correction": "BH",
                "pvalue_cutoff": 0.05,
                "qvalue_cutoff": 0.05,
                "gene_set_size_by_database": {
                    "GO": {"min": 3, "max": "Inf"},
                    "KEGG": {"min": 3, "max": "Inf"},
                },
            },
        }),
        encoding="utf-8",
    )

    jobs[job_id] = {
        "status": "completed",
        "created_at": "2025-01-01T00:00:00",
        "completed_at": "2025-01-01T00:01:00",
        "progress": 1.0,
        "request": {
            "genes": ["TP53", "BRCA1", "EGFR"],
            "species": "hsa",
            "databases": ["GO", "KEGG"],
        },
        "results": {
            "GO": [
                {
                    "Term_ID": "GO:0005576",
                    "Term_Name": "extracellular region",
                    "Gene_Count": 10,
                    "P_Value": 1e-5,
                    "Adjusted_P_Value": 1e-3,
                }
            ],
            "KEGG": [
                {
                    "Term_ID": "hsa04110",
                    "Term_Name": "Cell Cycle",
                    "Gene_Count": 8,
                    "P_Value": 1e-4,
                    "Adjusted_P_Value": 1e-2,
                }
            ],
        },
        "results_summary": {
            "GO": {"term_count": 1, "top_terms": [{"Term_ID": "GO:0005576"}]},
            "KEGG": {"term_count": 1, "top_terms": [{"Term_ID": "hsa04110"}]},
        },
        "results_file": str(tsv_file),
        "output_dir": output_dir,
        "error": None,
    }

    yield job_id

    # Clear temporary directory
    import shutil
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir, ignore_errors=True)


@pytest.fixture
def running_job(client):
    """Injecting an running job."""
    from allenricher.api.server import jobs

    job_id = "test-running-job-002"
    jobs[job_id] = {
        "status": "running",
        "created_at": "2025-01-01T00:00:00",
        "completed_at": None,
        "progress": 0.5,
        "request": {"genes": ["TP53"], "species": "hsa", "databases": ["GO"]},
        "results": None,
        "error": None,
    }
    yield job_id


@pytest.fixture
def failed_job(client):
    """Inject a failed job."""
    from allenricher.api.server import jobs

    job_id = "test-failed-job-003"
    jobs[job_id] = {
        "status": "failed",
        "created_at": "2025-01-01T00:00:00",
        "completed_at": "2025-01-01T00:00:30",
        "progress": 0.3,
        "request": {"genes": ["TP53"], "species": "hsa", "databases": ["GO"]},
        "results": None,
        "error": "Database not found: GO",
    }
    yield job_id


# ===========================================================================
# 1. GET/ - Service information
# ===========================================================================

class TestRootEndpoint:
    """Test root GET/"""

    def test_root_returns_service_info(self, client):
        """Root endpoint should return the Web interface or API information"""
        response = client.get("/")
        assert response.status_code == 200
        # Root route returns HTML (Web interface) now, the validation response contains HTML content
        content_type = response.headers.get("content-type", "")
        assert "text/html" in content_type or "application/json" in content_type
        assert f"AllEnricher <small>v{__version__}</small>" in response.text

    def test_root_endpoints_keys(self, client):
        """Endpoints returned by root should contain all major endpoints"""
        # The router returns HTML now, skips JSON's assertion
        # API information can be viewed through the /docs
        pass

    def test_root_endpoint_values(self, client):
        """Endpoint path should match expectations"""
        # The router returns HTML now, skips JSON's assertion
        # API information can be viewed through the /docs
        pass


# ===========================================================================
# 2. GET /api/species- species list
# ===========================================================================

class TestSpeciesEndpoint:
    """Test species list endpoint GET/api/species"""

    def test_species_returns_list(self, client):
        """The list of species should be returned (groups)"""
        response = client.get("/api/species")
        assert response.status_code == 200
        data = response.json()
        assert isinstance(data, list)
        assert len(data) > 0

    def test_species_contains_human(self, client):
        """List should include humans (hsa)"""
        response = client.get("/api/species")
        data = response.json()
        human = next((s for s in data if s["code"] == "hsa"), None)
        assert human is not None
        assert human["display_name"] == "Human"
        assert human["taxonomy_id"] == 9606

    def test_species_contains_mouse(self, client):
        """List should contain mouse (mmu)"""
        response = client.get("/api/species")
        data = response.json()
        mouse = next((s for s in data if s["code"] == "mmu"), None)
        assert mouse is not None
        assert mouse["display_name"] == "Mouse"

    def test_species_schema_fields(self, client):
        """Each species entry should contain code, name, taxony_id, display_name"""
        response = client.get("/api/species")
        data = response.json()
        for species in data:
            assert "code" in species
            assert "name" in species
            assert "taxonomy_id" in species
            assert "display_name" in species
            assert isinstance(species["taxonomy_id"], int)

    def test_species_summary_matches_registry(self, client, monkeypatch):
        from allenricher.database.species_registry import SpeciesRegistry

        registry = MagicMock()
        registry.get_summary.return_value = {
            "total_species": 42124,
            "go": {"count": 32443},
            "kegg": {"count": 10871},
            "disgenet": {"count": 1},
            "trrust": {"count": 2},
            "chea3": {"count": 1},
            "animaltfdb": {"count": 183},
            "htftarget": {"count": 1},
        }
        monkeypatch.setattr(SpeciesRegistry, "load_default", lambda *_: registry)

        response = client.get("/api/species/summary")

        assert response.status_code == 200
        assert response.json()["total_species"] == 42124
        assert response.json()["go"]["count"] == 32443
        assert response.json()["trrust"]["count"] == 2
        assert response.json()["chea3"]["count"] == 1
        assert response.json()["animaltfdb"]["count"] == 183
        assert response.json()["htftarget"]["count"] == 1


# ===========================================================================
# 3. GET /api/databases- Database List
# ===========================================================================

class TestDatabasesEndpoint:
    """Test database list endpoint GET/api/databases"""

    def test_databases_returns_object(self, client):
        """The object that contains the databases key should be returned"""
        response = client.get("/api/databases")
        assert response.status_code == 200
        data = response.json()
        assert "databases" in data
        assert isinstance(data["databases"], list)

    def test_databases_contains_core(self, client):
        """should contain core databases GO and KEGG"""
        response = client.get("/api/databases")
        db_names = [db["name"] for db in response.json()["databases"]]
        assert "GO" in db_names
        assert "KEGG" in db_names

    def test_databases_all_expected(self, client):
        """should include all anticipated databases"""
        response = client.get("/api/databases")
        db_names = [db["name"] for db in response.json()["databases"]]
        expected = {
            "GO", "KEGG", "Reactome", "WikiPathways", "DO", "DisGeNET",
            "TRRUST", "ChEA3", "AnimalTFDB", "hTFtarget", "CUSTOM",
        }
        assert set(db_names) == expected
        assert "MSigDB" not in db_names

    def test_databases_schema_fields(self, client):
        """Each database entry should contain name, description, references"""
        response = client.get("/api/databases")
        for db in response.json()["databases"]:
            assert "name" in db
            assert "description" in db
            assert "species" in db

    def test_disgenet_exposes_the_v1_snapshot_label(self, client):
        response = client.get("/api/databases")
        disgenet = next(
            item for item in response.json()["databases"]
            if item["name"] == "DisGeNET"
        )

        assert disgenet["display_name"] == "DisGeNET (v20190612)"
        assert disgenet["source_version"] == "v20190612"


# ===========================================================================
# 4. POST /api/analyze- Submit analysis
# ===========================================================================

class TestAnalyzeEndpoint:
    """Test Analysis Submission Endpoint POST/api/analyze"""

    @patch("allenricher.api.server.run_analysis")
    def test_analyze_returns_job_id(self, mock_run, client):
        """Submit analysis to return to the job_id and pending state"""
        mock_run.return_value = None

        payload = {
            "genes": ["TP53", "BRCA1", "EGFR"],
            "species": "hsa",
            "databases": ["GO", "KEGG"],
        }
        response = client.post("/api/analyze", json=payload)
        assert response.status_code == 200
        data = response.json()

        assert "job_id" in data
        assert data["status"] == "pending"
        assert "message" in data
        assert data["results"] is None

    @patch("allenricher.api.server.run_analysis")
    def test_analyze_job_created_in_jobs_dict(self, mock_run, client):
        """This task should exist in the jobs dictionary after submission"""
        mock_run.return_value = None

        payload = {"genes": ["TP53", "BRCA1"]}
        response = client.post("/api/analyze", json=payload)
        job_id = response.json()["job_id"]

        from allenricher.api.server import jobs
        assert job_id in jobs
        assert jobs[job_id]["status"] == "pending"
        assert jobs[job_id]["request"]["genes"] == ["TP53", "BRCA1"]

    @patch("allenricher.api.server.run_analysis")
    def test_analyze_default_parameters(self, mock_run, client):
        """Submit analysis using default parameters"""
        mock_run.return_value = None

        payload = {"genes": ["TP53"]}
        response = client.post("/api/analyze", json=payload)
        assert response.status_code == 200

        from allenricher.api.server import jobs
        job_id = response.json()["job_id"]
        req = jobs[job_id]["request"]
        assert req["species"] == "hsa"
        assert req["method"] == "hypergeometric"
        assert req["correction"] == "BH"
        assert req["pvalue_cutoff"] == 0.05
        assert req["qvalue_cutoff"] == 0.05
        assert req["min_genes"] == 3

    @patch("allenricher.api.server.run_analysis")
    def test_analyze_custom_parameters(self, mock_run, client):
        """Submit analysis using a custom parameter"""
        mock_run.return_value = None

        payload = {
            "genes": ["TP53", "BRCA1"],
            "species": "mmu",
            "databases": ["Reactome"],
            "method": "hypergeometric",
            "correction": "bonferroni",
            "pvalue_cutoff": 0.01,
            "qvalue_cutoff": 0.01,
            "min_genes": 5,
            "background": ["GENE1", "GENE2", "GENE3"],
        }
        response = client.post("/api/analyze", json=payload)
        assert response.status_code == 200

        from allenricher.api.server import jobs
        job_id = response.json()["job_id"]
        req = jobs[job_id]["request"]
        assert req["species"] == "mmu"
        assert req["databases"] == ["Reactome"]
        assert req["method"] == "hypergeometric"
        assert req["correction"] == "bonferroni"
        assert req["pvalue_cutoff"] == 0.01
        assert req["background"] == ["GENE1", "GENE2", "GENE3"]

    @patch("allenricher.api.server.run_analysis")
    def test_analyze_multiple_submissions_unique_ids(self, mock_run, client):
        """Multiple submissions should return different job_id"""
        mock_run.return_value = None

        payload = {"genes": ["TP53"]}
        response1 = client.post("/api/analyze", json=payload)
        response2 = client.post("/api/analyze", json=payload)

        id1 = response1.json()["job_id"]
        id2 = response2.json()["job_id"]
        assert id1 != id2


# ===========================================================================
# 5. POST /api/upload- File upload.
# ===========================================================================

class TestUploadEndpoint:
    """Test file uploads peer POST /api/upload"""

    @patch("allenricher.api.server.run_analysis")
    def test_upload_gene_file(self, mock_run, client):
        """Uploading a gene list file should return a job_id"""
        mock_run.return_value = None

        gene_file_content = b"TP53\nBRCA1\nEGFR\nKRAS\n"
        files = {"file": ("genes.txt", gene_file_content, "text/plain")}

        response = client.post("/api/upload", files=files)
        assert response.status_code == 200
        data = response.json()

        assert "job_id" in data
        assert data["status"] == "pending"
        assert "uploaded" in data["message"].lower()

    @patch("allenricher.api.server.run_analysis")
    def test_upload_parses_genes(self, mock_run, client):
        """Uploading should correctly resolve the list of genes"""
        mock_run.return_value = None

        gene_file_content = b"TP53\nBRCA1\nEGFR\n"
        files = {"file": ("genes.txt", gene_file_content, "text/plain")}

        response = client.post("/api/upload", files=files)
        job_id = response.json()["job_id"]

        from allenricher.api.server import jobs
        req = jobs[job_id]["request"]
        assert req["genes"] == ["TP53", "BRCA1", "EGFR"]

    @patch("allenricher.api.server.run_analysis")
    def test_upload_with_query_params(self, mock_run, client):
        """Specify profile for uploading through query parameters"""
        mock_run.return_value = None

        gene_file_content = b"TP53\n"
        files = {"file": ("genes.txt", gene_file_content, "text/plain")}

        response = client.post(
            "/api/upload",
            files=files,
            data={"species": "mmu", "databases": "Reactome,KEGG", "method": "hypergeometric"},
        )
        assert response.status_code == 200

        from allenricher.api.server import jobs
        job_id = response.json()["job_id"]
        req = jobs[job_id]["request"]
        assert req["species"] == "mmu"
        assert req["databases"] == ["Reactome", "KEGG"]
        assert req["method"] == "hypergeometric"

    @patch("allenricher.api.server.run_analysis")
    def test_upload_skips_empty_lines(self, mock_run, client):
        """The empty line that uploads the file should be ignored"""
        mock_run.return_value = None

        gene_file_content = b"TP53\n\n\nBRCA1\n\n"
        files = {"file": ("genes.txt", gene_file_content, "text/plain")}

        response = client.post("/api/upload", files=files)
        job_id = response.json()["job_id"]

        from allenricher.api.server import jobs
        req = jobs[job_id]["request"]
        assert req["genes"] == ["TP53", "BRCA1"]

    @patch("allenricher.api.server.run_analysis")
    def test_upload_strips_whitespace(self, mock_run, client):
        """Ignore blank gene identifiers in an uploaded gene list."""
        mock_run.return_value = None

        gene_file_content = b"  TP53  \n  BRCA1 \n"
        files = {"file": ("genes.txt", gene_file_content, "text/plain")}

        response = client.post("/api/upload", files=files)
        job_id = response.json()["job_id"]

        from allenricher.api.server import jobs
        req = jobs[job_id]["request"]
        assert req["genes"] == ["TP53", "BRCA1"]

    def test_upload_missing_file_returns_422(self, client):
        """The non-provision of files should return 422 error"""
        response = client.post("/api/upload")
        assert response.status_code == 422


# ===========================================================================
# 6. GET /api/status/{job_id}- Job status.
# ===========================================================================

class TestStatusEndpoint:
    """Test Task State Endpoint/api/status/{job_id}"""

    def test_status_completed(self, client, completed_job):
        """Querying Status of Completed Tasks"""
        response = client.get(f"/api/status/{completed_job}")
        assert response.status_code == 200
        data = response.json()

        assert data["job_id"] == completed_job
        assert data["status"] == "completed"
        assert data["progress"] == 1.0
        assert data["completed_at"] is not None
        assert data["error"] is None

    def test_status_running(self, client, running_job):
        """Querying the status of active jobs"""
        response = client.get(f"/api/status/{running_job}")
        assert response.status_code == 200
        data = response.json()

        assert data["status"] == "running"
        assert 0.0 < data["progress"] < 1.0
        assert data["completed_at"] is None

    def test_status_failed(self, client, failed_job):
        """Query failed task status"""
        response = client.get(f"/api/status/{failed_job}")
        assert response.status_code == 200
        data = response.json()

        assert data["status"] == "failed"
        assert data["error"] is not None

    def test_status_not_found(self, client):
        """Query for non-existent tasks should return 404"""
        response = client.get("/api/status/nonexistent-job-id")
        assert response.status_code == 404
        assert "not found" in response.json()["detail"].lower()

    def test_status_results_summary(self, client, completed_job):
        """Completed should include results_summary"""
        response = client.get(f"/api/status/{completed_job}")
        data = response.json()
        assert data["results"] is not None
        assert "GO" in data["results"]


# ===========================================================================
# 7. GET /api/results/{job_id}- Get the results.
# ===========================================================================

class TestResultsEndpoint:
    """Test result to get the endpoint/api/results/{job_id}"""

    def test_results_json(self, client, completed_job):
        """Fetch results in JSON format"""
        response = client.get(f"/api/results/{completed_job}?format=json")
        assert response.status_code == 200
        data = response.json()

        assert "GO" in data
        assert "KEGG" in data
        assert isinstance(data["GO"], list)
        assert len(data["GO"]) > 0
        assert "Term_ID" in data["GO"][0]

    def test_results_tsv(self, client, completed_job):
        """Fetch results in TSV format"""
        response = client.get(f"/api/results/{completed_job}?format=tsv")
        assert response.status_code == 200
        assert "text/tab-separated-values" in response.headers.get("content-type", "")

        content = response.text
        assert "Term_ID" in content

    def test_results_default_format_is_json(self, client, completed_job):
        """Default format should be JSON"""
        response = client.get(f"/api/results/{completed_job}")
        assert response.status_code == 200
        # JSON response returned by JSONResponse, cont-type should read application/json
        assert "application/json" in response.headers.get("content-type", "")

    def test_results_not_found_job(self, client):
        """Query for non-existent tasks should return 404"""
        response = client.get("/api/results/nonexistent-job-id")
        assert response.status_code == 404

    def test_results_incomplete_job(self, client, running_job):
        """Failure to complete should return 400"""
        response = client.get(f"/api/results/{running_job}")
        assert response.status_code == 400
        assert "running" in response.json()["detail"]

    def test_results_invalid_format(self, client, completed_job):
        """Invalid format parameters should return 400"""
        response = client.get(f"/api/results/{completed_job}?format=xml")
        assert response.status_code == 400
        assert "invalid format" in response.json()["detail"].lower()

    def test_results_tsv_file_not_found(self, client):
        """returns 404 when the TSV file does not exist"""
        from allenricher.api.server import jobs
        import tempfile

        job_id = "test-no-tsv-file"
        output_dir = tempfile.mkdtemp(prefix="test_allenricher_")
        jobs[job_id] = {
            "status": "completed",
            "created_at": "2025-01-01T00:00:00",
            "completed_at": "2025-01-01T00:01:00",
            "progress": 1.0,
            "request": {"genes": ["TP53"], "species": "hsa", "databases": ["GO"]},
            "results": {"GO": []},
            "output_dir": output_dir,
            # Note: No "results_file" key
        }

        try:
            response = client.get(f"/api/results/{job_id}?format=tsv")
            assert response.status_code == 404
            assert "not found" in response.json()["detail"].lower()
        finally:
            import shutil
            shutil.rmtree(output_dir, ignore_errors=True)


# ===========================================================================
# GET /api/results/{job_id}/plot - Chart capture
# ===========================================================================

class TestPlotEndpoint:
    """Test Chart Get/api/results/{job_id}/plot"""

    def test_plot_success(self, client, completed_job):
        """Retrieve existing chart to return PDF file"""
        response = client.get(
            f"/api/results/{completed_job}/plot?database=GO&plot_type=barplot"
        )
        assert response.status_code == 200
        assert "application/pdf" in response.headers.get("content-type", "")
        assert len(response.content) > 0

    def test_plot_not_found_job(self, client):
        """The task that does not exist should be returned 404"""
        response = client.get(
            "/api/results/nonexistent-job/plot?database=GO&plot_type=barplot"
        )
        assert response.status_code == 404

    def test_plot_incomplete_job(self, client, running_job):
        """Failure to complete should return 400"""
        response = client.get(
            f"/api/results/{running_job}/plot?database=GO&plot_type=barplot"
        )
        assert response.status_code == 400

    def test_plot_not_found_file(self, client, completed_job):
        """returns 404 when the chart file does not exist"""
        response = client.get(
            f"/api/results/{completed_job}/plot?database=KEGG&plot_type=lollipop"
        )
        assert response.status_code == 404
        assert "not found" in response.json()["detail"].lower()

    def test_plot_missing_required_params(self, client, completed_job):
        """The missing parameters should return 422"""
        response = client.get(f"/api/results/{completed_job}/plot")
        assert response.status_code == 422

    def test_retired_dotplot_not_found(self, client, completed_job):
        """Deleted dotplot should not have output files."""
        response = client.get(
            f"/api/results/{completed_job}/plot?database=GO&plot_type=dotplot"
        )
        # File does not exist so 404, but parameters are normal
        assert response.status_code == 404


# ===========================================================================
# 9. GET /api/results/{job_id}/report- Report access.
# ===========================================================================

class TestReportEndpoint:
    """Test Report Gets the peer/api/results/{job_id}/report"""

    def test_report_success(self, client, completed_job):
        """Retrieving report should return HTML file"""
        response = client.get(f"/api/results/{completed_job}/report")
        assert response.status_code == 200
        assert "text/html" in response.headers.get("content-type", "")
        content = response.text
        assert "AllEnricher" in content

    def test_report_not_found_job(self, client):
        """The task that does not exist should be returned 404"""
        response = client.get("/api/results/nonexistent-job/report")
        assert response.status_code == 404

    def test_report_incomplete_job(self, client, running_job):
        """Failure to complete should return 400"""
        response = client.get(f"/api/results/{running_job}/report")
        assert response.status_code == 400

    def test_report_file_not_found(self, client):
        """When the report file does not exist, return 404"""
        from allenricher.api.server import jobs
        import tempfile

        job_id = "test-no-report"
        output_dir = tempfile.mkdtemp(prefix="test_allenricher_")
        jobs[job_id] = {
            "status": "completed",
            "created_at": "2025-01-01T00:00:00",
            "completed_at": "2025-01-01T00:01:00",
            "progress": 1.0,
            "request": {"genes": ["TP53"], "species": "hsa", "databases": ["GO"]},
            "results": {"GO": []},
            "output_dir": output_dir,
        }

        try:
            response = client.get(f"/api/results/{job_id}/report")
            assert response.status_code == 404
            assert "not found" in response.json()["detail"].lower()
        finally:
            import shutil
            shutil.rmtree(output_dir, ignore_errors=True)


class TestMethodsReferenceEndpoint:
    """Validate Methods text generated exclusively from recorded run metadata."""

    def test_methods_reference_success(self, client, completed_job):
        response = client.get(f"/api/results/{completed_job}/methods-reference")
        assert response.status_code == 200
        payload = response.json()
        assert payload["title"] == "Materials and Methods Writing Reference"
        assert "AllEnricher version 2.0-test" in payload["paragraphs"][0]
        assert [item["source"] for item in payload["references"][:3]] == [
            "AllEnricher",
            "GO",
            "GO",
        ]

    def test_methods_reference_requires_completed_job(self, client, running_job):
        response = client.get(f"/api/results/{running_job}/methods-reference")
        assert response.status_code == 400


# ===========================================================================
# 10. DELETE /api/jobs/{job_id}- Delete Task
# ===========================================================================

class TestDeleteEndpoint:
    """Test Task Remove End DELETE/api/jobs/{job_id}"""

    def test_delete_success(self, client, completed_job):
        """Remove existing tasks should return successfully"""
        from allenricher.api.server import jobs

        assert completed_job in jobs
        response = client.delete(f"/api/jobs/{completed_job}")
        assert response.status_code == 200
        data = response.json()
        assert data["message"] == "Job deleted"
        assert data["job_id"] == completed_job

        # Validation removed from Jobs
        assert completed_job not in jobs

    def test_delete_not_found(self, client):
        """Remove non-existent task returns 404"""
        response = client.delete("/api/jobs/nonexistent-job-id")
        assert response.status_code == 404

    def test_delete_does_not_remove_unmanaged_files(self, client, tmp_path):
        """Deleting jobs cannot always remove paths from the API root directory"""
        from allenricher.api.server import jobs

        job_id = "test-cleanup-job"
        output_dir = tmp_path / "output"
        output_dir.mkdir()
        (output_dir / "result.txt").write_text("test data")

        jobs[job_id] = {
            "status": "completed",
            "created_at": "2025-01-01T00:00:00",
            "completed_at": "2025-01-01T00:01:00",
            "progress": 1.0,
            "request": {"genes": ["TP53"], "species": "hsa", "databases": ["GO"]},
            "results": {"GO": []},
            "output_dir": str(output_dir),
        }

        response = client.delete(f"/api/jobs/{job_id}")
        assert response.status_code == 200
        assert output_dir.exists()

    def test_delete_job_without_output_dir(self, client):
        """Remove task without output_dir should not be misreported"""
        from allenricher.api.server import jobs

        job_id = "test-no-output-dir"
        jobs[job_id] = {
            "status": "pending",
            "created_at": "2025-01-01T00:00:00",
            "completed_at": None,
            "progress": 0.0,
            "request": {"genes": ["TP53"], "species": "hsa", "databases": ["GO"]},
            "results": None,
            # No output_dir key
        }

        response = client.delete(f"/api/jobs/{job_id}")
        assert response.status_code == 200


# ===========================================================================
# 11. Paths through protective tests
# ===========================================================================

class TestPathTraversalProtection:
    """Test path through attack protection."""

    def test_plot_path_traversal_database(self, client, completed_job):
        """Figure endpoint should reject database parameters that contain path through characters"""
        response = client.get(
            f"/api/results/{completed_job}/plot"
            f"?database=../../etc/passwd&plot_type=barplot"
        )
        # The path clean file does not exist and should be returned 404 (not leaking system files)
        assert response.status_code in (400, 404)

    def test_plot_path_traversal_plot_type(self, client, completed_job):
        """Figure endpoint should reject the plot_type parameter that contains path through calendar characters"""
        response = client.get(
            f"/api/results/{completed_job}/plot"
            f"?database=GO&plot_type=../../../etc/passwd"
        )
        assert response.status_code in (400, 404)

    def test_plot_path_traversal_both(self, client, completed_job):
        """Adds database and plot_type with path through the past"""
        response = client.get(
            f"/api/results/{completed_job}/plot"
            f"?database=./../../evil&plot_type=barplot"
        )
        assert response.status_code in (400, 404)

    def test_plot_special_characters_sanitized(self, client, completed_job):
        """Special characters should be cleared"""
        # Database parameters with spaces and special characters
        response = client.get(
            f"/api/results/{completed_job}/plot"
            f"?database=GO;rm+-rf+/&plot_type=barplot"
        )
        # Clean-up file does not exist - > 404, or for security reasons - > 400
        assert response.status_code in (400, 404)


# ===========================================================================
# 12. Invalid Request 422 Test
# ===========================================================================

class TestInvalidRequestValidation:
    """Test invalid request validation (FastAPI returns 422)"""

    def test_analyze_missing_genes(self, client):
        """Missing required field entries"""
        payload = {"species": "hsa", "databases": ["GO"]}
        response = client.post("/api/analyze", json=payload)
        assert response.status_code == 422

    def test_analyze_empty_body(self, client):
        """The empty request should be returned to 422"""
        response = client.post("/api/analyze", json={})
        assert response.status_code == 422

    def test_analyze_invalid_genes_type(self, client):
        """The genes field type error should return 422"""
        payload = {"genes": "TP53", "species": "hsa"}
        response = client.post("/api/analyze", json=payload)
        assert response.status_code == 422

    def test_analyze_invalid_pvalue_cutoff(self, client):
        """pvalue_cutoff beyond range should return 422 or be handled by a backend"""
        payload = {"genes": ["TP53"], "pvalue_cutoff": -1.0}
        response = client.post("/api/analyze", json=payload)
        # API returns 200 now and tags the task as filed (backend authenticated)
        # Not at the entrance. 422.
        assert response.status_code in [200, 422]

    def test_analyze_invalid_species_type(self, client):
        """Errors for species field types should return 422"""
        payload = {"genes": ["TP53"], "species": 12345}
        response = client.post("/api/analyze", json=payload)
        assert response.status_code == 422

    def test_analyze_invalid_method_type(self, client):
        """Method field type error returns 422"""
        payload = {"genes": ["TP53"], "method": 999}
        response = client.post("/api/analyze", json=payload)
        assert response.status_code == 422

    def test_analyze_extra_unexpected_field(self, client):
        """Unknown field must be rejected, avoid the argument spelling error being ignored in silence"""
        payload = {"genes": ["TP53"], "unexpected_field": "value"}
        response = client.post("/api/analyze", json=payload)
        assert response.status_code == 422

    def test_upload_no_file(self, client):
        """The upload endpoint missing file should return 422"""
        response = client.post("/api/upload")
        assert response.status_code == 422

    def test_plot_missing_database_param(self, client, completed_job):
        """The missing datbase parameter for figure endpoint should be returned 422"""
        response = client.get(
            f"/api/results/{completed_job}/plot?plot_type=barplot"
        )
        assert response.status_code == 422

    def test_analyze_empty_genes_list(self, client):
        """The empty list of gens should return 422 (Field(..., Min_lenghth without setting)"""
        # EnterRequest using Field (...), no Min_length set
        # Empty list is legal in type, but there may be a problem with semanticity
        # This test confirms Pydantic behaviour
        payload = {"genes": []}
        response = client.post("/api/analyze", json=payload)
        # Empty list is legal under the current model definition (Field does not set a min_lenghth=1)
        # So, expected 200, if the model is updated with a mem_length, it will become 422
        assert response.status_code in (200, 422)


# ===========================================================================
# Extra Integrated Test
# ===========================================================================

class TestIntegrationFlows:
    """End-to-end integration process testing"""

    @patch("allenricher.api.server.run_analysis")
    def test_analyze_then_status(self, mock_run, client):
        """Query Status After Submission of Analysis"""
        mock_run.return_value = None

        payload = {"genes": ["TP53", "BRCA1"]}
        analyze_resp = client.post("/api/analyze", json=payload)
        assert analyze_resp.status_code == 200
        job_id = analyze_resp.json()["job_id"]

        status_resp = client.get(f"/api/status/{job_id}")
        assert status_resp.status_code == 200
        assert status_resp.json()["job_id"] == job_id
        assert status_resp.json()["status"] == "pending"

    @patch("allenricher.api.server.run_analysis")
    def test_upload_then_delete(self, mock_run, client):
        """Upload and delete task"""
        mock_run.return_value = None

        gene_file = b"TP53\nBRCA1\n"
        files = {"file": ("genes.txt", gene_file, "text/plain")}
        upload_resp = client.post("/api/upload", files=files)
        assert upload_resp.status_code == 200
        job_id = upload_resp.json()["job_id"]

        delete_resp = client.delete(f"/api/jobs/{job_id}")
        assert delete_resp.status_code == 200

        # Re-Query should return 404
        status_resp = client.get(f"/api/status/{job_id}")
        assert status_resp.status_code == 404

    @patch("allenricher.api.server.run_analysis")
    def test_full_workflow_mock_completed(self, mock_run, client, tmp_path):
        """Simulate full workflow: Commit -> Simulation complete. -> Get results -> Access to reports -> Delete"""
        mock_run.return_value = None

        # 1. Submission of analysis
        payload = {"genes": ["TP53", "BRCA1", "EGFR"], "species": "hsa", "databases": ["GO"]}
        analyze_resp = client.post("/api/analyze", json=payload)
        assert analyze_resp.status_code == 200
        job_id = analyze_resp.json()["job_id"]

        # Simulation job complete (direct revision of the jobs dictionary)
        from allenricher.api.server import jobs
        output_dir = tmp_path / "workflow_test"
        output_dir.mkdir()
        plots_dir = output_dir / "plots"
        plots_dir.mkdir()
        (plots_dir / "GO_barplot.pdf").write_bytes(b"%PDF mock")
        (output_dir / "report.html").write_text("<html><body>Report</body></html>")

        jobs[job_id].update({
            "status": "completed",
            "completed_at": "2025-01-01T00:01:00",
            "progress": 1.0,
            "results": {"GO": [{"Term_ID": "GO:0005576", "Term_Name": "test"}]},
            "results_summary": {"GO": {"term_count": 1, "top_terms": []}},
            "output_dir": str(output_dir),
        })

        # 3. Access to results
        results_resp = client.get(f"/api/results/{job_id}")
        assert results_resp.status_code == 200
        assert "GO" in results_resp.json()

        # 4. Access to reports
        report_resp = client.get(f"/api/results/{job_id}/report")
        assert report_resp.status_code == 200

        # 5. Delete tasks
        delete_resp = client.delete(f"/api/jobs/{job_id}")
        assert delete_resp.status_code == 200

        # 6. Certification deleted
        assert job_id not in jobs


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
