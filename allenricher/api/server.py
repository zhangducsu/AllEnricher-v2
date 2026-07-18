"""Serve the local AllEnricher web application and REST API through the canonical CLI workflow."""

from __future__ import annotations

import json
import logging
import mimetypes
import os
import re
import shutil
import subprocess
import sys
import threading
import uuid
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional

import pandas as pd
from fastapi import BackgroundTasks, FastAPI, File, Form, HTTPException, Query, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, ConfigDict, Field

from allenricher import __version__
from allenricher.core.config import Config, DATABASE_CATALOG, SPECIES_CONFIGS
from allenricher.database.manager import DatabaseManager
from allenricher.report.methods_reference import build_methods_reference


logger = logging.getLogger(__name__)

MAX_UPLOAD_SIZE = 100 * 1024 * 1024
SAFE_JOB_ID = re.compile(r"^[A-Za-z0-9-]+$")
SAFE_PLOT_TOKEN = re.compile(r"^[A-Za-z0-9_.-]+$")
SAFE_SPECIES_CODE = re.compile(r"^[A-Za-z0-9_-]+$")
RESULT_SUFFIX = "_enrichment.tsv"
LOG_ARTIFACTS = {"command.txt", "stdout.log", "stderr.log"}
DEFAULT_JOB_ROOT = Path.home() / ".allenricher" / "api_jobs"
AI_BACKEND_LABELS = {
    "openai": "OpenAI",
    "claude": "Claude",
    "deepseek": "DeepSeek",
    "glm": "GLM",
    "minimax": "MiniMax",
    "ollama": "Ollama",
    "mock": "Mock (validation only)",
}
AI_BACKEND_NAMES = tuple(AI_BACKEND_LABELS)


app = FastAPI(
    title="AllEnricher API",
    description="Local REST API for enrichment analysis, visualization and HTML reports",
    version=__version__,
    docs_url="/docs",
    redoc_url="/redoc",
)

allowed_origins = [
    value.strip()
    for value in os.getenv("ALLENRICHER_CORS_ORIGINS", "").split(",")
    if value.strip()
]
if allowed_origins:
    app.add_middleware(
        CORSMiddleware,
        allow_origins=allowed_origins,
        allow_credentials=False,
        allow_methods=["GET", "POST", "DELETE"],
        allow_headers=["Content-Type"],
    )

static_dir = Path(__file__).parent / "static"
if static_dir.exists():
    app.mount("/static", StaticFiles(directory=str(static_dir)), name="static")


jobs: Dict[str, Dict[str, Any]] = {}
_jobs_lock = threading.RLock()


class RankedGene(BaseModel):
    model_config = ConfigDict(extra="forbid")

    gene: str = Field(min_length=1)
    weight: float


class EnrichmentRequest(BaseModel):
    """Define an analysis request using the same parameters as `allenricher analyze`."""

    model_config = ConfigDict(extra="forbid")

    genes: List[str] = Field(default_factory=list)
    ranked_genes: Optional[List[RankedGene]] = None
    expression_matrix: Optional[Dict[str, Dict[str, float]]] = None
    gmt_lines: Optional[List[str]] = None
    background: Optional[List[str]] = None

    species: str = "hsa"
    databases: List[str] = Field(default_factory=lambda: ["GO", "KEGG"], min_length=1)
    method: Literal["hypergeometric", "gsea", "ssgsea", "gsva"] = "hypergeometric"
    correction: Literal["BH", "BY", "bonferroni", "holm", "none"] = "BH"
    pvalue_cutoff: float = Field(default=0.05, gt=0, le=1)
    qvalue_cutoff: float = Field(default=0.05, gt=0, le=1)
    min_genes: int = Field(default=3, ge=1)
    background_mode: Literal["annotated", "genome", "custom"] = "annotated"
    jobs: int = Field(default=1, ge=1, le=128)
    database_dir: Optional[str] = None
    use_version: Optional[str] = None
    groups: Optional[str] = None

    no_plot: bool = False
    no_report: bool = False
    methods_language: Literal["zh", "en"] = "en"
    only_significant: bool = False
    plot_types: Optional[str] = None
    plot_format: Literal["png", "pdf", "svg"] = "png"
    plot_dpi: int = Field(default=300, ge=72, le=1200)
    style: Literal["nature", "science", "presentation", "cell", "omicshare"] = "nature"
    palette: Optional[str] = None
    categorical_palette: Optional[str] = None
    sequential_palette: Optional[str] = None
    diverging_palette: Optional[str] = None
    use_r_plots: bool = False
    ai_backend: Optional[Literal["openai", "claude", "deepseek", "glm", "minimax", "ollama", "mock"]] = None
    ai_mode: Literal["summary", "reviewer", "caption"] = "summary"
    ai_top_n: Optional[int] = Field(default=None, ge=1)

    tf_database: Optional[Literal["trrust", "chea3", "animaltfdb", "htftarget", "both"]] = None
    tf_library: Optional[str] = None
    tf_tissue: Optional[str] = None
    tf_regulation: Literal["all", "activation", "repression", "unknown"] = "all"
    tf_min_size: Optional[int] = Field(default=None, ge=1)
    tf_max_size: Optional[int] = Field(default=None, ge=1)
    tf_combine: Literal["none", "meanrank", "toprank"] = "none"
    tf_only: bool = False

    emapplot_qvalue: float = Field(default=0.05, gt=0, le=1)
    emapplot_min_count: int = Field(default=3, ge=1)
    emapplot_top_n: int = Field(default=30, ge=1)
    gsea_enrichment_top_up: int = Field(default=5, ge=0)
    gsea_enrichment_top_down: int = Field(default=5, ge=0)
    gsea_multi_top_up: int = Field(default=3, ge=0)
    gsea_multi_top_down: int = Field(default=3, ge=0)
    verbose: bool = False


class EnrichmentResponse(BaseModel):
    job_id: str
    status: str
    message: str
    results: Optional[Dict[str, Any]] = None


class JobStatusResponse(BaseModel):
    job_id: str
    status: str
    created_at: str
    completed_at: Optional[str] = None
    progress: float
    results: Optional[Dict[str, Any]] = None
    error: Optional[str] = None
    ai_interpretation_error: Optional[Dict[str, Any]] = None


class DatabaseInfoResponse(BaseModel):
    databases: List[Dict[str, Any]]


class SpeciesInfoResponse(BaseModel):
    code: str
    name: str
    taxonomy_id: int
    display_name: str
    databases: List[str] = Field(default_factory=list)


class SpeciesDatabaseSupportResponse(BaseModel):
    species: str
    databases: List[str]


class AIBackendStatus(BaseModel):
    name: str
    label: str
    configured: bool


class AIBackendStatusResponse(BaseModel):
    backends: List[AIBackendStatus]
    configuration_error: Optional[str] = None


def _job_root() -> Path:
    root = Path(os.getenv("ALLENRICHER_API_JOB_DIR", str(DEFAULT_JOB_ROOT))).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    return root


def _managed_task_dir(job_id: str) -> Path:
    if not SAFE_JOB_ID.fullmatch(job_id):
        raise HTTPException(status_code=400, detail="Invalid job ID")
    return (_job_root() / job_id).resolve()


def _disk_safe_job(job: Dict[str, Any]) -> Dict[str, Any]:
    payload = dict(job)
    request = dict(payload.get("request") or {})
    for key in ("genes", "ranked_genes", "expression_matrix", "gmt_lines", "background"):
        value = request.pop(key, None)
        if value is not None:
            request[f"{key}_count"] = len(value)
    payload["request"] = request
    return payload


def _persist_job(job_id: str) -> None:
    with _jobs_lock:
        job = jobs.get(job_id)
        if not job or not job.get("task_dir"):
            return
        task_dir = Path(job["task_dir"])
        task_dir.mkdir(parents=True, exist_ok=True)
        target = task_dir / "job.json"
        temporary = task_dir / "job.json.tmp"
        temporary.write_text(
            json.dumps(_disk_safe_job(job), ensure_ascii=False, indent=2, default=str),
            encoding="utf-8",
        )
        temporary.replace(target)


def _server_config() -> tuple[Config, Optional[str]]:
    """Load optional server-side settings without exposing credentials."""
    config_path = os.getenv("ALLENRICHER_CONFIG")
    if not config_path:
        return Config(), None
    try:
        return Config.from_file(config_path), None
    except (OSError, ValueError, TypeError) as exc:
        return Config(), f"Could not load server configuration: {exc}"


def _ai_backend_statuses() -> tuple[List[AIBackendStatus], Optional[str]]:
    config, configuration_error = _server_config()
    statuses = []
    for name in AI_BACKEND_NAMES:
        backend_config = config.get_ai_backend_config(name)
        enabled = backend_config.enabled if backend_config else True
        if name == "mock":
            configured = True
        elif name == "ollama":
            configured = enabled and (backend_config is not None or bool(os.getenv("OLLAMA_BASE_URL")))
        else:
            configured = enabled and bool(config.get_ai_api_key(name))
        statuses.append(AIBackendStatus(name=name, label=AI_BACKEND_LABELS[name], configured=configured))
    return statuses, configuration_error


def _validate_ai_backend(request: EnrichmentRequest) -> None:
    if not request.ai_backend:
        return
    statuses, configuration_error = _ai_backend_statuses()
    if configuration_error:
        raise HTTPException(status_code=422, detail=configuration_error)
    status = next(item for item in statuses if item.name == request.ai_backend)
    if not status.configured:
        raise HTTPException(
            status_code=422,
            detail=(
                f"AI backend '{request.ai_backend}' is not configured on this server. "
                "Set ALLENRICHER_CONFIG or the backend environment variable, then restart the server; "
                "use 'mock' for offline validation."
            ),
        )


def _get_job(job_id: str) -> Optional[Dict[str, Any]]:
    with _jobs_lock:
        if job_id in jobs:
            return jobs[job_id]
    if not SAFE_JOB_ID.fullmatch(job_id):
        return None
    job_file = _managed_task_dir(job_id) / "job.json"
    if not job_file.exists():
        return None
    try:
        job = json.loads(job_file.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        logger.exception("Failed to load API job metadata: %s", job_file)
        return None
    if not job.get("ai_interpretation_error"):
        error_file = Path(job.get("output_dir", "")) / "ai_interpretation_error.json"
        if error_file.is_file():
            try:
                job["ai_interpretation_error"] = json.loads(error_file.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                pass
    interpretation = job.get("ai_interpretation")
    if isinstance(interpretation, dict) and not interpretation.get("backend"):
        request_backend = (job.get("request") or {}).get("ai_backend")
        if request_backend:
            interpretation["backend"] = request_backend
    with _jobs_lock:
        jobs[job_id] = job
    return job


def _create_job(request: EnrichmentRequest) -> tuple[str, Dict[str, Any]]:
    job_id = str(uuid.uuid4())
    task_dir = _managed_task_dir(job_id)
    (task_dir / "input").mkdir(parents=True, exist_ok=False)
    (task_dir / "output").mkdir(parents=True, exist_ok=False)
    job = {
        "status": "pending",
        "created_at": datetime.now().isoformat(),
        "completed_at": None,
        "progress": 0.0,
        "request": request.model_dump(),
        "input_files": {},
        "results": None,
        "results_summary": None,
        "artifacts": [],
        "error": None,
        "task_dir": str(task_dir),
        "output_dir": str(task_dir / "output"),
    }
    with _jobs_lock:
        jobs[job_id] = job
    _persist_job(job_id)
    return job_id, job


def _write_gene_list(path: Path, genes: List[str]) -> None:
    cleaned = list(dict.fromkeys(str(gene).strip() for gene in genes if str(gene).strip()))
    path.write_text("\n".join(cleaned) + ("\n" if cleaned else ""), encoding="utf-8")


def _prepare_json_inputs(job: Dict[str, Any], request: EnrichmentRequest) -> None:
    input_dir = Path(job["task_dir"]) / "input"
    files: Dict[str, str] = {}

    if request.genes:
        gene_path = input_dir / "genes.txt"
        _write_gene_list(gene_path, request.genes)
        files["input"] = str(gene_path)
    if request.background:
        background_path = input_dir / "background.txt"
        _write_gene_list(background_path, request.background)
        files["background"] = str(background_path)
    if request.ranked_genes:
        ranked_path = input_dir / "ranked_genes.tsv"
        pd.DataFrame([item.model_dump() for item in request.ranked_genes]).to_csv(
            ranked_path, sep="\t", index=False
        )
        files["ranked"] = str(ranked_path)
    if request.expression_matrix:
        expression_path = input_dir / "expression.tsv"
        expression = pd.DataFrame.from_dict(request.expression_matrix, orient="index")
        expression.index.name = "gene"
        expression.to_csv(expression_path, sep="\t")
        files["expression"] = str(expression_path)
    if request.gmt_lines:
        gmt_path = input_dir / "gene_sets.gmt"
        gmt_path.write_text("\n".join(request.gmt_lines) + "\n", encoding="utf-8")
        files["gmt"] = str(gmt_path)

    job["input_files"] = files


async def _save_upload(upload: UploadFile, target: Path) -> bytes:
    content = await upload.read(MAX_UPLOAD_SIZE + 1)
    if len(content) > MAX_UPLOAD_SIZE:
        raise HTTPException(status_code=413, detail="Uploaded file exceeds 100 MB")
    try:
        content.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise HTTPException(status_code=422, detail=f"{upload.filename or 'file'} must be UTF-8") from exc
    target.write_bytes(content)
    return content


def _derive_input_genes(files: Dict[str, str], input_dir: Path) -> None:
    if files.get("input"):
        return
    genes: List[str] = []
    if files.get("ranked"):
        frame = pd.read_csv(files["ranked"], sep="\t")
        if "gene" not in frame.columns:
            raise HTTPException(status_code=422, detail="Ranked file requires gene and weight columns")
        genes = frame["gene"].dropna().astype(str).tolist()
    elif files.get("expression"):
        path = Path(files["expression"])
        sep = "," if path.suffix.lower() == ".csv" else "\t"
        frame = pd.read_csv(path, sep=sep, index_col=0)
        genes = frame.index.dropna().astype(str).tolist()
    if genes:
        input_path = input_dir / "genes.txt"
        _write_gene_list(input_path, genes)
        files["input"] = str(input_path)


def _validate_inputs(request: EnrichmentRequest, files: Dict[str, str]) -> None:
    if not files.get("input"):
        raise HTTPException(status_code=422, detail="A gene list or a method-specific input file is required")
    if request.method == "gsea" and not files.get("ranked"):
        raise HTTPException(status_code=422, detail="GSEA requires a ranked gene TSV file")
    if request.method in {"ssgsea", "gsva"} and not files.get("expression"):
        raise HTTPException(status_code=422, detail=f"{request.method} requires an expression matrix")
    if request.background_mode == "custom" and not files.get("background"):
        raise HTTPException(status_code=422, detail="Custom background mode requires a background file")
    if request.tf_only and not request.tf_database:
        raise HTTPException(status_code=422, detail="TF-only mode requires a TF database")
    if request.tf_database and request.method != "hypergeometric":
        raise HTTPException(
            status_code=422,
            detail=(
                "The legacy TF-only endpoint supports ORA only. For GSEA, ssGSEA, or GSVA, "
                "select a TF database in the main database list."
            ),
        )


def _append_option(command: List[str], flag: str, value: Any) -> None:
    if value is not None and value != "":
        command.extend([flag, str(value)])


def build_cli_command(request: EnrichmentRequest, files: Dict[str, str], output_dir: Path) -> List[str]:
    """Translate a validated API request into the canonical CLI command."""

    command = [
        sys.executable,
        "-m",
        "allenricher",
        "analyze",
        "-i",
        files["input"],
        "-s",
        request.species,
        "-d",
        ",".join(request.databases),
        "-o",
        str(output_dir),
        "-m",
        request.method,
        "-c",
        request.correction,
        "-p",
        str(request.pvalue_cutoff),
        "-q",
        str(request.qvalue_cutoff),
        "-n",
        str(request.min_genes),
        "-j",
        str(request.jobs),
        "--background-mode",
        request.background_mode,
        "--plot-format",
        request.plot_format,
        "--plot-dpi",
        str(request.plot_dpi),
        "--style",
        request.style,
        "--methods-language",
        request.methods_language,
        "--emapplot-qvalue",
        str(request.emapplot_qvalue),
        "--emapplot-min-count",
        str(request.emapplot_min_count),
        "--emapplot-top-n",
        str(request.emapplot_top_n),
        "--gsea-enrichment-top-up",
        str(request.gsea_enrichment_top_up),
        "--gsea-enrichment-top-down",
        str(request.gsea_enrichment_top_down),
        "--gsea-multi-top-up",
        str(request.gsea_multi_top_up),
        "--gsea-multi-top-down",
        str(request.gsea_multi_top_down),
    ]
    file_flags = {
        "background": "--background",
        "expression": "--expression-matrix",
        "ranked": "--ranked-genes",
        "gmt": "--gmt",
    }
    for key, flag in file_flags.items():
        if files.get(key):
            command.extend([flag, files[key]])

    _append_option(command, "--database-dir", request.database_dir)
    _append_option(command, "--config", os.getenv("ALLENRICHER_CONFIG"))
    _append_option(command, "--use-version", request.use_version)
    _append_option(command, "--groups", request.groups)
    _append_option(command, "--plot-types", request.plot_types)
    _append_option(command, "--palette", request.palette)
    _append_option(command, "--categorical-palette", request.categorical_palette)
    _append_option(command, "--sequential-palette", request.sequential_palette)
    _append_option(command, "--diverging-palette", request.diverging_palette)
    _append_option(command, "--tf-database", request.tf_database)
    _append_option(command, "--tf-library", request.tf_library)
    _append_option(command, "--tf-tissue", request.tf_tissue)
    _append_option(command, "--tf-regulation", request.tf_regulation)
    _append_option(command, "--tf-min-size", request.tf_min_size)
    _append_option(command, "--tf-max-size", request.tf_max_size)
    _append_option(command, "--tf-combine", request.tf_combine)
    if request.ai_backend:
        command.extend(["--ai", request.ai_backend, "--ai-mode", request.ai_mode])
        _append_option(command, "--ai-top-n", request.ai_top_n)

    for enabled, flag in (
        (request.no_plot, "--no-plot"),
        (request.no_report, "--no-report"),
        (request.only_significant, "--only-significant"),
        (request.tf_only, "--tf-only"),
        (request.use_r_plots, "--use-r-plots"),
        (request.verbose, "--verbose"),
    ):
        if enabled:
            command.append(flag)
    return command


def _read_result_tables(output_dir: Path) -> tuple[Dict[str, List[Dict[str, Any]]], Optional[Path]]:
    tables: Dict[str, List[Dict[str, Any]]] = {}
    combined: List[pd.DataFrame] = []
    for path in sorted(output_dir.glob(f"*{RESULT_SUFFIX}")):
        database = path.name[: -len(RESULT_SUFFIX)]
        frame = pd.read_csv(path, sep="\t")
        for column in ("Term_ID", "Pathway_ID", "DO_ID", "Disease_ID", "TF_ID"):
            if column in frame.columns:
                frame[column] = frame[column].astype("string")
        tables[database] = json.loads(frame.to_json(orient="records", force_ascii=False))
        export_frame = frame.copy()
        if "Database" not in export_frame.columns:
            export_frame.insert(0, "Database", database)
        combined.append(export_frame)

    if not combined:
        return tables, None
    results_file = output_dir / "enrichment_results.tsv"
    pd.concat(combined, ignore_index=True, sort=False).to_csv(results_file, sep="\t", index=False)
    return tables, results_file


def _artifact_category(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix in {".png", ".pdf", ".svg"}:
        return "figure"
    if suffix in {".html", ".htm"}:
        return "report"
    if suffix in {".tsv", ".csv", ".xlsx", ".xls"}:
        return "table"
    return "log" if path.name in LOG_ARTIFACTS else "other"


def _collect_artifacts(job: Dict[str, Any]) -> List[Dict[str, Any]]:
    task_dir = Path(job["task_dir"])
    output_dir = Path(job["output_dir"])
    paths = [path for path in output_dir.rglob("*") if path.is_file()]
    paths.extend(task_dir / name for name in sorted(LOG_ARTIFACTS) if (task_dir / name).is_file())
    artifacts = []
    for path in sorted(paths, key=lambda item: str(item).lower()):
        relative = path.relative_to(task_dir).as_posix()
        artifacts.append(
            {
                "name": path.name,
                "path": relative,
                "size": path.stat().st_size,
                "category": _artifact_category(path),
                "media_type": mimetypes.guess_type(path.name)[0] or "application/octet-stream",
            }
        )
    return artifacts


def run_analysis(job_id: str, request: EnrichmentRequest) -> None:
    """Execute one analysis job through the CLI and record its artifacts."""

    job = _get_job(job_id)
    if job is None:
        logger.error("API job disappeared before execution: %s", job_id)
        return
    try:
        job["status"] = "running"
        job["progress"] = 0.1
        _persist_job(job_id)

        files = dict(job.get("input_files") or {})
        input_dir = Path(job["task_dir"]) / "input"
        _derive_input_genes(files, input_dir)
        _validate_inputs(request, files)
        job["input_files"] = files

        output_dir = Path(job["output_dir"])
        command = build_cli_command(request, files, output_dir)
        display_command = subprocess.list2cmdline(command)
        (Path(job["task_dir"]) / "command.txt").write_text(display_command + "\n", encoding="utf-8")
        job["command"] = command
        job["progress"] = 0.2
        _persist_job(job_id)

        environment = os.environ.copy()
        environment.setdefault("PYTHONUTF8", "1")
        completed = subprocess.run(
            command,
            cwd=str(Path(__file__).resolve().parents[2]),
            env=environment,
            text=True,
            encoding="utf-8",
            errors="replace",
            capture_output=True,
            check=False,
        )
        task_dir = Path(job["task_dir"])
        (task_dir / "stdout.log").write_text(completed.stdout or "", encoding="utf-8")
        (task_dir / "stderr.log").write_text(completed.stderr or "", encoding="utf-8")
        job["return_code"] = completed.returncode
        if completed.returncode != 0:
            error_lines = [line for line in (completed.stderr or completed.stdout).splitlines() if line.strip()]
            detail = error_lines[-1] if error_lines else f"CLI exited with code {completed.returncode}"
            raise RuntimeError(detail)

        results, results_file = _read_result_tables(output_dir)
        job["results"] = results
        job["results_summary"] = {
            database: {"term_count": len(rows), "top_terms": rows[:10]}
            for database, rows in results.items()
        }
        if results_file:
            job["results_file"] = str(results_file)
        reports = sorted(output_dir.rglob("*report*.html"))
        if reports:
            job["report_file"] = str(reports[0])
        interpretation_file = output_dir / "ai_interpretation.json"
        if interpretation_file.is_file():
            try:
                job["ai_interpretation"] = json.loads(interpretation_file.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as exc:
                raise RuntimeError(f"AI interpretation output is invalid: {exc}") from exc
        interpretation_error_file = output_dir / "ai_interpretation_error.json"
        if interpretation_error_file.is_file():
            try:
                job["ai_interpretation_error"] = json.loads(
                    interpretation_error_file.read_text(encoding="utf-8")
                )
            except (OSError, json.JSONDecodeError) as exc:
                job["ai_interpretation_error"] = {
                    "error_code": "AI_INTERPRETATION_ERROR_RECORD_INVALID",
                    "message": str(exc),
                }
        job["artifacts"] = _collect_artifacts(job)
        job["status"] = "completed"
        job["progress"] = 1.0
        job["completed_at"] = datetime.now().isoformat()
    except Exception as exc:
        logger.exception("Analysis job %s failed", job_id)
        job["status"] = "failed"
        job["error"] = str(exc)
        job["completed_at"] = datetime.now().isoformat()
        job["artifacts"] = _collect_artifacts(job)
    finally:
        _persist_job(job_id)


def _submit(request: EnrichmentRequest, background_tasks: BackgroundTasks) -> EnrichmentResponse:
    _validate_ai_backend(request)
    job_id, job = _create_job(request)
    try:
        _prepare_json_inputs(job, request)
        _derive_input_genes(job["input_files"], Path(job["task_dir"]) / "input")
        _validate_inputs(request, job["input_files"])
        _persist_job(job_id)
    except Exception:
        shutil.rmtree(job["task_dir"], ignore_errors=True)
        with _jobs_lock:
            jobs.pop(job_id, None)
        raise
    background_tasks.add_task(run_analysis, job_id, request)
    return EnrichmentResponse(
        job_id=job_id,
        status="pending",
        message="Analysis started. Use /api/status/{job_id} to check progress.",
    )


@app.get("/", response_class=HTMLResponse)
async def root() -> Any:
    index_file = static_dir / "index.html"
    if index_file.exists():
        html = index_file.read_text(encoding="utf-8").replace(
            "{{ ALLENRICHER_VERSION }}", __version__
        )
        return HTMLResponse(html)
    return {"name": "AllEnricher API", "version": __version__, "docs": "/docs"}


def _available_databases(species: str) -> List[str]:
    manager = DatabaseManager(os.getenv("ALLENRICHER_DATABASE_DIR", "./database"), species)
    return [item["name"] for item in DATABASE_CATALOG if manager.has_database(item["name"])]


@app.get("/api/species", response_model=List[SpeciesInfoResponse])
async def get_species() -> List[SpeciesInfoResponse]:
    database_dir = Path(os.getenv("ALLENRICHER_DATABASE_DIR", "./database"))
    organism_dir = database_dir / "organism"
    built_codes = {
        species_dir.name
        for version_dir in organism_dir.iterdir()
        if version_dir.is_dir()
        for species_dir in version_dir.iterdir()
        if species_dir.is_dir()
    } if organism_dir.is_dir() else set()
    species: Dict[str, SpeciesInfoResponse] = {
        code: SpeciesInfoResponse(
            code=code,
            name=config.name,
            taxonomy_id=config.taxonomy_id,
            display_name=config.display_name,
            databases=_available_databases(code),
        )
        for code, config in SPECIES_CONFIGS.items()
    }
    try:
        from allenricher.database.species_registry import SpeciesRegistry

        registry = SpeciesRegistry.load_default(database_dir)
        for entry in registry.entries.values():
            code = entry.kegg_code if entry.has_kegg and entry.kegg_code else str(entry.taxid)
            if code in species or code not in built_codes:
                continue
            available = _available_databases(code)
            if not available:
                continue
            display_name = (entry.common_name or "").strip()
            if len(display_name) < 3 or display_name.startswith("%"):
                display_name = entry.latin_name
            registry_item = SpeciesInfoResponse(
                    code=code,
                    name=entry.latin_name,
                    taxonomy_id=entry.taxid,
                    display_name=display_name,
                    databases=available,
                )
            species[code] = registry_item
    except Exception:
        logger.debug("SpeciesRegistry unavailable; using built-in species", exc_info=True)
    for code in built_codes - set(species):
        species[code] = SpeciesInfoResponse(
            code=code,
            name=code,
            taxonomy_id=0,
            display_name=code,
            databases=_available_databases(code),
        )
    return sorted(species.values(), key=lambda item: (item.code != "hsa", item.display_name.lower()))


@app.get("/api/species/summary")
async def get_species_summary() -> Dict[str, Any]:
    try:
        from allenricher.database.species_registry import SpeciesRegistry

        registry = SpeciesRegistry.load_default(Path(os.getenv("ALLENRICHER_DATABASE_DIR", "./database")))
        return registry.get_summary()
    except Exception as exc:
        logger.exception("Failed to load species registry summary")
        raise HTTPException(status_code=503, detail="Species registry summary is unavailable") from exc


@app.get("/api/species/{species}/databases", response_model=SpeciesDatabaseSupportResponse)
async def get_species_databases(species: str) -> SpeciesDatabaseSupportResponse:
    code = species.strip().lower()
    if not SAFE_SPECIES_CODE.fullmatch(code):
        raise HTTPException(status_code=400, detail="Invalid species code")
    return SpeciesDatabaseSupportResponse(species=code, databases=_available_databases(code))


@app.get("/api/databases", response_model=DatabaseInfoResponse)
async def get_databases() -> DatabaseInfoResponse:
    return DatabaseInfoResponse(databases=[dict(item) for item in DATABASE_CATALOG])


@app.get("/api/ai/backends", response_model=AIBackendStatusResponse)
async def get_ai_backends() -> AIBackendStatusResponse:
    statuses, configuration_error = _ai_backend_statuses()
    return AIBackendStatusResponse(backends=statuses, configuration_error=configuration_error)


@app.post("/api/analyze", response_model=EnrichmentResponse)
async def analyze_genes(request: EnrichmentRequest, background_tasks: BackgroundTasks) -> EnrichmentResponse:
    return _submit(request, background_tasks)


@app.post("/api/upload", response_model=EnrichmentResponse)
async def upload_analysis(
    background_tasks: BackgroundTasks,
    file: Optional[UploadFile] = File(None),
    gene_file: Optional[UploadFile] = File(None),
    background_file: Optional[UploadFile] = File(None),
    ranked_file: Optional[UploadFile] = File(None),
    expression_file: Optional[UploadFile] = File(None),
    gmt_file: Optional[UploadFile] = File(None),
    species: str = Form("hsa"),
    databases: str = Form("GO,KEGG"),
    method: Literal["hypergeometric", "gsea", "ssgsea", "gsva"] = Form("hypergeometric"),
    correction: Literal["BH", "BY", "bonferroni", "holm", "none"] = Form("BH"),
    pvalue_cutoff: float = Form(0.05),
    qvalue_cutoff: float = Form(0.05),
    min_genes: int = Form(3),
    background_mode: Literal["annotated", "genome", "custom"] = Form("annotated"),
    jobs_count: int = Form(1, alias="jobs"),
    database_dir: Optional[str] = Form(None),
    use_version: Optional[str] = Form(None),
    groups: Optional[str] = Form(None),
    plot_types: Optional[str] = Form(None),
    plot_format: Literal["png", "pdf", "svg"] = Form("png"),
    plot_dpi: int = Form(300),
    style: Literal["nature", "science", "presentation", "cell", "omicshare"] = Form("nature"),
    categorical_palette: Optional[str] = Form(None),
    sequential_palette: Optional[str] = Form(None),
    diverging_palette: Optional[str] = Form(None),
    only_significant: bool = Form(False),
    no_plot: bool = Form(False),
    no_report: bool = Form(False),
    methods_language: Literal["zh", "en"] = Form("en"),
    use_r_plots: bool = Form(False),
    ai_backend: Optional[Literal["openai", "claude", "deepseek", "glm", "minimax", "ollama", "mock"]] = Form(None),
    ai_mode: Literal["summary", "reviewer", "caption"] = Form("summary"),
    ai_top_n: Optional[int] = Form(None),
    tf_database: Optional[Literal["trrust", "chea3", "animaltfdb", "htftarget", "both"]] = Form(None),
    tf_library: Optional[str] = Form(None),
    tf_tissue: Optional[str] = Form(None),
    tf_regulation: Literal["all", "activation", "repression", "unknown"] = Form("all"),
    tf_min_size: Optional[int] = Form(None),
    tf_max_size: Optional[int] = Form(None),
    tf_combine: Literal["none", "meanrank", "toprank"] = Form("none"),
    tf_only: bool = Form(False),
    emapplot_qvalue: float = Form(0.05),
    emapplot_min_count: int = Form(3),
    emapplot_top_n: int = Form(30),
    gsea_enrichment_top_up: int = Form(5),
    gsea_enrichment_top_down: int = Form(5),
    gsea_multi_top_up: int = Form(3),
    gsea_multi_top_down: int = Form(3),
    verbose: bool = Form(False),
) -> EnrichmentResponse:
    request = EnrichmentRequest(
        species=species,
        databases=[item.strip() for item in databases.split(",") if item.strip()],
        method=method,
        correction=correction,
        pvalue_cutoff=pvalue_cutoff,
        qvalue_cutoff=qvalue_cutoff,
        min_genes=min_genes,
        background_mode=background_mode,
        jobs=jobs_count,
        database_dir=database_dir,
        use_version=use_version,
        groups=groups,
        plot_types=plot_types,
        plot_format=plot_format,
        plot_dpi=plot_dpi,
        style=style,
        categorical_palette=categorical_palette,
        sequential_palette=sequential_palette,
        diverging_palette=diverging_palette,
        only_significant=only_significant,
        no_plot=no_plot,
        no_report=no_report,
        methods_language=methods_language,
        use_r_plots=use_r_plots,
        ai_backend=ai_backend,
        ai_mode=ai_mode,
        ai_top_n=ai_top_n,
        tf_database=tf_database,
        tf_library=tf_library,
        tf_tissue=tf_tissue,
        tf_regulation=tf_regulation,
        tf_min_size=tf_min_size,
        tf_max_size=tf_max_size,
        tf_combine=tf_combine,
        tf_only=tf_only,
        emapplot_qvalue=emapplot_qvalue,
        emapplot_min_count=emapplot_min_count,
        emapplot_top_n=emapplot_top_n,
        gsea_enrichment_top_up=gsea_enrichment_top_up,
        gsea_enrichment_top_down=gsea_enrichment_top_down,
        gsea_multi_top_up=gsea_multi_top_up,
        gsea_multi_top_down=gsea_multi_top_down,
        verbose=verbose,
    )
    _validate_ai_backend(request)
    job_id, job = _create_job(request)
    input_dir = Path(job["task_dir"]) / "input"
    files: Dict[str, str] = {}
    try:
        primary_gene_file = gene_file or file
        if primary_gene_file:
            target = input_dir / "genes.txt"
            content = await _save_upload(primary_gene_file, target)
            genes = [line.strip() for line in content.decode("utf-8-sig").splitlines() if line.strip()]
            request.genes = genes
            files["input"] = str(target)
        if background_file:
            target = input_dir / "background.txt"
            await _save_upload(background_file, target)
            files["background"] = str(target)
        if ranked_file:
            raw_target = input_dir / "ranked_upload.txt"
            await _save_upload(ranked_file, raw_target)
            ranked = pd.read_csv(raw_target, sep=None, engine="python")
            if not {"gene", "weight"}.issubset(ranked.columns):
                raise HTTPException(status_code=422, detail="Ranked file requires gene and weight columns")
            target = input_dir / "ranked_genes.tsv"
            ranked.to_csv(target, sep="\t", index=False)
            raw_target.unlink(missing_ok=True)
            files["ranked"] = str(target)
        if expression_file:
            suffix = ".csv" if (expression_file.filename or "").lower().endswith(".csv") else ".tsv"
            target = input_dir / f"expression{suffix}"
            await _save_upload(expression_file, target)
            files["expression"] = str(target)
        if gmt_file:
            target = input_dir / "gene_sets.gmt"
            await _save_upload(gmt_file, target)
            files["gmt"] = str(target)

        job["request"] = request.model_dump()
        job["input_files"] = files
        _derive_input_genes(files, input_dir)
        _validate_inputs(request, files)
        _persist_job(job_id)
    except Exception:
        shutil.rmtree(job["task_dir"], ignore_errors=True)
        with _jobs_lock:
            jobs.pop(job_id, None)
        raise

    background_tasks.add_task(run_analysis, job_id, request)
    return EnrichmentResponse(job_id=job_id, status="pending", message="Files uploaded and analysis started.")


@app.get("/api/status/{job_id}", response_model=JobStatusResponse)
async def get_job_status(job_id: str) -> JobStatusResponse:
    job = _get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    return JobStatusResponse(
        job_id=job_id,
        status=job["status"],
        created_at=job["created_at"],
        completed_at=job.get("completed_at"),
        progress=job.get("progress", 0.0),
        results=job.get("results_summary"),
        error=job.get("error"),
        ai_interpretation_error=job.get("ai_interpretation_error"),
    )


@app.get("/api/results/{job_id}")
async def get_results(job_id: str, format: str = Query("json")) -> Any:
    job = _get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    if job["status"] != "completed":
        raise HTTPException(status_code=400, detail=f"Job status: {job['status']}")
    if format not in {"json", "tsv"}:
        raise HTTPException(status_code=400, detail="Invalid format. Use 'json' or 'tsv'")
    if format == "json":
        return JSONResponse(content=job.get("results") or {})
    results_file = Path(job.get("results_file", ""))
    if not results_file.is_file():
        raise HTTPException(status_code=404, detail="Results file not found")
    return FileResponse(
        results_file,
        media_type="text/tab-separated-values",
        filename=f"enrichment_results_{job_id}.tsv",
    )


@app.get("/api/results/{job_id}/ai-interpretation")
async def get_ai_interpretation(job_id: str) -> Dict[str, Any]:
    """Return the validated structured AI interpretation for one completed job."""
    job = _get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    if job["status"] != "completed":
        raise HTTPException(status_code=400, detail=f"Job status: {job['status']}")
    interpretation = job.get("ai_interpretation")
    if not isinstance(interpretation, dict):
        raise HTTPException(status_code=404, detail="AI interpretation was not generated for this analysis")
    return interpretation


def _find_plot(job: Dict[str, Any], database: str, plot_type: str, file_format: Optional[str]) -> Optional[Path]:
    if not SAFE_PLOT_TOKEN.fullmatch(database) or not SAFE_PLOT_TOKEN.fullmatch(plot_type):
        raise HTTPException(status_code=400, detail="Invalid database or plot type")
    extensions = [file_format] if file_format else ["png", "pdf", "svg"]
    output_dir = Path(job.get("output_dir", ""))
    candidates = [path for path in output_dir.rglob("*") if path.is_file()]
    for extension in extensions:
        expected = f"{database}_{plot_type}.{extension}".lower()
        exact = next((path for path in candidates if path.name.lower() == expected), None)
        if exact:
            return exact
    database_token = database.lower()
    plot_token = plot_type.lower()
    for path in candidates:
        if path.suffix.lower().lstrip(".") not in extensions:
            continue
        name = path.stem.lower()
        if database_token in name and plot_token in name:
            return path
    return None


@app.get("/api/results/{job_id}/plot")
async def get_plot(
    job_id: str,
    database: str = Query(...),
    plot_type: str = Query("barplot"),
    format: Optional[Literal["png", "pdf", "svg"]] = Query(None),
) -> FileResponse:
    job = _get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    if job["status"] != "completed":
        raise HTTPException(status_code=400, detail=f"Job status: {job['status']}")
    plot_file = _find_plot(job, database, plot_type, format)
    if plot_file is None:
        raise HTTPException(status_code=404, detail="Plot not found")
    return FileResponse(
        plot_file,
        media_type=mimetypes.guess_type(plot_file.name)[0] or "application/octet-stream",
        headers={"Content-Disposition": f'inline; filename="{plot_file.name}"'},
    )


@app.get("/api/results/{job_id}/report")
async def get_report(job_id: str) -> FileResponse:
    job = _get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    if job["status"] != "completed":
        raise HTTPException(status_code=400, detail=f"Job status: {job['status']}")
    report_file = Path(job.get("report_file") or Path(job.get("output_dir", "")) / "report.html")
    if not report_file.is_file():
        raise HTTPException(status_code=404, detail="Report not found")
    return FileResponse(
        report_file,
        media_type="text/html",
        headers={"Content-Disposition": f'inline; filename="{report_file.name}"'},
    )


@app.get("/api/results/{job_id}/methods-reference")
async def get_methods_reference(job_id: str) -> Dict[str, Any]:
    job = _get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    if job["status"] != "completed":
        raise HTTPException(status_code=400, detail=f"Job status: {job['status']}")
    metadata_file = Path(job.get("output_dir", "")) / "analysis_metadata.json"
    if not metadata_file.is_file():
        raise HTTPException(status_code=404, detail="Analysis metadata not found")
    try:
        metadata = json.loads(metadata_file.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise HTTPException(status_code=500, detail="Analysis metadata is invalid") from exc
    return build_methods_reference(metadata)


@app.get("/api/results/{job_id}/artifacts")
async def get_artifacts(job_id: str) -> Dict[str, Any]:
    job = _get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    artifacts = _collect_artifacts(job) if job.get("task_dir") else job.get("artifacts", [])
    return {"job_id": job_id, "status": job["status"], "artifacts": artifacts}


def _resolve_artifact(job: Dict[str, Any], relative_path: str) -> Path:
    task_dir = Path(job.get("task_dir", "")).resolve()
    output_dir = Path(job.get("output_dir", "")).resolve()
    candidate = (task_dir / relative_path).resolve()
    allowed = candidate.name in LOG_ARTIFACTS and candidate.parent == task_dir
    allowed = allowed or candidate == output_dir or output_dir in candidate.parents
    if not allowed:
        raise HTTPException(status_code=400, detail="Invalid artifact path")
    return candidate


@app.get("/api/results/{job_id}/files/{file_path:path}")
async def get_artifact_file(job_id: str, file_path: str) -> FileResponse:
    job = _get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    artifact = _resolve_artifact(job, file_path)
    if not artifact.is_file():
        raise HTTPException(status_code=404, detail="Artifact not found")
    return FileResponse(artifact, media_type=mimetypes.guess_type(artifact.name)[0] or "application/octet-stream")


@app.delete("/api/jobs/{job_id}")
async def delete_job(job_id: str) -> Dict[str, str]:
    job = _get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    if job.get("status") == "running":
        raise HTTPException(status_code=409, detail="Running jobs cannot be deleted")
    managed_dir = _managed_task_dir(job_id)
    task_dir = Path(job.get("task_dir", managed_dir)).resolve()
    if task_dir == managed_dir and task_dir.exists():
        shutil.rmtree(task_dir)
    with _jobs_lock:
        jobs.pop(job_id, None)
    return {"message": "Job deleted", "job_id": job_id}


def start_api(host: str = "127.0.0.1", port: int = 8000, reload: bool = False) -> None:
    import uvicorn

    target: Any = "allenricher.api.server:app" if reload else app
    uvicorn.run(target, host=host, port=port, reload=reload)


if __name__ == "__main__":
    start_api()
