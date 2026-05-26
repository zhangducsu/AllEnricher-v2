"""
REST API module for AllEnricher v2.0 using FastAPI

AllEnricher v2.0 REST API 服务模块
====================================

本模块基于 FastAPI 框架实现了基因集富集分析（Gene Set Enrichment Analysis）的 REST API 服务。
主要功能包括：

1. 基因列表的提交与文件上传
2. 异步执行富集分析任务（支持 GO、KEGG、Reactome 等多种数据库）
3. 任务状态查询与进度跟踪
4. 分析结果的获取（支持 JSON 和 TSV 格式）
5. 可视化图表的生成与下载（柱状图、气泡图、点图等）
6. HTML 交互式报告的生成与下载
7. 任务的删除与资源清理

本服务采用异步后台任务机制，提交分析请求后会立即返回任务ID，
客户端可通过轮询任务状态端点来获取分析进度和结果。

依赖模块：
    - allenricher.core.config: 配置管理
    - allenricher.core.enrichment: 富集分析核心引擎
    - allenricher.database.manager: 数据库管理器
    - allenricher.visualization.plotter: 可视化绑图工具
    - allenricher.report.generator: HTML 报告生成器
"""

import os
import logging
import tempfile
import uuid
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any

from fastapi import FastAPI, HTTPException, UploadFile, File, BackgroundTasks, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse, HTMLResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field
from pathlib import Path

from allenricher.core.config import Config
from allenricher.core.enrichment import EnrichmentAnalyzer
from allenricher.database.manager import DatabaseManager
from allenricher.visualization.plotter import Plotter
from allenricher.report.generator import ReportGenerator

logger = logging.getLogger(__name__)

# 创建 FastAPI 应用实例
# 配置 API 的基本信息，包括标题、描述、版本号，
# 以及 Swagger 文档（/docs）和 ReDoc 文档（/redoc）的访问路径
app = FastAPI(
    title="AllEnricher API",
    description="REST API for gene set enrichment analysis",
    version="2.0.0",
    docs_url="/docs",
    redoc_url="/redoc"
)

# 配置 CORS（跨域资源共享）中间件
# 允许所有来源（allow_origins=["*"]）的跨域请求，
# 支持携带凭证（cookies）、所有 HTTP 方法和所有请求头。
# 这使得前端应用可以从不同的域名/端口访问本 API 服务。
# CORS 配置：控制跨域请求策略
# 注意：allow_origins=["*"] 与 allow_credentials=True 同时使用时，
# 部分浏览器会拒绝该配置（CORS 规范不允许通配符来源携带凭据）。
# 生产环境中建议：1) 将 allow_origins 限制为具体前端域名，或
# 2) 设置 allow_credentials=False，或
# 3) 通过环境变量 ALLOWED_ORIGINS 动态配置（如 os.getenv("ALLOWED_ORIGINS", "*").split(",")）
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # 生产环境应限制为具体域名
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# 挂载静态文件目录（用于 Web 界面）
static_dir = Path(__file__).parent / "static"
if static_dir.exists():
    app.mount("/static", StaticFiles(directory=str(static_dir)), name="static")

# 任务存储字典（内存中存储所有分析任务的状态和结果）
# 键为任务ID（UUID字符串），值为包含任务详细信息的字典
# 生产环境中应替换为持久化数据库（如 Redis、PostgreSQL 等）
jobs: Dict[str, Dict[str, Any]] = {}


# ===================== 请求/响应数据模型 =====================

class EnrichmentRequest(BaseModel):
    """
    富集分析请求模型

    用于接收客户端提交的富集分析参数，包括基因列表、物种、
    目标数据库、统计方法、多重检验校正方法及各种阈值设置。

    Attributes:
        genes: 待分析的基因符号列表（必填）
        species: 物种代码，默认为 "hsa"（人类），也支持 "mmu"（小鼠）等
        databases: 要分析的数据库列表，默认为 GO 和 KEGG
        method: 富集分析方法，默认为 Fisher 精确检验（"fisher"）
        correction: 多重检验校正方法，默认为 BH（Benjamini-Hochberg）
        pvalue_cutoff: p 值显著性阈值，默认为 0.05
        qvalue_cutoff: q 值（校正后 p 值）阈值，默认为 0.05
        min_genes: 每个功能条目中最少需要的基因数，默认为 2
        background: 自定义背景基因列表，可选；若不提供则使用数据库默认背景
    """
    genes: List[str] = Field(..., description="List of gene symbols")
    species: str = Field(default="hsa", description="Species code (e.g., hsa, mmu)")
    databases: List[str] = Field(default=["GO", "KEGG"], description="Databases to analyze")
    method: str = Field(default="fisher", description="Enrichment method")
    correction: str = Field(default="BH", description="Multiple testing correction")
    pvalue_cutoff: float = Field(default=0.05, description="P-value cutoff")
    qvalue_cutoff: float = Field(default=0.05, description="Q-value cutoff")
    min_genes: int = Field(default=2, description="Minimum genes per term")
    background: Optional[List[str]] = Field(default=None, description="Background gene list")


class EnrichmentResponse(BaseModel):
    """
    富集分析响应模型

    当客户端提交分析请求后，返回此响应模型，包含任务ID、
    当前状态、提示信息以及可选的分析结果。

    Attributes:
        job_id: 分析任务的唯一标识符（UUID 格式）
        status: 任务当前状态（如 "pending"、"running"、"completed"、"failed"）
        message: 状态描述信息，用于提示用户下一步操作
        results: 分析结果数据，仅在任务完成时包含，默认为 None
    """
    job_id: str
    status: str
    message: str
    results: Optional[Dict[str, Any]] = None


class JobStatusResponse(BaseModel):
    """
    任务状态响应模型

    用于查询分析任务的详细状态信息，包括创建时间、完成时间、
    执行进度、结果摘要以及错误信息。

    Attributes:
        job_id: 分析任务的唯一标识符
        status: 任务当前状态（"pending" / "running" / "completed" / "failed"）
        created_at: 任务创建时间（ISO 8601 格式）
        completed_at: 任务完成时间，任务未完成时为 None
        progress: 任务执行进度（0.0 到 1.0 之间的浮点数）
        results: 分析结果摘要信息，仅在任务完成时包含
        error: 错误信息，仅在任务失败时包含
    """
    job_id: str
    status: str
    created_at: str
    completed_at: Optional[str] = None
    progress: float
    results: Optional[Dict[str, Any]] = None
    error: Optional[str] = None


class DatabaseInfoResponse(BaseModel):
    """
    数据库信息响应模型

    返回系统中可用的富集分析数据库列表及其描述信息。

    Attributes:
        databases: 数据库信息列表，每个元素为包含数据库名称、描述和支持物种的字典
    """
    databases: List[Dict[str, Any]]


class SpeciesInfoResponse(BaseModel):
    """
    物种信息 API 响应模型

    描述系统支持的物种信息，包括物种代码、名称、分类学ID和显示名称。
    用于 /api/species 端点的响应序列化。

    注意：内部物种数据结构使用 allenricher.database.species_lookup.SpeciesInfo（dataclass），
    本类仅用于 API 层的请求/响应序列化（Pydantic BaseModel）。

    Attributes:
        code: 物种代码（如 "hsa" 表示人类，"mmu" 表示小鼠）
        name: 物种拉丁学名（如 "Homo sapiens"）
        taxonomy_id: NCBI 分类学ID（如人类为 9606）
        display_name: 物种的常用显示名称（如 "Human"、"Mouse"）
    """
    code: str
    name: str
    taxonomy_id: int
    display_name: str


# ===================== API 端点 =====================

@app.get("/", response_class=HTMLResponse)
async def root():
    """
    根端点 (GET /)

    返回 Web 分析界面。如果 static/index.html 存在则返回该页面，
    否则返回 API 基本信息（JSON 格式）。
    """
    index_file = static_dir / "index.html"
    if index_file.exists():
        with open(index_file, "r", encoding="utf-8") as f:
            return f.read()
    
    # 如果没有静态文件，返回 API 信息
    return {
        "name": "AllEnricher API",
        "version": "2.0.0",
        "docs": "/docs",
        "webui": "/static/index.html (if available)",
        "endpoints": {
            "analyze": "/api/analyze",
            "upload": "/api/upload",
            "status": "/api/status/{job_id}",
            "results": "/api/results/{job_id}",
            "databases": "/api/databases",
            "species": "/api/species"
        }
    }


@app.get("/api/species", response_model=List[SpeciesInfoResponse])
async def get_species():
    """
    获取支持的物种列表 (GET /api/species)

    返回系统中所有支持进行富集分析的物种信息列表。
    每个物种包含代码、拉丁学名、NCBI 分类学 ID 和显示名称。
    客户端可使用返回的物种代码作为分析请求中的 species 参数。

    Returns:
        List[SpeciesInfoResponse]: 支持的物种信息列表
    """
    from allenricher.core.config import SPECIES_CONFIGS
    
    species_list = []
    for code, config in SPECIES_CONFIGS.items():
        species_list.append(SpeciesInfoResponse(
            code=code,
            name=config.name,
            taxonomy_id=config.taxonomy_id,
            display_name=config.display_name
        ))
    
    return species_list


@app.get("/api/databases", response_model=DatabaseInfoResponse)
async def get_databases():
    """
    获取可用数据库列表 (GET /api/databases)

    返回系统中所有可用的富集分析数据库信息，包括数据库名称、
    描述和支持的物种范围。当前支持的数据库包括：
    - GO（Gene Ontology）：基因本体数据库
    - KEGG：KEGG 通路数据库
    - Reactome：Reactome 通路数据库
    - WikiPathways：WikiPathways 通路数据库
    - MSigDB：分子特征数据库（仅支持人类）
    - DO（Disease Ontology）：疾病本体数据库（仅支持人类）
    - DisGeNET：疾病-基因关联数据库（仅支持人类）

    Returns:
        DatabaseInfoResponse: 包含所有可用数据库信息的响应对象
    """
    databases = [
        {"name": "GO", "description": "Gene Ontology", "species": "all"},
        {"name": "KEGG", "description": "KEGG Pathways", "species": "all"},
        {"name": "Reactome", "description": "Reactome Pathways", "species": "model_organisms"},
        {"name": "WikiPathways", "description": "WikiPathways", "species": "all"},
        {"name": "MSigDB", "description": "Molecular Signatures Database", "species": "hsa"},
        {"name": "DO", "description": "Disease Ontology", "species": "hsa"},
        {"name": "DisGeNET", "description": "Disease-Gene Associations", "species": "hsa"}
    ]
    
    return DatabaseInfoResponse(databases=databases)


@app.post("/api/analyze", response_model=EnrichmentResponse)
async def analyze_genes(request: EnrichmentRequest, background_tasks: BackgroundTasks):
    """
    提交富集分析任务 (POST /api/analyze)

    接收包含基因列表和分析参数的 JSON 请求体，创建异步后台分析任务。
    分析任务会在后台执行，本端点立即返回任务ID和状态信息。
    客户端可使用返回的 job_id 通过 GET /api/status/{job_id} 查询任务进度。

    请求体示例：
        {
            "genes": ["TP53", "BRCA1", "EGFR"],
            "species": "hsa",
            "databases": ["GO", "KEGG"],
            "method": "fisher",
            "correction": "BH",
            "pvalue_cutoff": 0.05
        }

    Args:
        request: 富集分析请求参数（EnrichmentRequest 模型）
        background_tasks: FastAPI 后台任务管理器，用于调度异步分析任务

    Returns:
        EnrichmentResponse: 包含任务ID、状态和提示信息的响应
    """
    # 生成唯一的任务ID（使用 UUID v4）
    job_id = str(uuid.uuid4())
    
    # 初始化任务记录，设置初始状态为 "pending"（等待中）
    jobs[job_id] = {
        "status": "pending",
        "created_at": datetime.now().isoformat(),
        "completed_at": None,
        "progress": 0.0,
        "request": request.model_dump(),
        "results": None,
        "error": None
    }
    
    # 将分析任务添加到后台任务队列中异步执行
    background_tasks.add_task(run_analysis, job_id, request)
    
    return EnrichmentResponse(
        job_id=job_id,
        status="pending",
        message="Analysis started. Use /api/status/{job_id} to check progress."
    )


@app.post("/api/upload", response_model=EnrichmentResponse)
async def upload_gene_list(
    background_tasks: BackgroundTasks,
    file: UploadFile = File(..., description="Gene list file (one gene per line)"),
    species: str = Query("hsa", description="Species code"),
    databases: str = Query("GO,KEGG", description="Comma-separated database list"),
    method: str = Query("fisher", description="Enrichment method")
):
    """
    上传基因列表文件并提交分析任务 (POST /api/upload)

    接收一个基因列表文本文件（每行一个基因符号），以及通过查询参数指定的
    分析配置（物种、数据库、分析方法等）。读取文件内容后自动创建
    富集分析任务并在后台异步执行。

    文件格式要求：纯文本文件，每行一个基因符号，支持 UTF-8 编码。

    Args:
        background_tasks: FastAPI 后台任务管理器
        file: 上传的基因列表文件（每行一个基因符号）
        species: 物种代码，默认为 "hsa"（人类）
        databases: 逗号分隔的数据库名称列表，默认为 "GO,KEGG"
        method: 富集分析方法，默认为 "fisher"

    Returns:
        EnrichmentResponse: 包含任务ID、状态和提示信息的响应
    """
    # 读取上传文件内容并解析基因列表（限制最大 10MB 防止内存耗尽）
    MAX_UPLOAD_SIZE = 10 * 1024 * 1024  # 10MB
    content = await file.read(MAX_UPLOAD_SIZE + 1)
    if len(content) > MAX_UPLOAD_SIZE:
        raise HTTPException(status_code=413, detail="上传文件过大，最大支持 10MB")
    genes = [g.strip() for g in content.decode('utf-8').split('\n') if g.strip()]
    
    # 根据解析的基因列表和查询参数构建富集分析请求对象
    request = EnrichmentRequest(
        genes=genes,
        species=species,
        databases=databases.split(','),
        method=method
    )
    
    # 生成唯一的任务ID
    job_id = str(uuid.uuid4())
    
    # 初始化任务记录
    jobs[job_id] = {
        "status": "pending",
        "created_at": datetime.now().isoformat(),
        "completed_at": None,
        "progress": 0.0,
        "request": request.model_dump(),
        "results": None,
        "error": None
    }
    
    # 将分析任务添加到后台任务队列中异步执行
    background_tasks.add_task(run_analysis, job_id, request)
    
    return EnrichmentResponse(
        job_id=job_id,
        status="pending",
        message="File uploaded and analysis started."
    )


@app.get("/api/status/{job_id}", response_model=JobStatusResponse)
async def get_job_status(job_id: str):
    """
    查询任务状态 (GET /api/status/{job_id})

    根据任务ID查询富集分析任务的当前状态，包括执行进度、
    创建时间、完成时间和错误信息（如果任务失败）。

    任务状态说明：
        - pending: 任务已创建，等待执行
        - running: 任务正在执行中
        - completed: 任务已完成，可获取结果
        - failed: 任务执行失败，可查看错误信息

    Args:
        job_id: 分析任务的唯一标识符（UUID 格式）

    Returns:
        JobStatusResponse: 包含任务详细状态信息的响应对象

    Raises:
        HTTPException 404: 当指定的任务ID不存在时
    """
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = jobs[job_id]
    
    return JobStatusResponse(
        job_id=job_id,
        status=job["status"],
        created_at=job["created_at"],
        completed_at=job.get("completed_at"),
        progress=job["progress"],
        results=job.get("results_summary"),
        error=job.get("error")
    )


@app.get("/api/results/{job_id}")
async def get_results(job_id: str, format: str = Query("json", description="Output format: json, tsv")):
    """
    获取分析结果 (GET /api/results/{job_id})

    获取已完成分析任务的详细结果数据。支持两种输出格式：
    - json: 返回 JSON 格式的完整结果数据（默认）
    - tsv: 返回 TSV（制表符分隔值）格式的结果文件下载

    注意：仅当任务状态为 "completed" 时才能获取结果。

    Args:
        job_id: 分析任务的唯一标识符
        format: 输出格式，可选 "json" 或 "tsv"，默认为 "json"

    Returns:
        JSON 格式的结果数据或 TSV 文件的 FileResponse

    Raises:
        HTTPException 404: 任务ID不存在或 TSV 结果文件未找到
        HTTPException 400: 任务尚未完成或格式参数无效
    """
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = jobs[job_id]
    
    if job["status"] != "completed":
        raise HTTPException(status_code=400, detail=f"Job status: {job['status']}")
    
    if format == "json":
        # 返回 JSON 格式的完整结果数据
        return JSONResponse(content=job["results"])
    elif format == "tsv":
        # 返回 TSV 格式的结果文件
        if "results_file" in job:
            return FileResponse(
                job["results_file"],
                media_type="text/tab-separated-values",
                filename=f"enrichment_results_{job_id}.tsv"
            )
        else:
            raise HTTPException(status_code=404, detail="Results file not found")
    else:
        raise HTTPException(status_code=400, detail="Invalid format. Use 'json' or 'tsv'")


@app.get("/api/results/{job_id}/plot")
async def get_plot(
    job_id: str,
    database: str = Query(..., description="Database name"),
    plot_type: str = Query("barplot", description="Plot type: barplot, bubble, dotplot")
):
    """
    获取可视化图表 (GET /api/results/{job_id}/plot)

    获取已完成分析任务的可视化图表文件（PDF 格式）。
    需要指定数据库名称和图表类型。

    支持的图表类型：
        - barplot: 柱状图（默认）
        - bubble: 气泡图
        - dotplot: 点图

    Args:
        job_id: 分析任务的唯一标识符
        database: 数据库名称（如 "GO"、"KEGG" 等）
        plot_type: 图表类型，可选 "barplot"、"bubble" 或 "dotplot"

    Returns:
        FileResponse: PDF 格式的图表文件

    Raises:
        HTTPException 404: 任务ID不存在或图表文件未找到
        HTTPException 400: 任务尚未完成
    """
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = jobs[job_id]
    
    if job["status"] != "completed":
        raise HTTPException(status_code=400, detail=f"Job status: {job['status']}")
    
    # 根据数据库名称和图表类型构建图表文件路径
    # 安全措施：清理用户输入，防止路径遍历攻击
    import re
    safe_database = re.sub(r'[^\w\-]', '', database)
    safe_plot_type = re.sub(r'[^\w\-]', '', plot_type)
    plot_dir = Path(job.get("output_dir", "")) / "plots"
    plot_file = plot_dir / f"{safe_database}_{safe_plot_type}.pdf"
    
    # 验证解析后的路径仍在预期的 plots 目录内，防止路径遍历
    resolved_plot_dir = plot_dir.resolve()
    resolved_plot_file = plot_file.resolve()
    if not str(resolved_plot_file).startswith(str(resolved_plot_dir) + os.sep) and resolved_plot_file != resolved_plot_dir:
        raise HTTPException(status_code=400, detail="Invalid database or plot type")
    
    if not plot_file.exists():
        raise HTTPException(status_code=404, detail="Plot not found")
    
    return FileResponse(
        plot_file,
        media_type="application/pdf",
        filename=f"{database}_{plot_type}.pdf"
    )


@app.get("/api/results/{job_id}/report")
async def get_report(job_id: str):
    """
    获取 HTML 报告 (GET /api/results/{job_id}/report)

    获取已完成分析任务的 HTML 交互式报告文件。
    报告包含所有数据库的富集分析结果汇总、可视化图表和详细数据表格。

    Args:
        job_id: 分析任务的唯一标识符

    Returns:
        FileResponse: HTML 格式的报告文件

    Raises:
        HTTPException 404: 任务ID不存在或报告文件未找到
        HTTPException 400: 任务尚未完成
    """
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = jobs[job_id]
    
    if job["status"] != "completed":
        raise HTTPException(status_code=400, detail=f"Job status: {job['status']}")
    
    report_file = Path(job.get("output_dir", "")) / "report.html"
    
    if not report_file.exists():
        raise HTTPException(status_code=404, detail="Report not found")
    
    return FileResponse(
        report_file,
        media_type="text/html",
        filename=f"enrichment_report_{job_id}.html"
    )


@app.delete("/api/jobs/{job_id}")
async def delete_job(job_id: str):
    """
    删除任务 (DELETE /api/jobs/{job_id})

    删除指定的分析任务及其所有相关文件（包括输出目录、
    图表文件、报告文件等）。删除后任务ID将不再可用。

    Args:
        job_id: 要删除的分析任务的唯一标识符

    Returns:
        dict: 包含删除确认信息的字典

    Raises:
        HTTPException 404: 当指定的任务ID不存在时
    """
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    # 清理任务相关的输出文件和目录
    job = jobs[job_id]
    if "output_dir" in job:
        output_dir = Path(job["output_dir"])
        if output_dir.exists():
            import shutil
            shutil.rmtree(output_dir)  # 递归删除整个输出目录及其内容
    
    del jobs[job_id]
    
    return {"message": "Job deleted", "job_id": job_id}


# ===================== 后台任务函数 =====================

def run_analysis(job_id: str, request: EnrichmentRequest):
    """
    后台执行富集分析任务

    该函数作为后台任务运行（通过 BackgroundTasks），执行完整的富集分析流程。
    注意：此函数为同步函数，因为内部调用的富集分析、数据库加载等均为同步阻塞操作。
    FastAPI 的 BackgroundTasks 会自动在线程池中执行同步后台任务。
    分析过程分为以下几个阶段：

    1. 初始化阶段（进度 0.0 - 0.1）：
       - 创建临时输出目录
       - 根据请求参数构建配置对象

    2. 数据加载阶段（进度 0.1 - 0.3）：
       - 加载指定的数据库文件
       - 准备背景基因集

    3. 分析执行阶段（进度 0.3 - 0.7）：
       - 对每个数据库执行富集分析
       - 应用统计检验和多重检验校正

    4. 可视化生成阶段（进度 0.7 - 0.85）：
       - 为每个数据库的结果生成多种图表

    5. 报告生成阶段（进度 0.85 - 0.95）：
       - 生成 HTML 交互式报告

    6. 结果保存阶段（进度 0.95 - 1.0）：
       - 保存结果文件和摘要信息

    如果任何阶段发生异常，任务状态将设置为 "failed" 并记录错误信息。

    Args:
        job_id: 分析任务的唯一标识符
        request: 富集分析请求参数对象
    """
    job = jobs[job_id]
    
    try:
        # ---- 阶段1：更新任务状态为"运行中" ----
        job["status"] = "running"
        job["progress"] = 0.1
        
        # 创建临时输出目录，目录名前缀包含任务ID以便识别
        output_dir = Path(tempfile.mkdtemp(prefix=f"allenricher_{job_id}_"))
        job["output_dir"] = str(output_dir)
        
        # 根据请求参数创建分析配置对象
        config = Config(
            species=request.species,
            databases=request.databases,
            method=request.method,
            correction=request.correction,
            pvalue_cutoff=request.pvalue_cutoff,
            qvalue_cutoff=request.qvalue_cutoff,
            min_genes=request.min_genes,
            output_dir=str(output_dir)
        )
        
        # ---- 阶段2：加载指定的数据库 ----
        job["progress"] = 0.2
        db_manager = DatabaseManager(config.database_dir, request.species)
        db_manager.load_databases(request.databases)
        
        # 获取背景基因集：如果用户提供了自定义背景则使用自定义的，否则使用数据库默认背景
        if request.background:
            background_set = set(request.background)
        else:
            background_set = db_manager.get_background_genes()
        
        # ---- 阶段3：执行富集分析 ----
        job["progress"] = 0.3
        analyzer = EnrichmentAnalyzer(config)
        
        # 将基因列表转换为集合去重，并获取所有数据库的功能条目数据
        gene_set = set(request.genes)
        database_data = db_manager.get_all_term_data()
        
        # 运行富集分析，返回各数据库的分析结果（DataFrame 字典）
        results = analyzer.run_analysis(gene_set, background_set, database_data)
        
        job["progress"] = 0.7
        
        # ---- 阶段4：生成可视化图表 ----
        plotter = Plotter(str(output_dir / "plots"), config)
        for db_name, df in results.items():
            plotter.plot_all(df, db_name)  # 为每个数据库生成所有类型的图表
        
        job["progress"] = 0.85
        
        # ---- 阶段5：生成 HTML 报告 ----
        report_gen = ReportGenerator(str(output_dir), config)
        report_gen.generate(results, str(output_dir / "report.html"))
        
        job["progress"] = 0.95
        
        # ---- 阶段6：保存结果文件 ----
        analyzer.save_results(str(output_dir))
        
        # 准备结果摘要信息（包含每个数据库的条目数量和前10个显著条目）
        results_summary = {}
        for db_name, df in results.items():
            results_summary[db_name] = {
                "term_count": len(df),
                "top_terms": df.head(10).to_dict(orient="records") if len(df) > 0 else []
            }
        
        # 更新任务状态为"已完成"，记录完成时间和完整结果
        job["status"] = "completed"
        job["completed_at"] = datetime.now().isoformat()
        job["progress"] = 1.0
        job["results"] = {k: v.to_dict(orient="records") for k, v in results.items()}
        job["results_summary"] = results_summary
        
    except Exception as e:
        # 捕获所有异常，记录错误日志并将任务状态设置为"失败"
        logger.error(f"Error in job {job_id}: {e}")
        job["status"] = "failed"
        job["error"] = str(e)
        job["completed_at"] = datetime.now().isoformat()


def start_api(host: str = "0.0.0.0", port: int = 8000):
    """
    启动 API 服务

    使用 uvicorn ASGI 服务器启动 FastAPI 应用。

    Args:
        host: 服务监听地址，默认为 "0.0.0.0"（所有网络接口）
        port: 服务监听端口，默认为 8000
    """
    import uvicorn
    uvicorn.run(app, host=host, port=port)


if __name__ == "__main__":
    start_api()

