"""Static checks for the reproducible container deployment contract."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_container_recipe_pins_runtime_and_runs_as_non_root() -> None:
    dockerfile = (ROOT / "Dockerfile").read_text(encoding="utf-8")

    assert "RELEASE_3_23@sha256:" in dockerfile
    assert 'as.character(getRversion()) == "4.6.1"' in dockerfile
    assert 'packageVersion("fgsea")) == "1.38.0"' in dockerfile
    assert "Archive/GSVA/GSVA_2.6.2.tar.gz" in dockerfile
    assert "133d1c3abb2bc886795d8ceb1a689d7a61e4cb4ee61f08a9ab0096a58af0064d" in dockerfile
    assert 'packageVersion("GSVA")) == "2.6.2"' in dockerfile
    assert "USER rstudio" in dockerfile
    assert 'ENTRYPOINT ["allenricher"]' in dockerfile
    assert ":latest" not in dockerfile


def test_container_constraints_cover_declared_runtime_dependencies() -> None:
    constraints = (ROOT / "docker" / "python-constraints.txt").read_text(
        encoding="utf-8"
    )
    required = {
        "fastapi",
        "Jinja2",
        "matplotlib",
        "numpy",
        "openpyxl",
        "pandas",
        "plotly",
        "python-multipart",
        "PyYAML",
        "requests",
        "scipy",
        "seaborn",
        "statsmodels",
        "tqdm",
        "uvicorn",
    }

    pinned = {line.split("==", 1)[0] for line in constraints.splitlines() if line}
    assert pinned == required


def test_readme_documents_container_build_and_batch_usage() -> None:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")

    assert "docker build --pull=false -t allenricher:2.1.1 ." in readme
    assert 'docker run --rm -v "${PWD}:/work"' in readme
    assert "serve --host 0.0.0.0 --port 8000" in readme
