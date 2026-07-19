from pathlib import Path


STATIC_HTML = Path(__file__).resolve().parents[1] / "allenricher" / "api" / "static" / "index.html"


def test_web_database_defaults_do_not_prefer_go_kegg():
    html = STATIC_HTML.read_text(encoding="utf-8")

    assert '["GO", "KEGG"].includes(database.name)' not in html
    assert "input.checked = supported && (previous.has(database.name) || !previous.size);" in html