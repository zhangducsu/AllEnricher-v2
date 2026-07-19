from pathlib import Path

from allenricher.database.animaltfdb_fetcher import AnimalTFDBFetcher
from allenricher.database.htftarget_fetcher import HTFtargetFetcher
from allenricher.database import animaltfdb_fetcher


def test_database_fetchers_reject_html_error_pages(tmp_path: Path):
    html = tmp_path / "error.txt"
    html.write_text("<!doctype html><html>405</html>", encoding="utf-8")

    assert not AnimalTFDBFetcher._is_valid_tabular_file(html)
    assert not HTFtargetFetcher._is_valid_tabular_file(html)


def test_database_fetchers_accept_tabular_data(tmp_path: Path):
    table = tmp_path / "data.txt"
    table.write_text("TF\ttarget\nA\tB\n", encoding="utf-8")

    assert AnimalTFDBFetcher._is_valid_tabular_file(table)
    assert HTFtargetFetcher._is_valid_tabular_file(table)


def test_animaltfdb_waf_cookie_matches_live_challenge_example():
    assert AnimalTFDBFetcher._compute_waf_cookie(
        "B2778FB68C745E541D4909A7A40059082DA2D1AA"
    ) == "6a5733d1c8c5c09ca3e2acc5d08740b73ff0de37"


def test_animaltfdb_download_uses_http_direct_session_and_waf_cookie(tmp_path, monkeypatch):
    challenge = "<html><script>var arg1='B2778FB68C745E541D4909A7A40059082DA2D1AA';</script></html>"
    payload = b"Species\tSymbol\nDrosophila_melanogaster\tfd59A\n"

    class Response:
        def __init__(self, content):
            self.content = content
            self.text = content.decode("utf-8")

        @staticmethod
        def raise_for_status():
            return None

        def iter_content(self, chunk_size=8192):
            yield self.content

    class Cookies:
        def __init__(self):
            self.values = {}

        def set(self, name, value, **kwargs):
            self.values[name] = value

    class Session:
        def __init__(self):
            self.trust_env = True
            self.cookies = Cookies()
            self.calls = []

        def get(self, url, **kwargs):
            self.calls.append((url, kwargs))
            return Response(challenge.encode() if len(self.calls) == 1 else payload)

    session = Session()
    monkeypatch.setattr(animaltfdb_fetcher.requests, "Session", lambda: session)
    fetcher = AnimalTFDBFetcher(str(tmp_path))
    output = fetcher.download_tf_list("Drosophila_melanogaster")

    assert str(session.calls[0][0]).startswith("http://")
    assert session.trust_env is False
    assert session.cookies.values["acw_sc__v2"] == "6a5733d1c8c5c09ca3e2acc5d08740b73ff0de37"
    assert output.read_bytes() == payload
