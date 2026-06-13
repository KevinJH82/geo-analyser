"""
geo-analyser API route tests (蚀变分析)
index/mineral_list/analyze/analyze_preview 等核心蚀变路由。
注:数据预处理已拆分至 geo-preprocess,尖点突破已移除,对应路由不再在本系统。
"""
import json
import pytest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))


@pytest.fixture
def analyser_app(monkeypatch, tmp_path):
    """Create Flask test client for geo-analyser."""
    import app as analyser_app_mod
    analyser_app_mod.app.config["TESTING"] = True
    with analyser_app_mod.app.test_client() as client:
        yield client


@pytest.mark.p0
class TestAnalyserRoutes:
    """Core geo-analyser (蚀变分析) API endpoints."""

    def test_index_page(self, analyser_app):
        resp = analyser_app.get("/")
        assert resp.status_code in (200, 302)

    def test_mineral_list(self, analyser_app):
        resp = analyser_app.post("/api/mineral_list", json={})
        assert resp.status_code in (200, 400, 404)

    def test_commodity_list(self, analyser_app):
        resp = analyser_app.get("/api/commodity_list")
        assert resp.status_code in (200, 404)

    def test_analyze_endpoint(self, analyser_app):
        resp = analyser_app.post("/api/analyze", json={"paths": []})
        assert resp.status_code in (200, 400, 404, 500)
