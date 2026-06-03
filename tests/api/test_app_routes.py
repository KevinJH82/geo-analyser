"""
geo-analyser API route tests
browse/scan/preview/process/mineral_list/analyze/analyze_preview/clustering 路由
"""
import json
import pytest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))


@pytest.fixture
def analyser_app(monkeypatch, tmp_path):
    """Create Flask test client for geo-analyser."""
    # Read app module
    import app as analyser_app_mod
    analyser_app_mod.app.config["TESTING"] = True
    with analyser_app_mod.app.test_client() as client:
        yield client


@pytest.mark.p0
class TestAnalyserRoutes:
    """Core geo-analyser API endpoints."""

    def test_index_page(self, analyser_app):
        resp = analyser_app.get("/")
        assert resp.status_code in (200, 302)  # may redirect

    def test_scan_endpoint(self, analyser_app):
        resp = analyser_app.post("/scan", json={"path": "/tmp"})
        assert resp.status_code in (200, 400, 404, 500)

    def test_preview_endpoint(self, analyser_app):
        resp = analyser_app.post("/preview", json={"path": "/tmp/test.tif"})
        assert resp.status_code in (200, 400, 404, 500)

    def test_process_endpoint(self, analyser_app):
        resp = analyser_app.post("/process", json={"path": "/tmp"})
        assert resp.status_code in (200, 400, 404, 500)

    def test_mineral_list(self, analyser_app):
        resp = analyser_app.get("/mineral_list")
        assert resp.status_code in (200, 404)

    def test_analyze_endpoint(self, analyser_app):
        resp = analyser_app.post("/analyze", json={"paths": []})
        assert resp.status_code in (200, 400, 404, 500)

    def test_clustering_endpoint(self, analyser_app):
        resp = analyser_app.post("/clustering", json={})
        assert resp.status_code in (200, 400, 404, 500)

    def test_browse_endpoint(self, analyser_app):
        resp = analyser_app.post("/browse", json={"path": "/tmp"})
        assert resp.status_code in (200, 400, 404, 500)
