from __future__ import annotations

import shutil
import tempfile
import uuid
from importlib.resources import files
from pathlib import Path

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse

from .models import AnalysisConfig
from .pipeline import run_analysis


ROOT = Path(__file__).resolve().parents[2]
RESULTS_ROOT = ROOT / "results"
REFERENCE = ROOT / "references" / "releases" / "nu3gen_animals_v1.fasta"
MAX_UPLOAD_BYTES = 2 * 1024 * 1024 * 1024

app = FastAPI(title="Nanopore Metabarcoding", version="0.1.0")


@app.get("/", response_class=HTMLResponse)
def index():
    return files("metabarcoder").joinpath("web/index.html").read_text(encoding="utf-8")


@app.post("/api/analyze")
def analyze(sample_name: str = Form("sample"), uploads: list[UploadFile] = File(...)):
    job_id = uuid.uuid4().hex[:12]
    result_dir = RESULTS_ROOT / job_id
    result_dir.mkdir(parents=True, exist_ok=False)
    try:
        with tempfile.TemporaryDirectory(prefix="metabarcoder-") as tmp:
            input_paths: list[Path] = []
            total = 0
            for upload in uploads:
                safe_name = Path(upload.filename or "upload.fastq").name
                if not (safe_name.lower().endswith(".fastq") or safe_name.lower().endswith(".fastq.gz")):
                    raise HTTPException(400, f"Nepodporovaný súbor: {safe_name}")
                target = Path(tmp) / f"{len(input_paths):04d}_{safe_name}"
                with target.open("wb") as handle:
                    while chunk := upload.file.read(1024 * 1024):
                        total += len(chunk)
                        if total > MAX_UPLOAD_BYTES:
                            raise HTTPException(413, "Upload prekročil limit 2 GB")
                        handle.write(chunk)
                input_paths.append(target)
            result = run_analysis(AnalysisConfig(
                sample_name=sample_name.strip() or "sample", inputs=input_paths,
                reference_path=REFERENCE, output_dir=result_dir,
            ))
        return JSONResponse({
            "job_id": job_id, "sample": result["sample"],
            "primer_detection": result["primer_detection"],
            "taxonomy": result["taxonomy"],
            "report_url": f"/api/jobs/{job_id}/report",
            "json_url": f"/api/jobs/{job_id}/json",
            "consensus_url": f"/api/jobs/{job_id}/consensus",
        })
    except HTTPException:
        shutil.rmtree(result_dir, ignore_errors=True)
        raise
    except Exception as exc:
        shutil.rmtree(result_dir, ignore_errors=True)
        raise HTTPException(422, str(exc)) from exc


@app.get("/api/jobs/{job_id}/report")
def report(job_id: str):
    path = RESULTS_ROOT / Path(job_id).name / "report.html"
    if not path.exists():
        raise HTTPException(404, "Report neexistuje")
    return FileResponse(path, media_type="text/html")


@app.get("/api/jobs/{job_id}/json")
def report_json(job_id: str):
    path = RESULTS_ROOT / Path(job_id).name / "report.json"
    if not path.exists():
        raise HTTPException(404, "Výsledok neexistuje")
    return FileResponse(path, media_type="application/json", filename="report.json")


@app.get("/api/jobs/{job_id}/consensus")
def consensus(job_id: str):
    path = RESULTS_ROOT / Path(job_id).name / "consensus.fasta"
    if not path.exists():
        raise HTTPException(404, "Konsenzuálna sekvencia neexistuje")
    return FileResponse(path, media_type="text/plain", filename="consensus.fasta")

