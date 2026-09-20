# app.py - FastAPI backend for the Genome Taxonomy Explorer web interface.
# Wraps the existing pipeline (via pipeline_web.py) with subprocess calls,
# it does not reimplement any bioinformatics logic.
import shutil
import uuid
from pathlib import Path
from typing import List

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from job_manager import ALLOWED_EXTENSIONS, MIN_GENOMES, STEP_KEYS, JobManager

BASE_DIR = Path(__file__).resolve().parent.parent      # ~/stage/webapp
JOBS_DIR = BASE_DIR / "jobs"
FRONTEND_DIR = BASE_DIR / "frontend"
JOBS_DIR.mkdir(exist_ok=True)

app = FastAPI(title="Genome Taxonomy Explorer")
manager = JobManager(JOBS_DIR)

# key -> path (relative to a job's directory) exposed for download/display
RESULT_FILES = {
    "qc": ("results/qc.txt", "text/plain"),
    "ani": ("results/ani_resultats.txt", "text/plain"),
    "dddh": ("results/dddh_out/dddh_distance.tsv", "text/tab-separated-values"),
    "aai": ("results/aai_resultats.tsv", "text/tab-separated-values"),
    "pocp": ("results/pocp_out/pocp_matrice.tsv", "text/tab-separated-values"),
    "tree_newick": ("results/arbre.tree", "text/plain"),
    "heatmap_ani": ("results/heatmap_ani.png", "image/png"),
    "heatmap_dddh": ("results/heatmap_dddh.png", "image/png"),
    "heatmap_aai": ("results/heatmap_aai.png", "image/png"),
    "heatmap_pocp": ("results/heatmap_pocp.png", "image/png"),
    "tree_png": ("results/tree.png", "image/png"),
}


def job_dir_or_404(job_id: str) -> Path:
    job_dir = JOBS_DIR / job_id
    if not job_dir.is_dir():
        raise HTTPException(404, "Unknown job id")
    return job_dir


@app.post("/api/jobs")
async def create_job(files: List[UploadFile] = File(...), steps: str = Form("all")):
    if len(files) < MIN_GENOMES:
        raise HTTPException(
            400,
            f"Please upload at least {MIN_GENOMES} genome files. The phylogenomic "
            f"tree step (EasyCGTree) cannot run with fewer. You uploaded {len(files)}.",
        )
    bad = [f.filename for f in files if Path(f.filename).suffix.lower() not in ALLOWED_EXTENSIONS]
    if bad:
        raise HTTPException(
            400,
            f"Unsupported file type(s): {', '.join(bad)}. Allowed extensions: "
            + ", ".join(sorted(ALLOWED_EXTENSIONS)),
        )

    if steps != "all":
        chosen = [s.strip() for s in steps.split(",") if s.strip()]
        unknown = [s for s in chosen if s not in STEP_KEYS]
        if not chosen or unknown:
            raise HTTPException(
                400,
                f"Invalid steps selection: {steps!r}. Choose at least one of: "
                + ", ".join(STEP_KEYS) + " (or 'all').",
            )

    job_id = uuid.uuid4().hex[:10]
    job_dir = JOBS_DIR / job_id
    genomes_dir = job_dir / "genomes"
    genomes_dir.mkdir(parents=True)
    for f in files:
        dest = genomes_dir / Path(f.filename).name
        with open(dest, "wb") as out:
            shutil.copyfileobj(f.file, out)

    manager.submit(job_id, job_dir, len(files), steps)
    return {"job_id": job_id}


@app.post("/api/jobs/{job_id}/cancel")
async def cancel_job(job_id: str):
    state = manager.status(job_id)
    if state is None:
        raise HTTPException(404, "Unknown job id")
    if not manager.cancel(job_id):
        raise HTTPException(409, "Job already finished, nothing to cancel")
    return {"status": "cancelling"}


@app.get("/api/jobs/{job_id}/status")
async def job_status(job_id: str):
    state = manager.status(job_id)
    if state is None:
        raise HTTPException(404, "Unknown job id")
    return state


@app.get("/api/jobs/{job_id}/results")
async def job_results(job_id: str):
    state = manager.status(job_id)
    if state is None:
        raise HTTPException(404, "Unknown job id")
    if state["status"] != "done":
        raise HTTPException(409, "Job is not finished yet")

    job_dir = job_dir_or_404(job_id)
    qc_rows = []
    qc_path = job_dir / "results" / "qc.txt"
    with open(qc_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            values = line.rstrip("\n").split("\t")
            qc_rows.append(dict(zip(header, values)))

    return {
        "qc": qc_rows,
        "n_genomes": state.get("n_genomes"),
        "files": {key: f"/api/jobs/{job_id}/files/{key}" for key in RESULT_FILES},
    }


@app.get("/api/jobs/{job_id}/files/{key}")
async def job_file(job_id: str, key: str):
    if key not in RESULT_FILES:
        raise HTTPException(404, "Unknown result file")
    job_dir = job_dir_or_404(job_id)
    relative_path, media_type = RESULT_FILES[key]
    path = job_dir / relative_path
    if not path.exists():
        raise HTTPException(404, "This result file was not produced")
    return FileResponse(path, media_type=media_type, filename=path.name)


# Serve the static frontend last, so it does not shadow the /api routes above.
app.mount("/", StaticFiles(directory=FRONTEND_DIR, html=True), name="frontend")
