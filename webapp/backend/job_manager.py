# job_manager.py - runs pipeline_web.py as a background subprocess per job,
# one job at a time (a queue), and keeps track of progress/log/status in memory.
import queue
import subprocess
import sys
import threading
from pathlib import Path

ALLOWED_EXTENSIONS = {".fas", ".fasta", ".fna"}
MIN_GENOMES = 5
TOTAL_STEPS = 7
STEP_KEYS = ["qc", "proteines", "ani", "aai", "pocp", "arbre", "figures"]

SCRIPT = Path(__file__).resolve().parent.parent.parent / "scripts" / "pipeline_web.py"


class JobManager:
    def __init__(self, jobs_dir: Path):
        self.jobs_dir = jobs_dir
        self.lock = threading.Lock()
        self.jobs = {}
        self.queue = queue.Queue()
        self.worker = threading.Thread(target=self._worker_loop, daemon=True)
        self.worker.start()

    def submit(self, job_id: str, job_dir: Path, n_genomes: int, steps: str = "all"):
        with self.lock:
            self.jobs[job_id] = {
                "status": "queued",
                "step": 0,
                "total": TOTAL_STEPS,
                "label": "Waiting in queue...",
                "log": [],
                "error": None,
                "n_genomes": n_genomes,
                "steps": steps,
            }
        self.queue.put((job_id, job_dir, steps))

    def status(self, job_id: str):
        with self.lock:
            job = self.jobs.get(job_id)
            return dict(job) if job else None

    def _update(self, job_id: str, **kwargs):
        with self.lock:
            if job_id in self.jobs:
                self.jobs[job_id].update(kwargs)

    def _append_log(self, job_id: str, line: str):
        with self.lock:
            log = self.jobs[job_id]["log"]
            log.append(line)
            if len(log) > 500:
                del log[0]

    def _worker_loop(self):
        while True:
            job_id, job_dir, steps = self.queue.get()
            self._run_job(job_id, job_dir, steps)

    def _run_job(self, job_id: str, job_dir: Path, steps: str):
        self._update(job_id, status="running", label="Starting...")
        genomes_dir = job_dir / "genomes"
        cmd = [sys.executable, "-u", str(SCRIPT),
               "--genomes-dir", str(genomes_dir),
               "--output-dir", str(job_dir),
               "--job-id", job_id,
               "--threads", "4",
               "--steps", steps]
        try:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                     text=True, bufsize=1)
            for raw_line in proc.stdout:
                line = raw_line.rstrip("\n")
                self._append_log(job_id, line)
                if line.startswith("PROGRESS"):
                    parts = line.split(" ", 2)
                    n, total = parts[1].split("/")
                    label = parts[2] if len(parts) > 2 else ""
                    self._update(job_id, step=int(n), total=int(total), label=label)
                elif line.startswith("STEP FAILED"):
                    self._update(job_id, error=line[len("STEP FAILED: "):])
            proc.wait()
            if proc.returncode == 0:
                self._update(job_id, status="done", label="Completed")
            else:
                with self.lock:
                    current_error = self.jobs[job_id].get("error")
                self._update(job_id, status="error", label="Failed",
                             error=current_error or "The pipeline failed. See the log for details.")
        except Exception as exc:
            self._update(job_id, status="error", label="Failed", error=str(exc))
