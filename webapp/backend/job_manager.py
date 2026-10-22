# job_manager.py - runs pipeline_web.py as a background subprocess per job,
# one job at a time (a queue), and keeps track of progress/log/status in memory.
import json
import os
import queue
import signal
import subprocess
import sys
import threading
from pathlib import Path

ALLOWED_EXTENSIONS = {".fas", ".fasta", ".fna"}
MIN_GENOMES = 2          # pairwise comparisons need two genomes...
MIN_GENOMES_TREE = 5     # ...and EasyCGTree cannot build a tree with fewer than five
THREADS = os.cpu_count() or 4
STEP_KEYS = ["qc", "checkm2", "proteines", "ani", "dddh", "aai", "pocp", "arbre", "figures"]

SCRIPT = Path(__file__).resolve().parent.parent.parent / "scripts" / "pipeline_web.py"


class JobManager:
    def __init__(self, jobs_dir: Path):
        self.jobs_dir = jobs_dir
        self.lock = threading.Lock()
        self.jobs = {}
        self.running_procs = {}   # job_id -> Popen, only while status == "running"
        self.cancelled = set()    # job_id's for which cancel() was called
        self.queue = queue.Queue()
        self.worker = threading.Thread(target=self._worker_loop, daemon=True)
        self.worker.start()

    def submit(self, job_id: str, job_dir: Path, n_genomes: int, steps: str = "all"):
        with self.lock:
            self.jobs[job_id] = {
                "status": "queued",
                "running_steps": [],   # steps in progress (Snakemake may run several at once)...
                "done_steps": [],      # ...and steps finished, read from the PROGRESS lines of pipeline_web.py
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
            if job is None:
                job = self._from_disk(job_id)
            return dict(job) if job else None

    def _from_disk(self, job_id: str):
        """Jobs live in memory, so a restart forgets them. Their outcome is also written to status.json in the
        job folder (see _mark), so after a restart a finished job keeps its results, a failed one stays failed,
        and one that was running when the server stopped is reported as interrupted."""
        if not job_id.isalnum():
            return None
        job_dir = self.jobs_dir / job_id
        results = job_dir / "results"
        marker = job_dir / "status.json"
        n_genomes = sum(1 for f in (job_dir / "genomes").iterdir()
                        if f.suffix.lower() in ALLOWED_EXTENSIONS) if (job_dir / "genomes").is_dir() else None
        job = {"running_steps": [], "done_steps": [], "log": [], "error": None,
               "n_genomes": n_genomes, "steps": "all"}
        if marker.exists():
            saved = json.loads(marker.read_text())
            job.update(steps=saved.get("steps", "all"), error=saved.get("error"))
            status = saved.get("status")
            if status in ("queued", "running"):
                status = "error"
                job["error"] = "The run was interrupted because the server stopped. Start the analysis again."
            job.update(status=status, label={"done": "Completed", "cancelled": "Cancelled"}.get(status, "Failed"))
        elif results.is_dir() and any(results.iterdir()):
            job.update(status="done", label="Completed")      # job run before status.json existed
        else:
            return None
        self.jobs[job_id] = job
        return job

    def _mark(self, job_dir: Path, status: str, steps: str, error=None):
        """Write the job's outcome next to its results, so it survives a server restart."""
        (job_dir / "status.json").write_text(json.dumps({"status": status, "steps": steps, "error": error}))

    def cancel(self, job_id: str) -> bool:
        """Cancel a queued or running job. Returns False if the job is
        unknown or already finished (nothing to cancel)."""
        with self.lock:
            job = self.jobs.get(job_id)
            if job is None or job["status"] not in ("queued", "running"):
                return False
            self.cancelled.add(job_id)
            proc = self.running_procs.get(job_id)
            if job["status"] == "queued":
                job["status"] = "cancelled"
                job["label"] = "Cancelled"

        if proc is not None:
            # Kill the whole process group: pipeline_web.py itself spawns
            # the real tools (fastANI, CheckM2, ...) as its own children,
            # so terminating just the wrapper would leave them running.
            try:
                os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
            except ProcessLookupError:
                pass  # already exited on its own
        return True

    def _update(self, job_id: str, **kwargs):
        with self.lock:
            if job_id in self.jobs:
                self.jobs[job_id].update(kwargs)

    def _progress(self, job_id: str, event: str, step: str):
        """'PROGRESS start <step>' and 'PROGRESS done <step>': a step can start several times (one job per genome
        for Prodigal and EzAAI), it is in progress until its 'done' line."""
        with self.lock:
            job = self.jobs[job_id]
            if event == "done":
                job["running_steps"] = [s for s in job["running_steps"] if s != step]
                if step not in job["done_steps"]:
                    job["done_steps"] = job["done_steps"] + [step]
            elif step not in job["running_steps"] and step not in job["done_steps"]:
                job["running_steps"] = job["running_steps"] + [step]

    def _append_log(self, job_id: str, line: str):
        with self.lock:
            log = self.jobs[job_id]["log"]
            log.append(line)
            if len(log) > 500:
                del log[0]

    def _worker_loop(self):
        while True:
            job_id, job_dir, steps = self.queue.get()
            with self.lock:
                already_cancelled = job_id in self.cancelled
            if already_cancelled:
                continue  # cancel() already marked it "cancelled" while queued
            self._run_job(job_id, job_dir, steps)

    def _run_job(self, job_id: str, job_dir: Path, steps: str):
        self._update(job_id, status="running", label="Starting...")
        self._mark(job_dir, "running", steps)
        genomes_dir = job_dir / "genomes"
        cmd = [sys.executable, "-u", str(SCRIPT),
               "--genomes-dir", str(genomes_dir),
               "--output-dir", str(job_dir),
               "--job-id", job_id,
               "--threads", str(THREADS),
               "--steps", steps]
        try:
            # start_new_session=True puts the process (and every tool it
            # launches below it: fastANI, CheckM2...) in its own process
            # group, so cancel() can kill the whole tree at once.
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                     text=True, bufsize=1, start_new_session=True)
            with self.lock:
                self.running_procs[job_id] = proc
            for raw_line in proc.stdout:
                line = raw_line.rstrip("\n")
                self._append_log(job_id, line)
                if line.startswith("STEPS "):
                    # the steps that will run: the selected ones, plus those they need (see etapes_prevues)
                    self._update(job_id, steps=line[len("STEPS "):].strip())
                elif line.startswith("PROGRESS ") and len(line.split()) >= 3:
                    self._progress(job_id, *line.split()[1:3])
                elif line.startswith("STEP FAILED"):
                    self._update(job_id, error=line[len("STEP FAILED: "):])
            proc.wait()

            with self.lock:
                self.running_procs.pop(job_id, None)
                was_cancelled = job_id in self.cancelled

            if was_cancelled:
                self._update(job_id, status="cancelled", label="Cancelled")
            elif proc.returncode == 0:
                self._update(job_id, status="done", label="Completed")
            else:
                with self.lock:
                    current_error = self.jobs[job_id].get("error")
                self._update(job_id, status="error", label="Failed",
                             error=current_error or "The pipeline failed. See the log for details.")
        except Exception as exc:
            self._update(job_id, status="error", label="Failed", error=str(exc))
        with self.lock:
            final = dict(self.jobs[job_id])
        self._mark(job_dir, final["status"], final["steps"], final.get("error"))
