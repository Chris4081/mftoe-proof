#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MFToE GUI (v2, clean) – simple Tkinter cockpit for common tasks:
Tabs: Info, Scan Relaxion, Run Relaxion, Run Baselines, BAO Quickcheck, Cov Sweep, Run All

- No dark/light mode.
- No fetch-DESI tab (as requested).
- Runs in threads; live log in a bottom console.
- Assumes you run this from the repo root (or places itself there).
"""

import os
import sys
import threading
import subprocess
import queue
import shlex
import webbrowser
from pathlib import Path
import tkinter as tk
from tkinter import ttk, messagebox, filedialog

# ---------- Repo root detection ----------
HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE  # file should live in the repo root
os.chdir(REPO_ROOT)

# ---------- Small helpers ----------
def which_python() -> str:
    """Return the python executable to use."""
    return sys.executable or "python3"

def open_path(p: Path):
    try:
        if sys.platform.startswith("darwin"):
            subprocess.run(["open", str(p)], check=False)
        elif os.name == "nt":
            os.startfile(str(p))  # type: ignore
        else:
            subprocess.run(["xdg-open", str(p)], check=False)
    except Exception as e:
        messagebox.showerror("Open path failed", f"{p}\n\n{e}")

# ---------- Command runner with live log ----------
class CmdRunner:
    def __init__(self, text_widget: tk.Text, start_btn: ttk.Button, stop_btn: ttk.Button | None = None):
        self.text = text_widget
        self.start_btn = start_btn
        self.stop_btn = stop_btn
        self.proc: subprocess.Popen | None = None
        self.th: threading.Thread | None = None
        self.q = queue.Queue()

    def log(self, s: str):
        self.text.configure(state="normal")
        self.text.insert("end", s)
        self.text.see("end")
        self.text.configure(state="disabled")

    def _enqueue_output(self, pipe, tag):
        try:
            for line in iter(pipe.readline, b""):
                self.q.put((tag, line.decode(errors="replace")))
        finally:
            pipe.close()

    def _drain_queue(self):
        try:
            while True:
                tag, line = self.q.get_nowait()
                self.log(line)
        except queue.Empty:
            pass
        if self.proc and self.proc.poll() is None:
            # still running
            self.text.after(100, self._drain_queue)
        else:
            # finished
            rc = self.proc.returncode if self.proc else None
            self.log(f"\n— process finished (returncode={rc}) —\n")
            self.start_btn.configure(state="normal")
            if self.stop_btn:
                self.stop_btn.configure(state="disabled")
            self.proc = None

    def run(self, cmd: list[str] | str, cwd: Path | None = None, env: dict | None = None, shell: bool = False):
        if self.proc:
            messagebox.showwarning("Busy", "Another process is still running.")
            return
        self.start_btn.configure(state="disabled")
        if self.stop_btn:
            self.stop_btn.configure(state="normal")

        self.log(f"\n$ {' '.join(cmd) if isinstance(cmd, (list, tuple)) else cmd}\n")

        try:
            self.proc = subprocess.Popen(
                cmd if not shell else " ".join(shlex.quote(c) for c in cmd),  # noqa
                cwd=str(cwd) if cwd else None,
                env=env or os.environ.copy(),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=False,
                shell=shell
            )
        except FileNotFoundError as e:
            self.log(f"ERROR: {e}\n")
            self.start_btn.configure(state="normal")
            if self.stop_btn:
                self.stop_btn.configure(state="disabled")
            self.proc = None
            return

        # threads to read stdout/stderr
        t_out = threading.Thread(target=self._enqueue_output, args=(self.proc.stdout, "OUT"), daemon=True)  # type: ignore
        t_err = threading.Thread(target=self._enqueue_output, args=(self.proc.stderr, "ERR"), daemon=True)  # type: ignore
        t_out.start(); t_err.start()
        self._drain_queue()

    def stop(self):
        if not self.proc:
            return
        try:
            self.proc.terminate()
        except Exception:
            pass

# ---------- GUI ----------
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("MFToE Proof – GUI v2")
        self.geometry("1024x720")

        # Notebook
        self.nb = ttk.Notebook(self)
        self.nb.pack(fill="both", expand=True)

        # bottom console
        self.console = tk.Text(self, height=14, wrap="word", state="disabled")
        self.console.pack(fill="both", side="bottom")
        self.console.configure(font=("Menlo", 11) if sys.platform == "darwin" else ("Courier New", 10))

        # runner control
        self.start_btn_dummy = ttk.Button(self, text="")  # placeholder for binding
        self.runner = CmdRunner(self.console, self.start_btn_dummy)

        # tabs
        self._build_info()
        self._build_scan_relaxion()
        self._build_run_relaxion()
        self._build_run_baselines()
        self._build_bao_quickcheck()
        self._build_cov_sweep()
        self._build_run_all()

    # --- Info tab ---
    def _build_info(self):
        frm = ttk.Frame(self.nb)
        self.nb.add(frm, text="Info")

        title = ttk.Label(frm, text="MFToE Proof — Dark Energy Reconstruction with DESI DR2 (2025)", font=("Helvetica", 14, "bold"))
        title.pack(anchor="w", padx=12, pady=(12, 6))

        descr = (
            "This GUI is a small cockpit to launch common MFToE tasks.\n"
            "Each tab includes a short description and a Start button.\n\n"
            "Docs (PDF):\n"
            " • docs/mftoe.pdf\n"
            " • docs/mftoe_proof.pdf\n\n"
            "Licensing & Data:\n"
            " • Code: AGPL-3.0\n"
            " • DESI DR2 data © DESI Collaboration (2025), CC BY 4.0\n"
            "   See: https://data.desi.lbl.gov/doc/releases/\n"
        )
        ttk.Label(frm, text=descr, justify="left").pack(anchor="w", padx=12)

        link_frame = ttk.Frame(frm)
        link_frame.pack(anchor="w", padx=12, pady=8)

        def link_btn(label, path):
            b = ttk.Button(link_frame, text=label, command=lambda: open_path(REPO_ROOT / path))
            b.pack(side="left", padx=(0, 10))

        link_btn("Open mftoe.pdf", "docs/mftoe.pdf")
        link_btn("Open mftoe_proof.pdf", "docs/mftoe_proof.pdf")
        ttk.Button(link_frame, text="Open runs/ folder", command=lambda: open_path(REPO_ROOT / "runs")).pack(side="left")

    # --- Scan Relaxion tab (scripts/scan_relaxion.sh) ---
    def _build_scan_relaxion(self):
        frm = ttk.Frame(self.nb)
        self.nb.add(frm, text="Scan Relaxion")

        ttk.Label(frm, text="Runs a small (gamma, sigma_noise) grid via scripts/scan_relaxion.sh", font=("Helvetica", 11, "bold")).pack(anchor="w", padx=12, pady=(12, 6))
        ttk.Label(frm, text="Outputs summary CSV + PNGs into runs/.").pack(anchor="w", padx=12, pady=(0, 12))

        row = ttk.Frame(frm); row.pack(anchor="w", padx=12, pady=(4, 6))
        start = ttk.Button(row, text="Start Scan")
        start.pack(side="left")
        self.runner_scan = CmdRunner(self.console, start)

        def run_scan():
            script = REPO_ROOT / "scripts" / "scan_relaxion.sh"
            if not script.exists():
                messagebox.showerror("Missing script", f"{script} not found.")
                return
            self.runner_scan.run(["bash", str(script)], cwd=REPO_ROOT)

        start.configure(command=run_scan)

        ttk.Button(frm, text="Open runs/", command=lambda: open_path(REPO_ROOT / "runs")).pack(anchor="w", padx=12, pady=(6, 6))

    # --- Run Relaxion (single example) ---
    def _build_run_relaxion(self):
        frm = ttk.Frame(self.nb)
        self.nb.add(frm, text="Run Relaxion")

        ttk.Label(frm, text="Single MFToE relaxion run (with RG+noise); saves CSV/PNG in runs/.", font=("Helvetica", 11, "bold")).pack(anchor="w", padx=12, pady=(12, 6))

        # Params file (optional)
        box = ttk.Labelframe(frm, text="Parameters (optional JSON overrides)")
        box.pack(fill="x", padx=12, pady=(0, 8))
        pth_var = tk.StringVar(value="")
        e = ttk.Entry(box, textvariable=pth_var, width=70)
        e.pack(side="left", padx=(8, 6), pady=8)
        def choose():
            p = filedialog.askopenfilename(initialdir=str(REPO_ROOT), filetypes=[("JSON", "*.json"), ("All", "*.*")])
            if p:
                pth_var.set(p)
        ttk.Button(box, text="Browse…", command=choose).pack(side="left", padx=6)

        # Start
        row = ttk.Frame(frm); row.pack(anchor="w", padx=12, pady=(4, 6))
        start = ttk.Button(row, text="Run Relaxion")
        start.pack(side="left")
        self.runner_rel = CmdRunner(self.console, start)

        def run_relaxion():
            py = which_python()
            cmd = [
                py, "mftoe_vacuum_astropy.py",
                "--mode", "relaxion",
                "--rg", "on",
                "--noise", "on",
                "--astropy", "on",
                "--desi-bestfit", "data/desi_dr2/iminuit/base/desi-bao-all/bestfit.minimum",
                "--out", "runs/relaxion_gui"
            ]
            if pth_var.get().strip():
                cmd += ["--params", pth_var.get().strip()]
            self.runner_rel.run(cmd, cwd=REPO_ROOT)

        start.configure(command=run_relaxion)

    # --- Run Baselines (scripts/run_baselines.sh) ---
    def _build_run_baselines(self):
        frm = ttk.Frame(self.nb)
        self.nb.add(frm, text="Run Baselines")

        ttk.Label(frm, text="Runs repository baseline pipeline via scripts/run_baselines.sh", font=("Helvetica", 11, "bold")).pack(anchor="w", padx=12, pady=(12, 6))
        ttk.Label(frm, text="Generates baseline CSV/PNG and BAO check in runs/.").pack(anchor="w", padx=12)

        row = ttk.Frame(frm); row.pack(anchor="w", padx=12, pady=(8, 6))
        start = ttk.Button(row, text="Start Baselines")
        start.pack(side="left")
        self.runner_base = CmdRunner(self.console, start)

        def run_base():
            script = REPO_ROOT / "scripts" / "run_baselines.sh"
            if not script.exists():
                messagebox.showerror("Missing script", f"{script} not found.")
                return
            self.runner_base.run(["bash", str(script)], cwd=REPO_ROOT)

        start.configure(command=run_base)

    # --- BAO Quickcheck (analysis/bao_compare.py) ---
    def _build_bao_quickcheck(self):
        frm = ttk.Frame(self.nb)
        self.nb.add(frm, text="BAO Quickcheck")

        ttk.Label(frm, text="Run a quick BAO χ² comparison with adjustable inputs.", font=("Helvetica", 11, "bold")).pack(anchor="w", padx=12, pady=(12, 6))

        # Inputs with sensible defaults
        defaults = {
            "model_csv": REPO_ROOT / "runs" / "mftoe_vacuum_astropy.csv",
            "bao_csv":   REPO_ROOT / "data" / "desi_dr2" / "bao_summary.csv",
            "H0phys":    "67.36",
            "rd":        "150.754",
            "out":       REPO_ROOT / "runs" / "bao_quickcheck_gui"
        }

        def row_entry(parent, label, var, browse=False, is_file=True):
            fr = ttk.Frame(parent); fr.pack(fill="x", padx=12, pady=4)
            ttk.Label(fr, text=label, width=18).pack(side="left")
            ent = ttk.Entry(fr, textvariable=var, width=70); ent.pack(side="left", padx=(6, 6))
            if browse:
                def pick():
                    if is_file:
                        p = filedialog.askopenfilename(initialdir=str(REPO_ROOT))
                    else:
                        p = filedialog.asksaveasfilename(initialdir=str(REPO_ROOT))
                    if p:
                        var.set(p)
                ttk.Button(fr, text="Browse…", command=pick).pack(side="left")
            return ent

        var_model = tk.StringVar(value=str(defaults["model_csv"]))
        var_bao   = tk.StringVar(value=str(defaults["bao_csv"]))
        var_H0    = tk.StringVar(value=defaults["H0phys"])
        var_rd    = tk.StringVar(value=defaults["rd"])
        var_out   = tk.StringVar(value=str(defaults["out"]))

        row_entry(frm, "Model CSV", var_model, browse=True, is_file=True)
        row_entry(frm, "BAO CSV",   var_bao,   browse=True, is_file=True)
        row_entry(frm, "H0phys",    var_H0)
        row_entry(frm, "r_d [Mpc]", var_rd)
        row_entry(frm, "Output base", var_out, browse=True, is_file=False)

        # Start
        row = ttk.Frame(frm); row.pack(anchor="w", padx=12, pady=(8, 6))
        start = ttk.Button(row, text="Run Quickcheck")
        start.pack(side="left")
        self.runner_bao = CmdRunner(self.console, start)

        def run_bao():
            py = which_python()
            model_csv = var_model.get().strip()
            bao_csv = var_bao.get().strip()
            H0 = var_H0.get().strip()
            rd = var_rd.get().strip()
            out = var_out.get().strip()

            if not Path(model_csv).exists():
                messagebox.showerror("Missing file", f"Model CSV not found:\n{model_csv}")
                return
            if not Path(bao_csv).exists():
                messagebox.showerror("Missing file", f"BAO CSV not found:\n{bao_csv}")
                return

            cmd = [
                py, "analysis/bao_compare.py",
                "--model-csv", model_csv,
                "--bao-csv", bao_csv,
                "--H0phys", H0,
                "--rd", rd,
                "--out", out
            ]
            self.runner_bao.run(cmd, cwd=REPO_ROOT)

        start.configure(command=run_bao)

    # --- Cov Sweep (scripts/cov_sweep.sh) ---
    def _build_cov_sweep(self):
        frm = ttk.Frame(self.nb)
        self.nb.add(frm, text="Cov Sweep")

        ttk.Label(frm, text="Run synthetic covariance rho-sweep via scripts/cov_sweep.sh", font=("Helvetica", 11, "bold")).pack(anchor="w", padx=12, pady=(12, 6))
        ttk.Label(frm, text="Generates covariance and BAO χ² for several ρ values.").pack(anchor="w", padx=12)

        row = ttk.Frame(frm); row.pack(anchor="w", padx=12, pady=(8, 6))
        start = ttk.Button(row, text="Start Cov Sweep")
        start.pack(side="left")
        self.runner_cov = CmdRunner(self.console, start)

        def run_cov():
            script = REPO_ROOT / "scripts" / "cov_sweep.sh"
            if not script.exists():
                messagebox.showerror("Missing script", f"{script} not found.")
                return
            self.runner_cov.run(["bash", str(script)], cwd=REPO_ROOT)

        start.configure(command=run_cov)

    # --- Run All (baseline → scan → cov sweep) ---
    def _build_run_all(self):
        frm = ttk.Frame(self.nb)
        self.nb.add(frm, text="Run All")

        ttk.Label(frm, text="Runs baseline, relaxion scan and covariance sweep sequentially.", font=("Helvetica", 11, "bold")).pack(anchor="w", padx=12, pady=(12, 6))
        ttk.Label(frm, text="Use this for a quick end-to-end reproducibility check (writes to runs/).").pack(anchor="w", padx=12)

        row = ttk.Frame(frm); row.pack(anchor="w", padx=12, pady=(8, 6))
        start = ttk.Button(row, text="Start Run-All")
        start.pack(side="left")
        stop = ttk.Button(row, text="Stop", state="disabled")
        stop.pack(side="left", padx=(8, 0))

        self.runner_all = CmdRunner(self.console, start, stop_btn=stop)

        def run_all():
            bash = "bash"
            scripts = [
                REPO_ROOT / "scripts" / "run_baselines.sh",
                REPO_ROOT / "scripts" / "scan_relaxion.sh",
                REPO_ROOT / "scripts" / "cov_sweep.sh",
            ]
            for s in scripts:
                if not s.exists():
                    messagebox.showerror("Missing script", f"{s} not found.")
                    return

            # chain them in one shell call so we keep a single process/log stream
            chain = " && ".join(shlex.quote(str(s)) for s in scripts)
            cmd = [bash, "-lc", chain]
            self.runner_all.run(cmd, cwd=REPO_ROOT, shell=False)

        def stop_all():
            self.runner_all.stop()

        start.configure(command=run_all)
        stop.configure(command=stop_all)

# ---------- main ----------
if __name__ == "__main__":
    App().mainloop()