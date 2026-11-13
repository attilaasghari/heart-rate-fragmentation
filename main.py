#!/usr/bin/env python3
"""
Heart Rate Fragmentation (HRF) Analyzer — Python/Tkinter rewrite of the
PhysioNet HRF.awk reference implementation.

Reference:
  Costa MD, Davis RB, Goldberger AL. Heart rate fragmentation: a new
  approach to the analysis of cardiac interbeat interval dynamics.
  Front Physiol. 2017;8:255. doi:10.3389/fphys.2017.00255
  Costa M. Software for computing Heart Rate Fragmentation (v1.0.0).
  PhysioNet. 2024. doi:10.13026/0mzj-gn98
"""

import os
import webbrowser          

import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk,
)

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "legend.fontsize": 8,
})


# =====================================================================
#  Core algorithm
# =====================================================================

def _signs(dnn, valid, r):
    """Return sign array (-1, 0, +1) after applying noise threshold r."""
    s = np.zeros(dnn.size, dtype=int)
    s[valid] = np.where(dnn[valid] >  r,  1,
                 np.where(dnn[valid] < -r, -1, 0))
    return s


def compute_hrf(times, labels, r=0.004, min_segment=3,
                count_zero_transitions=True):
    """
    Compute HRF metrics.

    Parameters
    ----------
    times : array-like of float
        R-peak times in seconds.
    labels : array-like of int
        Beat labels: 1 = normal sinus, 2 = PAC, 3 = PVC. Any label != 1
        excludes the beat and its adjacent NN intervals.
    r : float
        Noise tolerance in seconds (default 0.004 s = 1 sample at 250 Hz).
        |dNN| < r is treated as zero.
    min_segment : int
        Minimum dNN count for a "long" segment (default 3, per the paper).
    count_zero_transitions : bool
        If True, a transition between zero and nonzero dNN counts as an
        inflection (matches the PhysioNet description). If False, only
        genuine sign reversals (both nonzero) count.

    Returns
    -------
    dict with metrics and all intermediate arrays needed for plotting.
    """
    times  = np.asarray(times,  dtype=float)
    labels = np.asarray(labels, dtype=int)
    if times.size != labels.size:
        raise ValueError("times and labels must have the same length")
    if times.size < 4:
        raise ValueError("Need at least 4 beats")

    n = times.size

    # 1) Raw NN intervals and validity (both endpoints normal)
    nn       = np.diff(times)                                   # n-1
    valid_nn = (labels[:-1] == 1) & (labels[1:] == 1)           # n-1

    # 2) dNN intervals and validity (both NN intervals valid)
    dnn       = nn[1:] - nn[:-1]                                # n-2
    valid_dnn = valid_nn[1:] & valid_nn[:-1]                    # n-2

    # 3) Threshold
    dnn_thr = dnn.copy()
    dnn_thr[np.abs(dnn_thr) < r] = 0.0

    # 4) Signs and inflection points
    signs = _signs(dnn_thr, valid_dnn, r)
    valid_pairs = valid_dnn[:-1] & valid_dnn[1:]
    sign_diff   = (signs[:-1] != signs[1:])

    if count_zero_transitions:
        infl_mask = valid_pairs & sign_diff
    else:
        infl_mask = (valid_pairs & sign_diff
                     & (signs[:-1] != 0) & (signs[1:] != 0))

    n_inflections = int(np.sum(infl_mask))
    num_nn        = int(np.sum(valid_nn))
    PIP = 100.0 * n_inflections / num_nn if num_nn else float("nan")

    # 5) Segment grouping (for PNNSS / PNNLS)
    #    A segment is a maximal run of consecutive valid dNN with the
    #    same nonzero sign.  Zeroed dNN are neutral (skipped, do not
    #    break a segment).  Invalid dNN break segments.
    seg_info  = []       # list of dicts {start_j, length, sign}
    cur_sign  = 0
    cur_len   = 0
    cur_start = 0
    for j in range(dnn.size):
        if not valid_dnn[j]:
            if cur_len:
                seg_info.append({"start_j": cur_start,
                                 "length": cur_len, "sign": cur_sign})
            cur_sign, cur_len = 0, 0
            continue
        s = signs[j]
        if s == 0:
            continue
        if s == cur_sign:
            cur_len += 1
        else:
            if cur_len:
                seg_info.append({"start_j": cur_start,
                                 "length": cur_len, "sign": cur_sign})
            cur_sign, cur_len, cur_start = s, 1, j
    if cur_len:
        seg_info.append({"start_j": cur_start,
                         "length": cur_len, "sign": cur_sign})

    seg_lengths = [s["length"] for s in seg_info]
    total_seg   = sum(seg_lengths)
    short_dnn   = sum(L for L in seg_lengths if L <  min_segment)
    long_dnn    = sum(L for L in seg_lengths if L >= min_segment)

    if total_seg:
        PNNSS = 100.0 * short_dnn / total_seg
        PNNLS = 100.0 * long_dnn  / total_seg
    else:
        PNNSS = PNNLS = float("nan")

    return {
        "PIP": PIP, "PNNSS": PNNSS, "PNNLS": PNNLS,
        "n_inflections": n_inflections,
        "num_nn": num_nn,
        "num_valid_dnn": int(valid_dnn.sum()),
        "nn": nn, "valid_nn": valid_nn,
        "dnn": dnn, "dnn_thr": dnn_thr, "valid_dnn": valid_dnn,
        "signs": signs, "infl_mask": infl_mask,
        "seg_info": seg_info, "seg_lengths": seg_lengths,
        "r": r, "min_segment": min_segment,
    }


def load_input_file(path):
    """Load a 2-column file: R-peak time (s) and beat label."""
    times, labels = [], []
    with open(path, "r") as fh:
        for ln, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.replace(",", " ").split()
            if len(parts) < 2:
                raise ValueError(
                    f"Line {ln}: expected 2 columns, got {len(parts)}")
            try:
                times.append(float(parts[0]))
                labels.append(int(float(parts[1])))
            except ValueError:
                raise ValueError(f"Line {ln}: cannot parse '{line}'")
    if len(times) < 4:
        raise ValueError("File contains fewer than 4 beats")
    return np.array(times), np.array(labels)


# =====================================================================
#  Tkinter GUI
# =====================================================================

class HRFApp(tk.Tk):

    def __init__(self):
        super().__init__()
        self.title("Heart Rate Fragmentation Analyzer")
        self.geometry("1200x820")
        self.minsize(900, 600)

        self.times = self.labels = None
        self.file_path = None
        self.result    = None

        self._build_menu()
        self._build_layout()

    # ---------- menu -------------------------------------------------
    def _build_menu(self):
        mb = tk.Menu(self)
        fm = tk.Menu(mb, tearoff=0)
        fm.add_command(label="Open…",   command=self.open_file,
                       accelerator="Ctrl+O")
        fm.add_command(label="Export…", command=self.export_results)
        fm.add_separator()
        fm.add_command(label="Quit", command=self.destroy)
        mb.add_cascade(label="File", menu=fm)

        hm = tk.Menu(mb, tearoff=0)
        hm.add_command(label="About…", command=self.show_about)
        mb.add_cascade(label="Help", menu=hm)
        self.config(menu=mb)
        self.bind_all("<Control-o>", lambda e: self.open_file())

    # ---------- layout -----------------------------------------------
    def _build_layout(self):
        top = ttk.Frame(self, padding=8)
        top.pack(side=tk.TOP, fill=tk.X)

        ttk.Label(top, text="File:").grid(row=0, column=0, sticky="w")
        self.file_var = tk.StringVar(value="(none)")
        ttk.Label(top, textvariable=self.file_var, foreground="gray"
                  ).grid(row=0, column=1, columnspan=6, sticky="w",
                         padx=(4, 8))
        ttk.Button(top, text="Open…", command=self.open_file
                   ).grid(row=0, column=7, sticky="e")

        ttk.Label(top, text="r (s):").grid(row=1, column=0,
                                           sticky="w", pady=(6, 0))
        self.r_var = tk.StringVar(value="0.004")
        e_r = ttk.Entry(top, textvariable=self.r_var, width=10)
        e_r.grid(row=1, column=1, sticky="w", padx=(4, 12), pady=(6, 0))
        e_r.bind("<Return>", lambda e: self.compute())

        ttk.Label(top, text="min segment:").grid(row=1, column=2,
                                                 sticky="w", pady=(6, 0))
        self.minseg_var = tk.StringVar(value="3")
        e_m = ttk.Entry(top, textvariable=self.minseg_var, width=6)
        e_m.grid(row=1, column=3, sticky="w", padx=(4, 12), pady=(6, 0))
        e_m.bind("<Return>", lambda e: self.compute())

        self.zero_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            top, text="Count zero↔nonzero transitions",
            variable=self.zero_var
        ).grid(row=1, column=4, columnspan=2, sticky="w", pady=(6, 0))

        ttk.Button(top, text="Compute", command=self.compute
                   ).grid(row=1, column=7, sticky="e", pady=(6, 0))

        nb = ttk.Notebook(self)
        nb.pack(side=tk.TOP, fill=tk.BOTH, expand=True,
                padx=8, pady=(0, 8))

        # ----- Results tab
        rf = ttk.Frame(nb, padding=8)
        nb.add(rf, text="Results")

        box = ttk.LabelFrame(rf, text="HRF metrics", padding=8)
        box.pack(side=tk.TOP, fill=tk.X)
        self.mvars = {k: tk.StringVar(value="—") for k in
                      ("PIP", "PNNSS", "PNNLS", "infl", "nn",
                       "vdnn", "seg")}
        rows = [
            ("PIP (%)",              "PIP"),
            ("PNNSS (%)",            "PNNSS"),
            ("PNNLS (%)",            "PNNLS"),
            ("Inflection points",    "infl"),
            ("Usable NN intervals",  "nn"),
            ("Usable ΔNN intervals", "vdnn"),
            ("Segments",             "seg"),
        ]
        for i, (lbl, key) in enumerate(rows):
            ttk.Label(box, text=lbl + ":").grid(
                row=i, column=0, sticky="w", pady=2)
            ttk.Label(box, textvariable=self.mvars[key],
                      font=("TkFixedFont",)).grid(
                row=i, column=1, sticky="w", padx=(8, 0), pady=2)

        sf = ttk.LabelFrame(rf, text="Summary", padding=8)
        sf.pack(side=tk.TOP, fill=tk.BOTH, expand=True, pady=(8, 0))
        self.summary = tk.Text(sf, height=12, wrap="word")
        sb = ttk.Scrollbar(sf, orient="vertical",
                           command=self.summary.yview)
        self.summary.configure(yscrollcommand=sb.set)
        self.summary.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb.pack(side=tk.RIGHT, fill=tk.Y)

        # ----- Plots tab
        pf = ttk.Frame(nb, padding=4)
        nb.add(pf, text="Plots")
        self.fig    = Figure(figsize=(10, 7.5), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=pf)
        NavigationToolbar2Tk(self.canvas, pf).update()
        self.canvas.get_tk_widget().pack(side=tk.TOP,
                                         fill=tk.BOTH, expand=True)
        self._placeholder()

    def _placeholder(self):
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        ax.text(0.5, 0.5, "Load a file and press Compute.",
                ha="center", va="center", color="gray", fontsize=12)
        ax.set_axis_off()
        self.canvas.draw_idle()

    # ---------- file I/O ---------------------------------------------
    def open_file(self):
        path = filedialog.askopenfilename(
            title="Open R-peak file",
            filetypes=[("Text/CSV", "*.txt *.dat *.csv"),
                       ("All files", "*.*")])
        if not path:
            return
        try:
            self.times, self.labels = load_input_file(path)
        except Exception as exc:
            messagebox.showerror("Load error", str(exc))
            return
        self.file_path = path
        self.file_var.set(f"{os.path.basename(path)} "
                          f"({self.times.size} beats)")
        self.compute()

    def export_results(self):
        if self.result is None:
            messagebox.showinfo("Nothing to export", "Compute first.")
            return
        path = filedialog.asksaveasfilename(
            defaultextension=".tsv",
            filetypes=[("TSV", "*.tsv"), ("Text", "*.txt"),
                       ("All files", "*.*")])
        if not path:
            return
        r = self.result
        with open(path, "w") as fh:
            fh.write("# HRF results\n")
            if self.file_path:
                fh.write(f"# input: {self.file_path}\n")
            fh.write(f"# r = {r['r']} s\n")
            fh.write(f"# min_segment = {r['min_segment']}\n")
            fh.write(f"# zero_transitions = {self.zero_var.get()}\n")
            fh.write(f"PIP\t{r['PIP']:.6f}\n")
            fh.write(f"PNNSS\t{r['PNNSS']:.6f}\n")
            fh.write(f"PNNLS\t{r['PNNLS']:.6f}\n")
            fh.write(f"inflections\t{r['n_inflections']}\n")
            fh.write(f"num_nn\t{r['num_nn']}\n")
            fh.write(f"num_valid_dnn\t{r['num_valid_dnn']}\n")
            fh.write("segments\t" +
                     ",".join(map(str, r["seg_lengths"])) + "\n")

    # ---------- compute ----------------------------------------------
    def compute(self):
        if self.times is None:
            messagebox.showinfo("No data", "Open a file first.")
            return
        try:
            r_val   = float(self.r_var.get())
            min_seg = int(self.minseg_var.get())
        except ValueError:
            messagebox.showerror("Invalid parameter",
                                 "r must be a float, min_segment an int.")
            return
        try:
            self.result = compute_hrf(
                self.times, self.labels,
                r=r_val, min_segment=min_seg,
                count_zero_transitions=self.zero_var.get())
        except Exception as exc:
            messagebox.showerror("Compute error", str(exc))
            return
        self._update_text()
        self._draw_plots()

    def _update_text(self):
        r = self.result
        fmt = lambda x: "—" if (isinstance(x, float) and np.isnan(x)) \
                          else f"{x:.4f}"
        self.mvars["PIP"].set(fmt(r["PIP"]))
        self.mvars["PNNSS"].set(fmt(r["PNNSS"]))
        self.mvars["PNNLS"].set(fmt(r["PNNLS"]))
        self.mvars["infl"].set(str(r["n_inflections"]))
        self.mvars["nn"].set(str(r["num_nn"]))
        self.mvars["vdnn"].set(str(r["num_valid_dnn"]))
        sl = r["seg_lengths"]
        self.mvars["seg"].set(f"{len(sl)} segments; {sl[:15]}"
                              + ("…" if len(sl) > 15 else ""))

        n_non_normal = int(np.sum(self.labels != 1))
        lines = [
            f"Total beats:              {self.times.size}",
            f"Non-normal beats:         {n_non_normal}",
            f"Usable NN intervals:      {r['num_nn']}",
            f"Usable ΔNN intervals:     {r['num_valid_dnn']}",
            "",
            f"PIP   = {fmt(r['PIP'])} %   "
            f"({r['n_inflections']} inflections)",
            f"PNNSS = {fmt(r['PNNSS'])} %   "
            f"(short ΔNN, < {r['min_segment']})",
            f"PNNLS = {fmt(r['PNNLS'])} %   "
            f"(long ΔNN, ≥ {r['min_segment']})",
            "",
            "Interpretation guide:",
            "  • Low PIP + high PNNLS → smooth, vagally mediated rhythm.",
            "  • High PIP + high PNNSS → fragmented rhythm (aging, CAD,",
            "    autonomic dysfunction).",
            "  • HRF values are window-dependent — do not compare across",
            "    very different recording durations.",
            "",
            "Reference: Costa MD, Davis RB, Goldberger AL. "
            "Front Physiol. 2017;8:255.",
        ]
        self.summary.delete("1.0", tk.END)
        self.summary.insert("1.0", "\n".join(lines))

    # ---------- plots ------------------------------------------------
    def _draw_plots(self):
        r = self.result
        nn, dnn_thr, valid_nn, valid_dnn = (
            r["nn"], r["dnn_thr"], r["valid_nn"], r["valid_dnn"])
        r_val   = r["r"]
        min_seg = r["min_segment"]
        signs   = r["signs"]
        infl_m  = r["infl_mask"]

        t_nn  = np.cumsum(nn)                  # end time of each NN
        t_dnn = t_nn[:-1]                      # shared-beat time of dNN

        self.fig.clear()

        # Panel 1 — tachogram
        ax1 = self.fig.add_subplot(3, 1, 1)
        ax1.plot(t_nn[valid_nn],  nn[valid_nn] * 1000,
                 "o-", ms=3, lw=0.8, color="tab:blue",
                 label="valid NN")
        if (~valid_nn).any():
            ax1.plot(t_nn[~valid_nn], nn[~valid_nn] * 1000,
                     "x", ms=6, color="tab:red", label="excluded")
        ax1.set_ylabel("NN (ms)")
        ax1.set_title("Tachogram")
        ax1.legend(loc="upper right")
        ax1.grid(alpha=0.3)

        # Panel 2 — dNN, threshold band, inflections
        ax2 = self.fig.add_subplot(3, 1, 2, sharex=ax1)
        ax2.axhspan(-r_val * 1000, r_val * 1000, color="gray",
                    alpha=0.15, label=f"±r = ±{r_val*1000:.1f} ms")
        ax2.plot(t_dnn[valid_dnn], dnn_thr[valid_dnn] * 1000,
                 "-", lw=0.5, color="0.6")
        pos = valid_dnn & (dnn_thr >  0)
        neg = valid_dnn & (dnn_thr <  0)
        ax2.plot(t_dnn[pos], dnn_thr[pos] * 1000, ".",
                 color="tab:green", ms=6)
        ax2.plot(t_dnn[neg], dnn_thr[neg] * 1000, ".",
                 color="tab:red",   ms=6)
        if infl_m.any():
            infl_t = (t_dnn[:-1] + t_dnn[1:]) / 2.0
            ax2.plot(infl_t[infl_m], np.zeros(infl_m.sum()),
                     "v", color="magenta", ms=6,
                     label=f"inflection (n={infl_m.sum()})")
        ax2.set_ylabel("ΔNN (ms)")
        ax2.set_title("ΔNN sequence, threshold band, inflections")
        ax2.legend(loc="upper right")
        ax2.grid(alpha=0.3)

        # Panel 3 — segment bars
        ax3 = self.fig.add_subplot(3, 1, 3, sharex=ax1)
        seen = {"short": False, "long": False}
        for seg in r["seg_info"]:
            is_short = seg["length"] < min_seg
            key  = "short" if is_short else "long"
            col  = "tab:orange" if is_short else "tab:blue"
            lbl  = (f"short (<{min_seg})" if is_short
                    else f"long (≥{min_seg})")
            ax3.bar(t_dnn[seg["start_j"]], seg["length"],
                    width=2.0, color=col, alpha=0.8,
                    label=lbl if not seen[key] else None)
            seen[key] = True
        ax3.set_ylabel("segment length")
        ax3.set_xlabel("time (s)")
        ax3.set_title("Accelerative / decelerative segments")
        ax3.legend(loc="upper right")
        ax3.grid(alpha=0.3, axis="y")

        self.fig.tight_layout()
        self.canvas.draw_idle()

        # ---------- about ------------------------------------------------
    def show_about(self):
        """Custom About window with clickable author / reference links."""
        win = tk.Toplevel(self)
        win.title("About — Heart Rate Fragmentation Analyzer")
        win.resizable(False, False)
        win.transient(self)
        win.grab_set()

        pad = ttk.Frame(win, padding=16)
        pad.pack(fill=tk.BOTH, expand=True)

        ttk.Label(pad, text="Heart Rate Fragmentation Analyzer",
                  font=("TkDefaultFont", 13, "bold")).pack(anchor="w")
        ttk.Label(pad, text="Python / Tkinter port of the PhysioNet HRF.awk reference software",
                  foreground="gray").pack(anchor="w", pady=(0, 10))

        ttk.Separator(pad, orient="horizontal").pack(fill=tk.X, pady=8)

        ttk.Label(pad, text="Created by",
                  font=("TkDefaultFont", 10, "bold")).pack(anchor="w")
        for name, url in (
            ("Asal Karimi",    "https://asal.archtek.ir"),
            ("Attila Asghari", "https://attila.archtek.ir"),
            ("Archtek",        "https://archtek.ir"),
        ):
            row = ttk.Frame(pad)
            row.pack(anchor="w", fill=tk.X, pady=1)
            ttk.Label(row, text=f"  {name}  —  ").pack(side=tk.LEFT)
            link = tk.Label(row, text=url, fg="#1a73e8", cursor="hand2",
                            font=("TkDefaultFont", 9, "underline"))
            link.pack(side=tk.LEFT)
            link.bind("<Button-1>", lambda e, u=url: self._open_url(u))

        ttk.Separator(pad, orient="horizontal").pack(fill=tk.X, pady=8)

        ttk.Label(pad, text="References",
                  font=("TkDefaultFont", 10, "bold")).pack(anchor="w")
        tk.Label(
            pad, justify="left", wraplength=520,
            text=("Method: Costa MD, Davis RB, Goldberger AL. "
                  "Heart rate fragmentation: a new approach to the analysis "
                  "of cardiac interbeat interval dynamics. "
                  "Front Physiol. 2017;8:255.\n\n"
                  "Software: Costa M. Software for computing Heart Rate "
                  "Fragmentation (v1.0.0). PhysioNet. 2024. "
                  "doi:10.13026/0mzj-gn98"),
        ).pack(anchor="w")

        ttk.Separator(pad, orient="horizontal").pack(fill=tk.X, pady=8)
        tk.Label(
            pad, justify="left", wraplength=520, fg="#a04000",
            text=("NOTE: HRF values are window-dependent. "
                  "Validate against the PhysioNet reference test files "
                  "before any clinical use."),
        ).pack(anchor="w", pady=(0, 8))

        ttk.Button(pad, text="Close", command=win.destroy
                   ).pack(anchor="e", pady=(8, 0))

    def _open_url(self, url):
        """Open a URL in the system's default browser."""
        try:
            webbrowser.open(url, new=2)
        except Exception:
            pass


def main():
    HRFApp().mainloop()


if __name__ == "__main__":
    main()
