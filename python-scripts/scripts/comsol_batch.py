#!/usr/bin/env python3
"""Solve one or more studies of a COMSOL template with `comsol batch` (no GUI,
no MPh, no LiveLink) and write the eigenfrequencies to CSV / JSON.

Normally started through scripts/run_comsol_batch.sh; see that file for the
usage examples and for which licence features a batch run checks out.

Steps (each is one `comsolbatch` process; a seat is held only while it runs):
  1. solve     comsolbatch -inputfile <template> -outputfile <out>/solved_<study>.mph
                           -study <tag> [-pname a,w -plist 600[nm],1700[nm]]
                           -batchlog <out>/solve_<study>.log [-usebatchlic]
  2. extract   a small Java class, generated and compiled into <out>/, is run
               once with comsolbatch on all the solved files. It evaluates
               `freq` (and, for the optical studies, energy integrals) on the
               dataset of each study and writes <out>/eig_<study>.csv.
               This step cannot use the batch licence pool (COMSOL: "Input
               file does not support batch licenses" for a class file), so it
               runs on the ordinary features for its few seconds. If it fails
               for want of a seat the solved .mph is kept: rerun with
               --extract-only --outdir <out>.
  3. summary   <out>/summary.json plus a table on stdout: Re/Im of every
               eigenfrequency and whether it lies below the light line
               c*kF/(2*pi) (= c/2a at the X point).

Why a Java class for step 2: the batch log does not list eigenvalues, the
solution inside an .mph is an undocumented binary, and MPh would start a
second, interactive COMSOL session (COMSOL + module seat) just to read
numbers. The class uses the same documented API calls MPh makes
(result().numerical().create(.., "EvalGlobal") + computeResult()), holds its
seats for seconds instead of a session, and is compiled with the javac that
ships with COMSOL, which needs no licence.

Only the standard library is used, so any Python 3.8+ works: WSL/Linux
python3 (Windows COMSOL is called through WSL interop) or Windows python.
"""
from __future__ import annotations

import argparse
import csv
import datetime as _dt
import glob
import json
import math
import os
import re
import subprocess
import sys
import threading
import time
import zipfile

C0 = 299_792_458.0
F_TARGET = C0 / 1550e-9                 # 193.4 THz, as src/comsol_optics.py
PSEUDO_MODE_HZ = 10e12                  # below this: scattering-BC pseudo-modes

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DEFAULT_TEMPLATE = os.path.join(ROOT, "comsol", "omc_unitcell_rib.mph")

# Study labels of comsol/build_omc_unitcell_rib.m -> tags (-study wants a tag).
STUDY_ALIASES = {"mech sym": "std1", "mech antisym": "std3",
                 "opt te": "std2", "opt tm": "std4",
                 "te": "std2", "tm": "std4"}
STUDY_LABELS = {"std1": "mech sym", "std3": "mech antisym",
                "std2": "opt TE", "std4": "opt TM"}
OPTICAL_STUDIES = ("std2", "std4")

# Volume integrals evaluated per mode for the optical studies:
# (column, expression, named selection or "" for all domains).
# geom1_ribcell_dom is the beam domain of the rib template.
OPTICAL_INTEGRALS = [
    ("We_all",  "ewfd.Weav",       ""),
    ("We_beam", "ewfd.Weav",       "geom1_ribcell_dom"),
    ("Ex2_beam", "abs(ewfd.Ex)^2", "geom1_ribcell_dom"),
    ("Ey2_beam", "abs(ewfd.Ey)^2", "geom1_ribcell_dom"),
    ("Ez2_beam", "abs(ewfd.Ez)^2", "geom1_ribcell_dom"),
]

COMSOL_ROOT_GLOBS = [
    "/mnt/c/Program Files/COMSOL/COMSOL*/Multiphysics",
    "/c/Program Files/COMSOL/COMSOL*/Multiphysics",
    "C:/Program Files/COMSOL/COMSOL*/Multiphysics",
    "/usr/local/comsol*/multiphysics",
    "/opt/comsol*/multiphysics",
    "/Applications/COMSOL*/Multiphysics",
]

LICENCE_ERROR_RE = re.compile(
    r"licen[sc]e error|licensed number of users already reached|"
    r"flexnet licensing error|\(-4,\d+|could not (obtain|check ?out) licen[sc]e",
    re.I)


def log(msg):
    print(f"[{_dt.datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", file=sys.stderr,
          flush=True)


class BatchError(RuntimeError):
    def __init__(self, msg, code=1):
        super().__init__(msg)
        self.code = code


# ---- paths -----------------------------------------------------------------

def local_path(p):
    """Accept /mnt/c/... and /c/... spellings when running under Windows python."""
    if os.name == "nt":
        m = re.match(r"^/(?:mnt/)?([A-Za-z])/(.*)$", p.replace("\\", "/"))
        if m:
            return f"{m.group(1).upper()}:/{m.group(2)}"
    return p


def exe_path(p, windows_exe):
    """Absolute path in the form the COMSOL executable understands.

    A Windows executable started from WSL needs C:\\... instead of /mnt/c/...
    Forward slashes are kept out of it (COMSOL takes both, cmd tools do not).
    """
    p = os.path.abspath(local_path(p))
    if os.name == "nt" or not windows_exe:
        return p
    m = re.match(r"^/mnt/([A-Za-z])/(.*)$", p)
    if m:
        return f"{m.group(1).upper()}:\\" + m.group(2).replace("/", "\\")
    m = re.match(r"^/([A-Za-z])/(.*)$", p)            # Git Bash spelling
    if m and os.path.isdir(f"/{m.group(1)}/"):
        return f"{m.group(1).upper()}:\\" + m.group(2).replace("/", "\\")
    raise BatchError(
        f"{p} is not on a Windows drive, so Windows COMSOL cannot read or "
        "write it. Use a path under /mnt/<drive>/.")


def find_comsol_root(arg):
    cand = arg or os.environ.get("COMSOL_ROOT")
    if cand:
        cand = local_path(cand)
        if not os.path.isdir(cand):
            raise BatchError(f"COMSOL root not found: {cand}")
        return cand
    hits = []
    for pat in COMSOL_ROOT_GLOBS:
        hits.extend(h for h in glob.glob(pat) if os.path.isdir(h))
    if not hits:
        raise BatchError(
            "COMSOL installation not found. Pass --comsol-root or set "
            "COMSOL_ROOT to the .../COMSOL63/Multiphysics directory.")
    return sorted(hits)[-1]                              # newest version


def find_batch_command(root):
    """[executable, ...leading args] of the batch launcher under `root`."""
    env = os.environ.get("COMSOL_BATCH")
    if env:
        return [local_path(env)]
    for exe in glob.glob(os.path.join(root, "bin", "win*", "comsolbatch.exe")):
        return [exe]
    exe = os.path.join(root, "bin", "comsol")            # Linux / macOS
    if os.path.exists(exe):
        return [exe, "batch"]
    raise BatchError(f"no comsolbatch(.exe) or bin/comsol under {root}; "
                     "set COMSOL_BATCH to the batch executable")


def find_javac(root):
    env = os.environ.get("COMSOL_JAVAC")
    if env:
        return local_path(env)
    for pat in ("java/*/jre/bin/javac.exe", "java/*/jre/bin/javac",
                "java/*/bin/javac.exe", "java/*/bin/javac"):
        hits = glob.glob(os.path.join(root, pat))
        if hits:
            return hits[0]
    raise BatchError(f"javac not found under {root}/java; set COMSOL_JAVAC")


# ---- template --------------------------------------------------------------

def template_parameters(mph):
    """{name: (expression, SI value)} from the .mph, read without COMSOL.

    An .mph is a zip; smodel.json lists the global parameters with their
    evaluated values. Best effort: returns {} if the layout is not recognised.
    """
    out = {}
    try:
        with zipfile.ZipFile(mph) as z:
            text = z.read("smodel.json").decode("utf-8", "replace")
        for m in re.finditer(r"\{[^{}]*\"name\":\"([^\"]+)\"[^{}]*\}", text):
            blob = m.group(0)
            val = re.search(r"\"value\":\"([^\"]*)\"", blob)
            num = re.search(r"\"scalarReal\":\"([^\"]*)\"", blob)
            if val and num and '"description"' in blob:
                try:
                    out.setdefault(m.group(1), (val.group(1), float(num.group(1))))
                except ValueError:
                    pass
    except Exception:
        return {}
    return out


def repo_diamond_index():
    """refractive_index_1550nm of configs/materials.yaml (standard library only)."""
    path = os.path.join(ROOT, "configs", "materials.yaml")
    with open(path) as fh:
        m = re.search(r"^\s*refractive_index_1550nm:\s*([0-9.eE+-]+)", fh.read(), re.M)
    if not m:
        raise BatchError(f"refractive_index_1550nm not found in {path}")
    return float(m.group(1))


def template_diamond_eps(mph):
    """Largest relpermittivity stored in the template's materials (the
    diamond one; vacuum is 1), or None if it cannot be read."""
    try:
        with zipfile.ZipFile(mph) as z:
            text = z.read("smodel.json").decode("utf-8", "replace")
        vals = [float(v) for v in re.findall(
            r"\"name\":\"relpermittivity\"[^{}]*?\"scalarReal\":\"([^\"]+)\"", text)]
        vals += [float(v) for v in re.findall(
            r"\"scalarReal\":\"([^\"]+)\"[^{}]*?\"name\":\"relpermittivity\"", text)]
        return max(vals) if vals else None
    except Exception:
        return None


def diamond_index_override(template, tp, studies, overrides, n, keep):
    """Make the optical studies solve with diamond index n (materials.yaml).

    A batch run can only change GLOBAL PARAMETERS (-pname/-plist). Templates
    built by comsol/build_omc_unitcell_rib.m since 2026-11-23, and hole
    templates patched with comsol/patch_hole_templates_index.m, have the
    parameter n_diamond (relpermittivity = n_diamond^2): it is added to the
    overrides. An older template stores a number; if that is not n^2 the run
    is refused (exit 2) unless `keep`, because the index cannot be changed
    from here -- rebuild the template, or use the MPh drivers, which set the
    material property at run time (comsol_client.set_diamond_index).
    Returns a note for the log / summary.
    """
    if not any(s in OPTICAL_STUDIES for s in studies):
        return "no optical study: diamond index not used"
    if "n_diamond" in overrides:
        return f"n_diamond = {overrides['n_diamond']} (given with -p)"
    if "n_diamond" in tp:
        overrides["n_diamond"] = repr(float(n))
        return f"n_diamond = {n} (configs/materials.yaml, passed with -pname)"
    eps = template_diamond_eps(template)
    if eps is not None and abs(eps - n * n) < 1e-6 * n * n:
        return f"template relpermittivity {eps:g} = {n}^2"
    have = (f"relpermittivity {eps:g} (n = {math.sqrt(eps):.4f})"
            if eps else "an unknown permittivity")
    msg = (f"the template has no parameter n_diamond and stores {have}, but "
           f"the diamond index is {n} (configs/materials.yaml). comsol batch "
           "can only set global parameters")
    if keep:
        log("WARNING: " + msg + "; solving with the template's value "
            "(--keep-template-index)")
        return f"TEMPLATE VALUE KEPT: {have}"
    raise BatchError(
        msg + ". Give the template the parameter: rib templates, rebuild "
        "with build_omc_unitcell_rib (MATLAB); hole templates, run "
        "comsol/patch_hole_templates_index.m (MATLAB, patches in place). Or "
        "pass --keep-template-index to solve with the template's value.",
        code=2)


def air_region_note(tp, studies, allow_old=False):
    """Describe the air region the optical studies solve with.

    Rib templates built by comsol/build_omc_unitcell_rib.m since the air
    half-disk have the parameter air_r (= farthest beam corner + air_gap,
    both tunable with -p). An older rib template (parameter th, no air_r)
    has an air box whose outer faces sit inside the evanescent tail of the
    guided modes: it is refused for an optical study (exit 2) unless
    `allow_old` (--allow-old-air, e.g. for a deliberate box-vs-disk Q
    comparison). Mechanics does not use the air. Returns a note for the
    log / summary.
    """
    if not any(s in OPTICAL_STUDIES for s in studies):
        return "no optical study: air region not used"
    if "air_r" in tp:
        note = f"half-disk, air_r = {tp['air_r'][1] * 1e9:.1f} nm"
        if "air_gap" in tp:
            note += f" (air_gap = {tp['air_gap'][1] * 1e9:.1f} nm)"
        return note
    if "th" in tp:
        note = ("OLD AIR BOX (rib template built before the air half-disk, no "
                "parameter air_r; outer faces at x = +-2.5*th, inside the "
                "evanescent tail of the guided modes)")
        if allow_old:
            log("WARNING: " + note + "; solving it anyway (--allow-old-air)")
            return note
        raise BatchError(
            "the template is a rib template with the old air box (no "
            "parameter air_r): its outer faces sit inside the evanescent tail "
            "of the guided modes, so optical Q and frequencies are not to be "
            "trusted. Delete it and rebuild it in MATLAB with LiveLink "
            "(build_omc_unitcell_rib, air half-disk), or pass --allow-old-air "
            "to solve it deliberately (e.g. to compare box and disk).",
            code=2)
    return "as stored in the template (not a rib template)"


def parse_params(items):
    """['a=600[nm]', ...] -> ordered {name: value}; values are COMSOL expressions."""
    out = {}
    for item in items or []:
        if "=" not in item:
            raise BatchError(f"parameter override {item!r} is not NAME=VALUE", 64)
        name, val = (s.strip() for s in item.split("=", 1))
        if not re.fullmatch(r"[A-Za-z_]\w*", name):
            raise BatchError(f"bad parameter name {name!r}", 64)
        if not val or re.search(r"[,\s\"]", val):
            raise BatchError(
                f"value {val!r} of {name}: -plist is comma-separated, so a "
                "value cannot be empty or contain commas, spaces or quotes", 64)
        out[name] = val
    return out


def resolve_studies(items):
    out = []
    for item in items or ["std2"]:
        for s in item.split(","):
            s = s.strip()
            if s:
                tag = STUDY_ALIASES.get(s.lower(), s)
                if tag not in out:
                    out.append(tag)
    return out


# ---- licence watch (query only) --------------------------------------------

def find_lmutil(root):
    env = os.environ.get("LMUTIL")
    if env:
        return local_path(env)
    for pat in ("license/win64/lmutil.exe", "license/*/lmutil.exe",
                "license/*/lmutil"):
        hits = glob.glob(os.path.join(root, pat))
        if hits:
            return hits[0]
    return None


def lm_user():
    """User name as the licence server sees it (the Windows account under WSL)."""
    env = os.environ.get("COMSOL_LM_USER")
    if env:
        return env
    if os.name != "nt" and os.path.exists("/proc/sys/fs/binfmt_misc/WSLInterop"):
        try:
            out = subprocess.run(["cmd.exe", "/c", "echo %USERNAME%"],
                                 capture_output=True, text=True, timeout=20,
                                 cwd="/mnt/c").stdout.strip()
            if out and "%" not in out:
                return out
        except Exception:
            pass
    return os.environ.get("USERNAME") or os.environ.get("USER") or ""


def lmstat_snapshot(lmutil, server, user):
    """{feature: (seats in use, seats held by `user`)} from lmutil lmstat -a.

    Query only: nothing is checked out. Other users' names are not kept.
    """
    out = subprocess.run([lmutil, "lmstat", "-a", "-c", server],
                         capture_output=True, text=True, timeout=60).stdout
    snap, feat = {}, None
    holder = re.compile(r"^\s+(\S+) .* \(v[^)]*\) \([^)]*\), start (.*)$")
    for line in out.replace("\r", "").split("\n"):
        m = re.match(r"^Users of (\S+):\s+\(Total of (\d+) licenses? issued;"
                     r"\s+Total of (\d+) licenses? in use", line)
        if m:
            feat = m.group(1)
            snap.setdefault(feat, [int(m.group(3)), 0])
            continue
        if line.startswith("Users of "):
            feat = None
            continue
        m = holder.match(line) if feat else None
        if m and user and m.group(1).lower() == user.lower():
            n = re.search(r", (\d+) licenses?$", m.group(2))
            snap[feat][1] += int(n.group(1)) if n else 1
    return {k: tuple(v) for k, v in snap.items()}


class LicenceWatch:
    """Samples lmstat in a thread while the batch processes run.

    Reports, per feature, how many more seats this user held than in the
    baseline sample taken before the first process started (a LiveLink or
    Desktop session of the same user is in the baseline, so it is not counted).
    """

    def __init__(self, lmutil, server, user, interval):
        self.lmutil, self.server, self.user = lmutil, server, user
        self.interval = max(2.0, float(interval))
        self.baseline = lmstat_snapshot(lmutil, server, user)
        self.samples = []            # (time, phase, {feature: extra seats of user})
        self.peak_used = {}          # feature -> max total seats in use
        self.phase = "start"
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def _run(self):
        while not self._stop.is_set():
            try:
                snap = lmstat_snapshot(self.lmutil, self.server, self.user)
                extra = {}
                for feat, (used, mine) in snap.items():
                    d = mine - self.baseline.get(feat, (0, 0))[1]
                    if d > 0:
                        extra[feat] = d
                    self.peak_used[feat] = max(self.peak_used.get(feat, 0), used)
                self.samples.append((time.strftime("%H:%M:%S"), self.phase, extra))
            except Exception as e:                       # server hiccup: keep going
                self.samples.append((time.strftime("%H:%M:%S"), self.phase,
                                     {"_error": str(e)}))
            self._stop.wait(self.interval)

    def _extra(self, snap):
        return {f: m - self.baseline.get(f, (0, 0))[1] for f, (_, m) in snap.items()
                if m - self.baseline.get(f, (0, 0))[1] > 0}

    def stop(self, settle_s=30.0):
        """Stop sampling. The server can list a seat for a few seconds after
        the process has gone, so "still held" is only reported if it lasts
        longer than settle_s."""
        self._stop.set()
        self._thread.join(timeout=70)
        after, t0 = {}, time.time()
        while True:
            try:
                after = lmstat_snapshot(self.lmutil, self.server, self.user)
            except Exception:
                after = {}
            if not self._extra(after) or time.time() - t0 > settle_s:
                break
            time.sleep(3.0)
        released_after = round(time.time() - t0, 1)
        by_phase = {}
        for _, phase, extra in self.samples:
            feats = by_phase.setdefault(phase, {})
            for feat, n in extra.items():
                if not feat.startswith("_"):
                    feats[feat] = max(feats.get(feat, 0), n)
        still = self._extra(after)
        return dict(server=self.server, user_matched=bool(self.user),
                    interval_s=self.interval, n_samples=len(self.samples),
                    features_checked_out_by_phase=by_phase,
                    still_held_after_exit=still,
                    released_within_s=None if still else released_after,
                    samples=[dict(time=t, phase=p, extra=e)
                             for t, p, e in self.samples])


# ---- running COMSOL --------------------------------------------------------

def tail(path, n=25):
    try:
        with open(path, errors="replace") as fh:
            return "".join(fh.readlines()[-n:])
    except OSError:
        return ""


def run_logged(cmd, logfile, what, verbose, capture=False):
    """Run one COMSOL process; on failure show the end of its log and raise.

    Returns (seconds, stdout text). stdout is captured unless `verbose` shows
    it on the terminal; capture=True captures it in both cases (and echoes it
    to stderr when verbose).
    """
    log(f"{what}: " + " ".join(f'"{c}"' if " " in c else c for c in cmd))
    t0 = time.time()
    piped = capture or not verbose
    proc = subprocess.run(cmd, stdout=subprocess.PIPE if piped else None,
                          stderr=subprocess.STDOUT if piped else None,
                          text=True, errors="replace")
    dt = time.time() - t0
    if verbose and piped and proc.stdout:
        sys.stderr.write(proc.stdout)
    text = (proc.stdout or "") + tail(logfile, 400)
    # COMSOL batch can exit 0 after printing an error banner to the log
    # (seen with "Input file does not support batch licenses").
    banner = "*****Error" in text
    if proc.returncode != 0 or banner or LICENCE_ERROR_RE.search(text):
        # The log text goes to stderr so wait_for_comsol_seat.sh -r can see a
        # licence error and retry.
        sys.stderr.write((proc.stdout or "")[-4000:])
        sys.stderr.write(tail(logfile))
        licence = bool(LICENCE_ERROR_RE.search(text))
        raise BatchError(
            f"{what} failed after {dt:.0f} s (exit {proc.returncode})"
            + ("; licence error: no seat was available" if licence else "")
            + f". Log: {logfile}", code=3 if licence else 1)
    log(f"{what}: done in {dt:.0f} s")
    return dt, (proc.stdout or "").replace("\r", "")


JAVA_TEMPLATE = r"""// Generated by scripts/comsol_batch.py -- do not edit.
import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.*;
import java.util.*;

public class OmcEigExtract {
  // {study tag, solved .mph}
  static final String[][] JOBS = {
@JOBS@
  };
  // {column, global expression}: evaluated per mode with EvalGlobal
  static final String[][] GLOBALS = {
@GLOBALS@
  };
  // {study tag, column, expression, named selection or "" for all domains}
  static final String[][] INTEGRALS = {
@INTEGRALS@
  };
  static int counter = 0;

  public static void main(String[] args) { run(); }

  public static Model run() {
    Model model = null;
    for (int j = 0; j < JOBS.length; j++) {
      // No file is written from here: COMSOL's default security settings
      // ("File system access") forbid it for class files. The CSV text goes
      // to standard output between markers and comsol_batch.py saves it.
      try {
        model = ModelUtil.load("omcx" + (counter++), JOBS[j][1]);
        String csv = extract(model, JOBS[j][0]);
        System.out.println("@@OMC-CSV-BEGIN " + JOBS[j][0] + "@@");
        System.out.print(csv);
        System.out.println("@@OMC-CSV-END " + JOBS[j][0] + "@@");
      } catch (Throwable e) {
        StringWriter sw = new StringWriter();
        e.printStackTrace(new PrintWriter(sw));
        System.out.println("@@OMC-ERROR " + JOBS[j][0] + "@@ " + e);
        System.out.println(sw.toString());
      }
      System.out.flush();
    }
    return model;
  }

  static double[] flat(double[][] m) {
    int n = 0;
    for (double[] r : m) n += r.length;
    double[] out = new double[n];
    int k = 0;
    for (double[] r : m) for (double v : r) out[k++] = v;
    return out;
  }

  // Global expression on a solution dataset: {real[], imag[]}, one per mode.
  static double[][] evalGlobal(Model model, String dset, String expr) {
    String tag = "omcg" + (counter++);
    NumericalFeature f = model.result().numerical().create(tag, "EvalGlobal");
    try {
      f.set("data", dset);
      f.set("expr", new String[]{expr});
      double[][][] r = f.computeResult();
      double[] re = flat(r[0]);
      double[] im = (r.length > 1 && f.isComplex()) ? flat(r[1]) : new double[re.length];
      return new double[][]{re, im};
    } finally {
      model.result().numerical().remove(tag);
    }
  }

  // Volume integral of expr per mode (real part).
  static double[] integrate(Model model, String dset, String expr, String sel) {
    String tag = "omci" + (counter++);
    NumericalFeature f = model.result().numerical().create(tag, "IntVolume");
    try {
      f.set("data", dset);
      if (sel.length() == 0) f.selection().all();
      else f.selection().named(sel);
      f.set("expr", new String[]{expr});
      return flat(f.getReal());
    } finally {
      model.result().numerical().remove(tag);
    }
  }

  static String num(double v) { return String.format(Locale.ROOT, "%.15g", v); }

  static String extract(Model model, String study) throws Exception {
    // Dataset of the study: the one whose solution was computed by it.
    String dset = null;
    for (String t : model.result().dataset().tags()) {
      try {
        String sol = model.result().dataset(t).getString("solution");
        if (sol != null && study.equals(model.sol(sol).study())
            && !model.sol(sol).isEmpty()) { dset = t; break; }
      } catch (Throwable e) { }
    }
    if (dset == null)
      throw new RuntimeException("no computed solution dataset for study " + study);

    double[][] f = evalGlobal(model, dset, "freq");
    int n = f[0].length;

    List<String> cols = new ArrayList<String>();
    List<double[]> data = new ArrayList<double[]>();
    List<String> notes = new ArrayList<String>();
    for (String[] g : GLOBALS) {
      try {
        double[] v = evalGlobal(model, dset, g[1])[0];
        if (v.length != n) throw new RuntimeException("got " + v.length + " values for " + n + " modes");
        cols.add(g[0]); data.add(v);
      } catch (Throwable e) { notes.add("column " + g[0] + " (" + g[1] + ") failed: " + e); }
    }
    for (String[] g : INTEGRALS) {
      if (!g[0].equals(study)) continue;
      try {
        double[] v = integrate(model, dset, g[2], g[3]);
        if (v.length != n) throw new RuntimeException("got " + v.length + " values for " + n + " modes");
        cols.add(g[1]); data.add(v);
      } catch (Throwable e) { notes.add("column " + g[1] + " (" + g[2] + ") failed: " + e); }
    }

    StringWriter sw = new StringWriter();
    PrintWriter w = new PrintWriter(sw);
    w.println("# study=" + study);
    w.println("# study_label=" + model.study(study).label());
    w.println("# dataset=" + dset);
    for (String p : model.param().varnames()) {
      try { w.println("# param " + p + "=" + num(model.param().evaluate(p))
                      + " expr=" + model.param().get(p)); }
      catch (Throwable e) { notes.add("parameter " + p + " failed: " + e); }
    }
    for (String s : notes) w.println("# note " + s.replace('\n', ' ').replace('\r', ' '));
    StringBuilder h = new StringBuilder("mode,f_real_Hz,f_imag_Hz");
    for (String c : cols) h.append(",").append(c);
    w.println(h);
    for (int i = 0; i < n; i++) {
      StringBuilder r = new StringBuilder();
      r.append(i + 1).append(",").append(num(f[0][i])).append(",").append(num(f[1][i]));
      for (double[] d : data) r.append(",").append(num(d[i]));
      w.println(r);
    }
    w.flush();
    return sw.toString();
  }
}
"""


def java_string(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def java_rows(rows):
    return ",\n".join("    {" + ", ".join(java_string(x) for x in row) + "}"
                      for row in rows)


def write_extractor(outdir, jobs, global_exprs, integrals):
    src = (JAVA_TEMPLATE
           .replace("@JOBS@", java_rows(jobs))
           .replace("@GLOBALS@", java_rows(global_exprs))
           .replace("@INTEGRALS@", java_rows(integrals)))
    path = os.path.join(outdir, "OmcEigExtract.java")
    with open(path, "w", newline="\n") as fh:
        fh.write(src)
    return path


def compile_extractor(javac, root, java_file, outdir):
    win = javac.lower().endswith(".exe")
    sep = "\\" if (win and os.name != "nt") or os.name == "nt" else "/"
    cp = exe_path(os.path.join(root, "plugins"), win) + sep + "*"
    cmd = [javac, "-nowarn", "-cp", cp, "-d", exe_path(outdir, win),
           exe_path(java_file, win)]
    proc = subprocess.run(cmd, capture_output=True, text=True, errors="replace")
    if proc.returncode != 0:
        raise BatchError("javac failed for the extraction class:\n"
                         + proc.stdout + proc.stderr)
    return os.path.join(outdir, "OmcEigExtract.class")


# ---- results ---------------------------------------------------------------

def read_eig_csv(path):
    meta, params, notes, rows = {}, {}, [], []
    with open(path, newline="") as fh:
        lines = fh.read().splitlines()
    body = []
    for line in lines:
        if line.startswith("# param "):
            m = re.match(r"# param (\w+)=(\S+) expr=(.*)$", line)
            if m:
                params[m.group(1)] = dict(value=float(m.group(2)), expr=m.group(3))
        elif line.startswith("# note "):
            notes.append(line[7:])
        elif line.startswith("# ") and "=" in line:
            k, v = line[2:].split("=", 1)
            meta[k] = v
        elif line.strip():
            body.append(line)
    for row in csv.DictReader(body):
        rows.append({k: (int(v) if k == "mode" else float(v)) for k, v in row.items()})
    return meta, params, notes, rows


def summarise(study, meta, params, notes, rows):
    a = params.get("a", {}).get("value", float("nan"))
    kF = params.get("kF", {}).get("value", float("nan"))
    f_ll = C0 * kF / (2.0 * math.pi)                 # light line at this k
    modes = []
    for r in sorted(rows, key=lambda r: r["f_real_Hz"]):
        fr, fi = r["f_real_Hz"], r["f_imag_Hz"]
        m = dict(mode=r["mode"], f_real_Hz=fr, f_imag_Hz=fi,
                 f_real_THz=fr / 1e12, f_imag_THz=fi / 1e12,
                 wavelength_nm=(C0 / fr * 1e9) if fr > 0 else None,
                 Q=(fr / (2.0 * abs(fi))) if fi != 0 else None,
                 below_light_line=bool(fr < f_ll) if f_ll == f_ll else None,
                 pseudo_mode=bool(fr < PSEUDO_MODE_HZ),
                 rel_to_target=(fr - F_TARGET) / F_TARGET)
        if r.get("We_all"):
            m["energy_fraction_in_beam"] = r.get("We_beam", float("nan")) / r["We_all"]
        tot = sum(r.get(k, 0.0) for k in ("Ex2_beam", "Ey2_beam", "Ez2_beam"))
        if tot > 0:
            for c in "xyz":
                m[f"E{c}_fraction_in_beam"] = r[f"E{c}2_beam"] / tot
        modes.append(m)
    return dict(study=study, study_label=meta.get("study_label"),
                dataset=meta.get("dataset"), a_m=a, kF_per_m=kF,
                k_frac=(kF * a / math.pi) if a == a and kF == kF else None,
                light_line_Hz=f_ll, light_line_zone_edge_Hz=C0 / (2.0 * a),
                f_target_Hz=F_TARGET, params=params, notes=notes, modes=modes)


def print_summary(s, optical):
    k = s["k_frac"]
    print(f"\n{s['study']} '{s['study_label']}'  dataset {s['dataset']}")
    cell = ", ".join(f"{n}={s['params'][n]['value'] * 1e9:.4g} nm"
                     for n in ("a", "w", "t", "s", "th", "r1", "r2", "hx", "hy")
                     if n in s["params"])
    print(f"  cell: {cell}")
    print(f"  kF = {s['kF_per_m']:.6g} 1/m  (k*a/pi = {k:.4f})" if k is not None
          else "  kF: not available")
    if not optical:
        print("   mode   Re f [GHz]    Im f [GHz]")
        for m in s["modes"]:
            print(f"  {m['mode']:5d} {m['f_real_Hz'] / 1e9:12.5f} "
                  f"{m['f_imag_Hz'] / 1e9:13.3e}")
        return
    print(f"  light line at this k: {s['light_line_Hz'] / 1e12:.2f} THz"
          f"   (c/2a = {s['light_line_zone_edge_Hz'] / 1e12:.2f} THz)"
          f"   target {F_TARGET / 1e12:.2f} THz")
    print("   mode  Re f [THz]  Im f [THz]   lambda [nm]          Q  light line"
          "  W_beam   Ex/Ey/Ez share in beam")
    for m in s["modes"]:
        ll = ("pseudo" if m["pseudo_mode"]
              else "below" if m["below_light_line"] else "ABOVE")
        q = f"{m['Q']:10.3g}" if m["Q"] is not None else "       inf"
        wb = (f"{m['energy_fraction_in_beam']:6.3f}"
              if "energy_fraction_in_beam" in m else "     -")
        pol = ("/".join(f"{m[f'E{c}_fraction_in_beam']:.2f}" for c in "xyz")
               if "Ex_fraction_in_beam" in m else "-")
        lam = f"{m['wavelength_nm']:11.1f}" if m["wavelength_nm"] else "          -"
        near = " *" if abs(m["rel_to_target"]) < 0.30 and not m["pseudo_mode"] else ""
        print(f"  {m['mode']:5d} {m['f_real_THz']:11.3f} {m['f_imag_THz']:11.3e} "
              f"{lam} {q}  {ll:>10}  {wb}   {pol}{near}")
    guided = [m for m in s["modes"] if m["below_light_line"] and not m["pseudo_mode"]]
    print("  * within +/-30 % of the target. Guided (below the light line): "
          + (", ".join(f"{m['f_real_THz']:.2f}" for m in guided) + " THz"
             if guided else "none"))
    for n in s["notes"]:
        print(f"  note: {n}")


# ---- main ------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Solve studies of a COMSOL template with `comsol batch` "
                    "and extract the eigenfrequencies.",
        epilog="See scripts/run_comsol_batch.sh for examples and the licence notes.")
    ap.add_argument("--template", default=DEFAULT_TEMPLATE,
                    help="input .mph (default: comsol/omc_unitcell_rib.mph)")
    ap.add_argument("--study", action="append", metavar="TAG",
                    help="study tag, repeatable or comma-separated (default: "
                         "std2). The rib labels 'opt TE', 'opt TM', 'mech sym', "
                         "'mech antisym' are accepted too.")
    ap.add_argument("-p", "--param", action="append", metavar="NAME=VALUE",
                    help="global parameter override, e.g. -p a=600[nm] "
                         "(repeatable; passed as -pname/-plist). kF is left at "
                         "the template's pi/a (X point) unless given.")
    ap.add_argument("--outdir", help="output folder (default: a new "
                                     "results/batch_<time stamp>/)")
    ap.add_argument("--license", choices=("batch", "interactive"),
                    default=os.environ.get("COMSOL_BATCH_LICENSE", "batch"),
                    help="batch (default): pass -usebatchlic, so the job uses "
                         "the *BATCH features. interactive: the ordinary "
                         "features (COMSOL, WAVEOPTICS/RF, ...).")
    ap.add_argument("--comsol-root", help="the .../COMSOLxx/Multiphysics folder "
                                          "(default: env COMSOL_ROOT, then the "
                                          "usual install locations)")
    ap.add_argument("--expr", action="append", metavar="NAME=EXPR", default=[],
                    help="extra global expression per mode, e.g. "
                         "--expr 'fy=intop1(solid.uy^2)' (repeatable)")
    ap.add_argument("--no-integrals", action="store_true",
                    help="skip the optical energy integrals (frequencies only)")
    ap.add_argument("--watch-license", action="store_true",
                    help="poll lmstat while COMSOL runs and record which "
                         "features this user's job checked out "
                         "(license_watch.json)")
    ap.add_argument("--watch-interval", type=float, default=5.0,
                    help="seconds between lmstat polls with --watch-license "
                         "(default 5; only while the job runs)")
    ap.add_argument("--check-license", action="store_true",
                    help="print the licence features the template needs "
                         "(comsolbatch -checklicense) and exit; solves nothing")
    ap.add_argument("--extract-only", action="store_true",
                    help="skip the solve; extract from solved_<study>.mph "
                         "already in --outdir")
    ap.add_argument("--comsol-arg", action="append", default=[], metavar="ARG",
                    help="extra argument for the solve command (repeatable), "
                         "e.g. --comsol-arg=-np --comsol-arg=4")
    ap.add_argument("--n-diamond", type=float, default=None,
                    help="diamond refractive index of the optical studies "
                         "(default: configs/materials.yaml, 2.386); set through "
                         "the template parameter n_diamond")
    ap.add_argument("--keep-template-index", action="store_true",
                    help="solve an old template (no n_diamond parameter) with "
                         "the permittivity stored in it instead of refusing")
    ap.add_argument("--allow-old-air", action="store_true",
                    help="solve the optical studies of an old rib template "
                         "(air box, no parameter air_r) instead of refusing, "
                         "e.g. to compare its Q with the rebuilt half-disk "
                         "template")
    ap.add_argument("--verbose", action="store_true",
                    help="show COMSOL's own log on the terminal as well")
    args = ap.parse_args(argv)

    template = local_path(args.template)
    if not os.path.exists(template):
        raise BatchError(
            f"template not found: {template}\n"
            "The rib template is built once in MATLAB with LiveLink "
            "(comsol/build_omc_unitcell_rib.m); it cannot be built from here.",
            code=2)
    studies = resolve_studies(args.study)
    overrides = parse_params(args.param)
    global_exprs = []
    for item in args.expr:
        if "=" not in item:
            raise BatchError(f"--expr {item!r} is not NAME=EXPR", 64)
        global_exprs.append([s.strip() for s in item.split("=", 1)])

    root = find_comsol_root(args.comsol_root)
    batch = find_batch_command(root)
    win = batch[0].lower().endswith(".exe")
    lic = ["-usebatchlic"] if args.license == "batch" else []

    if args.check_license:
        proc = subprocess.run(batch + ["-checklicense", exe_path(template, win)],
                              capture_output=True, text=True, errors="replace")
        print(proc.stdout.replace("\r", "").strip())
        return proc.returncode

    tp = template_parameters(template)
    if tp:
        log("template parameters: " + ", ".join(
            f"{n}={tp[n][0]}" for n in ("a", "w", "th", "t", "s", "r1", "r2", "kF",
                                        "air_gap", "air_r")
            if n in tp))
        if "kF" in tp and "kF" not in overrides and tp["kF"][0].replace(" ", "") != "pi/a":
            log(f"WARNING: the template's kF is {tp['kF'][0]!r}, not pi/a; "
                "pass -p kF=... to choose the k-point")
        unknown = [n for n in overrides if n not in tp]
        if unknown:
            raise BatchError(f"template has no parameter(s) {unknown}; it has "
                             f"{sorted(tp)}", 64)
    index_note = "not solved here (--extract-only)"
    if not args.extract_only:
        n_diamond = (args.n_diamond if args.n_diamond is not None
                     else repo_diamond_index())
        index_note = diamond_index_override(template, tp, studies, overrides,
                                            n_diamond, args.keep_template_index)
        log(f"diamond index: {index_note}")
    # refused for an old rib air box only when something is solved here
    air_note = air_region_note(tp, studies,
                               allow_old=args.allow_old_air or args.extract_only)
    log(f"air region: {air_note}")

    outdir = local_path(args.outdir) if args.outdir else os.path.join(
        ROOT, "results", f"batch_{_dt.datetime.now():%Y%m%d_%H%M%S}")
    os.makedirs(outdir, exist_ok=True)
    outdir = os.path.abspath(outdir)
    javac = find_javac(root)
    log(f"COMSOL: {root}")
    log(f"output: {outdir}")
    log(f"licence pool: {args.license}" + (" (-usebatchlic)" if lic else ""))

    watch = None
    if args.watch_license:
        lmutil = find_lmutil(root)
        server = os.environ.get("COMSOL_LICENSE_SERVER",
                                "27005@license.seas.harvard.edu")
        user = lm_user()
        if lmutil is None:
            log("WARNING: lmutil not found; --watch-license is off")
        else:
            if not user:
                log("WARNING: user name unknown (set COMSOL_LM_USER); the watch "
                    "cannot tell this job's seats from other users'")
            watch = LicenceWatch(lmutil, server, user, args.watch_interval)

    timings, rc, extract_errors = {}, 0, {}
    try:
        # ---- 1. solve --------------------------------------------------------
        jobs = []
        for study in studies:
            solved = os.path.join(outdir, f"solved_{study}.mph")
            if not args.extract_only:
                before = set(glob.glob(os.path.join(outdir, "*.mph")))
                logfile = os.path.join(outdir, f"solve_{study}.log")
                cmd = batch + ["-inputfile", exe_path(template, win),
                               "-outputfile", exe_path(solved, win),
                               "-study", study,
                               "-batchlog", exe_path(logfile, win)] + lic
                if overrides:
                    cmd += ["-pname", ",".join(overrides),
                            "-plist", ",".join(overrides.values())]
                cmd += args.comsol_arg
                if args.verbose:
                    cmd.append("-batchlogout")
                if watch:
                    watch.phase = f"solve {study}"
                timings[f"solve_{study}_s"] = run_logged(
                    cmd, logfile, f"solve {study}", args.verbose)[0]
                if not os.path.exists(solved):
                    # With -pname/-plist COMSOL may append the parameter tuple
                    # to the output name; take the one new .mph.
                    new = sorted(set(glob.glob(os.path.join(outdir, "*.mph"))) - before)
                    if len(new) != 1:
                        raise BatchError(
                            f"solve {study} finished but {solved} was not "
                            f"written (new .mph files: {new}). Log: {logfile}")
                    log(f"solve {study}: COMSOL named the output "
                        f"{os.path.basename(new[0])}")
                    solved = new[0]
            elif not os.path.exists(solved):
                raise BatchError(f"--extract-only: {solved} not found")
            jobs.append([study, solved, os.path.join(outdir, f"eig_{study}.csv")])

        # ---- 2. extract ------------------------------------------------------
        integrals = [] if args.no_integrals else [
            [study, col, expr, sel] for study in studies if study in OPTICAL_STUDIES
            for col, expr, sel in OPTICAL_INTEGRALS]
        java_jobs = [[s, exe_path(m, win).replace("\\", "/")] for s, m, c in jobs]
        java_file = write_extractor(outdir, java_jobs, global_exprs, integrals)
        cls = compile_extractor(javac, root, java_file, outdir)
        logfile = os.path.join(outdir, "extract.log")
        # No -usebatchlic here: COMSOL refuses it for a class file ("Input
        # file does not support batch licenses"), so this short step uses the
        # ordinary features whatever --license says.
        cmd = batch + ["-inputfile", exe_path(cls, win), "-nosave",
                       "-batchlog", exe_path(logfile, win)]
        if args.verbose:
            cmd.append("-batchlogout")
        if watch:
            watch.phase = "extract"
        timings["extract_s"], text = run_logged(cmd, logfile, "extract",
                                                args.verbose, capture=True)
        with open(os.path.join(outdir, "extract_stdout.txt"), "w") as fh:
            fh.write(text)
        for study, solved, csv_path in jobs:
            m = re.search(rf"@@OMC-CSV-BEGIN {re.escape(study)}@@\n(.*?)"
                          rf"@@OMC-CSV-END {re.escape(study)}@@", text, re.S)
            if m:
                with open(csv_path, "w", newline="\n") as fh:
                    fh.write(m.group(1))
            else:
                e = re.search(rf"@@OMC-ERROR {re.escape(study)}@@.*", text, re.S)
                extract_errors[study] = e.group(0)[:3000] if e else "no output"
    finally:
        licence = None
        if watch:
            watch.phase = "after exit"
            time.sleep(min(watch.interval, 5.0))
            licence = watch.stop()
            with open(os.path.join(outdir, "license_watch.json"), "w") as fh:
                json.dump(licence, fh, indent=2)

    # ---- 3. summary ----------------------------------------------------------
    summaries = []
    for study, solved, csv_path in jobs:
        if study in extract_errors:
            err = extract_errors[study]
            sys.stderr.write(err + "\n")
            log(f"extract {study}: failed; the solved model is kept, rerun "
                f"with --extract-only --outdir {outdir}. Log: {logfile}")
            rc = 3 if LICENCE_ERROR_RE.search(err) else 1
            continue
        meta, params, notes, rows = read_eig_csv(csv_path)
        s = summarise(study, meta, params, notes, rows)
        s.update(solved_mph=solved, csv=csv_path)
        summaries.append(s)
        print_summary(s, optical=study in OPTICAL_STUDIES)

    summary = dict(template=os.path.abspath(template), comsol_root=root,
                   licence_pool=args.license, overrides=overrides,
                   diamond_index=index_note, air_region=air_note,
                   timings=timings, studies=summaries, licence_watch=licence)
    with open(os.path.join(outdir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    if licence:
        print("\nLicence features checked out by this job (lmstat, this user, "
              "above the baseline before the job):")
        for phase, feats in licence["features_checked_out_by_phase"].items():
            if phase in ("start", "after exit"):     # see the last line instead
                continue
            print(f"  {phase:14s} " + (", ".join(
                f"{f} x{n}" if n > 1 else f for f, n in sorted(feats.items()))
                or "(none seen)"))
        still = licence["still_held_after_exit"]
        print("  after exit     " + (", ".join(sorted(still)) + "  <-- STILL HELD"
                                     if still else "all released (within "
                                     f"{licence['released_within_s']:.0f} s)"))
    print(f"\nOutput folder: {outdir}")
    return rc


if __name__ == "__main__":
    try:
        sys.exit(main())
    except BatchError as e:
        print(f"comsol_batch: {e}", file=sys.stderr)
        sys.exit(e.code)
    except KeyboardInterrupt:
        print("comsol_batch: interrupted", file=sys.stderr)
        sys.exit(130)
