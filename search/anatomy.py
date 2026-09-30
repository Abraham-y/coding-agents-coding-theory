"""Produces the cost and cause-of-death statistics in the paper's anatomy section.

Exists because an earlier draft reported those numbers from an ad-hoc filter with
no script behind it -- a table of statistics resting on an unaudited summary,
which is the exact failure mode the paper is about. The filter is now explicit:
every results/*.json carrying an `elapsed_seconds` or `solve_seconds` field.

Run:  python search/anatomy.py
"""
from __future__ import annotations
import json, glob, os, re, statistics, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def timed_runs():
    """Every timed run we can account for, across ALL record formats.

    An earlier version globbed results/*.json for two key names only. Adversarial
    review showed that was a filter artifact, not a census: it missed 22 runs
    logging the same quantity under `wall_seconds`/`seconds`, and all 35 KaMIS
    runs, which record `Time:` in .log files. It also globbed the working tree, so
    the numbers moved when a verification run wrote a file. This version reads
    git-tracked files only and every format we actually emit.
    """
    import subprocess
    tracked = subprocess.run(["git", "ls-files", "results"], cwd=ROOT,
                             capture_output=True, text=True).stdout.split()
    out = []
    for rel in tracked:
        f = os.path.join(ROOT, rel)
        if rel.endswith(".json"):
            try:
                d = json.load(open(f))
            except Exception:
                continue
            if not isinstance(d, dict):
                continue
            for k in ("elapsed_seconds", "solve_seconds", "wall_seconds", "seconds"):
                t = d.get(k)
                if isinstance(t, (int, float)) and t > 0:
                    out.append((float(t), os.path.basename(rel)))
                    break
        elif rel.endswith(".log"):
            try:
                txt = open(f, errors="replace").read()
            except Exception:
                continue
            m = re.findall(r"^Time:\s*([0-9.]+)", txt, re.M)
            if m:
                out.append((max(float(x) for x in m), os.path.basename(rel)))
    return sorted(out, reverse=True)


# Runs that happened but leave no machine-readable record in results/. Counted
# explicitly rather than silently dropped, because the largest single run of the
# campaign is one of them and an earlier draft quoted a "share of all compute"
# that excluded it.
UNLOGGED_HOURS = {"tensor-rank F_3 search (tensor-rank-f3/README.md)": 25.5 * 10}


def entry_census():
    txt = open(os.path.join(ROOT, "notes", "failed-ideas.md")).read()
    return len(re.findall(r"^\*\*(.+?)\*\*", txt, re.M))

if __name__ == "__main__":
    runs = timed_runs()
    ts = [t for t, _ in runs]
    buckets = {"<10s": 0, "10s-1m": 0, "1m-10m": 0, "10m-1h": 0, "1h-6h": 0, ">6h": 0}
    for t in ts:
        k = ("<10s" if t < 10 else "10s-1m" if t < 60 else "1m-10m" if t < 600
             else "10m-1h" if t < 3600 else "1h-6h" if t < 21600 else ">6h")
        buckets[k] += 1
    ext = sum(UNLOGGED_HOURS.values())
    print(f"filter          : git-tracked results/*.json (any timing key) + *.log Time:")
    print(f"timed runs      : {len(ts)}")
    print(f"total compute   : {sum(ts)/3600:.1f} h")
    print(f"median run      : {statistics.median(ts):.0f} s")
    print(f"over one hour   : {sum(1 for t in ts if t >= 3600)}")
    print(f"top three       : {sum(ts[:3])/3600:.1f} h "
          f"({100*sum(ts[:3])/sum(ts):.0f}% of logged)")
    print(f"unlogged (core-h): {ext:.0f}  ({', '.join(UNLOGGED_HOURS)})")
    print(f"top three share of logged+unlogged: "
          f"{100*sum(ts[:3])/3600/(sum(ts)/3600+ext):.0f}%")
    for k, v in buckets.items():
        print(f"  {k:>7}: {v}")
    print(f"bolded log entries: {entry_census()}")
