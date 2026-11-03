#!/usr/bin/env python3
"""Two-directional test suite for verify_big2 (C, bit-parallel Levenshtein).

verify_big2.c is one of two independent verifiers for quaternary edit-metric
codes at n = 12 (CLAUDE.md rule 3).  This suite runs the compiled binary and
checks it against a reference Levenshtein DP written here (textbook
Wagner-Fischer on Python strings; nothing shared with the C code).

ACCEPT (must exit 0, 'valid: True'):
  * every best/edit4_*.txt and best/editq_*.txt.  d comes from '# d', or is 3
    (via --d 3) for the edit4 files that lack it; the two legacy edit4 files
    with no '# n' header get --n from their file name.  The printed sha256
    must equal the header's 'sha256(codewords)' where present, and always
    hashlib's sha256 of the sorted codewords each followed by '\\n'.
    Files of <= 400 words are also re-checked pair by pair with the reference
    DP (true minimum distance computed, must be >= d);
  * small codes built greedily here, with the true minimum distance D found
    by the reference DP over all pairs (accepted at D and at D-1);
  * random small word sets that the reference DP says are valid;
  * formatting variants (CRLF, blank lines, trailing spaces, no final newline).
REJECT (must exit 1, 'valid: False'), built from accepted codes:
  * one codeword replaced by a word within distance d-1 of another codeword:
    (a) d-1 substitutions; (b) a deletion + insertion (+ d-3 substitutions)
    chosen so the Hamming distance is >= d, i.e. a Hamming-only verifier
    would wrongly accept it;
  * a duplicated word (appended, and replacing another word);
  * a word one symbol short and one symbol long;
  * a symbol outside 0-3 ('4', 'x', an inner space);
  * a code whose true minimum distance is exactly D checked with D+1 (via the
    header and via --d); the reported first violating pair must be the
    reference's first pair in file order;
  * random small word sets that the reference DP says are invalid, with the
    same exact-first-pair check.
REFUSE (must exit 2): missing '# d' and no --d (including the real edit4
  files that lack '# d'), missing '# n', header/option conflicts, q != 4,
  d = 0, n out of range, malformed headers, bad options, missing file.
OTHER: verify_big2 --selftest; identical stdout with 1 and 2 threads.

SCALING (skip with --no-scaling): 100,000 seeded random distinct words of
length 12 at d = 3 (time to the first violation), and full all-pairs runs at
d = 1 (always valid for distinct words) for 20k, 50k and 100k words, giving
pairs/second and extrapolated times for 100k and 200k words.
--scaling-200k additionally runs the full 200k-word scan.

Run:   python tests/test_verify_big2.py [--no-scaling] [--scaling-200k]
Log:   results/verify_big2_tests.json
Exit:  0 iff every expectation holds.
"""
from __future__ import annotations

import datetime
import hashlib
import json
import os
import platform
import random
import re
import resource
import subprocess
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "verify_big2.c"
BIN = ROOT / "verify_big2"
BUILD_CMD = ["cc", "-O3", "-o", "verify_big2", "verify_big2.c", "-lpthread"]
LOG = ROOT / "results" / "verify_big2_tests.json"
SEED = 20261009
SCALE_SEED = 20261010
ALPHABET = "0123"


# --------------------------------------------------------------------------
# reference Levenshtein distance (boring on purpose)
# --------------------------------------------------------------------------

def lev(a: str, b: str) -> int:
    """Textbook Wagner-Fischer: unit-cost insert / delete / substitute."""
    prev = list(range(len(b) + 1))
    for i in range(1, len(a) + 1):
        cur = [i] + [0] * len(b)
        for j in range(1, len(b) + 1):
            cur[j] = min(prev[j] + 1,
                         cur[j - 1] + 1,
                         prev[j - 1] + (a[i - 1] != b[j - 1]))
        prev = cur
    return prev[len(b)]


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def ref_first_violation(words, d):
    """First pair (i, j), i < j, in file order with lev < d; None if valid."""
    for i in range(len(words)):
        for j in range(i + 1, len(words)):
            if lev(words[i], words[j]) < d:
                return i, j
    return None


def ref_min_distance(words):
    return min(lev(words[i], words[j])
               for i in range(len(words)) for j in range(i + 1, len(words)))


def sha_of(words) -> str:
    return hashlib.sha256("".join(w + "\n" for w in sorted(words)).encode()).hexdigest()


# --------------------------------------------------------------------------
# running the binary
# --------------------------------------------------------------------------

def build():
    if not BIN.exists() or BIN.stat().st_mtime < SRC.stat().st_mtime:
        print("building:", " ".join(BUILD_CMD))
        subprocess.run(BUILD_CMD, cwd=ROOT, check=True)


def child_cpu_seconds():
    r = resource.getrusage(resource.RUSAGE_CHILDREN)
    return r.ru_utime + r.ru_stime


def run_bin(path, opts=(), timeout=7200):
    """Run the binary once; returns (exit, parsed stdout, stdout, stderr, stats).

    stats has wall time and the child's CPU time (user + sys).  CPU time is
    the robust measure on a laptop that may sleep mid-run: wall time then
    includes the time asleep, CPU time does not."""
    c0 = child_cpu_seconds()
    t0 = time.perf_counter()
    p = subprocess.run([str(BIN), *map(str, opts), str(path)],
                       capture_output=True, text=True, timeout=timeout)
    wall = time.perf_counter() - t0
    cpu = child_cpu_seconds() - c0
    out = {}
    for line in p.stdout.splitlines():
        k, sep, v = line.partition(": ")
        if sep:
            out[k] = v
    stats = {}
    m = re.search(r"pairs_checked=(\d+) seconds=([0-9.]+)", p.stderr)
    if m:
        stats = {"pairs_checked": int(m.group(1)), "binary_seconds": float(m.group(2))}
    stats["wall_seconds"] = wall
    stats["cpu_seconds"] = cpu
    return p.returncode, out, p.stdout, p.stderr, stats


def write_code(tmp: Path, name, words, n=None, d=None, extra=(), newline="\n",
               final_newline=True):
    """Write an artifact; return (path, line number of words[0])."""
    header = ["# test artifact generated by tests/test_verify_big2.py"]
    if n is not None:
        header.append(f"# n = {n}")
    if d is not None:
        header.append(f"# d = {d}")
    header += list(extra)
    text = newline.join(header + list(words)) + (newline if final_newline else "")
    path = tmp / name
    with open(path, "w", newline="") as f:
        f.write(text)
    return path, len(header) + 1


class Suite:
    def __init__(self):
        self.counts = Counter()   # accept / reject / refuse / other
        self.kinds = Counter()    # finer labels
        self.failures = []

    def _record(self, cat, kind, label, problems, raw=""):
        self.counts[cat] += 1
        self.kinds[f"{cat}:{kind}"] += 1
        if problems:
            msg = f"[{cat}:{kind}] {label}: {'; '.join(problems)}"
            self.failures.append(msg)
            print("FAIL", msg, "|", raw.strip().replace("\n", " / ")[:400])
        return not problems

    def accept(self, kind, label, path, opts=(), size=None, sha=None, res=None):
        res = res or run_bin(path, opts)
        rc, out, raw, err, _ = res
        p = []
        if rc != 0:
            p.append(f"exit {rc}")
        if out.get("valid") != "True":
            p.append(f"valid={out.get('valid')}")
        if size is not None and out.get("size") != str(size):
            p.append(f"size {out.get('size')} != {size}")
        if sha is not None and out.get("sha256") != sha:
            p.append(f"sha256 {out.get('sha256')} != {sha}")
        return self._record("accept", kind, label, p, raw + err), res

    def reject(self, kind, label, path, opts=(), violation=None, words=None,
               first_line=None, d=None, pair=None, res=None):
        res = res or run_bin(path, opts)
        rc, out, raw, err, _ = res
        p = []
        if rc != 1:
            p.append(f"exit {rc}")
        if out.get("valid") != "False":
            p.append(f"valid={out.get('valid')}")
        if violation is not None and out.get("violation") != violation:
            p.append(f"violation={out.get('violation')} expected {violation}")
        if not p and violation == "distance" and words is not None:
            try:
                li, lj = map(int, out["first_violating_lines"].split())
                i, j = li - first_line, lj - first_line
                dist = int(out["first_violating_distance"])
                wi, wj = out["first_violating_pair"].split()
                ref = lev(words[i], words[j])
                if (wi, wj) != (words[i], words[j]):
                    p.append("printed pair does not match its line numbers")
                if ref != dist:
                    p.append(f"reported distance {dist} != reference {ref}")
                if not ref < d:
                    p.append(f"reported pair has reference distance {ref} >= d={d}")
                if pair is not None and (i, j) != tuple(pair):
                    p.append(f"first pair {(i, j)} != reference {tuple(pair)}")
            except (KeyError, ValueError, IndexError) as e:
                p.append(f"could not parse violation: {e!r}")
        return self._record("reject", kind, label, p, raw + err), res

    def refuse(self, kind, label, path, opts=()):
        rc, out, raw, err, _ = run_bin(path, opts)
        p = []
        if rc != 2:
            p.append(f"exit {rc}")
        if "valid" in out:
            p.append(f"printed a verdict valid={out['valid']}")
        return self._record("refuse", kind, label, p, raw + err)

    def other(self, kind, label, ok, detail=""):
        return self._record("other", kind, label, [] if ok else [detail])


# --------------------------------------------------------------------------
# corruptions of an accepted code
# --------------------------------------------------------------------------

def substitute(word, positions, rng):
    w = list(word)
    for p in positions:
        w[p] = rng.choice([c for c in ALPHABET if c != w[p]])
    return "".join(w)


def indel(word, rng):
    """Delete one symbol and insert one (length preserved): <= 2 edits."""
    n = len(word)
    p = rng.randrange(n)
    rest = word[:p] + word[p + 1:]
    r = rng.randrange(n)
    return rest[:r] + rng.choice(ALPHABET) + rest[r:]


def near_codeword(words, n, d, rng, mode):
    """Replace one codeword b by a word within distance d-1 of codeword a."""
    present = set(words)
    for _ in range(2000):
        a, b = rng.sample(range(len(words)), 2)
        if mode == "sub":
            if not 1 <= d - 1 <= n:
                return None
            new = substitute(words[a], rng.sample(range(n), d - 1), rng)
        else:  # deletion + insertion, then d-3 substitutions; Hamming >= d
            if d < 3 or d - 3 > n:
                return None
            new = indel(words[a], rng)
            if d > 3:
                new = substitute(new, rng.sample(range(n), d - 3), rng)
            if hamming(words[a], new) < d:
                continue
        if new in present:
            continue
        dist = lev(words[a], new)
        assert 1 <= dist <= d - 1, (words[a], new, dist, d)
        out = list(words)
        out[b] = new
        return out
    return None


def reject_battery(S, tmp, tag, words, n, d, rng):
    """All single-corruption rejects for one accepted code."""
    if len(words) >= 2:
        for mode in ("sub", "indel"):
            bad = near_codeword(words, n, d, rng, mode)
            if bad is None:
                continue
            path, fl = write_code(tmp, f"{tag}_near_{mode}.txt", bad, n=n, d=d)
            S.reject(f"near_{mode}", tag, path, violation="distance",
                     words=bad, first_line=fl, d=d)
        a, b = rng.sample(range(len(words)), 2)
        dup = list(words)
        dup[b] = words[a]
        path, _ = write_code(tmp, f"{tag}_dup_replace.txt", dup, n=n, d=d)
        S.reject("duplicate_replace", tag, path, violation="duplicate")
    k = rng.randrange(len(words))
    path, _ = write_code(tmp, f"{tag}_dup_append.txt", list(words) + [words[k]], n=n, d=d)
    S.reject("duplicate_append", tag, path, violation="duplicate")

    for kind, new in (("short", words[k][:-1]),
                      ("long", words[k] + rng.choice(ALPHABET))):
        bad = list(words)
        bad[k] = new
        path, _ = write_code(tmp, f"{tag}_{kind}.txt", bad, n=n, d=d)
        S.reject(f"length_{kind}", tag, path, violation="length")

    for kind, sym in (("4", "4"), ("x", "x"), ("space", " ")):
        p = rng.randrange(max(1, n - 1))  # never the last position
        bad = list(words)
        bad[k] = words[k][:p] + sym + words[k][p + 1:]
        path, _ = write_code(tmp, f"{tag}_sym_{kind}.txt", bad, n=n, d=d)
        S.reject(f"symbol_{kind}", tag, path, violation="symbol")


def d_plus_one(S, tmp, tag, words, n, dmin, exact_pair=True):
    """Code with true minimum distance exactly dmin, checked with dmin+1."""
    pair = ref_first_violation(words, dmin + 1) if exact_pair else None
    path, fl = write_code(tmp, f"{tag}_dplus1_hdr.txt", words, n=n, d=dmin + 1)
    S.reject("d_plus_1_header", tag, path, violation="distance", words=words,
             first_line=fl, d=dmin + 1, pair=pair)
    path, fl = write_code(tmp, f"{tag}_dplus1_opt.txt", words, n=n)
    S.reject("d_plus_1_option", tag, path, opts=("--d", dmin + 1),
             violation="distance", words=words, first_line=fl, d=dmin + 1, pair=pair)


# --------------------------------------------------------------------------
# test sections
# --------------------------------------------------------------------------

def load_artifact(path: Path):
    lines = path.read_text().splitlines()
    header = [l for l in lines if l.startswith("#")]
    words = [l.rstrip() for l in lines if l.strip() and not l.startswith("#")]
    hdr = {}
    for l in header:
        m = re.match(r"#\s*([ndq])\s*=\s*(\d+)", l)
        if m:
            hdr[m.group(1)] = int(m.group(2))
    m = re.search(r"sha256\(codewords\)\s*=\s*([0-9a-f]{64})", "\n".join(header))
    return words, hdr, (m.group(1) if m else None)


def section_best(S, tmp, rng):
    files = sorted((ROOT / "best").glob("edit4_*.txt")) + sorted((ROOT / "best").glob("editq_*.txt"))
    print(f"[best] {len(files)} artifacts")
    for f in files:
        words, hdr, hsha = load_artifact(f)
        tag = f.stem
        opts = []
        if "n" in hdr:
            n = hdr["n"]
        else:
            n = int(re.search(r"_n(\d+)_", f.name).group(1))
            opts += ["--n", n]
            S.refuse("missing_n_real", tag, f, opts=("--d", hdr.get("d", 3)))
        if "d" in hdr:
            d = hdr["d"]
        else:
            assert f.name.startswith("edit4_"), f
            d = 3
            S.refuse("missing_d_real", tag, f, opts=opts)  # no --d: must refuse
            opts += ["--d", 3]
        _, res = S.accept("best", tag, f, opts=opts, size=len(words), sha=sha_of(words))
        if hsha is not None:
            got = res[1].get("sha256")
            S.other("header_sha", tag, got == hsha, f"sha256 {got} != header {hsha}")

        # reference DP: true minimum on small artifacts, a witness at d on big ones
        if len(words) <= 400:
            dmin = ref_min_distance(words)
            S.other("reference_confirms_valid", tag, dmin >= d,
                    f"reference minimum distance {dmin} < d={d}")
            d_plus_one(S, tmp, tag, words, n, dmin, exact_pair=True)
        else:
            witness = next(((i, j) for i in range(min(len(words), 40))
                            for j in range(i + 1, len(words))
                            if lev(words[i], words[j]) == d), None)
            if witness is not None:
                d_plus_one(S, tmp, tag, words, n, d, exact_pair=False)
        reject_battery(S, tmp, tag, words, n, d, rng)


def greedy_code(n, d, rng, max_size, max_tries):
    code, seen = [], set()
    for _ in range(max_tries):
        w = "".join(rng.choice(ALPHABET) for _ in range(n))
        if w in seen:
            continue
        seen.add(w)
        if all(lev(w, c) >= d for c in code):
            code.append(w)
            if len(code) >= max_size:
                break
    return code


def section_bruteforce(S, tmp, rng):
    grid = [(3, 2, 12), (4, 2, 30), (4, 3, 12), (5, 3, 25), (6, 3, 40), (6, 4, 15),
            (7, 4, 30), (8, 5, 20), (10, 4, 60), (10, 6, 20), (12, 3, 120),
            (12, 4, 100), (12, 5, 80), (12, 7, 25), (12, 9, 8), (16, 6, 40), (16, 9, 15)]
    print(f"[bruteforce] {len(grid)} greedy codes")
    for n, d, size in grid:
        words = greedy_code(n, d, rng, size, 20 * size + 200)
        dmin = ref_min_distance(words)
        tag = f"greedy_n{n}_d{d}_m{len(words)}"
        assert dmin >= d, tag
        path, fl = write_code(tmp, f"{tag}.txt", words, n=n, d=dmin)
        S.accept("bruteforce_at_dmin", tag, path, size=len(words), sha=sha_of(words))
        if dmin > 1:
            path_lo, _ = write_code(tmp, f"{tag}_dm1.txt", words, n=n, d=dmin - 1)
            S.accept("bruteforce_below_dmin", tag, path_lo, size=len(words))
        d_plus_one(S, tmp, tag, words, n, dmin, exact_pair=True)
        reject_battery(S, tmp, tag, words, n, dmin, rng)


def random_set(n, size, rng):
    words, seen = [], set()
    while len(words) < size:
        if words and rng.random() < 0.6:
            base = rng.choice(words)
            w = base
            for _ in range(rng.randint(1, 3)):
                w = indel(w, rng) if rng.random() < 0.5 else substitute(w, [rng.randrange(n)], rng)
        else:
            w = "".join(rng.choice(ALPHABET) for _ in range(n))
        if w not in seen:
            seen.add(w)
            words.append(w)
    return words


def section_random(S, tmp, rng, trials=400):
    print(f"[random] {trials} differential sets")
    for t in range(trials):
        n = rng.choice([2, 3, 4, 5, 6, 8, 12, 16])
        size = rng.randint(2, min(40, 4 ** n))
        words = random_set(n, size, rng)
        d = rng.randint(1, n + 1)
        use_header = rng.random() < 0.5
        path, fl = write_code(tmp, f"rand_{t}.txt", words, n=n, d=d if use_header else None)
        opts = () if use_header else ("--d", d)
        ref = ref_first_violation(words, d)
        tag = f"rand{t}_n{n}_m{size}_d{d}"
        if ref is None:
            S.accept("random_valid", tag, path, opts=opts, size=size, sha=sha_of(words))
        else:
            S.reject("random_invalid", tag, path, opts=opts, violation="distance",
                     words=words, first_line=fl, d=d, pair=ref)
        if t % 3 == 0:  # output must not depend on the thread count
            r1 = run_bin(path, list(opts) + ["--threads", "1"])
            r2 = run_bin(path, list(opts) + ["--threads", "2"])
            S.other("threads_1_vs_2", tag, (r1[0], r1[2]) == (r2[0], r2[2]),
                    "stdout or exit code differs between 1 and 2 threads")


def section_format(S, tmp, rng):
    words = greedy_code(12, 5, rng, 30, 400)
    dmin = ref_min_distance(words)
    sha = sha_of(words)
    # CRLF, blank lines, trailing whitespace, no final newline
    path, _ = write_code(tmp, "fmt_crlf.txt", words, n=12, d=dmin, newline="\r\n")
    S.accept("format_crlf", "crlf", path, size=len(words), sha=sha)
    spaced = []
    for w in words:
        spaced += [w + "  \t", ""]
    path, _ = write_code(tmp, "fmt_blank.txt", spaced, n=12, d=dmin)
    S.accept("format_blank_trailing", "blank lines + trailing whitespace", path,
             size=len(words), sha=sha)
    path, _ = write_code(tmp, "fmt_nofinal.txt", words, n=12, d=dmin, final_newline=False)
    S.accept("format_no_final_newline", "no final newline", path, size=len(words), sha=sha)
    path, _ = write_code(tmp, "fmt_hdr_style.txt", words,
                         extra=[f"# n = 12, alphabet = 0123, size = {len(words)}",
                                f"# d = {dmin}", "# q = 4"])
    S.accept("format_header_style", "'# n = 12, alphabet...' header", path,
             size=len(words), sha=sha)
    path, _ = write_code(tmp, "fmt_single.txt", words[:1], n=12, d=dmin)
    S.accept("format_single_word", "single codeword", path, size=1, sha=sha_of(words[:1]))
    # leading whitespace is not stripped: ' ' + 11 symbols is a bad symbol,
    # ' ' + 12 symbols a bad length (length is checked first)
    path, _ = write_code(tmp, "fmt_indented.txt", [" " + words[0][1:]] + words[1:], n=12, d=dmin)
    S.reject("format_leading_space", "leading space, length n", path, violation="symbol")
    path, _ = write_code(tmp, "fmt_indented13.txt", [" " + words[0]] + words[1:], n=12, d=dmin)
    S.reject("format_leading_space_long", "leading space, length n+1", path, violation="length")
    return words, dmin


def section_refuse(S, tmp, words, dmin):
    n = 12
    cases = [
        ("missing_d", dict(n=n), ()),
        ("missing_n", dict(d=dmin), ()),
        ("missing_n_and_d", dict(), ()),
        ("conflict_d_option", dict(n=n, d=dmin), ("--d", dmin + 1)),
        ("conflict_n_option", dict(n=n, d=dmin), ("--n", n - 1)),
        ("conflict_d_headers", dict(n=n, d=dmin, extra=[f"# d = {dmin + 1}"]), ()),
        ("conflict_n_headers", dict(n=n, d=dmin, extra=["# n = 11"]), ()),
        ("q_not_4", dict(n=n, d=dmin, extra=["# q = 5"]), ()),
        ("d_zero_header", dict(n=n, d=0), ()),
        ("d_zero_option", dict(n=n), ("--d", 0)),
        ("n_zero", dict(n=0, d=dmin), ()),
        ("n_too_large", dict(n=33, d=dmin), ()),
        ("d_not_integer", dict(n=n, extra=["# d = three"]), ()),
        ("n_trailing_junk", dict(d=dmin, extra=["# n = 12x"]), ()),
        ("d_negative_option", dict(n=n), ("--d", -3)),
        ("d_option_no_value", dict(n=n), ("--d",)),
        ("threads_3", dict(n=n, d=dmin), ("--threads", 3)),
        ("unknown_option", dict(n=n, d=dmin), ("--bogus",)),
    ]
    for kind, kw, opts in cases:
        path, _ = write_code(tmp, f"refuse_{kind}.txt", words, **kw)
        S.refuse(kind, kind, path, opts=opts)
    S.refuse("missing_file", "missing file", tmp / "does_not_exist.txt", opts=("--d", 3))
    path, _ = write_code(tmp, "refuse_two_paths.txt", words, n=n, d=dmin)
    S.refuse("two_paths", "two positional arguments", path, opts=(str(path),))
    # a '# d' header equal to --d is not a conflict
    S.accept("header_equals_option", "'# d' equal to --d", path, opts=("--d", dmin))


def random_words_file(tmp, name, count, d, seed):
    rng = random.Random(seed)
    words = ["".join(ALPHABET[(x >> (2 * k)) & 3] for k in range(12))
             for x in rng.sample(range(4 ** 12), count)]
    return write_code(tmp, name, words, n=12, d=d), words


def power_state():
    """Best-effort record of the power state (macOS pmset); for the log only."""
    info = {}
    try:
        b = subprocess.run(["pmset", "-g", "batt"], capture_output=True, text=True,
                           timeout=10).stdout
        info["source"] = b.splitlines()[0].strip() if b else None
        g = subprocess.run(["pmset", "-g"], capture_output=True, text=True, timeout=10).stdout
        m = re.search(r"lowpowermode\s+(\d)", g)
        info["lowpowermode"] = int(m.group(1)) if m else None
    except (OSError, subprocess.SubprocessError, IndexError):
        pass
    return info


def timing_entry(pairs, st):
    cpu, wall = st["cpu_seconds"], st["wall_seconds"]
    return {**st,
            "pairs_per_wall_second": pairs / wall,
            "pairs_per_cpu_second": pairs / cpu if cpu > 0 else None,
            "cpu_over_wall": cpu / wall}


def section_scaling(S, tmp, run_200k):
    sc = {"seed": SCALE_SEED, "n": 12, "threads": 2,
          "load_avg_before": os.getloadavg(), "power_before": power_state(),
          "note": ("wall_seconds is time.perf_counter() around the subprocess (on macOS it "
                   "does not advance while the system sleeps); cpu_seconds is the child's "
                   "user+sys time summed over both threads.  A 2-thread wall-clock rate is "
                   "trusted only when cpu/wall >= 1.8.")}
    print("[scaling] 100,000 random words, d = 3", flush=True)
    (path, fl), words = random_words_file(tmp, "scale_100k_d3.txt", 100_000, 3, SCALE_SEED)
    res = run_bin(path)
    ok, _ = S.reject("scaling_random_100k_d3", "100k random words, d=3", path,
                     violation="distance", words=words, first_line=fl, d=3, res=res)
    rc, out, raw, err, st = res
    sc["random_100k_d3"] = {"exit": rc, "rejected_as_expected": ok,
                            "first_violating_lines": out.get("first_violating_lines"),
                            "first_violating_distance": out.get("first_violating_distance"),
                            **st}
    full = {}
    sizes = [20_000, 50_000, 100_000] + ([200_000] if run_200k else [])
    for count in sizes:
        print(f"[scaling] full all-pairs run, {count} words, d = 1", flush=True)
        (path, _), words = random_words_file(tmp, f"scale_{count}_d1.txt", count, 1,
                                             SCALE_SEED + count)
        res = run_bin(path)
        S.accept(f"scaling_full_{count}", f"{count} random words, d=1", path,
                 size=count, sha=sha_of(words), res=res)
        rc, out, raw, err, st = res
        pairs = count * (count - 1) // 2
        entry = {"words": count, "pairs": pairs, "exit": rc, **timing_entry(pairs, st)}
        entry["wall_reliable"] = entry["cpu_over_wall"] >= 1.8
        if count == 20_000:
            res1 = run_bin(path, ["--threads", "1"])
            S.other("scaling_threads_1_vs_2", "20k d=1",
                    (res1[0], res1[2]) == (rc, raw), "1-thread output differs")
            entry["one_thread"] = timing_entry(pairs, res1[4])
        full[str(count)] = entry
    sc["full_runs_d1"] = full
    largest = full[str(max(sizes))]
    rate_cpu2 = 2 * largest["pairs_per_cpu_second"]
    rate_wall = largest["pairs_per_wall_second"] if largest["wall_reliable"] else None
    sc["extrapolation"] = {
        "from_words": largest["words"],
        "two_thread_rate_from_cpu_time": rate_cpu2,
        "two_thread_rate_from_wall_time": rate_wall,
        "assumption": "2-thread rate = 2 x pairs per CPU-second (perfect 2-way split)",
        "seconds": {str(N): {"from_cpu_rate": (N * (N - 1) // 2) / rate_cpu2,
                             "from_wall_rate": ((N * (N - 1) // 2) / rate_wall) if rate_wall else None}
                    for N in (100_000, 200_000)},
    }
    sc["load_avg_after"] = os.getloadavg()
    sc["power_after"] = power_state()
    return sc


# --------------------------------------------------------------------------

def main():
    args = set(sys.argv[1:])
    unknown = args - {"--no-scaling", "--scaling-200k"}
    if unknown:
        sys.exit(f"unknown arguments: {sorted(unknown)}")
    t_start = time.perf_counter()
    build()
    S = Suite()
    rng = random.Random(SEED)

    # the reference DP itself, on hand-checked values
    for a, b, want in [("kitten", "sitting", 3), ("", "012", 3), ("0123", "1230", 2),
                       ("0000", "1111", 4), ("0101", "1010", 2), ("0123", "0123", 0),
                       ("012301", "123012", 2), ("000111", "111000", 6)]:
        S.other("reference_dp", f"lev({a!r},{b!r})", lev(a, b) == want,
                f"reference lev gives {lev(a, b)}, expected {want}")

    p = subprocess.run([str(BIN), "--selftest"], capture_output=True, text=True)
    print(p.stdout.strip())
    S.other("selftest", "verify_big2 --selftest", p.returncode == 0, p.stdout[-300:])

    timings = {}
    with tempfile.TemporaryDirectory(prefix="verify_big2_") as td:
        tmp = Path(td)
        for name, fn in (("best", lambda: section_best(S, tmp, rng)),
                         ("bruteforce", lambda: section_bruteforce(S, tmp, rng)),
                         ("random", lambda: section_random(S, tmp, rng)),
                         ("format_and_refuse", lambda: section_refuse(S, tmp, *section_format(S, tmp, rng)))):
            t0 = time.perf_counter()
            fn()
            timings[name] = time.perf_counter() - t0
        scaling = None
        if "--no-scaling" not in args:
            t0 = time.perf_counter()
            scaling = section_scaling(S, tmp, "--scaling-200k" in args)
            timings["scaling"] = time.perf_counter() - t0

    total = time.perf_counter() - t_start
    print()
    print("=" * 70)
    for cat in ("accept", "reject", "refuse", "other"):
        print(f"  {cat:7s}: {S.counts[cat]}")
    for k in sorted(S.kinds):
        print(f"      {k}: {S.kinds[k]}")
    if scaling:
        r = scaling["random_100k_d3"]
        print(f"  100k random words, d=3: exit {r['exit']}, wall {r['wall_seconds']:.3f} s, "
              f"cpu {r['cpu_seconds']:.3f} s (pair at lines {r['first_violating_lines']}, "
              f"distance {r['first_violating_distance']})")
        for k, e in scaling["full_runs_d1"].items():
            print(f"  full scan {k} words: wall {e['wall_seconds']:.2f} s, cpu {e['cpu_seconds']:.2f} s, "
                  f"{e['pairs_per_cpu_second']:.3e} pairs/cpu-s, {e['pairs_per_wall_second']:.3e} "
                  f"pairs/wall-s (wall reliable: {e['wall_reliable']})")
        for N, x in scaling["extrapolation"]["seconds"].items():
            print(f"  extrapolated full scan {N} words (2 threads): {x['from_cpu_rate']:.0f} s from "
                  f"CPU rate; from wall rate: {x['from_wall_rate']}")
    print(f"  failures: {len(S.failures)}   total wall time {total:.1f} s")
    print("PASS" if not S.failures else "FAIL")

    try:
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True,
                              text=True).stdout.strip()
    except OSError:
        head = None
    log = {
        "suite": "tests/test_verify_big2.py",
        "date_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "git_head": head,
        "verify_big2_c_sha256": hashlib.sha256(SRC.read_bytes()).hexdigest(),
        "build": " ".join(BUILD_CMD),
        "machine": {"platform": platform.platform(), "machine": platform.machine(),
                    "python": platform.python_version()},
        "seed": SEED,
        "counts": dict(S.counts),
        "counts_by_kind": dict(sorted(S.kinds.items())),
        "failures": S.failures,
        "passed": not S.failures,
        "section_seconds": timings,
        "total_seconds": total,
        "scaling": scaling,
    }
    LOG.parent.mkdir(exist_ok=True)
    LOG.write_text(json.dumps(log, indent=2) + "\n")
    print(f"log: {LOG.relative_to(ROOT)}")
    return 0 if not S.failures else 1


if __name__ == "__main__":
    sys.exit(main())
