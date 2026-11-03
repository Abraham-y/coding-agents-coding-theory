#!/usr/bin/env bash
# Every verification in the repo, in one command. Exits non-zero if anything fails.
#
#   ./run_tests.sh          fast suites only (~1 min)
#   ./run_tests.sh --full   adds the slow ones (~15 min)
set -u

# Interpreter: prefer the repo venv, fall back to system python3 so the suite
# runs from a clean checkout (CLAUDE.md rule 5) and inside the anonymised
# supplementary bundle, which ships no venv.
PY=${PY:-./.venv/bin/python}
[ -x "$PY" ] || PY=$(command -v python3)
[ -n "$PY" ] || { echo "no python found"; exit 1; }
FAILED=0
FULL=0
[ "${1:-}" = "--full" ] && FULL=1

run() {
  local label="$1"; shift
  echo "----------------------------------------------------------------------"
  echo "== $label"
  if "$@"; then
    echo "   PASS"
  else
    echo "   *** FAIL: $label"
    FAILED=1
  fi
}

echo "======================================================================"
echo " Track A -- binary s-deletion-correcting codes"
echo "======================================================================"
run "verify.py selftest"                 $PY verify.py --selftest
run "verify2.py selftest"                $PY verify2.py --selftest
run "two-directional suite (accept + reject sets)" $PY tests/test_verify.py
run "differential agreement verify vs verify2"     $PY tests/test_agreement.py

echo
echo "== certified LP upper bounds: re-check every stored certificate"
for f in results/cert_n*.json; do
  [ -e "$f" ] || continue
  if ! $PY search/lp_bound.py --verify-cert "$f" | grep -q '^PASS'; then
    echo "   *** FAIL: certificate $f"
    FAILED=1
  fi
done
echo "   re-checked $(ls results/cert_n*.json 2>/dev/null | wc -l | tr -d ' ') certificates"

echo
echo "======================================================================"
echo " Track B -- recursive 3-majority, MNSX constant alpha_k"
echo "======================================================================"
run "alpha_k structure + published alpha_1, alpha_2" $PY tests/test_maj3.py
run "three-solver agreement + Rule A/B premises"     $PY tests/test_maj3_agreement.py
run "alpha_1 == 2         (forced-rule solver)"      $PY search/maj3_fast.py --k 1 --quiet
run "alpha_2 == 24/7      (forced-rule solver)"      $PY search/maj3_fast.py --k 2 --quiet
run "alpha_3 == 12231/2203 (forced-rule solver)"     $PY search/maj3_fast.py --k 3 --quiet
run "alpha_1..3 (MNSX Appendix A transcription)"     bash -c "$PY search/maj3_mnsx.py --k 1 >/dev/null && $PY search/maj3_mnsx.py --k 2 >/dev/null && $PY search/maj3_mnsx.py --k 3 >/dev/null"
run "A_2(6,7,8,9) exact via two independent methods"  $PY tests/test_exact_cells.py

echo
echo "======================================================================"
echo " Lovasz theta (SDP) bounds"
echo "======================================================================"
run "theta vs known values (K_n, empty, C_5, C_7, Petersen) + sandwich" $PY tests/test_theta.py
run "edit-metric verifier: accept/reject + differential Levenshtein" $PY tests/test_edit4_verify.py
run "GROUND-TRUTH edit verifiers: accept/reject/agreement" \
    "$PY" tests/test_groundtruth_edit.py
run "general-d edit verifiers (verify_editq, verify_editq2): accept/reject/agreement" $PY tests/test_editq_verify.py
run "q = 2, 3 edit artifacts: verify.py vs verify_editq vs verify_editq2" $PY tests/test_agreement_q.py

if [ "$FULL" = "1" ]; then
  echo
  echo "== slow suites (--full)"
  run "two-directional suite, n up to 20"  $PY tests/test_verify.py --full
  run "differential agreement, extended"   $PY tests/test_agreement.py --full
  run "alpha_3 via structured solver (no forced rules, ~8 min)" \
      $PY search/maj3_struct.py --k 3
  run "alpha_4 == 2027349/216164 (MNSX transcription, ~3 min)" \
      $PY search/maj3_mnsx.py --k 4
fi

echo "======================================================================"
if [ "$FAILED" = "0" ]; then
  echo " ALL CHECKS PASSED"
else
  echo " SOMETHING FAILED -- see above"
fi
echo "======================================================================"
exit $FAILED
