"""Top-fill fix step 2 done-when check (decision 0114): rest relaxation
(materialProperties restRelax) removes the impact-locked lateral stress of a
box filled from the top, and leaves a box placed with gravity on as it was.

Runs tests/compaction/fill_box.py (cap on, lambda* 0.02) three times:
  [1] top fill, no relaxation: K mid-deposit above 0.7 (the effect is there)
  [2] top fill, restRelax time 0.05 s: K within 0.35-0.50 (Jaky K0 0.43)
  [3] gravity on, restRelax time 0.05 s: K within 0.35-0.50 (0.42 without)
About 6 minutes.

    python3 relax_checks.py WORKDIR
"""
import os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
work = os.path.abspath(sys.argv[1])


def K(mode, relax):
    env = dict(os.environ)
    env.pop("RELAX", None)
    if relax:
        env["RELAX"] = relax
    out = subprocess.run([sys.executable, f"{HERE}/fill_box.py", work, mode, "cap"],
                         env=env, capture_output=True, text=True)
    if out.returncode:
        print(out.stdout[-2000:], out.stderr[-3000:])
        sys.exit(1)
    line = [l for l in out.stdout.splitlines() if "K mid" in l][-1]
    print("   ", line, flush=True)
    return float(re.search(r"K mid ([-0-9.]+)", line).group(1))


ok = True
for i, (mode, relax, lo, hi, what) in enumerate([
        ("top", None, 0.7, 9.0, "top fill, no relaxation (impact-locked)"),
        ("top", "0.05", 0.35, 0.50, "top fill, restRelax 0.05 s"),
        ("gravity", "0.05", 0.35, 0.50, "gravity on, restRelax 0.05 s")], 1):
    k = K(mode, relax)
    good = lo <= k <= hi
    ok &= good
    print(f"[{i}] {what}: K {k:.2f}, band {lo}-{hi}: {'PASS' if good else 'FAIL'}", flush=True)
print("PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
