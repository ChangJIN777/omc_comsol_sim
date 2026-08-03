#!/usr/bin/env python3
"""COMSOL connectivity diagnostic -- RUN ON YOUR MAC (not in the sandbox).

Verifies, step by step, that Python can drive COMSOL via MPh:
  1. Python interpreter architecture (must match the COMSOL build).
  2. MPh + JPype import.
  3. COMSOL discovery + license (mph.start()).
  4. A trivial model: set a parameter and evaluate an expression.
  5. (If comsol/omc_unitcell_iso.mph exists) load it and report studies.

    python scripts/check_comsol.py
"""
import os
import platform
import sys

OK, BAD, WARN = "[ OK ]", "[FAIL]", "[warn]"


def main():
    print("=" * 64)
    print("COMSOL <-> Python (MPh) diagnostic")
    print("=" * 64)

    # 1. architecture
    arch = platform.machine()
    print(f"{OK} Python {sys.version.split()[0]}  arch={arch}  ({sys.executable})")
    print("      -> This arch MUST match your COMSOL build:")
    print("         native Apple Silicon COMSOL  -> arm64 Python")
    print("         Intel COMSOL (Rosetta)        -> x86_64 Python")

    # 2. imports
    try:
        import jpype  # noqa: F401
        import mph
        print(f"{OK} MPh {mph.__version__} and JPype import fine")
    except Exception as e:
        print(f"{BAD} cannot import MPh/JPype: {e}")
        print("      fix: pip install MPh")
        return 1

    # 3. start COMSOL
    try:
        client = mph.start()
        print(f"{OK} COMSOL started: version {client.version}")
        try:
            print(f"      install: {client.cores} core(s) available")
        except Exception:
            pass
    except Exception as e:
        print(f"{BAD} mph.start() failed: {e}")
        print("      checks:")
        print("       - COMSOL 6.2 installed under /Applications/COMSOL62/Multiphysics?")
        print("       - license available? (close other COMSOL sessions)")
        print("       - Python arch matches COMSOL build (see step 1)?")
        print("       - try: export COMSOL_HOME=/Applications/COMSOL62/Multiphysics")
        return 1

    # 4. trivial model: set a parameter and read it back via the Java layer
    try:
        model = client.create("diagnostic")
        model.parameter("a", "500[nm]")
        a_val = float(model.java.param().evaluate("a"))
        print(f"{OK} model parameter API works: a = {a_val*1e9:.0f} nm (expect 500 nm)")
        client.remove(model)
    except Exception as e:
        print(f"{WARN} model/eval step had an issue: {e}")
        print("      (non-fatal; proceed if mph.start() passed)")

    # 5. template, if present
    tpl = os.path.join(os.path.dirname(__file__), "..", "comsol",
                       "omc_unitcell_iso.mph")
    if os.path.exists(tpl):
        try:
            m = client.load(tpl)
            print(f"{OK} loaded template: studies = {list(m.studies())}")
        except Exception as e:
            print(f"{WARN} template present but failed to load: {e}")
    else:
        print(f"{WARN} no comsol/omc_unitcell_iso.mph yet "
              f"(build it per comsol/README_template.md)")

    print("\nAll critical steps passed -> you can run:")
    print("   python scripts/run_loop.py --config configs/run_loop.yaml"
          "   # set optical: mpb, mech: comsol in the config")
    return 0


if __name__ == "__main__":
    sys.exit(main())
