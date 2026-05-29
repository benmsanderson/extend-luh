#!/usr/bin/env python
"""
Run the extend-LUH notebook pipeline for one or more scenarios using papermill.

Usage:
    python run_notebooks.py [SCENARIO_KEY ...] [options]

Examples:
    python run_notebooks.py VL H          # run both scenarios end-to-end
    python run_notebooks.py VL --skip-calibration   # skip notebook 02
    python run_notebooks.py --no-save     # all scenarios, discard executed notebooks
"""

import argparse
import sys
import time
from pathlib import Path

NB_DIR     = Path(__file__).parent / "notebooks"
OUTPUT_DIR = Path(__file__).parent / "output"

NOTEBOOKS = {
    "calibration": NB_DIR / "02_afolu_calibration_multi.ipynb",
    "extension":   NB_DIR / "03_gridded_extension.ipynb",
    "verify":      NB_DIR / "04_verify_afolu.ipynb",
}

STEP_ORDER = ["calibration", "extension", "verify"]


def run_notebook(step: str, scenario_key: str, save: bool) -> tuple[bool, float, str]:
    """Execute one notebook for one scenario. Returns (success, elapsed_s, error_msg)."""
    import papermill as pm

    nb_in = NOTEBOOKS[step]
    if save:
        out_dir = OUTPUT_DIR / scenario_key / "notebooks"
        out_dir.mkdir(parents=True, exist_ok=True)
        nb_out = out_dir / nb_in.name
    else:
        nb_out = "/dev/null"

    t0 = time.monotonic()
    try:
        pm.execute_notebook(
            str(nb_in),
            str(nb_out),
            parameters={"SCENARIO_KEY": scenario_key},
            kernel_name="python3",
            progress_bar=False,
        )
        return True, time.monotonic() - t0, ""
    except Exception as exc:
        return False, time.monotonic() - t0, str(exc).splitlines()[0]


def main():
    import src.config as cfg

    all_keys = list(cfg.SCENARIOS.keys())

    parser = argparse.ArgumentParser(
        description="Run the extend-LUH notebook pipeline via papermill.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"Available scenarios: {', '.join(all_keys)}",
    )
    parser.add_argument(
        "scenarios",
        nargs="*",
        default=all_keys,
        help="Scenario keys to process (default: all)",
    )
    parser.add_argument("--skip-calibration", action="store_true",
                        help="Skip notebook 02 (AFOLU calibration)")
    parser.add_argument("--skip-extension",   action="store_true",
                        help="Skip notebook 03 (gridded extension)")
    parser.add_argument("--skip-verify",      action="store_true",
                        help="Skip notebook 04 (verification)")
    parser.add_argument("--no-save", action="store_true",
                        help="Discard executed notebooks (faster; no audit trail)")
    args = parser.parse_args()

    # Validate scenario keys
    bad = [k for k in args.scenarios if k not in all_keys]
    if bad:
        parser.error(f"Unknown scenario(s): {bad}. Available: {all_keys}")

    skip = {
        "calibration": args.skip_calibration,
        "extension":   args.skip_extension,
        "verify":      args.skip_verify,
    }
    save = not args.no_save

    print(f"\nextend-LUH notebook pipeline")
    print(f"  Scenarios : {args.scenarios}")
    print(f"  Steps     : {[s for s in STEP_ORDER if not skip[s]]}")
    print(f"  Save output notebooks: {save}")
    print()

    results = {}   # {(key, step): (success, elapsed, error)}

    for key in args.scenarios:
        print(f"── {key} ─────────────────────────────────────────────")
        for step in STEP_ORDER:
            if skip[step]:
                continue
            print(f"   {step:14s}  ...", end="", flush=True)
            ok, dt, err = run_notebook(step, key, save)
            status = "OK" if ok else "FAILED"
            print(f"\r   {step:14s}  {status}  ({dt:.0f}s)")
            if err:
                print(f"              {err}")
            results[(key, step)] = (ok, dt, err)

    # Summary table
    print()
    print("─" * 62)
    print(f"{'Scenario':<8}  {'Step':<14}  {'Status':<7}  {'Time':>6}  Note")
    print("─" * 62)
    for (key, step), (ok, dt, err) in results.items():
        status = "OK" if ok else "FAILED"
        note   = err[:30] if err else ("saved" if save else "discarded")
        print(f"{key:<8}  {step:<14}  {status:<7}  {dt:5.0f}s  {note}")
    print("─" * 62)

    n_fail = sum(1 for ok, _, _ in results.values() if not ok)
    if n_fail:
        print(f"\n{n_fail} step(s) failed.")
        sys.exit(1)
    else:
        print(f"\nAll steps completed successfully.")


if __name__ == "__main__":
    main()
