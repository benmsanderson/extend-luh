#!/usr/bin/env python3
"""
Driver: calibrate, extend, and verify all scenarios.

Usage:
    python run_all.py                    # run all scenarios in SCENARIOS dict
    python run_all.py VL H               # run only VL and H
    python run_all.py --skip-extension   # calibrate only (fast, no NetCDF I/O)
    python run_all.py --skip-verify      # skip verification step
"""

import argparse
import sys
import time

sys.path.insert(0, ".")

from src import config as cfg
from src.pipeline import calibrate_scenario, extend_scenario, verify_scenario
from src.diagnostics import plot_diagnostics


def main():
    parser = argparse.ArgumentParser(description="Run extend-LUH pipeline for all scenarios")
    parser.add_argument("scenarios", nargs="*", default=list(cfg.SCENARIOS.keys()),
                        help="Scenario keys to process (default: all)")
    parser.add_argument("--skip-calibration", action="store_true",
                        help="Skip calibration (use existing .npz files)")
    parser.add_argument("--skip-extension", action="store_true",
                        help="Skip gridded extension (calibrate only)")
    parser.add_argument("--skip-plots", action="store_true",
                        help="Skip diagnostic plot generation")
    parser.add_argument("--skip-verify", action="store_true",
                        help="Skip verification step")
    parser.add_argument("-q", "--quiet", action="store_true",
                        help="Suppress per-step logging")
    args = parser.parse_args()

    verbose = not args.quiet
    scenarios = args.scenarios

    # Validate scenario keys
    for key in scenarios:
        if key not in cfg.SCENARIOS:
            print(f"ERROR: unknown scenario '{key}'. "
                  f"Available: {list(cfg.SCENARIOS.keys())}")
            sys.exit(1)

    print(f"{'=' * 70}")
    print(f"extend-LUH pipeline — {len(scenarios)} scenario(s): {scenarios}")
    print(f"{'=' * 70}\n")

    cal_results = {}
    ext_results = {}
    ver_results = {}

    for key in scenarios:
        sc = cfg.get_scenario(key)
        t0 = time.time()
        print(f"── {key}: {sc.model} / {sc.scenario} ──")

        # Step 1: Calibration
        if not args.skip_calibration:
            print(f"  [1/4] Calibrating...")
            try:
                cal_results[key] = calibrate_scenario(sc, verbose=verbose)
            except Exception as e:
                print(f"  ✗ Calibration failed: {e}")
                cal_results[key] = {"key": key, "status": f"error: {e}"}
                continue
        else:
            print(f"  [1/4] Skipped (using existing calibration)")
            cal_results[key] = {"key": key, "status": "skipped"}

        # Step 2: Gridded extension
        if not args.skip_extension:
            print(f"  [2/4] Extending...")
            try:
                ext_results[key] = extend_scenario(sc, verbose=verbose)
            except Exception as e:
                print(f"  ✗ Extension failed: {e}")
                ext_results[key] = {"key": key, "status": f"error: {e}"}
                continue
        else:
            print(f"  [2/4] Skipped")
            ext_results[key] = {"key": key, "status": "skipped"}

        # Step 3: Diagnostic plots
        if not args.skip_plots and not args.skip_extension:
            print(f"  [3/4] Plots...")
            try:
                plot_diagnostics(sc, verbose=verbose)
            except Exception as e:
                print(f"  ✗ Diagnostic plots failed: {e}")
        else:
            print(f"  [3/4] Skipped")

        # Step 4: Verification
        if not args.skip_verify and not args.skip_extension:
            print(f"  [4/4] Verifying...")
            try:
                ver_results[key] = verify_scenario(sc, verbose=verbose)
            except Exception as e:
                print(f"  ✗ Verification failed: {e}")
                ver_results[key] = {"key": key, "status": f"error: {e}"}
        else:
            print(f"  [4/4] Skipped")
            ver_results[key] = {"key": key, "status": "skipped"}

        elapsed = time.time() - t0
        print(f"  Done ({elapsed:.0f}s)\n")

    # ═════════════════════════════════════════════════════════════════
    # Summary table
    # ═════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print(f"{'=' * 70}\n")

    # Calibration summary
    print("Calibration:")
    hdr = f"{'Scen':>5s}  {'Model':>12s}  {'τ':>4s}  {'R²':>6s}  {'RMSE_cal':>8s}  {'RMSE_ramp':>9s}  {'r(2101)':>7s}  {'r→0':>5s}  {'Status'}"
    print(hdr)
    print("-" * len(hdr))
    for key in scenarios:
        r = cal_results.get(key, {})
        if r.get("status") == "ok":
            print(f"{key:>5s}  {r['model_type']:>12s}  {r['tau']:>4.0f}  "
                  f"{r['r2']:>6.4f}  {r['rmse_cal']:>8.1f}  {r['rmse_ramp']:>9.1f}  "
                  f"{r['r_2101']:>7.4f}  {r.get('ramp_zero_yr', ''):>5}  ok")
        else:
            print(f"{key:>5s}  {'':>12s}  {'':>4s}  {'':>6s}  {'':>8s}  {'':>9s}  "
                  f"{'':>7s}  {'':>5s}  {r.get('status', '?')}")

    # Extension summary
    if ext_results:
        print(f"\nExtension:")
        hdr = f"{'Scen':>5s}  {'Ramp period':>14s}  {'Files':>5s}  {'Total MB':>8s}  {'Status'}"
        print(hdr)
        print("-" * len(hdr))
        for key in scenarios:
            r = ext_results.get(key, {})
            if r.get("status") == "ok":
                print(f"{key:>5s}  {r['ramp_years']:>14s}  {len(r['files']):>5d}  "
                      f"{r['total_mb']:>8.1f}  ok")
            else:
                print(f"{key:>5s}  {'':>14s}  {'':>5s}  {'':>8s}  {r.get('status', '?')}")

    # Verification summary
    if any(r.get("status") == "ok" for r in ver_results.values()):
        print(f"\nVerification (post-2100 vs IAM):")
        hdr = f"{'Scen':>5s}  {'RMSE':>8s}  {'Jump@2100':>10s}  {'Max|resid|':>10s}  {'Status'}"
        print(hdr)
        print("-" * len(hdr))
        for key in scenarios:
            r = ver_results.get(key, {})
            if r.get("status") == "ok":
                rmse = r["rmse"]
                flag = "WARN" if rmse > 100 else "ok"
                print(f"{key:>5s}  {rmse:>8.1f}  {r['jump_2100']:>+10.1f}  "
                      f"{r['max_residual']:>10.1f}  {flag}")
            else:
                print(f"{key:>5s}  {'':>8s}  {'':>10s}  {'':>10s}  {r.get('status', '?')}")

    print()


if __name__ == "__main__":
    main()
