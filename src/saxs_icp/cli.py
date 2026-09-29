"""
Command-line interface:

    saxs-icp align     MODEL ENVELOPE [-o OUTDIR]
    saxs-icp pack      PDB_FOLDER ENVELOPE [-o OUTDIR]
    saxs-icp benchmark --base FOLDER --out results.csv
"""
import argparse
import json
import os
import sys
import time
import traceback

import numpy as np

from . import __version__, core
from .io import (count_models, is_mmcif, read_coords, write_points_pdb,
                 write_transformed_cif, write_transformed_pdb)


def _nonneg_float(s):
    v = float(s)
    if not v >= 0:
        raise argparse.ArgumentTypeError(f"must be >= 0, got {s}")
    return v


def _int_at_least(n):
    def check(s):
        v = int(s)
        if v < n:
            raise argparse.ArgumentTypeError(f"must be >= {n}, got {s}")
        return v
    return check


def _existing_file(s):
    if not os.path.isfile(s):
        raise argparse.ArgumentTypeError(f"file not found: {s}")
    return s


def _add_icp_options(p):
    g = p.add_argument_group("alignment options (defaults = published method)")
    g.add_argument("--penalty", type=_nonneg_float, default=core.DEFAULTS["penalty"],
                   help="weight lambda of the leak term (default: 0.2; 0 = plain ICP)")
    g.add_argument("--mode", default=core.DEFAULTS["mode"],
                   choices=["author", "bidir", "uni"],
                   help="cost function (default: author = fill + lambda*leak)")
    g.add_argument("--restarts", type=_int_at_least(1), default=core.DEFAULTS["restarts"],
                   help="initial orientations per hand (default: 3)")
    g.add_argument("--max-iter", type=_int_at_least(1), default=core.DEFAULTS["max_iter"],
                   help="ICP iterations (default: 50)")
    g.add_argument("--max-points", type=_int_at_least(4),
                   default=core.DEFAULTS["max_points"],
                   help="C-alpha atoms used; larger models are subsampled (default: 3000)")
    g.add_argument("--sample-env", type=_int_at_least(0),
                   default=core.DEFAULTS["sample_env"],
                   help="envelope points used, 0 = all (default: 5000)")
    g.add_argument("--seed", type=int, default=core.DEFAULTS["seed"],
                   help="random seed (default: 42)")
    g.add_argument("--no-enantiomorphs", action="store_true",
                   help="do not test the mirror image")


def _icp_kwargs(a):
    return dict(mode=a.mode, penalty=a.penalty, restarts=a.restarts,
                max_iter=a.max_iter, max_points=a.max_points,
                sample_env=a.sample_env, seed=a.seed,
                enantiomorphs=not a.no_enantiomorphs)


def _clean(v):
    """JSON-safe values (NaN -> null, numpy -> python)."""
    if isinstance(v, dict):
        return {k: _clean(x) for k, x in v.items()}
    if isinstance(v, (list, tuple, np.ndarray)):
        return [_clean(x) for x in v]
    if isinstance(v, (np.floating, float)):
        v = float(v)
        return None if v != v else v
    if isinstance(v, np.integer):
        return int(v)
    if isinstance(v, np.bool_):
        return bool(v)
    return v


def cmd_align(a):
    t0 = time.time()
    model_ca = read_coords(a.model, ca_only=True)
    env = read_coords(a.envelope)
    if len(model_ca) < 4:
        raise ValueError(f"only {len(model_ca)} C-alpha atoms read from {a.model} "
                         f"(need at least 4); is it a PDB or mmCIF file?")
    if len(env) < 10:
        raise ValueError(f"no envelope points read from {a.envelope} "
                         f"({len(env)} found, need at least 10)")
    if count_models(a.model) > 1:
        print(f"Warning: {a.model} has several MODEL records; all of them are "
              f"aligned together as a single rigid body.")

    stem = os.path.splitext(os.path.basename(a.model))[0]
    out = a.output or f"{stem}_aligned"
    if os.path.exists(out) and not os.path.isdir(out):
        raise ValueError(f"output path exists and is not a folder: {out}")
    os.makedirs(out, exist_ok=True)

    print(f"Model:    {a.model} ({len(model_ca)} C-alpha)")
    print(f"Envelope: {a.envelope} ({len(env)} points)")
    print(f"Aligning (mode={a.mode}, lambda={a.penalty}, restarts={a.restarts}, "
          f"enantiomorphs={'no' if a.no_enantiomorphs else 'yes'}) ...")

    r = core.align(model_ca, env, **_icp_kwargs(a))
    T_model, M_env = core.proper_transform(r.transform)

    # aligned model, all atoms, same format as the input
    model_cif = is_mmcif(a.model)
    model_out = os.path.join(out, f"{stem}_aligned.{'cif' if model_cif else 'pdb'}")
    try:
        if model_cif:
            write_transformed_cif(a.model, model_out, T_model)
        else:
            write_transformed_pdb(a.model, model_out, T_model)
    except ValueError as e:
        model_out = os.path.join(out, f"{stem}_aligned_CA.pdb")
        print(f"Warning: could not rewrite {a.model} ({e}); "
              f"writing the aligned C-alpha trace instead.")
        write_points_pdb(core.apply_T(model_ca, T_model), model_out)

    env_out = None
    if M_env is not None:
        env_out = os.path.join(out, "envelope_mirrored.pdb")
        write_points_pdb(core.apply_T(env, M_env), env_out, atom="C", res="DUM")

    m = r.metrics
    report = dict(
        saxs_icp_version=__version__,
        model=os.path.abspath(a.model),
        envelope=os.path.abspath(a.envelope),
        parameters=_icp_kwargs(a),
        metrics=dict(chamfer_A=m["chamfer"], hausdorff_A=m["hausdorff"],
                     fraction_outside=m["frac_outside"], coverage=m["coverage"],
                     centroid_separation_A=m["centroid_sep"],
                     median_distance_A=m["median_dist"], dice=m["dice"],
                     cost=r.cost),
        envelope_mirrored=bool(r.mirrored),
        transform_model=T_model,
        transform_envelope=M_env,
        n_model_ca=r.n_model_full, n_model_used=r.n_model,
        n_envelope=len(env), n_envelope_used=r.n_env, grid_A=r.grid,
        outputs=dict(model=model_out, envelope_mirrored=env_out),
        seconds=round(time.time() - t0, 2),
    )
    rep_path = os.path.join(out, "report.json")
    with open(rep_path, "w", encoding="utf-8") as fh:
        json.dump(_clean(report), fh, indent=2)

    print()
    print(f"  Chamfer distance   {m['chamfer']:8.3f} A")
    print(f"  Fraction outside   {m['frac_outside']:8.3f}")
    print(f"  Coverage           {m['coverage']:8.3f}")
    print(f"  Centroid sep.      {m['centroid_sep']:8.3f} A")
    print()
    print(f"Aligned model: {model_out}")
    if env_out:
        print(f"The mirror image fits better. SAXS envelopes do not define "
              f"chirality, so the\nenvelope was reflected instead of the protein: "
              f"view the model with\n{env_out}")
    print(f"Report:        {rep_path}")
    return 0


def cmd_pack(a):
    from .packing import run_packing
    run_packing(a.input, a.envelope, a.output, copies=a.copies, penalty=a.penalty,
                max_iter=a.max_iter, sample_env=a.sample_env, workers=a.workers,
                seed=a.seed)
    return 0


def cmd_benchmark(a):
    from .benchmark import run
    run(a.base, a.out, workers=a.workers, save_aligned=a.save_aligned, limit=a.limit,
        w_leak=a.w_leak, w_cover=a.w_cover, **_icp_kwargs(a))
    return 0


def build_parser():
    p = argparse.ArgumentParser(
        prog="saxs-icp",
        description="Rigid-body fitting of atomic structures into SAXS envelopes "
                    "(ICP with bidirectional penalty and enantiomorph search).")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    p.add_argument("--debug", action="store_true", help="show full tracebacks")
    sub = p.add_subparsers(dest="command", metavar="COMMAND")

    s = sub.add_parser("align", help="fit one model into one envelope",
                       description="Fit one atomic model (PDB or mmCIF) into a SAXS "
                                   "envelope (dummy-atom PDB or mmCIF).")
    s.add_argument("model", type=_existing_file, help="atomic model (.pdb or .cif)")
    s.add_argument("envelope", type=_existing_file,
                   help="envelope, e.g. DAMMIF/DAMMIN/GASBOR output (.pdb or .cif)")
    s.add_argument("-o", "--output", default="",
                   help="output folder (default: <model>_aligned)")
    _add_icp_options(s)
    s.set_defaults(func=cmd_align)

    s = sub.add_parser("pack", help="rank a folder of conformations and pack the best",
                       description="Align every .pdb in a folder to the envelope "
                                   "and keep the best-scoring copies.")
    s.add_argument("input", help="folder with .pdb files")
    s.add_argument("envelope", type=_existing_file, help="envelope (.pdb or .cif)")
    s.add_argument("-o", "--output", default="packing_out",
                   help="output folder (default: packing_out)")
    s.add_argument("--copies", type=_int_at_least(1), default=20,
                   help="number of copies kept (default: 20)")
    s.add_argument("--penalty", type=_nonneg_float, default=0.2,
                   help="weight lambda of the leak term (default: 0.2)")
    s.add_argument("--max-iter", type=_int_at_least(1), default=30,
                   help="ICP iterations (default: 30)")
    s.add_argument("--sample-env", type=_int_at_least(0), default=2000,
                   help="envelope points used, 0 = all (default: 2000)")
    s.add_argument("--workers", type=_int_at_least(0), default=0,
                   help="processes, 0 = all CPUs (default: 0)")
    s.add_argument("--seed", type=int, default=None,
                   help="random seed for reproducible runs (default: none)")
    s.set_defaults(func=cmd_pack)

    s = sub.add_parser("benchmark", help="batch run over a SASBDB-style folder",
                       description="Align every <ACC>_model.pdb to "
                                   "<ACC>_envelope.cif under BASE/data/sasbdb/ "
                                   "and write one CSV row per entry.")
    s.add_argument("--base", required=True, help="folder containing data/sasbdb/")
    s.add_argument("--out", required=True, help="output CSV")
    s.add_argument("--workers", type=_int_at_least(0), default=0,
                   help="processes, 0 = all CPUs (default: 0)")
    s.add_argument("--limit", type=_int_at_least(0), default=0,
                   help="only the first N entries")
    s.add_argument("--save-aligned", default="",
                   help="folder to save the final C-alpha poses")
    s.add_argument("--w-leak", type=_nonneg_float, default=1.0, help=argparse.SUPPRESS)
    s.add_argument("--w-cover", type=_nonneg_float, default=1.0, help=argparse.SUPPRESS)
    _add_icp_options(s)
    s.set_defaults(func=cmd_benchmark)
    return p


def main(argv=None):
    parser = build_parser()
    a = parser.parse_args(argv)
    if not a.command:
        parser.print_help()
        return 1
    try:
        return a.func(a)
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        return 130
    except Exception as e:
        if a.debug:
            traceback.print_exc()
        print(f"saxs-icp: error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
