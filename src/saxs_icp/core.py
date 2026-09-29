"""
Rigid-body alignment of an atomic model into a SAXS envelope by ICP with a
bidirectional cost, C-alpha representation and enantiomorph search.

Modes:
    author  fill + lambda * leak, plain means (the published method; default)
    bidir   thresholded penalty on both terms
    uni     leak-only penalty (original version)
"""
from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree, Delaunay
from scipy.spatial.transform import Rotation

from .io import apply_T

DEFAULTS = dict(mode="author", penalty=0.2, restarts=3, max_iter=50,
                max_points=3000, sample_env=5000, seed=42, enantiomorphs=True,
                w_leak=1.0, w_cover=1.0)

MIRROR_Z = np.diag([1.0, 1.0, -1.0])


def random_rotation(seed):
    """
    Uniform random rotation from an integer seed.

    Same result, bit for bit, as Rotation.random(random_state=seed) (a normal
    quaternion drawn from RandomState(seed)), but independent of the SciPy
    random_state -> rng API transition.
    """
    q = np.random.RandomState(seed).normal(size=4)
    return Rotation.from_quat(q).as_matrix()


def grid_spacing(env):
    """Typical spacing between neighbouring dummy atoms."""
    if len(env) < 2:
        return 3.0
    d, _ = cKDTree(env).query(env, k=2)
    return float(np.median(d[:, 1]))


# cost and ICP

def bidirectional_cost(src, env_tree, src_tree_pts, env_pts, penalty, threshold,
                       w_leak=1.0, w_cover=1.0):
    """
    Thresholded bidirectional cost:

      leak term  : atoms far from the envelope (protein leaking out)
      cover term : envelope points far from any atom (empty space inside)
    """
    d_src, _ = env_tree.query(src, k=1)
    base = float(d_src.mean())
    if penalty <= 0:
        return base, base, 0.0, 0.0

    over = d_src[d_src > threshold]
    leak = float((over - threshold).sum()) / len(src) if len(over) else 0.0

    src_tree = cKDTree(src)
    d_env, _ = src_tree.query(env_pts, k=1)
    unc = d_env[d_env > threshold]
    cover = float((unc - threshold).sum()) / len(env_pts) if len(unc) else 0.0

    total = base + penalty * (w_leak * leak + w_cover * cover)
    return total, base, leak, cover


def author_cost(src, env_pts, env_tree):
    """
    Published formulation:

        score = mean(dist envelope->protein)            # fill
              + penalty * mean(dist protein->envelope)  # leak

    Plain means, no threshold. With penalty = 0.2 each angstrom of leak costs
    1/5 of an angstrom of empty space.
    """
    d_fill, _ = cKDTree(src).query(env_pts, k=1)
    d_leak, _ = env_tree.query(src, k=1)
    return float(d_fill.mean()), float(d_leak.mean())


def asymmetric_cost(src, env_tree, penalty, threshold):
    """Original one-sided version: leak only."""
    d, _ = env_tree.query(src, k=1)
    base = float(d.mean())
    if penalty <= 0:
        return base
    over = d[d > threshold]
    if len(over) == 0:
        return base
    return base + penalty * float((over - threshold).sum()) / len(src)


def icp(src_pts, env_pts, env_tree, max_iter, penalty, threshold,
        mode="author", w_leak=1.0, w_cover=1.0):
    """
    ICP with optional bidirectional cost.

    The weights enter the Procrustes step itself, not only the choice of the
    best iteration: an unweighted Procrustes produces exactly the same motion
    with or without a penalty.
    """
    src = src_pts.copy()
    T_total = np.eye(4)
    best_cost, best_T = float("inf"), np.eye(4)

    for _ in range(max_iter):
        dist, idx = env_tree.query(src, k=1)
        closest = env_pts[idx]

        if mode == "author":
            fill, leak = author_cost(src, env_pts, env_tree)
            cost = fill + penalty * leak
        elif mode == "bidir":
            cost, _, _, _ = bidirectional_cost(src, env_tree, None, env_pts,
                                               penalty, threshold, w_leak, w_cover)
        else:
            cost = asymmetric_cost(src, env_tree, penalty, threshold)
        if cost < best_cost:
            best_cost, best_T = cost, T_total.copy()

        if mode == "author" and penalty > 0:
            # no threshold: the weight grows with the distance itself
            w = 1.0 + penalty * dist / max(float(dist.mean()), 1e-9)
        elif penalty > 0:
            w = 1.0 + penalty * np.maximum(dist - threshold, 0.0) / max(threshold, 1e-9)
        else:
            w = np.ones(len(src))

        srcA, tgtA, wA = src, closest, w

        if mode == "author" and penalty > 0:
            # fill term: EVERY envelope point is a target, weighted by its
            # distance to the nearest atom
            src_tree = cKDTree(src)
            d_env, i_env = src_tree.query(env_pts, k=1)
            srcA = np.vstack([src, src[i_env]])
            tgtA = np.vstack([closest, env_pts])
            w_extra = d_env / max(float(d_env.mean()), 1e-9)
            wA = np.concatenate([w, w_extra])
        elif mode == "bidir" and penalty > 0 and w_cover > 0:
            # uncovered envelope points: extra targets that pull the
            # structure towards the empty regions
            src_tree = cKDTree(src)
            d_env, i_env = src_tree.query(env_pts, k=1)
            unc = d_env > threshold
            if unc.any():
                src_extra = src[i_env[unc]]
                tgt_extra = env_pts[unc]
                w_extra = (w_cover * penalty *
                           (d_env[unc] - threshold) / max(threshold, 1e-9))
                srcA = np.vstack([src, src_extra])
                tgtA = np.vstack([closest, tgt_extra])
                wA = np.concatenate([w, w_extra])

        wA = wA / wA.sum()
        sm = (wA[:, None] * srcA).sum(axis=0)
        tm = (wA[:, None] * tgtA).sum(axis=0)
        H = (srcA - sm).T @ (wA[:, None] * (tgtA - tm))
        U, _, Vt = np.linalg.svd(H)
        R = Vt.T @ U.T
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T
        t = tm - R @ sm

        src = (R @ src.T).T + t
        T = np.eye(4)
        T[:3, :3], T[:3, 3] = R, t
        T_total = T @ T_total

    if mode == "author":
        fill, leak = author_cost(src, env_pts, env_tree)
        cost = fill + penalty * leak
    elif mode == "bidir":
        cost, _, _, _ = bidirectional_cost(src, env_tree, None, env_pts,
                                           penalty, threshold, w_leak, w_cover)
    else:
        cost = asymmetric_cost(src, env_tree, penalty, threshold)
    if cost < best_cost:
        best_cost, best_T = cost, T_total.copy()
    return best_T, best_cost


# evaluation metrics

def chamfer(a, b):
    return float((cKDTree(b).query(a)[0].mean() + cKDTree(a).query(b)[0].mean()) / 2.0)


def hausdorff(a, b):
    return float(max(cKDTree(b).query(a)[0].max(), cKDTree(a).query(b)[0].max()))


def dice(a, b, voxel=1.0):
    """Volumetric overlap by voxelisation on a common grid."""
    allp = np.vstack([a, b])
    origin = allp.min(axis=0)
    ka = set(map(tuple, np.floor((a - origin) / voxel).astype(int)))
    kb = set(map(tuple, np.floor((b - origin) / voxel).astype(int)))
    if not ka or not kb:
        return 0.0
    return 2.0 * len(ka & kb) / (len(ka) + len(kb))


def frac_outside(pts, env):
    """Fraction outside the convex hull of the envelope. A lower bound of the
    real violation, since the hull is larger than a concave envelope."""
    try:
        return float((Delaunay(env).find_simplex(pts) < 0).mean())
    except Exception:
        return float("nan")


# one model / envelope pair

@dataclass
class AlignmentResult:
    transform: np.ndarray        # 4x4, original model frame -> envelope frame
    cost: float
    mirrored: bool               # True if the mirror image fitted better
    fitted: np.ndarray           # aligned (possibly downsampled) C-alpha trace
    metrics: dict = field(default_factory=dict)
    n_model: int = 0
    n_model_full: int = 0
    n_env: int = 0
    grid: float = 0.0


def align(model_ca, env_full, mode="author", penalty=0.2, restarts=3, max_iter=50,
          max_points=3000, sample_env=5000, seed=42, enantiomorphs=True,
          w_leak=1.0, w_cover=1.0):
    """
    Aligns a C-alpha trace to an envelope point cloud.

    If the result is mirrored, `transform` contains a reflection (det = -1):
    SAXS envelopes do not define chirality, so the envelope, not the protein,
    has the wrong hand. See proper_transform().
    """
    mod = np.asarray(model_ca, dtype=float)
    env_full = np.asarray(env_full, dtype=float)
    if len(env_full) < 10:
        raise ValueError(f"envelope has {len(env_full)} points (need at least 10)")
    if len(mod) < 4:
        raise ValueError(f"model has {len(mod)} atoms (need at least 4)")
    if restarts < 1:
        raise ValueError("restarts must be >= 1")

    rng = np.random.default_rng(seed)
    best_mirror = False

    # model downsampling (max_points; use >= n to disable)
    n_model_full = len(mod)
    if len(mod) > max_points:
        mod = mod[np.sort(rng.choice(len(mod), max_points, replace=False))]

    # envelope sampling
    env = env_full
    if sample_env and len(env) > sample_env:
        env = env[np.sort(rng.choice(len(env), sample_env, replace=False))]

    thr = max(3.0, grid_spacing(env))
    env_tree = cKDTree(env)
    center = mod.mean(axis=0)

    # SAXS envelopes do not define chirality: ab initio reconstruction gives a
    # shape and its mirror image with equal probability. Both hands are tested.
    mirrors = [np.eye(3)]
    if enantiomorphs:
        mirrors.append(MIRROR_Z.copy())        # reflection in the xy plane

    best = None
    for M in mirrors:
        if np.allclose(M, np.eye(3)):
            base, T_mir = mod, np.eye(4)
        else:
            base = (M @ (mod - center).T).T + center
            T_mir = np.eye(4)
            T_mir[:3, :3] = M
            T_mir[:3, 3] = center - M @ center

        for k in range(restarts):
            if k == 0:
                start, T_init = base, T_mir
            else:
                Rm = random_rotation(int(rng.integers(1 << 31)))
                start = (Rm @ (base - center).T).T + center
                T_rot = np.eye(4)
                T_rot[:3, :3] = Rm
                T_rot[:3, 3] = center - Rm @ center
                T_init = T_rot @ T_mir

            T_icp, cost = icp(start, env, env_tree, max_iter, penalty, thr,
                              mode=mode, w_leak=w_leak, w_cover=w_cover)
            if best is None or cost < best[0]:
                best = (cost, T_icp @ T_init)
                best_mirror = not np.allclose(M, np.eye(3))

    cost, T = best
    fit = apply_T(mod, T)

    metrics = dict(
        chamfer=chamfer(fit, env),
        hausdorff=hausdorff(fit, env),
        dice=dice(fit, env),
        frac_outside=frac_outside(fit, env_full),
        coverage=float((cKDTree(fit).query(env_full)[0] <= thr).mean()),
        centroid_sep=float(np.linalg.norm(fit.mean(0) - env_full.mean(0))),
        median_dist=float(np.median(cKDTree(env_full).query(fit)[0])),
    )
    return AlignmentResult(transform=T, cost=cost, mirrored=best_mirror, fitted=fit,
                           metrics=metrics, n_model=len(mod),
                           n_model_full=n_model_full, n_env=len(env), grid=thr)


def proper_transform(T):
    """
    Splits a transform that may contain a reflection into a proper rigid
    motion for the model plus a reflection for the envelope.

    If det(R) = -1, then  M @ (R p + t) = (M R) p + M t  with M R a rotation:
    the model moved by (M R, M t) fits the envelope reflected by M exactly as
    well as the mirrored model fits the original envelope (M is an isometry),
    without turning the protein into its D-amino-acid mirror image.

    Returns (T_model, M_envelope), M_envelope = None if no reflection.
    """
    if np.linalg.det(T[:3, :3]) > 0:
        return T, None
    M = np.eye(4)
    M[:3, :3] = MIRROR_Z
    return M @ T, M
