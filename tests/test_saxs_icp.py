import json
import os
import shutil

import numpy as np
import pytest

from saxs_icp import core
from saxs_icp.cli import main
from saxs_icp.io import (is_mmcif, read_coords, write_transformed_cif,
                         write_transformed_pdb)

EX = os.path.join(os.path.dirname(__file__), os.pardir, "examples")

# Rows of benchmark/final_lam02.csv (lambda = 0.2, author mode, 3000 points).
# The seed is 42 + the entry's index in the 50-entry benchmark.
PAPER = {
    "SASDTS8": dict(seed=44, cost=4.1854, chamfer=3.0231, hausdorff=17.3935,
                    dice=0.0052, frac_outside=0.2336, coverage=0.4274,
                    centroid_sep=1.041, median_dist=1.491, mirrored=False),
    "SASDXB2": dict(seed=83, cost=5.8884, chamfer=4.1722, hausdorff=24.5102,
                    dice=0.00168, frac_outside=0.0144, coverage=0.4422,
                    centroid_sep=3.106, median_dist=2.969, mirrored=True),
}
DIGITS = dict(cost=4, chamfer=4, hausdorff=4, dice=5, frac_outside=4,
              coverage=4, centroid_sep=3, median_dist=3)


def ex(acc, kind):
    return os.path.join(EX, f"{acc}_{kind}")


@pytest.mark.parametrize("acc", sorted(PAPER))
def test_reproduces_published_numbers(acc):
    exp = PAPER[acc]
    r = core.align(read_coords(ex(acc, "model.pdb"), ca_only=True),
                   read_coords(ex(acc, "envelope.cif")), seed=exp["seed"])
    got = dict(r.metrics, cost=r.cost)
    for k, nd in DIGITS.items():
        assert round(got[k], nd) == exp[k], k
    assert r.mirrored == exp["mirrored"]


def test_random_rotation_matches_scipy():
    from scipy.spatial.transform import Rotation
    for seed in (0, 1, 12345, 2**31 - 1):
        try:
            ref = Rotation.random(random_state=seed).as_matrix()
        except TypeError:
            pytest.skip("SciPy without the random_state argument")
        assert np.array_equal(core.random_rotation(seed), ref)


def test_format_detection():
    assert not is_mmcif(ex("SASDTS8", "envelope.cif"))    # PDB content, .cif name
    assert is_mmcif(ex("SASDXB2", "envelope.cif"))
    assert is_mmcif(ex("SASDXB2", "model.pdb"))           # mmCIF content, .pdb name
    assert len(read_coords(ex("SASDTS8", "envelope.cif"))) == 2033
    assert len(read_coords(ex("SASDXB2", "envelope.cif"))) == 778
    assert len(read_coords(ex("SASDXB2", "model.pdb"), ca_only=True)) == 416


def _rigid_T(seed):
    T = np.eye(4)
    T[:3, :3] = core.random_rotation(seed)
    T[:3, 3] = [10.0, -5.0, 3.0]
    return T


@pytest.mark.parametrize("acc,writer", [("SASDTS8", write_transformed_pdb),
                                        ("SASDXB2", write_transformed_cif)])
def test_rewrite_applies_transform(tmp_path, acc, writer):
    src = ex(acc, "model.pdb")
    T = _rigid_T(3)
    dst = str(tmp_path / "out")
    writer(src, dst, T)
    a, b = read_coords(src), read_coords(dst)
    assert a.shape == b.shape
    assert np.abs(core.apply_T(a, T) - b).max() < 1e-3


def test_rewrite_cif_quoted_tokens(tmp_path):
    cif = tmp_path / "q.cif"
    cif.write_text(
        "data_q\nloop_\n_atom_site.group_PDB\n_atom_site.id\n"
        "_atom_site.label_atom_id\n_atom_site.label_comp_id\n"
        "_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"
        "ATOM 1 \"O5'\" DA 1.000 2.000 3.000\n"
        "ATOM 2 CA 'A B' 4.000 5.000 6.000\n#\n")
    T = _rigid_T(5)
    write_transformed_cif(str(cif), str(tmp_path / "o.cif"), T)
    out = (tmp_path / "o.cif").read_text()
    assert "\"O5'\"" in out and "'A B'" in out
    got = read_coords(str(tmp_path / "o.cif"))
    assert np.abs(core.apply_T(np.array([[1., 2, 3], [4, 5, 6]]), T) - got).max() < 1e-3


def test_proper_transform_keeps_chirality():
    T = np.eye(4)
    T[:3, :3] = core.random_rotation(1) @ core.MIRROR_Z
    T[:3, 3] = [1.0, 2.0, 3.0]
    Tm, M = core.proper_transform(T)
    assert np.linalg.det(Tm[:3, :3]) > 0
    p = np.random.default_rng(0).normal(size=(20, 3))
    assert np.allclose(core.apply_T(p, Tm), core.apply_T(core.apply_T(p, T), M))
    assert core.proper_transform(np.eye(4))[1] is None


@pytest.mark.parametrize("acc", sorted(PAPER))
def test_cli_align(tmp_path, acc):
    exp = PAPER[acc]
    out = tmp_path / "out"
    rc = main(["align", ex(acc, "model.pdb"), ex(acc, "envelope.cif"),
               "-o", str(out), "--seed", str(exp["seed"])])
    assert rc == 0
    rep = json.loads((out / "report.json").read_text())
    assert round(rep["metrics"]["chamfer_A"], 4) == exp["chamfer"]
    assert rep["envelope_mirrored"] == exp["mirrored"]

    # the written model fits the (possibly mirrored) envelope as reported,
    # and is never itself mirrored
    env_file = rep["outputs"]["envelope_mirrored"] or ex(acc, "envelope.cif")
    fitted = read_coords(rep["outputs"]["model"], ca_only=True)
    assert abs(core.chamfer(fitted, read_coords(env_file))
               - core.chamfer(core.apply_T(read_coords(ex(acc, "model.pdb"), True),
                                           np.array(rep["transform_model"])),
                              read_coords(env_file))) < 1e-2
    assert np.linalg.det(np.array(rep["transform_model"])[:3, :3]) > 0


@pytest.mark.parametrize("args", [
    ["align", "missing.pdb", "missing.cif"],
    ["align", "EX_MODEL", "EX_ENV", "--penalty", "-1"],
    ["align", "EX_MODEL", "EX_ENV", "--restarts", "0"],
])
def test_cli_rejects_bad_arguments(args):
    args = [a.replace("EX_MODEL", ex("SASDTS8", "model.pdb"))
             .replace("EX_ENV", ex("SASDTS8", "envelope.cif")) for a in args]
    with pytest.raises(SystemExit) as e:
        main(args)
    assert e.value.code == 2


def test_cli_reports_unreadable_input(tmp_path, capsys):
    bad = tmp_path / "bad.pdb"
    bad.write_text("not a structure\n")
    rc = main(["align", str(bad), ex("SASDTS8", "envelope.cif"), "-o", str(tmp_path / "o")])
    assert rc == 1
    assert "error" in capsys.readouterr().err


def test_cli_pack(tmp_path):
    inp = tmp_path / "in"
    inp.mkdir()
    for i in range(3):
        shutil.copy(ex("SASDTS8", "model.pdb"), inp / f"conf{i}.pdb")
    (inp / "broken.pdb").write_text("garbage\n")
    out = tmp_path / "out"
    rc = main(["pack", str(inp), ex("SASDTS8", "envelope.cif"), "-o", str(out),
               "--copies", "2", "--workers", "1", "--seed", "1"])
    assert rc == 0
    assert (out / "ALL_TOP_ALIGNED.pdb").exists()
    assert "Number of Copies:         2" in (out / "report.txt").read_text()


def test_cli_benchmark(tmp_path):
    sas = tmp_path / "data" / "sasbdb"
    sas.mkdir(parents=True)
    for acc in PAPER:
        shutil.copy(ex(acc, "model.pdb"), sas / f"{acc}_model.pdb")
        shutil.copy(ex(acc, "envelope.cif"), sas / f"{acc}_envelope.cif")
    csv = tmp_path / "r.csv"
    rc = main(["benchmark", "--base", str(tmp_path), "--out", str(csv), "--workers", "1"])
    assert rc == 0
    lines = csv.read_text().splitlines()
    assert len(lines) == 3 and lines[0].startswith("entry,status,custo,chamfer")


def test_hybrid_mmcif_header_with_pdb_lines(tmp_path):
    # mmCIF header declaring 16 columns followed by PDB fixed-column lines
    # (found in SASBDB fit files): must fall back to the PDB reader
    head = "data_x\nloop_\n" + "".join(
        f"_atom_site.{c}\n" for c in (
            "group_PDB id type_symbol label_atom_id label_alt_id label_comp_id "
            "label_asym_id label_seq_id pdbx_PDB_ins_code Cartn_x Cartn_y Cartn_z "
            "occupancy B_iso_or_equiv pdbx_formal_charge pdbx_PDB_model_num").split())
    lines = [f"ATOM  {i:5d}  CA  ASP A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n"
             for i, (x, y, z) in enumerate([(-4.334, -2.157, 30.781),
                                            (-4.334, 3.843, 30.781),
                                            (16.666, 12.843, 68.965),
                                            (1.0, 2.0, 3.0)], 1)]
    f = tmp_path / "hybrid.pdb"
    f.write_text(head + "".join(lines))
    got = read_coords(str(f))
    assert got.shape == (4, 3)
    assert np.allclose(got[0], [-4.334, -2.157, 30.781])
