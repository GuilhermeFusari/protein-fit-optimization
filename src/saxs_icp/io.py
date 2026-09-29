"""
Reading and writing of atomic models and SAXS envelopes.

Format is detected from the file CONTENT, not the extension: SASBDB envelopes
carry a .cif extension but some are PDB-formatted, and reading one format with
the other parser returns zero coordinates silently.
"""
import re

import numpy as np

_ENC = dict(encoding="utf-8", errors="replace")

# mmCIF tokens: quoted strings (a quote only closes when followed by
# whitespace or end of line) or bare words
_CIF_TOKEN = re.compile(r"""'(?:[^']|'(?=\S))*'|"(?:[^"]|"(?=\S))*"|\S+""")


def is_mmcif(path):
    """True if the file declares an mmCIF _atom_site category."""
    with open(path, "r", **_ENC) as f:
        return any(line.startswith("_atom_site.") for line in f)


def _read_pdb_format(path, ca_only):
    pts = []
    with open(path, "r", **_ENC) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if ca_only and line[12:16].strip() != "CA":
                continue
            try:
                pts.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                continue
    return pts


def _atom_site_loop(path):
    """
    Returns (cols, rows) of the mmCIF _atom_site loop, reading columns by the
    positions declared in the header instead of assuming a fixed order.
    """
    cols, rows, in_loop = {}, [], False
    with open(path, "r", **_ENC) as f:
        for line in f:
            s = line.strip()
            if s.startswith("_atom_site."):
                cols[s.split(".", 1)[1].split()[0]] = len(cols)
                in_loop = True
                continue
            if in_loop:
                if s.startswith("#") or s.startswith("loop_") or not s:
                    if rows:
                        break
                    in_loop = bool(cols)
                    continue
                if s.startswith("_"):
                    break
                rows.append(_CIF_TOKEN.findall(s))
    return cols, rows


def _read_cif_format(path, ca_only):
    cols, rows = _atom_site_loop(path)

    need = ("Cartn_x", "Cartn_y", "Cartn_z")
    if not all(k in cols for k in need):
        return []
    ix, iy, iz = (cols[k] for k in need)
    ia = cols.get("label_atom_id", cols.get("auth_atom_id"))

    pts = []
    for parts in rows:
        # a row must have exactly one value per declared column; anything
        # else (e.g. PDB-formatted lines under an mmCIF header) would shift
        # the columns and give wrong coordinates
        if len(parts) != len(cols):
            continue
        if ca_only and ia is not None and len(parts) > ia:
            if parts[ia].strip('"').strip("'") != "CA":
                continue
        try:
            pts.append((float(parts[ix]), float(parts[iy]), float(parts[iz])))
        except ValueError:
            continue
    return pts


def read_coords(path, ca_only=False):
    """
    Reads coordinates from a PDB or mmCIF file as an (N, 3) float array.

    With ca_only=True returns the C-alpha atoms; if none can be identified
    (e.g. a coarse-grained model) returns all atoms.
    """
    if is_mmcif(path):
        # hybrid files (mmCIF header, PDB-formatted atom lines) fall back to
        # the fixed-column reader
        readers = (_read_cif_format, _read_pdb_format)
    else:
        readers = (_read_pdb_format,)
    for reader in readers:
        pts = reader(path, ca_only)
        if not pts and ca_only:
            pts = reader(path, False)
        if pts:
            break
    return np.asarray(pts, dtype=float).reshape(-1, 3)


def count_models(path):
    """Number of MODEL records in a PDB file (0 if none)."""
    with open(path, "r", **_ENC) as f:
        return sum(1 for line in f if line.startswith("MODEL "))


def apply_T(pts, T):
    return (T[:3, :3] @ pts.T).T + T[:3, 3]


def write_points_pdb(pts, path, atom="CA", res="ALA"):
    """Writes a point set as a PDB of single atoms (C-alpha trace, dummy beads)."""
    with open(path, "w", encoding="utf-8") as fh:
        for i, (x, y, z) in enumerate(pts, 1):
            fh.write(f"ATOM  {i % 100000:5d}  {atom:<3s} {res:3s} A{i % 10000:4d}    "
                     f"{_fmt_coord(x)}{_fmt_coord(y)}{_fmt_coord(z)}"
                     f"  1.00  0.00           C\n")
        fh.write("END\n")


def _fmt_coord(v):
    s = f"{v:8.3f}"
    if len(s) != 8:
        raise ValueError(f"coordinate {v:.3f} does not fit the PDB format")
    return s


def write_transformed_pdb(src, dst, T):
    """
    Applies the rigid transform T to every ATOM/HETATM record of a PDB file,
    keeping all other records. ANISOU records are dropped, since the
    anisotropic factors would no longer match the rotated frame.
    """
    R, t = T[:3, :3], T[:3, 3]
    out = []
    with open(src, "r", **_ENC) as f:
        for line in f:
            if line.startswith("ANISOU"):
                continue
            if line.startswith(("ATOM", "HETATM")):
                try:
                    p = np.array([float(line[30:38]), float(line[38:46]),
                                  float(line[46:54])])
                except ValueError:
                    out.append(line)
                    continue
                x, y, z = R @ p + t
                body = line.rstrip("\r\n")
                line = (body[:30] + _fmt_coord(x) + _fmt_coord(y) + _fmt_coord(z)
                        + body[54:] + "\n")
            out.append(line)
    with open(dst, "w", encoding="utf-8") as fh:
        fh.writelines(out)


def write_transformed_cif(src, dst, T):
    """
    Applies the rigid transform T to the _atom_site coordinates of an mmCIF
    file, keeping everything else. Raises ValueError if a row cannot be
    tokenised unambiguously, instead of writing untransformed coordinates.
    """
    R, t = T[:3, :3], T[:3, 3]
    with open(src, "r", **_ENC) as f:
        lines = f.readlines()

    cols, in_loop, out, n_done = {}, False, [], 0
    for line in lines:
        s = line.strip()
        if s.startswith("_atom_site."):
            cols[s.split(".", 1)[1].split()[0]] = len(cols)
            in_loop = True
            out.append(line)
            continue
        if in_loop and cols and s and not s.startswith(("#", "loop_", "_")):
            if not all(k in cols for k in ("Cartn_x", "Cartn_y", "Cartn_z")):
                raise ValueError("mmCIF _atom_site has no Cartn_x/y/z columns")
            tok = _CIF_TOKEN.findall(s)
            if len(tok) != len(cols):
                raise ValueError(f"cannot tokenise mmCIF row: {s[:60]}")
            ix, iy, iz = cols["Cartn_x"], cols["Cartn_y"], cols["Cartn_z"]
            p = np.array([float(tok[ix]), float(tok[iy]), float(tok[iz])])
            x, y, z = R @ p + t
            tok[ix], tok[iy], tok[iz] = f"{x:.3f}", f"{y:.3f}", f"{z:.3f}"
            out.append(" ".join(tok) + "\n")
            n_done += 1
            continue
        if in_loop and (s.startswith(("#", "loop_")) or (s.startswith("_") and n_done)):
            in_loop = False
            if n_done:
                cols = {}
        out.append(line)

    if not n_done:
        raise ValueError("no _atom_site coordinates found to transform")
    with open(dst, "w", encoding="utf-8") as fh:
        fh.writelines(out)
