"""Render the Pareto dominance and epsilon-dominance figure used in Section 2.7 of the guide.

The figure is drawn in TikZ, compiled with pdflatex and rasterised with pdftoppm
(both ship with MiKTeX and TeX Live).

    python tools/dominance-figure/make_dominance_figure.py

Output: assets/images/dominance-and-epsilon-dominance.png

The epsilon archive shown in panel (d) is computed here with the same rules the
guide describes (box dominance, then the closest point to the box's preferred
corner), so the panel always agrees with the points drawn.
"""

import math
import random
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
BUILD = HERE / "build"
OUT = HERE.parents[1] / "assets" / "images"
NAME = "dominance-and-epsilon-dominance"

DPI = 300
EPS = 0.8                  # box size in both objectives
SIZE = 6 * EPS             # extent of each plot
PARETO = [(0.40, 3.85), (0.95, 2.70), (1.40, 1.95), (1.72, 1.52),
          (2.15, 1.12), (2.85, 0.90), (3.85, 0.62)]
Z = (2.25, 2.25)           # reference solution for the relation panels
RED = "red!80!black"
DOMINATED = "black!12"     # fill for dominated regions
MUTED = "black!45"         # dominated solutions

PREAMBLE = r"""\documentclass[tikz,border=4pt]{standalone}
\usetikzlibrary{patterns.meta, shapes.misc}
\tikzset{
  every picture/.style={line width=0.5pt, line join=round, >=stealth},
  cross/.style={cross out, draw, minimum size=5pt, inner sep=0pt, outer sep=0pt},
  hatch/.style={pattern={Lines[angle=45, distance=4pt, line width=0.35pt]}, pattern color=black!40},
  plabel/.style={font=\small, inner sep=1pt, fill=white},
}
\begin{document}
"""


def pt(p):
    return f"({p[0]:.3f},{p[1]:.3f})"


def polyline(points):
    return " -- ".join(pt(p) for p in points)


def smooth(points, n=16):
    """Catmull-Rom spline through the points, returned as a dense polyline."""
    padded = [points[0], *points, points[-1]]
    out = []
    for i in range(1, len(padded) - 2):
        p0, p1, p2, p3 = padded[i - 1 : i + 3]
        for k in range(n):
            t = k / n
            out.append(tuple(
                0.5 * (2 * p1[j] + (p2[j] - p0[j]) * t
                       + (2 * p0[j] - 5 * p1[j] + 4 * p2[j] - p3[j]) * t * t
                       + (3 * p1[j] - p0[j] - 3 * p2[j] + p3[j]) * t ** 3)
                for j in range(2)))
    out.append(points[-1])
    return out


FRONT = smooth([(0.12, SIZE)] + PARETO + [(SIZE, 0.45)])


def front_y(x):
    for (x0, y0), (x1, y1) in zip(FRONT, FRONT[1:]):
        if x0 <= x <= x1:
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
    return FRONT[-1][1]


# --------------------------------------------------------------------------- epsilon archive

def box(p):
    return (math.floor(p[0] / EPS), math.floor(p[1] / EPS))


def box_dominates(a, b):
    return a != b and a[0] <= b[0] and a[1] <= b[1]


def corner_distance(p):
    b = box(p)
    return math.dist(p, (b[0] * EPS, b[1] * EPS))


def epsilon_archive(points):
    """Indices kept by the epsilon archive: one per box, never a box-dominated box."""
    kept = [i for i, p in enumerate(points)
            if not any(box_dominates(box(q), box(p)) for q in points)]
    best = {}
    for i in kept:
        b = box(points[i])
        if b not in best or corner_distance(points[i]) < corner_distance(points[best[b]]):
            best[b] = i
    return sorted(best.values())


RETAINED = epsilon_archive(PARETO)
assert RETAINED == [0, 2, 4, 6], RETAINED
RETAINED_BOXES = [box(PARETO[i]) for i in RETAINED]


def dominated_candidates(n=22, seed=3):
    """Feasible, Pareto-dominated solutions whose boxes are box-dominated by the archive."""
    rng = random.Random(seed)
    points = []
    while len(points) < n:
        p = (rng.uniform(0.3, SIZE - 0.25), rng.uniform(0.3, SIZE - 0.25))
        near_grid = min(abs(v / EPS - round(v / EPS)) for v in p) < 0.12
        if (p[1] < front_y(p[0]) + 0.25 or near_grid
                or not any(q[0] + 0.12 <= p[0] and q[1] + 0.12 <= p[1] for q in PARETO)
                or not any(box_dominates(b, box(p)) for b in RETAINED_BOXES)
                or min(math.dist(p, q) for q in points + PARETO) < 0.42):
            continue
        points.append(p)
    return points


CANDIDATES = dominated_candidates()


# --------------------------------------------------------------------------- drawing

def axes():
    return [rf"\draw[->] (0,0) -- ({SIZE + 0.45},0) node[below left, inner sep=2pt] {{$f_1$}};",
            rf"\draw[->] (0,0) -- (0,{SIZE + 0.45}) node[below left, inner sep=2pt] {{$f_2$}};"]


def epsilon_grid():
    lines = []
    for k in range(1, 7):
        v = k * EPS
        lines.append(rf"\draw[dashed, black!50, line width=0.4pt] ({v:.2f},0) -- ({v:.2f},{SIZE});")
        lines.append(rf"\draw[dashed, black!50, line width=0.4pt] (0,{v:.2f}) -- ({SIZE},{v:.2f});")
    lines.append(rf"\draw[|<->|, line width=0.4pt] (0,-0.3) -- ({EPS},-0.3) node[midway, below, font=\small] {{$\epsilon_1$}};")
    lines.append(rf"\draw[|<->|, line width=0.4pt] (-0.3,0) -- (-0.3,{EPS}) node[midway, left, font=\small] {{$\epsilon_2$}};")
    return lines


def infeasible_region():
    return [rf"\fill[hatch] (0,0) -- (0,{SIZE}) -- {polyline(FRONT)} -- ({SIZE},0) -- cycle;",
            rf"\draw[black!60] {polyline(FRONT)};",
            rf"\node[plabel, font=\footnotesize, text=black!70] at (0.85,0.55) {{Infeasible}};"]


def cross(p, colour, label=None, thick=False, place=("north east", -0.05, -0.05, "white")):
    width = "0.9pt" if thick else "0.6pt"
    lines = [rf"\node[cross, draw={colour}, line width={width}] at {pt(p)} {{}};"]
    if label:
        anchor, dx, dy, fill = place
        lines.append(rf"\node[plabel, fill={fill}, text={colour}, anchor={anchor}] at {pt((p[0] + dx, p[1] + dy))} {{{label}}};")
    return lines


# Label positions in panel (d), chosen to keep clear of the archive-box outlines
EPSILON_LABELS = {0: ("north", 0.0, -0.07, "white"), 3: ("north east", -0.16, -0.05, "white"),
                  4: ("west", 0.31, 0.06, "none")}


def candidates():
    lines = []
    for p in CANDIDATES:
        lines += cross(p, MUTED)
    return lines


def caption(text, y=-0.95):
    return [rf"\node[anchor=north] at ({SIZE / 2},{y}) {{{text}}};"]


def panel_dominance():
    zx, zy = Z
    lines = [rf"\fill[{DOMINATED}] {pt(Z)} rectangle ({SIZE},{SIZE});",
             rf"\draw[dashed, black!55, line width=0.4pt] ({zx},0) -- ({zx},{SIZE});",
             rf"\draw[dashed, black!55, line width=0.4pt] (0,{zy}) -- ({SIZE},{zy});"]
    lines += axes()
    text = r"\node[font=\small, align=center, text=black!75] at "
    lines += [text + rf"({(zx + SIZE) / 2:.2f},{(zy + SIZE) / 2:.2f}) {{dominated\\by $\mathbf{{z}}$}};",
              text + rf"({zx / 2:.2f},{zy / 2:.2f}) {{dominates\\$\mathbf{{z}}$}};",
              text + rf"({zx / 2:.2f},{(zy + SIZE) / 2:.2f}) {{incomparable}};",
              text + rf"({(zx + SIZE) / 2:.2f},{zy / 2:.2f}) {{incomparable}};"]
    lines += cross(Z, "black", thick=True)
    lines.append(rf"\node[anchor=north east, inner sep=1.5pt] at {pt(Z)} {{$\mathbf{{z}}$}};")
    return lines + caption("(a) Pareto dominance", -0.55)


def panel_pareto_set():
    stairs = [(PARETO[0][0], SIZE)]
    for p, q in zip(PARETO, PARETO[1:]):
        stairs += [p, (q[0], p[1])]
    stairs += [PARETO[-1], (SIZE, PARETO[-1][1])]
    lines = [rf"\fill[{DOMINATED}] {polyline(stairs)} -- ({SIZE},{SIZE}) -- cycle;"]
    lines += infeasible_region()
    lines.append(rf"\draw[{RED}, line width=0.8pt] {polyline(stairs)};")
    lines += axes() + candidates()
    for i, p in enumerate(PARETO):
        lines += cross(p, RED, rf"$p_{i + 1}$", thick=True)
    return lines + caption("(b) Pareto set", -0.55)


def panel_box_dominance():
    b = box(Z)
    x0, y0 = b[0] * EPS, b[1] * EPS
    lines = [rf"\fill[{DOMINATED}] ({x0},{y0}) rectangle ({SIZE},{SIZE});",
             rf"\fill[white] ({x0},{y0}) rectangle ({x0 + EPS},{y0 + EPS});"]
    lines += epsilon_grid() + axes()
    lines.append(rf"\draw[line width=0.8pt] ({x0},{y0}) rectangle ({x0 + EPS},{y0 + EPS});")
    lines.append(rf"\node[font=\small, align=center, text=black!75] at ({(x0 + EPS + SIZE) / 2:.2f},{(y0 + EPS + SIZE) / 2:.2f}) "
                 r"{box-dominated\\by $\mathbf{z}$};")
    lines += cross(Z, "black", thick=True)
    lines.append(rf"\node[anchor=north east, inner sep=1.5pt] at {pt(Z)} {{$\mathbf{{z}}$}};")
    return lines + caption(r"(c) $\epsilon$-box dominance")


def panel_epsilon_set():
    feasible = FRONT + [(SIZE, SIZE)]
    lines = [r"\begin{scope}", rf"\clip {polyline(feasible)} -- cycle;"]
    for i in range(6):
        for j in range(6):
            if any(box_dominates(b, (i, j)) for b in RETAINED_BOXES):
                lines.append(rf"\fill[{DOMINATED}] ({i * EPS:.2f},{j * EPS:.2f}) rectangle ++({EPS},{EPS});")
    lines.append(r"\end{scope}")
    lines += infeasible_region() + epsilon_grid() + axes()
    # tie-break between solutions sharing a box: the one nearer the preferred corner is kept
    shared = [i for i, p in enumerate(PARETO)
              if i not in RETAINED and box(p) in RETAINED_BOXES]
    for i in shared:
        b = box(PARETO[i])
        corner = (b[0] * EPS, b[1] * EPS)
        rival = PARETO[RETAINED[RETAINED_BOXES.index(b)]]
        for p in (PARETO[i], rival):
            lines.append(rf"\draw[densely dotted, line width=0.7pt] {pt(corner)} -- {pt(p)};")
        lines.append(rf"\fill {pt(corner)} circle (1.6pt);")
    lines += candidates()
    for i, p in enumerate(PARETO):
        colour = RED if i in RETAINED else MUTED
        place = EPSILON_LABELS.get(i, ("north east", -0.05, -0.05, "white"))
        lines += cross(p, colour, rf"$p_{i + 1}$", thick=i in RETAINED, place=place)
    for bx, by in RETAINED_BOXES:  # drawn last so no label covers an outline
        lines.append(rf"\draw[{RED}, line width=0.8pt] ({bx * EPS:.2f},{by * EPS:.2f}) rectangle ++({EPS},{EPS});")
    return lines + caption(r"(d) $\epsilon$-Pareto set")


def legend(y):
    swatch = r"\tikz[baseline=-0.6ex]{\fill[%s] (0,-0.12) rectangle (0.34,0.12); \draw[black!40, line width=0.3pt] (0,-0.12) rectangle (0.34,0.12);}"
    items = [
        (rf"\tikz[baseline=-0.6ex]\node[cross, draw={RED}, line width=0.9pt] {{}};", "Retained solution"),
        (rf"\tikz[baseline=-0.6ex]\node[cross, draw={MUTED}, line width=0.6pt] {{}};", "Dominated solution"),
        (swatch % DOMINATED, "Dominated region"),
        (swatch % "hatch", "Infeasible region"),
    ]
    row = r"\hspace{1.4em}".join(rf"{symbol}\hspace{{0.5em}}{text}" for symbol, text in items)
    return [rf"\node[font=\small] at ({(6.7 + SIZE) / 2:.2f},{y}) {{{row}}};"]

def figure():
    panels = [((0, 6.75), panel_dominance()), ((6.7, 6.75), panel_pareto_set()),
              ((0, 0), panel_box_dominance()), ((6.7, 0), panel_epsilon_set())]
    lines = [r"\begin{tikzpicture}[x=1cm, y=1cm]"]
    for (x, y), body in panels:
        lines.append(rf"\begin{{scope}}[shift={{({x},{y})}}]")
        lines += body
        lines.append(r"\end{scope}")
    lines += legend(-2.0)
    lines.append(r"\end{tikzpicture}")
    return "\n".join(lines)


if __name__ == "__main__":
    BUILD.mkdir(exist_ok=True)
    tex = BUILD / f"{NAME}.tex"
    tex.write_text(PREAMBLE + figure() + "\n\\end{document}\n", encoding="utf-8")
    subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
                    f"-output-directory={BUILD}", str(tex)], check=True, stdout=subprocess.DEVNULL)
    OUT.mkdir(parents=True, exist_ok=True)
    subprocess.run(["pdftoppm", "-r", str(DPI), "-png", "-singlefile",
                    str(BUILD / f"{NAME}.pdf"), str(OUT / NAME)], check=True)
    print(f"Wrote {OUT / NAME}.png")
