"""Render the animated genotype-to-phenotype figures used in Section 2.2 of the guide.

Every frame is drawn as a TikZ picture so the animations match the LaTeX figures
used in the journal articles. The frames are compiled with pdflatex, rasterised
with pdftoppm (both ship with MiKTeX and TeX Live) and assembled into looping
GIFs with Pillow.

    python tools/genotype-animations/make_animations.py

Outputs (written to assets/images/):
    genotype-zonal.gif / .png   zonal include-or-exclude encoding
    genotype-basin.gif / .png   local detention-basin encoding (Gray-coded)
The PNGs are the opening frame, used as a still for readers who prefer reduced motion.
"""

import math
import random
import subprocess
from pathlib import Path

from PIL import Image

HERE = Path(__file__).resolve().parent
BUILD = HERE / "build"
OUT = HERE.parents[1] / "assets" / "images"

DPI = 250                # renders at roughly twice the displayed size
GREY_LEVELS = 64
SHADE = 25               # fill for set bits and included zones, as black!SHADE
HIGHLIGHT = 12           # fill for the active rows of the mapping table
TRANSITION_FRAMES = 8
TRANSITION_MS = 40
HOLD_MS = 1100
GRAY_CODE = ["000", "001", "011", "010", "110", "111", "101", "100"]

PREAMBLE = r"""\documentclass[tikz,border=0pt]{standalone}
\usetikzlibrary{decorations.pathreplacing}
\tikzset{every picture/.style={line width=0.5pt, line join=round},
         strut/.style={text height=1.6ex, text depth=0.5ex}}
\begin{document}
"""


def ease(t):
    return t * t * (3 - 2 * t)


def mix(a, b, t):
    return a + (b - a) * t


def pt(p):
    return f"({p[0]:.3f},{p[1]:.3f})"


def polyline(points):
    return " -- ".join(pt(p) for p in points)


def smooth(points, n=12):
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


def closed_flip_sequence(bits, rng):
    """Flip every chosen bit twice, in random order, so the loop ends where it began."""
    while True:
        seq = bits * 2
        rng.shuffle(seq)
        if all(seq[i] not in seq[i + 1 : i + 3] for i in range(len(seq))):
            return seq


def timeline(genotype, flips):
    """Yield (old, new, flipped bit, t, duration) for every frame of the loop."""
    yield genotype, genotype, None, 1.0, HOLD_MS
    for step, bit in enumerate(flips):
        new = genotype.copy()
        new[bit] ^= 1
        for f in range(1, TRANSITION_FRAMES + 1):
            last = f == TRANSITION_FRAMES
            # The final hold is supplied by the opening frame when the GIF loops.
            hold = HOLD_MS if step < len(flips) - 1 else TRANSITION_MS
            yield genotype, new, bit, ease(f / TRANSITION_FRAMES), hold if last else TRANSITION_MS
        genotype = new


def bit_cells(x0, y0, w, h, old, new, active, t, index_labels=False):
    """Genotype cells: shaded when set, with the flipped bit outlined.

    Index labels are placed above the cells in the zonal figure's y-down coordinates.
    """
    lines = []
    for i, (a, b) in enumerate(zip(old, new)):
        level = mix(a, b, t) * SHADE
        x = x0 + i * w
        box = f"{pt((x, y0))} rectangle {pt((x + w, y0 + h))}"
        lines.append(rf"\filldraw[fill=black!{level:.1f}] {box};")
        digit = b if t >= 0.5 else a
        lines.append(rf"\node at {pt((x + w / 2, y0 + h / 2))} {{{digit}}};")
        if index_labels:
            lines.append(rf"\node[font=\tiny, text=black!55] at {pt((x + w / 2, y0 - 7))} {{{i + 1}}};")
    if active is not None:
        x = x0 + active * w
        lines.append(rf"\draw[line width=1.3pt] {pt((x, y0))} rectangle {pt((x + w, y0 + h))};")
    return lines


# --------------------------------------------------------------------------- zonal

def zonal_geometry():
    """Zone outlines traced from the journal figure, in its pixel coordinates."""
    A, B, C, D = (135, 192), (205, 192), (135, 237), (205, 182)
    E, F, G, H = (275, 182), (275, 122), (205, 245), (275, 252)
    I, J, K, L = (372, 122), (372, 204), (275, 204), (475, 122)
    M, N, O, P = (475, 204), (372, 282), (300, 282), (475, 282)
    Q, TIP, Z1, Z2, T4 = (475, 305), (32, 213), (135, 160), (140, 143), (205, 126)

    z1_top = smooth([TIP, (68, 199), (117, 192), (128, 181), Z1])
    z1_bottom = smooth([C, (122, 245), (108, 245), (70, 230), TIP])
    z2_top = smooth([Z2, (170, 135), (196, 128), T4])
    z4_top = smooth([T4, (240, 123), (262, 121), F])
    z6_edge = smooth([F, (300, 118), (333, 110), (350, 86), (386, 49)]) + [(465, 40), (500, 75), L]
    z9_edge = smooth([H, (290, 256), (297, 266), O])
    z11_edge = ([L] + smooth([(518, 213), (534, 232), (550, 243)]) + [(583, 240)]
                + smooth([(628, 278), (545, 287), (495, 297), Q]))
    z12_edge = [Q, (446, 350), (412, 361)] + smooth([(370, 376), (350, 337), (325, 305), (305, 290), O])

    zones = {
        1: z1_top + [A, C] + z1_bottom[1:],
        2: [Z2, Z1, A, B] + z2_top[::-1][:-1],
        3: [A, B, G, C],
        4: [T4, D, E] + z4_top[::-1][:-1],
        5: [D, E, H, G, B],
        6: z6_edge + [I],
        7: [F, I, J, K],
        8: [I, L, M, J],
        9: [K, J, N] + z9_edge[::-1],
        10: [J, M, P, N],
        11: z11_edge + [P, M],
        12: [O, N, P] + z12_edge[:-1],
    }
    labels = {1: (100, 217), 2: (170, 167), 3: (170, 214), 4: (240, 157), 5: (240, 217),
              6: (413, 93), 7: (323, 165), 8: (423, 165), 9: (323, 242), 10: (423, 242),
              11: (518, 256), 12: (395, 318)}
    return zones, labels


def zonal_frame(old, new, active, t, zones, labels):
    lines = [r"\begin{tikzpicture}[x=0.025cm, y=-0.025cm]",
             r"\useasboundingbox (22,30) rectangle (638,386);"]
    for z, outline in zones.items():
        level = mix(old[z - 1], new[z - 1], t) * SHADE
        lines.append(rf"\filldraw[fill=black!{level:.1f}] {polyline(outline)} -- cycle;")
    if active is not None:
        lines.append(rf"\draw[line width=1.3pt] {polyline(zones[active + 1])} -- cycle;")
    for z, xy in labels.items():
        lines.append(rf"\node at {pt(xy)} {{Zone {z}}};")
    lines += bit_cells(32, 65, 17.25, 25, old, new, active, t, index_labels=True)
    lines.append(r"\node at (135,102) {Genotype};")
    lines.append(r"\node at (240,282) {Phenotype};")
    lines.append(r"\end{tikzpicture}")
    return "\n".join(lines)


def zonal_frames(rng):
    start = [int(c) for c in "101010110100"]  # the example in the journal figure
    flips = closed_flip_sequence(rng.sample(range(12), 8), rng)
    zones, labels = zonal_geometry()
    for old, new, bit, t, ms in timeline(start, flips):
        yield zonal_frame(old, new, bit, t, zones, labels), ms


# --------------------------------------------------------------------------- basin

BASIN_CX, BASIN_CY = 8.3, 0.55   # centre of the basin's top face
RX, RY = 2.6, 0.82               # semi-axes of the largest top face (projected)
BASIN_LEFT = -0.25               # left edge of the figure's bounding box


def snap(x):
    """Centre a vertical line on a pixel boundary so it stays crisp as the basin resizes."""
    px = 2.54 / DPI
    return BASIN_LEFT + round((x - BASIN_LEFT) / px) * px


def depth_of(index):
    return 0.45 + 0.2 * index    # drawn depth, d0 ... dmax


def scale_of(index):
    return math.sqrt(0.42 + index * 0.58 / 7)  # linear scale giving evenly spaced areas


def blob(scale, n=240):
    """Top-face outline: a gently irregular kidney shape, seen obliquely."""
    points = []
    for k in range(n):
        th = 2 * math.pi * k / n
        dent = math.atan2(math.sin(th - 1.25), math.cos(th - 1.25))
        f = (1 + 0.07 * math.cos(2 * th) + 0.05 * math.sin(3 * th + 0.5)
             - 0.14 * math.exp(-dent * dent / 0.10))
        points.append((BASIN_CX + RX * scale * math.cos(th) * f,
                       BASIN_CY + RY * scale * math.sin(th) * f))
    return points


def basin_phenotype(depth, scale):
    top = blob(scale)
    bottom = [(x, y - depth) for x, y in top]
    n = len(top)
    left = min(range(n), key=lambda i: top[i][0])
    right = max(range(n), key=lambda i: top[i][0])
    lines = [rf"\draw[dashed] {polyline(bottom)} -- cycle;"]
    for i in (left, right, int(n * 0.77), int(n * 0.88)):  # silhouettes, then two front edges
        x = snap(top[i][0])
        lines.append(rf"\draw {pt((x, top[i][1]))} -- {pt((x, bottom[i][1]))};")
    lines.append(rf"\filldraw[fill=black!{SHADE}] {polyline(top)} -- cycle;")
    lines.append(rf"\node at {pt((BASIN_CX + 0.1, BASIN_CY))} {{$a$}};")
    # depth dimension, measured from the left-hand rim
    x, y = snap(top[left][0] - 0.4), top[left][1]
    lines.append(rf"\draw {pt((x, y))} -- {pt((x, y - depth))};")
    for yy in (y, y - depth):
        lines.append(rf"\draw ({x - 0.13:.2f},{yy:.2f}) -- ({x + 0.13:.2f},{yy:.2f});")
    lines.append(rf"\node[left] at ({x - 0.08:.2f},{y - depth / 2:.2f}) {{$d$}};")
    return lines


def mapping_table(old_idx, new_idx, t):
    cols = (12.95, 14.85, 16.75)
    rules = (13.9, 15.8)
    header_y, row0, dy = 1.95, 1.4, 0.5
    lines = []
    # highlight the decoded depth and area rows, cross-fading as they change
    for col, key in ((1, 0), (2, 1)):
        for idx, weight in ((old_idx[key], 1 - t), (new_idx[key], t)):
            if old_idx[key] == new_idx[key]:
                weight = 1.0
            y = row0 - idx * dy
            x0 = rules[0] if col == 1 else rules[1]
            x1 = rules[1] if col == 1 else 17.7
            lines.append(rf"\fill[black!{HIGHLIGHT * weight:.1f}] ({x0 + 0.04:.2f},{y - dy / 2 + 0.03:.2f})"
                         rf" rectangle ({x1 - 0.04:.2f},{y + dy / 2 - 0.03:.2f});")
    for x, text in zip(cols, ("Bit-string", "Depth", "Area")):
        lines.append(rf"\node[strut] at ({x},{header_y}) {{{text}}};")
    lines.append(rf"\draw (12.0,{header_y - 0.3:.2f}) -- (17.7,{header_y - 0.3:.2f});")
    for x in rules:
        lines.append(rf"\draw ({x},{header_y + 0.3:.2f}) -- ({x},{row0 - 7 * dy - 0.3:.2f});")
    for i, code in enumerate(GRAY_CODE):
        y = row0 - i * dy
        if i == 0:
            depth, area = r"d_0", r"a_0"
        elif i == 7:
            depth, area = r"d_{\max}", r"a_{\max}"
        else:
            step = "" if i == 1 else str(i)
            depth, area = rf"d_0 + {step}\epsilon_d", rf"a_0 + {step}\epsilon_a"
        lines.append(rf"\node[strut] at ({cols[0]},{y:.2f}) {{{code}}};")
        lines.append(rf"\node[strut] at ({cols[1]},{y:.2f}) {{${depth}$}};")
        lines.append(rf"\node[strut] at ({cols[2]},{y:.2f}) {{${area}$}};")
    return lines


def basin_frame(old, new, active, t):
    decode = lambda g: (GRAY_CODE.index("".join(map(str, g[:3]))),
                        GRAY_CODE.index("".join(map(str, g[3:]))))
    old_idx, new_idx = decode(old), decode(new)
    depth = mix(depth_of(old_idx[0]), depth_of(new_idx[0]), t)
    scale = mix(scale_of(old_idx[1]), scale_of(new_idx[1]), t)
    w, cy = 0.62, 0.0
    brace = r"decorate, decoration={brace, amplitude=5pt, raise=3pt}"
    lines = [r"\begin{tikzpicture}",
             rf"\useasboundingbox ({BASIN_LEFT},-3.45) rectangle (17.95,2.45);"]
    lines += bit_cells(0, cy - w / 2, w, w, old, new, active, t)
    lines.append(rf"\draw[{brace}] (0,{cy + w / 2}) -- ({3 * w},{cy + w / 2}) node[midway, above=9pt, strut] {{Depth}};")
    lines.append(rf"\draw[{brace}] ({3 * w},{cy + w / 2}) -- ({6 * w},{cy + w / 2}) node[midway, above=9pt, strut] {{Area}};")
    lines.append(rf"\draw[{brace}] ({6 * w},{cy - w / 2}) -- (0,{cy - w / 2}) node[midway, below=9pt, font=\large] {{Genotype}};")
    lines.append(rf"\node at ({BASIN_CX},2.0) {{Detention Basin}};")
    lines += basin_phenotype(depth, scale)
    lines.append(rf"\draw[{brace}] ({BASIN_CX + 3.1},-2.55) -- ({BASIN_CX - 3.1},-2.55) node[midway, below=9pt, font=\large] {{Phenotype}};")
    lines += mapping_table(old_idx, new_idx, t)
    lines.append(rf"\draw[{brace}] (17.7,-2.55) -- (12.0,-2.55) node[midway, below=9pt, font=\large] {{Binary--Discrete Mapping}};")
    lines.append(r"\end{tikzpicture}")
    return "\n".join(lines)


def basin_frames(rng):
    start = [int(c) for c in "101011"]  # depth d0 + 6 eps_d, area a0 + 2 eps_a
    flips = closed_flip_sequence(list(range(6)), rng)
    for old, new, bit, t, ms in timeline(start, flips):
        yield basin_frame(old, new, bit, t), ms


# --------------------------------------------------------------------------- build

def render(name, frames):
    BUILD.mkdir(exist_ok=True)
    pictures, durations = zip(*frames)
    tex = BUILD / f"{name}.tex"
    tex.write_text(PREAMBLE + "\n\n".join(pictures) + "\n\\end{document}\n", encoding="utf-8")
    subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
                    f"-output-directory={BUILD}", str(tex)], check=True, stdout=subprocess.DEVNULL)
    for old in BUILD.glob(f"{name}-*.png"):
        old.unlink()
    subprocess.run(["pdftoppm", "-r", str(DPI), "-gray", "-png",
                    str(BUILD / f"{name}.pdf"), str(BUILD / name)], check=True)
    pages = sorted(BUILD.glob(f"{name}-*.png"))
    images = [Image.open(p).convert("L") for p in pages]
    OUT.mkdir(parents=True, exist_ok=True)
    images[0].save(OUT / f"{name}.png", optimize=True)
    # 64 grey levels keep the anti-aliasing smooth at roughly a third of the file size
    step = 255 / (GREY_LEVELS - 1)
    images = [im.point(lambda v: round(round(v / step) * step)) for im in images]
    images[0].save(OUT / f"{name}.gif", save_all=True, append_images=images[1:],
                   duration=list(durations), loop=0, disposal=1)
    size = (OUT / f"{name}.gif").stat().st_size / 1024
    print(f"{name}: {len(images)} frames, {images[0].size[0]}x{images[0].size[1]} px, "
          f"{sum(durations) / 1000:.1f} s loop, {size:.0f} KB")


if __name__ == "__main__":
    render("genotype-zonal", zonal_frames(random.Random(7)))
    render("genotype-basin", basin_frames(random.Random(11)))
