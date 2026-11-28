"""Redraw Figure 2 (Refine): B&W circle style, paired-move direction matched
to the seven-value example.

Writes:
  paper/figs/v5_seven.pdf

With --preview, also writes paper/figs/fig2_preview_2026-11-27.png
(figure plus caption, about paper text width, 200 dpi). Does not edit
any manuscript .tex.
"""
from __future__ import annotations

from fractions import Fraction as F
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
from matplotlib.path import Path as MplPath

# Paper uses Latin Modern. Mathtext 'cm' is Computer Modern, the same family,
# so $L$, $X$, $R$ and the SSE fractions match the other serif labels.
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Latin Modern Roman"],
    "mathtext.fontset": "cm",
})

ROOT = Path(__file__).resolve().parents[1]
FIGDIR = ROOT / "paper" / "figs"
PREVIEW = ROOT / "paper" / "figs" / "fig2_preview_2026-11-27.png"

# ---------------------------------------------------------------------------
# Arithmetic (exact fractions). Fail the build if the example is wrong.
# ---------------------------------------------------------------------------
def sse(parts):
    t = F(0)
    for b in parts:
        m = F(sum(b), len(b))
        t += sum((F(v) - m) ** 2 for v in b)
    return t

BEFORE = [[3, 4], [6, 7, 9], [14, 22]]
AFTER = [[3, 4, 6], [7, 9, 14], [22]]
SINGLES = {
    "give left":  [[3, 4, 6], [7, 9], [14, 22]],
    "take left":  [[3], [4, 6, 7, 9], [14, 22]],
    "give right": [[3, 4], [6, 7], [9, 14, 22]],
    "take right": [[3, 4], [6, 7, 9, 14], [22]],
}

def check_arithmetic():
    s0, s1 = sse(BEFORE), sse(AFTER)
    assert s0 == F(223, 6), s0
    assert s1 == F(92, 3), s1
    assert s0 - s1 == F(13, 2)
    for name, parts in SINGLES.items():
        s = sse(parts)
        assert s > s0, f"{name} does not raise SSE: {s} vs {s0}"
    return s0, s1

# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------
R = 0.23
STEP = 0.58
BIN_PAD_X = 0.11
BIN_PAD_Y = 0.12
BIN_GAP = 0.30
ROW_DY = 1.62
STROKE = 1.15
SHADE = "#d0d0d0"
EDGE = "black"
# Horizontal gap from the left panel's right box edge to the right panel's
# first box edge. "Before" sits in this gap.
PANEL_GAP = 1.70
ARROW_KW = dict(
    arrowstyle="-|>,head_length=0.38,head_width=0.28",
    mutation_scale=9,
    lw=1.2,
    color="black",
    shrinkA=0.0,
    shrinkB=0.0,
    zorder=5,
    clip_on=False,
)
LABEL_FS = 7.5
SSE_FS = 9.5
BINNAME_FS = 8


def bin_box(ax, xs, y, fc="none"):
    x0, x1 = min(xs) - R - BIN_PAD_X, max(xs) + R + BIN_PAD_X
    y0, y1 = y - R - BIN_PAD_Y, y + R + BIN_PAD_Y
    box = FancyBboxPatch(
        (x0, y0), x1 - x0, y1 - y0,
        boxstyle="round,pad=0.02,rounding_size=0.08",
        facecolor=fc, edgecolor=EDGE, linewidth=STROKE,
        mutation_aspect=1, zorder=1, clip_on=False,
    )
    ax.add_patch(box)
    return x0, x1, y0, y1


def unit(ax, x, y, shade=False, label=None):
    circ = Circle(
        (x, y), R,
        facecolor=SHADE if shade else "white",
        edgecolor=EDGE,
        linewidth=1.7 if shade else STROKE,
        zorder=3, clip_on=False,
    )
    ax.add_patch(circ)
    if label is not None:
        ax.text(x, y, str(label), ha="center", va="center",
                fontsize=8, color="black", zorder=4, fontfamily="serif",
                clip_on=False)
    return x, y


def transfer_arrow(ax, x_src, y, dest_box, going):
    """Arrow from a moving unit to just inside the receiving bin's wall.

    going = +1: the unit moves right and enters the dest bin at its left wall.
    going = -1: the unit moves left and enters the dest bin at its right wall.

    The tip stops in the pad just past the wall, not on an occupied circle,
    so the move cannot be read as a swap.
    """
    x0, x1, y0, y1 = dest_box
    # Stop in the receiving bin's pad, just inside the wall and short of
    # the first occupied circle (BIN_PAD_X = 0.11).
    inset = 0.10
    x_end = (x0 + inset) if going > 0 else (x1 - inset)
    y_s = y + R + 0.04
    y_e = y + R + 0.02
    h = 0.30
    # Cubic bezier with a nearly horizontal arrival so the head points
    # into the bin. (A circular arc on this short gap either leaves the
    # head pointing at the wall or is too short to hold an arrowhead.)
    c1 = (x_src + 0.18 * going, y_s + h)
    c2 = (x_end - 0.22 * going, y_e + 0.05)
    path = MplPath(
        [(x_src, y_s), c1, c2, (x_end, y_e)],
        [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4],
    )
    ax.add_patch(FancyArrowPatch(path=path, **ARROW_KW))


def group_xs(n, x_left):
    return [x_left + i * STEP for i in range(n)]


def draw_schematic_row(ax, y, nL=3, nX=4, nR=3, x_left=0.0):
    xL = group_xs(nL, x_left)
    xX = group_xs(nX, xL[-1] + STEP + BIN_GAP)
    xR = group_xs(nR, xX[-1] + STEP + BIN_GAP)
    boxes = {
        "L": bin_box(ax, xL, y),
        "X": bin_box(ax, xX, y),
        "R": bin_box(ax, xR, y),
    }
    return xL, xX, xR, boxes


def draw_left(ax):
    """Three schematic rows. Bottom row is give-left / take-right."""
    y0, y1, y2 = 0.0, -ROW_DY, -2 * ROW_DY

    # (top) one transfer: X gives its right-edge unit to R
    L, X, Rgt, boxes = draw_schematic_row(ax, y0)
    for x in (*L, *X, *Rgt):
        unit(ax, x, y0, shade=(x == X[-1]))
    transfer_arrow(ax, X[-1], y0, boxes["R"], going=+1)
    # L, X, R under this row only
    y_lab = boxes["L"][2] - 0.11  # just below the box
    for name, xs, box in (("L", L, boxes["L"]), ("X", X, boxes["X"]), ("R", Rgt, boxes["R"])):
        xm = 0.5 * (box[0] + box[1])
        ax.text(xm, y_lab, rf"${name}$", ha="center", va="top",
                fontsize=BINNAME_FS, color="black", fontfamily="serif",
                clip_on=False)

    # (middle) two inward transfers: X takes from L and from R
    L, X, Rgt, boxes = draw_schematic_row(ax, y1)
    movers = {L[-1], Rgt[0]}
    for x in (*L, *X, *Rgt):
        unit(ax, x, y1, shade=(x in movers))
    transfer_arrow(ax, L[-1], y1, boxes["X"], going=+1)
    transfer_arrow(ax, Rgt[0], y1, boxes["X"], going=-1)

    # (bottom) paired give-left / take-right, matching the seven-value example
    L, X, Rgt, boxes = draw_schematic_row(ax, y2)
    movers = {X[0], Rgt[0]}
    for x in (*L, *X, *Rgt):
        unit(ax, x, y2, shade=(x in movers))
    transfer_arrow(ax, X[0], y2, boxes["L"], going=-1)   # 6: X's left unit joins L
    transfer_arrow(ax, Rgt[0], y2, boxes["X"], going=-1) # 14: R's left unit joins X

    xs = [0.0 - R - BIN_PAD_X, Rgt[-1] + R + BIN_PAD_X]
    # Room above for curved arrows; bottom is the last row of boxes
    ys = [y2 - R - BIN_PAD_Y, y0 + R + BIN_PAD_Y + 0.42]
    return xs, ys


def pack_row(ax, parts, x_left, y, movers):
    """Pack circles within bins, with a visible gap between bins. Return {value: x}, last box."""
    pos = {}
    x = x_left
    last_box = None
    first_box = None
    for b in parts:
        xs = [x + i * STEP for i in range(len(b))]
        last_box = bin_box(ax, xs, y)
        if first_box is None:
            first_box = last_box
        for xv, v in zip(xs, b):
            unit(ax, xv, y, shade=(v in movers), label=v)
            pos[v] = xv
        x = xs[-1] + STEP + BIN_GAP
    return pos, first_box, last_box


def draw_right(ax, x_left, y_top):
    """Before / after of the seven-value window, with small row labels and SSE."""
    movers = {6, 14}
    y_before = y_top
    y_after = y_top - ROW_DY

    pos_b, box_b0, box_b1 = pack_row(ax, BEFORE, x_left, y_before, movers)
    pos_a, box_a0, box_a1 = pack_row(ax, AFTER, x_left, y_after, movers)

    for v in (6, 14):
        ax.add_patch(FancyArrowPatch(
            (pos_b[v], y_before - R - 0.05),
            (pos_a[v], y_after + R + 0.05),
            **ARROW_KW,
        ))

    # Small black labels: Before / After at the left, SSE at the right.
    x_lab_l = box_b0[0] - 0.12
    x_lab_r = max(box_b1[1], box_a1[1]) + 0.12
    ax.text(x_lab_l, y_before, "Before", ha="right", va="center",
            fontsize=LABEL_FS, color="black", fontfamily="serif", clip_on=False)
    ax.text(x_lab_l, y_after, "After", ha="right", va="center",
            fontsize=LABEL_FS, color="black", fontfamily="serif", clip_on=False)
    ax.text(x_lab_r, y_before, r"SSE $= \frac{223}{6}$", ha="left", va="center",
            fontsize=SSE_FS, color="black", fontfamily="serif", clip_on=False)
    ax.text(x_lab_r, y_after, r"SSE $= \frac{92}{3}$", ha="left", va="center",
            fontsize=SSE_FS, color="black", fontfamily="serif", clip_on=False)

    xs_all = list(pos_b.values()) + list(pos_a.values())
    # SSE text extends ~1.35 data units to the right; Before ~0.70 to the left.
    xs_ext = [x_lab_l - 0.70, x_lab_r + 1.65]
    ys_ext = [y_after - R - BIN_PAD_Y, y_before + R + BIN_PAD_Y]
    return xs_ext, ys_ext


# Fractions now appear on the figure, so the caption no longer repeats them.
CAPTION = (
    "Refine scores eight moves of the two cuts in a three-bin window. "
    "The rows at the left show a single transfer, two inward transfers, "
    "and a paired give-and-take, with shaded units moving. "
    "On the seven values at the right, that paired move lowers the SSE, "
    "while each of the four single transfers raises it."
)


def build_figure(with_caption=False):
    check_arithmetic()

    fig_w = 6.50
    fig_h = 2.95 if not with_caption else 4.05
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_aspect("equal")
    ax.axis("off")

    left_xs, left_ys = draw_left(ax)

    left_h = left_ys[1] - left_ys[0]
    right_h = ROW_DY + 2 * (R + BIN_PAD_Y)
    y_top_right = left_ys[1] - (left_h - right_h) / 2 - (R + BIN_PAD_Y)

    x_right = left_xs[1] + PANEL_GAP
    right_xs, right_ys = draw_right(ax, x_right, y_top_right)

    # Faint divider in the leftover space between the left boxes and "Before".
    # "Before" is right-aligned at (first right box left - 0.12).
    # Midpoint of the whitespace between the left panel and the Before/After column.
    x_lab_l = x_right - R - BIN_PAD_X - 0.12
    before_w = 0.62
    x_div = 0.5 * (left_xs[1] + (x_lab_l - before_w))
    y_mid = 0.5 * (left_ys[0] + left_ys[1])
    half = 0.38 * (left_ys[1] - left_ys[0])
    ax.plot(
        [x_div, x_div], [y_mid - half, y_mid + half],
        color="#c8c8c8", lw=0.6, solid_capstyle="round", zorder=0, clip_on=False,
    )

    xmin = min(left_xs[0], right_xs[0]) - 0.28
    xmax = max(left_xs[1], right_xs[1]) + 0.28
    ymin = min(left_ys[0], right_ys[0]) - 0.32
    ymax = max(left_ys[1], right_ys[1]) + 0.50

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    fig.subplots_adjust(
        left=0.03, right=0.97,
        top=0.96 if not with_caption else 0.97,
        bottom=0.04 if not with_caption else 0.30,
    )
    return fig, ax


def main():
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--preview", action="store_true",
                   help="also write the captioned PNG preview")
    args = p.parse_args()

    FIGDIR.mkdir(parents=True, exist_ok=True)

    fig, _ = build_figure(with_caption=False)
    pdf_path = FIGDIR / "v5_seven.pdf"
    fig.savefig(pdf_path, format="pdf", dpi=300, bbox_inches="tight", pad_inches=0.20)
    plt.close(fig)
    print("wrote", pdf_path)

    if args.preview:
        fig, ax = build_figure(with_caption=True)
        fig.text(
            0.05, 0.03, r"$\mathbf{Figure\ 2.}$  " + CAPTION,
            ha="left", va="bottom",
            fontsize=9, fontfamily="serif",
            wrap=True,
            linespacing=1.4,
        )
        png_path = PREVIEW
        fig.savefig(png_path, format="png", dpi=200, bbox_inches="tight", pad_inches=0.20,
                    facecolor="white")
        plt.close(fig)
        print("wrote", png_path)

    print("caption:")
    print(CAPTION)


if __name__ == "__main__":
    main()
