"""Build the three TGA concepts (black first): symbol + horizontal lockup for each."""
import math
from glyphs import text_path, bounds

INK = "#111"
NAME = "TAXO GENOMIC ANALYSES"
T15 = math.tan(math.radians(15))


def svg(w, h, body, title):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w:g} {h:g}" width="{w:g}" height="{h:g}" '
            f'role="img" aria-labelledby="title"><title id="title">{title}</title>{body}</svg>\n')


def rect(x0, y0, x1, y1):
    return f"M{x0:g} {y0:g}H{x1:g}V{y1:g}H{x0:g}Z"


def poly(*pts):
    return "M" + "L".join(f"{x:.2f} {y:.2f}" for x, y in pts) + "Z"


def arc_band(cx, cy, r0, r1, a0, a1):
    """Annular sector between radii r0<r1 from angle a0 to a1 (degrees, 0 = up, clockwise)."""
    def pt(r, a):
        t = math.radians(a)
        return cx + r * math.sin(t), cy - r * math.cos(t)
    large = 1 if (a1 - a0) % 360 > 180 else 0
    (x0, y0), (x1, y1) = pt(r1, a0), pt(r1, a1)
    (x2, y2), (x3, y3) = pt(r0, a1), pt(r0, a0)
    return (f"M{x0:.2f} {y0:.2f}A{r1} {r1} 0 {large} 1 {x1:.2f} {y1:.2f}"
            f"L{x2:.2f} {y2:.2f}A{r0} {r0} 0 {large} 0 {x3:.2f} {y3:.2f}Z")


def centred_text(text, font, size, cx, baseline, **kw):
    b = bounds(text, font, size, kw.get("wght"), kw.get("opsz"), kw.get("track", 0))
    d, _ = text_path(text, font, size, cx - (b[0] + b[2]) / 2, baseline, **kw)
    return d


def lockup(symbol_body, fname, title, sym_w=256):
    """Symbol (256 high) on the left, TGA + name on the right, height 256."""
    x = sym_w + 40
    tb = bounds("TGA", "bodoni", 150, 600, 36, 6)
    nb = bounds(NAME, "spectral-600", 21, track=4.2)
    d_tga, _ = text_path("TGA", "bodoni", 150, x - tb[0], 160, wght=600, opsz=36, track=6)
    d_name, _ = text_path(NAME, "spectral-600", 21, x - nb[0], 206, track=4.2)
    width = x + max(tb[2] - tb[0], nb[2] - nb[0]) + 8
    body = symbol_body + f'<path fill="{INK}" d="{d_tga}{d_name}"/>'
    open(fname, "w").write(svg(width, 256, body, title))


# ---------- A: Divergence ----------
# A Didone A drawn as a split: both strokes leave one apex at exactly 15 deg from vertical (two lineages from a
# common ancestor). Thin left stroke 14, thick right stroke 42 (unequal rates), hairline foot serifs. The crossbar
# is the species threshold: it crosses both lineages and, in v2, runs past them as a measuring line.
def divergence(extend=False, dx=-10):
    top, base = 24, 232
    run = (base - top) * T15                                  # horizontal travel of each stroke
    thick = poly((116, top), (158, top), (158 + run, base), (116 + run, base))
    thin = poly((116, top), (130, top), (130 - run, base), (116 - run, base))
    feet = rect(116 - run - 16, base - 6, 130 - run + 16, base) + rect(116 + run - 14, base - 6, 158 + run + 18, base)
    yb0, yb1 = 150, 164
    if extend:
        bar = rect(116 - run - 16, yb0, 158 + run + 18, yb1)
    else:
        xl = 130 - T15 * (yb1 - top) - 2
        xr = 116 + T15 * (yb0 - top) + 2
        bar = rect(xl, yb0, xr, yb1)
    return f'<g transform="translate({dx} 0)"><path fill="{INK}" d="{thick}{thin}{feet}"/><path fill="{INK}" d="{bar}"/></g>'


open("a-symbol.svg", "w").write(svg(256, 256, divergence(), "TGA — Divergence"))
open("a-symbol-ext.svg", "w").write(svg(256, 256, divergence(True), "TGA — Divergence, open threshold"))
lockup(divergence(), "a-lockup.svg", "TGA — Divergence lockup")


# ---------- B: Chromosome seal ----------
# Circular chromosome: a heavy outer ring (r 120/106) broken at 12 o'clock where a wedge marks the origin of
# replication; an inner track of irregular arcs (genes on the two strands) at r 92/100; TGA in the middle.
GENES = [(18, 52), (60, 74), (96, 141), (150, 166), (188, 236), (246, 262), (278, 300), (318, 344)]


def chromosome():
    cx = cy = 128
    ring = arc_band(cx, cy, 106, 120, 7, 353)
    ori = poly((cx - 10, cy - 124), (cx + 10, cy - 124), (cx, cy - 100))
    genes = "".join(arc_band(cx, cy, 91, 99, a0, a1) for a0, a1 in GENES)
    d_tga = centred_text("TGA", "bodoni", 66, cx, cy + 23, wght=700, opsz=28, track=1)
    return f'<path fill="{INK}" d="{ring}{ori}{genes}{d_tga}"/>'


open("b-symbol.svg", "w").write(svg(256, 256, chromosome(), "TGA — Chromosome"))
lockup(chromosome(), "b-lockup.svg", "TGA — Chromosome lockup")


# ---------- C: Dendrogram ----------
# The wordmark is the leaf row of a dendrogram, as in the pipeline's figures: G and A are sisters (joined low),
# T joins them higher up, a root stub rises from the top bar. Hairlines 6 wide, letters hang from the tips.
def dendrogram(size=150, x0=0, base=230, line=6, gap=14):
    kw = dict(wght=600, opsz=36)
    track = 14
    tb = bounds("TGA", "bodoni", size, track=track, **kw)
    d_tga, _ = text_path("TGA", "bodoni", size, x0 - tb[0], base, track=track, **kw)
    # centre x of each letter's ink
    centres, cx = [], x0 - tb[0]
    for ch in "TGA":
        b = bounds(ch, "bodoni", size, **kw)
        adv = text_path(ch, "bodoni", size, 0, 0, **kw)[1]
        centres.append(cx + (b[0] + b[2]) / 2)
        cx += adv + track
    cap = base + tb[1]                       # top of the caps (tb[1] is negative)
    tip = cap - gap                          # where the branches stop
    h_ga, h_root = tip - 26, tip - 58        # node heights
    t, g, a = centres
    m_ga = (g + a) / 2
    m_root = (t + m_ga) / 2
    hl = line / 2
    tree = (rect(g - hl, h_ga, g + hl, tip) + rect(a - hl, h_ga, a + hl, tip) + rect(g - hl, h_ga - hl, a + hl, h_ga + hl) +
            rect(t - hl, h_root, t + hl, tip) + rect(m_ga - hl, h_root, m_ga + hl, h_ga) +
            rect(t - hl, h_root - hl, m_ga + hl, h_root + hl) + rect(m_root - hl, h_root - 22, m_root + hl, h_root))
    return f'<path fill="{INK}" d="{tree}{d_tga}"/>', tb[2] - tb[0], h_root - 22


body, w, top = dendrogram()
open("c-wordmark.svg", "w").write(svg(w, 256, body, "TGA — Dendrogram"))
# square symbol: the same construction scaled into 256
b2, w2, top2 = dendrogram(size=78, x0=0, base=196, line=6, gap=10)
open("c-symbol.svg", "w").write(svg(256, 256, f'<g transform="translate({(256 - w2) / 2:.1f} 6)">{b2}</g>', "TGA — Dendrogram symbol"))
# lockup: dendrogram wordmark + name beneath, centred
body3, w3, _ = dendrogram(size=120, base=196)
nb = bounds(NAME, "spectral-600", 19, track=3.8)
d_name, _ = text_path(NAME, "spectral-600", 19, (w3 - (nb[2] - nb[0])) / 2 - nb[0], 240, track=3.8)
nw = nb[2] - nb[0]
W = max(w3, nw) + 8
d_name, _ = text_path(NAME, "spectral-600", 19, (W - nw) / 2 - nb[0], 240, track=3.8)
body3 = f'<g transform="translate({(W - w3) / 2:.1f} 0)">{body3}</g>'
open("c-lockup.svg", "w").write(svg(W, 256, body3 + f'<path fill="{INK}" d="{d_name}"/>', "TGA — Dendrogram lockup"))
print("built")
