"""TGA — Dendrogram: final masters (colour), lockups and small-size cuts, written to final/."""
import os
from fontTools.pens.recordingPen import RecordingPen
from fontTools.pens.transformPen import TransformPen
from glyphs import text_path, bounds, _font

OUT = "final"
os.makedirs(OUT, exist_ok=True)

INK = "#2A211B"       # encre (interface)
BLUE = "#2F4A63"      # bleu de Prusse (interface: "same taxon")
PAPER = "#F6EEE3"     # papier de planche
NAME = "TAXO GENOMIC ANALYSES"
WORD = dict(wght=600, opsz=36)


def svg(w, h, body, title, bg=None):
    back = f'<rect width="{w:g}" height="{h:g}" fill="{bg}"/>' if bg else ""
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w:g} {h:g}" width="{w:g}" height="{h:g}" '
            f'role="img" aria-labelledby="title"><title id="title">{title}</title>{back}{body}</svg>\n')


def rect(x0, y0, x1, y1):
    f = lambda v: f"{round(v, 2):g}"
    return f"M{f(x0)} {f(y0)}H{f(x1)}V{f(y1)}H{f(x0)}Z"


def apex_x(ch, size, **kw):
    """x of the topmost point of a glyph's outline (the A's apex is left of its ink centre in Bodoni)."""
    font = _font("bodoni", kw.get("wght"), kw.get("opsz"))
    s = size / font["head"].unitsPerEm
    rec = RecordingPen()
    font.getGlyphSet()[font.getBestCmap()[ord(ch)]].draw(TransformPen(rec, (s, 0, 0, -s, 0, 0)))
    pts = [p for _, args in rec.value for p in args]
    top = min(y for _, y in pts)
    xs = [x for x, y in pts if y <= top + 0.5]
    return (min(xs) + max(xs)) / 2


def dendrogram(size, line, x0=0, base=0, track=None, gap=None, word=WORD, tree_fill=BLUE, ink=INK, tk=None):
    """Leaves = the letters T, G, A. Returns (body, ink_width, top_y_of_root, cap_y).

    Proportions scale with the type size so every lockup has the same drawing."""
    k = size / 150
    tk = k if tk is None else tk
    track = 14 * k if track is None else track
    gap = 14 * k if gap is None else gap
    tb = bounds("TGA", "bodoni", size, track=track, **word)
    d_tga, _ = text_path("TGA", "bodoni", size, x0 - tb[0], base, track=track, **word)
    tips, cx = [], x0 - tb[0]
    for ch in "TGA":
        b = bounds(ch, "bodoni", size, **word)
        tips.append(cx + (apex_x(ch, size, **word) if ch == "A" else (b[0] + b[2]) / 2))
        cx += text_path(ch, "bodoni", size, 0, 0, **word)[1] + track
    cap = base + tb[1]
    tip = cap - gap
    h_ga, h_root = tip - 26 * tk, tip - 58 * tk
    t, g, a = tips
    m_ga = (g + a) / 2
    m_root = (t + m_ga) / 2
    hl = line / 2
    root_top = h_root - 22 * tk
    tree = (rect(g - hl, h_ga - hl, g + hl, tip) + rect(a - hl, h_ga - hl, a + hl, tip) +
            rect(g - hl, h_ga - hl, a + hl, h_ga + hl) +
            rect(t - hl, h_root - hl, t + hl, tip) + rect(m_ga - hl, h_root - hl, m_ga + hl, h_ga - hl) +
            rect(t - hl, h_root - hl, m_ga + hl, h_root + hl) + rect(m_root - hl, root_top, m_root + hl, h_root - hl))
    body = f'<path fill="{tree_fill}" d="{tree}"/><path fill="{ink}" d="{d_tga}"/>'
    return body, tb[2] - tb[0], root_top, cap


def write(name, w, h, body, title, bg=None):
    open(os.path.join(OUT, name), "w").write(svg(w, h, body, title, bg))


# ---- 1. Wordmark (primary) ------------------------------------------------------------------------------------
PAD = 16
size, line = 150, 6
probe = dendrogram(size, line)
top = probe[2]
W = round(probe[1] + 2 * PAD)
H = round(-top + 2 * PAD)
body, *_ = dendrogram(size, line, x0=PAD, base=PAD - top)
write("tga-wordmark-color.svg", W, H, body, "TGA — Taxo Genomic Analyses")
body_rev, *_ = dendrogram(size, line, x0=PAD, base=PAD - top, word=dict(wght=500, opsz=36), tree_fill="#C9D3DB", ink=PAPER)
write("tga-wordmark-reversed.svg", W, H, body_rev, "TGA — reversed", bg=INK)

# ---- 2. Stacked: wordmark + name ------------------------------------------------------------------------------
name_size, name_track = 21, 4.6
nb = bounds(NAME, "spectral-600", name_size, track=name_track)
nw = nb[2] - nb[0]
SW = round(max(probe[1], nw) + 2 * PAD)
wx = (SW - probe[1]) / 2
body, *_ = dendrogram(size, line, x0=wx, base=PAD - top)
name_base = PAD - top + 52
d_name, _ = text_path(NAME, "spectral-600", name_size, (SW - nw) / 2 - nb[0], name_base, track=name_track)
SH = round(name_base + PAD)
write("tga-stacked-color.svg", SW, SH, body + f'<path fill="{INK}" d="{d_name}"/>', "TGA — Taxo Genomic Analyses, stacked")

# ---- 3. Horizontal: wordmark | name on three lines ------------------------------------------------------------
base = PAD - top
cap = probe[3] + base
lines = ["TAXO", "GENOMIC", "ANALYSES"]
ls, lt = 25, 5
rule_x = PAD + probe[1] + 34
tx = rule_x + 30
lead = (base - cap) / 3
parts, tw = "", 0
for i, word in enumerate(lines):
    b = bounds(word, "spectral-600", ls, track=lt)
    d, _ = text_path(word, "spectral-600", ls, tx - b[0], cap + lead * (i + 1) - 4 + (i == 0) * 0, track=lt)
    parts += d
    tw = max(tw, b[2] - b[0])
body, *_ = dendrogram(size, line, x0=PAD, base=base)
rule = rect(rule_x, cap, rule_x + 1.5, base)
HW = round(tx + tw + PAD)
write("tga-horizontal-color.svg", HW, H, body + f'<path fill="{INK}" d="{rule}{parts}"/>', "TGA — Taxo Genomic Analyses, horizontal")

# ---- 4. Small-size cuts (256 square) --------------------------------------------------------------------------
# a) symbol for 32-64 px: same drawing, sturdier type (opsz 6, heavier) and thicker branches
sw = dict(wght=800, opsz=6)
p = dendrogram(96, 12, word=sw, track=6, gap=10, tk=1.0)
k = 256
x = (k - p[1]) / 2
y_base = 124 - p[2] / 2                         # centre root-top..baseline, 4 px above geometric centre
body, *_ = dendrogram(96, 12, x0=x, base=y_base, word=sw, track=6, gap=10, tk=1.0)
write("tga-symbol-color.svg", 256, 256, body, "TGA — symbol")


# b) favicon for 16-32 px: the tree alone on a paper tile (no ground line: with it the tree read as a building)
def favicon(L, leaves=False, tile=True):
    hl = L / 2
    t, g, a = 58, 150, 206
    tip, h_ga, h_root, root_top = 196 if not leaves else 168, 128, 84, 44
    m_ga = (g + a) / 2
    m_root = (t + m_ga) / 2
    tree = (rect(g - hl, h_ga - hl, g + hl, tip) + rect(a - hl, h_ga - hl, a + hl, tip) + rect(g - hl, h_ga - hl, a + hl, h_ga + hl) +
            rect(t - hl, h_root - hl, t + hl, tip) + rect(m_ga - hl, h_root - hl, m_ga + hl, h_ga - hl) +
            rect(t - hl, h_root - hl, m_ga + hl, h_root + hl) + rect(m_root - hl, root_top, m_root + hl, h_root - hl))
    leaf = "".join(rect(x - 18, 184, x + 18, 212) for x in (t, g, a)) if leaves else ""
    back = f'<rect width="256" height="256" rx="48" fill="{PAPER}"/>' if tile else ""
    return back + f'<path fill="{BLUE}" d="{tree}"/>' + (f'<path fill="{INK}" d="{leaf}"/>' if leaf else "")


write("tga-favicon-color.svg", 256, 256, favicon(24, leaves=True), "TGA — favicon")
write("tga-favicon-mark.svg", 256, 256, favicon(24, leaves=True, tile=False), "TGA — favicon mark (no tile)")
print("final masters written:", sorted(os.listdir(OUT)))
