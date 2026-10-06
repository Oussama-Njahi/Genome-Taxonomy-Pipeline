"""Turn text set in the interface's fonts into SVG path data (no live <text>).

    from glyphs import text_path
    d, width = text_path("TGA", "bodoni", size=120, x=0, y=200, wght=500, opsz=96, track=0)
"""
from functools import lru_cache
from fontTools.ttLib import TTFont
from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.pens.transformPen import TransformPen
from fontTools.varLib.instancer import instantiateVariableFont

FONTS = "/home/zamme/stage/webapp/frontend/fonts/"
FILES = {"bodoni": "bodoni-moda-var.woff2", "spectral": "spectral-500.woff2",
         "spectral-600": "spectral-600.woff2", "fraktur": "unifraktur-maguntia.woff2"}


@lru_cache(None)
def _font(name, wght=None, opsz=None):
    f = TTFont(FONTS + FILES[name])
    if "fvar" in f:
        axes = {"wght": wght or 400, "opsz": opsz or 96}
        f = instantiateVariableFont(f, axes)
    return f


def text_path(text, font="bodoni", size=100, x=0, y=0, wght=None, opsz=None, track=0, kern=None):
    """Return (path d, advance width). y is the baseline; track is extra space between letters in px.
    kern: optional list of per-gap adjustments in px (len(text)-1)."""
    f = _font(font, wght, opsz)
    upm = f["head"].unitsPerEm
    s = size / upm
    gs = f.getGlyphSet()
    cmap = f.getBestCmap()
    hmtx = f["hmtx"]
    pen = SVGPathPen(gs, ntos=lambda v: ("%.2f" % v).rstrip("0").rstrip("."))
    cx = x
    for i, ch in enumerate(text):
        g = cmap[ord(ch)]
        tp = TransformPen(pen, (s, 0, 0, -s, cx, y))
        gs[g].draw(tp)
        cx += hmtx[g][0] * s + track
        if kern and i < len(kern):
            cx += kern[i]
    return pen.getCommands(), cx - x - track


def bounds(text, font="bodoni", size=100, wght=None, opsz=None, track=0):
    """Ink bounds (xmin, ymin, xmax, ymax) relative to origin/baseline, y down."""
    from fontTools.pens.boundsPen import BoundsPen
    f = _font(font, wght, opsz)
    s = size / f["head"].unitsPerEm
    gs, cmap, hmtx = f.getGlyphSet(), f.getBestCmap(), f["hmtx"]
    bp = BoundsPen(gs)
    cx = 0
    for ch in text:
        g = cmap[ord(ch)]
        gs[g].draw(TransformPen(bp, (s, 0, 0, -s, cx, 0)))
        cx += hmtx[g][0] * s + track
    return bp.bounds
