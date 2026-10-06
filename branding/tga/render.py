#!/home/zamme/.venvs/browser/bin/python
"""render.py OUT_DIR HEIGHT file.svg...  -> PNG of each SVG at HEIGHT px tall on white, in one browser session."""
import sys, os, re, base64
from playwright.sync_api import sync_playwright

out, h, files = sys.argv[1], int(sys.argv[2]), sys.argv[3:]
os.makedirs(out, exist_ok=True)
with sync_playwright() as p:
    b = p.chromium.launch(env={**os.environ, "LD_LIBRARY_PATH": os.path.expanduser("~/.venvs/browser-libs/lib")})
    pg = b.new_page()
    for f in files:
        src = open(f).read()
        vw, vh = map(float, re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', src).groups())
        w = round(h * vw / vh)
        pg.set_viewport_size({"width": w, "height": h})
        data = base64.b64encode(src.encode()).decode()
        pg.set_content(f"<body style='margin:0;background:#fff'><img style='display:block;width:{w}px;"
                       f"height:{h}px' src='data:image/svg+xml;base64,{data}'></body>")
        pg.wait_for_timeout(50)
        o = os.path.join(out, os.path.basename(f).replace(".svg", ".png"))
        pg.screenshot(path=o)
        print(o, w, "x", h)
    b.close()
