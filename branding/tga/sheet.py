#!/home/zamme/.venvs/browser/bin/python
"""sheet.py OUT.png file.svg... -> one contact sheet: each SVG large, plus symbols at 64/32/16 px."""
import sys, os, base64
from playwright.sync_api import sync_playwright

out, files = sys.argv[1], sys.argv[2:]
cells = ""
for f in files:
    data = base64.b64encode(open(f, "rb").read()).decode()
    src = f"data:image/svg+xml;base64,{data}"
    small = ""
    if "symbol" in f:
        small = "".join(f"<img src='{src}' style='height:{s}px;margin-left:10px'>" for s in (64, 32, 16))
    cells += (f"<figure style='margin:0'><div style='display:flex;align-items:flex-end;background:#fff;padding:10px;"
              f"border:1px solid #999'><img src='{src}' style='height:220px'>{small}</div>"
              f"<figcaption>{os.path.basename(f)}</figcaption></figure>")
html = ("<body style='margin:0;background:#ddd;font:14px sans-serif;display:flex;flex-wrap:wrap;gap:14px;"
        f"padding:14px;width:1560px'>{cells}</body>")
with sync_playwright() as p:
    b = p.chromium.launch(env={**os.environ, "LD_LIBRARY_PATH": os.path.expanduser("~/.venvs/browser-libs/lib")})
    pg = b.new_page(viewport={"width": 1590, "height": 600})
    pg.set_content(html)
    pg.wait_for_timeout(100)
    pg.screenshot(path=out, full_page=True)
    b.close()
print(out)
