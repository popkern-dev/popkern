"""Render the popkern mark at high resolution (PNG + SVG)."""
import math, sys
from PIL import Image, ImageDraw
out = sys.argv[1]
AMBER = (245, 165, 36)
DARK, LIGHT = (10, 12, 15), (236, 237, 245)
INK_ON_LIGHT = (20, 24, 29)

def render(size, bg, stroke, path):
    SS = 4; N = size * SS; k = N / 32
    img = Image.new("RGBA", (N, N), bg if bg else (0, 0, 0, 0))
    d = ImageDraw.Draw(img)
    P = lambda x, y: (x * k, y * k)
    arc = lambda cx, cy, r, a0, a1: [P(cx + r * math.cos(math.radians(a)), cy + r * math.sin(math.radians(a)))
                                     for a in [a0 + (a1 - a0) * i / 180 for i in range(181)]]

    def contour(o):
        # Offset outline of the mark at signed distance o (sharp corners mitred).
        return (arc(14, 13.5, 7 + o, 180, 270) + [P(25 + o, 6.5 - o)]
                + arc(18, 17.5, 7 + o, 0, 90) + [P(7 - o, 24.5 + o)])

    h = 1.2                                   # half the 2.4 stroke
    mask = Image.new("L", img.size, 0)
    md = ImageDraw.Draw(mask)
    md.polygon(contour(h), fill=255)
    md.polygon(contour(-h), fill=0)
    img.paste(Image.new("RGBA", img.size, stroke), (0, 0), mask)

    ah = 1.9                                  # half the 3.8 amber stroke
    d = ImageDraw.Draw(img)
    d.polygon(arc(14, 13.5, 7 + ah, 180, 270) + arc(14, 13.5, 7 - ah, 270, 180), fill=AMBER)
    for cx, cy in ((7, 13.5), (14, 6.5)):     # round caps
        c = P(cx, cy); r = ah * k
        d.ellipse([c[0] - r, c[1] - r, c[0] + r, c[1] + r], fill=AMBER)
    img.resize((size, size), Image.LANCZOS).save(f"{out}/{path}", optimize=True)

for s in (1024, 512):
    render(s, DARK + (255,), LIGHT, f"popkern-logo-{s}.png")
render(1024, (255, 255, 255, 255), INK_ON_LIGHT, "popkern-logo-light-1024.png")
render(1024, None, LIGHT, "popkern-logo-transparent-dark-bg-1024.png")
render(1024, None, INK_ON_LIGHT, "popkern-logo-transparent-light-bg-1024.png")

svg = lambda bg, stroke: f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 32 32" width="1024" height="1024">
{f'  <rect width="32" height="32" fill="{bg}"/>' if bg else ''}
  <path d="M7 13.5a7 7 0 0 1 7-7h11v11a7 7 0 0 1-7 7H7V13.5Z" fill="none" stroke="{stroke}" stroke-width="2.4"/>
  <path d="M7 13.5a7 7 0 0 1 7-7" fill="none" stroke="#f5a524" stroke-width="3.8" stroke-linecap="round"/>
</svg>
'''
open(f"{out}/popkern-logo.svg", "w").write(svg("#0a0c0f", "#ecedf5"))
open(f"{out}/popkern-logo-transparent-dark-bg.svg", "w").write(svg(None, "#ecedf5"))
open(f"{out}/popkern-logo-transparent-light-bg.svg", "w").write(svg(None, "#14181d"))
print("ok")
