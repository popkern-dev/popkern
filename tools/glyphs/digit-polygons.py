"""Flatten a font's digit glyphs into polygon data for Rhai scripts.

    uv run --with fonttools python tools/glyphs/digit-polygons.py <font.ttf> [chars]

Prints a Rhai map literal: for each character, its advance and its
contours as [x, y] point lists (cap height 1.0, baseline y = 0), outer
contours counter-clockwise, counters (holes) clockwise, each flagged.
The polyhedral dice example embeds this table (generated from Inter
Bold, SIL Open Font License 1.1, https://rsms.me/inter/), so a script
engraves real numerals without a font loader. Quadratic and cubic
segments are flattened to a chord deviation of FLATNESS cap heights;
points closer than MIN_EDGE are merged (the kernel's text.rs bound).
"""
import sys
from fontTools.ttLib import TTFont
from fontTools.pens.recordingPen import DecomposingRecordingPen

FLATNESS = float(__import__("os").environ.get("FLATNESS", "0.02"))
MIN_EDGE = float(__import__("os").environ.get("MIN_EDGE", "0.03"))

def flatten(ops, scale):
    contours, cur = [], []
    def quad(p0, p1, p2):
        n = max(2, int(((abs(p0[0] - 2 * p1[0] + p2[0]) + abs(p0[1] - 2 * p1[1] + p2[1])) * scale / FLATNESS) ** 0.5))
        for i in range(1, n + 1):
            t = i / n
            cur.append(((1 - t) ** 2 * p0[0] + 2 * (1 - t) * t * p1[0] + t * t * p2[0],
                        (1 - t) ** 2 * p0[1] + 2 * (1 - t) * t * p1[1] + t * t * p2[1]))
    def cubic(p0, p1, p2, p3):
        d = max(abs(p0[0] - 2 * p1[0] + p2[0]) + abs(p0[1] - 2 * p1[1] + p2[1]),
                abs(p1[0] - 2 * p2[0] + p3[0]) + abs(p1[1] - 2 * p2[1] + p3[1]))
        n = max(2, int((1.5 * d * scale / FLATNESS) ** 0.5))
        for i in range(1, n + 1):
            t = i / n; s = 1 - t
            cur.append((s**3 * p0[0] + 3 * s * s * t * p1[0] + 3 * s * t * t * p2[0] + t**3 * p3[0],
                        s**3 * p0[1] + 3 * s * s * t * p1[1] + 3 * s * t * t * p2[1] + t**3 * p3[1]))
    for op, args in ops:
        if op == "moveTo":
            cur = [args[0]]
        elif op == "lineTo":
            cur.append(args[0])
        elif op == "qCurveTo":
            pts = list(args); p0 = cur[-1]
            # TrueType implied on-curve points between consecutive offs.
            offs, end = pts[:-1], pts[-1]
            for k, off in enumerate(offs):
                nxt = end if k == len(offs) - 1 else ((off[0] + offs[k + 1][0]) / 2, (off[1] + offs[k + 1][1]) / 2)
                quad(p0, off, nxt); p0 = nxt
        elif op == "curveTo":
            cubic(cur[-1], *args)
        elif op in ("closePath", "endPath"):
            if len(cur) > 2:
                contours.append(cur)
            cur = []
    return contours

def area(c):
    return sum(c[i][0] * c[(i + 1) % len(c)][1] - c[(i + 1) % len(c)][0] * c[i][1] for i in range(len(c))) / 2

def clean(c):
    out = []
    for p in c:
        if not out or ((p[0] - out[-1][0]) ** 2 + (p[1] - out[-1][1]) ** 2) ** 0.5 >= MIN_EDGE:
            out.append(p)
    while len(out) > 3 and ((out[0][0] - out[-1][0]) ** 2 + (out[0][1] - out[-1][1]) ** 2) ** 0.5 < MIN_EDGE:
        out.pop()
    return out if len(out) >= 3 else None

def inside(p, c):
    x, y = p; hit = False
    for i in range(len(c)):
        (x0, y0), (x1, y1) = c[i], c[(i + 1) % len(c)]
        if (y0 > y) != (y1 > y) and x < x0 + (y - y0) * (x1 - x0) / (y1 - y0):
            hit = not hit
    return hit

font = TTFont(sys.argv[1])
chars = sys.argv[2] if len(sys.argv) > 2 else "0123456789"
cap = font["OS/2"].sCapHeight
scale = 1.0 / cap
gs = font.getGlyphSet(); cmap = font.getBestCmap()
rows = []
for ch in chars:
    name = cmap[ord(ch)]
    pen = DecomposingRecordingPen(gs); gs[name].draw(pen)
    contours = [c for c in (clean([(x * scale, y * scale) for x, y in c]) for c in flatten(pen.value, 1.0)) if c]
    entries = []
    for c in contours:
        # Even-odd depth: a contour inside an odd number of others is a counter.
        depth = sum(1 for o in contours if o is not c and inside(c[0], o))
        hole = depth % 2 == 1
        if (area(c) > 0) == hole:
            c = c[::-1]
        pts = ", ".join(f"[{x:.4f}, {y:.4f}]" for x, y in c)
        entries.append(f"#{{ hole: {'true' if hole else 'false'}, pts: [{pts}] }}")
    adv = font["hmtx"][name][0] * scale
    rows.append(f'    "{ch}": #{{ adv: {adv:.4f}, contours: [\n        ' + ",\n        ".join(entries) + "] }")
print("#{\n" + ",\n".join(rows) + "\n}")
