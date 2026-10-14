"""Crop a transparent render to its subject (2% padding) and downscale it.

    python3 cutout.py <render.png> <out.png> <max width>   (needs Pillow)
"""
import sys
from PIL import Image
src, out, width = sys.argv[1], sys.argv[2], int(sys.argv[3])
im = Image.open(src).convert("RGBA")
a = im.getchannel("A").point(lambda v: 255 if v > 8 else 0)
x0, y0, x1, y1 = a.getbbox()
pad = int(0.02 * max(x1 - x0, y1 - y0))
im = im.crop((max(0, x0 - pad), max(0, y0 - pad), min(im.width, x1 + pad), min(im.height, y1 + pad)))
if im.width > width:
    im = im.resize((width, round(im.height * width / im.width)), Image.LANCZOS)
im.save(out, optimize=True)
print(out, im.size)
