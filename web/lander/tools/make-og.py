"""Regenerate assets/img/og.png (1200x630) in the teal + amber palette."""
import math, sys
from PIL import Image, ImageDraw, ImageFont, ImageFilter

site, fonts, out = sys.argv[1], sys.argv[2], sys.argv[3]
W, H, SS = 1200, 630, 3            # supersample for crisp vector bits
BG = (10, 12, 15)
TEAL, BLUE, AMBER = (46, 196, 182), (92, 157, 247), (245, 165, 36)
TEXT, DIM, FAINT = (230, 236, 245), (151, 163, 179), (109, 120, 133)

inter_b = lambda s: ImageFont.truetype(f"{fonts}/inter/extras/ttf/Inter-Bold.ttf", s * SS)
inter_r = lambda s: ImageFont.truetype(f"{fonts}/inter/extras/ttf/Inter-Regular.ttf", s * SS)
mono = lambda s: ImageFont.truetype(f"{fonts}/jb/fonts/ttf/JetBrainsMono-Regular.ttf", s * SS)

img = Image.new("RGB", (W * SS, H * SS), BG)

def glow(cx, cy, rx, ry, color, alpha):
    layer = Image.new("L", (W, H), 0)
    d = ImageDraw.Draw(layer)
    d.ellipse([cx - rx, cy - ry, cx + rx, cy + ry], fill=int(255 * alpha))
    layer = layer.filter(ImageFilter.GaussianBlur(min(rx, ry) * 0.55)).resize((W * SS, H * SS))
    img.paste(Image.new("RGB", img.size, color), (0, 0), layer)

glow(930, 320, 330, 280, TEAL, 0.17)
# (no amber glow: teal + amber mix to olive)

# Torus cutout (the same render the hero uses).
torus = Image.open(f"{site}/assets/img/torus-boolean.png").convert("RGBA")
tw = 470
torus = torus.resize((tw * SS, int(torus.height * tw / torus.width) * SS), Image.LANCZOS)
img.paste(torus, ((W - tw - 40) * SS, (H - torus.height // SS) // 2 * SS - 20 * SS), torus)

d = ImageDraw.Draw(img)

# Logo: the brand mark, white outline with the amber arc.
def mark(x, y, size):
    k = size / 32 * SS
    P = lambda px, py: (x * SS + px * k, y * SS + py * k)
    pts = []
    for t in range(0, 91, 5):        # top-left rounded corner, centre (12,12) r 8
        a = math.radians(180 + t); pts.append(P(12 + 8 * math.cos(a), 12 + 8 * math.sin(a)))
    pts += [P(28, 4), P(28, 20)]
    for t in range(0, 91, 5):        # bottom-right rounded corner, centre (20,20) r 8
        a = math.radians(t); pts.append(P(20 + 8 * math.cos(a), 20 + 8 * math.sin(a)))
    pts += [P(4, 28), P(4, 12)]
    d.line(pts + [pts[0]], fill=TEXT, width=int(2.2 * k), joint="curve")
    arc = [P(12 + 8 * math.cos(math.radians(180 + t)), 12 + 8 * math.sin(math.radians(180 + t))) for t in range(0, 91, 3)]
    w = int(3.6 * k)
    d.line(arc, fill=AMBER, width=w, joint="curve")
    for e in (arc[0], arc[-1]):
        d.ellipse([e[0] - w / 2, e[1] - w / 2, e[0] + w / 2, e[1] + w / 2], fill=AMBER)

mark(80, 74, 46)
d.text((142 * SS, 97 * SS), "popkern", font=inter_b(34), fill=TEXT, anchor="lm")

# Headline: white first line, teal->blue gradient after.
def gradient_text(x, y, text, font):
    box = d.textbbox((x, y), text, font=font, anchor="ls")
    mask = Image.new("L", img.size, 0)
    ImageDraw.Draw(mask).text((x, y), text, font=font, fill=255, anchor="ls")
    grad = Image.new("RGB", img.size)
    gd = ImageDraw.Draw(grad)
    x0, x1 = box[0], box[0] + 640 * SS
    for gx in range(box[0], box[2] + 1):
        f = max(0.0, min(1.0, (gx - x0) / (x1 - x0)))
        gd.line([(gx, box[1]), (gx, box[3])], fill=tuple(int(TEAL[i] + (BLUE[i] - TEAL[i]) * f) for i in range(3)))
    img.paste(grad, (0, 0), mask)

hf = inter_b(58)
d.text((80 * SS, 238 * SS), "A CAD kernel for", font=hf, fill=TEXT, anchor="ls")
gradient_text(80 * SS, 306 * SS, "everyone who", hf)
gradient_text(80 * SS, 374 * SS, "makes things.", hf)

sf = inter_r(23)
d.text((80 * SS, 428 * SS), "Free and open source, written in Rust.", font=sf, fill=DIM, anchor="ls")
d.text((80 * SS, 462 * SS), "Booleans, fillets, sweeps, shells, STEP and IGES.", font=sf, fill=DIM, anchor="ls")

d.line([(80 * SS, 540 * SS), (1120 * SS, 540 * SS)], fill=(38, 44, 51), width=SS)
mf = mono(19)
d.text((80 * SS, 574 * SS), "AGPL-3.0  ·  Rhai scripts  ·  built by agents", font=mf, fill=FAINT, anchor="lm")
d.text((1120 * SS, 574 * SS), "github.com/graemedouglas/popkern", font=mf, fill=FAINT, anchor="rm")

img.resize((W, H), Image.LANCZOS).save(out, optimize=True)
print("wrote", out)
