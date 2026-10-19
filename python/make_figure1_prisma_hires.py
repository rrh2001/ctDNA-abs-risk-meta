# -*- coding: utf-8 -*-
"""Figure 1: PRISMA 2020 flow diagram, double-column 7.2in @600dpi (4320 px wide).
All geometry/fonts/strokes are scaled uniformly from the 2600-wide master (S=4320/2600)."""
import os
from PIL import Image, ImageDraw, ImageFont

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "figures")
os.makedirs(OUT, exist_ok=True)

DPI = 600
BASE_W = 2600.0
S = 4320.0 / BASE_W            # 1.6615 -> 7.2in wide @600dpi
def px(v):
    return int(round(v * S))

W, H = 4320, px(2010)
img = Image.new("RGB", (W, H), "white")
d = ImageDraw.Draw(img)

def font(sz, bold=False):
    candidates = [
        r"C:\Windows\Fonts\arialbd.ttf" if bold else r"C:\Windows\Fonts\arial.ttf",
        "/Library/Fonts/Arial Bold.ttf" if bold else "/Library/Fonts/Arial.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if bold else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    ]
    for c in candidates:
        if os.path.exists(c):
            return ImageFont.truetype(c, px(sz))
    return ImageFont.load_default()

F_BODY = font(40); F_BODYB = font(40, bold=True)
F_HEAD = font(44, bold=True)
F_SIDE = font(32); F_SIDEB = font(32, bold=True)
INK = (20, 20, 20); BLUE = (44, 62, 80)
GREY = (245, 246, 248); EDGE = (60, 60, 60)

def wrap(text, fnt, max_w):
    words, lines, cur = text.split(), [], ""
    for w in words:
        t = (cur + " " + w).strip()
        if d.textlength(t, font=fnt) <= max_w:
            cur = t
        else:
            if cur: lines.append(cur)
            cur = w
    if cur: lines.append(cur)
    return lines

def box(x, y, w, h, fill, edge=EDGE, lw=3):
    d.rectangle([px(x), px(y), px(x) + px(w), px(y) + px(h)], fill=fill, outline=edge, width=px(lw))

def para(x, y, max_w, segments, leading=10, gap=6):
    cy = px(y)
    for text, f, color in segments:
        for line in wrap(text, f, px(max_w)):
            d.text((px(x), cy), line, font=f, fill=color)
            cy += f.size + px(leading)
        cy += px(gap)
    return cy

def arrow_down(x, y1, y2, color=INK):
    x=px(x); y1=px(y1); y2=px(y2)
    d.line([x, y1, x, y2 - px(18)], fill=color, width=px(5))
    d.polygon([(x - px(12), y2 - px(18)), (x + px(12), y2 - px(18)), (x, y2)], fill=color)

def arrow_right(x1, x2, y, color=INK):
    x1=px(x1); x2=px(x2); y=px(y)
    d.line([x1, y, x2 - px(18), y], fill=color, width=px(5))
    d.polygon([(x2 - px(18), y - px(12)), (x2 - px(18), y + px(12)), (x2, y)], fill=color)

LX, LW = 110, 1480
RX, RW = 1720, 770
FULLW = LW + (RX - LX - LW) + RW
cw = 40

def phase(y, label):
    d.text((px(20), px(y)), label, font=F_HEAD, fill=BLUE)

phase(40, "IDENTIFICATION")
phase(396, "SCREENING")
phase(1166, "INCLUDED")

# Identification
b1_y, b1_h = 120, 250
box(LX, b1_y, LW, b1_h, "white")
para(LX + cw, b1_y + 26, LW - 2 * cw, [
    ("Records identified from databases:", F_BODYB, INK),
    ("PubMed/MEDLINE, Embase, Cochrane CENTRAL, Web of Science  (n = 1,963)", F_BODY, INK),
    ("Records identified through other sources (reference lists)  (n = 0)", F_BODY, INK),
])
box(RX, b1_y, RW, 250, GREY)
para(RX + cw, b1_y + 26, RW - 2 * cw, [
    ("Records removed before screening:", F_SIDEB, INK),
    ("Duplicate records removed  (n = 413)", F_SIDE, INK),
    ("Records marked ineligible by automation  (n = 0)", F_SIDE, INK),
])

# Screening
b2_y, b2_h = 470, 190
box(LX, b2_y, LW, b2_h, "white")
para(LX + cw, b2_y + 52, LW - 2 * cw, [
    ("Records screened by title and abstract  (n = 1,550)", F_BODY, INK),
])
box(RX, b2_y, RW, 190, GREY)
para(RX + cw, b2_y + 30, RW - 2 * cw, [
    ("Records excluded after title/abstract screening (reviews, case reports, meta-analyses, RCTs, single-arm or non-relevant studies)  (n = 1,464)", F_SIDE, INK),
])

b3_y, b3_h = 760, 340
box(LX, b3_y, LW, b3_h, "white")
para(LX + cw, b3_y + 110, LW - 2 * cw, [
    ("Full-text reports sought and assessed for eligibility  (n = 86)", F_BODY, INK),
])
box(RX, b3_y, RW, 340, GREY)
para(RX + cw, b3_y + 26, RW - 2 * cw, [
    ("Full-text reports excluded  (n = 62):", F_SIDEB, INK),
    ("\u2022 no usable Kaplan\u2013Meier curve or numbers-at-risk table  (n = 31)", F_SIDE, INK),
    ("\u2022 no ctDNA \u00d7 recurrence time-to-event comparison  (n = 16)", F_SIDE, INK),
    ("\u2022 overlapping or duplicate cohort  (n = 8)", F_SIDE, INK),
    ("\u2022 wrong setting or endpoint  (n = 7)", F_SIDE, INK),
], leading=8, gap=4)

# Included
b4_y, b4_h = 1240, 300
box(LX, b4_y, FULLW, b4_h, (235, 240, 245), edge=BLUE, lw=4)
para(LX + cw, b4_y + 34, FULLW - 2 * cw, [
    ("Studies included in the quantitative synthesis  (n = 24 independent studies)", F_BODYB, BLUE),
    ("190 reconstructed Kaplan\u2013Meier arms; 11,297 reconstructed patients; 3,041 recurrence events (1,765 ctDNA-positive; 1,276 ctDNA-negative).", F_BODY, INK),
])
b5_y, b5_h = 1620, 170
box(LX, b5_y, FULLW, b5_h, "white", edge=BLUE, lw=3)
para(LX + cw, b5_y + 48, FULLW - 2 * cw, [
    ("Of these, 8 studies also contributed serial / multi-state curves to the longitudinal trajectory analysis.", F_BODY, INK),
])

# arrows
arrow_down(LX + LW // 2, b1_y + b1_h, b2_y)
arrow_down(LX + LW // 2, b2_y + b2_h, b3_y)
arrow_down(LX + LW // 2, b3_y + b3_h, b4_y)
arrow_down(LX + LW // 2, b4_y + b4_h, b5_y)
arrow_right(LX + LW, RX, b2_y + b2_h // 2)
arrow_right(LX + LW, RX, b3_y + b3_h // 2)

png = os.path.join(OUT, "Figure_1_PRISMA_flow.png")
tif = os.path.join(OUT, "Figure_1_PRISMA_flow.tif")
img.save(png, dpi=(DPI, DPI))
img.save(tif, format="TIFF", dpi=(DPI, DPI), compression="tiff_lzw")
print("saved hires Figure 1", img.size)
