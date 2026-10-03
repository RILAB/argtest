"""Generate docs/graphical_abstract.svg. Usage: python docs/make_graphical_abstract.py docs/graphical_abstract.svg"""
import sys

W, H = 1600, 780
INK, MUTED, LABEL, LINE, PANEL = "#1A1D21", "#3F4650", "#5A626C", "#C9CED4", "#F3F5F7"
REC, ACC, LOAD, TEAL = "#2B5FA8", "#6B7280", "#C25A12", "#0E6B6A"
BASE, KEEP = "#E6E9ED", "#CFE5E3"
SANS = "'IBM Plex Sans', 'Helvetica Neue', Helvetica, Arial, sans-serif"
MONO = "'IBM Plex Mono', Menlo, Consolas, monospace"

N = 24
REC_COLS = [2, 3, 11, 12, 19]
ACC_COLS = [5, 15, 22]
LOAD_COLS = [7, 8, 9, 16, 17, 20]
PRUNED = {1: [20], 2: [7, 8, 9], 5: [16, 17]}

out = []
def e(s):
    out.append(s)

def text(x, y, s, size=14, weight=400, fill=MUTED, anchor="start", font=SANS, spacing=None):
    ls = f' letter-spacing="{spacing}"' if spacing else ""
    e(f'<text x="{x}" y="{y}" font-family="{font}" font-size="{size}" font-weight="{weight}" '
      f'fill="{fill}" text-anchor="{anchor}"{ls}>{s}</text>')

def rect(x, y, w, h, fill, rx=0, stroke=None, sw=1.5):
    st = f' stroke="{stroke}" stroke-width="{sw}"' if stroke else ""
    e(f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" rx="{rx}" fill="{fill}"{st}/>')

def badge(cx, cy, n, color):
    e(f'<circle cx="{cx}" cy="{cy}" r="17" fill="{color}"/>')
    text(cx, cy + 6, n, size=17, weight=700, fill="#FFFFFF", anchor="middle")

def cells(x0, y, width, h, colors, gap=3):
    cw = (width - gap * (len(colors) - 1)) / len(colors)
    for i, c in enumerate(colors):
        rect(x0 + i * (cw + gap), y, cw, h, c, rx=2)

def arrow(x, y):
    e(f'<path d="M{x} {y}h34M{x+24} {y-9}l10 9-10 9" fill="none" stroke="{INK}" '
      f'stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round"/>')

e(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" role="img" '
  f'aria-label="ARGtest pipeline overview: inferred ARGs are screened for low recombination, poor accessibility '
  f'and aberrant mutation load; flagged regions are trimmed from all samples and outlier individuals are pruned '
  f'per window; outputs are cleaned ARGs, validation plots and a summary report.">')
rect(0, 0, W, H, "#FFFFFF")

# Header
text(56, 92, "ARGtest", size=52, weight=700, fill=INK, spacing=-1)
text(56, 130, "Masking, trimming and validating ancestral recombination graphs", size=22, fill=MUTED)
text(1544, 104, "Snakemake pipeline · local or SLURM", size=15, anchor="end", font=MONO)
text(1544, 128, "input: SINGER-style .trees / .tsz", size=15, anchor="end", font=MONO)
rect(56, 150, 1488, 2, INK)

# Inputs
text(56, 196, "INPUT", size=14, weight=600, fill=LABEL, spacing=2)
rect(56, 212, 270, 170, PANEL, rx=10)
trees = [
    "M27 4v10M27 14L12 26M27 14l15 12M12 26L5 42M12 26l7 16M42 26l-6 16M42 26l8 16",
    "M27 4v8M27 12L8 30M27 12l16 10M8 30L4 42M8 30l6 12M43 22L26 42M43 22l4 20",
    "M27 4v12M27 16L16 22M27 16l14 18M16 22L5 42M16 22l8 20M41 34l-5 8M41 34l8 8",
]
for i, d in enumerate(trees):
    e(f'<g transform="translate({74 + i * 60} 228)"><path d="{d}" fill="none" stroke="{INK}" '
      f'stroke-width="2" stroke-linecap="round" stroke-linejoin="round"/></g>')
text(74, 306, "Inferred ARGs", size=19, weight=600, fill=INK)
text(74, 332, "Tree sequences, one per", size=15)
text(74, 352, "chromosome × replicate", size=15)
for i, (t, s) in enumerate([("Recombination map", "HapMap format"),
                            ("Mutation-rate map", "embedded or scalar µ"),
                            ("Reference index", ".fai chromosome lengths")]):
    y = 396 + i * 80
    rect(56, y, 270, 66, PANEL, rx=10)
    text(74, y + 28, t, size=17, weight=600, fill=INK)
    text(74, y + 50, s, size=14)

arrow(338, 470)
arrow(1170, 470)

# Detect
PX, TX, TRX, TRW = 399, 430, 694, 464
text(382, 196, "DETECT", size=14, weight=600, fill=LABEL, spacing=2)
text(1158, 196, "genomic windows →", size=13, fill=LABEL, anchor="end", font=MONO)
rows = [("1", REC, "Low recombination", "from the recombination map", REC_COLS),
        ("2", ACC, "Poor accessibility", "windows mostly masked in input", ACC_COLS),
        ("3", LOAD, "Aberrant mutation load", "per individual, vs. simulated", LOAD_COLS)]
for i, (n, col, t, s, cols) in enumerate(rows):
    cy = 234 + i * 50
    badge(PX, cy, n, col)
    text(TX, cy - 3, t, size=18, weight=600, fill=INK)
    text(TX, cy + 16, s, size=14)
    cells(TRX, cy - 11, TRW, 22, [col if c in cols else BASE for c in range(N)])

# Remove
text(382, 388, "REMOVE", size=14, weight=600, fill=LABEL, spacing=2)
rect(464, 382, 694, 1, LINE)
badge(PX, 424, "4", INK)
text(TX, 430, "Trim regions", size=18, weight=600, fill=INK)
text(TX, 452, "Masks 1 + 2 cut from every", size=14)
text(TX, 470, "sample; coordinates preserved", size=14)
badge(PX, 540, "5", INK)
text(TX, 546, "Trim samples", size=18, weight=600, fill=INK)
text(TX, 568, "Outlier individuals pruned only", size=14)
text(TX, 586, "in the windows where they fail", size=14)

rect(TRX, 406, TRW, 222, PANEL, rx=10)
for r in range(8):
    colors = []
    for c in range(N):
        if c in REC_COLS or c in ACC_COLS:
            colors.append("#FFFFFF")
        elif c in PRUNED.get(r, []):
            colors.append(LOAD)
        else:
            colors.append(KEEP)
    cells(TRX + 14, 420 + r * 23, TRW - 28, 19, colors)
text(TRX + 14, 618, "rows = individuals", size=12, fill=LABEL, font=MONO)
text(TRX + TRW - 14, 618, "columns = windows", size=12, fill=LABEL, anchor="end", font=MONO)

# Legend
lx = TX
for label, fill, stroke in [("retained", KEEP, None), ("region removed", "#FFFFFF", LINE),
                            ("sample pruned", LOAD, None)]:
    rect(lx, 652, 14, 14, fill, rx=2, stroke=stroke, sw=1)
    text(lx + 22, 664, label, size=14)
    lx += 22 + len(label) * 8 + 26
text(TX, 692, "Optional: drop windows with too few samples left", size=14, fill=LABEL)

# Outputs
OX = 1214
text(OX, 196, "OUTPUT", size=14, weight=600, fill=LABEL, spacing=2)
rect(OX, 212, 330, 112, TEAL, rx=10)
text(OX + 18, 240, "Cleaned ARGs", size=19, weight=600, fill="#FFFFFF")
text(OX + 18, 264, "Per chromosome, plus one merged", size=14, fill="#FFFFFF")
text(OX + 18, 284, "genome-wide tree sequence per", size=14, fill="#FFFFFF")
text(OX + 18, 304, "replicate · optional VCF", size=14, fill="#FFFFFF")

rect(OX, 338, 330, 224, "#FFFFFF", rx=10, stroke=LINE)
text(OX + 18, 366, "Validation", size=18, weight=600, fill=INK)
e(f'<g transform="translate({OX + 18} 380)">'
  f'<path d="M6 4v60h76" fill="none" stroke="{INK}" stroke-width="1.5"/>'
  f'<path d="M8 62L80 8" stroke="#8A929C" stroke-width="1.2" stroke-dasharray="3 3"/>')
for cx, cy, f in [(20, 50, REC), (25, 52, REC), (30, 44, REC), (40, 33, "#6E93C8"), (46, 35, "#6E93C8"),
                  (50, 29, "#6E93C8"), (58, 20, "#B7C9E3"), (68, 16, "#B7C9E3")]:
    e(f'<circle cx="{cx}" cy="{cy}" r="3" fill="{f}"/>')
e('</g>')
e(f'<g transform="translate({OX + 158} 380)"><path d="M6 4v60h76" fill="none" stroke="{INK}" stroke-width="1.5"/>')
for x, y in [(11, 14), (23, 36), (35, 45), (47, 50), (59, 53), (71, 55)]:
    e(f'<rect x="{x}" y="{y}" width="9" height="{64 - y}" fill="{TEAL}"/>')
e('</g>')
text(OX + 61, 464, "expected vs", size=12, anchor="middle")
text(OX + 61, 479, "observed", size=12, anchor="middle")
text(OX + 201, 464, "site frequency", size=12, anchor="middle")
text(OX + 201, 479, "spectrum", size=12, anchor="middle")
text(OX + 18, 506, "Load, π, Tajima's D and SFS,", size=14)
text(OX + 18, 526, "compared with mutations", size=14)
text(OX + 18, 546, "simulated on the ARG", size=14)

rect(OX, 576, 330, 90, "#FFFFFF", rx=10, stroke=LINE)
text(OX + 18, 604, "Summary report", size=18, weight=600, fill=INK)
text(OX + 18, 628, "One HTML page: genome retained at", size=14)
text(OX + 18, 648, "each step, outlier counts, all plots", size=14)

# Footer
rect(56, 724, 1488, 1, LINE)
text(56, 752, "github.com/RILAB/argtest", size=14, font=MONO)
text(1544, 752, "Ross-Ibarra 2026 · doi:10.5281/zenodo.19698118", size=14, anchor="end", font=MONO)

e("</svg>")
open(sys.argv[1], "w").write("\n".join(out) + "\n")
