"""Figures of the comment.

  Fig. 1  the first three components follow two drug contrasts and report length
          a  LC1 receptor loadings vs DMT - MDMA (standardised affinities)
          b  LC2 drug means vs sentences per report
          c  LC3 receptor loadings vs LSD - DMT
  Fig. 2  drug-level t-u scores for LC1 and LC3, with and without the two dominant drugs

Run:  python figures.py      (writes results/figure1.svg, results/figure2.svg, and PNGs if cairosvg is installed)
"""
from pathlib import Path

import numpy as np

from common import load, standardise, r
from svg import SVG, INK, INK2, MUTED, RULE, ACCENT, GREY_DRUG, DRUG_COLOR

OUT = Path(__file__).parent / "results"
REDUCED = "#c0392b"
NICE = {"Imidazoline1": "I1", "Alpha2A": "α2A", "Alpha2B": "α2B", "Alpha2C": "α2C"}
PW, PH, PY = 330, 300, 120          # plot width, height, top

d = load()
n = d.n
Z = standardise(d.npki.to_numpy(), n)
z = dict(zip(d.drugs, Z))
marker_size = 3 + 9 * np.sqrt(n / n.max())
drug_colours = [DRUG_COLOR.get(x, GREY_DRUG) for x in d.drugs]


def fmt(v):
    return f"{v:g}".replace("-", "−")


def ticks(lo, hi, step):
    start = np.ceil(lo / step) * step
    return [round(v, 6) for v in np.arange(start, hi + step * 1e-6, step)]


class Panel:
    """A framed scatter panel with gridlines, ticks and axis labels."""

    def __init__(self, s, x0, letter, title, subtitle, xs, ys, xlabel, ylabel, xstep, ystep, pad=0.08):
        self.s, self.px, self.py = s, x0 + 62, PY
        s.label(x0, 40, letter, title)
        s.text(x0 + 26, 62, subtitle, 12, INK2)
        dx, dy = (xs.max() - xs.min()) * pad, (ys.max() - ys.min()) * pad
        self.dx = dx
        self.xlim = xs.min() - dx, xs.max() + dx
        self.ylim = ys.min() - dy, ys.max() + dy
        for t in ticks(*self.xlim, xstep):
            s.line(self.X(t), PY, self.X(t), PY + PH, RULE, 0.8, opacity=0.6)
            s.text(self.X(t), PY + PH + 18, fmt(t), 11, INK2, anchor="middle")
        for t in ticks(*self.ylim, ystep):
            s.line(self.px, self.Y(t), self.px + PW, self.Y(t), RULE, 0.8, opacity=0.6)
            s.text(self.px - 8, self.Y(t) + 4, fmt(t), 11, INK2, anchor="end")
        s.rect(self.px, PY, PW, PH, "none", stroke=RULE)
        s.text(self.px + PW / 2, PY + PH + 40, xlabel, 11.5, INK2, anchor="middle")
        s.vtext(self.px - 40, PY + PH / 2, ylabel)

    def X(self, v):
        lo, hi = self.xlim
        return self.px + PW * (v - lo) / (hi - lo)

    def Y(self, v):
        lo, hi = self.ylim
        return self.py + PH - PH * (v - lo) / (hi - lo)

    def fit(self, xs, ys, colour=MUTED, width=1.0):
        b1, b0 = np.polyfit(xs, ys, 1)
        a, b = xs.min() - self.dx / 2, xs.max() + self.dx / 2
        self.s.line(self.X(a), self.Y(b0 + b1 * a), self.X(b), self.Y(b0 + b1 * b), colour, width, dash="4 3")

    def points(self, xs, ys, sizes=None, colours=None):
        sizes = np.full(len(xs), 5.0) if sizes is None else sizes
        for i in np.argsort(-sizes):                     # large markers first
            colour = ACCENT if colours is None else colours[i]
            self.s.circle(self.X(xs[i]), self.Y(ys[i]), sizes[i], colour, 0.85, stroke="#ffffff")

    def labels(self, xs, ys, names, offsets):
        for name, (dx, dy, anchor) in offsets.items():
            i = names.index(name)
            self.s.text(self.X(xs[i]) + dx, self.Y(ys[i]) + dy, NICE.get(name, name), 10.5, INK2, anchor=anchor)

    def note(self, text, dy=66):
        self.s.text(self.px, PY + PH + dy, text, 10.5, MUTED)


# ---------------------------------------------------------------- Figure 1

def figure1():
    s = SVG(1400, 520)
    contrast = z["DMT"] - z["MDMA"]
    q = d.fig5.LC1.to_numpy()
    p = Panel(s, 20, "a", "LC1 follows DMT − MDMA", "receptor loadings vs difference of affinity profiles",
              contrast, q, "DMT − MDMA (standardised npKi)", "LC1 receptor loading", 1, 0.1, pad=0.14)
    p.fit(contrast, q)
    p.points(contrast, q)
    p.labels(contrast, q, d.receptors, {"5-HT7a": (-2, -10, "start"), "5-HT1E": (-6, -9, "end"), "M3": (8, 4, "start")})
    s.text(p.px + 10, PY + 20, f"r = {r(contrast, q):.2f}", 13, INK, weight="700")
    p.note("points: 19 receptors")

    length, lc2 = d.table1.sentences_per_report.to_numpy(), d.fig4.LC2.to_numpy()
    p = Panel(s, 480, "b", "LC2 follows report length", "drug means vs mean sentences per report",
              length, lc2, "sentences per report", "LC2 drug mean", 10, 2)
    p.fit(length, lc2)
    p.points(length, lc2, marker_size, drug_colours)
    p.labels(length, lc2, d.drugs, {"MDMA": (12, 4, "start"), "DMT": (12, 4, "start"), "LSD": (-12, -6, "end"),
                                    "psilocin": (0, -16, "middle"), "MEM": (9, 15, "start"), "DOET": (9, -7, "start")})
    s.text(p.px + 10, PY + 20, f"r = {r(length, lc2):.2f}", 13, INK, weight="700")
    p.note("points: 21 drugs; marker area proportional to number of reports")

    contrast = z["LSD"] - z["DMT"]
    q = d.fig5.LC3.to_numpy()
    p = Panel(s, 940, "c", "LC3 follows LSD − DMT", "receptor loadings vs difference of affinity profiles",
              contrast, q, "LSD − DMT (standardised npKi)", "LC3 receptor loading", 1, 0.05, pad=0.14)
    p.fit(contrast, q)
    p.points(contrast, q)
    p.labels(contrast, q, d.receptors, {"5-HT1B": (-8, -6, "end")})
    s.text(p.px + 10, PY + 20, f"r = {r(contrast, q):.2f}", 13, INK, weight="700")
    p.note("points: 19 receptors")
    return s


# ---------------------------------------------------------------- Figure 2

def figure2():
    s = SVG(940, 560)
    panels = [("a", "LC1", ("DMT", "MDMA"), 1,
               {"DMT": (-14, 4, "end"), "MDMA": (14, 4, "start"), "LSD": (-12, -10, "end"), "psilocin": (8, 20, "start")}),
              ("b", "LC3", ("LSD", "DMT"), 0.5,
               {"DMT": (14, 4, "start"), "MDMA": (14, 4, "start"), "LSD": (-14, 4, "end"), "psilocin": (12, 4, "start")})]
    for k, (letter, comp, pair, xstep, offsets) in enumerate(panels):
        u = Z @ d.fig5[comp].to_numpy()        # drug affinity profiles projected on the loadings
        t = d.fig4[comp].to_numpy()            # drug mean scores
        keep = ~np.isin(d.drugs, pair)
        p = Panel(s, 20 + 460 * k, letter, f"{comp} at the drug level",
                  "t–u scores: drug means vs affinity profiles projected on the loadings",
                  u, t, "u (affinity profile · receptor loadings)", f"t ({comp} drug mean)", xstep, 0.5, pad=0.1)
        p.fit(u, t, MUTED, 1.2)
        p.fit(u[keep], t[keep], REDUCED, 1.6)
        p.points(u, t, marker_size, drug_colours)
        p.labels(u, t, d.drugs, offsets)
        ly = PY + PH + 66
        s.line(p.px, ly - 4, p.px + 22, ly - 4, MUTED, 1.2, dash="4 3")
        s.text(p.px + 28, ly, f"all 21 drugs: r = {r(u, t):.2f}", 12, INK, weight="700")
        s.line(p.px, ly + 16, p.px + 22, ly + 16, REDUCED, 1.6, dash="4 3")
        s.text(p.px + 28, ly + 20, f"without {pair[0]} and {pair[1]}: r = {r(u[keep], t[keep]):.2f}",
               12, REDUCED, weight="700")
        p.note("points: 21 drugs; marker area proportional to number of reports", dy=110)
    return s


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    for name, make in [("figure1", figure1), ("figure2", figure2)]:
        make().save(OUT / f"{name}.svg")
        try:
            import cairosvg
            cairosvg.svg2png(url=str(OUT / f"{name}.svg"), write_to=str(OUT / f"{name}.png"), output_width=2800)
        except ImportError:
            pass
        print(f"results/{name}.svg")
