"""Minimal SVG writer and palette, shared with the methods-paper figures."""
INK, INK2, MUTED, RULE, ACCENT = "#1b2128", "#4d5966", "#8a96a3", "#d7dde3", "#1d5d7f"
SANS = "Helvetica, Arial, sans-serif"
SERIF = "Georgia, serif"
GREY_DRUG = "#a7b0ba"
DRUG_COLOR = {"LSD": "#0072B2", "psilocin": "#E69F00", "MDMA": "#CC79A7", "DMT": "#009E73"}


class SVG:
    def __init__(self, w, h):
        self.w, self.h, self.o = w, h, []
        self.o.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">')
        self.o.append('<defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" '
                      'markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" '
                      f'fill="{INK2}"/></marker></defs>')
        self.o.append(f'<rect width="{w}" height="{h}" fill="#ffffff"/>')

    def text(self, x, y, s, size=12, fill=INK, anchor="start", weight="400", style="normal", fam=SANS):
        self.o.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{fill}" text-anchor="{anchor}" '
                      f'font-weight="{weight}" font-style="{style}" font-family="{fam}">{s}</text>')

    def rect(self, x, y, w, h, fill, opacity=1, stroke="none", sw=1, rx=0):
        if h < 0:
            y, h = y + h, -h
        self.o.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{max(w,0):.1f}" height="{h:.1f}" rx="{rx}" fill="{fill}" '
                      f'fill-opacity="{opacity:.2f}" stroke="{stroke}" stroke-width="{sw}"/>')

    def line(self, x1, y1, x2, y2, stroke=INK2, sw=1, dash=None, arrow=False, opacity=1):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        m = ' marker-end="url(#ah)"' if arrow else ""
        self.o.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{stroke}" '
                      f'stroke-width="{sw}" stroke-opacity="{opacity}"{d}{m}/>')

    def path(self, d, stroke=INK2, sw=1, fill="none", opacity=1):
        self.o.append(f'<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" stroke-opacity="{opacity}"/>')

    def circle(self, x, y, r, fill, opacity=1, stroke="none"):
        self.o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{fill}" fill-opacity="{opacity}" stroke="{stroke}"/>')

    def vtext(self, x, y, s, size=11.5, fill=INK2):
        """Text rotated to read bottom-to-top, centred on (x, y)."""
        self.o.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{fill}" text-anchor="middle" '
                      f'font-family="{SANS}" transform="rotate(-90 {x:.1f} {y:.1f})">{s}</text>')

    def label(self, x, y, letter, title):
        self.text(x, y, letter, 20, weight="700")
        self.text(x + 26, y, title, 16, weight="700")

    def save(self, path):
        self.o.append("</svg>")
        with open(path, "w") as f:
            f.write("\n".join(self.o))
