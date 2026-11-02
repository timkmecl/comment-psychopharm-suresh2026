"""Read numbers out of the published figures of Suresh et al. (2026).

Three kinds of figure are read:
  annotated heatmaps (Fig. 1, Fig. 5, Fig. S3): the number printed in each cell, by OCR
  scatter plots (Fig. 4, Fig. S2): marker centres from pixel positions, LC3 from marker colour
  colourbars (Fig. 4, Fig. S2): colour along the bar, calibrated by its tick labels

Pixel boxes and gridline positions for each figure are given in extract.py.
"""
import io
import re
import zipfile

import numpy as np
import pymupdf
import pytesseract
from PIL import Image
from scipy import ndimage as ndi
from scipy.cluster.vq import kmeans2


# ---------------------------------------------------------------- images from the source files

def pdf_image(pdf_path, page, index=0):
    """The index-th embedded image on a (1-based) page of a PDF, as an RGB float array."""
    doc = pymupdf.open(pdf_path)
    xref = doc[page - 1].get_images()[index][0]
    return as_array(doc.extract_image(xref)["image"])


def docx_image(docx_path, name):
    """An image stored in a .docx file, e.g. name='image2.png'."""
    with zipfile.ZipFile(docx_path) as z:
        return as_array(z.read(f"word/media/{name}"))


def as_array(data):
    return np.asarray(Image.open(io.BytesIO(data)).convert("RGB")).astype(float)


# ---------------------------------------------------------------- annotated heatmaps

def heatmap_text(img, box, nrow, ncol, scale=4):
    """OCR the number printed in each cell of a heatmap.

    box = (left, top, right, bottom) pixel edges of the cell grid. Each cell is turned into
    dark-on-white "ink" relative to its own background colour, so light text on dark cells
    reads the same as dark text on light cells. All cells are pasted onto one canvas and read
    in a single OCR call. Minus signs are small and OCR drops them, so the sign is read from
    the glyph shape instead (see is_negative).

    Returns an array of OCR strings (digits and '.') and a boolean array of minus signs.
    """
    x0, y0, x1, y1 = box
    cw, ch = (x1 - x0) / ncol, (y1 - y0) / nrow
    inks, negative = {}, np.zeros((nrow, ncol), bool)
    for i in range(nrow):
        for j in range(ncol):
            cell = img[int(y0 + (i + .25) * ch):int(y0 + (i + .75) * ch),
                       int(x0 + (j + .08) * cw):int(x0 + (j + .92) * cw)]
            inks[i, j] = ink(cell)
            negative[i, j] = is_negative(inks[i, j] < 128)

    pad = 30
    h = max(m.shape[0] for m in inks.values()) + 2 * pad
    w = max(m.shape[1] for m in inks.values()) + 2 * pad
    canvas = np.full((nrow * h, ncol * w), 255, np.uint8)
    for (i, j), m in inks.items():
        canvas[i * h + pad:i * h + pad + m.shape[0], j * w + pad:j * w + pad + m.shape[1]] = m
    big = Image.fromarray(canvas).resize((canvas.shape[1] * scale, canvas.shape[0] * scale), Image.BICUBIC)
    words = pytesseract.image_to_data(big, config="--psm 11 -c tessedit_char_whitelist=0123456789.",
                                      output_type=pytesseract.Output.DICT, timeout=120)

    text = np.full((nrow, ncol), "", object)
    for t, left, top, ww, hh in zip(words["text"], words["left"], words["top"], words["width"], words["height"]):
        if t.strip():
            i = int((top + hh / 2) / scale // h)
            j = int((left + ww / 2) / scale // w)
            text[i, j] += t.strip()
    return text, negative


def ink(cell):
    """Greyscale image: 255 where the cell has its background colour, darker where it has text."""
    edge = np.concatenate([cell[0], cell[-1], cell[:, 0], cell[:, -1]])
    background = np.median(edge, 0)
    distance = np.sqrt(((cell - background) ** 2).sum(2))
    return (255 - np.clip(distance / distance.max() * 255 * 1.3, 0, 255)).astype(np.uint8)


def is_negative(mask):
    """True if the leftmost glyph is a minus: its ink lies only in the middle of the text height."""
    rows, cols = np.where(mask)
    if len(rows) == 0:
        return False
    top, bottom, left = rows.min(), rows.max(), cols.min()
    height = bottom - top + 1
    lead = np.where(mask[:, left:left + 3].any(1))[0]
    return lead.min() > top + 0.25 * height and lead.max() < bottom - 0.25 * height


def parse_loadings(text, negative):
    """Heatmap cells printed as -0.dd / 0.dd. All loadings are below 1 in magnitude, so only the
    two decimals are taken from OCR (the leading 0 is sometimes misread)."""
    out = np.full(text.shape, np.nan)
    for idx, t in np.ndenumerate(text):
        m = re.search(r"\.(\d\d)$", t)
        if m:
            out[idx] = (-1 if negative[idx] else 1) * int(m.group(1)) / 100
    return out


def colour_check(img, box, values, bar_box):
    """Independent check of OCR'd heatmap values against cell colours.

    Each cell's background colour is located on the vertical colourbar (bar_box); position on the
    bar is linear in value, so regressing the OCR values on bar position should leave only small
    residuals. Returns the largest absolute residual.
    """
    nrow, ncol = values.shape
    bx0, by0, bx1, by1 = bar_box
    xc = (bx0 + bx1) // 2
    bar = img[by0 + 2:by1 - 1, xc - 3:xc + 4].mean(1)
    x0, y0, x1, y1 = box
    cw, ch = (x1 - x0) / ncol, (y1 - y0) / nrow
    position = np.zeros(values.shape)
    for i in range(nrow):
        for j in range(ncol):
            strip = img[int(y0 + (i + .05) * ch):int(y0 + (i + .2) * ch), int(x0 + (j + .1) * cw):int(x0 + (j + .9) * cw)]
            colour = np.median(strip.reshape(-1, 3), 0)
            position[i, j] = np.argmin(((bar - colour) ** 2).sum(1))
    b, a = np.polyfit(position.ravel(), values.ravel(), 1)
    return float(np.abs(values - (a + b * position)).max())


def parse_percent(text):
    """Heatmap cells printed as dd.d (no negative values)."""
    out = np.full(text.shape, np.nan)
    for idx, t in np.ndenumerate(text):
        m = re.search(r"(\d{1,2}\.\d)$", t)
        if m:
            out[idx] = float(m.group(1))
    return out


# ---------------------------------------------------------------- scatter plots

def linear_map(pixels, values):
    """Pixel -> data value, fitted through gridline positions."""
    b, a = np.polyfit(pixels, values, 1)
    return lambda p: a + b * np.asarray(p, float)


def markers(img, plot_box, approx, to_x, to_y, size=13):
    """Centre (in data units) and colour of each filled marker.

    Markers are the saturated regions that survive an opening with a size x size square,
    which removes the coloured text labels. Each marker is named by the drug whose
    approximate position (approx: name -> (x, y), read by eye) lies nearest. Touching markers
    form one region; it is split by k-means started from the approximate positions.
    """
    x0, y0, x1, y1 = plot_box
    saturation = img.max(2) - img.min(2)
    mask = np.zeros(saturation.shape, bool)
    mask[y0:y1, x0:x1] = saturation[y0:y1, x0:x1] > 80
    mask = ndi.binary_opening(mask, np.ones((size, size)))
    inner = ndi.binary_erosion(mask, np.ones((5, 5)))          # avoids marker outlines
    labels, _ = ndi.label(mask)

    to_px = inverse(to_x, x0, x1)
    to_py = inverse(to_y, y0, y1)
    guess = {d: (to_px(x), to_py(y)) for d, (x, y) in approx.items()}
    ys, xs = np.where(labels > 0)
    members = {}
    for d, (px, py) in guess.items():
        k = np.argmin((xs - px) ** 2 + (ys - py) ** 2)
        members.setdefault(labels[ys[k], xs[k]], []).append(d)

    out = {}
    for region, drugs in members.items():
        ry, rx = np.where(labels == region)
        points = np.c_[rx, ry].astype(float)
        if len(drugs) == 1:
            groups = {drugs[0]: np.ones(len(points), bool)}
        else:
            _, assign = kmeans2(points, np.array([guess[d] for d in drugs]), minit="matrix", iter=50)
            groups = {d: assign == k for k, d in enumerate(drugs)}
        for d, sel in groups.items():
            cx, cy = points[sel].mean(0)
            core = sel & inner[ry, rx]
            out[d] = dict(x=float(to_x(cx)), y=float(to_y(cy)),
                          colour=np.median(img[ry[core], rx[core]], 0))
    return out


def inverse(f, lo, hi):
    """Data value -> pixel for a linear map f."""
    p = np.array([lo, hi], float)
    v = f(p)
    return lambda value: float(np.interp(value, *(v, p) if v[0] < v[1] else (v[::-1], p[::-1])))


# ---------------------------------------------------------------- colourbars

def colourbar(img, bar_box, tick_values):
    """Colour and value at each pixel column of a horizontal colourbar.

    The tick labels below the bar are found as groups of dark pixels; their centres, in order,
    are matched to tick_values to calibrate pixel -> value.
    """
    x0, y0, x1, y1 = bar_box
    colours = ndi.uniform_filter1d(img[(y0 + y1) // 2 - 2:(y0 + y1) // 2 + 3, x0:x1].mean(0), 9, axis=0)
    below = img[y1 + 3:y1 + 50, x0 - 40:x1 + 40].mean(2) < 128
    cols = np.where(below.any(0))[0]
    groups = np.split(cols, np.where(np.diff(cols) > 15)[0] + 1)
    centres = [x0 - 40 + (g[0] + g[-1]) / 2 for g in groups]
    centres = [c for c in centres if x0 - 10 <= c <= x1 + 10]
    assert len(centres) == len(tick_values), f"found {len(centres)} tick labels, expected {len(tick_values)}"
    to_value = linear_map(centres, tick_values)
    return colours, to_value(np.arange(x0, x1))


def colour_to_value(colour, colours, values):
    return float(values[np.argmin(((colours - colour) ** 2).sum(1))])
