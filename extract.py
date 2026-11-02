"""Extract the published values used in the comment into data/*.csv.

Sources (put them in sources/ under the names below; see README.md, "Rerun from scratch"):
  the article PDF        Table 1 (report counts), Fig. 1 (cannabis co-mentions),
                         Fig. 4 (drug means), Fig. 5 (receptor loadings)
  ESM1 (.xlsx)           receptor affinities (Ki, npKi)
  ESM2 (.xlsx)           Table S2: word lists of each component end (summarised only, see below)
  ESM4 (.docx)           Fig. S2 (cannabis-negative drug means), Fig. S3 (cannabis-negative loadings)
  en_full.txt            OpenSubtitles 2018 English word counts (downloaded if missing)

The word lists themselves and their word frequencies are written to local/, which is not published;
data/ receives only per-component summary statistics.

Run:  python extract.py
"""
import re
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import pymupdf
from scipy.stats import mannwhitneyu

import digitise as dg

SOURCES = Path(__file__).parent / "sources"
DATA = Path(__file__).parent / "data"
LOCAL = Path(__file__).parent / "local"
PAPER = SOURCES / "Suresh_2026_Psychopharmacology.pdf"
ESM1 = SOURCES / "213_2026_7170_MOESM1_ESM.xlsx"
ESM2 = SOURCES / "213_2026_7170_MOESM2_ESM.xlsx"
ESM4 = SOURCES / "213_2026_7170_MOESM4_ESM.docx"
FREQUENCIES = SOURCES / "en_full.txt"
FREQUENCIES_URL = "https://raw.githubusercontent.com/hermitdave/FrequencyWords/master/content/2018/en/en_full.txt"

ARTICLE = ("Suresh P, Kins J, Levkovskyy L, et al. (2026) Moving beyond 5-HT2A: language models map "
           "neurotransmitter receptors to experiential semantics. Psychopharmacology. "
           "https://doi.org/10.1007/s00213-026-07170-0")
RAY = ("Ray TS (2010) Psychedelics and the human receptorome. PLoS ONE 5:e9019. "
       "https://doi.org/10.1371/journal.pone.0009019 (CC BY)")
OPENSUBTITLES = ("Dave H (2018) FrequencyWords: OpenSubtitles 2018 English word counts. "
                 "https://github.com/hermitdave/FrequencyWords (CC BY-SA 4.0)")


def save(df, name, source, **kwargs):
    """Write a CSV to data/ with comment lines naming its source (read back with comment='#')."""
    with open(DATA / name, "w") as f:
        f.write(f"# Source: {source}\n# Values extracted by extract.py; see README.md.\n")
        df.to_csv(f, **kwargs)


# Drug names as in Table 1 (the canonical spelling used throughout this repository).
DRUGS = ["MEM", "DOET", "DiPT", "DOM", "2C-B-fly", "5-MeO-MiPT", "DOI", "TMA-2", "MDA", "DOB",
         "5-MeO-DiPT", "5-MeO-DMT", "2C-B", "DPT", "2C-T-2", "mescaline", "DMT", "2C-E", "LSD",
         "MDMA", "psilocin"]
BY_LOWER = {d.lower(): d for d in DRUGS}


# ---------------------------------------------------------------- Table 1

def table1():
    """Rows of Table 1 end in: n reports, n sentences, sentences/report."""
    page = pymupdf.open(PAPER)[3]
    lines = {}
    for w in page.get_text("words"):
        lines.setdefault(round(w[3]), []).append(w)
    rows = []
    for y in sorted(lines):
        tokens = [w[4] for w in sorted(lines[y], key=lambda w: w[0])]
        if (len(tokens) >= 5 and re.fullmatch(r"[\d,]+", tokens[-3]) and re.fullmatch(r"[\d,]+", tokens[-2])
                and re.fullmatch(r"\d+\.\d", tokens[-1])):
            k = next(i for i, t in enumerate(tokens) if t in ("PEA", "TRYPT", "ERG", "N/A"))
            name = "".join(tokens[:k]).rstrip("*")
            rows.append(dict(drug=name, drug_class=tokens[k],
                             n_reports=int(tokens[-3].replace(",", "")),
                             n_sentences=int(tokens[-2].replace(",", "")),
                             sentences_per_report_printed=float(tokens[-1])))
    t = pd.DataFrame(rows).set_index("drug")
    total = t.loc["total"]
    t = t.drop("total")
    assert sorted(t.index) == sorted(DRUGS), t.index
    assert t.n_reports.sum() == total.n_reports and t.n_sentences.sum() == total.n_sentences
    t["sentences_per_report"] = t.n_sentences / t.n_reports
    return t.loc[DRUGS]


# ---------------------------------------------------------------- Fig. 1: cannabis co-mentions

FIG1_COLUMNS = ["Overall", "MDMA", "5-MeO-DiPT", "psilocin", "2C-B", "2C-E", "2C-T-2", "2C-B-fly", "LSD",
                "TMA-2", "mescaline", "MDA", "DMT", "DOI", "DOB", "5-MeO-MiPT", "DOM", "DPT", "5-MeO-DMT",
                "DiPT", "MEM", "DOET"]
FIG1_ROWS = ["Any tracked", "Cannabis", "Alcohol", "Nicotine", "Caffeine", "Medication terms"]


def fig1(t1):
    img = dg.pdf_image(PAPER, page=7)
    text, _ = dg.heatmap_text(img, box=(300, 18, 1860, 677), nrow=6, ncol=22)
    pct = pd.DataFrame(dg.parse_percent(text), index=FIG1_ROWS, columns=FIG1_COLUMNS)
    assert not pct.isna().any().any(), pct
    out = pct.T.loc[DRUGS, ["Cannabis"]].rename(columns={"Cannabis": "cannabis_pct"})
    out["n_reports"] = t1.n_reports
    out["n_cannabis_negative"] = (out.n_reports * (1 - out.cannabis_pct / 100)).round().astype(int)
    # the article reports 3,368 cannabis-negative reports
    assert out.n_cannabis_negative.sum() == 3368, out.n_cannabis_negative.sum()
    return out


# ---------------------------------------------------------------- receptor loadings (Fig. 5, Fig. S3)

# Row (component) and column (receptor) order of the clustered heatmaps, transcribed from their labels.
FIG5_COMPONENTS = [9, 1, 8, 4, 16, 7, 18, 11, 3, 6, 13, 14, 12, 17, 10, 19, 5, 2, 15]
FIG5_RECEPTORS = ["Imidazoline1", "Alpha2B", "Alpha2C", "Alpha2A", "M3", "5-HT1A", "D2", "5-HT1B", "D3",
                  "SERT", "H1", "Sigma1", "5-HT2A", "5-HT2C", "5-HT1D", "5-HT7a", "5-HT5A", "5-HT1E", "5-HT6"]
FIGS3_COMPONENTS = [8, 17, 7, 19, 6, 13, 1, 18, 5, 14, 15, 3, 9, 2, 16, 12, 4, 10, 11]
FIGS3_RECEPTORS = ["Alpha2A", "M3", "H1", "SERT", "Sigma1", "Imidazoline1", "Alpha2B", "Alpha2C", "5-HT1A",
                   "D2", "5-HT1B", "D3", "5-HT2A", "5-HT1E", "5-HT2C", "5-HT1D", "5-HT7a", "5-HT5A", "5-HT6"]


def loadings(img, box, bar_box, components, receptors):
    text, negative = dg.heatmap_text(img, box, nrow=len(components), ncol=len(receptors))
    values = dg.parse_loadings(text, negative)
    assert not np.isnan(values).any(), "unread cells"
    worst = dg.colour_check(img, box, values, bar_box)
    print(f"  largest difference between OCR value and cell colour: {worst:.3f}")
    assert worst < 0.03, "an OCR value disagrees with its cell colour"
    out = pd.DataFrame(values.T, index=receptors, columns=[f"LC{c}" for c in components])
    return out[[f"LC{c}" for c in range(1, len(components) + 1)]].sort_index()


# ---------------------------------------------------------------- drug means (Fig. 4, Fig. S2)

# Approximate marker positions read by eye from the labelled markers; used only to name markers.
FIG4_APPROX = {"MEM": (0.6, -7.3), "DOET": (0.1, -7.1), "DOM": (-0.2, -4.5), "DiPT": (-0.1, -4.5),
               "2C-B-fly": (-0.5, -4.3), "5-MeO-MiPT": (-0.5, -3.9), "MDA": (-0.8, -3.9), "TMA-2": (-0.2, -3.4),
               "DOI": (-0.3, -3.2), "5-MeO-DiPT": (-0.8, -2.5), "DOB": (-0.3, -2.4), "2C-T-2": (-0.6, -1.9),
               "5-MeO-DMT": (0.3, -1.2), "mescaline": (0.1, -1.1), "2C-B": (0.1, -0.7), "DPT": (0.5, -0.6),
               "MDMA": (-2.1, -0.1), "DMT": (1.6, 1.6), "2C-E": (0.2, 1.7), "LSD": (0.6, 2.3), "psilocin": (0.3, 2.4)}
FIGS2_APPROX = {"MEM": (-5.2, 1.0), "DOET": (-5.4, 0.0), "DOM": (-3.8, -0.5), "DiPT": (-3.9, 0.3),
                "2C-B-fly": (-4.7, -0.1), "5-MeO-MiPT": (-3.3, -0.8), "MDA": (-2.8, -1.3), "TMA-2": (-2.8, 0.0),
                "DOI": (-2.6, -0.2), "5-MeO-DiPT": (-1.9, -1.5), "DOB": (-2.2, -0.7), "2C-T-2": (-2.1, -1.4),
                "5-MeO-DMT": (-0.8, 0.0), "mescaline": (-1.2, 0.6), "2C-B": (-0.4, -0.7), "DPT": (-0.2, 0.0),
                "MDMA": (-0.2, -2.5), "DMT": (1.3, 2.1), "2C-E": (1.3, -0.1), "LSD": (2.5, 1.4), "psilocin": (2.1, 0.9)}


def drug_means(img, plot_box, x_grid, y_grid, approx, bar_box, ticks):
    """x_grid, y_grid: (pixel, value) pairs of labelled gridlines, measured on the image."""
    to_x = dg.linear_map(*zip(*x_grid))
    to_y = dg.linear_map(*zip(*y_grid))
    found = dg.markers(img, plot_box, approx, to_x, to_y)
    colours, values = dg.colourbar(img, bar_box, ticks)
    rows = {d: dict(LC1=m["x"], LC2=m["y"], LC3=dg.colour_to_value(m["colour"], colours, values))
            for d, m in found.items()}
    return pd.DataFrame(rows).T.loc[DRUGS].round(2)


# ---------------------------------------------------------------- word lists (ESM2, Table S2)

def word_lists():
    """One row per component end: its LLM label and ordered word list.

    The sheet has the positive ends in columns 0-3 and the negative ends in columns 4-7
    (end, LLM label, expert label, words)."""
    sheet = pd.read_excel(ESM2, header=None).iloc[1:]
    rows = []
    for _, row in sheet.iterrows():
        for c in (0, 4):
            end = row.iloc[c]
            if isinstance(end, str):                                  # e.g. "L2+"
                words = [w.strip().lower() for w in row.iloc[c + 3].split(";") if w.strip()]
                rows.append(dict(component=int(end[1:-1]), end=end[-1], llm_label=str(row.iloc[c + 1]).strip(),
                                 words=";".join(words)))
    out = pd.DataFrame(rows).sort_values(["component", "end"])
    assert len(out) == 38
    return out


def word_frequencies(lists):
    """General-English frequency of every listed word, in log10 occurrences per billion words
    (missing if the word is not in the frequency list)."""
    if not FREQUENCIES.exists():
        print(f"downloading {FREQUENCIES_URL}")
        urllib.request.urlretrieve(FREQUENCIES_URL, FREQUENCIES)
    counts = pd.read_csv(FREQUENCIES, sep=" ", header=None, names=["word", "count"],
                         keep_default_na=False, quoting=3)
    total = counts["count"].sum()
    counts = counts.set_index("word")["count"]
    words = sorted({w for ws in lists.words for w in ws.split(";")})
    out = pd.DataFrame({"word": words})
    out["count"] = [int(counts.get(w, 0)) for w in words]
    out["log10_per_billion"] = np.where(out["count"] > 0, np.log10(out["count"].clip(lower=1) / total * 1e9), np.nan)
    return out.set_index("word"), int(total)


def word_list_summary(lists, freq):
    """Per component: list lengths, median word frequency of each end, the AUC, and Erowid page metadata.

    AUC = P(a word of the + list is more frequent than a word of the - list), i.e. the Mann-Whitney U
    statistic divided by the product of the list lengths. Words absent from the frequency list are
    ranked below all others."""
    rows = []
    for comp, group in lists.groupby("component"):
        row = dict(component=comp)
        logf = {}
        for end, name in (("+", "plus"), ("-", "minus")):
            entry = group[group.end == end].iloc[0]
            words = entry.words.split(";")
            logf[end] = freq.log10_per_billion[words].fillna(-np.inf).to_numpy()
            row[f"llm_label_{name}"] = entry.llm_label
            row[f"n_words_{name}"] = len(words)
            row[f"median_log10_freq_{name}"] = round(float(np.median(logf[end])), 3)
            row[f"n_not_in_frequency_list_{name}"] = int(np.isinf(logf[end]).sum())
            row[f"page_metadata_{name}"] = any(re.fullmatch(r"\d{4}(expid|views)", w) for w in words)
        a, b = logf["+"], logf["-"]
        row["auc"] = round(float(mannwhitneyu(a, b).statistic / (len(a) * len(b))), 4)
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- run

def main():
    DATA.mkdir(exist_ok=True)

    t1 = table1()
    save(t1, "table1.csv", f"Table 1 of {ARTICLE}", index_label="drug")
    print(f"Table 1: {len(t1)} drugs, {t1.n_reports.sum()} reports")

    aff = pd.read_excel(ESM1)
    aff["drug"] = aff["drug"].map(BY_LOWER)
    assert aff.drug.notna().all()
    save(aff, "affinities_esm1.csv", f"ESM1 of {ARTICLE}; values from {RAY}", index=False)
    print(f"ESM1: {len(aff)} drug-receptor pairs")

    c = fig1(t1)
    save(c, "fig1_cannabis.csv", f"Fig. 1 (cannabis row, digitised) and Table 1 of {ARTICLE}", index_label="drug")
    print(f"Fig. 1: {c.n_cannabis_negative.sum()} cannabis-negative reports")

    fig5 = loadings(dg.pdf_image(PAPER, page=11), (132, 128, 1373, 2214), (16, 169, 50, 448),
                    FIG5_COMPONENTS, FIG5_RECEPTORS)
    save(fig5, "fig5_loadings.csv", f"Fig. 5 (digitised) of {ARTICLE}", index_label="receptor")
    figS3 = loadings(dg.docx_image(ESM4, "image3.png"), (112, 106, 1042, 1641), (22, 134, 47, 340),
                     FIGS3_COMPONENTS, FIGS3_RECEPTORS)
    save(figS3, "figS3_loadings_cannabis_negative.csv", f"ESM4, Fig. S3 (digitised) of {ARTICLE}",
         index_label="receptor")
    print("Fig. 5 and Fig. S3: 19 x 19 loadings each")

    fig4 = drug_means(dg.pdf_image(PAPER, page=10), plot_box=(140, 15, 2025, 920),
                      x_grid=[(209.5, -2), (454, -1.5), (700, -1), (945, -.5), (1191, 0), (1437.5, .5),
                              (1682.5, 1), (1928, 1.5)],
                      y_grid=[(67, 2), (247, 0), (426, -2), (604.5, -4), (784.5, -6)],
                      approx=FIG4_APPROX, bar_box=(218, 1108, 1113, 1125), ticks=[-1, -.5, 0, .5, 1, 1.5])
    save(fig4, "fig4_drug_means.csv", f"Fig. 4 (digitised) of {ARTICLE}", index_label="drug")
    figS2 = drug_means(dg.docx_image(ESM4, "image2.png"), plot_box=(190, 16, 2215, 1005),
                       x_grid=[(333, -5), (578, -4), (821, -3), (1067, -2), (1311, -1), (1556.5, 0),
                               (1802, 1), (2046, 2)],
                       y_grid=[(54.5, 2), (158, 1.5), (262, 1), (367, .5), (470, 0), (574, -.5), (679, -1),
                               (782, -1.5), (889, -2), (993, -2.5)],
                       approx=FIGS2_APPROX, bar_box=(247, 1218, 1230, 1234), ticks=[-1, -.5, 0, .5, 1, 1.5])
    save(figS2, "figS2_drug_means_cannabis_negative.csv", f"ESM4, Fig. S2 (digitised) of {ARTICLE}",
         index_label="drug")
    print("Fig. 4 and Fig. S2: 21 drug means each")

    LOCAL.mkdir(exist_ok=True)
    lists = word_lists()
    lists.to_csv(LOCAL / "esm2_word_lists.csv", index=False)               # not published
    freq, total = word_frequencies(lists)
    freq.to_csv(LOCAL / "word_frequencies.csv")                             # not published
    summary = word_list_summary(lists, freq)
    save(summary, "word_list_summary.csv",
         f"ESM2, Table S2 of {ARTICLE}; word frequencies from {OPENSUBTITLES}", index=False)
    print(f"ESM2: {len(lists)} word lists summarised; {freq['count'].eq(0).sum()} of {len(freq)} distinct words "
          f"not in the frequency list ({total:,} words)")

if __name__ == "__main__":
    main()
