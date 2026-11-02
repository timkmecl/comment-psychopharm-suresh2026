"""Recompute every number stated in the comment from the extracted data.

Run:  python check_numbers.py     (writes results/numbers.md)
"""
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from common import load, standardise, kish, loadings_from_drug_means, r

d = load()
n, t1 = d.n, d.table1
Z = standardise(d.npki.to_numpy(), n)                       # drugs x 19 receptors
z = dict(zip(d.drugs, Z))
rows = []


def claim(section, text, stated, value):
    rows.append((section, text, stated, value))


# ---------------------------------------------------------------- sample
claim("Opening", "reports", "4,496", f"{int(n.sum()):,}")
claim("Opening", "drugs", "21", f"{len(d.drugs)}")
claim("Opening", "receptors", "19", f"{len(d.receptors)}")

# ---------------------------------------------------------------- loadings follow from drug means
for comp in ["LC1", "LC2", "LC3"]:
    q = loadings_from_drug_means(d.fig4[comp], Z, n)
    stated = "r ≥ 0.998" if comp in ("LC1", "LC3") else "(not stated)"
    claim("Drug means", f"published {comp} loadings (Fig. 5) vs loadings implied by drug means (Fig. 4)",
          stated, f"r = {r(q, d.fig5[comp]):.4f}")

n_neg = d.cannabis.n_cannabis_negative.to_numpy(float)
Z_neg = standardise(d.npki.to_numpy(), n_neg)
for comp in ["LC1", "LC2", "LC3"]:
    q = loadings_from_drug_means(d.figS2[comp], Z_neg, n_neg)
    claim("Drug means", f"cannabis-negative {comp} loadings (Fig. S3) vs implied by drug means (Fig. S2)",
          "r ≥ 0.987", f"r = {r(q, d.figS3[comp]):.4f}")

# ---------------------------------------------------------------- effective sample size
top3 = ["psilocin", "MDMA", "DMT"]
k = t1.loc[top3, "n_reports"].sum()
claim("Sample size", "psilocin, MDMA and DMT: share of reports", "42% (1,893 of 4,496)",
      f"{k / n.sum():.1%} ({k:,} of {int(n.sum()):,})")
for drug in ["MEM", "2C-B-fly"]:
    claim("Sample size", f"{drug}: share of reports", "under 1%", f"{t1.loc[drug, 'n_reports'] / n.sum():.2%}")
claim("Sample size", "Kish effective number of drugs", "10.8", f"{kish(n):.2f}")

# ---------------------------------------------------------------- saturation
rng = np.random.default_rng(1)
r2 = []
for _ in range(20000):
    X = np.c_[np.ones(21), rng.standard_normal((21, 19))]
    y = rng.standard_normal(21)
    resid = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
    r2.append(1 - resid @ resid / ((y - y.mean()) @ (y - y.mean())))
claim("Saturation", "expected R² of an unrelated response, 19 predictors, 21 observations",
      "19/20 = 0.95", f"{19 / 20:.3f} (simulated: {np.mean(r2):.3f})")

# ---------------------------------------------------------------- report length
s = t1.n_sentences
top4 = ["psilocin", "MDMA", "DMT", "LSD"]
claim("Embeddings", "psilocin, MDMA, DMT and LSD: share of sentences", "60%", f"{s[top4].sum() / s.sum():.1%}")
claim("Embeddings", "MEM and DOET: share of sentences", "0.3%", f"{s[['MEM', 'DOET']].sum() / s.sum():.2%}")

length = t1.sentences_per_report
claim("Embeddings", "LC2 drug means (Fig. 4) vs sentences per report", "r = 0.96",
      f"r = {r(d.fig4.LC2, length):.3f} (printed ratios: {r(d.fig4.LC2, t1.sentences_per_report_printed):.3f})")
claim("Embeddings", "cannabis-negative LC1 drug means (Fig. S2) vs sentences per report", "r = 0.92",
      f"r = {r(d.figS2.LC1, length):.3f}")
claim("Embeddings", "sentences per report vs number of reports (Spearman)", "ρ = 0.83",
      f"ρ = {spearmanr(length, t1.n_reports)[0]:.3f}")

# ---------------------------------------------------------------- word lists (ESM2; summarised by extract.py)
w = d.word_lists
claim("Embeddings", "length of the LC2 word lists", "153 words each",
      f"{w.loc[2, 'n_words_plus']} and {w.loc[2, 'n_words_minus']}")
claim("Embeddings", "LC2 word lists: median log10 frequency per billion words (+ vs −)", "5.2 vs 3.0",
      f"{w.loc[2, 'median_log10_freq_plus']:.2f} vs {w.loc[2, 'median_log10_freq_minus']:.2f}")
claim("Embeddings", "LC2 word lists: AUC", "0.96", f"{w.loc[2, 'auc']:.3f}")
others = w.auc.drop(2)
outside = others[(others.round(2) < 0.41) | (others.round(2) > 0.66)]
claim("Embeddings", "AUC of the other components", "most between 0.41 and 0.66",
      f"{len(others) - len(outside)} of {len(others)} in [0.41, 0.66]; outside: "
      + ", ".join(f"LC{c} {a:.2f}" for c, a in outside.items()))
claim("Embeddings", "LC2 labels (LLM)", "‘Perceptual Distortion’, ‘Physical Sensation and Emotional Response’",
      f"‘{w.loc[2, 'llm_label_plus']}’ / ‘{w.loc[2, 'llm_label_minus']}’")
metadata = [f"LC{c}{sign}" for c, row in w.iterrows()
            for sign, col in (("+", "page_metadata_plus"), ("−", "page_metadata_minus")) if row[col]]
claim("Embeddings", "word lists containing Erowid page metadata", "LC3−, LC5+, LC9+, LC12−, LC18+",
      ", ".join(metadata))

# ---------------------------------------------------------------- two-drug contrasts
claim("Contrasts", "LC1 drug means of MDMA and DMT (Fig. 4)", "about −2.1 and +1.6",
      f"{d.fig4.loc['MDMA', 'LC1']:.2f} and {d.fig4.loc['DMT', 'LC1']:+.2f}")
claim("Contrasts", "LC1 loadings vs DMT − MDMA standardised affinities", "r = 0.93",
      f"r = {r(d.fig5.LC1, z['DMT'] - z['MDMA']):.3f}")
claim("Contrasts", "LC3 loadings vs LSD − DMT standardised affinities", "r = 0.97",
      f"r = {r(d.fig5.LC3, z['LSD'] - z['DMT']):.3f}")
mdma = d.npki.loc["MDMA"]
claim("Contrasts", "receptors (of 19) with non-zero MDMA affinity", "five: α2A, α2B, α2C, I1, M3",
      f"{(mdma > 0).sum()}: {', '.join(mdma[mdma > 0].index)}")
sert = d.affinities.query("drug == 'MDMA' and receptor == 'SERT'").iloc[0]
claim("Contrasts", "MDMA at SERT in ESM1", "no affinity (Ki > 10 µM)", f"Ki = {sert.Ki} nM, npKi = {sert.npKi}")

# ---------------------------------------------------------------- figures
for comp, pair in [("LC1", ("DMT", "MDMA")), ("LC3", ("LSD", "DMT"))]:
    u = Z @ d.fig5[comp].to_numpy()
    t = d.fig4[comp].to_numpy()
    keep = ~np.isin(d.drugs, pair)
    claim("Fig. 2", f"{comp} t–u correlation, all 21 drugs / without {pair[0]} and {pair[1]}", "(figure)",
          f"r = {r(u, t):.2f} / {r(u[keep], t[keep]):.2f}")

# ---------------------------------------------------------------- output
lines = ["| Section | Quantity | In the comment | Recomputed |", "|---|---|---|---|"]
lines += [f"| {a} | {b} | {c} | {v} |" for a, b, c, v in rows]
out = Path(__file__).parent / "results" / "numbers.md"
out.parent.mkdir(exist_ok=True)
out.write_text("\n".join(lines) + "\n")
print("\n".join(lines))
