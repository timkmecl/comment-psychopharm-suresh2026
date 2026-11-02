# Comment on Suresh et al. (2026): code and data

Code and data for the comment *"Comment on Suresh et al. (2026): effective sample size and validation of
receptor–semantic mappings"* (Tim Kmecl). Every number in the comment is computed from values published in
the article below and its supplementary material; no report-level data were used.

> Suresh P, Kins J, Levkovskyy L, Thieme A, Wren P, Ray T, Barrett FS, Carhart-Harris R, Bzdok D (2026)
> Moving beyond 5-HT2A: language models map neurotransmitter receptors to experiential semantics.
> *Psychopharmacology*. https://doi.org/10.1007/s00213-026-07170-0

## Quick check (no source files needed)

The extracted values are included in `data/`, so the numbers and figures can be reproduced directly:

```
pip install numpy pandas scipy cairosvg
python check_numbers.py      # -> results/numbers.md: each number in the comment next to its recomputed value
python figures.py            # -> results/figure1.svg, figure2.svg (and PNGs)
```

## Rerun from scratch

This repeats the extraction from the published files.

**1. Install.** Python 3.9 or later, then

```
pip install -r requirements.txt
```

and the Tesseract OCR program (version 5), which `pytesseract` calls:
`apt install tesseract-ocr` (Debian/Ubuntu), `brew install tesseract` (macOS), or the installer from
https://github.com/UB-Mannheim/tesseract/wiki (Windows).

**2. Download the source files** from the article page (https://doi.org/10.1007/s00213-026-07170-0; the
supplementary files are under "Supplementary Information") and save them in `sources/` under these names:

| Save as | What it is |
|---|---|
| `sources/Suresh_2026_Psychopharmacology.pdf` | the article PDF |
| `sources/213_2026_7170_MOESM1_ESM.xlsx` | ESM1: receptor affinities |
| `sources/213_2026_7170_MOESM2_ESM.xlsx` | ESM2: Table S2, component word lists and labels |
| `sources/213_2026_7170_MOESM4_ESM.docx` | ESM4: Supplementary Figures S1–S3 |
| `sources/en_full.txt` | OpenSubtitles 2018 English word counts. Optional: downloaded automatically if missing, from https://raw.githubusercontent.com/hermitdave/FrequencyWords/master/content/2018/en/en_full.txt |

The ESM files keep the names Springer gives them on download. The article PDF must be the published
version: figures are located by fixed pixel positions on its embedded images.

**3. Run**

```
python extract.py            # sources/ -> data/*.csv (overwrites them) and local/  (about 15 s)
python check_numbers.py      # data/ -> results/numbers.md
python figures.py            # data/ -> results/figure1.svg, figure2.svg (+ PNG)
```

`extract.py` stops with an error if a check fails (see "How the figures are read").

Resulting layout:

```
sources/   the files above                               not published (see sources/README.md)
local/     word lists and word frequencies from ESM2     not published, written by extract.py
data/      extracted values, each CSV headed by its source
results/   numbers.md and the figures
```

## Files

| File | Content |
|---|---|
| `extract.py` | Reads the published values: Table 1, ESM1, ESM2, and Figs. 1, 4, 5, S2, S3 |
| `digitise.py` | Functions for reading figures: OCR of annotated heatmaps, marker positions, colourbars |
| `common.py` | Loads `data/`; weighted standardisation, Kish's effective sample size, loadings implied by drug means |
| `check_numbers.py` | Recomputes each number in the comment next to the value stated there |
| `figures.py`, `svg.py` | The two figures of the comment |

## Data (`data/`)

| File | Source | Content |
|---|---|---|
| `table1.csv` | Table 1 | Reports, sentences and sentences per report for each drug (PDF text) |
| `affinities_esm1.csv` | ESM1; values from Ray (2010) | Ki and npKi per drug and receptor; drug names spelled as in Table 1 |
| `fig1_cannabis.csv` | Fig. 1, Table 1 | % of reports with cannabis keywords (OCR); cannabis-negative reports = reports × (1 − %) |
| `fig5_loadings.csv` | Fig. 5 | Printed receptor loadings, 19 receptors × 19 components (OCR) |
| `figS3_loadings_cannabis_negative.csv` | ESM4, Fig. S3 | As Fig. 5, cannabis-negative model |
| `fig4_drug_means.csv` | Fig. 4 | Drug means on LC1 and LC2 (marker centres) and LC3 (marker colour) |
| `figS2_drug_means_cannabis_negative.csv` | ESM4, Fig. S2 | As Fig. 4, cannabis-negative model |
| `word_list_summary.csv` | ESM2, OpenSubtitles 2018 | Per component: LLM labels, list lengths, median word frequency of each end, AUC, Erowid page metadata |

The word lists of ESM2 are not redistributed; `extract.py` writes them to `local/` and publishes only the
per-component summary.

## How the figures are read

- **Heatmaps** (Figs. 1, 5, S3). Each cell is converted to dark-on-white text relative to its own background
  colour, all cells are read in one OCR call, and minus signs are read from the glyph shape (OCR tends to drop
  them). Loadings are then checked against their cell colours on the colourbar; extraction stops if any value
  disagrees by more than 0.03 (largest differences: 0.008 for Fig. 5, 0.007 for Fig. S3). The cannabis
  percentages are checked by the cannabis-negative total, which must equal the 3,368 reports stated in the
  article.
- **Scatter plots** (Figs. 4, S2). Axes are calibrated from the pixel positions of labelled gridlines.
  Markers are the saturated regions left after removing thin coloured text; each is named from an approximate
  position read by eye (`FIG4_APPROX`, `FIGS2_APPROX` in `extract.py`), and touching markers are split by
  k-means. LC3 is the colourbar value whose colour is nearest to the marker's, with the colourbar calibrated
  from the positions of its tick labels.

Values read from figures are approximate: loadings to the two printed decimals, drug means to about ±0.02
on LC1 and LC3 and ±0.05 on LC2.

## Conventions

- Affinities are the npKi values of ESM1 at the 19 receptors of the published model; missing drug–receptor
  pairs are 0, as in the article.
- Affinities are standardised per receptor across drugs, weighting each drug by its number of reports. This
  equals report-level standardisation, where each report carries its drug's profile.
- Sentences per report are sentences / reports from Table 1 (the printed ratios give the same results to two
  decimals).
- Word frequencies are log10 occurrences per billion words in the full OpenSubtitles 2018 English list
  (734,777,659 words). For each component, the AUC is the probability that a word from the + list is more
  frequent than a word from the − list (Mann–Whitney U divided by the product of the list lengths). The 16
  listed words absent from the frequency list (mostly Erowid page metadata such as `2004expid`, and rare
  terms such as `cyanescens`) are ranked below all others.

## Sources and citation

Data in `data/` are values extracted from the following works, for verification of the comment. Please cite
the original sources when using them.

- Suresh P, Kins J, Levkovskyy L, Thieme A, Wren P, Ray T, Barrett FS, Carhart-Harris R, Bzdok D (2026)
  Moving beyond 5-HT2A: language models map neurotransmitter receptors to experiential semantics.
  *Psychopharmacology*. https://doi.org/10.1007/s00213-026-07170-0 — Table 1, Figs. 1, 4, 5, ESM1, ESM2,
  ESM4 (Figs. S2, S3).
- Ray TS (2010) Psychedelics and the human receptorome. *PLoS ONE* 5:e9019.
  https://doi.org/10.1371/journal.pone.0009019 — original source of the affinities in ESM1 (CC BY).
- Dave H (2018) FrequencyWords: frequency word lists generated from OpenSubtitles 2018.
  https://github.com/hermitdave/FrequencyWords — word counts (CC BY-SA 4.0).

The article and its supplementary files are not included in this repository.

## Licence

- **Code** (`*.py`): MIT, see `LICENSE`.
- **Figures, `results/numbers.md` and this README**: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/),
  © 2026 Tim Kmecl.
- **Values in `data/`** are extracted from the works listed under "Sources and citation" and are not covered by
  these licences. They are provided to verify the comment; cite the original sources when using them.
