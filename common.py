"""Load the extracted data and the few quantities shared by check_numbers.py and figures.py."""
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

DATA = Path(__file__).parent / "data"


def load():
    t1 = pd.read_csv(DATA / "table1.csv", index_col="drug", comment="#")
    drugs = list(t1.index)
    fig5 = pd.read_csv(DATA / "fig5_loadings.csv", index_col="receptor", comment="#")
    receptors = list(fig5.index)

    aff = pd.read_csv(DATA / "affinities_esm1.csv", comment="#")
    npki = aff.pivot(index="drug", columns="receptor", values="npKi").loc[drugs, receptors]

    return SimpleNamespace(
        drugs=drugs,
        receptors=receptors,
        table1=t1,
        n=t1.n_reports.to_numpy(float),
        affinities=aff,
        npki=npki,
        cannabis=pd.read_csv(DATA / "fig1_cannabis.csv", index_col="drug", comment="#").loc[drugs],
        fig4=pd.read_csv(DATA / "fig4_drug_means.csv", index_col="drug", comment="#").loc[drugs],
        fig5=fig5,
        figS2=pd.read_csv(DATA / "figS2_drug_means_cannabis_negative.csv", index_col="drug", comment="#").loc[drugs],
        word_lists=pd.read_csv(DATA / "word_list_summary.csv", index_col="component", comment="#"),
        figS3=pd.read_csv(DATA / "figS3_loadings_cannabis_negative.csv", index_col="receptor", comment="#").loc[receptors],
    )


def standardise(Y, w):
    """Standardise each column across drugs, weighting drugs by w (report counts).

    This equals report-level standardisation, where each drug's profile is repeated once per report,
    as in the published PLSR."""
    Y = np.asarray(Y, float)
    mean = (w[:, None] * Y).sum(0) / w.sum()
    sd = np.sqrt((w[:, None] * (Y - mean) ** 2).sum(0) / w.sum())
    return (Y - mean) / sd


def kish(w):
    """Kish's effective sample size."""
    return w.sum() ** 2 / (w ** 2).sum()


def loadings_from_drug_means(t, Z, w):
    """Receptor loading direction implied by drug means alone: q = sum_m n_m t_m z_m.

    For report-level PLS, t'Y = sum_m n_m * mean score of drug m * profile of drug m, so the y-loadings
    are proportional to this (any component, after any deflation of X)."""
    return Z.T @ (w * np.asarray(t, float))


def r(a, b):
    return float(np.corrcoef(a, b)[0, 1])
