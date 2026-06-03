import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA
from sklearn.metrics import auc, roc_curve
from statsmodels.stats.multitest import multipletests


CONTROL_GROUP = "Control"
TREATMENT_GROUP = "Hypoxia"
MIN_PRESENT = 0.5
P_CUTOFF = 0.05
LOG2FC_CUTOFF = 1
VIP_CUTOFF = 1


def load_data(raw_file, meta_file):
    raw = pd.read_csv(raw_file)
    meta = pd.read_csv(meta_file)

    raw = raw.rename(columns={raw.columns[0]: "Name"})
    raw = raw.set_index("Name").T
    raw.index = raw.index.astype(str).str.strip()

    meta = meta.rename(
        columns={
            meta.columns[0]: "SampleID",
            meta.columns[1]: "Group",
        }
    )
    meta["SampleID"] = meta["SampleID"].astype(str).str.strip()
    meta["Group"] = meta["Group"].astype(str).str.strip()
    meta = meta.set_index("SampleID")

    df = meta.join(raw, how="inner")
    metabolites = df.columns.drop("Group")
    df[metabolites] = df[metabolites].apply(pd.to_numeric, errors="coerce")

    if df.empty:
        raise ValueError("No matching sample IDs found between raw data and metadata.")

    missing_groups = {CONTROL_GROUP, TREATMENT_GROUP} - set(df["Group"].unique())
    if missing_groups:
        raise ValueError(f"Missing expected group labels: {sorted(missing_groups)}")

    return df


def filter_data(df):
    metabolites = df.drop(columns="Group")

    present_control = metabolites[df["Group"] == CONTROL_GROUP].notna().mean()
    present_treatment = metabolites[df["Group"] == TREATMENT_GROUP].notna().mean()
    keep = (present_control >= MIN_PRESENT) & (present_treatment >= MIN_PRESENT)

    return pd.concat([df["Group"], metabolites.loc[:, keep]], axis=1)


def fill_missing_values(df):
    x = df.drop(columns="Group")
    return x.fillna(x.mean(numeric_only=True))


def plot_pca(df, outdir):
    x = fill_missing_values(df)
    y = df["Group"]

    pca = PCA(n_components=2)
    scores = pca.fit_transform(x)

    pca_df = pd.DataFrame(scores, columns=["PC1", "PC2"])
    pca_df["Group"] = y.values
    pca_df.to_csv(outdir / "PCA_scores.csv", index=False)

    for group in [CONTROL_GROUP, TREATMENT_GROUP]:
        temp = pca_df[pca_df["Group"] == group]
        plt.scatter(temp["PC1"], temp["PC2"], label=group)

    plt.title("PCA - simulated demo data")
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "PCA.png", dpi=300)
    plt.close()


def ttest_and_volcano(df, outdir):
    control = df[df["Group"] == CONTROL_GROUP].drop(columns="Group")
    treatment = df[df["Group"] == TREATMENT_GROUP].drop(columns="Group")

    rows = []
    for metabolite in control.columns:
        x = control[metabolite].dropna()
        y = treatment[metabolite].dropna()

        if len(x) >= 2 and len(y) >= 2 and x.mean() > 0 and y.mean() > 0:
            p_value = ttest_ind(x, y, equal_var=False).pvalue
            log2fc = np.log2(y.mean() / x.mean())
            rows.append([metabolite, p_value, log2fc])

    result = pd.DataFrame(rows, columns=["Metabolite", "p_value", "log2FC"])
    if result.empty:
        raise ValueError("No metabolites passed the basic t-test input checks.")

    result["FDR"] = multipletests(result["p_value"], method="fdr_bh")[1]
    result.to_csv(outdir / "t_test.csv", index=False)

    plt.scatter(result["log2FC"], -np.log10(result["p_value"]), alpha=0.7)
    plt.axhline(-np.log10(P_CUTOFF), linestyle="--")
    plt.axvline(LOG2FC_CUTOFF, linestyle="--")
    plt.axvline(-LOG2FC_CUTOFF, linestyle="--")
    plt.title("Volcano Plot - simulated demo data")
    plt.xlabel(f"log2FC ({TREATMENT_GROUP}/{CONTROL_GROUP})")
    plt.ylabel("-log10 p-value")
    plt.tight_layout()
    plt.savefig(outdir / "Volcano_plot.png", dpi=300)
    plt.close()

    return result


def calculate_vip(pls):
    t = pls.x_scores_
    w = pls.x_weights_
    q = pls.y_loadings_

    p, _ = w.shape
    ss = np.sum(t**2, axis=0) * np.sum(q**2, axis=0)

    vip = []
    for i in range(p):
        weight = (w[i, :] ** 2) / np.sum(w**2, axis=0)
        score = np.sqrt(p * np.sum(ss * weight) / np.sum(ss))
        vip.append(score)

    return np.array(vip)


def plsda_and_vip(df, outdir):
    x = fill_missing_values(df)
    y = df["Group"].map({CONTROL_GROUP: 0, TREATMENT_GROUP: 1})

    pls = PLSRegression(n_components=2)
    scores = pls.fit_transform(x, y)[0]

    score_df = pd.DataFrame(scores, columns=["Comp1", "Comp2"])
    score_df["Group"] = df["Group"].values
    score_df.to_csv(outdir / "PLS_DA_scores.csv", index=False)

    for group in [CONTROL_GROUP, TREATMENT_GROUP]:
        temp = score_df[score_df["Group"] == group]
        plt.scatter(temp["Comp1"], temp["Comp2"], label=group)

    plt.title("Exploratory PLS scores - simulated demo data")
    plt.xlabel("Component 1")
    plt.ylabel("Component 2")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "PLS_DA_plot.png", dpi=300)
    plt.close()

    vip_df = pd.DataFrame(
        {
            "Metabolite": x.columns,
            "VIP": calculate_vip(pls),
        }
    ).sort_values("VIP", ascending=False)
    vip_df.to_csv(outdir / "VIP_scores.csv", index=False)

    top = vip_df[vip_df["VIP"] > VIP_CUTOFF].head(20)
    if not top.empty:
        plt.barh(top["Metabolite"], top["VIP"])
        plt.gca().invert_yaxis()
        plt.title("Top VIP Features - simulated demo data")
        plt.xlabel("VIP")
        plt.tight_layout()
        plt.savefig(outdir / "VIP_plot.png", dpi=300)
        plt.close()

    return vip_df


def find_biomarkers(ttest_df, vip_df, outdir):
    df = ttest_df.merge(vip_df, on="Metabolite")
    df["FC"] = 2**df["log2FC"]

    biomarkers = df[
        (df["p_value"] < P_CUTOFF)
        & (df["log2FC"].abs() >= LOG2FC_CUTOFF)
        & (df["VIP"] > VIP_CUTOFF)
    ]

    biomarkers.to_csv(outdir / "biomarkers.csv", index=False)
    biomarkers[["Metabolite", "log2FC"]].to_csv(
        outdir / "pathway_input_demo.csv",
        index=False,
    )

    return biomarkers


def plot_roc(df, biomarkers, outdir):
    if biomarkers.empty:
        print("No demo biomarkers found. ROC skipped.")
        return

    top_metabolite = biomarkers.iloc[0]["Metabolite"]
    x = df[top_metabolite].fillna(df[top_metabolite].mean())
    y = df["Group"].map({CONTROL_GROUP: 0, TREATMENT_GROUP: 1})

    fpr, tpr, _ = roc_curve(y, x)
    roc_auc = auc(fpr, tpr)

    plt.plot(fpr, tpr, label=f"AUC = {roc_auc:.2f}")
    plt.plot([0, 1], [0, 1], linestyle="--")
    plt.title(f"ROC for {top_metabolite} - simulated demo data")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "ROC_top_biomarker.png", dpi=300)
    plt.close()


def write_summary(df, filtered, biomarkers, outdir):
    summary = [
        "Demo run summary",
        "================",
        "This output was generated from simulated demo data, not thesis data.",
        f"Samples analyzed: {len(df)}",
        f"Features before filtering: {df.shape[1] - 1}",
        f"Features after filtering: {filtered.shape[1] - 1}",
        f"Features passing demo biomarker criteria: {len(biomarkers)}",
        "",
        "Private thesis values should not be copied into this public output.",
    ]
    (outdir / "run_summary.txt").write_text("\n".join(summary), encoding="utf-8")


def analyze_folder(folder, results_root):
    raw_file = folder / "raw_data.csv"
    meta_file = folder / "metadata.csv"

    if not raw_file.exists() or not meta_file.exists():
        print(f"Skipping {folder.name}: expected raw_data.csv and metadata.csv")
        return

    outdir = results_root / folder.name
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Analyzing {folder.name}")

    df = load_data(raw_file, meta_file)
    filtered = filter_data(df)
    filtered.to_csv(outdir / "filtered_data.csv")

    plot_pca(filtered, outdir)
    ttest_df = ttest_and_volcano(filtered, outdir)
    vip_df = plsda_and_vip(filtered, outdir)
    biomarkers = find_biomarkers(ttest_df, vip_df, outdir)
    plot_roc(filtered, biomarkers, outdir)
    write_summary(df, filtered, biomarkers, outdir)

    print(f"Done: {outdir}")


def main():
    parser = argparse.ArgumentParser(description="Beginner metabolomics demo pipeline.")
    parser.add_argument("--data-root", default="data/demo", help="Folder containing comparison folders.")
    parser.add_argument("--results-root", default="results_demo", help="Folder for generated demo outputs.")
    args = parser.parse_args()

    data_root = Path(args.data_root)
    results_root = Path(args.results_root)
    results_root.mkdir(parents=True, exist_ok=True)

    for folder in sorted(data_root.iterdir()):
        if folder.is_dir():
            analyze_folder(folder, results_root)


if __name__ == "__main__":
    main()
