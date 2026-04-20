import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt

# script to summarize results from the open chromatin pipeline
# reads in all the output files and makes plots and a summary text file

# defining expected command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="summarize pipeline results")
    parser.add_argument("--openChromDir", required=True, help="folder with BED files from bedtools")
    parser.add_argument("--goDir", required=True, help="folder with GO CSV files from rGREAT")
    parser.add_argument("--homerDir", required=True, help="folder with HOMER motif results")
    parser.add_argument("--outputDir", required=True, help="folder to save summary files")
    return parser.parse_args()

def count_peaks(bed_dir):
    # count how many peaks are in each BED file
    # each line in a BED file is one peak
    files = {
        "human_atac": "human_atac.bed",
        "mouse_atac": "mouse_atac.bed",
        "mouse_shared_in_human": "mouse_shared_in_human.humanCoords.bed",
        "human_shared_in_mouse": "human_shared_in_mouse.mouseCoords.bed",
        "mouse_open_human_closed": "mouse_open_human_closed.humanCoords.bed",
        "human_open_mouse_closed": "human_open_mouse_closed.mouseCoords.bed",
        "mouse_no_ortholog": "mouse_no_ortholog.mouseCoords.bed",
        "human_no_ortholog": "human_no_ortholog.humanCoords.bed"
    }

    # getting counts for each bed file in open chrom directory
    counts = {}
    for label, filename in files.items():
        filepath = os.path.join(bed_dir, filename)
        if os.path.exists(filepath):
            with open(filepath) as f:
                counts[label] = sum(1 for line in f if line.strip())
        else:
            counts[label] = 0
            print("WARNING: could not find " + filepath)
    return counts

def count_promoters_enhancers(homer_dir):
    # read each promoter and enhancer output file from homer
    # count number of promoters vs enhancers for each group
    categories = ["shared_human", "shared_mouse", "unique_human", "unique_mouse"]
    counts = {}

    # getting counts of peaks in each type of promoter/enhancer file in the homer directory
    for cat in categories:
        for ptype in ("promoter", "enhancer"):
            bed = os.path.join(homer_dir, cat + "_" + ptype + "_peaks.bed")
            if os.path.exists(bed):
                with open(bed) as f:
                    counts[cat + "_" + ptype] = sum(1 for line in f if line.strip())
            else:
                counts[cat + "_" + ptype] = 0
                print("WARNING: could not find " + bed)

    return counts

def get_top_go_terms(go_dir, n=5):
    # read each rGREAT output CSV and get the top significant GO terms
    # only keep terms with adjusted p-value below 0.05
    results = {}
    files = {
        "mouse_shared_in_human": "mouse_shared_in_human_GREAT_GO_BP.csv",
        "human_shared_in_mouse": "human_shared_in_mouse_GREAT_GO_BP.csv",
        "mouse_open_human_closed": "mouse_open_human_closed_GREAT_GO_BP.csv",
        "human_open_mouse_closed": "human_open_mouse_closed_GREAT_GO_BP.csv",
        "mouse_no_ortholog": "mouse_no_ortholog_GREAT_GO_BP.csv",
        "human_no_ortholog": "human_no_ortholog_GREAT_GO_BP.csv",
        "mouse_allpeaks": "liver_mouse_GO_BP_allpeaks.csv",
        "human_allpeaks": "liver_human_GO_BP_allpeaks.csv"
    }

    # filtering out csv files for only significant results (p val < 0.05
    for label, filename in files.items():
        filepath = os.path.join(go_dir, filename)
        if os.path.exists(filepath):
            df = pd.read_csv(filepath)
            # only keeping top 5 significant terms
            sig = df[df["p_adjust"] < 0.05].head(n)
            results[label] = sig
        else:
            results[label] = pd.DataFrame()

    return results

def get_top_motifs(homer_dir, n=5):
    # read the top known motifs from each HOMER knownResults.txt file
    results = {}
    subdirs = [
        "shared_human_enhancer_motifs",
        "shared_mouse_enhancer_motifs",
        "shared_human_promoter_motifs",
        "shared_mouse_promoter_motifs",
        "unique_human_enhancer_motifs",
        "unique_mouse_enhancer_motifs",
        "unique_human_promoter_motifs",
        "unique_mouse_promoter_motifs"
    ]

    # reading known motif results and only keeping top 5
    for subdir in subdirs:
        known_results = os.path.join(homer_dir, subdir, "knownResults.txt")
        if os.path.exists(known_results):
            df = pd.read_csv(known_results, sep="\t")
            df.columns = df.columns.str.strip()
            results[subdir] = df.head(n)
        else:
            results[subdir] = pd.DataFrame()

    return results

def get_top_denovo_motifs(homer_dir, n=5):
    # parsing de novo motifs from homerMotifs.all.motifs output
    # each motif header line starts with > and contains
    # headers: >consensus  name  log-odds  log-p-value  0,  T:...,B:...,P:raw_pval
    results = {}
    subdirs = [
        "shared_human_enhancer_motifs",
        "shared_mouse_enhancer_motifs",
        "shared_human_promoter_motifs",
        "shared_mouse_promoter_motifs",
        "unique_human_enhancer_motifs",
        "unique_mouse_enhancer_motifs",
        "unique_human_promoter_motifs",
        "unique_mouse_promoter_motifs"
    ]

    # getting de novo motifs for each OCR category
    for subdir in subdirs:
        motif_file = os.path.join(homer_dir, subdir, "homerMotifs.all.motifs")
        if not os.path.exists(motif_file):
            results[subdir] = pd.DataFrame()
            continue

        # making list of all motifs
        motifs = []
        with open(motif_file) as f:
            for line in f:
                # checking that line represents a new motif and has all info for motif
                if not line.startswith(">"):
                    continue
                parts = line[1:].strip().split("\t")
                if len(parts) < 6:
                    continue
                consensus = parts[0]
                name = parts[1]
                log_pval = parts[3]
                # raw p-value is at the end of the last field after "P:"
                raw_pval = ""
                for field in parts[5].split(","):
                    if field.startswith("P:"):
                        raw_pval = field[2:]
                motifs.append({
                    "name": name,
                    "consensus": consensus,
                    "log_p_value": log_pval,
                    "p_value": raw_pval
                })

        # only keeping top 5 motifs
        results[subdir] = pd.DataFrame(motifs).head(n) if motifs else pd.DataFrame()

    return results

def plot_peak_counts(counts, output_dir):
    # make a bar chart showing how many peaks are in each category
    labels = [
        "Mouse shared\n(human coords)",
        "Human shared\n(mouse coords)",
        "Mouse open\nhuman closed",
        "Human open\nmouse closed",
        "Mouse no\northolog",
        "Human no\northolog"
    ]
    keys = [
        "mouse_shared_in_human",
        "human_shared_in_mouse",
        "mouse_open_human_closed",
        "human_open_mouse_closed",
        "mouse_no_ortholog",
        "human_no_ortholog"
    ]

    # blue for shared, orange for species-specific, gray for no ortholog
    colors = [
        "blue", "blue",
        "orange", "orange",
        "gray", "gray"
    ]

    # using values in count dict from count_peaks function results
    values = [counts.get(k, 0) for k in keys]

    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.bar(labels, values, color=colors, edgecolor="white")

    # add the count number on top of each bar
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 200, str(val),
            ha="center", va="bottom", fontsize=9)

    # add legend since bars of same color share a label
    bars[0].set_label("shared peaks")
    bars[2].set_label("species-specific peaks")
    bars[4].set_label("no ortholog")
    ax.legend()

    ax.set_ylabel("number of peaks")
    ax.set_title("open chromatin peak counts by category")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "peak_counts.png"), dpi=150)
    plt.close()
    print("saved peak_counts.png")

def plot_promoter_enhancer(pe_counts, output_dir):
    # grouped bar chart comparing promoter vs enhancer counts for each category
    categories = ["shared_human", "shared_mouse", "unique_human", "unique_mouse"]
    labels = ["shared (human)", "shared (mouse)", "unique (human)", "unique (mouse)"]

    # getting count of promoters / enhancers from count promoters and enhancers function results
    promoter_counts = [pe_counts.get(c + "_promoter", 0) for c in categories]
    enhancer_counts = [pe_counts.get(c + "_enhancer", 0) for c in categories]


    x = range(len(categories))
    width = 0.35
    fig, ax = plt.subplots(figsize=(10, 6))

    # plotting promoters and enhancers side by side for each category
    ax.bar([i - width/2 for i in x], promoter_counts, width, label="promoters", color="green")
    ax.bar([i + width/2 for i in x], enhancer_counts, width, label="enhancers", color="purple")

    ax.set_xticks(list(x))
    ax.set_xticklabels(labels)
    ax.set_ylabel("number of peaks")
    ax.set_title("promoter vs enhancer peaks by category")
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "promoter_vs_enhancer.png"), dpi=150)
    plt.close()
    print("saved promoter_vs_enhancer.png")

def plot_go_terms(go_results, output_dir, n=5):
    import math
    plot_order = [
        ("mouse_shared_in_human", "mouse shared (human coords)"),
        ("human_shared_in_mouse", "human shared (mouse coords)"),
        ("mouse_open_human_closed", "mouse open / human closed"),
        ("human_open_mouse_closed", "human open / mouse closed"),
        ("mouse_no_ortholog", "mouse no ortholog"),
        ("human_no_ortholog", "human no ortholog"),
        ("mouse_allpeaks", "all mouse peaks"),
        ("human_allpeaks", "all human peaks"),
    ]

    # getting data for each cateogory of GO plots
    panels = [(key, label, go_results[key]) for key, label in plot_order
              if key in go_results and go_results[key] is not None and not go_results[key].empty]

    if not panels:
        print("no significant GO terms to plot")
        return

    # calculating the grid dimensions of subplots so they fit into 2 cols
    ncols = 2
    nrows = math.ceil(len(panels) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 4 * nrows))
    axes = axes.flatten() if len(panels) > 1 else [axes]

    # only taking top 5 enriched pathways
    for ax, (key, label, df) in zip(axes, panels):
        top = df.head(n).copy()
        # shortening GO pathway description
        top["desc_short"] = top["description"].str.slice(0, 50)
        # reversing order so highest fold enrichment is first
        top = top.iloc[::-1]

        ax.barh(top["desc_short"], top["fold_enrichment"], color="steelblue")
        ax.set_xlabel("fold enrichment")
        ax.set_title(label, fontsize=10)
        ax.tick_params(axis="y", labelsize=8)

    for ax in axes[len(panels):]:
        ax.set_visible(False)

    plt.suptitle("GO Biological Process enrichment - fold enrichment (top " + str(n) + " per category)", fontsize=12, y=1.01)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "go_enrichment.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print("saved go_enrichment.png")

def write_summary(counts, pe_counts, go_results, motif_results, denovo_motif_results, output_dir):
    # write all results to a single text file pipeline_summary.txt

    # initializing where summary text file will be (output dir in summary dir)
    summary_path = os.path.join(output_dir, "pipeline_summary.txt")

    with open(summary_path, "w") as f:

        f.write("pipeline summary - comparative open chromatin analysis\n")
        f.write("-" * 60 + "\n\n")

        ### section 1 - peak counts per category
        f.write("1. peak mapping results\n")
        f.write("-" * 60 + "\n\n")

        # getting counts of peaks for each category
        total_human = counts.get("human_atac", 0)
        total_mouse = counts.get("mouse_atac", 0)
        mouse_shared = counts.get("mouse_shared_in_human", 0)
        human_shared = counts.get("human_shared_in_mouse", 0)
        mouse_unique = counts.get("mouse_open_human_closed", 0)
        human_unique = counts.get("human_open_mouse_closed", 0)
        mouse_no_ortho = counts.get("mouse_no_ortholog", 0)
        human_no_ortho = counts.get("human_no_ortholog", 0)

        # calculate percentages of shared peaks, avoiding zero division
        if total_mouse > 0:
            mouse_shared_pct = mouse_shared / total_mouse * 100
        else:
            mouse_shared_pct = 0

        if total_human > 0:
            human_shared_pct = human_shared / total_human * 100
        else:
            human_shared_pct = 0

        f.write("total human peaks: " + str(total_human) + "\n")
        f.write("total mouse peaks: " + str(total_mouse) + "\n\n")
        f.write("mouse peaks with human ortholog: " + str(mouse_shared + mouse_unique) + "\n")
        f.write("human peaks with mouse ortholog: " + str(human_shared + human_unique) + "\n\n")
        f.write("mouse peaks shared (open in both): " + str(mouse_shared) +
                " (" + str(round(mouse_shared_pct, 1)) + "% of mouse peaks)\n")
        f.write("human peaks shared (open in both): " + str(human_shared) +
                " (" + str(round(human_shared_pct, 1)) + "% of human peaks)\n")
        f.write("mouse open / human closed: " + str(mouse_unique) + "\n")
        f.write("human open / mouse closed: " + str(human_unique) + "\n")
        f.write("mouse peaks with no ortholog: " + str(mouse_no_ortho) + "\n")
        f.write("human peaks with no ortholog: " + str(human_no_ortho) + "\n\n")

        ### section 2 - GO terms (categorized peaks)
        f.write("2. top enriched biological processes - categorized peaks (rGREAT GO:BP)\n")
        f.write("-" * 60 + "\n")
        f.write("note: p_adj = 0.00e+00 indicates floating point underflow (p-value too small to represent -- extremely significant)\n\n")

        go_labels = {
            "mouse_shared_in_human": "mouse shared in human (human coords)",
            "human_shared_in_mouse": "human shared in mouse (mouse coords)",
            "mouse_open_human_closed": "mouse open / human closed",
            "human_open_mouse_closed": "human open / mouse closed",
            "mouse_no_ortholog": "mouse no ortholog",
            "human_no_ortholog": "human no ortholog"
        }

        # getting go results for each cateogry and writing top 5 with fold and adjusted p val
        for key, label in go_labels.items():
            f.write(label + "\n")
            df = go_results.get(key, pd.DataFrame())
            if df is not None and not df.empty:
                for _, row in df.iterrows():
                    desc = str(row["description"])[:45]
                    f.write("  " + row["id"] + "  " + desc +
                            "  fold=" + str(round(row["fold_enrichment"], 2)) +
                            "  p_adj=" + format(row["p_adjust"], ".2e") + "\n")
            else:
                f.write("  no significant GO terms found (p_adjust > 0.05)\n")
            f.write("\n")

        ### section 2b - GO terms (full raw peak sets)
        f.write("2b. top enriched biological processes - all peaks (rGREAT GO:BP)\n")
        f.write("-" * 60 + "\n")
        f.write("note: p_adj = 0.00e+00 indicates floating point underflow (p-value too small to represent -- extremely significant)\n\n")

        allpeaks_labels = {
            "mouse_allpeaks": "all mouse ATAC peaks",
            "human_allpeaks": "all human ATAC peaks"
        }

        # same thing as 2a
        for key, label in allpeaks_labels.items():
            f.write(label + "\n")
            df = go_results.get(key, pd.DataFrame())
            if df is not None and not df.empty:
                for _, row in df.iterrows():
                    desc = str(row["description"])[:45]
                    f.write("  " + row["id"] + "  " + desc +
                            "  fold=" + str(round(row["fold_enrichment"], 2)) +
                            "  p_adj=" + format(row["p_adjust"], ".2e") + "\n")
            else:
                f.write("  no significant GO terms found (p_adjust > 0.05)\n")
            f.write("\n")

        ### section 3 - promoter vs enhancer counts
        f.write("3. promoters vs enhancers\n")
        f.write("-" * 60 + "\n\n")

        pe_labels = [
            ("shared_human", "shared peaks (human coords)"),
            ("shared_mouse", "shared peaks (mouse coords)"),
            ("unique_human", "human-unique peaks"),
            ("unique_mouse", "mouse-unique peaks")
        ]

        # getting count of promoters and enhancers for each category and total identified
        for key, label in pe_labels:
            prom = pe_counts.get(key + "_promoter", 0)
            enh = pe_counts.get(key + "_enhancer", 0)
            total = prom + enh

            # finding percentage of promoters vs enhancers
            if total > 0:
                prom_pct = round(prom / total * 100, 1)
                enh_pct = round(enh / total * 100, 1)
            else:
                prom_pct = 0
                enh_pct = 0

            # writing in results to summary file
            f.write(label + "\n")
            f.write("  promoters: " + str(prom) + " (" + str(prom_pct) + "%)\n")
            f.write("  enhancers: " + str(enh) + " (" + str(enh_pct) + "%)\n")
            f.write("  total: " + str(total) + "\n\n")

        ### section 4 - top known motifs
        f.write("4. top enriched known motifs in promoters and enhancers\n")
        f.write("-" * 60 + "\n\n")

        motif_labels = [
            ("shared_human_promoter_motifs", "shared human promoters"),
            ("shared_human_enhancer_motifs", "shared human enhancers"),
            ("shared_mouse_promoter_motifs", "shared mouse promoters"),
            ("shared_mouse_enhancer_motifs", "shared mouse enhancers"),
            ("unique_human_promoter_motifs", "human-unique promoters"),
            ("unique_human_enhancer_motifs", "human-unique enhancers"),
            ("unique_mouse_promoter_motifs", "mouse-unique promoters"),
            ("unique_mouse_enhancer_motifs", "mouse-unique enhancers")
        ]

        # getting top 5 known motifs and writing in motif name and p val
        for key, label in motif_labels:
            f.write(label + "\n")
            df = motif_results.get(key, pd.DataFrame())
            if df is not None and not df.empty:
                for _, row in df.iterrows():
                    name = str(row.get("Motif Name", ""))[:45]
                    pval = row.get("P-value", "N/A")
                    f.write("  " + name + "  p=" + str(pval) + "\n")
            else:
                f.write("  no motif results found\n")
            f.write("\n")

        ### section 5 - top de novo motifs
        f.write("5. top de novo motifs in promoters and enhancers (HOMER)\n")
        f.write("-" * 60 + "\n\n")

        # getting top 5 de novo discovered motifs and writing in name, consensus seq and p val
        for key, label in motif_labels:
            f.write(label + "\n")
            df = denovo_motif_results.get(key, pd.DataFrame())
            if df is not None and not df.empty:
                for _, row in df.iterrows():
                    name = str(row.get("name", ""))[:45]
                    consensus = str(row.get("consensus", ""))
                    pval = row.get("p_value", "N/A")
                    f.write("  " + name + "  consensus=" + consensus + "  p=" + str(pval) + "\n")
            else:
                f.write("  no de novo motif results found\n")
            f.write("\n")

        f.write("-" * 60 + "\n")
        f.write("end of summary\n")

    print("saved pipeline_summary.txt")

# main function that performs full summary workflow
def main():
    options = parse_arguments()

    # create output folder if it doesn't exist
    os.makedirs(options.outputDir, exist_ok=True)

    # getting peak count summary
    print("counting peaks...")
    counts = count_peaks(options.openChromDir)

    # getting promoter vs enhancer summary
    print("counting promoters and enhancers...")
    pe_counts = count_promoters_enhancers(options.homerDir)

    # getting top 5 go results
    print("reading GO results...")
    go_results = get_top_go_terms(options.goDir)

    # getting top 5 known motifs results
    print("reading HOMER motifs...")
    motif_results = get_top_motifs(options.homerDir)
    denovo_motif_results = get_top_denovo_motifs(options.homerDir)

    # plotting peaks, promoters vs. enhancers and go results
    print("making plots...")
    plot_peak_counts(counts, options.outputDir)
    plot_promoter_enhancer(pe_counts, options.outputDir)
    plot_go_terms(go_results, options.outputDir)

    # writing summary into text file
    print("writing summary...")
    write_summary(counts, pe_counts, go_results, motif_results, denovo_motif_results, options.outputDir)

    print("done! files saved to " + options.outputDir)

if __name__ == "__main__":
    main()