import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from matplotlib.lines import Line2D

# --- SETUP PATHS ---
# Defining paths relative to the project root
input_folder = "results/enrichment"
output_folder = "results/plots"

# Create output folder if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# List all .tsv files in the enrichment results folder
files_to_plot = [f for f in os.listdir(input_folder) if f.endswith(".tsv")]

if not files_to_plot:
    print(f"No enrichment results found in {input_folder}. Please run the R enrichment script first.")

# --- PLOTTING LOOP ---
for filename in files_to_plot:
    print(f"Generating bubble plot for: {filename}...")
    
    # Load data
    file_path = os.path.join(input_folder, filename)
    df = pd.read_csv(file_path, sep="\t")

    if df.empty:
        print(f"Skipping {filename}: File is empty.")
        continue

    # 1. CALCULATE METRICS
    # -log10 of adjusted p-value for significance visualization
    df["log10_pval"] = -np.log10(df["p_value"])
    
    # Gene Ratio: fraction of the user's query genes found in the term
    df["gene_ratio"] = df["intersection_size"] / df["query_size"]
    
    # Term Ratio: fraction of the total genes in the GO term found in the query
    df["term_ratio"] = df["intersection_size"] / df["term_size"]

    # 2. DATA PREPARATION
    # Sort by p-value (most significant first) and take top 20 terms
    df = df.sort_values("p_value").head(20)
    
    # Shorten long term names if they exceed 50 characters (better for plot labels)
    df["term_name"] = df["term_name"].apply(lambda x: x[:47] + '...' if len(x) > 50 else x)

    # 3. CREATE FIGURE
    fig, ax = plt.subplots(figsize=(10, 8))

    # Create the bubble plot
    # x = significance, y = term name, size = term ratio, hue = gene ratio
    bubble = sns.scatterplot(
        data=df,
        x="log10_pval",
        y="term_name",
        size="term_ratio",
        hue="gene_ratio",
        palette="viridis",
        sizes=(100, 1000), # Range of bubble sizes
        alpha=0.8,
        edgecolor="black",
        ax=ax,
        legend=False  # We will define custom legends below
    )

    # Labels and Title
    ax.set_xlabel("-log10(Adjusted P-Value)", fontweight='bold')
    ax.set_ylabel("Enriched Terms", fontweight='bold')
    
    # Clean title from filename (e.g., 'degs_BP_results.tsv' -> 'degs BP Enrichment')
    clean_title = filename.replace("_results.tsv", "").replace("_", " ")
    ax.set_title(f"Enrichment Analysis: {clean_title}", fontsize=14, fontweight='bold', pad=20)

    # 4. COLORBAR (Gene Ratio)
    norm = plt.Normalize(df["gene_ratio"].min(), df["gene_ratio"].max())
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("Gene Ratio", fontweight='bold')

    # 5. CUSTOM SIZE LEGEND (Term Ratio)
    # Define example sizes for the legend
    size_min, size_max = df["term_ratio"].min(), df["term_ratio"].max()
    example_sizes = np.linspace(size_min, size_max, 3)
    
    markers = [
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor='gray', markersize=np.sqrt(s * 1000), 
               label=f"{s:.2f}")
        for s in example_sizes
    ]

    ax.legend(
        handles=markers,
        title="Term Ratio",
        loc='upper left',
        bbox_to_anchor=(1.05, 0.5), # Place outside the plot
        frameon=True
    )

    ax.grid(True, linestyle='--', alpha=0.5)

    # 6. SAVE PLOT
    output_filename = filename.replace(".tsv", "_bubble_plot.png")
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, output_filename), dpi=300, bbox_inches='tight')
    plt.close()

print(f"\nAll plots have been saved in: {output_folder}")