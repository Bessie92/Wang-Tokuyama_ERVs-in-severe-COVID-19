# Import packages 
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import decoupler as dc
import os
from scipy.stats import mannwhitneyu
from scipy.sparse import issparse
from scipy import sparse


# Configurate PDF saving
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
# Set font to Arial
matplotlib.rcParams['font.family'] = 'Arial'
matplotlib.rcParams['text.usetex'] = False

# ========= Create merged anndata containing severe COVID-19 and healthy individuals ========= 
Healthy = "Healthy_h5"
COVID_severe = "COVID_h5"

## Function to load .h5 files and create adata containing sample + condition information 
def load_and_label_data(folder, condition):
    adatas = {}
    for filename in os.listdir(folder):
        if filename.endswith("_filtered_feature_bc_matrix.h5"):
            sample_id = filename.split("_filtered_feature_bc_matrix.h5")[0]
            h5_file = os.path.join(folder, filename)
            sample_adata = sc.read_10x_h5(h5_file)
            sample_adata.var_names_make_unique()
            sample_adata.obs['sample'] = sample_id
            sample_adata.obs['condition'] = condition
            adatas[sample_id] = sample_adata
    return adatas

## Load data from folders
Healthy_adatas = load_and_label_data(Healthy, "healthy")
COVID_severe_adatas = load_and_label_data(COVID_severe, "severe COVID")

## Merge the dictionaries
all_adatas = {**Healthy_adatas, **COVID_severe_adatas}  

## Concatenate all the AnnData objects into a single AnnData object
adata = ad.concat(all_adatas, label="sample")
adata.obs_names_make_unique()

## Print the number of observations per sample and condition to verify the concatenation
print(adata.obs["sample"].value_counts())
print(adata.obs["condition"].value_counts())

# Reorder and save adata 
adata.obs['condition'] = pd.Categorical(adata.obs['condition'], categories=["healthy", "severe COVID"], ordered=True)
adata.obs_names_make_unique()
# Save the raw AnnData object
adata.write("./AnnData/merged_adata.h5ad")

#  ========= Quality control for cells ========= 

## mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")

num_mt_genes = adata.var["mt"].sum()
num_mt_genes

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], inplace=True, log1p=True
)

## Plot a histogram for the total detected genes
plt.figure(figsize=(10, 6))
plt.hist(adata.obs['n_genes_by_counts'], bins=100, edgecolor='black')
plt.xlabel('Total Detected Genes')
plt.ylabel('Number of Cells')
plt.title('Histogram of Total Detected Genes')
plt.grid(False)
plt.axvline(x=200, color='red', linestyle='--', linewidth=2)

## Save as PDF
plt.savefig("./Cell_QC/total_detected_genes_histogram.pdf", format='pdf')
plt.show()

## Plot a histogram for the percentage of mitochondrial genes detected per cell
plt.figure(figsize=(10, 6))
plt.hist(adata.obs['pct_counts_mt'], bins=100, edgecolor='black')
plt.axvline(x=13, color='red', linestyle='--', linewidth=2)  # Example threshold at 5%
plt.xlabel('Percentage of Mitochondrial Genes')
plt.ylabel('Number of Cells (log scale)')
plt.title('Histogram of Percentage of Mitochondrial Genes Detected')
plt.grid(True)
## Save as PDF
plt.savefig("./Cell_QC/percentage_mt_histogram.pdf", format='pdf')
plt.show()

## Define the filtering criteria
mt_threshold = 13
gene_threshold = 200

mt_filter = adata.obs['pct_counts_mt'] > mt_threshold
gene_filter = adata.obs['n_genes_by_counts'] < gene_threshold
combined_filter = mt_filter | gene_filter

num_cells_mt = mt_filter.sum()
num_cells_gene = gene_filter.sum()

overlap_filter = mt_filter & gene_filter
num_cells_overlap = overlap_filter.sum()

print(f"Number of cells with >{mt_threshold}% mitochondrial genes: {num_cells_mt}")
print(f"Number of cells with <{gene_threshold} genes detected: {num_cells_gene}")
print(f"Number of cells with both criteria: {num_cells_overlap}")

plt.figure(figsize=(10, 6))
plt.scatter(adata.obs['n_genes_by_counts'], adata.obs['pct_counts_mt'], c='black', label='Retained cells', alpha=0.5)
plt.scatter(adata.obs.loc[combined_filter, 'n_genes_by_counts'], adata.obs.loc[combined_filter, 'pct_counts_mt'], c='grey', label='Filtered out cells', alpha=0.5)
plt.axvline(x=gene_threshold, color='blue', linestyle='--', linewidth=2, label=f'<{gene_threshold} genes')
plt.axhline(y=mt_threshold, color='green', linestyle='--', linewidth=2, label=f'>{mt_threshold}% mt genes')
plt.xlabel('Number of Genes Detected')
plt.ylabel('Percentage of Mitochondrial Genes')
plt.title('Scatter Plot of Cells by Gene Count and Mitochondrial Gene Percentage')
plt.legend()
plt.grid(True)
plt.savefig("./Cell_QC/scatter_plot.pdf", format='pdf')
plt.show()

# Filter the cells
adata = adata[~combined_filter].copy()
print(f"Number of remaining cells after filtering: {adata.n_obs}")

# Save the filtered Anndata
adata.write("./AnnData/merged_adata_filtered.h5ad")

# ========= Normalization ========= 
## Import proviral ERV names 
ERVs_list = pd.read_csv("proviral_ERV_list.csv")

## Filter adata so that it only contains non-ERV genes
non_ERV_genes = list(set(adata.var_names) - set(ERVs_list.erv_id))
non_ERV_adata = adata[:, non_ERV_genes]
non_ERV_adata

sc.experimental.pp.highly_variable_genes(
        non_ERV_adata, flavor="pearson_residuals", n_top_genes=3000
    )

non_ERV_adata = non_ERV_adata[:, non_ERV_adata.var["highly_variable"]]

non_ERV_adata.layers["raw"] = non_ERV_adata.X.copy()
non_ERV_adata.layers["sqrt_norm"] = np.sqrt(
    sc.pp.normalize_total(non_ERV_adata, inplace=False)["X"]
)

sc.experimental.pp.normalize_pearson_residuals(non_ERV_adata)

# ========= Dimensionality reduction ========= 
sc.pp.pca(non_ERV_adata, n_comps=50)
n_cells = len(non_ERV_adata)

sc.pl.pca(
    non_ERV_adata,
    color=["sample", "sample", "pct_counts_mt", "pct_counts_mt"],
    dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
    ncols=2,
    size=2,
)

sc.pp.neighbors(non_ERV_adata, n_neighbors=50, n_pcs=50)
sc.tl.umap(non_ERV_adata, random_state=123)

sc.tl.leiden(non_ERV_adata, flavor="igraph", n_iterations=2)
for res in [0.02, 0.1, 0.3, 0.5, 0.8, 1.5]:
    sc.tl.leiden(
        non_ERV_adata, flavor="igraph", key_added=f"leiden_res_{res:4.2f}", resolution=res
    )

sc.pl.umap(
    non_ERV_adata,
    color=["leiden_res_0.02", "leiden_res_0.10", "leiden_res_0.30", "leiden_res_0.50", "leiden_res_0.80", "leiden_res_1.50"],
    legend_loc="on data",
)

# ========= Cell type annotation ========= 
markers = dc.op.resource("PanglaoDB", organism="human")
markers

## Filter by canonical_marker and human
markers = markers[
    markers["human"].astype(bool)
    & markers["canonical_marker"].astype(bool)
    & (markers["human_sensitivity"].astype(float) > 0.5)
    & (markers['organ'].isin(['Blood', 'Immune system']))
]

markers = markers[~markers.duplicated(["cell_type", "genesymbol"])]
markers = markers.rename(columns={"cell_type": "source", "genesymbol": "target"})
markers = markers[["source", "target", "human_sensitivity", "human_specificity", "ubiquitiousness"]]
markers
# Save markers being used to annotate
markers.to_csv("markers_used_for_cell_type_annotation.csv", index=False)

dc.mt.ulm(data=non_ERV_adata, net=markers, tmin=3)
non_ERV_adata.obsm["score_ulm"]

score = dc.pp.get_obsm(non_ERV_adata, key="score_ulm")
score

df = dc.tl.rankby_group(adata=score, groupby="leiden_res_0.50", reference="rest", method="t-test_overestim_var")
df = df[df["stat"] > 0]
df

n_ctypes = 3
ctypes_dict = df.groupby("group").head(n_ctypes).groupby("group")["name"].apply(lambda x: list(x)).to_dict()
ctypes_dict

sc.pl.matrixplot(
    adata=score,
    var_names=ctypes_dict,
    groupby="leiden",
    dendrogram=True,
    standard_scale="var",
    colorbar_title="Z-scaled scores",
    cmap="RdBu_r",
)

dict_ann = df[df["stat"] > 0].groupby("group").head(1).set_index("group")["name"].to_dict()
dict_ann

non_ERV_adata.obs["cell_type"] = (
    non_ERV_adata.obs["leiden_res_0.50"]
      .map(dict_ann)          # uses the same cluster labels
      .fillna("Unannotated")  # clusters with no positive markers
      .astype("category")
)

# ========= Figure S6A =========
sc.pl.umap(non_ERV_adata, color='cell_type')
sc.pl.umap(non_ERV_adata, color='condition', save='_condition.pdf')

# ========= Figure 6B =========
obs_df = non_ERV_adata.obs[['condition', 'cell_type']]
cell_type_counts = obs_df.groupby(['condition', 'cell_type']).size().unstack(fill_value=0)
cell_type_proportions = cell_type_counts.div(cell_type_counts.sum(axis=1), axis=0)
num_cell_types = cell_type_proportions.shape[1]
colors = sns.color_palette("tab20", num_cell_types) 

cell_type_proportions.plot(kind='bar', stacked=True, figsize=(10, 6), color=colors)

plt.xlabel('Condition')
plt.ylabel('Proportion of Cell Types')
plt.title('Proportion of Different Cell Types within Each Condition')
plt.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

plt.savefig('cell_type_bar_plot_leiden_0.5.pdf', format='pdf')
plt.show()

# ========= Map UMAP coordinates to complete adata =========
adata = adata[:, np.array(adata.X.sum(axis=0)).flatten() > 0]

adata.layers["raw"] = adata.X.copy()
adata.layers["sqrt_norm"] = np.sqrt(
    sc.pp.normalize_total(adata, inplace=False)["X"]
)

cluster_to_ct = non_ERV_adata.obs['leiden_res_0.50'].map(dict_ann).fillna("unannotated")
adata.obs['cell_type'] = (
    cluster_to_ct
    .reindex(adata.obs_names)
    .fillna("unannotated")
)

adata.obsm['X_umap'] = non_ERV_adata.obsm['X_umap']
adata.write("./AnnData/merged_adata_filtered_with_cell_types_and_clustering.h5ad")

adata = adata[adata.obs['cell_type'] != 'Plasma cells']
adata = adata[adata.obs['cell_type'] != 'unannotated']

common_ERVs = list(set(ERVs_list.erv_id) & set(adata.var_names))
ERV_adata = adata[:, common_ERVs] 
ERV_adata

severe_covid_ERV_adata = ERV_adata[ERV_adata.obs['condition'] == 'severe COVID']
healthy_ERV_adata = ERV_adata[ERV_adata.obs['condition'] == 'healthy']
combined_ERV_adata = ad.concat([healthy_ERV_adata, severe_covid_ERV_adata])

# ========= Figure 6C =========
dotplot = sc.pl.DotPlot(
    combined_ERV_adata,
    ["2637", "1874", "2124",
        "3052", "4830", "2095", "4896","4174", "565", "882", "W-8", "3659", "606", "4745", 
        "4426", "909", "1379", 
        "942", "1179", "2724", "5302", "3743", "814", "1879", "3687", "3104", "3508", "4657", "6080", "1043"],
    log=True,
    #layer="sqrt_norm",
    groupby=['condition', 'cell_type'],
).style(
    cmap='inferno', 
    dot_min=0.01,
    dot_edge_color='black',
    dot_edge_lw=1,
).swap_axes(False)

dotplot.show()
dotplot.savefig("bulk_signatures_by_condition_pearson_residual.pdf")

# ========= Figure 6E =========
adata_E = adata[adata.obs['cell_type'] == 'Erythroid-like and erythroid precursor cells']
adata_E_covid = adata_E[adata_E.obs['condition'] == 'severe COVID']
adata_E_covid

# Define ERV
## Change to 2124 for the bottom row in the figure
gene_A = "2637"

# Extract ERV expression values and ensure they are in a dense array format
gene_A_expression = adata_E_covid[:, gene_A].X
if issparse(gene_A_expression):
    gene_A_expression = gene_A_expression.toarray().flatten()

# Add binary column indicating whether each cell expresses the ERV
adata_E_covid.obs["GeneA_status"] = (gene_A_expression > 0).astype(int)

# Convert to categorical so it will display properly on the dotplot
adata_E_covid.obs["GeneA_status"] = adata_E_covid.obs["GeneA_status"].astype(str)

# Define gene (set) of interest
genes_of_interest = ["TLR1","TLR2","TLR6","JAK2","ANXA3","IL18RAP","ITGAM","CD177","BATF","DYSF","NDRG1","LRRK2","PLA2G4A","CLEC4D","NECTIN2","SYK","PREX1","FES","PLCG2","LDLR","APP","CLU","SNCA","MMP9","BPGM"] 

# Visualize using dotplot
sc.pl.dotplot(adata_E_covid, var_names=genes_of_interest, groupby="GeneA_status",
              cmap="inferno",  dot_min=0.1, dot_max=1, save= 'Erythroid_2637_GO_term_genes_more_than_10_per.pdf')

adata_mac = adata[adata.obs['cell_type'] == 'Macrophages']
adata_mac_covid = adata_mac[adata_mac.obs['condition'] == 'severe COVID']
adata_mac_covid

# Define ERV
## Change to 2124 for the bottom row in the figure
gene_A = "1874"

# Extract ERV expression values and ensure they are in a dense array format
gene_A_expression = adata_mac_covid[:, gene_A].X
if issparse(gene_A_expression):
    gene_A_expression = gene_A_expression.toarray().flatten()  # Convert sparse to dense and flatten to 1D

# Add binary column indicating whether each cell expresses the ERV
adata_mac_covid.obs["GeneA_status"] = (gene_A_expression > 0).astype(int)

# Convert to categorical so it will display properly on the dotplot
adata_mac_covid.obs["GeneA_status"] = adata_mac_covid.obs["GeneA_status"].astype(str)

# Define gene (set) of interest - epigenetic and post-transcriptional silencing of ERVs 
genes_of_interest = ["RNASE2","HMGB2","PLD1","RAB13","RNASE3","CNR2","SEMA3C","NRG1","CX3CR1","ITGA1","CYP19A1","C5","CCR2","KIT","CH25H","EFNA5","F3","CMTM8","ITGA9","CMKLR1","TNFRSF11A","CXCL3","FGFR1","SLAMF1","ITGB3","MET","GCNT1","CD99L2","ITGA7","TRIM55","SELP","MMP14","ADTRP","ITGA2B"] 

# Visualize using dotplot
sc.pl.dotplot(adata_mac_covid, var_names=genes_of_interest, groupby="GeneA_status",
              cmap="inferno",   dot_min=0.1, dot_max=1, save= 'Macrophage_1874_GO_genes_1_more_than_10_per.pdf')

# ========= Figure 6D =========

# Define ERV
gene_A = "2637"

# Extract ERV expression values and ensure they are in a dense array format
gene_A_expression = adata_E_covid[:, gene_A].X
if issparse(gene_A_expression):
    gene_A_expression = gene_A_expression.toarray().flatten()  # Convert sparse to dense and flatten to 1D

# Add binary column indicating whether each cell expresses the ERV
adata_E_covid.obs["GeneA_status"] = (gene_A_expression > 0).astype(int)

# Convert to categorical so it will display properly on the dotplot
adata_E_covid.obs["GeneA_status"] = adata_E_covid.obs["GeneA_status"].astype(str)

# Define gene (set) of interest - epigenetic and post-transcriptional silencing of ERVs 
genes_of_interest = ["IL13RA1", "IL1RAP", "PELI1", "IL18R1", "IL1RL1", "IL1R2", "IRAK3","IL1R1","IL18RAP"] 

# Visualize using dotplot
sc.pl.dotplot(adata_E_covid, var_names=genes_of_interest, groupby="GeneA_status",
              cmap="inferno",  dot_min=0.01, dot_max=1, save= 'Erythroid_2637_il1_genes.pdf')


# Define ERV
gene_A = "1874"

# Extract ERV expression values and ensure they are in a dense array format
gene_A_expression = adata_mac_covid[:, gene_A].X
if issparse(gene_A_expression):
    gene_A_expression = gene_A_expression.toarray().flatten()  # Convert sparse to dense and flatten to 1D

# Add binary column indicating whether each cell expresses the ERV
adata_mac_covid.obs["GeneA_status"] = (gene_A_expression > 0).astype(int)

# Convert to categorical so it will display properly on the dotplot
adata_mac_covid.obs["GeneA_status"] = adata_mac_covid.obs["GeneA_status"].astype(str)

# Define gene (set) of interest - epigenetic and post-transcriptional silencing of ERVs 
genes_of_interest = ["IL13RA1", "IL1RAP", "PELI1", "IL18R1", "IL1RL1", "IL1R2", "IRAK3","IL1R1","IL18RAP"] 

# Visualize using dotplot
sc.pl.dotplot(adata_mac_covid, var_names=genes_of_interest, groupby="GeneA_status",
              cmap="inferno",  dot_min=0.01, dot_max=1, save= 'Mac_1874_il1_genes.pdf')



## Figures S6C and S6D were generated using the same code as above, with the genes shown in the respective figures.
