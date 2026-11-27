# Paper Supporting Code and Data

R scripts and supporting datasets for reproducing the main figures (Figures 1–4) of the manuscript.

## Repository Structure

```text
.
├── 图1/                        # Figure 1: Soil properties, EV characterization, and community composition
│   ├── lefse图/                # Signed LEfSe LDA score barplot (LDA_score_LEfSe_signed.R)
│   ├── 囊泡电镜图/             # Extracellular vesicle TEM images and size distribution
│   ├── 土壤理化/               # Soil physicochemical properties radar chart (土壤理化雷达图.R)
│   └── 组成柱状图/             # Phylum-level community composition barplots (组成柱状图.R, 组成组装图（打断）.R)
├── 图2/                        # Figure 2: Carbon-cycling functional gene diversity and pathways
│   ├── α多样性图/              # Alpha diversity of carbon-cycling genes (α多样性.R)
│   ├── 碳循环β多样性图/        # Beta diversity / PCoA based on Bray-Curtis distance (功能基因beta多样性.R)
│   ├── 碳循环热图/             # Circular heatmap of key carbon-cycling genes (环形功能基因热图.r)
│   ├── 碳循环组成柱状图/       # Stacked barplot of carbon-cycling pathways (碳循环组成柱状图.R)
│   └── 功能基因通路图.*        # Functional gene pathway schematic (PPTX / PDF)
├── 图3/                        # Figure 3: Life-history trait & MCOA multi-omics integration
│   ├── annotation/             # Life-history trait annotation database
│   └── gene_base/              # Scripts (LH1–LH3) and outputs for EV, M, and combined 4-group analyses
└── 图4/                        # Figure 4: MAG phylogeny and global plastisphere distribution
    ├── MAG全球地图/            # Global distribution map by strict Genus/Family matching
    └── MAG进化树/              # Rectangular MAG phylogeny with genome size & carbon-gene tracks (MAG进化树.R)
```

## Requirements

- **R** (>= 4.0)
- Common R packages used across scripts: `ggplot2`, `dplyr`, `tidyr`, `readxl`, `openxlsx`, `vegan`, `ggpubr`, `fmsb`, `patchwork`, `ape`, `ade4`

## Usage

1. Open the target figure subfolder as your R working directory (e.g., `setwd("图2/α多样性图")`).
2. Run the corresponding `.R` script.
3. For **Figure 3** (`图3/gene_base/`), run scripts in order: `LH1.Life_history.R` → `LH2.MCOA.R` → `LH3.MCOA_BAR.R` (see `图3/使用说明.txt` for details).
