# Fibroblast Gene Regulatory Network Inference

Scripts infer networks of transcriptional regulation from gene expression data using the [GRNBoost2 algorithm](https://doi.org/10.1093/bioinformatics/bty916). Original networks are refined to match annotated transcription factor (TF)-target databases as well as genes related to cardiac fibrosis. Additional capabilities include k-fold cross validation and conversion to a database for dynamic simulation using [Netflux](https://github.com/saucermanlab/Netflux).

## Required dependencies

- Python: 3.8 or later
- Pandas: 1.0 or later
- Numpy: 1.17 or later
- Dask distributed: 1.28 or later
- [Arboreto](https://arboreto.readthedocs.io/) containing GRNBoost2 algorithm: 0.1.15 or later

## Data requirements

- Gene expression data: CSV file containing genes (rows) x samples (columns)
  - Currently utilizes CPM data for inference
- Library data: TXT files of annotated TF-target databases as obtained via [Harmonizome](https://maayanlab.cloud/Harmonizome/download)
  - Current databases: [CHEA](https://maayanlab.cloud/Harmonizome/dataset/CHEA+Transcription+Factor+Targets), [TRANSFAC-predicted](https://maayanlab.cloud/Harmonizome/dataset/TRANSFAC+Predicted+Transcription+Factor+Targets), [TRANSFAC-curated](https://maayanlab.cloud/Harmonizome/dataset/TRANSFAC+Curated+Transcription+Factor+Targets), [ENCODE](https://maayanlab.cloud/Harmonizome/dataset/ENCODE+Transcription+Factor+Targets)
