#!/usr/bin/env bash

# Crear ambiente en conda
conda create -n drexml python=3.10

# Activar conda
conda activate drexml

# Instalar la versión de Drexml
pip install --no-cache-dir --no-binary=shap  -i https://test.pypi.org/simple/     drexml --extra-index-url=https://pypi.org/simple

mkdir -p results

# Ejecutar drexml
drexml run disease.env


# plot

${CONDA_RUN} drexml plot \
 results/shap_selection_symbol.tsv.gz \
 results/shap_summary_symbol.tsv.gz \
 results/stability_results_symbol.tsv.gz \
 results/
