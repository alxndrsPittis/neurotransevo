# Reproduce the single-cell and Boltz-2 analyses and figures.
#   conda env create -f environment.yml && conda activate neurotransevo
#   make all            # everything (matrix building needs a large-memory machine for the mouse atlas)
#   make figures        # re-plot from existing results
RSCRIPT ?= Rscript
PYTHON  ?= python

.PHONY: all matrices coregulons figures coregulon-figures family-figures boltz-figures report clean

all: matrices coregulons figures report

matrices:
	$(RSCRIPT) scripts/01_build_matrices.R

coregulons: matrices
	$(RSCRIPT) scripts/02_compute_coregulons.R

figures: coregulon-figures family-figures boltz-figures

coregulon-figures:
	$(RSCRIPT) scripts/03_plot_coregulons.R

family-figures:
	$(RSCRIPT) scripts/04_plot_family_heatmaps.R

boltz-figures:
	$(PYTHON) python/boltz/plot_panels.py

report:
	$(PYTHON) scripts/05_summary_tables.py

clean:
	rm -rf results
