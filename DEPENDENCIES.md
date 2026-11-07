# Software dependencies

The repository does not contain source datasets. Install dependencies before placing licensed or controlled data in each project's local `Data` directory.

## Python projects

Projects 2 and 5 require Python 3.10 or newer. Create an isolated environment and install the documented packages:

```text
python -m venv .venv
.venv/Scripts/python -m pip install -r requirements-python.txt
```

On macOS or Linux, use `.venv/bin/python` instead. Package versions are lower compatibility bounds, not claims about the versions used to generate archived scientific results. Record the resolved environment with `python -m pip freeze` when reproducing an analysis.

## R projects

Projects 1, 3 and 4 require R 4.3 or newer. CRAN dependencies are listed in `requirements-r.txt`. Bioconductor packages must be installed through `BiocManager`:

```r
install.packages("BiocManager")
BiocManager::install(c("AnnotationDbi", "BPCells", "ComplexHeatmap", "edgeR", "fgsea", "org.Hs.eg.db", "speckle"))
```

Project 1 also requires Seurat. Project 4's complete figure atlas uses `ComplexHeatmap`, `circlize` and `RColorBrewer`. Use `renv::snapshot()` in a reproduction environment to record exact resolved package versions; no lockfile is supplied because the original environment was not captured and an invented lockfile would be misleading.

## Lightweight validation

The automated tests do not run the scientific analyses or require source datasets:

```text
python -m unittest discover -s tests -v
```

R syntax validation requires a local R installation and can be performed with `Rscript -e "parse(file='path/to/script.R')"`. Full numerical reproduction requires the documented datasets, compatible software environments and an explicitly initiated analysis run.
