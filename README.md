# metacleaner <img src="man/figures/logo.svg" align="right" height="139" alt="metacleaner hex logo" />

<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/sbresnahan/metacleaner/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sbresnahan/metacleaner/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/sbresnahan/metacleaner/graph/badge.svg)](https://app.codecov.io/gh/sbresnahan/metacleaner)
[![pkgdown](https://github.com/sbresnahan/metacleaner/actions/workflows/pkgdown.yaml/badge.svg)](https://sbresnahan.github.io/metacleaner/)
[![DOI](https://zenodo.org/badge/572172411.svg)](https://zenodo.org/doi/10.5281/zenodo.7920682)
<!-- badges: end -->

### Automated curation of DNA barcode reference databases for metabarcoding and metagenomics

DNA barcode reference databases built from NCBI GenBank (for example with
[MetaCurator](https://github.com/RTRichar/MetaCurator)) contain mislabeled
sequences. In pollen metabarcoding with plant ITS1/ITS2 databases, this
produces reads assigned to fungi, to the wrong plant genus, or to several
genera at the same identity.

metacleaner checks each sequence in your database. It removes the sequence
if it

1. matches a database of unwanted sequences (e.g. fungal ITS),
2. has no match in a database of wanted sequences (e.g. plant ITS), or
3. matches wanted sequences whose taxonomy disagrees with its own at the rank
   you choose (genus by default).

Everything runs inside R. The sequence search is built in (C++,
multithreaded), so you do not need to install BLAST, Python or command-line
tools.

## Installation

Install from CRAN:

```r
install.packages("metacleaner", dependencies = TRUE)
```

If R reports that package 'metacleaner' is not available, it is still under
review at CRAN. Install it from GitHub instead:

```r
install.packages("remotes")
remotes::install_github("sbresnahan/metacleaner", dependencies = TRUE)
```

Installing from GitHub compiles the package on your computer, so it needs
build tools: [Rtools](https://cran.r-project.org/bin/windows/Rtools/) on
Windows, or the Xcode command line tools on macOS (run
`xcode-select --install` in Terminal). Most Linux systems already have them.

## Quick start

```r
library(metacleaner)

result <- clean_reference(
  query        = "ITS2_database.fasta",   # the database to clean
  reference    = "ITS_plants.fasta",      # sequences you want to keep (e.g. plants)
  contaminants = "ITS_nonplants.fasta",   # sequences you want to remove (e.g. fungi)
  taxonomy     = "ITS2_database.tax",     # or "ncbi" to look taxonomy up online
  level        = "genus",
  threads      = "auto"
)
result                                    # summary
write_results(result, outdir = "cleaned") # cleaned FASTA + taxonomy + report
```

To try it on the example data shipped with the package:

```r
ex <- function(f) system.file("extdata", f, package = "metacleaner")
clean_reference(ex("query.fasta"), ex("reference.fasta"), ex("contaminants.fasta"),
                taxonomy = ex("taxonomy.tax"))
```

To build the reference and contaminant databases from NCBI:

```r
build_reference_db('((((ITS1) OR 5.8S) OR 28S) OR ITS2) AND Embryophyta[Organism]',
                   "ITS_plants.fasta")
build_reference_db('((((ITS1) OR 5.8S) OR 28S) OR ITS2) NOT Embryophyta[Organism]',
                   "ITS_nonplants.fasta", dereplicate = TRUE)
```

## Documentation

Full documentation, including all help pages and guides, is on the
package website: <https://sbresnahan.github.io/metacleaner/>.

* `vignette("metacleaner")`: step-by-step guide, explanation of the
  settings and outputs, and how the old command-line options map onto the R
  arguments.
* `vignette("building-databases")`: downloading reference and contaminant
  databases from NCBI.
* `?clean_reference`, `?classify_queries` (exact decision rules) and
  `?search_sequences` (search engine details).

## Output

| File | Contents |
|------|----------|
| `<prefix>_clean.fasta` | kept sequences |
| `<prefix>_clean.tax` | their taxonomy (`taxid,accession,taxid,domain,...,species`; or QIIME-style with `tax_format = "lineage"`) |
| `<prefix>_flagged.tsv` | removed sequences: query taxonomy, the hit that caused removal, its taxonomy, and the reason |
| `<prefix>_summary.txt` | settings and counts per reason |

## Citation

If you use metacleaner in your research, please cite:

Crone, M., Boyle, N., Bresnahan, S.T., Biddinger, D., Grozinger, C.M. (2023).
More than mesolectic: Characterizing the nutritional niche of *Osmia
cornifrons*. Ecology and Evolution 13, e10640.
<https://doi.org/10.1002/ece3.10640>
