[README.md](https://github.com/user-attachments/files/33468189/README.md)
# CelluloScope

R script for filtering cellulosome-associated domain hits and annotating domain overlaps in prepared HMMER results.

## Overview

CelluloScope groups domain hits by protein, filters proteins using predefined domain-combination rules, and annotates overlapping hits within each domain class. Results are exported to Excel for further inspection.

The rules are specific to this analysis and can be edited for other datasets. Retained proteins are candidates for further evaluation, not experimentally validated cellulosome components.

## Requirements

Install R and the packages loaded by the script:

```r
install.packages(c("readxl", "dplyr", "stringr", "writexl", "tidyr", "purrr", "readr"))
```

Use a version of `dplyr` that supports `relationship = "many-to-many"` in `inner_join()` (1.1.1 or later).

## Input

The script reads **`hmmer_minor_results+1.xlsx`** from the current R working directory. It reads the first worksheet with a header row and expects one row per domain hit.

**This is a prepared Excel table, not a raw HMMER `.domtblout` file.** Raw HMMER output must be converted and arranged into the expected format before running the script.

The table must contain **exactly 28 columns in the order below**. The script replaces the original column names by position, so column order matters even if the input headers have different names.

| Position | Assigned column name |
| --- | --- |
| 1 | `target_name` |
| 2 | `accession` |
| 3 | `acc` |
| 4 | `acc_clean` |
| 5 | `name` |
| 6 | `tlen` |
| 7 | `query_name` |
| 8 | `acce2` |
| 9 | `qlen` |
| 10 | `Evalue_full` |
| 11 | `score_full` |
| 12 | `bias_full` |
| 13 | `instance` |
| 14 | `of` |
| 15 | `cEvalue_dom` |
| 16 | `iEvalue_dom` |
| 17 | `score_dom` |
| 18 | `bias_dom` |
| 19 | `hmm_from` |
| 20 | `hmm_to` |
| 21 | `ali_from` |
| 22 | `ali_to` |
| 23 | `env_from` |
| 24 | `env_to` |
| 25 | `acc2` |
| 26 | `description_target` |
| 27 | `DESCRIPTION_UNIPROT` |
| 28 | `Domain` |

Key fields used in the analysis:

- `query_name`: protein identifier used to group hits.
- `name`: domain label used to assign a domain class; matching is case-sensitive.
- `target_name`: hit/model identifier used in overlap annotations.
- `ali_from` and `ali_to`: alignment start and end coordinates; both must contain valid numeric values without missing entries.

Other fields are retained in the exported table. In this script, classification uses `name`, not `Domain`.

## Domain classification

| Input value in `name` | Assigned `domain_class` |
| --- | --- |
| `cohesin` | `cohesin` |
| `dockerin` | `dockerin` |
| `cbm` | `cbm` |
| `slh` | `slh` |
| `alt_anchor` | `alt_anchor` |
| `gh`, `ce`, `pl` | `enzyme` |
| `rdgl` | `regulator` |
| `x` | `x_domain` |
| Any other value | `other` |

## How to run

1. Place `celluloscope.R` and the input Excel file in the same directory.
2. Set that directory as your R working directory, for example through an RStudio project.
3. Run the full script:

```r
source("celluloscope.R")
```

Alternatively, from a terminal opened in that directory:

```bash
Rscript celluloscope.R
```

To use another input filename, edit this line in the script:

```r
df <- readxl::read_excel("hmmer_minor_results+1.xlsx", col_names = TRUE)
```

The script has no command-line arguments. Input and output filenames are configured directly in the code.

## Filtering logic

For each protein, the script builds a set of unique domain classes and compares it with `remove_patterns` and `keep_patterns` using exact set matching. Domain order and repeated copies do not affect the match.

A protein is retained when its combination is absent from `remove_patterns`, or when it also matches `keep_patterns`. Thus, **keep rules take precedence**, and combinations absent from both lists are retained. Filtering removes all rows belonging to an excluded protein.

For example, an `enzyme`-only protein is removed, while an `enzyme` + `dockerin` protein is retained. The `regulator`-only combination occurs in both lists and is retained because the keep rule overrides removal.

## Overlap annotation

Hits are compared only when they share both `query_name` and `domain_class`, using `ali_from` and `ali_to` coordinates.

| Annotation | Meaning |
| --- | --- |
| `partial_overlap` | Intervals overlap without either fully containing the other. |
| `contained_in` | The current interval lies entirely within another interval. |
| `contains_other` | The current interval fully contains another interval. |
| `none` | No other hit exists in the same protein/domain-class group. |

Annotations include the related hit identifier and coordinates, for example `contained_in:model_B(20-150)`. Multiple relationships are separated by semicolons; multiple hits within a relationship are separated by commas.

Current implementation details:

- Hits with identical intervals receive both containment labels.
- Partial overlap uses strict inequalities, so intervals that only touch at an endpoint are not labeled as partially overlapping.
- If other hits exist in the group but none overlaps, the output is an empty string rather than `none`.
- Overlapping hits are annotated, not removed or ranked.

## Output

The script creates the following file in the working directory:

```text
Filtered_Domains_With_Overlap+Evalue_RF1.xlsx
```

It contains retained input rows with the assigned column names and three additional columns:

- `domain_class`: domain category assigned from `name`.
- `id`: local hit number within each protein/domain-class group, ordered by alignment start.
- `overlap`: overlap relationship annotations.

The console also reports row counts at loading, classification, and filtering stages.

**Despite the output filename, this version does not filter by E-value.** E-value columns are retained as input metadata. Change the filename in `write_xlsx()` and the final console message if desired.

## Limitations

- The input schema and filtering rules are fixed in the script and must be adapted manually for different table structures or biological criteria.
- The script does not run HMMER, parse raw HMMER output, or validate protein function.
- Overlaps between different domain classes are not evaluated.
- Pairwise comparisons can consume substantial memory for proteins with many hits in the same domain class.
- This README documents the supplied code; execution requires a compatible input Excel file.
