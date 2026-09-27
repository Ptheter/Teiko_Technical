# Teiko Take-Home Technical Assignment

## Overview

This project analyzes immune-cell count data from biological samples using a reproducible Python and SQLite pipeline.

The analysis includes:

* Data validation and cleaning
* Relational database design and SQLite loading
* Relative immune-cell frequency calculations
* Statistical comparison of responders vs. nonresponders
* Baseline melanoma sample summaries
* Interactive visualization with Streamlit and Plotly

The complete workflow can be executed through the provided `Makefile`.

---

## Project Structure

```text
.
├── data/
│   ├── cell-count.csv
│   └── cell-count-cleaned.csv
│
├── output/
│   ├── cleaning_report.txt
│   ├── relative_frequencies.csv
│   ├── part3_data.csv
│   ├── part3_statistics.csv
│   ├── part4_samples_by_project.csv
│   ├── part4_subjects_by_response.csv
│   ├── part4_subjects_by_sex.csv
│   └── part4_average_b_cells.csv
│
├── clean_data.py
├── load_data.py
├── analysis.py
├── dashboard.py
├── schema.sql
├── Makefile
├── requirements.txt
└── README.md
```

### Key files

| File               | Description                                                   |
| ------------------ | ------------------------------------------------------------- |
| `clean_data.py`    | Validates and cleans the input CSV                            |
| `schema.sql`       | Defines the normalized SQLite database schema                 |
| `load_data.py`     | Creates and populates the SQLite database                     |
| `analysis.py`      | Performs Parts 2–4 of the analysis                            |
| `dashboard.py`     | Runs the interactive Streamlit dashboard                      |
| `Makefile`         | Provides reproducible setup, pipeline, and dashboard commands |
| `requirements.txt` | Python dependencies                                           |

---

## Requirements

The project was developed using:

* Python 3
* SQLite
* pandas
* NumPy
* SciPy
* statsmodels
* Plotly
* Streamlit

The project can be run locally or through GitHub Codespaces.

---

## Setup

Clone the repository and open the project in GitHub Codespaces or a local environment.

Install the required Python dependencies:

```bash
make setup
```

---

## Running the Analysis Pipeline

The complete analysis pipeline can be executed with:

```bash
make pipeline
```

This runs the following steps sequentially:

```text
clean_data.py
      ↓
load_data.py
      ↓
analysis.py
```

### Step 1 — Data cleaning

`clean_data.py` validates the input dataset and creates:

```text
data/cell-count-cleaned.csv
```

The cleaning process checks:

* Required columns
* Missing values
* Categorical values
* Numeric values
* Negative cell counts
* Invalid ages and treatment timepoints
* Duplicate sample identifiers
* Completely duplicated rows

Missing response values are retained because response is an outcome variable and should not be imputed or inferred.

A cleaning report is also generated at:

```text
output/cleaning_report.txt
```

### Step 2 — Database loading

`load_data.py` creates:

```text
cell-count.db
```

The database is rebuilt from the cleaned dataset each time the pipeline runs.

The relational schema consists of three tables:

```text
subjects
   │
   └──< samples
           │
           └──< cell_counts
```

#### `subjects`

Contains subject-level information:

* subject
* project
* condition
* age
* sex

#### `samples`

Contains sample-level information:

* sample
* subject
* treatment
* response
* sample_type
* time_from_treatment_start

#### `cell_counts`

Stores cell-population measurements in long format:

* sample
* population
* count

This structure separates subject-level metadata, sample-level metadata, and cell-population measurements while maintaining foreign-key relationships.

### Step 3 — Analysis

`analysis.py` generates the requested Part 2–4 results and saves them under `output/`.

---

# Analysis

## Part 2 — Relative Cell Frequencies

For each sample, the five immune-cell populations are summed to calculate the total cell count:

```text
total_count =
    B cell
    + CD8 T cell
    + CD4 T cell
    + NK cell
    + monocyte
```

The relative frequency of each population is then calculated as:

```text
percentage = population count / total count × 100
```

The resulting dataset is saved to:

```text
output/relative_frequencies.csv
```

The output contains:

```text
sample
total_count
population
count
percentage
```

The analysis produced:

* **10,500 samples**
* **52,500 population-level records**
* **5 immune-cell populations per sample**

---

## Part 3 — Responders vs. Nonresponders

The comparison was restricted to:

* Condition: melanoma
* Sample type: PBMC
* Treatment: miraclib
* Response: yes or no

This produced:

* **1,968 samples**
* **656 subjects**
* **3 timepoints per subject:** 0, 7, and 14
* **993 responder samples**
* **975 nonresponder samples**

Relative frequencies were compared separately for:

* B cells
* CD4 T cells
* CD8 T cells
* NK cells
* Monocytes

### Statistical analysis

A two-sided **Mann–Whitney U test** was used to compare responders and nonresponders for each cell population.

Because five populations were tested, p-values were adjusted using the **Benjamini–Hochberg false discovery rate (FDR)** procedure.

The statistical results are saved to:

```text
output/part3_statistics.csv
```

### Results

| Population | Responder Median (%) | Nonresponder Median (%) | Raw p-value | Adjusted p-value |
| ---------- | -------------------: | ----------------------: | ----------: | ---------------: |
| B cell     |                 9.43 |                    9.79 |      0.0557 |           0.1393 |
| CD4 T cell |                30.22 |                   29.66 |      0.0133 |           0.0667 |
| CD8 T cell |                24.73 |                   24.60 |      0.6391 |           0.6391 |
| Monocyte   |                19.61 |                   19.94 |      0.1632 |           0.2039 |
| NK cell    |                14.51 |                   14.80 |      0.1211 |           0.2018 |

Using an FDR threshold of 0.05, **none of the five populations reached statistical significance after multiple-testing correction**.

The CD4 T-cell comparison had a raw p-value below 0.05, but its adjusted p-value was 0.0667.

### Repeated measurements

The Part 3 dataset contains three measurements for each of the 656 subjects.

The requested comparison was performed at the sample level. This treats samples as the unit of comparison and follows the sample-level framing of the assignment.

For subject-level inference accounting for repeated measurements, a mixed-effects or other repeated-measures statistical model would be more appropriate.

---

## Part 4 — Baseline Melanoma Analysis

Baseline was defined as:

```text
time_from_treatment_start = 0
```

For Part 4A, the analysis was restricted to:

* Condition: melanoma
* Sample type: PBMC
* Treatment: miraclib
* Timepoint: 0

### Samples by project

| Project | Samples |
| ------- | ------: |
| prj1    |     384 |
| prj3    |     272 |

Total baseline samples:

**656**

### Subjects by response

| Response | Subjects |
| -------- | -------: |
| No       |      325 |
| Yes      |      331 |

### Subjects by sex

| Sex    | Subjects |
| ------ | -------: |
| Female |      312 |
| Male   |      344 |

---

## Part 4B — Average B Cells

The requested calculation was performed using:

* Condition: melanoma
* Sex: male
* Response: yes
* Timepoint: 0
* Population: B cell

The calculation includes **all sample and treatment types**, as specified in the assignment.

### Result

**Average B-cell count: 10206.15**

The result is saved to:

```text
output/part4_average_b_cells.csv
```

---

# Interactive Dashboard

The project includes an interactive Streamlit dashboard using Plotly.

Launch it with:

```bash
make dashboard
```

The dashboard provides interactive views of:

* Part 2 relative cell frequencies
* Part 3 responder vs. nonresponder distributions
* Part 3 statistical results
* Part 4 baseline sample summaries
* Part 4 average B-cell calculation

The Part 3 dashboard allows filtering the visualization by cell population and timepoint.

The statistical results displayed in the dashboard correspond to the overall sample-level comparison across all three timepoints.

---

# Reproducibility

The project is designed so that the main workflow does not require manual execution of individual analysis steps.

Run:

```bash
make setup
make pipeline
```

to install dependencies and reproduce the analysis outputs.

Then run:

```bash
make dashboard
```

to launch the interactive dashboard.

To view the dashboard, go to https://teikotechnical-kfayrucj3gzdoa22b5wvwg.streamlit.app/

The database is regenerated from the cleaned CSV during `make pipeline`, ensuring that downstream analysis is based on the current cleaned input data.

---

# Output Files

| Output                           | Description                                 |
| -------------------------------- | ------------------------------------------- |
| `cleaning_report.txt`            | Data validation and cleaning summary        |
| `relative_frequencies.csv`       | Part 2 cell-population relative frequencies |
| `part3_data.csv`                 | Part 3 analysis dataset                     |
| `part3_statistics.csv`           | Mann–Whitney U and FDR results              |
| `part4_samples_by_project.csv`   | Baseline samples by project                 |
| `part4_subjects_by_response.csv` | Baseline subjects by response               |
| `part4_subjects_by_sex.csv`      | Baseline subjects by sex                    |
| `part4_average_b_cells.csv`      | Part 4B average B-cell count                |

---

# Notes

The analysis preserves missing response values in the cleaned dataset rather than imputing them. Response-specific analyses exclude samples where response is unavailable.

Statistical conclusions are based on the specified sample-level analysis and should be interpreted in the context of repeated measurements from the same subjects.