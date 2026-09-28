from pathlib import Path
import sqlite3
import csv
import pandas as pd
import plotly.express as px
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


# ---------------------------------------------------------
# Paths
# ---------------------------------------------------------

BASE_DIR = Path(__file__).resolve().parent
DB_PATH = BASE_DIR / "cell-count.db"
OUTPUT_DIR = BASE_DIR / "output"
OUTPUT_PATH = OUTPUT_DIR / "relative_frequencies.csv"


# ---------------------------------------------------------
# Part 2: Calculate relative cell frequencies
# ---------------------------------------------------------

def calculate_relative_frequencies(conn):
    """
    Calculate total cell count and relative frequency for
    each cell population within every sample.
    """

    query = """
        WITH sample_totals AS (
            SELECT
                sample,
                SUM(count) AS total_count
            FROM cell_counts
            GROUP BY sample
        )
        SELECT
            cc.sample,
            st.total_count,
            cc.population,
            cc.count,
            (cc.count * 100.0 / st.total_count) AS percentage
        FROM cell_counts AS cc
        JOIN sample_totals AS st
            ON cc.sample = st.sample
        ORDER BY
            cc.sample,
            cc.population;
    """

    cursor = conn.execute(query)
    rows = cursor.fetchall()

    return rows


# ---------------------------------------------------------
# Part 3: Responders vs nonresponders
# ---------------------------------------------------------

def analyze_responders_vs_nonresponders(conn):
    """
    Prepare the sample-level dataset for comparing responders and nonresponders 
    among melanoma PBMC samples treated with miraclib. Each subject contributes 
    samples at timepoints 0, 7, and 14.
    """

    query = """
        WITH sample_totals AS (
            SELECT
                sample,
                SUM(count) AS total_count
            FROM cell_counts
            GROUP BY sample
        )
        SELECT
            sm.sample,
            sm.subject,
            sm.time_from_treatment_start,
            sm.response,
            cc.population,
            cc.count,
            st.total_count,
            (cc.count * 100.0 / st.total_count) AS percentage
        FROM samples AS sm
        JOIN cell_counts AS cc
            ON sm.sample = cc.sample
        JOIN sample_totals AS st
            ON sm.sample = st.sample
        JOIN subjects AS s
            ON sm.subject = s.subject
        WHERE LOWER(s.condition) = 'melanoma'
          AND LOWER(sm.sample_type) = 'pbmc'
          AND LOWER(sm.treatment) = 'miraclib'
          AND sm.response IN ('yes', 'no')
        ORDER BY
            sm.time_from_treatment_start,
            cc.population,
            sm.response,
            sm.sample;
    """

    df = pd.read_sql_query(query, conn)

    if df.empty:
        raise ValueError("No melanoma PBMC miraclib samples found.")

    print("\nPart 3 dataset")
    print("----------------")
    print(f"Samples: {df['sample'].nunique()}")
    print(f"Subjects: {df['subject'].nunique()}")
    print(f"Rows: {len(df)}")

    print("\nSamples by response:")
    print(
        df[['sample', 'response']]
        .drop_duplicates()
        ['response']
        .value_counts()
    )

    print("\nSamples by timepoint:") 
    print( 
        df[["sample", "time_from_treatment_start"]].drop_duplicates() 
        ["time_from_treatment_start"] 
        .value_counts() 
        .sort_index() )

    return df


# ---------------------------------------------------------
# Part 4: Baseline melanoma PBMC + miraclib analysis
# ---------------------------------------------------------

def analyze_part4_baseline(conn):
    """
    Analyze baseline (time = 0) melanoma PBMC samples
    treated with miraclib.

    Returns:
        project-level sample counts
        responder/nonresponder subject counts
        male/female subject counts
    """

    # ---------------------------------------------
    # 4A. Number of samples per project
    # ---------------------------------------------

    project_query = """
        SELECT
            s.project,
            COUNT(*) AS sample_count
        FROM samples AS sm
        JOIN subjects AS s
            ON sm.subject = s.subject
        WHERE LOWER(s.condition) = 'melanoma'
          AND LOWER(sm.sample_type) = 'pbmc'
          AND LOWER(sm.treatment) = 'miraclib'
          AND sm.time_from_treatment_start = 0
        GROUP BY s.project
        ORDER BY s.project;
    """

    project_df = pd.read_sql_query(
        project_query,
        conn
    )

    # ---------------------------------------------
    # 4A. Responder / nonresponder subjects
    # ---------------------------------------------

    response_query = """
        SELECT
            sm.response,
            COUNT(DISTINCT sm.subject) AS subject_count
        FROM samples AS sm
        JOIN subjects AS s
            ON sm.subject = s.subject
        WHERE LOWER(s.condition) = 'melanoma'
          AND LOWER(sm.sample_type) = 'pbmc'
          AND LOWER(sm.treatment) = 'miraclib'
          AND sm.time_from_treatment_start = 0
          AND sm.response IN ('yes', 'no')
        GROUP BY sm.response
        ORDER BY sm.response;
    """

    response_df = pd.read_sql_query(
        response_query,
        conn
    )

    # ---------------------------------------------
    # 4A. Male / female subjects
    # ---------------------------------------------

    sex_query = """
        SELECT
            s.sex,
            COUNT(DISTINCT s.subject) AS subject_count
        FROM samples AS sm
        JOIN subjects AS s
            ON sm.subject = s.subject
        WHERE LOWER(s.condition) = 'melanoma'
          AND LOWER(sm.sample_type) = 'pbmc'
          AND LOWER(sm.treatment) = 'miraclib'
          AND sm.time_from_treatment_start = 0
        GROUP BY s.sex
        ORDER BY s.sex;
    """

    sex_df = pd.read_sql_query(
        sex_query,
        conn
    )

    # ---------------------------------------------
    # Print results
    # ---------------------------------------------

    print("\nPart 4A: Baseline melanoma PBMC + miraclib")
    print("--------------------------------------------")

    print("\nSamples per project:")
    print(project_df.to_string(index=False))

    print("\nSubjects by response:")
    print(response_df.to_string(index=False))

    print("\nSubjects by sex:")
    print(sex_df.to_string(index=False))

    return project_df, response_df, sex_df


def calculate_average_b_cells(conn):
    """
    Calculate the average number of B cells for
    male melanoma responders at time = 0.

    The assignment specifies all sample and treatment
    types, so treatment and sample_type are NOT filtered.
    """

    query = """
        SELECT
            ROUND(AVG(cc.count), 2) AS average_b_cells
        FROM cell_counts AS cc
        JOIN samples AS sm
            ON cc.sample = sm.sample
        JOIN subjects AS s
            ON sm.subject = s.subject
        WHERE LOWER(s.condition) = 'melanoma'
          AND LOWER(s.sex) = 'm'
          AND LOWER(sm.response) = 'yes'
          AND sm.time_from_treatment_start = 0
          AND cc.population = 'b_cell';
    """

    result = pd.read_sql_query(
        query,
        conn
    )

    average_b_cells = result.loc[0, "average_b_cells"]

    print("\nPart 4B: Average B cells")
    print("------------------------")
    print(
        f"Average B cells for male melanoma responders "
        f"at time 0: {average_b_cells:.2f}"
    )

    return average_b_cells


def save_part4_results(
    project_df,
    response_df,
    sex_df,
    average_b_cells
):
    """
    Save Part 4 results for reproducibility and dashboard use.
    """

    project_path = OUTPUT_DIR / "part4_samples_by_project.csv"
    response_path = OUTPUT_DIR / "part4_subjects_by_response.csv"
    sex_path = OUTPUT_DIR / "part4_subjects_by_sex.csv"
    average_path = OUTPUT_DIR / "part4_average_b_cells.csv"

    project_df.to_csv(
        project_path,
        index=False
    )

    response_df.to_csv(
        response_path,
        index=False
    )

    sex_df.to_csv(
        sex_path,
        index=False
    )

    average_df = pd.DataFrame({
        "average_b_cells": [average_b_cells]
    })

    average_df.to_csv(
        average_path,
        index=False,
        float_format="%.2f"
    )

    print(f"\nSaved: {project_path}")
    print(f"Saved: {response_path}")
    print(f"Saved: {sex_path}")
    print(f"Saved: {average_path}")


def run_statistical_tests(df):
    """
    Run Mann-Whitney U tests comparing responders and
    nonresponders for each cell population.

    Benjamini-Hochberg correction is applied across
    the five population-level tests.
    """

    results = []

    populations = sorted(df["population"].unique())

    for population in populations:

        responder_values = df[
            (df["population"] == population)
            & (df["response"] == "yes")
        ]["percentage"]

        nonresponder_values = df[
            (df["population"] == population)
            & (df["response"] == "no")
        ]["percentage"]

        statistic, p_value = mannwhitneyu(
            responder_values,
            nonresponder_values,
            alternative="two-sided"
        )

        results.append({
            "population": population,
            "n_responder": len(responder_values),
            "n_nonresponder": len(nonresponder_values),
            "responder_median": responder_values.median(),
            "nonresponder_median": nonresponder_values.median(),
            "u_statistic": statistic,
            "p_value": p_value
        })

    results_df = pd.DataFrame(results)

    # Benjamini-Hochberg false discovery rate correction
    rejected, adjusted_p_values, _, _ = multipletests(
        results_df["p_value"],
        method="fdr_bh"
    )

    results_df["adjusted_p_value"] = adjusted_p_values
    results_df["significant_fdr_0.05"] = rejected

    return results_df

'''
def make_part3_boxplots(df):
    """
    Create interactive boxplots for each immune-cell population.
    """

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    populations = sorted(df["population"].unique())

    for population in populations:

        population_df = df[
            df["population"] == population
        ].copy()

        fig = px.box(
            population_df,
            x="response",
            y="percentage",
            points="all",
            category_orders={
                "response": ["yes", "no"]
            },
            labels={
                "response": "Response",
                "percentage": "Relative frequency (%)"
            },
            title=f"{population}: Melanoma PBMC + Miraclib"
        )

        fig.update_xaxes(
            ticktext=["Responder", "Nonresponder"],
            tickvals=["yes", "no"]
        )

        fig.update_layout(
            width=800,
            height=600
        )

        output_path = OUTPUT_DIR / f"part3_{population}_boxplot.html"

        fig.write_html(output_path)

        print(f"Saved: {output_path}")
'''

def save_part3_dataset(df): 
    """Save the Part 3 sample-level dataset.""" 
    output_path = OUTPUT_DIR / "part3_data.csv" 
    df.to_csv( output_path, index=False, float_format="%.6f" ) 
    print(f"Saved: {output_path}")


def save_part3_results(results_df):
    """Save statistical test results."""

    output_path = OUTPUT_DIR / "part3_statistics.csv"

    results_df.to_csv(
        output_path,
        index=False,
        float_format="%.6g"
    )

    print(f"Saved: {output_path}")


# ---------------------------------------------------------
# Save results
# ---------------------------------------------------------

def save_results(rows):
    """Save relative-frequency results to a CSV file."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    with open(OUTPUT_PATH, "w", newline="") as file:
        writer = csv.writer(file)

        writer.writerow([
            "sample",
            "total_count",
            "population",
            "count",
            "percentage"
        ])

        for row in rows:
            writer.writerow([
                row[0],
                int(row[1]),
                row[2],
                int(row[3]),
                f"{row[4]:.2f}%"
            ])


# ---------------------------------------------------------
# Main
# ---------------------------------------------------------

def main():
    if not DB_PATH.exists():
        raise FileNotFoundError(
            f"Database not found: {DB_PATH}\n"
            "Run 'python3 load_data.py' first."
        )

    print(f"Reading database: {DB_PATH}")

    conn = sqlite3.connect(DB_PATH)

    try:
        rows = calculate_relative_frequencies(conn)

        print(f"Relative-frequency rows calculated: {len(rows)}")

        # Each sample should have 5 populations.
        expected_rows = 10500 * 5

        if len(rows) != expected_rows:
            raise ValueError(
                f"Expected {expected_rows} rows, but found {len(rows)}."
            )

        save_results(rows)

        print()
        print(f"Output saved to: {OUTPUT_PATH}")

        # Part 3
        df_part3 = analyze_responders_vs_nonresponders(conn)

        results_df = run_statistical_tests(df_part3)

        print("\nStatistical results")
        print("-------------------")
        print(results_df.to_string(index=False))

        save_part3_dataset(df_part3)

        save_part3_results(results_df)

        # Part 4
        project_df, response_df, sex_df = analyze_part4_baseline(conn)

        average_b_cells = calculate_average_b_cells(conn)

        save_part4_results(
            project_df,
            response_df,
            sex_df,
            average_b_cells
        )

    
    finally:
        conn.close()


if __name__ == "__main__":
    main()