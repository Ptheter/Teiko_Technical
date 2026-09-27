from pathlib import Path
import sys
import sqlite3
import subprocess

import pandas as pd
import plotly.express as px
import streamlit as st


# ---------------------------------------------------------
# Database setup
# ---------------------------------------------------------

BASE_DIR = Path(__file__).resolve().parent
DB_PATH = BASE_DIR / "cell-count.db"
OUTPUT_DIR = BASE_DIR / "output"


def database_is_ready():
    if not DB_PATH.exists():
        return False

    connection = sqlite3.connect(DB_PATH)
    try:
        tables = pd.read_sql_query(
            """
            SELECT name
            FROM sqlite_master
            WHERE type = 'table'
              AND name IN ('subjects', 'samples', 'cell_counts')
            """,
            connection,
        )
        return set(tables["name"]) == {"subjects", "samples", "cell_counts"}
    finally:
        connection.close()


if not database_is_ready():
    if DB_PATH.exists():
        DB_PATH.unlink()

    subprocess.run(
        [sys.executable, str(BASE_DIR / "load_data.py")],
        check=True,
        cwd=BASE_DIR,
    )


@st.cache_resource
def get_connection():
    return sqlite3.connect(DB_PATH, check_same_thread=False)


conn = get_connection()


# ---------------------------------------------------------
# Database queries
# ---------------------------------------------------------

@st.cache_data
def load_relative_frequencies():
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

    return pd.read_sql_query(query, conn)


@st.cache_data
def load_part3_data():
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

    return pd.read_sql_query(query, conn)


@st.cache_data
def load_part3_statistics():
    return pd.read_csv(OUTPUT_DIR / "part3_statistics.csv")


@st.cache_data
def load_part4_project():
    query = """
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

    return pd.read_sql_query(query, conn)


@st.cache_data
def load_part4_response():
    query = """
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

    return pd.read_sql_query(query, conn)


@st.cache_data
def load_part4_sex():
    query = """
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

    return pd.read_sql_query(query, conn)


@st.cache_data
def load_part4_average():
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

    return pd.read_sql_query(query, conn)


relative_df = load_relative_frequencies()
part3_df = load_part3_data()
part3_stats = load_part3_statistics()
project_df = load_part4_project()
response_df = load_part4_response()
sex_df = load_part4_sex()
average_df = load_part4_average()

# ---------------------------------------------------------
# Header
# ---------------------------------------------------------

st.title("🧬 Cell Count Analysis Dashboard")

st.markdown(
    """
    Interactive analysis of immune-cell composition across biological samples.

    The dashboard presents results from Parts 2–4 of the analysis pipeline,
    using the SQLite database and generated analysis outputs.
    """
)


# ---------------------------------------------------------
# Dataset overview
# ---------------------------------------------------------

st.header("Dataset Overview")

col1, col2, col3 = st.columns(3)

with col1:
    st.metric("Samples", f"{relative_df['sample'].nunique():,}")

with col2:
    st.metric("Cell populations", f"{relative_df['population'].nunique()}")

with col3:
    st.metric("Relative-frequency rows", f"{len(relative_df):,}")


# ---------------------------------------------------------
# Part 2
# ---------------------------------------------------------

st.header("Part 2 — Relative Cell Frequencies")

st.markdown(
    """
    For each sample, total cell count is calculated across the five immune
    populations. Each population is then expressed as a percentage of the
    sample total.
    """
)

sample_options = sorted(relative_df["sample"].unique())

selected_sample = st.selectbox(
    "Select a sample",
    sample_options,
)

sample_df = relative_df[
    relative_df["sample"] == selected_sample
].copy()

total_count = int(sample_df["total_count"].iloc[0])

st.metric(
    "Total cell count",
    f"{total_count:,}",
)

fig_part2 = px.bar(
    sample_df,
    x="population",
    y="percentage",
    text="percentage",
    labels={
        "population": "Cell population",
        "percentage": "Relative frequency (%)",
    },
    title=f"Cell Population Composition — {selected_sample}",
)

fig_part2.update_traces(
    texttemplate="%{text:.2f}%",
    textposition="outside",
)

fig_part2.update_layout(
    yaxis_range=[0, max(sample_df["percentage"]) * 1.2],
)

st.plotly_chart(
    fig_part2,
    use_container_width=True,
)


with st.expander("View Part 2 data"):
    display_df = sample_df.copy()
    display_df["count"] = display_df["count"].astype(int)
    display_df["total_count"] = display_df["total_count"].astype(int)
    display_df["percentage"] = display_df["percentage"].round(2)

    st.dataframe(
        display_df[
            [
                "sample",
                "total_count",
                "population",
                "count",
                "percentage",
            ]
        ],
        use_container_width=True,
        hide_index=True,
    )


# ---------------------------------------------------------
# Part 3
# ---------------------------------------------------------

st.header("Part 3 — Responders vs Nonresponders")

st.markdown(
    """
    Relative immune-cell frequencies are compared between responders and
    nonresponders among melanoma PBMC samples receiving miraclib.

    Statistical comparisons use a two-sided Mann–Whitney U test with
    Benjamini–Hochberg false-discovery-rate correction.
    """
)

col1, col2 = st.columns(2)

with col1:
    population_options = sorted(part3_df["population"].unique())

    selected_population = st.selectbox(
        "Cell population",
        population_options,
    )

with col2:
    timepoint_options = ["All"] + sorted(
        part3_df["time_from_treatment_start"].unique().tolist()
    )

    selected_timepoint = st.selectbox(
        "Timepoint",
        timepoint_options,
    )


plot_df = part3_df[
    part3_df["population"] == selected_population
].copy()

if selected_timepoint != "All":
    plot_df = plot_df[
        plot_df["time_from_treatment_start"] == selected_timepoint
    ]


fig_part3 = px.box(
    plot_df,
    x="response",
    y="percentage",
    color="response",
    points="outliers",
    labels={
        "response": "Response",
        "percentage": "Relative frequency (%)",
    },
    title=(
        f"{selected_population} Relative Frequency — "
        f"{selected_timepoint if selected_timepoint != 'All' else 'All Timepoints'}"
    ),
)

fig_part3.update_layout(
    showlegend=False,
)

st.plotly_chart(
    fig_part3,
    use_container_width=True,
)


# ---------------------------------------------------------
# Part 3 statistics
# ---------------------------------------------------------

selected_stats = part3_stats[
    part3_stats["population"] == selected_population
].iloc[0]

col1, col2, col3, col4 = st.columns(4)

with col1:
    st.metric(
        "Responder median",
        f"{selected_stats['responder_median']:.2f}%",
    )

with col2:
    st.metric(
        "Nonresponder median",
        f"{selected_stats['nonresponder_median']:.2f}%",
    )

with col3:
    st.metric(
        "Raw p-value",
        f"{selected_stats['p_value']:.4f}",
    )

with col4:
    st.metric(
        "Adjusted p-value",
        f"{selected_stats['adjusted_p_value']:.4f}",
    )

if selected_stats["significant_fdr_0.05"]:
    st.success(
        "This population shows a statistically significant difference "
        "after Benjamini–Hochberg FDR correction (adjusted p < 0.05)."
    )
else:
    st.info(
        "No statistically significant difference was detected after "
        "Benjamini–Hochberg FDR correction (adjusted p ≥ 0.05)."
    )

st.caption(
    "Statistical results shown above are based on all timepoints and the "
    "sample-level comparison used in the analysis pipeline."
)

with st.expander("View statistical results for all populations"):
    stats_display = part3_stats.copy()

    stats_display["responder_median"] = (
        stats_display["responder_median"].round(2)
    )

    stats_display["nonresponder_median"] = (
        stats_display["nonresponder_median"].round(2)
    )

    stats_display["p_value"] = (
        stats_display["p_value"].round(4)
    )

    stats_display["adjusted_p_value"] = (
        stats_display["adjusted_p_value"].round(4)
    )

    st.dataframe(
        stats_display,
        use_container_width=True,
        hide_index=True,
    )


# ---------------------------------------------------------
# Part 4
# ---------------------------------------------------------

st.header("Part 4 — Baseline Melanoma Analysis")

st.markdown(
    """
    Baseline melanoma PBMC samples receiving miraclib are summarized by
    project, response status, and sex.
    """
)


col1, col2 = st.columns(2)

with col1:
    st.subheader("Samples by Project")

    st.dataframe(
        project_df,
        use_container_width=True,
        hide_index=True,
    )

with col2:
    st.subheader("Subjects by Response")

    response_display = response_df.copy()
    response_display["response"] = (
        response_display["response"].str.capitalize()
    )

    st.dataframe(
        response_display,
        use_container_width=True,
        hide_index=True,
    )


st.subheader("Subjects by Sex")

sex_display = sex_df.copy()
sex_display["sex"] = sex_display["sex"].str.upper()

st.dataframe(
    sex_display,
    use_container_width=True,
    hide_index=True,
)


# ---------------------------------------------------------
# Part 4B
# ---------------------------------------------------------

st.subheader("Average B Cells")

average_b_cells = float(
    average_df["average_b_cells"].iloc[0]
)

st.metric(
    "Male melanoma responders at time 0",
    f"{average_b_cells:,.2f}",
)

st.caption(
    "Calculated for melanoma males with response = yes at time 0, "
    "across all sample and treatment types."
)


# ---------------------------------------------------------
# Footer
# ---------------------------------------------------------

st.divider()