from pathlib import Path
import sqlite3
import os

import pandas as pd


# -------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------

BASE_DIR = Path(__file__).resolve().parent

INPUT_FILE = BASE_DIR / "data" / "cell-count-cleaned.csv"
SCHEMA_FILE = BASE_DIR / "schema.sql"
DATABASE_FILE = Path(
    os.environ.get(
        "TEIKO_DB_PATH",
        BASE_DIR / "cell-count.db",
    )
)

CELL_COUNT_COLUMNS = [
    "b_cell",
    "cd8_t_cell",
    "cd4_t_cell",
    "nk_cell",
    "monocyte",
]


# -------------------------------------------------------------------
# Database setup
# -------------------------------------------------------------------

def create_database(connection):
    """Create database tables using schema.sql."""

    schema = SCHEMA_FILE.read_text(encoding="utf-8")

    connection.executescript(schema)

    connection.commit()


# -------------------------------------------------------------------
# Load subjects
# -------------------------------------------------------------------

def load_subjects(connection, df):
    """Load unique subject-level information."""

    subjects = (
        df[
            [
                "subject",
                "project",
                "condition",
                "age",
                "sex",
            ]
        ]
        .drop_duplicates(subset=["subject"])
    )

    records = subjects.itertuples(index=False, name=None)

    connection.executemany(
        """
        INSERT INTO subjects (
            subject,
            project,
            condition,
            age,
            sex
        )
        VALUES (?, ?, ?, ?, ?)
        """,
        records,
    )


# -------------------------------------------------------------------
# Load samples
# -------------------------------------------------------------------

def load_samples(connection, df):
    """Load sample-level metadata."""

    samples = df[
        [
            "sample",
            "subject",
            "treatment",
            "response",
            "sample_type",
            "time_from_treatment_start",
        ]
    ].drop_duplicates(subset=["sample"])

    records = samples.itertuples(index=False, name=None)

    connection.executemany(
        """
        INSERT INTO samples (
            sample,
            subject,
            treatment,
            response,
            sample_type,
            time_from_treatment_start
        )
        VALUES (?, ?, ?, ?, ?, ?)
        """,
        records,
    )


# -------------------------------------------------------------------
# Load cell counts
# -------------------------------------------------------------------

def load_cell_counts(connection, df):
    """Convert wide cell-count columns to long format and load them."""

    cell_counts = df[
        ["sample"] + CELL_COUNT_COLUMNS
    ].melt(
        id_vars="sample",
        value_vars=CELL_COUNT_COLUMNS,
        var_name="population",
        value_name="count",
    )

    records = cell_counts.itertuples(index=False, name=None)

    connection.executemany(
        """
        INSERT INTO cell_counts (
            sample,
            population,
            count
        )
        VALUES (?, ?, ?)
        """,
        records,
    )


# -------------------------------------------------------------------
# Verification
# -------------------------------------------------------------------

def verify_database(connection, expected_samples):
    """Verify that the database contains the expected number of records."""

    subject_count = connection.execute(
        "SELECT COUNT(*) FROM subjects"
    ).fetchone()[0]

    sample_count = connection.execute(
        "SELECT COUNT(*) FROM samples"
    ).fetchone()[0]

    cell_count = connection.execute(
        "SELECT COUNT(*) FROM cell_counts"
    ).fetchone()[0]

    expected_cell_count = expected_samples * len(CELL_COUNT_COLUMNS)

    print("\nDatabase verification")
    print("---------------------")
    print(f"Subjects:     {subject_count}")
    print(f"Samples:      {sample_count}")
    print(f"Cell counts:  {cell_count}")

    if sample_count != expected_samples:
        raise ValueError(
            f"Expected {expected_samples} samples, "
            f"but loaded {sample_count}."
        )

    if cell_count != expected_cell_count:
        raise ValueError(
            f"Expected {expected_cell_count} cell-count records, "
            f"but loaded {cell_count}."
        )

    print("\nDatabase loaded successfully.")


# -------------------------------------------------------------------
# Main
# -------------------------------------------------------------------

def main():

    # ---------------------------------------------------------------
    # Check input files
    # ---------------------------------------------------------------

    if not INPUT_FILE.exists():
        raise FileNotFoundError(
            f"Could not find cleaned dataset:\n{INPUT_FILE}\n\n"
            "Run clean_data.py first."
        )

    if not SCHEMA_FILE.exists():
        raise FileNotFoundError(
            f"Could not find schema file:\n{SCHEMA_FILE}"
        )

    # ---------------------------------------------------------------
    # Read cleaned dataset
    # ---------------------------------------------------------------

    print(f"Reading: {INPUT_FILE}")

    df = pd.read_csv(INPUT_FILE)

    print(f"Rows loaded from CSV: {len(df)}")

    # ---------------------------------------------------------------
    # Validate required columns
    # ---------------------------------------------------------------

    required_columns = [
        "project",
        "subject",
        "condition",
        "age",
        "sex",
        "treatment",
        "response",
        "sample",
        "sample_type",
        "time_from_treatment_start",
        *CELL_COUNT_COLUMNS,
    ]

    missing_columns = [
        column
        for column in required_columns
        if column not in df.columns
    ]

    if missing_columns:
        raise ValueError(
            f"Missing required columns: {missing_columns}"
        )

    # ---------------------------------------------------------------
    # Convert missing response values to Python None
    # ---------------------------------------------------------------

    # SQLite stores Python None as SQL NULL.
    df["response"] = df["response"].where(
        df["response"].notna(),
        None,
    )

    # ---------------------------------------------------------------
    # Create database
    # ---------------------------------------------------------------

    # Delete the existing database so every run starts clean.
    if DATABASE_FILE.exists():
        DATABASE_FILE.unlink()

    print(f"Creating database: {DATABASE_FILE}")

    connection = sqlite3.connect(DATABASE_FILE)

    try:
        # Enable foreign-key enforcement.
        connection.execute("PRAGMA foreign_keys = ON")

        # Create tables.
        create_database(connection)

        # Load data.
        print("Loading subjects...")
        load_subjects(connection, df)

        print("Loading samples...")
        load_samples(connection, df)

        print("Loading cell counts...")
        load_cell_counts(connection, df)

        # Commit all inserts.
        connection.commit()

        # Verify.
        verify_database(
            connection,
            expected_samples=len(df),
        )

    except Exception:
        connection.rollback()
        raise

    finally:
        connection.close()

    print(f"\nDatabase created at:")
    print(DATABASE_FILE)


if __name__ == "__main__":
    main()