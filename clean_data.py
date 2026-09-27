from pathlib import Path
import pandas as pd


# -------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------

INPUT_FILE = Path("data/cell-count.csv")
OUTPUT_FILE = Path("data/cell-count-cleaned.csv")
REPORT_FILE = Path("output/cleaning_report.txt")

REQUIRED_COLUMNS = [
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
    "b_cell",
    "cd8_t_cell",
    "cd4_t_cell",
    "nk_cell",
    "monocyte",
]

CELL_COUNT_COLUMNS = [
    "b_cell",
    "cd8_t_cell",
    "cd4_t_cell",
    "nk_cell",
    "monocyte",
]

# Values that are explicitly allowed to be missing.
# Response is allowed to be missing because it is an outcome variable.
OPTIONAL_COLUMNS = {"response"}


# -------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------

def add_report(report, message=""):
    report.append(message)


def clean_string_columns(df):
    """Strip leading/trailing whitespace from string columns."""
    string_columns = df.select_dtypes(include=["object", "string"]).columns

    for column in string_columns:
        df[column] = df[column].astype("string").str.strip()

    return df


def normalize_categorical_values(df):
    """
    Normalize categorical values where case should not matter.

    This does NOT invent values. It only makes common formatting
    differences consistent.
    """
    for column in ["sex", "response"]:
        if column in df.columns:
            df[column] = df[column].str.lower()

    return df


def validate_required_columns(df):
    missing_columns = [
        column for column in REQUIRED_COLUMNS
        if column not in df.columns
    ]

    if missing_columns:
        raise ValueError(
            f"Missing required columns: {missing_columns}"
        )


def validate_missing_values(df, report):
    add_report(report, "MISSING VALUES")
    add_report(report, "--------------")

    for column in df.columns:
        missing_count = df[column].isna().sum()

        if missing_count > 0:
            add_report(
                report,
                f"{column}: {missing_count}"
            )

    # Check required fields that should not be missing.
    required_non_null = [
        column
        for column in REQUIRED_COLUMNS
        if column not in OPTIONAL_COLUMNS
    ]

    invalid_required = {}

    for column in required_non_null:
        count = df[column].isna().sum()

        if count > 0:
            invalid_required[column] = count

    if invalid_required:
        raise ValueError(
            "Missing values found in required columns: "
            f"{invalid_required}"
        )

    add_report(report)


def validate_categorical_values(df, report):
    add_report(report, "CATEGORICAL VALUES")
    add_report(report, "------------------")

    # Values are case-normalized above.
    expected_values = {
        "sex": {"m", "f"},
        "response": {"yes", "no"},
    }

    for column, allowed_values in expected_values.items():

        # Exclude missing values from unexpected-value checking.
        observed = set(
            df[column]
            .dropna()
            .astype(str)
            .str.lower()
        )

        unexpected = observed - allowed_values

        if unexpected:
            add_report(
                report,
                f"{column}: unexpected values = {sorted(unexpected)}"
            )
        else:
            add_report(
                report,
                f"{column}: no unexpected values"
            )

    add_report(report)

    # Fail the pipeline if categorical values are invalid.
    problems = []

    for column, allowed_values in expected_values.items():

        observed = set(
            df[column]
            .dropna()
            .astype(str)
            .str.lower()
        )

        unexpected = observed - allowed_values

        if unexpected:
            problems.append(
                f"{column}: {sorted(unexpected)}"
            )

    if problems:
        raise ValueError(
            "Unexpected categorical values found:\n"
            + "\n".join(problems)
        )


def validate_numeric_values(df, report):
    add_report(report, "NUMERIC VALIDATION")
    add_report(report, "------------------")

    numeric_columns = [
        "age",
        "time_from_treatment_start",
        *CELL_COUNT_COLUMNS,
    ]

    conversion_errors = []

    for column in numeric_columns:
        original = df[column]

        converted = pd.to_numeric(
            original,
            errors="coerce"
        )

        # Values that were non-null but couldn't be converted.
        invalid = original.notna() & converted.isna()

        if invalid.any():
            conversion_errors.append(
                f"{column}: {invalid.sum()} non-numeric values"
            )

        df[column] = converted

    if conversion_errors:
        raise ValueError(
            "Invalid numeric values found:\n"
            + "\n".join(conversion_errors)
        )

    # Age should be positive.
    invalid_age = df["age"].notna() & (df["age"] <= 0)

    if invalid_age.any():
        raise ValueError(
            f"Found {invalid_age.sum()} invalid age values."
        )

    # Time can be zero or positive.
    invalid_time = (
        df["time_from_treatment_start"].notna()
        & (df["time_from_treatment_start"] < 0)
    )

    if invalid_time.any():
        raise ValueError(
            f"Found {invalid_time.sum()} negative time values."
        )

    # Cell counts cannot be negative.
    for column in CELL_COUNT_COLUMNS:
        invalid_counts = (
            df[column].notna()
            & (df[column] < 0)
        )

        if invalid_counts.any():
            raise ValueError(
                f"Found {invalid_counts.sum()} negative "
                f"values in {column}."
            )

    add_report(report, "All numeric values passed validation.")
    add_report(report)


def validate_identifiers(df, report):
    add_report(report, "IDENTIFIER VALIDATION")
    add_report(report, "---------------------")

    identifier_columns = [
        "project",
        "subject",
        "sample",
    ]

    for column in identifier_columns:

        missing = df[column].isna().sum()

        if missing > 0:
            raise ValueError(
                f"{column} contains {missing} missing values."
            )

        empty = (
            df[column]
            .astype(str)
            .str.strip()
            .eq("")
            .sum()
        )

        if empty > 0:
            raise ValueError(
                f"{column} contains {empty} empty values."
            )

    # A biological sample should be uniquely identifiable.
    duplicate_samples = df["sample"].duplicated().sum()

    if duplicate_samples > 0:
        duplicate_ids = (
            df.loc[
                df["sample"].duplicated(keep=False),
                "sample"
            ]
            .unique()
            .tolist()
        )

        add_report(
            report,
            f"Duplicate sample IDs: {duplicate_samples}"
        )
        add_report(
            report,
            f"Duplicate IDs: {duplicate_ids[:20]}"
        )

        raise ValueError(
            f"Found {duplicate_samples} duplicate sample IDs."
        )

    add_report(report, "No duplicate sample IDs found.")
    add_report(report)


def validate_cell_counts(df, report):
    add_report(report, "CELL COUNT VALIDATION")
    add_report(report, "---------------------")

    for column in CELL_COUNT_COLUMNS:

        missing = df[column].isna().sum()

        if missing > 0:
            add_report(
                report,
                f"{column}: {missing} missing values"
            )

    # For cell-count analysis, missing cell counts are problematic
    # because total_count cannot be calculated reliably.
    missing_counts = df[CELL_COUNT_COLUMNS].isna().any(axis=1)

    if missing_counts.any():
        add_report(
            report,
            f"Rows with at least one missing cell count: "
            f"{missing_counts.sum()}"
        )
    else:
        add_report(
            report,
            "No missing cell counts."
        )

    add_report(report)


# -------------------------------------------------------------------
# Main cleaning pipeline
# -------------------------------------------------------------------

def main():

    report = []

    add_report(report, "CELL COUNT DATA CLEANING REPORT")
    add_report(report, "==============================")
    add_report(report)

    # ---------------------------------------------------------------
    # Load
    # ---------------------------------------------------------------

    if not INPUT_FILE.exists():
        raise FileNotFoundError(
            f"Could not find input file: {INPUT_FILE}"
        )

    df = pd.read_csv(INPUT_FILE, sep=None, engine="python")

    original_rows = len(df)

    add_report(
        report,
        f"Input file: {INPUT_FILE}"
    )
    add_report(
        report,
        f"Rows loaded: {original_rows}"
    )
    add_report(
        report,
        f"Columns loaded: {len(df.columns)}"
    )
    add_report(report)

    # ---------------------------------------------------------------
    # Column validation
    # ---------------------------------------------------------------

    validate_required_columns(df)

    # ---------------------------------------------------------------
    # Basic formatting
    # ---------------------------------------------------------------

    df = clean_string_columns(df)
    df = normalize_categorical_values(df)

    # ---------------------------------------------------------------
    # Validation
    # ---------------------------------------------------------------

    validate_missing_values(df, report)
    validate_categorical_values(df, report)
    validate_numeric_values(df, report)
    validate_identifiers(df, report)
    validate_cell_counts(df, report)

    # ---------------------------------------------------------------
    # Remove completely duplicated rows
    # ---------------------------------------------------------------

    duplicate_rows = df.duplicated().sum()

    if duplicate_rows > 0:
        add_report(
            report,
            f"Completely duplicated rows removed: {duplicate_rows}"
        )

        df = df.drop_duplicates()

    else:
        add_report(
            report,
            "No completely duplicated rows found."
        )

    add_report(report)

    # ---------------------------------------------------------------
    # Final summary
    # ---------------------------------------------------------------

    final_rows = len(df)

    add_report(report, "FINAL SUMMARY")
    add_report(report, "-------------")
    add_report(report, f"Rows before cleaning: {original_rows}")
    add_report(report, f"Rows after cleaning:  {final_rows}")
    add_report(
        report,
        f"Rows removed:         {original_rows - final_rows}"
    )

    missing_response = df["response"].isna().sum()

    add_report(
        report,
        f"Missing response values retained: {missing_response}"
    )

    add_report(report)
    add_report(
        report,
        "Cleaning completed successfully."
    )

    # ---------------------------------------------------------------
    # Save
    # ---------------------------------------------------------------

    df.to_csv(
        OUTPUT_FILE,
        index=False
    )

    REPORT_FILE.write_text(
        "\n".join(report),
        encoding="utf-8"
    )

    print("\n".join(report))

    print(
        f"\nCleaned dataset saved to: {OUTPUT_FILE}"
    )

    print(
        f"Cleaning report saved to: {REPORT_FILE}"
    )


if __name__ == "__main__":
    main()