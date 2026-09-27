PRAGMA foreign_keys = ON;

DROP TABLE IF EXISTS cell_counts;
DROP TABLE IF EXISTS samples;
DROP TABLE IF EXISTS subjects;

CREATE TABLE subjects (
    subject TEXT PRIMARY KEY,
    project TEXT NOT NULL,
    condition TEXT NOT NULL,
    age INTEGER NOT NULL,
    sex TEXT NOT NULL
);

CREATE TABLE samples (
    sample TEXT PRIMARY KEY,
    subject TEXT NOT NULL,
    treatment TEXT NOT NULL,
    response TEXT,
    sample_type TEXT NOT NULL,
    time_from_treatment_start INTEGER NOT NULL,

    FOREIGN KEY (subject)
        REFERENCES subjects(subject)
);

CREATE TABLE cell_counts (
    sample TEXT NOT NULL,
    population TEXT NOT NULL,
    count REAL NOT NULL,

    PRIMARY KEY (sample, population),

    FOREIGN KEY (sample)
        REFERENCES samples(sample)
);