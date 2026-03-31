#!/usr/bin/env python3

import re
from pathlib import Path
import pandas as pd
from pdfminer.high_level import extract_text
from pdfminer.layout import LAParams

# -----------------------------
# Config
# -----------------------------
PDF_FILE = "amjuel.pdf"
OUTPUT_DIR = "amjuel_data"

SECTION_MATRIX = {
    "H.0": (1, 9),
    "H.1": (1, 9),
    "H.2": (1, 9),
    "H.3": (9, 9),
    "H.4": (9, 9),
    "H.5": (1, 9),
    "H.6": (9, 9),
    "H.7": (9, 9),
    "H.8": (1, 9),
    "H.9": (9, 9),
    "H.10": (9, 9),
    "H.11": (1, 9),
    "H.12": (9, 9),
}

# -----------------------------
# Regex patterns
# -----------------------------
section_pattern = re.compile(r"H\.\d+")
reaction_pattern = re.compile(r"Reaction\s+(\S+)")
coef_pattern = re.compile(r"[-+]?\d+\.\d+[eEdD][-+]\d+")

# -----------------------------
# Helpers
# -----------------------------
def convert_fortran(n):
    """Convert Fortran D exponent to float."""
    return float(n.replace("D", "E"))

def build_matrix(values, rows, cols):
    """Convert flat list into matrix."""
    return [[v for c in range(3) for v in values[c*27+i*3:c*27+i*3+3]] for i in range(rows)]

def build_vector(values, rows):
    """Convert flat list into matrix."""
    return [values[:rows]]

def safe_filename(name):
    """Make filename filesystem-safe."""
    return re.sub(r"[^\w.-]", "_", name)

def save_matrix(section, reaction, matrix, outdir):
    """Save matrix to CSV according to AMJUEL rules."""
    filename = outdir / f"{section}_{safe_filename(reaction)}.csv"
    df = pd.DataFrame(matrix)

    rows, cols = SECTION_MATRIX[section]

    if rows == 1:
        df = df.T
        df.to_csv(filename, header=False, index=False, float_format="%.12e")

    else:
        # transpose + no headers
        df.to_csv(filename, header=False, index=False, float_format="%.12e")

# -----------------------------
# Parse PDF
# -----------------------------
print("Reading AMJUEL PDF...")
laparams = LAParams(
    line_margin=0.2,
    char_margin=3.0,
    word_margin=0.1,
    boxes_flow=None  # controls column detection
)
text = extract_text(PDF_FILE, laparams= laparams)
lines = text.splitlines()

current_section = ""
current_reaction = ""
coeffs = []

database = {}

def store_reaction():
    """Store reaction if enough coefficients collected."""

    if current_section not in SECTION_MATRIX:
        return

    rows, cols = SECTION_MATRIX[current_section]
    needed = rows * cols

    # allow extra coefficients for 1-row sections
    if rows == 1 and len(coeffs) >= cols:
    # if rows == 1:
        matrix = build_vector(coeffs, len(coeffs))
        database.setdefault(current_section, {}).setdefault(current_reaction, []).append(matrix)
    # elif len(coeffs) >= needed:
    else:
        matrix = build_matrix(coeffs[:needed], rows, cols)
        database.setdefault(current_section, {}).setdefault(current_reaction, []).append(matrix)

for idx, line in enumerate(lines):
    line = line.strip()
    if not line:
        continue


    # detect section
    sec = section_pattern.search(line)
    if sec:
        if current_reaction!="":
            store_reaction()
        current_reaction=""
        coeffs = []
        current_section = sec.group()
        continue

    # detect reaction
    reac = reaction_pattern.search(line)
    if reac:
        if current_reaction!="":
            store_reaction()
        current_reaction = reac.group(1)
        coeffs = []
        continue

    # extract coefficients
    nums = coef_pattern.findall(line)
    if nums:
        coeffs.extend(convert_fortran(n) for n in nums)

# store last reaction
store_reaction()

# -----------------------------
# Write output
# -----------------------------
output = Path(OUTPUT_DIR)
output.mkdir(exist_ok=True)

print("Saving CSV files...")

for section, reactions in database.items():
    section_dir = output / section
    section_dir.mkdir(exist_ok=True)

    for reaction, matrices in reactions.items():
        for matrix in matrices:
            save_matrix(
                section,
                reaction,
                matrix,
                section_dir,
            )

print("\nFinished.")
print("Output directory:", output)
print("Sections parsed:", len(database))