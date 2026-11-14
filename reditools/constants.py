"""Constants."""

"""Complementary base dictionary."""
comp_map = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
    "N": "N",
    "-": "-",
}

"""Nucleotides in output order."""
bases = ("A", "C", "G", "T")

"""Output strand symbols."""
forward_strand_symbol = "+"
reverse_strand_symbol = "-"
undetermined_strand_symbol = "*"

"""Output strand numbers."""
forward_strand_number = "1"
reverse_strand_number = "0"
undetermined_strand_number = "2"

"""Empty output field string."""
empty_field = "-"

"""Analyze output fieldnames."""
analyze_output_fieldnames = [
    "Region",
    "Position",
    "Reference",
    "Strand",
    "Coverage",
    "MeanQ",
    "BaseCount[A,C,G,T]",
    "AllSubs",
    "Frequency",
    "gCoverage",
    "gMeanQ",
    "gBaseCount[A,C,G,T]",
    "gAllSubs",
    "gFrequency",
]

