"""Class to store compiled information for a specific genomic position."""
from __future__ import annotations

from dataclasses import dataclass, field

from reditools.constants import bases as base_order
from reditools.constants import (
    comp_map,
    empty_field,
    forward_strand_number,
    forward_strand_symbol,
    reverse_strand_number,
    reverse_strand_symbol,
    undetermined_strand_number,
    undetermined_strand_symbol,
)


def _calc_edit_ratio(counter: dict[str, int], ref: str) -> float:
    """
    Caclulate the edit ratio.

    Parameters
    ----------
    counter : dict[str, int]
        Base counts
    ref : str
        Base from genomic reference

    Returns
    -------
    float
        Edit ratio, or zero if no edits detected.
    """
    max_edits = 0
    for base, count in counter.items():
        if base != ref and count > max_edits:
            max_edits = count
    try:
        return max_edits  / (counter[ref] + max_edits)
    except ZeroDivisionError:
        return 0

def _strand_from_string(strand_value: str) -> str:
    """
    Get the strand symbol from string input.

    Depending on the parameters used from the analyze tool, the strand may
    be represented as a number ("0", "1", or "2") or as a string
    ("+", "-", "*"). This function converts the former to the latter, or just
    reports the number if that was used instead.

    Parameters
    ----------
    strand_value : str
        Value from the Strand column in analyze output

    Returns
    -------
    str
        Strand number (as a string)
    """
    if strand_value == forward_strand_number:
        return forward_strand_symbol
    if strand_value == reverse_strand_number:
        return reverse_strand_symbol
    if strand_value == undetermined_strand_number:
        return undetermined_strand_symbol
    return strand_value

@dataclass
class CompiledPosition:
    """Class to store compiled information for a specific genomic position.

    Attributes
    ----------
    ref : str
        The reference nucleotide at this position.
    position : int
        The genomic position (0-based).
    contig : str
        The name of the contig/chromosome.
    qualities : list[int]
        List of Phred quality scores for each base at this position.
    strands : list[str]
        List of strand orientations ('+', '-', or '*') for each base.
    bases : list[str]
        List of nucleotides observed at this position.
    """

    ref: str
    position: int
    contig: str
    qualities: list[int] = field(default_factory=list)
    strands: list[str] = field(default_factory=list)
    bases: list[str] = field(default_factory=list)

    def __len__(self) -> int:
        """Return the number of bases at this position.

        Returns
        -------
        int
            Number of bases.
        """
        return len(self.bases)

    def add_base(self, quality: int, strand: str, base: str) -> None:
        """Add a base observation to this position.

        Parameters
        ----------
        quality : int
            Phred quality score of the base.
        strand : str
            Strand orientation of the read containing the base
            ('+', '-', or '*').
        base : str
            The observed nucleotide.
        """
        self.bases.append(base)
        self.strands.append(strand)
        self.qualities.append(quality)

    def calculate_strand(self, threshold: float = 0) -> str:
        """Calculate the consensus strand based on observations.

        Parameters
        ----------
        threshold : float, default 0
            The fraction of observations required to assign a strand.

        Returns
        -------
        str
            The calculated strand ('+', '-', or '*').
        """
        pos_count = 0
        neg_count = 0
        for strand in self.strands:
            if strand == forward_strand_symbol:
                pos_count += 1
            elif strand == reverse_strand_symbol:
                neg_count += 1
        if pos_count == neg_count:
            return undetermined_strand_symbol
        if pos_count / (pos_count + neg_count) >= threshold:
            return forward_strand_symbol
        if neg_count / (pos_count + neg_count) >= threshold:
            return reverse_strand_symbol
        return undetermined_strand_symbol

    def filter_by_strand(self, strand: str) -> None:
        """Filter observations to keep only those from a specific strand.

        Parameters
        ----------
        strand : str
            The strand to keep ('+', '-', or '*'). If '*', no filtering is done.
        """
        if strand == undetermined_strand_symbol:
            return
        keep = [
            idx for idx in range(len(self.bases))
            if self.strands[idx] == strand
        ]
        self.qualities = [self.qualities[_] for _ in keep]
        self.strands = [self.strands[_] for _ in keep]
        self.bases = [self.bases[_] for _ in keep]

    def split_by_strand(self) -> tuple[CompiledPosition, ...]:
        """Split the CompiledPosition into new objects for each strand.

        Returns
        -------
        tuple[CompiledPosition]
            CompiledPositions for each strand.

        """
        idx_by_strand: dict[str, list[int]] = {
            "+": [],
            "-": [],
            "*": [],
        }
        for idx in range(len(self.bases)):  # noqa: WPS518
            idx_by_strand[self.strands[idx]].append(idx)

        return tuple(
            CompiledPosition(
                self.ref,
                self.position,
                self.contig,
                [self.qualities[_] for _ in strand_sublist],
                [self.strands[_] for _ in strand_sublist],
                [self.bases[_] for _ in strand_sublist],
            )
            for strand_sublist in idx_by_strand.values()
            if len(strand_sublist) > 0
        )

    def complement(self) -> None:
        """Replace all bases and the reference with their complements."""
        self.bases = [comp_map[base] for base in self.bases]
        self.ref = comp_map[self.ref]

@dataclass(slots=True)
class RTResult:
    """Class to represent the results of REDItools analysis at a position.

    Attributes
    ----------
    cp : CompiledPosition
        The compiled position data.
    strand : str
        The strand assigned to this position.
    reference : str
        The reference nucleotide.
    position : int
        The genomic position.
    contig : str
        The contig name.
    counter : dict
        A dictionary counting occurrences of each nucleotide.
    variants : list
        A list of observed alternate bases ('A', 'C', 'G', and/or 'T').
    """

    strand: str
    reference: str
    position: int
    contig: str
    base_counts: dict[str, int]
    variants: list[str]
    depth: int
    edit_ratio: float
    mean_quality: float
    dna: RTResult | None = None

    def complement(self) -> None:
        """Compute base complement.

        Updates the reference, variants, and base_counts properties.
        """
        self.reference = comp_map[self.reference]
        self.base_counts = {
            base: self.base_counts[comp_map[base]]
            for base in self.base_counts
        }
        self.variants = [comp_map[_] for _ in self.variants]

    def output_list(self, strand_as_numbers: bool = False) -> list[str]:
        """Generate list for Analayze tool output.

        Parameters
        ----------
        strand_as_numbers : bool
            If True, uses numbers for the Strand field. Otherwise uses symbols.

        Returns
        -------
        list[str]
            Output fields for the Analyze tool.
        """
        if strand_as_numbers:
            if self.strand == reverse_strand_symbol:
                strand = reverse_strand_number
            elif self.strand == forward_strand_symbol:
                strand = forward_strand_number
            else:
                strand = undetermined_strand_number
        else:
            strand = self.strand
        if self.variants:
            variants = " ".join([
                f"{self.reference}{_}" for _ in self.variants
            ])
        else:
            variants = empty_field
        output_list = [
            self.contig,
            str(self.position + 1),
            self.reference,
            strand,
            str(self.depth),
            f"{self.mean_quality:.2f}",
            str(list(self.base_counts.values())),
            variants,
            f"{self.edit_ratio:.2f}",
        ]
        if self.dna is None:
            output_list.extend([empty_field for _ in range(5)])
        else:
            output_list.extend(self.dna.output_list()[4:9])
        return output_list

def rtresult_from_dict(src_dict: dict[str, str]) -> RTResult:
    """Create an RTResult object from a dictionary.

    This will NOT fill in annotation column values.

    Parameters
    ----------
    src_dict : dict[str, str]
        Dictionary data from a REDItools analyze file.

    Returns
    -------
    RTResult
    """
    strand = _strand_from_string(src_dict["Strand"])
    if "Coverage" in src_dict:
        depth = src_dict["Coverage"]  # noqa: WPS529
    else:
        depth = src_dict["Coverage-q30"]
    counter = {
        base: int(count)
        for base, count in zip(
            base_order,
            src_dict["BaseCount[A,C,G,T]"][1:-1].split(", "),
        )
    }
    if src_dict["AllSubs"] == "-":
        variants = []
    else:
        variants=[_[1] for _ in src_dict["AllSubs"].split(" ")]
    return RTResult(
        strand=strand,
        reference=src_dict["Reference"],
        position=int(src_dict["Position"]) - 1,
        contig=src_dict["Region"],
        base_counts=counter,
        variants=variants,
        depth=int(depth),
        edit_ratio=float(src_dict["Frequency"]),
        mean_quality=float(src_dict["MeanQ"]),
    )


def rtresult_from_compiled_position(
    cp: CompiledPosition,
    strand: str,
    noise_threshold: float = 0,
) -> RTResult:
    """Initialize RTResult.

    Parameters
    ----------
    compiled_position : CompiledPosition
        The compiled position data.
    strand : str
        The strand assigned to this position.
    noise_threshold : float
        Noise threshold for populating variants list.
    """
    counter = dict.fromkeys(base_order, 0)
    for base in cp.bases:
        counter[base] += 1

    edit_ratio = _calc_edit_ratio(counter, cp.ref)

    mean_quality = (
        sum(cp.qualities) / len(cp)
        if len(cp) else 0
    )

    ref_reads = counter[cp.ref]
    variants = [
        _ for _ in sorted(
            base_order,
            key=counter.__getitem__,
            reverse=True,
        )
        if counter[_] and \
            counter[_] / (counter[_] + ref_reads) > noise_threshold and \
            _ != cp.ref
    ]

    return RTResult(
        strand=strand,
        reference=cp.ref,
        position=cp.position,
        contig=cp.contig,
        base_counts=counter,
        variants=variants,
        depth=len(cp),
        edit_ratio=edit_ratio,
        mean_quality=mean_quality,
    )

