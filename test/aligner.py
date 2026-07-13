"""Needleman-Wunsch sequence aligner."""
from __future__ import annotations


class Aligner:
    """Needleman-Wunsch sequence aligner."""

    _trace_match = frozenset((2, 5, 6, 9))
    _trace_del = frozenset((3, 7))
    _trace_ins = 4

    def __init__(self, match: int=1, mismatch: int=1, gap: int=1) -> None:
        """Initialize self.

        Parameters
        ----------
        match : int
            Match score.
        mismatch : int
            Mismatch penalty.
        gap : int
            Gap penalty.
        """
        self.match = match
        self.mismatch = mismatch
        self.gap = gap

    def align(self, ref_seq: str, query_seq: str) -> tuple[str, str]:
        """Perform alignment.

        Parameters
        ----------
        ref_seq : str
            Reference sequence.
        query_seq : str
            Query sequence.

        Returns
        -------
        tuple[str, str]
            Reference and query alignment strings, respectively.
        """
        matrix = NWMatrix(
            ref_seq,
            query_seq,
            self.match,
            self.mismatch,
            self.gap,
        )
        matrix.run_dp()
        return self.trace_matrix(
            ref_seq,
            query_seq,
            matrix.trace_matrix,
        )

    def trace_matrix(
        self,
        ref_seq: str,
        qry_seq: str,
        trace_mat: list[list[int]],
    ) -> tuple[str, str]:
        """Trace path from matrix to create alignment strings.

        Parameters
        ----------
        ref_seq : str
            Reference sequence.
        qry_seq : str
            Query sequence.
        trace_mat : list[list[int]]
            Trace matrix.

        Returns
        -------
        tuple[str, str]
            Reference and query alignment strings, respectively.
        """
        row_idx = len(ref_seq)
        col_idx = len(qry_seq)

        ref_align: list[str] = []
        qry_align: list[str] = []
        while row_idx > 0 or col_idx > 0:
            trace_val = trace_mat[row_idx][col_idx]

            if trace_val in self._trace_match:
                row_idx -= 1
                col_idx -= 1
                ref_align.insert(0, ref_seq[row_idx])
                qry_align.insert(0, qry_seq[col_idx])
            elif trace_val in self._trace_del:
                row_idx -= 1
                ref_align.insert(0, ref_seq[row_idx])
                qry_align.insert(0, "-")
            elif trace_val == self._trace_ins:
                col_idx -= 1
                ref_align.insert(0, "-")
                qry_align.insert(0, qry_seq[col_idx])

        return (
            "".join(ref_align),
            "".join(qry_align),
        )

class NWMatrix:
    """Needleman-Wunsch matrix and trace matrix object."""

    def __init__(
        self,
        ref_seq: str,
        query_seq: str,
        match: int,
        mismatch: int,
        gap: int,
    ) -> None:
        """Initialize self.

        Parameters
        ----------
        ref_seq : str
            Reference sequence.
        query_seq : str
            Query sequence.
        match : int
            Match score.
        mismatch : int
            Mismatch penalty.
        gap : int
            Gap penalty.
        """
        self.ref_seq = ref_seq
        self.qry_seq = query_seq
        self.gap = gap
        self.match = match
        self.mismatch = mismatch

        self.nw_matrix = self.init_nw_matrix(ref_seq, query_seq)
        self.trace_matrix = self.init_trace_matrix(ref_seq, query_seq)

    def assess_cell(
        self,
        col_idx: int,
        ref_base: str,
        row_idx: int,
        query_base: str,
    ) -> None:
        """Fill in a cell of the matrices.

        Parameters
        ----------
        col_idx : int
            Column index - 1.
        ref_base : str
            Reference nucleotide.
        row_idx : int
            Row index - 1.
        query_base : str
            Query nucleotide.
        """
        align_val = self.match if ref_base == query_base else -self.mismatch
        t_list = [
            self.nw_matrix[row_idx][col_idx] + align_val,
            self.nw_matrix[row_idx][col_idx + 1] - self.gap,
            self.nw_matrix[row_idx + 1][col_idx] - self.gap,
        ]
        t_max = max(t_list)
        self.nw_matrix[row_idx + 1][col_idx + 1] = t_max
        self.trace_matrix[row_idx + 1][col_idx + 1] += sum((
            idx + 2 for idx, tv in enumerate(t_list) if tv == t_max
        ))

    def run_dp(self) -> None:
        """Run dynamic programming."""
        for col_idx, ref_base in enumerate(self.qry_seq):
            for row_idx, query_base in enumerate(self.ref_seq):
                self.assess_cell(col_idx, ref_base, row_idx, query_base)

    def init_nw_matrix(self, ref_seq: str, query_seq: str) -> list[list[int]]:
        """Create a starting Needleman-Wunsch matrix.

        Parameters
        ----------
        ref_seq : str
            Reference sequence.
        query_seq : str
            Query sequence.

        Returns
        -------
        list[list[int]]
            Initialized matrix.
        """
        nw_mat = [
            [-self.gap * (row_idx + 1)] + \
                [0 for _ in range(len(query_seq))]
            for row_idx in range(len(ref_seq))
        ]
        nw_mat.insert(
            0, [
                -self.gap * _
                for _ in range(len(query_seq) + 1)
            ],
        )
        return nw_mat

    def init_trace_matrix(
        self,
        ref_seq: str,
        query_seq: str,
    ) -> list[list[int]]:
        """Create an initialized trace matrix.

        Parameters
        ----------
        ref_seq : str
            Reference sequence.
        query_seq : str
            Query sequence.

        Returns
        -------
        list[list[int]]
            Initialized matrix.
        """
        trace_mat = [
            [3] + [0 for _ in range(len(query_seq))]
            for _ in range(len(ref_seq))
        ]
        trace_mat.insert(
            0,
            [4 for _ in range(len(query_seq) + 1)],
        )
        return trace_mat
