class Aligner:
    def __init__(self, match=1, mismatch=1, gap=1):
        self.match = 1
        self.mismatch = 1
        self.gap = 1

    def align(self, ref_seq, query_seq):
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

    def trace_matrix(self, ref_seq, qry_seq, trace_mat):
        row_idx = len(ref_seq)
        col_idx = len(qry_seq)

        ref_align = []
        qry_align = []
        while row_idx > 0 or col_idx > 0:
            trace_val = trace_mat[row_idx][col_idx]

            if trace_val in [2, 5, 6, 9]:
                row_idx -= 1
                col_idx -= 1
                ref_align.insert(0, ref_seq[row_idx])
                qry_align.insert(0, qry_seq[col_idx])
            elif trace_val in [3, 7]:
                row_idx -= 1
                ref_align.insert(0, ref_seq[row_idx])
                qry_align.insert(0, '-')
            elif trace_val == 4:
                col_idx -= 1
                ref_align.insert(0, '-')
                qry_align.insert(0, qry_seq[col_idx])

        return (
            ''.join(ref_align),
            ''.join(qry_align),
        )

class NWMatrix:
    def __init__(self, ref_seq, query_seq, match, mismatch, gap):
        self.ref_seq = ref_seq
        self.qry_seq = query_seq
        self.gap = gap
        self.match = match
        self.mismatch = mismatch

        self.nw_matrix = self.init_nw_matrix(ref_seq, query_seq)
        self.trace_matrix = self.init_trace_matrix(ref_seq, query_seq)

    def assess_cell(self, col_idx, ref_base, row_idx, query_base):
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

    def run_dp(self):
        for col_idx, ref_base in enumerate(self.qry_seq):
            for row_idx, query_base in enumerate(self.ref_seq):
                self.assess_cell(col_idx, ref_base, row_idx, query_base)

    def init_nw_matrix(self, ref_seq, query_seq):
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

    def init_trace_matrix(self, ref_seq, query_seq):
        trace_mat = [
            [3] + [0 for _ in range(len(query_seq))]
            for _ in range(len(ref_seq))
        ]
        trace_mat.insert(
            0,
            [4 for _ in range(len(query_seq) + 1)],
        )
        return trace_mat
