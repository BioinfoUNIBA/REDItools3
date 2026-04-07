import csv
from tempfile import NamedTemporaryFile


def write_results(rtools, sam_manager, file_name, region, output_format,
                  temp_dir):
    """
    Write the results from a REDItools analysis to a temporary file.

    Parameters:
        rtools (REDItools): REDItools instance
        sam_manager (AlignmentManager): Source of reads
        file_name (string): Input file name for analysis
        region: Region to analyze
        output_format (dict): keyword arguments for csv.writer constructor.

    Returns:
        string: Name of the temporary file.
    """
    with NamedTemporaryFile(mode='w', delete=False, dir=temp_dir) as stream:
        writer = csv.writer(stream, **output_format)
        for rt_result in rtools.analyze(sam_manager, region):
            variants = rt_result.variants
            writer.writerow([
                rt_result.contig,
                rt_result.position + 1,
                rt_result.reference,
                rt_result.strand,
                rt_result.depth,
                f'{rt_result.mean_quality:.2f}',
                rt_result.per_base_depth,
                ' '.join(sorted(variants)) if variants else '-',
                f'{rt_result.edit_ratio:.2f}',
                '-', '-', '-', '-', '-',
            ])
        return stream.name
