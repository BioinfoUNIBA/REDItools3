from reditools.region import Region
from pysam import AlignmentFile


def region_args(bam_fname, region, window):
    """
    Split a region into segments for paralllel processing.

    Parameters:
        bam_fname (str): BAM file to collect contig info from
        region (Region): Genomic region to split
        window (int): How large the sub regions should be.

    Returns:
        (list): Sub regions
    """
    if region is not None:
        if window:
            return region.split(window)
        return [region]

    sub_regions = []
    with AlignmentFile(bam_fname, ignore_truncation=True) as bam:
        for contig, size in zip(bam.references, bam.lengths):
            region = Region(contig=contig, start=0, stop=size)
            if window:
                sub_regions.extend(region.split(window))
            else:
                sub_regions.append(region)
    return sub_regions
