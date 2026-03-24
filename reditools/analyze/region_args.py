from reditools import utils
from reditools.region import Region


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

    args = []
    for contig, size in utils.get_contigs(bam_fname):
        region = Region(contig=contig, start=1, stop=size+1)
        if window:
            args.extend(region.split(window))
        else:
            args.append(region)
    return args
