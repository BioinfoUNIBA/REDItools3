from pysam import AlignmentFile

from reditools.region import Region


def region_args(options):
    """
    Split a region into segments for paralllel processing.

    Parameters:
        options (namespace): analyze tool options

    Returns:
        list: region windows
    """
    if options.region is not None:
        region = Region.from_string(options.region, options.file[0])
        if options.window:
            return region.split(options.window)
        return [region]

    sub_regions = []
    with AlignmentFile(options.file[0], ignore_truncation=True) as bam:
        for contig, size in zip(bam.references, bam.lengths):
            region = Region(contig=contig, start=0, stop=size)
            if options.window:
                sub_regions.extend(region.split(options.window))
            else:
                sub_regions.append(region)
    return sub_regions
