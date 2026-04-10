from reditools.logger import Logger


def check_min_read_depth(namespace, bases):
    if len(bases) < namespace.min_read_depth:
        return (
            'DISCARDING COLUMN {} [MIN_READ_DEPTH={}]',
            len(bases),
            namespace.min_read_depth,
        )
