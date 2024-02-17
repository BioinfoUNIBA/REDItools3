'''Representation of genomic regions'''
import re

class Region:
    '''Representation of genomic regions'''
    @staticmethod
    def _to_int(text):
        return int(text.replace(",", ""))
    def __init__(self, text):
        match = re.match(
            r'(?P<contig>[^:]+)(:(?P<start>\d+(,\d{3})*)-(?P<stop>\d+(,\d{3})*))?',
            text)
        if match is None:
            raise ValueError(f"Cannot parse region: {text}")

        self.contig = match.group("contig")
        self.start = match.group("start")
        self.stop = match.group("stop")

        if self.start is not None:
            try:
                self.start = Region._to_int(self.start)
                self.stop = Region._to_int(self.stop)
            except Exception as exc:
                raise ValueError(f"Cannot parse region: {text}") from exc

    def args(self):
        '''Reformats the Region into a tuple.

        Returns:
            tuple (contig, start, stop)
        '''
        return (self.contig, self.start, self.stop)
    def kwargs(self):
        '''Reformates the Region into a dictionary.

        Returns:
            dict {"contig", "start", "stop"}
        '''
        return {"contig": self.contig, "start": self.start, "stop": self.stop}
