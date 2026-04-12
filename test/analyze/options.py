from dataclasses import InitVar, dataclass


@dataclass
class Options:
    bamfile: InitVar[str]
    region: str | None = None
    window: int = 0

    def __post_init__(self, bamfile):
        self.file = [bamfile]  # noqa: WPS110
