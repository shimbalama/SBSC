from dataclasses import dataclass, field

import pandas as pd
from Bio import SeqIO


@dataclass(frozen=True, slots=True)
class GenomicRegion:
    """Data pertaining to a chromosomal region"""

    chrom: str
    start: int
    end: int
    seq: str = field(repr=False)
    seq_padding_start: str = field(repr=False)
    seq_padding_end: str = field(repr=False)
    homopolymer_positions: set[int] = field(init=False, repr=False)

    def __post_init__(self):
        """set attribute homopolymer_positions"""

        object.__setattr__(self, "homopolymer_positions", self.find_hom_pol_positions())

    def __str__(self) -> str:
        return f"{self.chrom}:{str(self.start)}-{str(self.end)}"

    def find_hom_pol_positions(self, length: int = 3) -> set[int]:
        """Returns all zero based genomic positions that are in or
        adjacent to homoploymers of given length.

        Example:
        self.seq_padding_start = 'GTA'
        self.seq_padding_end = 'GTA'
        self.seq = 'AATATGGGTGAT'
        self.start = 10

        Returns:
        {9, 10, 11, 12, 13, 15, 16, 17, 18, 19}
        """
        positions = set([])
        seq = self.seq_padding_start + self.seq + self.seq_padding_end
        start = self.start - len(self.seq_padding_start)
        for pos, nuc in enumerate(seq):
            if nuc.upper() != "N":
                if (
                    len(seq[pos : pos + length]) == length
                    and len(set(seq[pos : pos + length])) == 1
                ):
                    for hom_pol_pos in range(
                        pos - 1, pos + length + 1
                    ):
                        positions.add(start + hom_pol_pos + 1)
        return positions

    def get_hom_pol_lengths(self) -> dict[int, int]:
        """
        Adds the length of the associated homolymer to to each zero based position
        in homopolymer_positions.

        Example:
        self.homopolymer_positions = {9, 10, 11, 12, 13, 15, 16, 17, 18, 19}

        Returns:
        {9: 3, 10: 3, 11: 3, 12: 3, 13: 3, 15: 3, 16: 3, 17: 3, 18: 3, 19: 3}"""

        df = pd.DataFrame({"pos": sorted(self.homopolymer_positions)})
        df["rank"] = df["pos"].rank(method="first")
        df["key"] = df["pos"] - df["rank"]
        result = df.groupby("key")["pos"].agg(["min", "max"])
        result["length"] = (result.iloc[:, 1] - result.iloc[:, 0]) - 1
        d = {}
        df = df.reset_index()  # make sure indexes pair with number of rows
        for _, row in result.iterrows():
            start: int = row["min"]
            stop: int = row["max"] + 1
            for pos in range(start, stop):
                d[pos] = row["length"]
        return d

    @property
    def homopolymer_lengths(self) -> dict[int, int]:
        return self.get_hom_pol_lengths()


def chunk_ref(size: int, reference: str, chroms: list[str]):
    """
    Split ref into chunks for parallel processing
    """
    for record in SeqIO.parse(reference, "fasta"):
        if record.id in chroms:
            end_of_chrom = False
            seqlen = len(record.seq)
            seq = str(record.seq)
            for start in range(0, seqlen, size):
                assert not end_of_chrom
                if seqlen > start + size:
                    end = start + size
                else:
                    end_of_chrom = True
                    end = seqlen
                seq_chunk = seq[start:end]
                chars = set(seq_chunk)
                if len(chars) == 1 and chars.pop().upper() == "N":
                    continue
                else:
                    pad_start = seq[start - 1000 : start]
                    pad_end = seq[end : end + 1000] if not end_of_chrom else ""
                    yield GenomicRegion(
                        str(record.id).replace("chr", ""),
                        start,
                        end - 1,
                        seq_chunk,
                        pad_start,
                        pad_end,
                    )
