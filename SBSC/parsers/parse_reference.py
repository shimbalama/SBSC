from dataclasses import dataclass, field
from typing import Set
from Bio import SeqIO
import pandas as pd

@dataclass(frozen=True, slots=True)
class GenomicRegion:
    chrom: str
    start: int
    end: int
    seq: str = field(repr=False)
    seq_padding_start: str = field(repr=False)
    seq_padding_end: str = field(repr=False)
    homopolymer_positions: Set = field(init=False, repr=False)

    def __post_init__(self):
        '''set attribute 'homopolymer_positions'''

        object.__setattr__(self, 'homopolymer_positions', self.find_hom_pol_positions())

    def __str__(self) -> str:
        return f'{self.chrom}:{str(self.start)}-{str(self.end)}'

    def find_hom_pol_positions(self, length=3):
     
        s = set([])
        seq = self.seq_padding_start + self.seq + self.seq_padding_end
        start = self.start - len(self.seq_padding_start)
        for pos, nuc in enumerate(seq):
            if nuc.upper() != 'N':
                if len(seq[pos:pos+length]) == length and len(set(seq[pos:pos+length])) == 1:
                    # hom pols len 3 + get flagged
                    for hom_pol_pos in range(pos - 1, pos + length + 1):  # in or adjacent to
                        s.add(start + hom_pol_pos + 1)
        return s
        
    def get_hom_pol_lengths(self):

        df = pd.DataFrame({'pos': sorted(self.homopolymer_positions)})
        df['rank'] = df['pos'].rank(method='first')
        df['key'] = df['pos'] - df['rank']
        result = (df.groupby('key')['pos'].agg(['min', 'max']))
        result['length'] = (result.iloc[:,1] - result.iloc[:,0]) -1
        d= {}
        df = df.reset_index()  # make sure indexes pair with number of rows
        for _, row in result.iterrows():
            start = row['min']
            stop = row['max'] + 1
            for pos in range(start, stop):
                d[pos] = row['length']
        return d

    @property
    def homopolymer_lengths(self) -> list[str]:
        return self.get_hom_pol_lengths()

    
def chunk_ref(args, chroms):
    '''
    Split ref into chunks for parallel processing
    '''
    size = args.window_size #- dont go lower than 100k, or super slow
    for record in SeqIO.parse(args.ref, 'fasta'):
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
                if len(chars) == 1 and chars.pop().upper() == 'N':
                    continue
                else:
                    pad_start = seq[start - 1000: start]
                    pad_end = seq[end: end + 1000] if not end_of_chrom else ''
                    yield GenomicRegion(
                        str(record.id).replace('chr', ''),#TODO - why rm chr?
                        start,
                        end -1,
                        seq_chunk,
                        pad_start,
                        pad_end
                    )
