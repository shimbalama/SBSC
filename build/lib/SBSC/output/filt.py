from collections import defaultdict
from itertools import product
import pandas as pd
import numpy as np


def filt(d, args):

    d2 = defaultdict(lambda: defaultdict(list))
    for pos, calls in d.items():
        for k, v in calls.items():
            if type(v) != list:
                d2[pos][k] = v
            else:
                for i, option in enumerate(v):
                    p = args.P_value
                    alt = calls.get('tumour_alts')[i]
                    percent_norm = calls.get('read_count_normal')[i][0] /\
                        calls.get('read_count_normal')[i][1]
                    if percent_norm > 0.07:
                        continue
                    if alt[0] in ['-', '+']:
                        if percent_norm > 0.00:
                            continue
                    if alt == calls.get('ref'):  # any LOH that slip through that arent interesting
                        continue
                    if calls.get('read_count_tumour')[i][1] > 120 and \
                            calls.get('read_count_normal')[i][1] > 60:
                        p *= 0.001  # if in CNV requires lower P
                    if calls.get('read_count_tumour')[i][1] > 480 and \
                            calls.get('read_count_normal')[i][1] > 240:
                        p *= 0.00001  # if in CNV requires lower P
                    if 0.07 > calls.get('strand')[i] or \
                            calls.get('strand')[i] > 0.93:
                        continue
                    if float(calls.get('tumour_P')[i]) > p:
                        continue
                    if calls.get('read_count_tumour')[i][1] < args.min_depth_tumour or \
                        calls.get('read_count_normal')[i][1] < args.min_depth_normal:
                        break  # HOM dels that actually line up properly should have read depth zero
                    if all(map(lambda x: x > args.min_mean_base_quality, [
                            float(calls.get('tumour_qual')[i]),
                            calls.get('tumour_qual_all'),
                            calls.get('normal_qual_all')])):
                        d2[pos][k].append(option)
        if not d2[pos].get('tumour_alts'):
            del d2[pos]
    df = pd.DataFrame.from_dict(d2, orient='index')

    if not len(df):
        return

    # make pseudo-VCF
    # TODO make a proper VCF
    df['tmp_index_col'] = df.index
    df[['chrom', 'pos']] = df.tmp_index_col.str.split(':', expand=True)
    df.drop('tmp_index_col', axis=1, inplace=True)
    cols = df.columns.tolist()
    cols = cols[-2:] + cols[:-2]
    df = df[cols]

    if args.MNVs:
        df = MNV(df)

    # multiple alts as separate records
    df = df.explode(['tumour_alts',
                     'tumour_P', 'tumour_qual', 'normal_qual', 'strand',
                     'read_count_tumour', 'read_count_normal'])

    fout_name = '_'.join([
        args.output,
        str(args.P_value).split('.')[1],
        str(args.min_mean_base_quality),
        str(args.min_depth_tumour),
        str(args.min_depth_normal)])

    df.to_csv(fout_name + '.tsv', index=False, sep='\t')


def distance(r1, r2):
    """
    Return r2['pos'] - r1['pos'] if r1['chrom'] == r2['chrom']
    (corollary, return 0 iff r1 and r2 have same position)
    Return np.inf if r1['chrom'] != r2['chrom']
    """
    if r1.get('chrom') is None or r1.get('chrom') != r2.get('chrom'):
        return np.inf
    return r2['pos'] - r1['pos']


def mergepos(altalleles):
    """
    Given a list-of-lists representation of (possibly poly-allelic)
    contiguous SNVs, return a list-of-sequences representation of the MNV
    alleles at the first position

    Args:
        altalleles: list of lists, the inner lists representing allele(s)
            at first and subsequent contiguous positions of a variant.
            For example:
            [['A']]                         # mono-allelic SNV
            [['A', 'C']]                    # poly-allelic SNV
            [['A'], ['C']]                  # mono-allelic MNV
            [['A', 'C'], ['A'], ['C', 'G']] # poly-allelic MNV
    Returns:
        list of MNV alleles at the first position
        >>> mergepos([['A']])
        ['A']
        >>> mergepos([['A']])
        ['A', 'C']
        >>> mergepos([['A'], ['C']])
        ['AC']
        >>> # note: the result is Cartesian product
        >>> mergepos([['A', 'C'], ['A'], ['C', 'G']])
        ['AAC', 'AAG', 'CAC', 'CAG']
    """
    if not altalleles:
        return []
    return [''.join(bases) for bases in product(*altalleles)]


def mergesnvs(snvs):
    """
    Takes a list of contiguous SNV records with a list-of-alleles
    representation in tumour_alts and returns MNV starting at the first
    position.

    Returns:
        Empty list if snvs is empty, otherwise a single-element list
        containing the MNV record

    Raises:
        ValueError if SNV positions are not contiguous
    """
    if not snvs:
        return []

    posx = [snv['pos'] for snv in snvs]
    posx0 = [pos - posx[0] for pos in posx]
    if not posx0 == list(range(len(snvs))):
        raise ValueError(f'{posx} are not contiguous')

    tumour_alts = mergepos([snv['tumour_alts'] for snv in snvs])
    n_alleles = len(tumour_alts)
    # simplification: in combining contiguous SNV calls into MNVs it's not
    # immediately obvious what to do with the different values for tumour_P,
    # tumour_qual, normal_qual, strand, read_count_tumour, read_count_normal.
    # Currently we take the simplest possible approach by repeating the first
    # value from the first position only. Other simple approaches include
    # taking average values or worst values across positions and alleles.
    first = snvs[0]
    return [{
        **snvs[0],
        'tumour_alts':       tumour_alts,
        'ref':               ''.join(snv['ref'] for snv in snvs),
        'normal_alts':       ''.join(snv['normal_alts'] for snv in snvs),
        'tumour_P':          [first['tumour_P'][0]] * n_alleles,
        'tumour_qual':       [first['tumour_qual'][0]] * n_alleles,
        'normal_qual':       [first['normal_qual'][0]] * n_alleles,
        'strand':            [first['strand'][0]] * n_alleles,
        'read_count_tumour': [first['read_count_tumour'][0]] * n_alleles,
        'read_count_normal': [first['read_count_normal'][0]] * n_alleles
    }]


def MNV(df):
    """
    Convert any contiguous SNVs to MNVs.
    """
    df.pos = df.pos.astype(int)
    df.sort_values(by=['chrom', 'pos'])
    records = df.to_dict('records')

    merged = []
    prev = {}
    contiguous = []
    for record in records:
        d = distance(prev, record)
        if d == 1:
            contiguous.append(record)
        else:
            assert d > 1, 'records must be sorted by (chrom, pos) ascending'
            merged.extend(mergesnvs(contiguous))
            contiguous = [record]
        prev = record
    merged.extend(mergesnvs(contiguous))

    return df.from_records(
        merged, index=[f'{x["chrom"]}:{x["pos"]}' for x in merged])
