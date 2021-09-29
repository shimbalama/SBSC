import pysam
from Bio import SeqIO
from collections import Counter, defaultdict
import scipy.stats as stats
import operator


def doit(tup):
    args, genomic_region = tup
    return genomic_region.run(args)


class GenomicRegion:

    def __init__(self, chrom, start, end, seq, tumour, normal):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.seq = seq
        self.tumour = tumour
        self.normal = normal
        self.homopolymer_positions = set([])
        self.var_d = {}

    def name(self):

        print(f'{str(self.chrom)}:{str(self.start)}:{str(self.end)}')

    def run(self, args):

        self.find_hompols()

        self.parse_tab(pysam.TabixFile(self.tumour), 't', args)

        self.parse_tab(pysam.TabixFile(self.normal), 'n', args)

        vars_caller = Vars()
        vars_caller.get_vars(self.var_d, self.homopolymer_positions)

        return vars_caller.variants  # annotated

    def annotate(self, d, args):

        tbx = pysam.TabixFile(args.gnomad)
        # UCSC
        repeats = self.parse_UCSC(args.low_complexity, 'genoName', 'genoStart', 'genoEnd')
        centromeres = self.parse_UCSC(args.meres, 'chrom', 'chromStart', 'chromEnd')
        segmental_dups = self.parse_UCSC(args.seg_dups, 'chrom', 'chromStart', 'chromEnd')
        print('ress', list(repeats)[:4])
        for pos, calls in d.items():
            chrom, position = pos.split(':')
            if pos in repeats:
                d[pos]['repeats'] = 'Y'
            else:
                d[pos]['repeats'] = 'N'
            if pos in centromeres:
                d[pos]['centromeres'] = 'Y'
            else:
                d[pos]['centromeres'] = 'N'
            if pos in segmental_dups:
                d[pos]['segmental_dups'] = 'Y'
            else:
                d[pos]['segmental_dups'] = 'N'
            d[pos]['gnomad'] = 'N'
            d[pos]['dbSNP'] = 'N'
            # this is to determine if region dark in short reads, not var specific.
            # very permissive. counts even if just indel overlap and +-5bp
            if args.gnomad:
                for i, row in enumerate(tbx.fetch(chrom, int(position) - 101, int(position) + 100)):
                    row = str(row).split()
                    # if row[1] == position:
                    d[pos]['gnomad'] = 'Y'
                    if row[1] == position:
                        if row[2].startswith('rs'):
                            d[pos]['dbSNP'] = 'Y'
        return d

    def parse_UCSC(self, ucsc, chrom, start, end):

        s = set([])
        with open(ucsc, 'r') as fin:
            header = fin.readline().strip().split()
            for line in fin:
                tmp_d = dict(zip(header, line.strip().split()))
                tmp_s = int(tmp_d.get(start))
                tmp_e = int(tmp_d.get(end))
                if tmp_d.get(chrom) == 'chr' + self.chrom and \
                        (self.start < tmp_s < self.end or
                         self.start < tmp_e < self.end):
                    for i in range(tmp_s, tmp_e):
                        s.add(tmp_d.get(chrom) + ':' + str(i))
        return s

    def extract_ins(self, res, val):

        # not currently used which sux as we lose seq info
        # problem is that len wobble means that the exact seq won't also be in
        # norm so germ gets through
        indels = []
        bits = res.strip().split(val)
        if len(bits) > 1:
            res = ''
            for bit in bits:
                if bit[0].isdigit():
                    indel_len = []
                    for char in bit:
                        if char.isdigit():
                            indel_len.append(char)
                        else:
                            break
                    indel_len = ''.join(indel_len)
                    indel = bit[len(indel_len):int(indel_len) + len(indel_len)]
                    indels.append(val + ':' + indel)
                    res += bit[len(indel_len) + int(indel_len):]
                else:
                    res += bit

        return indels, res

    def extract_indel(self, res, tmp_d, i, args):

        SVs = ['+-:break_point' for i in range(res.count('^') + res.count('$'))]
        # ends = ['+-:end' for i in range(res.count('$'))]
        if '^' in res:
            # ^ (caret) marks the start of a read segment and the ASCII
            # of the character following `^' minus 33 gives the mapping quality
            if res.startswith('^'):
                res = ''.join([bit[1:] for bit in res.split('^')])
            else:
                res = ''.join([bit[1:] if j != 0 else bit for j, bit in enumerate(res.split('^'))])
        ins, res = self.extract_ins(res, '+')
        dels, res = self.extract_ins(res, '-')

        ref_up = tmp_d.get('ref').upper()
        to_keep = ['.', ',', 'A', 'T', 'C', 'G', '*']  # keep * cuz have phred score...
        if i == 0:
            res = ''.join([nuc for nuc in res if nuc.upper() in to_keep])
            res = res.replace('.', ref_up).replace(',', tmp_d.get('ref').lower())
        else:
            res = ''.join([nuc for nuc in res if nuc.upper() in to_keep])
            res = res.replace('.', ref_up).replace(',', ref_up)
        if not len(res) == len(tmp_d.get('qual')):
            print('not_same_len', len(res), len(tmp_d.get('qual')),
                  res, tmp_d, i)
        res_qual = list(zip(res, tmp_d.get('qual')))
        # if (ord(quality) - 33) > args.base_qual]
        # this is realling killing depth - makes sense for SNVs but not for SVs/indels
        res = ''.join([base for base, quality in res_qual if (ord(quality) - 33) > args.base_qual])

        # for base QC reporting and len becomes the new read depth
        quals = [(read, quality) for read, quality in res_qual if (ord(quality) - 33) > args.base_qual]
        quals_unfiltered = [(read, quality) for read, quality in res_qual]  # for indel SV depth... bit messy tidy up logic one day
        vars_dict = Counter(SVs + ins + dels + list(res))  # only qual > args.base_qual are counted
        skip = ['>', '<', '^']  # '*'
        for i, var in enumerate(vars_dict.copy()):
            if var in skip:
                del vars_dict[var]
        return vars_dict, quals, quals_unfiltered

    def parse_tab(self, tabixfile, tn, args):

        # new pileup
        # chrm, pos, ref, depth, calls, phreds, position in read, read names, flags, map

        # keys = ['seq', 'pos', 'ref', 'reads', 'res', 'qual']#old
        keys = ['seq', 'pos', 'ref', 'reads', 'res', 'qual', 'pos_in_read', 'read_names', 'flags', 'map']
        # d={}
        for pos in tabixfile.fetch('chr' + str(self.chrom), self.start, self.end):
            tmp_d = dict(zip(keys, pos.split('\t')))
            if tn == 't':
                pos_d = defaultdict(lambda: defaultdict((dict)))#hack to deal with pickling
            chrom_pos = tmp_d.get('seq')+':'+tmp_d.get('pos')
            res = str(tmp_d.get('res'))
            res_ss = res.upper().replace(',', '.')
            for i, version in enumerate([res, res_ss]):  # i==0 is raw, 1 is all upper
                vars_dict, res_mod, res_mod_unfiltered = self.extract_indel(version, tmp_d, i, args)
                for base, count in vars_dict.items():
                    if tn == 't':
                        pos_d[i][base][tn] = count
                        pos_d[i][base]['n'] = 0
                    else:
                        self.var_d[chrom_pos][i][base][tn] = count
                tmp_d['res_mod'] = res_mod  # this should always be the upper version as its 2nd
                tmp_d['res_mod_unfiltered'] = res_mod_unfiltered
            if tn == 't':
                self.var_d[chrom_pos] = pos_d
            self.var_d[chrom_pos]['meta'][tn] = tmp_d

    def find_hompols(self):

        for i, nuc in enumerate(self.seq):
            if nuc.upper() != 'N':
                if len(self.seq[i:i+3]) == 3 and len(set(self.seq[i:i+3])) == 1:
                    # hom pols len 3 + get flagged
                    for j in range(i, i+5):  # in or adjacent to
                        self.homopolymer_positions.add(self.start+j)


class Vars:

    def __init__(self, variants=None):
        self.variants = variants
        if self.variants:
            self.del_positions = set([])
        else:
            self.variants = {}
            self.del_positions = None
        if self.del_positions:
            print(f'total len dels{len(self.del_positions)}')

    def get_proportion(self, percent_tumour, percent_normal):

        if all([percent_normal, percent_tumour]):
            return percent_normal/percent_tumour
        else:
            return 0.0

    def fishing(self,
                tumour_count,
                normal_count,
                non_base_tumour,
                non_base_normal,
                calls,
                read_count_tumour,
                read_count_normal,
                pos,
                base,
                info,
                homopolymer_positions):

        assert calls is not None
        if base[0] not in ['-', '+']:
            if not all(map(lambda x: x > -1, [
                    tumour_count,
                    normal_count,
                    non_base_tumour,
                    non_base_normal])):
                print('wat da!', tumour_count, normal_count, non_base_tumour,
                      non_base_normal, calls, read_count_tumour,
                      read_count_normal, pos, base, info)
        vals = [tumour_count, normal_count]
        for i, val in enumerate(vals):
            if val < 0:
                vals[i] = 0
        nons = [non_base_tumour, non_base_normal]
        for i, val in enumerate(nons):
            if val < 0:
                nons[i] = 0

        oddsratio, pvalue1 = stats.fisher_exact([vals, nons],
                                                alternative='greater')
        if pvalue1 < 0.1:
            # let everything through, then filter later
            meta = info.get('meta').get('t')
            calls['ref'] = meta.get('ref')
            calls['tumour_alts'] = [base]
            calls['normal_alts'] = self.get_max_val(info, pos, 'n')
            calls['tumour_P'] = [pvalue1]
            chrom, base_pos = pos.split(':')
            if int(base_pos) in homopolymer_positions:
                calls['hom_pol'] = 'Y'
            else:
                calls['hom_pol'] = 'N'
            calls = self.quals(calls, meta, 'tumour_qual', base, info)  # del info once testing done
            calls = self.quals(calls, info.get('meta').get('n'), 'normal_qual',
                               base, info)
            # check strand
            if info.get(1).get(base.upper()):
                base_upper = info.get(1).get(base.upper()).get('t', 0)
            else:
                base_upper = 0
            if info.get(0).get(base.lower()):
                base_lower = info.get(0).get(base.lower()).get('t', 0)
            else:
                base_lower = 0
            proportion = self.get_proportion(base_upper, base_lower)
            calls['strand'] = [proportion]
            calls['read_count_tumour'] = [[tumour_count, read_count_tumour]]
            calls['read_count_normal'] = [[normal_count, read_count_normal]]

        return calls, pvalue1

    def quals(self, calls, meta, key, base, info):

        quals = [(ord(quality) - 33) for read, quality in meta.get('res_mod')
                 if read == base]
        if len(quals) == 0:
            calls[key] = [0]
        else:
            calls[key] = [sum(quals)/len(quals)]
        quals = [(ord(quality) - 33) for read, quality in meta.get('res_mod')]
        if len(quals) == 0:
            calls[key+'_all'] = 0
        else:
            calls[key+'_all'] = sum(quals)/len(quals)
        return calls

    def test_base(self, info, base, calls, pos, homopolymer_positions):

        assert calls is not None
        tumour_count = int(info.get(1).get(base).get('t'))
        normal_count = int(info.get(1).get(base).get('n'))
        if base[0] not in ['+', '-']:  # arguments can be made either way but I think its better to have indels count for depth
            read_count_tumour = len(info.get('meta').get('t').get('res_mod'))
            read_count_normal = len(info.get('meta').get('n').get('res_mod'))
        else:
            read_count_tumour = len(info.get('meta').get('t').get('res_mod_unfiltered'))
            read_count_normal = len(info.get('meta').get('n').get('res_mod_unfiltered'))
        pvalue1 = 1
        if tumour_count > 3 and all([read_count_tumour, read_count_normal]):
            proportion = self.get_proportion(tumour_count/read_count_tumour,
                                             normal_count/read_count_normal)
            if proportion < 0.5:
                non_base_tumour = read_count_tumour - tumour_count
                non_base_normal = read_count_normal - normal_count
                calls, pvalue1 = self.fishing(tumour_count,
                                              normal_count,
                                              non_base_tumour,
                                              non_base_normal,
                                              calls,
                                              read_count_tumour,
                                              read_count_normal,
                                              pos,
                                              base,
                                              info,
                                              homopolymer_positions)
        return calls

    def get_vars(self, var_d, homopolymer_positions):

        for pos in list(var_d.keys()):
            chrom, position = pos.split(':')
            info = var_d.get(pos)
            if info.get(1):  # just fluf flike *^ etc
                for i, base in enumerate(info.get(1)):  # 1==uppers, 0==raw
                    calls = {}
                    if not base[0] in ['-', '+'] and base != '*':
                        if info.get(1).get(base).get('t'):
                            if info.get(1).get(base).get('t') > 3:
                                assert calls is not None
                                calls = self.test_base(info, base, calls, pos,
                                                       homopolymer_positions)
                    if calls:
                        if pos in self.variants:
                            self.add_more_calls(calls, pos)
                        else:
                            self.variants[pos] = calls
            del var_d[pos]  # free ram?

    def get_max_val(self, info_copy, pos, tn):

        d = {base: count.get(tn) if count.get(tn) else 0
             for base, count in info_copy.get(1).items() if base != '*'}

        return max(d.items(), key=operator.itemgetter(1))[0]

    def add_more_calls(self, calls, pos):

        cols = ['tumour_alts', 'tumour_P', 'tumour_qual', 'normal_qual',
                'strand', 'read_count_tumour', 'read_count_normal']
        for key, val in calls.items():
            if key in cols:
                try:
                    self.variants[pos][key] += val
                except Exception as e:
                    print('wtf', e, key, val, self.variants[pos])


def chunk_ref(args, chroms):
    '''
    Split ref into chunks for parallel processing
    '''
    chunks = []
    size = args.window_size
    total_len = 0
    for record in SeqIO.parse(args.ref, 'fasta'):
        if record.id in chroms:
            seqlen = len(record.seq)
            for start in range(0, seqlen, size):
                if seqlen > start + size:
                    end = start + size
                else:  # end of chrom
                    end = seqlen
                chunk = record.seq[start:end]
                chunks.append(GenomicRegion(
                    str(record.id).replace('chr', ''),
                    start,
                    end,
                    chunk,
                    args.cancer_pile,
                    args.normal_pile))
                total_len += len(chunk)

    assert total_len == sum([genomic_region.end - genomic_region.start for genomic_region in chunks])
    print('split', ' '.join(chroms), 'into', len(chunks), 'chunks')
    return chunks
