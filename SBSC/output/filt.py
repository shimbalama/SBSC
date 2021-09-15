from collections import defaultdict
import pandas as pd


def filt(d, args, chroms):

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

    if len(df):
        if args.MNVs:
            df = MNVs(df)
        fout_name = '_'.join([
            args.output,
            str(args.P_value).split('.')[1],
            str(args.min_mean_base_quality),
            str(args.min_depth_tumour),
            str(args.min_depth_normal)])
        # get rid of lsits
        df = df.explode('tumour_alts')
        df = df.explode('tumour_P')
        df = df.explode('tumour_qual')
        df = df.explode('normal_qual')
        df = df.explode('strand')
        df = df.explode('read_count_tumour')
        df = df.explode('read_count_normal')

        df.to_csv(fout_name + '.tsv', sep='\t')


def MNVs(df):

    df['index_col'] = df.index
    df[['chrom', 'pos']] = df.index_col.str.split(':', expand=True)
    df.pos = df.pos.astype(int)
    already_updated = set([])
    for chrom, df_tmp in df.groupby('chrom'):
        df_tmp = df_tmp.sort_values(by=['pos'])
        for var_len in [5, 4, 3, 2, 1]:
            df_tmp['dif'] = df_tmp.pos.diff(var_len)
            hits = df_tmp[df_tmp['dif'] == var_len]
            if len(hits):  # not empty
                for hit in hits.pos:
                    rows = df_tmp[(df_tmp.pos <= hit) & (df_tmp.pos >= (hit - var_len))]
                    update = dict(rows.iloc[0])
                    if len(update.get('tumour_alts')) == 1 and len(update.get('tumour_alts')[0]) == 1:
                        for i in range(var_len):
                            i += 1
                            update_tmp = dict(rows.iloc[i])
                            if len(update.get('tumour_alts')) == 1 and len(update.get('tumour_alts')[0]) == 1:
                                key = update_tmp.get('chrom') + ':'+str(update_tmp.get('pos'))
                                if key not in already_updated:
                                    update['tumour_alts'][0] += update_tmp.get('tumour_alts')[0]
                                    update['tumour_P'] += update_tmp.get('tumour_P')
                                    update['tumour_qual'] += update_tmp.get('tumour_qual')
                                    update['normal_qual'] += update_tmp.get('normal_qual')
                                    update['strand'] += update_tmp.get('strand')
                                    update['read_count_tumour'] += update_tmp.get('read_count_tumour')
                                    update['read_count_normal'] += update_tmp.get('read_count_normal')
                                    df = df.drop(index=(key))
                                    already_updated.add(key)
                            else:
                                break
    return df
