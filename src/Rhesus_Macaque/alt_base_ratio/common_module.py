from common_parameter import *


def get_df(file_path, mode_ , v_type):
    fpd = pd.read_csv(file_path, index_col=0).astype(str)
    fpd_snp = fpd[fpd['TYPE'].str.contains(v_type + '.snp')].reset_index(drop=True)
    filter_refcall(fpd_snp)
    # print(fpd)

    if mode_ == 'filter_sampled':
        fpd_target = fpd_snp[(fpd_snp['MATCH'] == '0') & (fpd_snp['result'] == '1')].reset_index(drop=True)
    else:
        fpd_target = fpd_snp

    # print(mode_, fpd_target)

    return fpd_target


def check_alt(base, ref, alt):
    alt_list = alt.split(',')
    alt_list_ = list()
    for r in alt_list:
        alt_list_.append(r.strip())
    if '*' in alt_list_:
        if base != ref:
            return True
        else:
            return False
    else:
        if base in alt_list_:
            return True
        else:
            return False


# get alt base ratio of rhesus macaque individual by open .bam file
def get_dis(fpd_target, file_path_dedup, save_path, sampling_n):
    index_list = list()
    dis_list = []
    bamfile = pysam.AlignmentFile(file_path_dedup, "rb")
    reffile = pysam.FastaFile(ref_file_path)

    while len(index_list) < sampling_n:
        sample_index = np.random.choice(fpd_target.index, 1)[0]
        if sample_index in index_list:
            continue
        sample_info = {'CHR': fpd_target.loc[sample_index, 'CHROM'], 'POS': int(fpd_target.loc[sample_index, 'POS']),
                       'REF': fpd_target.loc[sample_index, 'REF'], 'ALT': fpd_target.loc[sample_index, 'ALT']}

        ref_base = reffile.fetch(rhe_mac_chr_dict_10_reverse[sample_info['CHR']],
                                 sample_info['POS'] - 1,sample_info['POS']).upper()
        if ref_base == sample_info['REF']:
            pass
        else:
            print('error in', sample_info)
            continue

        base_counter = Counter()
        total_read = 0

        for pileupcolumn in bamfile.pileup(rhe_mac_chr_dict_10_reverse[sample_info['CHR']], sample_info['POS'] - 1,
                                           sample_info['POS'], min_base_quality=0, min_mapping_quality=0):
            if pileupcolumn.pos == sample_info['POS'] - 1:  # 0-based indexing
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        read_base = pileupread.alignment.query_sequence[pileupread.query_position].upper()
                        base_counter[read_base] += 1
                        total_read += 1

        var_count = 0
        for base, count in base_counter.items():
            if check_alt(base, ref_base, sample_info['ALT']):
                var_count += count

        if var_count == 0 or var_count > total_read:
            continue
        else:
            dis_list.append(var_count/total_read * 100)
            index_list.append(sample_index)

    sample_index_ = np.sort(index_list)

    fpd_sampled = fpd_target.loc[sample_index_]
    fpd_sampled = fpd_sampled.reset_index(drop=True)
    fpd_sampled['distribution'] = dis_list

    fpd_sampled.to_csv(save_path)

    print('distribution')
    print('mean :', np.mean(dis_list))
    print('var :', np.var(dis_list))
    print(len(dis_list) )

    bamfile.close()
    reffile.close()

    return dis_list


# sampling the variants and get the alt base ratio from all individuals
def sample_fpd(fpd, save_path, sampling_n):
    sample_index = np.random.choice(fpd.index, size=sampling_n, replace=False)
    sample_index_ = np.sort(sample_index)

    fpd_sampled = fpd.loc[sample_index_]
    fpd_sampled = fpd_sampled.reset_index(drop=True)

    fpd_sampled.to_csv(save_path)

    return fpd_sampled['distribution'].to_list()


# plot KDE distribution of alt base ratio
def plot_dis(dis_list, plot_path, title_):
    plt.figure(figsize=(10, 6))

    color_l = ['red', 'blue']

    for i, dis_dict in enumerate(dis_list):
        dis = dis_dict['distribution']
        dis_label = dis_dict['label']
        sns.kdeplot(dis, color=color_l[i], fill=True, label=dis_label)

    plt.xticks(np.arange(0, 101, 10))
    plt.xlim(0, 100)

    plt.title(title_, fontsize=18)
    plt.xlabel('ABR (%)', fontsize=16)
    plt.ylabel('ABR Density', fontsize=16)
    plt.legend(loc='upper left', fontsize=14)
    plt.savefig(plot_path, dpi=300)


# Single sample t-test
def statistic_normal_check(dis, alpha):
    t_statistic, p_value = stats.ttest_1samp(dis, alpha)

    if p_value < 0.05:
        return t_statistic, p_value, 'False'
    else:
        return t_statistic, p_value, 'True'


def statistic_dis_compare(dis1, dis2, is_t_test=False, is_wt_test=False):
    t1_stat = None
    p1_value = None
    t_result = None
    t2_stat = None
    p2_value = None
    wt_result = None

    if is_t_test:
        # Independent sample t-test
        t1_stat, p1_value = stats.ttest_ind(dis1, dis2)

        if p1_value > 0.05:
            t_result = 'similar'
        else:
            t_result = 'not similar'
    if is_wt_test:
        # Welch's t-test
        t2_stat, p2_value = stats.ttest_ind(dis1, dis2, equal_var=False)

        if p2_value > 0.05:
            wt_result = 'similar'
        else:
            wt_result = 'not similar'

    return t1_stat, p1_value, t_result, t2_stat, p2_value, wt_result


def filter_refcall(df):
    return df[(df['FILTER'] != 'RefCall') & (df['TYPE'].str.contains('.snp'))].reset_index(drop=True)
