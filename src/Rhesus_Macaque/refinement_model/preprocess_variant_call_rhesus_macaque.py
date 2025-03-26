import pandas as pd

from common_parameter import *


def get_data_path(origin, ind):
    if origin == 'Indian':
        gt_path = path_preprocess_data + 'vcf/indian_origin/indian_origin_ground_truth.vcf'
        dv_path = path_preprocess_data + 'vcf/indian_origin/' + str(ind) + '_output.vcf'
        gt_csv_path = path_preprocess_data + 'csv/indian_origin/gt_labeled_' + str(ind) + '.csv'
        dv_csv_path = path_preprocess_data + 'csv/indian_origin/dv_labeled_' + str(ind) + '.csv'
        match_csv_path = path_preprocess_data + 'csv/indian_origin/matched_csv_' + str(ind) + '.csv'
        learning_csv_path = path_preprocess_data + 'csv/indian_origin/learning_csv_' + str(ind) + '.csv'
        result_csv_path = path_result_data + 'csv/indian_origin/result_csv_' + str(ind) + '.csv'
        confusion_csv_path = path_result_data + 'csv/indian_origin/confusion_' + str(ind) + '.csv'
        score_csv_path = path_result_data + 'csv/indian_origin/score_' + str(ind) + '.csv'
        alignment_path = path_preprocess_data + 'sequence/indian_origin/' + str(ind) + '_dedup.bam'
        alignment_csv_path = path_preprocess_data + 'csv/indian_origin/aligned_csv_' + str(ind) + '.csv'
    else:
        gt_path = path_preprocess_data + 'vcf/chinese_origin/china_origin_ground_truth.vcf'
        dv_path = path_preprocess_data + 'vcf/chinese_origin/' + str(ind) + '_output.vcf'
        gt_csv_path = path_preprocess_data + 'csv/chinese_origin/gt_labeled_' + str(ind) + '.csv'
        dv_csv_path = path_preprocess_data + 'csv/chinese_origin/dv_labeled_' + str(ind) + '.csv'
        match_csv_path = path_preprocess_data + 'csv/chinese_origin/matched_csv_' + str(ind) + '.csv'
        learning_csv_path = path_preprocess_data + 'csv/chinese_origin/learning_csv_' + str(ind) + '.csv'
        result_csv_path = path_result_data + 'csv/chinese_origin/result_csv_' + str(ind) + '.csv'
        confusion_csv_path = path_result_data + 'csv/chinese_origin/confusion_' + str(ind) + '.csv'
        score_csv_path = path_result_data + 'csv/chinese_origin/score_' + str(ind) + '.csv'
        alignment_path = path_preprocess_data + 'sequence/chinese_origin/' + str(ind) + '_dedup.bam'
        alignment_csv_path = path_preprocess_data + 'csv/chinese_origin/aligned_csv_' + str(ind) + '.csv'
    return gt_path, dv_path, gt_csv_path, dv_csv_path, match_csv_path, \
        learning_csv_path, result_csv_path, confusion_csv_path,\
        score_csv_path, alignment_path, alignment_csv_path

# vcf to csv for rhesus macaque individual
def read_vcf(path):
    with open(path, 'r') as f:
        lines_info = [l for l in tqdm(f) if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines_info)),
        dtype={'#CHROM': str, 'POS': str, 'ID': str, 'REF': str, 'ALT': str, 'QUAL': float, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})


# vcf to csv for rhesus macaque label
def read_vcf_(path, ind_list, id=None):
    if id is not None:
        idx = ind_list.index(id) + 8
        with open(path, 'r') as f:
            lines_info = [l.split(' ')[0:8] + [l.split(' ')[idx]] for l in tqdm(f) if not l.startswith('#')]
        column_ = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "FORMAT", "default"]
        df_ = pd.DataFrame(lines_info, columns=column_)
        return df_
    else:
        with open(path, 'r') as f:
            lines_info = [l.split(' ') for l in tqdm(f) if not l.startswith('#')]
        column_ = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "FORMAT", "default"]
        df_ = pd.DataFrame(lines_info, columns=column_)
        return df_


def check_snp(ref, alt_set):
    if len(ref) > 1:
        return False

    for l_ in alt_set:
        if l_ != '*' and len(ref) != len(l_):
            return False
    else:
        return True


def return_type_(info, ref, alt):
    # print(info)
    if pd.isna(info):
        return 'WRONG'
    if len(info.split(':')[0].split('/')) == 1 or \
            (info.split(':')[0].split('/')[0] == '.' and info.split(':')[0].split('/')[1] == '.'):
        o1 = 'DROP'
        o2 = ''
    else:
        if info.split(':')[0].split('/')[0] == info.split(':')[0].split('/')[1]:
            o1 = 'hm'
        else:
            o1 = "ht"

        if len(ref) == 0:
            o2 = '.wrong'
        else:
            alt_set = set()
            for l_ in alt.split(','):
                alt_set.add(l_)

            if check_snp(ref, alt_set):
                o2 = '.snp'
            else:
                o2 = '.indel'
    if 'DROP' in o1 + o2 or 'wrong' in o1 + o2:
        return 'WRONG'
    else:
        return o1 + o2


# convert vcf to csv and labeling the variant types
def labeling_M(label_data, idx_number, gt_path, dv_path, gt_csv_path, dv_csv_path, ind_list):
    if label_data == 'gt':
        label_path = gt_path
        csv_path = gt_csv_path

        vcf_df = read_vcf_(label_path, ind_list, idx_number)
    else:
        label_path = dv_path
        csv_path = dv_csv_path

        vcf_df = read_vcf(label_path)
        vcf_df = vcf_df[vcf_df['CHROM'].isin(rhe_mac_chr_dict_10.keys())]
        vcf_df['CHROM'] = vcf_df['CHROM'].replace(rhe_mac_chr_dict_10)
        vcf_df = vcf_df.reset_index(drop=True)

    vcf_df['QUAL'] = vcf_df['QUAL'].astype(float)
    vcf_df['CHROM'] = vcf_df['CHROM'].astype(str)

    vcf_df = vcf_df[vcf_df['CHROM'].isin(chrom_values)]
    vcf_df = vcf_df.reset_index(drop=True)

    print(vcf_df)
    vcf_df['TYPE'] = ''
    labels = dict()
    for i in tqdm(range(len(vcf_df.index))):
        labels[i] = return_type_(vcf_df.loc[i, 'default'], vcf_df.loc[i, 'REF'], vcf_df.loc[i, 'ALT'])

    label_df = pd.DataFrame.from_dict(labels, orient='index', columns=['label'])
    vcf_df['TYPE'] = label_df['label']
    if 'INFO' in vcf_df.columns:
        vcf_df = vcf_df.drop('INFO', axis=1)
    vcf_df = vcf_df.loc[vcf_df['TYPE'] != 'WRONG']
    vcf_df = vcf_df.reset_index(drop=True)
    vcf_df.to_csv(csv_path)


# comparing and labeling the label between ground truth and rhesus macaque individuals
def match_L(gt_csv_path, dv_csv_path, match_csv_path):
    gt_df = pd.read_csv(gt_csv_path, index_col=0)
    print(gt_df)

    gt_dict = dict()
    for i in tqdm(range(len(gt_df.index))):
        gt_dict[str(gt_df.loc[i, 'CHROM']) + '_' + str(gt_df.loc[i, 'POS'])] = {'REF': gt_df.loc[i, 'REF'],
                                                                                'ALT': gt_df.loc[i, 'ALT'],
                                                                                'TYPE': gt_df.loc[i, 'TYPE'],
                                                                                'default': gt_df.loc[i, 'default']}

    count_dict = {'NOT_MATCH': 0, 'MATCH': 0, 'MATCH_BUT_WRONG': 0, 'MATCH_RIGHT_BUT_DIFF': 0, 'MATCH_ALL_RIGHT': 0,
                  'DIFFERENT_ALLELE': 0}
    labels = dict()
    dv_df = pd.read_csv(dv_csv_path, index_col=0)
    print(dv_df)
    for i in tqdm(range(len(dv_df.index))):
        key_ = str(dv_df.loc[i, 'CHROM']) + '_' + str(dv_df.loc[i, 'POS'])
        ref = dv_df.loc[i, 'REF']
        alt = dv_df.loc[i, 'ALT']
        _type = dv_df.loc[i, 'TYPE']
        fflag = None
        if key_ in gt_dict:
            ori_ref = gt_dict[key_]['REF']
            ori_ALT = gt_dict[key_]['ALT']
            ori_TYPE = gt_dict[key_]['TYPE']
            ori_INFO = gt_dict[key_]['default']
            count_dict['MATCH'] += 1
            if _type == ori_TYPE:
                if ref == ori_ref and (alt == ori_ALT or alt in ori_ALT):
                    count_dict['MATCH_ALL_RIGHT'] += 1
                    fflag = 'MATCH_ALL_RIGHT'
                else:
                    count_dict['MATCH_RIGHT_BUT_DIFF'] += 1
                    fflag = 'MATCH_RIGHT_BUT_DIFF'
                    if ref != ori_ref:
                        count_dict['DIFFERENT_ALLELE'] += 1
                        fflag += '_diff_allele'
            else:
                count_dict['MATCH_BUT_WRONG'] += 1
                fflag = 'MATCH_BUT_WRONG'
        else:
            count_dict['NOT_MATCH'] += 1
            ori_ref = None
            ori_ALT = None
            ori_TYPE = None
            ori_INFO = None
            fflag = 'NOT_MATCH'
        labels[i] = [ori_ref, ori_ALT, ori_TYPE, ori_INFO, fflag]

    print(count_dict)
    label_df = pd.DataFrame.from_dict(labels, orient='index',
                                      columns=['ORI_REF', 'ORI_ALT', 'ORI_TYPE', 'ORI_INFO', 'MATCH'])
    dv_df = dv_df.join(label_df[['ORI_REF', 'ORI_ALT', 'ORI_TYPE', 'ORI_INFO', 'MATCH']])
    dv_df.to_csv(match_csv_path)


def extract_alignment_info(bam, chrom, pos):
    cigar_features = []
    mapq_values = []
    chr_ = rhe_mac_chr_dict_10_reverse[str(chrom)]
    strand_counts = {"forward": 0, "reverse": 0}  # Strand count 초기화

    for read in bam.fetch(chr_, int(pos) - 1, int(pos)):  # 0-based 좌표
        # Strand 정보 계산
        if read.flag & 16:  # Reverse strand
            strand_counts["reverse"] += 1
        else:  # Forward strand
            strand_counts["forward"] += 1

        # CIGAR 분석
        cigar_string = read.cigarstring
        if cigar_string is None:
            continue

        # Count CIGAR operations
        operations = re.findall(r'(\d+)([MIDNSHP=X])', cigar_string)
        feature = {"M": 0, "S": 0, "I": 0, "D": 0}
        for length, op in operations:
            if op in feature:
                feature[op] += int(length)

        feature["Total_length"] = sum(feature.values())  # 리드 총 길이
        feature["M_ratio"] = feature["M"] / feature["Total_length"] if feature["Total_length"] > 0 else 0
        feature["S_ratio"] = feature["S"] / feature["Total_length"] if feature["Total_length"] > 0 else 0
        feature["CIGAR"] = cigar_string

        # MAPQ
        mapq_values.append(read.mapping_quality)

        # Merge features
        cigar_features.append(feature)

    # 평균 CIGAR feature 계산
    avg_cigar = {
        "M_ratio": np.mean([feat["M_ratio"] for feat in cigar_features]) if cigar_features else None,
        "S_ratio": np.mean([feat["S_ratio"] for feat in cigar_features]) if cigar_features else None,
        "Total_reads": len(cigar_features),
    }

    # MAPQ 통계 계산
    mapq_stats = {
        "mean_mapq": np.mean(mapq_values) if mapq_values else None,
        "low_mapq_ratio": sum(1 for mq in mapq_values if mq < 20) / len(mapq_values) if mapq_values else None,
    }

    # Strand 비율 계산
    total_strands = strand_counts["forward"] + strand_counts["reverse"]
    strand_ratios = {
        "forward_strand_ratio": strand_counts["forward"] / total_strands if total_strands > 0 else 0,
        "reverse_strand_ratio": strand_counts["reverse"] / total_strands if total_strands > 0 else 0,
    }

    return {**avg_cigar, **mapq_stats, **strand_ratios}


def process_variant(args):
    """병렬 처리를 위한 함수."""
    bam_path, chrom, pos = args
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        align_info = extract_alignment_info(bam, chrom, pos)

        return {"CHROM": chrom, "POS": pos, **align_info}


def add_align(match_csv_path, alignment_path, alignment_csv_path, num_processes=32):
    match_df = pd.read_csv(match_csv_path, index_col=0)
    match_df = match_df[match_df['TYPE'].str.contains('.snp')].reset_index(drop=True)

    # 병렬 처리를 위한 인자 준비 (벡터화 사용)
    args_list = list(zip([alignment_path] * len(match_df), match_df['CHROM'], match_df['POS']))

    # 병렬 처리 실행
    with Pool(num_processes) as pool:
        align_info_results = list(tqdm(pool.imap(process_variant, args_list), total=len(args_list)))

    # 병렬 처리 결과를 DataFrame으로 변환
    align_info_df = pd.DataFrame(align_info_results)

    # match_df와 병합
    align_df = match_df.merge(align_info_df, on=["CHROM", "POS"], how="left")
    align_df.fillna({
        "M_ratio": 0, "S_ratio": 0, "Total_reads": 0,
        "mean_mapq": 0, "low_mapq_ratio": 0,
        "forward_strand_ratio": 0, "reverse_strand_ratio": 0
    }, inplace=True)
    align_df.to_csv(alignment_csv_path)
    return align_df


# extract the features
def make_learning_data(data_csv, column_name, learning_csv_path):
    data_csv['QUAL'] = data_csv['QUAL'].astype(float)
    data_csv['M_ratio'] = data_csv['M_ratio'].astype(float)
    data_csv['S_ratio'] = data_csv['S_ratio'].astype(float)
    data_csv['Total_reads'] = data_csv['Total_reads'].astype(float)
    data_csv['mean_mapq'] = data_csv['mean_mapq'].astype(float)
    data_csv['low_mapq_ratio'] = data_csv['low_mapq_ratio'].astype(float)
    data_csv['forward_strand_ratio'] = data_csv['forward_strand_ratio'].astype(float)
    data_csv['reverse_strand_ratio'] = data_csv['reverse_strand_ratio'].astype(float)
    data_csv['FORMAT'] = data_csv['FORMAT'].astype(str)
    data_csv[column_name] = data_csv[column_name].astype(str)

    print("Convert to learning data...")

    # parsing features
    format_split = data_csv['FORMAT'].str.split(':').apply(lambda x: x[1:])
    column_split = data_csv[column_name].str.split(':').apply(lambda x: x[1:])

    print('calculate meand diff')

    # calculate mean and diff for each feature
    mean_values = column_split.apply(lambda x: [np.mean(eval(i)) for i in x])
    diff_values = column_split.apply(lambda x: [np.max(eval(i)) - np.min(eval(i)) for i in x])

    print('make new column')

    # 새로운 컬럼 생성
    feature_columns = pd.DataFrame({
        f"{c}_mean": mean_values.str[idx] for idx, c in enumerate(format_split.iloc[0])
    }).join(
        pd.DataFrame({
            f"{c}_diff": diff_values.str[idx] for idx, c in enumerate(format_split.iloc[0])
        })
    )

    # 최종 데이터 생성
    feature_columns['QUAL'] = data_csv['QUAL']
    feature_columns['M_ratio'] = data_csv['M_ratio']
    feature_columns['S_ratio'] = data_csv['S_ratio']
    feature_columns['Total_reads'] = data_csv['Total_reads']
    feature_columns['mean_mapq'] = data_csv['mean_mapq']
    feature_columns['low_mapq_ratio'] = data_csv['low_mapq_ratio']
    feature_columns['forward_strand_ratio'] = data_csv['forward_strand_ratio']
    feature_columns['reverse_strand_ratio'] = data_csv['reverse_strand_ratio']
    feature_columns['MATCH'] = data_csv['MATCH']

    feature_columns.to_csv(learning_csv_path, index=False)
    return feature_columns


# apply the GVRP refinement model
def apply_model(input_csv, data_csv, model_dir):
    print('\nmodel inference step')
    feature_order = ['AD_diff', 'DP_mean', 'GQ_mean', 'PL_diff', 'VAF_mean', 'QUAL', 'M_ratio', 'S_ratio',
                     'Total_reads', 'mean_mapq', 'low_mapq_ratio', 'forward_strand_ratio']  # 학습할 때 사용했던 순서
    m = 'mix_all'
    tag = '_all_info'

    data_type_ = input_csv[feature_order]
    print(len(data_type_))
    print(data_type_)

    model = pickle.load(open(model_dir + m + '/' + tag + '_LGBM_model.pkl', 'rb'))

    pred_list = model.predict(data_type_)
    pred_p_list = model.predict_proba(data_type_)[:, 1]
    print(Counter(pred_list))

    data_csv['result'] = pred_list
    data_csv['result_p'] = pred_p_list


# apply the GVRP refinement model
def rhe_L(match_csv_path, learning_csv_path, result_csv_path, model_dir):
    match_csv = pd.read_csv(match_csv_path, index_col=0)
    match_csv['MATCH'] = match_csv['MATCH'].apply(lambda x: 1 if x == 'MATCH_ALL_RIGHT' else 0)
    match_csv = filter_refcall(match_csv)
    print(match_csv['MATCH'].value_counts())
    print(match_csv)

    input_csv = make_learning_data(match_csv, 'default', learning_csv_path)
    # input_csv = pd.read_csv(learning_csv_path)

    apply_model(input_csv, match_csv, model_dir)
    print(match_csv)
    print(match_csv['result'].value_counts())
    match_csv.to_csv(result_csv_path)


# calculate the confusion matrix
def return_metric(df_t):
    tp = len(df_t[df_t['confusion_M'] == 'TP'].reset_index(drop=True).index)
    tn = len(df_t[df_t['confusion_M'] == 'TN'].reset_index(drop=True).index)
    fp = len(df_t[df_t['confusion_M'] == 'FP'].reset_index(drop=True).index)
    fn = len(df_t[df_t['confusion_M'] == 'FN'].reset_index(drop=True).index)
    gt_1 = len(df_t[df_t['MATCH'] == 1].reset_index(drop=True).index)
    gt_0 = len(df_t[df_t['MATCH'] == 0].reset_index(drop=True).index)
    y_ = df_t['MATCH']
    predict_proba = df_t['result_p']
    if tp == 0 and fp == 0:
        precision = None
    else:
        precision = round(tp / (tp + fp), 3)
    if tp == 0 and fn == 0:
        recall = None
    else:
        recall = round(tp / (tp + fn), 3)
    if gt_1 == 0 and gt_0 == 0:
        miscalling_dv = None
    else:
        miscalling_dv = round(gt_0 / (gt_1 + gt_0), 3)
    if tp == 0 and fp == 0:
        miscalling_gvrp = None
    else:
        miscalling_gvrp = round(fp / (tp + fp), 3)
    if tp == 0 and tn == 0 and fp == 0 and fn == 0:
        accuracy = None
    else:
        accuracy = round((tp + tn) / (tp + tn + fp + fn), 3)
    if precision == 0 or precision is None or recall == 0 or recall is None:
        f1_score = None
    else:
        f1_score = round(2 / (1 / precision + 1 / recall), 3)

    fpr, tpr, threshold_ = roc_curve(y_, predict_proba)
    auc_ = auc(fpr, tpr)

    return tp, tn, fp, fn, precision, recall, f1_score, miscalling_dv, miscalling_gvrp, accuracy, auc_


# save the result of GVRP refinement model appliance
def confusion_rhe(result_csv_path, confusion_csv_path):
    result_df = pd.read_csv(result_csv_path, index_col=0)

    print(result_df)
    confusion_M = {}
    for i in tqdm(range(len(result_df.index))):
        if result_df.loc[i, 'MATCH'] == 0 and result_df.loc[i, 'result'] == 0:
            flag_ = 'TN'
        elif result_df.loc[i, 'MATCH'] == 0 and result_df.loc[i, 'result'] == 1:
            flag_ = 'FP'
        elif result_df.loc[i, 'MATCH'] == 1 and result_df.loc[i, 'result'] == 0:
            flag_ = 'FN'
        elif result_df.loc[i, 'MATCH'] == 1 and result_df.loc[i, 'result'] == 1:
            flag_ = 'TP'
        else:
            flag_ = None
        confusion_M[i] = flag_

    confusion_df = pd.DataFrame.from_dict(confusion_M, orient='index', columns=['confusion_M'])

    result_df['confusion_M'] = confusion_df['confusion_M']
    result_df.to_csv(confusion_csv_path)
    return confusion_df


def score_rhe(confusion_csv_path, score_csv_path):
    confusion_df = pd.read_csv(confusion_csv_path, index_col=0)

    type_list = ['ht.snp', 'hm.snp']
    column_list = ['TYPE', 'total_#', 'sum of all', '# of TP', '# of TN', '# of FP', '# of FN', 'precision',
                   'recall', 'f1 score', 'miscalling_dv', 'miscalling_gvrp', 'variance rate', 'accuracy', 'AUC-ROC']
    result_list = []
    for type_ in type_list:
        df_t = confusion_df[confusion_df['TYPE'] == type_].reset_index(drop=True)
        tp, tn, fp, fn, precision, recall, f1_score, miscalling_dv, miscalling_gvrp, accuracy, auc_ = return_metric(df_t)
        v_rate = round((miscalling_gvrp / miscalling_dv) * 100, 3)
        print('TYPE :', type_, 'total_# :', len(df_t.index), 'sum of all :', tp + tn + fp + fn,
              '# of TP :', tp, '# of TN :', tn, '# of FP :', fp, '# of FN :', fn,
              'precision :', precision, 'recall :', recall, 'f1 score :', f1_score,
              'miscalling_dv :', miscalling_dv, 'miscalling_gvrp :', miscalling_gvrp,
              'variance rate :', v_rate, 'accuracy :', accuracy, 'AUC-ROC :', auc_)
        t_list = [type_, len(df_t.index), tp + tn + fp + fn, tp, tn, fp, fn, precision, recall,
                  f1_score, miscalling_dv, miscalling_gvrp, v_rate, accuracy, auc_]
        result_list.append(t_list)

    tp, tn, fp, fn, precision, recall, f1_score, miscalling_dv, miscalling_gvrp, accuracy, auc_ = return_metric(confusion_df)
    v_rate = round((miscalling_gvrp / miscalling_dv) * 100, 3)
    print('TYPE : ALL', 'total_# :', len(confusion_df.index), 'sum of all :', tp + tn + fp + fn,
          '# of TP :', tp, '# of TN :', tn, '# of FP :', fp, '# of FN :', fn,
          'precision :', precision, 'recall :', recall, 'f1 score :', f1_score,
          'miscalling_dv :', miscalling_dv, 'miscalling_gvrp :', miscalling_gvrp,
          'variance rate :', v_rate, 'accuracy :', accuracy, 'AUC-ROC :', auc_)
    t_list = ['ALL', len(confusion_df.index), tp + tn + fp + fn, tp, tn, fp, fn, precision, recall,
              f1_score, miscalling_dv, miscalling_gvrp, v_rate, accuracy, auc_]
    result_list.append(t_list)

    df_r = pd.DataFrame(result_list, columns=column_list)
    df_r.to_csv(score_csv_path)
    return df_r


def filter_refcall(df):
    return df[(df['FILTER'] != 'RefCall') & (df['TYPE'].str.contains('.snp'))].reset_index(drop=True)
