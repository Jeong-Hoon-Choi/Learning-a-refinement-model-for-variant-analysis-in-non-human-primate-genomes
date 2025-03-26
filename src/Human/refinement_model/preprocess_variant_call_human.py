import pandas as pd

from common_module import *


# vcf to csv
def read_vcf(path):
    with open(path, 'r') as f:
        lines_info = [l for l in tqdm(f) if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines_info)),
        dtype={'#CHROM': str, 'POS': str, 'ID': str, 'REF': str, 'ALT': str, 'QUAL': float, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})


# setting csv path for preprocessing data
def get_parameters(data_dir, output_dir, target_data, target_type):
    if target_data == 'HG001':
        label_path = data_dir + 'HG001_005-6_dedup_WGS.deepvar.giab_benchmark.GRCh38_concordance.vcf'
        data_path = data_dir + 'HG001_005-6_dedup_WGS.deepvar.GRCh38.PASS.' + target_type + '.vcf'
        match_data_path = output_dir + target_data + target_type + '_match.csv'
        align_data_path = output_dir + target_data + target_type + '_align.csv'
        learn_data_path = output_dir + target_data + target_type + '.csv'
        column_name = 'HG001_005'

    else:
        label_path = data_dir + 'HG002.giab_0028-9_dedup_WGS.deepvar.giab_benchmark.GRCh38_concordance.vcf'
        data_path = data_dir + 'HG002.giab_0028-9_WGS.GRCh38.deepvar.PASS.' + target_type + '.vcf'
        match_data_path = output_dir + target_data + target_type + '_match.csv'
        align_data_path = output_dir + target_data + target_type + '_align.csv'
        learn_data_path = output_dir + target_data + target_type + '.csv'
        column_name = 'HG002.giab_0029'

    return label_path, data_path, match_data_path, learn_data_path, align_data_path, column_name


# labeling the human variant based on GIAB data
def labeling_vcf(label_path, data_path, match_data_path, target_type):
    true_list = ['CONC_ST=TP', 'CONC_ST=TP,TP', 'CONC_ST=TP,TN', 'CONC_ST=TN,TP']
    non_list = ['CONC_ST=EMPTY']
    chrom_values = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']

    label_key_dict = {}
    vcf_df = read_vcf(label_path)
    vcf_df = vcf_df[vcf_df['CHROM'].isin(chrom_values)]
    vcf_df = vcf_df.reset_index(drop=True)

    print('convert vcf file to csv file')
    for i in tqdm(range(len(vcf_df.index))):
        key_ = str(vcf_df.loc[i, 'CHROM']) + '_' + str(vcf_df.loc[i, 'POS'])
        if vcf_df.loc[i, 'INFO'] in non_list:
            label_key_dict[key_] = 'NONE'
        elif vcf_df.loc[i, 'INFO'] in true_list:
            label_key_dict[key_] = 1
        else:
            label_key_dict[key_] = 0

    label_dict = {}
    data_df = read_vcf(data_path)
    data_df = data_df[data_df['CHROM'].isin(chrom_values)]
    data_df = data_df.reset_index(drop=True)
    print('labeling data')
    for i in tqdm(range(len(data_df.index))):
        key_ = str(data_df.loc[i, 'CHROM']) + '_' + str(data_df.loc[i, 'POS'])
        label_dict[i] = label_key_dict[key_]

    label_df = pd.DataFrame.from_dict(label_dict, orient='index', columns=['label'])
    data_df['TYPE'] = target_type
    data_df['label'] = label_df['label']
    data_df = data_df.loc[(data_df['label'] != 'NONE') & (data_df['label'] != 'MISS')]
    data_df = data_df.reset_index(drop=True)
    data_df.to_csv(match_data_path)
    print(data_df)
    return data_df


def extract_alignment_info(bam, chrom, pos):
    cigar_features = []
    mapq_values = []
    strand_counts = {"forward": 0, "reverse": 0}  # Strand count 초기화

    for read in bam.fetch(str(chrom), int(pos) - 1, int(pos)):  # 0-based 좌표
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


def get_alignment_information(seq_dir, labeled_df, align_data_path):
    align_info_results = []
    chroms = labeled_df['CHROM'].values
    positions = labeled_df['POS'].values

    with pysam.AlignmentFile(seq_dir, "rb") as bam:
        for i in tqdm(range(len(labeled_df.index))):
            chrom = chroms[i]
            pos = positions[i]
            align_info = extract_alignment_info(bam, chrom, pos)
            align_info_results.append({"CHROM": chrom, "POS": pos, **align_info})

    # align_info_results를 DataFrame으로 변환
    align_info_df = pd.DataFrame(align_info_results)

    # labeled_df와 병합
    align_df = labeled_df.merge(align_info_df, on=["CHROM", "POS"], how="left")
    align_df.fillna({
        "M_ratio": 0, "S_ratio": 0, "Total_reads": 0,
        "mean_mapq": 0, "low_mapq_ratio": 0,
        "forward_strand_ratio": 0, "reverse_strand_ratio": 0
    }, inplace=True)
    align_df.to_csv(align_data_path)
    return align_df


# make learning data based on preprocessing labeling data
def make_learning_data(data_csv, column_name, learn_data_path):
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
    feature_columns['label'] = data_csv['label']

    feature_columns.to_csv(learn_data_path, index=False)
    return feature_columns


if __name__ == '__main__':
    # data_list = ['HG001', 'HG002']
    data_list = ['HG001']
    # data_list = ['HG002']
    type_list = ['hm.indel', 'hm.snp', 'ht.indel', 'ht.snp']

    for target_data in data_list:
        print(target_data)
        concat_df = pd.DataFrame()
        concat_match_df = pd.DataFrame()
        for target_type in type_list:
            print('---------------------------------------------------------')
            print(target_type)

            # 1. first set parameters
            print('1. parameter setting')
            label_path, data_path, match_data_path, learn_data_path, align_data_path, column_name = get_parameters(vcf_dir, data_dir, target_data, target_type)

            # 2. labeling data
            print('2. labeling data')
            # labeled_df = labeling_vcf(label_path, data_path, match_data_path, target_type)
            labeled_df = pd.read_csv(match_data_path, index_col=0)
            print(labeled_df)

            # 3. get alignment information
            print('3. get alignment information')
            # align_df = get_alignment_information(seq_dir, labeled_df, align_data_path)
            align_df = pd.read_csv(align_data_path, index_col=0)
            # print(align_df)

            # 3-1. concat alignment information data
            concat_match_df = pd.concat([concat_match_df, align_df]).reset_index(drop=True)

            # 4. convert to input csv data
            print('3. convert to input csv data')
            input_csv = make_learning_data(align_df, column_name, learn_data_path)

            # 4-1. concat learning data
            concat_df = pd.concat([concat_df, input_csv]).reset_index(drop=True)

        concat_df.to_csv(data_dir + target_data + '_concat_all_.csv')
        concat_match_df.to_csv(data_dir + target_data + '_concat_match_all_.csv')
