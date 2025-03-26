from common_parameter import *

df_m = {'concat_data_train_HG001_test_HG002': 'HG002_concat_match_all_',
        'concat_data_train_HG002_test_HG001': 'HG001_concat_match_all_',
        'mix_all': 'mixed_match_test'}

dict_feature = {'_all_info': ['AD_diff', 'DP_mean', 'GQ_mean', 'PL_diff', 'VAF_mean', 'QUAL', 'M_ratio',
                              'S_ratio', 'Total_reads', 'mean_mapq', 'low_mapq_ratio', 'forward_strand_ratio']}

list_type = ['.snp', '.indel', 'all']

fptp_dict = dict()
fptp_df = pd.DataFrame()

for data_ in df_m:
    print(data_)
    for v_ in dict_feature:
        print(v_)

        df_confusion = pd.read_csv(performance_dir + df_m[data_] + v_ + '_confusion.csv', index_col=0)

        for v_type in list_type:
            if v_type == '.snp':
                df_ = df_confusion[df_confusion['TYPE'].str.contains('.snp')].reset_index(drop=True)
            elif v_type == '.indel':
                df_ = df_confusion[df_confusion['TYPE'].str.contains('.indel')].reset_index(drop=True)
            else:
                df_ = df_confusion


            s_ = data_ + '_' + v_
            fptp_dict[s_] = {'data_type': [], 'v_type': [],
                             'FP_DV': [], 'TP_DV': [], 'FP_GVRP': [], 'TP_GVRP': [],
                             'filtered_FP': [], 'filtered_TP': [],
                             'FP_filter_rate': [], 'TP_filter_rate': [],
                             'TP_rate_DV': [], 'TP_rate_GVRP': [], 'Enhance_rate': []}

            fptp_dict[s_]['data_type'].append(data_)
            fptp_dict[s_]['v_type'].append(v_type)

            fp_dv = len(df_[df_['label'] == 0].reset_index(drop=True).index)
            tp_dv = len(df_[df_['label'] == 1].reset_index(drop=True).index)
            fptp_dict[s_]['FP_DV'].append(fp_dv)
            fptp_dict[s_]['TP_DV'].append(tp_dv)

            df_GVRP = df_[(df_['confusion_M'] == 'TP') | (df_['confusion_M'] == 'FP')].reset_index(drop=True)

            fp_gvrp = len(df_GVRP[df_GVRP['label'] == 0].reset_index(drop=True).index)
            tp_gvrp = len(df_GVRP[df_GVRP['label'] == 1].reset_index(drop=True).index)
            fptp_dict[s_]['FP_GVRP'].append(fp_gvrp)
            fptp_dict[s_]['TP_GVRP'].append(tp_gvrp)

            filtered_fp = fp_dv - fp_gvrp
            filtered_tp = tp_dv - tp_gvrp
            fptp_dict[s_]['filtered_FP'].append(filtered_fp)
            fptp_dict[s_]['filtered_TP'].append(filtered_tp)

            fp_filter_rate = round(filtered_fp / fp_dv, 3)
            tp_filter_rate = round(filtered_tp / tp_dv, 3)
            fptp_dict[s_]['FP_filter_rate'].append(fp_filter_rate)
            fptp_dict[s_]['TP_filter_rate'].append(tp_filter_rate)

            tp_rate_dv = round(tp_dv / (fp_dv + tp_dv), 3)
            tp_rate_gvrp = round(tp_gvrp / (fp_gvrp + tp_gvrp), 3)
            enhance_rate = round(tp_rate_gvrp / tp_rate_dv, 3)
            fptp_dict[s_]['TP_rate_DV'].append(tp_rate_dv)
            fptp_dict[s_]['TP_rate_GVRP'].append(tp_rate_gvrp)
            fptp_dict[s_]['Enhance_rate'].append(enhance_rate)

            print(fptp_dict[s_])

            temp_df = pd.DataFrame.from_dict(fptp_dict[s_])
            fptp_df = pd.concat([fptp_df, temp_df]).reset_index(drop=True)

fptp_df.to_csv(filtering_dir + 'filtering_result_0313.csv')
