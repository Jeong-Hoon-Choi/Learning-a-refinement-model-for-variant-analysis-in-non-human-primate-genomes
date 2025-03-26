from common_parameter import *

fptp_dict = dict()
df_dict = {'hm': pd.DataFrame(), 'ht': pd.DataFrame(), 'all': pd.DataFrame()}

for index_list in ind_list_all:
    for i_ in index_list:
        if i_ in ind_list_china:
            GVRP_result_p = path_result_data + 'csv/chinese_origin/confusion_' + str(i_) +'.csv'
        else:
            GVRP_result_p = path_result_data + 'csv/indian_origin/confusion_' + str(i_) + '.csv'

        df_p = pd.read_csv(GVRP_result_p, index_col=0)

        for l_ in ['hm', 'ht', 'all']:
            s_ = str(i_) + '_' + l_
            fptp_dict[s_] = {'IND': [], 'v_type': [], 'FP_DV': [],
                             'TP_DV': [], 'FP_GVRP': [], 'TP_GVRP': [],
                             'filtered_FP': [], 'filtered_TP': [],
                             'FP_filter_rate': [], 'TP_filter_rate': [],
                             'MR_DV': [], 'MR_GVRP': [], 'decrease_rate': []}

            print(s_)
            if l_ != 'all':
                df_ = df_p[df_p['TYPE'].str.contains(l_ + '.snp')].reset_index(drop=True)
            else:
                df_ = df_p

            fptp_dict[s_]['IND'].append(i_)
            fptp_dict[s_]['v_type'].append(l_)

            fp_dv = len(df_[df_['MATCH'] == 0].reset_index(drop=True).index)
            tp_dv = len(df_[df_['MATCH'] == 1].reset_index(drop=True).index)
            fptp_dict[s_]['FP_DV'].append(fp_dv)
            fptp_dict[s_]['TP_DV'].append(tp_dv)

            df_GVRP = df_[(df_['confusion_M'] == 'TP') | (df_['confusion_M'] == 'FP')].reset_index(drop=True)

            fp_gvrp = len(df_GVRP[df_GVRP['MATCH'] == 0].reset_index(drop=True).index)
            tp_gvrp = len(df_GVRP[df_GVRP['MATCH'] == 1].reset_index(drop=True).index)
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

            mr_dv = round(fp_dv / (fp_dv + tp_dv), 3)
            mr_gvrp = round(fp_gvrp / (fp_gvrp + tp_gvrp), 3)
            decrease_rate = round(mr_gvrp / mr_dv, 3)
            fptp_dict[s_]['MR_DV'].append(mr_dv)
            fptp_dict[s_]['MR_GVRP'].append(mr_gvrp)
            fptp_dict[s_]['decrease_rate'].append(decrease_rate)

            print(fptp_dict[s_])

            temp_df = pd.DataFrame.from_dict(fptp_dict[s_])
            df_dict[l_] = pd.concat([df_dict[l_], temp_df]).reset_index(drop=True)

for l_ in df_dict:
    df_dict[l_].to_csv(path_result_data + 'filtering/filtering_result_' + l_ + '.csv')
