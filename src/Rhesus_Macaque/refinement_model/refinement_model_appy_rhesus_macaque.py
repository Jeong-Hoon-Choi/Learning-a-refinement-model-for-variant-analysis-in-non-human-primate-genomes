import pandas as pd

from preprocess_variant_call_rhesus_macaque import *


if __name__ == '__main__':
    df_dict = {}
    for ind_list in ind_list_all:
        if ind_list == ind_list_china:
            origin = 'Chinese'
        else:
            origin = 'Indian'
        for ind in ind_list:
            print('\n', ind)

            gt_path, dv_path, gt_csv_path, dv_csv_path, match_csv_path, \
                learning_csv_path, result_csv_path, confusion_csv_path,\
                score_csv_path, alignment_path, alignment_csv_path = get_data_path(origin, ind)

            labeling_M('gt', ind, gt_path, dv_path, gt_csv_path, dv_csv_path, ind_list)
            labeling_M('dv', ind, gt_path, dv_path, gt_csv_path, dv_csv_path, ind_list)
            match_L(gt_csv_path, dv_csv_path, match_csv_path)
            add_align(match_csv_path, alignment_path, alignment_csv_path)
            rhe_L(alignment_csv_path, learning_csv_path, result_csv_path, model_dir)
            confusion_rhe(result_csv_path, confusion_csv_path)
            score_df = score_rhe(confusion_csv_path, score_csv_path)

            # score_df = pd.read_csv(score_csv_path, index_col=0)
            df_dict[ind] = {'ht snp': score_df[score_df['TYPE'] == 'ht.snp'].drop(['TYPE'], axis=1),
                            'hm snp': score_df[score_df['TYPE'] == 'hm.snp'].drop(['TYPE'], axis=1),
                            'all snp': score_df[score_df['TYPE'] == 'ALL'].drop(['TYPE'], axis=1)}

    concat_ht = pd.DataFrame()
    concat_hm = pd.DataFrame()
    concat_all = pd.DataFrame()
    for key_ in df_dict:
        concat_ht = pd.concat([concat_ht, df_dict[key_]['ht snp']]).reset_index(drop=True)
        concat_hm = pd.concat([concat_hm, df_dict[key_]['hm snp']]).reset_index(drop=True)
        concat_all = pd.concat([concat_all, df_dict[key_]['all snp']]).reset_index(drop=True)

    concat_ht['ind'] = ind_list_china + ind_list_indian
    concat_hm['ind'] = ind_list_china + ind_list_indian
    concat_all['ind'] = ind_list_china + ind_list_indian

    concat_ht.to_csv(path_result_data + 'csv/' + 'concat_ht_score.csv')
    concat_hm.to_csv(path_result_data + 'csv/' + 'concat_hm_score.csv')
    concat_all.to_csv(path_result_data + 'csv/' + 'concat_all_score.csv')
