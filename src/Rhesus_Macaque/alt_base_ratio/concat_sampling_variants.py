from common_module import *


if __name__ == '__main__':
    ind_list_all = [ind_list_china, ind_list_30]
    sampling_n = 100
    df_concat = pd.DataFrame()
    for v_type in ['hm', 'ht']:
        total_path = path_alt_base_ratio_result + v_type + '/total/'
        save_path_label_sampled_t = total_path + 'label_total_sampled.csv'
        save_path_GVRP_sampled_t = total_path + 'GVRP_total_sampled.csv'
        save_path_GVRP_new_sampled_t = total_path + 'GVRP_new_total_sampled.csv'
        save_path_raw_dv_sampled_t = total_path + 'raw_dv_sampled.csv'

        if v_type == 'ht':
            alpha = 50
        else:
            alpha = 100
        print(v_type)

        fpd_label_concat = pd.DataFrame()
        fpd_GVRP_concat = pd.DataFrame()
        fpd_GVRP_new_concat = pd.DataFrame()
        fpd_raw_dv_concat = pd.DataFrame()
        for ind_list in ind_list_all:
            print(ind_list)
            if ind_list == ind_list_china:
                path_result_csv = path_refinement_result + 'chinese_origin/'
                path_pre_csv = path_preprocess_csv + 'chinese_origin/'
                path_pre_seq = path_preprocess_sequence + 'chinese_origin/'
                ind_tag = 'MMUL.IN-'
            else:
                path_result_csv = path_refinement_result + 'indian_origin/'
                path_pre_csv = path_preprocess_csv + 'indian_origin/'
                path_pre_seq = path_preprocess_sequence + 'indian_origin/'
                ind_tag = 'MMUL.CH-'
            for ind in ind_list:
                print(ind)
                ind_name_ = ind_tag + str(ind)
                path_ind = path_alt_base_ratio_result + v_type + '/' + str(ind)
                save_path_label_sampled = path_ind + '/label_' + str(ind) + '_sampled.csv'
                save_path_GVRP_sampled = path_ind + '/GVRP_' + str(ind) + '_sampled.csv'
                save_path_GVRP_new_sampled = path_ind + '/GVRP_new_' + str(ind) + '_sampled.csv'
                save_path_raw_dv_sampled = path_ind + '/raw_dv_' + str(ind) + '_sampled.csv'

                if not os.path.exists(path_ind):
                    os.makedirs(path_ind)

                print('get dataframe')

                fpd_label = pd.read_csv(save_path_label_sampled, index_col=0)
                fpd_label['ind_tag'] = ind
                fpd_GVRP = pd.read_csv(save_path_GVRP_sampled, index_col=0)
                fpd_GVRP['ind_tag'] = ind
                fpd_GVRP_new = pd.read_csv(save_path_GVRP_new_sampled, index_col=0)
                fpd_GVRP_new['ind_tag'] = ind
                fpd_raw_dv = pd.read_csv(save_path_raw_dv_sampled, index_col=0)
                fpd_raw_dv['ind_tag'] = ind

                fpd_label_concat = pd.concat([fpd_label_concat, fpd_label]).reset_index(drop=True)
                fpd_GVRP_concat = pd.concat([fpd_GVRP_concat, fpd_GVRP]).reset_index(drop=True)
                fpd_GVRP_new_concat = pd.concat([fpd_GVRP_new_concat, fpd_GVRP_new]).reset_index(drop=True)
                fpd_raw_dv_concat = pd.concat([fpd_raw_dv_concat, fpd_raw_dv]).reset_index(drop=True)

        fpd_label_concat.to_csv(save_path_label_sampled_t)
        fpd_GVRP_concat.to_csv(save_path_GVRP_sampled_t)
        fpd_GVRP_new_concat.to_csv(save_path_GVRP_new_sampled_t)
        fpd_raw_dv_concat.to_csv(save_path_raw_dv_sampled_t)

        fpd_label_dis = sample_fpd(fpd_label_concat, save_path_label_sampled_t, sampling_n)
        fpd_GVRP_dis = sample_fpd(fpd_GVRP_concat, save_path_GVRP_sampled_t, sampling_n)
        fpd_GVRP_new_dis = sample_fpd(fpd_GVRP_new_concat, save_path_GVRP_new_sampled_t, sampling_n)
        fpd_raw_dv_dis = sample_fpd(fpd_raw_dv_concat, save_path_raw_dv_sampled_t, sampling_n)
