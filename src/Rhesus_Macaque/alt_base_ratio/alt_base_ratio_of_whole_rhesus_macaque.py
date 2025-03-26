from common_module import *


# get the alt base ratio distribution of SNPs in whole rhesus macaque
# For this process you need to concat all the sampled variants first by 'concat_sampling_variants.py'
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
        plot_path_1 = total_path + 'GVRP_newly_detected.png'
        plot_path_2 = total_path + 'compare_between_GVRP_n_label.png'
        plot_title_ = v_type.upper() + '-' + 'SNVs' + ' ABR Density Distribution'

        if v_type == 'ht':
            alpha = 50
        else:
            alpha = 100
        print(v_type)
        fpd_label_concat = pd.read_csv(save_path_label_sampled_t)
        fpd_GVRP_concat = pd.read_csv(save_path_GVRP_sampled_t)
        fpd_GVRP_new_concat = pd.read_csv(save_path_GVRP_new_sampled_t)
        fpd_raw_dv_concat = pd.read_csv(save_path_raw_dv_sampled_t)

        fpd_label_dis = fpd_label_concat['distribution'].to_list()
        fpd_GVRP_dis = fpd_GVRP_concat['distribution'].to_list()
        fpd_GVRP_new_dis = fpd_GVRP_new_concat['distribution'].to_list()
        fpd_raw_dv_dis = fpd_raw_dv_concat['distribution'].to_list()

        print('get distribution')
        print('part1')
        dis_list_1 = list()
        dis_list_1.append({'distribution': fpd_GVRP_new_dis, 'label': 'DeepVariant-Refinement Model Newly Detected SNVs'})
        dis_list_1.append({'distribution': fpd_raw_dv_dis, 'label': 'DeepVariant Called Raw SNVs'})

        print('part2')
        dis_list_2 = list()
        dis_list_2.append({'distribution': fpd_GVRP_dis, 'label': 'After Applying Refinement Model'})
        dis_list_2.append({'distribution': fpd_label_dis, 'label': 'Ground Truth SNVs'})

        print('plot distribution')
        plot_dis(dis_list_1, plot_path_1, plot_title_)
        plot_dis(dis_list_2, plot_path_2, plot_title_)

        print('statistic tests')
        # part1_t, part1_p, part1_result = statistic_normal_check(dis_list_1[0]['distribution'], alpha)
        part1_t_t, part1_p_t, part1_result_t, part1_t_w, part1_p_w, part1_result_w \
            = statistic_dis_compare(dis_list_1[0]['distribution'], dis_list_2[1]['distribution'], is_t_test=True)
        part2_t_t, part2_p_t, part2_result_t, part2_t_w, part2_p_w, part2_result_w \
            = statistic_dis_compare(dis_list_2[0]['distribution'], dis_list_2[1]['distribution'], is_t_test=True)

        # part3_t_t, part3_p_t, part3_result_t, part3_t_w, part3_p_w, part3_result_w \
        #     = statistic_dis_compare(dis_list_1[0]['distribution'], dis_list_1[1]['distribution'], is_t_test=True)

        ind_dict_ = {'v_type': [v_type], }
        GVRP_new_dict = {'gvrp new mean': [np.mean(dis_list_1[0]['distribution'])],
                         'gvrp new var': [np.var(dis_list_1[0]['distribution'])]}
        raw_dv_dict = {'raw dv mean': [np.mean(dis_list_1[1]['distribution'])],
                      'raw dv var': [np.var(dis_list_1[1]['distribution'])]}
        GVRP_filtered_dict = {'gvrp filtered mean': [np.mean(dis_list_2[0]['distribution'])],
                              'gvrp filtered var': [np.var(dis_list_2[0]['distribution'])]}
        label_dict = {'label mean': [np.mean(dis_list_2[1]['distribution'])],
                      'label var': [np.var(dis_list_2[1]['distribution'])]}
        stat_dict = {'t statistics for part1 t-test': [part1_t_t], 'p value for part1 t-test': [part1_p_t], 'result for part1 t-test': [part1_result_t],
                     't statistics for part2 t-test': [part2_t_t], 'p value for part2 t-test': [part2_p_t], 'result for part2 t-test': [part2_result_t]}

        result_dict = ind_dict_ | raw_dv_dict | GVRP_new_dict | GVRP_filtered_dict | label_dict | stat_dict

        df_r = pd.DataFrame.from_dict(result_dict)
        df_concat = pd.concat([df_concat, df_r]).reset_index(drop=True)

    df_concat.to_csv(path_alt_base_ratio_result + 'total_alt_base_result.csv')
