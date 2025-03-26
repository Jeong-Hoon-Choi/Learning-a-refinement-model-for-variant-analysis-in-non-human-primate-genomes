from common_module import *


# get alt base ratio distribution of SNPs in each rhesus macaque individuals
if __name__ == '__main__':
    ind_list_all = [ind_list_china, ind_list_30]
    sampling_n = 100
    # v_list = ['hm', 'ht']
    # v_list = ['ht']
    v_list = ['hm']
    for v_type in v_list:
        df_concat = pd.DataFrame()
        if v_type == 'ht':
            alpha = 50
        else:
            alpha = 100
        print(v_type)
        for ind_list in ind_list_all:
            print(ind_list)
            # input data path
            if ind_list == ind_list_china:
                path_result_csv = path_refinement_result + 'chinese_origin/'
                path_pre_csv = path_preprocess_csv + 'chinese_origin/'
                path_pre_seq = path_preprocess_sequence + 'chinese_origin/'
                ind_tag = 'MMUL.CH-'
            else:
                path_result_csv = path_refinement_result + 'indian_origin/'
                path_pre_csv = path_preprocess_csv + 'indian_origin/'
                path_pre_seq = path_preprocess_sequence + 'indian_origin/'
                ind_tag = 'MMUL.IN-'
            for ind in ind_list:
                print(ind)
                ind_name_ = ind_tag + str(ind)
                path_ind = path_alt_base_ratio_result + v_type + '/' + str(ind)
                file_path_GVRP = path_result_csv + 'result_csv_' + str(ind) + '.csv'
                file_path_labeled = path_pre_csv + 'gt_labeled_' + str(ind) + '.csv'
                file_path_raw_dv = path_pre_csv + 'dv_labeled_' + str(ind) + '.csv'
                file_path_dedup = path_pre_seq + str(ind) + '_dedup.bam'
                plot_path_1 = path_ind + '/GVRP_newly_detected.png'
                plot_path_2 = path_ind + '/compare_between_GVRP_n_label.png'
                save_path_label_sampled = path_ind + '/label_' + str(ind) + '_sampled.csv'
                save_path_raw_dv_sampled = path_ind + '/raw_dv_' + str(ind) + '_sampled.csv'
                save_path_GVRP_sampled = path_ind + '/GVRP_' + str(ind) + '_sampled.csv'
                save_path_GVRP_new_sampled = path_ind + '/GVRP_new_' + str(ind) + '_sampled.csv'
                plot_title_ = ind_name_ + ' ' + v_type.upper() + 'SNVs' + ' ABR Density Distribution'

                if not os.path.exists(path_ind):
                    os.makedirs(path_ind)

                print('get dataframe')
                label_df = get_df(file_path_labeled, 'label', v_type)
                raw_dv_df = get_df(file_path_raw_dv, 'raw_dv', v_type)
                GVRP_new_df = get_df(file_path_GVRP, 'filter_sampled', v_type)
                GVRP_df = get_df(file_path_GVRP, 'GVRP', v_type)

                print('get distribution')
                print('part1') # the alt base ratio of GVRP newly detected variants
                dis_list_1 = list()
                dis_list_1.append({'distribution': get_dis(GVRP_new_df, file_path_dedup, save_path_GVRP_new_sampled, sampling_n),
                                   'label': 'DeepVariant-Refinement Model Newly Detected SNVs'})
                dis_list_1.append({'distribution': get_dis(raw_dv_df, file_path_dedup, save_path_raw_dv_sampled, sampling_n),
                                   'label': 'DeepVariant Called Raw SNVs'})

                print('part2') # comparing the alt base ratio between the label and GVRP variatns
                dis_list_2 = list()
                dis_list_2.append({'distribution': get_dis(GVRP_df, file_path_dedup, save_path_GVRP_sampled, sampling_n),
                                   'label': 'After Applying Refinement Model'})
                dis_list_2.append({'distribution': get_dis(label_df, file_path_dedup, save_path_label_sampled, sampling_n),
                                   'label': 'Ground Truth SNVs'})

                print('plot distribution')  # plot KDE distribution based on alt base ratio we got before
                plot_dis(dis_list_1, plot_path_1, plot_title_)
                plot_dis(dis_list_2, plot_path_2, plot_title_)

                print('statistic tests') # conduct the statistical analysis
                # part1_t, part1_p, part1_result = statistic_normal_check(dis_list_1[0]['distribution'], alpha)
                part1_t_t, part1_p_t, part1_result_t, part1_t_w, part1_p_w, part1_result_w \
                    = statistic_dis_compare(dis_list_1[0]['distribution'], dis_list_2[1]['distribution'], is_t_test=True)
                part2_t_t, part2_p_t, part2_result_t, part2_t_w, part2_p_w, part2_result_w\
                    = statistic_dis_compare(dis_list_2[0]['distribution'], dis_list_2[1]['distribution'], is_t_test=True)
                # part3_t_t, part3_p_t, part3_result_t, part3_t_w, part3_p_w, part3_result_w \
                #     = statistic_dis_compare(dis_list_1[0]['distribution'], dis_list_1[1]['distribution'], is_t_test=True)

                ind_dict = {'ind': [ind],}
                GVRP_new_dict = {'gvrp new mean': [np.mean(dis_list_1[0]['distribution'])],
                                 'gvrp new var': [np.var(dis_list_1[0]['distribution'])]}
                GVRP_filtered_dict = {'gvrp filtered mean': [np.mean(dis_list_2[0]['distribution'])],
                                      'gvrp filtered var': [np.var(dis_list_2[0]['distribution'])]}
                label_dict = {'label mean': [np.mean(dis_list_2[1]['distribution'])],
                              'label var': [np.var(dis_list_2[1]['distribution'])]}
                stat_dict = {'t statistics for part1 t-test': [part1_t_t], 'p value for part1 t-test': [part1_p_t], 'result for part1 t-test': [part1_result_t],
                             't statistics for part2 t-test': [part2_t_t], 'p value for part2 t-test': [part2_p_t], 'result for part2 t-test': [part2_result_t]}

                result_dict = ind_dict | GVRP_new_dict | GVRP_filtered_dict | label_dict | stat_dict

                df_r = pd.DataFrame.from_dict(result_dict)
                df_concat = pd.concat([df_concat, df_r]).reset_index(drop=True)

        df_concat.to_csv(path_alt_base_ratio_result + v_type + '_alt_base_result.csv')

    print('done')
