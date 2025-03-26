from common_module import *


# get alt base ratio distribution of SNPs in each rhesus macaque individuals
if __name__ == '__main__':
    ind_list_all = [ind_list_china, ind_list_30]
    sampling_n = 100
    for v_type in ['hm', 'ht']:
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

                plot_path_1_ = path_alt_base_ratio_result + v_type + '/plot/' + 'compared_between_rm-dv ' + ind_name_ + '_' + v_type + '.png'
                plot_path_2_ = path_alt_base_ratio_result + v_type + '/plot/' + 'compared_between_rm_gt ' + ind_name_ + '_' + v_type + '.png'

                save_path_label_sampled = path_ind + '/label_' + str(ind) + '_sampled.csv'
                save_path_raw_dv_sampled = path_ind + '/raw_dv_' + str(ind) + '_sampled.csv'
                save_path_GVRP_sampled = path_ind + '/GVRP_' + str(ind) + '_sampled.csv'
                save_path_GVRP_new_sampled = path_ind + '/GVRP_new_' + str(ind) + '_sampled.csv'
                plot_title_ = ind_name_ + ' ' + v_type.upper() + '-SNVs' + ' ABR Density Distribution'

                if not os.path.exists(path_ind):
                    os.makedirs(path_ind)

                print('get dataframe')
                fpd_label_concat = pd.read_csv(save_path_label_sampled)
                fpd_raw_dv_concat = pd.read_csv(save_path_raw_dv_sampled)
                fpd_GVRP_concat = pd.read_csv(save_path_GVRP_sampled)
                fpd_GVRP_new_concat = pd.read_csv(save_path_GVRP_new_sampled)

                fpd_label_dis = fpd_label_concat['distribution'].to_list()
                fpd_raw_dv_dis = fpd_raw_dv_concat['distribution'].to_list()
                fpd_GVRP_dis = fpd_GVRP_concat['distribution'].to_list()
                fpd_GVRP_new_dis = fpd_GVRP_new_concat['distribution'].to_list()

                print('get distribution')
                print('part1')
                dis_list_1 = list()
                dis_list_1.append({'distribution': fpd_GVRP_new_dis, 'label': 'DeepVariant-Refinement Model Newly Detected SNVs'})
                dis_list_1.append({'distribution': fpd_raw_dv_dis, 'label': 'DeepVariant Called Raw SNVs'})

                print('part2')
                dis_list_2 = list()
                dis_list_2.append({'distribution': fpd_GVRP_dis, 'label': 'After Applying Refinement Model'})
                dis_list_2.append({'distribution': fpd_label_dis, 'label': 'Ground Truth SNVs'})

                print('plot distribution')  # plot KDE distribution based on alt base ratio we got before
                plot_dis(dis_list_1, plot_path_1_, plot_title_)
                plot_dis(dis_list_2, plot_path_2_, plot_title_)

                print('statistic tests') # conduct the statistical analysis
                part1_t_t, part1_p_t, part1_result_t, part1_t_w, part1_p_w, part1_result_w \
                    = statistic_dis_compare(dis_list_1[0]['distribution'], dis_list_2[1]['distribution'], is_t_test=True)
                part2_t_t, part2_p_t, part2_result_t, part2_t_w, part2_p_w, part2_result_w \
                    = statistic_dis_compare(dis_list_2[0]['distribution'], dis_list_2[1]['distribution'], is_t_test=True)

                ind_dict = {'ind': [ind_name_]}
                GVRP_new_dict = {'gvrp new mean': [str(round(np.mean(dis_list_1[0]['distribution']), 3)) + '/' +\
                                                   str(round(np.var(dis_list_1[0]['distribution']), 3))]}
                raw_dv_dict = {'raw dv mean': [str(round(np.mean(dis_list_1[1]['distribution']), 3)) + '/' +\
                                                   str(round(np.var(dis_list_1[1]['distribution']), 3))]}
                GVRP_filtered_dict = {'gvrp filtered mean': [str(round(np.mean(dis_list_2[0]['distribution']), 3)) + '/' +\
                                                   str(round(np.var(dis_list_2[0]['distribution']), 3))]}
                label_dict = {'label mean': [str(round(np.mean(dis_list_2[1]['distribution']), 3)) + '/' +\
                                                   str(round(np.var(dis_list_2[1]['distribution']), 3))]}
                stat_dict = {'p value for part1 t-test': [part1_p_t], 'result for part1 t-test': [part1_result_t],
                             'p value for part2 t-test': [part2_p_t], 'result for part2 t-test': [part2_result_t]}

                result_dict = ind_dict | GVRP_new_dict | raw_dv_dict | GVRP_filtered_dict | label_dict | stat_dict

                df_r = pd.DataFrame.from_dict(result_dict)
                df_concat = pd.concat([df_concat, df_r]).reset_index(drop=True)

        df_concat.to_csv(path_alt_base_ratio_result + v_type + '_alt_base_result.csv')

    print('done')
