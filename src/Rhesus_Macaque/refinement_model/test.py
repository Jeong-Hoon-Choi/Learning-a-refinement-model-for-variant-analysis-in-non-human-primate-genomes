from common_parameter import *
from preprocess_variant_call_rhesus_macaque import *


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

        df_l = pd.read_csv(gt_csv_path)
        print(len(df_l[df_l['FILTER'] == 'RefCall'].reset_index(drop=True).index))
