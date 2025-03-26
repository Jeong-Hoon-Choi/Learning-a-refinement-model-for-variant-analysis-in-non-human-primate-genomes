from common_module import *


def confusion_(result_csv_path, confusion_csv_path):
    result_df = pd.read_csv(result_csv_path, index_col=0)

    confusion_M = {}
    for i in tqdm(range(len(result_df.index))):
        if result_df.loc[i, 'label'] == 0 and result_df.loc[i, 'result'] == 0:
            flag_ = 'TN'
        elif result_df.loc[i, 'label'] == 0 and result_df.loc[i, 'result'] == 1:
            flag_ = 'FP'
        elif result_df.loc[i, 'label'] == 1 and result_df.loc[i, 'result'] == 0:
            flag_ = 'FN'
        elif result_df.loc[i, 'label'] == 1 and result_df.loc[i, 'result'] == 1:
            flag_ = 'TP'
        else:
            flag_ = None
        confusion_M[i] = flag_

    confusion_df = pd.DataFrame.from_dict(confusion_M, orient='index', columns=['confusion_M'])

    result_df['confusion_M'] = confusion_df['confusion_M']
    result_df.to_csv(confusion_csv_path)
    return result_df


def return_metric(df_t):
    tp = len(df_t[df_t['confusion_M'] == 'TP'].reset_index(drop=True).index)
    tn = len(df_t[df_t['confusion_M'] == 'TN'].reset_index(drop=True).index)
    fp = len(df_t[df_t['confusion_M'] == 'FP'].reset_index(drop=True).index)
    fn = len(df_t[df_t['confusion_M'] == 'FN'].reset_index(drop=True).index)
    gt_1 = len(df_t[df_t['label'] == 1].reset_index(drop=True).index)
    gt_0 = len(df_t[df_t['label'] == 0].reset_index(drop=True).index)
    y_ = df_t['label']
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


def score_rhe(confusion_csv_path, score_csv_path):
    confusion_df = pd.read_csv(confusion_csv_path, index_col=0)

    type_list = ['snp', 'indel']
    column_list = ['TYPE', 'total_#', 'sum of all', '# of TP', '# of TN', '# of FP', '# of FN', 'precision',
                   'recall', 'f1 score', 'miscalling_dv', 'miscalling_gvrp', 'variance rate', 'accuracy', 'AUC-ROC']
    result_list = []
    for type_ in type_list:
        df_t = confusion_df[confusion_df['TYPE'].str.contains(type_)].reset_index(drop=True)
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


if __name__ == '__main__':
    for m in mode_all:
        print('-----------------------------------------------------------------------------')
        print('mode :', m)

        result_dict_ = {'f1': [], 'precision': [], 'recall': [], 'Accuracy':[]}
        pred_dict_ = {}

        model_name = 'LGBM'

        print('1. load train test data')
        train_X, train_y, test_X, test_y = load_data(m, data_dir, 'without_types')

        dict_feature = {'_all_info': ['AD_diff', 'DP_mean', 'GQ_mean', 'PL_diff', 'VAF_mean', 'QUAL', 'M_ratio',
                                      'S_ratio', 'Total_reads', 'mean_mapq', 'low_mapq_ratio', 'forward_strand_ratio']}

        df_m = {'concat_data_train_HG001_test_HG002': 'HG002_concat_match_all_',
                'concat_data_train_HG002_test_HG001': 'HG001_concat_match_all_',
                'mix_all': 'mixed_match_test'}

        for tag in dict_feature:
            match_df = pd.read_csv(data_dir + df_m[m] + '.csv', index_col=0)

            feature_less = dict_feature[tag]

            train_X_ = train_X[feature_less]
            test_X_ = test_X[feature_less]

            print(tag)
            print(train_X_.columns)

            print('2. model inference performance')
            file_name = model_dir + m + '/' + tag + '_' + model_name + '_model.pkl'
            model = pickle.load(open(file_name, 'rb'))
            predict_ = model.predict(test_X_)
            predict_pro_ = model.predict_proba(test_X_)[:, 1]

            match_df['result'] = predict_
            match_df['result_p'] = predict_pro_

            print('match')
            match_df.to_csv(performance_dir + df_m[m] + tag + '.csv')

            print('confusion')
            con_df = confusion_(performance_dir + df_m[m] + tag + '.csv', performance_dir + df_m[m] + tag + '_confusion.csv')

            print('score')
            sc_df = score_rhe(performance_dir + df_m[m] + tag + '_confusion.csv', performance_dir + df_m[m] + tag + '_score.csv')
