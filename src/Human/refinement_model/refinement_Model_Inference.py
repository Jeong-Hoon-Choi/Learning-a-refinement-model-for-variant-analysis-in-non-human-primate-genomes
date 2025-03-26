from common_module import *


if __name__ == '__main__':
    for m in mode_all:
        print('-----------------------------------------------------------------------------')
        print('mode :', m)

        result_dict_ = {'f1': [], 'precision': [], 'recall': [], 'Accuracy':[]}
        pred_dict_ = {}

        # models = ['XG_Boost', 'LGBM', 'RF', 'LR', 'KNN', 'NB', 'MLP', 'ft_transformer']
        models = ['ft_transformer']
        # models = ['XG_Boost']

        print('1. load train test data')
        train_X, train_y, test_X, test_y = load_data(m, data_dir, 'without_types')

        tag = '_only_likelihood'
        feature_less = ['QUAL', 'GQ_mean', 'PL_diff']

        # tag = '_add_DeepVariant_alignment_info'
        # feature_less = ['AD_diff', 'DP_mean', 'GQ_mean', 'PL_diff', 'VAF_mean', 'QUAL']

        # tag = '_add_all_alignment_info'
        # feature_less = ['AD_diff', 'DP_mean', 'GQ_mean', 'PL_diff', 'VAF_mean', 'QUAL', 'M_ratio', 'S_ratio', 'Total_reads', 'mean_mapq', 'low_mapq_ratio']

        train_X = train_X[feature_less]
        test_X = test_X[feature_less]

        print(train_X.columns)
        print(tag)

        result_dict = {'f1': [], 'precision': [], 'recall': [], 'Accuracy':[], 'TNR': [], 'ROC-AUC': []}
        print('2. model inference performance')
        for i, model_name in enumerate(models):
            print('2.' + str(i + 1) + ' | model:', model_name)
            if model_name == 'MLP':
                file_name = model_dir + m + '/' + tag + '_' + model_name + '_model.h5'
                model = tf.keras.models.load_model(file_name)
                predict_ = model.predict(test_X)
                predict_ = np.where(predict_ > 0.5, 1, 0)
                predict_pro_ = model.predict(test_X)
            elif model_name == 'ft_transformer':
                file_name = model_dir + m + '/' + tag + '_' + model_name + '_model.pth'
                predict_pro_a, predict_ = ft_transformer_inf(file_name, test_X, test_y)
                predict_pro_ = [pre[1] for pre in predict_pro_a]
            else:
                file_name = model_dir + m + '/' + tag + '_' + model_name + '_model.pkl'
                model = pickle.load(open(file_name, 'rb'))
                predict_ = model.predict(test_X)
                predict_pro_ = model.predict_proba(test_X)[:, 1]
            score_f(test_y, predict_, result_dict)

            fpr, tpr, threshold_ = roc_curve(test_y, predict_pro_)
            auc_ = auc(fpr, tpr)
            result_dict['ROC-AUC'].append(auc_)

        label_df = pd.DataFrame.from_dict(result_dict, orient='index', columns=models)
        label_df = label_df.transpose()
        print()
        print('tag :', tag)
        print('#train :', len(train_y), '#test :', len(test_y), '\n')
        print(tabulate(label_df, headers='keys', tablefmt='psql'))
        label_df.to_csv(performance_dir + m + '_' + tag + '_result.csv')
