from common_module import *

title_dict = {'mix_all': 'HG001, HG002 Mixed Data',
              'concat_data_train_HG001_test_HG002': 'Train HG001, Test HG002',
              'concat_data_train_HG002_test_HG001': 'Train HG002, Test HG001'}

if __name__ == '__main__':
    for m in mode_all:
        print('-----------------------------------------------------------------------------')
        print('mode :', m)

        result_dict_ = {'f1': [], 'precision': [], 'recall': [], 'Accuracy':[]}
        models = ['LGBM', 'XG_Boost', 'RF', 'LR', 'KNN', 'NB', 'MLP', 'ft_transformer']
        # models = ['LGBM']
        m_dict_ = {}
        l_dict_ = {'label': []}
        for mo in models:
            m_dict_[mo] = list()

        print('1. load train test data')
        train_X, train_y, test_X, test_y = load_data(m, data_dir, 'without_types')

        tag = '_all_info'
        feature_less = ['AD_diff', 'DP_mean', 'GQ_mean', 'PL_diff', 'VAF_mean', 'QUAL', 'M_ratio', 'S_ratio',
                        'Total_reads', 'mean_mapq', 'low_mapq_ratio', 'forward_strand_ratio']

        train_X = train_X[feature_less]
        test_X = test_X[feature_less]

        print('2. model inference ROC')
        pred_dict = {}
        for i, model_name in enumerate(models):
            print('2.' + str(i + 1) + ' | model:', model_name)
            if model_name == 'MLP':
                file_name = model_dir + m + '/' + tag + '_' + model_name + '_model.h5'
                model = tf.keras.models.load_model(file_name)
                predict_pro_ = model.predict(test_X)
            elif model_name == 'ft_transformer':
                file_name = model_dir + m + '/' + tag + '_' + model_name + '_model.pth'
                predict_pro_, predict_ = ft_transformer_inf(file_name, test_X, test_y)
                predict_pro_ = np.array(predict_pro_)[:, 1]
            else:
                file_name = model_dir + m + '/' + tag + '_' + model_name + '_model.pkl'
                model = pickle.load(open(file_name, 'rb'))
                predict_pro_ = model.predict_proba(test_X)[:, 1]

            fpr, tpr, threshold_ = roc_curve(test_y, predict_pro_)
            auc_ = auc(fpr, tpr)
            pred_dict[model_name] = {'fpr': fpr, 'tpr': tpr, 'auc_': auc_}

        print(m)
        for m_ in pred_dict:
            print(m_, pred_dict[m_]['auc_'])

        plt.figure(figsize=(8, 8))
        for m_ in pred_dict:
            plt.plot(pred_dict[m_]['fpr'], pred_dict[m_]['tpr'], lw=2,
                     label=m_ + ' (AUC = %0.4f)' % pred_dict[m_]['auc_'])
        plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate', fontsize=16)
        plt.ylabel('True Positive Rate', fontsize=16)
        plt.title('ROC curve in ' + title_dict[m], fontsize=18)
        plt.legend(loc="lower right", fontsize=12)
        # plt.savefig(roc_dir + m + '_roc_curve.png', dpi=600)
