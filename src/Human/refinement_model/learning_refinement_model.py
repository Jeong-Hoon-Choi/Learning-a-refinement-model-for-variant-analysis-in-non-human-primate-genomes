from common_module import *


if __name__ == '__main__':
    for m in mode_all:
        print('-----------------------------------------------------------------------------')
        print('mode :', m)

        result_dict_ = {'f1': [], 'precision': [], 'recall': [], 'Accuracy':[]}
        pred_dict_ = {}

        print('1. load train test data')
        train_X, train_y, test_X, test_y = load_data(m, data_dir, 'without_types')

        tag = 'all_info'
        feature_less = dv_feature + dv_result_feature + alignment_feature

        train_X = train_X[feature_less]
        test_X = test_X[feature_less]

        print(train_X.columns)
        print(tag)

        result_dict = {'F1 Score': [], 'Precision': [], 'Recall': [], 'Accuracy':[], 'AUC': []}
        print('2. model learning')
        for i, model_name in enumerate(models):
            print('2.' + str(i + 1) + ' | model:', model_name)
            model, pred, pred_proba = eval(model_name)(train_X, test_X, train_y, test_y)
            temp_result = score_f(test_y, pred, pred_proba, result_dict)
            print(model_name, temp_result)
            if model_name == 'MLP':
                file_name = model_dir + m + '/' + tag + '_' + model_name + '_model.h5'
                model.save(file_name)
            elif model_name == 'ft_transformer':
                file_name = model_dir + m + '/' + tag + '_' + model_name + '_model.pth'
                torch.save(model.state_dict(), file_name)
            else:
                file_name = model_dir + m + '/' + tag + '_' + model_name + '_model.pkl'
                pickle.dump(model, open(file_name, 'wb'))

        label_df = pd.DataFrame.from_dict(result_dict, orient='index', columns=models)
        label_df = label_df.transpose()
        print()
        print('task :', m)
        print('#train :', len(train_y), '#test :', len(test_y), '\n')
        print(tabulate(label_df, headers='keys', tablefmt='psql'))
        label_df.to_csv(performance_dir + m + '_' + tag + '_result.csv')
