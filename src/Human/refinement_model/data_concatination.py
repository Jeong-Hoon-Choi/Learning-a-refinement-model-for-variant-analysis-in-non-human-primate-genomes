import pandas as pd
from common_module import *
from common_parameter import *


def return_1_0(df):
    num1 = len(df.loc[df['label'] == 1].reset_index(drop=True).index)
    num2 = len(df.loc[df['label'] == 0].reset_index(drop=True).index)
    return num1, num2

data_2 = pd.read_csv(data_dir + 'HG001_concat_match_all_.csv', index_col=0)
data_1 = pd.read_csv(data_dir + 'HG002_concat_match_all_.csv', index_col=0)

data_2_X, data_2_Y = train_test_split(data_2, test_size=0.5, random_state=1024, stratify=data_2[label])
data_1_X, data_1_Y = train_test_split(data_1, test_size=0.5, random_state=1024, stratify=data_1[label])
data_train = pd.concat([data_2_X, data_1_X]).reset_index(drop=True)
data_test = pd.concat([data_2_Y, data_1_Y]).reset_index(drop=True)

data_train.to_csv(data_dir + 'mixed_match_train.csv')
data_test.to_csv(data_dir + 'mixed_match_test.csv')

train_1, train_0 = return_1_0(data_train)
train_snv_1, train_snv_0 = return_1_0(data_train.loc[data_train['TYPE'].str.contains('.snp')].reset_index(drop=True))
train_indel_1, train_indel_0 = return_1_0(data_train.loc[data_train['TYPE'].str.contains('.indel')].reset_index(drop=True))
print(train_1 + train_0, train_1, train_0, train_0 / (train_1 + train_0))
print(train_snv_1 + train_snv_0, train_snv_1, train_snv_0, train_snv_0 / (train_snv_1 + train_snv_0))
print(train_indel_1 + train_indel_0, train_indel_1, train_indel_0, train_indel_0 / (train_indel_1 + train_indel_0))

print()

test_1, test_0 = return_1_0(data_test)
test_snv_1, test_snv_0 = return_1_0(data_test.loc[data_test['TYPE'].str.contains('.snp')].reset_index(drop=True))
test_indel_1, test_indel_0 = return_1_0(data_test.loc[data_test['TYPE'].str.contains('.indel')].reset_index(drop=True))
print(test_1 + test_0, test_1, test_0, test_0 / (test_1 + test_0))
print(test_snv_1 + test_snv_0, test_snv_1, test_snv_0, test_snv_0 / (test_snv_1 + test_snv_0))
print(test_indel_1 + test_indel_0, test_indel_1, test_indel_0, test_indel_0 / (test_indel_1 + test_indel_0))
