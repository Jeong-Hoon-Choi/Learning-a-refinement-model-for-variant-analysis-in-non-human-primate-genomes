import pandas as pd

from common_parameter import *
from common_module import *

def return_1_0(df):
    num1 = len(df.loc[df['label'] == 1].reset_index(drop=True).index)
    num2 = len(df.loc[df['label'] == 0].reset_index(drop=True).index)
    return num1, num2

path_hg001 = data_dir + 'HG001' + '_concat_match_all_.csv'
hg001 = pd.read_csv(path_hg001)
hg001_snp = hg001[hg001['TYPE'].str.contains('.snp')].reset_index(drop=True)
hg001_indel = hg001[hg001['TYPE'].str.contains('.indel')].reset_index(drop=True)
tp, fp = return_1_0(hg001_snp)
print('HG001_snp', tp, fp, tp+fp, fp/(tp+fp))
tp, fp = return_1_0(hg001_indel)
print('HG001_indel', tp, fp, tp+fp, fp/(tp+fp))
tp, fp = return_1_0(hg001)
print('HG001', tp, fp, tp+fp, fp/(tp+fp))
print()

path_hg002 = data_dir + 'HG002' + '_concat_match_all_.csv'
hg002 = pd.read_csv(path_hg002)
hg002_snp = hg002[hg002['TYPE'].str.contains('.snp')].reset_index(drop=True)
hg002_indel = hg002[hg002['TYPE'].str.contains('.indel')].reset_index(drop=True)
tp, fp = return_1_0(hg002_snp)
print('HG002_snp', tp, fp, tp+fp, fp/(tp+fp))
tp, fp = return_1_0(hg002_indel)
print('HG002_indel', tp, fp, tp+fp, fp/(tp+fp))
tp, fp = return_1_0(hg002)
print('HG002', tp, fp, tp+fp, fp/(tp+fp))
