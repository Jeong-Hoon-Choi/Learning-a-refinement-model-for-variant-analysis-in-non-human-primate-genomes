import pysam
import pandas as pd
import io
from tqdm import tqdm
from collections import Counter
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import os

# indian origin rhesus macaque individual id
ind_list_30 = [35044, 35045, 35046, 35048, 35049, 35051, 35055, 35059, 35060, 35061]

# chinese origin rhesus macaque individual id
ind_list_china = [35082, 35086, 36013, 36390, 36394, 37854, 37945, 37950]
ind_list_china_high_x = [35082, 35086, 36013, 36390, 36394, 37854, 37945, 37950]

# all rhesus macaque individual
ind_list_all = [ind_list_china, ind_list_30]

# chromosome index to name dictionary
rhe_mac_chr_dict_10_reverse = {'1': 'NC_041754.1', '2': 'NC_041755.1', '3': 'NC_041756.1', '4': 'NC_041757.1',
                               '5': 'NC_041758.1', '6': 'NC_041759.1', '7': 'NC_041760.1', '8': 'NC_041761.1',
                               '9': 'NC_041762.1', '10': 'NC_041763.1', '11': 'NC_041764.1', '12': 'NC_041765.1',
                               '13': 'NC_041766.1', '14': 'NC_041767.1', '15': 'NC_041768.1', '16': 'NC_041769.1',
                               '17': 'NC_041770.1', '18': 'NC_041771.1', '19': 'NC_041772.1', '20': 'NC_041773.1',
                               'X': 'NC_041774.1', 'Y': 'NC_027914.1', 'M': 'NC_005943.1'}

# rhesus macaque reference genome sequence file
ref_file_path = '../../../data/Rhesus_Macaque/sequence/reference/GCF_003339765.1_Mmul_10_genomic.fna'


path_refinement_result = '../../../result/Rhesus_Macaque/csv/'
path_preprocess_csv = '../../../data/Rhesus_Macaque/csv/'
path_preprocess_sequence = '../../../data/Rhesus_Macaque/sequence/'
path_alt_base_ratio_result = '../../../result/Rhesus_Macaque/alt_base_ratio/'
