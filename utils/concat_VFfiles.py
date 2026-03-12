#!/usr/bin/python3


"""

"""


import pandas as pd
import os
import numpy as np
import getopt
import sys
import glob
from pandas import Series, DataFrame

insep='\t'

opts, args = getopt.getopt(sys.argv[1:], 
                           '-i:-o:-s:-h', 
                           ['--input_dir', '--out_dir', '--sep', '--help'])

def interpret_escape_sequences(value):
  """Convert escape sequences in a string to actual characters."""
  return value.encode('utf-8').decode('unicode_escape')
  
for opt_name,opt_value in opts:
  if opt_name in ('-i'):
    ind=opt_value
    #print('input path =',ind)
  elif opt_name in ('-o'):
    oud=opt_value
    #print('output path =',oud)
  elif opt_name in ('-s'):
    outsep=opt_value
    outsep=interpret_escape_sequences(outsep)
    #print("output delitmmer",outsep)
  elif opt_name in ('-h','--help'):
    print(' '.join([
        'python3 concat_VFfiles.py \\\n',
        '-i input_path\\\n',
        '-o output_path\\\n',
        '-s $out_sep\\\n'
        ]))


# 1. 合并所有文件
all_files = glob.glob(ind + "/*.addVFDBinfo.tsv")
#print(all_files)


combined_df = pd.concat([pd.read_csv(f, sep=outsep) for f in all_files], ignore_index=True)
#print(combined_df)
ouf=oud+'/merge_VFDB_long.tsv'
combined_df.to_csv(ouf, sep="\t",)

# 2. 生成矩阵（自动处理重复值）
matrix = combined_df.pivot_table(
        index="filename",
        columns="VF_gene_symbol",
        values="fix_iden",
        aggfunc="first"  # 保留第一个出现的值（无重复时等效于 pivot）
        )

 
# 3. 填充缺失值并保存
# matrix.fillna(0).to_csv("gene_matrix.csv", sep = '\t')
#matrix.fillna(0).to_csv(otbl, sep = '\t')
#[]matrix.fillna(0).to_csv(otbl, sep="\t", quoting=csv.QUOTE_NONE)
ouf=oud+"/merge_VFDB_fix_iden_matrix.tsv"
matrix.fillna(0).to_csv(ouf, sep="\t", float_format="%.6g")

print('Merge all VFDB output to ', ouf )

