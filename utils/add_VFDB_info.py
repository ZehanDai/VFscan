#!/usr/bin/python3


"""
update 
v 0.1 与VFDB数据库自带注释信息表格对比，增加到blast结果表格
v 0.2 0818 修正相似度不超100%；
    一个菌株匹配多个subject，选择最高分的保留
        计算方法 加权相似度 * subject基因全局覆盖率
    （注意基于现有最优判断脚本无法判断多拷贝）
"""

import pandas as pd
import os
import numpy as np
import getopt
import sys
from pandas import Series, DataFrame

insep='\t'

opts, args = getopt.getopt(sys.argv[1:],
                           '-i:-o:-a:-s:-S:-h',                        
                           ['--infile', '--out_tbl', '--anno_tbl', 
                            '--insep', '--outsep2', '--help'])

#print('opts',opts)
#print('args',args)

def interpret_escape_sequences(value):
    """Convert escape sequences in a string to actual characters."""
    return value.encode('utf-8').decode('unicode_escape')

for opt_name,opt_value in opts:
    if opt_name in ('-i'):
        inf=opt_value
        #print('input file =',inf)
    elif opt_name in ('-o'):
        ouf=opt_value
        #print('outfile table =',ouf)
    elif opt_name in ('-s'):
        #insep=opt_value
        insep=interpret_escape_sequences(opt_value)
        #print("input delitmmer",insep)
    elif opt_name in ('-S'):
        outsep=opt_value
        outsep=interpret_escape_sequences(outsep)
        #print("output delitmmer",outsep)
    elif opt_name in ('-S'):
        reftbl=opt_value
        # print("annotation file")
    elif opt_name in ('-h','--help'):
        print(' '.join([
            'python3 add_VFDB_info_v0.3.py \\\n',
            '-i $input_file\\\n',
            '-o $out_file_tsv\\\n',
            '-s $in_sep\\\n',
            '-S $out_sep\\\n'
            ]))


anno_file = "/mnt/d/database/VFDB/SetB_info.tsv"

anno = pd.read_csv(anno_file, sep='\t')
anno.columns = anno.columns.str.replace(r'\r|\^M', '', regex=True)
#anno['Reference'] = anno['Reference'].str.replace(r'\r|\^M', '', regex=True)
#print(anno)


def map_info(odf, anno):
  '''
  匹配blast结果和VF注释文件的结果，增加gene_symbol信息
  '''
  odf['saccver_id'] = odf['saccver'].str.split('(').str[0]  # 获取括号前的部分作为基因 ID
  #print('匹配前odf',odf)
  #print(odf.columns)
  
  merged_df = pd.merge(
    odf, 
    anno[['VF_gene_id', 'VF_gene_symbol', 'Gene_description', 'VF_name', 'taxonomy']], 
    left_on='saccver_id', 
    right_on='VF_gene_id', 
    how='left')  # 使用左连接以保留 odf 中的所有行
  merged_df.drop(columns=['saccver_id'], inplace=True)
  #print('匹配后odf',merged_df)
  #print(merged_df.columns)
  return merged_df

def pick_top1_symbol(odf):
  """
  按基因symbol选择最优结果
  """
  odf = odf.copy()
  genes = odf['VF_gene_symbol'].unique()
  #print(genes)
  odf['fix_iden'] = odf['pident'] * odf['scov']
  
  top1_rows = odf.loc[odf.groupby('VF_gene_symbol')['fix_iden'].idxmax()]
  #for sbj in genes:
  #  hit_rows = odf[odf['VF_gene_symbol'] == sbj ]
  #  print(hit_rows)
  top1_rows.loc[ top1_rows['fix_iden'] > 1, 'fix+iden' ] = 1
#    exit()
  #print(top1_rows)
  #print(top1_rows['VF_gene_symbol'])
  return top1_rows

odf = pd.read_csv(inf, sep=insep)
#print(odf)
odf2 = map_info(odf, anno) # 匹配VF内容内容
odf2 = pick_top1_symbol(odf2)
#print(odf2)
odf2.to_csv(ouf, sep = outsep, index=False)

