#!/usr/bin/python3


"""

"""


import pandas as pd
import os
import numpy as np
import getopt
import sys
from pandas import Series, DataFrame


insep='\t'

opts, args = getopt.getopt(sys.argv[1:],
    '-i:-o:-s:-S:-h',
    ['--infile', '--out_tbl', '--insep', '--outsep2', '--help'])
threshold_global = 0.95
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
  elif opt_name in ('-h','--help'):
    print(' '.join([
        'python3 blast-best_hit_calling2.py \\\n',
        '-i $input_file\\\n',
        '-o $out_file_tsv\\\n',
        '-s $in_sep\\\n',
        '-S $out_sep\\\n'
        ]))

label = 'original' # 'bn' or 'original'

def check_global_alignment(hit_rows): 
    #if len(pd.unique(df['#bn_query'])) != 1:
    #    raise ValueError("bn_query列的值不唯一")

    #if len(pd.unique(df['saccver'])) != 1 :
    #    raise ValueError("saccver列有多个不同值")
    max_score = hit_rows['scov'].max()
    max_rows = hit_rows[hit_rows['scov'] == max_score].iloc[[0]]
    return max_rows

def merge_regions(df):

  # 预处理：确保每个区间的起始<=终止
  df1 = df.copy()
  df1.loc[:, 'start'] = df[['sstart', 'send']].min(axis=1)
  df1.loc[:, 'end'] = df[['sstart', 'send']].max(axis=1) 
  #df1['start'] = df[['sstart', 'send']].min(axis=1) # 两种方法都可以
  #df1['end'] = df[['sstart', 'send']].max(axis=1) # 两种方法都可以
  
  # 按起始坐标排序
  df1 = df1.sort_values('start').reset_index(drop=True)
  #print(df1)

  # 初始化结果存储
  merged_intervals = [] 

  # 处理第一个区间
  if len(df1) > 0:
    current_start = df1.at[0, 'start']
    current_end = df1.at[0, 'end']
    
    # 遍历后续区间
    for i in range(1, len(df1)):
      next_start = df1.at[i, 'start']
      next_end = df1.at[i, 'end']
      
      # 检查重叠
      if next_start <= current_end:
        # 有重叠，合并区间\
        current_end = max(current_end, next_end)
      else:
        # 无重叠，保存当前合并后的区间
        merged_intervals.append([current_start, current_end])
        # 开始新的区间
        current_start = next_start
        current_end = next_end
        
    # 添加最后一个合并后的区间
    merged_intervals.append([current_start, current_end])
    
  # 转换为DataFrame
  output_df = pd.DataFrame(merged_intervals, columns=['sstart', 'send'])
  
  #print("合并后的不重叠区间：")
  #print(output_df)
  
  # 计算总长度
  total_length = (output_df['send'] - output_df['sstart'] + 1).sum()  
  #print(f"\n合并后区间总长度：{total_length}")
  return output_df

def cal_coverage(data, gene_size):
    # 计算每个区间的长度并求和  
    data['length'] = data['send'] - data['sstart'] + 1
    total_covered = data['length'].sum()
    
    # 计算覆盖度
    coverage = total_covered / gene_size
    #print('coverage', coverage)

    return coverage, total_covered

def test_merge_region():
    data = { 'sstart': [4953, 4971, 5391, 889, 3656, 3929, 3892, 759, 3620, 3778, 4965, 5388, 3926, 870, 759, 3889, 3617, 3775],
             'send': [1, 1, 1, 796, 3571, 3844, 3817, 667, 3571, 3730, 1, 1, 3841, 796, 667, 3837, 3568, 3727] }  
    #print(pd.DataFrame(data))
    merge_regions(pd.DataFrame(data))
    
    data = { 'sstart':[889, 3656, 759],
        'send': [ 796,3571,667]}
    data = merge_regions(pd.DataFrame(data))
    scov = cal_coverage(data, 5391)
    #print('参考基因覆盖度', scov)

    interval_strings = data.apply(lambda row: f"{row['sstart']}-{row['send']}", axis=1)
    result = ";".join(interval_strings)
    #print('result', result)

#test_merge_region()
#exit()


def main_bn(df):
  
  # main run base name mode
  cols = list(df.columns)
  #cols.append('PI')
  #cols.append('label')
  #print(cols)
  
  df = df.copy()
  #print(df)
  df = df[df['length'] >= 100]
  df['scov'] = df['length'] / df['slen']
  df = df[df['scov'] >= 0.05 ]
  #print(df)
  
  ncols = [ 'filename', 'queries', 'saccver', 'stitle', 'pident', 'scov', 's_cov_len', 'slen', 'label' ]
  results_df = pd.DataFrame(columns = ncols)
  
  results = []
  bn_lst = pd.unique(df['#bn_query'])
  
  count_bn = len(bn_lst)
  
  i = 0
  for bn in bn_lst:
    i += 1
    #print('  [main_bn] '+ str(i) + ' /' + str(count_bn) )
    bn_que_hits = df[ df['#bn_query'] == bn ] # all rows of the same basename

    #que_lst = pd.unique(bn_que_hits['qaccver']).tolist() 
    sbj_lst = pd.unique(bn_que_hits['saccver']).tolist()
    
    #print('sbj_lst', sbj_lst)
    
    for sbj in sbj_lst:
      #print('subject', sbj)
      hit_rows = bn_que_hits[ bn_que_hits['saccver'] == sbj ]
      
      # 1） 先判断有无全局匹配，有则跳到下一行
      max_row = check_global_alignment(hit_rows)

      #print(max_row['scov'].item())
      #print('threshold_global', threshold_global)
      #print('max_row', max_row)
      #print("max_row['scov'].iloc[0]", max_row['scov'].iloc[0])
      #print("max_row['scov'].item()", max_row['scov'].item())
      #exit()

      if max_row['scov'].item() >= threshold_global:
        olst = [ max_row['#bn_query'].item(), max_row['qaccver'].item(), max_row['saccver'].item(), max_row['stitle'].item(), 
                max_row['pident'].item(), max_row['scov'].item(), max_row['scov'].item(), max_row['slen'].item(), 'global' ]
        orow = pd.DataFrame([olst], columns = ncols )
        results_df = pd.concat([results_df, orow], axis = 0)
        #continue 
      else: # 2) 
        orow = pd.DataFrame(columns = ncols)
        merge_regs = merge_regions( hit_rows[['sstart','send']] )
        #print('合并区间merge_regs', merge_regs, '\n')
        
        # 用分号连接所有区间字符串
        interval_strings = merge_regs.apply(lambda row: f"{row['sstart']}-{row['send']}", axis=1)
        result = ";".join(interval_strings)
        #print('result', result,'\n')

        sbj_gene_size = max_row['slen'].item()
        #print('sbj_gene_size',sbj_gene_size)

        scov,s_cov_len = cal_coverage(merge_regs, sbj_gene_size)
        #print('scov', scov, '\n')
        #print('s_cov_len', s_cov_len, '\n')
        

        queries = hit_rows['qaccver'].drop_duplicates().str.cat(sep=';')
        #print(hit_rows[['pident','scov']])
        #print('query序列', queries, '\n')
        
        # 计算加权平均值（pident 是数据，scov 是权重）
        weighted_avg = np.average(hit_rows['pident'], weights=hit_rows['scov'])
        #print(f"加权相似度平均值: {weighted_avg:.4f}",'\n')

        refgene_desp = max_row['stitle'].item()

        # ncols = [ 'filename', 'queries', 'saccver', 'stitle', 'pident', 'scov', 's_cov_len', 'slen', 'label' ]  
        olst = [ bn, queries, sbj, refgene_desp,
                 weighted_avg, scov, s_cov_len, sbj_gene_size, 
                 'global' ]
        #print('olst', olst)

        orow = pd.DataFrame( [olst], columns = ncols )  
        #print('orow', orow)

        results_df = pd.concat([results_df, orow], axis = 0)
        #print('results_df', results_df)
        #results_df.to_csv('test_out.csv', index=False, encoding='utf-8')

      #print(hit_rows)
  #print(results_df)
 
  return results_df

df = pd.read_csv(inf, sep = insep)
results_df = main_bn(df)
results_df.to_csv(ouf, sep=outsep, index=False)

