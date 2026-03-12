#!/usr/bin/python3

"""
将VFDB的setX.fas的header提取处理，梳理为一个表格
"""

# setB

import re
import pandas as pd
import sys

args = sys.argv[1:]


if len(args) < 2:
    print("用法: python3 process.py <input> <output>")
    sys.exit(1)

input_file = args[0]
output_file = args[1]


with open(input_file, "r") as f:
    data=f.readlines()
anno_list=[]
for line in data:
    info = re.findall('>(\S+)\(gb\|\S+\) \((.*?)\) (.*) \[(.*) \((.*)\) - (.*) \((.*)\)\] \[(.*)\]',line)
    if len(info)==0:
        info = re.findall('>(\S+) \((.*?)\) (.*) \[(.*) \((.*)\) - (.*) \((.*)\)\] \[(.*)\]',line)
    tmp = pd.DataFrame(info)
    tmp.columns=['VF_gene_id','VF_gene_symbol','Gene_description','VF_name',
                 'VF_id','VF_category_level1','VF_category_id','taxonomy']
    anno_list.append(tmp)
anno = pd.concat(anno_list)
f.close()
anno.to_csv(output_file, sep='\t',index=0)

