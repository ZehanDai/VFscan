#!/bin/bash

set -e 

# 参数说明函数
usage() {
    echo "Usage: $0 -i <input_path> -r <reference_db_path> -o <output_path>"
    echo "Options:"
    echo "  -i  Input directory path (required)"
    echo "  -r  Reference database directory path (required)"
    echo "  -o  Output directory path (required)"
    exit 1
}

# 初始化变量
INPUT_PATH=""
REFERENCE_PATH=""
OUTPUT_PATH=""

# 解析命令行参数
while getopts ":i:r:o:" opt; do
    case $opt in
        i)
            INPUT_PATH="$OPTARG"
            ;;
        r)
            REFERENCE_PATH="$OPTARG"
            ;;
        o)
            OUTPUT_PATH="$OPTARG"
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            ;;
    esac
done

# 检查必需参数
if [ -z "$INPUT_PATH" ] || [ -z "$REFERENCE_PATH" ] || [ -z "$OUTPUT_PATH" ]; then
    echo "Error: Missing required arguments." >&2
    usage
fi


# 打印参数校验结果（调试用）
echo "=== Parameters ==="
echo "Input Path:      $INPUT_PATH"
echo "Reference Path:  $REFERENCE_PATH"
echo "Output Path:     $OUTPUT_PATH"




 

# 建数据库index
echo "Checking VFDB reference database index for blastn"
[ ! -f $REFERENCE_PATH/VFDB_setB_nt.ndb ] && makeblastdb -in $REFERENCE_PATH/VFDB_setB_nt.fas -out $REFERENCE_PATH/VFDB_setB_nt -dbtype nucl 

mkdir -p $OUTPUT_PATH

# 判断输入文件后缀-输出列表存到一个变量里
fasta_extensions=("fasta" "fa" "fna" "fas")  # 支持的FASTA后缀
found_extensions=()  # 存储实际找到的后缀

# [测试用] echo real files `ls $INPUT_PATH/*.fna`
# 检查每种后缀是否存在文件
# [测试用] echo fasta file extension
for ext in "${fasta_extensions[@]}"; do
    if ls "$INPUT_PATH"/*."$ext" 1> /dev/null 2>&1; then
        found_extensions+=("$ext")
    fi
done
# [测试用] echo foudn extensions: $found_extensions


# 判断结果
if [ ${#found_extensions[@]} -eq 0 ]; then
    echo "Error: No FASTA files found in $INPUT_PATH with supported extensions (.fasta, .fa, .fna, .fas)"
    exit 1
elif [ ${#found_extensions[@]} -gt 1 ]; then
    echo "Error: Multiple FASTA file extensions found in $INPUT_PATH: ${found_extensions[*]}"
    echo "Please use only one type of FASTA file extension in the input directory."
    exit 1
else
    # 只有一种后缀，获取文件列表
    ext="${found_extensions[0]}"
    input_fasta=($(ls "$INPUT_PATH"/*."$ext"))
    #echo "Found ${#input_fasta[@]} .$ext files:"
    #printf '%s\n' "${input_fasta[@]}"
    
    # 文件列表已存储在 $input_fasta 数组中
    # 后续可以通过 ${input_fasta[@]} 访问所有文件
fi

threads=$(( $(nproc) - 2 ))
count=0
echo Homolog scanning for virulence factor gene using blastn
for fa in `echo ${input_fasta[@]}`; do
  count=$(( $count + 1 ))
  #echo fa $fa
  
  bn=`basename $fa .fasta`
  # echo $bn
  bn=${bn%%.*}
  echo "  [$count/${#input_fasta[@]}]" aligning file $bn
   
  # blast比对VFDB
  [ ! -f $OUTPUT_PATH/$bn.fmt6 ] && blastn -db $REFERENCE_PATH/VFDB_setB_nt \
    -query $fa \
    -out $OUTPUT_PATH/$bn.fmt6 \
    -word_size 21 \
    -outfmt '6 qaccver saccver pident length qstart qend sstart send qlen slen stitle' \
    -num_threads $threads
  
  # 给fmt6增加header
  hder=`echo '#bn_query qaccver saccver pident length qstart qend sstart send qlen slen stitle' | sed "s/ /\t/g"`
  [ ! -f $OUTPUT_PATH/$bn.add_hder.fmt6 ] && sed "s/^/$bn\t/g" $OUTPUT_PATH/$bn.fmt6 | sed "1i $hder" > $OUTPUT_PATH/$bn.add_hder.fmt6 # 给format 6格式输出增加header 

  # 最优判断
  echo "  "Running besthit calling...
  [ ! -f $OUTPUT_PATH/$bn.tsv ] && python3 utils/best_hit_calling.py -i $OUTPUT_PATH/$bn.add_hder.fmt6 -o $OUTPUT_PATH/$bn.tsv \
    -s "\t" -S "\t"


  # 增加VFDB信息
  echo "  "Add VFDB information...
  [ ! -f $OUTPUT_PATH/$bn.addVFDBinfo.tsv ] && python3 utils/add_VFDB_info.py -i $OUTPUT_PATH/$bn.tsv \
      -o $OUTPUT_PATH/$bn.addVFDBinfo.tsv  \
      -a $REFERENCE_PATH/SetB_info.tsv \
      -s "\t" -S "\t"
  
done

python3 utils/concat_VFfiles.py -i $OUTPUT_PATH -o $OUTPUT_PATH -s "\t"


 

  
