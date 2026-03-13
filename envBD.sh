#!/bin/bash
# script for setting up VFDB environment
# usage: bash ./envBD_v0.2.sh or ./envBD.sh -d 新VFSCAN路径

# version 2026-03-13 15:56


# color code for checking and alarming 
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

# 1 check whether the module path has been set or not
# 检查全局变量是否设置了模块路径
if [[ -n "$VFSCAN" ]]; then
    echo -e "${GREEN}[SUCCESS]${NC} VFSCAN has been set: $VFSCAN"
else
    echo -e "${RED}[ERROR]${NC} Modude path not found in global environment"
    echo "Example for setting: VFSCAN=/mnt/d/database/VFDB/VFscan_module/VFscan_v0.2; export VFSCAN"
    echo ""
    exit 1 
fi





# 2 检查数据库路径及文件
# 默认VFDB文件在模块路径下的database
# 如果有传参，则设置传入的数据库路径
DB_PATH="${VFSCAN}/database" # default database path 
# intepret parameters
while getopts ":d:" opt; do
    case $opt in
        d)
            DB_PATH="$OPTARG"
            ;;
        \?)
            echo -e "${RED}[ERROR]${NC} invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo -e "${RED}[ERROR]${NC} option -$OPTARG require assignment" >&2
            exit 1
            ;;
    esac
done


wd=`pwd`; cd $DB_PATH
[ ! -r $DB_PATH/Comparative_tables_from_VFDB ] && unzip Comparative_tables_from_VFDB.zip
[ ! -f $DB_PATH/VFDB_setA_nt/VFDB_setA_nt.fas ] && unzip VFDB_setA_nt.zip
[ ! -f $DB_PATH/VFDB_setA_pro/VFDB_setA_pro.fas ] && unzip VFDB_setA_pro.zip
[ ! -f $DB_PATH/VFDB_setB_nt/VFDB_setB_nt.fas ] && unzip VFDB_setB_nt.zip
[ ! -f $DB_PATH/VFDB_setB_pro/VFDB_setB_pro.fas ] && unzip VFDB_setB_pro.zip                             
cd $wd


# 3 fasta header信息转表格
[ ! -f $DB_PATH/SetA_anno.txt ] && zcat $DB_PATH/VFDB_setA_nt.fas.gz | grep '^>' > $DB_PATH/SetA_anno.txt
[ ! -f $DB_PATH/SetB_anno.txt ] && zcat $DB_PATH/VFDB_setB_nt.fas.gz | grep '^>' > $DB_PATH/SetB_anno.txt
[ ! -f $DB_PATH/SetB_info.tsv ] && python3 $VFSCAN/utils/process.py $DB_PATH/SetB_anno.txt $DB_PATH/SetB_info.tsv
[ ! -f $DB_PATH/SetA_info.tsv ] && python3 $VFSCAN/utils/process.py $DB_PATH/SetA_anno.txt $DB_PATH/SetA_info.tsv
[ ! -f $DB_PATH/VFDB_setA_nt/SetA_info.tsv ] && ln -s $DB_PATH/SetA_info.tsv $DB_PATH/VFDB_setA_nt/SetA_info.tsv
[ ! -f $DB_PATH/VFDB_setB_nt/SetB_info.tsv ] && ln -s $DB_PATH/SetB_info.tsv $DB_PATH/VFDB_setB_nt/SetB_info.tsv
 

# 4 检查 blastn
echo -e "\n=== 检查 BLAST+ 环境 ==="
if command -v blastn &> /dev/null; then
    blast_version=$(blastn -version | head -n 1)
    echo -e "${GREEN}[SUCCESS]${NC} BLAST+ 已安装: $blast_version"
else
    echo -e "${RED}[FAIL]${NC} 未找到 blastn 命令"
 fi

# 5 构建blast数据库索引

[ ! -f $DB_PATH/VFDB_setB_nt/VFDB_setB_nt.nin ] && makeblastdb -in $DB_PATH/VFDB_setB_nt/VFDB_setB_nt.fas -out $DB_PATH/VFDB_setB_nt/VFDB_setB_nt -dbtype nucl 
[ ! -f $DB_PATH/VFDB_setA_nt/VFDB_setA_nt.nin ] && makeblastdb -in $DB_PATH/VFDB_setA_nt/VFDB_setA_nt.fas -out $DB_PATH/VFDB_setA_nt/VFDB_setA_nt -dbtype nucl


# 6 检查 python3
echo -e "\n=== 检查 Python3 环境 ==="
if command -v python3 &> /dev/null; then
    py_version=$(python3 --version)
    echo -e "${GREEN}[SUCCESS]${NC} Python3 已安装: $py_version"
    
    # 检查 Python 依赖
    echo -e "\n=== 检查 Python 依赖包 ==="
    dependencies=("numpy" "pandas")
    
    for pkg in "${dependencies[@]}"; do
        if python3 -c "import $pkg" &> /dev/null; then
            pkg_version=$(python3 -c "import $pkg; print($pkg.__version__)")
            echo -e "${GREEN}[SUCCESS]${NC} ${pkg} 已安装 (版本: $pkg_version)"
        else
            echo -e "${RED}[FAIL]${NC} 未找到 Python 包: ${YELLOW}$pkg${NC}"
            echo -e "安装命令: ${YELLOW}python3 -m pip install $pkg${NC}"
        fi
    done
else
    echo -e "${RED}[FAIL]${NC} 未找到 python3 命令"
    echo -e "请先安装 Python3:"
    echo -e "1. Linux: ${YELLOW}sudo apt-get install python3${NC}"
    echo -e "2. Mac: ${YELLOW}brew install python${NC}"
    echo -e "3. Windows: 从官网下载安装 ${YELLOW}https://www.python.org/downloads/${NC}"
fi

echo -e "\n=== 环境检查完成 ===\n"
