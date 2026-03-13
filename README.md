# VFscan
VFscan is a lightweight pipeline for identifying homologs of virulence factor (VF) genes in genomic assemblies. It integrates BLAST searches against the VFDB database and generates a summarized annotation table.

## 📌 Version
* Tool version: 0.2.1
* built-in VFDB database version: 2024-12-27

## ✨ Features
* Accepts FASTA files (assemblies from WGS or metagenomic data)
* Automatically runs BLAST against the VFDB core dataset (SetB)
* Merges BLAST results with VFDB metadata to produce a comprehensive annotation table
* Supports custom databases with the same directory structure

## ⚙️ Installation
### 1. Clone the repository
```
bash

git clone https://github.com/yourusername/VFscan.git
cd VFscan
```
unzip all zip files in the database path

### 2. Set the VFSCAN environment variable
Add the following line to your ~/.bashrc (or equivalent shell config file):
```
bash

export VFSCAN=/path/to/VFscan   # replace with actual path
```

Then apply the changes:
```
bash

source ~/.bashrc
```
*Optional*: If you prefer not to modify your shell config, you can set the variable manually before each run:
```
bash

VFSCAN=/path/to/VFscan; export VFSCAN
```

3. (Optional) Update the VFDB database
The repository includes the VFDB database at the time of release. To use an updated version, replace the files inside the database/ directory with the new ones, preserving the same directory structure. See the database/ tree below:
```
> tree database/

database/
├── Comparative_tables_from_VFDB
├── SetA_anno.txt
├── SetB_anno.txt
├── SetB_info.tsv
├── VFDB_setA_nt
│   └── VFDB_setA_nt.fas
├── VFDB_setA_pro
│   └── VFDB_setA_pro.fas
├── VFDB_setB_nt
│   ├── VFDB_setB_nt.fas
├── VFDB_setB_pro
│   └── VFDB_setB_pro.fas
├── VFs.xls
└── current_version.TXT
```

## 🧪 Dependencies
Operating system: Linux (tested on Ubuntu 22.04 )
* BLAST+ (teseted in 2.15.0+)
* Python3 3.8+ (tested in 3.10.6)
* Python packages:
  * numpy  (tested in 1.23.4)
  * pandas  (tested in 1.5.1)

## 📥 Input / 📤 Output
### Input
One or more FASTA files containing genomic assemblies (contigs/scaffolds).
Example:
 text
/path/to/input/Rick22/assembly.fasta

### Output
A tab-separated summary table (.tsv) with the following columns:

| Column | Description |
|--------|-------------|
| `filename` | Name of the input file |
| `queries` | Query sequence identifier |
| `saccver` | Subject accession (VFDB gene identifier) |
| `stitle` | Subject title (full VFDB annotation) |
| `pident` | Percentage of identical matches |
| `scov` | Subject coverage (aligned length / subject length) |
| `s_cov_len` | Subject coverage length (aligned length) |
| `slen` | Subject sequence length |
| `label` | Global or specific hit classification (currently always "global") |
| `VF_gene_id` | VFDB gene identifier |
| `VF_gene_symbol` | Gene symbol (e.g., HP0256) |
| `Gene_description` | Functional description of the gene |
| `VF_name` | Virulence factor name (e.g., Flagella) |
| `taxonomy` | Source organism of the VFDB entry |
| `fix_iden` | Fixed identity (average pident weighted by scov) |
| `fix+iden` | fix_iden converted to positive/negative |
 
A sample row:
```
text

1200045C5  1200045C5_00001  VFG043369  VFG043369 (HP0256) ...  95.571  1.0  1.0  429  global  VFG043369  HP0256  involved in motility ...  Flagella  Helicobacter pylori 26695  95.571  1.0
```

## usage 
### 1. Run the helper script to verify your setup and reference paths:
`bash envBD_v0.2.sh`

### 2. Basic run
```
bash

INDIR=/path/to/input/fasta_files
OUTDIR=/path/to/output

bash batch_blastn_setB.sh -i "$INDIR" -o "$OUTDIR"
```
This will use the default VFDB database (SetB) included in the package.

### 3. Run with custom database 
```
bash

REFD=/path/to/custom/VFDB_setB_nt
bash batch_blastn_setB.sh -i "$INDIR" -o "$OUTDIR" -r "$REFD"
```
Arguments for batch_blastn_setB.sh:
* -i : Input directory containing FASTA files
* -o : Output directory (will be created if not exists)
* -r : (Optional) Path to VFDB_setB_nt directory; defaults to $VFSCAN/database/VFDB_setB_nt

