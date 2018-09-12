# ChIP-Seq Workflow Tutorial

### 1. Make project directory
```
# the project directory contains specific GROUP name and current DATE
GROUP=XXX
DATE=`date +"%Y%m%d"`
mkdir ~/Project/${GROUP}_${DATE}
```

### 2. Clone the repository
```
cd ~/Project/${GROUP}_${DATE}
git clone https://github.com/bioxfu/ChIP-Seq
cd ChIP-Seq
```

### 3. Copy/Download the raw data
```
DATAPATH=/the/path/of/the/raw/data/on/HPC

# if you are working on the HPC, copy the raw data
./script/copy_rawdata.sh $DATAPATH

# if you are working on the local machine, download the raw data
./script/copy_rawdata.sh $DATAPATH --download

# check the sample (folder) names
# the pattern of the sample names should be:
# IP_NNN.1, IP_NNN.2, IP_NNN.3
# Input_NNN.1, Input_NNN.2, Input_NNN.3
# NNN should NOT start with a number
```

### 4. Rename the raw data
```
# if one sample has one run, rename the run
# dry run to check if mv command is correct
./script/rename_rawdata.sh --dry-run
# then do it 
./script/rename_rawdata.sh

# if one sample has multiple runs, combine the runs
./script/rename_and_combine_rawdata.sh
```

### 5. Create *config.yaml* and *Snakefile* based on the examples
```
cp example/example_config.yaml config.yaml
cp example/example_Snakefile Snakefile
cp example/example_sample_sheet.csv sample_sheet.csv

# edit config.yaml and sample_sheet.csv 
```

### 6. Initiate the project
```
source init.sh
```

### 7. Dry run the workflow to check any mistakes
```
./dry_run.sh
```

### 8. Run the workflow
```
# if you are working on the HPC
./run_HPC.sh

# if you are working on the local machine
./run.sh

# check the workflow progress in nohup.out file
tail nohup.log 

# check the jobs on HPC
qstat

# if you get the error: Directory cannot be locked.
snakemake --unlock 
```

### 9. Remove the temporary files
```
./clean.sh
```

