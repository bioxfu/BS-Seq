# BS-Seq workflow

### 1. Clone the repository
```
git clone https://github.com/bioxfu/BS-Seq
```

### 2. Initiate the project
```
source init.sh
```

### 3. Download SRA data
```
bash download_sra.sh
```

### 4. Create *config.yaml* and *Snakefile* based on the examples

### 5. Dry run the workflow to check any mistakes
```
./dry_run.sh
```

### 6. Start the workflow
```
./run.sh
```

### 7. Check the workflow progress in *nohup.out* 

### 8. Remove the temporary files
```
./clean.sh
```

### install HOMER tair10 annotation
perl /home/xfu/miniconda2/envs/gmatic/share/homer-4.9.1-5/configureHomer.pl -install tair10
