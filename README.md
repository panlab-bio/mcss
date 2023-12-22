# mcss
MCSS: microbial community simulator based on structure
## installation

1. clone repository
```
git clone https://github.com/panlab-bio/mcss
```

2.  enter mcss path
```
cd mcss
```
3. create python environment
```
conda env create -f env_mcss.yaml 
conda create -n env_mcss python=3.8.13
source activate env_mcss
conda install -c bioconda pbccs
conda install -c bioconda biopython=1.79
conda install -c bioconda samtools
conda install -c bioconda kraken2
conda install -c bioconda pbsim3
conda install pandas
conda install requests
conda install scipy




```
    
4. activate environment
```
conda activate env_mcss / source activate env_mcss
```

## usage

1. simulate based on environment

```
./mcss.sh env -e gut -o sim_env pbsim 
   
       # -e choice environment 
   
       # -o output file 
   
       # call pbsim to generate fastq reads
```
   
   
2. simulate based on GTDB species (User-specified )
```
./mcss.sh community -f example/example_sp/sp_ex_1.txt -a example/example_abu/abu_ex_1.txt -o sim_com_f 
   
       #  file with GTDB species
       #  User-specified abundance
./mcss.sh community -d example/example_sp/ -a example/example_abu/ -o sim_com_d
   
       #  dir with fastq files
       #  User-specified abundance
       
./mcss.sh community -d example/example_sp/ -a gut -o sim_com_g
   
       #  dir with fastq files
       #  gut environment abundance
```

3. simulate based on fastq reads
```
./mcss.sh sample -i kraken_fastq/kraken_oral/ -o sim_sample_paired p 
   
       #  paired fastq reads
./mcss.sh sample -i kraken_fastq/kraken_rhi/ -o sim_sample_single
   
       #  single fastq reads
```

4. call pbsim3 ( simulate community and call pbsim3 can be done separately )
```
./mcss.sh pbsim -i sim_com_g/
   
       #  call pbsim3 to generate reads of communities in sim_com_g
       
   
```
5. help doc
```
./mcss.sh -h
./mcss.sh env -h
./mcss.sh sample -h
./mcss.sh community -h
./mcss.sh pbsim -h

```

### citation
if you use MCSS, please cite:





