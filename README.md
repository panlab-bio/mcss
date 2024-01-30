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
# The tool is currently developed based on amd64 POSIX and has undergone testing only on the Linux operating system.
# We are dedicated to resolving any issues encountered by users on different operating systems and have plans to introduce support for additional systems in future updates.
conda env create -f env_mcss.yaml
or:
conda create -n env_mcss python=3.8.13 -c conda-forge
source activate env_mcss
conda install -c bioconda pbccs=6.4.0
conda install biopython=1.78 -c conda-forge
conda install -c bioconda samtools=0.1
conda install -c bioconda kraken2=2.0.7_beta
conda install -c bioconda pbsim3=3.0.1
conda install pandas=1.2.5 -c conda-forge
conda install requests=2.31.0 -c conda-forge
conda install scipy=1.10.1 -c conda-forge
```
    
4. activate environment
```
conda activate env_mcss / source activate env_mcss
```
5. download strain genomes
```
tar xvzf data/strain_file/down_strain.tar.gz -C data/strain_file/
bash tools/strain_download.sh
```

## usage

1. simulate based on environment

```
./mcss.sh env -e gut -o sim_env 
   
       # -e choice environment 
   
       # -o output file 
   
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
# Recommended for execution exclusively on high-performance servers
./mcss.sh sample -i kraken_fastq/kraken_oral/ -o sim_sample_paired p 
   
       #  paired fastq reads
./mcss.sh sample -i kraken_fastq/kraken_rhi/ -o sim_sample_single
   
       #  single fastq reads
```

4. call pbsim3 ( simulate community and call pbsim3 can be done separately )
```
./mcss.sh pbsim -i sim_com_f/ --mean_depth 5
   
       #  call pbsim3 to generate reads of communities in sim_com_f
       #  sim_com_f/pbsim/sim_concat_1/sim.fastq
```

5. generate strains ( only strain genomes )
```
./mcss.sh strain -i sim_com_f -s 2
   
       #  strain genomes and abu
       #  sim_com_f/strains.csv

```

6. help doc
```
./mcss.sh -h
./mcss.sh env -h
./mcss.sh sample -h
./mcss.sh community -h
./mcss.sh pbsim -h

```

### citation
if you use MCSS, please cite:





