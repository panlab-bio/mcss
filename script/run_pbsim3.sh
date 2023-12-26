path_out=$1
method=$2
model=$3
depth=$4
genome=$5
prefix=$6
pass_num=$7
sample_n=$8

sim_res=sim_res_${sample_n}
sim_concat=sim_concat_${sample_n}
# sim_strain=sim_strain${sample_n}
if [ ! -d $path_out/$sim_res ];then
    echo "mkdir $path_out/$sim_res"
    mkdir -p $path_out/$sim_res
fi

# if [ ! -d $path_out/$sim_strain ];then
#     echo "mkdir $path_out/$sim_strain"
#     mkdir -p $path_out/$sim_strain
# fi

if [ ! -d $path_out/$sim_concat ];then
    echo "mkdir $path_out/$sim_concat"
    mkdir -p $path_out/$sim_concat
fi

if [ ! -f $path_out/$sim_concat/sim.fastq ];then
    echo "mkdir $path_out/$sim_concat/sim.fastq "
    > $path_out/$sim_concat/sim.fastq
fi

# echo $genome
genome_name=$(basename $genome)
# echo $genome_name
# genome_new=$path_out/$sim_strain/$genome_name
genome_new=$genome
# echo "gzip -dc ${genome}.gz > $genome_new"
# gzip -dc ${genome}.gz > $genome_new

cd $path_out/$sim_res

echo "pbsim --strategy wgs --method $method --$method $model --depth $depth --genome $genome_new --prefix $prefix"

pbsim --strategy wgs --method $method --$method $model --depth $depth --genome $genome_new --prefix $prefix --pass-num $pass_num

samtools view -bS ${prefix}_0001.sam -o ${prefix}_0001.subreads.bam
ccs ${prefix}_0001.subreads.bam ${prefix}_0001.ccs.fastq.gz
gunzip ${prefix}_0001.ccs.fastq.gz
cat ${prefix}_0001.ccs.fastq >> $path_out/$sim_concat/sim.fastq

