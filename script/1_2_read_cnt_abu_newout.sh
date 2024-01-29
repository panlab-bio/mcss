# Read the srr.out file to obtain the number of reads for sample.
#bash 1_2_get_read_cnt.sh test batch
out_file_root=$1
batch=$2
path_out=$out_file_root/kraken_res
path_data=$out_file_root/data/kraken
file_path="$path_data/read_cnt/read_cnt.txt"
path_out_new="$path_data/out_new"

#abu
path_abu="$path_data/abundance/"

while read -r srr; do
    name=$path_out/${srr}.out
    rep=$path_out/${srr}_rep.txt
    if [ -s $rep ]; then
        
        wc -l $name >>$file_path
        
        awk -F '[\t ]' '{print $2"\t"$3}' $name>$path_out_new/${srr}.out
        
        grep "s__" $rep |awk -F '[| \t]' '{print$(NF-2)"_"$(NF-1)"\t"$NF}'   >$path_abu/${srr}_abundance.txt 
        
    else
        if [ -s $name ]; then
            rm $name
        fi
        
    fi
    
done < $batch

