#读取srr.out,得到srr的read数
#bash 1_2_get_read_cnt.sh test batch
#把get_newout整合到这里面来了
#把get_read_abu也整合进来吧
out_file_root=$1
# echo $out_file_root
batch=$2
path_out=$out_file_root/kraken_res
#这个其实是输入文件，只是输入的.out的结果
path_data=$out_file_root/data/kraken

file_path="$path_data/read_cnt/read_cnt.txt"
#get_newout的文件
path_out_new="$path_data/out_new"

#abu
path_abu="$path_data/abundance/"

while read -r srr; do
    name=$path_out/${srr}.out
    rep=$path_out/${srr}_rep.txt
    if [ -s $rep ]; then
        # echo $name
        wc -l $name >>$file_path
        # echo "wc -l $name >>$file_path"
        awk -F '[\t ]' '{print $2"\t"$3}' $name>$path_out_new/${srr}.out
        # echo "awk -F '[\t ]' '{print $2"\t"$3}' $name>$out_file/${srr}.out"
        grep "s__" $rep |awk -F '[| \t]' '{print$(NF-2)"_"$(NF-1)"\t"$NF}'   >$path_abu/${srr}_abundance.txt 
        # echo "grep "s__" $rep |awk -F '[| \t]' '{print$(NF-2)"_"$(NF-1)"\t"$NF}'   >$path_abu/${srr}_abundance.txt"
    else
        if [ -s $name ]; then
            rm $name
        fi
        
    fi
    
done < $batch

