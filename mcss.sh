#一个脚本将所有的联系起来
#!/bin/bash
#total是要根据fastq个数来的
#确定输入的参数
env_array=("gut" "oral" "skin" "marine" "soil" "rhizosphere")

current_dir=$(dirname $(realpath $0))
# path_pbsim=$current_dir/tools/pbsim3/src/pbsim
# echo $path_pbsim
# echo "路径" $current_dir
# 默认值
sim_param=""
get_param=""
date=$(date +'%m_%d')

# 检查是否提供了足够的参数
if [ $# -lt 1 ]; then
    echo "bash $0 kraken or sim"
    exit 1
fi

# 通过命令行参数确定要执行的操作
command="$1"
shift # 移除第一个参数，剩下的参数将由 getopts 处理

case "$command" in
    "env")
        current_dir=$(dirname $(realpath $0))
        path_tax="$current_dir/data/gtdb_taxonomy.tsv"
        output=$current_dir/sim_${date}
        n_samples=1
        flag_in=false
        flag_sl=0
        # flag_im=0
        flag_acc=0
        flag_acc_v=false
        # 解析cmss sim的参数
        while [[ $# -gt 0 ]]; do
            case "$1" in
                "-e" | "--env")
                    env="$2"
                    flag_in=true
                    shift 2
                    ;;
                -n|--n_samples)
                    n_samples="$2"
                    shift 2
                    ;;
                # sl|st_len)
                #     flag_sl=1
                #     shift 
                #     ;;
                acc|accurate)
                    if [ "$flag_acc_v" = true ]; then
                        echo "these two options (acc, prol) are in conflict"
                        exit 1
                    fi
                    flag_acc=1
                    shift 
                    ;;
                prol|prolific)
                    if [ "$flag_acc_v" = true ]; then
                        echo "these two options (acc, prol) are in conflict"
                        exit 1
                    fi
                    flag_acc=0
                    shift 
                    ;;
                -o|--output)
                    
                    output="$2"
                    output=$(realpath $output)
                    shift 2
                    ;;
                "pbsim")
                    flag_genome=1
                    flag_pbsim=true
                    cnt_strain=1
                    depth=$cnt_strain
                    type_depth=0
                    flag_depth=false
                    method="qshmm"
                    model="QSHMM-RSII.model"
                    pass_num=10
                    # if which pbsim3 &>/dev/null; then
#                     if [ -f $path_pbsim ]; then
#                         # path_pbsim=$(which pbsim3)
#                         echo " pbsim3 in tools "
                        
#                     elif which pbsim3 &>/dev/null; then
#                         echo "pbsim3 in path"
#                         path_pbsim=$(which pbsim3)
#                     else
#                         echo "please install pbsim3"
#                         exit 1
                    
#                     fi
                    
                    shift
                    ;;
                -s|--strain)
                    if [ "$flag_pbsim" = false ]; then
                        echo "only meaningful with pbsim"
                        exit 1
                    fi
                    cnt_strain="$2"
                    depth=$cnt_strain
                    
                    shift 2
                    ;;
                g|genome)
                    if [ "$flag_pbsim" = false ]; then
                        echo "only meaningful with pbsim"
                        exit 1
                    fi
                    flag_genome=0
                    shift 2
                    ;;
                --min_depth)
                    if [ "$flag_pbsim" = false ]; then
                        echo "only meaningful with pbsim"
                        exit 1
                    fi
                    
                    if [ "$flag_depth" = true ]; then
                        echo "these two options (--min_depth, --mean_depth) are in conflict"
                        exit 1
                    fi
                    flag_depth=true
                    type_depth=0
                    depth_pb="$2"
                    shift 2
                    ;;
                --mean_depth)
                    if [ "$flag_pbsim" = false ]; then
                        echo "only meaningful with pbsim"
                        exit 1
                    fi
                    if [ "$flag_depth" = true ]; then
                        echo "these two options (--min_depth, --mean_depth) are in conflict"
                        exit 1
                    fi
                    flag_depth=true
                    type_depth=1
                    depth_pb="$2"
                    shift 2
                    ;;

                --method)
                    if [ "$flag_pbsim" = false ]; then
                        echo "only meaningful with pbsim"
                        exit 1
                    fi
                    method="$2"
                    shift 2
                    ;;
                --model)
                    if [ "$flag_pbsim" = false ]; then
                        echo "only meaningful with pbsim"
                        exit 1
                    fi
                    model="$2"
                    shift 2
                    ;;
                --pass-num)
                    if [ "$flag_pbsim" = false ]; then
                        echo "only meaningful with pbsim"
                        exit 1
                    fi
                    pass_num="$2"
                    shift 2
                    ;;
                    
                -h|--help)
                    echo "mcss env"
                    
                    echo "-e,--env: choice an environment. gut,oral,skin,marine,soil,rhizosphere"
                    echo "-o,--output: output dir " 
                    echo "-n,--n_samples: count of simulate samples. default:1 "
                    # echo "sl,st_len: exclude significant species count bias"
                    echo "acc,accurate: accurate mode to seek out the subtree"
                    echo "prol,prolific: prolific mode to seek out the subtree "
                    echo "Choose one between acc and prol, with prol as the default value in env mode"
                    
                    echo "pbsim：execute  pbsim3"
                    echo "options for pbsim:"
                    echo "-s,--strain：number of strain"
                    echo "--min_depth,--mean_depth：min or average depth of species, min "
                    echo "--method: qshmm, errhmm. default:qshmm"
                    echo "--model QSHMM-ONT.model, QSHMM-RSII.model, ERRHMM-ONT.model, ERRHMM-RSII.model, ERRHMM-SEQUEL.model. default:QSHMM-RSII.model"
                    echo "--pass-num: if number of passes(pass-num)>1,multi-pass sequencing is performed to get hifi reads. default:10"
                    # echo "im,intermediate: keep intermediate results"
                    
                    exit 0
                    ;;
                *)
                    echo "error: -$1"
                    echo "mcss env -h for usage"
                    exit 1
                    ;;
            esac
        done
        # echo "执行cmss sim的操作，参数 -n $sim_param"
        # 在这里添加您的cmss sim代码，可以使用$sim_param
        #env模式下只需要abu这个文件夹，sp_cut啥的都不需要
        if [ "$flag_in" = false ]; then
            echo "error : choice an env"
            exit 1
        fi
        if [ "$flag_depth" = true ]; then
                $depth=$depth_pb
        fi
        
        if [ ! -d $output/community/abu ]; then
            mkdir -p $output/community/abu
        fi
        if [ ! -d $output/species ]; then
            mkdir -p $output/species
        fi
        
        path_abu_new=$current_dir/data/env/$env/abu.pkl
        cp $path_abu_new $output/community/abu/abu_env.pkl
        path_script=$current_dir/script
        python $path_script/sim_tree_env.py $output $current_dir $n_samples $env $flag_sl $path_tax $flag_acc
        if [ "$flag_pbsim" = true ]; then
            abu_name=$(ls $output/community/abu/)
            echo $abu_name
            if [ "$abu_name" = "abu_env.pkl" ]; then
                sim_mode=0
            elif [ "$abu_name" = "abu_sample.pkl" ]; then
                sim_mode=1
            elif [ "$abu_name" = "abu_user.pkl" ]; then
                sim_mode=2
            fi
            echo "sim_mode "$sim_mode
            path_abu_new=$output/community/abu/$abu_name
            echo $path_abu_new
            model_new=$current_dir/data/pbsim_model/$model
            
            if [ ! -f $model_new ]; then
                echo "no model named $model"
                exit 1
            fi
            if [ ! -d $output/pbsim ]; then
                mkdir -p $output/pbsim
            fi
            #这里有两种思路，一组就是直接调用mcss pbsim 
            #一种就是python 
            #pyhon需要把pbsim的东西全都复刻
            #mcss需要判断输入了哪些
            python $path_script/get_strain.py $current_dir $output $cnt_strain $depth $type_depth $path_abu_new $sim_mode $method $model_new $pass_num $flag_genome
            #这个的第二个参数和pbsim单独选项里的参数不一样，因为那个直接输入的input就是sim文件，我们的输出才是
        fi

        ;;
############################################################################################# sample
    "sample")
        normal_exit=true
        flag_pbsim=false
        kill_background_commands() {
            for pid in "${bg_pids[@]}"; do
                kill "$pid"
            done
        }

        run_background_commands(){
            local_arr=("$@")
            bg_pids=()
            normal_exit=false
            for localcmd in "${local_arr[@]}"; do
                # echo "$localcmd"
                $localcmd &
                bg_pids+=($!)  
            done
            wait
            
            normal_exit=true
        }

        get_part(){
            local path_dir=$1
            local mode_flag=$2
            local batch=$3
            if [ "$mode_flag" = true ]; then
                # echo "2"
                local total=$(($(ls $path_dir |wc -l)/2))

            else
                # echo "1"
                local total=$(ls $path_dir |wc -l)
            fi

            if [ $batch -gt $total ]; then
                # echo "warning: the number of batch exceeds the total count"
                
                batch=$total
            fi
            # echo "batch" $batch
            local part=$((total/batch))
            
            if [ $((part*batch)) != $total ]; then
                part=$(((total+10)/$batch))
            fi
            
            local parts=()
            local sum=0
            for ((i=1; i<=$batch; i++)); do

                sum=$((sum+part))
                # echo "sum"$sum
                parts+=("$sum") 
            done
            parts_len=${#parts[@]}
            parts[parts_len-1]=$total
            echo "${parts[@]}"
        }
        #得到kraken的目录
        # if which kraken2 &>/dev/null; then
        #     path_kraken=$(which kraken2)
        # else
        #     echo "install kraken2 first"
        #     exit 1
        # fi
        #path_kraken是执行文件，dir_kraken是执行目录，高一级
        # dir_kraken=${path_kraken%/*}
        
        dir_kraken=$current_dir/tools
        # echo $dir_kraken
        # kraken数据库地址
        path_db=$dir_kraken/kraken_db
        # echo $path_db
        # kraken
        #用来判断是不是pair模式
        mode_flag=false
        #是否从断点重新开始
        flag_in=false
        current_dir=$(dirname $(realpath $0))
        output=$current_dir/sim_${date}
        path_script=$current_dir/script
        path_tax="$current_dir/data/gtdb_taxonomy.tsv"
        # flag_len=0
        batch=1
        same_env=0
        n_samples=1
        flag_acc=1
        flag_acc_v=false   
        abu_sp="0.001"
        suf_pre="_1"
        suf="fastq.gz"
        while [[ $# -gt 0 ]]; do
            case "$1" in
            #i必须是绝对路径,实现了将相对路径转换为绝对路径
                "-i" | "--input")
                    flag_in=true
                    input="$2"
                    input=$(realpath $input)
                    shift 2
                    ;;
                -n|--n_samples)
                    n_samples="$2"
                    shift 2
                    ;;
                -asp|--abu_sp)
                    abu_sp="$2"
                    shift 2
                    ;;
                -suf|--suffix)
                    suf_pre="$2"
                    shift 2
                    ;;
                -suf_pair)
                    suf="$2"
                    shift 2
                    ;;
                acc|accurate)
                    if [ "$flag_acc_v" = true ]; then
                        echo "these two options (acc, prol) are in conflict"
                        exit 1
                    fi
                    flag_acc=1
                    flag_acc_v=true
                    shift 
                    ;;
                prol|prolific)
                    if [ "$flag_acc_v" = true ]; then
                        echo "these two options (acc, prol) are in conflict"
                        exit 1
                    fi
                    flag_acc=0
                    flag_acc_v=true
                    shift 
                    ;;
                -o|--output)
                    
                    output="$2"
                    output=$(realpath $output)
                    shift 2
                    ;;
                "p"|"paired")
                    # echo "p"
                    mode_flag=true
                    shift 
                    ;;
                -b|--batch)
                    batch="$2"
                    shift 2
                    ;;
                    
                -db|--kraken_db)
                    path_db="$2"
                    shift 2
                    ;;
                # "nl"|"new_len")
                #     flag_len=1
                #     shift
                #     ;;
                "ml"|"multi_sample")
                    same_env=1
                    shift
                    ;;
                "pbsim")
                    flag_genome=1
                    flag_pbsim=true
                    cnt_strain=1
                    depth=$cnt_strain
                    type_depth=0
                    flag_depth=false
                    method="qshmm"
                    model="QSHMM-RSII.model"
                    pass_num=10
#                     if [ -f $path_pbsim ]; then
#                         # path_pbsim=$(which pbsim3)
#                         echo " pbsim3 in tools "
                        
#                     elif which pbsim3 &>/dev/null; then
#                         echo "pbsim3 in path"
#                         path_pbsim=$(which pbsim3)
#                     else
#                         echo "please install pbsim3"
#                         exit 1
                    
#                     fi
                    
                    shift
                    ;;
                -s|--strain)
                    if [ "$flag_pbsim" = false ]; then
                        echo "only meaningful with pbsim"
                        exit 1
                    fi
                    cnt_strain="$2"
                    depth=$cnt_strain
                    
                    shift 2
                    ;;
                g|genome)
                    if [ "$flag_pbsim" = false ]; then
                        echo "only meaningful with pbsim"
                        exit 1
                    fi
                    flag_genome=0
                    shift 2
                    ;;
                --min_depth)
                    if [ "$flag_pbsim" = false ]; then
                        echo "only meaningful with pbsim"
                        exit 1
                    fi
                    
                    if [ "$flag_depth" = true ]; then
                        echo "these two options (--min_depth, --mean_depth) are in conflict"
                        exit 1
                    fi
                    flag_depth=true
                    type_depth=0
                    depth_pb="$2"
                    shift 2
                    ;;
                --mean_depth)
                    if [ "$flag_pbsim" = false ]; then
                        echo "only meaningful with pbsim"
                        exit 1
                    fi
                    if [ "$flag_depth" = true ]; then
                        echo "these two options (--min_depth, --mean_depth) are in conflict"
                        exit 1
                    fi
                    flag_depth=true
                    type_depth=1
                    depth_pb="$2"
                    shift 2
                    ;;

                --method)
                    if [ "$flag_pbsim" = false ]; then
                        echo "only meaningful with pbsim"
                        exit 1
                    fi
                    method="$2"
                    shift 2
                    ;;
                --model)
                    if [ "$flag_pbsim" = false ]; then
                        echo "only meaningful with pbsim"
                        exit 1
                    fi
                    model="$2"
                    shift 2
                    ;;
                --pass-num)
                    if [ "$flag_pbsim" = false ]; then
                        echo "only meaningful with pbsim"
                        exit 1
                    fi
                    pass_num="$2"
                    shift 2
                    ;;
                -h|--help)
                    echo "mcss sample"
                    echo "-i,--input: a file with fastq reads"
                    echo "-o,--output: output dir " 
                    echo "-n,--n_samples: count of simulate samples. default:1 "
                    echo "p,paired: paired fastq reads. default:single "
                    echo -e "-b, --batch: Divide into \"batch\" parts and run simultaneously"
                    echo "-db,--kraken_db: path of kraken database. defalut kraken/db"
                    echo "-asp,--abu_sp: minimum species abundance. default:0.001"
                    # echo "nl,new_len: change count of sample species."
                    echo "ml,multi_sample: more than one samples from the same env"
                    echo "acc,accurate: accurate mode to seek out the subtree"
                    echo "prol,prolific: prolific mode to seek out the subtree "
                    echo "Choose one between acc and prol, with acc as the default value in sample mode"
                    
                    echo "pbsim：execute  pbsim3"
                    echo "options for pbsim:"
                    echo "-s,--strain：number of strain"
                    echo "--min_depth,--mean_depth：min or average depth of species"
                    echo "--method: qshmm, errhmm. default:qshmm"
                    echo "--model QSHMM-ONT.model, QSHMM-RSII.model, ERRHMM-ONT.model, ERRHMM-RSII.model, ERRHMM-SEQUEL.model. default:QSHMM-RSII.model"
                    echo "--pass-num: if number of passes(pass-num)>1,multi-pass sequencing is performed to get hifi reads. default:10"
                    exit 0
                    ;;
                *)
                    echo "error: $1"
                    echo "mcss sample -h for usage"
                    exit 1
                    ;;
            esac
        done
        # echo "执行cmss sample的操作，参数 -i "
        # 在这里添加您的cmss get代码，可以使用$get_param
        if [ "$flag_in" = false ]; then
            echo "error : no input"
            exit 1
        fi
          
        if [ "$mode_flag" = true ]; then
            mode=2    
        else  
            mode=1
        fi
        
        if [ ! -d $output ]; then
            mkdir -p $output/community
        fi
        path_com=$output/community
        echo -e "\033[32minput file\033[0m :$input"
        echo -e "\033[32moutput file\033[0m :$output"
        echo $path_db
        
        if [ "$flag_depth" = true ]; then
                $depth=$depth_pb
        fi
        
        if [ ! -d $path_db ]; then
            cleanup() {
                echo "deleting the unfinished database..."
                rm -rf "$path_db"

                exit 1
            }
            trap cleanup INT
            read -p "download kraken ref database (GTDB 207) ? (y/n): " choice
            if [ "$choice" = "y" ] || [ "$choice" = "Y" ]; then
                # 用户输入是 "y" 或 "Y"，执行程序
                bash ./tools/install_kraken_db.sh
                # 在这里添加你要执行的命令
            elif [ "$choice" = "n" ] || [ "$choice" = "N" ]; then
                
                # 用户输入是 "n" 或 "N"，不执行程序
                exit 1
            else
                # 用户输入既不是 "y" 也不是 "n"
                exit 1
            fi
        fi
        # echo "batch" $batch
        if [ ! -d /dev/shm/kraken_db ]; then
            echo "cp db"
            cp -r $path_db /dev/shm/kraken_db
        fi
        path_map_db=/dev/shm/kraken_db
        echo "batch" $batch
        parts=($(get_part $input $mode_flag $batch))
        parts_len=${#parts[@]}
        # echo ${parts[*]}
        # echo $parts_len
        arry_1=()
        start_l=0 
        for ((i=0;i<parts_len;i++)); do
            # echo ${parts[i]}
            cmd_tmp="python $path_script/get_sample_sp.py $input $path_com $mode $dir_kraken $path_map_db $start_l ${parts[i]} $suf $suf_pre"
            #需要给cmd加一个括号，不然不行
            # arry_1+=("$cmd_tmp")
            arry_1=("${arry_1[@]}" "$cmd_tmp")
            start_l=${parts[i]}
        done
        for cmd in "${arry_1[@]}"; do
          echo "$cmd"
        done
        if [ ! -d $output/community/batch ]; then
            mkdir -p $output/community/batch
        fi
        path_data=$output/community/data
        
        if [ ! -d $output/community/kraken_res ]; then
            mkdir $output/community/kraken_res
        fi
        if [ ! -d $output/community/abu ]; then
            mkdir $output/community/abu
        fi
        if [ ! -d $output/species ]; then
            mkdir $output/species
        fi
        if [ ! -d $output/community/data/kraken ]; then
            mkdir -p $output/community/data/kraken
            mkdir $output/community/data/kraken/read_cnt
            mkdir $output/community/data/kraken/out_new
            mkdir $output/community/data/kraken/abundance
            
        fi
        

        
        run_background_commands "${arry_1[@]}"
        trap 'if [ "$normal_exit" = false ]; then kill_background_commands; fi' EXIT
        
        path_cut=$path_data/sp_cut
        if [ ! -d $path_cut ]; then
            mkdir $path_cut
        fi
        path_s1="$current_dir/data/s1_ancs_database.tsv"
        path_abu="$path_data/kraken/abundance/"
        path_out_new="$path_data/kraken/out_new"
        mode_2=false
        parts=($(get_part $path_abu $mode_2 $batch))
        parts_len=${#parts[@]}
        

        arry_2=()
        start_l=0
        for ((i=0;i<parts_len;i++)); do
            # echo ${parts[i]}
            cmd_tmp="python $path_script/2_get_sp.py $path_data $path_s1 $abu_sp $start_l ${parts[i]} "
            #需要给cmd加一个括号，不然不行
            arry_2+=("$cmd_tmp")
            start_l=${parts[i]}
        done

        run_background_commands "${arry_2[@]}"
        trap 'if [ "$normal_exit" = false ]; then kill_background_commands; fi' EXIT
        path_nw_label=$(which nw_labels)
        path_nw=$(dirname $path_nw_label)
        if [ ! -d $output/community/data/sample ]; then
            mkdir $output/community/data/sample
        fi
        if [ ! -d $output/community/tmp_data ]; then
            mkdir $output/community/tmp_data
        fi
        
            
        python $path_script/sim_tree_sample.py $path_com $path_tax $path_nw $current_dir $flag_acc $n_samples $same_env
        
        if [ "$flag_pbsim" = true ]; then
            abu_name=$(ls $output/community/abu/)
            echo $abu_name
            if [ "$abu_name" = "abu_env.pkl" ]; then
                sim_mode=0
            elif [ "$abu_name" = "abu_sample.pkl" ]; then
                sim_mode=1
            elif [ "$abu_name" = "abu_user.pkl" ]; then
                sim_mode=2
            fi
            echo "sim_mode "$sim_mode
            path_abu_new=$output/community/abu/$abu_name
            echo $path_abu_new
            model_new=$current_dir/data/pbsim_model/$model
       
            if [ ! -f $model_new ]; then
                echo "no model named $model"
                exit 1
            fi
            if [ ! -d $output/pbsim ]; then
                mkdir -p $output/pbsim
            fi
            #这里有两种思路，一组就是直接调用mcss pbsim 
            #一种就是python 
            #pyhon需要把pbsim的东西全都复刻
            #mcss需要判断输入了哪些
            python $path_script/get_strain.py $current_dir $output $cnt_strain $depth $type_depth $path_abu_new $sim_mode $method $model_new $pass_num $flag_genome
            #这个的第二个参数和pbsim单独选项里的参数不一样，因为那个直接输入的input就是sim文件，我们的输出才是
        fi
        ;;
########################################################################################### community        
    "community")
        # 解析cmss get的参数
        flag_in=false
        current_dir=$(dirname $(realpath $0))
        path_script=$current_dir/script
        output=$current_dir/sim_${date}
        flag_f=false
        flag_d=false
        flag_pbsim=false
        flag_a=false
        while [[ $# -gt 0 ]]; do
            case "$1" in
                "-f" | "--file")
                    if [ "$flag_d" = true ]; then
                        echo "these two options (-d, -f) are in conflict"
                        exit 1
                    fi
                    
                    flag_in=true
                    flag_f=true
                    input="$2"
                    input=$(realpath $input)
                    shift 2
                    ;;
                "-d" | "--dir")
                    if [ "$flag_f" = true ]; then
                        echo "these two options (-d, -f) are in conflict"
                        exit 1
                    fi
                    flag_in=true
                    flag_d=true
                    input="$2"
                    input=$(realpath $input)
                    shift 2
                    ;;
                
                -o|--output)
                
                    output="$2"
                    output=$(realpath $output)
                    shift 2
                    ;;

                -a|--abu)
                    path_abu="$2"
                    flag_a=true

                    shift 2
                    ;;
                "pbsim")
                    flag_genome=1
                    flag_pbsim=true
                    cnt_strain=1
                    depth=$cnt_strain
                    type_depth=0
                    flag_depth=false
                    method="qshmm"
                    model="QSHMM-RSII.model"
                    pass_num=10
#                     if [ -f $path_pbsim ]; then
#                         # path_pbsim=$(which pbsim3)
#                         echo " pbsim3 in tools "
                        
#                     elif which pbsim3 &>/dev/null; then
#                         echo "pbsim3 in path"
#                         path_pbsim=$(which pbsim3)
#                     else
#                         echo "please install pbsim3"
#                         exit 1
                    
#                     fi
                    
                    shift
                    ;;
                -s|--strain)
                    if [ "$flag_pbsim" = false ]; then
                        echo "only meaningful with pbsim"
                        exit 1
                    fi
                    cnt_strain="$2"
                    depth=$cnt_strain
                    
                    shift 2
                    ;;
                g|genome)
                    if [ "$flag_pbsim" = false ]; then
                        echo "only meaningful with pbsim"
                        exit 1
                    fi
                    flag_genome=0
                    shift 2
                    ;;
                --min_depth)
                    if [ "$flag_pbsim" = false ]; then
                        echo "only meaningful with pbsim"
                        exit 1
                    fi
                    
                    if [ "$flag_depth" = true ]; then
                        echo "these two options (--min_depth, --mean_depth) are in conflict"
                        exit 1
                    fi
                    flag_depth=true
                    type_depth=0
                    depth_pb="$2"
                    shift 2
                    ;;
                --mean_depth)
                    if [ "$flag_pbsim" = false ]; then
                        echo "only meaningful with pbsim"
                        exit 1
                    fi
                    if [ "$flag_depth" = true ]; then
                        echo "these two options (--min_depth, --mean_depth) are in conflict"
                        exit 1
                    fi
                    flag_depth=true
                    type_depth=1
                    depth_pb="$2"
                    shift 2
                    ;;

                --method)
                    if [ "$flag_pbsim" = false ]; then
                        echo "only meaningful with pbsim"
                        exit 1
                    fi
                    method="$2"
                    shift 2
                    ;;
                --model)
                    if [ "$flag_pbsim" = false ]; then
                        echo "only meaningful with pbsim"
                        exit 1
                    fi
                    model="$2"
                    shift 2
                    ;;
                --pass-num)
                    if [ "$flag_pbsim" = false ]; then
                        echo "only meaningful with pbsim"
                        exit 1
                    fi
                    pass_num="$2"
                    shift 2
                    ;;
                -h|--help)
                    echo "mcss community"
                    echo "-f,--file: a file with species"
                    echo "-d,--dir: path of species files"
                    echo "-o,--output: output dir "  
                    echo "-a,--abu：abundance file or abundance of a env"
                    echo "pbsim：execute  pbsim3"
                    echo "options for pbsim:"
                    echo "-s,--strain：number of strain"
                    echo "--min_depth,--mean_depth：min or average depth of species"
                    echo "--method: qshmm, errhmm. default:qshmm"
                    echo "--model QSHMM-ONT.model, QSHMM-RSII.model, ERRHMM-ONT.model, ERRHMM-RSII.model, ERRHMM-SEQUEL.model. default:QSHMM-RSII.model"
                    echo "--pass-num: if number of passes(pass-num)>1,multi-pass sequencing is performed to get hifi reads. default:10"
                    exit 0
                    ;;

                *)
                    echo "error: -$1"
                    echo " mcss community -h for usage"
                    exit 1
                    ;;
            esac
        done
        
        if [ "$flag_f" = true ]; then
            #n_samples=1
            multi=0
        else
            #n_samples=$(ls $input |wc -l )
            multi=1
        fi
        if [ ! -d $output ]; then
            mkdir -p $output/community
            mkdir -p $output/community/abu
        fi
        #这个还是要跑一下的，万一community有，而abu没有
        if [ ! -d $output/community/abu ]; then
            mkdir -p $output/community/abu
        else
            rm -r $output/community/abu
            mkdir -p $output/community/abu
        fi
        
        
        if [ "$flag_a" = false ]; then
            echo "abundance file missing"
            exit  1
        else
            if [[ " ${env_array[*]} " = *" $path_abu "* ]]; then
                flag_abu=0
                sim_mode=0
                path_abu_new=$current_dir/data/env/$path_abu/abu.pkl

                cp $path_abu_new $output/community/abu/abu_env.pkl
            else
                flag_abu=1
                sim_mode=2
                path_abu_new=$(realpath $path_abu)
                python $path_script/save_community_abu.py $input $path_abu_new $multi $output
                #执行完上面的代码就已经保存了结果
                
            fi
        fi
        
        
        if [ "$flag_in" = false ]; then
            echo "error : no input"
            exit 1
        fi
        
        if [ "$flag_depth" = true ]; then
                $depth=$depth_pb
        fi
        # echo $flag_f $multi

        path_com=$output/community
        echo -e "\033[32minput file\033[0m :$input"
        echo -e "\033[32moutput file\033[0m :$output"
        if [ ! -d $output/species ]; then
            mkdir $output/species
        fi
        
        python $path_script/sim_tree_community.py $input $output $multi
        
        if [ "$flag_pbsim" = true ]; then
            abu_name=$(ls $output/community/abu/)
            echo $abu_name
            if [ "$abu_name" = "abu_env.pkl" ]; then
                sim_mode=0
            elif [ "$abu_name" = "abu_sample.pkl" ]; then
                sim_mode=1
            elif [ "$abu_name" = "abu_user.pkl" ]; then
                sim_mode=2
            fi
            echo "sim_mode "$sim_mode
            path_abu_new=$output/community/abu/$abu_name
            echo $path_abu_new
            model_new=$current_dir/data/pbsim_model/$model
            
            if [ ! -f $model_new ]; then
                echo "no model named $model"
                exit 1
            fi
            if [ ! -d $output/pbsim ]; then
                mkdir -p $output/pbsim
            fi
            #这里有两种思路，一组就是直接调用mcss pbsim 
            #一种就是python 
            #pyhon需要把pbsim的东西全都复刻
            #mcss需要判断输入了哪些
            python $path_script/get_strain.py $current_dir $output $cnt_strain $depth $type_depth $path_abu_new $sim_mode  $method $model_new $pass_num $flag_genome
            #这个的第二个参数和pbsim单独选项里的参数不一样，因为那个直接输入的input就是sim文件，我们的输出才是
        fi
        
        
        # 在这里添加您的cmss get代码，可以使用$get_param
        ;;
############################################################################################# pbsim   
    "pbsim")
        flag_genome=1
        cnt_strain=1
        depth=$cnt_strain
        type_depth=0
        flag_depth=false
        flag_in=false
        method="qshmm"
        model="QSHMM-RSII.model"
        pass_num=10
        current_dir=$(dirname $(realpath $0))
        path_script=$current_dir/script
#         if [ -f $path_pbsim ]; then
#             # path_pbsim=$(which pbsim3)
#             echo " pbsim3 in tools "

#         elif which pbsim3 &>/dev/null; then
#             echo "pbsim3 in path"
#             path_pbsim=$(which pbsim3)
#         else
#             echo "please install pbsim3"
#             exit 1

#         fi
        # 解析cmss sim的参数
        while [[ $# -gt 0 ]]; do
            case "$1" in
                -i|--input)
                    input="$2"
                    input=$(realpath $input)
                    flag_in=true
                    shift 2
                    ;;
                -s|--strain)
                    cnt_strain="$2"
                    depth=$cnt_strain
                    shift 2
                    ;;
                g|genome)
                    if [ "$flag_pbsim" = false ]; then
                        echo "only meaningful with pbsim"
                        exit 1
                    fi
                    flag_genome=0
                    shift 2
                    ;;
                --min_depth)
                    if [ "$flag_depth" = true ]; then
                        echo "these two options (--min_depth, --mean_depth) are in conflict"
                        exit 1
                    fi
                    flag_depth=true
                    type_depth=0
                    depth_pb="$2"
                    shift 2
                    ;;
                --mean_depth)
                    if [ "$flag_depth" = true ]; then
                        echo "these two options (--min_depth, --mean_depth) are in conflict"
                        exit 1
                    fi
                    flag_depth=true
                    type_depth=1
                    depth_pb="$2"
                    shift 2
                    ;;
                # -d|--depth)
                #     depth="$2"
                #     shift 2
                #     ;;
                --method)
                    method="$2"
                    shift 2
                    ;;
                --model)
                    model="$2"
                    shift 2
                    ;;
                --pass-num)
                    pass_num="$2"
                    shift 2
                    ;;
                -h|--help)
                    echo "-i,--input: input dir " 
                    echo "-s,--strain：number of strain"
                    # echo "-a,--abu：abundance file or abundance of a env"
                    echo "--min_depth, --mean_depth：min or average depth of species"
                    # echo "-d,--depth: depth"
                    echo "--method: qshmm, errhmm. default:qshmm"
                    echo "--model QSHMM-ONT.model, QSHMM-RSII.model, ERRHMM-ONT.model, ERRHMM-RSII.model, ERRHMM-SEQUEL.model. default:QSHMM-RSII.model"
                    echo "--pass-num: if number of passes(pass-num)>1,multi-pass sequencing is performed to get hifi reads. default:10"
                    exit 0
                    ;;
                *)
                    echo "error: -$1"
                    echo " mcss pbsim -h for usage"
                    exit 1
                    ;;
            esac
        done
        if [ "$flag_in" = false ]; then
            echo "error : no input"
            exit 1
        fi
        if [ "$flag_depth" = true ]; then
                $depth=$depth_pb
        fi
        abu_name=$(ls $input/community/abu/)
        if [ "$abu_name" = "abu_env.pkl" ]; then
            sim_mode=0
        elif [ "$abu_name" = "abu_sample.pkl" ]; then
            sim_mode=1
        elif [ "$abu_name" = "abu_user.pkl" ]; then
            sim_mode=2
        fi
        echo "sim_mode "$sim_mode
        path_abu_new=$input/community/abu/$abu_name
        model_new=$current_dir/data/pbsim_model/$model
        if [ ! -f $model_new ]; then
            echo "no model named $model"
            exit 1
        fi
        if [ ! -d $output/pbsim ]; then
            mkdir -p $output/pbsim
        fi
        # echo $type_depth
        # if [ ! -d $input/sim_concat ]; then
        #     mkdir $input/sim_concat 
        # fi
        python $path_script/get_strain.py $current_dir $input $cnt_strain $depth $type_depth $path_abu_new $sim_mode $method $model_new $pass_num $flag_genome
        ;;
    "-h")
        echo "mcss.sh env: Environment-specific simulated data generation "
        echo "mcss.sh sample: Generating simulated data based on real data "
        echo "mcss.sh community: Commuinty specified by the user"
        echo "mcss.sh pbsim:  execute pbsim3"
        exit 0
        ;;
    *)
        echo "error: -$1"
        echo " mcss -h for usage"
        exit 1
        ;;
esac

exit 0
