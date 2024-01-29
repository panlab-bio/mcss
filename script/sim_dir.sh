# Create a folder to store the results.
output=$1

if [ ! -d $output ]; then
    mkdir -p $output/community
fi

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
path_cut=$path_data/sp_cut
if [ ! -d $path_cut ]; then
    mkdir $path_cut
fi
if [ ! -d $output/community/data/sample ]; then
    mkdir $output/community/data/sample
fi
if [ ! -d $output/community/tmp_data ]; then
    mkdir $output/community/tmp_data
fi
