cd data
wget -O strain_download.tar.gz 'http://storage.live.com/items/B2A5926AB0601A6E!3084:/strain_download.tar.gz?authkey=!AIA94oJyoNOd104'
if [ -f strain_download.tar.gz ]; then
	tar xvzf strain_download.tar.gz
fi
