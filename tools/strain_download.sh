cd ../data
wget -O strain_download.tar.gz https://public.dm.files.1drv.com/y4mTUGiDqbUSdjUQkJjK1ZTC1R7WmnTLEThlWvjZv8vK9tjaFAf2i9U4ULrwM_RxHCsGVpXlsp25FiQbO5_tthZz2MNtlJoQlosZjZEc1eLPhZYH5LlOj3vLVRwL72qCdh9E_6qtumPWOHfQgOtTso-SbGDKlBwk9rZ4NaxOz8TcmECSw1I0Ogyv6P5pRMjiTf23MWFwWeaO4E54QXp_Z3rviju9lC83zd-mB7tgR5MJ3k?AVOverride=1
if [ -f strain_download.tar.gz ]; then
	tar xvzf strain_download.tar.gz
fi