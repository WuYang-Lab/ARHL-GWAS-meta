metal="/public/home/shilulu/software/METAL/build/bin/metal"

# EAS meta
qsubshcom "$metal EAS_MVP_BBJ_metal.conf" 1 100G METAL 1:00:00 ""

# EUR Meta
qsubshcom "$metal EUR_MVP_Trpchevska_De-Angelis_metal.conf" 1 100G METAL 1:00:00 ""

# all meta
qsubshcom "$metal All_MVP_Trpchevska_De-Angelis_BBJ_metal.conf" 1 100G METAL 90:00:00 ""