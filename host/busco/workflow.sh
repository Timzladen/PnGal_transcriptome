####BUSCO workflow for SCG completeness estimation ####
#setup env
conda create -n busco -c conda-forge -c bioconda busco=5.7.1
conda activate busco

bash run_busco_all.sh

source busco_scg.R
