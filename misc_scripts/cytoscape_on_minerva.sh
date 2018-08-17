
## running cytoscape from Minerva
## open xquartz.app

## option-CLICK to paste
ssh -X richtf01@minerva.hpc.mssm.edu


echo $DISPLAY
# localhost:10.0

cd /sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq/wgcna/nuclear_results/cytoscape_modules
wget --no-check-certificate https://www.dropbox.com/s/j70ayepibhl5lc6/aracne_fromR_turquoise.txt?dl=0


# running as an interactive job
bsub -XF -P acc_chdiTrios -q expressalloc -n 1 -W 1:00 -m manda -Ip /bin/bash
module load cytoscape
bash $CYTOSCAPE_HOME/cytoscape.sh --network /sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq/wgcna/nuclear_results/cytoscape_modules/aracne_fromR_turquoise.txt

## option-click to PASTE into XQuartz
cp -r $CYTOSCAPE_HOME/.* /sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq/wgcna/cytoscape_cp


#################################

ssh -X interactive2
screen -R -D test_cytoscape
cd /sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq/wgcna/cytoscape_cp
export HOME=/sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq/wgcna/cytoscape_cp
module load java/1.8.0_141
module load java/1.7.0_60
bash cytoscape.sh --network /sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq/wgcna/nuclear_results/cytoscape_modules/aracne_fromR_turquoise.txt

## alternatively, interactive1

ssh -X interactive1



# module load cytoscape
# bash $CYTOSCAPE_HOME/cytoscape.sh --network /sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq/wgcna/nuclear_results/cytoscape_modules/cs_edges_black.txt
# echo $DISPLAY

