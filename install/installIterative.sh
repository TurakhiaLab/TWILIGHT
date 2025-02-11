# Install packages
pip install snakemake
pip install numpy

# Install tree inference tools with conda
conda install bioconda::fasttree -y # FastTree
conda install bioconda::mashtree -y # MashTree
conda install bioconda::mafft -y    # MAFFT
conda install bioconda::raxml -y    # RAxML
conda install bioconda::iqtree -y   # IQ-Tree