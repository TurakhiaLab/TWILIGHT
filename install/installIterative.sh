# Install packages
conda install pip -y
pip install snakemake
pip install numpy

# Set up channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# Get system architecture
ARCH=$(uname -m)
OS=$(uname -s)

# Install tree inference tools with conda
conda install bioconda::fasttree -y # FastTree
conda install mafft -y    # MAFFT
conda install bioconda::raxml -y    # RAxML
conda install bioconda::iqtree -y   # IQ-Tree

if [[ "$OS" == "Linux" && "$ARCH" == "x86_64" ]]; then
    conda install bioconda::mashtree -y # MashTree
else
    echo "Skipping mashtree installation: not supported on $OS/$ARCH."
fi

