# Install packages
conda install pip -y

# Set up channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

conda install -y snakemake ete3 numpy numba


# Get system architecture
ARCH=$(uname -m)
OS=$(uname -s)

# Iterative mode
conda install bioconda::fasttree -y # FastTree
conda install mafft -y              # MAFFT
conda install bioconda::raxml -y    # RAxML
conda install bioconda::iqtree -y   # IQ-Tree
git clone https://github.com/somme89/rapidNJ.git ## Rapid NJ
make -C rapidNJ/
# Placement mode
conda install bioconda::epa-ng -y   # EPA-NG
conda install bioconda::gappa -y    # GAPPA

if [[ "$OS" == "Linux" && "$ARCH" == "x86_64" ]]; then
    conda install bioconda::mashtree -y # MashTree
else
    echo "Skipping mashtree installation: not supported on $OS/$ARCH."
fi

