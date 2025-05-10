FROM --platform=linux/amd64 condaforge/miniforge3
ENV DEBIAN_FRONTEND noninteractive

# BASE
RUN apt-get update && \
    apt-get install -y \
       	wget \
	    git \
        vim \
	    build-essential \
	    libboost-all-dev \
		libtbb-dev \
	    cmake && \
	apt-get clean

RUN conda create -n twilight -y
RUN echo "source activate twilight" > ~/.bashrc
SHELL ["conda", "run", "-n", "twilight", "/bin/bash", "-c"]

RUN pip install snakemake
RUN pip install numpy

RUN conda install bioconda::fasttree -y --quiet # FastTree
RUN conda install bioconda::mashtree -y --quiet # MashTree
RUN conda install bioconda::mafft -y --quiet    # MAFFT
RUN conda install bioconda::raxml -y --quiet    # RAxML
RUN conda install bioconda::iqtree -y --quiet   # IQ-Tree

WORKDIR /home/

# TWILIGHT
RUN git clone https://github.com/TurakhiaLab/TWILIGHT.git && \
    cd TWILIGHT && \
    bash ./install/buildTWILIGHT.sh

WORKDIR /home/TWILIGHT

CMD ["bash"]
