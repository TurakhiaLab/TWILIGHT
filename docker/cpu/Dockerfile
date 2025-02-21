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
	libtbb2 \
	cmake \
	protobuf-compiler && \
    apt-get clean

RUN conda create -n twilight -y
RUN echo "source activate twilight" > ~/.bashrc
SHELL ["conda", "run", "-n", "twilight", "/bin/bash", "-c"]

RUN pip install snakemake
RUN pip install numpy

RUN conda install bioconda::fasttree -y # FastTree
RUN conda install bioconda::mashtree -y # MashTree
RUN conda install bioconda::mafft -y    # MAFFT
RUN conda install bioconda::raxml -y    # RAxML
RUN conda install bioconda::iqtree -y   # IQ-Tree

WORKDIR /home/

# TWILIGHT
RUN git clone https://github.com/TurakhiaLab/TWILIGHT.git && \
    mkdir ./TWILIGHT/build && \
    cd ./TWILIGHT/build && \ 
    wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz && \
    tar -xvzf 2019_U9.tar.gz && \
    cmake  -DTBB_DIR=${PWD}/oneTBB-2019_U9  -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake  .. && \
    make -j

WORKDIR /home/TWILIGHT

CMD ["bash"]
