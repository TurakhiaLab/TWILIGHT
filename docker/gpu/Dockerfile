FROM --platform=linux/amd64 nvidia/cuda:12.2.0-devel-ubuntu20.04
ENV DEBIAN_FRONTEND noninteractive
RUN chmod 1777 /tmp

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

WORKDIR /home/

ENV CONDA_DIR /opt/conda

# Install conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
   /bin/bash ~/miniconda.sh -b -p /opt/conda

ENV PATH=$CONDA_DIR/bin:$PATH

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

# TWILIGHT
RUN git clone https://github.com/TurakhiaLab/TWILIGHT.git && \
    cd TWILIGHT && \
    bash ./install/buildTWILIGHT.sh cuda

WORKDIR /home/TWILIGHT

CMD ["bash"]
