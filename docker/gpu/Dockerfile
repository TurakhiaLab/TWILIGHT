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

RUN conda create -n twilight python=3.11 -y
RUN echo "source activate twilight" > ~/.bashrc
SHELL ["conda", "run", "-n", "twilight", "/bin/bash", "-c"]

# TWILIGHT
RUN git clone https://github.com/TurakhiaLab/TWILIGHT.git && \
    cd TWILIGHT && \
    bash ./install/buildTWILIGHT.sh cuda && \
	bash ./install/installIterative.sh

WORKDIR /home/TWILIGHT

CMD ["bash"]
