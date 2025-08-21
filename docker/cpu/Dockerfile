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

RUN conda create -n twilight python=3.11 -y
RUN echo "source activate twilight" > ~/.bashrc
SHELL ["conda", "run", "-n", "twilight", "/bin/bash", "-c"]

WORKDIR /home/

# TWILIGHT
RUN git clone https://github.com/TurakhiaLab/TWILIGHT.git && \
    cd TWILIGHT && \
    bash ./install/buildTWILIGHT.sh && \
	bash ./install/installIterative.sh

WORKDIR /home/TWILIGHT

CMD ["bash"]
