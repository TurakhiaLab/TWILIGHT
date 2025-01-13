FROM ubuntu:20.04

ENV DEBIAN_FRONTEND noninteractive

# BASE
RUN apt-get update && \
    apt-get install -y wget git build-essential libboost-all-dev libtbb2 cmake protobuf-compiler

WORKDIR /app
RUN git clone https://github.com/TurakhiaLab/TWILIGHT.git && mkdir /app/TWILIGHT/build
WORKDIR /app/TWILIGHT/build
RUN wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz && \
    tar -xvzf 2019_U9.tar.gz && \
    cmake  -DTBB_DIR=${PWD}/oneTBB-2019_U9  -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake  .. && \
    make
ENV PATH=/app/TWILIGHT/build:${PATH}

CMD ["bash"]