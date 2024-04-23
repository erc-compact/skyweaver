FROM nvidia/cuda:12.4.1-devel-ubuntu22.04

LABEL maintainer="compact-erc"

RUN apt-get update && \
    apt-get install -y libboost-all-dev

# need to install PSRDADA_CPP

WORKDIR /usr/src/skyweaver

COPY . . 

