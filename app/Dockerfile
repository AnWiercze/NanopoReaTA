FROM ubuntu:22.04


# Install base utilities
RUN apt-get update
RUN apt-get install -y build-essential autoconf libtool
RUN apt-get install -y wget
RUN apt-get install -y zip 
RUN apt-get install -y apt-utils
RUN apt-get install -y firefox
RUN apt-get update
RUN apt-get install -y gcc-12 gcc-12-base gcc-12-doc g++-12
RUN apt-get install -y libstdc++-12-dev libstdc++-12-doc
RUN apt-get update
RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/*

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
     /bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

RUN wget https://github.com/AnWiercze/NanopoReaTA/archive/refs/heads/master.zip

RUN unzip master.zip

WORKDIR ./NanopoReaTA-master/app/

RUN conda env create -f ./requirements_nanoporeata.yml

RUN conda run -n nanoporeata Rscript install.R

RUN conda init bash

ENV PORT=8080

EXPOSE 8080

ENV SHINY_LOG_STDERR=1

ENTRYPOINT ["conda","run", "--no-capture-output","-n","nanoporeata", "Rscript", "app_docker.R"]

