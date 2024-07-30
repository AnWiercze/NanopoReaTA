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
RUN apt-get install -y bedops
RUN apt-get update
RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/*

# Install Human reference
RUN mkdir Human_reference_data
RUN wget https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_GENCODE_V42_Comprehensive.bed.gz/download -O hg38_GENCODE_V42_Comprehensive.bed.gz
RUN gunzip hg38_GENCODE_V42_Comprehensive.bed.gz

RUN wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz
RUN gunzip GRCh38.primary_assembly.genome.fa.gz

RUN wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.transcripts.fa.gz
RUN gunzip gencode.v43.transcripts.fa.gz

RUN wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.primary_assembly.annotation.gtf.gz
RUN gunzip gencode.v43.primary_assembly.annotation.gtf.gz

RUN mv hg38_GENCODE_V42_Comprehensive.bed Human_reference_data
RUN mv GRCh38.primary_assembly.genome.fa Human_reference_data
RUN mv gencode.v43.transcripts.fa Human_reference_data
RUN mv gencode.v43.primary_assembly.annotation.gtf Human_reference_data

# Install Mouse reference
RUN mkdir Mouse_reference_data
RUN wget https://sourceforge.net/projects/rseqc/files/BED/Mouse_Mus_musculus/GRCm39_GENCODE_VM27.bed.gz/download -O GRCm39_GENCODE_VM27.bed.gz
RUN gunzip GRCm39_GENCODE_VM27.bed.gz

RUN wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/GRCm39.primary_assembly.genome.fa.gz
RUN gunzip GRCm39.primary_assembly.genome.fa.gz

RUN wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.transcripts.fa.gz
RUN gunzip gencode.vM32.transcripts.fa.gz

RUN wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.primary_assembly.annotation.gtf.gz
RUN gunzip gencode.vM32.primary_assembly.annotation.gtf.gz

RUN mv GRCm39_GENCODE_VM27.bed Mouse_reference_data 
RUN mv GRCm39.primary_assembly.genome.fa Mouse_reference_data 
RUN mv gencode.vM32.transcripts.fa Mouse_reference_data 
RUN mv gencode.vM32.primary_assembly.annotation.gtf Mouse_reference_data 

# Install Yeast reference
RUN mkdir Yeast_reference_data
RUN wget https://ftp.ensembl.org/pub/release-111/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz 
RUN gunzip Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
RUN wget https://ftp.ensembl.org/pub/release-111/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz
RUN gunzip Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz
RUN wget https://ftp.ensembl.org/pub/release-111/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.111.gtf.gz
RUN gunzip Saccharomyces_cerevisiae.R64-1-1.111.gtf.gz
RUN mv Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa Yeast_reference_data
RUN mv Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa Yeast_reference_data
RUN cat Saccharomyces_cerevisiae.R64-1-1.111.gtf |  sed s/gene_biotype/gene_type/g | sed s/transcript_biotype/transcript_type/g > ./Yeast_reference_data/Saccharomyces_cerevisiae.R64-1-1.111.gtf
RUN rm Saccharomyces_cerevisiae.R64-1-1.111.gtf

RUN mkdir Reference_data
RUN mv Human_reference_data Reference_data
RUN mv Mouse_reference_data Reference_data 
RUN mv Yeast_reference_data Reference_data

RUN apt-get update && apt-get install -y curl
RUN curl https://sh.rustup.rs -sSf | bash -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"
RUN cargo install gxf2bed
RUN gxf2bed --input /Reference_data/Yeast_reference_data/Saccharomyces_cerevisiae.R64-1-1.111.gtf --output /Reference_data/Yeast_reference_data/Saccharomyces_cerevisiae.R64-1-1.111.bed

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

RUN conda init bash
RUN source /root/.bashrc
#RUN conda run -n nanoporeata Rscript install.R
SHELL ["conda", "run", "--no-capture-output", "-n", "nanoporeata", "Rscript", "install.R"]


ENV PORT=8080

EXPOSE 8080

ENV SHINY_LOG_STDERR=1


ENTRYPOINT ["conda","run", "--no-capture-output","-n","nanoporeata", "Rscript", "app_docker.R"]


