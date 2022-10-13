FROM ubuntu:20.04

LABEL Maintainer="fabio.malta.de.sa.patroni@uni-hamburg.de"

USER root

EXPOSE 80

# NOTE basic set up
ENV DEBIAN_FRONTEND=noninteractive
ENV LANG "en_US.UTF-8"
ENV LC_ALL "en_US.UTF-8"
ENV LC_CTYPE "en_US.UTF-8"
ENV PYTHONDONTWRITEBYTECODE "1"
ENV HOME /root

SHELL ["/bin/bash", "--login", "-c"]

RUN apt update && apt install -y locales tzdata && \
    apt upgrade -y && locale-gen "en_US.UTF-8" && \
    dpkg-reconfigure --frontend noninteractive tzdata locales && \
    apt update --fix-missing

# NOTE installing dependencies
RUN apt install -y \
    htop git make build-essential libssl-dev zlib1g-dev libbz2-dev nano \
    libreadline-dev libsqlite3-dev wget curl llvm libncursesw5-dev \
    xz-utils tk-dev libxml2-dev libxmlsec1-dev libffi-dev liblzma-dev \
    bzip2 g++ gcc gfortran libcairo2-dev libcurl4-openssl-dev libjpeg-dev \
    libpango1.0-dev libpcre++-dev libpcre2-dev libpng-dev libreadline6-dev \
    libx11-dev libxt-dev libzstd-dev screen texinfo texlive xvfb \
    texlive-fonts-extra xorg-dev software-properties-common dirmngr \
    apt-transport-https --fix-missing

RUN wget https://packages.microsoft.com/config/ubuntu/20.04/packages-microsoft-prod.deb -O packages-microsoft-prod.deb && \
    dpkg -i packages-microsoft-prod.deb && \
    rm packages-microsoft-prod.deb && \
    apt update && apt install -y dotnet-sdk-6.0

# NOTE create and set workdir
RUN mkdir -p $HOME/projects/data/real_data
RUN mkdir -p $HOME/projects/scripts
RUN mkdir -p $HOME/projects/outputs

WORKDIR $HOME/projects

RUN git clone https://github.com/ozolotareva/DESMOND2.git && \
    git clone https://github.com/ozolotareva/DESMOND.git

WORKDIR $HOME

# NOTE preparing miniconda3
RUN apt update -qq && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN chmod +x Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3

RUN echo "export PATH=$HOME/miniconda3/bin:$PATH" >> $HOME/.bashrc && \
    echo "conda init bash" >> $HOME/.bashrc && \
    echo "export R_HOME=/root/miniconda3/envs/renv/lib/R" >> $HOME/.bashrc && \
    env R_HOME=/root/miniconda3/envs/renv/lib/R && \
    echo ".libPaths(c(.libPaths()[grep('/envs/', .libPaths())], .libPaths()[-grep('/envs/', .libPaths())]))" > $HOME/.Rprofile

RUN $HOME/miniconda3/bin/conda update -n base -y -c defaults conda && \
    $HOME/miniconda3/bin/conda config --add channels defaults && \
    $HOME/miniconda3/bin/conda config --add channels bioconda && \
    $HOME/miniconda3/bin/conda config --add channels conda-forge

# NOTE python and R instalation
COPY desmond2_Python.yml .
COPY desmond2_R.yml .

RUN $HOME/miniconda3/bin/conda install -n base -y mamba

RUN $HOME/miniconda3/bin/mamba env create -f desmond2_Python.yml && \
    $HOME/miniconda3/bin/mamba env create -f desmond2_R.yml && \
    $HOME/miniconda3/bin/mamba install -n renv bioconductor-genomeinfodbdata \
    bioconductor-org.hs.eg.db bioconductor-simpintlists

# NOTE post instalation setup
RUN rm Miniconda3-latest-Linux-x86_64.sh && \
    echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc
RUN rm -rf /var/lib/apt/lists/*

CMD ["bash"]
