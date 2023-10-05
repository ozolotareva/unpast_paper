# Use the latest version of the Ubuntu base image
FROM ubuntu:latest

# Set the working directory
WORKDIR /data

# Set the environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Add the deadsnakes PPA and install Python 3.8 and pip
RUN apt-get update -y && \
    apt-get install -y software-properties-common && \
    apt-get install -y libcurl4-openssl-dev libssl-dev && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt-get update -y && \
    apt-get install -y python3.8 python3.8-venv python3.8-distutils python3.8-dev python3-pip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Update the default Python and pip symlinks
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.8 1 && \
    update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 1
    
# Add the R repository
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" >> /etc/apt/sources.list

# Import the public key for the R repository
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9

# Install the latest version of R
RUN apt-get update -y && \
    apt-get install -y r-base && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Upgrade pip
RUN pip install --upgrade pip

# Install WGCNA
RUN R -e "install.packages('BiocManager'); BiocManager::install('WGCNA')"

# Install unpast
RUN pip install unpast==0.1.5

# Create a directory for user data
RUN mkdir /user_data

# Create a new user 'myuser' and set it as the default user
RUN useradd -m myuser
USER myuser

# Set the default command to "run_unpast"
ENTRYPOINT ["run_unpast"]
CMD ["-h"]
