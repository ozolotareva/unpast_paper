# Use an official Python runtime as a base image
FROM python:3.8.16

# Set the working directory in the container to /app
WORKDIR /app

# Copy the current directory contents into the container at /app
COPY ./dist/ /app

# Install R and necessary libraries
RUN apt-get update && \
    apt-get install -y software-properties-common && \
    apt-get install -y libcurl4-openssl-dev libssl-dev && \
    apt-get install -y r-base && \
    apt-get update -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Upgrade pip
RUN pip install --upgrade pip

# Install the package using pip
RUN pip install numpy==1.22.3
RUN pip install /app/unpast-0.1.8.tar.gz

# Run the script to install R dependencies
RUN python -m unpast.install_r_dependencies > /app/install.log 2>&1

# Create a new user 'myuser' and set it as the default user
RUN useradd -m user
USER user

# Set the default command to "run_unpast"
ENTRYPOINT ["run_unpast"]
CMD ["-h"]
