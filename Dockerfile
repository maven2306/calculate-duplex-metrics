# syntax=docker/dockerfile:1
ARG R_VER=4.4.1
FROM rocker/r-ver:${R_VER}

# Install system dependencies  + Python (argparse/json are stdlib modules)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    python3 python3-pip python3-venv python3-distutils \
    zlib1g-dev \
    liblzma-dev \ 
    libbz2-dev \
    && rm -rf /var/lib/apt/lists/*

# Hint to findpython (optional but helpful)
ENV PYTHON=/usr/bin/python3
ENV PYTHON3=/usr/bin/python3

# Set working directory
WORKDIR /app

# Keep renv's library outside the project so mounts don't hide it
ENV RENV_PATHS_LIBRARY=/opt/renv/library
RUN mkdir -p /opt/renv/library

# Copy renv infra FIRST
COPY renv.lock ./
COPY renv/ ./renv/
COPY .Rprofile ./

# Restore R packages; fail build if argparse didn't install
RUN R -q -e "options(repos=c(CRAN='https://cloud.r-project.org')); install.packages('renv')" \
 && R -q -e "renv::restore(prompt = FALSE)"

# Copy the rest
COPY . .

# Set default command to bash for flexibility
CMD ["/bin/bash"]
