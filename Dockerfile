# Using the ubuntu:20.04 base image for compatibility with the ICESat-2 SlideRule software
# Some code to add conda from https://github.com/conda-forge/miniforge-images/blob/master/ubuntu/Dockerfile
# This file builds the docker image for the ICESat-2 tracks project. Includes the ICESat-2 SlideRule tool, conda, and the python3.11 environment for the project with dependencies provided in the environment_small.yml file.
FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive
ARG MINIFORGE_NAME=Miniforge3
ARG MINIFORGE_VERSION=23.3.1-0
ARG TARGETPLATFORM

ENV CONDA_DIR=/opt/conda
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH=${CONDA_DIR}/bin:${PATH}
ENV ENVNAME=2021-icesat2-tracks
RUN apt-get update && \
    apt-get install -y git curl wget

# Install sliderule dependencies
RUN apt-get install -y build-essential libreadline-dev liblua5.3-dev zlib1g-dev cmake
RUN git clone https://github.com/ICESat2-SlideRule/sliderule.git
WORKDIR /sliderule/
RUN make config && \
    make && \
    make install

# Install conda
RUN apt-get update > /dev/null && \
    apt-get install --no-install-recommends --yes \
    bzip2 ca-certificates \
    tini \
    > /dev/null && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    wget --no-hsts --quiet https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/${MINIFORGE_NAME}-${MINIFORGE_VERSION}-Linux-$(uname -m).sh -O /tmp/miniforge.sh && \
    /bin/bash /tmp/miniforge.sh -b -p ${CONDA_DIR} && \
    rm /tmp/miniforge.sh && \
    conda clean --tarballs --index-cache --packages --yes && \
    find ${CONDA_DIR} -follow -type f -name '*.a' -delete && \
    find ${CONDA_DIR} -follow -type f -name '*.pyc' -delete && \
    conda clean --force-pkgs-dirs --all --yes && \
    echo ". ${CONDA_DIR}/etc/profile.d/conda.sh && conda activate base" >> /etc/skel/.bashrc && \
    echo ". ${CONDA_DIR}/etc/profile.d/conda.sh && conda activate base" >> ~/.bashrc

# Install conda environment
COPY environment_py37_small.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml

ENTRYPOINT ["tini", "--"]
CMD ["/bin/bash"]
