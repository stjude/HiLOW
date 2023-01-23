#################################################################
# Dockerfile
#
# Software:         HiLOW tool dependencies
# Software Version: 1.0
# Website:          https://github.com/stjude/HiLOW
# Provides:         All dependencies needed to run HiLOW
# Base Image:       ghcr.io/stjude/abralab/binf-base:1.1.0
# Build Cmd:        docker build -t ghcr.io/stjude/abralab/hilow:latest -t ghcr.io/stjude/abralab/hilow:v1.0 -f hilow-dockerfile .
# Pull Cmd:         docker pull ghcr.io/stjude/abralab/hilow:v1.0
# Run Cmd:          docker run --rm -ti ghcr.io/stjude/abralab/hilow:v1.0
#################################################################

## TOOLS INCLUDED
#	HiCPro dependencies
#	HiCPro v3.1.0 (with minor edit)
#	JUICEBOX (Juicer Tools) v1.22.01
#	HICHIP-PEAKS v0.1.2 (with minor edit)
#	BEDTOOLS v2.30.0
#	KENTUTIL bedGraphToBigWig
#	FitHiChIP v10
#	MACS2 v2.2.7.1
#

FROM ghcr.io/stjude/abralab/bedtools:v2.30.0 as frombed
FROM ghcr.io/stjude/abralab/kentutils:latest as fromkent

FROM ghcr.io/stjude/abralab/binf-base:1.1.0

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

LABEL description="Docker image containing all requirements for HiLOW"

MAINTAINER Modupeore Adetunji "modupeore.adetunji@stjude.org"

## Install system tools
RUN apt-get update \
  && apt-get -y upgrade \
  && apt-get -y install pkg-config \
  libssl-dev curl \
  libcurl4-openssl-dev \
  bc \
  bzip2 \
  gcc \
  g++ && apt-get clean


## Install miniconda.
RUN wget https://repo.continuum.io/miniconda/Miniconda3-py37_4.8.2-Linux-x86_64.sh -O ~/anaconda.sh
RUN bash ~/anaconda.sh -b -p /usr/local/anaconda
RUN rm ~/anaconda.sh
ENV PATH /usr/local/anaconda/bin:$PATH

## Install all dependencies using conda
COPY environment.yml /
RUN conda env create -f /environment.yml && \
    conda install -y -c bioconda iced && \
    conda clean -a
ENV PATH /usr/local/anaconda/envs/HiC-Pro_v3.1.0/bin:$PATH

## Install HiCPro
ENV H_VERSION 3.1.0
RUN cd /tmp && \
    echo "master.zip" | wget https://github.com/nservant/HiC-Pro/archive/refs/tags/v${H_VERSION}.zip && \
    unzip v${H_VERSION}.zip && \
    cd HiC-Pro-${H_VERSION}  && \ 
    sed -i "s/TORQUE/LSF/" config-install.txt && \
    wget https://raw.githubusercontent.com/stjude/HiLOW/main/modifications/scripts/hic.inc.sh && \
    cp hic.inc.sh scripts && \
    make configure prefix=/ && \
    make install && \
    cd .. && \
    rm -fr HiC-Pro*
ENV PATH="/HiC-Pro-3.1.0/bin/:${PATH}"
ENV PATH="/HiC-Pro-3.1.0/bin/utils/:${PATH}"

## Install JuicerTools
RUN cd /tmp && \
    wget https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar && \
    mv juicer_tools_1.22.01.jar /usr/local/bin
    
## Get configuration file
RUN mkdir /data && cd /data && wget https://raw.githubusercontent.com/stjude/HiLOW/main/modifications/config-hicpro.txt

## ADD DEPENDENCIES TO BASE IMAGE
COPY --from=frombed /opt/bedtools2/bin /usr/local/bin
COPY --from=fromkent /usr/local/bin/bedGraphToBigWig /usr/local/bin

## INSTALL HICHIP-PEAKS
RUN pip install hichip-peaks
RUN cd /usr/local/anaconda/envs/HiC-Pro_v3.1.0/lib/python3.8/site-packages/hichip_peaks/ && \
    rm sparse_to_peaks.py && \
    wget https://raw.githubusercontent.com/stjude/HiLOW/main/modifications/sparse_to_peaks.py
RUN /HiC-Pro_3.1.0/bin/HiC-Pro -h

# Install required R packages
RUN R -e "install.packages(c('BiocManager', 'optparse', 'ggplot2', 'data.table', 'splines', 'fdrtool', 'parallel', 'tools', 'plyr', 'dplyr'), quitely=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e 'BiocManager::install()'
RUN R -e "BiocManager::install('GenomicRanges');"
RUN R -e "BiocManager::install('edgeR');"

# Install MACS2
RUN pip install macs2

# Install FitHiChIP v10
RUN cd / && wget https://github.com/ay-lab/FitHiChIP/archive/refs/tags/10.0.zip && \
    unzip 10.0.zip && \
    mv FitHiChIP-10.0 FitHiChIP 
RUN cd /FitHiChIP && sed -i 's|HiCProExec=`which HiC-Pro`|HiCProExec=`which /HiC-Pro_3.1.0/bin/HiC-Pro`|g' FitHiChIP_HiCPro.sh
RUN pip install networkx

RUN chmod -R 777 /root

# Add HILOW

