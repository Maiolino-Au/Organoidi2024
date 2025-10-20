FROM satijalab/seurat:5.0.0

RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common dirmngr gpg curl build-essential \
    libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libharfbuzz-dev \
    libfribidi-dev make cmake gfortran libxt-dev liblapack-dev libblas-dev \
    sudo wget zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev pandoc git && \
    rm -rf /var/lib/apt/lists/*

# Python
RUN sudo apt update && sudo apt install -y python3 python3-pip python3-venv
RUN pip3 install anndata h5py numpy scipy pandas scanpy scib scvi muon

# R Packages
RUN R -e "install.packages('devtools')"
RUN R -e "devtools::install_github('dviraran/SingleR')"
RUN R -e "install.packages('tictoc')"
RUN R -e "BiocManager::install(c('zellkonverter', 'scuttle'))"
RUN R -e "remotes::install_github('mojaveazure/seurat-disk')"

RUN mkdir /scripts
# HNOA conversion
COPY SETTEMBRE/hnoa_conversion /scripts/

# Download the scripts for sep/oct_2025 analyses
WORKDIR /scripts_bersia_plots
COPY SETTEMBRE/6_run.r .
COPY SETTEMBRE/6_functions.r .


ENV SHELL=/bin/bash
CMD ["/bin/bash"]
#CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=9999", "--no-browser", "--allow-root", "--ServerApp.allow_origin='*'", "--ServerApp.token=''"]
