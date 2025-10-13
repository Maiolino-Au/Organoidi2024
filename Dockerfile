FROM ubuntu
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get -y update && apt-get -y upgrade && apt-get -y install r-base
RUN apt-get -y install gdebi-core
RUN apt-get -y install wget
RUN wget https://download2.rstudio.org/server/jammy/amd64/rstudio-server-2023.12.1-402-amd64.deb
RUN gdebi -n rstudio-server-2023.12.1-402-amd64.deb
RUN useradd rstudio -p "\$y\$j9T\$/.6YKeUOB4ifaPjuG/xaC1\$0162SW98NtTo5c6I7uXbwlNlKGuu9LTcUanCzz6DF/C" -d /home/rstudio -m
RUN apt-get update
RUN apt-get install -y \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg \
    lsb-release
RUN apt-get update && apt install -y libudunits2-dev libgdal-dev
RUN apt-get update
RUN apt-get -y install gfortran
RUN apt-get -y install build-essential
RUN apt-get -y install fort77
RUN apt-get -y install xorg-dev
RUN apt-get -y install liblzma-dev  libblas-dev gfortran
RUN apt-get -y install gobjc++
RUN apt-get -y install aptitude
RUN apt-get -y install libbz2-dev
RUN apt-get -y install libpcre3-dev
RUN aptitude -y install libreadline-dev
RUN apt-get -y install libcurl4-openssl-dev
RUN apt install -y build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev
RUN apt-get install -y libcurl4-openssl-dev libssl-dev
RUN apt-get install -y libgit2-dev
RUN apt-get install -y libharfbuzz-dev
RUN apt-get install -y libfribidi-dev
RUN apt-get install -y cmake
RUN apt-get install -y libcairo2-dev
RUN Rscript -e 'install.packages(c("devtools", "dplyr", "ggplot2", "tidyr", "stringr", "viridis", "ggthemes", "tidyverse", "ggsignif", "umap", "heatmap3", "plyr", "compareGroups", "dbscan", "reshape2", "msigdbr", "BiocManager", "tibble", "purrr", "magrittr", "ggplotify", "ggrepel", "igraph", "readxl", "pathviewr", "Seurat", "SeuratObject"), dependencies = TRUE)'
RUN Rscript -e 'BiocManager::install(c("limma", "Glimma", "edgeR", "scran", "fgsea", "DESeq2", "ensembldb", "BiocGenerics", "DelayedArray", "DelayedMatrixStats", "lme4", "S4Vectors", "SingleCellExperiment", "SummarizedExperiment", "MAST", "batchelor", "HDF5Array", "terra", "ggrastr", "topGO", "org.Hs.eg.db", "clusterProfiler", "enrichplot", "sf", "AnnotationHub", "vsn", "ComplexHeatmap", "ADImpute", "RRHO"), update = TRUE, ask = FALSE)'
RUN Rscript -e 'install.packages(c("RNAseqQC", "misgdbr", "fcors", "FactoMineR", "factoextra", "corrplot", "dendextend", "WGCNA", "pals"), dependencies = TRUE)'
RUN Rscript -e 'devtools::install_github("Biometeor/monocle3", dependencies = TRUE)'
RUN Rscript -e 'devtools::install_github("RRHO2/RRHO2", build_opts = c("--no-resave-data", "--no-manual"))'
RUN Rscript -e 'devtools::install_version("matrixStats", version="1.1.0")'
ENV TZ=UTC
RUN chmod +x /home/*

# Install JupyterLab
RUN apt update && apt install -y python3 python3-pip python3-venv
# create a virtual environment in which JupyterLab can be installed
RUN python3 -m venv /opt/venv
# Activate virtual environment and install JupyterLab
RUN /opt/venv/bin/pip install --upgrade pip && /opt/venv/bin/pip install jupyterlab
# Set the virtual environment as the default Python path
ENV PATH="/opt/venv/bin:$PATH"
# Make R visible to jupyter
RUN R -e "install.packages('IRkernel')" \
    R -e "IRkernel::installspec(user = FALSE)"

# Install the monocle3 branch of garnett
RUN R -e "BiocManager::install(c('org.Mm.eg.db', 'org.Hs.eg.db'))" \
    R -e "devtools::install_github('cole-trapnell-lab/garnett', ref='monocle3')"

RUN R -e "install.packages(c('openai', 'enrichR', 'pachwork', 'rtracklayer', 'tinytex'))" 
# RUN R -e "remotes::install_github('Winnie09/GPTCelltype')" 
RUN R -e "remotes::install_github('immunogenomics/presto')"
# To add: SingleR
RUN R -e "devtools::install_github('dviraran/SingleR')"

ENV SHELL=/bin/bash
CMD ["/bin/bash"]
#CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=9999", "--no-browser", "--allow-root", "--ServerApp.allow_origin='*'", "--ServerApp.token=''"]
