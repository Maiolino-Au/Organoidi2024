FROM ghcr.io/maiolino-au/monocle:latest
# To add: SingleR
RUN R -e "devtools::install_github('dviraran/SingleR')"
RUN R -e "install.packages('tictoc')"
RUN R -e "BiocManager::install(c('zellkonverter', 'scuttle'))"
#; library(reticulate); reticulate::install_miniconda()"
# RUN /root/.local/share/r-miniconda/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
# RUN /root/.local/share/r-miniconda/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
# RUN R -e "library(reticulate); reticulate::py_install(c('anndata', 'h5py', 'numpy', 'scipy'))"


# Download the scripts for sep/oct_2025 analyses
WORKDIR /scripts_bersia_plots
COPY SETTEMBRE/6_run.r .
COPY SETTEMBRE/6_functions.r .


ENV SHELL=/bin/bash
CMD ["/bin/bash"]
#CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=9999", "--no-browser", "--allow-root", "--ServerApp.allow_origin='*'", "--ServerApp.token=''"]
