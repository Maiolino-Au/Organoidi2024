FROM ghcr.io/maiolino-au/monocle:latest
# To add: SingleR
RUN R -e "devtools::install_github('dviraran/SingleR')"
RUN R -e "install.packages('tictoc')"
RUN R -e "BiocManager::install('zellkonverter')"

# Download the scripts for sep/oct_2025 analyses
RUN mkdir -p /scripts_bersia_plots && \
  cd /scripts_bersia_plots && \
  curl -O https://raw.githubusercontent.com/Maiolino-Au/ProgrAppBioinfo_Exam/main/Scripts/Script_PAB_exam_Maiolino.R && \
  curl -O https://raw.githubusercontent.com/Maiolino-Au/ProgrAppBioinfo_Exam/main/Scripts/Maiolino_Au.Rmd 

ENV SHELL=/bin/bash
CMD ["/bin/bash"]
#CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=9999", "--no-browser", "--allow-root", "--ServerApp.allow_origin='*'", "--ServerApp.token=''"]
