FROM ghcr.io/maiolino-au/monocle:latest
# To add: SingleR
RUN R -e "devtools::install_github('dviraran/SingleR')"

ENV SHELL=/bin/bash
CMD ["/bin/bash"]
#CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=9999", "--no-browser", "--allow-root", "--ServerApp.allow_origin='*'", "--ServerApp.token=''"]
