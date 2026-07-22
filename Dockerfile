# syntax=docker/dockerfile:1.7

FROM bioconductor/bioconductor_docker:RELEASE_3_23@sha256:1d871e1ca9cca76b220eb16e22677e728f4352f81a9ee91aaf29e24aea43e624

LABEL org.opencontainers.image.title="AllEnricher" \
      org.opencontainers.image.version="2.1.1" \
      org.opencontainers.image.source="https://github.com/zhangducsu/AllEnricher-v2" \
      org.opencontainers.image.licenses="MIT"

USER root

ENV VIRTUAL_ENV=/opt/allenricher-venv \
    PATH=/opt/allenricher-venv/bin:${PATH} \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    MPLCONFIGDIR=/tmp/matplotlib \
    XDG_CACHE_HOME=/tmp/.cache

RUN apt-get update \
    && apt-get install -y --no-install-recommends python3-venv \
    && rm -rf /var/lib/apt/lists/*

ADD --checksum=sha256:133d1c3abb2bc886795d8ceb1a689d7a61e4cb4ee61f08a9ab0096a58af0064d \
    https://bioconductor.org/packages/3.23/bioc/src/contrib/Archive/GSVA/GSVA_2.6.2.tar.gz \
    /tmp/GSVA_2.6.2.tar.gz

RUN Rscript -e 'BiocManager::install(c("fgsea", "GSVA", "ComplexHeatmap"), ask = FALSE, update = FALSE)'

RUN R CMD INSTALL /tmp/GSVA_2.6.2.tar.gz \
    && Rscript -e 'install.packages(c("circlize", "dplyr", "ggplot2", "scales", "tidyr"), repos = Sys.getenv("CRAN"))' \
    && Rscript -e 'stopifnot(as.character(getRversion()) == "4.6.1", as.character(BiocManager::version()) == "3.23", as.character(packageVersion("fgsea")) == "1.38.0", as.character(packageVersion("GSVA")) == "2.6.2")'

WORKDIR /opt/allenricher
COPY pyproject.toml README.md LICENSE ./
COPY allenricher ./allenricher
COPY docker/python-constraints.txt /tmp/python-constraints.txt

RUN python3 -m venv "${VIRTUAL_ENV}" \
    && python -m pip install --no-cache-dir --constraint /tmp/python-constraints.txt ".[api,visualization]" \
    && allenricher --version \
    && python -m pip freeze > /opt/allenricher-python-freeze.txt \
    && Rscript -e 'keys <- c("R", "Bioconductor", "fgsea", "GSVA"); values <- c(as.character(getRversion()), as.character(BiocManager::version()), as.character(packageVersion("fgsea")), as.character(packageVersion("GSVA"))); writeLines(paste(keys, values, sep = "\t"), "/opt/allenricher-r-versions.tsv")'

RUN mkdir -p /work /tmp/matplotlib /tmp/.cache \
    && chown -R rstudio:rstudio /work /tmp/matplotlib /tmp/.cache

WORKDIR /work
USER rstudio

EXPOSE 8000

ENTRYPOINT ["allenricher"]
CMD ["--help"]
