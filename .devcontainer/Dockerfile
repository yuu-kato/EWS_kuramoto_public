FROM quay.io/jupyter/datascience-notebook:julia-1.10.5

USER root

WORKDIR /workdir
EXPOSE 8888

RUN apt-get update && apt-get install -y \
    texlive-publishers \
    latexmk \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

USER ${NB_UID}

# set default to the current project
ENV JULIA_PROJECT=@.
# enable multi-threading
ENV JULIA_NUM_THREADS=4

ENV JULIA_DEPOT_PATH=/home/${NB_USER}/.julia

# R packages needed to (re)generate networks/nets/*.csv from the `networkdata`
# package (networks/export_empirical_networks.R). `remotes` already ships
# with this base image; `igraph` and `networkdata` do not. Versions pinned to
# what this pipeline was validated against.
RUN Rscript -e 'remotes::install_version("igraph", version = "2.3.3", repos = "https://cloud.r-project.org")' \
    && Rscript -e 'remotes::install_github("schochastics/networkdata@521594049f725c3962c11496d5b0f48386ee84b2", upgrade = "never")'

# Python packages needed for TIPMOC detection (simulate/03_tipmoc_detect.py,
# calling tipmoc_lib/delta_aicc.py and fit_powerlaw.py directly). Versions
# pinned to what this pipeline was validated against.
RUN pip install --no-cache-dir \
    numpy==2.0.2 \
    scipy==1.14.1 \
    pandas==2.2.3
