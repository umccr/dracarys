FROM condaforge/mambaforge:22.9.0-2 as conda
LABEL maintainer="https://github.com/pdiakumis"

# install conda-lock
RUN mamba config \
      --set always_yes yes \
      --set always_copy yes && \
    mamba install \
      -c conda-forge \
      -c nodefaults \
      conda-lock && \
    mamba clean --all --force-pkgs-dirs

ARG DRACARYS_LOCK="dracarys-linux-64.lock"
ARG ENV_NAME="dracarys_env"
COPY ./conda/env/lock/${DRACARYS_LOCK} .
RUN conda-lock install --name ${ENV_NAME} --file ${DRACARYS_LOCK} && \
    mamba clean --all --force-pkgs-dirs

ARG MAMBA_PREFIX="/opt/conda"
ENV PATH="${MAMBA_PREFIX}/envs/${ENV_NAME}/bin:${PATH}"
ENV CONDA_PREFIX="${MAMBA_PREFIX}/envs/${ENV_NAME}"

CMD [ "dracarys.R" ]
