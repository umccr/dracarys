FROM condaforge/mambaforge:4.12.0-2 as conda
LABEL maintainer="https://github.com/pdiakumis"

# install conda-lock
RUN mamba config \
      --set always_yes yes \
      --set always_copy yes && \
    mamba install \
      -c conda-forge \
      -c nodefaults \
      conda-lock==1.0.5 && \
    mamba clean --all --force-pkgs-dirs

ARG ENV_NAME="dracarys_env"
COPY ./conda/env/lock/conda-lock.yml .
RUN conda-lock install --name ${ENV_NAME} conda-lock.yml && \
    mamba clean --all --force-pkgs-dirs && \
    rm conda-lock.yml


# CLEANUP

ARG D1="/opt/conda"
ARG D2="/opt/conda/envs/${ENV_NAME}"
RUN echo "HELLO"
RUN rm -fr /usr/bin/git-shell \
     /usr/bin/git \
     /usr/bin/perl \
     /usr/share/perl/ \
     /usr/share/perl5/ \
     /usr/share/doc/ \
     /opt/conda/lib/libicudata.so.70.1 \
     /opt/conda/lib/libstdc++.so.6.0.30 \
     /opt/conda/lib/libmamba.so.2.0.0 \
     /opt/conda/share/doc/ \
     /opt/conda/envs/dracarys_env/share/doc/ \
     /opt/conda/bin/nghttpx \
     /opt/conda/bin/openssl \
     /opt/conda/bin/mamba-package \
     /opt/conda/bin/bsdtar \
     /opt/conda/conda-meta \
     /opt/conda/envs/dracarys_env/conda-meta \
     /opt/conda/include \
     /opt/conda/envs/dracarys_env/include \
     /opt/conda/lib/libpython3.9.so.1.0 \
     /opt/conda/envs/dracarys_env/lib/libpython3.10.so.1.0 \
     /opt/conda/bin/x86_64-conda-linux-gnu-ld \
     /opt/conda/envs/dracarys_env/bin/x86_64-conda-linux-gnu-ld \
     /opt/conda/bin/sqlite3 \
     /opt/conda/envs/dracarys_env/bin/sqlite3 \
     /opt/conda/bin/openssl \
     /opt/conda/envs/dracarys_env/bin/openssl \
     /opt/conda/lib/python3.9/site-packages/pip \
     /opt/conda/lib/python3.9/idlelib \
     /opt/conda/lib/python3.9/ensurepip \
     /opt/conda/envs/dracarys_env/lib/python3.10/site-packages \
     /opt/conda/envs/dracarys_env/lib/python3.10/idlelib \
     /opt/conda/envs/dracarys_env/lib/python3.10/ensurepip \
     /opt/conda/share/terminfo \
     /opt/conda/envs/dracarys_env/share/terminfo
#RUN find -name '__pycache__' -type d -exec rm -rf '{}' '+'

ARG MAMBA_PREFIX="/opt/conda"
ENV PATH="${MAMBA_PREFIX}/envs/${ENV_NAME}/bin:${PATH}"
ENV CONDA_PREFIX="${MAMBA_PREFIX}/envs/${ENV_NAME}"

CMD [ "dracarys.R" ]
