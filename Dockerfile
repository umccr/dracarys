FROM ubuntu:20.04
LABEL maintainer="https://github.com/pdiakumis"

RUN apt-get update && \
    apt-get install --yes --no-install-recommends \
    bash bzip2 curl git less vim wget zip ca-certificates && \
    apt-get clean && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*
