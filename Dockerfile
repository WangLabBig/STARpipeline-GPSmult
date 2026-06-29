# ============================================================
# Dockerfile for Multi-Ancestry Multi-Trait PRS Pipeline
# ============================================================

FROM ubuntu:20.04

LABEL maintainer="wangshaoqi@cncb.ac.cn"
LABEL description="Multi-ancestry, multi-trait polygenic risk score (PRS) pipeline"

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

# ============================================================
# 1. Switch apt to Alibaba mirror
# ============================================================
RUN sed -i 's|http://archive.ubuntu.com|http://mirrors.aliyun.com|g' /etc/apt/sources.list && \
    sed -i 's|http://security.ubuntu.com|http://mirrors.aliyun.com|g' /etc/apt/sources.list

# ============================================================
# 2. System dependencies + R in one layer
#    gpg and dirmngr are installed first so the CRAN key step works.
#    Ubuntu 20.04 focal-cran40 provides R >= 4.0
# ============================================================
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    curl \
    unzip \
    gzip \
    bzip2 \
    zstd \
    parallel \
    build-essential \
    cmake \
    gfortran \
    libgomp1 \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses5-dev \
    zlib1g-dev \
    libdeflate-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    python3 \
    python3-pip \
    python3-dev \
    ca-certificates \
    git \
    less \
    vim \
    datamash \
    gnupg \
    dirmngr \
    apt-transport-https \
    software-properties-common \
 && rm -rf /var/lib/apt/lists/*

RUN ln -s /usr/bin/python3 /usr/bin/python
# ============================================================
# 3. Install R 4.x from CRAN (focal-cran40)
#    gpg and dirmngr are now available from step 2.
#    Using Tsinghua mirror for the repo, key from keyserver.
# ============================================================
RUN gpg --keyserver keyserver.ubuntu.com \
        --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    gpg --export E298A3A825C0D65DFD57CBB651716619E084DAB9 \
        | gpg --dearmor -o /usr/share/keyrings/r-cran.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/r-cran.gpg] https://mirrors.tuna.tsinghua.edu.cn/CRAN/bin/linux/ubuntu focal-cran40/" \
        > /etc/apt/sources.list.d/r-cran.list && \
    apt-get update && \
    apt-get install -y --no-install-recommends r-base r-base-dev && \
    rm -rf /var/lib/apt/lists/*

# ============================================================
# 4. Install htslib (bgzip, tabix)
# ============================================================
ARG HTSLIB_VERSION=1.19
RUN cd /tmp && \
    wget -q https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    tar -xjf htslib-${HTSLIB_VERSION}.tar.bz2 && \
    cd htslib-${HTSLIB_VERSION} && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    rm -rf /tmp/htslib-${HTSLIB_VERSION}*

# ============================================================
# 5. Install PLINK2
# ============================================================
RUN cd /tmp && \
    wget -q https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64_20240105.zip && \
    unzip plink2_linux_x86_64_20240105.zip -d /usr/local/bin/ && \
    chmod +x /usr/local/bin/plink2 && \
    rm /tmp/plink2_linux_x86_64_20240105.zip

# ============================================================
# 6. Install csvtk
# ============================================================
ARG CSVTK_VERSION=0.30.0
RUN cd /tmp && \
    wget -q https://github.com/shenwei356/csvtk/releases/download/v${CSVTK_VERSION}/csvtk_linux_amd64.tar.gz && \
    tar -xzf csvtk_linux_amd64.tar.gz && \
    mv csvtk /usr/local/bin/ && \
    chmod +x /usr/local/bin/csvtk && \
    rm /tmp/csvtk_linux_amd64.tar.gz

# ============================================================
# 7. Install R packages
# ============================================================
COPY install_r_packages.R /tmp/install_r_packages.R
RUN Rscript /tmp/install_r_packages.R

COPY install_bigsnpr.R /tmp/install_bigsnpr.R
RUN Rscript /tmp/install_bigsnpr.R

# ============================================================
# 8. Install Python packages
# ============================================================
RUN pip3 install --no-cache-dir \
    pandas \
    numpy \
    scipy \
    docopt

# ============================================================
# 9. Copy pipeline scripts
# ============================================================
COPY scripts/ /opt/pipeline/scripts/
RUN chmod +x /opt/pipeline/scripts/* && \
    ln -sf /opt/pipeline/scripts/* /usr/local/bin/

COPY *.sh /opt/pipeline/

# ============================================================
# 10. Environment
# ============================================================
ENV PATH="/opt/pipeline/scripts:/usr/local/bin:${PATH}"

WORKDIR /data

CMD ["/bin/bash"]
