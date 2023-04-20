FROM ubuntu:20.04

# Ensure /bin/bash is the default shell
SHELL ["/bin/bash", "-c"]

RUN apt -y update \
    && apt -y install gcc make autoconf git zip python3 python3-pip libbz2-dev liblzma-dev curl wget

RUN git clone https://github.com/mrcepid-rap/general_utilities.git \
    && cd general_utilities \
    && pip3 install .

ADD https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_x86_64_20220603.zip plink2.zip

RUN mkdir plink2 \
    && unzip plink2.zip -d plink2/ \
    && ln plink2/plink2 /usr/bin/ \
    && rm plink2.zip

RUN pip3 install simuPOP pysam

ADD https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed bigBedToBed

RUN chmod +x bigBedToBed \
    && mv bigBedToBed /usr/bin/

