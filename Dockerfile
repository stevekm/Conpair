# replicate continuumio/miniconda2:4.7.12 but need to install everything in a different location because /opt is not available for us
FROM debian:stretch
# https://hub.docker.com/r/continuumio/miniconda3/dockerfile
# they use 'debian:latest' which is currently 'buster' but the old container was built when 'stretch' was latest so use that instead

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# ENV PATH /usr/conda/bin:$PATH
 ENV PATH=/usr/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

RUN apt-get update --fix-missing && \
apt-get install -y wget bzip2 ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 git mercurial subversion procps make && \
apt-get clean

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda2-4.7.12-Linux-x86_64.sh -O ~/miniconda.sh && \
/bin/bash ~/miniconda.sh -b -p /usr/conda && \
rm ~/miniconda.sh && \
ln -s /usr/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
echo ". /usr/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
echo "conda activate base" >> ~/.bashrc && \
find /usr/conda/ -follow -type f -name '*.a' -delete && \
find /usr/conda/ -follow -type f -name '*.js.map' -delete && \
/usr/conda/bin/conda clean -afy


# add Conpair
RUN git clone https://github.com/mskcc/Conpair.git && \
    cd Conpair && \
    pip install -r requirements.txt
ENV PATH=/Conpair/:/Conpair/scripts/:$PATH
