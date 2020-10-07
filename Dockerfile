FROM continuumio/miniconda2:4.7.12
RUN git clone https://github.com/mskcc/Conpair.git && \
    cd Conpair && \
    pip install -r requirements.txt
