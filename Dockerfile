FROM friend1ws/annot_utils 
MAINTAINER Yuichi Shiraishi <friend1ws@gmail.com> 

WORKDIR /usr/local/bin

RUN pip install --upgrade pip
RUN pip install pysam

RUN wget https://github.com/Genomon-Project/fusionfusion/archive/v0.3.0.tar.gz && \
    tar xzvf v0.3.0.tar.gz && \
    cd fusionfusion-0.3.0 && \
    python setup.py build install

RUN wget https://github.com/friend1ws/fusion_utils/archive/v0.2.0.tar.gz && \
    tar xzvf v0.2.0.tar.gz && \
    cd fusion_utils-0.2.0 && \
    python setup.py build install

WORKDIR /data

ENTRYPOINT ["/usr/local/bin/fusion_utils"]
CMD ["--help"]

