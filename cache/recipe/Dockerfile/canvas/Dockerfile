FROM ubuntu:22.04
WORKDIR /opt

ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    echo $TZ > /etc/timezone
RUN apt-get -y update && \
    apt-get -y install wget perl make gcc libicu70 libcurl4-openssl-dev
RUN mkdir -p /opt/dotnet/2.1.30 && cd /opt/dotnet/2.1.30 && \
    wget https://download.visualstudio.microsoft.com/download/pr/d6040f80-8343-4771-9c02-dbc9a35ac88a/68e74e6e46cf36fa1a50f68af6831d6d/aspnetcore-runtime-2.1.30-linux-x64.tar.gz && \
    tar zxvf aspnetcore-runtime-2.1.30-linux-x64.tar.gz && rm aspnetcore-runtime-2.1.30-linux-x64.tar.gz && mkdir lib && cd lib && wget https://www.openssl.org/source/openssl-1.1.1c.tar.gz && \
    tar zxvf openssl-1.1.1c.tar.gz && rm openssl-1.1.1c.tar.gz && cd openssl-1.1.1c && ./config --prefix=/opt/dotnet/lib/openssl --openssldir=/opt2/dotnet/lib/openssl && make && make install
ENV PATH="${PATH}:/opt/dotnet/2.1.30"
RUN mkdir -p /opt/canvas && cd /opt/canvas && \
    wget https://github.com/Illumina/canvas/releases/download/1.40.0.1613%2Bmaster/Canvas-1.40.0.1613.master_x64.tar.gz && \
    tar zxvf Canvas-1.40.0.1613.master_x64.tar.gz && rm -rf 1.40.0.1613+master Canvas-1.40.0.1613.master_x64.tar.gz && \
    echo '#!/bin/bash' > /opt/canvas/Canvas.sh && \
    echo 'export LD_LIBRARY_PATH=/opt/dotnet/lib/openssl/lib' >> /opt/canvas/Canvas.sh && \
    echo '/opt/dotnet/2.1.30/dotnet /opt/canvas/Canvas-1.40.0.1613+master_x64/Canvas.dll "$@"' >> /opt/canvas/Canvas.sh && \
    chmod 755 /opt/canvas/Canvas.sh /opt/canvas/Canvas-1.40.0.1613+master_x64/* /opt/canvas/Canvas-1.40.0.1613+master_x64/*/* 
RUN find /opt/canvas |grep bedGraph|while read n; do   wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig -O $n; done
ENV PATH="${PATH}:/opt/canvas"
CMD Canvas.sh
