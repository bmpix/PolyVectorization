FROM ubuntu:16.04

RUN apt update && \
    apt install -y libopencv-dev libomp-dev libboost-all-dev wget git software-properties-common

WORKDIR /tmp/cmake
RUN wget https://github.com/Kitware/CMake/releases/download/v3.30.3/cmake-3.30.3-linux-x86_64.sh && \
    sh cmake-3.30.3-linux-x86_64.sh --skip-license && rm -rf cmake-3.30.3-linux-x86_64.sh && \
    cp -r * /usr/ && rm -rf /tmp/cmake

WORKDIR /tmp
RUN git clone --branch 3.3 https://gitlab.com/libeigen/eigen.git && \
    cd eigen && mkdir build && cd build && cmake .. && make install && \
    cd /usr/local/include && \
    ln -sf eigen3/Eigen Eigen && \
    ln -sf eigen3/unsupported unsupported && \
    rm -rf /tmp/eigen

RUN add-apt-repository -y ppa:beineri/opt-qt571-xenial && \
    apt update && apt install -y qt57-meta-minimal mesa-common-dev libglu1-mesa-dev && \
    apt clean && apt autoremove

WORKDIR /app

COPY cmake /app/cmake
COPY src /app/src
COPY CMakeLists.txt /app/

WORKDIR /app/build_cli
RUN cmake .. -D CMAKE_BUILD_TYPE=Release && \
    make -j8 && mv polyvector_thing /usr/bin/polyvector_thing_cli && \
    rm -rf /app/build_cli

WORKDIR /app/build_cli_qt
RUN cmake .. -D CMAKE_BUILD_TYPE=Release -D WITH_QT=1 && \
    make -j8 && mv polyvector_thing /usr/bin/polyvector_thing_cli_qt && \
    rm -rf /app/build_cli_qt

WORKDIR /app/build_gui
RUN cmake .. -D CMAKE_BUILD_TYPE=Release -D WITH_GUI=1 -D WITH_QT=1 && \
    make -j8 && mv polyvector_thing /usr/bin/polyvector_thing && \
    rm -rf /app/build_gui
        
WORKDIR /app/
RUN rm -rf *