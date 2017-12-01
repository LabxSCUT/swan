#NOTE: developer tested the commands by folloing actions:
#      docker pull charade/xlibbox:latest
#      docker images
#      docker run --memory=2g -i -t charade/xlibbox:latest /bin/bash
#NOTE: merging RUN as file become stable as every RUN creates a commit which has a limit
FROM charade/xlibbox:latest

MAINTAINER Charlie Xia <lxia.usc@gmail.com>

WORKDIR $HOME/
ENV PATH "PATH=$PATH:$HOME/bin"
RUN ulimit -s unlimited

### install swan ###
RUN cd $HOME/setup && git clone https://bitbucket.org/charade/swan.git
RUN cd $HOME/setup/swan && R CMD INSTALL .

### run test scripts ### 
RUN cd $HOME/setup/swan/test && wget https://s3-us-west-2.amazonaws.com/lixiabucket/swan_test.tgz -O swan_test.tgz && tar -zxvf swan_test.tgz
RUN export SWAN_BIN=$HOME/setup/swan/inst && cd $HOME/setup/swan/test && ./swan_test.sh

### now mount your server dir for data analysis
#docker run -d --name="foo" -v "/home/lixia:/home/lixia" ubuntu
