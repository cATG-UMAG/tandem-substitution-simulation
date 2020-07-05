FROM debian:10-slim
MAINTAINER Diego Alvarez S. (dialvarezs@gmail.com)
LABEL description="to use with https://github.com/CATG-UMAG/tandem-substitution-simulation"

RUN apt update && apt upgrade -y
RUN apt install -y python3 python3-pip wget
RUN wget https://raw.githubusercontent.com/CATG-UMAG/tandem-substitution-simulation/master/env/requirements.txt
RUN pip3 install -r requirements.txt
