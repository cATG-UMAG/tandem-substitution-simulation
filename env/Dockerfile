FROM python:3.8-slim-buster
MAINTAINER Diego Alvarez S. (dialvarezs@gmail.com)
LABEL description="to use with https://github.com/CATG-UMAG/tandem-substitution-simulation"

RUN apt update && apt upgrade -y
RUN pip3 install -r https://raw.githubusercontent.com/CATG-UMAG/tandem-substitution-simulation/master/env/requirements.txt
