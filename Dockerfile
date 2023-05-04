FROM ubuntu:bionic
#i.e. unbuntu 18.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update 

RUN apt-get install gnupg -y \
  && apt-get clean
  
RUN apt-get install wget -y \
  && apt-get clean
 
#this is needed for the R install, but turns out it also installs many things, including python3 
RUN apt-get install --no-install-recommends software-properties-common dirmngr -y \
  && apt-get clean

#symbolic path to python3.6 from python
RUN ln -s /usr/bin/python3.6 /usr/bin/python

#this is just to update the cran links in the apt-get calls so that it points to the correct R version
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
  && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran35/"
  
#finally install R 3.6
RUN apt-get install r-base -y \
  && apt-get install r-base-dev -y \
  && apt-get clean

#now R packages
RUN R -e "install.packages(c('argparse', 'stringr', 'purrr', 'dplyr', 'multidplyr', 'tidyr', 'data.table', 'parallel', 'rcompanion'), rdependencies=TRUE, lib='/usr/lib/R/site-library')"

#install pandas
RUN apt-get -y install --no-install-recommends python3-pandas \
  && apt-get clean

#install java
RUN apt-get install default-jre -y --no-install-recommends \
  && apt-get clean

#copy the main branch of the git (alternatively could also git clone it with 'git clone https://github.com/DrGBL/snp2hla_redux --single-branch main')
#copy snp2hla_redux
COPY main main

#make your way there
WORKDIR /main/

ENTRYPOINT ["/bin/bash"]