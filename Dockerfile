#
# This container will allow you to run RAILS and Cobbler
#
FROM debian:testing

#
# Authorship
#
MAINTAINER rwarren@bcgsc.ca

#
# Update and Install dependencies
#
RUN apt-get update -qq && apt-get install -y bwa wget cpanminus 

#
# Download the software
#
RUN wget https://github.com/warrenlr/RAILS/raw/master/rails_v1-1.tar.gz && tar xvfz rails_v1-1.tar.gz && rm rails_v1-1.tar.gz 

#
# Set the default working directory
#
WORKDIR /RAILS_v1.1
