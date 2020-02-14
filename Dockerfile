FROM ubuntu:18.04

ARG   DEBIAN_FRONTEND=noninteractive
ARG   BUILD_DIR=/opt/fastaq

# Install the dependancies
RUN   apt-get update && \
      apt-get install --yes apt-utils && \
      apt-get --yes upgrade && \
      apt-get install --yes python3 python3-pip samtools

RUN   apt-get install -y locales && \
      sed -i -e 's/# \(en_GB\.UTF-8 .*\)/\1/' /etc/locale.gen && \
      touch /usr/share/locale/locale.alias && \
      locale-gen
ENV   LANG     en_GB.UTF-8
ENV   LANGUAGE en_GB:en
ENV   LC_ALL   en_GB.UTF-8

COPY  . ${BUILD_DIR}

RUN   cd ${BUILD_DIR} && \
      python3 setup.py test && \
      python3 setup.py install && \
      echo 'Installed '`which fastaq`

# check fastaq is an executable in the path
RUN   fastaq 2> /tmp/fastaq.out; \
      EXIT=$?; \
      if [ 1 -ne ${EXIT} ]; then \
         cat /tmp/fastaq.out && exit ${EXIT}; \
      fi;

CMD   echo "Usage:  docker run -v \`pwd\`:/var/data -it <IMAGE_NAME> bash" && \
      echo "" && \
      echo "This will place you in a shell with your current working directory accessible as /var/data." && \
      echo "You can then run commands like:" && \
      echo "   fastaq reverse_complement /var/data/in.fastq /var/data/out.fastq" && \
      echo "" && \
      fastaq | tail +3
