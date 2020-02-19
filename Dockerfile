FROM ubuntu:18.04

ARG   DEBIAN_FRONTEND=noninteractive
ARG   BUILD_DIR=/opt/fastaq

ARG   SAMTOOLS_VERSION=0.1.19
ENV   SAMTOOLS_URL=https://github.com/samtools/samtools/archive/${SAMTOOLS_VERSION}.tar.gz
ENV   SAMTOOLS_DIR=/usr/local/samtools-${SAMTOOLS_VERSION}

# Install the dependancies
RUN   apt-get update && \
      apt-get install --yes apt-utils && \
      apt-get --yes upgrade && \
      apt-get install --yes curl python3 python3-pip

RUN   apt-get install -y locales && \
      sed -i -e 's/# \(en_GB\.UTF-8 .*\)/\1/' /etc/locale.gen && \
      touch /usr/share/locale/locale.alias && \
      locale-gen
ENV   LANG     en_GB.UTF-8
ENV   LANGUAGE en_GB:en
ENV   LC_ALL   en_GB.UTF-8

# samtools build
# nore specific version requirement
RUN   apt-get install --yes libncurses5-dev libz-dev
RUN   mkdir -p ${SAMTOOLS_DIR} && \
      cd ${SAMTOOLS_DIR} && \
      curl -fsSL ${SAMTOOLS_URL} | tar xzf - --strip-components=1 && \
      echo -n "Building samtools... " && \
      make > make.log 2>&1; \
      EXIT=$?; \
      if [ $? -ne ${EXIT} ]; then echo "FAILED:" && cat make.log && exit ${EXIT}; else echo "OK."; fi
ENV   PATH ${SAMTOOLS_DIR}:${SAMTOOLS_DIR}/misc:${PATH}

# check samtools is an executable in the path
RUN   samtools 2> /tmp/samtools.out; \
      EXIT=$?; \
      if [ 1 -ne ${EXIT} ]; then \
         cat /tmp/samtools.out && exit ${EXIT}; \
      fi;

# Fastaq build
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
