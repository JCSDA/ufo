FROM  jcsda/docker:latest

#ENV OOPS_CMAKE_OPTIONS="-DCMAKE_BUILD_TYPE=DEBUG -DENABLE_GPROF=ON"
#ENV COVERAGE=ON
RUN touch /env.txt
RUN printenv > /env.txt
  
RUN mkdir -p /var/run/sshd \
    && ssh-keygen -A \
    && sed -i 's/#PermitRootLogin yes/PermitRootLogin yes/g' /etc/ssh/sshd_config \
    && sed -i 's/#RSAAuthentication yes/RSAAuthentication yes/g' /etc/ssh/sshd_config \
    && sed -i 's/#PubkeyAuthentication yes/PubkeyAuthentication yes/g' /etc/ssh/sshd_config

RUN groupadd jcsda -g 9999
RUN adduser jcsdauser

RUN mkdir -p  /jcsda /build_repo \
 &&  chown -R jcsdauser:jcsda /jcsda /usr/local /build_repo \
 &&  chmod 6755 /jcsda /build_repo /usr/local

RUN mkdir /jcsda/.ssh ; echo "StrictHostKeyChecking no" > /jcsda/.ssh/config
COPY default-mca-params.conf /jcsda/.openmpi/mca-params.conf
RUN mkdir -p /jcsda/.openmpi
RUN chown -R jcsdauser:jcsda /jcsda/
RUN chown -R jcsdauser:jcsda /build_repo/

USER jcsdauser
WORKDIR /jcsda

RUN ssh-keygen -f /jcsda/.ssh/id_rsa -t rsa -N '' \
    && chmod 600 /jcsda/.ssh/config \
    && chmod 700 /jcsda/.ssh \
    && cp /jcsda/.ssh/id_rsa.pub /jcsda/.ssh/authorized_keys

VOLUME /jcsda
VOLUME /build_repo

CMD ["/bin/bash"]
