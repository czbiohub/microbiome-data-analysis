FROM biocontainers/biocontainers:debian-stretch-backports

USER root:root

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get clean && \
    apt-get update && \
    apt-get install -y grinder awscli pigz vim htop && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/*

COPY run_grinder.sh ./
RUN chmod +rx run_grinder.sh

LABEL    software="grinder" \
    maintainer="Sunit Jain" \
    container="grinder" \ 
    about.summary="Versatile omics shotgun and amplicon sequencing read simulator" \ 
    about.home="http://sourceforge.net/projects/biogrinder/" \ 
    software.version="0.5.4-1-deb" \ 
    version="1" \ 
    about.copyright="2009-2011, Florent Angly <florent.angly@gmail.com>" \ 
    about.license="GPL-3+" \ 
    about.license_file="/usr/share/doc/grinder/copyright" \ 
    extra.binaries="/usr/bin/average_genome_size,/usr/bin/change_paired_read_orientation,/usr/bin/grinder" \ 
    container.category="read_simulator"
