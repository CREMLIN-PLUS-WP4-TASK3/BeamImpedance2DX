from dolfinx/lab

RUN apt-get update \
	&& DEBIAN_FRONTEND=noninteractive \
	apt-get install -y python3-meshio \
	&& apt-get autoremove \
	&& apt-get autoclean \
	&& rm -rf /var/lib/apt/lists/*

RUN pip3 install git+https://gitlab.com/matsievskiysv/bi2d.git@master \
    && wget https://gitlab.com/matsievskiysv/bi2d/-/raw/dev/tools/convert_msh.py -O /usr/bin/convert_msh.py \
    && chmod 555 /usr/bin/convert_msh.py
