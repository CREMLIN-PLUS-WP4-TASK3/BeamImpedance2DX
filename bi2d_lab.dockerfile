from dolfinx/lab

RUN apt-get update \
	&& DEBIAN_FRONTEND=noninteractive \
	apt install -y python3-meshio \
	&& apt autoremove \
	&& apt autoclean \
	&& rm -rf /var/lib/apt/lists/*

RUN pip3 install git+https://gitlab.com/matsievskiysv/bi2d.git@dev \
    && wget https://gitlab.com/matsievskiysv/bi2d/-/raw/dev/tools/convert_msh.py -O /usr/bin/convert_msh.py \
    && chmod 555 /usr/bin/convert_msh.py
