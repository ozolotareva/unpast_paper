FROM conda/miniconda3
WORKDIR /usr/src/desmond2/

ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

RUN apt-get update && apt-get upgrade -y ; apt-get -y install ssh git build-essential wget unzip

RUN mkdir ~/.ssh/
RUN ssh-keyscan github.com >> ~/.ssh/known_hosts

RUN conda install -c conda-forge -y r-base=4.2.1


RUN R -e "chooseCRANmirror(ind=1) ; install.packages('BiocManager') ; BiocManager::install('mosbi')"

RUN git clone https://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@github.com/ozolotareva/DESMOND2.git

WORKDIR "DESMOND2"
COPY eval.py ./utils/
RUN python setup.py install

WORKDIR "evaluation/biclustering/simulated"
RUN pip install -r requirements.txt

RUN wget https://github.com/maqin2001/qubic2/archive/refs/heads/master.zip
RUN unzip master.zip
RUN rm master.zip
WORKDIR "qubic2-master"
RUN make
WORKDIR "../"

RUN apt install -y openjdk-8-jdk

COPY . ./

RUN mkdir -p "/local/DESMOND2_data_simulated/"

ENTRYPOINT ["tail","-f","/dev/null"]
