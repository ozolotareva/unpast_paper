FROM conda/miniconda3
WORKDIR /usr/src/desmond2/

ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

RUN apt-get update && apt-get upgrade -y ; apt-get -y install ssh git build-essential

RUN mkdir ~/.ssh/
RUN ssh-keyscan github.com >> ~/.ssh/known_hosts
#RUN apt-get update && apt-get install -y build-essential
#RUN apt-get install -y --no-install-recommends libxt6

#RUN conda update conda
#RUN conda install python=3.9
RUN conda install -c conda-forge -y r-base=4.2.1

RUN git clone https://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@github.com/ozolotareva/DESMOND2.git
RUN R -e "chooseCRANmirror(ind=1) ; install.packages('BiocManager') ;  BiocManager::install('WGCNA')"
WORKDIR "DESMOND2"

#RUN apt install -y python3.9
RUN python setup.py install

WORKDIR "evaluation/clustering"
COPY . ./

RUN pip install -r requirements.txt
RUN mkdir -p "/local/DESMOND2_data_simulated/"
ENTRYPOINT ["tail","-f","/dev/null"]
