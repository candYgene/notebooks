FROM python:2-slim

# RUN apk update
RUN apt-get update
# RUN apk add git gcc python-dev musl-dev libzmq zeromq-dev libpng freetype-dev
RUN apt-get install -y git gcc

RUN pip install grlc jupyter ipywidgets matplotlib pandas
RUN jupyter nbextension enable --py widgetsnbextension

RUN mkdir /queries
ADD config.docker.ini /queries/config.ini
ADD gunicorn_config.py /queries/
ADD *.rq /queries/
ADD *.ipynb /queries/

WORKDIR /queries

CMD gunicorn -c gunicorn_config.py grlc.server:app & jupyter notebook --port=8888 --ip=0.0.0.0 --allow-root --NotebookApp.token=''
