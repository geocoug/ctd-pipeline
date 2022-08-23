FROM python:3.10

RUN apt-get update && \
    apt-get install -y binutils libproj-dev libgdal-dev gdal-bin

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

WORKDIR /app

ENV PYTHON_VENV=/opt/.venv

RUN python -m venv $PYTHON_VENV

ENV PATH="$PYTHON_VENV/bin:$PATH"

RUN python -m pip install --upgrade pip

COPY requirements.txt .

RUN python -m pip install --no-cache-dir -r requirements.txt

COPY config.json .

COPY *.py .
