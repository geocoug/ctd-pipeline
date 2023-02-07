# Base image
FROM python:3.10-slim

# Update outdated system packages, install additional dependencies
RUN apt-get update && \
    apt-get install -y \
    binutils \
    # GDAL dependencies
    libproj-dev \
    libgdal-dev \
    gdal-bin \
    # Udunits dependencies
    udunits-bin \
    libudunits2-dev && \
    apt-get update -y && \
    pip install --no-cache-dir --upgrade pip==23.0 && \
    rm -rf /var/lib/apt/lists/*

# Set env variables
ENV HOME=/usr/local/app
ENV VENV=/opt/.venv
ENV PATH="$VENV/bin:$PATH"
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1
ENV UDUNITS2_XML_PATH=/usr/share/xml/udunits/udunits2-common.xml

# Set the current working directory
WORKDIR $HOME

# Create a Python virtual environment
RUN python -m venv ${VENV}

# Copy Python requirements to the container
COPY ./requirements.txt .

# Install Python dependencies
RUN python -m pip install --no-cache-dir -r requirements.txt

# Copy codebase to the container
COPY . ${HOME}
