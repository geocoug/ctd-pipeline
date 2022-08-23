# Base image
FROM python:3.10

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
    libudunits2-dev

# Environment variables
ENV UDUNITS2_XML_PATH=/usr/share/xml/udunits/udunits2-common.xml
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# App working directory
WORKDIR /app

# Python virtual environment path
ENV PYTHON_VENV=/opt/.venv

# Create Python virtual environment
RUN python -m venv $PYTHON_VENV

# Default Python path
ENV PATH="$PYTHON_VENV/bin:$PATH"

# Upgrade pip
RUN python -m pip install --upgrade pip

# Copy Python requirements to the container
COPY requirements.txt .

# Install Python dependencies
RUN python -m pip install --no-cache-dir -r requirements.txt

# Copy codebase to the container
COPY . .
