#!/usr/bin/python
import logging
import sys

APP_DIR = "/usr/local/bin/ctd-pipeline"

logging.basicConfig(stream=sys.stderr)
sys.path.insert(0, APP_DIR)

from viewer.app import dashboard as application  # noqa
