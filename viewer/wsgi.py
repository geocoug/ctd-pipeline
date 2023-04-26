#!/usr/bin/python
import logging
import os
import sys

APP_DIR = "/usr/local/webs/asv"
os.environ["APP_DIR"] = APP_DIR

logging.basicConfig(stream=sys.stderr)
sys.path.insert(0, APP_DIR)

from viewer.app import dashboard as application  # noqa
