#!/usr/bin/python
import logging
import os
import sys

from dotenv import load_dotenv

load_dotenv()

APP_DIR = os.getenv("APP_DIR")

logging.basicConfig(stream=sys.stderr)
sys.path.insert(0, APP_DIR)

from viewer.app import dashboard as application  # noqa
