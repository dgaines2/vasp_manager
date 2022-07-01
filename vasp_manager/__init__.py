"""
vasp_manager

A python package to run and analyze vasp calculations
"""


__package__ = "vasp_manager"
__version__ = "0.1.0"
__author__ = "Dale Gaines II"

import logging

from vasp_manager.main_manager import VaspManager

logging.getLogger(__name__).addHandler(logging.NullHandler())
