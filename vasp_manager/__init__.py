"""
vasp_manager

A python package to run and analyze VASP calculations
"""

__package__ = "vasp_manager"
__version__ = "1.1.4"
__author__ = "Dale Gaines II"

import logging

from vasp_manager.vasp_manager import VaspManager

logging.getLogger(__name__).addHandler(logging.NullHandler())
