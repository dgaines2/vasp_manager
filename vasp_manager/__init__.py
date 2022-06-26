"""
vasp_manager

A python package to run and analyze vasp calculations
"""


__package__ = "vasp_manager"
__version__ = "0.1.0"
__author__ = "Dale Gaines II"

import logging

from .calculation_manager import (
    manage_calculation,
    manage_calculations,
    setup_bulkmod,
    setup_coarse_relax,
    setup_elastic,
    setup_relax,
)

logging.getLogger(__name__).addHandler(logging.NullHandler())
