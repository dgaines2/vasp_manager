"""
vasp_manager

A python package to run and analyze vasp calculations
"""


__version__ = "1.0.0"
__author__ = "Dale Gaines"


from .calculation_manager import (
    analyze_bulkmod,
    check_elastic_calculation,
    check_relax,
    check_volume_difference,
    do_analyze_elastic,
    manage_calculations,
    setup_coarse_relax,
    setup_elastic,
    setup_relax,
    setup_simple_bulkmod,
)
from .elastic_analysis import (
    analyze_elastic,
    get_primitive_structure_from_poscar,
    make_elastic_constants,
)
