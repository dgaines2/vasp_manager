# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import logging
from pathlib import Path

import numpy as np
from pymatgen.analysis.eos import BirchMurnaghan
from pymatgen.core import Structure
from pymatgen.io.vasp import Vasprun

logger = logging.getLogger(__name__)


class BulkmodAnalyzer:
    def __init__(self, calc_dir, rounding_precision=3):
        """
        Args:
            calc_dir (str): path to bulkmod calculation directory
            rounding_precision (int): precision to round calculated quantities
        """
        self.calc_dir = Path(calc_dir)
        self.rounding_precision = rounding_precision
        self._results = None

    @property
    def calc_dir(self):
        return self._calc_dir

    @calc_dir.setter
    def calc_dir(self, value):
        if not value.exists():
            raise ValueError(f"Could not set calc_dir to {value} as it does not exist")
        self._calc_dir = value

    @property
    def rounding_precision(self):
        return self._rounding_precision

    @rounding_precision.setter
    def rounding_precision(self, value):
        if not isinstance(value, int):
            raise ValueError(
                f"Could not set rounding_precision to {value}, rounding precision must"
                " be an int"
            )
        self._rounding_precision = value

    @staticmethod
    def analyze_bulkmod(calc_dir, rounding_precision):
        """
        Fit an EOS to calculate the bulk modulus from a finished bulkmod calculation
        """
        strain_dirs = [path for path in calc_dir.glob("strain*") if path.is_dir()]
        strain_dirs = sorted(strain_dirs, key=lambda d: int(d.name.split("_")[-1]))
        volumes = []
        final_energies = []
        for i, strain_dir in enumerate(strain_dirs):
            poscar_path = strain_dir / "POSCAR"
            # search for vasprun.xml or vasprun.xml.gz
            vasprun_glob = list(strain_dir.glob("vasprun.xml*"))
            if len(vasprun_glob) == 0:
                raise Exception(f"No vasprun.xml available at {strain_dir}")
            vasprun_path = vasprun_glob[0]
            volume = Structure.from_file(poscar_path).volume
            vasprun = Vasprun(
                filename=vasprun_path,
                parse_dos=False,
                parse_eigen=False,
                parse_potcar_file=False,
            )
            final_energy = vasprun.final_energy
            volumes.append(volume)
            final_energies.append(final_energy)
        logger.debug(f"Volumes:\n\t{volumes}")
        logger.debug(f"Final Energies:\n\t{final_energies}")
        eos_analyzer = BirchMurnaghan(volumes, final_energies)
        eos_analyzer.fit()
        bulk_modulus = np.round(eos_analyzer.b0_GPa, rounding_precision)
        logger.debug(f"BULK MODULUS: {bulk_modulus}")
        b_dict = {"B": bulk_modulus}
        return b_dict

    @property
    def results(self):
        if self._results is None:
            self._results = self.analyze_bulkmod(self.calc_dir, self.rounding_precision)
        return self._results
