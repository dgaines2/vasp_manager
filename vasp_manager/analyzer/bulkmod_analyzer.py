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
    def __init__(self, calc_path=None, rounding_precision=3):
        """
        Args:
            calc_path (str): path to bulkmod calculation folder
            rounding_precision (int): precision to round calculated quantities
        """
        self._calc_path = Path(calc_path) if calc_path else calc_path
        self._rounding_precision = rounding_precision
        self._results = None

    @property
    def calc_path(self):
        if self._calc_path is not None:
            self.calc_path = self._calc_path
        return self._calc_path

    @calc_path.setter
    def calc_path(self, value):
        # make sure cij wasn't already defined
        if not value.exists():
            raise ValueError(f"Could not set calc_path to {value} as it does not exist")
        self._calc_path = value

    @property
    def rounding_precision(self):
        if self._rounding_precision is not None:
            self.rounding_precision = self._rounding_precision
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
    def analyze_bulkmod(calc_path, rounding_precision):
        """
        Fit an EOS to calculate the bulk modulus from a finished bulkmod calculation
        """
        strain_paths = [path for path in calc_path.glob("strain*") if path.is_dir()]
        strain_paths = sorted(strain_paths, key=lambda d: int(d.name.split("_")[-1]))
        volumes = []
        final_energies = []
        for i, strain_path in enumerate(strain_paths):
            poscar_path = strain_path / "POSCAR"
            # search for vasprun.xml or vasprun.xml.gz
            vasprun_glob = list(strain_path.glob("vasprun.xml*"))
            if len(vasprun_glob) == 0:
                raise Exception(f"No vasprun.xml available at {strain_path}")
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
        logger.debug("Volumes")
        logger.debug(f"{volumes}")
        logger.debug("Final Energies")
        logger.debug(f"{final_energies}")
        eos_analyzer = BirchMurnaghan(volumes, final_energies)
        eos_analyzer.fit()
        bulk_modulus = np.round(eos_analyzer.b0_GPa, rounding_precision)
        logger.info(f"BULK MODULUS: {bulk_modulus}")
        b_dict = {"B": bulk_modulus}
        return b_dict

    @property
    def results(self):
        if self._results is None:
            self._results = self.analyze_bulkmod(self.calc_path, self.rounding_precision)
        return self._results
