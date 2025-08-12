# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

from __future__ import annotations

from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray
from pymatgen.core import Structure

from vasp_manager.utils import pgrep

if TYPE_CHECKING:
    from vasp_manager.types import Floating, WorkingDirectory


class ElasticAnalyzer:
    def __init__(
        self,
        cij: NDArray,
        structure: Structure,
        rounding_precision: int = 3,
    ) -> None:
        """
        Args:
            cij: stiffness tensor in Voigt notation in units of GPa
            structure: the associated pymatgen structure
            rounding_precision: precision to round calculated quantities

        Note:
            Be careful with cij! VASP does not output stiffness tensors in
            typical Voigt notation. See change_elastic_constants_from_vasp
            method.
        """
        self.cij = cij
        self.structure = structure
        self.rounding_precision = rounding_precision
        self._results: dict

    @property
    def cij(self) -> NDArray:
        return self._cij

    @cij.setter
    def cij(self, values: NDArray) -> None:
        values = np.asarray(values)
        if values.shape != (6, 6):
            raise ValueError(
                "cij must be a 6x6 array," f"found shape={values.shape} instead"
            )
        if values.dtype != float:
            raise TypeError(
                "cij must be a 6x6 array of floats," f" found type={values.dtype} instead"
            )
        self._cij = values

    @property
    def structure(self) -> Structure:
        return self._structure

    @structure.setter
    def structure(self, value: Structure) -> None:
        if not isinstance(value, Structure):
            raise TypeError(
                "Structure must be pymatgen Structure object,"
                f" found type={type(value)} instead"
            )
        self._structure = value

    @property
    def density(self) -> float:
        return self.structure.density

    @property
    def rounding_precision(self) -> int:
        return self._rounding_precision

    @rounding_precision.setter
    def rounding_precision(self, value: int) -> None:
        if not isinstance(value, int):
            raise TypeError(
                f"Could not set rounding_precision to {value}. rounding_precision must"
                " be an int, found type={type(value)} instead)"
            )
        self._rounding_precision = value

    @staticmethod
    def change_elastic_constants_from_vasp(vasp_elastic_tensor: NDArray) -> NDArray:
        """
        Args:
            vasp_elastic_tensor: elastic tensor from VASP OUTCAR

        Returns:
            elastic_tensor: reordered to match Voigt notation

        Note:
            VASP ordering of elastic constants is not the same as those expected by
            my equations below

            We should expect from Voigt notation:
            1: 11 or xx
            2: 22 or yy
            3: 33 or zz
            4: 23 or yz
            5: 13 or xz
            6: 12 or xy
            but VASP differs as it presents 1, 2, 3, 6 (xy), 4 (yz), 5 (xz)
            This function performs swapping to match expectations
        """
        elastic_tensor = np.asarray(vasp_elastic_tensor).copy()
        for j in range(6):
            elastic_tensor[3, j], elastic_tensor[4, j], elastic_tensor[5, j] = (
                elastic_tensor[4, j],
                elastic_tensor[5, j],
                elastic_tensor[3, j],
            )
        for i in range(6):
            elastic_tensor[i, 3], elastic_tensor[i, 4], elastic_tensor[i, 5] = (
                elastic_tensor[i, 4],
                elastic_tensor[i, 5],
                elastic_tensor[i, 3],
            )
        return elastic_tensor

    @staticmethod
    def get_compliance_tensor(cij: NDArray) -> NDArray:
        """
        Args:
            cij: stiffness tensor

        Returns:
            sij: compliance tensor
        """
        sij = np.linalg.inv(cij)
        return sij

    @staticmethod
    def get_VRH_average(mod1: Floating, mod2: Floating) -> Floating:
        """
        Args:
            mod1: B or G Voigt
            mod2: B or G Reuss

        Returns:
            VRH_average
        """
        return (mod1 + mod2) / 2.0

    @staticmethod
    def get_B_Reuss(sij: NDArray) -> float:
        """
        Args:
            sij: compliance tensor

        Returns:
            B_Reuss: Reuss bulk modulus
        """
        B_Reuss = 1 / (
            (sij[0, 0] + sij[1, 1] + sij[2, 2]) + 2 * (sij[0, 1] + sij[1, 2] + sij[2, 0])
        )
        return B_Reuss

    @staticmethod
    def get_B_Voigt(cij: NDArray) -> float:
        """
        Args:
            cij: stiffness tensor

        Returns:
            B_Voigt: Voigt bulk modulus
        """
        B_Voigt = (
            (cij[0, 0] + cij[1, 1] + cij[2, 2]) + 2 * (cij[0, 1] + cij[1, 2] + cij[2, 0])
        ) / 9
        return B_Voigt

    @staticmethod
    def get_G_Reuss(sij: NDArray) -> float:
        """
        Args:
            sij: compliance tensor

        Returns:
            G_Reuss: Reuss shear modulus
        """
        G_Reuss = 15 / (
            4 * (sij[0, 0] + sij[1, 1] + sij[2, 2])
            - 4 * (sij[0, 1] + sij[1, 2] + sij[2, 0])
            + 3 * (sij[3, 3] + sij[4, 4] + sij[5, 5])
        )
        return G_Reuss

    @staticmethod
    def get_G_Voigt(cij: NDArray) -> float:
        """
        Args:
            cij: stiffness tensor

        Returns:
            G_Voigt: Voigt shear modulus
        """
        G_Voigt = (
            (cij[0, 0] + cij[1, 1] + cij[2, 2])
            - (cij[0, 1] + cij[1, 2] + cij[2, 0])
            + 3 * (cij[3, 3] + cij[4, 4] + cij[5, 5])
        ) / 15
        return G_Voigt

    @staticmethod
    def get_vl(
        bulkmod: Floating,
        shearmod: Floating,
        density: Floating,
    ) -> Floating:
        """
        Returns:
            vl: longitudinal speed of sound
        """
        vl = np.sqrt((bulkmod + 4 / 3 * shearmod) / density)
        return vl

    @staticmethod
    def get_vt(
        shearmod: Floating,
        density: Floating,
    ) -> Floating:
        """
        Returns:
            vl: transverse speed of sound
        """
        vt = np.sqrt(shearmod / density)
        return vt

    @staticmethod
    def get_vs(
        bulkmod: Floating,
        shearmod: Floating,
        density: Floating,
    ) -> Floating:
        """
        Returns:
            vs: speed of sound
        """
        vl = np.sqrt((bulkmod + 4 / 3 * shearmod) / density)
        vt = np.sqrt(shearmod / density)
        vs = (1 / 3 * (1 / vl**3 + 2 / vt**3)) ** (-1 / 3)
        return vs

    @staticmethod
    def check_elastically_unstable(cij: NDArray) -> bool:
        """
        Returns:
            True if compound is elastically unstable
        """
        eigenvalues, eigenvectors = np.linalg.eig(cij)
        born_critera_satisfied = np.all(eigenvalues > 0)
        elastically_unstable = not born_critera_satisfied
        return elastically_unstable

    @cached_property
    def sij(self) -> NDArray:
        # DO NOT ROUND
        return self.get_compliance_tensor(self.cij)

    @cached_property
    def b_reuss(self) -> Floating:
        return np.round(self.get_B_Reuss(self.sij), self.rounding_precision)

    @cached_property
    def b_voigt(self) -> Floating:
        return np.round(self.get_B_Voigt(self.cij), self.rounding_precision)

    @cached_property
    def b_vrh(self) -> Floating:
        return np.round(
            self.get_VRH_average(self.b_reuss, self.b_voigt), self.rounding_precision
        )

    @cached_property
    def g_reuss(self) -> Floating:
        return np.round(self.get_G_Reuss(self.sij), self.rounding_precision)

    @cached_property
    def g_voigt(self) -> Floating:
        return np.round(self.get_G_Voigt(self.cij), self.rounding_precision)

    @cached_property
    def g_vrh(self) -> Floating:
        return np.round(
            self.get_VRH_average(self.g_reuss, self.g_voigt), self.rounding_precision
        )

    @cached_property
    def elastically_unstable(self) -> bool:
        return self.check_elastically_unstable(self.cij)

    @cached_property
    def vt(self) -> Floating:
        return np.round(self.get_vt(self.g_vrh, self.density), self.rounding_precision)

    @cached_property
    def vl(self) -> Floating:
        return np.round(
            self.get_vl(self.b_vrh, self.g_vrh, self.density), self.rounding_precision
        )

    @cached_property
    def vs(self) -> Floating:
        return np.round(
            self.get_vs(self.b_vrh, self.g_vrh, self.density), self.rounding_precision
        )

    def _analyze_elastic(self, properties: list | None = None) -> dict:
        """
        Grabs important quantities from the elastic calculation results

        Args:
            properties (list): names of properties to calculate

        Returns:
            elastic_dict (dict): dict of extracted info from
                elastic constants
        """
        properties_map = {
            "B_Reuss": "b_reuss",
            "B_Voigt": "b_voigt",
            "B_VRH": "b_vrh",
            "G_Reuss": "g_reuss",
            "G_Voigt": "g_voigt",
            "G_VRH": "g_vrh",
            "unstable": "elastically_unstable",
            "elastic_tensor": "cij",
            "vl": "vl",
            "vt": "vt",
            "vs": "vs",
        }
        if properties is None:
            properties = list(properties_map.keys())

        elastic_dict = {}
        for property in properties:
            elastic_dict[property] = getattr(self, properties_map[property])
        return elastic_dict

    @property
    def results(self) -> dict:
        if getattr(self, "_results", None) is None:
            self._results = self._analyze_elastic()
        return self._results

    @classmethod
    def from_calc_dir(
        cls,
        calc_dir: WorkingDirectory,
        **kwargs,
    ) -> ElasticAnalyzer:
        """
        Method to construct ElasticAnalyzer from a calculation directory.
        The directory should contain a POSCAR and an OUTCAR with elastic constants

        Args:
            calc_dir:

        Returns:
            Instance of ElasticAnalyzer
        """
        calc_dir = Path(calc_dir)
        if not calc_dir.exists():
            raise ValueError(f"Could not set calc_dir to {calc_dir} as it does not exist")

        structure = Structure.from_file(calc_dir / "POSCAR")
        outcar_glob = list(calc_dir.glob("OUTCAR*"))
        if len(outcar_glob) == 0:
            raise FileNotFoundError(f"No OUTCAR available in {calc_dir}")
        outcar_filepath = outcar_glob[0]

        elastic_filepath = calc_dir / "elastic_constants.txt"
        if not elastic_filepath.exists():
            # Scrapes OUTCAR for elastic constants and write to elastic_constants.txt
            elastic_table = pgrep(
                outcar_filepath,
                str_to_grep="TOTAL ELASTIC MOD",
                stop_after_first_match=True,
                after=8,
                as_string=True,
            )
            assert isinstance(elastic_table, str)
            with open(elastic_filepath, "w+") as fw:
                fw.write(elastic_table)

        # Read the elastic constants and construct the stiffness tensor
        with open(elastic_filepath) as fr:
            raw_elastic_data = fr.readlines()
        try:
            # Skip first 3 rows as they are just header
            # Skip the first column as it contains xx, xy, etc
            elastic_data = [line.strip().split()[1:] for line in raw_elastic_data[3:]]
            # Divide by 10 to get GPa instead of kBar
            elastic_tensor = np.array(elastic_data, dtype=float) / 10.0
            # VASP only outputs the stiffness tensor to 4 places
            elastic_tensor = np.round(elastic_tensor, 4)
            elastic_tensor = ElasticAnalyzer.change_elastic_constants_from_vasp(
                elastic_tensor
            )
        except Exception as e:
            raise RuntimeError(
                "Could not construct elastic tensor from OUTCAR. Please verify the"
                " OUTCAR contains a full elastic tensor.\n"
                f"{e}"
            )
        return cls(cij=elastic_tensor, structure=structure, **kwargs)
