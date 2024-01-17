# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import logging
from functools import cached_property

import numpy as np

from vasp_manager.calculation_manager.base import BaseCalculationManager
from vasp_manager.utils import get_pmg_structure_from_poscar, pgrep, ptail
from vasp_manager.vasp_input_creator import VaspInputCreator

logger = logging.getLogger(__name__)


class RlxCalculationManager(BaseCalculationManager):
    """
    Runs relaxation job workflow for a single material
    """

    def __init__(
        self,
        material_path,
        to_rerun,
        to_submit,
        primitive=True,
        ignore_personal_errors=True,
        from_coarse_relax=True,
        from_scratch=False,
        tail=5,
        max_reruns=3,
        magmom_per_atom_cutoff=0.0,
    ):
        """
        For material_path, to_rerun, to_submit, ignore_personal_errors, and from_scratch,
        see BaseCalculationManager

        Args:
            from_coarse_relax (bool): if True, use CONTCAR from coarse relax
            tail (int): number of last lines to log in debugging if job failed
        """
        self.from_coarse_relax = from_coarse_relax
        self.tail = tail
        self.max_reruns = max_reruns
        self.magmom_per_atom_cutoff = magmom_per_atom_cutoff
        super().__init__(
            material_path=material_path,
            to_rerun=to_rerun,
            to_submit=to_submit,
            primitive=primitive,
            ignore_personal_errors=ignore_personal_errors,
            from_scratch=from_scratch,
        )
        self._is_done = None
        self._results = None

    @cached_property
    def mode(self):
        return "rlx"

    @cached_property
    def poscar_source_path(self):
        if self.from_coarse_relax:
            poscar_source_path = self.material_path / "rlx-coarse" / "CONTCAR"
        else:
            poscar_source_path = self.material_path / "POSCAR"
        return poscar_source_path

    @cached_property
    def vasp_input_creator(self):
        return VaspInputCreator(
            self.calc_path,
            mode=self.mode,
            poscar_source_path=self.poscar_source_path,
            primitive=self.primitive,
            name=self.material_name,
        )

    def setup_calc(
        self,
        increase_nodes_by_factor=1,
        increase_walltime_by_factor=1,
        make_archive=False,
        use_spin=True,
    ):
        """
        Sets up a fine relaxation
        """
        self.vasp_input_creator.increase_nodes_by_factor = increase_nodes_by_factor
        self.vasp_input_creator.increase_walltime_by_factor = increase_walltime_by_factor
        self.vasp_input_creator.use_spin = use_spin

        if make_archive:
            self.vasp_input_creator.make_archive_and_repopulate()
        else:
            self.vasp_input_creator.create()

        if self.to_submit:
            job_submitted = self.submit_job()
            if not job_submitted:
                self.setup_calc()

    def check_calc(self):
        """
        Checks if calculation has finished and reached required accuracy

        Returns:
            relaxation_successful (bool): if True, relaxation completed successfully
        """
        if not self.job_complete:
            logger.info(f"{self.mode.upper()} not finished")
            return False

        stdout_path = self.calc_path / "stdout.txt"
        if not stdout_path.exists():
            # calculation never actually ran
            # shouldn't get here unless function was called with submit=False
            # or job was manually cancelled
            logger.info(f"{self.mode.upper()} Calculation: No stdout.txt available")
            if self.to_rerun:
                self._cancel_previous_job()
                self.setup_calc()
            return False

        vasp_errors = self._check_vasp_errors()
        if len(vasp_errors) > 0:
            all_errors_addressed = self._address_vasp_errors(vasp_errors)
            if all_errors_addressed:
                if self.to_rerun:
                    logger.info(f"Rerunning {self.calc_path}")
                    self.setup_calc(make_archive=True)
            else:
                msg = (
                    f"{self.mode.upper()} Calculation: ",
                    "Couldn't address all VASP Errors\n",
                    "\tRefusing to continue...\n",
                    f"\tVasp Errors: {vasp_errors}\n",
                )
                logger.error(msg)
                self.stop()
            return False

        previous_magmom_per_atom = self._parse_magmom_per_atom()
        if previous_magmom_per_atom is None:
            use_spin = False
        else:
            use_spin = np.abs(previous_magmom_per_atom) >= self.magmom_per_atom_cutoff

        tail_output = ptail(stdout_path, n_tail=self.tail, as_string=True)
        grep_output = pgrep(
            stdout_path, "reached required accuracy", stop_after_first_match=True
        )
        if len(grep_output) == 0:
            archive_dirs = list(self.calc_path.glob("archive*"))
            if len(archive_dirs) >= self.max_reruns - 1:
                msg = (
                    "Many archives exist, calculations may not be converging\n"
                    "\t Refusing to continue..."
                )
                logger.error(msg)
                return False

            logger.warning(f"{self.mode.upper()} FAILED")
            logger.debug(tail_output)
            if self.to_rerun:
                logger.info(f"Rerunning {self.calc_path}")
                # increase nodes as its likely the calculation failed
                self.setup_calc(
                    increase_walltime_by_factor=2, make_archive=True, use_spin=use_spin
                )
            return False

        # in the case that the calculation finishes but results in a spin value lower
        # than the cutoff, rerun it without spin
        if not use_spin and previous_magmom_per_atom is not None and self.to_rerun:
            logger.info(f"Rerunning {self.calc_path}")
            self.setup_calc(make_archive=True, use_spin=False)
            return False

        logger.info(f"{self.mode.upper()} Calculation: reached required accuracy")
        logger.debug(tail_output)
        return True

    def check_volume_difference(self):
        """
        Checks relaxation runs for volume difference

        if abs(volume difference) is >= 5%, reruns relaxation
        only checks for mode='rlx' as that's the structure for further calculation

        Returns:
            volume_converged (bool): if True, relaxation completed successfully
        """
        original_poscar_path = self.material_path / "POSCAR"
        poscar_path = self.calc_path / "POSCAR"
        contcar_path = self.calc_path / "CONTCAR"
        try:
            orig_structure, orig_spacegroup = get_pmg_structure_from_poscar(
                original_poscar_path, return_spacegroup=True
            )
            p_structure, p_spacegroup = get_pmg_structure_from_poscar(
                poscar_path, return_spacegroup=True
            )
            c_structure, c_spacegroup = get_pmg_structure_from_poscar(
                contcar_path, return_spacegroup=True
            )
        except Exception as e:
            logger.error(f"  RLX CONTCAR doesn't exist or is empty: {e}")
            return False

        if orig_spacegroup == c_spacegroup:
            logger.debug(
                f"  Spacegroups match orig-{orig_spacegroup} == rlx-{c_spacegroup}"
            )
        else:
            logger.warning(
                "   Warning: spacegroups do not match "
                + f"orig-{orig_spacegroup} != rlx-{c_spacegroup}"
            )

        volume_diff = (c_structure.volume - p_structure.volume) / p_structure.volume
        if np.abs(volume_diff) >= 0.05:
            logger.warning(f"  NEED TO RE-RELAX: dV = {volume_diff:.4f}")
            volume_converged = False
            previous_magmom_per_atom = self._parse_magmom_per_atom()
            if previous_magmom_per_atom is None:
                use_spin = False
            else:
                use_spin = np.abs(previous_magmom_per_atom) >= self.magmom_per_atom_cutoff
            if self.to_rerun:
                self.setup_calc(make_archive=True, use_spin=use_spin)
        else:
            logger.info("  RLX volume converged")
            logger.debug(f"  dV = {volume_diff:.4f}")
            volume_converged = True
        orig_volume_diff = (
            c_structure.volume - orig_structure.volume
        ) / orig_structure.volume
        self._results = {}
        self._results["initial_spacegroup"] = orig_spacegroup
        self._results["relaxed_spacegroup"] = c_spacegroup
        self._results["total_dV"] = np.round(orig_volume_diff, 4)
        return volume_converged

    @property
    def is_done(self):
        if self._is_done is None:
            self._is_done = False
            calc_done = self.check_calc()
            if calc_done:
                volume_converged = self.check_volume_difference()
                if volume_converged:
                    self._is_done = True
        return self._is_done

    @property
    def results(self):
        if not self.is_done:
            if self.stopped:
                return "STOPPED"
            else:
                return None
        return self._results
