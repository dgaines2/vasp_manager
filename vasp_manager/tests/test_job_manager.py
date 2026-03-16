# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import json

from vasp_manager.job_manager import JobManager


def _make_job_manager(tmp_path, computer="cluster"):
    """
    Helper: write a minimal computing_config.json two levels up from calc_dir
    and return a JobManager pointed at calc_dir.

    Directory layout:
        tmp_path/
            computing_config.json
            material/
                rlx-coarse/      <- calc_dir
    """
    material_dir = tmp_path / "material"
    calc_dir = material_dir / "rlx-coarse"
    calc_dir.mkdir(parents=True)

    # JobManager reads config_dir = calc_dir.parents[1] = tmp_path
    computing_config = {
        "computer": computer,
        computer: {
            "user_id": "testuser",
            "potcar_dir": "/fake/potcars",
            "queuetype": "regular",
            "allocation": "fake_alloc",
            "constraint": "cpu",
            "vasp_module": "vasp/6.0",
            "ncore_per_node": 128,
            "ncore": 16,
        },
    }
    config_path = tmp_path / "computing_config.json"
    config_path.write_text(json.dumps(computing_config))

    # Create a dummy vasp.q so submit_job can find it
    vaspq_path = calc_dir / "vasp.q"
    vaspq_path.write_text("#!/bin/bash\n#SBATCH -n 1\nsrun vasp_std\n")

    return JobManager(calc_dir=calc_dir)


# ---------------------------------------------------------------------------
# Personal-computer bypass (no subprocess calls made)
# ---------------------------------------------------------------------------


def test_submit_job_personal_computer_returns_true(tmp_path):
    """
    submit_job on a personal computer skips sbatch and returns True
    """
    jm = _make_job_manager(tmp_path, computer="personal")
    result = jm.submit_job()
    assert result is True


def test_check_job_complete_personal_computer_returns_true(tmp_path):
    """
    _check_job_complete on a personal computer skips squeue and returns True
    """
    jm = _make_job_manager(tmp_path, computer="personal")
    # Manually write a jobid so _check_job_complete is reachable
    jobid_path = jm.calc_dir / "jobid"
    jobid_path.write_text("99999\n")
    result = jm._check_job_complete()
    assert result is True


# ---------------------------------------------------------------------------
# submit_job: no vasp.q file → returns False (signals need for resubmission)
# ---------------------------------------------------------------------------


def test_submit_job_missing_vaspq_returns_false(tmp_path):
    """
    submit_job returns False when vasp.q is absent
    """
    jm = _make_job_manager(tmp_path, computer="cluster")
    vaspq_path = jm.calc_dir / "vasp.q"
    vaspq_path.unlink()
    result = jm.submit_job()
    assert result is False


# ---------------------------------------------------------------------------
# submit_job: job already exists → returns True without sbatch
# ---------------------------------------------------------------------------


def test_submit_job_already_exists_returns_true(tmp_path):
    """
    submit_job returns True immediately if a jobid file already exists
    """
    jm = _make_job_manager(tmp_path, computer="cluster")
    jobid_path = jm.calc_dir / "jobid"
    jobid_path.write_text("12345\n")
    # If sbatch were called here it would fail (no SLURM); the test would error
    result = jm.submit_job()
    assert result is True


# ---------------------------------------------------------------------------
# submit_job: monkeypatched sbatch → tests jobid parsing
# ---------------------------------------------------------------------------


def test_submit_job_parses_jobid(tmp_path, monkeypatch):
    """
    submit_job calls sbatch, parses the jobid from stdout, and writes it
    to the jobid file
    """
    jm = _make_job_manager(tmp_path, computer="cluster")

    monkeypatch.setattr(
        "vasp_manager.job_manager.subprocess.check_output",
        lambda cmd, **kw: b"Submitted batch job 42001",
    )

    result = jm.submit_job()
    assert result is True
    assert jm.jobid == 42001
    jobid_file = jm.calc_dir / "jobid"
    assert jobid_file.exists()
    assert jobid_file.read_text().strip() == "42001"


# ---------------------------------------------------------------------------
# _check_job_complete: monkeypatched squeue
# ---------------------------------------------------------------------------


def test_check_job_complete_job_still_running(tmp_path, monkeypatch):
    """
    _check_job_complete returns False when the jobid is still in the squeue output
    """
    jm = _make_job_manager(tmp_path, computer="cluster")
    jobid_path = jm.calc_dir / "jobid"
    jobid_path.write_text("55555\n")

    # squeue output that includes our jobid
    squeue_output = (
        b"JOBID PARTITION  NAME  USER ST  TIME NODES\n"
        b" 55555  regular  rcNaCl  testuser  R 0:10 1\n"
    )
    monkeypatch.setattr(
        "vasp_manager.job_manager.subprocess.check_output",
        lambda cmd, **kw: squeue_output,
    )

    result = jm._check_job_complete()
    assert result is False


def test_check_job_complete_job_finished(tmp_path, monkeypatch):
    """
    _check_job_complete returns True when the jobid is absent from squeue output
    (job has completed or been cancelled)
    """
    jm = _make_job_manager(tmp_path, computer="cluster")
    jobid_path = jm.calc_dir / "jobid"
    jobid_path.write_text("55555\n")

    # squeue output that does NOT include our jobid
    squeue_output = (
        b"JOBID PARTITION  NAME  USER ST  TIME NODES\n"
        b" 99999  regular  other  testuser  R 0:05 1\n"
    )
    monkeypatch.setattr(
        "vasp_manager.job_manager.subprocess.check_output",
        lambda cmd, **kw: squeue_output,
    )

    result = jm._check_job_complete()
    assert result is True


# ---------------------------------------------------------------------------
# job_exists property
# ---------------------------------------------------------------------------


def test_job_exists_no_file(tmp_path):
    """
    job_exists returns False when no jobid file is present
    """
    jm = _make_job_manager(tmp_path)
    assert jm.job_exists is False


def test_job_exists_empty_file(tmp_path):
    """
    job_exists returns False when the jobid file is empty
    """
    jm = _make_job_manager(tmp_path)
    jobid_path = jm.calc_dir / "jobid"
    jobid_path.write_text("")
    assert jm.job_exists is False


def test_job_exists_with_jobid(tmp_path):
    """
    job_exists returns True and sets jobid when a valid jobid file is present
    """
    jm = _make_job_manager(tmp_path)
    jobid_path = jm.calc_dir / "jobid"
    jobid_path.write_text("12345\n")
    assert jm.job_exists is True
    assert jm.jobid == 12345
