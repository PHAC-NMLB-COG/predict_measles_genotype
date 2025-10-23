"""
Testing different utility aspects of the tool
"""

import os
from collections import defaultdict

import pytest

from src.predict_genotype.predict import check_external_dependencies, find_fastq_files


@pytest.fixture
def script_path():
    return os.path.dirname(os.path.abspath(__file__))


# Tests for the dependency check
def test_dependency_fail():
    deps = ["this_dep_doesn't_exist"]
    with pytest.raises(SystemExit):
        check_external_dependencies(deps)


def test_dependency_pass():
    deps = ["python"]
    check_external_dependencies(deps)


# Fastq file pair check
def test_find_fastq_files_pass(script_path):
    indir = os.path.join(script_path, "data")
    samples = find_fastq_files(indir)

    # Structure
    assert isinstance(samples, defaultdict)
    assert len(samples) == 3

    # Files found
    assert "data/sample2_R1_001.fq" in samples.get("sample2").get("R1")
    assert "data/test_R2.fastq" in samples.get("test").get("R2")
    assert "data/ont.fastq.gz" in samples.get("ont").get("single")

    # Files not found
    assert samples.get("test").get("single") is None
    assert samples.get("ont").get("R1") is None
    assert samples.get("DNE") is None


def test_find_fastq_files_no_dir():
    indir = "DNE"
    with pytest.raises(SystemExit):
        _ = find_fastq_files(indir)


def test_find_fastq_no_samples(tmp_path):
    d = tmp_path / "no_data"
    d.mkdir()

    with pytest.raises(SystemExit):
        _ = find_fastq_files(d)
