"""
Testing the analysis aspects of the tool

To run make sure to have the external dependencies installed
- Minimap2
- Samtools
"""

import os
from argparse import Namespace
from collections import defaultdict
from importlib.resources import files

import pytest

import src.predict_genotype.reference_data
from src.predict_genotype.predict import map_reads_to_ref, process_sample, summarize_genotype_hits


@pytest.fixture
def fastq_file():
    basepath = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(basepath, "data", "ont.fastq.gz")


@pytest.fixture
def ref_file():
    ref_dir = files(src.predict_genotype.reference_data)
    return ref_dir.joinpath("measles_N450_genotypes.mmi")


@pytest.fixture
def bam_file(ref_file, fastq_file, tmp_path):
    # Sorted_bam is the tmp output path
    sorted_bam = tmp_path / "test1-ont.sorted.bam"
    map_reads_to_ref(ref_file, fastq_file, sorted_bam)
    return sorted_bam


def test_genotype_summary(bam_file):
    hits, total = summarize_genotype_hits(bam_file)
    assert isinstance(hits, defaultdict)
    assert total == 125

    sorted_hits = sorted(hits.items(), key=lambda x: x[1], reverse=True)
    top_ref, top_count = sorted_hits[0]
    sec_ref, sec_count = sorted_hits[1]

    assert top_ref == "D8"
    assert top_count == 100
    assert sec_ref == "D4"
    assert sec_count == 25


@pytest.fixture
def simple_namespace(ref_file):
    return Namespace(reference=ref_file, threads=1, majority_threshold=60, min_reads=100)


@pytest.fixture
def test_dir(tmp_path):
    test_dir = tmp_path / "tmp"
    test_dir.mkdir()
    return test_dir


def test_process_sample_pass(simple_namespace, fastq_file, test_dir):
    sample_name = "test2-ont"
    sample_dict = (sample_name, {"single": fastq_file})

    d = process_sample(sample_dict, simple_namespace, test_dir)

    # Asserts all based on the test fastq file made with
    #  125 'reads', 100 D8, 25 D4
    assert isinstance(d, dict)
    assert d["sample"] == sample_name
    assert d["predicted_genotype"] == "D8"
    assert d["percent_supporting"] == 80
    assert d["supporting_count"] == 100
    assert d["total_count"] == 125
    assert d["fastq_1"] == str(fastq_file)


def test_process_sample_mixed(simple_namespace, fastq_file, test_dir):
    sample_name = "test3-ont"
    sample_dict = (sample_name, {"single": fastq_file})
    simple_namespace.majority_threshold = 90

    d = process_sample(sample_dict, simple_namespace, test_dir)

    assert isinstance(d, dict)
    assert isinstance(d, dict)
    assert d["sample"] == sample_name
    assert d["predicted_genotype"] == "Mixed"
    assert d["percent_supporting"] == 80
    assert d["supporting_count"] == 100
    assert d["total_count"] == 125
    assert d["fastq_1"] == str(fastq_file)
