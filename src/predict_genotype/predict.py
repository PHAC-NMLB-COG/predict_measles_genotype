#!/usr/bin/env python3

"""
This is a quick script to process either a single paired (or single) fastq file or
take a directory of fastq files and then run them against the reference database
to check what measles genotype the sample likely is
"""

import argparse
import csv
import glob
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from multiprocessing import Pool
from typing import Tuple

from . import __version__, program_description, program_name

# Get ref file from script path
script_dir = os.path.dirname(os.path.abspath(__file__))
default_ref = os.path.join(script_dir, "reference_data/measles_N450_genotypes.mmi")

# Illumina REGEX
ILLUMINA_FASTQ_REGEX = re.compile(r"^(.*?)(?:_S\d+)?(_R|\.pair)([12])(?:.+)?\.(fastq|fq)(\.gz)?$")
SIMPLE_FASTQ_REGEX = re.compile(r"^(.*?)\.(fastq|fq)(\.gz)?$")

# Really basic logging setup
logging.basicConfig(
    format="[%(asctime)s] %(levelname)s: %(message)s",
    level=logging.INFO,
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def check_external_dependencies(dependencies: list) -> None:
    missing = []
    for dep in dependencies:
        try:
            subprocess.check_output(["which", dep])
        except subprocess.CalledProcessError:
            missing.append(dep)

    if missing != []:
        logger.error(f"Missing external dependencies {missing}")
        sys.exit(1)


def find_fastq_files(input_dir: str) -> defaultdict:
    """Find and try to pair FASTQ files within the input directory

    Parameters
    ----------
        input_dir (str): Directory containing FASTQ files

    Returns
    -------
        defaultdict: Dictionary mapping sample names to their read file paths
    """
    ext_patterns = ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
    fastqs_found = []

    # Check that it is actually a directory
    if not os.path.isdir(input_dir):
        logger.error(f"Input Fastq Dir {input_dir} cannot be found")
        sys.exit(1)

    for pattern in ext_patterns:
        fastqs_found.extend(glob.glob(os.path.join(input_dir, pattern)))

    fastqs_found = sorted(fastqs_found)

    samples = defaultdict(dict)

    # Organize found fastqs into pairs
    for fq in fastqs_found:
        basename = os.path.basename(fq)
        fq_match = ILLUMINA_FASTQ_REGEX.match(basename)

        if fq_match:
            sample = fq_match.group(1)
            read_num = fq_match.group(3)
            if read_num == "1":
                samples[sample]["R1"] = fq
            elif read_num == "2":
                samples[sample]["R2"] = fq
            else:
                samples[sample]["single"] = fq
            continue

        # Check if we have ONT data or single data
        fq_simple_match = SIMPLE_FASTQ_REGEX.match(basename)
        if fq_simple_match:
            sample = fq_simple_match.group(1)
            samples[sample]["single"] = fq

    # Check that we actually had fastqs
    if samples == {}:
        logger.error(f"No fastq files found in {input_dir}")
        sys.exit(1)

    return samples


def map_reads_to_ref(ref: str, read1: str, bam_sorted: str, read2=None, threads=1) -> None:
    """Map and sort input reads and read pairs to the multifasta reference with minimap2 and samtools

    Parameters
    ----------
        ref (str): Path to multifasta reference
        read1 (str): Path to first input fastq file
        bam_sorted (str): Output sorted BAM file path
        read2 (str): Path to second fastq file if sample was paired
        threads (int): Number of threads to run minimap2 with
    """
    mode = "sr" if read2 else "map-ont"
    cmd_minimap = ["minimap2", "-ax", mode, "-t", str(threads), ref, read1]
    if read2:
        cmd_minimap.append(read2)

    # Running minimap2 and samtools sort back to back to save IO
    with open(bam_sorted, "wb") as bam_out:
        minimap2_proc = subprocess.run(
            cmd_minimap, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=True
        )
        # For removing supplementary and secondary
        samtools_view_proc = subprocess.run(
            ["samtools", "view", "-bS", "-F", "0x900"],
            input=minimap2_proc.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            check=True,
        )
        # Sorting and writing to final file
        subprocess.run(
            ["samtools", "sort"],
            input=samtools_view_proc.stdout,
            stdout=bam_out,
            stderr=subprocess.DEVNULL,
            check=True,
        )

    # Index final file
    subprocess.run(
        ["samtools", "index", bam_sorted],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True,
    )


def summarize_genotype_hits(bam_sorted: str) -> Tuple[defaultdict, int]:
    """Summarize how many hits each reference genotype got along with compressing genotypes with multiple references together

    Parameters
    ----------
        bam_sorted (str): Sorted BAM file path

    Returns
    -------
        Tuple
            defaultdict: A dictionary of summarized genotypes

            int: Total reads
    """
    cmd = ["samtools", "idxstats", bam_sorted]
    result = subprocess.check_output(cmd).decode().strip()
    hits = defaultdict(int)
    total = 0

    for line in result.split("\n"):
        ref, ref_len, num_mapped, num_unmapped = line.split("\t")

        # Multifasta ref IDs are formatted as genotype_seq_name
        #  So we can split on the '_' to merge them together
        genotype = ref.split("_")[0]
        hits[genotype] += int(num_mapped)
        total += int(num_mapped)
    return hits, total


def process_sample(sample_tuple: Tuple[str, dict], args: argparse.Namespace, tmpdir: str) -> dict:
    """Process sample to determine and return predicted genotype and information

    Parameters
    ----------
        sample_tuple (Tuple[str, dict]): String is the sample name, dict is
            {'R1': '<PATH>', 'R2': '<PATH>'} or {'single': '<PATH>'}
        args (argparse.Namespace): Parsed args
        tmpdir (str): Path to tmpdir to write data to

    Returns
    -------
        dict: A dictionary of the sample, predicted genotype, file paths, and supporting data
    """
    sample, files = sample_tuple

    bam_sorted = os.path.join(tmpdir, f"{sample}.sorted.bam")
    read1 = files.get("R1", files.get("single"))
    read2 = files.get("R2")

    map_reads_to_ref(args.reference, read1, bam_sorted, read2, args.threads)

    hits, total = summarize_genotype_hits(bam_sorted)
    sorted_hits = sorted(hits.items(), key=lambda x: x[1], reverse=True)

    top_5_info = ";".join([f"{ref}:{count}" for ref, count in sorted_hits[:5]])
    top_ref, top_count = sorted_hits[0]

    percent_support = 0
    if total > 0:
        percent_support = round((top_count / total * 100), 2)
        if percent_support < args.majority_threshold:
            top_ref = "Mixed"

    if total < args.min_reads:
        top_ref = "NA"

    if read2:
        read2 = os.path.abspath(read2)

    # CL statement for prediction
    predict_statement = f"Predicted: {top_ref} for sample {sample}"
    if top_ref == "NA":
        predict_statement = f"Not enough reads mapped to N450 Database to predict sample {sample}"
    logger.info(predict_statement)

    return {
        "sample": sample,
        "predicted_genotype": top_ref,
        "percent_supporting": percent_support,
        "supporting_count": top_count,
        "total_count": total,
        "top_5_info": top_5_info,
        "fastq_1": os.path.abspath(read1),
        "fastq_2": read2,
    }


def main():
    """Main entry point"""
    # Slurm check for max threads
    if "SLURM_CPUS_PER_TASK" in os.environ:
        max_threads = int(os.environ["SLURM_CPUS_PER_TASK"])
    else:
        max_threads = os.cpu_count()

    # Generate parser and parse args
    parser = argparse.ArgumentParser(prog=program_name, description=program_description)
    req_group = parser.add_argument_group(
        "Required Input Data", "How to find/input fastq data into the pipeline"
    )
    exclusive_group = req_group.add_mutually_exclusive_group(required=True)
    exclusive_group.add_argument("-d", "--directory", help="Input directory with FASTQ files")
    exclusive_group.add_argument("-1", "--read1", help="Illumina R1 or ONT Single-end FASTQ file")

    opt_group = parser.add_argument_group("Optional Args")
    opt_group.add_argument("-2", "--read2", help="Illumina R2 read file")
    opt_group.add_argument(
        "-r",
        "--reference",
        default=default_ref,
        help="Multifasta reference file or minimap2 index file (default: %(default)s)",
    )
    opt_group.add_argument(
        "-o",
        "--output",
        default="predictions.csv",
        help="Output filename for genotype predicitons (default: predictions.csv)",
    )
    opt_group.add_argument(
        "-m",
        "--majority-threshold",
        default=60,
        help="Threshold for marking the input as mixed (default: 60)",
    )
    opt_group.add_argument(
        "-l",
        "--min-reads",
        default=50,
        help="Minimum total reads mapping to any N450 reference to call a genotype (default: 50)",
    )
    opt_group.add_argument(
        "-s",
        "--samplesheet-split",
        action="store_true",
        help="Create and split samplesheet by genotype",
    )
    opt_group.add_argument(
        "-t",
        "--threads",
        type=int,
        default=max_threads,
        help="Number of concurrent threads to use when processing (default: {})".format(
            max_threads
        ),
    )
    parser.add_argument(
        "-v",
        "--version",
        help="Outputs the current version and exits",
        action="version",
        version="%(prog)s {}".format(__version__),
    )
    args = parser.parse_args()

    # Dependency Checks
    external_deps = ["minimap2", "samtools"]
    check_external_dependencies(external_deps)

    # Threading checks and calcs after args
    if args.threads > max_threads:
        logger.warning(
            f"Requested {args.threads} threads but only {max_threads} available. Using {max_threads}."
        )
        args.threads = max_threads
    elif args.threads < 1:
        logger.warning(f"Invalid thread count ({args.threads}). Using 1.")
        args.threads = 1

    # Minimap2 max threads = 2 for now, so use that to determine concurrent processes
    max_concurrent = int(max(args.threads / 2, 1))  # Number to run at once
    minimap2_threads = int(max(args.threads / max_concurrent, 1))
    args.threads = minimap2_threads  # To easily pass to minimap2 command

    logger.info(
        f"Max Processes set to: {max_concurrent} using a total of {args.threads} thread(s) each"
    )

    # Tmpdir for writing data to remove later
    tmpdir = tempfile.mkdtemp()

    # Based on input, determine samples and what type of sample they are
    if args.directory:
        samples = find_fastq_files(args.directory)
    else:
        basename = os.path.basename(args.read1)
        fq_match = ILLUMINA_FASTQ_REGEX.match(basename)
        if fq_match:
            sample_name = fq_match.group(1)
            samples = {
                sample_name: {"R1": args.read1, "R2": args.read2}
                if args.read2
                else {"single": args.read1}
            }
        else:
            logger.error(f"Input fastq file {args.read1} did not match expected Illumina pattern")
            sys.exit(1)

    # MultiProcess and Track Running
    sample_items = list(samples.items())
    with Pool(processes=max_concurrent) as pool:
        predictions = list(pool.starmap(process_sample, [(s, args, tmpdir) for s in sample_items]))

    # Don't write nothing
    if len(predictions) == 0:
        logger.error("No data found in input args")
        sys.exit(1)

    # Output to file
    csv_out = args.output
    with open(csv_out, "w", newline="") as csvfile:
        fieldnames = list(predictions[0].keys())
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(predictions)

    # Split samplesheet output if wanted
    if args.samplesheet_split:
        # Reformat to group by genotype
        #  Each key is the genotype and the value is the list of matching dicts
        grouped_predictions = defaultdict(list)
        keep_keys = ["sample", "fastq_1", "fastq_2"]
        for d in predictions:
            key = d["predicted_genotype"]
            filtered_d = {k: d[k] for k in keep_keys if k in d}
            grouped_predictions[key].append(filtered_d)

        # Create our outputs
        for prediction, group_list in grouped_predictions.items():
            csv_out = f"{prediction}_samplesheet.csv"
            with open(csv_out, "w", newline="") as csvfile:
                fieldnames = list(group_list[0].keys())
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(group_list)

    # Remove intermediates
    shutil.rmtree(tmpdir)


# Run script
if __name__ == "__main__":
    main()
