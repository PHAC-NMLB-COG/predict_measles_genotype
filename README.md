# Predict Measles Genotype

A really quick script that takes in either an input FASTQ file (or an Illumina FASTQ pair) or a directory of FASTQ files and predicts which genotype they are based on mapping with [`minimap2`](https://github.com/lh3/minimap2) to the measles WHO N450 reference database with some supplemented additional N450 diversity.

The end output is a file called `predicitions.csv` that has the sample name, predicted genotype, stats on the genotype call, and the path to the used FASTQ file(s). See [output](#output) for examples

## Index

- [Installing](#installing)
  - [Conda](#conda)
  - [Just the Tool](#just-the-tool)
- [Running](#running)
  - [Basic](#basic)
  - [Generate Split Samplesheets](#generate-split-samplesheets)
  - [All Parameters](#all-parameters)
- [Outputs](#outputs)
  - [`predictions.csv`](#predicitonscsv)
  - [`<GENOTYPE>_samplesheet.csv`](#genotype_samplesheetcsv)
- [Database Info](#database-info)
- [Development](#development)
- [Contributing](#contributing)
- [Legal](#legal)

## Installing

### Conda

To install, clone the repo and create the conda enviroment with the environment file:

```bash
# 1. Clone
git clone https://github.com/PHAC-NMLB-COG/predict_measles_genotype.git

# 2. cd
cd predict_measles_genotype

# 3. Env create
conda env create -f environment.yml
```

### Just the Tool

If you already have the needed dependencies below listed installed you can clone the repo and just pip install.

Required dependencies:
- `minimap2`
- `samtools`
- `python3`

Install:

```bash
# 1. Clone
git clone https://github.com/PHAC-NMLB-COG/predict_measles_genotype.git

# 2. cd
cd predict_measles_genotype

# 3. Pip install
pip install .
```

## Running

### Basic

To run on a directory of FASTQ files:

```bash
predict_genotype -d <DIRECTORY>
```

The output from here can be used to run our measles analysis pipeline, the [`MeaSeq Pipeline`](https://github.com/phac-nml/measeq), or you can split the samplesheet by genotype and use that as input

### Generate Split Samplesheets

To generate additional split samplesheets called `<GENOTYPE>_samplesheet.csv` you can run:

```bash
predict_genotype -d <DIRECTORY> --samplesheet-split
```

### All Parameters

```bash
usage: predict_genotype [-h] (-d DIRECTORY | -1 READ1) [-2 READ2] [-r REFERENCE] [-o OUTPUT] [-m MAJORITY_THRESHOLD] [-l MIN_READS] [-s] [-t THREADS] [-v]

Measles N450 Genotyping Predictions from FASTQ data using the WHO N450 genotyping samples and minimap2

options:
  -h, --help            show this help message and exit
  -v, --version         Outputs the current version and exits

Required Input Data:
  How to find/input fastq data into the pipeline

  -d, --directory DIRECTORY
                        Input directory with FASTQ files
  -1, --read1 READ1     Illumina R1 or ONT Single-end FASTQ file

Optional Args:
  -2, --read2 READ2     Illumina R2 read file
  -r, --reference REFERENCE
                        Multifasta reference file or minimap2 index file (default:
                        <INSTALL_ENV>/lib/python3.X/site-
                        packages/predict_genotype/reference_data/measles_N450_genotypes.mmi)
  -o, --output OUTPUT   Output filename for genotype predicitons (default: predictions.csv)
  -m, --majority-threshold MAJORITY_THRESHOLD
                        Threshold for marking the input as mixed (default: 60)
  -l, --min-reads MIN_READS
                        Minimum total reads mapping to any N450 reference to call a genotype (default: 50)
  -s, --samplesheet-split
                        Create and split samplesheet by genotype
  -t, --threads THREADS
                        Number of concurrent threads to use when processing (default: os.cpu_count())
```

## Outputs

### `Predicitons.csv`

Format:

| sample | predicted_genotype | percent_supporting | supporting_count | total_count | top_5_info | fastq_1 | fastq_2 |
| - | - | - | - | - | - | - | - |
| sample1 | B3    | 84.20 | 1610 | 1912 | B3:1610;C2:136;D8:38;B2:33;D4:29 | <ABS PATH> | <ABS PATH> |
| sample2 | D8    | 86.65 | 2110 | 2435 | D8:2110;E:88;D4:56;B3:52;B2:42   | <ABS PATH> | <ABS PATH> |
| sample3 | A     | 88.19 | 3749 | 4251 | B3:1610;C2:136;D8:38;B2:33;D4:29 | <ABS PATH> | <ABS PATH> |
| sample3 | NA    | 100   | 10   | 10   | D3:10;A:0;C2:0;D4:0;B2:0         | <ABS PATH> | <ABS PATH> |
| sample3 | Mixed | 55.20 | 552  | 1000 | B3:552;D8:448;A:0;B2:0;D4:0      | <ABS PATH> |            |

*Note* that NA can have a 100% supporting because it has such few reads

### `<GENOTYPE>_samplesheet.csv`

Format:

| sample | fastq_1 | fastq_2 |
| - | - | - |
| sample1 | <ABS PATH> | <ABS PATH> |
| sample2 | <ABS PATH> |            |
| sample3 | <ABS PATH> | <ABS PATH> |

This can go right into the [MeaSeq Pipeline](https://github.com/phac-nml/measeq) as well

## Database Info

The database used was based on this [WHO Reference Table](https://cdn.who.int/media/docs/default-source/immunization/vpd_surveillance/lab_networks/measles_rubella/manual/table_7.pdf?sfvrsn=13c20a88_5) along with some additional sequences added to capture more diversity in some of the genotypes aiding in predictions and preventing some mixed calls.

## Development

Using [Ruff](https://github.com/astral-sh/ruff) for code linting and formatting

Tests currently written with pytest. To run all tests, external dependencies need to be available on the path

## Contributing

Contributions are welcome through creating PRs or Issues

## Legal

Copyright 2025 Government of Canada

Licensed under the MIT License (the "License"); you may not use this work except in compliance with the License. You may obtain a copy of the License at:

https://opensource.org/license/mit/

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
