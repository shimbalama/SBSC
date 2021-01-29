# Somatic base substitution caller for LRS

Python tool that calls somatic base substitutions from matched tumour normal whole genome pileups.

## Getting Started

Pip install

### Prerequisites

samtools/1.10
htslib/1.10.2
tabix/0.2.6

Pileups from samtools
bgzip pileup_name
tabix -b 2 -e 2 pileup_name.gz

#### Required

Python 3

#### Optional

NA

### Installing

pip install SBSC

## Running the program

SBSC.py -h

## Authors

* **Liam McIntyre** - https://github.com/shimbalama/

## License

This project is licensed under the MIT License. See LICENSE

## Acknowledgments

QMIR medical genomics team
