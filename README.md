# Somatic base substitution caller for LRS

Python tool that calls somatic base substitutions from matched tumour-normal whole genome pileups.

#### Required

Python 3.6+

## Install

Install into a virtualenv is recommended, _e.g._:

    python3 -m venv SBSC
    source SBSC/bin/activate
    tempdir=$(mktemp -d)
    git clone https://github.com/shimbalama/SBSC $tempdir
    cd $tempdir
    tox # run tests: this is optional and requires tox be installed
    pip3 install -r requirements.txt .

### Prerequisites

Tabix-indexed pileups from samtools. `samtools mpileup` should be run with the following options:

    samtools mpileup \
        --fasta-ref [reference.fa] \
        --count-orphans \
        --no-BAQ \
        --max-depth 500 \
        --min-BQ 1 \
        --ignore-RG \
        --ignore-overlaps \
        -a -a \
        --output-BP \
        --output-MQ \
        --output-extra FLAG,QNAME \
        [bam] | bgzip > [pileup.gz]

`tabix` should be run with the following options:

    tabix -b 2 -e 2 [pileup.gz]

When the tool runs `.tbi` index files are presumed to exist adjacent to input pileup files.

## Running the program

    # help
    SBSCall.py --help

    # basic example
    SBSCall.py call \
        -r /path/to/reference \
        -c /path/to/cancer_pileup.gz \
        -n /path/to/normal_pileup.gz 

## Authors

* **Liam McIntyre** - https://github.com/shimbalama/

## License

This project is licensed under the MIT License. See LICENSE

## Acknowledgments

QIMR Berghofer medical genomics team
