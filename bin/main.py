#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
#import liams_simple_scaffold_pkg.prefil.prefil_files as fil

def parse_args(args):

    parser = argparse.ArgumentParser(
            description=".")
    required = parser.add_argument_group(
            'Required',
            'Project name')
    required.add_argument(
            "-f",
            "--foo",
            type=str,
            help="")

    optional = parser.add_argument_group(
            'Optional parameters',
            'All have defaults')
    optional.add_argument(
            "-b",
            "--bar",
            type=str,
            help="",
            default = None)

    return parser.parse_args(args)

def main():

    '''
    Make a simple python package scaffold
    '''

    args = parse_args(sys.argv[1:])
    fil.main(args)


if __name__ == "__main__":
    main()


