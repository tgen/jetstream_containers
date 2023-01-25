#!/usr/bin/env python3

import os
import sys
import argparse
from SigProfilerAssignment import Analyzer as Analyze

currentDirectory = os.getcwd()

def parse_arguments():
    """Parse arguments, validate and return the args"""

    parser = argparse.ArgumentParser(
        description='Compute mutational signatures.',
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-o', '--output', default=currentDirectory,
                        help='Path to and name of output directory - Can be relative or full path')

    parser.add_argument('-i', '--vcfpath',
                        help='Path to directory containing vcf file(s) - Can be relative or full path')

    parser.add_argument('-g', '--genome', default='GRCh38',
                        help='Optional definition of genome, defaults to GRCh38')

    parser.add_argument('-p', '--project', metavar="<arg>",
                        help='Name of the project')

    parser.add_argument('-e', '--exome', default=False, action="store_true",
                        help='Set if input is from exome')

    # prints help message when 0 arguments are entered
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    parser_args = parser.parse_args()

    return parser_args


def main():
    # Parse and validate arguments
    args = parse_arguments()

    #Getting SBS96 signatures
    signature_SBS96_DB="/sigprofiler_resources/COSMIC_v3.3.1_SBS_GRCh38.txt"
    try:
        Analyze.cosmic_fit(samples=args.vcfpath, output=args.output+"/SBS", input_type="vcf", context_type="96", signatures=None, signature_database=signature_SBS96_DB, genome_build="GRCh38", exome=args.exome)
    except KeyError as e:
        print(f'No {e} signatures found!')
    
    #Getting DBS signautures
    signature_DBS_DB="/sigprofiler_resources/COSMIC_v3.3_DBS_GRCh38.txt"
    try:
        Analyze.cosmic_fit(samples=args.vcfpath, output=args.output+"/DBS", input_type="vcf", context_type="DINUC", signatures=None, signature_database=signature_DBS_DB, genome_build="GRCh38", exome=args.exome, collapse_to_SBS96 = False)
    except KeyError as e:
        print(f'No {e} signatures found!')

if __name__ == '__main__':
    main()
