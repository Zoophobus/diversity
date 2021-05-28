#!/usr/bin/env python3.7

import zoo_diversity_analysis

import re,glob,argparse

parser = argparse.ArgumentParser(description="Process a shorah read file to filter and calculate the average number of polymorphic sites\n")
essentialArgs = parser.add_argument_group('Files', 'Required arguments')
essentialArgs.add_argument('-i','--input',metavar='input',required=True,type=str,dest='shorah_input',help="The reads.fas output file from ShoRAH")
essentialArgs.add_argument('-o','--output',metavar='output',required=True,type=str,dest='output_file',help="The file to print the results to")
#essentialArgs.add_argument('-r','--reference',metavar='reference',required=True,type=str,dest='reference_file',help="The file containing the reference sequence used in ShoRAH")
nonEssentialArgs = parser.add_argument_group('Conditions','Optional arguments')

nonEssentialArgs.add_argument('-l','--pipeline-shorah',metavar='shorah',required=False,type=bool,dest='shorah',help="A boolean to indicate if the input is from the shorah pipeline, by default this is False.")
nonEssentialArgs.add_argument('-d','--directory',metavar='directory',required=False,type=bool,dest='input_is_directory',help="A boolean to indicate if the input is a directory or not, by default this is False.")
nonEssentialArgs.add_argument('-c','--continuous-mismatches',metavar='continuous',required=False,type=int,dest='contiguous_filter',help="The number of consecutive mismatches after which a particular read should be dropped.")
nonEssentialArgs.add_argument('-C','--minimum-predicted-coverage',metavar='coverage',required=False,type=float,dest='coverage',help="The average number of expected read coverages for the predicted local haplotypes.")
nonEssentialArgs.add_argument('-p','--minimum-sequence-posterior-probability',metavar='min_prob',required=False,type=float,dest='min_prob',help="The required posterior probability for the predicted local haplotypes.")
nonEssentialArgs.add_argument('-f','--identity-filter',metavar='identity',required=False,type=float,dest='identity',help="The required number of matches for a read to be included.")
nonEssentialArgs.add_argument('-M','--mega_output',metavar='MEGA',required=False,type=str,dest='MEGA',help="A prefix for outputing reads in short local alignments for subsequent analysis.")
# NOTE the following two arguments are dealt with here
nonEssentialArgs.add_argument('-s','--input-files-suffix',metavar='input_suffix',required=False,type=str,dest='input_suffix',help="The suffix, or extension for the files that are to be used as input, by default .fas.")
nonEssentialArgs.add_argument('-O','--merge-output',metavar='merge_output',required=False,type=bool,dest='merge_output',help="A boolean to indicate whether the output should be merged or not, either True or False (default).")
nonEssentialArgs.add_argument('-w','--window-length',metavar='window_length',required=False,type=int,dest='window_length',help="The length of the windows to be used (by default 66).")
nonEssentialArgs.add_argument('-W','--wattersons-statistic',metavar='wattersons_statistic',required=False,type=bool,dest='wattersons_statistic',help="A boolean to indicate whether to use wattersons theta (True), or pi (if -W is not provided, default)")
nonEssentialArgs.add_argument('-N','--minimum-no-haplotypes',metavar='minimum_haplotypes',required=False,type=int,dest='minimum_haplotypes',help="The minimum required number of haplotypes for the diversity analysis to ve performed (default 2)")

arguments = parser.parse_args()

# use the input_is_directory argument to either find the files, or to use the input argument directly

if arguments.input_suffix is None:
    arguments.input_suffix = ".fas"
if arguments.merge_output is None:
    arguments.merge_output = False
if arguments.input_is_directory is None:
    arguments.input_is_directory = False
if arguments.shorah is None:
    arguments.shorah = False
if arguments.wattersons_statistic is None:
    arguments.wattersons_statistic = False

if arguments.minimum_haplotypes is None:
    arguments.minimum_haplotypes = 2


if arguments.contiguous_filter is None:
    arguments.contiguous_filter = 10
if arguments.min_prob is None:
    arguments.min_prob = 0.95
if arguments.identity is None:
    arguments.identity = 0.95
if arguments.coverage is None:
    arguments.coverage = 0
if arguments.input_suffix is None:
    arguments.input_suffix = ""
if arguments.merge_output is None:
    arguments.merge_output = False
if arguments.window_length is None:
    arguments.window_length = 66




if arguments.input_is_directory:
    file_list = glob.glob(arguments.shorah_input + "/*" + arguments.input_suffix)
    for aFile in file_list:
        if arguments.merge_output:
            output = open(arguments.output_file,"a")
            zoo_diversity_analysis.main(aFile,output,arguments.contiguous_filter,arguments.identity,arguments.MEGA,arguments.window_length,arguments.wattersons_statistic,arguments.coverage,arguments.min_prob,arguments.minimum_haplotypes,arguments.shorah)
        else:
            exp = re.compile(r".*/(.*)\.\w+$")
            excludingExtension = re.match(exp,aFile)
            theFile = excludingExtension.group(1) + arguments.output_file
            output = open(theFile,"w")
            zoo_diversity_analysis.main(aFile,output,arguments.contiguous_filter,arguments.identity,arguments.MEGA,arguments.window_length,arguments.wattersons_statistic,arguments.coverage,arguments.min_prob,arguments.minimum_haplotypes,arguments.shorah)
else:
    output = open(arguments.output_file,"a")
    zoo_diversity_analysis.main(arguments.shorah_input,output,arguments.contiguous_filter,arguments.identity,arguments.MEGA,arguments.window_length,arguments.wattersons_statistic,arguments.coverage,arguments.min_prob,arguments.minimum_haplotypes,arguments.shorah)
