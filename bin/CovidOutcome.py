#!/usr/bin/env python
# -*-coding: utf8 -*-

# Created on 2022

import os
from os.path import join, isfile, abspath, dirname
import inspect
import yaml
import argparse
import pandas
from uuid import UUID
from uuid import uuid4

def command_line_parser_pipeline_mngt():
    """
    Parsing input arguments for covid outcome.
    :return:
    """
    parser = argparse.ArgumentParser(description="Covid Outcome application\n"
                                                 "For the detailed description about the outputs and format, please visit our github page: ...\n",
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("input", help="Path to the sequences. The sequence should be in fasta format", type=str)
    parser.add_argument("output", help="Path to the ", type=str, default='./outcomes')

    parser.add_argument("--config_file",
                        help="Path to the config file, that contains the neccessery paremeters for the analysis and for running the pipeline properly.",
                        type=str,
                        default='data/PipelineConfig.yml')
    parser.add_argument("--age data",
                        help="OPTIONAL. Path to the age data. ")
    args = parser.parse_args()
    return args, parser

def get_absolut_file_path(input_args):
    """
    Get the absolut path to the config file
    :param input_args:   input arguments for the CovidOutcome
    :return: the absoutl path, absolut path to the config file
    """
    abs_path = abspath((inspect.stack()[0])[1])
    base_dir = join(dirname(abs_path), '..')

    if input_args.config_file == 'data/PipelineConfig.yml':
        covid_config_path = join(base_dir, input_args.config_file)
    else:
        covid_config_path = input_args.config_file

    return base_dir, covid_config_path






def set_envirement_variables(input_args):
    """
    Setting up the envirement variables needed for CovidOutcome
    :param input_args:   set of input arguments for the caller
    """

    base_dir,covid_config_path = get_absolut_file_path(input_args)
    os.environ["BASEPATH"] = base_dir
    os.environ["COVIDCONFIG"] = covid_config_path

def check_input_files_and_output(input_args):
    """
    Checking the input, output paramaters
    :param input_args:  input arguments provided by the caller
    :return:
    """
    base_dir, covid_config_path = get_absolut_file_path(input_args)
    output_dir = input_args.output
    input_sequence_file = input_args.input

    if not os.path.exists(covid_config_path):
        raise argparse.ArgumentTypeError("Configuration file is missing. {0} does not exist".format(covid_config_path))

    if not os.path.exists(output_dir):
        try:
            print('Output directory does not exists.')
            os.mkdir(output_dir)
            print('The directory: {0} has been created successfully.'.format(output_dir))
        except PermissionError:
            raise argparse.ArgumentTypeError("PermissionError: [Errno 13] Permission denied: {0}. Check the permission".format(output_dir))
    # Check input file:
    if not os.path.exists(input_sequence_file):
        raise argparse.ArgumentTypeError("Input sequence file is missing. {0} does not exist".format(input_sequence_file))


def read_sequence_file(input_args, input_id, CovidOutcome):
    print('Reading input file.')

    input_file = input_args.input
    max_file_size = CovidOutcome.CONFIG['MaxFileSize']
    file_size = os.path.getsize(input_file)
    if file_size > max_file_size:
        print('Too large file, interrupt, modify the config file if needed. ')
        raise Exception("Too large input sequence file. The maximum allowed file size is %d." %
                   CovidOutcome.CONFIG['MaxFileSize'],
        )
    with open(input_file) as in_file:
        sequences = in_file.read()

    number_of_headers =sequences.count('>')
    number_of_lines = len(sequences.split())

    CovidOutcome.set_session_data_null(input_id, sequences)

    act_session_data = CovidOutcome.sessions[input_id]
    act_session_data['raw_sequence_data'] = sequences
    sequence_response = CovidOutcome.check_input_sequences(input_id)
    N_valid_seqs = CovidOutcome.get_seqqc_number_of_valid_sequences(input_id)

    if not N_valid_seqs > 0:
        print('The number of valid sequences after QC is zero, interrupt. \
        Please check the output QC, modify the parameters, or the input if neccessary.')
        raise Exception('Input error.')

def run_ml_evaluation(input_id, CovidOutcome):
    ml_results = CovidOutcome.evaluate_ml_model(input_id, 'deep')
    ml_results_file = CovidOutcome.sessions[input_id]['ml_results_file']
    if not (os.path.exists(ml_results_file) and os.path.getsize(ml_results_file) > 0):
        ml_results.to_csv(ml_results_file, sep='\t', index=False)


def run_pipeline(input_args):
    '''
    Pipeline for covid outcome. Main steps: i) create a database record; ii) QC of the input sequence;
    iii) mutation calling and annotation; iv) evaluation ml results model.
    :param input_args:  input arguments provdided by the caller
    :return:
    '''
    print('Start analysis')
    import CovidOutcome
    input_id = str(uuid4())
    CovidOutcome.add_new_session(input_id)
    read_sequence_file(input_args,input_id, CovidOutcome)
    CovidOutcome.mutation_calling_procedure(input_id)
    run_ml_evaluation(input_id, CovidOutcome)


def main():
    print("Starting COVID OUTCOME standalone application! \n")
    args, parser = command_line_parser_pipeline_mngt()
    check_input_files_and_output(args)
    set_envirement_variables(args)
    run_pipeline(args)


if __name__ == "__main__":
    main()

