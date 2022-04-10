#!/usr/bin/env python
# -*-coding: utf8 -*-

"""
Library for Pipeline config.
Copyright 2021 by Regina Kalcsevszki (kalcsevszkiregi@gmail.com), Balázs Ligeti (obalasz@gmail.com). All rights reserved.

"""

import os
from os.path import join
import pandas as pd
import json
import gzip
import sys
from sys import exit

import yaml


def yaml_parsers(yaml_config_file):
    """
    Parse the configuration file!
    :param yaml_config_file:
    :return:
    """
    with open(yaml_config_file, 'r') as ymlfile:
        cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)
        return cfg

class Config():

    max_fasta_header_size = 2000
    offset = 20
    char_byte = 2
    valid_qc_steps = ["check_ids",
                     'check_uppercase',
                     'check_sequence_length',
                     'check_nucleotid_proportions',
                     'check_congruence']

    def __init__(self, config_path, base_path):
        #try:
            self.CONFIG = yaml_parsers(config_path)
            self.set_base_path(base_path)
            self.creating_directory_structure()
            self.get_max_file_input_sizes()
            self.get_max_input_char_count()
            self.check_qc_steps()
            self.load_database_info()
        #except:
        #    print("Config file not found or can not be parsed:", config_path)
        #    exit(2)
    
    def get(self):
        return self.CONFIG

    def check_qc_steps(self):

        if self.CONFIG['InputQC']['QCSteps'] and (len(set(self.CONFIG['InputQC']['QCSteps']) - set(Config.valid_qc_steps)) > 0):
            raise ValueError('Invalid QC step in the config file.')



    # Not tested
    # Calculating the limit for the input file size
    def get_max_file_input_sizes(self):
        Nseq =  self.CONFIG['MaxSampleSize']
        max_seq_length = self.CONFIG['InputQC']['Sequence_length']['MaxSequenceLength']
        self.CONFIG["MaxFileSize"] = (Nseq*max_seq_length + Nseq*Config.max_fasta_header_size + Nseq + Config.offset)*Config.char_byte

    # Setting upper limit for the input size
    def get_max_input_char_count(self):
        Nseq =  self.CONFIG['MaxSampleSize']
        max_seq_length = self.CONFIG['InputQC']['Sequence_length']['MaxSequenceLength']
        self.CONFIG["MaxSequenceSequenceChars"] = (Nseq*max_seq_length + Nseq*Config.max_fasta_header_size + Nseq + Config.offset)


    def set_base_path(self, base_path):
        self.CONFIG["BasePath"] = base_path

    def set_mummer_path(self, mummer_path):
        self.CONFIG["InputQC"]["MummerPath"] = mummer_path
        
    def set_mafft_path(self, mafft_path):
        self.CONFIG["MulipleAlignment"]["MAFFTPath"] = mafft_path

    def creating_directory_structure(self):
        '''
        Creating the directory structure to hold the temporary results!
        :param cfg:
        :return:
        '''
        act_output_folder = self.CONFIG['OutputFolders']['output_folder']
        new_output_folder = act_output_folder
        if not act_output_folder.startswith('/'):
            new_output_folder = os.path.join(self.CONFIG['BasePath'], act_output_folder)

        self.CONFIG['OutputFolders']['output_folder'] =new_output_folder
        os.makedirs(new_output_folder, exist_ok=True)

        # Testing write permission?
    def load_database_info(self):
        base_path = self.CONFIG['BasePath']
        suggestion_db_file = join(base_path, self.CONFIG['Database']['mutation_db_suggestions_file'])
        mut_suggestion_db = pd.read_csv(suggestion_db_file)
        mut_suggestion_db['mut_short_lowercase'] = mut_suggestion_db.apply(lambda x: str(x['mut_short']).lower(),
                                                                           axis=1)
        self.CONFIG['Database']['mutation_db_suggestions'] = mut_suggestion_db
        # Loading the mutation mapping information
        print('Loading the mutation mapping information into the memory!')
        mutation_name_mapping_file = join(base_path, self.CONFIG['Database']['mutation_db_mapping_file'])
        with gzip.open(mutation_name_mapping_file, 'rt', encoding='UTF-8') as zipfile:
            mutation_name_mapping = json.load(zipfile)
        self.CONFIG['Database']['mutation_name_mapping'] = mutation_name_mapping

        # Loading the mutation annotations
        print('Loading the current mutation annotation information into the memory!')
        annotations_db_file = join(base_path, self.CONFIG['Database']['mutation_deduplicated_annotation_file'])
        mut_annotation_db = pd.read_csv(annotations_db_file, sep='\t')
        self.CONFIG['Database']['mutation_annotation'] = mut_annotation_db


