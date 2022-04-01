#!/usr/bin/env python
# -*-coding: utf8 -*-
from fastapi.middleware.cors import CORSMiddleware
"""
Library for input fasta file quality control in hCov-19 analysis.
Copyright 2021 by Regina Kalcsevszki (kalcsevszkiregi@gmail.com), Balázs Ligeti (obalasz@gmail.com). All rights reserved.

"""

import os
from os.path import join, isfile
import sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import date

def run_input_qc(fasta_file_path, output_fasta_path, id_mapper_file, CONFIG):

    print('Running QC!')

    max_sample_size = CONFIG["InputQC"]["MaxSampleSize"]
    sample_id_prefix = CONFIG["InputQC"]["SampleIDPrefix"]
    min_seq_len = CONFIG["InputQC"]["MinSequenceLength"]
    max_seq_len = CONFIG["InputQC"]["MaxSequenceLength"]
    min_acgt_ratio = CONFIG["InputQC"]["MinACGTRatio"]
    
    id_mapper = []
    with open(fasta_file_path) as handle:
        fasta = list(SeqIO.parse(handle,'fasta'))
        if len(fasta) > CONFIG["InputQC"]["MaxSampleSize"]:
            raise(Exception("There are more than "+str(max_sample_size)+" samples in the fasta file!"))
        sample_number = 1
        valid_records = []
        for record in fasta:
            new_id = get_id_in_gisaid_format(record, sample_number, sample_id_prefix)
            old_id = record.id
            record.id = new_id
            record.description = ''
            id_mapper.append([old_id, new_id])
            if is_length_valid(record, min_seq_len, max_seq_len) and is_quality_valid(record, min_acgt_ratio):
                valid_records.append(record)
                sample_number = sample_number + 1
    if len(fasta) < 1:
        with open(fasta_file_path) as handle:
            seq = Seq(handle.read().replace('\n', ''))
            if len(seq) > 0: 
                desc = os.path.splitext(os.path.basename(fasta_file_path))[0]
                record = SeqIO.SeqRecord(seq, id=desc)
                path_to_new_fasta = os.path.splitext(fasta_file_path)[0] + '.fasta'
                SeqIO.write(record, path_to_new_fasta, "fasta")
                return run_input_qc(path_to_new_fasta, output_fasta_path,id_mapper_file, CONFIG)
            else:
                raise(Exception("Input file format is invalid!")) 
    else:
        pd.DataFrame(id_mapper, columns =['Old ID', 'New ID']).to_csv(id_mapper_file, sep='\t')
        SeqIO.write(valid_records, output_fasta_path, "fasta")
        return valid_records
    
def change_spec_characters_in_id(input_str):
    return input_str.replace('|', '_').replace(' ', '')

def get_id_in_gisaid_format(record, sample_number, sample_id_prefix):
    today = date.today().strftime("%Y-%m-%d")   
    return '|'.join([change_spec_characters_in_id(record.description), sample_id_prefix + "_" + str(sample_number).zfill(7), today])

def is_length_valid(record, min_seq_len, max_seq_len):
    length_of_sequence = len(record.seq)
    return length_of_sequence > min_seq_len and length_of_sequence < max_seq_len

def is_quality_valid(record, min_acgt_ratio):
    seq = record.seq.upper()
    A = seq.count('A')
    T = seq.count('T')
    C = seq.count('C')
    G = seq.count('G')
    valid_letter_count = A + T + C + G
    gc_content = (G+C) / (valid_letter_count)

    return valid_letter_count/len(seq) >= min_acgt_ratio

