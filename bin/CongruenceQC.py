#!/usr/bin/env python
# -*-coding: utf8 -*-

"""
Library for congruence with reference quality control in hCov-19 analysis.
Copyright 2021 by Regina Kalcsevszki (kalcsevszkiregi@gmail.com), Balázs Ligeti (obalasz@gmail.com). All rights reserved.

"""

import os
from os.path import join, isfile
import pandas as pd
from Bio import SeqIO

def run_congruence_qc(fasta_file_path, mummer_output_folder_path, output_fasta_path, output_congruance_table_path, CONFIG, print_mummer_cmd = True):
    congruence_threshold = CONFIG['InputQC']["CongruenceQC"]["CongruenceThreshold"]
    reference_length = CONFIG["Reference"]["Length"]
    reference_path = join(CONFIG["BasePath"], CONFIG["Reference"]["ReferencePathFromBase"])
    mummer_path = CONFIG["InputQC"]["CongruenceQC"]["MummerPath"]
    
    try:
        os.mkdir(mummer_output_folder_path)
        print("Directory " , output_folder_path ,  " created")
    except FileExistsError:
        pass
    
    with open(fasta_file_path) as handle:
        records = list(SeqIO.parse(handle,'fasta'))
    
        delta_output_path = run_mummer(mummer_path, reference_path, fasta_file_path, mummer_output_folder_path, print_mummer_cmd)

        congruence_table = process_congruence_table(delta_output_path, records, reference_length, congruence_threshold)

        valid_records = select_valid_records(congruence_table, records)
        
        congruence_table.to_csv(output_congruance_table_path, sep='\t')
        SeqIO.write(valid_records, output_fasta_path, "fasta")
        
        return valid_records, congruence_table

def run_mummer(mummer_path, reference_path, fasta_path, mummer_output_folder_path, print_mummer_cmd=True):
    fasta_file_name = os.path.basename(fasta_path)
    output_pre = join(mummer_output_folder_path, fasta_file_name)
    nucmer_output_path = join(mummer_output_folder_path, output_pre + ".delta")
    filtered_nucmer_output_path = join(mummer_output_folder_path, output_pre + ".filtered.delta")
    show_coords_output_path = join(mummer_output_folder_path, output_pre + ".coords")
    
    # nucmer <referemce fasta file> <subjects fasta file> -p <prefix of the output delta file>
    nucmer_cmd = '{0} {1} {2} -p {3}'.format(join(mummer_path, 'nucmer'), reference_path, fasta_path, output_pre) 
    if print_mummer_cmd: 
        print(nucmer_cmd)
    return_value = os.system(nucmer_cmd)
    if return_value != 0:
        raise(Exception("Something went wrong while running nucmer:" + nucmer_cmd))
   
    # delta-filter -g <delta file> > <output filtered delata file> 
    # -g: 1-to-1 global alignment not allowing rearrangements
    filter_cmd = '{0} -g {1}>{2}'.format(join(mummer_path, 'delta-filter'), nucmer_output_path, filtered_nucmer_output_path) 
    if print_mummer_cmd:
        print(filter_cmd)
    return_value = os.system(filter_cmd)
    if return_value != 0:
        raise(Exception("Something went wrong while running delta-filter:" + filter_cmd))
    
    # show-coords <filtered delata file> > <show-coords out file>
    show_coords_cmd = '{0} {1}>{2}'.format(join(mummer_path, 'show-coords'), filtered_nucmer_output_path, show_coords_output_path) 
    if print_mummer_cmd:
        print(show_coords_cmd)
    return_value = os.system(show_coords_cmd)
    if return_value != 0:
        raise(Exception("Something went wrong while running show-coords:" + show_coords_cmd))
    
    return show_coords_output_path

def process_congruence_table(coords_output_path, records, reference_length, congruence_threshold):
    congruence_table = process_coords_output(coords_output_path)
    congruence_table["Differences"] = None
    congruence_table["Flanks"] = None
    for record in records:
        calculate_further_columns_for_record(congruence_table, record, reference_length)

    congruence_table["DifferenceRate"] = congruence_table["Differences"]/reference_length
    congruence_table["Score"] = 1 - (congruence_table["Differences"] + congruence_table["Flanks"])/reference_length
    congruence_table["Accepted"] = congruence_table["Score"] > congruence_threshold
    
    return congruence_table

def process_coords_output(path_to_coords_file):
    df = pd.read_csv(path_to_coords_file, skiprows=5, sep='\s\|\s', engine='python', names=['S1 E1', 'S2 E2', 'LEN1 LEN2', '% IDY', 'TAGS'])

    df[['S1', 'E1']] = df['S1 E1'].apply(lambda x: x.strip()).str.split('\s', 1, expand=True)
    df[['S2', 'E2']] = df['S2 E2'].apply(lambda x: x.strip()).str.split('\s', 1, expand=True)
    df[['LEN1', 'LEN2']] = df['LEN1 LEN2'].apply(lambda x: x.strip()).str.split('\s', 1, expand=True)
    df[['REFERENCE', 'SUBJECT']] = df['TAGS'].apply(lambda x: x.strip()).str.split('\t', 1, expand=True)
    df = df[['S1', 'E1', 'S2', 'E2', 'LEN1', 'LEN2', '% IDY', 'REFERENCE', 'SUBJECT']]
    
    df['AlignedLength'] = df['E1'].astype(int) - df['S1'].astype(int)
    df['CongurentCharacters'] = round((df['E1'].astype(int) - df['S1'].astype(int))*df['% IDY'].astype(int)/100)

    df_grouped = df[['SUBJECT', 'AlignedLength', 'CongurentCharacters']].groupby('SUBJECT').sum()
    return df_grouped

def calculate_further_columns_for_record(congruence_table, record, reference_length):
    if record.name in congruence_table.index:
        congruence_table.loc[record.name, "Differences"] = reference_length - congruence_table.loc[record.name, 'CongurentCharacters']
        congruence_table.loc[record.name, "Flanks"] = 0 if reference_length > len(record.seq) else len(record.seq)-reference_length

def select_valid_records(congruence_table, records):
    congruence_table_valid = congruence_table.loc[congruence_table["Accepted"]]
    valid_ids = congruence_table_valid.index.to_list()

    valid_records = []
    for record in records:
        if record.name in valid_ids:
            valid_records.append(record)
    return valid_records