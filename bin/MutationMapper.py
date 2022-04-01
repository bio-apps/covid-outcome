#!/usr/bin/env python
# -*-coding: utf8 -*-

"""
Library for multiple alignment and mutation calling in hCov-19 analysis.
Copyright 2021 by Regina Kalcsevszki (kalcsevszkiregi@gmail.com), Balázs Ligeti (obalasz@gmail.com). All rights reserved.

"""

import os
from os.path import join, isfile
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'libs'))

import covlib
from covlib import *
#from data_preprocessing import *

import InputQC

def map_mutations(agg_mutations_file, output_flat_table_path, CONFIG):
    mutation_mapper_path = join(CONFIG["BasePath"], CONFIG["MutationCalling"]["MutationMapperPathFromBase"])
    
    mapper = pd.read_csv(mutation_mapper_path, sep='\t')
    mapper = mapper.set_index("Mutation")
    
    mutations = pd.read_csv(agg_mutations_file, sep='\t')
    mutations = mutations.join(mapper, on="Mutation", how='inner')[["Target ID", "mutation_name"]]
    mutations["Value"] = True
    
    mutations_flat = mutations.pivot(index="Target ID", columns='mutation_name', values='Value')
    
    missing_mutations = pd.DataFrame(columns=set(mapper['mutation_name']).difference(set(mutations['mutation_name'])))
    mutations_flat = mutations_flat.join(missing_mutations).fillna(False)
    
    mutations_flat.to_csv(output_flat_table_path, sep="\t")
    
    return mutations_flat

# Metadata file should be a tsv file, with \t separator with at least two columns: 'Sample' and 'Patient age', where 'Sample ID' is the sample id in the original fasta file for each sequence, and 'Patient age' is the age of the patient!
def join_metadata(mutations_flat_path, metadata_path, output_flat_table_with_age_path, CONFIG):
    flat_table = pd.read_csv(mutations_flat_path, sep='\t')
    target_id_mapper = flat_table[["Target ID"]]
    target_id_mapper[['Description', 'Mapped ID', 'Date']] = target_id_mapper["Target ID"].str.split('|', 3, expand=True)

    metadata = pd.read_csv(metadata_path, sep='\t')
    metadata["Description"] = metadata["Sample"].apply(InputQC.change_spec_characters_in_id)
    
    age_mapper = target_id_mapper.set_index('Description').join(metadata.set_index("Description")).set_index("Target ID")[['Patient age']]

    flat_table_with_age = flat_table.set_index('Target ID').join(age_mapper)
    
    flat_table_with_age.to_csv(output_flat_table_with_age_path)
    
    return flat_table_with_age
    