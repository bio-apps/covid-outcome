import os
from os.path import join
import sys

import numpy as np
import pandas as pd

import Config
from InputQC import *
from CongruenceQC import *
from MutationCalling import *
import CovidMLApply
import datetime
import SequenceQC

# expected to have the following enviremtn variables:
# BASEPATH #Ebben van az útvonal a projekthez
# COVIDCONFIG #Ebben van a config file, amiben a paramétereket tesszük


try:
    BASEPATH = os.environ['BASEPATH']
    COVIDCONFIG = os.environ['COVIDCONFIG']
    print(BASEPATH, COVIDCONFIG)
except KeyError as e:
    print('One of the envirement variable is not avaialbe: %s' % str(e))

CONFIG = Config.Config(COVIDCONFIG, BASEPATH).get()
sessions = {}


def add_new_session(session_id, raw_sequence_data=None):
    ct = datetime.datetime.now()
    new_int_id = len(sessions) + 1
    new_id = 'S' + f'{new_int_id:05}'
    new_id = str(session_id)

    output_folder = join(CONFIG['OutputFolders']['output_folder'], new_id)
    sessions[new_id] = {'create_date': ct,
                        'output_folder': output_folder,
                        'raw_sequence_data': raw_sequence_data,
                        'raw_sequence_file': join(output_folder, new_id + '.fasta'),
                        'input_sequence': None,
                        'input_sequence_qc_table' : None,
                        'input_sequence_qc': None,
                        'sample_mutations': None,
                        'input_sequence_status': None,  # If it is not OK, then do not procced with the analysis
                        'mutation_calculation_status': None, #If it is ok, then do not proceed with ML
                        'ml_status': None,
                        'age_status': None,
                        'prediction_score_status': None,
                        'age_template.xlsx': None,
                        'input_cleaned_fasta_path': join(output_folder, new_id + '.filtered.fasta'),
                        'outout_congr_cleaned_fasta_path': join(output_folder, new_id + '.qc.fasta'),
                        'output_congr_table_path': join(output_folder, new_id + ".congruence.tsv"),
                        'output_alignment_path': join(output_folder, new_id + ".ma.fasta"),
                        'raw_mutations_output_file': join(output_folder, new_id + ".mutations.tsv"),
                        'agg_raw_mutations_output_file': join(output_folder, new_id + ".agg_mutations.tsv"),
                        'mutation_vcf_file': join(output_folder, new_id + ".mutations.vcf"),
                        'mutation_annotated_vcf_file': join(output_folder, new_id + ".mutations.ann.vcf"),
                        'mutations_samples': join(output_folder, new_id + ".mutations.filt.ann.tsv"),


                        'status': {'created': {'timestamp': ct, 'status': 'OK'},
                                   'qc': {},  # keeps track on the internal states,
                                   'adding_age' : {},
                                   'mutation_calling': {} }
                        }

    return new_id


def initiate_processing_input_data(input_id):
    print('Session id:   %s' % input_id)
    print('Creating output directory!')
    output_folder = join(CONFIG['OutputFolders']['output_folder'], input_id)
    os.makedirs(output_folder, exist_ok=True)
    sessions[input_id]['status']['init'] = {'timestamp': datetime.datetime.now(), 'status': 'OK'}

def adding_age_information_check_pre_req(input_id):

    print('input_sequence_status:   ', sessions[input_id]['input_sequence_status'])
    if sessions[input_id]['input_sequence_status']:
        return True
    else:
        raise ValueError('Sequence input error.')


def adding_age_information(input_id, age_list):
    print('Adding age information!')
    # We should have valid sequences here.
    adding_age_information_check_pre_req(input_id)
    sequence_info = sessions[input_id]['input_sequence_qc_table']
    N = len(sequence_info)
    # Parsing the id-s 1 by 1
    ages_year = []
    for age in age_list:
        try:
            act_age = float(age)
            if act_age >= CONFIG['InputQC']['Age_min'] and act_age <= CONFIG['InputQC']['Age_max']:
                ages_year.append(act_age)
            else:
                print('Invalid age: ', str(act_age))
                ages_year.append(np.NaN)
        except ValueError as VE:
            print('Invalid age: ', str(age), 'Skipping ...')
            ages_year.append(np.NaN)
    if len(ages_year) > N:
        message = 'Too many age had been provided. The expected number of patent age is {0}, here we received: {1}'.format(N, len(ages_year))
        print(message)
        # Set failed flag
        raise ValueError(message)
    # Merging in the the age information
    sequence_info['age'] = np.NaN
    ext_ages = [np.NaN for i in range(N)]
    for i in range(len(ages_year)):
        ext_ages[i] = ages_year[i]
    sequence_info['age'] = ext_ages
    sessions[input_id]['age_status'] = 'OK'

    return sequence_info[['user_seq_id', 'sequence_id','age']]


def check_input_sequence_set_status(input_id, status_flag):
    sessions[input_id]['status']['qc']['status'] = status_flag


def check_input_sequences(input_id):
    print('Starting checking the input sequences!')
    sequences = SequenceQC.parse_sequence_from_string(sessions[input_id]['raw_sequence_data'])
    sessions[input_id]['input_sequence'] = sequences
    sequence_qc = pd.DataFrame()
    qc_fasta_path = sessions[input_id]['outout_congr_cleaned_fasta_path']

    if len(sequences) > CONFIG['MaxSampleSize']:
        check_input_sequence_set_status(input_id, 'Failed')
        raise ValueError('Too many sequences. Maximum number of supported sequences: %d.' % CONFIG['MaxSampleSize'])
    if len(sequences) == 0:
        check_input_sequence_set_status(input_id, 'Failed')
        raise ValueError('No valid sequence in the input data. Please check the input!')
    sequence_qc, N_fix = SequenceQC.quality_control_sequences(sequences, CONFIG, input_id)
    # Dumping the outputs, setting the flags etc...
    sessions[input_id]['input_sequence_qc_table'] = sequence_qc
    ok_seq_data = sequence_qc[sequence_qc.qc_summary == 'OK']
    if len(ok_seq_data) == 0:
        print('No valid sequence avalialbe for the analysis.')
        check_input_sequence_set_status(input_id, 'Failed')
        sessions[input_id]['input_sequence_status'] = False
        raise ValueError('No valid sequences are available')
    else:
        valid_sequence_ids = {seq_id: 0 for seq_id in list(ok_seq_data[CONFIG['DataSchema']['sequence_id']])}
        qc_sequences = [seq_record for seq_record in sequences if seq_record.id in valid_sequence_ids]
        sessions[input_id]['input_sequence_qc'] = qc_sequences
        with open(qc_fasta_path, 'w') as fout:
            for qc_sequence in qc_sequences:
                SeqIO.write(qc_sequence, fout, 'fasta')
    check_input_sequence_set_status(input_id, 'OK')
    sessions[input_id]['input_sequence_status'] = 'OK'

    return sequence_qc


def mutation_calling_procedure_check_prereq(input_id):
    print('Checking prereqestics for the mutation calling')
    print('input_sequence_status:   ', sessions[input_id]['input_sequence_status'])

    if sessions[input_id]['input_sequence_status']:
        return True
    else:
        raise ValueError('Sequence input error.')


def mutation_calling_procedure(input_id):
    print('Identifing mutations for session data: {0}'.format(input_id))
    mutation_calling_procedure_check_prereq(input_id)

    output_folder = join(CONFIG['OutputFolders']['output_folder'], input_id)
    snp_eff_path = CONFIG['MutationCalling']['snp_eff_command']
    snp_eff_reference_genome_id = CONFIG['MutationCalling']['snp_eff_reference_genome_id']
    annotation_file = join(CONFIG['BasePath'], CONFIG['Reference']['annotation'])

    mafft_alignment_input_fasta = sessions[input_id]['outout_congr_cleaned_fasta_path']

    output_alignment_path = sessions[input_id]['output_alignment_path']
    raw_mutations_output_file = sessions[input_id]['raw_mutations_output_file']

    agg_raw_mutations_output_file = sessions[input_id]['agg_raw_mutations_output_file']
    mutation_vcf_file = sessions[input_id]['mutation_vcf_file']
    mutation_annotated_vcf_file = sessions[input_id]['mutation_annotated_vcf_file']
    mutations_samples = sessions[input_id]['mutations_samples']

    run_mafft(mafft_alignment_input_fasta, output_alignment_path, CONFIG)
    mutation_calling(output_alignment_path, raw_mutations_output_file, agg_raw_mutations_output_file, CONFIG)
    raw_agg_mutations = pd.read_csv(agg_raw_mutations_output_file, sep='\t')
    raw_agg_mutations_cleaned, mutations = filter_and_preprocess_aggr_mutation_table(raw_agg_mutations)
    vcf_out = toVCF(mutations)
    vcf_out.to_csv(mutation_vcf_file, sep='\t', index=False)
    annotate_mutations_with_vcf(snp_eff_path, snp_eff_reference_genome_id, mutation_vcf_file,
                                mutation_annotated_vcf_file)
    annotation_snpeff = pd.read_csv(mutation_annotated_vcf_file, sep='\t', header=5)
    annotations = get_annotation_df(annotation_snpeff)
    genome_annot, genome_annot_proteins = get_protein_annotation_ncbi(annotation_file)
    mutation_annotations = mutation_annotation_preprocessing(annotations, genome_annot_proteins)
    raw_agg_mutations_cleaned_annot = raw_agg_mutations_cleaned.merge(mutation_annotations, how='left',
                                                                      right_on='VCF_MUTATION_ID', left_on='Mutation')

    raw_agg_mutations_cleaned_annot.to_csv(mutations_samples, sep='\t', index=False)
    sessions[input_id]['sample_mutations'] = raw_agg_mutations_cleaned_annot
    sessions[input_id]['mutation_calculation_status'] = 'OK'
    return raw_agg_mutations_cleaned_annot


def evaluate_ml_model_build_mutation_feature_table_annotated(input_id=None, mutation_table = None):
    print('Building features for mutation data')

    return None



def evaluate_ml_model_with_age(input_id):
    print('Evaluating the model with the provided features!')


def evaluate_ml_model_withoud_age(input_id):
    print('Evaluating ML model without age')



def evaluate_ml_model_check_prerequisite(input_id):
    print('Checking input for the analysis.')




def writing_results_into_db(input_id, table_name, data):
    print('Writing the results into a database!')


def evaluate_ml_model(input_id):
    import random
    print('Evaluating the defined ML model for the uploaded data')
    ml_results = sessions[input_id]['input_sequence_qc_table']
    ml_results['prediction_score'] = np.NaN
    ml_results.loc[:, 'prediction_score'] = ml_results.apply(lambda x: random.random(), axis=1)
    evaluate_ml_model_check_prerequisite(input_id)
    feature_table_for_sequences = evaluate_ml_model_build_mutation_feature_table_annotated(input_id)
    ml_results_age = ml_results[~ml_results.age.isnull()]
    print(ml_results_age)
    ml_results_without_age = ml_results[ml_results.age.isnull()]

    sessions[input_id]['ml_status'] = 'OK'
    print('bla')
    return ml_results


def evaluate_ml_model_generate_distribution(input_id):

    print('Generating background distribution!')



def run_pipeline(input_id, input_fasta, CONFIG, input_ages=None):
    '''
    Running the mutation calling and identification pipeline!
    :param input_id: input_id
    :param input_fasta:  input fasta to process
    :param input_ages:   inp
    :return:
    '''

    print('Running the pipeline!')
    print('Input fasta:   %s' % input_fasta)
    print('Input ages:    %s' % str(input_ages))

    print('Creating output directory!')
    output_folder = join(CONFIG['OutputFolders']['output_folder'], input_id)
    os.makedirs(output_folder, exist_ok=True)

    # Ide kellene egy QC kontroll a szekvenciákra
    # Ide kellene egy QC kontroll a korra

    # raw_agg_mutations_cleaned_annot = run_mutation_finder_pipeline(input_id, input_fasta, CONFIG)
    mutations_samples = join(output_folder, input_id + ".mutations.filt.ann.tsv")
    mutations_samples_file = join(output_folder, input_id + ".mutations.filt.ann.tsv")
    raw_agg_mutations_cleaned_annot = pd.read_csv(mutations_samples_file, sep='\t')
    raw_agg_mutations_cleaned_annot.rename({'gisaid_epi_isl': 'sequence_id'}, inplace=True, axis=1)
    prediction_scores = CovidMLApply.apply_ML_model(input_id, input_fasta, CONFIG, raw_agg_mutations_cleaned_annot,
                                                    input_ages)
    # print(prediction_scores)
    return raw_agg_mutations_cleaned_annot, prediction_scores
    # print(raw_agg_mutations_cleaned_annot)

def run_pipeline_call_sample(CONFIG, input_fasta=None, input_ages=None):
    '''

    :return:
    '''
    input_id = 'T0001'
    print('Initiating the pipeline on the query: %s' % input_id)
    if ~input_fasta:
        input_fasta = '/home/ligeti/gitrepos/hCov-19/bin/Pipeline/test/TestData/TestSequences.fasta'
    mutation_profile, prediction_scores = run_pipeline(input_id, input_fasta, CONFIG, input_ages)
    return mutation_profile, prediction_scores

# run_pipeline_call_sample(CONFIG, 'TEST0001')

def generate_age_template_from_qc_table(input_id):

    print('Generating an excel file with the sequence ids for the user.')
    if sessions[input_id]['input_sequence_status']:
        pass
    else:
        print('Some error!')

