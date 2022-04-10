import os
from os.path import join
import sys

import numpy as np
import pandas as pd

import Config
from InputQC import *
from CongruenceQC import *
from MutationCalling import *
from MutationVis import *
from File_handling import *
import CovidMLApply
from CovidML_DL import *
import datetime
import SequenceQC
import shutil


try:
    BASEPATH = os.environ['BASEPATH']
    COVIDCONFIG = os.environ['COVIDCONFIG']
    #print(BASEPATH, COVIDCONFIG)
except KeyError as e:
    print('One of the envirement variable is not avaialbe: %s' % str(e))

CONFIG = Config.Config(COVIDCONFIG, BASEPATH).get()
sessions = {}


def set_session_data_null(session_id, raw_sequence_data):
    print('resetting the session data!')
    ct = datetime.datetime.now()
    new_id = session_id
    output_folder = join(CONFIG['OutputFolders']['output_folder'], new_id)
    sessions[session_id] = {'create_date': ct,
                            'output_folder': output_folder,
                            'raw_sequence_data': raw_sequence_data,
                            'raw_sequence_file': join(output_folder, new_id + '.fasta'),
                            'input_sequence': None,
                            'input_sequence_qc_table': None,
                            'input_sequence_qc': None,
                            'sample_mutations': None,
                            'sample_ml_results': None,
                            'applied_ml_model': None,
                            'input_sequence_status': None,  # If it is not OK, then do not procced with the analysis
                            'mutation_calculation_status': None,  # If it is ok, then do not proceed with ML
                            'ml_status': None,
                            'age_status': None,
                            'mutation_fig_status': None,
                            'mutations_anal_fig_status': None,
                            'ml_results_summary_fig_status': None,
                            'prediction_score_status': None,
                            'age_template.xlsx': None,
                            'qc_results_file': join(output_folder, new_id + '.qc_results.tsv'),
                            'input_cleaned_fasta_path': join(output_folder, new_id + '.filtered.fasta'),
                            'outout_congr_cleaned_fasta_path': join(output_folder, new_id + '.qc.fasta'),
                            'output_congr_table_path': join(output_folder, new_id + ".congruence.tsv"),
                            'output_alignment_path': join(output_folder, new_id + ".ma.fasta"),
                            'raw_mutations_output_file': join(output_folder, new_id + ".mutations.tsv"),
                            'agg_raw_mutations_output_file': join(output_folder, new_id + ".agg_mutations.tsv"),
                            'mutation_vcf_file': join(output_folder, new_id + ".mutations.vcf"),
                            'mutation_annotated_vcf_file': join(output_folder, new_id + ".mutations.ann.vcf"),
                            'mutations_samples': join(output_folder, new_id + ".mutations.filt.ann.tsv"),
                            'mutations_pic_file': join(output_folder, new_id + ".mutations.pic.svg"),
                            'mutations_anal_pic_file': join(output_folder, new_id + ".mutations.pic.svg"),
                            'protein_bar_plot_file': join(output_folder, new_id + ".mutations.protein_bar_plot.svg"),
                            'mutation_histogram_file': join(output_folder,
                                                            new_id + ".mutations.mutation_histogram_plot.svg"),
                            'ml_results_summary_file': join(output_folder, new_id + ".ml_results_summary_plot.svg"),
                            'zip_file_path': join(output_folder, new_id + ".compressed.zip"),
                            'ml_results_file': join(output_folder, new_id + ".ml_results.tsv"),
                            'status': {'created': {'timestamp': ct, 'status': 'OK'},
                                       'qc': {},  # keeps track on the internal states,
                                       'adding_age': {},
                                       'mutation_calling': {'status': None,
                                                            'calculation_start': None,
                                                            'end_time': None
                                                            }}
                            }

    # print (deleting the files)
    print('Deleting the content of the folder:')
    for filename in os.listdir(output_folder):
        file_path = os.path.join(output_folder, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except OSError as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))


def add_new_session(session_id, raw_sequence_data=None):
    print('Session id:   %s' % (session_id))
    print('Creating output directory!')
    new_id = str(session_id)
    output_folder = join(CONFIG['OutputFolders']['output_folder'], new_id)
    os.makedirs(output_folder, exist_ok=True)
    set_session_data_null(new_id, raw_sequence_data)
    sessions[new_id]['status']['init'] = {'timestamp': datetime.datetime.now(), 'status': 'OK'}
    return new_id


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
    sequence_info['age'] = np.NaN
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
        message = 'Too many age had been provided. The expected number of patent age is {0}, here we received: {1}'.format(
            N, len(ages_year))
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
    return sequence_info[['user_seq_id', 'sequence_id', 'age']]


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
    sequence_qc.to_csv(sessions[input_id]['qc_results_file'], index=False, sep='\t')
    check_input_sequence_set_status(input_id, 'OK')
    sessions[input_id]['input_sequence_status'] = 'OK'

    return sequence_qc


def get_seqqc_number_of_valid_sequences(input_id):
    '''
    Checks how many valid sequence record we have in the input qc. It assumes that there is QC record
    '''
    print('Checking the number of valid sequences after QC for batch: {0}. It should be non-zero'.format(input_id))

    try:
        seq_stat = sessions[input_id]['input_sequence_qc_table']
        N = len(seq_stat[seq_stat.qc_summary == 'OK'])
    except KeyError as ke:
        print('Key error. There is no such input id in this database. ID={0}'.format(input_id))
        N=0
    except AttributeError as ae:
        print('Attribute error. There is no qc_summary in the QC table. It is probably empty, abort. ID={0}'.format(input_id))
        N = 0

    print('Number of valid sequences after QC for batch: {0}; N={1}. It should be non-zero'.format(input_id,N ))
    return N


def mutation_calling_procedure_check_prereq(input_id):
    print('Checking prereqestics for the mutation calling')
    print('input_sequence_status:   ', sessions[input_id]['input_sequence_status'])
    if sessions[input_id]['input_sequence_status']:
        return True
    else:
        raise ValueError('Sequence input error.')


def get_mutation_call_status(input_id):
    """
    Query the actual mutation status
    """
    ct = datetime.datetime.now()
    mutation_calling_start_time = sessions[input_id]['status']['mutation_calling']['calculation_start']
    mutation_calling_status_flag = sessions[input_id]['status']['mutation_calling']['status']
    status_calc_var = 0
    if sessions[input_id]['status']['mutation_calling']['status'] and sessions[input_id]['status']['mutation_calling'][
        'end_time'] is not None \
            and mutation_calling_start_time is not None:
        running_time = (sessions[input_id]['status']['mutation_calling'][
                            'end_time'] - mutation_calling_start_time).total_seconds()
    elif mutation_calling_start_time is not None:
        running_time = (ct - mutation_calling_start_time).total_seconds()
    else:
        running_time = 0

    if mutation_calling_status_flag == 'OK':
        status_calc_var = 100
    elif mutation_calling_status_flag == 'MAFFT':
        status_calc_var = 5
    elif mutation_calling_status_flag == 'DNAMUTCALL':
        status_calc_var = 20
    elif mutation_calling_status_flag == 'DNAMUTCALLCLEAN':
        status_calc_var = 30
    elif mutation_calling_status_flag == 'DNAVCF':
        status_calc_var = 32
    elif mutation_calling_status_flag == 'SNPEFF':
        status_calc_var = 40
    elif mutation_calling_status_flag == 'ANNOT':
        status_calc_var = 90
    elif mutation_calling_status_flag == 'ERROR':
        status_calc_var = 0
    return mutation_calling_status_flag, status_calc_var, running_time




def mutation_calling_procedure(input_id):
    print('Identifing mutations for session data: {0}'.format(input_id))
    ct = datetime.datetime.now()

    sessions[input_id]['status']['mutation_calling']['calculation_start'] = ct
    sessions[input_id]['status']['mutation_calling']['status'] = 'INIT'

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

    sessions[input_id]['status']['mutation_calling']['status'] = 'MAFFT'
    run_mafft(mafft_alignment_input_fasta, output_alignment_path, CONFIG)
    sessions[input_id]['status']['mutation_calling']['status'] = 'DNAMUTCALL'
    mutation_calling(output_alignment_path, raw_mutations_output_file, agg_raw_mutations_output_file, CONFIG)
    raw_agg_mutations = pd.read_csv(agg_raw_mutations_output_file, sep='\t')
    sessions[input_id]['status']['mutation_calling']['status'] = 'DNAMUTCALLCLEAN'
    raw_agg_mutations_cleaned, mutations = filter_and_preprocess_aggr_mutation_table(raw_agg_mutations)
    sessions[input_id]['status']['mutation_calling']['status'] = 'DNAVCF'
    vcf_out = toVCF(mutations)
    vcf_out.to_csv(mutation_vcf_file, sep='\t', index=False)
    sessions[input_id]['status']['mutation_calling']['status'] = 'SNPEFF'
    annotate_mutations_with_vcf(snp_eff_path, snp_eff_reference_genome_id, mutation_vcf_file,
                                mutation_annotated_vcf_file)
    sessions[input_id]['status']['mutation_calling']['status'] = 'ANNOT'
    annotation_snpeff = pd.read_csv(mutation_annotated_vcf_file, sep='\t', header=5)
    annotations = get_annotation_df(annotation_snpeff)
    genome_annot, genome_annot_proteins = get_protein_annotation_ncbi(annotation_file)
    mutation_annotations = mutation_annotation_preprocessing(annotations, genome_annot_proteins)
    raw_agg_mutations_cleaned_annot = raw_agg_mutations_cleaned.merge(mutation_annotations, how='left',
                                                                      right_on='VCF_MUTATION_ID', left_on='Mutation')
    raw_agg_mutations_cleaned_annot.to_csv(mutations_samples, sep='\t', index=False)
    sessions[input_id]['sample_mutations'] = raw_agg_mutations_cleaned_annot
    sessions[input_id]['mutation_calculation_status'] = 'OK'
    sessions[input_id]['status']['mutation_calling']['status'] = 'OK'
    end_time = datetime.datetime.now()
    sessions[input_id]['status']['mutation_calling']['end_time'] = end_time
    return raw_agg_mutations_cleaned_annot


def evaluate_ml_model_check_prerequisite(input_id):
    '''
    We can only do prediction if we have properly formatted input data.
    '''
    import time
    print('Checking input for the analysis.')

    timout_limit = CONFIG['MutationCalling']['tim_out_sec']
    print('Time out limit: ' + str(timout_limit))

    ml_prereq_qc_status = False
    ml_prereq_mutation_status = False
    ml_prereq_mutation_results = False

    print('Checking QC call status.')
    if sessions[input_id]['input_sequence_status'] == 'OK':
        print('The input query has been processed successfuly.')
        ml_prereq_qc_status = True
    else:
        print('The input query has some preprocessing issue. Raising error')
        raise ValueError(
            "Error during the prediction. There are quality control problems with the input sequences. Please check the QC results first. ")
    #  sessions[input_id]['input_sequence_status']
    #
    print('Checking mutation call status: ')
    mutation_calling_status_flag, status_calc_var, running_time = get_mutation_call_status(input_id)
    # If the mutation calculation is not ready. let's just wait.
    if (not sessions[input_id]['mutation_calculation_status'] == 'ERROR') and sessions[input_id][
        'mutation_calculation_status'] != 'OK' and running_time <= timout_limit:
        print('The mutation calling is still running lets wait ...')
        for i in range(timout_limit * 10):
            print('The mutation calling is still running lets wait ...')
            time.sleep(0.1)
            mutation_calling_status_flag, status_calc_var, running_time = get_mutation_call_status(input_id)

            if running_time > timout_limit:
                break
            elif mutation_calling_status_flag == 'ERROR':
                break
            elif mutation_calling_status_flag == 'OK':
                break

    mutation_calling_status_flag, status_calc_var, running_time = get_mutation_call_status(input_id)
    if sessions[input_id]['mutation_calculation_status'] == 'OK':
        print('The mutation calculation status is OK.')
        ml_prereq_mutation_status = True
    else:
        print('The mutation calculation went wrong. Raising error')
        raise ValueError(
            "Error during the prediction. Mutation idenfitication went wrong, please check inputs and logs.")

    print('Checking number of mutations - for empty input')
    if len(sessions[input_id]['sample_mutations']) > 0:
        ml_prereq_mutation_results = True

    print('Checking the available age columns')
    sequence_metadata = sessions[input_id]['input_sequence_qc_table']
    if ('age' not in list(sequence_metadata.columns)):
        print('Adding age columns')
        sequence_metadata['age'] = np.NaN

    ml_prereq = ml_prereq_qc_status and ml_prereq_mutation_status and ml_prereq_mutation_results

    return [ml_prereq, ml_prereq_qc_status, ml_prereq_mutation_status, ml_prereq_mutation_results]


def writing_results_into_db(input_id, table_name, data):
    print('Writing the results into a database!')


def evaluate_ml_model_build_mutation_feature_table_annotated(input_id):
    print('Building features for mutation data')

    sequence_metadata = sessions[input_id]['input_sequence_qc_table']
    valid_sequences_with_age = sequence_metadata[(sequence_metadata.age > 0) &
                                                 (sequence_metadata.qc_summary == 'OK')][['sequence_id', 'age']].copy()
    valid_sequences_without_age = \
        sequence_metadata[(sequence_metadata.age.isna()) & (sequence_metadata.qc_summary == 'OK')][
            ['sequence_id']].copy()

    return [valid_sequences_with_age, valid_sequences_without_age]


def evaluate_ml_model_run_ml(input_id, model_type, input_mut_data, ml_model_type='automl'):
    # print('input_mut_data')
    # print(input_mut_data)

    base_path = CONFIG['BasePath']
    samples_ids = input_mut_data['sequence_id']

    input_mut_data.rename({'sequence_id': 'Accession ID'}, axis=1, inplace=True)
    input_mut_data.rename({'age': 'Patient age'}, axis=1, inplace=True, errors='ignore')
    # print(input_mut_data)

    model_dir = CONFIG['MachineLearning'][model_type]['model_dir']
    template_file = join(base_path, model_dir, CONFIG['MachineLearning'][model_type]['flat_table_template'])

    output_file = join(sessions[input_id]['output_folder'], 'ml_flat_table_' + model_type + '.csv')
    predicted_data_file = join(sessions[input_id]['output_folder'], 'ml_precitions_' + model_type + '.csv')

    headers = open(template_file).read().strip().split(',')
    template = pd.DataFrame([['template_00', col, 0] for col in headers[1:]],
                            columns=[headers[0], 'mutation_id', 'mutstatus'])

    sequence_mutations = sessions[input_id]['sample_mutations']
    sample_mutations = sequence_mutations[sequence_mutations.gisaid_epi_isl.isin(samples_ids)][
        ['gisaid_epi_isl', 'mutation_name']]
    sample_mutations['mutstatus'] = 1
    sample_mutations.columns = list(template.columns)

    sample_mutations.rename({'Accession ID': 'Accession ID'}, axis=1, inplace=True)
    pred_matrix_piv = pd.concat([sample_mutations, template]).pivot(index=headers[0], columns='mutation_id',
                                                                    values='mutstatus').fillna(0)
    pred_matrix_piv.drop(labels=['Patient age'], inplace=True, errors='ignore', axis=1)

    # print('Checking the input data')
    # print(input_mut_data)
    pred_matrix = input_mut_data.merge(pred_matrix_piv, how='inner', left_on='Accession ID', right_index=True)[headers]

    mut_col_list = list(set(pred_matrix.columns).difference(set(['Accession ID', 'Cohort', 'Patient age'])))
    type_list = [int for _ in range(len(mut_col_list))]
    zip_iterator = zip(mut_col_list, type_list)
    type_dictionary = dict(zip_iterator)
    pred_matrix = pred_matrix.astype(type_dictionary)
    print('Precition matrix after the conversion')
    #print(pred_matrix)
    pred_matrix.to_csv(output_file, sep=',', index=False)

    if ml_model_type == 'automl':
        jadbio_cmd = '''{0} -m {1} \
        -i {2} \
        -o {3} \
        '''.format(CONFIG['MachineLearning']['jadbio_executor'],
                   join(BASEPATH, model_dir, CONFIG['MachineLearning'][model_type]['model_file']),
                   output_file,
                   predicted_data_file)
        print(jadbio_cmd)
        os.system(jadbio_cmd)
        ml_results = pd.read_csv(predicted_data_file)
    elif ml_model_type == 'deep':
        print('Evaluating the deep learning model model')
        print('Using CPU as primary model')
        print('Using the BCEWithLogitsLoss as criterion')
        UseCuda = False
        criterion = torch.nn.BCEWithLogitsLoss()

        if model_type == 'deep_model_without_age':
            print('Appling the deep learning model withouth age')
            pred_matrix.iloc[:, 1:] = pred_matrix.iloc[:, 1:] * 2 - 1
            model_without_age_file = join(BASEPATH, model_dir, CONFIG['MachineLearning'][model_type]['model_file'])
            device = torch.device("cuda" if UseCuda else "cpu")
            model = C_FCClassifier().to(device)
            model.load_state_dict(torch.load(model_without_age_file, map_location='cpu'))
        elif model_type == 'deep_model_with_age':
            print('Appling the deep learning model with age')
            pred_matrix.iloc[:, 2:] = pred_matrix.iloc[:, 2:] * 2 - 1
            model_with_age_file = join(BASEPATH, model_dir, CONFIG['MachineLearning'][model_type]['model_file'])
            device = torch.device("cuda" if UseCuda else "cpu")
            model = CA_FCClassifier().to(device)
            model.load_state_dict(torch.load(model_with_age_file, map_location='cpu'))

        #print('The created prediction matrix:')
        #print(pred_matrix)

        with torch.no_grad():
            TensorX = torch.Tensor(pred_matrix.iloc[:, 1:].to_numpy())
            model.eval()
            prediction = torch.sigmoid(model(TensorX))
            prediction = prediction.cpu().detach().numpy()
            # Creating a new ml_results dataframe
            #print(pred_matrix_piv)
            ml_results = pred_matrix[['Accession ID']].copy()
            # print(ml_results)
            # print(prediction)
            # print(prediction.shape)
            # print(ml_results.info())
            # print(1-prediction)
            ml_results['P=Mild'] = 1 - prediction
            ml_results['P=Severe'] = prediction

    # print(ml_results)
    return ml_results


def evaluate_ml_model_get_confidence(score):
    confidence_flag = 'Undefined'

    if score > 0.8:
        confidence_flag = 'High'
    elif score > 0.6:
        confidence_flag = 'Low'
    elif score > 0.4:
        confidence_flag = 'Undefined'
    elif score > 0.2:
        confidence_flag = 'Low'
    elif score > 0.0:
        confidence_flag = 'High'

    return confidence_flag


def evaluate_ml_model_get_label(score):
    label = 'Undefined'
    if score > 0.5:
        label = 'Severe'
    else:
        label = 'Mild'
    return label


def evaluate_ml_model(input_id, ml_model_type='automl'):
    import random
    print('Evaluating the defined ML model for the uploaded data')
    ml_cols = ['sequence_id', 'prediction_score_mild', 'prediction_score', 'applied_model_type']
    [ml_prereq, ml_prereq_qc_status, ml_prereq_mutation_status,
     ml_prereq_mutation_results] = evaluate_ml_model_check_prerequisite(input_id)
    feature_table_for_sequences = evaluate_ml_model_build_mutation_feature_table_annotated(input_id)
    [valid_sequences_with_age, valid_sequences_without_age] = evaluate_ml_model_build_mutation_feature_table_annotated(
        input_id)

    print('Set analysis type')
    sessions[input_id]['applied_ml_model'] = ml_model_type

    if len(valid_sequences_with_age) > 0:
        print('Calculating the scores for the covid sequences with age!')
        if ml_model_type == 'automl':
            jadbio_results_age = evaluate_ml_model_run_ml(input_id, 'model_with_age', valid_sequences_with_age.copy(),
                                                          ml_model_type)
            jadbio_results_age['applied_model_type'] = 'model_with_age'
        elif ml_model_type == 'deep':
            jadbio_results_age = evaluate_ml_model_run_ml(input_id, 'deep_model_with_age',
                                                          valid_sequences_with_age.copy(),
                                                          ml_model_type)
            jadbio_results_age['applied_model_type'] = 'deep_model_with_age'
            print('Working with the deep model! (with age)')
    else:
        jadbio_results_age = pd.DataFrame()
    if len(valid_sequences_without_age) > 0:
        print('Calculating the scores for the covid sequences withOUT age!')
        if ml_model_type == 'automl':
            jadbio_results_wage = evaluate_ml_model_run_ml(input_id, 'model_without_age',
                                                           valid_sequences_without_age.copy(), ml_model_type)
            jadbio_results_wage['applied_model_type'] = 'model_without_age'
        elif ml_model_type == 'deep':
            print('Working with the deep model!')
            jadbio_results_wage = evaluate_ml_model_run_ml(input_id, 'deep_model_without_age',
                                                           valid_sequences_without_age.copy(), ml_model_type)
            jadbio_results_wage['applied_model_type'] = 'deep_model_without_age'


    else:
        jadbio_results_wage = pd.DataFrame()
    jadbio_results = pd.concat([jadbio_results_age, jadbio_results_wage])
    # print('jadbio_results')
    # print(jadbio_results)
    jadbio_results.columns = ml_cols
    jadbio_results['prediction_confidence'] = jadbio_results.apply(
        lambda x: evaluate_ml_model_get_confidence(x['prediction_score']), axis=1)
    jadbio_results['Label'] = jadbio_results.apply(lambda x: evaluate_ml_model_get_label(x['prediction_score']), axis=1)

    # print(jadbio_results)
    # merging in the prediction results
    current_cols = list(sessions[input_id]['input_sequence_qc_table'].columns)
    diff_cols = set(current_cols) - set(ml_cols[1:])
    keep_cols = [act_colname for act_colname in current_cols if act_colname in diff_cols]
    # print('keep_cols  ', str(keep_cols))

    act_seq_data = sessions[input_id]['input_sequence_qc_table'][diff_cols]

    # print(act_seq_data)
    # print(act_seq_data.columns)
    # print(jadbio_results)
    # print(jadbio_results.columns)

    ml_results = act_seq_data.merge(jadbio_results, how='left', left_on='sequence_id', right_on='sequence_id')
    sessions[input_id]['sample_ml_results'] = jadbio_results
    # print('______________')
    # print(jadbio_results)
    # print('______________')
    # print(ml_results)
    # print('Writning the input qs results')
    # sessions[input_id]['input_sequence_qc_table'] = ml_results
    # print('ok')
    sessions[input_id]['ml_status'] = 'OK'
    # print('sffsddf')
    filtered_ml_results = ml_results[~ml_results['prediction_score'].isnull()]
    # print(filtered_ml_results)
    return filtered_ml_results


def get_ml_results(input_id):
    print('Query the ml results only')
    if sessions[input_id]['ml_status'] == 'OK':
        raw_ml_results = sessions[input_id]['sample_ml_results']
        act_seq_data = sessions[input_id]['input_sequence_qc_table']
        ml_results = act_seq_data.merge(raw_ml_results, how='inner', left_on='sequence_id', right_on='sequence_id')
        # print(ml_results)
    else:
        print('No available ML calculation')
        ml_results = None
        raise ValueError()
    return ml_results


def get_suggestion_for_string(query_mut, max_hits=10):
    # Let's iterate throug the list
    best_hits = []
    query_mut_mod = query_mut.lower()

    for record in CONFIG['Database']['mutation_db_suggestions'].itertuples():
        try:
            if query_mut_mod in record[3]:
                best_hits.append(record[1])
        except TypeError as te:
            print(te)
        if len(best_hits) > max_hits:
            break
    return best_hits


def evaluate_ml_model_generate_distribution(input_id):
    print('Generating background distribution!')


def generate_age_template_from_qc_table(input_id):
    print('Generating an excel file with the sequence ids for the user.')
    if sessions[input_id]['input_sequence_status']:
        pass
    else:
        print('Some error!')


def get_static_gene_picture(input_id):
    print('Getting a static gene picture')
    print('Here it comes the figure creating scripts')
    print('Chenking ML prereqeistics')

    if sessions[input_id]['mutation_fig_status'] == 'OK' and \
            sessions[input_id]['mutations_anal_fig_status'] == 'OK' and \
            sessions[input_id]['prediction_score_status'] == 'OK':
        return None

    [ml_prereq, ml_prereq_qc_status, ml_prereq_mutation_status,
     ml_prereq_mutation_results] = evaluate_ml_model_check_prerequisite(input_id)
    base_path = CONFIG['BasePath']
    sessions[input_id]['mutations_pic_file'] = join(base_path, 'data/pictures/mutation_hist.svg')
    sessions[input_id]['mutations_anal_pic_file'] = join(base_path, 'data/pictures/genes.png')

    all_sample_data = get_mutation_data(sessions[input_id]['sample_mutations'])
    prediction_table = get_ml_results(input_id)
    modified_prediction_table = get_prediction_table(prediction_table)

    protein_bars = make_protein_bar_plot(all_sample_data, sessions[input_id]['protein_bar_plot_file'])
    histograms = make_mutation_histogram(all_sample_data, sessions[input_id]['mutation_histogram_file'])
    ml_plot = make_prediction_histogram(modified_prediction_table, sessions[input_id]['ml_results_summary_file'])

    sessions[input_id]['mutation_fig_status'] = 'OK'
    sessions[input_id]['mutations_anal_fig_status'] = 'OK'
    sessions[input_id]['ml_results_summary_fig_status'] = 'OK'


def test_stuff_static_pictures(input_id):
    all_sample_data = get_mutation_data(sessions[input_id]['sample_mutations'])

    protein_bars = make_protein_bar_plot(all_sample_data, sessions[input_id]['protein_bar_plot_file'])
    histograms = make_mutation_histogram(all_sample_data, sessions[input_id]['mutation_histogram_file'])

    print(all_sample_data)


def build_artificial_genome_data(session_id, mutation_lists):
    '''
    Building the artificial genome representation for the data we have
    Steps:
      - generating IDS,
      - transforming the received mutation data
      - mapping to annototation
      - setting the proper flags
    '''
    print('Creating database record!')
    create_artifical_genomes_session(session_id)
    # Check for the valid mutations and some report would be needed ...
    updated_normalized_mutlists, empty_genomes_row_id, unmapped_mutations = build_artificial_genome_data_check_and_map_the_mutation_names(
        mutation_lists)
    # Generating IDs and QC table
    print('Generating QC table and IDs for the records!')
    qc_results, id2mutations = get_qc_table_for_artifical_genomes(session_id, updated_normalized_mutlists)
    print('Building the mutation data tables!')
    create_artifical_genomes_mutation(session_id, id2mutations)

    return updated_normalized_mutlists, empty_genomes_row_id, unmapped_mutations


def build_artificial_genome_data_check_and_map_the_mutation_names(mutation_lists):
    ''' Checking the mutation names in the uploaded '''

    print('Checking the uploaded mutations strings')
    updated_normalized_mutlists = []
    empty_genomes_row_id = []
    unmapped_mutations = []
    for i, mutlists in enumerate(mutation_lists):
        # print('Genome {0}'.format(i))
        # print(' mutations {0}'.format(mutlists))
        mapped_list = []
        unmapped_list = []
        for act_mutation in mutlists:
            act_mutation_q = act_mutation.strip().lower()
            if act_mutation_q in CONFIG['Database']['mutation_name_mapping']:
                mapped_list.append(CONFIG['Database']['mutation_name_mapping'][act_mutation_q])
            else:
                unmapped_list.append(act_mutation)
        mapped_list = list(set(mapped_list))
        if len(mapped_list) > 0:
            updated_normalized_mutlists.append(mapped_list)
        else:
            empty_genomes_row_id.append(i)
        unmapped_mutations.extend(unmapped_list)
        # print(mapped_list)
        # print(unmapped_list)
    unmapped_mutations = list(set(unmapped_mutations))
    print('Unmapped mutations {0}'.format(unmapped_mutations))
    print('Empty genomes mutations {0}'.format(empty_genomes_row_id))

    return updated_normalized_mutlists, empty_genomes_row_id, unmapped_mutations


def create_artifical_genomes_mutation(session_id, id2mutations):
    ct = datetime.datetime.now()
    sessions[session_id]['status']['mutation_calling']['calculation_start'] = ct
    sessions[session_id]['status']['mutation_calling']['status'] = 'INIT'
    output_folder = join(CONFIG['OutputFolders']['output_folder'], session_id)
    mutations_samples_file = sessions[session_id]['mutations_samples']
    mutation_data = artifical_genomes_get_id2mutation_map(id2mutations)
    filtered_mutations_dedup = CONFIG['Database']['mutation_annotation']
    annotaded_mutations = mutation_data.merge(filtered_mutations_dedup, how='inner', left_on='mutation_id',
                                              right_on='mutation_name')
    sessions[session_id]['sample_mutations'] = annotaded_mutations
    annotaded_mutations.to_csv(mutations_samples_file, sep='\t', index=False)
    sessions[session_id]['mutation_calculation_status'] = 'OK'
    sessions[session_id]['status']['mutation_calling']['status'] = 'OK'
    sessions[session_id]['status']['mutation_calling']['end_time'] = datetime.datetime.now()


def create_artifical_genomes_session(session_id):
    set_session_data_null(session_id, None)
    act_session_data = sessions[session_id]
    act_session_data['raw_sequence_data'] = None
    check_input_sequence_set_status(session_id, 'OK')
    sessions[session_id]['input_sequence_status'] = 'OK'


def artifical_genomes_get_id2mutation_map(id2mutations):
    seq_id_mutations_list = []
    for seq_id in id2mutations.keys():
        for mutation in id2mutations[seq_id]:
            seq_id_mutations_list.append([seq_id, mutation])
    mutation_data = pd.DataFrame(seq_id_mutations_list,
                                 columns=['gisaid_epi_isl', 'mutation_id'])
    return mutation_data


def get_qc_table_for_artifical_genomes(session_id, mutation_lists):
    N_seqs = len(mutation_lists)
    new_qc_data_list = []
    id2mutations = {}
    for i in range(N_seqs):
        new_id = '{0}_{1}'.format(session_id, str(i).rjust(3, '0'))
        new_qc_data_list.append([new_id, new_id, 'OK', 'OK'])
        id2mutations[new_id] = mutation_lists[i]
    if N_seqs > 0:
        qc_results = pd.DataFrame(new_qc_data_list)
        qc_results.columns = ['user_seq_id', 'sequence_id', 'check_ids', 'qc_summary']
        sessions[session_id]['input_sequence_qc_table'] = qc_results
    else:
        sessions[session_id]['input_sequence_status'] = None
        qc_results = None

    return qc_results, id2mutations


def get_all_zipped_files(session_id):
    print('Start zipping the available files')
    output_folder = sessions[session_id]['output_folder']
    expected_zip_file = sessions[session_id]['zip_file_path']
    readme_file = join(CONFIG['BasePath'], 'data/Readme')
    ml_results_file = sessions[session_id]['ml_results_file']
    # Creating Ml results if not exists

    try:
        if not (os.path.exists(ml_results_file) and os.path.getsize(ml_results_file) > 0):
            ml_results = get_ml_results(session_id)
            ml_results.to_csv(ml_results_file, sep='\t', index=False)
    except IOError:
        print('IOError during the creation of ML results file')
    except ValueError:
        print('Value error during the process of dumping the ML file')
    files_to_compress = [readme_file,
                         sessions[session_id]['qc_results_file'],
                         sessions[session_id]['outout_congr_cleaned_fasta_path'],
                         sessions[session_id]['mutation_annotated_vcf_file'],
                         sessions[session_id]['mutations_samples'],
                         ml_results_file
                         ]
    print(files_to_compress)

    try:
        if os.path.exists(expected_zip_file) and os.path.getsize(expected_zip_file) > 0:
            print('The zipped file exists and it is not null! Rewrting it!')
            valid_files = check_uploadable_files(files_to_compress)
            generate_zipped_files(valid_files, expected_zip_file)
        else:
            print('Check ...')
            valid_files = check_uploadable_files(files_to_compress)
            generate_zipped_files(valid_files, expected_zip_file)
    except IOError as ioerror:
        print('Some IOError happened')

    print(expected_zip_file)

    return expected_zip_file
