'''
Python module for testing the basic quality of sequences
'''
import pandas as pd
from Bio import SeqIO
from io import StringIO
import CongruenceQC
from os.path import join
import sys
import traceback
import numpy as np


def parse_sequence_from_string(sequences):
    '''
    Parsing a string if it is looks like a fasta file
    :param CONFIG:
    :param sequences:
    :return: set of fasta records
    '''

    seqs = [seq for seq in SeqIO.parse(StringIO(sequences), 'fasta')]

    return seqs


def quality_control_sequences_check_ids(sequences, qc_records, main_cols, CONFIG, input_id):
    print('Checking the sequence IDS and generating new ones!')
    id_prefix = CONFIG['InputQC']['SampleIDPrefix']
    isSuccess = False
    for i, record in enumerate(qc_records):
        new_seq_id = id_prefix + '_' + str(input_id) + '_' + str(i)
        sequences[i].id = new_seq_id
        sequences[i].name = ''
        sequences[i].description = ''
        sequences[i].dbxrefs = []

        record.append(new_seq_id)
        record.append(1)
        isSuccess = True
    if isSuccess:
        main_cols.append('sequence_id')
        main_cols.append('check_ids')


def quality_control_sequences_check_seq_lengths(qc_records, main_cols, CONFIG):
    print('Checking input sequence sizes!')
    isSuccess = False
    min_length = CONFIG['InputQC']['Sequence_length']['MinSequenceLength']
    max_length = CONFIG['InputQC']['Sequence_length']['MaxSequenceLength']

    for i, record in enumerate(qc_records):
        if record[1] > max_length or record[1] < min_length:
            record.append(0)
        else:
            record.append(1)
        isSuccess = True

    if isSuccess:
        main_cols.append('check_sequence_length')


def quality_control_sequences_check_nucleotid_props(sequences, qc_records, main_cols, CONFIG):
    print('Checking GC content and ACTG content')
    isSuccess = False
    min_gc = CONFIG['InputQC']['GCContent']['GCmin']
    max_gc = CONFIG['InputQC']['GCContent']['GCmax']
    actg_prop = CONFIG['InputQC']['MinACGTRatio']

    for i, record in enumerate(qc_records):
        seq_length = record[1]
        A = sequences[i].seq.count('A')
        T = sequences[i].seq.count('T')
        C = sequences[i].seq.count('C')
        G = sequences[i].seq.count('G')
        valid_letter_count = A + T + C + G
        # Non valid nucleotides
        if seq_length > 0 and valid_letter_count / seq_length > actg_prop:
            record.append(1)
        else:
            record.append(0)
        # GC content
        if valid_letter_count > 0 and ( \
                        min_gc < (G + C) / (valid_letter_count) < max_gc):
            record.append(1)
        else:
            record.append(0)
        isSuccess = True
    if isSuccess:
        main_cols.append('check_ACTG_proportion')
        main_cols.append('check_GC_content')


def quality_control_sequences_check_uppercase(sequences, qc_records, main_cols):
    print('Converting sequence characters into all uppercase!')
    isSuccess = False
    for i, record in enumerate(qc_records):
        sequences[i].seq = sequences[i].seq.upper()
        record.append(1)
        isSuccess = True
    if isSuccess:
        main_cols.append('check_uppercase')


def quality_control_sequences_add_sequence_id(qc_records, main_cols):
    isSuccess = False
    for i, record in enumerate(qc_records):
        record.append(record[0])
        isSuccess = True

    if isSuccess:
        main_cols.append('sequence_id')


def quality_control_sequences_congruence_check(sequences, qc_records_df, CONFIG, N_fix, qc_steps,
                                               input_cleaned_fasta_path,
                                               output_folder,
                                               outout_congr_cleaned_fasta_path,
                                               output_congr_table_path):
    print('Here comes the congruence checking procedure!')
    N = len(qc_records_df.columns)
    qc_sum = qc_records_df.iloc[:, N_fix:].sum(axis=1)
    N_ok = len(qc_sum[qc_sum == N - N_fix])
    qc_checked = list(qc_sum)
    with open(input_cleaned_fasta_path, 'w') as fout:
        for i, qc_value in enumerate(qc_checked):
            if qc_value == N - N_fix:
                SeqIO.write(sequences[i], fout, 'fasta')
    if N_ok == 0:
        print('Empty input for congruence calculation.')
        qc_records_df['check_congruence'] = np.NaN
        return qc_records_df

    try:
        valid_records, congurence_table = CongruenceQC.run_congruence_qc(input_cleaned_fasta_path,
                                                                         output_folder,
                                                                         outout_congr_cleaned_fasta_path,
                                                                         output_congr_table_path, CONFIG)
        congurence_table.loc[congurence_table.Accepted, 'check_congruence'] = 1
        qc_records_df = qc_records_df.merge(congurence_table[congurence_table.Accepted][['check_congruence']],
                                            how='left',
                                            left_on='sequence_id', right_index=True).fillna(0)
    except ValueError as VE:
        print('Error occured during the congruence value calculation. Check the input.')
        qc_records_df['check_congruence'] = np.NaN
        print(sys.exc_info())
        print(traceback.format_exc())

    return qc_records_df


def quality_control_sequences(sequences, CONFIG, input_id):
    # set processing status?
    qc_steps = CONFIG['InputQC']['QCSteps']
    if qc_steps is None:
        print('No QC is required!')
        qc_steps = []
    else:
        print('Checking sequences!\nQC steps: ')
        for qc_step in qc_steps:
            print('  %s' % qc_step)

    main_cols = ['user_seq_id', 'sequence_length', 'user_fasta_header']
    N_fix = len(main_cols) + 1
    qc_records = [[seq_record.description, len(seq_record.seq), seq_record.format('fasta').split('\n')[0]]
                  for seq_record in sequences]

    output_folder = join(CONFIG['OutputFolders']['output_folder'], input_id)
    input_cleaned_fasta_path = join(output_folder, input_id + '.filtered.fasta')
    outout_congr_cleaned_fasta_path = join(output_folder, input_id + ".qc.fasta")
    output_congr_table_path = join(output_folder, input_id + ".congruence.tsv")

    if 'check_ids' in qc_steps:
        quality_control_sequences_check_ids(sequences, qc_records, main_cols, CONFIG, input_id)
    else:
        # Adding sequnce id attributre
        quality_control_sequences_add_sequence_id(qc_records, main_cols)
    if 'check_uppercase' in qc_steps:
        quality_control_sequences_check_uppercase(sequences, qc_records, main_cols)
    if 'check_sequence_length' in qc_steps:
        quality_control_sequences_check_seq_lengths(qc_records, main_cols, CONFIG)
    if 'check_nucleotid_proportions' in qc_steps:
        quality_control_sequences_check_nucleotid_props(sequences, qc_records, main_cols, CONFIG)

    # One should create a dataframe for simple QC checks
    qc_records_df = pd.DataFrame(qc_records)
    qc_records_df.columns = main_cols

    if 'check_congruence' in qc_steps:
        qc_records_df = quality_control_sequences_congruence_check(sequences, qc_records_df, CONFIG, N_fix, qc_steps,
                                                                   input_cleaned_fasta_path,
                                                                   output_folder,
                                                                   outout_congr_cleaned_fasta_path,
                                                                   output_congr_table_path)
    # Summarizing the results into a QC summary columns.
    N = len(qc_records_df.columns)
    qc_sum = qc_records_df.iloc[:, N_fix:].sum(axis=1)
    qc_records_df.loc[qc_sum == N - N_fix,'qc_summary'] = 'OK'
    qc_records_df.loc[qc_sum != N - N_fix, 'qc_summary'] = 'Failed'

    qc_records_df.fillna(0, inplace=True)
    qc_records_df.replace(0, 'Failed', inplace=True)
    qc_records_df.replace(1, 'OK', inplace=True)


    return qc_records_df, N_fix
