#!/usr/bin/env python
# -*-coding: utf8 -*-

"""
Library for multiple alignment and mutation calling in hCov-19 analysis.
Copyright 2021 by Regina Kalcsevszki (kalcsevszkiregi@gmail.com), Balázs Ligeti (obalasz@gmail.com). All rights reserved.

"""

import os
from os.path import join
from Bio import SeqIO
import re
import numpy as np
import pandas as pd


def run_mafft(fasta_file_path, output_alignment_path, CONFIG):
    reference_path = join(CONFIG["BasePath"], CONFIG["Reference"]["ReferencePathFromBase"])
    mafft_path = CONFIG["MulipleAlignment"]["MAFFTPath"]
    mafft_threads = CONFIG["MulipleAlignment"]["CPU_cores"]
    
    mafft_cmd = '{0} --6merpair --thread {4} --addfragments {1} {2}>{3}'.format(mafft_path, fasta_file_path,
                                                                                reference_path, output_alignment_path,
                                                                                mafft_threads)
    print(mafft_cmd)
    return_value = os.system(mafft_cmd)
    if return_value != 0:
        raise(Exception("Something went wrong while running mafft:" + mafft_cmd))
    
def mutation_calling(mafft_alignment_file, raw_mutations_output_file, agg_raw_mutations_output_file, CONFIG):
    reference_file = join(CONFIG["BasePath"], CONFIG["Reference"]["ReferencePathFromBase"])
    reference_genome_id = CONFIG["Reference"]["GenomeID"]
    subject_genome_id_gisaid = CONFIG["MutationCalling"]["add_gisaid_id"]

    original_ref_genome = SeqIO.read(reference_file, 'fasta')
    all_records = [record for record in SeqIO.parse(mafft_alignment_file, "fasta")]

    genome_ids =[record.id for record in all_records]
    id2genome = {record.id : record for record in all_records}
    #charset = get_char_stat_from_alignment(all_records)

    ref_genome = id2genome[reference_genome_id]
    ref_genome_raw = re.sub(r"-+", "", str(ref_genome.seq)).upper()


    # Running the mutation finding process

    header_aggragated_mutatations = ['Reference ID', 'Target ID', 'Mut_tpye', 'Mutation', 'Variation_type', 'REF_START', 'REF_END', 'TAR_START', 'TAR_END',
                                    'REF_STR', 'TAR_STR','Mutation_ids', 'gisaid_epi_isl', 'REF_STR_SNPEFF', 'REF_STR_END_SNPEFF']
    header_mutations = ['Reference ID', 'Target ID', 'Mutation', 'Mut_tpye', 'MA_pos', 'POS_REF', 'POS_TAR', 'REF_STRING', 'TAR_STRING', 'gisaid_epi_isl']

    with open(raw_mutations_output_file,'w') as fout,  open(agg_raw_mutations_output_file, 'w') as fout_agg:
        fout.write('\t'.join(header_mutations) + '\n')
        fout_agg.write('\t'.join(header_aggragated_mutatations) + '\n')
        for i, act_genome in enumerate(all_records):        
            act_results = mutation_calling_compare_to_ref(act_genome, ref_genome)
            if len(act_results) > 0:
                # Some basic filtering? 
                act_results = filtering_internal_results(act_results)
                #print(act_results)
                # Adding target id
                if subject_genome_id_gisaid:
                    adding_gisaid_id(act_results)
                else:
                    for record in act_results:
                        record.append(record[1])
                aggr_mutation_list = aggregate_mutation_profile(act_results)
                agg_mutation_records = [create_mutation_record(r, ref_genome_raw) for r in aggr_mutation_list]

                for record in act_results:
                    a = fout.write('\t'.join([str(act_rec) for act_rec in record]) + '\n')
                for record in agg_mutation_records:
                    a = fout_agg.write('\t'.join([str(act_rec) for act_rec in record]) + '\n')
            if i % 100 == 0:
                print(i/len(all_records)*100)
    fout.close()


# Mutation calling
# meg kell találni a referencia genome első nem gap karakterét
def get_first_nucleotid_in_reference_genome(ref_genome):
    '''
    Get the first non-gap nucleotide sequence
    '''
    for i, act_char in enumerate(ref_genome.seq):
        if act_char != '-':
            return i
    return len(ref_genome.seq)


def mutation_calling_compare_to_ref(act_genome, ref_genome, gap='-'):
    ref_genome_seq = str(ref_genome.seq).upper()
    Nref = len(ref_genome_seq)
    not_gap_pointer_ref = get_first_nucleotid_in_reference_genome(ref_genome)  # első nem gap pozició koordinátája
    ref_id = ref_genome.description
    tar_id = act_genome.description

    # analysis
    act_genome_seq = str(act_genome.seq).upper()
    Nq = len(act_genome_seq)

    # ref_genome_seq
    # act_genome_seq
    if Nref != Nq:
        print(
            'Warning the sequence lenghts are not match! The length of the query %d, but the reference is %d na long.',
            (Nq, Nref))

    results = []

    ref_counter = 0
    ref_gap_counter = 0
    ref_internal_gap_counter = 1

    query_counter = 0

    for i in range(Nref):
        ref_char = ref_genome_seq[i]
        tar_char = act_genome_seq[i]

        # setting the internal gap counter
        if ref_genome_seq[i - 1] != gap:
            ref_internal_gap_counter = 1
        elif ref_genome_seq[i - 1] == gap:
            ref_internal_gap_counter += 1

        # increasing query counter
        if tar_char != gap:
            query_counter += 1

        if ref_char != gap:
            ref_counter += 1

        if ref_char == gap and tar_char == gap:
            ref_gap_counter += 1

        # print('Status: ref_counter, query_counter: ', ref_counter, query_counter)

        if ref_char != tar_char:
            # print('Mutation: %s vs %s' % (ref_char, tar_char))
            if ref_char == gap and tar_char != gap:
                # print('  GAP IN REF (insertion) ANALYSIS')
                # print('  i', i, ref_counter, ref_gap_counter, (i - ref_counter + ref_gap_counter))
                mut_type = 'I'
                if i < not_gap_pointer_ref:
                    mut_ref_string = 'G' + str(-1 * not_gap_pointer_ref + i) + 'REF1' + tar_char
                else:
                    # print('  ref_gap_counter, ', ref_gap_counter, ref_internal_gap_counter)
                    mut_ref_string = 'G' + str(
                        i - (ref_counter + ref_gap_counter) + ref_internal_gap_counter) + 'REF' + str(
                        ref_counter) + tar_char
                    # print('   prev: ', ref_genome_seq[i-1])
                    if ref_genome_seq[i - 1] == gap:
                        # print('INCRESING THE INTERNAL GAP COUNTER')
                        ref_internal_gap_counter += 1
                ref_gap_counter += 1
                # print(' ', mut_ref_string)
                record = [ref_id, tar_id, mut_ref_string, 'insertion', i, ref_counter, query_counter, ref_char,
                          tar_char]

            else:
                # print('  Subsitition analysis')
                if tar_char == gap:
                    mut_string = 'deletion'
                else:
                    mut_string = 'substitution'

                mut_type = 'S'
                mut_ref_string = 'S' + ref_char + str(ref_counter) + tar_char
                # print('  ', mut_ref_string)
                record = [ref_id, tar_id, mut_ref_string, mut_string, i, ref_counter, query_counter, ref_char, tar_char]

                # ref_counter +=1
            # print(results)
            results.append(record)

        # else:
        #    ref_counter +=1

        # print('Match: %s vs %s', (ref_char, tar_char))
    return results


def filtering_internal_results(mutation_records):
    valid_chars = {c: 0 for c in ['A', 'T', 'C', 'G', 'N', '-']}
    new_records = []
    for record in mutation_records:
        tar_char = record[8]
        ref_char = record[7]
        if tar_char == 'N' or ref_char == 'N':
            continue
        if tar_char not in valid_chars:
            continue

        new_records.append(record)
    return new_records


def adding_gisaid_id(mutation_records):
    for record in mutation_records:
        record.append(record[1].split('|')[1])


def set_prev_values(record):
    ref_id, target_id, mutation, mut_type, pos_ma, pos_ref, pos_tar, ref_char, tar_char, gisaid_id = record

    act_mut_type = mut_type
    act_target_id = target_id
    act_pos_ref = pos_ref
    act_pos_target = pos_tar

    return act_mut_type, act_target_id, act_pos_ref, act_pos_target


def create_new_mutation_record(record):
    ref_id, target_id, mutation, mut_type, pos_ma, pos_ref, pos_tar, ref_char, tar_char, gisaid_id = record
    act_mutation = {}
    act_mutation['pos_tar_end'] = pos_tar
    act_mutation['pos_tar_start'] = pos_tar
    act_mutation['pos_ref_start'] = pos_ref
    act_mutation['pos_ref_end'] = pos_ref
    act_mutation['tar_chars'] = [tar_char]
    act_mutation['ref_chars'] = [ref_char]  # Ez itt nem számít, mivel a referencia elvileg gap

    act_mutation['mut_ids'] = [mutation]

    act_mutation['mut_type'] = mut_type
    act_mutation['ref_id'] = ref_id
    act_mutation['target_id'] = target_id
    act_mutation['gisaid_id'] = gisaid_id
    return act_mutation


def create_mutation_record(mutation_data, reference_genome=None, target_genome=None, nucl_context=False):
    # Print creating a new record from the data object
    # Record structure:
    # ['Reference ID', 'Target ID', 'Mutation', 'Mut_tpye', 'Variation_type', 'POS_REF_START', 'POS_REF_END', 'POS_TAR_START', 'POS_TAR_END',
    # 'STRING_REF', 'STRING_TAR', 'covered_mutatations', 'gisaid_epi_isl']

    if len(mutation_data['ref_chars']) > 1 or len(mutation_data['tar_chars']) > 1:
        Variation_type = 'non_single'
    else:
        Variation_type = 'single'

    Mutation = 'temp'

    ref_string = ''.join(mutation_data['ref_chars'])
    tar_string = ''.join(mutation_data['tar_chars'])

    # Creating the new ID string
    if mutation_data['mut_type'] == 'substitution':
        Mutation = 'S' + ref_string + str(mutation_data['pos_ref_start']) + tar_string

    elif mutation_data['mut_type'] == 'insertion':
        # G1REF11A
        # Itt előfordulhat, hogy duplikálódik az ID, pl. ha egymást követő inzerciók vannak, ezért inklúdolni kellene
        # a target koordinátákat is

        Mutation = 'G' + str(len(tar_string)) + 'REF' + str(mutation_data['pos_ref_start']) + tar_string


    elif mutation_data['mut_type'] == 'deletion':
        Mutation = 'S' + ref_string + str(mutation_data['pos_ref_start']) + tar_string
    else:
        Mutation = ';'.join(mutation_data['mut_ids'])
    # Adding subsequence char in the reference if possible
    prev_ref_char = ''
    next_ref_char = ''
    if reference_genome is not None:
        if mutation_data['pos_ref_start'] == 0:
            prev_ref_char = reference_genome[0]
        elif mutation_data['mut_type'] == 'deletion' or mutation_data['mut_type'] == 'substitution':
            prev_ref_char = reference_genome[mutation_data['pos_ref_start'] - 2]
        elif mutation_data['mut_type'] == 'insertion':
            if mutation_data['pos_ref_start'] > 1:
                prev_ref_char = reference_genome[mutation_data['pos_ref_start'] - 1]
            else:
                prev_ref_char = ''
    else:
        prev_ref_char = ''
    if reference_genome is not None:
        if mutation_data['pos_ref_end'] < len(reference_genome):
            next_ref_char = reference_genome[mutation_data['pos_ref_end']]
        else:
            next_ref_char = reference_genome[-1]

    new_record = [mutation_data['ref_id'],
                  mutation_data['target_id'],
                  mutation_data['mut_type'],
                  Mutation,
                  Variation_type,
                  mutation_data['pos_ref_start'],
                  mutation_data['pos_ref_end'],
                  mutation_data['pos_tar_start'],
                  mutation_data['pos_tar_end'],

                  ref_string,
                  tar_string,
                  ','.join(mutation_data['mut_ids']),
                  mutation_data['gisaid_id'],
                  prev_ref_char,
                  next_ref_char
                  ]

    return new_record


def aggregate_mutation_profile(mutation_records):
    act_mut_type = None
    act_target_id = None
    act_pos_ref = None
    act_pos_target = None

    act_mutation = {}
    act_mutations = []

    for record in mutation_records:
        ref_id, target_id, mutation, mut_type, pos_ma, pos_ref, pos_tar, ref_char, tar_char, gisaid_id = record
        # Ha a mostani mutáció az előző folytatás, akkor következőket kellene csinálni
        # Ugyanaz target genom?
        if act_target_id == target_id:
            if mut_type == act_mut_type:
                if mut_type == 'insertion':
                    # print('Ez egy inzerció')
                    # Referenciában ugyanaz a mutató van és poziciók egymást követik, akkor összelehet vonni az előzővel
                    if pos_ref == act_pos_ref and pos_tar - 1 == act_pos_target:
                        # print('Extending previous insertion')
                        act_mutation['pos_tar_end'] = pos_tar
                        act_mutation['tar_chars'].append(tar_char)
                        act_mutation['mut_ids'].append(mutation)

                        # Be kell állítani az act-okat
                        act_mut_type, act_target_id, act_pos_ref, act_pos_target = set_prev_values(record)

                    else:
                        # print('New insertion have been detected!')
                        # print('Dupming the previous one')
                        act_mutations.append(act_mutation)
                        # print('Create a new record! ')
                        act_mutation = create_new_mutation_record(record)
                        # print('Default értékek beállítása')
                        act_mut_type, act_target_id, act_pos_ref, act_pos_target = set_prev_values(record)


                elif mut_type == 'deletion':
                    # print('Deléció')
                    if pos_ref - 1 == act_pos_ref and pos_tar == act_pos_target:
                        # print('cont Deletion detected')
                        act_mutation['pos_ref_end'] = pos_ref
                        act_mutation['ref_chars'].append(ref_char)
                        act_mutation['mut_ids'].append(mutation)
                        act_mut_type, act_target_id, act_pos_ref, act_pos_target = set_prev_values(record)

                    else:
                        act_mutations.append(act_mutation)
                        act_mutation = create_new_mutation_record(record)
                        act_mut_type, act_target_id, act_pos_ref, act_pos_target = set_prev_values(record)



                elif mut_type == 'substitution':
                    # print('subsztitució')
                    if pos_ref - 1 == act_pos_ref and pos_tar - 1 == act_pos_target:
                        # print('Egymás követő szubsztuciók')
                        act_mutation['pos_tar_end'] = pos_tar
                        act_mutation['tar_chars'].append(tar_char)
                        act_mutation['pos_ref_end'] = pos_ref
                        act_mutation['ref_chars'].append(ref_char)
                        act_mutation['mut_ids'].append(mutation)
                        act_mut_type, act_target_id, act_pos_ref, act_pos_target = set_prev_values(record)

                    else:
                        # print('Nem egymást követő szubsztituciók')
                        act_mutations.append(act_mutation)
                        act_mutation = create_new_mutation_record(record)
                        act_mut_type, act_target_id, act_pos_ref, act_pos_target = set_prev_values(record)

            else:
                # print('mutáció kiírása')
                act_mutations.append(act_mutation)
                act_mutation = create_new_mutation_record(record)
                act_mut_type, act_target_id, act_pos_ref, act_pos_target = set_prev_values(record)

        else:
            # print('Uj genom, az előző mutáció kiírása')
            if act_target_id:
                act_mutations.append(act_mutation)
            # print('Default értékek beállítása, új mutációs objektum létrehozása')
            act_mut_type, act_target_id, act_pos_ref, act_pos_target = set_prev_values(record)
            act_mutation = create_new_mutation_record(record)
            # előző mutáció kiírása
    act_mutations.append(act_mutation)
    return act_mutations


######### BUILDING VCF FILES #############

def setMutationCode(row):
    mut_code = None
    if row["Mut_tpye"] == "substitution":
        mut_code = "c." + str(row["POS"]) + row["REF"] + ">" + row["ALT"]
    if row["Mut_tpye"] == "deletition":
        mut_code = "c." + str(row["POS"]) + "del" + row["REF"]
    if row["Mut_tpye"] == "insertion":
        mut_code = "c." + str(row["POS"]) + "ins" + row["ALT"]

    return mut_code


def get_annotation_df_mutation_id_str(mutation_id_string):
    mut_string_id_l = mutation_id_string.split('=')
    if len(mut_string_id_l) > 1:
        return mut_string_id_l[-1]
    else:
        return mutation_id_string


def get_annotation_df(df):
    annotations_list = []

    for row in df.iterrows():
        CHROM = row[1][0]
        POS = row[1][1]
        ID = row[1][2]
        REF = row[1][3]
        ALT = row[1][4]
        QUAL = row[1][5]
        FILTER = row[1][6]
        INFO_splitted = row[1][7].split(';')
        INFO = row[1][7].split(';')[3].split('=')[1]
        MUTATATION_ID = get_annotation_df_mutation_id_str(row[1][7].split(';')[2])

        # Parse mutation id if possible

        annotations = INFO.split(',')
        j = 1
        # print(row[1][7])
        # print(INFO_splitted)
        # print(annotations)

        for a in annotations:
            attributes = a.split('|')
            names = ["Allele", "Effect", "Putative impact", "Gene Name", 'Gene ID',
                     'Feature type', 'Feature ID', 'Transcript biotype', 'Rank / total',
                     'HGVS.c', 'HGVS.p', 'cDNA_position / cDNA_len', 'CDS_position / CDS_len',
                     'Protein_position / Protein_len', 'Distance to feature', 'Messages']

            ann_dict = dict(zip(names, attributes))
            ann_dict['#CHROM'] = CHROM
            ann_dict['POS'] = POS
            ann_dict['ID'] = ID
            ann_dict['REF'] = REF
            ann_dict['ALT'] = ALT
            ann_dict['QUAL'] = QUAL
            ann_dict['FILTER'] = FILTER
            ann_dict['VCF_MUTATION_ID'] = MUTATATION_ID
            annotations_list.append(ann_dict)

    annotations_df = pd.DataFrame(annotations_list)
    annotations_df = annotations_df.reset_index()
    annotations_df = annotations_df[annotations_df["ALT"].str.contains("A|C|G|T")]

    return annotations_df


def get_protein_annotation_ncbi(annotation_file):
    # annotation_file = join(base_path, 'data/NCBI/NC_045512.2_annot.xlsx')
    genome_annot = pd.read_excel(annotation_file, engine='openpyxl')

    genome_annot_proteins_proteins_with_id = genome_annot[genome_annot.protein_id.notnull()][
        ['protein_id', 'product', 'display_name', 'note']]
    genome_annot_proteins_locus = genome_annot[(genome_annot.start_pos > 21552) &
                                               (genome_annot.protein_id.notnull())][
        ['locus_tag', 'product', 'display_name', 'note']]
    genome_annot_proteins_locus.columns = ['protein_id', 'product', 'display_name', 'note']

    genome_annot_proteins = pd.concat([genome_annot_proteins_locus, genome_annot_proteins_proteins_with_id])

    # genome_annot_proteins['protein_name'] = genome_annot_proteins.apply(lambda x: x['product'].replace(' ', '_'), axis=1)
    genome_annot_proteins['protein_name'] = genome_annot_proteins['display_name']

    genome_annot_proteins['Feature ID'] = genome_annot_proteins['protein_id']

    return genome_annot, genome_annot_proteins


def mutation_annotation_preprocessing_get_protein_annotations(protein_annotations, genome_annot_proteins):
    protein_annotations_non_duplicated = protein_annotations[~protein_annotations.duplicated(subset=['ID'], keep=False)]
    protein_annotations_duplicated = protein_annotations[protein_annotations.duplicated(subset=['ID'], keep=False)]

    protein_annotations_duplicated_preffered_protein_id = protein_annotations_duplicated.merge(
        genome_annot_proteins[['Feature ID']], how='inner', left_on='Feature ID', right_on='Feature ID')
    protein_annotations_duplicated_wiout_protein_id = protein_annotations_duplicated[
        ~protein_annotations_duplicated['ID'].isin(protein_annotations_duplicated_preffered_protein_id['ID'])]

    return protein_annotations_non_duplicated, protein_annotations_duplicated_preffered_protein_id, protein_annotations_duplicated_wiout_protein_id


def mutation_annotation_preprocessing_get_mature_and_non_integenic_annotations(annotations, genome_annot_proteins):
    effect_proteins = ['conservative_inframe_insertion&splice_region_variant',
                       'start_lost&disruptive_inframe_insertion',
                       'start_lost&disruptive_inframe_deletion',
                       'disruptive_inframe_deletion&splice_region_variant',
                       'stop_lost&disruptive_inframe_insertion&splice_region_variant',
                       'stop_lost&disruptive_inframe_deletion&splice_region_variant',
                       'stop_lost&conservative_inframe_deletion&splice_region_variant',
                       'conservative_inframe_deletion&splice_region_variant',
                       'start_lost&conservative_inframe_deletion',
                       'stop_gained&conservative_inframe_insertion',
                       'initiator_codon_variant',
                       'stop_gained&splice_region_variant',
                       'stop_gained&disruptive_inframe_insertion',
                       'frameshift_variant&splice_region_variant',
                       'stop_gained&disruptive_inframe_deletion',
                       'splice_region_variant&stop_retained_variant',
                       'frameshift_variant&start_lost',
                       'frameshift_variant&stop_lost&splice_region_variant',
                       'stop_lost&splice_region_variant',
                       'splice_region_variant&synonymous_variant',
                       'start_lost',
                       'frameshift_variant&stop_gained',
                       'missense_variant&splice_region_variant',
                       'conservative_inframe_insertion',
                       'disruptive_inframe_insertion',
                       'conservative_inframe_deletion',
                       'disruptive_inframe_deletion',
                       'stop_gained',
                       'frameshift_variant',
                       'synonymous_variant',
                       'missense_variant']

    effect_nonprotein = ['intergenic_region']

    non_protein_annotations = annotations[annotations['Effect'].isin(effect_nonprotein)]
    protein_annotations = annotations[annotations['Effect'].isin(effect_proteins)]

    prot_uq, prot_mature, prot_nonuq_non_mature = mutation_annotation_preprocessing_get_protein_annotations(
        protein_annotations, genome_annot_proteins)

    mutation_annotations = pd.concat([non_protein_annotations,
                                      prot_uq,
                                      prot_mature,
                                      prot_nonuq_non_mature])

    return mutation_annotations


def map_to_gisaid_AA_notation(code):
    if code != None:
        code = str(code)
        AAdict = {"Ala": "A", "Arg": "R", "Asn": 'N', "Asp": "D", "Cys": 'C', 'Gln': "Q", "Glu": 'E',
                  'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F',
                  'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'}
        for key, value in AAdict.items():
            code = code.replace(key, value)
    return code


def mutation_annotation_add_protein_data_mut_name(mut_record):
    new_name = 'temp'

    if isinstance(mut_record['protein_id'], str):
        if len(mut_record['protein_id']) > 0:
            new_name = mut_record['protein_name'] + '_' + map_to_gisaid_AA_notation(mut_record['HGVS.p'])
        else:
            new_name = mut_record['ID'] + '_INVALID'
    elif np.isnan(mut_record['protein_id']):
        new_name = mut_record['Gene Name'] + '_' + mut_record['HGVS.c']

    return new_name


def mutation_annotation_add_protein_data(mutation_annotations, genome_annot_proteins):
    mutation_annotations_name = mutation_annotations.merge(genome_annot_proteins, how='left', left_on='Feature ID',
                                                           right_on='Feature ID')

    mutation_annotations_name['mutation_name'] = mutation_annotations_name.apply(
        lambda x: mutation_annotation_add_protein_data_mut_name(x), axis=1)

    return mutation_annotations_name


def mutation_annotation_preprocessing(annotations, genome_annot_proteins):
    ''' Preprocessing the raw output coming from SNPEFF'''

    # Filtering the annotations!
    mutation_annotations = mutation_annotation_preprocessing_get_mature_and_non_integenic_annotations(annotations,
                                                                                                      genome_annot_proteins)

    # Adding amibous flag
    mutation_annotations.loc[
        mutation_annotations.duplicated(subset=['ID'], keep=False), 'Mutation_annotation_ambiguity'] = 'ambiguous'
    mutation_annotations.loc[
        ~mutation_annotations.duplicated(subset=['ID'], keep=False), 'Mutation_annotation_ambiguity'] = 'non_ambiguous'

    # Adding mutation name?
    mutation_annotations = mutation_annotation_add_protein_data(mutation_annotations, genome_annot_proteins)

    return mutation_annotations


def filter_and_preprocess_aggr_mutation_table(raw_agg_mutations):
    columns_to_drop = ['Reference ID', 'Target ID', 'Mutation_ids']
    vcf_cols_to_keep = ['Mutation', 'Mut_tpye', 'REF_START', 'REF_END', 'REF_STR', 'TAR_STR', 'REF_STR_SNPEFF',
                        'REF_STR_END_SNPEFF']

    print('Filtering the large insertions and deletions at beginning and at the and of the genome!')

    raw_agg_mutations_cleaned = raw_agg_mutations[~(
            (raw_agg_mutations['Mut_tpye'] == 'deletion') & (raw_agg_mutations['REF_START'] == 1) & (
            np.abs(raw_agg_mutations['REF_START'] - raw_agg_mutations['REF_END']) > 1)) &
                                                  ~((raw_agg_mutations['Mut_tpye'] == 'deletion') & (
                                                          raw_agg_mutations['REF_END'] == 29891) & (np.abs(
                                                      raw_agg_mutations['REF_START'] - raw_agg_mutations[
                                                          'REF_END']) > 1)) &
                                                  ~((raw_agg_mutations['Mut_tpye'] == 'insertion') & (
                                                          raw_agg_mutations['REF_END'] == 29891) & (np.abs(
                                                      raw_agg_mutations['TAR_START'] - raw_agg_mutations[
                                                          'TAR_END']) > 1)) &
                                                  ~((raw_agg_mutations['Mut_tpye'] == 'insertion') & (
                                                          raw_agg_mutations['REF_START'] == 1) & (np.abs(
                                                      raw_agg_mutations['TAR_START'] - raw_agg_mutations[
                                                          'TAR_END']) > 1))
                                                  ]
    print('Finished')
    len(raw_agg_mutations_cleaned)
    # adding mutation_lenght
    print('Calculating the mutation lengths')
    raw_agg_mutations_cleaned['mutation_length'] = raw_agg_mutations_cleaned.apply(
        lambda x: max(np.abs(x['TAR_END'] - x['TAR_START']),
                      np.abs(x['REF_END'] - x['REF_START'])) + 1, axis=1)
    print('Finished!')
    # Túl hosszú mutációk kihagyása
    print('Removing the large mutations!')
    raw_agg_mutations_cleaned = raw_agg_mutations_cleaned[raw_agg_mutations_cleaned.mutation_length < 50]
    print('Finished')
    len(raw_agg_mutations_cleaned)

    # Duplikálódott ID cseréje unique ID-ra
    raw_agg_mutations_cleaned.loc[
        raw_agg_mutations_cleaned.duplicated(subset=['gisaid_epi_isl', 'Mutation'], keep=False), 'Mutation'] = \
        raw_agg_mutations_cleaned.loc[
            raw_agg_mutations_cleaned.duplicated(subset=['gisaid_epi_isl', 'Mutation'], keep=False), 'Mutation_ids']

    print('Dropping columns: ' + str(columns_to_drop))
    raw_agg_mutations_cleaned.drop(columns_to_drop, axis=1, inplace=True)
    print('Finished!')
    mutations = raw_agg_mutations_cleaned[vcf_cols_to_keep].drop_duplicates()

    return raw_agg_mutations_cleaned, mutations


def annotate_mutations_with_vcf(snp_eff_path, snp_eff_reference_genome_id, mutation_vcf_file,
                                mutation_annotated_vcf_file):
    annotation_cmd = '{0} {1} {2} > {3}'.format(snp_eff_path, snp_eff_reference_genome_id, mutation_vcf_file,
                                                mutation_annotated_vcf_file)
    print('Running the annotation!')
    print(annotation_cmd)
    os.system(annotation_cmd)
    print('Finished!')



"""
Converts aggragated mutation profiles to a VCF file, that is approproate for SNPEff input

Input: df: path to aggregated mutations dataframe

Output: vcf dataframe
"""


def toVCF(df):

    vcf_dfs = []
    try:
        vcf_subs = get_vcf_for_substitutions(df)
        vcf_dfs.append(vcf_subs)
    except:
        print('Error in substituion calculation!')
    try:
        vcf_ins = get_vcf_for_insertions(df)
        vcf_dfs.append(vcf_ins)
    except:
        print('Error in insertion calculation!')
    try:
        vcf_del = get_vcf_for_deletions(df)
        vcf_dfs.append(vcf_del)
    except:
        print('Error in deletion. VCF')
    extended_vcf_df = pd.concat(vcf_dfs)
    vcf_df = add_general_fields(extended_vcf_df)

    return select_vcf_columns(vcf_df)


def get_vcf_for_substitutions(df):
    out_df_subs = df.loc[df["Mut_tpye"] == "substitution"].copy()

    out_df_subs["POS"] = out_df_subs["REF_START"]
    out_df_subs["REF"] = out_df_subs["REF_STR"]
    out_df_subs["ALT"] = out_df_subs["TAR_STR"]

    out_df_subs["MUT_CODE"] = "c." + out_df_subs["POS"].apply(str) + out_df_subs["REF"] + ">" + out_df_subs["ALT"]

    return out_df_subs


def get_vcf_for_insertions(df):
    out_df_ins = df.loc[df["Mut_tpye"] == "insertion"].copy()

    out_df_ins["POS"] = out_df_ins["REF_START"]
    out_df_ins.loc[out_df_ins["POS"] == 0, "POS"] = 1

    out_df_ins["REF"] = out_df_ins["REF_STR_SNPEFF"]
    out_df_ins["ALT"] = out_df_ins["TAR_STR"]

    # (https://samtools.github.io/hts-specs/VCFv4.2.pdf 1.4.1.4) the REF and ALT Strings must include the base before the event (which must be reflected in the POS field),
    out_df_ins.loc[out_df_ins["REF_START"] != 0, "ALT"] = out_df_ins.loc[out_df_ins["REF_START"] != 0, "REF"] + \
                                                          out_df_ins.loc[out_df_ins["REF_START"] != 0, "ALT"]
    # unless the event occurs at position 1 on the contig in which case it must include the base after the event
    out_df_ins.loc[out_df_ins["REF_START"] == 0, "ALT"] = out_df_ins.loc[out_df_ins["REF_START"] == 0, "ALT"] + \
                                                          out_df_ins.loc[out_df_ins["REF_START"] == 0, "REF"]

    out_df_ins["MUT_CODE"] = "c." + out_df_ins["REF_START"].apply(str) + "_" + (out_df_ins["REF_START"] + 1).apply(
        str) + "ins" + out_df_ins["TAR_STR"]

    return out_df_ins


def get_vcf_for_deletions(df):
    out_df_del = df.loc[df["Mut_tpye"] == "deletion"].copy()

    out_df_del["POS"] = out_df_del["REF_START"]
    out_df_del.loc[out_df_del["POS"] != 1, "POS"] = out_df_del.loc[out_df_del["POS"] != 1, "POS"] - 1

    out_df_del["REF"] = out_df_del["REF_STR"]

    # (https://samtools.github.io/hts-specs/VCFv4.2.pdf 1.4.1.4) the REF and ALT Strings must include the base before the event (which must be reflected in the POS field),
    out_df_del.loc[out_df_del["REF_START"] != 1, "ALT"] = out_df_del.loc[out_df_del["REF_START"] != 1, "REF_STR_SNPEFF"]
    out_df_del.loc[out_df_del["REF_START"] != 1, "REF"] = out_df_del.loc[out_df_del["REF_START"] != 1, "ALT"] + \
                                                          out_df_del.loc[out_df_del["REF_START"] != 1, "REF"]
    # unless the event occurs at position 1 on the contig in which case it must include the base after the event
    out_df_del.loc[out_df_del["REF_START"] == 1, "ALT"] = out_df_del.loc[
        out_df_del["REF_START"] == 1, "REF_STR_END_SNPEFF"]
    out_df_del.loc[out_df_del["REF_START"] == 1, "REF"] = out_df_del.loc[out_df_del["REF_START"] == 1, "REF"] + \
                                                          out_df_del.loc[out_df_del["REF_START"] == 1, "ALT"]

    out_df_del["MUT_CODE"] = "c." + out_df_del["REF_START"].apply(str) + "_" + out_df_del["REF_END"].apply(str) + "del"
    out_df_del.loc[out_df_del["REF_START"] == out_df_del["REF_END"], "MUT_CODE"] = out_df_del.loc[
        out_df_del["REF_START"] == out_df_del["REF_END"], "REF_START"].apply(lambda x: "c." + str(x) + "del")

    return out_df_del


def add_general_fields(df):
    # REFERENCE_GENOME = "hCoV-19/Wuhan/WIV04/2019" # MN996528.1
    SNPEFF_DB_ID = "NC_045512.2"

    df["#CHROM"] = SNPEFF_DB_ID
    df["QUAL"] = 100
    df["FILTER"] = "PASS"
    df["ID"] = df["MUT_CODE"]

    return df


def select_vcf_columns(df):
    df = df[["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "MUT_CODE", "Mutation"]]
    df["CNT"] = 1
    df = df[["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "MUT_CODE", "CNT", "Mutation"]].drop_duplicates()
    df = df.reset_index()
    df["INFO"] = "NS=" + df["CNT"].apply(str) + ';' "MC=" + df["MUT_CODE"] + ';' + 'Mutation=' + df["Mutation"]

    return df[["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]]
