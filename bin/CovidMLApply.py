# Apply ML model for the COVID application

import random

import pandas as pd


def get_random_scores(seq_ids):

    records = []

    for seq_id in seq_ids:
        records.append([seq_id, random.random()])
    scores = pd.DataFrame(records)
    scores.columns = ['sequence_id', 'prediction_score']

    scores.set_index('sequence_id', inplace=True)

    return scores



# Ehhez kell majd egy szofisztikáltabb választó szkript
def apply_ML_model(input_id, input_fasta, CONFIG, mutation_profile, input_ages):
    # Mathcing age
    print('Applying ML model on the dataset')

    seq_id = 'sequence_id'
    seq_ids = list(mutation_profile[seq_id].unique())
    prediction_scores = get_random_scores(seq_ids)

    return prediction_scores


