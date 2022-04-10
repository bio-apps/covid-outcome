import pygal
import math
import pandas as pd

class MutationDataException(Exception):
    def __init__(self, message="Mutation data is invalid"):
        self.message = message
        super().__init__(self.message)

class MLResultsDataException(Exception):
    def __init__(self, message="MLResults data is invalid"):
        self.message = message
        super().__init__(self.message)

def get_protein_bar_plot_data(mutation_data):
    '''
    Check if the mutation calculation has run successfully. Also should work valid, but empty data.
    Input is sessions[input_id]['sample_mutations']
    Output: dataframe that contains all the data necessary to make_protein_bar_plot()
    '''

    if mutation_data is None:
        raise MutationDataException("No available mutation data calculation")
    columns = mutation_data.columns
    for required_column in ["Effect", "protein_name", "mutation_name"]:
        if required_column not in columns:
            raise MutationDataException("Missing column in mutation data:", required_column)

    data = mutation_data.copy()
    data = data.loc[data["Effect"] != "synonymous_variant"]
    data["protein_name"] = data["protein_name"].fillna("UTR")

    data_grouped = data[["mutation_name", "protein_name"]].drop_duplicates().groupby(
        "protein_name").count().reset_index().sort_values(by="mutation_name", ascending=False)

    proteins = data_grouped["protein_name"].to_list()
    bar_list = []
    for protein in proteins:
        height = data_grouped.loc[data_grouped["protein_name"] == protein]["mutation_name"].values[0]
        mutation_names = data.loc[data["protein_name"] == protein]["mutation_name"].drop_duplicates()
        mutation_names = [m.split('.')[1] for m in mutation_names.to_list()]
        bar = pd.DataFrame(data={'Protein': [protein], 'Occurrences': [height], 'Mutations': [mutation_names]})
        bar_list.append(bar)
    protein_bar_plot_data = pd.concat(bar_list)

    return protein_bar_plot_data


def make_protein_bar_plot(protein_bar_plot_data, output_svg_file=None):
    bar = pygal.Bar(title="Mutations by proteins", x_title='Protein', y_title="Mutation occurrences")
    for index, row in protein_bar_plot_data.iterrows():
        bar.add(row.Protein, [{"value": row.Height, "label": row.Label}])

    if output_svg_file:
        bar.render_to_file(output_svg_file)
    return bar.render()

def get_mutation_histogram_plot_data(mutation_data):
    '''
    Check if the mutation calculation has run successfully
    Input is sessions[input_id]['sample_mutations']
    Output: dataframe that contains all the data necessary to mutation_histogram()
    '''

    if mutation_data is None:
        raise MutationDataException("No available mutation data calculation")
    columns = mutation_data.columns
    for required_column in ["Effect", "protein_name", "mutation_name", "POS", "mutation_length", "gisaid_epi_isl"]:
        if required_column not in columns:
            raise MutationDataException("Missing column in mutation data:", required_column)

    data = mutation_data.copy()
    data = data.loc[data["Effect"] != "synonymous_variant"]
    data["protein_name"] = data["protein_name"].fillna("UTR")

    protein_record_list = []

    for protein in data.sort_values(by="POS")["protein_name"].drop_duplicates().to_list():
        for mutation in data.loc[data["protein_name"] == protein]["mutation_name"]:
            mutation_records = data.loc[data["mutation_name"] == mutation]
            mutation_records_grouped = mutation_records.groupby(
                ["mutation_name", "POS", "protein_name", "mutation_length"]).count().reset_index()

            height = mutation_records_grouped["gisaid_epi_isl"][0]
            start_pos = mutation_records_grouped["POS"][0] + math.floor(
                mutation_records_grouped["mutation_length"][0] / 2)
            length = mutation_records_grouped["mutation_length"][0]
            end_pos = start_pos + length

            protein_record = pd.DataFrame(data={'Protein': [protein], 'Mutation': [mutation], 'Occurrences': [height],
                                                'Start_POS': [start_pos], 'End_POS': [end_pos]})
            protein_record_list.append(protein_record)

    mutation_histogram_plot_data = pd.concat(protein_record_list)

    return mutation_histogram_plot_data


def make_mutation_histogram(mutation_histogram_plot_data, output_svg_file):
    hist = pygal.Histogram(title="Mutations by positions", x_title='SARS-CoV-2 reference position',
                           y_title="Mutation occurrences")
    for index, row in mutation_histogram_plot_data.iterrows():
        hist.add(row.Protein, row.Records)
    if output_svg_file:
        hist.render_to_file(output_svg_file)

    return hist.render()

def get_color_for_label(label):
    if "High confidence Severe" == label:
        return "red"
    if "Low confidence Severe" == label:
        return "orange"
    if "Undefined" == label:
        return "grey"
    if "Low confidence Mild" == label:
        return "yellow"
    if "High confidence Mild" == label:
        return "green"
    raise MLResultsDataException("Color is undefined for label:", label)

def get_prediction_histogram_plot_data(ml_results):
    '''
    Check if the mutation calculation and ML has run successfully
    Input is get_ml_results(input_id)
    Output: dataframe that contains all the data necessary to mutation_histogram()
    '''

    if ml_results is None:
        raise MLResultsDataException("No available ML calculation")
    columns = ml_results.columns
    for required_column in ["prediction_confidence", "Label", "prediction_score", "user_seq_id"]:
        if required_column not in columns:
            raise MLResultsDataException("Missing column in ML result:", required_column)

    prediction_table = ml_results.copy()
    prediction_table["Display"] = prediction_table["prediction_confidence"] + " confidence " + prediction_table["Label"]
    prediction_table.loc[prediction_table["prediction_confidence"] == "Undefined", "Display"] = "Undefined"

    prediction_table_grouped = prediction_table[["prediction_score", "Display",
                                                 "user_seq_id"]].groupby(["prediction_score",
                                                                         "Display"]).count().reset_index()

    records_list = []

    for index, row in prediction_table_grouped.iterrows():
        number_of_samples = row["user_seq_id"]
        sample_ids = prediction_table.loc[prediction_table["prediction_score"] ==
                                          row['prediction_score']]["user_seq_id"].to_list()

        record = pd.DataFrame(data={"Label": [row["Display"]], "Occurrences": [number_of_samples],
                                    "Prediction_score": [row["prediction_score"]],
                                    "Color": [get_color_for_label(row["Display"])], "SampleIDs": [sample_ids]})

        records_list.append(record)

    prediction_histogram_plot_data = pd.concat(records_list)

    return prediction_histogram_plot_data


def make_prediction_histogram(prediction_histogram_plot_data, output_svg_file):

    colors = prediction_histogram_plot_data["Color"].to_list()

    hist = pygal.Histogram(title="Predictions", x_title='Predicted',
                           style=pygal.style.styles['default'](colors=colors), xrange=(0, 1))

    for index, row in prediction_histogram_plot_data.iterrows():
        hist.add(row.Label, row.Records)

    if output_svg_file:
        hist.render_to_file(output_svg_file)

    return hist.render()
