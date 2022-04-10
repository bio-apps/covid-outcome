# Ebben lesz az API
import os.path

from fastapi import FastAPI, HTTPException, Response, File, UploadFile, BackgroundTasks
from pydantic import BaseModel
from typing import Optional, List
from fastapi.encoders import jsonable_encoder
# from fastapi.middleware.cors import CORSMiddleware
from starlette.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi_sessions.frontends.implementations import SessionCookie, CookieParameters
from fastapi_sessions.session_verifier import SessionVerifier
from uuid import UUID
from uuid import uuid4
from fastapi import Depends
from fastapi_sessions.backends.implementations import InMemoryBackend
import aiofiles
import io
import chardet
import pandas as pd
import time
import CovidOutcome

app = FastAPI()

class FormData(BaseModel):
    chars: str

class CovData(BaseModel):
    sequences: str
    age: float


class CovidSequenceData(BaseModel):
    sequences: str

class AgeData(BaseModel):
    age_list: str


class SessionData(BaseModel):
    session_id: str

class ProteinBarPlotData(BaseModel):
    proteinName: str
    mutations : List[str]
    occurences: int

class ProteinBarPlotDataRecord(BaseModel):
    proteinName: str
    mutations : List[str]
    occurrences: int

class ProteinBarPlotData(BaseModel):
    records: List[ProteinBarPlotDataRecord]

class MutationHistogramPlotDataRecord(BaseModel):
    proteinName: str
    mutation : str
    occurrences: int
    startPos: int
    endPos: int

class MutationHistogramPlotData(BaseModel):
    records: List[MutationHistogramPlotDataRecord]

class PredictionHistogramPlotDataRecord(BaseModel):
    label: str
    occurrences: int
    predictionScore: float
    color: str
    sampleIds: List[str]

class PredictionHistogramPlotData(BaseModel):
    records: List[PredictionHistogramPlotDataRecord]

class ArtificialCovidGenomes(BaseModel):
    '''
    Storing the mutation data

    '''
    genomes : List[List[str]]


class BasicVerifier(SessionVerifier[UUID, SessionData]):
    def __init__(
            self,
            *,
            identifier: str,
            auto_error: bool,
            backend: InMemoryBackend[UUID, SessionData],
            auth_http_exception: HTTPException,
    ):
        self._identifier = identifier
        self._auto_error = auto_error
        self._backend = backend
        self._auth_http_exception = auth_http_exception

    @property
    def identifier(self):
        return self._identifier

    @property
    def backend(self):
        return self._backend

    @property
    def auto_error(self):
        return self._auto_error

    @property
    def auth_http_exception(self):
        return self._auth_http_exception

    def verify_session(self, model: SessionData) -> bool:
        """If the session exists, it is valid"""
        return True


backend = InMemoryBackend[UUID, SessionData]()
verifier = BasicVerifier(
    identifier="general_verifier",
    auto_error=True,
    backend=backend,
    auth_http_exception=HTTPException(status_code=403, detail="invalid session"),
)

# Session handling
cookie_params = CookieParameters()
cookie_params.samesite = "none"
#cookie_params.samesite = "Strict"
cookie_params.secure="True"

print(cookie_params)

# Uses UUID
cookie = SessionCookie(
    cookie_name="covid-cookie",
    identifier="general_verifier",
    auto_error=True,
    secret_key="DONOTUSE",
    cookie_params=cookie_params,
)



origins = [
    "covidoutcome.bio-ml.com",
    "wwww.covidoutcome.bio-ml.com",
    "https://covidoutcome.bio-ml.com",
    "https://www.covidoutcome.bio-ml.com",
    "https://covidoutcome.bio-apps.itk.ppke.hu",
    "https://covidoutcome.com",
    "https://www.covidoutcome.com",
    "http://localhost",
    "http://localhost:8080",
    "https://localhost",
    "https://localhost:8080",
    "https://lively-pebble-092280103.azurestaticapps.net",
    "https://lively-pebble-092280103.azurestaticapps.net",
    "http://localhost:3000",
    "https://localhost:3000",
    "https://bio-apps.itk.ppke.hu",
    "bio-apps.itk.ppke.hu",
    "https://covid-outcome-api-2m7pws7dyq-ew.a.run.app",
    "https://covid-outcome-front-end-2m7pws7dyq-ew.a.run.app",
    "covid-outcome-front-end-2m7pws7dyq-ew.a.run.app"
]
origins = list(set(origins))

def background_calculation_of_mutations(session_id):
    """
    Calculating the mutation profile in the background
    """
    print('Identify mutation in the sequence!')
    # session_id = session_data.session_id
    print(session_id)
    try:
        CovidOutcome.mutation_calling_procedure(session_id)
    except ValueError as VE:
        print('Value error during the mutation call process: {0}'.format(str(VE)))
        print('Set the mutation calling status flag to ERROR')
        CovidOutcome.sessions[session_id]['status']['mutation_calling']['status'] = 'ERROR'
        raise HTTPException(
            status_code=500,
            detail=str(VE),
            headers={"Mutation calling error!": "Check the logs!"},
        )
    # If status flag is not okay, then set to ERROR!
    if CovidOutcome.sessions[session_id]['status']['mutation_calling']['status'] != 'OK':
        print(CovidOutcome.sessions[session_id]['status']['mutation_calling']['status'])
        print('Set the mutation calling status flag to ERROR')
        CovidOutcome.sessions[session_id]['status']['mutation_calling']['status'] = 'ERROR'
    print('done')

@app.get("/")
async def root():
    '''
    Health check ...
    '''
    print(app.routes)
    return {"message": "Köszi!"}


@app.get("/health_check")
async def health_check():
    ''' Just for checking if the API is available '''
    return {'Health': 'OK'}


@app.get("/maxsequencesize")
async def get_max_file_size():
    return {'MaxSequenceSequenceChars': CovidOutcome.CONFIG['MaxSequenceSequenceChars']}


@app.get("/maxfilesize")
async def get_max_file_size():
    return {'MaxFileSize': CovidOutcome.CONFIG['MaxFileSize']}

@app.get("/get_mutation_name_suggestions")
async def get_mutation_name_suggestion(q: Optional[str] = None):
    import json
    print('Query the query:', q)
    suggested_muts = CovidOutcome.get_suggestion_for_string(q, max_hits=10)
    print('Suggestions: ', str(suggested_muts))

    return json.dumps(suggested_muts)



@app.get("/get_session_id", dependencies=[Depends(cookie)])
async def get_session_id(session_data: SessionData = Depends(verifier)):
    """
    Query session ID stored in the covid cookie
    # A Header
    **bold**
    *italic*
    ## Smaller Header
   `` `
   a multiline code block
   `` `
   `a inline codeblock` Lorem Ipsum
    """

    return {'session_id': session_data.session_id}

@app.get("/get_mutation_information_table", dependencies=[Depends(cookie)])
async def get_mutation_information_table(session_data: SessionData = Depends(verifier)):
    """
   Get the stored and preprocessed mutation data table.
    ## Smaller Header
   `` `
   a multiline code block
   `` `
   `a inline codeblock` Lorem Ipsum
    """
    print('Query the mutation table if available')
    try:
        act_mutation_table = CovidOutcome.sessions[session_data.session_id]['sample_mutations']
    except KeyError as KE:
        raise HTTPException(
            status_code=500,
            detail=str('Key error. The session or the mutation table does not exist for that session, please recalculate first.'),
        )

    if act_mutation_table is None:
        act_mutation_table = pd.DataFrame()
        raise HTTPException(
            status_code=500,
            detail=str('The mutation table does not exist for that session, please recalculate first.'),
        )
    return act_mutation_table.to_json(orient="records")

@app.get("/get_filtered_mutation_information_table", dependencies=[Depends(cookie)])
async def get_filtered_mutation_information_table(session_data: SessionData = Depends(verifier)):
    '''
    Get the filtered mutation table.
    '''

    act_mutation_table = CovidOutcome.sessions[session_data.session_id]['sample_mutations']
    act_sequence_data =  CovidOutcome.sessions[session_data.session_id]['input_sequence_qc_table']

    mod_seq_table = act_sequence_data[['sequence_id', 'user_seq_id']]

    new_list = [
        'sequence_id',
        'user_seq_id',
        'ID',
        'HGVS.p',
        'mutation_name',
        'protein_name',
        'Variation_type',
        'Transcript biotype',
        'mutation_length',
        'Effect',
        'Putative impact',
        'Gene Name',
        'cDNA_position / cDNA_len',
        'CDS_position / CDS_len',
        'Protein_position / Protein_len', 'Distance to feature',
        'note'
    ]
    final_mutation_table = act_mutation_table.merge(mod_seq_table, how='left', left_on='gisaid_epi_isl', right_on='sequence_id')[new_list]

    return final_mutation_table.to_json(orient="records")


@app.get("/mutation_call_status", dependencies=[Depends(cookie)])
async def get_mutation_call_status(session_data: SessionData = Depends(verifier)):
    """
    Get mutation call status!
    """
    session_id = session_data.session_id
    mutation_calling_status_flag = None
    status_calc_var = None
    running_time = None
    timout_limit = CovidOutcome.CONFIG['MutationCalling']['tim_out_sec']
    print('Time out limit: ' + str(timout_limit))
    try:
        [mutation_calling_status_flag, status_calc_var, running_time] = CovidOutcome.get_mutation_call_status(
            session_id)
        if mutation_calling_status_flag != 'ERROR' and (mutation_calling_status_flag != 'OK') and running_time > timout_limit:
            print('Time out error!')
            raise HTTPException(
                status_code=500,
                detail=str(
                    'Time out error in mutation calling process for session {0}.\n  timeout limit: {1}. \n  running_time: '.format(
                        session_id,
                        timout_limit, running_time)),
            )
        if mutation_calling_status_flag == 'ERROR':
            print('Error happened during the mutation calculation process in session: {0}. Please retry.'.format(session_id))
            raise HTTPException(
                status_code=500,
                detail=str('Error happened during the mutation calculation process.'),
            )
    except KeyError as KE:
        raise HTTPException(
            status_code=500,
            detail=str(KE),
        )

    return {'mutation_calling_status_flag': mutation_calling_status_flag,
            'status_calc_var': status_calc_var,
            'mutation_call_elapsed_time': running_time}

@app.get("/get_sequence_information_table", dependencies=[Depends(cookie)])
async def get_sequence_information_table(session_data: SessionData = Depends(verifier)):
    print('Get the sequence information for a session.')
    session_id = session_data.session_id
    sequence_info = CovidOutcome.sessions[session_id]['input_sequence_qc_table']

    return {'sequence_info': sequence_info.to_json(orient="records")}

@app.get("/get_ml_results_table", dependencies=[Depends(cookie)])
async def get_ml_results_table(session_data: SessionData = Depends(verifier)):
    print('Get the Ml Results table!')
    session_id = session_data.session_id

    try:
        ml_results = CovidOutcome.get_ml_results(session_id)
    except ValueError:
        raise HTTPException(
            status_code=500,
            detail=str('Error in quering the ML results. Please check the output of the ML results, or run it first.'),
        )
    return ml_results.to_json(orient="records")

@app.get("/get_qc_controlled_sequence_file", dependencies=[Depends(cookie)])
async def get_qc_controlled_sequence_file(session_data: SessionData = Depends(verifier)):
    '''
    The output is in fasta format. The extension should be .fasta
    '''
    print('Returning with the QC controlled file if exists.')
    session_id = session_data.session_id
    qc_controlled_file = CovidOutcome.sessions[session_id]['outout_congr_cleaned_fasta_path']

    if not (qc_controlled_file is not None and os.path.isfile(qc_controlled_file)):
        raise HTTPException(
            status_code=500,
            detail="FileNotFound. " + "The quality controlled sequence file does not exist",
        )
    return FileResponse(qc_controlled_file)


@app.get("/get_annotated_vcf_file", dependencies=[Depends(cookie)])
async def get_annotated_vcf_file(session_data: SessionData = Depends(verifier)):
    print('Returning with the annotated VCF file if exists.')
    session_id = session_data.session_id
    print('Checking mutation calculation prereqestics')
    [ml_prereq, ml_prereq_qc_status, ml_prereq_mutation_status,
     ml_prereq_mutation_results] = CovidOutcome.evaluate_ml_model_check_prerequisite(session_id)
    if not ml_prereq:
        raise HTTPException(
            status_code=500,
            detail="Mutation calculation went wrong. " + "The output does not exist",
        )
    return FileResponse(CovidOutcome.sessions[session_id]['mutation_annotated_vcf_file'])

@app.get("/get_annotated_mutations_file", dependencies=[Depends(cookie)])
async def get_annotated_mutations_file(session_data: SessionData = Depends(verifier)):
    print('Returning with the filtered and annotated mutation file.')
    session_id = session_data.session_id
    print('Checking mutation calculation prereqestics')
    [ml_prereq, ml_prereq_qc_status, ml_prereq_mutation_status,
     ml_prereq_mutation_results] = CovidOutcome.evaluate_ml_model_check_prerequisite(session_id)
    if not ml_prereq:
        raise HTTPException(
            status_code=500,
            detail="Mutation calculation went wrong. " + "The output does not exist",
        )
    return FileResponse(CovidOutcome.sessions[session_id]['mutations_samples'])

@app.get("/get_protein_bar_plot_data", dependencies=[Depends(cookie)])
async def get_protein_bar_plot_data(session_data: SessionData = Depends(verifier)):
    '''
    Returns with the dataframe containing the data necessary to visualize the protein bar plot.

    '''
    print('Returning get_protein_bar_plot_data if exists.')
    session_id = session_data.session_id
    try:
        data_for_protein_bar_plot = CovidOutcome.get_protein_bar_plot_data(CovidOutcome.sessions[session_id]['sample_mutations'])
    except CovidOutcome.MutationDataException as err:
        raise HTTPException(
            status_code=404,
            detail=err.message,
        )
    except:
        raise HTTPException(
            status_code=500,
            detail='Something went wrong in the protein bar plot data generation. Are the mutations calculated successfully?',
        )

    records = []
    for  index, row  in data_for_protein_bar_plot.iterrows():
        record_data = {
            'proteinName': row.Protein,
            'mutations': row.Mutations,
            'occurrences': row.Occurrences
        }
        records.append(ProteinBarPlotDataRecord(**record_data))
    protein_bar_plot_data = {
        'records': records
    }
    proteinBarPlotData = ProteinBarPlotData(**protein_bar_plot_data)

    return proteinBarPlotData

@app.get("/get_mutation_histogram_plot_data", dependencies=[Depends(cookie)])
async def get_mutation_histogram_plot_data(session_data: SessionData = Depends(verifier)):
    '''
    Returns with the dataframe containing the data necessary to visualize the mutations histogram plot.

    '''
    print('Returning with the get_mutation_histogram_plot_data if exists.')
    session_id = session_data.session_id
    try:
        data_for_mutation_histogram = CovidOutcome.get_mutation_histogram_plot_data(CovidOutcome.sessions[session_id]['sample_mutations'])
    except CovidOutcome.MutationDataException as err:
        raise HTTPException(
            status_code=404,
            detail=err.message,
        )
    except:
        raise HTTPException(
            status_code=500,
            detail='Something went wrong in the mutation histogram data generation. Are the mutations calculated successfully?',
        )


    data_for_mutation_histogram.sort_values(by='Start_POS', inplace=True)
    print(data_for_mutation_histogram)
    records = []
    for index, row in data_for_mutation_histogram.iterrows():
        record_data = {
            'proteinName': row['Protein'],
            'mutation': row['Mutation'],
            'occurrences': row['Occurrences'],
            'startPos': row['Start_POS'],
            'endPos': row['End_POS']
        }
        records.append(MutationHistogramPlotDataRecord(**record_data))
    mutation_historgam_plot_data = {
        'records': records
    }
    mutationHistogramPlotData = MutationHistogramPlotData(**mutation_historgam_plot_data)

    return mutationHistogramPlotData

@app.get("/get_download_all_files_zipped", dependencies=[Depends(cookie)])
async def get_download_all_files_zipped(session_data: SessionData = Depends(verifier)):
    '''
    Returns with the dataframe containing the data necessary to visualize the ML results.

    '''
    print('get_download_all_files_zipped')
    session_id = session_data.session_id
    zipped_file_path = CovidOutcome.get_all_zipped_files(session_id)
    return FileResponse(zipped_file_path)

@app.get("/get_prediction_histogram_data", dependencies=[Depends(cookie)])
async def get_prediction_histogram_data(session_data: SessionData = Depends(verifier)):
    '''
    Returns with the dataframe containing the data necessary to visualize the ML results.

    '''
    print('Returning with the get_prediction_histogram_data if exists.')
    session_id = session_data.session_id
    try:
        data_for_prediction_histogram = CovidOutcome.get_prediction_histogram_plot_data(CovidOutcome.get_ml_results(session_id))
    except CovidOutcome.MLResultsDataException as err:
        raise HTTPException(
            status_code=404,
            detail=err.message,
        )
    except ValueError as err:
        raise HTTPException(
            status_code=404,
            detail="No available ML calculation",
        )
    except:
        raise HTTPException(
            status_code=500,
            detail='Something went wrong in the prediction histogram data generation. Are the predictions calculated successfully?',
        )

    records = []
    #print(data_for_prediction_histogram)

    for index, row in data_for_prediction_histogram.iterrows():
        record_data = {
            'label': row['Label'],
            'occurrences': row['Occurrences'],
            'predictionScore': row['Prediction_score'],
            'color': row['Color'],
            'sampleIds': row['SampleIDs']
        }
        records.append(PredictionHistogramPlotDataRecord(**record_data))
    prediction_historgam_plot_data = {
        'records': records
    }
    predictionHistogramPlotData = PredictionHistogramPlotData(**prediction_historgam_plot_data)

    return predictionHistogramPlotData

@app.post('/start_session')
async def create_session(response: Response):
    print('Creating a session!')
    session = uuid4()
    data = SessionData(session_id=str(session))
    a = await backend.create(session, data)
    cookie.attach_to_response(response, session)
    new_session_id = CovidOutcome.add_new_session(session)

    return {'session_id': new_session_id}

@app.post('/upload_and_qc_check_char_sequence_data', dependencies=[Depends(cookie)])
async def upload_char_sequence_data(covidsequencedata: CovidSequenceData, background_tasks: BackgroundTasks,
                               session_data: SessionData = Depends(verifier)):
    session_id = session_data.session_id

    number_of_headers =covidsequencedata.sequences.count('>')
    number_of_lines = len(covidsequencedata.sequences.split())
    to_log = """
    Some statistics: \nLenght of the sequence data: {0} 
    number of headers:  {1},
    number of lines: {2}
    """.format(len(covidsequencedata.sequences), number_of_headers,number_of_lines)
    #print(to_log)
    QC_RESULTS_FLAG = 'FAILED'
    if len(covidsequencedata.sequences) > CovidOutcome.CONFIG['MaxSequenceSequenceChars']:
        raise HTTPException(
            status_code=500,
            detail="Too long input sequence(s). The maximum length is %d, try smaller seqs." % CovidOutcome.CONFIG[
                'MaxSequenceSequenceChars'],
            headers={"SequenceInputError": "Too long input sequence(s)"},
        )
    try:
        CovidOutcome.set_session_data_null(session_id, covidsequencedata.sequences)
        act_session_data = CovidOutcome.sessions[session_id]
        act_session_data['raw_sequence_data'] = covidsequencedata.sequences
    except KeyError as KE:
        raise HTTPException(
            status_code=500,
            detail="SessionError" + "Invalid session, unregistered session. Reload the page." + str(KE),
        )
    try:
        sequence_response = CovidOutcome.check_input_sequences(session_id)
        background_tasks.add_task(background_calculation_of_mutations, session_id)
        N_valid_seqs = CovidOutcome.get_seqqc_number_of_valid_sequences(session_id)
        if N_valid_seqs > 0:
            QC_RESULTS_FLAG = 'OK'


    except ValueError as VE:
        CovidOutcome.check_input_sequence_set_status(session_id, 'Failed')
        raise HTTPException(
            status_code=500,
            detail="SequenceInputError" + " Invalid input sequence data." + str(VE),
        )
    return {'session_id': session_id,
            'sequence_profile': sequence_response.to_json(orient="records"),
            'QC_RESULTS_FLAG': QC_RESULTS_FLAG}


@app.post('/upload_and_qc_check_file_sequence_data', dependencies=[Depends(cookie)])
async def upload_sequence_file(background_tasks: BackgroundTasks, file: UploadFile = File(...), session_data: SessionData = Depends(verifier)):
    print('Uploading sequence file!')
    session_id = session_data.session_id
    filename = CovidOutcome.sessions[session_id]['raw_sequence_file']
    byteio = io.BytesIO()
    max_file_size = CovidOutcome.CONFIG['MaxFileSize']
    act_file_size = 0

    async with aiofiles.open(filename, 'wb') as out_file:
        while content := await file.read(10240):  # async read chunk
            act_file_size += 10240
            byteio.write(content)
            await out_file.write(content)
            if act_file_size > max_file_size:
                print('Too large file, interrupt')
                raise HTTPException(
                    status_code=500,
                    detail="Too large input sequence file. The maximum allowed file size is %d." %
                           CovidOutcome.CONFIG['MaxFileSize'],
                )
    raw_bytes = byteio.getvalue()
    encoding_info = chardet.detect(raw_bytes)
    #print(chardet.detect(raw_bytes))
    raw_sequence_chars = raw_bytes.decode(encoding_info['encoding'])
    covid_sequence = CovidSequenceData(sequences=raw_sequence_chars)
    upload_results = await upload_char_sequence_data(covid_sequence, background_tasks, session_data)
    return upload_results


@app.post('/create_artificial_covid_genome_data', dependencies=[Depends(cookie)])
async def create_artificial_covid_genome_data(covid_genomes: ArtificialCovidGenomes, session_data: SessionData = Depends(verifier)):
    '''
    Uploading a set of artificial covid genomes. No limits regarding the samples as of now.

    '''

    print('Creating artifical genome data')
    session_id = session_data.session_id
    updated_normalized_mutlists, empty_genomes_row_id, unmapped_mutations = CovidOutcome.build_artificial_genome_data(session_id, covid_genomes.genomes)

    return {'empty_genomes_row_id': empty_genomes_row_id,
            'unmapped_mutations': unmapped_mutations,
            'accepted_mapped_mutations': updated_normalized_mutlists}





@app.post('/mutation_call', dependencies=[Depends(cookie)])
async def mutation_calling(session_data: SessionData = Depends(verifier)):
    session_id = session_data.session_id
    print(session_id)
    try:
        CovidOutcome.mutation_calling_procedure(session_id)
    except ValueError as VE:
        raise HTTPException(
            status_code=500,
            detail="Mutation calling error!" + "Check the logs! " + str(VE),
        )


@app.post('/upload_age_file', dependencies=[Depends(cookie)])
async def upload_age_file(file: UploadFile = File(...), session_data: SessionData = Depends(verifier)):
    print('Uploading age data file!')
    session_id = session_data.session_id
    byteio = io.BytesIO()
    max_file_size = CovidOutcome.CONFIG['MaxFileSize']
    act_file_size = 0
    while content := await file.read(10240):  # async read chunk
        act_file_size += 10240
        byteio.write(content)
        if act_file_size > max_file_size:
            print('Too large file, interrupt')
            raise HTTPException(
                status_code=500,
                detail="Too large input sequence file. The maximum allowed file size is %d." %
                       CovidOutcome.CONFIG['MaxFileSize'],
                headers={"SequenceInputError": "Too large file ..."},
            )
    print('Checking filetype!')
    act_filename = str(file.filename)
    if act_filename.endswith('xlsx') or act_filename.endswith('xls'):
        print('This is probably an excel file. Trying to read it with pandas, looking for a column named age')
        data = pd.read_excel(byteio, engine='openpyxl')
        age_data = data['age']
        age_string = '\n'.join(list([str(age) for age in age_data]))
        age_data = AgeData(age_list=age_string)
        age_results = await set_sample_ages_char(age_data, session_data)
    else:
        print('It is probably a text file. ')
        raw_bytes = byteio.getvalue()
        encoding_info = chardet.detect(raw_bytes)
        print(chardet.detect(raw_bytes))
        raw_age_data_chars = raw_bytes.decode(encoding_info['encoding'])
        age_data = AgeData(age_list=raw_age_data_chars)
        age_results = await set_sample_ages_char(age_data, session_data)
    return age_results


@app.post('/set_sample_ages_char', dependencies=[Depends(cookie)])
async def set_sample_ages_char(age_data: AgeData, session_data: SessionData = Depends(verifier)):
    session_id = session_data.session_id
    print('Adding age info')
    if ';' or ',' in age_data.age_list:
        pass
    age_list = age_data.age_list.strip().splitlines()
    age_info = CovidOutcome.adding_age_information(session_id, age_list)
    return {"age": age_info.to_json(orient="records")}


@app.post('/generate_age_excel_form', dependencies=[Depends(cookie)])
async def generate_age_excel_form(session_data: SessionData = Depends(verifier)):
    session_id = session_data.session_id
    sequence_info = CovidOutcome.sessions[session_id]['input_sequence_qc_table']
    return {'filename': 'random_file',
            'form_data': sequence_info}


@app.post("/run_ml", dependencies=[Depends(cookie)])
async def run_ml(machine_learning_methods_type: Optional[str] = 'deep', session_data: SessionData = Depends(verifier)):
    '''
    Evaluate ML model for the data
    '''
    import traceback
    print('Evaluating the machine learning model on the provided data with the ml type %s. ' % machine_learning_methods_type)
    session_id = session_data.session_id

    if machine_learning_methods_type == 'automl':
        ml_model_type = 'automl'
    elif machine_learning_methods_type == 'deep':
        ml_model_type = 'deep'
    else:
        ml_model_type = 'automl'

    try:
        ml_results = CovidOutcome.evaluate_ml_model(session_id, ml_model_type)
    except ValueError as VE:
        print('Value error!!!')
        traceback.print_exc()

        raise HTTPException(
            status_code=500,
            detail=str(VE)
        )

    return (ml_results.to_json(orient="records"))



app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["GET", "POST", "OPTIONS", "PUT"],
    allow_headers=["*"],
)
