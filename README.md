# Overview
[CovidOutcome](https://covidoutcome.bio-ml.com) is a SARS-CoV-2 mutation identification pipeline and a corresponding machine learning tool that is capable of predicting the outcome of the disease using the genome information and the age of a given patient.

Our research goal was to apply state-of-the-art machine learning techniques to reveal and predict such possible links between the mutation status and the outcome. 

The model is based on 67708 SARS-CoV-2 genomes and corresponding patient data from the GISAID database. After rigorous data cleaning and preprocessing the machine learning models were trained with not only the single nucleotide substitutions, but mutations affecting UTR regions as well. The training set was further stratified to time-periods and age groups.

It also provides a prediction pipeline, where one can predict the outcome of the disease from a genome sample. The uploaded genome is analyzed and a prediction is made by one of the suitable models based on the user’s choice. Next to the prediction, we also output the found annotated mutations for the sample.

__TOC__


# Documentation

What does CovidOutcome do? Check [CovidOutcome wiki](https://github.com/bio-apps/covid-outcome/wiki/CovidOutcome)!

How to use CovidOutcome? Check [our documentation](https://github.com/bio-apps/covid-outcome/wiki/Documentation)!

# Installation
The covidoutcome available as a standalone version as well. 

## Install using conda
Clone the repository and install the dependencies using the conda package manager.

```
git clone https://github.com/bio-apps/covid-outcome.git
cd covid-outcome
conda env create -f environment_cov_env.yml
conda activate cov_env
```
The usage:
```
usage: covidoutcome.py [-h] [--config_file CONFIG_FILE] [--age data AGE DATA] input output
 ```
 Using the example:
 ```
 covid-outcome/bin/covidoutcome.py --config covid-outcome/data/PipelineConfig.yml \
  covid-server-app/data/test_data/samples_for_prediction_viz.fasta ./outcome
```

# Citation

Please kindly cite our papers to support further development:
```
[1]    Nagy, Á., Ligeti, B., Szebeni, J., Pongor, S., & Gyrffy, B. (2021).
       COVIDOUTCOME—estimating COVID severity based on mutation signatures in 
       the SARS-CoV-2 genome. Database, 2021. 
```

