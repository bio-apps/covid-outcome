# Overview
[CovidOutcome](https://covidoutcome.bio-ml.com) is a SARS-CoV-2 mutation identification pipeline and a corresponding machine learning tool that is capable of predicting the outcome of the disease using the genome information and the age of a given patient.

Our research goal was to apply state-of-the-art machine learning techniques to reveal and predict such possible links between the mutation status and the outcome. The tool provides the followings:
  * Applying state-of-the art machine learning techniques, autoML and deep-learning models for predicting the expected clinical outcome with expected accuracy of 0.795 CI: [0.762, 0.831]
  * Predicting novel-unseen mutation combination effects based on the existing data using ML models
  * Providing fast and high quality annotated mutation calling for SARS-COV2 genomes in standard (VCF) format

The model is based on 67708 SARS-CoV-2 genomes and corresponding patient data from the GISAID database. After rigorous data cleaning and preprocessing the machine learning models were trained with not only the single nucleotide substitutions, indels but mutations affecting UTR regions as well. The training set was further stratified to time-periods and age groups.



# Documentation

What does CovidOutcome do? Check [CovidOutcome wiki](https://github.com/bio-apps/covid-outcome/wiki/CovidOutcome).

How to use CovidOutcome? Check [the documentation](https://github.com/bio-apps/covid-outcome/wiki/Documentation).

# Installation
The covidoutcome is available as a standalone version as well. 

## Install with conda
Clone the repository and install the dependencies using the conda package manager.

```
git clone https://github.com/bio-apps/covid-outcome.git
cd covid-outcome
conda env create -f environment_cov_env.yml
conda activate cov_env
```
Basic usage:
```
usage: covidoutcome.py input output
 ```
``` input ``` = the covid sequences in fasta format

``` output ``` = the output directory
 
 For a detailed description of the outputs, please see the [documentation](https://github.com/bio-apps/covid-outcome/wiki/Documentation). 
 
 
 Using the example:
 ```
 covid-outcome/bin/covidoutcome.py --config covid-outcome/data/PipelineConfig.yml \
  covid-server-app/data/test_data/samples_for_prediction_viz.fasta ./outcome
```


# Citation

Please kindly cite our paper to support further development:

```
[1]    Ádám Nagy, Balázs Ligeti, János Szebeni, Sándor Pongor, Balázs Győrffy, COVIDOUTCOME—estimating 
       COVID severity based on mutation signatures in the SARS-CoV-2 genome, Database, Volume 2021,
       2021, baab020, https://doi.org/10.1093/database/baab020
```

```
@article{10.1093/database/baab020,
    author = {Nagy, Ádám and Ligeti, Balázs and Szebeni, János and Pongor, Sándor and Győrffy, Balázs},
    title = "{COVIDOUTCOME—estimating COVID severity based on mutation signatures in the SARS-CoV-2 genome}",
    journal = {Database},
    volume = {2021},
    year = {2021},
    month = {06},
    issn = {1758-0463},
    doi = {10.1093/database/baab020},
    url = {https://doi.org/10.1093/database/baab020},
    note = {baab020},
    eprint = {https://academic.oup.com/database/article-pdf/doi/10.1093/database/baab020/38828465/baab020.pdf},
}
```

