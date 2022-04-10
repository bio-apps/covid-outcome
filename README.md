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

## Install using conda
Clone the repository

```
git clone https://github.com/bio-apps/covid-outcome.git
cd covid-outcome
conda env create -f environment_cov_env.yml
conda activate cov_env
```


# Citation

Please kindly cite our papers to support further development:

Nagy A. et al: [COVIDOUTCOME—estimating COVID severity based on mutation signatures in the SARS-CoV-2 genome Database](https://academic.oup.com/database/article-abstract/doi/10.1093/database/baab020/6272506) (2021)
