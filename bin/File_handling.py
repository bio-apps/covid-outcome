#!/usr/bin/env python
# -*-coding: utf8 -*-
import zipfile
import os
'''
File handling parts
'''

def check_uploadable_files(filenames):
    '''
    checking the that the files are available, existing and such
    '''

    print('Checking files:')

    valid_files = []
    for filename in filenames:
        if os.path.exists(filename) and os.path.getsize(filename) > 0:
            valid_files.append(filename)

    return valid_files



def generate_zipped_files(filenames, output_file_name):
    '''
    Print generating zipped files
    '''

    print('Generating a zipped file')

    with zipfile.ZipFile(output_file_name, 'w') as zipF:
        for file in filenames:
            _, tail = os.path.split(file)
            zipF.write(file, arcname=tail, compress_type=zipfile.ZIP_DEFLATED)

    print('All the files were compressed successfully!')


