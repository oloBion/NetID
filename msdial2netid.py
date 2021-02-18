#!/usr/bin/env python

import argparse
import csv
import itertools
import pandas as pd
import pathlib
import tqdm


def parse_alignment_table(filename: str):
    N_HEADER_ROWS = 4

    # load sample properties
    with open(filename) as fin:
        # read the header from the MS-DIAL output
        reader = csv.reader(fin, delimiter='\t', quotechar='"')
        data = list(itertools.islice(reader, N_HEADER_ROWS + 1))
        
        # count the number of metadata columns
        metadata_col_count = len(list(itertools.takewhile(lambda x: x == '', data[0])))
        
        # drop the unwanted metadata
        data = [x[metadata_col_count:] for x in data][::-1]
        
        # extract column and row names
        sample_names = data.pop(0)[1:]
        property_names = [x[0] for x in data]
        data = [x[1:] for x in data]
        
        # build data frame and drop statistics columns
        samples = pd.DataFrame(data, columns=sample_names, index=property_names).transpose()
        samples = samples[samples.Class != 'NA']
        
        # rename duplicate columns
        samples.columns = pd.io.parsers.ParserBase({'names': samples.columns})._maybe_dedup_names(samples.columns)


    # load alignment table
    df = pd.read_csv(filename, delimiter='\t', skiprows=N_HEADER_ROWS)

    print(f'Loaded {filename}: {df.shape}')

    return samples, df


def msdial2netid(alignment_table_file, output_directory, msms_per_excel=100):
    output_directory.mkdir(exist_ok=True)
    (output_directory / 'msms').mkdir(exist_ok=True)

    samples, df = parse_alignment_table(alignment_table_file)

    netid_df = df[['Alignment ID', 'Average Mz', 'Average Rt(min)']]
    netid_df.columns = ['groupId', 'medMz', 'medRt']

    columns = [
        'label', 'metaGoupId', 'groupId', 'goodPeakCount', 'medMz',
        'medRt', 'maxQuality', 'isotopeLabel', 'compound', 'compoundId',
        'formula', 'isotopeLabel', 'expectedRtDiff', 'ppmDiff', 'parent'
    ]

    netid_df = netid_df.reindex(columns=columns)
    netid_df.to_csv(output_directory / 'raw_data.csv', index=False)


    # export MS/MS
    msms = df.dropna(subset=['MS/MS spectrum'])

    counter = 0
    metadata = []
    spectra = []

    for _, row in tqdm.tqdm(msms.iterrows(), total=len(msms)):
        metadata.append({
            'Mass_m_z_': row['Average Mz'],
            'Start_min_': row['Average Rt(min)'],
            'Comment': f'ID={row["Alignment ID"]}'
        })

        spectra.append(row['MS/MS spectrum'])
        counter += 1

        if counter > 0 and counter % msms_per_excel == 0 or counter == len(msms) - 1:
            writer = pd.ExcelWriter(output_directory / 'msms' / (alignment_table_file.stem + f'_msms_{counter // msms_per_excel:03d}.xlsx'))

            metadata = pd.DataFrame(metadata)
            metadata = metadata.reindex(columns=['Mass_m_z_', 'Formula_M_', 'FormulaType', 'Species', 'CS_z_', 'Polarity', 'Start_min_', 'End_min_', 'x_N_CEType', 'MSXID', 'Comment'])
            metadata.to_excel(writer, index=False)

            for i, s in enumerate(spectra):
                s = pd.DataFrame([x.split(':') for x in s.strip().split()])
                s.to_excel(writer, index=False, header=False, sheet_name=str(i + 2))

            writer.save()

            metadata = []
            spectra = []


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('alignment_table_file', type=pathlib.Path)
    parser.add_argument('output_directory', type=pathlib.Path)
    args = parser.parse_args()

    msdial2netid(args.alignment_table_file, args.output_directory)
