import numpy
import pandas
import os.path
import shutil
import yaml

# define functions (probably this is a later step once the script works)
def extract_and_save_cols(out_dir, subj_str_value, row_value, column_to_extract):
    fname_out = "exclude_runs.txt"
    subj_dir_out = os.path.join(out_dir, subj_str_value)
    fpath_out = os.path.join(subj_dir_out, fname_out)
    extracted_col = row_value[column_to_extract]

    if not os.path.exists(subj_dir_out):
        print("making subject directory " + subj_dir_out)
        os.makedirs(subj_dir_out)

    # based on: https://stackoverflow.com/questions/5214578/python-print-string-to-text-file
    text_file = open(fpath_out, "w")
    text_file.write(extracted_col)
    text_file.close()

if __name__ == '__main__':
    config = yaml.load(open('config.yml', 'r'))
    directories = config['directories']
    base_dir = directories['base_dir']
    raw_behavioral_dir = directories['raw_behavioral']

    # read in the csv file with subject information
    info_fname = os.path.join(base_dir, 'subject_exclusion_info.csv')
    if not os.path.exists(info_fname):
        print(info_fname + " does not exist!")

    # read in using pandas because it plays nicely with csv files
    df = pandas.read_csv(info_fname)
    print(type(df))
    print("column names are: " + str(list(df)))

    # based on: https://stackoverflow.com/questions/16476924/how-to-iterate-over-rows-in-a-dataframe-in-pandas/16476974
    for index, row in df.iterrows():
        subj_str = row['subject']
        print("working on subject: " + subj_str)

        extract_and_save_cols(raw_behavioral_dir, subj_str, row, 'exclude_runs')
