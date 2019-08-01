import sys
import os

print "os.system('pwd') is "
print os.system('pwd')

print "sys.path is: "
print sys.path

# need to add cwd to path (according to https://rcc.uchicago.edu/docs/tutorials/kicp-tutorials/running-jobs.html)
sys.path.append(os.getcwd())
import yaml
import os.path
import shutil
import zipfile
import os

if __name__ == '__main__':
    config = yaml.load(open('config.yml', 'r'))
    directories = config['directories']

    # define some variables that get re-used
    unzipped_mri_dir = directories['local_unzipped_dir']
    raw_mri_dir = directories['raw_mri']
    # range based on: http://stackoverflow.com/questions/18265935/python-create-list-with-numbers-between-2-values
    # adding string based on: http://stackoverflow.com/questions/2847386/python-string-and-integer-concatenation
    missing_zip_subjects = ['s' + `i` for i in range(1,9+1)]
    print 'missing_zip_subjects: [' + ', '.join(missing_zip_subjects) + ']'

    # loop across subject (s*) directories w/in raw_mri
    for d in os.listdir(raw_mri_dir):
        subj_dirname = os.path.relpath(d)

        print '--------'
        print 'Working on ' + subj_dirname + ':'
        if subj_dirname.startswith('s'):
        # check if the current subject's data needs to be copied from Maria Montchal's analyses
            if subj_dirname in missing_zip_subjects:
                print '- ' + subj_dirname + ' is contained in "missing_zip_subjects" - skipping.'
                continue

            cur_subj_raw_mri_dir = os.path.join(raw_mri_dir,subj_dirname)

            # figure out subject's ds*.zip filename
            for d in os.listdir(cur_subj_raw_mri_dir):
                if not d.endswith('.zip'):
                    print '- skipping non-zip file: "' + d + '"'
                    continue

                zip_file_name = d
                zip_file_path = os.path.join(cur_subj_raw_mri_dir,zip_file_name)
                print('- zip_file_path: "' + zip_file_path) + '"'

                # check to see if the subject already has an analyzed directory
                # if not, then create one
                subj_unzipped_mri_dir = os.path.join(unzipped_mri_dir,subj_dirname)
		print '- subj_unzipped_mri_dir is: " ' + subj_unzipped_mri_dir + '"'
                if not os.path.exists(subj_unzipped_mri_dir):
                    print '- creating directory: "' + subj_unzipped_mri_dir + '"'
                    os.makedirs(subj_unzipped_mri_dir)

                # loop through files in zip file
                # unzip into subj_unzipped_mri_dir if doesn't exist
                fh = open(zip_file_path,'r')
                z = zipfile.ZipFile(fh)

                # list all of the files in the *.zip file
                for f in z.namelist():
                    # limit to directories within the *.zip file
                    if not f.endswith('/'):
                        continue

                    zip_out_file = os.path.join(subj_unzipped_mri_dir,f)
                    # check to ensure
                    # that folder hasn't already been unzipped
                    if os.path.exists(zip_out_file) or os.path.isfile(zip_out_file):
                        print '- directory already exists: "' + zip_out_file + '" - skipping.'
                        continue

                    print '- extracting: "' + f + '" to: "' + subj_unzipped_mri_dir + '"'
                    z.extractall(subj_unzipped_mri_dir)

                fh.close()
