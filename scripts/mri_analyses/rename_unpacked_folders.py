import yaml
import os.path
import shutil

def move_if_exists(source, destination):
    if os.path.exists(destination):
        print 'Skipping move of ' + source + ' to ' + destination + ' - destination already exists.'
        return

    if not os.path.exists(source):
        print 'Cannot move ' + source + ' to ' + destination + ' - source does not exist.'
        return

    print 'Moving ' + source + ' to ' + destination + '.'
    shutil.move(source, destination)


if __name__ == '__main__':
    config = yaml.load(open('config.yml', 'r'))
    directories = config['directories']

    default_analyzed_dirs_path = os.path.join(directories['mri_analysis_scripts'], 'default_dirs.yml')
    default_analyzed_dirs = yaml.load(open(default_analyzed_dirs_path, 'r'))


    dir_map = {
        'localizerDir': 'localizer',
        'mprageDir': 'mprage',
        't219Dir': 't2_19',
        't215Dir': 't2_15',
        'run1Dir': 'run1',
        'run2Dir': 'run2',
        'run3Dir': 'run3',
        'run4Dir': 'run4',
        'run5Dir': 'run5',
        'run6Dir': 'run6',
        'fieldmap1Dir': 'fieldmap1',
        'fieldmap2Dir': 'fieldmap2',
        'segmentedpartialDir': 'ep_seg_partial',
        'segmentedwholeDir': 'ep_seg_wholebrain'
    }

    for d in os.listdir(directories['analyzed_mri']):
        dirname = os.path.relpath(d)
        if dirname.startswith('s'):
            print("executing " + dirname)
            analyzed_dirs = default_analyzed_dirs.copy()

            # If there is a yaml file with overwrite directories,
            # we will merge anything that exists there with the existing defaults
            # The overwrite file will only contain entries that have changed,
            # so merging is necessary - otherwise we will clobber any defaults
            # that haven't changed.
            overwrite_analyzed_dirs_path = os.path.join(directories['raw_behavioral'], dirname, dirname + '.yml')
            if os.path.exists(overwrite_analyzed_dirs_path):
                print("found overwrite file " + overwrite_analyzed_dirs_path)
                overwrite_analyzed_dirs = yaml.load(open(overwrite_analyzed_dirs_path, 'r'))
                for k in overwrite_analyzed_dirs:
                    analyzed_dirs[k] = overwrite_analyzed_dirs[k]

            # Move all of the files if they exist, whether they were defaults or overrides
            for k in dir_map:
                move_if_exists(
                    os.path.join(directories['analyzed_mri'],dirname,analyzed_dirs[k]),
                    os.path.join(directories['analyzed_mri'],dirname,dir_map[k]))
