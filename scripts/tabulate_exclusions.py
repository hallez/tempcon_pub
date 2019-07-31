import numpy
import pandas
import os.path
import shutil
import yaml

if __name__ == '__main__':
    config = yaml.load(open('config.yml', 'r'))
    directories = config['directories']

    base_dir = directories['base_dir']

    info_fname = os.path.join(base_dir, 'subject_exclusion_info.csv')
    if not os.path.exists(info_fname):
        print(info_fname + " does not exist!")

    # read in using pandas because it plays nicely with csv files
    df = pandas.read_csv(info_fname)

    # we need to convert the csv object into a numeric so that we can deal with it
    # this gets split up into multiple steps - first, removing the comma
    df['exclude_runs_split'] = df.exclude_runs.str.split(',')

    # count the number of excluded runs
    # this for loop feels like a hack. seems like should be able to use "series" or "apply" instead
    tmp = None
    for r in range(len(df.exclude_runs_split)):
        print('--- working on row %d') % r
        cur_val = df.exclude_runs_split[r]

        if cur_val != ['0']:
            tmp_count = len(cur_val)
        else:
            tmp_count = 0
        print('%d row has %s excluded runs, which is length %s') % (r, cur_val, tmp_count)

        # add these values into tmp for adding to dataframe later
        if not bool(tmp):
            tmp = tmp_count
        else:
            if isinstance(tmp, int):
                tmp = [tmp, tmp_count]
            else:
                tmp.append(tmp_count)

    if len(df) == len(tmp):
        df['exclude_runs_count'] = tmp

    # create a mask for whether subjects had any excluded runs and whether they were excluded
    subj_w_excl_runs = df['exclude_runs_count'] > 0
    excl_subj = df['exclude_subject'] == 1

    print('\nacross %d subjects, %.02f runs [SD = %.02f] were excluded') % (len(df), df.exclude_runs_count.mean(),\
                                                                            df.exclude_runs_count.std())

    print('\nsplit up by excluded (1) and included (0) subjects')
    print(df.groupby('exclude_subject', as_index=False)['exclude_runs_count'].agg([numpy.mean, numpy.std, numpy.amin,\
                                                                                   numpy.amax, numpy.size]))
    print('\n%d out of %d included subjects had excluded runs') %\
         (df.exclude_runs_count.loc[subj_w_excl_runs][~excl_subj].size,\
                                                        len(df[~excl_subj]))