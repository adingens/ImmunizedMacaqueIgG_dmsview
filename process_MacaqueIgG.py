"""
Retrieve, process, and output data file for IDC561 and mAbs (1-18 and 2-12).

This script assumes the following:
- the dms data is in BG505 numbring. #NOTE - I do not think this is true! As is, this produces the correct visualized data, with input MAP data in HSB2 numbering already
- the pdb file is in HXB2 numbering.
"""

import pandas as pd


def site_long_form(HXB2_map, antibodies, output=None):
    'Create site metric long form data file'
    #site_metrics = {'mediansitefracsurvive': 'median'}
    site_metrics = {'mediansitediffsel': 'median'}
    print('Processing site data\n')
    # process the site metrics
    df = []
    for antibody in antibodies:
        for metric in site_metrics.keys():
            fname = (f'../results/diffsel/summary_{antibody}-{metric}.csv')
            print(f'reading ... {fname}')
            sitemetricdf = pd.read_csv(fname)
            sitemetricdf = pd.melt(sitemetricdf, id_vars='site',
                                   var_name='metric')
            sitemetricdf['metric'] = (sitemetricdf['metric']
                                      .apply(lambda x:
                                             f'{x} ({site_metrics[metric]} '
                                             'of reps)'))
            sitemetricdf['condition'] = antibody
            df.append(sitemetricdf)
    df = pd.concat(df).rename(columns={'site': 'label_site'})

    # process the site numbering
    HXB2_map['site'] = [x+1 for x in range(len(HXB2_map))]
    df = pd.merge(df, HXB2_map, on=["label_site"]).astype({'site': 'int32'})
    df = df.sort_values("site", ascending=True)

    # output
    if output:
        df.to_csv(output, index=False)
    print('\nDone processing site data\n')
    return df


def mut_long_form(HXB2_map, antibodies, output=None):
    mut_metrics = {'medianmutdiffsel':
                   'differential selection (median of reps)'}

    print('Processing mut data\n')
    # process the mut metrics
    df = []
    for antibody in antibodies:
        for mut_metric in mut_metrics.keys():
            fname = (f'../results/diffsel/summary_{antibody}-{mut_metric}.csv')
            print(f'reading ... {fname}')
            mutmetricdf = pd.read_csv(fname)
            mutmetricdf = pd.melt(mutmetricdf,
                                  id_vars=['site', 'wildtype', 'mutation'],
                                  var_name='metric')
            mutmetricdf['metric'] = mut_metrics[mut_metric]
            mutmetricdf['value'] = mutmetricdf['value'].round(3)
            mutmetricdf['condition'] = antibody
            df.append(mutmetricdf)
    df = pd.concat(df).rename(columns={'site': 'label_site'})

    # process the site numbering
    HXB2_map['site'] = [x+1 for x in range(len(HXB2_map))]
    df = pd.merge(df, HXB2_map, on=["label_site"]).astype({'site': 'int32'})
    df = df.sort_values("site", ascending=True)

    # output
    if output:
        df.to_csv(output, index=False)
    print('\nDone processing mut data\n')
    return df


def final(site, mut):
    '''Combine site and mut long form into dms-view format.'''
    print('\nCombining mut and site data\n')
    # process the mutation data
    mut['metric'] = mut['metric'].apply(lambda x: f'mut_{x}')
    mut = pd.pivot_table(mut, index=['site', 'label_site', 'condition',
                                     'wildtype', 'mutation'],
                         columns='metric', values='value').reset_index()

    # process the site data
    site['metric'] = site['metric'].apply(lambda x: f'site_{x}')
    site = pd.pivot_table(site, index=['site', 'label_site', 'condition'],
                          columns='metric', values='value').reset_index()

    # combine the two dataframes
    df = pd.merge(mut, site, on=['site', 'label_site', 'condition'])
    assert len(mut) == len(df)
    print('Done!')
    return df


def main():
    'Create site long form, mut long form, and final data file.'
    # input files
    protein_fname = '5fyl_trimer_renumber.pdb'
    # set up
    HXB2_fname = ('https://raw.githubusercontent.com/jbloomlab/'
                  'EnvsAntigenicAtlas/master/results/HXB2_numbering/'
                  'BG505_to_HXB2.csv')
    HXB2_map = (pd.read_csv(HXB2_fname).rename(columns={"new": "label_site"})
                [["original", "label_site"]]).sort_values(by='original')
    antibodies = ['33203', '33311', '33433', '34945', 'CG41', 'MA255']
    chains = ["G", "B"]

    # create the site metric long form datafile
    site = site_long_form(HXB2_map, antibodies)
    mut = mut_long_form(HXB2_map, antibodies)
    df = final(site, mut)
    df['protein_site'] = df['label_site']
    df['protein_chain'] = " ".join(chains)
    df.to_csv('process_ImmunoMacaque.csv', index=False)


if __name__ == '__main__':
    main()
