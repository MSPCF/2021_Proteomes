# %% IMPORTING PACKAGES =======================================================
# =============================================================================

# Importing supporting packages
import re

# Importing main packages
import numpy as np

# Importing Scikit Learn packages for incomming DBSCAN clustering
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import DBSCAN

# %% ==========================================================================
# =============================================================================


def start_finder(seq, fasta_dict={}):
    '''TODO'''
    if seq in fasta_dict['P1']:
        return fasta_dict['P1'].find(seq) + 1
    elif seq in fasta_dict['pre-P2']:
        return fasta_dict['pre-P2'].find(seq) + 1
    elif seq in fasta_dict['P2-2']:
        return fasta_dict['P2-2'].find(seq) + 1
    else:
        return fasta_dict['P2-3'].find(seq) + 1


def stop_finder(seq, fasta_dict={}):
    '''TODO'''
    if seq in fasta_dict['P1']:
        return fasta_dict['P1'].find(seq) + len(seq)
    elif seq in fasta_dict['pre-P2']:
        return fasta_dict['pre-P2'].find(seq) + len(seq)
    elif seq in fasta_dict['P2-2']:
        return fasta_dict['P2-2'].find(seq) + len(seq)
    else:
        return fasta_dict['P2-3'].find(seq) + len(seq)

# %% ==========================================================================
# =============================================================================


def PS_parser(PS_df):

    '''Parses ProSigth 'Modifications' column to get 'Targets', 'Positions',
    and 'Positions (abs.)'.'''

    # Filling nans in 'Modifications' column with empty string
    PS_df['Modifications'].fillna(value='', inplace=True)

    # Parsing 'Modifications' to get phosphorylation 'Targets' ([A-Z]+)
    ph_targets_λ = lambda c: re.findall(r'([A-Z]+)\d+\(.{0,3}\-phospho', c)
    PS_df['Targets'] = PS_df['Modifications'].apply(ph_targets_λ)

    # Parsing 'Modifications' to get phosphorylation 'Positions' (\d+)
    ph_positions_λ = lambda c: re.findall(r'[A-Z]+(\d+)\(.{0,3}\-phospho', c)
    PS_df['Positions'] = PS_df['Modifications'].apply(ph_positions_λ)

    # Converting 'Positions' list of strings into list of ints
    s_to_i_λ = lambda c: [int(p) for p in c]  # c --> cell p --> position
    PS_df['Positions'] = PS_df['Positions'].apply(s_to_i_λ)

    # Converting 'Positions' list into numpy arrays
    l_to_a_λ = lambda c: np.asarray(c)
    PS_df['Positions'] = PS_df['Positions'].apply(l_to_a_λ)

    # Converting relative 'Positions' into absolute 'Positions'
    PS_df['Positions (abs.)'] = PS_df['Positions'] + PS_df['Start'] - 1

    # Converting absolute and relative 'Positions' lists into numpy arrays
    a_to_l_λ = lambda c: list(c)
    PS_df['Positions (abs.)'] = PS_df['Positions (abs.)'].apply(a_to_l_λ)
    PS_df['Positions'] = PS_df['Positions'].apply(a_to_l_λ)

    # Add column with empty lists to uniformize 'Probabilities'
    PS_df['Probabilities'] = [[]] * len(PS_df.index)

    return PS_df

# %% ==========================================================================
# =============================================================================


def TP_parser(TP_df):

    '''Parses TopPic 'MIScore' column to get 'Targets', 'Positions',
    'Probabilities' and 'Positions (abs.)'.'''

    import re
    import numpy as np
    import pandas as pd

    # Transforms retention times from seconds to minuts
    TP_df['Retention time'] = TP_df['Retention time'] / 60

    # Parsing 'Protein accession'  to get the clenaed AC
    ac_regex = r'\|(.*)\|'
    ac_ser = TP_df['Protein accession'].apply(lambda i: re.findall(ac_regex, i))
    TP_df['Protein'] = ac_ser.str[0]
    del ac_regex, ac_ser

    # Filling '-' in 'MIScore' column  with empty string
    TP_df['MIScore'] = TP_df['MIScore'].replace({'-': ''}, regex=True)

    # Parsing 'MIScore' to get modification 'Targets'
    t_regex = r'\[(\w{1})'
    t_regex = '(\w{1})\d*\:'
    t_ser = TP_df['MIScore'].apply(lambda i: re.findall(t_regex, i))
    TP_df['Targets'] = t_ser
    del t_regex, t_ser

    # Parsing 'MIScore' to get modification 'Positions'
    p_regex = r'(\d.)\:'
    p_ser = TP_df['MIScore'].apply(lambda i: re.findall(p_regex, i))
    TP_df['Positions'] = p_ser
    del p_regex, p_ser

    # Converting 'Positions' list of strings into list of ints
    s_to_i_λ = lambda c: [int(p) for p in c]  # c --> cell p --> position
    TP_df['Positions'] = TP_df['Positions'].apply(s_to_i_λ)

    # Converting 'Positions' list into numpy arrays
    l_to_a_λ = lambda c: np.asarray(c)
    TP_df['Positions'] = TP_df['Positions'].apply(l_to_a_λ)

    # TODO
    μ = TP_df['Start'] == TP_df['First residue']

    # TODO
    TP_df.loc[μ, 'Positions (abs.)'] = TP_df.loc[μ, 'Positions']
    TP_df.loc[~μ, 'Positions (abs.)'] = TP_df.loc[~μ, 'Start'] + TP_df.loc[~μ, 'Positions'] -1

    # Converting 'Positions' arrays into numpy list
    a_to_l_λ = lambda c: list(c)
    TP_df['Positions'] = TP_df['Positions'].apply(a_to_l_λ)
    TP_df['Positions (abs.)'] = TP_df['Positions (abs.)'].apply(a_to_l_λ)

    # Parsing 'MIScore' to get modification 'Probabilities'
    prob_regex = r'\:(.*?)\%'
    prob_ser = TP_df['MIScore'].apply(lambda i: re.findall(prob_regex, i))
    TP_df['Probabilities'] = prob_ser
    del prob_regex, prob_ser

    # Converting 'Probabilities' list of strings into list of rounded floats
    s_to_f_λ = lambda c: [round(float(p), 1) for p in c]  # c -> cell p -> prob
    TP_df['Probabilities'] = TP_df['Probabilities'].apply(s_to_f_λ)

    return TP_df

# %% ==========================================================================
# =============================================================================


def col_renamer(df):
    '''TODO'''

    # Defining standardizing renaming dictionary
    uniform_dict = {
                    'Spectrum File': 'File',             # PS
                    'Data file name': 'File',            # TP

                    # 'Modifications': 'Modifications',  # PS
                    'MIScore': 'Modifications',          # TP

                    'MH [Da]': 'm (Da)',                 # PS
                    'Precursor mass': 'm (Da)',          # TP

                    'RT [min]': 'RT (min)',              # PS
                    'Retention time': 'RT (min)',        # TP

                    'Intensity': 'Int',                  # PS
                    'Feature intensity': 'Int',          # TP

                    'Fragmentation Scan(s)': 'Scan',     # PS
                    'Scan(s)': 'Scan',                   # TP

                    'Protein Descriptions': 'Protein',   # PS
                    # 'Protein': 'Protein',              # TP

                    'Annotated Sequence': 'Ann. Seq.',   # PS
                    'Proteoform': 'Ann. Seq.',           # TP

                    '-Log E-Value': '-log E',            # PS
                    # '-log E': '-log E',                # TP
                    }

    df.rename(columns=uniform_dict, inplace=True)

# %% ==========================================================================
# =============================================================================


# Defining supporting UDF for the incomming main Bokeh CallBack
def DBSCAN_clustering(data, eps, min_samples):
    '''Now we will perform the DBSCAN clustering. This time we only clustered
    the 'm (Da)', but the idea is to cluster both 'm (Da)' and 'RT (min)' in
    the future.
    '''

    # Standarizing & scaling the mass for the clustering
    scaler = StandardScaler()
    mass_scaled = scaler.fit_transform(data[['m (Da)']])

    # TODO
    scale = scaler.scale_[0]
    print(f'The clustering size is {(eps * scale) / 2} Da.')

    # Applying the DBSCAN clustering on the standarized & scaled masses
    clustering = DBSCAN(eps=eps, min_samples=min_samples,
                        # metric='euclidean', metric_params=None,
                        # algorithm='auto', leaf_size=30, p=None, n_jobs=1
                        ).fit_predict(mass_scaled)

    # Adding the clusters to our data
    data['Cluster'] = clustering

    # Preparing two lists to map cluster numbers with their cluster tags
    cluster_numbers = list(set(clustering))
    cluster_tags = ['P|' + str(number).zfill(3) for number in cluster_numbers]

    # Renaming unclustered (cluster numner -1) with tag 'Unclustered'
    cluster_tags = ['Unclustered' if t == 'P|-01' else t for t in cluster_tags]

    # Making a dictionary with the cluster numbers and the cluster tags
    cluster_dict = dict(zip(cluster_numbers, cluster_tags))

    # Replacing cluster numbers with our brand new cluster tags
    data['Cluster'] = data['Cluster'].map(cluster_dict)
