import pandas as pd
import numpy as np
import requests
import requests_ftp
import os
import re
import csv
from unidecode import unidecode


def EFO_parent_mapper(Cat_Stud, Cat_Anc_byN):
    """A function to map the parents and traits linkage file to the remaining
    catalogue datasets. It also does a bit of cleaning of terms to make strings
    which appear long on graphs in the notebook a bit cleaner
    """
    with open(os.path.abspath(os.path.join('__file__',
                                           '../..',
                                           'data',
                                           'Catalogue',
                                           'Synthetic',
                                           'Mapped_EFO.csv')), 'w') as fileout:
        efo_out = csv.writer(fileout, delimiter=',', lineterminator='\n')
        efo_out.writerow(['EFO URI', 'STUDY ACCESSION',
                          'PUBMEDID', 'ASSOCIATION COUNT'])
        for index, row in Cat_Stud.iterrows():
            listoftraits = row['MAPPED_TRAIT_URI'].split(',')
            for trait in listoftraits:
                efo_out.writerow([trait.lower().strip(),
                                  row['STUDY ACCESSION'],
                                  str(row['PUBMEDID']),
                                  str(row['ASSOCIATION COUNT'])])
    EFOsPerPaper = pd.read_csv(os.path.abspath(
                               os.path.join('__file__',
                                            '../..',
                                            'data',
                                            'Catalogue',
                                            'Synthetic',
                                            'Mapped_EFO.csv')), sep=',')
    EFO_Parent_Map = pd.read_csv(os.path.abspath(
                                 os.path.join('__file__',
                                              '../..',
                                              'data',
                                              'Catalogue',
                                              'Raw',
                                              'Cat_Map.tsv')), sep='\t')
    EFO_Parent_Map['EFO URI'] = EFO_Parent_Map['EFO URI'].str.lower(
    ).str.strip()
    EFO_Parent_Map = EFO_Parent_Map[[
        'EFO URI', 'Parent term', 'EFO term']].drop_duplicates()
    EFO_Parent_Paper_Merged = pd.merge(
        EFOsPerPaper, EFO_Parent_Map, on='EFO URI', how='left')
    Mapped = pd.merge(EFO_Parent_Paper_Merged, Cat_Anc_byN,
                      how='left', on='STUDY ACCESSION')
    Mapped = Mapped[~Mapped['N'].isnull(
    )]

    Mapped['EFO term'] = Mapped['EFO term'].str.replace('measurement', 'meas.')
    Mapped['EFO term'] = Mapped['EFO term'].str.title()
    Mapped['EFO term'] = Mapped['EFO term'].str.replace(
        'Hd Lc', 'HD LC')
    Mapped['EFO term'] = Mapped['EFO term'].str.replace(
        'Ld Lc', 'LD LC')
    Mapped['EFO term'] = Mapped['EFO term'].replace(
        'High Density Lipoprotein Cholesterol Meas.', 'HD LC measurement')
    Mapped['EFO term'] = Mapped['EFO term'].replace(
        'Low Density Lipoprotein Cholesterol Meas.', 'LD LC measurement')
    Mapped['EFO term'] = Mapped['EFO term'].str.replace(' Ii ', ' II ')
    Mapped['Parent term'] = Mapped['Parent term'].str.replace('measurement',
                                                              'meas.')
    Mapped['Parent term'] = Mapped['Parent term'].str.replace(' or ', '/')
    Mapped['Parent term'] = Mapped['Parent term'].str.title()

    with open(os.path.abspath(
              os.path.join('__file__',
                           '../..',
                           'data',
                           'Catalogue',
                           'Synthetic',
                           'Map_NoComs.csv')), 'w') as fileout:
        efo_out = csv.writer(fileout, delimiter=',', lineterminator='\n')
        efo_out.writerow(['EFO URI', 'STUDY ACCESSION',
                          'PUBMEDID', 'ASSOCIATION COUNT'])
        for index, row in Cat_Stud.iterrows():
            if ',' not in row['MAPPED_TRAIT_URI']:
                efo_out.writerow([row['MAPPED_TRAIT_URI'].lower().strip(),
                                  row['STUDY ACCESSION'],
                                  str(row['PUBMEDID']),
                                  str(row['ASSOCIATION COUNT'])])
    EFOsPer_NoComs = pd.read_csv(os.path.abspath(
                                 os.path.join('__file__', '../..', 'data',
                                              'Catalogue', 'Synthetic',
                                              'Map_NoComs.csv')), sep=',')
    EFO_Parent_Paper_Merged = pd.merge(
        EFOsPer_NoComs, EFO_Parent_Map, on='EFO URI', how='left')
    Map_NoComs = pd.merge(EFO_Parent_Paper_Merged, Cat_Anc_byN,
                          how='left', on='STUDY ACCESSION')
    Map_NoComs = Map_NoComs[~Map_NoComs['N'].isnull()]
    Map_NoComs['EFO term'] = Map_NoComs['EFO term'].str.replace('measurement',
                                                                'meas.')
    Map_NoComs['EFO term'] = Map_NoComs['EFO term'].str.title()
    Map_NoComs['EFO term'] = Map_NoComs['EFO term'].str.replace('Hd Lc',
                                                                'HD LC')
    Map_NoComs['EFO term'] = Map_NoComs['EFO term'].str.replace('Ld Lc',
                                                                'LD LC')
    Map_NoComs['EFO term'] = Map_NoComs['EFO term'].replace(
        'High Density Lipoprotein Cholesterol Meas.', 'HD LC measurement')
    Map_NoComs['EFO term'] = Map_NoComs['EFO term'].replace(
        'Low Density Lipoprotein Cholesterol Meas.', 'LD LC measurement')
    Map_NoComs['EFO term'] = Map_NoComs['EFO term'].str.replace(' Ii ', ' II ')
    Map_NoComs['Parent term'] = Map_NoComs['Parent term'].str.replace(
        'measurement', 'meas.')
    Map_NoComs['Parent term'] = Map_NoComs['Parent term'].str.replace(
        ' or ', '/').str.title()
    return Mapped, Map_NoComs


def load_gwas_cat():
    """A function which cleans and loads all the main GWAS Catalog files into
    the notebook workspace. It renames a couple of fields, and creates
    one groupedby.sum for return
    """
    Cat_Stud = pd.read_csv(os.path.abspath(
                           os.path.join('__file__',
                                        '../..',
                                        'data',
                                        'Catalogue',
                                        'Raw',
                                        'Cat_Stud.tsv')),
                           header=0, sep='\t', encoding='utf-8',
                           index_col=False)
    Cat_Stud.fillna('N/A', inplace=True)
    Cat_Anc = pd.read_csv(os.path.abspath(
                          os.path.join('__file__', '../..',
                                       'data',
                                       'Catalogue',
                                       'Raw',
                                       'Cat_Anc.tsv')),
                          header=0, sep='\t', encoding='utf-8',
                          index_col=False)
    Cat_Anc.rename(columns={'BROAD ANCESTRAL CATEGORY': 'BROAD ANCESTRAL',
                            'NUMBER OF INDIVDUALS': 'N'}, inplace=True)
    Cat_Anc = Cat_Anc[~Cat_Anc['BROAD ANCESTRAL'].isnull()]
    Cat_Anc.columns = Cat_Anc.columns.str.replace('ACCCESSION', 'ACCESSION')
    Cat_Anc_byN = Cat_Anc[['STUDY ACCESSION', 'N',
                           'DATE']].groupby(by='STUDY ACCESSION').sum()
    Cat_Anc_byN = Cat_Anc_byN.reset_index()
    Cat_Anc_byN = pd.merge(Cat_Anc_byN, Cat_Stud[[
        'STUDY ACCESSION', 'DATE']], how='left', on='STUDY ACCESSION')
    cleaner_broad = pd.read_csv(os.path.abspath(
        os.path.join('__file__',
                     '../..',
                     'data',
                     'Support',
                     'dict_replacer_broad.tsv')),
        sep='\t', header=0, index_col=False)
    Cat_Anc = pd.merge(Cat_Anc, cleaner_broad, how='left',
                       on='BROAD ANCESTRAL')
    Cat_Anc['Dates'] = [pd.to_datetime(d) for d in Cat_Anc['DATE']]
    Cat_Anc['N'] = pd.to_numeric(Cat_Anc['N'], errors='coerce')
    Cat_Anc = Cat_Anc[Cat_Anc['N'].notnull()]
    Cat_Anc['N'] = Cat_Anc['N'].astype(int)
    Cat_Anc = Cat_Anc.sort_values(by='Dates')
    Cat_Anc['Broader']
    Cat_Anc['Broader'] = Cat_Anc['Broader'].str.replace(
        'African American or Afro-Caribbean', 'African Am./Caribbean')
    Cat_Anc['Broader'] = Cat_Anc['Broader'].str.replace(
        'Hispanic or Latin American', 'Hispanic/Latin American')
    Cat_Full = pd.read_csv(os.path.abspath(os.path.join('__file__',
                                                        '../..',
                                                        'data',
                                                        'Catalogue',
                                                        'Raw',
                                                        'Cat_Full.tsv')),
                           header=0, sep='\t', encoding='utf-8',
                           index_col=False, low_memory=False)

    Cat_Anc.to_csv(os.path.abspath(
        os.path.join('__file__',
                     '../..',
                     'data',
                     'Catalogue',
                     'Synthetic',
                     'Cat_Anc_withBroader.tsv')),
                   sep='\t', index=False)
    return Cat_Stud, Cat_Anc, Cat_Anc_byN, Cat_Full


def load_pubmed_data():
    """" Load the PUBMED data into the workspace"""
    FunderInfo = pd.read_csv(os.path.abspath(
                             os.path.join('__file__',
                                          '../..',
                                          'data',
                                          'PUBMED',
                                          'Pubmed_FunderInfo.csv')))
    FunderInfo = FunderInfo[FunderInfo['Agency'].notnull()]
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('National Natural Science Foundation of China \(National Science Foundation of China\)',
                                       'National Natural Science Foundation of China')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('Medical Research Council', 'MRC')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('National Research Foundation of Korea \(KR\)',
                                       'National Research Foundation of Korea')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('National Human Genome Research Institute',
                                       'NHGRI')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('NHGRI NIH HHS',
                                       'NHGRI')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('NCI NIH HHS',
                                       'NCI')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('National Cancer Institute',
                                       'NCI')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('National Institute of Arthritis and Musculoskeletal and Skin Diseases',
                                       'NIAMS NIH HHS')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('National Institute on Aging',
                                       'NIA NIH HHS')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('NIDA NIH HHS',
                                       'NIDA')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('National Institute of Diabetes and Digestive and Kidney Diseases',
                                       'NIDDK NIH HHS')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('National Institutes of Health',
                                       'NIH HHS')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('NIMH NIH HHS',
                                       'NIMH')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('National Institute of Neurological Disorders and Stroke',
                                       'NINDS')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('NINDS NIH HHS',
                                       'NINDS')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('National Natural Science Foundation of China (National Science Foundation of China)',
                                       'National Natural Science Foundation of China')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('National Research Foundation of Korea (KR)',
                                       'National Research Foundation of Korea')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('National Human Genome Research Institute',
                                       'NHGRI')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('NHGRI NIH HHS',
                                       'NHGRI')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('National Heart, Lung, and Blood Institute',
                                       'NHLBI')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('NHLBI NIH HHS',
                                       'NHLBI')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('National Cancer Institute',
                                       'NCI')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('NCI NIH HHS',
                                       'NCI')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('National Institute of Arthritis and Musculoskeletal and Skin Diseases',
                                       'NIAMS NIH HHS')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('National Institute on Aging',
                                       'NIA NIH HHS')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('National Institute of Arthritis and Musculoskeletal and Skin Diseases',
                                       'NIAMS NIH HHS')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('NIDA NIH HHS',
                                       'NIDA')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('National Institute of Diabetes and Digestive and Kidney Diseases',
                                       'NIDDK NIH HHS')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('National Institutes of Health',
                                       'NIH HHS')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('NIMH NIH HHS',
                                       'NIMH')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('National Institute of Neurological Disorders and Stroke',
                                       'NINDS')
    FunderInfo['Agency'] = FunderInfo['Agency'].\
                           str.replace('NINDS NIH HHS',
                                       'NINDS')
    AbstractInfo = pd.read_csv(os.path.abspath(
                               os.path.join('__file__',
                                            '../..',
                                            'data',
                                            'PUBMED',
                                            'Pubmed_AbstractInfo.csv')))
    AbstractInfo['Abstracts'] = AbstractInfo['Abstract'].apply(
        lambda x: re.sub('[^a-zA-Z ]+', '', x)).str.lower()

    AbstractCount = pd.DataFrame(AbstractInfo.Abstract.apply(lambda x:
                                                             pd.value_counts(x.split(' '))).sum(axis=0),
                                 columns=['num_words'])
    AbstractCount.index.name = 'word'
    AbstractCount.to_csv(os.path.abspath(
                         os.path.join('__file__',
                                      '../..',
                                      'data',
                                      'PUBMED',
                                      'Pubmed_AbstractCount.csv')))

    AuthorMaster = pd.read_csv(os.path.abspath(
                               os.path.join('__file__',
                                            '../..',
                                            'data',
                                            'PUBMED',
                                            'Pubmed_AuthorInfo.csv')))
    AuthorMaster = AuthorMaster.drop_duplicates(
        subset=['PUBMEDID', 'AUTHORNAME'])
    for col in ['FORENAME', 'LASTNAME', 'AUTHORNAME']:
        AuthorMaster[col] = AuthorMaster[col].apply(unidecode)
    return FunderInfo, AbstractInfo, AuthorMaster


def make_timely(variables, yearlist, yearquarterlist, Cat_Stud,
                Cat_Anc, Cat_Anc_byN):
    """ build the dataframe which plots the longitudinal growth of GWAS over
    time (figure 1a)
    """
    df_years = pd.DataFrame(columns=variables, index=yearlist)
    df_yearquarters = pd.DataFrame(
        columns=variables, index=yearquarterlist)
    for year_ in range(2007, 2018):
        df_years['N ≤ 5,000'][str(year_)] = \
            len(Cat_Anc_byN[(Cat_Anc_byN['N'] <= 5000) &
                            (Cat_Anc_byN['DATE'].str.contains(str(year_)))])
        df_years['5,001 ≤ N ≤ 50,000'][str(year_)] = \
            len(Cat_Anc_byN[(Cat_Anc_byN['N'] >= 5001) &
                            (Cat_Anc_byN['N'] <= 50000) &
                            (Cat_Anc_byN['DATE'].str.contains(str(year_)))])
        df_years['50,001 ≤ N ≤ 100,000'][str(year_)] = \
            len(Cat_Anc_byN[(Cat_Anc_byN['N'] > 50000) &
                            (Cat_Anc_byN['N'] <= 100000) &
                            (Cat_Anc_byN['DATE'].str.contains(str(year_)))])
        df_years['100,001 ≤ N'][str(year_)] = \
            len(Cat_Anc_byN[(Cat_Anc_byN['N'] > 100000) &
                            (Cat_Anc_byN['DATE'].str.contains(str(year_)))])
        df_years['N'][str(year_)] = Cat_Anc[Cat_Anc['DATE'].str.contains(
            str(year_))]['N'].sum()
        df_years['Associations'][str(year_)] = \
            Cat_Stud[Cat_Stud['DATE'].str.contains(
                str(year_))]['ASSOCIATION COUNT'].sum()
        df_years['Journals Printing GWAS'][str(year_)] = \
            len(Cat_Stud[Cat_Stud['DATE'].str.contains(
                str(year_))]['JOURNAL'].unique())
        df_years['# Diseases Studied'][str(year_)] = \
            len(Cat_Stud[Cat_Stud['DATE'].str.contains(
                str(year_))]['DISEASE/TRAIT'].unique())
        df_yearquarters['N ≤ 5,000'][str(year_) + 'Q1'] = \
            len(Cat_Anc_byN[(Cat_Anc_byN['N'] <= 5000) &
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-01-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-02-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-03-'))])
        df_yearquarters['5,001 ≤ N ≤ 50,000'][str(year_) + 'Q1'] = \
            len(Cat_Anc_byN[(Cat_Anc_byN['N'] >= 5001) &
                            (Cat_Anc_byN['N'] <= 50000) &
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-01-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-02-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-03-'))])
        df_yearquarters['50,001 ≤ N ≤ 100,000'][str(year_) + 'Q1'] = \
            len(Cat_Anc_byN[(Cat_Anc_byN['N'] <= 5000) &
                            (Cat_Anc_byN['N'] <= 100000) &
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-01-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-02-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-03-'))])
        df_yearquarters['100,001 ≤ N'][str(year_) + 'Q1'] = \
            len(Cat_Anc_byN[(Cat_Anc_byN['N'] > 100000) &
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-01-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-02-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-03-'))])
        df_yearquarters['N ≤ 5,000'][str(year_) + 'Q2'] = \
            len(Cat_Anc_byN[(Cat_Anc_byN['N'] <= 5000) &
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-04-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-05-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-06-'))])
        df_yearquarters['5,001 ≤ N ≤ 50,000'][str(year_) + 'Q2'] = \
            len(Cat_Anc_byN[(Cat_Anc_byN['N'] >= 5001) &
                            (Cat_Anc_byN['N'] <= 50000) &
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-04-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-05-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-06-'))])
        df_yearquarters['50,001 ≤ N ≤ 100,000'][str(year_) + 'Q2'] = \
            len(Cat_Anc_byN[(Cat_Anc_byN['N'] <= 5000) &
                            (Cat_Anc_byN['N'] <= 100000) &
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-04-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-05-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-06-'))])
        df_yearquarters['100,001 ≤ N'][str(year_) + 'Q2'] = \
            len(Cat_Anc_byN[(Cat_Anc_byN['N'] > 100000) &
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-04-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-05-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-06-'))])
        df_yearquarters['N ≤ 5,000'][str(year_) + 'Q3'] = \
            len(Cat_Anc_byN[(Cat_Anc_byN['N'] <= 5000) &
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-07-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-08-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-09-'))])
        df_yearquarters['5,001 ≤ N ≤ 50,000'][str(year_) + 'Q3'] = \
            len(Cat_Anc_byN[(Cat_Anc_byN['N'] >= 5001) &
                            (Cat_Anc_byN['N'] <= 50000) &
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-07-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-08-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-09-'))])
        df_yearquarters['50,001 ≤ N ≤ 100,000'][str(year_) + 'Q3'] = \
            len(Cat_Anc_byN[(Cat_Anc_byN['N'] <= 5000) &
                            (Cat_Anc_byN['N'] <= 100000) &
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-07-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-08-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-09-'))])
        df_yearquarters['100,001 ≤ N'][str(year_) + 'Q3'] = \
            len(Cat_Anc_byN[(Cat_Anc_byN['N'] > 100000) &
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-07-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-08-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-09-'))])
        df_yearquarters['N ≤ 5,000'][str(year_) + 'Q4'] = \
            len(Cat_Anc_byN[(Cat_Anc_byN['N'] <= 5000) &
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-10-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-11-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-12-'))])
        df_yearquarters['5,001 ≤ N ≤ 50,000'][str(year_) + 'Q4'] = \
            len(Cat_Anc_byN[(Cat_Anc_byN['N'] >= 5001) &
                            (Cat_Anc_byN['N'] <= 50000) &
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-10-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-11-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-12-'))])
        df_yearquarters['50,001 ≤ N ≤ 100,000'][str(year_) + 'Q4'] = \
            len(Cat_Anc_byN[(Cat_Anc_byN['N'] <= 5000) &
                            (Cat_Anc_byN['N'] <= 100000) &
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-10-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-11-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-12-'))])
        df_yearquarters['100,001 ≤ N'][str(year_) + 'Q4'] = \
            len(Cat_Anc_byN[(Cat_Anc_byN['N'] > 100000) &
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-10-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-11-')) |
                            (Cat_Anc_byN['DATE'].str.contains(
                                str(year_) + '-12-'))])
    return df_years, df_yearquarters


def make_clean_CoR(Cat_Anc):
    """ clean the country of recruitment field for the geospatial analysis
    """
    with open(os.path.abspath(
              os.path.join('__file__',
                           '../..',
                           'data',
                           'Catalogue',
                           'Synthetic',
                           'ancestry_CoR.csv')), 'w') as fileout:
        rec_out = csv.writer(fileout, delimiter=',', lineterminator='\n')
        rec_out .writerow(['Date', 'PUBMEDID', 'N', 'Cleaned Country'])
        for index, row in Cat_Anc.iterrows():
            if len(row['COUNTRY OF RECRUITMENT'].split(',')) == 1:
                rec_out .writerow([row['DATE'],
                                   str(row['PUBMEDID']),
                                   str(row['N']),
                                   row['COUNTRY OF RECRUITMENT']])

    Clean_CoR = pd.read_csv(os.path.abspath(
        os.path.join('__file__',
                     '../..',
                     'data',
                     'Catalogue',
                     'Synthetic',
                     'ancestry_CoR.csv')))
    Clean_CoR['Cleaned Country'] = Clean_CoR['Cleaned Country'].str.replace(
        'U.S.', 'United States')
    Clean_CoR['Cleaned Country'] = Clean_CoR['Cleaned Country'].str.replace(
        'Gambia', 'Gambia, The')
    Clean_CoR['Cleaned Country'] = Clean_CoR['Cleaned Country'].str.replace(
        'U.K.', 'United Kingdom')
    Clean_CoR['Cleaned Country'] = Clean_CoR['Cleaned Country'].str.replace(
        'Republic of Korea', 'Korea, South')
    Clean_CoR['Cleaned Country'] = Clean_CoR['Cleaned Country'].str.replace(
        'Czech Republic', 'Czechia')
    Clean_CoR['Cleaned Country'] = Clean_CoR['Cleaned Country'].str.replace(
        'Russian Federation', 'Russia')
    Clean_CoR['Cleaned Country'] = Clean_CoR['Cleaned Country'].str.replace(
        'Iran \(Islamic Republic of\)', 'Iran')
    Clean_CoR['Cleaned Country'] = Clean_CoR['Cleaned Country'].str.replace(
        'Viet Nam', 'Vietnam')
    Clean_CoR['Cleaned Country'] = Clean_CoR['Cleaned Country'].str.replace(
        'United Republic of Tanzania', 'Tanzania')
    Clean_CoR['Cleaned Country'] = Clean_CoR['Cleaned Country'].str.replace(
        'Republic of Ireland', 'Ireland')
    Clean_CoR['Cleaned Country'] = Clean_CoR['Cleaned Country'].str.replace(
        'Micronesia \(Federated States of\)',
        'Micronesia, Federated States of')
    Clean_CoR = Clean_CoR[Clean_CoR['Cleaned Country'] != 'NR']
    print('Cleaning for single country of recruitment:\n' +
          str(round((len(Clean_CoR) / len(Cat_Anc)) * 100, 2)) +
          '% of the rows remain.')
    print(str(round((Clean_CoR['N'].sum() / Cat_Anc['N'].sum()) * 100, 2)) +
          '% of the N remains.')
    Clean_CoR.to_csv(os.path.abspath(
                     os.path.join('__file__',
                                  '../..',
                                  'data',
                                  'Catalogue',
                                  'Synthetic',
                                  'GWASCatalogue_CleanedCountry.tsv')),
                     sep='\t', index=False)
    return Clean_CoR


def download_cat(path, ebi_download):
    """ download the data from the ebi main site and ftp"""
    r = requests.get(ebi_download + 'studies_alternative')
    with open(os.path.join(path, 'Cat_Stud.tsv'), 'wb') as tsvfile:
        tsvfile.write(r.content)
    r = requests.get(ebi_download + 'ancestry')
    with open(os.path.join(path, 'Cat_Anc.tsv'), 'wb') as tsvfile:
        tsvfile.write(r.content)
    r = requests.get(ebi_download + 'full')
    with open(os.path.join(path, 'Cat_Full.tsv'), 'wb') as tsvfile:
        tsvfile.write(r.content)
    requests_ftp.monkeypatch_session()
    s = requests.Session()
    ftpsite = 'ftp://ftp.ebi.ac.uk/'
    subdom = '/pub/databases/gwas/releases/latest/'
    file = 'gwas-efo-trait-mappings.tsv'
    r = s.get(ftpsite+subdom+file)
    with open(os.path.join(path, 'Cat_Map.tsv'), 'wb') as tsvfile:
        tsvfile.write(r.content)
