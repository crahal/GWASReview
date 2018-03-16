import pandas as pd
import os
import re


def EFO_parent_mapper(Cat_Stud, Cat_Anc_byN):
    """A function to map the parents and traits linkage file to the remaining
    catalogue datasets. It also does a bit of cleaning of terms to make strings
    which appear long on graphs in the notebook a bit cleaner
    """
    fileout = open(os.path.abspath(os.path.join('__file__', '../..', 'Data',
                                                'Catalogue', 'Synthetic',
                                                'Mapped_EFO.csv')), 'w')
    fileout.write('EFO URI,STUDY ACCESSION,PUBMEDID,ASSOCIATION COUNT\n')
    for index, row in Cat_Stud.iterrows():
        listoftraits = row['MAPPED_TRAIT_URI'].split(',')
        for trait in listoftraits:
            fileout.write(trait.lower().strip() + ',' + row['STUDY ACCESSION']
                          + ',' + str(row['PUBMEDID']) + ','
                          + str(row['ASSOCIATION COUNT']) + '\n')
    EFOsPerPaper = pd.read_csv(os.path.abspath(
                               os.path.join('__file__', '../..', 'Data',
                                            'Catalogue', 'Synthetic',
                                            'Mapped_EFO.csv')), sep=',')
    EFO_Parent_Map = pd.read_csv(os.path.abspath(
                                 os.path.join('__file__', '../..', 'Data',
                                              'Catalogue', 'Raw',
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

    fileout = open(os.path.abspath(
                   os.path.join('__file__', '../..', 'Data', 'Catalogue',
                                'Synthetic', 'Map_NoComs.csv')), 'w')
    fileout.write('EFO URI,STUDY ACCESSION,PUBMEDID,ASSOCIATION COUNT\n')
    for index, row in Cat_Stud.iterrows():
        if ',' not in row['MAPPED_TRAIT_URI']:
            fileout.write(row['MAPPED_TRAIT_URI'].lower().strip() + ',' +
                          row['STUDY ACCESSION'] + ',' + str(row['PUBMEDID']) +
                          ',' + str(row['ASSOCIATION COUNT']) + '\n')
    EFOsPer_NoComs = pd.read_csv(os.path.abspath(
                                 os.path.join('__file__', '../..', 'Data',
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
                           os.path.join('__file__', '../..', 'Data',
                                        'Catalogue', 'Raw',
                                        'Cat_Stud.tsv')),
                           header=0, sep='\t', encoding='utf-8',
                           index_col=False)
    Cat_Stud.fillna('N/A', inplace=True)
    Cat_Anc = pd.read_csv(os.path.abspath(
                          os.path.join('__file__', '../..',
                                       'Data', 'Catalogue', 'Raw',
                                       'Cat_Anc.tsv')),
                          header=0, sep='\t', encoding='utf-8',
                          index_col=False)
    Cat_Anc.rename(columns={'BROAD ANCESTRAL CATEGORY': 'BROAD ANCESTRAL',
                            'NUMBER OF INDIVDUALS': 'N'}, inplace=True)
    Cat_Anc.columns = Cat_Anc.columns.str.replace('ACCCESSION', 'ACCESSION')
    Cat_Anc_byN = Cat_Anc[['STUDY ACCESSION', 'N',
                           'DATE']].groupby(by='STUDY ACCESSION').sum()
    Cat_Anc_byN = Cat_Anc_byN.reset_index()
    Cat_Anc_byN = pd.merge(Cat_Anc_byN, Cat_Stud[[
        'STUDY ACCESSION', 'DATE']], how='left', on='STUDY ACCESSION')
    cleaner_broad = pd.read_excel(os.path.abspath(
                                  os.path.join('__file__', '../..', 'Data',
                                               'Support',
                                               'dict_replacer_broad.tsv')),
                                  sep='\t', header=0, index_col=False)
    Cat_Anc = pd.merge(Cat_Anc, cleaner_broad, how='left',
                       on='BROAD ANCESTRAL')
    Cat_Anc = Cat_Anc[Cat_Anc['N'].notnull()]  # jankys?
    Cat_Anc['Dates'] = [pd.to_datetime(d) for d in Cat_Anc['DATE']]
    Cat_Anc = Cat_Anc[Cat_Anc['Broader Ancestral'].notnull()]  # jankys?
    Cat_Anc = Cat_Anc.sort_values(by='Dates')
    Cat_Anc['Broader Ancestral'] = Cat_Anc['Broader Ancestral'].str.replace(
        'African American or Afro-Caribbean', 'African Am.\Caribbean')
    Cat_Anc['Broader Ancestral'] = Cat_Anc['Broader Ancestral'].str.replace(
        'Hispanic or Latin American', 'Hispanic\Latin American')
    Cat_Anc['BROAD ANCESTRAL'] = Cat_Anc['BROAD ANCESTRAL'].str.replace(
        'unspecified', r'(NR)')
    Cat_Anc['BROAD ANCESTRAL'] = Cat_Anc['BROAD ANCESTRAL'].str.replace(
        'Hispanic or Latin American', 'Hispanic\Latin American')
    Cat_Anc['BROAD ANCESTRAL'] = Cat_Anc['BROAD ANCESTRAL'].str.replace(
        'African American or Afro-Caribbean', 'African Am.\Caribbean')
    Cat_Full = pd.read_csv(os.path.abspath(os.path.join('__file__', '../..',
                                                        'Data', 'Catalogue',
                                                        'Raw',
                                                        'Cat_Full.tsv')),
                           header=0, sep='\t', encoding='utf-8',
                           index_col=False, low_memory=False)
    return Cat_Stud, Cat_Anc, Cat_Anc_byN, Cat_Full


def load_pubmed_data():
    """" Load the PUBMED data into the workspace"""
    FunderInfo = pd.read_csv(os.path.abspath(
                             os.path.join('__file__', '../..', 'Data',
                                          'PUBMED', 'Pubmed_FunderInfo.csv')))
    AbstractInfo = pd.read_csv(os.path.abspath(
                               os.path.join('__file__', '../..', 'Data',
                                            'PUBMED',
                                            'Pubmed_AbstractInfo.csv')))
    AbstractInfo['Abstracts'] = AbstractInfo['Abstract'].apply(
        lambda x: re.sub('[^a-zA-Z ]+', '', x)).str.lower()

    AbstractCount = pd.DataFrame(AbstractInfo.Abstract.apply(lambda x:
                                                             pd.value_counts(x.split(' '))).sum(axis=0),
                                                             columns = ['num_words'])
    AbstractCount.index.name = 'word'
    AbstractCount.to_csv(os.path.abspath(os.path.join(
        '__file__', '../..', 'Data', 'PUBMED', 'Pubmed_AbstractCount.csv')))

    AuthorMaster = pd.read_csv(os.path.abspath(
                               os.path.join('__file__', '../..',
                                            'Data', 'PUBMED',
                                            'Pubmed_AuthorInfo.csv')))
    AuthorMaster = AuthorMaster.drop_duplicates(
        subset=['PUBMEDID', 'AUTHORNAME'])
    FunderInfo['Agency'] = FunderInfo['Agency'].str.replace(
        'Medical Research Council', 'MRC')
    return FunderInfo, AbstractInfo, AuthorMaster


def make_timely(variables, yearlist, yearquarterlist, Cat_Stud,
                Cat_Anc, Cat_Anc_byN):
    """ build the dataframe which plots the longitudinal growth of GWAS over
    time
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
    fileout = open(os.path.abspath(
                   os.path.join('__file__', '../..', 'Data',
                                'Catalogue', 'Synthetic',
                                'ancestry_CoR.csv')),
                   'w')
    fileout.write('Date,PUBMEDID,N,Cleaned Country\n')
    for index, row in Cat_Anc.iterrows():
        if len(row['COUNTRY OF RECRUITMENT'].split(',')) == 1:  # or contains ?
            fileout.write(row['DATE'] + ',' + str(row['PUBMEDID']) + ',' +
                          str(row['N']) + ',' + row['COUNTRY OF RECRUITMENT'] +
                          '\n')
    fileout.close()
    Clean_CoR = pd.read_csv(os.path.abspath(
                                  os.path.join('__file__', '../..',
                                               'Data', 'Catalogue',
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
    print('Cleaning for single country of recruitment field results in ' +
          str(round((len(Clean_CoR) / len(Cat_Anc)) * 100, 3)) +
          '% of ancestry observations remaining. \nThis represents about ' +
          str(round((Clean_CoR['N'].sum() / Cat_Anc['N'].sum()) * 100, 3)) +
          '% of the total GWAS N. \nWhen we drop for country==NR,' +
          ' we lose another: ' +
          str(round(len(Clean_CoR[Clean_CoR['Cleaned Country'] == 'NR']) /
                       (len(Clean_CoR)) * 100, 2)) + '% papers.')
    Clean_CoR = Clean_CoR[Clean_CoR['Cleaned Country'] != 'NR']
    print('In total, we have ' + str(len(Clean_CoR)) +
          ' observations for this field, out of a total of ' +
          str(len(Cat_Anc)) + ' rows of Cat_Anc data')
    Clean_CoR.to_csv(os.path.abspath(
                     os.path.join('__file__', '../..', 'Data', 'Catalogue',
                                  'Synthetic',
                                  'GWASCatalogue_CleanedCountry.tsv')),
                     sep='\t', index=False)
    return Clean_CoR
