import re
import pandas as pd


def simple_facts(Cat_Studies, Cat_Ancestry,
                 Cat_Ancestry_groupedbyN,
                 Cat_Full, AuthorMaster):
    """ Print out some simple facts from the GWAS Catalog """

    print('There are currently: ' +
          str(len(Cat_Studies['PUBMEDID'].unique())) +
          ' GWAS papers published (in terms of unique PUBMEDIDs)')
    print('The first observed GWAS in the catalogue was dated: ' +
          str(Cat_Studies['DATE'].min()))
    print('However, only: ' +
          str(len(Cat_Studies[Cat_Studies['DATE'].str.contains('2005|2006')])) +
          ' papers were published in 2005 and 2006 combined')
    print('There are currently: ' +
          str(len(Cat_Studies['STUDY ACCESSION'].unique())) +
          ' unique Study Accessions')
    print('There are currently: ' +
          str(len(Cat_Studies['DISEASE/TRAIT'].unique())) +
          ' unique Diseases/Traits studied')
    print('These correspond to: ' +
          str(len(Cat_Studies['MAPPED_TRAIT'].unique())) +
          ' unique EBI "Mapped Traits"')
    print('The total number of Associations found is currently: ' +
          str(round(Cat_Studies['ASSOCIATION COUNT'].sum(), 1)))
    print('The average number of Associations is currently: ' +
          str(round(Cat_Studies['ASSOCIATION COUNT'].mean(), 1)))
    print('Mean P-Value for the strongest SNP risk allele is currently: ' +
          str(round(Cat_Full['P-VALUE'].astype(float).mean(), 10)))
    print('The number of associations reaching the 5e-8 threshold is ' +
          str(len(Cat_Full[Cat_Full['P-VALUE'].astype(float) < 5.000000e-8])))
    print('The journal to feature the most GWAS studies since ' +
          str(Cat_Studies['DATE'].min()) +
          ' is: ' + Cat_Studies['JOURNAL'].mode()[0])
    print('However, in 2017, ' + Cat_Studies[Cat_Studies['DATE'].str.contains(
          '2017')]['JOURNAL'].mode()[0] +
          ' published the largest number of GWAS papers')
    print('Largest Accession to date: ' +
          str(Cat_Ancestry_groupedbyN.sort_values(by='N',
                                                  ascending=False)['N'].iloc[0]) +
          ' people.\nThis was published in ' +
          str(Cat_Studies[
              Cat_Studies['STUDY ACCESSION'] ==
              Cat_Ancestry_groupedbyN.sort_values(by='N',
                                                  ascending=False)
              ['STUDY ACCESSION'].iloc[0]]['JOURNAL'].iloc[0]) +
          '.\nThe first author was: ' +
          str(Cat_Studies[
              Cat_Studies['STUDY ACCESSION'] ==
              Cat_Ancestry_groupedbyN.sort_values(by='N',
                                                  ascending=False)
              ['STUDY ACCESSION'].iloc[0]]
              ['FIRST AUTHOR'].iloc[0]) + '.')
    print('Total number of SNP-trait associations is ' +
          str(Cat_Studies['ASSOCIATION COUNT'].sum()) + '.')
    print('Total number of journals publishing GWAS is ' +
          str(len(Cat_Studies['JOURNAL'].unique())) + '.')
    print('The study with the largest number of authors has: ' +
          str(AuthorMaster.groupby(['PUBMEDID'])['AUTHORNAME'].count().max()) +
          ' authors.')


def ancestry_parser(output_path, input_series, Cat_Studies):
    fileout = open(output_path, 'w', encoding='utf-8')
    fileout.write('STUDY ACCESSION,Cleaned Ancestry,Cleaned Ancestry Size\n')
    for index, row in Cat_Studies.iterrows():
        checksum = 0
        for ancestry in row[input_series].split(';'):
            number = re.findall(r'\d+', ancestry.strip())
            if (len(number) == 1):
                checksum += 1
        if checksum == len(row[input_series].split(';')):
            for ancestry in row[input_series].split(';'):
                number = re.findall(r'\d+', ancestry.strip())
                words = ''.join(i for i in ancestry.strip() if not i.isdigit())
                if (len(number) == 1) and (len(words.strip()) > 3) and \
                   (sum(1 for c in words if c.isupper()) > 0):
                    fileout.write(row['STUDY ACCESSION'] + ',' +
                                  words.strip() + ',' + str(number[0]) + '\n')
    fileout.close()


def make_meanfemale_andranks(AuthorMaster, EFO):
    """ make inputs for the box and swarm figure """
    AuthorMaster_merged = pd.merge(AuthorMaster[(
        AuthorMaster['MaleFemale'] == 'male') |
        (AuthorMaster['MaleFemale'] == 'female')],
        EFO, how='left', on='PUBMEDID')
    AuthorMaster_EFO = AuthorMaster_merged.groupby(
        ['EFO term'])['isfemale'].mean().to_frame()
    AuthorMaster_Parent = AuthorMaster_merged.groupby(
        ['Parent term'])['isfemale'].mean().to_frame()
    countstudiesperEFO = EFO.groupby(['EFO term'])['PUBMEDID'].count(
    ).sort_values(ascending=False)[0:14].index.tolist()
    meanfemalestudyaccession = AuthorMaster_merged.groupby(
        ['STUDY ACCESSION'])['isfemale'].mean().to_frame().reset_index()
    meanfemalestudyaccession_withparent = pd.merge(EFO[[
        'STUDY ACCESSION', 'Parent term']], meanfemalestudyaccession,
        on='STUDY ACCESSION', how='left')
    meanfemalestudyaccession_withparent = meanfemalestudyaccession_withparent[
        ~meanfemalestudyaccession_withparent['isfemale'].isnull()]
    meanfemalestudyaccession_withEFO = pd.merge(EFO[[
        'STUDY ACCESSION', 'EFO term']], meanfemalestudyaccession,
        on='STUDY ACCESSION', how='left')
    meanfemalestudyaccession_withEFO = meanfemalestudyaccession_withEFO[
        meanfemalestudyaccession_withEFO['EFO term'].isin(
            countstudiesperEFO)]
    meanfemalestudyaccession_withEFO = meanfemalestudyaccession_withEFO[
        ~meanfemalestudyaccession_withEFO['isfemale'].isnull()]
    AuthorMaster_EFO.reset_index(inplace=True)
    ranks = AuthorMaster_EFO[AuthorMaster_EFO['EFO term'].isin(
        countstudiesperEFO)].sort_values(by='isfemale',
                                         ascending=False)['EFO term'].tolist()
    return ranks, countstudiesperEFO, meanfemalestudyaccession_withEFO, AuthorMaster_EFO, AuthorMaster_Parent


def make_funders(FunderInfo_Parent, FunderInfo, Cat_Ancestry):
    """ make the dataframes for the heatmaps """
    funder_parent = pd.DataFrame(
        index=FunderInfo.groupby(
            ['Agency'])['Agency'].count().sort_values(ascending=False).
        index.values[0:10].tolist(),
        columns=FunderInfo_Parent.groupby(['Parent term'])['Parent term'].count().
        sort_values(ascending=False).index.values.tolist())

    for parent in FunderInfo_Parent.groupby(['Parent term'])['Parent term'].count().\
            sort_values(ascending=False).index.values.tolist():
        for funder in FunderInfo.groupby(['Agency'])['Agency'].count().sort_values(
                ascending=False).index.values[0:10].tolist():
            funder_parent[parent][funder] = len(
                FunderInfo_Parent[(FunderInfo_Parent['Agency'] == funder) & (
                    FunderInfo_Parent['Parent term'] == parent)])
    FunderInfo_Ancestry = pd.merge(FunderInfo, Cat_Ancestry, left_on='PUBMEDID',
                                   right_on='PUBMEDID', how='left')

    funder_ancestry = pd.DataFrame(index=FunderInfo.groupby(
        ['Agency'])['Agency'].count().sort_values(ascending=False).
        index.values[0:10].tolist(),
        columns=FunderInfo_Ancestry.groupby(
        ['Broader'])['Broader'].count().
        sort_values(ascending=False).index.values.tolist())

    for ancestry in FunderInfo_Ancestry.groupby(
        ['Broader'])['Broader'].count().\
            sort_values(ascending=False).index.values.tolist():
        for funder in FunderInfo.groupby(['Agency'])['Agency'].count().sort_values(
                ascending=False).index.values[0:15].tolist():
            funder_ancestry[ancestry][funder] = len(
                FunderInfo_Ancestry[(FunderInfo_Ancestry['Agency'] == funder) & (
                    FunderInfo_Ancestry['Broader'] == ancestry)])
    return funder_ancestry, funder_parent


def mapped_trait_summary(input_df, input_series):
    """ Make a summary of the mapped traits seen in the Catalog """
    df_Count = input_df.groupby(input_series)[input_series].count(
    ).to_frame().rename(columns={input_series: 'Number of Studies'}).reset_index()
    df_AncSum = input_df.groupby(input_series)['N'].sum(
    ).to_frame().rename(columns={'N': 'Total Sample'}).reset_index()
    df_AssSum = input_df.groupby(
        input_series)['ASSOCIATION COUNT'].sum().to_frame().rename(columns={
            'ASSOCIATION COUNT': 'Total Associations'}).reset_index()
    df_summary = pd.merge(
        df_AssSum, df_AncSum,
        how='left', on=input_series)
    df_summary = pd.merge(
        df_summary, df_Count, how='left', on=input_series)
    df_summary['Studies'] = (
        df_summary['Number of Studies'] /
        df_summary['Number of Studies'].sum()) * 100
    df_summary['Associations'] = (
        df_summary['Total Associations'] /
        df_summary['Total Associations'].sum()) * 100
    df_summary['Sample'] = (
        df_summary['Total Sample'] /
        df_summary['Total Sample'].sum()) * 100
    df_summary = df_summary.set_index(input_series)
    df_summary = df_summary.sort_values(
        by='Studies', ascending=False)[0:17]
    df_summary = df_summary[[
        'Sample', 'Studies', 'Associations']]
    print(df_summary)
