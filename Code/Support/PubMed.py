import csv
import os
from Bio import Entrez
import time


def build_author(papers):
    ''' Take the read return from the Entrez.efetch utility and
    parse out the list of authors which is then saved to the PUBMED
    subdirectory in the Data folder.
    '''
    noforename = 0
    yesforename = 0
    with open(os.path.abspath(os.path.join('__file__', '../..',
                                           'Data', 'PUBMED',
                                           'Pubmed_AuthorInfo.csv')),
              'w', newline='', encoding='utf-8') as csv_file:
        writer_author = csv.writer(csv_file, delimiter=',')
        writer_author.writerow(
            ['PUBMEDID', 'FORENAME', 'LASTNAME', 'AUTHORNAME'])
        for i, paper in enumerate(papers['PubmedArticle']):
            AuthorList = paper['MedlineCitation']['Article']['AuthorList']
            for author in range(0, len(AuthorList)):
                PUBMEDID = paper['MedlineCitation']['PMID']
                if 'LastName' in AuthorList[author]:
                    LASTNAME = AuthorList[author]['LastName']
                    if 'ForeName' in AuthorList[author]:
                        FORENAME = AuthorList[author]['ForeName']
                        yesforename += 1
                    else:
                        FORENAME = 'N/A'
                        noforename += 1
                    if ('ForeName' in AuthorList[author]) and \
                            ('LastName' in AuthorList[author]):
                        AUTHORNAME = FORENAME + ' ' + LASTNAME
                        writer_author.writerow(
                            [PUBMEDID, FORENAME, LASTNAME, AUTHORNAME])
        print('Authors with last names but no forenames: ' +
              str(noforename) + ' out of ' + str(yesforename))
    print('Built a database of Authors from list of PUBMEDID IDs')


def build_collective(papers):
    ''' Take the read return from the Entrez.efetch utility and
    parse out the list of collectives which is then saved to the
    PUBMED subdirectory in the Data folder.
    '''
    numberofcollectives = 0
    with open(os.path.abspath(os.path.join('__file__', '../..',
                                           'Data', 'PUBMED',
                                           'Pubmed_CollectiveInfo.csv')),
              'w', newline='', encoding='utf-8') as csv_file:
        writer_author = csv.writer(csv_file, delimiter=',')
        writer_author.writerow(['PUBMEDID', 'COLLECTIVE'])
        for i, paper in enumerate(papers['PubmedArticle']):
            AuthorList = paper['MedlineCitation']['Article']['AuthorList']
            for author in range(0, len(AuthorList)):
                PUBMEDID = paper['MedlineCitation']['PMID']
                if 'CollectiveName' in AuthorList[author]:
                    COLLECTIVE = AuthorList[author]['CollectiveName']
                    numberofcollectives += 1
                    writer_author.writerow([PUBMEDID, COLLECTIVE])
        print('Number of Collectives Found: ' + str(numberofcollectives))
    print('Built a database of Authors from list of PUBMEDID IDs')


def build_funder(papers):
    ''' Take the read return from the Entrez.efetch utility and
    parse out the list of collectives which is then saved to the
    PUBMED subdirectory in the Data folder.
    '''
    with open(os.path.abspath(os.path.join('__file__', '../..',
                                           'Data', 'PUBMED',
                                           'Pubmed_FunderInfo.csv')),
              'w', newline='', encoding='utf-8') as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow(['PUBMEDID', 'Agency', 'GrantCountry', 'GrantID'])
        for i, paper in enumerate(papers['PubmedArticle']):
            try:
                GrantList = paper['MedlineCitation']['Article']['GrantList']
                for grant in range(0, len(GrantList)):
                    try:
                        PUBMEDID = paper['MedlineCitation']['PMID']
                    except:
                        PUBMEDID = 'N/A'
                    try:
                        Agency = GrantList[grant]['Agency']
                    except:
                        Agency = 'N/A'
                    try:
                        GrantCountry = GrantList[grant]['Country']
                    except:
                        GrantCountry = 'N/A'
                    try:
                        GrantID = GrantList[grant]['GrantID']
                    except:
                        GrantID = 'N/A'
                    writer.writerow([PUBMEDID, Agency, GrantCountry, GrantID])
            except:
                pass
    print('Built a database of Funders from list of PUBMEDID IDs')


def build_abstract(papers):
    '''Take the read return from the Entrez.efetch utility and
    parse out the list of abstracts which is then saved to the
    PUBMED subdirectory in the Data folder.
    '''
    with open(os.path.abspath(os.path.join('__file__', '../..',
                                           'Data', 'PUBMED',
                                           'Pubmed_AbstractInfo.csv')),
              'w', newline='', encoding='utf-8') as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow(['PUBMEDID', 'Abstract'])
        for i, paper in enumerate(papers['PubmedArticle']):
            try:
                AbstractList = paper['MedlineCitation']['Article']['Abstract']
                PUBMEDID = paper['MedlineCitation']['PMID']
                Abstract = AbstractList['AbstractText']
            except KeyError:
                Abstract = 'N/A'
            writer.writerow([PUBMEDID, ', '.join(Abstract)])
    print('Built a database of Abstracts from list of PUBMEDID IDs')


def build_citation(id_list, emailaddress):
    ''' This function queries the PubMed Central database to return
    information on citation counts which is then saved to the
    PUBMED subdirectory in the Data folder.
    '''
    citationcountout = open(os.path.abspath(os.path.join('__file__',
                                                         '../..',
                                                         'Data',
                                                         'PUBMED',
                                                         'Pubmed_Cites.csv')),
                            'w', encoding='utf-8')
    citationcountout.write('PUBMEDID' + ',' + 'citedByCount\n')
    Entrez.email = emailaddress
    for pubmedid in id_list:
        results = Entrez.read(Entrez.elink(dbfrom='pubmed', db='pmc',
                                           LinkName='pubmed_pmc_refs',
                                           id=pubmedid))
        try:
            citationcountout.write(
                str(pubmedid) + ',' +
                str(len(results[0]['LinkSetDb'][0]['Link'])) + '\n')
        except IndexError:
            citationcountout.write(str(pubmedid) + ',' + '0\n')
        time.sleep(1.5)
    citationcountout.close()
