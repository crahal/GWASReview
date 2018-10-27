import re


def ancestry_cleaner(row, field):
    """ clean up the ancestry fields in GWASCatalogue_Ancestry.

    Keyword arguments:
    row: the row of the ancestry DataFrame
    field: the field of the ancestry dataframe ('initial' or 'replication')
    """
    free_text = re.sub(r'(\d+),?([\d+]?)', r'\1\2', row[field])
    free_text = re.sub(r'(\d+)', r'; \1', free_text)
    free_text = punctuation_cleaner(free_text)
    free_text = remove_lower(free_text)
    free_text = remove_lower(free_text)
    free_text = free_text.replace('  ', ' ')
    free_text = free_text.replace('  ', ' ')
    free_text = list_remover(free_text)
    free_text = dict_replace(free_text)
    try:
        if free_text[-1] == ';':
            free_text = free_text[:-1]
    except ValueError:
        pass
    cleaned = []
    for ancestry in free_text[1:].split(';'):
        if " and" in ancestry.strip()[-4:]:
            cleaned.append(ancestry.replace(' and', '').strip())
        elif " or" in ancestry.strip()[-4:]:
            cleaned.append(ancestry.replace(' or', '').strip())
        else:
            cleaned.append(ancestry.strip())
    cleaned = ';'.join(cleaned)
    cleaned = cleaned.replace(';', ' ; ')
    for word in cleaned.split(' '):
        if (word.isalpha()) and (len(word) < 3) and word != "or":
            cleaned = cleaned.replace(word, '')
    cleaned = re.sub(r';\s+;', ';', cleaned)
    return cleaned


def dict_replace(text):
    """ sanitize the free text strings from the initial/replication fields.

    Keyword arguements:
    text: the free text string prior to splitting
    """
    replacedict = {'Arabian': 'Arab',
                   'HIspanic': 'Hispanic',
                   'Korculan': 'Korcula',
                   'Hispaic': 'Hispanic',
                   'Hispanics': 'Hispanic',
                   'Chineses': 'Chinese',
                   'Europea ': 'European ',
                   'Finish': 'Finnish',
                   'Val Bbera': 'Val Borbera',
                   'Erasmus Rchen': 'Erasmus Rucphen',
                   'Erasmus Rupchen': 'Erasmus Rucphen',
                   'Erasmus Ruchpen': 'Erasmus Rucphen',
                   'Chinese Han': 'Han Chinese',
                   'Cilen ': 'Cilento',
                   'Clien': 'Cilento',
                   'Old Amish': 'Old Order Amish',
                   'Americans': 'American',
                   'Geman': 'German',
                   'Japnese': 'Japanese',
                   'Finland': 'Finnish',
                   'Eat Aian': 'East Asian',
                   'Hipanic': 'Hispanic',
                   'Sub African': 'Sub-saharan African',
                   'Erasmus Rucphen Family': 'Erasmus Rucphen',
                   'Nfolk Island': 'Norfolk Island',
                   'Israeli': 'Isreali',
                   'Sh Asian': 'Asian',
                   'Hispanic Latino': 'Hispanic/Latino',
                   'Hispanic Latino': 'Hispanic/Latino',
                   'Hispanic Latin ': 'Hispanic/Latino',
                   'European ad': 'European',
                   'European End': 'European',
                   'LatinoAmerican': 'Latino American',
                   'Val Bbera': 'Val Borbera',
                   'Oceanian': 'Oceania',
                   'Cilentoto': 'Cilento',
                   'Uighur': 'Uyghur',
                   'Friuli': 'Fruili',
                   'Giuli': 'Giulia',
                   'Giuliaa': 'Giulia',
                   'Rupchen': 'Rucphen',
                   'Korcula': 'Korkulan',
                   'Ruchpen': 'Rucphen',
                   'Brazillian': 'Brazilian',
                   'Sub-saharan': 'Sub Saharan',
                   'Tyrolian': 'Tyrolean',
                   'Seychelles': 'Seychellois',
                   'South Tyrolean': 'South Tyrol',
                   'Europen': 'European',
                   '≥': ''}
    for key in replacedict:
        if key in text:
            text = text.replace(key, replacedict[key])
    return text


def list_remover(text):
    """ titlecase/capitalised words to remove from the strings which
    are not associated with countries, races or ancestries.

    Keyword arguements:
    text: the free text string prior to splitting
    """
    removelist = ['AIS', 'APOE', 'HIV', 'Â', 'HER2-', '1000 Genomes',
                  'MYCN-amplification', 'Alzheimer', 'ASD', 'OCB', 'BD',
                  'Genetically', 'Homogenous', 'BRCA', 'ALL', 'Coronary',
                  'Amyotrophic', 'Large', 'anti-dsDNA', 'Up ', 'Biracial',
                  'Follicular', 'Hodgkin', 'Lymphoma', 'GI', 'Abstinent',
                  'Schizophrenia', 'Îµ', 'JAK', 'ADHD', 'Diabetes',
                  'Allogenic', 'BGPI', 'Ischemic', 'Chronic', 'Major',
                  'Diabetic', 'Microalbuminuria', 'Asthma', 'Individuals',
                  'At ', "Barrett's", 'Crohn', 'Bipolar', 'MMR', 'HBV', 'RA',
                  'Elated', 'Escitalpram', 'Irritable', 'Lymphoblastoid',
                  'ACPA', 'HCC', 'pPhe508del', 'Anti', 'B2GPI', 'Kashin Beck',
                  '(LDL-cholesterol)', 'TPO', 'OCD', 'CCT', 'FTD', 'CAPOX B',
                  'LAC', 'LOAD', ' So ', 'MYCN-amplification', 'Yang', 'Tae',
                  'Eum', 'Non-abstinent', 'EBWL', 'Semantic', 'General',
                  'Cluster', 'Frontremporal', 'Frontotremporal',
                  'Frontotemporal', 'Graves', 'Attention', 'Autism', 'Liu',
                  'High', 'Low', 'HCV', 'Citalopram', 'Haemophilia', ' III ',
                  ' II ', ' I ', 'NFT', 'Progressive', 'Ancestry', 'Parkinson',
                  'Lin', 'BMD', 'GBA', 'Traylor', 'Consortium', ' Torgerson',
                  'EVE', 'Germain', 'Boraska', 'Cases', 'HapMap', 'vWF', 'HDL',
                  'LDL', ' Mild', 'Cognitive', 'Impairment', 'Sarcoidosis',
                  'Yu Zhi', 'Lymphoma', 'Impairment', 'Type', 'Kuru',
                  'Frontemporal', 'Erasmus', 'Barrett', 'Lofgren', 'Hashimoto',
                  'Family', 'Multiple', 'Richardson', 'Metropolitan']
    for word in removelist:
        text = text.replace(word, '')
    return text


def remove_lower(free_text):
    """ remove lowercase letters (assumed to not be associated with
    countries, races or ancestries.)

    Keyword arguements:
    text: the free text string prior to splitting
    """
    free_text = free_text.replace('up to', '')
    for word in free_text.split(' '):
        if (word.title() != word.strip()):
            try:
                float(word)
            except ValueError:
                if ';' in word:
                    free_text = free_text.replace(word, ';').strip()
                elif (';' not in word) and (word != "and") and (word != "or"):
                    if free_text.find(word) == 0:
                        free_text = free_text.replace(word + ' ', ' ')
                    else:
                        free_text = free_text.replace(' ' + word, ' ')
    return free_text.strip()


def punctuation_cleaner(temp):
    """ remove various punctuation (assumed to not be associated with
    countries, races or ancestries.)

    Keyword arguements:
    text: the free text string prior to splitting
    """
    temp = temp.replace(',', ';')
    for pmark in ['-', '\'', '’', '?', '+']:
        temp = temp.replace(pmark, ' ')
    for pmark in ['(', ')', '.', '*', '~', '<', '>']:
        temp = temp.replace(pmark, '')
    return temp
