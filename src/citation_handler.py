__author__ = 'Rasmus Magnusson'
__COPYRIGHT__ = 'Copyright (C) 2020 Rasmus Magnusson'
__contact__ = 'rasma774@gmail.com'

def load_citations():
    path =  __file__.split('/src')[0] 
    
    with open(path + '/data/refs.txt') as f:
        cits = {}
        for line in f:
            line = line.split('%')
            cits[line[0]] = line[1][:-2]
    return cits

def write_citations(used_dbs):
    cits = load_citations()
    
    print('Thank you for using TFTenricher. To support our work, please cite:\n')
    print('Magnusson, R., Lubovac-Pilav, Z. TFTenricher: a python toolbox for annotation enrichment analysis of transcription factor target genes. BMC Bioinformatics 22, 440 (2021). https://doi.org/10.1186/s12859-021-04357-4')
    
    print('\n\nIn addition, the following third-party material should be referenced:')
    for db in used_dbs:
        print('\n')
        if db.upper() not in ['KEGG', 'REACTOME']:
            print(cits[db.lower()])
        else:
            print(cits['keggreactome'])
