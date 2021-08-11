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
    print('No citable material yet, all your base are belong to us')
    
    print('\n\nIn addition, the following third-party material should be referenced:')
    for db in used_dbs:
        print('\n')
        if db.upper() not in ['KEGG', 'REACTOME']:
            print(cits[db.lower()])
        else:
            print(cits['keggreactome'])
