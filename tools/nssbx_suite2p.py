"""
    nssbx_suite2p 
    A few functions to save dicts to json from matlab, then read those again 
    from python and then call suite2p. 

    This is done to avoid calling inprocess Python from Matlab which 
    on some installations leads to library conflicts.
    
"""
import json
from  suite2p import run_s2p
from  sys import argv

def save_dict_to_file(dictionary, filename):
    with open(filename, 'w') as file:
        json.dump(dictionary, file)

def load_dict_from_file(filename):
    with open(filename, 'r') as file:
        dictionary = json.load(file)
    return dictionary

if __name__ == '__main__':
    print('Hello world, reading ops and db, then starting suite2p')
    ops=load_dict_from_file(argv[1])
    db =load_dict_from_file(argv[2])
    print(ops)
    print(db)
    run_s2p(ops, db);

