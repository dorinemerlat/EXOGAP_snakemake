
import pandas as pd
import argparse
import yaml 
import os 
import requests
import sys
import json
from Bio import Entrez
Entrez.email = 'A.N.Other@example.com'
import time

def get_lbgi_infos(parent):
    return parent['name'] + '/' + str(parent['id'])


def get_ncbi_infos(parent):
    return parent['ScientificName'] + '/' + str(parent['TaxId'])


def extract_name(row):
    if row['specie'] != '':
        name = row['specie'].split('/')[0].split(' ')
    else:
        name = row['genome'].split('-')
    
    return name


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genomes", help="specie to add in the info_genomes.tsv", type=str, required=True, metavar="FILE")
    parser.add_argument("-o", "--output", help='Phylogeny file', type=str, default = os.getcwd(), metavar="FILE")
    args = parser.parse_args()

    # convert the string with genomes given to argument into a dictionnary
    genomes = json.loads(args.genomes.replace("'", '"'))

    # create dataframe with good column names
    ranks = ['superkingdom', 'kingdom', 'phylum', 'subphylum', 'class', 'subclass', 'order', 'family', 'genus', 'specie']
    column_names = ['genome', 'taxid'] + ranks + ['mnemonic']
    df = pd.DataFrame(columns = column_names)

    # fill dataframe
    index_genomes = 0 
    for genome, taxid in genomes.items():
        taxid = str(taxid)
        new_line = dict.fromkeys(column_names, '')
        new_line['genome'] = genome

        #if user gives a taxid
        if taxid != '0':
            new_line['taxid'] = taxid

            # try with LBGI API
            r = requests.get("https://lbgi.fr/api/taxonomy/lineage/" + str(taxid), headers={ "Accept" : "application/json"})
                    
            # if they are a result in taxonomy database -> extract informations
            if r.ok:
                new_line = dict.fromkeys(column_names, '')
                new_line['genome'] = genome
                new_line['taxid'] = taxid

                for parent in r.json()['data']:
                    if 'rank' in parent:
                        for rank in ranks:
                            if rank == parent['rank']:
                                new_line[rank] = get_lbgi_infos(parent)
                            elif parent['rank'] == 'species':
                                new_line['specie'] = get_lbgi_infos(parent)
                                break

                error = 'no'

            # if they are no results, try with NCBI efetch
            else:
                try:
                    handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
                    error = 0
                except:
                    try:
                        handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
                        error = 0
                    except:
                        pass
                        error = 1
            
            # if NCBI util gives a result
            if error == 0:
                records = Entrez.parse(handle)

                for record in records:
                    new_line['specie'] = record['ScientificName'] + '/' + record['TaxId']

                    for parent in record['LineageEx']:
                        for rank in ranks:
                            if rank == parent['Rank']:
                                new_line[rank] = get_ncbi_infos(parent)
                                break

        # add the mnemonic name
        if taxid != 0:
            # get the Mnemonic name from uniprot if it exists
            url = "https://rest.uniprot.org/taxonomy/stream?fields=id%2Cmnemonic%2Crank&format=json&query=%28{taxid}%29".format(taxid = taxid)
            uniprot = requests.get(url).text
            uniprot = json.loads(uniprot)['results']
            if len(uniprot) > 0:
                mnemonic = uniprot[0]['mnemonic']
                # check if the mnemonic belongs to a species
                url = "https://rest.uniprot.org/taxonomy/stream?fields=id%2Cmnemonic%2Crank&format=json&query=%28{taxid}%29".format(taxid = mnemonic)
                uniprot = requests.get(url).text
                uniprot = json.loads(uniprot)['results']
                if len(uniprot) == 1:
                    new_line['mnemonic'] = mnemonic


        # add the genome to the dataframe 
        new_line = pd.DataFrame(new_line, index = [index_genomes])
        index_genomes += 1
        df = pd.concat([df, new_line])

    df.sort_values(by = ['mnemonic', 'genome'], ascending = [False, True], inplace = True)
    
    # add a mnemonic name for genomes which don't have one 
    for index, line in df.iterrows():
        if len(line['mnemonic']) == 0:
            name = extract_name(line)
            mnemonic = name[0][0:3].upper() + name[1][0:2].upper()
            # check if another genome has the same mnemonic name
            for previous_index, previous_line in df.iterrows():
                previous_mnemonic = previous_line['mnemonic']
                if mnemonic == previous_mnemonic:
                    end_name = name[1]
                    previous_name = extract_name(previous_line)
                    previous_end_name = previous_name[1]
                    i = 0
                    letter = 0
                    previous_letter = 0
                    while i < len(end_name) and i < len(previous_end_name) and letter == 0:
                        if end_name[i] != previous_end_name[i]:
                            letter = end_name[i]
                            previous_letter = previous_end_name[i]
                            mnemonic = mnemonic[0:4] + letter.upper()
                            previous_mnemonic = previous_mnemonic[0:4] + previous_letter.upper()
                            df['mnemonic'][previous_index] = previous_mnemonic
                            df['mnemonic'][index] = mnemonic
                            break
                        i += 1
                else:
                    df['mnemonic'][index] = mnemonic
                
    # sort by taxonomy
    df.sort_values(by = ranks, ascending = [True] * len(ranks), inplace = True)
    
    # Save results
    df.to_csv(args.output, index = False, sep = '\t')

if __name__ == "__main__":
    main()
