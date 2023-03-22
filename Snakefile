"""
https://github.com/dorinemerlat/EXOGAP

Snakefile to run EXOGAP
-----------------------------------------

Requirements:
 - SnakeMake (http://snakemake.readthedocs.io/en/stable/)

Basic usage:
  snakemake -p \
    --configfile config.yaml \
    -j 64 

Author:
  Dorine MERLAT

Contact:
  dorine.merlat@etu.unistra.fr  

License:
  GNU General Public License
"""

import os
import sys


# general paths
exogap_dir = os.getcwd()
info_genomes = "config/info_genomes.tsv"

# define which module/pipeline has been launched
if sys.argv[-1] in ['annotate_repetitive_elements']:
  launched = sys.argv[-1]
else:
  launched = "pipeline"


# environment paths
envs_dir = exogap_dir + "/workflow/envs/"
env_annotation = "%sannotation.yaml" % envs_dir
env_repeatmodeler = "%srepeatmodeler.yaml" % envs_dir
env_toolbox = "%stoolbox.yaml" % envs_dir


# scripts
scripts_dir = exogap_dir + "/workflow/scripts/"


# rules
rules_dir = exogap_dir + "/workflow/rules/"
include: "%sinstall_environments.smk" % rules_dir
include: "%scommons.smk" % rules_dir
include: "%sannotate_repetitive_elements.smk" % rules_dir



# target rule
rule all:
    """
    Rule to define all output
    """
    output: "file.1"

def read_infos_genomes(specie, rank = "Specie"):
    """
    get information in the config/info_genomes.tsv in a dictionnary
    """
    if rank is not False:
        with open('config/info_genomes.tsv', 'r') as file:
            lines = file.read().splitlines()
            # get header as description of fields
            description = lines[0].split('\t')

            # read file if it is not empty
            if len(lines) > 1:
                INFOS_SPECIES = dict()
                for line in lines[1:]: 
                    line = line.split('\t') 

                    # create a dictionnary for each entry in the file
                    specie = dict() 
                    for i in list(range(1, len(line))):

                        # if only one information in the feature
                        if '/' not in line[i]:
                            specie[description[i]] = line[i]
                        
                        # for feature with taxon name and taxon taxid
                        else:
                            elt = line[i].split('/')
                            specie[description[i]] = {'genome': elt[0], 'taxid': elt[1]}
                    
                    # save dictionnary in the global doctionnary
                    INFOS_SPECIES[line[0]] = specie
 
            else: 
                print("error: phylogeny.csv is empty")

            return INFOS_SPECIES