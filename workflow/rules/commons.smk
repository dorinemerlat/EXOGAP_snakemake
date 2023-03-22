# define list of species to annotate
yaml_species = config['species']

rule create_info_genomes:
    """
    create a tsv file with all informations about genomes (phylogeny, taxid, 5 letters code)
    """
    input:
        env = rules.install_toolbox.output.env
    output:
        info_genomes = 'config/info_genomes.tsv'
    params: 
        yaml_species = yaml_species,
        script = "%screate_info_genomes.py" % scripts_dir
    conda:
        env_toolbox
    shell:
        """
        python {params.script} -g "{params.yaml_species}" -o {output.info_genomes}
        """
    

# read file if it is not empty
if os.path.isfile("config/info_genomes.tsv"):
    with open('config/info_genomes.tsv', 'r') as file:
        lines = file.read().splitlines()
        # get header as description of fields
        description = lines[0].split('\t')

        if len(lines) > 1:
            GENOMES = dict()
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
                        specie[description[i]] = {'name': elt[0], 'taxid': elt[1]}
                
                # save dictionnary in the global doctionnary
                GENOMES[line[0]] = specie


def get_species_in_group(GENOMES, clade, level):
    group = list()

    for genome, infos in GENOMES.items():
        if infos[level]['name'] == clade:
            group.append(genome)

    return group


rule reformat_genome:
    input:
        genome = "resources/genomes/{specie}.fa"
    output:
        reformated_genome = "results/{specie}/reformated_genome/{specie}.fa"
    params:
        script = "workflow/scripts/reformatGenome.py",
        outfiles = "results/{specie}/reformated_genome/{specie}",
        code = lambda wildcards : get_taxid_in_phylogeny('{name}'.format(name=wildcards.specie), "Code")
    threads:
        1
    conda:
        "../envs/toolbox.yaml"
    shell:
        "{params.script} -i {input.genome} -o {params.outfiles} -p {params.code}seq"
