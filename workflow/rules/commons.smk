
# define list of species to annotate
yaml_species = config['species']
print(yaml_species)

rule create_infos_genomes:
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
