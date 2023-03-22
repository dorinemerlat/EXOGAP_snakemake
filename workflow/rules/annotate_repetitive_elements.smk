def get_repetitive_elements_options(config, launched):
    opts_repeats = {'group_set': False, 'mites': False, 'helitrons': False, 'ltr': False}

    if launched in ['annotate_repetitive_elements', 'pipeline']:
        opts_repeats['group_set'] = config['repetitive_elements']['group_set']
        opts_repeats['mites'] = config['repetitive_elements']['mites']
        opts_repeats['helitrons'] = config['repetitive_elements']['helitrons']
        opts_repeats['ltr'] = config['repetitive_elements']['ltr']

    return opts_repeats
       

# define options
opts_repeats = get_repetitive_elements_options(config, launched)
        
# rules
rule run_repeatmodeler:
    input:
        genome = rules.reformat_genome.output.reformated_genome,
        env = rules.install_repeatmodeler.output.env
    output:
        repeatsLib = "results/data_sets/repetitive_elements/{specie}/repeatModeler/{specie}-families.fa"
    log:
        repeatModeler = 'logs/run_repeatmodeler/{specie}_repeatModeler.log'
    params:
        rmDir = "results/data_sets/repetitive_elements/{specie}/repeatModeler/",
        exogap = exogap_dir
    threads:
        64
    conda:
        env_repeatmodeler
    shell:
        """
        mkdir -p {params.rmDir}
        cd {params.rmDir}

        BuildDatabase -name {wildcards.specie} -engine ncbi {params.exogap}/{input.genome}

        RepeatModeler -pa {threads} -engine ncbi -database {wildcards.specie} -LTRStruct \
            -ninja_dir $CONDA_DEFAULT_ENV/bin 2>&1 {params.exogap}/{log.repeatModeler}
        """


rule reformat_repeatmodeler_library:
    input:
        unreformated = rules.run_repeatmodeler.output.repeatsLib,
        env = rules.install_toolbox.output.env,
    output:
        formated = "results/data_sets/repetitive_elements/{specie}/repeatModeler/{specie}_repetitive_elements_library.fa"
    conda: 
        env_toolbox
    params:
        script = "workflow/scripts/renameRepeats.py"
    shell:
        """
        python {params.script} -i {input.unreformated} -o {output.formated} -p repeats
        """


rule concatanate_repeats_libraries:
    input:
        repeatsLib = lambda wildcards : expand("results/data_sets/repetitive_elements/{specie}/repeatModeler/{specie}_repetitive_elements_library.fa", 
            specie = get_species_in_group(GENOMES, '{name}'.format(name=wildcards.repeatsLib), opts_repeats['group_set'])),
        phylogeny = info_genomes
    output:
        repeatsLibGrouped = "results/data_sets/repetitive_elements/{repeatsLib}/concatenation/{repeatsLib}_repetitive_elements_library.fa"
    shell:
        """
        cat {input.repeatsLib} > {output.repeatsLibGrouped}
        """


def choice_library_to_split(wildcards):
    if opts_repeats['group_set'] == False:
        library =  expand("results/data_sets/repetitive_elements/{repeatsLib}/repeatModeler/{specie}_repetitive_elements_library.fa",
            specie = wildcards.repeatsLib)
    else:
        library = "results/data_sets/repetitive_elements/{repeatsLib}/concatenation/{repeatsLib}_repetitive_elements_library.fa"
    return library



rule split_repetitive_elements:
    input: 
        repeatsLib = choice_library_to_split,
        env = rules.install_toolbox.output.env
    output:
        known = "results/data_sets/repetitive_elements/{repeatsLib}/split/{repeatsLib}_repetitive_elements_library_known.fa",
        unknown = "results/data_sets/repetitive_elements/{repeatsLib}/split/{repeatsLib}_repetitive_elements_library_unknown.fa"
    params:
        script = "workflow/scripts/sortingRepeats.py"
    conda:
        env_toolbox
    shell:
        """
        mkdir -p results/{wildcards.repeatsLib}/repetitive_elements/split/

        python {params.script} -i {input.repeatsLib} -k {output.known} -u {output.unknown}
        """



def choice_library_for_repeatmasker(wildcards):
    # user's choice to group libraries
    if opts_repeats['group_set'] == False:
        repeatsLib = wildcards.specie
    else:
        repeatsLib = GENOMES[wildcards.specie][opts_repeats['group_set']]['name']
        
    # choice library to use in fuction of the iteration: repbase, the known repetitive elements or the unknown repetitive elements
    if wildcards.repeatsIteration == 'repbaseMask':
        return config["repetitive_elements"]["repBase"]
    elif wildcards.repeatsIteration == 'knownMask':
        library = "results/data_sets/repetitive_elements/{repeatsLib}/split/{repeatsLib}_repetitive_elements_library_known.fa".format(
            repeatsLib = repeatsLib)
    elif wildcards.repeatsIteration == 'unknownmask':
        library = "results/data_sets/repetitive_elements/{repeatsLib}/split/{repeatsLib}_repetitive_elements_library_unknown.fa".format(
            repeatsLib = repeatsLib)

    return library


def choice_genome_for_repeatmasker(wildcards):
    # choice the genome to use in function of the iteration
    if wildcards.repeatsIteration == 'repbaseMask':
        genome =  "results/{specie}/reformated_genome/{specie}.fa"
    elif wildcards.repeatsIteration == 'knownMask':
        genome =  "results/{specie}/repetitive_elements/repeatMasker/repbaseMask/{specie}_repbaseMask.fa"
    elif wildcards.repeatsIteration == 'unknownmask':
        genome = "results/{specie}/repetitive_elements/repeatMasker/knownMask/{specie}_knownMask.fa"
        
    return genome


rule run_repeatmasker:
    input:
        phylogeny = info_genomes,
        genome = choice_genome_for_repeatmasker,
        library = choice_library_for_repeatmasker,
        env = rules.install_annotation.output.env
    output:
        maskedGenome = "results/{specie}/repetitive_elements/repeatMasker/{repeatsIteration}/{specie}_{repeatsIteration}.fa",
        outRM = "results/{specie}/repetitive_elements/repeatMasker/{repeatsIteration}/{specie}_{repeatsIteration}.out",
        catgz = "results/{specie}/repetitive_elements/repeatMasker/{repeatsIteration}/{specie}_{repeatsIteration}.cat"
    params:
        name = "{specie}_{repeatsIteration}",
        repeatmaskerDir = "results/{specie}/repetitive_elements/repeatMasker/{repeatsIteration}"
    threads:
        20
    conda:
        env_annotation
    shell:
        """
        RepeatMasker -pa {threads} -e ncbi -nolow -lib {input.library} -a -gccalc -norna -excln \
        -gff -s -xsmall -gccalc -excln -gff -s -dir {params.repeatmaskerDir} {input.genome}

        rename -v -E "s/[[:alnum:]_-]*\.fa/{params.name}/" {params.repeatmaskerDir}/*
        mv {params.repeatmaskerDir}/{params.name}.masked {params.repeatmaskerDir}/{params.name}.fa

        for i in {params.repeatmaskerDir}/*
        do 
            if [[ $i == *"cat.gz" ]]
            then 
                gunzip $i
            fi
        done
        """


repeatsIteration = ["unknownmask", "knownMask", "repbaseMask"]


rule process_repeatmasker:
    input: 
        maskedGenome = "results/{specie}/repetitive_elements/repeatMasker/unknownmask/{specie}_unknownmask.fa",
        outGenome = 'results/{specie}/repetitive_elements/repeatMasker/unknownmask/{specie}_unknownmask.out', 
        cat = expand("results/{specie}/repetitive_elements/repeatMasker/{repeatsIteration}/{specie}_{repeatsIteration}.cat",
            specie = "{specie}",
            repeatsIteration = repeatsIteration),
        library = config["repetitive_elements"]["repBase"],
        env = rules.install_annotation.output.env
    output:
        maskedGenomeFinal = 'results/{specie}/repetitive_elements/repeatMasker/finalMask/{specie}_repetitive_elements.fa',
        outFinal = 'results/{specie}/repetitive_elements/repeatMasker/finalMask/{specie}_repetitive_elements.out',
        catFinal = 'results/{specie}/repetitive_elements/repeatMasker/finalMask/{specie}_repetitive_elements.cat',
        tblFinal = 'results/{specie}/repetitive_elements/repeatMasker/finalMask/{specie}_repetitive_elements.tbl',
        alignFinal = 'results/{specie}/repetitive_elements/repeatMasker/finalMask/{specie}_repetitive_elements.align',
        gffAll = 'results/{specie}/repetitive_elements/repeatMasker/finalMask/{specie}_repetitive_elements_all.gff',
        gffComplex = 'results/{specie}/repetitive_elements/repeatMasker/finalMask/{specie}_repetitive_elements_complex.gff',
        gffFinal = "results/{specie}/repetitive_elements/repeatMasker/finalMask/{specie}_repetitive_elements.gff"
    params:
        finalDir = 'results/{specie}/repetitive_elements/repeatMasker/finalMask',
        name = "{specie}_repetitive_elements"
    conda:
        env_annotation
    shell:
        """
        exogap=$(pwd)

        cp {input.maskedGenome} {output.maskedGenomeFinal}
        cp {input.outGenome} {output.maskedGenomeFinal}

        cat {input.cat} > {output.catFinal}

        # To generate tbl file and align files
        cd {params.finalDir}
        ProcessRepeats -a -gff -lib $exogap/{input.library} $exogap/{output.catFinal} 

        cd $exogap
        mv {params.finalDir}/*out.gff {output.gffAll}

        # Isolation of complex repeats
        grep -v -e "Satellite" -e ")n" -e "-rich" {output.gffAll} > {output.gffComplex}

        # Reformatting to make it work with Maker
        cat {output.gffComplex} | \
            perl -ane '$id; if(!/^\#/){{@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ =join("\t", @F)."\n"}} print $_' \
            > {output.gffFinal}
        """


rule annotate_repetitive_elements:
    input: expand('results/{specie}/repetitive_elements/repeatMasker/finalMask/{specie}_repetitive_elements.fa', specie = GENOMES)