rule install_toolbox:
    output:
        env = 'results/envs_installation/install_toolbox.txt'
    conda: 
        "../envs/toolbox.yaml"
    shell:
        """
        echo $CONDA_PREFIX > {output.env}
        """


rule install_sratool:
    output:
        env = 'results/envs_installation/install_sratool.txt'
    conda: 
        "../envs/sratool.yaml"
    shell:
        """
        echo $CONDA_PREFIX > {output.env}
        """


rule install_trinity:
    output:
        env = 'results/envs_installation/install_trinity.txt'
    conda: 
        "../envs/trinity.yaml"
    shell:
        """
        echo $CONDA_PREFIX > {output.env}
        """


rule install_annotation:
    input:
        lib_embl = config['lib_embl']
    output:
        env = 'results/envs_installation/install_annotation.txt'
    conda: 
        "../envs/annotation.yaml"
    shell:
        """
        # Create symlinks for the RepeatMasker scripts
        for file in $CONDA_PREFIX/share/RepeatMasker/* ; do
            if [[ "$file" == *".pm"* ]] ; then
                file_name=$(echo $file| rev |cut -d '/' -f 1 |rev)
                link=$CONDA_PREFIX/bin/$file_name
                if [[ ! -L "$link" ]] ; then
                    ln -s $file $link
                fi
            fi
        done
        
        # Create a symlink for the repbase library
        lib_embl_symlink=$CONDA_PREFIX/share/RepeatMasker/Libraries/RepeatMaskerLib.embl
        lib_embl=$(pwd)/{input.lib_embl}

        if [[ ! -L "$lib_embl_symlink" ]] ; then
            ln -s $lib_embl $lib_embl_symlink
        fi
        
        # test if RepeatMasker works
        RepeatMasker -help |grep 'RepeatMasker version'
        rmOutToGFF3.pl -h

        # test if maker works
        maker --version

        echo $CONDA_PREFIX > {output.env}
        """

rule install_repeatmodeler:
    input:
        lib_embl = config['lib_embl']
    output:
        env = 'results/envs_installation/install_repeatmodeler.txt'
    conda: 
        "envs/repeatmodeler.yaml"
    shell:
        """
        exogap=$(pwd)

        # Create symlinks for the RepeatMasker scripts
        for file in $CONDA_PREFIX/share/RepeatMasker/* ; do
            if [[ "$file" == *".pm"* ]] ; then
                file_name=$(echo $file| rev |cut -d '/' -f 1 |rev)
                link=$CONDA_DEFACONDA_PREFIXULT_ENV/bin/$file_name
                if [[ ! -L "$link" ]] ; then
                    ln -s $file $link
                fi
            fi
        done

        #  configure repeatmodeler
        cd $CONDA_PREFIX/share/RepeatMasker
        BIN_DIR=$CONDA_PREFIX/bin

        perl ./configure -crossmatch_dir $BIN_DIR \
        -abblast_dir $BIN_DIR \
        -hmmer_dir $BIN_DIR \
        -libdir $CONDA_PREFIX/share/RepeatMasker/Libraries/ \
        -rmblast_dir $BIN_DIR \
        -trf_prgm $BIN_DIR/trf

        # test if RepeatModeler works
        cd $exogap
        # RepeatModeler |grep "RepeatModeler - Model repetitive DNA"
        echo $?

        echo $CONDA_PREFIX > {output.env}
        """