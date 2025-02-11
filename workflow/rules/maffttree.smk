
rule maffttree:
    input: config["sequences"],
    output: config["work_dir"]+"/tree_iter0.nwk"
    params:
        mafft_exe=config["mafft"],
        tempFile=config["sequences"]+".tree"
    threads: config["num_threads"]
    shell:
        '''
        {params.mafft_exe} --retree 0 --treeout --reorder --thread {threads} {input} > mafft.out
        python3 scripts/mafft2nwk.py {params.tempFile} {input} {output}
        rm mafft.out
        rm {params.tempFile}
        '''