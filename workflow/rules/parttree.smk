
rule parttree:
    input: config["sequences"],
    output: config["work_dir"]+"/tree_iter0.nwk"
    params:
        mafft_exe=config["mafft"],
        tempFile=config["sequences"]+".tree"
    threads: config["num_threads"]
    shell:
        '''
        {params.mafft_exe} --retree 0 --treeout --parttree --reorder --thread {threads} {input} > mafft.out
        python3 workflow/scripts/mafft2nwk.py {params.tempFile} {input} {output} --parttree
        rm mafft.out
        rm {params.tempFile}
        '''