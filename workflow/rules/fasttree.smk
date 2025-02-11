
rule fasttree_iter1:
    input: msa=config["work_dir"]+"/msa_iter1.fa",
    output: tree=config["work_dir"]+"/tree_iter1.nwk"
    params:
        fasttree_exe=config["fasttree"],
        model= "-gtr -nt" if config["ft_model"] == "GTR" else "-nt",
        threshold=config["mask_gappy"],
        tempFile=config["work_dir"]+"/msa_iter1.mask.fa"
    threads: config["num_threads"]
    shell:
        '''
        python3 scripts/reduceLen.py {input.msa} {params.tempFile} {params.threshold}
        export OMP_NUM_THREADS={threads}
        {params.fasttree_exe} {params.model} -fastest {params.tempFile} > {output.tree} 
        rm {params.tempFile}
        '''

rule fasttree_iter2:
    input: msa=config["work_dir"]+"/msa_iter2.fa",
    output: tree=config["work_dir"]+"/tree_iter2.nwk"
    params:
        fasttree_exe=config["fasttree"],
        model= "-gtr -nt" if config["ft_model"] == "GTR" else "-nt",
        threshold=config["mask_gappy"],
        tempFile=config["work_dir"]+"/msa_iter2.mask.fa"
    shell:
        '''
        python3 scripts/reduceLen.py {input.msa} {params.tempFile} {params.threshold}
        {params.fasttree_exe} {params.model} -fastest {params.tempFile} > {output.tree} 
        rm {params.tempFile}
        '''

rule fasttree_iter3:
    input: msa=config["work_dir"]+"/msa_iter3.fa",
    output: tree=config["work_dir"]+"/tree_iter3.nwk"
    params:
        fasttree_exe=config["fasttree"],
        model= "-gtr -nt" if config["ft_model"] == "GTR" else "-nt",
        threshold=config["mask_gappy"],
        tempFile=config["work_dir"]+"/msa_iter3.mask.fa"
    shell:
        '''
        python3 scripts/reduceLen.py {input.msa} {params.tempFile} {params.threshold}
        {params.fasttree_exe} {params.model} -fastest {params.tempFile} > {output.tree} 
        rm {params.tempFile}
        '''

rule fasttree_iter4:
    input: msa=config["work_dir"]+"/msa_iter4.fa",
    output: tree=config["work_dir"]+"/tree_iter4.nwk"
    params:
        fasttree_exe=config["fasttree"],
        model= "-gtr -nt" if config["ft_model"] == "GTR" else "-nt",
        threshold=config["mask_gappy"],
        tempFile=config["work_dir"]+"/msa_iter4.mask.fa"
    shell:
        '''
        python3 scripts/reduceLen.py {input.msa} {params.tempFile} {params.threshold}
        {params.fasttree_exe} {params.model} -fastest {params.tempFile} > {output.tree} 
        rm {params.tempFile}
        '''

rule fasttree_iter5:
    input: msa=config["work_dir"]+"/msa_iter5.fa",
    output: tree=config["work_dir"]+"/tree_iter5.nwk"
    params:
        fasttree_exe=config["fasttree"],
        model= "-gtr -nt" if config["ft_model"] == "GTR" else "-nt",
        threshold=config["mask_gappy"],
        tempFile=config["work_dir"]+"/msa_iter5.mask.fa"
    shell:
        '''
        python3 scripts/reduceLen.py {input.msa} {params.tempFile} {params.threshold}
        {params.fasttree_exe} {params.model} -fastest {params.tempFile} > {output.tree} 
        rm {params.tempFile}
        '''


  
        