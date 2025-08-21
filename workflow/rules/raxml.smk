
rule raxml:
    input: msa=config["work_dir"]+"/msa_iter" + str(config["iterations"]) + ".fa",
    output: tree=config["work_dir"]+"/tree_iter" + str(config["iterations"]) + ".nwk"
    params: 
        raxml_exe=config["raxml"],
        model=config["rx_model"],
        threshold=config["mask_gappy"],
        tempFile=config["work_dir"]+"/msa.mask.fa"
    threads: config["num_threads"]
    shell:
        '''
        python3 scripts/reduceLen.py {input.msa} {params.tempFile} {params.threshold}
        {params.raxml_exe} -s {params.tempFile} -m {params.model} -n raxml.tree -T {threads} -p 235813
        mv RAxML_bestTree.raxml.tree {output}
        rm *.raxml.tree
        rm {params.tempFile}
        '''