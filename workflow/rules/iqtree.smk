
rule iqtree:
    input: msa=config["work_dir"]+"/msa_iter" + str(config["iterations"]) + ".fa",
    output: tree=config["work_dir"]+"/tree_iter" + str(config["iterations"]) + ".nwk"
    params:
        iqtree_exe=config["iqtree"],
        temp=config["work_dir"],
        model="",
        threshold=config["mask_gappy"],
        tempFile=config["work_dir"]+"/msa.mask.fa"
    threads: config["num_threads"]
    shell:
        '''
        python3 scripts/reduceLen.py {input.msa} {params.tempFile} {params.threshold}
        {params.iqtree_exe} -s {params.tempFile} {params.model} --threads-max {threads}
        mv {params.temp}/msa_iter1.mask.fa.treefile {output}
        rm {params.temp}/msa_iter1.mask.fa.*
        rm {params.tempFile}
        '''