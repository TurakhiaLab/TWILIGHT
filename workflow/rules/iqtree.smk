
rule iqtree_iter1:
    input: msa=config["work_dir"]+"/msa_iter1.fa",
    output: tree=config["work_dir"]+"/tree_iter1.nwk"
    params:
        iqtree_exe=config["iqtree"],
        temp=config["work_dir"],
        model="" if config["iq_model"] == "AUTO" else "-m " + config["iq_model"],
        threshold=config["mask_gappy"],
        tempFile=config["work_dir"]+"/msa_iter1.mask.fa"
    threads: config["num_threads"]
    shell:
        '''
        python3 workflow/scripts/reduceLen.py {input.msa} {params.tempFile} {params.threshold}
        {params.iqtree_exe} -s {params.tempFile} {params.model} --threads-max {threads}
        mv {params.temp}/msa_iter1.mask.fa.treefile {output}
        rm {params.temp}/msa_iter1.mask.fa.*
        rm {params.tempFile}
        '''

rule iqtree_iter2:
    input: msa=config["work_dir"]+"/msa_iter2.fa",
    output: tree=config["work_dir"]+"/tree_iter2.nwk"
    params:
        iqtree_exe=config["iqtree"],
        temp=config["work_dir"],
        model="" if config["iq_model"] == "AUTO" else "-m " + config["iq_model"],
        threshold=config["mask_gappy"],
        tempFile=config["work_dir"]+"/msa_iter2.mask.fa"
    threads: config["num_threads"]
    shell:
        '''
        python3 workflow/scripts/reduceLen.py {input.msa} {params.tempFile} {params.threshold}
        {params.iqtree_exe} -s {params.tempFile} {params.model} --threads-max {threads}
        mv {params.temp}/msa_iter2.mask.fa.treefile {output}
        rm {params.temp}/msa_iter2.mask.fa.*
        rm {params.tempFile}
        '''

rule iqtree_iter3:
    input: msa=config["work_dir"]+"/msa_iter3.fa",
    output: tree=config["work_dir"]+"/tree_iter3.nwk"
    params:
        iqtree_exe=config["iqtree"],
        temp=config["work_dir"],
        model="" if config["iq_model"] == "AUTO" else "-m " + config["iq_model"],
        threshold=config["mask_gappy"],
        tempFile=config["work_dir"]+"/msa_iter3.mask.fa"
    threads: config["num_threads"]
    shell:
        '''
        python3 workflow/scripts/reduceLen.py {input.msa} {params.tempFile} {params.threshold}
        {params.iqtree_exe} -s {params.tempFile} {params.model} --threads-max {threads}
        mv {params.temp}/msa_iter3.mask.fa.treefile {output}
        rm {params.temp}/msa_iter3.mask.fa.*
        rm {params.tempFile}
        '''

rule iqtree_iter4:
    input: msa=config["work_dir"]+"/msa_iter4.fa",
    output: tree=config["work_dir"]+"/tree_iter4.nwk"
    params:
        iqtree_exe=config["iqtree"],
        temp=config["work_dir"],
        model="" if config["iq_model"] == "AUTO" else "-m " + config["iq_model"],
        threshold=config["mask_gappy"],
        tempFile=config["work_dir"]+"/msa_iter4.mask.fa"
    threads: config["num_threads"]
    shell:
        '''
        python3 workflow/scripts/reduceLen.py {input.msa} {params.tempFile} {params.threshold}
        {params.iqtree_exe} -s {params.tempFile} {params.model} --threads-max {threads}
        mv {params.temp}/msa_iter4.mask.fa.treefile {output}
        rm {params.temp}/msa_iter4.mask.fa.*
        rm {params.tempFile}
        '''

rule iqtree_iter5:
    input: msa=config["work_dir"]+"/msa_iter5.fa",
    output: tree=config["work_dir"]+"/tree_iter5.nwk"
    params:
        iqtree_exe=config["iqtree"],
        temp=config["work_dir"],
        model="" if config["iq_model"] == "AUTO" else "-m " + config["iq_model"],
        threshold=config["mask_gappy"],
        tempFile=config["work_dir"]+"/msa_iter5.mask.fa"
    threads: config["num_threads"]
    shell:
        '''
        python3 workflow/scripts/reduceLen.py {input.msa} {params.tempFile} {params.threshold}
        {params.iqtree_exe} -s {params.tempFile} {params.model} --threads-max {threads}
        mv {params.temp}/msa_iter5.mask.fa.treefile {output}
        rm {params.temp}/msa_iter5.mask.fa.*
        rm {params.tempFile}
        '''