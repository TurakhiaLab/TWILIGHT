
rule iqtree_iter1:
    input: msa=config["work_dir"]+"/msa_iter1.fa",
    output: tree=config["work_dir"]+"/tree_iter1.nwk"
    params:
        iqtree_exe=config["iqtree"],
        temp=config["work_dir"]
    threads: config["num_threads"]
    shell:
        '''
        {params.iqtree_exe} -s {input.msa} -T {threads}
        mv {params.temp}/msa_iter1.fa.treefile {output}
        rm {params.temp}/msa_iter1.fa.*
        '''

rule iqtree_iter2:
    input: msa=config["work_dir"]+"/msa_iter2.fa",
    output: tree=config["work_dir"]+"/tree_iter2.nwk"
    params:
        iqtree_exe=config["iqtree"],
        temp=config["work_dir"]
    threads: config["num_threads"]
    shell:
        '''
        {params.iqtree_exe} -s {input.msa} -T {threads}
        mv {params.temp}/msa_iter2.fa.treefile {output}
        rm {params.temp}/msa_iter2.fa.*
        '''

rule iqtree_iter3:
    input: msa=config["work_dir"]+"/msa_iter3.fa",
    output: tree=config["work_dir"]+"/tree_iter3.nwk"
    params:
        iqtree_exe=config["iqtree"],
        temp=config["work_dir"]
    threads: config["num_threads"]
    shell:
        '''
        {params.iqtree_exe} -s {input.msa} -T {threads}
        mv {params.temp}/msa_iter2.fa.treefile {output}
        rm {params.temp}/msa_iter2.fa.*
        '''

rule iqtree_iter4:
    input: msa=config["work_dir"]+"/msa_iter4.fa",
    output: tree=config["work_dir"]+"/tree_iter4.nwk"
    params:
        iqtree_exe=config["iqtree"],
        temp=config["work_dir"]
    threads: config["num_threads"]
    shell:
        '''
        {params.iqtree_exe} -s {input.msa} -T {threads}
        mv {params.temp}/msa_iter2.fa.treefile {output}
        rm {params.temp}/msa_iter2.fa.*
        '''

rule iqtree_iter5:
    input: msa=config["work_dir"]+"/msa_iter5.fa",
    output: tree=config["work_dir"]+"/tree_iter5.nwk"
    params:
        iqtree_exe=config["iqtree"],
        temp=config["work_dir"]
    threads: config["num_threads"]
    shell:
        '''
        {params.iqtree_exe} -s {input.msa} -T {threads}
        mv {params.temp}/msa_iter2.fa.treefile {output}
        rm {params.temp}/msa_iter2.fa.*
        '''