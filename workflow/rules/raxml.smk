
rule raxml_iter1:
    input: msa=config["work_dir"]+"/msa_iter1.fa",
    output: tree=config["work_dir"]+"/tree_iter1.nwk"
    params: 
        raxml_exe=config["raxml"],
        model=config["rx_model"],
        threshold=config["mask_gappy"],
        tempFile=config["work_dir"]+"/msa_iter1.mask.fa"
    threads: config["num_threads"]
    shell:
        '''
        python3 scripts/reduceLen.py {input.msa} {params.tempFile} {params.threshold}
        {params.raxml_exe} -s {params.tempFile} -m {params.model} -n raxml.tree -T {threads} -p 235813
        mv RAxML_bestTree.raxml.tree {output}
        rm *.raxml.tree
        rm {params.tempFile}
        '''

rule raxml_iter2:
    input: msa=config["work_dir"]+"/msa_iter2.fa",
    output: tree=config["work_dir"]+"/tree_iter2.nwk"
    params: 
        raxml_exe=config["raxml"],
        model=config["rx_model"],
        threshold=config["mask_gappy"],
        tempFile=config["work_dir"]+"/msa_iter2.mask.fa"
    threads: config["num_threads"]
    shell:
        '''
        python3 scripts/reduceLen.py {input.msa} {params.tempFile} {params.threshold}
        {params.raxml_exe} -s {params.tempFile} -m {params.model} -n raxml.tree -T {threads} -p 235813
        mv RAxML_bestTree.raxml.tree {output}
        rm *.raxml.tree
        rm {params.tempFile}
        '''

rule raxml_iter3:
    input: msa=config["work_dir"]+"/msa_iter3.fa",
    output: tree=config["work_dir"]+"/tree_iter3.nwk"
    params: 
        raxml_exe=config["raxml"],
        model=config["rx_model"],
        threshold=config["mask_gappy"],
        tempFile=config["work_dir"]+"/msa_iter3.mask.fa"
    threads: config["num_threads"]
    shell:
        '''
        python3 scripts/reduceLen.py {input.msa} {params.tempFile} {params.threshold}
        {params.raxml_exe} -s {params.tempFile} -m {params.model} -n raxml.tree -T {threads} -p 235813
        mv RAxML_bestTree.raxml.tree {output}
        rm *.raxml.tree
        rm {params.tempFile}
        '''

rule raxml_iter4:
    input: msa=config["work_dir"]+"/msa_iter4.fa",
    output: tree=config["work_dir"]+"/tree_iter4.nwk"
    params: 
        raxml_exe=config["raxml"],
        model=config["rx_model"],
        threshold=config["mask_gappy"],
        tempFile=config["work_dir"]+"/msa_iter4.mask.fa"
    threads: config["num_threads"]
    shell:
        '''
        python3 scripts/reduceLen.py {input.msa} {params.tempFile} {params.threshold}
        {params.raxml_exe} -s {params.tempFile} -m {params.model} -n raxml.tree -T {threads} -p 235813
        mv RAxML_bestTree.raxml.tree {output}
        rm *.raxml.tree
        rm {params.tempFile}
        '''

rule raxml_iter5:
    input: msa=config["work_dir"]+"/msa_iter5.fa",
    output: tree=config["work_dir"]+"/tree_iter5.nwk"
    params: 
        raxml_exe=config["raxml"],
        model=config["rx_model"],
        threshold=config["mask_gappy"],
        tempFile=config["work_dir"]+"/msa_iter5.mask.fa"
    threads: config["num_threads"]
    shell:
        '''
        python3 scripts/reduceLen.py {input.msa} {params.tempFile} {params.threshold}
        {params.raxml_exe} -s {params.tempFile} -m {params.model} -n raxml.tree -T {threads} -p 235813
        mv RAxML_bestTree.raxml.tree {output}
        rm *.raxml.tree
        rm {params.tempFile}
        '''