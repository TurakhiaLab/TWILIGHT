
rule raxml_iter1:
    input: msa=config["work_dir"]+"/msa_iter1.fa",
    output: tree=config["work_dir"]+"/tree_iter1.nwk"
    params: 
        raxml_exe=config["raxml"],
        model="" if config["rx_model"] != "JC69" and config["rx_model"] != "K80" and config["rx_model"] != "HKY85" else "--"+config["rx_model"]
    threads: config["num_threads"]
    shell:
        '''
        {params.raxml_exe} -s {input.msa} -m GTRCAT -n raxml.tree -T {threads} -p 235813 {params.model}
        mv RAxML_bestTree.raxml.tree {output}
        rm *.raxml.tree
        '''

rule raxml_iter2:
    input: msa=config["work_dir"]+"/msa_iter2.fa",
    output: tree=config["work_dir"]+"/tree_iter2.nwk"
    params: 
        raxml_exe=config["raxml"],
        model="" if config["rx_model"] != "JC69" and config["rx_model"] != "K80" and config["rx_model"] != "HKY85" else "--"+config["rx_model"]
    threads: config["num_threads"]
    shell:
        '''
        {params.raxml_exe} -s {input.msa} -m GTRCAT -n raxml.tree -T {threads} -p 235813 {params.model}
        mv RAxML_bestTree.raxml.tree {output}
        rm *.raxml.tree
        '''

rule raxml_iter3:
    input: msa=config["work_dir"]+"/msa_iter3.fa",
    output: tree=config["work_dir"]+"/tree_iter3.nwk"
    params: 
        raxml_exe=config["raxml"],
        model="" if config["rx_model"] != "JC69" and config["rx_model"] != "K80" and config["rx_model"] != "HKY85" else "--"+config["rx_model"]
    threads: config["num_threads"]
    shell:
        '''
        {params.raxml_exe} -s {input.msa} -m GTRCAT -n raxml.tree -T {threads} -p 235813 {params.model}
        mv RAxML_bestTree.raxml.tree {output}
        rm *.raxml.tree
        '''

rule raxml_iter4:
    input: msa=config["work_dir"]+"/msa_iter4.fa",
    output: tree=config["work_dir"]+"/tree_iter4.nwk"
    params: 
        raxml_exe=config["raxml"],
        model="" if config["rx_model"] != "JC69" and config["rx_model"] != "K80" and config["rx_model"] != "HKY85" else "--"+config["rx_model"]
    threads: config["num_threads"]
    shell:
        '''
        {params.raxml_exe} -s {input.msa} -m GTRCAT -n raxml.tree -T {threads} -p 235813 {params.model}
        mv RAxML_bestTree.raxml.tree {output}
        rm *.raxml.tree
        '''

rule raxml_iter5:
    input: msa=config["work_dir"]+"/msa_iter5.fa",
    output: tree=config["work_dir"]+"/tree_iter5.nwk"
    params: 
        raxml_exe=config["raxml"],
        model="" if config["rx_model"] != "JC69" and config["rx_model"] != "K80" and config["rx_model"] != "HKY85" else "--"+config["rx_model"]
    threads: config["num_threads"]
    shell:
        '''
        {params.raxml_exe} -s {input.msa} -m GTRCAT -n raxml.tree -T {threads} -p 235813 {params.model}
        mv RAxML_bestTree.raxml.tree {output}
        rm *.raxml.tree
        '''