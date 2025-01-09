
rule twilight_iter1:
    input:
        seq=config["sequences"],
        tree="output/tree_iter0.nwk"
    output:
        msa="output/msa_iter1.fa"
    params:
        twilight_exe=config["twilight"],
        psgop=config["psgop"],           
        max_group=config["max_group"],   
        max_subtree=config["max_subtree"],  
        rgc=config["rgc"]               
    threads: config["num_threads"]
    shell:
        '''
        {params.twilight_exe} -i {input.seq} -t {input.tree} -o {output.msa} -C {threads} -p {params.psgop} -g {params.max_group} -m {params.max_subtree} -r {params.rgc}
        '''

rule twilight_iter2:
    input:
        seq=config["sequences"],
        tree="output/tree_iter1.nwk"
    output:
        msa="output/msa_iter2.fa"
    params:
        twilight_exe=config["twilight"],
        psgop=config["psgop"],           
        max_group=config["max_group"],   
        max_subtree=config["max_subtree"],  
        rgc=config["rgc"]               
    threads: config["num_threads"]
    shell:
        '''
        {params.twilight_exe} -i {input.seq} -t {input.tree} -o {output.msa} -C {threads} -p {params.psgop} -g {params.max_group} -m {params.max_subtree} -r {params.rgc}
        '''
   
rule twilight_iter3:
    input:
        seq=config["sequences"],
        tree="output/tree_iter2.nwk"
    output:
        msa="output/msa_iter3.fa"
    params:
        twilight_exe=config["twilight"],
        psgop=config["psgop"],           
        max_group=config["max_group"],   
        max_subtree=config["max_subtree"],  
        rgc=config["rgc"]               
    threads: config["num_threads"]
    shell:
        '''
        {params.twilight_exe} -i {input.seq} -t {input.tree} -o {output.msa} -C {threads} -p {params.psgop} -g {params.max_group} -m {params.max_subtree} -r {params.rgc}
        '''

rule twilight_iter4:
    input:
        seq=config["sequences"],
        tree="output/tree_iter3.nwk"
    output:
        msa="output/msa_iter4.fa"
    params:
        twilight_exe=config["twilight"],
        psgop=config["psgop"],           
        max_group=config["max_group"],   
        max_subtree=config["max_subtree"],  
        rgc=config["rgc"]               
    threads: config["num_threads"]
    shell:
        '''
        {params.twilight_exe} -i {input.seq} -t {input.tree} -o {output.msa} -C {threads} -p {params.psgop} -g {params.max_group} -m {params.max_subtree} -r {params.rgc}
        '''

rule twilight_iter5:
    input:
        seq=config["sequences"],
        tree="output/tree_iter4.nwk"
    output:
        msa="output/msa_iter5.fa"
    params:
        twilight_exe=config["twilight"],
        psgop=config["psgop"],           
        max_group=config["max_group"],   
        max_subtree=config["max_subtree"],  
        rgc=config["rgc"]               
    threads: config["num_threads"]
    shell:
        '''
        {params.twilight_exe} -i {input.seq} -t {input.tree} -o {output.msa} -C {threads} -p {params.psgop} -g {params.max_group} -m {params.max_subtree} -r {params.rgc}
        '''