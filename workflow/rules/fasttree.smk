
rule fasttree_iter1:
    input: msa="output/msa_iter1.fa",
    output: tree="output/tree_iter1.nwk"
    params:
        fasttree_exe=config["fasttree"],
        model= "-gtr -nt" if config["ft_model"] == "GTR" else "-nt"
    threads: config["num_threads"]
    shell:
        '''
        export OMP_NUM_THREADS={threads}
        {params.fasttree_exe} {params.model} -fastest {input.msa} > {output.tree} 
        '''

rule fasttree_iter2:
    input: msa="output/msa_iter2.fa",
    output: tree="output/tree_iter2.nwk"
    params:
        fasttree_exe=config["fasttree"],
        model= "-gtr -nt" if config["ft_model"] == "GTR" else "-nt"
    shell:
        '''
        {params.fasttree_exe} {params.model} -fastest {input.msa} > {output.tree} 
        '''

rule fasttree_iter3:
    input: msa="output/msa_iter3.fa",
    output: tree="output/tree_iter3.nwk"
    params:
        fasttree_exe=config["fasttree"],
        model= "-gtr -nt" if config["ft_model"] == "GTR" else "-nt"
    shell:
        '''
        {params.fasttree_exe} {params.model} -fastest {input.msa} > {output.tree} 
        '''

rule fasttree_iter4:
    input: msa="output/msa_iter4.fa",
    output: tree="output/tree_iter4.nwk"
    params:
        fasttree_exe=config["fasttree"],
        model= "-gtr -nt" if config["ft_model"] == "GTR" else "-nt"
    shell:
        '''
        {params.fasttree_exe} {params.model} -fastest {input.msa} > {output.tree} 
        '''

rule fasttree_iter5:
    input: msa="output/msa_iter5.fa",
    output: tree="output/tree_iter5.nwk"
    params:
        fasttree_exe=config["fasttree"],
        model= "-gtr -nt" if config["ft_model"] == "GTR" else "-nt"
    shell:
        '''
        {params.fasttree_exe} {params.model} -fastest {input.msa} > {output.tree} 
        '''


  
        