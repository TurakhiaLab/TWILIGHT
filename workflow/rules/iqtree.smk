
rule iqtree_iter1:
    input: msa="output/msa_iter1.fa",
    output: tree="output/tree_iter1.nwk"
    params:
        iqtree_exe=config["iqtree"]
    threads: config["num_threads"]
    shell:
        '''
        {params.iqtree_exe} -s {input.msa} -T {threads}
        mv output/msa_iter1.fa.treefile {output}
        rm output/msa_iter1.fa.*
        '''

rule iqtree_iter2:
    input: msa="output/msa_iter2.fa",
    output: tree="output/tree_iter2.nwk"
    params:
        iqtree_exe=config["iqtree"]
    threads: config["num_threads"]
    shell:
        '''
        {params.iqtree_exe} -s {input.msa} -T {threads}
        mv output/msa_iter2.fa.treefile {output}
        rm output/msa_iter2.fa.*
        '''

rule iqtree_iter3:
    input: msa="output/msa_iter3.fa",
    output: tree="output/tree_iter3.nwk"
    params:
        iqtree_exe=config["iqtree"]
    threads: config["num_threads"]
    shell:
        '''
        {params.iqtree_exe} -s {input.msa} -T {threads}
        mv output/msa_iter2.fa.treefile {output}
        rm output/msa_iter2.fa.*
        '''

rule iqtree_iter4:
    input: msa="output/msa_iter4.fa",
    output: tree="output/tree_iter4.nwk"
    params:
        iqtree_exe=config["iqtree"]
    threads: config["num_threads"]
    shell:
        '''
        {params.iqtree_exe} -s {input.msa} -T {threads}
        mv output/msa_iter2.fa.treefile {output}
        rm output/msa_iter2.fa.*
        '''

rule iqtree_iter5:
    input: msa="output/msa_iter5.fa",
    output: tree="output/tree_iter5.nwk"
    params:
        iqtree_exe=config["iqtree"]
    threads: config["num_threads"]
    shell:
        '''
        {params.iqtree_exe} -s {input.msa} -T {threads}
        mv output/msa_iter2.fa.treefile {output}
        rm output/msa_iter2.fa.*
        '''