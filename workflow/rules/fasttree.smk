
rule fasttree:
    input: msa=config["work_dir"]+"/msa_template.fa"
    output: tree=config["work_dir"]+"/tree_template.nwk"
    params:
        fasttree_exe=config["fasttree"],
        model=config["ft_model"],
        threshold=config["mask_gappy"],
        tempFile=config["work_dir"]+"/msa.mask.fa",
        tempTree=config["work_dir"]+"/temp_tree.nwk"
    threads: config["num_threads"]
    shell:
        '''
        python3 scripts/reduceLen.py {input.msa} {params.tempFile} {params.threshold}
        export OMP_NUM_THREADS={threads}
        {params.fasttree_exe} {params.model} -fastest {params.tempFile} > {params.tempTree} 
        python3 scripts/resolveTree.py {params.tempTree} {output.tree}
        rm {params.tempFile} {params.tempTree}
        '''

itree = (config["iter_tree"] == "fasttree") and (config["backbone_aln"] == "")
ftree = (config["final_tree"] == "fasttree")
iters = int(config["iterations"])

if (itree and iters >= 2) or (ftree and iters == 1):
    use rule fasttree as fasttree_iter1 with:
        input: msa=config["work_dir"]+"/msa_iter1.fa",
        output: tree=config["work_dir"]+"/tree_iter1.nwk"
if (itree and iters >= 3) or (ftree and iters == 2):
    use rule fasttree as fasttree_iter2 with:
        input: msa=config["work_dir"]+"/msa_iter2.fa",
        output: tree=config["work_dir"]+"/tree_iter2.nwk"
if (itree and iters >= 4) or (ftree and iters == 3):
    use rule fasttree as fasttree_iter3 with:
        input: msa=config["work_dir"]+"/msa_iter3.fa",
        output: tree=config["work_dir"]+"/tree_iter3.nwk"
if (itree and iters >= 5) or (ftree and iters == 4):
    use rule fasttree as fasttree_iter4 with:
        input: msa=config["work_dir"]+"/msa_iter4.fa",
        output: tree=config["work_dir"]+"/tree_iter4.nwk"
if (ftree and iters == 5):
    use rule fasttree as fasttree_iter5 with:
        input: msa=config["work_dir"]+"/msa_iter5.fa",
        output: tree=config["work_dir"]+"/tree_iter5.nwk"

    

if config["backbone_aln"] != "" and config["backbone_tree"] == "":
    use rule fasttree as fasttree_init with:
        input: msa=config["backbone_aln"]
        output: tree=config["work_dir"]+"/backbone_tree.nwk"
