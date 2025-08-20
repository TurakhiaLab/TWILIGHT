
rule placement:
    input: 
        backbone_tree=config["work_dir"]+"/backbone_tree.nwk",
        seq_msa=config["work_dir"]+"/placed_template.fa",
        backbone_msa=config["work_dir"]+"/backbone_template.fa"
    output: 
        tree=config["work_dir"]+"/tree_template.nwk"
    params:
        model=config["epa_model"],
        epang_exe=config["epang"],
        gappa_exe=config["gappa"],
        dir=config["work_dir"],
    threads: config["num_threads"]
    shell:
        '''
        {params.epang_exe} --tree {input.backbone_tree} --ref-msa {input.backbone_msa} --query {input.seq_msa} --model {params.model} --outdir {params.dir}
        {params.gappa_exe} examine graft --jplace-path {params.dir}/epa_result.jplace --out-dir {params.dir}
        mv {params.dir}/epa_result.newick {output.tree}
        rm -f {params.dir}/epa_*
        '''

iters = int(config["iterations"])

if iters >= 2:
    use rule placement as placement_iter1 with:
        input: 
            backbone_tree=config["work_dir"]+"/backbone_tree.nwk",
            seq_msa=config["work_dir"]+"/placed_iter1.fa",
            backbone_msa=config["work_dir"]+"/backbone_iter1.fa"
        output: 
            tree=config["work_dir"]+"/tree_iter1.nwk"
if iters >= 3:
    use rule placement as placement_iter2 with:
        input: 
            backbone_tree=config["work_dir"]+"/backbone_tree.nwk",
            seq_msa=config["work_dir"]+"/placed_iter2.fa",
            backbone_msa=config["work_dir"]+"/backbone_iter2.fa"
        output: 
            tree=config["work_dir"]+"/tree_iter2.nwk"
if iters >= 4:
    use rule placement as placement_iter3 with:
        input: 
            backbone_tree=config["work_dir"]+"/backbone_tree.nwk",
            seq_msa=config["work_dir"]+"/placed_iter3.fa",
            backbone_msa=config["work_dir"]+"/backbone_iter3.fa"
        output: 
            tree=config["work_dir"]+"/tree_iter3.nwk"
if iters >= 5:
    use rule placement as placement_iter4 with:
        input: 
            backbone_tree=config["work_dir"]+"/backbone_tree.nwk",
            seq_msa=config["work_dir"]+"/placed_iter4.fa",
            backbone_msa=config["work_dir"]+"/backbone_iter4.fa"
        output: 
            tree=config["work_dir"]+"/tree_iter4.nwk"
