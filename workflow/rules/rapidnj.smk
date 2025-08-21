
rule rapidnj:
    input: msa=config["work_dir"]+"/msa_template.fa",
    output: tree=config["work_dir"]+"/tree_template.nwk"
    params:
        rapidnj_exe=config["rapidnj"],
        model=config["ft_model"],
        threshold=config["mask_gappy"],
        tempFile=config["work_dir"]+"/msa.mask.fa"
    threads: config["num_threads"]
    shell:
        '''
        {params.rapidnj_exe} {input.msa} -i fa -o t -x {output.tree} -c {threads}
        '''

iters = int(config["iterations"])

if iters >= 2:
    use rule rapidnj as rapidnj_iter1 with:
        input: msa=config["work_dir"]+"/msa_iter1.fa",
        output: tree=config["work_dir"]+"/tree_iter1.nwk"
if iters >= 3:
    use rule rapidnj as rapidnj_iter2 with:
        input: msa=config["work_dir"]+"/msa_iter2.fa",
        output: tree=config["work_dir"]+"/tree_iter2.nwk"
if iters >= 4:
    use rule rapidnj as rapidnj_iter3 with:
        input: msa=config["work_dir"]+"/msa_iter3.fa",
        output: tree=config["work_dir"]+"/tree_iter3.nwk"
if iters >= 5:
    use rule rapidnj as rapidnj_iter4 with:
        input: msa=config["work_dir"]+"/msa_iter4.fa",
        output: tree=config["work_dir"]+"/tree_iter4.nwk"

    
