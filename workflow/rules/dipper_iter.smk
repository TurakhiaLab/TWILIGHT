
rule dipper_iter:
    input: msa=config["work_dir"]+"/msa_template.fa",
    output: tree=config["work_dir"]+"/tree_template.nwk"
    params:
        dipper_exe=config["dipper"],
    threads: config["num_threads"]
    shell:
        '''
        {params.dipper_exe} -i m -o t -m 1 -I {input.msa} -O {output.tree}
        '''

iters = int(config["iterations"])

if iters >= 2:
    use rule dipper_iter as dipper_iter1 with:
        input: msa=config["work_dir"]+"/msa_iter1.fa",
        output: tree=config["work_dir"]+"/tree_iter1.nwk"
if iters >= 3:
    use rule dipper_iter as dipper_iter2 with:
        input: msa=config["work_dir"]+"/msa_iter2.fa",
        output: tree=config["work_dir"]+"/tree_iter2.nwk"
if iters >= 4:
    use rule dipper_iter as dipper_iter3 with:
        input: msa=config["work_dir"]+"/msa_iter3.fa",
        output: tree=config["work_dir"]+"/tree_iter3.nwk"
if iters >= 5:
    use rule dipper_iter as dipper_iter4 with:
        input: msa=config["work_dir"]+"/msa_iter4.fa",
        output: tree=config["work_dir"]+"/tree_iter4.nwk"

    
