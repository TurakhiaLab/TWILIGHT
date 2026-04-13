
if config["backbone_aln"] == "":
    rule dipper_init:
        input: config["sequences"],
        output: config["work_dir"]+"/tree_iter0.nwk"
        params:
            dipper_exe=config["dipper"],
            tempFile=config["sequences"]+".tree"
        threads: config["num_threads"]
        shell:
            '''
            {params.dipper_exe} -i r -o t -m 1 -I {input} -O {output}
            '''