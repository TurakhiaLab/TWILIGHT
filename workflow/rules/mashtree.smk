
rule mashtree:
    input: config["sequences"],
    output: config["work_dir"]+"/tree_iter0.nwk"
    params:
        mashtree_exe=config["mashtree"],
        tempDir=config["sequences"]+".tempDir"
    threads: config["num_threads"]
    shell:
        '''
        bash workflow/scripts/mashtree.sh {input} {params.tempDir}
		{params.mashtree_exe} --mindepth 0 --numcpus {threads} --outtree {output} {params.tempDir}/*.fa
        rm -rf {params.tempDir}
        '''
