
if config["backbone_aln"] == "":
    rule mashtree:
        input: config["sequences"],
        output: config["work_dir"]+"/tree_iter0.nwk"
        params:
            mashtree_exe=config["mashtree"],
            tempDir=config["sequences"]+".tempDir"
        threads: config["num_threads"]
        shell:
            '''
            bash scripts/mashtree.sh {input} {params.tempDir}
    		{params.mashtree_exe} --numcpus {threads} --outtree {output} {params.tempDir}/*.fa
            rm -rf {params.tempDir}
            '''
