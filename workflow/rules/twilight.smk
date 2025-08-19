if config["backbone_aln"] == "":
    rule twilight_denovo:
        input:
            seq=config["sequences"],
            tree=config["work_dir"]+"/tree_template.nwk"
        output:
            msa=config["work_dir"]+"/msa_template.fa"
        params:
            twilight_exe=config["twilight"],
            max_subtree=config["max_subtree"],  
            rgc=config["rgc"],
            gop=config["gapopen"],
            gex=config["gapextend"],
            matrix="" if config["matrix"]=="" else "-x "+config["matrix"]

        threads: config["num_threads"]
        shell:
            '''
            {params.twilight_exe} -i {input.seq} -t {input.tree} -o {output.msa} -C {threads} \
            -m {params.max_subtree} -r {params.rgc} --gap-open {params.gop} --gap-extend {params.gex} {params.matrix} -v
            '''
    use rule twilight_denovo as twilight_iter1 with:
        input:
            seq=config["sequences"],
            tree=config["work_dir"]+"/tree_iter0.nwk"
        output:
            msa=config["work_dir"]+"/msa_iter1.fa"
    if int(config["iterations"]) >= 2: 
        use rule twilight_denovo as twilight_iter2 with:
            input:
                seq=config["sequences"],
                tree=config["work_dir"]+"/tree_iter1.nwk"
            output:
                msa=config["work_dir"]+"/msa_iter2.fa"

    if int(config["iterations"]) >= 3: 
        use rule twilight_denovo as twilight_iter3 with:
            input:
                seq=config["sequences"],
                tree=config["work_dir"]+"/tree_iter2.nwk"
            output:
                msa=config["work_dir"]+"/msa_iter3.fa"

    if int(config["iterations"]) >= 4: 
        use rule twilight_denovo as twilight_iter4 with:
            input:
                seq=config["sequences"],
                tree=config["work_dir"]+"/tree_iter3.nwk"
            output:
                msa=config["work_dir"]+"/msa_iter4.fa"

    if int(config["iterations"]) >= 5: 
        use rule twilight_denovo as twilight_iter5 with:
            input:
                seq=config["sequences"],
                tree=config["work_dir"]+"/tree_iter4.nwk"
            output:
                msa=config["work_dir"]+"/msa_iter5.fa"


else:
    if config["backbone_tree"] != "":
        rule backbone_tree:
            input: tree_in=config["backbone_tree"]
            output: tree_out=config["work_dir"]+"/backbone_tree.nwk"
            params: temp_tree=config["work_dir"]+"/temp_tree.nwk"
            shell:
                '''
                python3 scripts/resolveTree.py {input.tree_in} {params.temp_tree}
                mv {params.temp_tree} {output.tree_out}
                '''

    rule generate_final_msa:
        input: 
            placed=config["work_dir"]+"/placed_iter"+ str(config["iterations"]) + ".fa",
            backbone=config["work_dir"]+"/backbone_iter"+ str(config["iterations"]) + ".fa",
        output:
            msa=config["work_dir"]+"/msa_iter" + str(config["iterations"]) + ".fa"
        shell:
            '''
            cat {input.placed} {input.backbone} > {output.msa}
            '''

    rule twilight_place_at_root:
        input:
            seq=config["sequences"],
            msa=config["backbone_aln"]
        output:
            placed=config["work_dir"]+"/placed_iter1.fa",
            backbone=config["work_dir"]+"/backbone_iter1.fa"
        params:
            twilight_exe=config["twilight"],
            max_subtree=config["max_subtree"],  
            rgc=config["rgc"],
            gop=config["gapopen"],
            gex=config["gapextend"],
            matrix="" if config["matrix"]=="" else "-x "+config["matrix"],
            tempFile=config["work_dir"]+"/placement_temp.fa",
            tempDir=config["work_dir"]+"/twilight_temp",
            place_msa=config["work_dir"]+"/twilight_temp/"+os.path.basename(config["sequences"])+".final.aln",
            backbone_msa=config["work_dir"]+"/twilight_temp/"+os.path.basename(config["backbone_aln"])+".final.aln"
        threads: config["num_threads"]
        shell:
            '''
            {params.twilight_exe} -i {input.seq} -a {input.msa} -o {params.tempFile} -C {threads} -k -d {params.tempDir} \
            -m {params.max_subtree} -r {params.rgc} --gap-open {params.gop} --gap-extend {params.gex} {params.matrix} -v
            mv {params.place_msa} {output.placed}
            mv {params.backbone_msa} {output.backbone}
            rm -rf {params.tempDir}
            rm {params.tempFile}
            '''

    rule twilight_place_at_tips:
        input:
            seq=config["sequences"],
            msa=config["backbone_aln"],
            tree=config["work_dir"]+"/tree_template.fa"
        output:
            placed=config["work_dir"]+"/placed_template.fa",
            backbone=config["work_dir"]+"/backbone_template.fa"
        params:
            twilight_exe=config["twilight"],
            max_subtree=config["max_subtree"],  
            rgc=config["rgc"],
            gop=config["gapopen"],
            gex=config["gapextend"],
            matrix="" if config["matrix"]=="" else "-x "+config["matrix"],
            tempFile=config["work_dir"]+"/placement_temp.fa",
            tempDir=config["work_dir"]+"/twilight_temp",
            place_msa=config["work_dir"]+"/twilight_temp/"+os.path.basename(config["sequences"])+".final.aln",
            backbone_msa=config["work_dir"]+"/twilight_temp/"+os.path.basename(config["backbone_aln"])+".final.aln"
        threads: config["num_threads"]
        shell:
            '''
            {params.twilight_exe} -i {input.seq} -a {input.msa} -t {input.tree} -o {params.tempFile} -C {threads} -k -d {params.tempDir} \
            -m {params.max_subtree} -r {params.rgc} --gap-open {params.gop} --gap-extend {params.gex} {params.matrix} -v
            mv {params.place_msa} {output.placed}
            mv {params.backbone_msa} {output.backbone}
            rm -rf {params.tempDir}
            rm {params.tempFile}
            '''
    
    if int(config["iterations"]) >= 2: 
        use rule twilight_place_at_tips as twilight_tips_iter2 with:
            input:
                seq=config["sequences"],
                msa=config["backbone_aln"],
                tree=config["work_dir"]+"/tree_iter1.nwk"
            output:
                placed=config["work_dir"]+"/placed_iter2.fa",
                backbone=config["work_dir"]+"/backbone_iter2.fa"
    if int(config["iterations"]) >= 3: 
        use rule twilight_place_at_tips as twilight_tips_iter3 with:
            input:
                seq=config["sequences"],
                msa=config["backbone_aln"],
                tree=config["work_dir"]+"/tree_iter2.nwk"
            output:
                placed=config["work_dir"]+"/placed_iter3.fa",
                backbone=config["work_dir"]+"/backbone_iter3.fa"
    if int(config["iterations"]) >= 4: 
        use rule twilight_place_at_tips as twilight_tips_iter4 with:
            input:
                seq=config["sequences"],
                msa=config["backbone_aln"],
                tree=config["work_dir"]+"/tree_iter3.nwk"
            output:
                placed=config["work_dir"]+"/placed_iter4.fa",
                backbone=config["work_dir"]+"/backbone_iter4.fa"
    if int(config["iterations"]) >= 5: 
        use rule twilight_place_at_tips as twilight_tips_iter5 with:
            input:
                seq=config["sequences"],
                msa=config["backbone_aln"],
                tree=config["work_dir"]+"/tree_iter4.nwk"
            output:
                placed=config["work_dir"]+"/placed_iter5.fa",
                backbone=config["work_dir"]+"/backbone_iter5.fa"