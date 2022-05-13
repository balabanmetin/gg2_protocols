
APPLES_F=0.2


# checkpoint here
#rule


rule create_apples_db:
    input: t="data/single_copy_seqs.nwk", a="data/single_copy_seqs.fa"
    output: "data/apples.dtb"
    params: f=APPLES_F
    shell: 
        '''
            build_applesdtb.py -o {output} -f {params.f} -T 1 -s {input.a} -t {input.t}
        '''


rule single_copy_bb:
    input: e="data/ent_rep.csv", t="data/backbone.nwk", a="data/backbone_align.fa"
    output: t="data/single_copy_seqs.nwk", a="data/single_copy_seqs.fa"
    shell: 
        '''
            awk '$2 == 0' {input.e} | cut -f1 > data/single_copy_seqs.txt
            seqkit grep -f <(sed "s/$/_1/g" data/single_copy_seqs.txt) -w 0 {input.a} | sed "s/_1//g" | seqkit rmdup -s -o {output.a} -i -w 0 -D data/single_copy_rm_map.txt
            grep ">" {output.a} | sed "s/>//g" > data/single_copy_seqs_dedupe.txt   
            nw_prune -v <(nw_topology -bI {input.t}) `cat data/single_copy_seqs_dedupe.txt` > {output.t}
        '''


rule entropy_report:
    input: "data/backbone_align.fa"
    output: "data/ent_rep.csv"
    shell: 
        '''
            python scripts/multicopy_consensus_renewed.py < {input} > {output}
        '''
