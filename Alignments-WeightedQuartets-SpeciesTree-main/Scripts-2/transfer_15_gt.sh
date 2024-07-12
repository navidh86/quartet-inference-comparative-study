for gene in 100 1000; do
    for bp in 100 1000; do
        for r in {1..10}; do
            from="/mnt/nvme/1605005-1605006/quertet_correction/15-taxon/${gene}gene-${bp}bp/R${r}/all_gt.tre"
            to="/mnt/nvme/1605005-1605006/msc_project/Alignments-WeightedQuartets-SpeciesTree-main/15-taxon/${gene}gene-${bp}bp/estimated-genetrees/R${r}/all_gt.tre"
            cp $from "$to"
            # echo $from
        done
    done
done