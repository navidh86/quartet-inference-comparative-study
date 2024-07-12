for gene in 100 1000; do
    for bp in 100 1000; do
        for r in {1..10}; do
            gt="/mnt/nvme/1605005-1605006/msc_project/Alignments-WeightedQuartets-SpeciesTree-main/15-taxon/${gene}gene-${bp}bp/estimated-genetrees/R${r}/all_gt.tre"
            out="/mnt/nvme/1605005-1605006/msc_project/Alignments-WeightedQuartets-SpeciesTree-main/15-taxon/Output/${gene}gene-${bp}bp-quartets/R${r}-GTF.wqrts"
            # run wqfm on the specified gene tree
            #wqfm_path=/mnt/nvme/1605005-1605006/msc_project/Alignments-WeightedQuartets-SpeciesTree-main/Scripts-2/wQFM-v1.3.jar
            #java -jar ${wqfm_path} -i ${gt} -o ${out} -im gene-trees

            # # find out rf-rates
            # rf_out=


            # java -jar phylonet_v2_4.jar rf -m {} -e {} >> {}
            bash /mnt/nvme/1605005-1605006/msc_project/Alignments-WeightedQuartets-SpeciesTree-main/Scripts-2/compute_weighted_quartet_from_geneTree/quartet-controller.sh ${gt} ${out}
        done
    done
done
