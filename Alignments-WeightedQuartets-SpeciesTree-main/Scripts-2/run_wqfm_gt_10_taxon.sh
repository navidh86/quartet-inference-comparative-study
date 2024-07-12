cd Scripts-2
for mode in higher-ILS lower-ILS; do
    for r in {1..20}; do
        gt="/mnt/nvme/1605005-1605006/msc_project/Alignments-WeightedQuartets-SpeciesTree-main/10-taxon/${mode}/estimated-genetrees/R${r}/all_gt.tre"
        out="/mnt/nvme/1605005-1605006/msc_project/Alignments-WeightedQuartets-SpeciesTree-main/Estimated-Species-Trees/10-taxa-species-trees-estimated/${mode}/R${r}-wQFM-GTF.tre"
        # run wqfm on the specified gene tree
        wqfm_path="/mnt/nvme/1605005-1605006/msc_project/Alignments-WeightedQuartets-SpeciesTree-main/Scripts-2/wQFM-v1.3.jar"
        java -jar ${wqfm_path} -i ${gt} -o ${out} -im gene-trees
    done
done
cd ..
