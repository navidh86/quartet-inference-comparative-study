# /mnt/nvme/1605005-1605006/msc_project/Alignments-WeightedQuartets-SpeciesTree-main/10-taxon/higher-ILS/estimated-genetrees/R1/1/RAxML_bipartitions.final.f200
modes = ['higher-ILS', 'lower-ILS']
for mode in modes:
    for r in range(1, 21):
        output_file = "/mnt/nvme/1605005-1605006/msc_project/Alignments-WeightedQuartets-SpeciesTree-main/10-taxon/{}/estimated-genetrees/R{}/all_gt.tre".format(mode, r)
        with open(output_file, "w") as fout:
            for gene in range(1, 201):
                input_file = "/mnt/nvme/1605005-1605006/msc_project/Alignments-WeightedQuartets-SpeciesTree-main/10-taxon/{}/estimated-genetrees/R{}/{}/RAxML_bipartitions.final.f200".format(mode, r, gene)
                with open(input_file, "r") as fin:
                    tree = fin.readline()
                    fout.write(tree)

