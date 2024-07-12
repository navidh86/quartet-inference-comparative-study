taxon = 15
if taxon == 10:
	n = 11
else:
	n = 15

if taxon == 15:
	modes = ["100gene-100bp", "100gene-1000bp", "1000gene-100bp", "1000gene-1000bp"]
else:
	modes = ["higher-ILS", "lower-ILS"]
for mode in modes:
	input_file="/mnt/nvme/1605005-1605006/msc_project/Alignments-WeightedQuartets-SpeciesTree-main/Estimated-Species-Trees/FP_FN/FP-FN-{}tax-{}.txt".format(taxon, mode)
	output_file="RF-rates/{}-taxon/RF-{}tax-{}.txt".format(taxon, taxon, mode)
	fout = open(output_file,"w")
	with open(input_file,"r") as fin:
		lines=fin.readlines()
		for line in lines:
			arr = line.split(" ")
			# print("arr", arr)
			fn=int(arr[1])
			fp=int(arr[2])
			rf=(fn+fp)/(2*n - 6)
			print(arr,fn,fp,rf)
			fout.write(arr[0]+" "+str(rf)+"\n")
			
	fout.close()
	fin.close()