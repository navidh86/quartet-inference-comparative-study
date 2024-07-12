fout=open('/home/navidh86/Desktop/comparative_study/RF_Rates_Folder/compute_rf/compute_RF-10tax-LowILS-11.txt',"a")
with open("/home/navidh86/Desktop/comparative_study/RF_Rates_Folder/fp,fn,branch/RF-10tax-LowILS.txt","r") as fin:
	lines=fin.readlines()
	for line in lines:
		arr = line.split(" ")
		fn=int(arr[1])
		fp=int(arr[2])
		n=11
		rf=(fn+fp)/(2*n - 6)
		print(arr,fn,fp,rf)
		fout.write(arr[0]+" "+str(rf)+"\n")
		
fout.close()
fin.close()