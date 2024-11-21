from sys import argv


def read_regions(infile, pstart,pend):
	f=open(infile)
	f.readline()
	dic={}
	for i in f:
		line=i.split("\t")
		if line[0][3:] not in dic:
			dic[line[0][3:]]=[[int(line[1])-int(pstart), int(line[2])+int(pend)]]
		else:
			dic[line[0][3:]].append([int(line[1])-int(pstart), int(line[2])+int(pend)])
	f.close()
	return dic

def get_Overlap(a,b):
	return max(0,min(a[1],b[1])-max(a[0],b[0]))

def read_vcf(infile, dic):
	f=open(infile)
	for i in f:
		if i.startswith("#"):
			print(i.strip())
		else:
			line=i.split("\t")
			if "chr" in line[0]:
				line[0]=line[0][3:]
			if line[0] in dic:
				for ele in dic[line[0]]:
					if get_Overlap(ele, [int(line[1])-1, int(line[1])+1])>0:
						print(i.strip())
						break
	f.close()

read_vcf(argv[1], read_regions(argv[2], argv[3], argv[4]))
#python3 genome2exome.py vcf regions start_padding end_padding
				
