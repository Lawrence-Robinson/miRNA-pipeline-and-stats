import re
import numpy

# Define desired location of outfile.
outfile = open('/media/Lawrence/C065-6FD2/rp2/parsed_mirdeep2_output/mirdeep2_data.csv','w') 

# Specify location of miRDeep2 '.mrd' file.
path = "/media/stephen/C065-6FD2/rp2/hsa_mirdeep2/mirdeep_runs/run_16_06_2019_t_03_39_21/output.mrd"

# Write header to outfile
outfile.write("id,sample,group,count,region,score,exp_seq,pri_seq,mirbase_accession,mirbase_id,related_accession,related_id\n")

# Iterate through '.mrd' file, line by line, and extract relevant information.
with open(path,"r") as file:
	for line in file:
		if line.startswith(">"):
			homolog_present = False
			novel = 1
			mirna_id = line[1:].strip('\n')
		elif line.startswith("score total"):		
			score = line.split()[2]
		elif line.startswith("exp"):
			exp = line.split()[1]
			star_start = exp.find('S')
			star_end = len(exp)-exp[::-1].find('S')
			loop_start = exp.find('l')
			loop_end = len(exp)-exp[::-1].find('l')
			mat_start = exp.find('M')
			mat_end = len(exp)-exp[::-1].find('M')
			mat_len = mat_end - mat_start		
		elif line.startswith("pri_seq"):
			pri_seq = line.split()[1]
			exp_star_seq = pri_seq[star_start:star_end]
			exp_loop_seq = pri_seq[loop_start:loop_end]
			exp_mat_seq = pri_seq[mat_start:mat_end]
		elif line.startswith("miRNA with same seed"):
			homolog_present = True
			homolog = line.split()[4]
			homolog = line.split("_")
			a = line.split("\t")[1]
			homolog_id = a[0:a.find("MIMA")-1]
			if (len(homolog) == 10):
				homolog_accession = homolog[4]
				
			elif (len(homolog) == 8):
				homolog_accession = homolog[3]
			else:
				homolog_accession = homolog[5]
		elif line.startswith("oar") and novel == 1:
			novel = 0
			mirbase_id = line[0:line.find("MIMA")-1]
			mirbase = line.split("_")
			if (len(mirbase) == 10):
				mirbase_accession = mirbase[4]
			else:
				mirbase_accession = mirbase[3]
			outfile.write("{0},na,na,0,na,0,na,na,{1},{2},na,na\n".format(mirna_id,mirbase_accession, mirbase_id))
		elif line.startswith("H") or line.startswith("F") or line.startswith("R") :
			if (homolog_present == False):
				homolog_accession = "none"
				homolog_id = "none"
			a = line.split()		
			sample = a[0][0:3]
			if sample.startswith("H"):
				group = "healthy"
			elif sample.startswith("F"):
				group = "failure"
			elif sample.startswith("R"):
				group = "recover"
			else:
				print("Unknown")
			read_start = re.search("\w", a[1]).start()
			read_end = len(a[1]) - re.search("\w", a[1][::-1]).start()
			mean_pos = (read_start + read_end)/2
			if (mean_pos >= star_start and mean_pos <= star_end):
				mirna = "star"
				exp_seq = exp_star_seq
			elif (mean_pos >= mat_start and mean_pos <= mat_end):
				mirna = "mature"
				exp_seq = exp_mat_seq
			elif (mean_pos >= loop_start and mean_pos <= loop_end):
				mirna = "loop"
				exp_seq = exp_loop_seq
			else:
				mirna = "unknown"
				exp_seq = "unknown"
			outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},na,na,{8},{9}\n".format(mirna_id,sample,group,a[0][line.find('x')+1:],mirna,score,exp_seq,pri_seq,homolog_accession, homolog_id))
outfile.close()

