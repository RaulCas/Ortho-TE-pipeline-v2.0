

inserts=open('insertsize.txt', 'r')
LTR=open('lengths.txt', 'r')
output=open('classified_events.txt', 'a')

loop=LTR.readlines()

output.write('Element'+'\t'+'ortho_distance'+'\t'+'Real_distance'+'\t'+'Ratio_True/ortho'+'\t'+'Classification'+'\n')
for element in inserts.readlines():
	element=element.strip()
	name_a=str(element.split('\t')[0])
	len_a=float(element.split('\t')[1])
	if len_a > -50 and len_a < 50:
		output.write(name_a+'\t'+str(len_a)+'\t'+'N/A'+'\t'+'N/A'+'\t'+'New insertion'+'\n')
	else:
		for item in loop:
			try:
				item=item.strip()
				name_b=str(item.split('\t')[0])
				len_b=float(item.split('\t')[1])
				if name_a == name_b:
					ratio=len_a/len_b
					if ratio < 0.66 and ratio > 0:
						output.write(name_a+'\t'+str(len_a)+'\t'+str(len_b)+'\t'+str(ratio)+'\t'+'Deletion'+'\n')
					elif ratio > 0.66 and ratio < 1.33:
						output.write(name_a+'\t'+str(len_a)+'\t'+str(len_b)+'\t'+str(ratio)+'\t'+'Ortholog'+'\n')	
					elif ratio > 1.33:
						output.write(name_a+'\t'+str(len_a)+'\t'+str(len_b)+'\t'+str(ratio)+'\t'+'Putative nested'+'\n')
					else:
						output.write(name_a+'\t'+str(len_a)+'\t'+str(len_b)+'\t'+str(ratio)+'\t'+'Other'+'\n')
			except:
				pass
