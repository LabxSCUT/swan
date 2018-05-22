#!/usr/bin/env python
import argparse, csv, gzip, sys
def main():

	parser = argparse.ArgumentParser()
	parser.add_argument("callFile", metavar="callFile", \
						help="the input call file")
	parser.add_argument("-o", dest="outPrefix", default=None, \
						help="the output prefix, default assume rstrip(callFile,'.txt')")
	parser.add_argument("-f", "--format", dest="inFormat", default="CN",\
						help="input file format, either CN (cnvnator), BD (breakdancer), DL (delly), PD (pindel) or LP (Lumpy)")
	parser.add_argument("-l", "--low", dest="lowBound", default=5, type=int,\
						help="lower limit of SV size")
	parser.add_argument("-u", "--up", dest="upBound", default=100000, type=int,\
						help="up limit of SV size")
	parser.add_argument("-a", "--adjust", dest="adjust", default=0, type=int,\
						help="adjust for reference index shift")
	parser.add_argument("-p", "--pad", dest="pad", default=10, type=int, \
						help="add a pad to compensate for break point inaccuracy")

	arg_namespace = parser.parse_args()
	callFile = vars(arg_namespace)['callFile']
	outPrefix = vars(arg_namespace)['outPrefix']
	if outPrefix == None: #assume 4 digit .txt 
		outPrefix=callFile[:-4]
	inFormat = vars(arg_namespace)['inFormat']
	lowBound = vars(arg_namespace)['lowBound']
	upBound = vars(arg_namespace)['upBound']
	adjust = vars(arg_namespace)['adjust']
	pad = vars(arg_namespace)['pad']

	inFile = open(callFile,'rU')
	inContent = csv.reader(filter(lambda row: row[0]!='#', inFile), delimiter='\t')
	delContent = [] #deletion content
	dupContent = [] #duplication content
	insContent = [] #insertion content
	invContent = [] #inversion content
	trdContent = [] #traslocation deletion content
	triContent = [] #traslocation insertion content
	othContent = [] #other content

	if inFormat == 'CN': #cnvnator, no translocation
		for row in inContent:
			bedRow = ['']*5
			tmp1 = row[1].split(':')
			tmp2 = tmp1[1].split('-')
			bedRow[0] = (tmp1[0]).lstrip("chr") 
			bedRow[1] = str(int(tmp2[0])-1+adjust)
			bedRow[2] = str(int(tmp2[1])-1+adjust)
			#if int(bedRow[2])>=int(bedRow[1]): #ensure correct formatting
			#	tmp=bedRow[1]; bedRow[1]=bedRow[2]; bedRow[2]=tmp 
			bedRow[3] = row[0]
			bedRow[4] = row[2]
			if bedRow[3] == 'deletion':
				bedRow[3] = 'DEL'
				delContent.append(bedRow) 
			elif bedRow[3] == 'duplication':
				bedRow[3] = 'DUP'
				dupContent.append(bedRow) 
			elif bedRow[3] == 'insertion':
				bedRow[3] = 'INS'
				insContent.append(bedRow) 
			elif bedRow[3] == 'inversion':
				bedRow[3] = 'INV'
				invContent.append(bedRow) 
			#elif bedRow[3] == 'translocation':
			#	bedRow[3] = 'TRA'
			#	trdContent.append(bedRow) 
			#	triContent.append(bedRow) 
			else:
				othContent.append(bedRow)
	elif inFormat == 'BD': #breakdancer
		for row in inContent:
			#print row
			bedRow = ['']*5
			bedRow[0] = row[0]
			bedRow[1] = str(int(row[1])-1+adjust)
			bedRow[2] = str(int(row[4])-1+adjust)
			bedRow[3] = row[6]
			bedRow[4] = row[7]
			if bedRow[3] == 'DEL':
				if int(bedRow[2])<=int(bedRow[1]): #ensure correct formatting, there is obvious bug in BD such that CTX was reported as DEL
					continue
				delContent.append(bedRow) 
			elif bedRow[3] == 'DUP':
				dupContent.append(bedRow) 
			elif bedRow[3] == 'INS':
				if int(bedRow[2])<=int(bedRow[1]): #ensure correct formatting, there is obvious bug in BD such that CTX was reported as DEL
					continue
				insContent.append(bedRow) 
			elif bedRow[3] == 'INV':
				invContent.append(bedRow) 
			elif bedRow[3] == 'ITX':
				donRow = bedRow
				donRow[3] = 'INS'; donRow[2]=str(int(bedRow[1])+1+adjust)
				insContent.append(donRow) 
				recRow = bedRow
				recRow[3] = 'INS'; recRow[1]=str(int(bedRow[2])-1+adjust)
				insContent.append(recRow) 
			elif bedRow[3] == 'CTX': #insert two entries
				donRow = bedRow
				donRow[3] = 'INS'; donRow[2]=str(int(bedRow[1])+1+adjust)
				insContent.append(donRow) 
				recRow = bedRow
				recRow[3] = 'INS'; recRow[1]=str(int(bedRow[2])-1+adjust); recRow[0]=row[3]
				insContent.append(recRow) 
			else:
				othContent.append(bedRow)
	elif inFormat == 'PD': #pindel
		for row in inContent:
			#print row
			tmp1 = row[1].split(' ')
			tmp2 = row[4].split(' ')
			size = abs(int(float(tmp1[1])))
			if size>upBound or size<lowBound:
				continue
			bedRow = ['']*5
			bedRow[0] = (row[3].split(' '))[1]
			bedRow[1] = str(int(tmp2[1])+adjust)
			bedRow[2] = str(int(row[5])+adjust)
			bedRow[3] = tmp1[0]
			bedRow[4] = tmp1[1] 
			if bedRow[3] == 'D':
				bedRow[3] = 'DEL'
				delContent.append(bedRow) 
			elif bedRow[3] == 'I':
				bedRow[3] = 'INS'
				insContent.append(bedRow) 
			elif bedRow[3] == 'INV':
				bedRow[3] = 'INV'
				invContent.append(bedRow) 
			elif bedRow[3] == 'TD': #tandom duplication
				bedRow[3] = 'DUP'
				dupContent.append(bedRow) 
			else:
				othContent.append(bedRow)
	elif inFormat == 'DL': #delly translocation and others have different format
		for row in inContent:
			#print row
			bedRow = ['']*5
			bedRow[0] = row[0]
			bedRow[1] = str(int(row[1])-1+adjust)
			bedRow[3] = ((row[6].split('_'))[0]).lstrip('>')
			bedRow[4] = row[3]
			if bedRow[3] == 'Translocation': #3rd col is the chr name
				bedRow[2] = str(int(row[3])-1+adjust)
				donRow = bedRow
				donRow[3] = 'INS'; donRow[2]=str(int(bedRow[1])+1)
				insContent.append(donRow) 
				recRow = bedRow
				recRow[3] = 'INS'; recRow[1]=str(int(bedRow[2])-1); recRow[0]=row[2]
				insContent.append(recRow) 
			elif bedRow[3] == 'Deletion': #3rd col is the end pos
				bedRow[2] = str(int(row[2])-1+adjust)
				bedRow[3] = 'DEL'
				delContent.append(bedRow) 
			elif bedRow[3] == 'Duplication':
				bedRow[2] = str(int(row[2])-1+adjust)
				bedRow[3] = 'DUP'
				dupContent.append(bedRow) 
			elif bedRow[3] == 'Insertion':
				bedRow[2] = str(int(row[2])-1+adjust)
				bedRow[3] = 'INS'
				insContent.append(bedRow) 
			elif bedRow[3] == 'Inversion':
				bedRow[2] = str(int(row[2])-1+adjust)
				bedRow[3] = 'INV'
				invContent.append(bedRow) 
			else:
				othContent.append(bedRow)
	elif inFormat == 'LP':
		for row in inContent:
			bedRow = ['']*5
			bedRow[0] = row[0]
			bedRow[1] = str(int(row[1])-1+adjust) #start
			bedRow[3] = row[10].split(':')[1]  #type
			bedRow[4] = "0" #size
			if bedRow[3] == 'INTERCHROM': #3rd col is the chr name
				donRow = bedRow
				donRow[3] = 'INS'; donRow[2]=str(int(row[2])+1+adjust) #end
				insContent.append(donRow) 
				recRow = bedRow
				recRow[3] = 'INS'; recRow[1]=str(int(row[4])+adjust); recRow[2]=str(int(row[4])+1+adjust); recRow[0]=row[3]
				insContent.append(recRow) 
			elif bedRow[3] == 'DELETION': #3rd col is the end pos
				bedRow[2] = str(int(row[4])-1+adjust)
				bedRow[3] = 'DEL'
				delContent.append(bedRow) 
			elif bedRow[3] == 'DUPLICATION':
				bedRow[2] = str(int(row[4])-1+adjust)
				bedRow[3] = 'DUP'
				dupContent.append(bedRow) 
			elif bedRow[3] == 'INSERTION':
				bedRow[2] = str(int(row[4])-1+adjust)
				bedRow[3] = 'INS'
				insContent.append(bedRow) 
			elif bedRow[3] == 'INVERSION':
				bedRow[2] = str(int(row[4])-1+adjust)
				bedRow[3] = 'INV'
				invContent.append(bedRow) 
			else:
				othContent.append(bedRow)
	elif inFormat == 'SW':
		for row in inContent:
			bedRow = ['']*5 #chr, start, end, type, size
			bedRow[0] = row[0] #chr
			bedRow[1] = str(int(row[1])-1+adjust) #start
			bedRow[2] = str(int(bedRow[1]) + 1) #end
			bedRow[3] = row[5] #type
			bedRow[4] = "0" # default size for unprocessd
			if bedRow[3] == 'DEL': #deletion
				bedRow[4] = str((row[8].split(';')[4]).split('=')[-1]) #size
				bedRow[2] = str(int(bedRow[1])+int(bedRow[4])+adjust) #end
				delContent.append(bedRow) 
			elif bedRow[3] == 'INV':
				bedRow[4] = str((row[8].split(';')[4]).split('=')[-1]) #size
				bedRow[2] = str(int(bedRow[1])+int(bedRow[4])+adjust) #end
				invContent.append(bedRow)
			elif bedRow[3] == 'DUP':
                                bedRow[2] = str(int(bedRow[2]) + adjust)
				bedRow[4] = str(int(bedRow[2]) - int(bedRow[1]) + 1)
				dupContent.append(bedRow)
                        elif bedRow[3] == 'INS':
				bedRow[1] = str(int(bedRow[1]) - pad)
				bedRow[2] = str(int(bedRow[1]) + pad)
				bedRow[4] = str(2*pad)
				insContent.append(bedRow) 
			else:
				othContent.append(bedRow)
	elif inFormat == 'LP-VCF':
		for row in inContent:
			bedRow = ['']*5 #chr, start, end, type, size
			bedRow[0] = row[0] #chr
			bedRow[1] = str(int(row[1])-1+adjust) #start
			bedRow[2] = str(int(bedRow[1]) + 1) #end
			bedRow[3] = (row[7].split(';')[0]).split('=')[-1]  #type
			bedRow[4] = "0" # default size for unprocessd
			if bedRow[3] == 'DEL': #deletion
				bedRow[4] = str(-int((row[7].split(';')[2]).split('=')[-1])) #size
				bedRow[2] = str(int(bedRow[1])+int(bedRow[4])+adjust) #end
				delContent.append(bedRow) 
			elif bedRow[3] == 'INV':
				bedRow[4] = (row[7].split(';')[2]).split('=')[-1] #size
				bedRow[2] = str(int(bedRow[1])+int(bedRow[4])+adjust) #end
				invContent.append(bedRow)
			elif bedRow[3] == 'BND':
				bedRow[1] = str(int(bedRow[1]) - pad)
				bedRow[2] = str(int(bedRow[1]) + pad)
				bedRow[4] = str(2*pad)
				bedRow[3] = 'UNK'
				othContent.append(bedRow) 
			else:
				othContent.append(bedRow)
	elif inFormat == 'MT-VCF': #manta VCF
		for row in inContent:
			bedRow = ['']*5 #chr, start, end, type, size
			#print(row[7])
			info = dict((k, v) for k,v in ((item.split('=') + [None]*99)[:2] for item in row[7].split(';'))) #hack for unknown count splitting
			#info = dict((k.strip(), v.strip()) for k,v in 
      #      (item.split('=') for item in row[7].split(';')))
			bedRow[0] = row[0] #chr
			bedRow[1] = str(int(row[1])-1+adjust) #start
			bedRow[2] = str(int(bedRow[1]) + 1) #end
			bedRow[3] = info['SVTYPE']
			bedRow[4] = "0" # default size for unprocessd
			if bedRow[3] == 'DEL': #deletion
				bedRow[4] = str(-int(info['SVLEN'])) #size
				bedRow[2] = str(int(bedRow[1])+int(bedRow[4])+adjust) #end
				delContent.append(bedRow) 
			elif bedRow[3] == 'INV':
				bedRow[4] = info['SVLEN']
				bedRow[2] = str(int(bedRow[1])+int(bedRow[4])+adjust) #end
				invContent.append(bedRow)
			elif bedRow[3] == 'INS':
				bedRow[1] = str(int(bedRow[1]) - pad)
				bedRow[2] = str(int(bedRow[1]) + pad)
				bedRow[4] = str(2*pad)
				insContent.append(bedRow)
			elif bedRow[3] == 'DUP':
				bedRow[4] = info['SVLEN']
				bedRow[2] = str(int(bedRow[1])+int(bedRow[4])+adjust) #end
				dupContent.append(bedRow)
			elif bedRow[3] == 'BND':
				bedRow[1] = str(int(bedRow[1]) - pad)
				bedRow[2] = str(int(bedRow[1]) + pad)
				bedRow[4] = str(2*pad)
				bedRow[3] = 'UNK'
				othContent.append(bedRow) 
			else:
				othContent.append(bedRow)
	else:
		print >>sys.stderr, "format not supported, quit!"
		pass

	print delContent
	delOut = gzip.open(outPrefix+".del.bed.gz",'wb') #can only use gzip file here
	if delContent != []:
		print >>sys.stderr, "has content"
		print >>delOut,	"\n".join(["\t".join(v) for v in delContent])
	delOut.close()
	print >>sys.stderr, "covert finished"

	print dupContent
	dupOut = gzip.open(outPrefix+".dup.bed.gz",'wb') #can only use gzip file here
	if dupContent != []:
		print >>sys.stderr, "has content"
		print >>dupOut,	"\n".join(["\t".join(v) for v in dupContent])
	dupOut.close()
	print >>sys.stderr, "covert finished"

	print insContent
	insOut = gzip.open(outPrefix+".ins.bed.gz",'wb') #can only use gzip file here
	if insContent != []:
		print >>sys.stderr, "has content"
		print >>insOut,	"\n".join(["\t".join(v) for v in insContent])
	insOut.close()
	print >>sys.stderr, "covert finished"

	print invContent
	invOut = gzip.open(outPrefix+".inv.bed.gz",'wb') #can only use gzip file here
	if invContent != []:
		print >>sys.stderr, "has content"
		print >>invOut,	"\n".join(["\t".join(v) for v in invContent])
	invOut.close()
	print >>sys.stderr, "covert finished"

	print othContent
	othOut = gzip.open(outPrefix+".oth.bed.gz",'wb') #can only use gzip file here
	if othContent != []:
		print >>sys.stderr, "has content"
		print >>othOut,	"\n".join(["\t".join(v) for v in othContent])
	othOut.close()
	print >>sys.stderr, "covert finished"

	allContent = delContent + dupContent + insContent + invContent + othContent
	print allContent
	allOut = gzip.open(outPrefix+".all.bed.gz",'wb') #can only use gzip file here
	if allContent != []:
		print >>sys.stderr, "has content"
		print >>allOut,	"\n".join(["\t".join(v) for v in allContent])
	allOut.close()
	print >>sys.stderr, "covert finished"

	print >>sys.stderr, "written results in", outPrefix, ".{del,dup,ins,inv,oth}.bed.gz"

if __name__=="__main__":
	main()
