#!/usr/bin/env python3


#python 3 standard library

import os
import sys
from datetime import datetime
import pickle

def hash_djb2(s):

	'''
	Convert string to int key
	'''
	
	hash = 5381
	
	for x in s:
		
		hash = (( hash << 5) + hash) + ord(x)
	
	return hash & 0xFFFFFFFF



def Index(REF,BIN,kmer):

	'''
	Indexing part of the reference genome
	'''

	hashtab=dict()

	#read

	genome = ''
	startcounter=0

	with open(REF, 'r') as f:
	
		for line in f:
		
			if line[0] == '>':

				startcounter+=1

				if startcounter > 1:


					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + ']' + '[Warning] SPARKLE can process just one chromosome for the time being')

					break

				else:

					continue

			else:

				genome+=line.rstrip()


	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + ']' + '[Warning] Subsetting chromosome to first 300000 nucleotides')

	for i in range(300000):

		k_=genome[i:i+kmer]

		if 'N' in k_:

			continue

		else:

			k_int=hash_djb2(k_)

			if k_int not in hashtab:

				hashtab[k_int] = []
				hashtab[k_int].append(i)

			else:

				hashtab[k_int].append(i)


	binout=open(BIN, 'wb')
	data=pickle.dumps(hashtab)
	binout.write(data)
	binout.close()


def run(parser,args):

	'''
	Execute the code and and dump binary output to file
	'''

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')

	REF=os.path.abspath(args.genome)

	try:

		with open(REF) as genome:

			assert(genome.readline().startswith('>'))

	except:

		print('[' + now + ']' + '[Error] Invalid reference FASTA file')
		sys.exit(1)

	BIN=os.path.abspath(REF + '.bin')

	if not os.access(os.path.dirname(REF),os.W_OK):

		print('[' + now + ']' + '[Error] Missing write permissions on the reference folder')
		sys.exit(1)

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + ']' + '[Message] Indexing')

	Index(REF,BIN,args.kmer)

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + ']' + '[Message] Done')

	sys.exit(0)