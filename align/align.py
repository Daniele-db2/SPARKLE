#!/usr/bin/env python3

#python 3 standard library

import os
import sys
from datetime import datetime
from collections import namedtuple
import pickle
import gzip

#spark

from pyspark.shell import sqlContext,spark
from pyspark import SparkContext
from pyspark.sql import Row
from pyspark.sql import functions as F
from pyspark.sql.window import Window
from pyspark.sql.types import StructType,StructField, StringType, IntegerType
from pyspark.sql.functions import *

#numpy and tabulate

import numpy as np
import tabulate as tb


def FastqREAD(sc,FQ):

	'''
	Parse FASTQ and create SPARK dataframe from this
	'''

	#parse FASTQ and store in table

	sc.setLogLevel("WARN")
	file = sc.textFile(FQ)
	ls_ = file.take(file.count())
	ks=namedtuple('SEQUENCE', ['SEQ'])
	lines=[]

	for v in ls_:

		lines.append(v.rstrip())
		
		if len(lines) == 4:

			meanqual=np.mean(np.array([ord(x) - 33 for x in lines[3]]))
			#print(meanqual)

			if meanqual < 7:

				lines = []
				continue

			else:

				df=ks(SEQ = lines[1])

			lines=[]
			yield df

	#rdd = sc.parallelize(DFs)
	#seqDF = rdd.map(lambda x: Row(ID=x[0], SEQ=x[1], OP=x[2], QUAL=x[3]))
	#schemaSeqDF = sqlContext.createDataFrame(seqDF)
	#schemaSeqDF.show()


def hash_djb2(s):

	'''
	Convert string to int key
	'''
	
	hash = 5381
	
	for x in s:
		
		hash = (( hash << 5) + hash) + ord(x)
	
	return hash & 0xFFFFFFFF



def best_choice(SEQ,reDF,genome,sc):

	'''
	Retain only hits with lower edit distance with respect to the reference
	'''


	PG = [x["POS_GEN"] for x in reDF.rdd.collect()]
	SA = [x["POS_SEQ"] for x in reDF.rdd.collect()]


	schema = StructType([
		StructField('segmentSeq', StringType(), True),
		StructField('segmentGen', StringType(), True),
		StructField('posSeq', IntegerType(), True),
		StructField('posGen', IntegerType(), True)
	])


	df = spark.createDataFrame(spark.sparkContext.emptyRDD(), schema)
	
	for z in range(len(PG)):

		for pos_gen in PG[z]:

			seq = [(SEQ, genome[pos_gen - SA[z]: pos_gen - SA[z] + len(SEQ)], SA[z], pos_gen)]
			rddSeq = sc.parallelize(seq)
			schemaSeqDF = rddSeq.map(lambda x: Row(SEQ=x[0], GEN=x[1], POS_SEQ=x[2], POS_GEN=x[3]))
			seq_compare = sqlContext.createDataFrame(schemaSeqDF)
			df = df.union(seq_compare)
	
	df = df.withColumn("dist", F.levenshtein(F.col("segmentSeq"), F.col("segmentGen")))
	val = (1 / float(len(SEQ))) * 100
	df = df.withColumn("percentage", val*F.col( "dist")).drop("dist")
	minDF = df.agg(min(col("percentage")).alias("percentage"))
	min_percentage = [x["percentage"] for x in minDF.rdd.collect()][0]
	df = df.filter(df.percentage == min_percentage)

	return df, min_percentage


def make_matrix(sizex, sizey):

	'''
	Initialize matrix
	'''

	return np.zeros((sizex,sizey))



def createB(word_1, word_2):

	'''
	Create M matrix to store operations
	'''

	n = len(word_1) + 1
	m = len(word_2) + 1
	D = np.zeros(shape=(n, m), dtype=np.int)
	D[:,0] = range(n)
	D[0,:] = range(m)
	B = np.zeros(shape=(n, m), dtype=[("del", 'b'), ("mis", 'b'), ("ins", 'b')])
	B[1:,0] = (1, 0, 0)
	B[0,1:] = (0, 0, 1)
	for i, l_1 in enumerate(word_1, start=1):
		for j, l_2 in enumerate(word_2, start=1):
			deletion = D[i-1,j] + 1
			insertion = D[i, j-1] + 1
			mismatch = D[i-1,j-1] + (0 if l_1==l_2 else 1)
			mo = np.min([deletion, insertion, mismatch])
			B[i,j] = (deletion==mo, mismatch==mo, insertion==mo)
			D[i,j] = mo
	return B


def backtrace(A, optloc,B):

	'''
	Backtracing to retrieve the optimal path in the alignment
	'''

	if optloc == None:
		
		i, j = B.shape[0] - 1, B.shape[1] - 1
		backtrace_idxs = [(i, j)]

		while A[i][j] != 0:

			if B[i, j][1]:

				i, j = i - 1, j - 1

			elif B[i, j][0]:

				i, j = i - 1, j

			elif B[i, j][2]:

				i, j = i, j - 1

			backtrace_idxs.append((i, j))
	
	else:

		i, j = optloc
		backtrace_idxs = [(i, j)]

		while A[i][j] != 0:

			if B[i,j][1]:

				i, j = i - 1, j - 1

			elif B[i, j][0]:

				i, j = i - 1, j

			elif B[i, j][2]:

				i, j = i, j - 1

			backtrace_idxs.append((i, j))

	return backtrace_idxs




def local_align(x, y, match, mismatch, gapopen):

	'''
	Smith-Waterman alignment
	'''

	A = make_matrix(len(x)+1, len(y)+1)
	B = createB(x,y)
	best = 0
	optloc = (0, 0)

	for i in range(1, len(x) + 1):

		for j in range(1, len(y) + 1):

			if x[i - 1] == y[j - 1]:

				score=match

			else:

				score=mismatch

			#print(A[i][j - 1]+ gapopen)
			#print(A[i - 1][j] + gapopen)
			#print( A[i - 1][j - 1] + score)

			M=sorted([A[i][j - 1] + gapopen, A[i - 1][j] + gapopen, A[i - 1][j - 1] + score, 0])[-1]

			A[i][j] = M

			if M >= best:

				best = M
				optloc = (i, j)

	backtrace_idxs=backtrace(A,optloc,B)

	return backtrace_idxs


def global_align(x, y, match, mismatch, gapopen,gapextend):

	'''
	Needleman-Wunsch alignment
	'''
	Infinity = float('inf')

	M = make_matrix(len(x)+1, len(y)+1)
	B = createB(x,y)
	
	for i in range(1, len(x) + 1):
		
		M[i][0] = -Infinity
	
	for i in range(1, len(y) + 1):
		
		M[0][i] = -Infinity
	
	for i in range(1, len(x) + 1):
		
		for j in range(1, len(y) + 1):

			if x[i - 1] == y[j - 1]:

				score=match

			else:

				score=mismatch
			
			M[i][j] = M[i - 1][j - 1] + score

	backtrace_idxs=backtrace(M,None,B)

	return backtrace_idxs

def align(word_1, word_2, bt):

	'''
	Visually align strings
	'''
	
	aligned_word_1 = []
	aligned_word_2 = []
	operations = []
	line = []
	backtrace = bt[::-1]
	
	for k in range(len(backtrace) - 1):
	
		i_0, j_0 = backtrace[k]
		i_1, j_1 = backtrace[k+1]
		w_1_letter = None
		w_2_letter = None
		op = None
	
		if i_1 > i_0 and j_1 > j_0:  # either substitution or no-op
	
			if word_1[i_0] == word_2[j_0]:  # no-op, same symbol
	
				w_1_letter = word_1[i_0]
				w_2_letter = word_2[j_0]
				op = "M"
	
			else:  # cost increased: substitution
	
				w_1_letter = word_1[i_0]
				w_2_letter = word_2[j_0]
				op = "X"
	
		elif i_0 == i_1:  # insertion
	
			w_1_letter = "-"
			w_2_letter = word_2[j_0]
			op = "I"
	
		else: #  j_0 == j_1,  deletion
			w_1_letter = word_1[i_0]
			w_2_letter = "-"
			op = "D"
	
		aligned_word_1.append(w_1_letter)
		aligned_word_2.append(w_2_letter)
		operations.append(op)
		line.append("|")
	
	operations = (''.join(operations))
	
	return aligned_word_1, aligned_word_2, operations, line



def Align(REF,BIN,FQ,OUT,match,mismatch,gapopen,gapextend,kmer,alignment_type,distance):

	'''
	Proof of concept: aligner
	'''

	result=[]
	sc=SparkContext.getOrCreate()

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + ']' + '[Message] Loading index file')

	binin=open(BIN,'rb')
	data=pickle.load(binin)
	binin.close()

	#convert to df for SPARK

	rdd = sc.parallelize(data.items())
	schemaHashDF = rdd.map(lambda x: Row(ID_GEN = x[0], POS_GEN = x[1]))
	hashDF = sqlContext.createDataFrame(schemaHashDF)
	hashDF.show()


	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + ']' + '[Message] Parsing reference FASTA file')

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
	print('[' + now + ']' + '[Message] Parsing FASTQ file')

	for i,el in enumerate(FastqREAD(sc,FQ)):
		
		seq=el.SEQ

		#create k-mers for each seq and compare to reference k-mers

		word = [(i,hash_djb2(seq[j:j+kmer]),j) for j in range(0, len(seq)-kmer)]
		rddW = sc.parallelize(word)
		schemaWordDF = rddW.map(lambda x: Row(NUM_SEQ=x[0],ID_SEQ=x[1], POS_SEQ=x[2]))
		df = sqlContext.createDataFrame(schemaWordDF)
		reDF = df.join(hashDF, df.ID_SEQ == hashDF.ID_GEN, how='inner')
		reDF = reDF.orderBy(reDF.POS_SEQ).select(reDF.NUM_SEQ, reDF.ID_SEQ, reDF.POS_SEQ, reDF.POS_GEN )
		my_window = Window.partitionBy(reDF.NUM_SEQ).orderBy(reDF.POS_SEQ)
		reDF = reDF.withColumn("prev_value", F.lag(reDF.POS_SEQ).over(my_window))
		reDF = reDF.withColumn("dist", F.when(F.isnull(reDF.POS_SEQ - reDF.prev_value), 0).otherwise(reDF.POS_SEQ - reDF.prev_value))
		reDF = reDF.select(reDF.NUM_SEQ, reDF.ID_SEQ, reDF.POS_SEQ, reDF.dist, reDF.POS_GEN)
		reDF = reDF.withColumn("dist0", F.lead(reDF.dist).over(my_window))
		elDF = reDF.filter(((reDF.dist == 0) | (reDF.dist >= distance)) & ((reDF.dist0.isNull()) | (reDF.dist0 >= distance)))
		reDF = reDF.subtract(elDF)
		reDF = reDF.orderBy(reDF.POS_SEQ).select(reDF.ID_SEQ, reDF.POS_SEQ, reDF.POS_GEN)
		
		if reDF.count() >= 3: #force to try alignment only if >=3 anchors. Cannot be changed hy the user for the time being

			df,min_percentage=best_choice(seq,reDF,genome,sc)
			s1=[x["segmentSeq"] for x in df.rdd.collect()][0]
			s2=[x["segmentGen"] for x in df.rdd.collect()][0]

			if alignment_type == 'guess':

				if 100-min_percentage >= 60.0:

					t=global_align(s1,s2, match, mismatch, gapopen,gapextend)

				else:

					t=local_align(s1,s2, match, mismatch, gapopen)

			elif alignment_type == 'global':

				t=global_align(s1,s2, match, mismatch, gapopen,gapextend)

			else:

				t=local_align(s1,s2, match, mismatch, gapopen)

			#print ("Allineamento sequenza n°: ", i+1)
			as2,as1,ops,line =align(s2,s1,t)
			#print("Lunghezza sequenze: ", len(s1), "| Numero operazioni: ", len(ops))
			alignment_table = [as2, line, ops, line, as1]
			result.append(tb.tabulate(alignment_table, tablefmt="orgtbl"))
			
			if len(result) >=3:

				return result

			else:

				continue


	#print(REF,BIN,FQ,match,mismatch,gapopen,gapextend,kmer,alignment_type,OUT)

def run(parser,args):

	'''
	Execute the code and and alignment to compressed text file
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

	if not os.path.isfile(BIN):

		print('[' + now + ']' + '[Error] Missing required index for the reference genome')
		sys.exit(1)

	FQ=os.path.abspath(args.reads)

	try:

		with open(FQ) as fastq:

			assert(fastq.readline().startswith('@'))

	except:

		print('[' + now + ']' + '[Error] Invalid reads FASTQ file')
		sys.exit(1)


	OUT=os.path.abspath(args.output)

	if not os.access(os.path.dirname(OUT),os.W_OK):

		print('[' + now + ']' + '[Error] Missing write permissions on the output folder')
		sys.exit(1)



	match=args.match
	mismatch=args.mismatch
	gapopen=args.gapopen
	gapextend=args.gapextend
	kmer=args.kmer
	alignment_type=args.alignment
	distance=args.distance

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + ']' + '[Message] Mapping')

	result=Align(REF,BIN,FQ,OUT,match,mismatch,gapopen,gapextend,kmer,alignment_type,distance)

	if not OUT.endswith('.gz'):

		OUTG = OUT + '.gz'

	with gzip.open(OUTG, 'wt') as fout:

		for el in result:

			fout.write(el + '\n')

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + ']' + '[Message] Done')

	sys.exit(0)
