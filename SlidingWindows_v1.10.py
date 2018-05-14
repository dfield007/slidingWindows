#!/usr/bin/env python
# Essential modules
import re, sys, getopt, itertools, os, math
import numpy as np
from collections import OrderedDict
from collections import defaultdict
from operator import truediv
#cwd=os.getcwd()
print "                                  "
print "----------------------------------------------"
print "-    Genome scan sliding window analyses     -"
print "-  Diversity and differentiation statistics  -"
print "-   For whole genome data (poolSeq)          -"
print "-    requires *.sync data file               -"
print "-      David L. Field 09/05/2018             -"
print "-      09/05/2018 (last update)              -"
print "-      david.field@univie.ac.at              -"
print "----------------------------------------------"
print "                                  "
#  to run, you need the following arguments: 
#  (1) list of scaffolds to analyse, 
#  (2) file with population names and sample size, 
#  (3) output file name,
#  (4) the min depth to include a site, 
#  (5) the max depth to include a site, 
#  (6) min copies for an allele call (i.e. > 1, no singletons)
#  (7) min number of populations (pools) with min depth to include a site, 
#  (8) window size (bp), 
#  (9) overlap between windows (bp),
#  (10) keep site specific file (1 = Y,0 = N), 
#  (11) window analysis only (1 = Y,0 = N),
#  (12) site analysis only (1 = Y, 0 = N) 
#  (13) keep all positions (1 = Y,0 = N),
#
#  min depth: if running on individual based data, make minimum depth the size of the smallest number of haploid genomes (i.e. 4x2 = 8),
#  min copies: if running on individual based data, make minimum copies = 1 (i.e. allow for singletons),
#  input: scaffold input list must have line breaks saved as Unix (LF). Use text wrangler or similar text editor,
#  window analysis only: only works if a site specific file(s) already created and present in the same folder,
#  min scaffold size: currently skips over any scaffold < 15kB. Ensures at least two windows minimum.
# 
# Example 1, Individual based WGS data (running site specific and window analyses):"
#  SlidingWindows_C.py scaffoldList_RosElcombined.txt popDetails.txt output_windows_ind_ROSassembly_ABC_run1.txt 8 8 8 1 8 10000 5000 1 0"
# Example 2, pooled data (Cython version) - site analysis and windows (latest version August 2016) :"
#   SlidingWindows_C.py scaffoldList_rosAssembly_new.txt popDetails_10pool.txt windowOutput_RosAssembly_test.txt 15 200 2 2 10000 5000 1 0"
# Example 3, pooled data (Python version) - site analysis and windows (latest version August 2016) :"
#   python SlidingWindowsEditsNew.py scaffoldList_rosAssembly_new.txt popDetails_10pool.txt windowOutput_RosAssembly_test.txt 15 200 2 2 10000 5000 1 0 0"
# Example 4, pooled data - window analysis ONLY (latest version August 2016) :"
#   python SlidingWindowsEditsNew.py scaffoldList_rosAssembly_new.txt popDetails_10pool.txt windowOutput_RosAssembly_test.txt 15 200 2 2 10000 5000 1 1 0"
#
print "Importing data... "
print 'Input file (scaffold list):', sys.argv[1]
print 'Input file (pop details):', sys.argv[2]
print 'Output file (window analyses):', sys.argv[3]
print 'Min depth:', sys.argv[4]
print 'Max depth:', sys.argv[5]
print 'Min copies for an allele call:', sys.argv[6]
print 'Min population number with min depth to include:', sys.argv[7]
print 'Window size (bp):', sys.argv[8]
print 'Window overlap (bp):', sys.argv[9]
print 'Keep site specific file:', sys.argv[10]
print 'Window analysis only:', sys.argv[11]
print 'Site analysis only:', sys.argv[12]
print 'Keep all sites temp output:', sys.argv[13]

minDepthSeqInfo = int(sys.argv[4])
maxDepthSeqInfo = int(sys.argv[5])
minAlleleCount = int(sys.argv[6])
minPoolNum = int(sys.argv[7])
windowSize = int(sys.argv[8])
overlapSize = int(sys.argv[9])
keepTemp = int(sys.argv[10])
windowOnly = int(sys.argv[11])
sitesOnly = int(sys.argv[12])
allSites = int(sys.argv[13])

# Gen stat functions
def adjustData(LineList2, Pool, minDepthSeqInfo, maxDepthSeqInfo, minAlleleCount):
	# David Field 12.08.16 
	LineData = {}
	thisPoolVector=(Pool+2)
	alleles=str(LineList2[thisPoolVector])
	alleles=alleles.split(":")
	alleles=alleles[0:4]
	alleles=map(float,alleles)
	LineData['pop']= Pool
	LineData['count']= alleles
	LineData['depth']= sum(LineData['count'])
	if (LineData['depth']>0):
		LineData['alleleFreq']=[x / LineData['depth'] for x in LineData['count']]
	if (LineData['depth']==0):
		LineData['alleleFreq']=[0,0,0,0]
	seqCode=["A","T","C","G"]			
	counter = 0
	for val in LineData['count']:
		LineData['count'][counter]=int(val)
		counter=counter+1
	LineData['base_present'] = [1 if x>=minAlleleCount else 0 for x in LineData['count']]
	counter=0
	LineData['count_adj']=LineData['count']
	LineData['bases_s']=['','','','']
	for allele in LineData['base_present']:
		LineData['bases_s'][counter]=int(allele)*seqCode[counter]
		counter=counter+1
	counter=0	
	LineData['bases_t']=''.join(LineData['bases_s'])
	LineData['count_adj'] = [q*w for q,w in zip(LineData['count_adj'],LineData['base_present'])]
	LineData['depth_adj']=float(sum(LineData['count_adj']))
	if (LineData['depth_adj']>0):
		LineData['alleleFreq_adj']=[x / LineData['depth_adj'] for x in LineData['count_adj']]
	if (LineData['depth_adj']==0):
		LineData['alleleFreq_adj']=[0,0,0,0]
	LineData['minDepth_pass'] = 0
	if (LineData['depth_adj']>=minDepthSeqInfo and LineData['depth_adj']<= maxDepthSeqInfo):
		LineData['minDepth_pass'] = 1
	LineData['poly'] = 0
	if (len(LineData['bases_t'])>1):
		LineData['poly']=1
	alleles=[]
	return LineData
def divergeStats(Set1, Set2, Set1Div, Set2Div, p_1, q_1, p_2, q_2, minC):
	# David Field 12.08.16 
	Pbr,PbrA,Pt,PtA,Pb_raw,Pb_fA,F_Pb,F_Pbr,F_PbrA,BP,counter,counter2=[0,0,0,0,0,0,0,0,0,0,0,1]
	Px,Qy,Py,Qx=[0,0,0,0]
	if (Set1[0]['minDepth_pass']==1 and Set2[0]['minDepth_pass']==1):
		BP = 1
		Pb_raw = (p_1*q_2)+(q_1*p_2)
		PbrA = ((Set1Div[0]['Pi_adj']+Set2Div[0]['Pi_adj'])/2)
		Pbr = ((Set1Div[0]['Pi']+Set2Div[0]['Pi'])/2)
		Pmean = (p_1+p_2)/2
		Qmean = (q_1+q_2)/2
		Pt = 2*(Pmean*Qmean)
		t_C_Sets = [q+w for q,w in zip(Set1[0]['count_adj'],Set2[0]['count_adj'])]
		Dp_t_C_Sets = sum(t_C_Sets)
		Dp_t_C_Sets = float(Dp_t_C_Sets)
		AVal_Total = ((Dp_t_C_Sets)-1)/((Dp_t_C_Sets)-(2*minC)+1)
		# or 
		dp1 = Set1[0]['depth_adj']
		dp2 = Set2[0]['depth_adj']
		dp_m = (dp1+dp2)/2
		Aval_TotalMean  = ((dp_m)-1)/((dp_m)-(2*minC)+1)
		PtA = Pt*Aval_TotalMean	
		num2 = (Pt-Pbr)
		den2 = (Pt)
		F_Pbr=0
		if (den2>0):
			F_Pbr = num2/den2		
		num3 = (PtA-PbrA)
		den3 = (PtA)
		F_PbrA=0
		if (den3>0):
			F_PbrA = num3/den3
		num = (Pb_raw-Pbr)
		den = (Pb_raw+Pbr)
		F_Pb = 0
		if (den>0):
			F_Pb = num/den
		num4 = (1+F_PbrA)
		den4 = (1-F_PbrA)
		if (den4>0):
			Pb_fA = PbrA*(num4/den4)
	PairwiseDivergence = {}
	PairwiseDivergence['Pi_bar'] = Pbr
	PairwiseDivergence['Pi_bar_adj'] = PbrA
	PairwiseDivergence['Pi_T'] = Pt
	PairwiseDivergence['Pi_T_adj'] = PtA
	PairwiseDivergence['d_xy_raw'] = Pb_raw
	PairwiseDivergence['d_xy_fromFstAdj'] = Pb_fA
	PairwiseDivergence['Fst_fromPi'] = F_Pbr
	PairwiseDivergence['Fst_fromPiAdj'] = F_PbrA
	PairwiseDivergence['Fst_fromDxy'] = F_Pb
	PairwiseDivergence['BothPassed'] = BP
	return PairwiseDivergence
def diverseStats(poolData,minAlleleCount):
	# David Field 12.08.16 
	P_diversity = {}
	Pi_raw = 0; Pi_adj = 0; AdjVal = 0
	if (poolData['minDepth_pass']==1):
		if (poolData['poly']==1):
			Vals = [x for x in poolData['alleleFreq_adj'] if x > 0]
			Pi_raw = 2*Vals[0]*Vals[1]
			AdjVal = (float(poolData['depth_adj'])-1)/(float(poolData['depth_adj'])-(2*minAlleleCount)+1)
			Pi_adj = Pi_raw*AdjVal
	P_diversity['Pi'] = Pi_raw
	P_diversity['Pi_adj'] = Pi_adj
	P_diversity['AdjVal'] = AdjVal
	return P_diversity
# random functions
def MyDivision(num, denom):
	if denom==0:
		return "NaN"
	if (denom=='NaN' or denom=='NaN'):	
		return "NaN"
	else:
		return (num/denom)
def MySubtraction(first,second):
	if (first=='NaN' or second=='NaN'):
		return 'NaN'
	if (first!='NaN' and second!='NaN'):
		return (first-second)
def MyMultiplication(first, second):
	if second=='NaN':
		return 'NaN'
	if second!='NaN':
		return (first * second)
def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False
def round_down(num, divisor):
	# test
	return (num - (num%divisor))	
InFilePop = open(sys.argv[2], 'r')
InFilePop2 = [i.split()[0:] for i in InFilePop.readlines()]
InFilePop.close() 
numPops = len(InFilePop2)
print 'numPops:', numPops

pops = range(1, numPops+1)
comboPops_nums = []
print '\n population names: ', pops
for subset in itertools.combinations(pops, 2):
    comboPops_nums.append(subset) 
print '\n population pairs: ' , comboPops_nums
comboPops_header_sites = []
for subset in itertools.combinations(pops, 2):
    comboPops_header_sites.append(str(subset[0]) + '_' + str(subset[1]))
# header site output
header_sites_start = ["scaffold", "position", "LG", "cM", "ref", "bases_t", "allele_num", "1or2bases", "minDepth_pools", "minDepth_pass"]
headers_sites_pops = ['p_{0}'.format(i) for i in pops] + ['q_{0}'.format(i) for i in pops] + ['poly_pool_{0}'.format(i) for i in pops] \
+ ['dpthAdj_{0}'.format(i) for i in pops] + ['dpthPass_{0}'.format(i) for i in pops] + ['pi_{0}'.format(i) for i in pops] \
+ ['piAdj_{0}'.format(i) for i in pops]
headers_sites_pairs = ['dpthPolyPass_{0}'.format(i) for i in comboPops_header_sites] + ['piBar_{0}'.format(i) for i in comboPops_header_sites] \
+ ['piBarAdj_{0}'.format(i) for i in comboPops_header_sites] + ['piT_{0}'.format(i) for i in comboPops_header_sites] \
+ ['piTadj_{0}'.format(i) for i in comboPops_header_sites] + ['dXYraw_{0}'.format(i) for i in comboPops_header_sites] \
+ ['dXYfromFstadj_{0}'.format(i) for i in comboPops_header_sites] + ['FstPi_{0}'.format(i) for i in comboPops_header_sites] \
+ ['FstPiAdj_{0}'.format(i) for i in comboPops_header_sites] + ['FstfromDxy_{0}'.format(i) for i in comboPops_header_sites]
OutputStringHeaderTempFile = '\t'.join(header_sites_start+headers_sites_pops+headers_sites_pairs)

# header depth output
header_sites_dpth_start = ["scaffold", "position", "LG", "cM", "ref", "bases_t", "allele_num", "1or2bases", "minDepth_pools", "minDepth_pass"]
headers_sites_dpth_pops = ['p_{0}'.format(i) for i in pops] + ['dpthAdj_{0}'.format(i) for i in pops] + ['dpthPass_{0}'.format(i) for i in pops]
headers_sites_dpth_pairs = ['dpthPolyPass_{0}'.format(i) for i in comboPops_header_sites]
OutputStringHeaderTempFile_dpth = '\t'.join(header_sites_dpth_start+headers_sites_dpth_pops+headers_sites_dpth_pairs)

# header window output
header_window_start = ["scaffold","LG","cM","windowStart","windowEnd","coverage"]
headers_window_pops = ['poly_pool_{0}'.format(i) for i in pops] \
+ ['dpthAdj_{0}'.format(i) for i in pops] + ['dpthPass_{0}'.format(i) for i in pops] + ['pi_{0}'.format(i) for i in pops] \
+ ['piAdj_{0}'.format(i) for i in pops]
# ['p_{0}'.format(i) for i in pops] + ['q_{0}'.format(i) for i in pops] + 
headers_window_pairs = ['dpthPolyPass_{0}'.format(i) for i in comboPops_header_sites] + ['piBar_{0}'.format(i) for i in comboPops_header_sites] \
+ ['piBarAdj_{0}'.format(i) for i in comboPops_header_sites] + ['piT_{0}'.format(i) for i in comboPops_header_sites] + ['piTadj_{0}'.format(i) for i in comboPops_header_sites] \
+ ['dXYraw_{0}'.format(i) for i in comboPops_header_sites] + ['dXYfromFstadj_{0}'.format(i) for i in comboPops_header_sites] \
+ ['FstPi_{0}'.format(i) for i in comboPops_header_sites] + ['FstPiAdj_{0}'.format(i) for i in comboPops_header_sites] \
+ ['FstfromDxy_{0}'.format(i) for i in comboPops_header_sites] + ['FstfromMeanPi_{0}'.format(i) for i in comboPops_header_sites] \
+ ['FstfromMeanPiAdj_{0}'.format(i) for i in comboPops_header_sites] + ['dXYfromMeanFst_{0}'.format(i) for i in comboPops_header_sites] \
+ ['dXYfromMeanFstAdj_{0}'.format(i) for i in comboPops_header_sites] + ['dAfromMeanFst_{0}'.format(i) for i in comboPops_header_sites] \
+ ['dAfromMeanFstAdj_{0}'.format(i) for i in comboPops_header_sites]
OutputStringHeaderSlidingWindows = '\t'.join(header_window_start+headers_window_pops+headers_window_pairs)

# tracking:
headerCombined = header_sites_start + headers_sites_pops + headers_sites_pairs
headerCombined_dpth = header_sites_dpth_start + headers_sites_dpth_pops + headers_sites_dpth_pairs
# indexing statistic locations
# site output file
# 7 single pop stats
#statsOrderSingle = ['poly_pool', 'dpthAdj', 'dpthPass', 'pi', 'piAdj']
statsOrderSingle = ['p', 'q', 'poly_pool', 'dpthAdj', 'dpthPass', 'pi', 'piAdj']

# 10 pairwise stats (+ 5 later from window averages)
statsOrderPairs = ['dpthPolyPass', 'piBar', 'piBarAdj', 'piT', 'piTadj', 'dXYraw', 'dXYfromFstadj', 'FstPi', 'FstPiAdj', 'FstfromDxy']
# dpth output file (3 single pop, 1 pairwise)
statsOrderSingle_dpth = ['p','dpthAdj', 'dpthPass']; statsOrderPairs_dpth = ['dpthPolyPass']

# single pop stats
# initial vals
# initial vals
statsColumnList = defaultdict(list)
thisStart = len(header_sites_start)
thisEnd = (thisStart+numPops)

singlePopStatsCounter = 1
counter = 1
for thisStat in statsOrderSingle: 
	if counter == 1:
		thisRange = headerCombined[thisStart:thisEnd]
		statsColumnList[thisStat].append(thisStart)
		statsColumnList[thisStat].append(thisEnd)
		thisStart = thisEnd
		thisEnd = thisStart+numPops
	if counter != 1:
		thisRange = headerCombined[thisStart:thisEnd]
		statsColumnList[thisStat].append(thisStart)
		statsColumnList[thisStat].append(thisEnd)
		thisStart = thisEnd
		thisEnd = thisStart+numPops

# pop pairs
# initial vals
statsColumnListPairs = defaultdict(list)
numPairs = len(comboPops_nums)
thisStart = thisEnd-numPops
thisEnd = (thisStart+numPairs)
counter = 1
for thisStat in statsOrderPairs: 
    if counter == 1:
        thisRange = headerCombined[thisStart:thisEnd]
        statsColumnListPairs[thisStat].append(thisStart)
        statsColumnListPairs[thisStat].append(thisEnd)
        thisStart = thisEnd
        thisEnd = thisStart+numPairs
    if counter != 1:
        thisRange = headerCombined[thisStart:thisEnd]
        statsColumnListPairs[thisStat].append(thisStart)
        statsColumnListPairs[thisStat].append(thisEnd)
        thisStart = thisEnd
        thisEnd = thisStart+numPairs	

# depth file
statsColumnList_dpth = defaultdict(list)
thisStart = len(header_sites_dpth_start)
thisEnd = (thisStart+numPops)

singlePopStatsCounter = 1
counter = 1
for thisStat in statsOrderSingle_dpth: 
	if counter == 1:
		thisRange = headerCombined_dpth[thisStart:thisEnd]
		statsColumnList_dpth[thisStat].append(thisStart)
		statsColumnList_dpth[thisStat].append(thisEnd)
		thisStart = thisEnd
		thisEnd = thisStart+numPops
	if counter != 1:
		thisRange = headerCombined_dpth[thisStart:thisEnd]
		statsColumnList_dpth[thisStat].append(thisStart)
		statsColumnList_dpth[thisStat].append(thisEnd)
		thisStart = thisEnd
		thisEnd = thisStart+numPops
				
# Depth pop pairs
# initial vals
statsColumnList_Pairs_dpth = defaultdict(list)
numPairs = len(comboPops_nums)
thisStart = thisEnd-numPops
thisEnd = (thisStart+numPairs)
counter = 1
for thisStat in statsOrderPairs_dpth: 
    if counter == 1:
        thisRange = headerCombined_dpth[thisStart:thisEnd]
        statsColumnList_Pairs_dpth[thisStat].append(thisStart)
        statsColumnList_Pairs_dpth[thisStat].append(thisEnd)
        thisStart = thisEnd
        thisEnd = thisStart+numPairs
    if counter != 1:
        thisRange = headerCombined_dpth[thisStart:thisEnd]
        statsColumnList_Pairs_dpth[thisStat].append(thisStart)
        statsColumnList_Pairs_dpth[thisStat].append(thisEnd)
        thisStart = thisEnd
        thisEnd = thisStart+numPairs	

# allele frequency function 
def alleleFreq(allPools,pops_allAdjustments):
	#print "allPools:", allPools
	#print "Counts all pools:", allPools['count_allPools']
	highestCount = max(x for x in allPools['count_allPools'] if x > 0)
	lowestCount = min(x for x in allPools['count_allPools'])
	#print "\nHighestCount:", highestCount
	#print "\nLowestCount:", lowestCount
	lowestCountNonZero = min(x for x in allPools['count_allPools'] if x > 0)
	if (highestCount!=lowestCountNonZero):
		#print "\nAlternative Allele found"
		lowestCount = lowestCountNonZero
		#print "\nLowestCount Non Zero:",lowestCountNonZero
		P_allele = [1 if x==highestCount else 0 for x in allPools['count_allPools']]
		#print "\n P allele:",P_allele
		Q_allele = [1 if x==lowestCount else 0 for x in allPools['count_allPools']]
		#print "\n Q allele:",Q_allele
		alleleFreqs = defaultdict(list)
		#print '\n pops: ', pops
		for thisPop in pops:
			#print '\n this pop: ', thisPop
			#print '\n allele freq adj - this pop: ', pops_allAdjustments[thisPop][0]['alleleFreq_adj']
			p_1 = max([q*w for q,w in zip(pops_allAdjustments[thisPop][0]['alleleFreq_adj'],P_allele)])
			q_1 = max([q*w for q,w in zip(pops_allAdjustments[thisPop][0]['alleleFreq_adj'],Q_allele)])
			alleleFreqs[thisPop].append(p_1)
			alleleFreqs[thisPop].append(q_1)
			#print '\n alleleFreqs: ', alleleFreqs
	if (highestCount==lowestCountNonZero):
		alleleFreqs = defaultdict(list)
		for thisPop in pops:
			#print '\n this pop: ', thisPop
			#print '\n pops_allAdjustments: ', pops_allAdjustments
			#print '\n pops_allAdjustments [1]: ', pops_allAdjustments[thisPop]
			#print '\n allele freq adj - this pop: ', pops_allAdjustments[thisPop][0]['alleleFreq_adj']
			p_1 = 1
			q_1 = 0
			alleleFreqs[thisPop].append(p_1)
			alleleFreqs[thisPop].append(q_1)
			#print '\n alleleFreqs: ', alleleFreqs
		minAllele  = min(x for x in allPools['count_allPools'] if x > 0)
	return alleleFreqs

# read 
InFileTest = open(sys.argv[1], 'r')
InFile  = [i.split()[0:] for i in InFileTest.readlines()]
InFileTest.close() 
lengthInFile=len(InFile)
numLines = float(sum(1 for line in InFile))
# Initialise 
LineNumber = 0; thisScaff = []; thisPos = []; theDataSummary = []; theDataFull = []

# Scaffold list
for Line in InFile:
	thisScaff = [];	thiscM = []; thisLG = []; theDataTemp = []
	PercentCompleteWindows = 0
	if LineNumber > 0:
		thisScaff = str(Line[2])
		thisLG = str(Line[0])
		thiscM = Line[3]
		print '\n '
		print "\n Running analyses... Scaffold ",int(round(LineNumber,0)),"of",int(round((numLines-1),0))
		if windowOnly == 0:
			InFile2 = open(thisScaff + '.sync', 'r')
			InFile2lines  = [i.split()[0:] for i in InFile2.readlines()]
			InFile2.close() 
			InFile2lineNum = len(InFile2lines)
			totalScaffSize = InFile2lines[InFile2lineNum-2][1]
			print '\n thisScaff: ', thisScaff
			print ' totalScaffSize: ',totalScaffSize
			print ' LG', thisLG
			print ' cM', thiscM
		    # Scan positions (only) to locate non-consecutive jumps
			with open(thisScaff + '.sync') as InFile2:
				lineNum_count = 0; lineNum_consec = []; position_consec = []
				print "\n Scanning for consecutive positions..."
				for line in InFile2:
					LineList2 = line.split()[0:]
					#print '\n LineList2: ', LineList2 
					currentScaff = str(LineList2[0])
					currentPos = int(LineList2[1])
					# print ",", currentPos
					lineNum_consec.append(lineNum_count)
					position_consec.append(currentPos)
					lineNum_count = lineNum_count + 1
			# calculating jumps
			#print "\n line numbers..."
			#print lineNum_consec
			#print "\n position_consec difference..."
			diff_consec = []
			#diff_consec.append([abs(y - x) for x,y in zip(position_consec,position_consec[1:])])
			diff_consec = ([abs(y - x) for x,y in zip(position_consec,position_consec[1:])])
			#print "\n length diff_consec:", len(diff_consec)
			#print 'diff_consec:', diff_consec
			count_consec = 0
			pass_consec = [0]*len(diff_consec)
			eval_flag = 0
			#print "\n pos 66356:66363 diff_consec:", diff_consec[66356:66363]
			pass_consec[0] = 1
			#print "\n all diff_consec:", diff_consec
			for i in diff_consec:
				#print 'count_consec:', count_consec
				if count_consec==len(diff_consec):
					next
				if (count_consec<(len(diff_consec)-1)):
					#print '\n'
					#print 'count_consec:', count_consec, ', i:', i, ', position_consec:', position_consec[count_consec], ', diff_consec+1:',diff_consec[count_consec+1]
					if i == 1:
						#print 'normal'
						pass_consec[count_consec] = 1 
						#print 'pass_consec[count_consec]:', pass_consec[count_consec]
						eval_flag = 1
					if i > 1 and diff_consec[count_consec-1]==1:
						#print 'one i > 1 : normal sequence gap'
						pass_consec[count_consec] = 1	
						#print 'pass_consec[count_consec]:', pass_consec[count_consec]
						eval_flag = 1
					if i > 1 and diff_consec[count_consec-1]>1:
						#print 'set of two i > 1 : non-consecutive error'
						pass_consec[count_consec] = 0
						#print 'pass_consec[count_consec]:', pass_consec[count_consec]
						eval_flag = 1
					if eval_flag == 0:
						print 'not evaluated!'
					eval_flag = 0
					count_consec = count_consec + 1
				sys.stdout.flush()
			#print "\n pass_consec:", pass_consec
			#print "\n length pass_consec:", len(pass_consec)
			#print "\n Python line num 66358:", pass_consec[66358], "position_consec[66358]:", position_consec[66358] ,"diff_consec[66358]:", diff_consec[66358]
			#print "\n Python line num 66359:", pass_consec[66359], "position_consec[66359]:", position_consec[66359], "diff_consec[66359]:", diff_consec[66359]
			#print "\n Python line num 66360:", pass_consec[66360], "position_consec[66360]:", position_consec[66360], "diff_consec[66360]:", diff_consec[66360]
			#print "\n Python line num 66361:", pass_consec[66361], "position_consec[66361]:", position_consec[66361], "diff_consec[66361]:", diff_consec[66361]
			#print "\n Python line num 66362:", pass_consec[66362], "position_consec[66362]:", position_consec[66362], "diff_consec[66362]:", diff_consec[66362]

			#for thisPos in position_consec:
			#	diffCounters = abs(int(position_consec)-int(position_consec[thisPos+1]))	
			#	print thisPos, diffCounters
			with open(thisScaff + '.sync') as InFile2:
				if keepTemp ==1:
					OutFileTemp = open(thisScaff + '_OutputSites.txt', 'w')
					OutFileTemp_dpth = open(thisScaff + '_OutputSites_dpth.txt', 'w')
					# collecting line and position indexes
					# depth file
					lineNum_dpth = 0
					lineNumCounter_dpth = []
					positionCounter_dpth = []
					# sites file
					lineNum_sites = 0
					lineNumCounter_sites = []
					positionCounter_sites = []
				if keepTemp == 0:
					OutFileTemp = open(thisLG + '_OutputSites.txt', 'w')
					OutFileTemp_dpth = open(thisLG + '_OutputSites_dpth.txt', 'w')
					# collecting line and position indexes
					# depth file
					lineNum_dpth = 0
					lineNumCounter_dpth = []
					positionCounter_dpth = []
					# sites file
					lineNum_sites = 0
					lineNumCounter_sites = []
					positionCounter_sites = []
				# Initialise
				startFlag=0; stopcounter=0; switchOn=0; numNonBiallele=0; dataOutputFile=0; lineNumberScaffold=0
				progDisp = range(1,100); progDispAdj = range(1,100); LineNum_sync = 0; posCounter_previous =0; posCounter=0; badSeq = 0
				countProg = int(0)
				print "\n Processing sites..."
				#print 'length pass_consec:', len(pass_consec)
				OutFileTemp.write(OutputStringHeaderTempFile+"\n")
				OutFileTemp_dpth.write(OutputStringHeaderTempFile_dpth+"\n")
				for line in InFile2:
					#print '\n'
					#print 'line:', line
					#print 'LineNum_sync:', LineNum_sync, 'pass_consec[LineNum_sync]:', pass_consec[LineNum_sync]
					OutputString = 1
					#if (LineNum_sync < len(pass_consec) and pass_consec[LineNum_sync]==0):
						#print 'skip non-consecutive:', line
					#	next
					if (LineNum_sync < len(pass_consec)-1):
						PercentCompleteSites = round((float(lineNumberScaffold)/(float(InFile2lineNum)-1))*100,2)
						PercentCompleteSites = float(PercentCompleteSites)
						lineNumberScaffold = lineNumberScaffold+1
						if PercentCompleteSites in progDispAdj:
							sys.stdout.write('%s' % ('.'))
							sys.stdout.flush()
							countProg=int(countProg)+1
							progDispAdj=progDisp[countProg:]	
						# this line extract
						LineList2 = line.split()[0:]
						if 	pass_consec[LineNum_sync]==1:
							#print 'pass'
							#print '\n LineList2: ', LineList2 
							currentScaff = str(LineList2[0])
							currentPos = int(LineList2[1])
							posCounter = currentPos
							refAllele = str(LineList2[2])
							position_passed = []
							pops_allAdjustments = defaultdict(list)
							for thisPop in pops:
								# Python version
								pops_allAdjustments[thisPop].append(adjustData(LineList2,thisPop,minDepthSeqInfo, maxDepthSeqInfo, minAlleleCount))
								# C version
								#pops_allAdjustments[thisPop].append(pool_adjustments.adjustData(LineList2,thisPop,minDepthSeqInfo, maxDepthSeqInfo, minAlleleCount))
							allPools={}
							minDepthPassed = 0
							for thisPop in pops:
								minDepthPassed = float(pops_allAdjustments[thisPop][0]['minDepth_pass']) + minDepthPassed
							allPools['minDepth_pools'] = minDepthPassed
							allPools['minDepth_thresh'] = 0
							if (allPools['minDepth_pools'] >= minPoolNum):
								allPools['minDepth_thresh'] = 1
							# Check if locus has >2 bases with a different set of most common alleles present in pools
							totalCounts = [0]*4			
							for thisPop in pops:
								totalCounts = [q+w for q,w in zip(totalCounts,pops_allAdjustments[thisPop][0]['count_adj'])]
							allPools['count_allPools'] = totalCounts
							allPools['base_present'] = [1 if x>=minAlleleCount else 0 for x in allPools['count_allPools']]
							# Alleles present
							counter = 0
							seqCode = ["A","T","C","G"]
							allPools['bases_s'] = ['','','','']
							for allele in allPools['base_present']:
								allPools['bases_s'][counter] = int(allele)*seqCode[counter]
								counter = counter+1
							counter = 0	
							allPools['bases_t'] = ''.join(allPools['bases_s'])
							allPools['numAlleles'] = len(allPools['bases_t'])
							allPools['OneOrTwoBp'] = 0		
							if (allPools['minDepth_thresh'] == 0):
								# skips output if below threshold depth across pools
								next
							if (sum(allPools['base_present'])<=2 and allPools['minDepth_thresh'] != 0 and allPools['numAlleles'] <= 2):
								p=[]; q=[]; poly_pool=[]; dpthAdj=[]; dpthPass=[]; pi=[]; piAdj=[]
								dpthPolyPass=[]; pi_bar=[]; pi_bar_adj=[]; pi_T=[]; pi_T_adj=[]; d_xy_raw=[]; d_xy_fromFstAdj=[] 
								Fst_fromPi=[]; Fst_fromPiAdj=[]; Fst_fromDxy=[]
								allPools['OneOrTwoBp'] = 1	
								# Diversity statistics 
								pops_diversityStats = defaultdict(list)
								for thisPop in pops:
									# Python version
									pops_diversityStats[thisPop].append(diverseStats(pops_allAdjustments[thisPop][0],minAlleleCount))
									# C version
									#pops_diversityStats[thisPop].append(pool_diversityStats.diverseStats(pops_allAdjustments[thisPop][0],minAlleleCount))
								# Allele frequencies
								alleleFreq_allPops = alleleFreq(allPools,pops_allAdjustments)
								for thisPop in pops:
									p.append(alleleFreq_allPops[thisPop][0])
									dpthAdj.append(pops_allAdjustments[thisPop][0]['depth_adj'])
									dpthPass.append(pops_allAdjustments[thisPop][0]['minDepth_pass'])
								# Pairwise divergence statistics
								pairs_divergenceStats = defaultdict(list)
								for thisPair in comboPops_nums:
									# Python version
									pairs_divergenceStats[thisPair].append(divergeStats(pops_allAdjustments[thisPair[0]],pops_allAdjustments[thisPair[1]], \
									pops_diversityStats[thisPair[0]],pops_diversityStats[thisPair[1]], \
									alleleFreq_allPops[thisPair[0]][0],alleleFreq_allPops[thisPair[0]][1], \
									alleleFreq_allPops[thisPair[1]][0],alleleFreq_allPops[thisPair[1]][1],minAlleleCount))
								for thisPair in comboPops_nums:
									dpthPolyPass.append(pairs_divergenceStats[thisPair][0]['BothPassed'])
								# Depth file write
								OutputStringStart_dpth = [thisScaff,currentPos,thisLG,thiscM,refAllele,allPools['bases_t'],allPools['numAlleles'],allPools['OneOrTwoBp'],allPools['minDepth_pools'],allPools['minDepth_thresh']]
								lineNum_dpth = lineNum_dpth + 1
								if OutputString==1 and lineNum_dpth!=0:
									OutputStringData_dpth = OutputStringStart_dpth + p + dpthAdj + dpthPass + dpthPolyPass
									OutFileTemp_dpth.write("\t".join(map(str, OutputStringData_dpth))+"\n")
									lineNumCounter_dpth.append(float(lineNum_dpth))
									positionCounter_dpth.append(float(currentPos))	
							# new feature - only records polymorphic sites (if allSites==0)
							if (sum(allPools['base_present'])<=2 and allPools['minDepth_thresh'] != 0 and allPools['numAlleles'] == 2 and allSites == 0):
								outputStringStart = [thisScaff,currentPos,thisLG,thiscM,refAllele,allPools['bases_t'],allPools['numAlleles'],allPools['OneOrTwoBp'],allPools['minDepth_pools'],allPools['minDepth_thresh']]		
								p=[]; q=[]; poly_pool=[]; dpthAdj=[]; dpthPass=[]; pi=[]; piAdj=[]
								dpthPolyPass=[]; pi_bar=[]; pi_bar_adj=[]; pi_T=[]; pi_T_adj=[]; d_xy_raw=[]; d_xy_fromFstAdj=[] 
								Fst_fromPi=[]; Fst_fromPiAdj=[]; Fst_fromDxy=[]
								for thisPop in pops:
									p.append(alleleFreq_allPops[thisPop][0])
									q.append(alleleFreq_allPops[thisPop][1])
									poly_pool.append(pops_allAdjustments[thisPop][0]['poly'])
									dpthAdj.append(pops_allAdjustments[thisPop][0]['depth_adj'])
									dpthPass.append(pops_allAdjustments[thisPop][0]['minDepth_pass'])
									pi.append(pops_diversityStats[thisPop][0]['Pi'])
									piAdj.append(pops_diversityStats[thisPop][0]['Pi_adj'])
								for thisPair in comboPops_nums:
									dpthPolyPass.append(pairs_divergenceStats[thisPair][0]['BothPassed'])
									pi_bar.append(pairs_divergenceStats[thisPair][0]['Pi_bar'])
									pi_bar_adj.append(pairs_divergenceStats[thisPair][0]['Pi_bar_adj'])
									pi_T.append(pairs_divergenceStats[thisPair][0]['Pi_T'])
									pi_T_adj.append(pairs_divergenceStats[thisPair][0]['Pi_T_adj'])
									d_xy_raw.append(pairs_divergenceStats[thisPair][0]['d_xy_raw'])
									d_xy_fromFstAdj.append(pairs_divergenceStats[thisPair][0]['d_xy_fromFstAdj'])
									Fst_fromPi.append(pairs_divergenceStats[thisPair][0]['Fst_fromPi'])
									Fst_fromPiAdj.append(pairs_divergenceStats[thisPair][0]['Fst_fromPiAdj'])
									Fst_fromDxy.append(pairs_divergenceStats[thisPair][0]['Fst_fromDxy'])
								# Main site output write
								lineNum_sites = lineNum_sites + 1
								if OutputString==1 and lineNum_sites!=0:
									OutputStringData = outputStringStart + p + q + poly_pool + dpthAdj + dpthPass + pi + piAdj + dpthPolyPass + pi_bar + pi_bar_adj + pi_T \
									+ pi_T_adj + d_xy_raw + d_xy_fromFstAdj + Fst_fromPi + Fst_fromPiAdj + Fst_fromDxy
									OutFileTemp.write("\t".join(map(str, OutputStringData))+"\n")
									lineNumCounter_sites.append(float(lineNum_sites))
									positionCounter_sites.append(float(currentPos))	
							# previous approach still available - records all sites (if allSites==1). Essentially becomes the length of output as dpth file
							if (sum(allPools['base_present'])<=2 and allPools['minDepth_thresh'] != 0 and allPools['numAlleles'] <= 2 and allSites == 1):
								outputStringStart = [thisScaff,currentPos,thisLG,thiscM,refAllele,allPools['bases_t'],allPools['numAlleles'],allPools['OneOrTwoBp'],allPools['minDepth_pools'],allPools['minDepth_thresh']]		
								p=[]; q=[]; poly_pool=[]; dpthAdj=[]; dpthPass=[]; pi=[]; piAdj=[]
								dpthPolyPass=[]; pi_bar=[]; pi_bar_adj=[]; pi_T=[]; pi_T_adj=[]; d_xy_raw=[]; d_xy_fromFstAdj=[] 
								Fst_fromPi=[]; Fst_fromPiAdj=[]; Fst_fromDxy=[]
								for thisPop in pops:
									p.append(alleleFreq_allPops[thisPop][0])
									q.append(alleleFreq_allPops[thisPop][1])
									poly_pool.append(pops_allAdjustments[thisPop][0]['poly'])
									dpthAdj.append(pops_allAdjustments[thisPop][0]['depth_adj'])
									dpthPass.append(pops_allAdjustments[thisPop][0]['minDepth_pass'])
									pi.append(pops_diversityStats[thisPop][0]['Pi'])
									piAdj.append(pops_diversityStats[thisPop][0]['Pi_adj'])
								for thisPair in comboPops_nums:
									dpthPolyPass.append(pairs_divergenceStats[thisPair][0]['BothPassed'])
									pi_bar.append(pairs_divergenceStats[thisPair][0]['Pi_bar'])
									pi_bar_adj.append(pairs_divergenceStats[thisPair][0]['Pi_bar_adj'])
									pi_T.append(pairs_divergenceStats[thisPair][0]['Pi_T'])
									pi_T_adj.append(pairs_divergenceStats[thisPair][0]['Pi_T_adj'])
									d_xy_raw.append(pairs_divergenceStats[thisPair][0]['d_xy_raw'])
									d_xy_fromFstAdj.append(pairs_divergenceStats[thisPair][0]['d_xy_fromFstAdj'])
									Fst_fromPi.append(pairs_divergenceStats[thisPair][0]['Fst_fromPi'])
									Fst_fromPiAdj.append(pairs_divergenceStats[thisPair][0]['Fst_fromPiAdj'])
									Fst_fromDxy.append(pairs_divergenceStats[thisPair][0]['Fst_fromDxy'])
								# Main site output write
								lineNum_sites = lineNum_sites + 1
								if OutputString==1 and lineNum_sites!=0:
									OutputStringData = outputStringStart + p + q + poly_pool + dpthAdj + dpthPass + pi + piAdj + dpthPolyPass + pi_bar + pi_bar_adj + pi_T \
									+ pi_T_adj + d_xy_raw + d_xy_fromFstAdj + Fst_fromPi + Fst_fromPiAdj + Fst_fromDxy
									OutFileTemp.write("\t".join(map(str, OutputStringData))+"\n")
									lineNumCounter_sites.append(float(lineNum_sites))
									positionCounter_sites.append(float(currentPos))	
						LineNum_sync = LineNum_sync + 1
				OutFileTemp.close()
				OutFileTemp_dpth.close()
			#InFile2.close()
		#############################
		# Calculate window averages #
		#############################
		if sitesOnly == 1:
			print "\n Site analyses only ... exiting "
			print " "
			print " "
			print "*** Analyses complete ***"
			print " "
  			quit()
		#InFile2 = open(thisScaff + '.sync', 'r')
		#InFile2lines  = [i.split()[0:] for i in InFile2.readlines()]
		#InFile2.close() 
		#InFile2lineNum = len(InFile2lines)
		#totalScaffSize = InFile2lines[InFile2lineNum-2][1]
		totalScaffSize = 2000000
		#print '\n totalScaffSize: ', totalScaffSize
		# New window output file name
		OutFileWindows = open(thisScaff + '_OutputWindows_' + sys.argv[3] + '.txt', 'w')
		OutFileWindows.write(OutputStringHeaderSlidingWindows+"\n")
		#print "\n positionCounter_dpth: ", positionCounter_dpth
		if float(totalScaffSize) > 15000:
			thisWindowCounter = 1; PercentCompleteWindows=0; lineNum_sites = 0
			# collecting line and position indexes
			# sites file
			# start and end of scaffold to calculate total window number
			bpStart_scaffold = positionCounter_dpth[0]
			bpEnd_scaffold = positionCounter_dpth[-1]
			firstWindowStart = round_down(bpStart_scaffold, windowSize)
			lastWindowEnd = int(-(-float(bpEnd_scaffold)//windowSize)*windowSize)
 			#print '\n bpStart_scaffold: ', bpStart_scaffold
 			#print '\n bpEnd_scaffold: ', bpEnd_scaffold
 			#print '\n firstWindowStart:', firstWindowStart
 			#print '\n lastWindowEnd:', lastWindowEnd
#			OutFileTempAgain.close() 
			# bpStart = 0
			thisStart = firstWindowStart 
			thisEnd = thisStart + windowSize
			totalFullWindows = (lastWindowEnd-firstWindowStart)/(windowSize-overlapSize)
			windowList = range(1,int(totalFullWindows))
			print '\n Total windows ~ ', totalFullWindows
			print "\n Processing windows..."
			for thisWindow in windowList:
 				startCounter = 0
				sys.stdout.write('%s' % ('.'+ str(thisWindow)))
				sys.stdout.flush()
				# for dpth file
				collect_p_dpthFile = [0]*numPops
				collect_dpthAdj_dpthFile = [0]*numPops
				collect_dpthPass_dpthFile = [0]*numPops
				collect_dpthPolyPass_dpthFile = [0]*len(comboPops_nums)
				count_dpthAdj_dpthFile = [0]*numPops
				count_p_dpthFile = [0]*numPops
				count_dpthPass_dpthFile = [0]*numPops
				count_dpthPolyPass_dpthFile = [0]*len(comboPops_nums)
				collect_dpthAdj_dpthFile = [float(x) for x in collect_dpthAdj_dpthFile]
				collect_dpthPass_dpthFile = [float(x) for x in collect_dpthPass_dpthFile]
				collect_dpthPolyPass_dpthFile = [float(x) for x in collect_dpthPolyPass_dpthFile]
				count_dpthAdj_dpthFile = [float(x) for x in count_dpthAdj_dpthFile]
				count_dpthPass_dpthFile = [float(x) for x in count_dpthPass_dpthFile]
				count_dpthPolyPass_dpthFile = [float(x) for x in count_dpthPolyPass_dpthFile]
				# for sites file
				collect_poly = [0]*numPops
				count_poly = [0]*numPops
				collect_pi = [0]*numPops
				count_pi = [0]*numPops
				collect_piAdj = [0]*numPops
				count_piAdj = [0]*numPops
				collect_PairPass = [0]*len(comboPops_nums)
				count_PairPass = [0]*len(comboPops_nums)
				collect_Pibar = [0]*len(comboPops_nums)
				count_Pibar = [0]*len(comboPops_nums)
				collect_PibarAdj = [0]*len(comboPops_nums)
				count_PibarAdj = [0]*len(comboPops_nums)
				collect_PiT = [0]*len(comboPops_nums)
				count_PiT = [0]*len(comboPops_nums)
				collect_PiTAdj = [0]*len(comboPops_nums)
				count_PiTAdj = [0]*len(comboPops_nums)
				collect_dXYraw = [0]*len(comboPops_nums)
				count_dXYraw = [0]*len(comboPops_nums)	
				collect_dXYfromFst = [0]*len(comboPops_nums)
				count_dXYfromFst = [0]*len(comboPops_nums)
				collect_FstPi = [0]*len(comboPops_nums)
				count_FstPi = [0]*len(comboPops_nums)
				collect_FstPiAdj = [0]*len(comboPops_nums)
				count_FstPiAdj = [0]*len(comboPops_nums)
				collect_FstDxy =[0]*len(comboPops_nums)
				count_FstDxy = [0]*len(comboPops_nums)
				collect_pi = [float(x) for x in collect_pi]
				count_pi = [float(x) for x in count_pi]
				collect_piAdj = [float(x) for x in collect_piAdj]
				count_piAdj = [float(x) for x in count_piAdj]
				collect_PairPass = [float(x) for x in collect_PairPass]
				count_PairPass = [float(x) for x in count_PairPass]
				collect_Pibar = [float(x) for x in collect_Pibar]
				count_Pibar = [float(x) for x in count_Pibar]
				collect_PibarAdj = [float(x) for x in collect_PibarAdj]
				count_PibarAdj = [float(x) for x in count_PibarAdj]
				collect_PiT = [float(x) for x in collect_PiT]
				count_PiT = [float(x) for x in count_PiT]
				collect_PiTAdj = [float(x) for x in collect_PiTAdj]
				count_PiTAdj = [float(x) for x in count_PiTAdj]
				collect_dXYraw = [float(x) for x in collect_dXYraw]
				count_dXYraw = [float(x) for x in count_dXYraw]
				collect_dXYfromFst = [float(x) for x in collect_dXYfromFst]
				count_dXYfromFst = [float(x) for x in count_dXYfromFst]
				collect_FstPi = [float(x) for x in collect_FstPi]
				count_FstPi = [float(x) for x in count_FstPi]
				collect_FstPiAdj = [float(x) for x in collect_FstPiAdj]
				count_FstPiAdj = [float(x) for x in count_FstPiAdj]
				collect_FstDxy = [float(x) for x in collect_FstDxy]
				count_FstDxy = [float(x) for x in count_FstDxy]

				if keepTemp == 1:
					indexRange = []
					indexRange_dpth=[]
					# indexRange collects all the positions which lie within the current window range 
					indexRange = [i for i, x in enumerate(positionCounter_sites) if x >= thisStart and x < thisEnd]
					indexRange_dpth = [i for i, x in enumerate(positionCounter_dpth) if x >= thisStart and x < thisEnd]					
					#print 'indexRange:', indexRange
					#print 'indexRange_dpth:', indexRange_dpth
					# these arguments designed to remove empty windows
					if len(indexRange_dpth)>100 and len(indexRange)>10:
						# dpth file
						lineNumsFile_dpth = lineNumCounter_dpth[min(indexRange_dpth):max(indexRange_dpth)]
						lineMin_dpth = int(min(lineNumsFile_dpth))
						lineMax_dpth = int(max(lineNumsFile_dpth))+1
						#print 'lineMin_dpth:', lineMin_dpth
						#print 'lineMax_dpth:', lineMax_dpth
						with open(thisScaff + '_OutputSites_dpth.txt') as OutFileTempAgain_dpth:
							for lineRaw in itertools.islice(OutFileTempAgain_dpth, lineMin_dpth, lineMax_dpth):
								line = lineRaw.split()[0:]
								currentPos = line[1]
								currentPos = float(currentPos)
								minDepthPools = float(line[8])
								Biallelic = float(line[7])
								#print '\n minDepthPools: ', minDepthPools
								# note depth data has to be multiplied by the number of stats
								if (minDepthPools >= minPoolNum and Biallelic==1):
									vals_p = line[statsColumnList_dpth['p'][0]:statsColumnList_dpth['p'][1]]
									vals_dpthAdj = line[statsColumnList_dpth['dpthAdj'][0]:statsColumnList_dpth['dpthAdj'][1]]
									vals_dpthPass = line[statsColumnList_dpth['dpthPass'][0]:statsColumnList_dpth['dpthPass'][1]]	
									vals_dpthPolyPass = line[statsColumnList_Pairs_dpth['dpthPolyPass'][0]:statsColumnList_Pairs_dpth['dpthPolyPass'][1]]									
									# single pop stats for this line
									thisVals_p_dpthFile = [float(x) for x in vals_p]
									thisVals_dpthAdj_dpthFile = [float(x) for x in vals_dpthAdj]
									thisVals_dpthPass_dpthFile = [float(x) for x in vals_dpthPass]
									# pop pair stats for this line
									thisVals_dpthPolyPass_dpthFile = [float(x) for x in vals_dpthPolyPass]
									for i,j in enumerate(thisVals_p_dpthFile):
										collect_p_dpthFile[i] += j
										count_p_dpthFile[i] += 1
									for i,j in enumerate(thisVals_dpthAdj_dpthFile):
										collect_dpthAdj_dpthFile[i] += j
										count_dpthAdj_dpthFile[i] += 1
									for i,j in enumerate(thisVals_dpthPass_dpthFile):
										collect_dpthPass_dpthFile[i] += j
										count_dpthPass_dpthFile[i] += 1
									for i,j in enumerate(thisVals_dpthPolyPass_dpthFile):
										collect_dpthPolyPass_dpthFile[i] += j
										count_dpthPolyPass_dpthFile[i] += 1
						# sites file
						lineNumsFile = lineNumCounter_sites[min(indexRange):max(indexRange)]
						lineMin = int(min(lineNumsFile))
						lineMax = int(max(lineNumsFile))+1		
						with open(thisScaff + '_OutputSites.txt') as OutFileTempAgain:
							for lineRaw in itertools.islice(OutFileTempAgain, lineMin, lineMax):
								line = lineRaw.split()[0:]
								currentPos = line[1]
								currentPos = float(currentPos)
								minDepthPools = float(line[8])
								Biallelic = float(line[7])
								if (minDepthPools >= minPoolNum and Biallelic==1):
									# Single pool - 7 stats 
									# 7 lots of single pop statistics 
									# but only 3 required
									#vals_P = line[statsColumnList['p'][0]:statsColumnList['p'][1]]
									#vals_Q = line[statsColumnList['q'][0]:statsColumnList['q'][1]]
									#vals_poly = line[statsColumnList['poly_pool'][0]:statsColumnList['poly_pool'][1]]
									#vals_dpthAdj = line[statsColumnList['dpthAdj'][0]:statsColumnList['dpthAdj'][1]]
									#vals_dpthPass = line[statsColumnList['dpthPass'][0]:statsColumnList['dpthPass'][1]]
									vals_poly = line[statsColumnList['poly_pool'][0]:statsColumnList['poly_pool'][1]]
									vals_pi = line[statsColumnList['pi'][0]:statsColumnList['pi'][1]]
									vals_piAdj = line[statsColumnList['piAdj'][0]:statsColumnList['piAdj'][1]]
									thisVals_poly = [float(x) for x in vals_poly]
									thisVals_pi = [float(x) for x in vals_pi]
									thisVals_piAdj = [float(x) for x in vals_piAdj]
									# Pair statistics - 10 lots of pair-wise statistics
									vals_PairPass = line[statsColumnListPairs['dpthPolyPass'][0]:statsColumnListPairs['dpthPolyPass'][1]]
									vals_Pibar = line[statsColumnListPairs['piBar'][0]:statsColumnListPairs['piBar'][1]]
									vals_PibarAdj = line[statsColumnListPairs['piBarAdj'][0]:statsColumnListPairs['piBarAdj'][1]]
									vals_PiT = line[statsColumnListPairs['piT'][0]:statsColumnListPairs['piT'][1]]
									vals_PiTAdj = line[statsColumnListPairs['piTadj'][0]:statsColumnListPairs['piTadj'][1]]
									vals_dXYraw = line[statsColumnListPairs['dXYraw'][0]:statsColumnListPairs['dXYraw'][1]]
									vals_dXYfromFst = line[statsColumnListPairs['dXYfromFstadj'][0]:statsColumnListPairs['dXYfromFstadj'][1]]
									vals_FstPi = line[statsColumnListPairs['FstPi'][0]:statsColumnListPairs['FstPi'][1]]
									vals_FstPiAdj = line[statsColumnListPairs['FstPiAdj'][0]:statsColumnListPairs['FstPiAdj'][1]]
									vals_FstDxy = line[statsColumnListPairs['FstfromDxy'][0]:statsColumnListPairs['FstfromDxy'][1]]								
									thisVals_PairPass = [float(x) for x in vals_PairPass]
									thisVals_Pibar = [float(x) for x in vals_Pibar]
									thisVals_PibarAdj = [float(x) for x in vals_PibarAdj]
									thisVals_PiT = [float(x) for x in vals_PiT]
									thisVals_PiTAdj = [float(x) for x in vals_PiTAdj]
									thisVals_dXYraw = [float(x) for x in vals_dXYraw]
									thisVals_dXYfromFst = [float(x) for x in vals_dXYfromFst]
									thisVals_FstPi = [float(x) for x in vals_FstPi]
									thisVals_FstPiAdj = [float(x) for x in vals_FstPiAdj]
									thisVals_FstDxy = [float(x) for x in vals_FstDxy]
									for i,j in enumerate(thisVals_pi):
										collect_pi[i] += j
										count_pi[i] += 1
									for i,j in enumerate(thisVals_poly):
										collect_poly[i] += j
										count_poly[i] += 1
									for i,j in enumerate(thisVals_piAdj):
										collect_piAdj[i] += j
										count_piAdj[i] += 1
									for i,j in enumerate(thisVals_PairPass):
										collect_PairPass[i] += j
										count_PairPass[i] += 1
									for i,j in enumerate(thisVals_Pibar):
										collect_Pibar[i] += j
										count_Pibar[i] += 1
									for i,j in enumerate(thisVals_PibarAdj):
										collect_PibarAdj[i] += j
										count_PibarAdj[i] += 1
									for i,j in enumerate(thisVals_PiT):
										collect_PiT[i] += j
										count_PiT[i] += 1
									for i,j in enumerate(thisVals_PiTAdj):
										collect_PiTAdj[i] += j
										count_PiTAdj[i] += 1	
									for i,j in enumerate(thisVals_dXYraw):
										collect_dXYraw[i] += j
										count_dXYraw[i] += 1	
									for i,j in enumerate(thisVals_dXYfromFst):
										collect_dXYfromFst[i] += j
										count_dXYfromFst[i] += 1
									for i,j in enumerate(thisVals_FstPi):
										collect_FstPi[i] += j
										count_FstPi[i] += 1
									for i,j in enumerate(thisVals_FstPiAdj):
										collect_FstPiAdj[i] += j
										count_FstPiAdj[i] += 1
									for i,j in enumerate(thisVals_FstDxy):
										collect_FstDxy[i] += j
										count_FstDxy[i] += 1
						OutFileTempAgain_dpth.close()
						OutFileTempAgain.close()
				if keepTemp == 0:
					indexRange = []
					indexRange_dpth=[]
					# indexRange collects all the positions which lie within the current window range 
					indexRange = [i for i, x in enumerate(positionCounter_sites) if x >= thisStart and x < thisEnd]
					indexRange_dpth = [i for i, x in enumerate(positionCounter_dpth) if x >= thisStart and x < thisEnd]

					# these arguments designed to remove empty windows
					if len(indexRange_dpth)>100 and len(indexRange)>10:
						# dpth file
						lineNumsFile_dpth = lineNumCounter_dpth[min(indexRange_dpth):max(indexRange_dpth)]
						lineMin_dpth = int(min(lineNumsFile_dpth))
						lineMax_dpth = int(max(lineNumsFile_dpth))+1
						with open(thisLG + '_OutputSites_dpth.txt') as OutFileTempAgain_dpth:
							for lineRaw in itertools.islice(OutFileTempAgain_dpth, lineMin_dpth, lineMax_dpth):
								#print '\n lineRaw :', lineRaw
								line = lineRaw.split()[0:]
								currentPos = line[1]
								currentPos = float(currentPos)
								minDepthPools = float(line[8])
								Biallelic = float(line[7])
								#print '\n minDepthPools: ', minDepthPools
								# note depth data has to be multiplied by the number of stats
								if (minDepthPools >= minPoolNum and Biallelic==1):
									vals_p = line[statsColumnList_dpth['p'][0]:statsColumnList_dpth['p'][1]]
									vals_dpthAdj = line[statsColumnList_dpth['dpthAdj'][0]:statsColumnList_dpth['dpthAdj'][1]]
									vals_dpthPass = line[statsColumnList_dpth['dpthPass'][0]:statsColumnList_dpth['dpthPass'][1]]	
									vals_dpthPolyPass = line[statsColumnList_Pairs_dpth['dpthPolyPass'][0]:statsColumnList_Pairs_dpth['dpthPolyPass'][1]]									
									# single pop stats for this line
									thisVals_p_dpthFile = [float(x) for x in vals_p]
									thisVals_dpthAdj_dpthFile = [float(x) for x in vals_dpthAdj]
									thisVals_dpthPass_dpthFile = [float(x) for x in vals_dpthPass]
									# pop pair stats for this line
									thisVals_dpthPolyPass_dpthFile = [float(x) for x in vals_dpthPolyPass]
									for i,j in enumerate(thisVals_p_dpthFile):
										collect_p_dpthFile[i] += j
										count_p_dpthFile[i] += 1
									for i,j in enumerate(thisVals_dpthAdj_dpthFile):
										collect_dpthAdj_dpthFile[i] += j
										count_dpthAdj_dpthFile[i] += 1
									for i,j in enumerate(thisVals_dpthPass_dpthFile):
										collect_dpthPass_dpthFile[i] += j
										count_dpthPass_dpthFile[i] += 1
									for i,j in enumerate(thisVals_dpthPolyPass_dpthFile):
										collect_dpthPolyPass_dpthFile[i] += j
										count_dpthPolyPass_dpthFile[i] += 1
						# sites file
						lineNumsFile = lineNumCounter_sites[min(indexRange):max(indexRange)]
						lineMin = int(min(lineNumsFile))
						lineMax = int(max(lineNumsFile))+1		
						with open(thisLG + '_OutputSites.txt') as OutFileTempAgain:
							for lineRaw in itertools.islice(OutFileTempAgain, lineMin, lineMax):
								line = lineRaw.split()[0:]
								currentPos = line[1]
								currentPos = float(currentPos)
								minDepthPools = float(line[8])
								Biallelic = float(line[7])
								if (minDepthPools >= minPoolNum and Biallelic==1):
									# Single pool - 7 stats 
									# 7 lots of single pop statistics 
									# but only 3 required
									#vals_P = line[statsColumnList['p'][0]:statsColumnList['p'][1]]
									#vals_Q = line[statsColumnList['q'][0]:statsColumnList['q'][1]]
									#vals_poly = line[statsColumnList['poly_pool'][0]:statsColumnList['poly_pool'][1]]
									#vals_dpthAdj = line[statsColumnList['dpthAdj'][0]:statsColumnList['dpthAdj'][1]]
									#vals_dpthPass = line[statsColumnList['dpthPass'][0]:statsColumnList['dpthPass'][1]]
									vals_poly = line[statsColumnList['poly_pool'][0]:statsColumnList['poly_pool'][1]]
									vals_pi = line[statsColumnList['pi'][0]:statsColumnList['pi'][1]]
									vals_piAdj = line[statsColumnList['piAdj'][0]:statsColumnList['piAdj'][1]]
									thisVals_poly = [float(x) for x in vals_poly]
									thisVals_pi = [float(x) for x in vals_pi]
									thisVals_piAdj = [float(x) for x in vals_piAdj]
									# Pair statistics - 10 lots of pair-wise statistics
									vals_PairPass = line[statsColumnListPairs['dpthPolyPass'][0]:statsColumnListPairs['dpthPolyPass'][1]]
									vals_Pibar = line[statsColumnListPairs['piBar'][0]:statsColumnListPairs['piBar'][1]]
									vals_PibarAdj = line[statsColumnListPairs['piBarAdj'][0]:statsColumnListPairs['piBarAdj'][1]]
									vals_PiT = line[statsColumnListPairs['piT'][0]:statsColumnListPairs['piT'][1]]
									vals_PiTAdj = line[statsColumnListPairs['piTadj'][0]:statsColumnListPairs['piTadj'][1]]
									vals_dXYraw = line[statsColumnListPairs['dXYraw'][0]:statsColumnListPairs['dXYraw'][1]]
									vals_dXYfromFst = line[statsColumnListPairs['dXYfromFstadj'][0]:statsColumnListPairs['dXYfromFstadj'][1]]
									vals_FstPi = line[statsColumnListPairs['FstPi'][0]:statsColumnListPairs['FstPi'][1]]
									vals_FstPiAdj = line[statsColumnListPairs['FstPiAdj'][0]:statsColumnListPairs['FstPiAdj'][1]]
									vals_FstDxy = line[statsColumnListPairs['FstfromDxy'][0]:statsColumnListPairs['FstfromDxy'][1]]
									thisVals_PairPass = [float(x) for x in vals_PairPass]
									thisVals_Pibar = [float(x) for x in vals_Pibar]
									thisVals_PibarAdj = [float(x) for x in vals_PibarAdj]
									thisVals_PiT = [float(x) for x in vals_PiT]
									thisVals_PiTAdj = [float(x) for x in vals_PiTAdj]
									thisVals_dXYraw = [float(x) for x in vals_dXYraw]
									thisVals_dXYfromFst = [float(x) for x in vals_dXYfromFst]
									thisVals_FstPi = [float(x) for x in vals_FstPi]
									thisVals_FstPiAdj = [float(x) for x in vals_FstPiAdj]
									thisVals_FstDxy = [float(x) for x in vals_FstDxy]
									for i,j in enumerate(thisVals_pi):
										collect_pi[i] += j
										count_pi[i] += 1
									for i,j in enumerate(thisVals_poly):
										collect_poly[i] += j
										count_poly[i] += 1	
									for i,j in enumerate(thisVals_piAdj):
										collect_piAdj[i] += j
										count_piAdj[i] += 1
									for i,j in enumerate(thisVals_PairPass):
										collect_PairPass[i] += j
										count_PairPass[i] += 1
									for i,j in enumerate(thisVals_Pibar):
										collect_Pibar[i] += j
										count_Pibar[i] += 1
									for i,j in enumerate(thisVals_PibarAdj):
										collect_PibarAdj[i] += j
										count_PibarAdj[i] += 1
									for i,j in enumerate(thisVals_PiT):
										collect_PiT[i] += j
										count_PiT[i] += 1
									for i,j in enumerate(thisVals_PiTAdj):
										collect_PiTAdj[i] += j
										count_PiTAdj[i] += 1	
									for i,j in enumerate(thisVals_dXYraw):
										collect_dXYraw[i] += j
										count_dXYraw[i] += 1	
									for i,j in enumerate(thisVals_dXYfromFst):
										collect_dXYfromFst[i] += j
										count_dXYfromFst[i] += 1
									for i,j in enumerate(thisVals_FstPi):
										collect_FstPi[i] += j
										count_FstPi[i] += 1
									for i,j in enumerate(thisVals_FstPiAdj):
										collect_FstPiAdj[i] += j
										count_FstPiAdj[i] += 1
									for i,j in enumerate(thisVals_FstDxy):
										collect_FstDxy[i] += j
										count_FstDxy[i] += 1
						#OutFileTempAgain_dpth.close()
						#OutFileTempAgain.close()
				######################
				# calculate averages #
				######################
			    # collect_dpthAdj_dpthFile/count_dpthPass_dpthFile = mean depth (dpthAdj in window output)
				# count_dpthPass_dpthFile/windowSize = coverage of window per pool (dpthPass in window output)
				# mean of this last one will be average window coverage (coverage in window output)
				# pair stats have to be divided by [collect_dpthPolyPass_dpthFile]*len(comboPops_nums))
				#if count_coverage_dpth==0:
					# write -9 to rows?
				if len(indexRange_dpth) != 0:
					windowSize_rep = [windowSize]*numPops
					windowSize_rep_pair = [windowSize]*len(comboPops_nums)
					windowSize_rep = [float(x) for x in windowSize_rep]
					thisWindow_meanDepth = [MyDivision(a,b) for a, b in zip(collect_dpthAdj_dpthFile,count_dpthPass_dpthFile)]
					thisWindow_coverage = [MyDivision(a,b) for a, b in zip(collect_dpthPass_dpthFile,windowSize_rep)]
					coverage = np.mean(thisWindow_coverage)
					#print 'collect_dpthPass_dpthFile:' ,collect_dpthPass_dpthFile
					#print 'windowSize_rep:' ,windowSize_rep
					#print 'thisWindow_coverage:' ,thisWindow_coverage
					#print 'coverage = np.mean(thisWindow_coverage):' ,coverage
					if coverage!=0:	
						# First set of columns
						firstComponent = []
						# note need to fix -9 value which was window coverage. Will this be the mean window coverage?
						firstComponent = [thisScaff,thisLG,thiscM,thisStart,thisEnd,coverage]
						# new means
						# single pop stats
						# will come from each statistic / number of sites with reads above threshold (count_dpthPass_dpthFile) 
						p_mean = [MyDivision(a,b) for a, b in zip(collect_p_dpthFile,collect_dpthPass_dpthFile)]
						polyPool_mean = [MyDivision(a,b) for a, b in zip(collect_poly,collect_dpthPass_dpthFile)]
						pi_mean = [MyDivision(a,b) for a, b in zip(collect_pi,collect_dpthPass_dpthFile)]
						piAdj_mean = [MyDivision(a,b) for a, b in zip(collect_piAdj,collect_dpthPass_dpthFile)]
# 						print '\n collect_p_dpthFile: ', collect_p_dpthFile
# 						print '\n collect_dpthPass_dpthFile: ', collect_dpthPass_dpthFile
# 						print 'get p _ mean'
# 						print '\n p_mean: ', p_mean
# 						print '\n polyPool_mean: ', polyPool_mean
# 						print '\n pi_mean: ', pi_mean
# 						print '\n piAdj_mean: ', piAdj_mean
						# pair stats averaged
							# first one a little different, mean number of sites passed for pair/ total sites window
						PairPass_mean = [MyDivision(a,b) for a, b in zip(collect_dpthPolyPass_dpthFile,windowSize_rep_pair)]
							# the rest are stat/total passed for pair (
						Pibar_mean = [MyDivision(a,b) for a, b in zip(collect_Pibar,collect_dpthPolyPass_dpthFile)]
						PibarAdj_mean = [MyDivision(a,b) for a, b in zip(collect_PibarAdj,collect_dpthPolyPass_dpthFile)]
						PiT_mean = [MyDivision(a,b) for a, b in zip(collect_PiT,collect_dpthPolyPass_dpthFile)]
						PiTAdj_mean = [MyDivision(a,b) for a, b in zip(collect_PiTAdj,collect_dpthPolyPass_dpthFile)]
						dXYraw_mean = [MyDivision(a,b) for a, b in zip(collect_dXYraw,collect_dpthPolyPass_dpthFile)]
						dXYfromFst_mean = [MyDivision(a,b) for a, b in zip(collect_dXYfromFst,collect_dpthPolyPass_dpthFile)]
						FstPi_mean = [MyDivision(a,b) for a, b in zip(collect_FstPi,collect_dpthPolyPass_dpthFile)]
						FstPiAdj_mean = [MyDivision(a,b) for a, b in zip(collect_FstPiAdj,collect_dpthPolyPass_dpthFile)]
						FstDxy_mean = [MyDivision(a,b) for a, b in zip(collect_FstDxy,collect_dpthPolyPass_dpthFile)]

						# stats calculated from window means
						# FstfromPiMean
						firstBit = [MySubtraction(a,b) for a, b in zip(PiT_mean,Pibar_mean)]
						FstfromPiMean = [MyDivision(a,b) for a, b in zip(firstBit,PiT_mean)]
						# FstfromPiAdjMean
						firstBitAdj = [MySubtraction(a,b) for a, b in zip(PiTAdj_mean,PibarAdj_mean)]
						FstfromPiAdjMean = [MyDivision(a,b) for a, b in zip(firstBitAdj,PiTAdj_mean)]
						#print '\n FstfromPiMean :' ,FstfromPiMean
						#print '\n FstfromPiAdjMean :' ,FstfromPiAdjMean
						# dXYfromFstMean
						numero = [1+x if x!='NaN' else (x) for x in FstfromPiMean]
						denomo = [1-x if x!='NaN' else (x) for x in FstfromPiMean]
						#denomo = [1-x for x in FstfromPiMean]
						rightTerm = [MyDivision(a,b) for a, b in zip(numero,denomo)]
						dXYfromFstMean = [MyMultiplication(a, b) for a,b in zip(Pibar_mean, rightTerm)]
						#print '\n dXYfromFstMean :' ,dXYfromFstMean
						# dXYfromFstAdjMean
						numero = [1+x if x!='NaN' else (x) for x in FstfromPiAdjMean]
						denomo = [1-x if x!='NaN' else (x) for x in FstfromPiAdjMean]
						rightTerm = [MyDivision(a,b) for a, b in zip(numero,denomo)]
						dXYfromFstAdjMean = [MyMultiplication(a, b) for a,b in zip(PibarAdj_mean, rightTerm)]
						#print '\n dXYfromFstAdjMean :' ,dXYfromFstAdjMean
						# dAfromFstMean
						dAfromFstMean = []
						dAfromFstMean = [MySubtraction(a,b) for a,b in zip(FstfromPiMean, Pibar_mean)]
						# dAfromFstAdjMean
						dAfromFstAdjMean = []
						dAfromFstAdjMean = [MySubtraction(a,b) for a,b in zip(FstfromPiAdjMean, PibarAdj_mean)]
						toExport = firstComponent + polyPool_mean + thisWindow_meanDepth + thisWindow_coverage + pi_mean + piAdj_mean
						toExport = toExport + PairPass_mean + Pibar_mean + PibarAdj_mean + PiT_mean + PiTAdj_mean + dXYraw_mean + dXYfromFst_mean
						toExport = toExport + FstPi_mean + FstPiAdj_mean + FstDxy_mean + FstfromPiMean + FstfromPiAdjMean + dXYfromFstMean + dXYfromFstAdjMean
						toExport = toExport + dAfromFstMean + dAfromFstAdjMean
						toExport_convert = [str(i) for i in toExport]
						OutFileWindows.write("\t".join(map(str, toExport_convert))+"\n")
					if coverage == 0:
						# First set of columns
						firstComponent = []
						firstComponent = [thisScaff,thisLG,thiscM,thisStart,thisEnd,coverage]
						toExport = firstComponent
						toExport_convert = [str(i) for i in toExport]
						OutFileWindows.write("\t".join(map(str, toExport_convert))+"\n")
				PercentCompleteWindows = round((thisWindowCounter/(totalFullWindows))*100,2)
				PercentCompleteWindows = str(float(PercentCompleteWindows))
				thisWindowCounter=thisWindowCounter+1
				thisStart = int(thisEnd-float(overlapSize))
				thisEnd = int(thisStart+float(windowSize))
	LineNumber = LineNumber + 1
	LineNumber = float(int(LineNumber))
	numLines = float(int(numLines))
OutFileWindows.close()
print " "
print " "
print "*** Analyses complete ***"
print " "
