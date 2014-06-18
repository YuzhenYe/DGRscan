#!/usr/bin/env python

#Package Name: DGRscan
#First Released on June 15, 2014
#Developer: Yuzhen Ye <yye@indiana.edu>
#Affiliation: School of Informatics and Computing, Indiana University, Bloomington

import time
import sys

minAmut = 7
minAmutPer = 0.7
#minAmutSeg = 5
#minAmutPerSeg = 0.5
minAmutSeg = 7
minAmutPerSeg = 0.7
semi = False 

def myzeros(dim):
        ar = []
        for i in range(dim[0]):
                ar.append([0] * dim[1])
        return ar
def scorefun():
    	t = ['A', 'T', 'C', 'G', 'X', 'K', 'N']
    	similarityMatrixMap = dict()
    	for i in range(len(t)):
        	base1 = t[i]
        	for j in range(len(t)):
            		base2 = t[j]
	    		if base1 == base2:
				similarityMatrixMap[base1 + base2] = 1
	    		else:
				similarityMatrixMap[base1 + base2] = -1 
    	return similarityMatrixMap

def align(seq1, seq2, gap=-10):
	global semi;
	if semi == True:
		alignedSeq1, alignedSeq2, alignPointer, alignedSeq1Full, alignedSeq2Full, alignPointerFull, alignstart = computeFMatrixSemi(seq1, seq2, gap)
	else:
		alignedSeq1, alignedSeq2, alignPointer, alignedSeq1Full, alignedSeq2Full, alignPointerFull, alignstart = computeFMatrix(seq1, seq2, gap)
	return (alignedSeq1, alignedSeq2, alignPointer, alignedSeq1Full, alignedSeq2Full, alignPointerFull, alignstart)

def computeFMatrix(seq1, seq2, gap=-10): #gap=-10, large gap penalty to disable gap introduction
    	rows = len(seq1) + 1
    	cols = len(seq2) + 1
   	fMatrix = myzeros([rows, cols])
    	pointers = myzeros([rows, cols])
    	maxScore, iOfMax, jOfMax = 0, 0, 0

    	similarityMatrixMap = scorefun() 

    	for i in range(1, rows):
        	for j in range(1, cols):
            		mtch = fMatrix[i - 1][j - 1] + similarityMatrixMap[seq1[i - 1] + seq2[j - 1]]
            		delete = fMatrix[i - 1][j] + gap
            		insert = fMatrix[i][j - 1] + gap
            		fMatrix[i][j] = max(0, mtch, delete, insert)

            		if(fMatrix[i][j] == 0):
                		pointers[i][j] = -1
            		elif(fMatrix[i][j] == delete):
                		pointers[i][j] = 1
            		elif(fMatrix[i][j] == insert):
                		pointers[i][j] = 2
            		elif(fMatrix[i][j] == mtch):
                		pointers[i][j] = 3
            		if fMatrix[i][j] > maxScore :
                		iOfMax = i
                		jOfMax = j
                		maxScore = fMatrix[i][j]
    	(aligned1, aligned2, startOfAlign1, startOfAlign2) = traceBack(pointers, seq1, seq2, gap, similarityMatrixMap, iOfMax, jOfMax)

    	'''Some formatting for displaying alignment'''
    	numOfSpacesToAdd1 = {True: 0, False: startOfAlign2 - startOfAlign1}[startOfAlign1 >= startOfAlign2]
   	numOfSpacesToAdd2 = {True: 0, False: startOfAlign1 - startOfAlign2}[startOfAlign2 >= startOfAlign1]
    	aligned1full = ' ' * numOfSpacesToAdd1 + seq1[:startOfAlign1] + aligned1 + seq1[iOfMax:]
    	aligned2full = ' ' * numOfSpacesToAdd2 + seq2[:startOfAlign2] + aligned2 + seq2[jOfMax:]
    	alignPointer = ''
    	for cSeq1,cSeq2 in zip(aligned1, aligned2):
        	alignPointer += {True: ':', False: ' '}[cSeq1 == cSeq2]
    	alignPointerfull = ''
    	for cSeq1,cSeq2 in zip(aligned1full, aligned2full):
        	alignPointerfull += {True: ':', False: ' '}[cSeq1 == cSeq2]

    	return (aligned1, aligned2, alignPointer, aligned1full, aligned2full, alignPointerfull, [startOfAlign1, startOfAlign2])

#semi global (ends no penalty)
def computeFMatrixSemi(seq1, seq2, gap=-10): #gap=-10, large gap penalty to disable gap introduction
    	rows = len(seq1) + 1
    	cols = len(seq2) + 1
   	fMatrix = myzeros([rows, cols])
    	pointers = myzeros([rows, cols])
    	maxScore, iOfMax, jOfMax = 0, 0, 0

    	similarityMatrixMap = scorefun() 

    	for i in range(1, rows):
        	for j in range(1, cols):
            		mtch = fMatrix[i - 1][j - 1] + similarityMatrixMap[seq1[i - 1] + seq2[j - 1]]
            		delete = fMatrix[i - 1][j] + gap
            		insert = fMatrix[i][j - 1] + gap
			#don't consider 0
            		fMatrix[i][j] = max(mtch, delete, insert)
            		if(fMatrix[i][j] == delete):
                		pointers[i][j] = 1
            		elif(fMatrix[i][j] == insert):
                		pointers[i][j] = 2
            		elif(fMatrix[i][j] == mtch):
                		pointers[i][j] = 3
            		if (fMatrix[i][j] > maxScore) and ((i == rows - 1) or (j == cols - 1)) :
                		iOfMax = i
                		jOfMax = j
                		maxScore = fMatrix[i][j]
    	(aligned1, aligned2, startOfAlign1, startOfAlign2) = traceBack(pointers, seq1, seq2, gap, similarityMatrixMap, iOfMax, jOfMax)

    	'''Some formatting for displaying alignment'''
    	numOfSpacesToAdd1 = {True: 0, False: startOfAlign2 - startOfAlign1}[startOfAlign1 >= startOfAlign2]
   	numOfSpacesToAdd2 = {True: 0, False: startOfAlign1 - startOfAlign2}[startOfAlign2 >= startOfAlign1]
    	aligned1full = ' ' * numOfSpacesToAdd1 + seq1[:startOfAlign1] + aligned1 + seq1[iOfMax:]
    	aligned2full = ' ' * numOfSpacesToAdd2 + seq2[:startOfAlign2] + aligned2 + seq2[jOfMax:]
    	alignPointer = ''
    	for cSeq1,cSeq2 in zip(aligned1, aligned2):
        	alignPointer += {True: ':', False: ' '}[cSeq1 == cSeq2]
    	alignPointerfull = ''
    	for cSeq1,cSeq2 in zip(aligned1full, aligned2full):
        	alignPointerfull += {True: ':', False: ' '}[cSeq1 == cSeq2]

    	return (aligned1, aligned2, alignPointer, aligned1full, aligned2full, alignPointerfull, [startOfAlign1, startOfAlign2])

def traceBack(pointers, seq1, seq2, gap, similarityMap, i, j):
    	alignedSeq1 = ''
    	alignedSeq2 = ''
    	while (pointers[i][j] != -1) and (i > 0 and j > 0):
        	#print "i", i, "j", j, "pointer", pointers[i][j]
        	if pointers[i][j] == 1:
            		alignedSeq1 = seq1[i - 1] + alignedSeq1
            		alignedSeq2 = '-' + alignedSeq2
            		i = i - 1
        	elif pointers[i][j] == 2:
            		alignedSeq1 = '-' + alignedSeq1
            		alignedSeq2 = seq2[j - 1] + alignedSeq2
            		j = j - 1
        	elif pointers[i][j] == 3:
            		alignedSeq1 = seq1[i - 1] + alignedSeq1
            		alignedSeq2 = seq2[j - 1] + alignedSeq2
            		i = i - 1
            		j = j - 1
    	return (alignedSeq1, alignedSeq2, i, j)

def get_align2ref(aln1, aln2):
	aln_2_ref = ""
	for idx in range(len(aln1)):
		if aln1[idx] != ' ':
			if idx >= len(aln2):
				aln_2_ref += "-"
			elif aln2[idx] == ' ':
				aln_2_ref += "-"
			else:
				aln_2_ref += aln2[idx]
	return aln_2_ref

def align_test():
    	seq1 = "GTCAGGCTCTAACCGTGTTAAACGCGGCGGCAGCTGGAACAACAACGCGAACAACTGCACTGTAGGCAAACGGAATAACAACAGTCCTGACAACAGGAACAACAATCTTGGCTTCCGCTTGGCTTGTCGGC"
    	seq2 = "ATCAGGCTCTCTCCGTGTTATGCGCGGCGGCAGCTGGATCAACAACGCGAACAACTGCTCTGTAGGCGTGCGGGGCAACAGCAGTCCTGGCCTCAGGAGCAACCGTCTTGGCTTCCGCTTGGCTTGTCGGC"
    	(alignedSeq1, alignedSeq2, alignPointer, alignedSeq1Full, alignedSeq2Full, alignPointerFull) = computeFMatrix(seq1, seq2)
    	print alignedSeq1
    	print alignPointer
    	print alignedSeq2
#functions for sw-align: end

#basic functions for sequence/file processing: begin
def seqretrieve(fastafile, idlist):
	print "load seq", fastafile
        idnew = []
        if idlist:
                seq = [""] * len(idlist)
        else:
                seq = []
        infile = open(fastafile, "r")
        for aline in infile:
                aline = aline.strip()
                if not aline:
                        continue
                if aline[0] == '>':
                        subs = aline[1:].split()
                        if idlist and (subs[0] in idlist):
                                qidx = idlist.index(subs[0])
                        elif (not idlist) and (subs[0] in idnew):
                                qidx = idnew.index(subs[0])
                        else:
                                if idlist:
                                        qidx = -1
                                else:
                                        qidx = len(idnew)
                                        idnew.append(subs[0])
                                        seq.append("")
                elif qidx != -1:
                        seq[qidx] += aline
        infile.close()
	if idlist:
		print "done.. total seq", len(idlist)
	else:
		print "done.. total seq", len(idnew)
        return (idnew, seq)

def complement(seq):
	code = {"A":"T", "T":"A", "C":"G", "G":"C", "X":"X"}
	comp = ""
	for base in seq:
		comp += code[base]	
	return comp
			
def revcomplement(seq):
	code = {"A":"T", "T":"A", "C":"G", "G":"C", "a":"t", "t":"a", "c":"g", "g":"c", "X":"X", "N":"N"}
	comp = ""
	for base in seq[::-1]:
		comp += code[base]	
	return comp

def writefasta(filename, seq_id, seq_seq):
	out = open(filename, "w")
	for idx in range(len(seq_id)):
		print >>out, ">" + seq_id[idx]
		for pos in range(0, len(seq_seq[idx]), 70):
			print >>out, seq_seq[idx][pos:pos+70]
	out.close()
#basic functions for sequence/file processing: end
			
#basic functions for DGR discovery: begin
def istemplate(alignedSeq1, alignedSeq2, target='A', speedup=False, minAmut = 7, minAmutPer = 0.7):
	'''target can be 'A' for forward strand, 'T' or reverse strand'''
	involveA, totMis, other = 0, 0, 0
	involveN = 0
	for pos in range(len(alignedSeq1)):
		if (alignedSeq1[pos] == 'N') or (alignedSeq2[pos] == 'N'):
			involveN += 1
		if alignedSeq1[pos] != alignedSeq2[pos]:
			if alignedSeq1[pos] == target:
				involveA += 1
			else:
				other += 1
			totMis += 1
			if speedup and (other > 6): #speedup to check further!!!
				return False, involveA, totMis
	if (not involveN) and (involveA >= totMis * minAmutPer) and (involveA >= minAmut): #7 is used in the DiGReF pipeline
		return True, involveA, totMis
	else:
		return False, involveA, totMis

def findpair(summaryfile, candidate_id, candidate_seq, candidate_beg, candidate_end, candidate_strand, paired=True, minLen=60, description="", loose=False):
	print "total candidate", len(candidate_seq)
	tot = len(candidate_seq)
	totpair, tottr = 0, 0
	step = 1
	if paired:
		step = 2
	summary = ""
	if loose:
		global minAmutSeg, minAmutPerSeg
		mina, minp = minAmutSeg, minAmutPerSeg
	else:
		global minAmut, minAmutPer
		mina, minp = minAmut, minAmutPer

	for tr in range(0, tot, step):
		if paired: #TR and VR candidates alternate
			v1, v2 = tr + 1, tr + 2
		else: #all against all
			v1, v2 = 0, tot
		iftr = False
		for vr in range(v1, v2):
			if (vr == tr) or (vr >= tot):
				continue
			seq1, seq2 = candidate_seq[tr], candidate_seq[vr]
    			(alignedSeq1, alignedSeq2, alignPointer, alignedSeq1F, alignedSeq2F, alignPointerF, alignstart) = computeFMatrix(seq1, seq2)
			if len(alignedSeq1) < minLen:
				continue
			iftr, involveA, totMis = istemplate(alignedSeq1, alignedSeq2, minAmut=mina, minAmutPer=minp)
			iftrt, involveT, totMis = istemplate(alignedSeq1, alignedSeq2, target="T", minAmut=mina, minAmutPer=minp)
			if not (iftr or iftrt):
				continue
			iftr = True
			tr_beg, vr_beg = alignstart[0] + candidate_beg[tr] + 1, alignstart[1] + candidate_beg[vr] + 1
			tr_end, vr_end = tr_beg + len(alignedSeq1), vr_beg + len(alignedSeq2)
			if not summary:
	        		summary = open(summaryfile, "w")
				if description:
					print >>summary, description
			print >>summary, "\nTemplate: ", candidate_id[tr], "strand:", candidate_strand[tr], str(tr_beg) + "-" + str(tr_end)
			print >>summary, "Variable: ", candidate_id[vr], "strand:", candidate_strand[vr], str(vr_beg) + "-" + str(vr_end), "alnlen", len(alignedSeq1), "mutations", totMis, "involve-A", involveA, "involve-T", involveT
			if candidate_strand[vr] != candidate_strand[tr]:
				print >>summary, "(on different strand)"
			else:
				print >>summary
			print >>summary, "%-10d %s" % (tr_beg, alignedSeq1)
			print >>summary, "%10s %s" % ("", alignPointer)
			print >>summary, "%-10d %s" % (vr_beg, alignedSeq2)
			totpair += 1
		if iftr:
			tottr += 1

	print "\nTotal tr-vr pair:", totpair
	if summary:
		print >>summary, "\nTotal tr-vr pair:", totpair
		summary.close()


def findpair_seg(summaryfile, seq_id, seq_seq, ref_tr_id, ref_tr_seq, minLen=40):

	global semi
	semi = True

	print "total seq", len(seq_id), "ref-TR", len(ref_tr_id), "semi", semi
	tot = len(seq_seq)
	global minAmutSeg, minAmutPerSeg
	mina, minp = minAmutSeg, minAmutPerSeg

	totpair, tottr = 0, 0
	summary = ""
	seq_rev = []
	for aseq in seq_seq:
		seq_rev.append(revcomplement(aseq))
	seq_tag = ['u'] * tot
	seq_strand = ['+'] * tot
	seq_des = [""] * tot
	seq_aln = [""] * tot
	tr_all, tr_id = ref_tr_seq[:], ref_tr_id[:]

	#assign sequences first according to given reference TRs (if given)
	tagged_tr = 0
	tagged_tr_list = [0] * len(ref_tr_id)
	tagged_vr = 0
	for s in range(tot):
		same, same_rev, vr, vr_rev, align_2_ref = False, False, False, False, ""
		whichtr = -1 
		for tr in range(len(ref_tr_id)):
			seq1, seq2, seq2r = ref_tr_seq[tr], seq_seq[s], seq_rev[s]
    			(alignedSeq1, alignedSeq2, alignPointer, alignedSeq1F, alignedSeq2F, alignPointerF, alignstart) = align(seq1, seq2)
			align_2_ref = get_align2ref(alignedSeq1F, alignedSeq2F)
			iftr, involveA, totMis = istemplate(alignedSeq1, alignedSeq2, minAmut=mina, minAmutPer=minp)
			whichtr = tr
			if iftr:
				vr = True
				break
			if (totMis <= 2) and (len(alignedSeq1) >= minLen):
				same = True
				break
    			(alignedSeq1, alignedSeq2, alignPointer, alignedSeq1F, alignedSeq2F, alignPointerF, alignstart) = align(seq1, seq2r)
			align_2_ref = get_align2ref(alignedSeq1F, alignedSeq2F)
			iftr_rev, involveA, totMis = istemplate(alignedSeq1, alignedSeq2, minAmut=mina, minAmutPer=minp)
			if iftr_rev:
				vr_rev = True
				break
			if (totMis <= 2) and (len(alignedSeq1) >= minLen):
				same_rev = True
				break
		if same or same_rev:
			tagged_tr += 1
			tagged_tr_list[whichtr] += 1
			seq_tag[s] = 't'
			seq_des[s] = "similar to reference TR " + ref_tr_id[whichtr] + " alnlen " + str(len(alignedSeq1)) + " misMatch " + str(totMis)
			seq_aln[s] = align_2_ref
			if same_rev:
				seq_strand[s] = '-'
		elif vr or vr_rev:
			tagged_vr += 1
			seq_tag[s] = 'v'
			seq_des[s] = "pair with reference TR " + ref_tr_id[whichtr] + " mutations " + str(totMis) + " involve-A " + str(involveA)
			seq_aln[s] = align_2_ref
			if vr_rev:
				seq_strand[s] = '-'
	print "Total tagged as TR that shares similarity with reference TR:", tagged_tr, tagged_tr_list
	print "Total tagged as VR that pairs with reference TR:", tagged_vr
	print "Total tagged after comparing to reference TR:", tagged_tr + tagged_vr

	#first check if an unclassified seq is a TR, by checking if it forms TR-VR pair with another unclassified seq
	tagged_pair = 0
	for tr in range(0, tot):
		if seq_tag[tr] != 'u': continue 
		v1, v2 = tr + 1, tot
		iftr, iftr_rev = False, False
		print "check", tr, seq_id[tr] 
		for vr in range(v1, v2):
			if seq_tag[vr] != 'u': continue
			this_vr = vr
			seq1, seq1r, seq2, seq2r = seq_seq[tr], seq_rev[tr], seq_seq[vr], seq_rev[vr]
    			(alignedSeq1, alignedSeq2, alignPointer, alignedSeq1F, alignedSeq2F, alignPointerF, alignstart) = align(seq1, seq2)
			iftr, iftr_r = 0, 0
			if(len(alignedSeq1) >= minLen):
				iftr, involveA, totMis = istemplate(alignedSeq1, alignedSeq2, minAmut=mina, minAmutPer=minp)
				info_f = " forward-forward alignment len %d, involveA %d toMis %d" % (len(alignedSeq1), involveA, totMis)
			if not iftr:
		    		(alignedSeq1, alignedSeq2, alignPointer, alignedSeq1F, alignedSeq2F, alignPointerF, alignstart) = align(seq1, seq2r)
				if(len(alignedSeq1) >= minLen):
					iftr, involveA, totMis = istemplate(alignedSeq1, alignedSeq2, minAmut=mina, minAmutPer=minp)
					info_f = " forward-reverse alignment len %d, involveA %d toMis %d" % (len(alignedSeq1), involveA, totMis)

			if iftr:
				break
    			(alignedSeq1, alignedSeq2, alignPointer, alignedSeq1F, alignedSeq2F, alignPointerF, alignstart) = align(seq1r, seq2)
			if(len(alignedSeq1) >= minLen):
				iftr_rev, involveA, totMis = istemplate(alignedSeq1, alignedSeq2, minAmut=mina, minAmutPer=minp)
				info_r = " reverse-forward alignment len %d, involveA %d toMis %d" % (len(alignedSeq1), involveA, totMis)
			if not iftr_rev:
		    		(alignedSeq1, alignedSeq2, alignPointer, alignedSeq1F, alignedSeq2F, alignPointerF, alignstart) = align(seq1r, seq2r)
				if(len(alignedSeq1) >= minLen):
					iftr_rev, involveA, totMis = istemplate(alignedSeq1, alignedSeq2, minAmut=mina, minAmutPer=minp)
					info_r = " reverse-reverse alignment len %d, involveA %d toMis %d" % (len(alignedSeq1), involveA, totMis)
			if iftr_rev:
				break
		if iftr or iftr_rev:
			seq_tag[tr] = 't' 
			seq_des[tr] = " predicted TR (VR %s alnlen %d mutations %d involve-A %d)" % (seq_id[this_vr], len(alignedSeq1), totMis, involveA)
			seq_tag[this_vr] = 'v'
			seq_des[this_vr] = " predicted VR (TR %s alnlen %d mutations %d involve-A %d)" % (seq_id[tr], len(alignedSeq1), totMis, involveA)
			tagged_pair += 1
			if iftr:
				tr_id.append(seq_id[tr])
				tr_all.append(seq_seq[tr])
				print seq_id[tr], "is a TR paired with", this_vr, seq_id[this_vr], info_f
			else:
				tr_id.append(seq_id[tr] + "-reverse")
				tr_all.append(seq_rev[tr])
				seq_strand[tr] = '-'
				print seq_id[tr], "is a TR", " (reverse) paired with", this_vr, seq_id[this_vr], info_r
	print "For unclassifed seq: Total tagged as TR through pairing with a potential VR: ", tagged_pair

	#check if the remaining seq is a TR by comparing it to known TRs
	print "now tr_id", len(tr_id)
	tagged_sim = 0
	for vr in range(tot):
		if seq_tag[vr] != 'u':
			continue
		same, same_rev = False, False
		whichtr = 0
		for tr in range(len(tr_id)):
			seq1, seq2, seq2r = tr_all[tr], seq_seq[vr], seq_rev[vr]
    			(alignedSeq1, alignedSeq2, alignPointer, alignedSeq1F, alignedSeq2F, alignPointerF, alignstart) = align(seq1, seq2)
			iftr, involveA, totMis = istemplate(alignedSeq1, alignedSeq2, minAmut=mina, minAmutPer=minp)
			whichtr = tr
			if (totMis <= 2) and (len(alignedSeq1) >= minLen):
				same = True
				break
    			(alignedSeq1, alignedSeq2, alignPointer, alignedSeq1F, alignedSeq2F, alignPointerF, alignstart) = align(seq1, seq2r)
			iftr_rev, involveA, totMis = istemplate(alignedSeq1, alignedSeq2, minAmut=mina, minAmutPer=minp)
			if (totMis <= 2) and (len(alignedSeq1) >= minLen):
				same_rev = True
				break
		if same or same_rev:
			tagged_sim += 1
			seq_tag[vr] = 't'
			seq_des[vr] = "similar to detected TR"
			if same_rev:
				seq_strand[vr] = '-'
	print "Total tagged as TR that share similarity with predicted TRs : ", tagged_sim

	#then check if a seq is a variable region
	for vr in range(tot):
		if seq_tag[vr] != 'u':
			continue
		iftr, iftr_rev = False, False
		#A specific (typical type)
		for tr in range(len(tr_id)):
			seq1, seq2, seq2r = tr_all[tr], seq_seq[vr], seq_rev[vr]
    			(alignedSeq1, alignedSeq2, alignPointer, alignedSeq1F, alignedSeq2F, alignPointerF, alignstart) = align(seq1, seq2)
			iftr, involveA, totMis = istemplate(alignedSeq1, alignedSeq2, minAmut=mina, minAmutPer=minp)
			if iftr:
				break
    			(alignedSeq1, alignedSeq2, alignPointer, alignedSeq1F, alignedSeq2F, alignPointerF, alignstart) = align(seq1, seq2r)
			iftr_rev, involveA, totMis = istemplate(alignedSeq1, alignedSeq2, minAmut=mina, minAmutPer=minp)
			if iftr_rev:
				break	
		if iftr or iftr_rev:
			seq_tag[vr] = 'v'
			bcount = {'A':0, 'T':0, 'C':0, 'G':0}
			for idx in range(len(alignedSeq1)):
				if alignedSeq1[idx] != alignedSeq2[idx]:
					bcount[alignedSeq1[idx]] += 1
			seq_des[vr] = "mutations " + str(totMis) + " involve-A " + str(involveA) + " T " + str(bcount['T']) + " C " + str(bcount['C']) + " G " + str(bcount['G'])
			#seq_des[vr] = "mutations " + str(totMis) + " involve-A " + str(involveA)
			if iftr_rev:
				seq_strand[vr] = '-'
			continue
		#other specific (untypical type)??
		for tr in range(len(tr_id)):
			seq1, seq2, seq2r = tr_all[tr], seq_seq[vr], seq_rev[vr]
    			(alignedSeq1, alignedSeq2, alignPointer, alignedSeq1F, alignedSeq2F, alignPointerF, alignstart) = align(seq1, seq2)
			iftr, involveA, totMis = istemplate(alignedSeq1, alignedSeq2, minAmut=mina, minAmutPer=minp)
			if (totMis >= 10) and (len(alignedSeq1) > minLen):
				iftr = True
				break
    			(alignedSeq1, alignedSeq2, alignPointer, alignedSeq1F, alignedSeq2F, alignPointerF, alignstart) = align(seq1, seq2r)
			iftr_rev, involveA, totMis = istemplate(alignedSeq1, alignedSeq2, minAmut=mina, minAmutPer=minp)
			if (totMis >= 10) and (len(alignedSeq1) > minLen):
				iftr_rev = True
				break	
		if iftr or iftr_rev:
			seq_tag[vr] = 's'
			bcount = {'A':0, 'T':0, 'C':0, 'G':0, '-':0}
			for idx in range(len(alignedSeq1)):
				if alignedSeq1[idx] != alignedSeq2[idx]:
					bcount[alignedSeq1[idx]] += 1
			seq_des[vr] = "mutations " + str(totMis) + " involve-A " + str(involveA) + " T " + str(bcount['T']) + " C " + str(bcount['C']) + " G " + str(bcount['G'])
			if iftr_rev:
				seq_strand[vr] = '-'
			continue
				

	# write to summary file
	descode = {'t':"TR", 'v':"VR", 's':"VR-other", 'u':"unassigned"}
	count = {'t':0, 'v':0, 's':0, 'u':0}
	out = open(summaryfile, "w")
	print "Summary"
	for tr in range(tot):
		print seq_id[tr], descode[seq_tag[tr]], seq_strand[tr], seq_des[tr]
		count[seq_tag[tr]] += 1
		print >>out, seq_id[tr], descode[seq_tag[tr]], seq_strand[tr], seq_des[tr]
	print >>out, "#", count
	out.close()

	# save multiple alignments (all aligned to the reference TRs)
	tralnfile = summaryfile + "-tr.aln"
	vralnfile = summaryfile + "-vr.aln"
	traln = open(tralnfile, "w")
	vraln = open(vralnfile, "w")
	for tr in range(len(ref_tr_id)):
		print >>traln, ">" + ref_tr_id[tr]
		print >>traln, ref_tr_seq[tr]
		print >>vraln, ">" + ref_tr_id[tr]
		print >>vraln, ref_tr_seq[tr]
	for tr in range(tot):
		if seq_aln[tr] == "":
			continue
		if seq_tag[tr] == 't':
			print >>traln, ">" + seq_id[tr] + " " + seq_des[tr]		
			print >>traln, seq_aln[tr]		
		else:
			print >>vraln, ">" + seq_id[tr] + " " + seq_des[tr]		
			print >>vraln, seq_aln[tr]		
	traln.close()
	vraln.close()

#basic functions for DGR discovery: begin
def sim2candidates(fastafile, m8file, alnlencut=60, tr_as_db=True, tr_db="", outseqfile="", flank=30): #alnlencut = 60? to check YY May 23
	#get hits
	print "m8file", m8file
	infile = open(m8file, "r")
	seq_id = []
	candidate_idx, candidate_id, candidate_beg, candidate_end, candidate_seq, candidate_strand = [], [], [], [], [], []
	reftr_id, reftr_seq, reftr_2_candidate, reftr_found_in_qry = [], [], [], []
	if tr_db:
		reftr_id, reftr_seq = seqretrieve(tr_db, [])
		reftr_2_candidate = [[]] * len(reftr_id)
		reftr_found_in_qry = [False] * len(reftr_id)
	for aline in infile:
		subs = aline.split()
		if int(subs[3]) > alnlencut:
			if tr_as_db:
				query, dgr, beg, end = subs[0], subs[1], 6, 7 #6 7?
			else:
				dgr, query, beg, end = subs[0], subs[1], 8, 9
			alnide = float(subs[2])
			if reftr_id:
				dgridx = reftr_id.index(dgr)	
			if query in seq_id:
				sidx = seq_id.index(query)
			else:
				sidx = len(seq_id)
				seq_id.append(query)
			if int(subs[beg]) < int(subs[end]):
				this_beg, this_end = int(subs[beg]), int(subs[end])
			else:
				this_beg, this_end = int(subs[end]), int(subs[beg])
			this_strand = '+'
			if int(subs[8]) > int(subs[9]):
				this_strand = '-'
			#check if it overlaps with previous candidates
			overlap = False 
			this_len = this_end - this_beg + 1
			for old in range(len(candidate_id)):
				if (query != candidate_id[old]) or (this_strand != candidate_strand[old]): #check query as well, YY May 24, 2014
					continue
				old_len = candidate_end[old] - candidate_beg[old] + 1
				minbeg = min(candidate_beg[old], this_beg)
				maxend = max(candidate_end[old], this_end)
				span = maxend - minbeg + 1
				if span < this_len + old_len:
					overlap = True
					break
			#merge overlap candidates
			if overlap: #update the 'beg' and 'end' information of the previous candidate
				cidx = old
				candidate_beg[old] = minbeg
				candidate_end[old] = maxend
			else: #no overlap; add as a new candidate
				cidx = len(candidate_id)
				candidate_idx.append(sidx)
				candidate_id.append(query)
				candidate_beg.append(this_beg)
				candidate_end.append(this_end)
				candidate_strand.append(this_strand)
			if reftr_id and alnide >= 97:
				reftr_found_in_qry[dgridx] = True
			if reftr_id and (cidx not in reftr_2_candidate[dgridx]):
				reftr_2_candidate[dgridx].append(cidx) 
				
	infile.close()	
	print "candidate-seq", len(candidate_idx)

	print "retrieve seq..", seq_id
	if seq_id:
		seq_id_new, seq_seq = seqretrieve(fastafile, seq_id)
		if outseqfile:
			writefasta(outseqfile, seq_id, seq_seq)

	#get candidate seq
	candidate_seq = []
	for idx in range(len(candidate_idx)):
		print "candidate", candidate_id[idx], candidate_beg[idx], candidate_end[idx]
		candidate_seq.append("")
		if candidate_strand[idx] == '+':
			beg = max(0, candidate_beg[idx]-1-flank)
			end = candidate_end[idx] + flank
			#candidate_seq[idx] = seq_seq[candidate_idx[idx]][candidate_beg[idx]-1:candidate_end[idx]] 	
			candidate_seq[idx] = seq_seq[candidate_idx[idx]][beg:end]  #to check YY May 23	
		else: #minus
			beg = max(0, candidate_beg[idx]-flank)
			end = candidate_end[idx] + 1 + flank
			#candidate_seq[idx] = revcomplement(seq_seq[candidate_idx[idx]][candidate_beg[idx]:candidate_end[idx] + 1])
			candidate_seq[idx] = revcomplement(seq_seq[candidate_idx[idx]][beg:end]) #to check YY May 23

	#complemented with reference DGR systems
	tot_add = 0
	if tr_db:
		add_dgr_list, add_dgr_seq = [], []
		candidate_with_tr = [False] * len(candidate_id)
		for idx in range(len(reftr_id)):
			if reftr_found_in_qry[idx]:
				for c in reftr_2_candidate[idx]:
					candidate_with_tr[c] = True
		for idx in range(len(reftr_id)):
			match = 0
			for c in reftr_2_candidate[idx]:
				if not candidate_with_tr[c]:
					candidate_with_tr[c] = True
					match += 1
			if match:
				candidate_id.append(reftr_id[idx]+ " (Ref)")	
				candidate_seq.append(reftr_seq[idx])
				candidate_beg.append(0)
				candidate_end.append(len(reftr_seq[idx]))
				candidate_strand.append("+")
				tot_add += 1
	print "Reference TR added: ", tot_add
	print "Done"
	return (candidate_id, candidate_beg, candidate_end, candidate_strand, candidate_seq) 

def findseed(id, seq, brange = 10000, seedlen = 60, skip = 150, flank = 100, minaay=3, window=10, targetregion=[], aayfilter=True):
	slen = len(seq)
	candidateseq, candidatebeg, candidateend, candidateid, candidatestrand = [], [], [], [], []

	aay = [0] * seedlen
	aay_add = [0] * slen 
	rtt = [0] * seedlen
	rtt_add = [0] * slen 
	if aayfilter: #if target region is given, don't use AAC (GTT) filter
		#count AAY (AAC or AAT) for all the windows
		#RTT (GTT or ATT) (AAY on reverse complementary)
		if (seq[:3] == 'AAC') or (seq[:3] == 'AAT'): 
			aay[2], aay_add[2] = 1, 1
		elif (seq[:3] == 'GTT') or (seq[:3] == 'ATT'): 
			rtt[2], rtt_add[2] = 1, 1

		for p in range(3, seedlen):
			if (seq[p-2:p+1] == "AAC") or (seq[p-2:p+1] == 'AAT'):
				aay[p] = 1
			elif (seq[p-2:p+1] == "GTT") or (seq[p-2:p+1] == "ATT"):
				rtt[p] = 1
			aay_add[p] = aay_add[p - 3] + aay[p]
			rtt_add[p] = rtt_add[p - 3] + rtt[p]
		for p in range(seedlen, slen, 1):
			last = aay.pop(0)
			if (seq[p-2:p+1] == "AAC") or (seq[p-2:p+1] == "AAT"):
				aay.append(1)
			else:
				aay.append(0)
			aay_add[p] = aay_add[p - 3] + aay[-1] - last 
			last = rtt.pop(0)
			if (seq[p-2:p+1] == "GTT") or (seq[p-2:p+1] == "ATT"):
				rtt.append(1)
			else:
				rtt.append(0)
			rtt_add[p] = rtt_add[p - 3] + rtt[-1] - last 
		
	print "targetregion", targetregion
	if targetregion:
		p1 = targetregion[0]
		slen_scan = targetregion[1]
	else:
		p1 = 0
		slen_scan = slen
	while(p1 < slen_scan - seedlen):
		#print "seed:", p1, "-", p1 + seedlen, "aay", aay_add[p1 + seedlen]
		if aayfilter and (aay_add[p1 + seedlen] < minaay) and (rtt_add[p1 + seedlen] < minaay):
			p1 += window
			continue
		p2 = max(0, p1 - brange)
		maxp2 = min(p1 + brange, slen_scan - seedlen)
		ifseed1 = False
		while(p2 < maxp2):
			#use small minAmut for seed discovery
			tr1, involveA1, totMis1 = istemplate(seq[p1:p1+seedlen], seq[p2:p2+seedlen], target='A', speedup=True, minAmut = 3)
			if not tr1:
				tr2, involveA2, totMis2 = istemplate(seq[p1:p1+seedlen], seq[p2:p2+seedlen], target='T', speedup=True, minAmut = 3)
			if tr1 or tr2:
				candidateid.append(id)
				candidatebeg.append(max(0, p1 - flank))
				candidateend.append(min(slen_scan, p1 + seedlen + flank))
				candidateseq.append(seq[candidatebeg[-1]:candidateend[-1]])
				candidateid.append(id)
				candidatebeg.append(max(0, p2 - flank))
				candidateend.append(min(slen_scan, p2 + seedlen + flank))
				candidateseq.append(seq[candidatebeg[-1]:candidateend[-1]])
				if tr1:
					candidatestrand.extend(["+", "+"])
				else:
					candidatestrand.extend(["-", "-"])
				print "find seed", len(candidateseq)
				p2 += skip
				ifseed1 = True
			else:
				p2 += 1
		if ifseed1:
			p1 += skip
		else:
			p1 += window 

	return candidateid, candidatebeg, candidateend, candidatestrand, candidateseq

def denovo_scan(inseqfile, revtm8file="", aayfilter=True, outseqfile=""):
	seqid, seqseq = [], []
	targetregion = []
	rtgene_des = ""
	if revtm8file:
		infile = open(revtm8file, "r")
		for aline in infile:
			if aline[0] == '#':
				continue
			subs = aline.split()
			query, alnide, alnlen, alnbeg, alnend = subs[0], float(subs[2]), int(subs[3]), int(subs[6]), int(subs[7])	
			if not ((alnide >= 30) and (alnlen >= 200)):
				continue
			if query in seqid:
				continue
			seqid.append(query)
			if alnbeg > alnend:
				scan_beg, scan_end = alnend, alnbeg
			else:
				scan_beg, scan_end = alnbeg, alnend
			targetregion.append([scan_beg, scan_end])
			print "Putative RT: ", aline,
			rtgene_des += "Putative RT: " + aline
			#print "query", query, alnbeg, "alnend", alnend, "targetregion", targetregion[-1]
		infile.close()
		if seqid:
			seq_id_new, seqseq = seqretrieve(inseqfile, seqid)
			for idx in range(len(seqid)): #add flanking region
				targetregion[idx][0] = max(0, targetregion[idx][0] - 10000) 
				targetregion[idx][1] = min(len(seqseq[idx]), targetregion[idx][1] + 10000)
			if outseqfile:
				writefasta(outseqfile, seqid, seqseq)
	else: #all--slow
		seqid, seqseq = seqretrieve(inseqfile, [])
		for idx in range(len(seqid)):
			targetregion.append([])

	candidate_id, candidate_seq, candidate_beg, candidate_end, candidate_strand = [], [], [], [], []
	for idx in range(len(seqid)):
		print "Now check: ", seqid[idx], len(seqseq[idx]), seqseq[idx][:100], "..."
		id0, beg0, end0, strand0, seq0 = findseed(seqid[idx], seqseq[idx], targetregion=targetregion[idx], aayfilter=aayfilter)
		if id0:
			candidate_id.extend(id0)
       		        candidate_beg.extend(beg0)
       		        candidate_end.extend(end0)
       		        candidate_seq.extend(seq0)
       		        candidate_strand.extend(strand0)
	return candidate_id, candidate_beg, candidate_end, candidate_strand, candidate_seq, rtgene_des

def main():
	inseqfile, trm8file, revtm8file = "", "", ""
	tr_db, tr_as_db, aayfilter, seg = "", True, True, False
	outseqfile, summaryfile = "", ""
	for idx in range(len(sys.argv)):
		if sys.argv[idx] == '-inseq' and len(sys.argv) > idx + 1:
			inseqfile = sys.argv[idx + 1]
		elif sys.argv[idx] == '-tr_hom' and len(sys.argv) > idx + 1:
			trm8file = sys.argv[idx + 1]
		elif sys.argv[idx] == '-rev_hom' and len(sys.argv) > idx + 1:
			revtm8file = sys.argv[idx + 1]
		elif sys.argv[idx] == "-tr_db" and len(sys.argv) > idx + 1:
			tr_db = sys.argv[idx + 1]
		elif sys.argv[idx] == "-tr_as_qry":
			tr_as_db = False 
		elif sys.argv[idx] == "-seg":
			seg = True
		elif sys.argv[idx] == '-NOaayfilter':
			aayfilter = False
		elif sys.argv[idx] == '-saveseq' and len(sys.argv) > idx + 1:
			outseqfile = sys.argv[idx + 1]
		elif sys.argv[idx] == '-summary' and len(sys.argv) > idx + 1:
			summaryfile = sys.argv[idx + 1]
	if not (inseqfile and summaryfile):
		print "Usage: ", sys.argv[0], "-inseq seq-file -summary summary-file"
		print " Options for similarity-based:"
		print "   -tr_hom file: TR (temperate)-region-search-result-in-m8file, if defined, will narrow down the search space" 
		print " Options for denovo search:"
		print "   -rev_hom file: Reverse-transcriptase-search-result-in-m8file; constain the scanning region"
		print "   -NOaayfilter: if specified, will disable the AAC based filtering"
		print " Options for finding TR-VR pairs in sequence fragments (e.g., metagenomic sequences)"
		print "   -seg"
		print " Other options:"
		print "   -tr_db file: augument TRs with the reference TR sequences; this option is used with -seg option"
		print "   -saveseq seq-file"
		sys.exit()

	paired = True
	rtgene_des = ""
	if seg:
		seq_id, seq_seq = seqretrieve(inseqfile, [])
		reftr_id, reftr_seq = [], []
		if tr_db:
			reftr_id, reftr_seq = seqretrieve(tr_db, [])
               	findpair_seg(summaryfile, seq_id, seq_seq, reftr_id, reftr_seq)
	else:
		if trm8file:
			id, beg, end, strand, seq = sim2candidates(inseqfile, trm8file, tr_as_db=tr_as_db, tr_db=tr_db, outseqfile=outseqfile)
			paired = False
		else:
			id, beg, end, strand, seq, rtgene_des = denovo_scan(inseqfile, revtm8file=revtm8file, aayfilter=aayfilter, outseqfile=outseqfile)
			paired = True
	
       	 	if seq:
			description = "Input-seq: " + inseqfile
			if rtgene_des:
				description += "\n" + rtgene_des
                	findpair(summaryfile, id, seq, beg, end, strand, paired=paired, description=description)

if __name__ == '__main__':
	time1 = time.time()
	main()
	time2 = time.time()
	print "total time used: ", (time2 - time1)/60, "min"
