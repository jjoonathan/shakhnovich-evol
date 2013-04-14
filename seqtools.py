import re, code, itertools, unittest
from math import *

def read_fasta(infile):
	""" Input is an iterator that returns lines of the fasta file.
	    Returns (headers, seqs). """
	if type(infile)==str:
		infile = open(infile)
	out_headers, out_prot_seqs = [], []
	curheader = None
	seqlines = []
	for line in infile:
		if line.startswith('>'):
			if seqlines:
				out_headers.append(curheader); curheader = None
				out_prot_seqs.append(''.join(seqlines)); seqlines = []
			curheader = line
		else:
			seqlines.append(line.strip())
	if curheader != None:
		out_headers.append(curheader)
		out_prot_seqs.append(''.join(seqlines))
	return out_headers, out_prot_seqs

def write_fasta((headers, seqs), outfile):
	if type(outfile)==str:
		outfile = open(outfile,'w')
	for hdr,seq in zip(headers,seqs):
		outfile.write(hdr.strip())
		outfile.write("\n")
		outfile.write(seq.strip())
		outfile.write("\n")

def viterbi(obs, states, pStartInState_, pXitionFrom_to_, pState_emits_):
	""" Credit to the wiki page for this especially clean formulation.
			obs: the array of observations (comparables)
			states: an array of all possible states (comparables)
			pStartInState_: a dictionary mapping state -> P(start in state)
			pXitionFrom_to_[state1][state2] -> P(state1 xitions to state2)
			pState_emits_[state][observation] -> P(state emits observation)
		Returns: (log(P(path)),[path[0],path[1],...]) where path is the maximum likelihood path"""
	for s in states:
		pStartInState_[s] = log(pStartInState_[s])
	for s1,s2 in itertools.product(states,states):
		pXitionFrom_to_[s1][s2] = log(pXitionFrom_to_[s1][s2])
	for s in states:
		for o in pState_emits_[s].keys():
			pState_emits_[s][o] = log(pState_emits_[s][o])
	pAtTime_inState_ = [{}]
	pathGivenFinalState_ = {}
	for s in states:
		pAtTime_inState_[0][s] = pStartInState_[s] * pState_emits_[s][obs[0]]
		pathGivenFinalState_[s] = [s]
	for t in range(1,len(obs)):
		pAtTime_inState_.append({})
		newpath = {}
		for cur_s in states:
			(p,prev_s) = max([(pAtTime_inState_[t-1][prev_s]+\
								pXitionFrom_to_[prev_s][cur_s]+\
								pState_emits_[cur_s][obs[t]],\
							   prev_s) for prev_s in states])
			pAtTime_inState_[t][cur_s] = p
			newpath[cur_s] = pathGivenFinalState_[prev_s]+[cur_s]
		pathGivenFinalState_ = newpath
	(p,final_s) = max([ (pAtTime_inState_[-1][final_s],final_s) for final_s in states])
	return (p, pathGivenFinalState_[final_s])

def fasta_viterbi(fasta_file_path):
	""" Appends an observation and max-likelihood sequence to fasta_file_path,
	puts it in fasta_file_path+'.v' """
	headers, seqs = read_fasta(open(fasta_file_path))
	observations, maximum_likelihood_path = viterbify_seqs(seqs)
	of = open(fasta_file_path+'.v',)
	of.write(open(fasta_file_path).read())
	of.write('>obs\n'); of.write(observations)
	of.write('\n>mlp\n'); of.write(maximum_likelihood_path)
	of.write('\n')
	of.close()

def viterbify_seqs(seqs):
	""" Input:
	       seqs: a list of protein sequences
	    Returns: [
	       observations: a string classifying each nucleotide as G=gap, D=disagree, A=agree
	       maximum_likelihood_path: a string of 'I', 'H' (indel / homomorphic)
		]"""
	fraction_identity = []
	if len(seqs)==0:
		print "No sequences to viterbify!"
		return "", ""
	out_seq = seqs[-1]  # Sequence of the outgroup
	best_match, best_match_fi = 0, 0
	for i in range(len(seqs)-1):
		in_seq = seqs[i]
		num_identical = 0
		for j in range(len(in_seq)):
			if in_seq[j] == out_seq[j]:
				num_identical += 1
		fi = float(num_identical)/len(in_seq)
		if fi >= best_match_fi:
			best_match, best_match_fi = i, fi
		fraction_identity.append(fi)
	# ---------------------------------------------
	in_seq = seqs[best_match]  # In-group sequence
	#out_seq = outgroup sequence (defined above)
	obs_gap, obs_disagree, obs_agree = 'GDA'  # Output state of HMM
	observations = []
	for i in range(len(in_seq)):
		AAs = {}
		has_gap = False
		for j in range(len(seqs)):
			aa = seqs[j][i]
			if aa == '-':
				has_gap = True
				break
			AAs[aa] = AAs.get(aa,0)+1
		if has_gap:
			observations.append(obs_gap)
		elif len(AAs) > 4:  # We could use entropy calculations if this proves to be a bad measure of disagreement
			observations.append(obs_disagree)
		else:
			observations.append(obs_agree)
	observations = ''.join(observations)
	state_indel, state_hom   = 0,1
	all_states = [state_indel, state_hom]
	pState_emits_         = {state_indel:{obs_gap:.1,  obs_disagree:.65, obs_agree:.25},
							 state_hom:  {obs_gap:.01, obs_disagree:.04, obs_agree:.95}}
	pXitionFrom_to_       = {state_indel: {state_indel:.99, state_hom:.01   },
							 state_hom:   {state_indel:.01, state_hom:.99 }}
	pStartInState_        = {state_indel:.3,  state_hom:.7}
	p, path = viterbi(observations, all_states, pStartInState_, pXitionFrom_to_, pState_emits_)
	maximum_likelihood_path = ''.join('IH'[i] for i in path)
	return observations, maximum_likelihood_path

def ranges_for_indel_calls(maximum_likelihood_path):
	""" maximum_likelihood_path is a string like 'HHHHHHHHHIIIHHHHIIIIIH'
	    containing the max-likelihood calls (I=indel, H=homologous).
	    ranges_for_indel_calls then identifies blocks like 'IIII' and
	    returns a list of their pythonic ranges. """
	# max_p_calls is a string of the same length as prot_seqs containing:
	max_p_calls = maximum_likelihood_path
	xitions = []  # A list of transitions from state I to H (even indexes) and H to I (odd indexes)
	if max_p_calls[0] == 'H':
		xitions.append(0)  # 0th transition is always I to H
		last_state = 'H'
	else:
		last_state = 'I'
	for i in range(1,len(max_p_calls)):
		if last_state != max_p_calls[i]:
			xitions.append(i)
			last_state = max_p_calls[i]
	if last_state == 'H':
		xitions.append(len(max_p_calls))
	homologous_ranges = [(xitions[i],xitions[i+1]) for i in range(0,len(xitions),2)]
	# --------- Check to see if we got it correct
	#max_p_calls_from_ranges = bytearray('I'*len(max_p_calls))
	#for (start, end) in homologous_ranges:
	#	for i in range(start, end):
	#		max_p_calls_from_ranges[i] = 'H'
	#assert str(max_p_calls_from_ranges) == max_p_calls
	# --------- End verification
	return homologous_ranges

default_prot_hdr_filt = lambda s: re.split('[:_]', s.strip('> \n').split(' ')[0])
default_nuc_hdr_filt = lambda s: re.split('[:_]', s.strip('> \n').split(' ')[0])

def pal2nal(aligned_scrambled_prot_tup, unaligned_nuc_tup, prot_hdr_filt=None, nuc_hdr_filt=None):
	""" Takes (hdrs,seqs) tuples for an aligned protein FASTA file and an unaligned nucleotide FASTA file.
	[prot,nuc]_hdr_filt are applied to headers to generate python objects whose equality is used to match
	protein and nucleotide sequences. Then, gaps are inserted into the nucleotide sequences and they are returned."""
	aligned_scrambled_prot_hdrs,aligned_scrambled_prot_seqs = aligned_scrambled_prot_tup
	unaligned_nuc_hdrs,unaligned_nuc_seqs = unaligned_nuc_tup
	if prot_hdr_filt == None:
		prot_hdr_filt = default_prot_hdr_filt
	if nuc_hdr_filt == None:
		nuc_hdr_filt = default_nuc_hdr_filt
	nuc_seq_ids = [nuc_hdr_filt(hdr) for hdr in unaligned_nuc_hdrs]
	prot_seq_ids = [prot_hdr_filt(hdr) for hdr in aligned_scrambled_prot_hdrs]
	try:
		idxs = [prot_seq_ids.index(nuc_id) for nuc_id in nuc_seq_ids]
	except ValueError as e:
		code.interact(None,None,locals())
	aligned_prot_seqs = [aligned_scrambled_prot_seqs[i] for i in idxs]
	aligned_prot_hdrs = [aligned_scrambled_prot_hdrs[i] for i in idxs]
	aligned_nuc_seqs = []
	for i in range(len(unaligned_nuc_seqs)):
		unaligned_nuc_seq = unaligned_nuc_seqs[i]
		aligned_prot_seq = aligned_prot_seqs[i]
		aligned_nuc_seq = bytearray(3*len(aligned_prot_seq))
		aligned_nuc_idx = 0
		unaligned_nuc_idx = 0
		for j in range(len(aligned_prot_seq)):
			if aligned_prot_seq[j] != '-':
				unaligned_codon = bytearray(unaligned_nuc_seq[unaligned_nuc_idx:unaligned_nuc_idx+3],'utf8')
				aligned_nuc_seq[(j*3):(j*3+3)] = unaligned_codon
				unaligned_nuc_idx += 3
			else:
				aligned_nuc_seq[(j*3):(j*3+3)] = '---'
		aligned_nuc_seqs.append(aligned_nuc_seq)
	aligned_nuc_seqs = [str(s) for s in aligned_nuc_seqs]
	aligned_nuc_tup = (unaligned_nuc_hdrs, aligned_nuc_seqs)
	unscrambled_prot_tup = (aligned_prot_hdrs, aligned_prot_seqs)
	return unscrambled_prot_tup, aligned_nuc_tup

def chop_indels(aligned_prot_tup, aligned_nuc_tup):
	""" Given an aligned protein FASTA file handle (open()'d file) and an unaligned
	nucleic acid FASTA file handle, chops out anything viterbi calls as an indel.
	Returns: ((headers, seq),good_prot_ranges) of the aligned nucleic acid sequences or False."""
	if len(aligned_prot_tup[0])==0:
		print "Can't chop indels -- no seqs!"
		return False, []
	(prot_hdrs,prot_seqs),(nuc_hdrs,nuc_seqs) = aligned_prot_tup, aligned_nuc_tup
	observations, max_p_calls = viterbify_seqs(prot_seqs)
	chopped_nuc_seqs = []
	homologous_ranges = ranges_for_indel_calls(max_p_calls)
	for i in range(len(nuc_seqs)):
		chopped_nuc_seq = ''.join(str(nuc_seqs[i][start*3:end*3]) for start,end in homologous_ranges)
		chopped_nuc_seqs.append(chopped_nuc_seq)
	num_sites_eq_to_outgrp = sum(1 if chopped_nuc_seqs[-1][i]==chopped_nuc_seqs[0][i] else 0 for i in range(len(chopped_nuc_seqs[0])))
	if len(chopped_nuc_seqs[0])==0:
		print "Entirely chopped."
		return False, []
	percent_identity = num_sites_eq_to_outgrp*1.0/len(chopped_nuc_seqs[0])
	if percent_identity<.35:
		print "Seq 0 is too dissimilar to the outgroup, aborting chop."
		return False, []
	return (nuc_hdrs, chopped_nuc_seqs), homologous_ranges

def drop_bad_seqs((headers,seqs), threshold=.5):
	""" Drops sequences from the input FASTA tuple that have
	a percent identity to the last sequence (outgroup) of less than
	threshold. Returns a FASTA tuple (headers,seqs)"""
	outgrp = seqs[-1]
	bad_seq_ids = []
	def should_keep_seq((header, seq, idx)):
		num_sites_eq_to_outgrp = sum(int(seq[i]==outgrp[i]) for i in range(len(outgrp)))
		percent_identity = num_sites_eq_to_outgrp*1.0/len(outgrp)
		if percent_identity>=threshold:
			return True
		else:
			bad_seq_ids.append(idx)
			return False
	headers, seqs, idxs = zip(*filter(should_keep_seq,zip(headers,seqs,itertools.count())))
	return (headers, seqs), bad_seq_ids

def num_seqs_in_file(file_path):
	n=0
	for line in open(file_path):
		n += int(line.startswith('>'))
	return n

def idx2aln_map(aln):
	idx_map = []
	for i in range(len(aln)):
		if aln[i]!='-':
			idx_map.append(i)
	idx_map.append(len(aln))
	return idx_map

def map_idxs_from_seq_to_aln(idxs, seq, aln):  # Compatibility
	idx_map = idx2aln_map(aln)
	return [idx_map[i] for i in idxs]

class Test(unittest.TestCase):
	def test_idx_map(self):
		sa = "ABC---D-E--FG"
		su = "ABCDEFG"
		for i0,i1 in itertools.combinations(range(len(su)+1),2):
			if i0 > i1:
				i0,i1 = i1,i0
			m0,m1 = map_idxs_from_seq_to_aln([i0,i1],su,sa)
			su_cut = su[i0:i1]
			sa_cut = ''.join([c for c in sa[m0:m1] if c!='-'])
			self.assertEqual(su_cut,sa_cut)
