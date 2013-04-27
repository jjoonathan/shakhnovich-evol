import seqtools as st
from shared import *

eh,es = st.read_fasta('eugene_seqs.faa')
eh = [h.strip()[1:].split(' ')[0] for h in eh]
eh2s = {}
reh,res = [],[]
beh,bes = [],[]
for h,s in zip(eh,es):
    if h.endswith('_RED'):
        reh.append(h[:-4])
        res.append(s)
        eh2s[h[:-4]]=s
    elif h.endswith('_BLK'):
        eh2s[h[:-4]]=s
        beh.append(h[:-4])
        bes.append(s)
    else:
        print "Couldn't classify seq "+repr(h)

print 'Red symmetric difference: ' + repr(list(set(reh).symmetric_difference(set(rah))))
print 'Black symmetric difference: ' + repr(list(set(beh).symmetric_difference(set(bah))))

for h in rah+bah:
    try:
        aln_seq = ''.join([c for c in asmap[h] if c!='-'])
        e_seq = ''.join([c for c in eh2s[h] if c!='-'])
    except KeyError:
        print "Couldn't find eugene seq corresponding to "+h
        continue
    if aln_seq!=e_seq:
        print 'Mismatched sequences for '+h+' (written to mismatch.faa)'
        mm=open('mismatch.faa','w')
        mm.write('>aln_seq '+h+'\n'+aln_seq+'\n')
        mm.write('>e_seq '+h+'\n'+e_seq+'\n')
