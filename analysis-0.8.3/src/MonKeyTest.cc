// vim: ts=6:sts=6:sw=6
#include <cassert>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <utility>
#include <numeric>
#include <limits>
#include <algorithm>
#include <iterator>
#include <functional>
#include <set>
#include <cctype>
#if defined (__SVR4) && defined (__sun)
#include <ieeefp.h>
#endif

#include <boost/tuple/tuple.hpp>

//libsequence stuff
#if defined(__GNUG__) && __GNUC__ >= 3
#include <Sequence/FastaExplicit.hpp>
#else
#include <Sequence/Fasta.hpp>
#include <Sequence/Alignment.hpp>
#endif
#include <Sequence/PolySites.hpp>
#include <Sequence/stateCounter.hpp>
#include <Sequence/shortestPath.hpp>

//base routines for dealing with coding sequence data
#include <polydNdSbase.hpp>
#include <codingRegionProcessor.hpp>
#include <MKtestOutput.hpp>

using namespace std;
using namespace Sequence;
using namespace Alignment;

using boost::tuple;
using boost::get;
using boost::make_tuple;

enum {FIXEDA=0,FIXEDS,POLYA,POLYS};

const double MAXDBL = numeric_limits<double>::max();


struct isPresent : public std::binary_function <
  iterator_traits<Sequence::PolySites::const_site_iterator>::value_type,
  double,
  bool >
{
  typedef iterator_traits<Sequence::PolyTable::const_site_iterator>::value_type vt;
  inline bool operator()(const vt & v, const double &d) const
  {
    return (v.first == d);
  }
};

vector<string> makePopCodons(const int codon_pos,
			     const int site_index,
			     const vector<Fasta> & alignment,
			     const vector<int> & intervals);

set<string> makeCodSet(vector<string>::const_iterator beg,
		       vector<string>::const_iterator end);

set<shortestPath::pathType> catalogueDiffs(const set<string> &cod1,
					   const set<string> &cod2,
					   const int & cPos);

set<shortestPath::pathType> catalogueDiffsComplex(const set<string> &cod1,
						  const set<string> &cod2,
						  const int & cPos);
unsigned countFixations(const PolySites & table,
			const size_t & n1);

size_t IncrementSNPitr(PolySites::const_site_iterator beg,
		       PolySites::const_site_iterator end,
		       const unsigned &n,
		       const vector<int> &intervals);

void updateFixCells( const shortestPath::pathType & type,
		     unsigned * fixedA,
		     unsigned * fixedS);

std::tuple<double,double,double> getPositions(PolySites::const_site_iterator beg,
						PolySites::const_site_iterator end,
						const unsigned &n,
						const vector<int> &intervals);

void usage(void);

void (*usg)(void) = &usage;

void print_coding_warnings(const codingRegionProcessor& crp) {
      Warnings w = crp.warnings();
      //cerr << "Warnings:";
      std::copy(w.begin(),w.end(),
		std::ostream_iterator<const string>(std::cerr,"\n"));
}

vector<Fasta> expand_multiseq(const vector<Fasta>& multiseqs) {
      // In this modification to MKtest, some Fasta sequences are
      // implicitly repeated a certain number of times. This function
      // makes that repitition explicit by turning multiplicity>1 seqs
      // into (multiplicity) sequences of multiplicity=1.
      vector<Fasta> ret;
      for (Fasta mseq : multiseqs) {
            // Don't want to mess with libsequence's "copy constructor"
            // so we do it manually :/
            Fasta m1seq(mseq.GetName(), mseq.GetSeq());
            m1seq.clone_id = mseq.clone_id;
            for (int i=0; i<mseq.multiplicity; i++) {
                  ret.push_back(m1seq);
            }
      }
      return ret;
}

int main (int argc, char *argv[])
{
  params args;
  //try{/
      parseargs(argc,argv,&args,usg);

      vector<Fasta> in_seqs, out_seqs;  // Sequences from ingroup, outgroup
      ifstream igf(args.i1);  // In Group File
      ifstream ogf(args.i2);  // Out Group File
      if (!igf) {cerr<<"Couldn't open "<<args.i1<<"\n"; exit(10);}
      if (!ogf) {cerr<<"Couldn't open "<<args.i2<<"\n"; exit(10);}
      while(!igf.eof()) {Fasta tmp; tmp.read_multiseq(igf); in_seqs.push_back(tmp);}
      while(!ogf.eof()) {Fasta tmp; tmp.read_multiseq(ogf); out_seqs.push_back(tmp);}

      map< int, vector<Fasta> >  in_cid2ms;  // clone id to array of seqs
      map< int, vector<Fasta> >  out_cid2ms;  // clone id to array of seqs

      // Bin the sequences based on clone id
      for (auto in_seq : in_seqs) {
            auto clone_it = in_cid2ms.find(in_seq.clone_id);
            if (clone_it==in_cid2ms.end())
                in_cid2ms[in_seq.clone_id] = vector<Fasta>();
            if (in_seq.length())
                  in_cid2ms[in_seq.clone_id].push_back(in_seq);
      }
      for (auto out_seq : out_seqs) {
            auto clone_it = out_cid2ms.find(out_seq.clone_id);
            if (clone_it==out_cid2ms.end())
                out_cid2ms[out_seq.clone_id] = vector<Fasta>();
            if (out_seq.length())
                  out_cid2ms[out_seq.clone_id].push_back(out_seq);
      }

      set<int> in_cids, out_cids;
      for (auto kv : in_cid2ms)
            in_cids.insert(kv.first);
      for (auto kv : out_cid2ms)
            out_cids.insert(kv.first);
      // Find clone IDs (cids) common to infile and outfile
      vector<int> common_cids;
      std::set_intersection(in_cids.begin(), in_cids.end(),
                            out_cids.begin(), out_cids.end(),
                            std::back_inserter(common_cids));

      // Ensure that all sequences in each bin have the same length
      for (int cid : common_cids) {
            bool all_good = true;
            unsigned common_len = in_cid2ms[cid][0].length();
            for (Fasta& seq : in_cid2ms[cid])
                  if (seq.length()!=common_len) all_good = false;
            for (Fasta& seq : out_cid2ms[cid])
                  if (seq.length()!=common_len) all_good = false;
            if (!all_good)
                  cerr<<"Clone "<<cid<<" (in /tmp/a.ffn) had uneven lengths!\n";
            if (!all_good || cid==0) {
                  ofstream fout("/tmp/a.ffn");
                  vector<Fasta>& ingroup = in_cid2ms[cid];
                  vector<Fasta>& outgroup = out_cid2ms[cid];
                  for (int i=0; i<ingroup.size(); i++) {
                        fout<<">in"<<i<<" "<<ingroup[i].length()<<"\n";
                        fout<<ingroup[i].GetSeq()<<"\n";
                  }
                  for (int i=0; i<outgroup.size(); i++) {
                        fout<<">out"<<i<<" "<<outgroup[i].length()<<"\n";
                        fout<<outgroup[i].GetSeq()<<"\n";
                  }
            }
            if (!all_good) exit(1);
      }

      // Count mutations according to the multiplicity of the underlying strand
      bool should_print_header = true;  // Set to false after first iteration
      for (int cid : common_cids) {
            // Count ingroup SNPs per aligned site
            vector<Fasta>& ingroup = in_cid2ms[cid];
            int m=ingroup.size(), n=ingroup[0].length();  // m=#rows, n=#cols
            vector<stateCounter> inctr(n);
            for (int i=0; i<m; i++) for (int j=0; j<n; j++) {
                  char nucleotide = ingroup[i][j];
                  inctr[j].count(nucleotide, ingroup[i].multiplicity);
            }

            // Count outgroup SNPs per aligned site
            vector<Fasta>& outgroup = out_cid2ms[cid];
            m = outgroup.size(), n=outgroup[0].length();
            vector<stateCounter> outctr(n);
            for (int i=0; i<m; i++) for (int j=0; j<n; j++) {
                  char nucleotide = outgroup[i][j];
                  outctr[j].count(nucleotide, outgroup[i].multiplicity);
            }

            // Find sites with simple SNPs (no gaps, only 2 different letters)
            vector<double> in_ss_locs, out_ss_locs;
            for (int j=0; j<n; j++) {
                  bool notgap = inctr[j].gap==0 && outctr[j].gap==0;
                  if (notgap && inctr[j].nStates()==2)
                        in_ss_locs.push_back(j+1);
                  if (notgap && outctr[j].nStates()==2)
                        out_ss_locs.push_back(j+1);
            }

            // Hack: expand sequences with multiplicity>1 into seqs with mult=1
            vector<Fasta> gene1 = expand_multiseq(ingroup);
            vector<Fasta> gene2 = expand_multiseq(outgroup);
            vector<Fasta> data = ingroup;
            data.insert(data.end(), outgroup.begin(), outgroup.end());
            args.n1 = ingroup.size();  // Yuck!

            // Call mutations as synonymous or nonsynonymous
            vector<int> fullseq = {0,n-1}; args.intervals = fullseq;
            codingRegionProcessor in_crp(ingroup, in_ss_locs, fullseq);
            codingRegionProcessor out_crp(ingroup, out_ss_locs, fullseq);
            print_coding_warnings(in_crp);
            print_coding_warnings(out_crp);
            PolySites *A1 = new PolySites(in_crp.replacementTable());
            PolySites *S1 = new PolySites(in_crp.synonymousTable());
            PolySites *A2 = new PolySites(out_crp.replacementTable());
            PolySites *S2 = new PolySites(out_crp.synonymousTable());
            //the following is necessary to handle sites polymorphic in
            //both genes/species
            set<double> uniqueA,uniqueS;
            uniqueA.insert(A1->pbegin(),A1->pend());
            uniqueA.insert(A2->pbegin(),A2->pend());
            uniqueS.insert(S1->pbegin(),S1->pend());
            uniqueS.insert(S2->pbegin(),S2->pend());

            //Process divergence
            vector<unsigned> contingency_table(4,0);   //store cell counts
            // ^ indexed by FIXEDA FIXEDS POLYA POLYS
            set<double> complex;  // Really?
            PolySites *allSNP = new PolySites(data);
            // Iterate over segregating sites, each of which looks like
            // (double position, string seqNumToCharMap)
            PolySites::const_site_iterator sbeg = allSNP->sbegin(),
            send = allSNP->send();
            contingency_table[POLYA] = uniqueA.size(); //repl Poly
            contingency_table[POLYS] = uniqueS.size(); //syn poly
            unsigned tot_fixations = 0;
            while (sbeg < send)
            {
              int indexed_pos = int(sbeg->first)-1;
              unsigned nPolyPerCodon = 0;
              if (InCoding(indexed_pos,args.intervals)==true)
                {
                  // need to assign a temp
                  // b/c erase-remove bit below
                  // doesn't work w/std::set
                  string temp(sbeg->second.begin(),
                          sbeg->second.begin()+args.n1);
                  //remove missing data
                  temp.erase( remove(temp.begin(),temp.end(),'N'),temp.end() );
                  set<char> one(temp.begin(),temp.end()); // uniq nucs in ingrp

                  //same for 2nd aligment partition
                  temp = string(sbeg->second.begin()+args.n1,
                            sbeg->second.end());
                  temp.erase( remove(temp.begin(),temp.end(),'N'),temp.end() );
                  set<char> two(temp.begin(),temp.end()); // uniq nucs in outgrp

                  //calculate overlap b/w partitions
                  vector<char> overlap(one.size()+two.size());

                  vector<char>::iterator itr = set_intersection(one.begin(),one.end(),
                                                    two.begin(),two.end(),
                                                    overlap.begin());


                  if ( itr-overlap.begin() == 0 && one.size()==1 && two.size()==1) //fixed diff
                  {
                  int cPos = GetCodonPos(indexed_pos,args.intervals);//codon pos of mut
                    //The set overlap is of size 0, therefore
                    //no states are shared b/w the two partitions.
                    //This is a fixed difference.
                    vector<string> codons = makePopCodons(cPos,indexed_pos,
                                                data,args.intervals);
                    PolySites codPol(codons);
                    nPolyPerCodon = codPol.numsites();
                    //number of fixations at this codon
                    unsigned nfixations = countFixations(codPol,args.n1);
                    tot_fixations += nfixations;
                    set<string> cod1 = makeCodSet(codons.begin(),
                                          codons.begin()+args.n1);
                    set<string> cod2 = makeCodSet(codons.begin()+args.n1,
                                          codons.end());

                    if (nfixations == 1)
                      {
                        //this is the easiest case
                        set<shortestPath::pathType> types = catalogueDiffs(cod1,cod2,cPos);
                        if(types.size()==1)
                        {
                          switch ( *(types.begin()) )
                            {
                            case shortestPath::S :
                              contingency_table[FIXEDS]++;
                              break;
                            case shortestPath::N :
                              contingency_table[FIXEDA]++;
                              break;
                            default:
                              complex.insert(sbeg->first);
                              break;
                            }
                        }
                        else
                        {
                          complex.insert(sbeg->first);
                        }
                      }
                    else if (nfixations > 1)
                      {
                        //harder case
                        set<shortestPath::pathType> types = catalogueDiffsComplex(cod1,cod2,cPos);
                        if (types.size()==1)
                        {
                          updateFixCells( *(types.begin()),
                                      &contingency_table[FIXEDA],
                                      &contingency_table[FIXEDS]);
                        }
                        else
                        {
                          //fixation occurs at a codon w/complex history,
                          //so we save the positions in the set "complex"
                              std::tuple<double,double,double> compPos =
                            getPositions(sbeg,send,nPolyPerCodon,args.intervals);
                          double x = get<0>(compPos);
                          if ( x != MAXDBL)
                            {
                              complex.insert(x);
                            }
                          x = get<1>(compPos);
                          if ( x != MAXDBL)
                            {
                              complex.insert(x);
                            }
                          x = get<2>(compPos);
                          if ( x != MAXDBL)
                            {
                              complex.insert(x);
                            }
                        }
                      }
                  }
                  else
                  {
                    //states are shared, therefore no fixation
                  }
                }
              sbeg += IncrementSNPitr(sbeg,send,nPolyPerCodon,
                                args.intervals);
            }
            if (should_print_header) {
                  should_print_header = false;
                  cout<<"clone_id\tdN\tdS\tpN\tpS"<<endl;
            }
            vector<unsigned>& res = contingency_table;
            cout<<cid<<"\t";
            cout<<res[FIXEDA]<<"\t"<<res[FIXEDS]<<"\t";
            cout<<res[POLYA]<<"\t"<<res[POLYS]<<endl;
      }
/*
    }
  catch (Sequence::SeqException &e)
    {
      cerr << "Processing file "
	   << args.infile
	   << "resulting in the following exception begin thrown:\n";
      cerr << e << endl;
      exit(10);
    }
  catch (std::exception &e)
    {
      cerr << "Processing file "
	   << args.infile
	   << "resulting in the following exception begin thrown:\n";
      cerr << e.what() << endl;
      exit(10);
    }*/
}

vector<string> makePopCodons(const int codon_pos,
			     const int site_index,
			     const vector<Fasta> & alignment,
			     const vector<int> & intervals)
/*!
  Reduces a data set to just the codon of interest.  Returns
  the codon in a vector of strings
*/
{
  vector<string> rv;

  for(unsigned i=0;i<alignment.size();++i)
    {
      string codon;
      MakeCodon(site_index,i,codon_pos,alignment,intervals,codon);
      rv.push_back(codon);
    }
  return rv;
}

set<string> makeCodSet(vector<string>::const_iterator beg,
		       vector<string>::const_iterator end)
/*!
  Makes a std::set of the strings in the range beg,end.
  Useful if you just want the unique codons in a range
*/
{
  set<string>rv;
  while(beg<end)
    {
      //skip codons w/missing data
      if (beg->find('N')==string::npos)
	rv.insert(*beg);
      ++beg;
    }
  return rv;
}

set<shortestPath::pathType> catalogueDiffs(const set<string> &cod1,
					   const set<string> &cod2,
					   const int & cPos)
/*!
  Calculates the set of shortestPath::pathType differences at
  codon position cPos between all codons in cod1 and all in cod2.
*/
{
  typedef boost::tuples::tuple<shortestPath::pathType,
    shortestPath::pathType,
    shortestPath::pathType> codonTuple;
  set<shortestPath::pathType> types;
  set<string>::const_iterator b1,b2;
  for(b1=cod1.begin() ; b1 != cod1.end() ; ++b1)
    {
      for(b2=cod2.begin() ; b2 != cod2.end() ; ++b2)
	{
	  codonTuple t = diffTypeMulti(*b1,*b2);
	  shortestPath::pathType type;

	  //retrieve type at the site of substitution
	  if (cPos == 0)
	    {
	      type = get<0>(t);
	    }
	  else if (cPos == 1)
	    {
	      type = get<1>(t);
	    }
	  else if (cPos == 2)
	    {
	      type = get<2>(t);
	    }
	  types.insert(type);
	}
    }
  return types;
}

set<shortestPath::pathType> catalogueDiffsComplex(const set<string> &cod1,
						  const set<string> &cod2,
						  const int & cPos)
{
  set<shortestPath::pathType> types;
  set<string>::const_iterator b1,b2;
  for(b1=cod1.begin() ; b1 != cod1.end() ; ++b1)
    {
      for(b2=cod2.begin() ; b2 != cod2.end() ; ++b2)
	{
	  shortestPath sp(*b1,*b2);
	  shortestPath::pathType type = sp.type();
 	  types.insert(type);
	}
    }
  return types;
}

unsigned countFixations(const PolySites & table,const size_t &n1)
/*!
  Counts all the fixations in a snp table
*/
{
  PolySites::const_site_iterator beg = table.sbegin(),
    end=table.send();
  unsigned nfix=0;
  while(beg<end)
    {
      set<char> _one(beg->second.begin(),
		     beg->second.begin()+n1);
      set<char> _two(beg->second.begin()+n1,
		     beg->second.end());
      vector<char> overlap(_one.size()+_two.size());
      vector<char>::iterator itr = set_intersection(_one.begin(),_one.end(),
						    _two.begin(),_two.end(),
						    overlap.begin());
      if (itr - overlap.begin() == 0)
	++nfix;
      ++beg;
    }
  return nfix;
}

void updateFixCells( const shortestPath::pathType & type,
		     unsigned * fixedA,
		     unsigned * fixedS )
/*!
  Updates cell entries in MK tables.
*/
{
  switch (type)
    {
    case shortestPath::S :
      ++(*fixedS);
      break;
    case shortestPath::N :
      ++(*fixedA);
      break;
    case shortestPath::SN :
      (*fixedS)++;
      (*fixedA)++;
      break;
    case shortestPath::SS :
      (*fixedS)+=2;
      break;
    case shortestPath::NN :
      (*fixedA)+=2;

      break;
    case shortestPath::SSS :
            (*fixedS)+=3;

      break;
    case shortestPath::SSN :
            (*fixedS)+=2;

            (*fixedA)++;
      break;
    case shortestPath::SNN :
      (*fixedS)++;
            (*fixedA)+=2;
      break;
    case shortestPath::NNN :
            (*fixedA)+=3;
	    //(*fixedA)+=1;
      break;
    default:
      break;
    }
}

size_t IncrementSNPitr(PolySites::const_site_iterator beg,
		       PolySites::const_site_iterator end,
		       const unsigned &n,
		       const vector<int> &intervals)
/*!
  Figures out how much to increment the pointer to site columns.  Used in
  processing fixed differences.
*/
{
  size_t k=1;
  while(beg<end && k<n)
    {
      int index = int(beg->first)-1;
      if (InCoding(index,intervals))
	{
	  ++k;
	}
      ++beg;
    }
  return k;
}

std::tuple<double,double,double> getPositions(PolySites::const_site_iterator beg,
						PolySites::const_site_iterator end,
						const unsigned &n,
						const vector<int> &intervals)
{
  double rv[3];
  rv[0]=rv[1]=rv[2]=MAXDBL;
  unsigned k = 0;
  while(beg<end && k<n)
    {
      int index = int(beg->first)-1;
      if (InCoding(index,intervals))
	{
	  rv[k] = beg->first;
	  ++k;
	}
      ++beg;
    }
  return make_tuple(rv[0],rv[1],rv[2]);
}

void usage(void)
{
  cerr << "MKtest -i <infile> -I nint i j k l ..."<<endl;
  cerr << "see man MKtest for details"<<endl;
}
