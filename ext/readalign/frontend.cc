#include "readalign.h"

#include <string.h>
#include <stdio.h>

// defined in alignutils.c
extern "C"
void findband( char* iSeq1, char* iSeq2, int iLen1, int iLen2, char* pQuals1, int iInsNucPrior, int iGapOpen, int iGapExtend, int band, 
	       int global_maxband, int global_minband, int* yposmin, int* iNewlen2, int* newWidth );


double posterior( const string state, int i, int j, AlignDPTable* pFW, AlignDPTable* pBW ) {
  logspace total = pBW->getProb("start",0,0);
  logspace fw = pFW->getProb( state, i, j );
  logspace bw = pBW->getProb( state, i, j );
  return (fw * bw)/total;
}

double posterior( const string state, int i, int j, AlignBandingDPTable* pFW, AlignBandingDPTable* pBW ) {
  logspace total = pBW->getProb("start",0,0);
  logspace fw = pFW->getProb( state, i, j );
  logspace bw = pBW->getProb( state, i, j );
  return (fw * bw)/total;
}


// AlignDPTable version
double get_alignment( Path* path, char* iSeq1, char* iSeq2, int pad1id, int* seq2locus, int* y, int iLen1, char* line1, char* line2, char* line3, char* baq, char* pQuals, AlignDPTable* pFW, AlignDPTable* pBW) {

    int x = 0, aln=0;
    int i, j;
    *y = 0;
    char baqvalue;
    vector<double> minima;
    vector<int> indices;
    const int window = 15;
    double maximum = 0.0;

    // Loop over all transitions, except the last one (to the end state),
    // and the first one (into the delete state)
    j = (int)path->size()-1;
    for (i=1; i < j; i++) {
      if (path->toState(i) != pad1id) break;
      (*y)++;
    }

    *seq2locus = i-1;

    while (x < iLen1) {

      const vector<int>& v = path->emission(i++);

      if (line3 != NULL || baq != NULL) {
	int stateId = path->toState(i-1);
	string state = pFW->stateId[stateId];
	double pdouble = posterior(state, x, *y, pFW, pBW);
	if (line3) {
	  int p = (int)(0.5 + 25.0*pdouble);
	  line3[aln] = 'a' + (p<0 ? 0 : (p>25 ? 25 : p));
	}
	if (baq) {
	  baqvalue = pQuals[x];
	  if ((pdouble < 1.0) && (pdouble > 0.0))
	    baqvalue = min( baqvalue, char(33 + 0.5 + min(127-33,int(log(1.0-pdouble)*(-4.34294481903252)))) );
	  else if (pdouble <= 0.0)
	    baqvalue = char(33);
	  baq[x] = char(64 + pQuals[x] - baqvalue);
	}
	// implement maximum minimum over window
	if (v[1] == 0) pdouble=0;
	j = minima.size();
	while (j>0) {
	  j--;
	  if (minima[j] < pdouble) 
	    j=0;
	  else {
	    minima.pop_back();
	    indices.pop_back();
	  }
	  minima.push_back(pdouble);
	  indices.push_back(x);
	  while (indices[0] < x - window) {
	    minima.erase( minima.begin() );
	    indices.erase( indices.begin() );
	  }
	  if (x >= window && (maximum < minima[0]))
	    maximum = minima[0];
	}
      }

      if (v[0] == 1) {
	line1[aln] = iSeq1[x++];
      } else {
	line1[aln] = '-';
      }

      if (v[1] == 1) {
	line2[aln++] = iSeq2[(*y)++];
      } else {
	line2[aln++] = '-';
      }

    }

    line1[aln] = 0;
    line2[aln] = 0;
    if (line3 != NULL) line3[aln] = 0;
    if (baq != NULL) baq[x] = 0;
    
    return maximum;
}


// Align*Banding*DPTable version
double get_alignment( Path* path, char* iSeq1, char* iSeq2, int pad1id, int* seq2locus, int* y, int iLen1, char* line1, char* line2, char* line3, char* baq, char* pQuals, AlignBandingDPTable* pFW, AlignBandingDPTable* pBW) {

    int x = 0, aln=0;
    int i, j;
    *y = 0;
    char baqvalue;
    vector<double> minima;
    vector<int> indices;
    const int window = 15;
    double maximum = 0.0;
    int padding = 1;   // flag to indicate initial padding insertion/deletion

    // Loop over all transitions, except the last one (to the end state),
    // and the first one (into the delete/padding state)
    j = (int)path->size()-1;
    for (i=1; i < j; i++) {
      if (path->toState(i) != pad1id) break;
      (*y)++;
    }

    *seq2locus = i-1;
    
    while (x < iLen1) {

      const vector<int>& v = path->emission(i++);

      if (line3 != NULL || baq != NULL) {
	int stateId = path->toState(i-1);
	string state = pFW->stateId[stateId];
	double pdouble = posterior(state, x, *y, pFW, pBW);
	if (line3) {
	  int p = (int)(0.5 + 25.0*pdouble);
	  line3[aln] = 'a' + (p<0 ? 0 : (p>25 ? 25 : p));
	}
	if (baq) {
	  baqvalue = pQuals[x];
	  if ((pdouble < 1.0) && (pdouble > 0.0))
	    baqvalue = min( baqvalue, char(33 + 0.5 + min(127-33,int(log(1.0-pdouble)*(-4.34294481903252)))) );
	  else if (pdouble <= 0.0)
	    baqvalue = char(33);
	  baq[x] = char(64 + pQuals[x] - baqvalue);
	}
	// implement maximum minimum over window
	j = minima.size();
	if (v[1] == 0) pdouble=0;
	while (j>0) {
	  j--;
	  if (minima[j] < pdouble) 
	    j=0;
	  else {
	    minima.pop_back();
	    indices.pop_back();
	  }
	}
	minima.push_back(pdouble);
	indices.push_back(x);
	while (indices[0] < x - window) {
	  minima.erase( minima.begin() );
	  indices.erase( indices.begin() );
	}
	if (x >= window && (maximum < minima[0])) {
	  maximum = minima[0];
	}
      }

      if (v[0] == 1 && v[1] == 1) padding = 0;

      if (v[0] == 1) {
	line1[aln] = iSeq1[x++];
      } else {
	line1[aln] = '-';
      }

      if (v[1] == 1) {
	if (padding) {
	  // a deletion in the padding section (possibly following an insertion) -- exclude deletion from alignment
	  (*y)++;
	  //Bugfix: if this happened, *seq2locus would be wrong, leading to wrong mapping positions!!
	  (*seq2locus)++;
	} else {
	  line2[aln++] = iSeq2[(*y)++];
	}
      } else {
	line2[aln++] = '-';
      }

    }

    line1[aln] = 0;
    line2[aln] = 0;
    if (line3 != NULL) line3[aln] = 0;
    if (baq != NULL) baq[x] = 0;

    return maximum;
}


// Unbanded version
unsigned int sensitiveAlign( char* iSeq1, char* iSeq2, int iLen1, int iLen2, char* pQuals1, int iInsNucPrior, int iGapOpen, int iGapExtend, 
			     char* line1, char* line2, char* line1q, char* baq1, int maxlikelihood, 
			     char* line3, char* line4, char* line3q, char* baq3, int iStartMean, int iStartSD, double* amapq ) {

  // Returns an alignment phred score, computed in the obvious way, except that a penalty of 6 is levied for every unaligned
  // nucleotide within the read footprint, representing the equilibrium distribution.  Effectively this increases the gap extension penalty by 6.

  AlignDPTable* pViterbiTable;
  AlignDPTable* pFWTable;
  AlignDPTable* pBWTable;
  const double db = 10.0/log(0.1);
  int i;
  
  int score = int (0.5 + log(Viterbi_recurse(&pViterbiTable, pQuals1, iSeq1, iSeq2, iGapExtend, iGapOpen,iInsNucPrior, iLen1, iLen2, iStartMean, iStartSD)) * db);

  if (score > maxlikelihood) {
    // bail out, with correct score but bogus alignment
    delete pViterbiTable;
    for (i=0;i<iLen1;i++) {
      if (line1 != NULL) {
	line1[i] = iSeq1[i];
	line2[i] = iSeq2[i];
      }
      if (line3 != NULL) {
	line3[i] = iSeq1[i];
	line4[i] = iSeq2[i];
      }
      if (line1q != NULL) line1q[i] = 'z';
      if (line3q != NULL) line3q[i] = 'z';
      if (baq1 != NULL) baq1[i] = pQuals1[i];
      if (baq3 != NULL) baq3[i] = pQuals1[i];
    }
    if (line1 != NULL) {
      line1[iLen1] = 0;
      line2[iLen1] = 0;
    }
    if (line3 != NULL) {
      line3[iLen1] = 0;
      line4[iLen1] = 0;
    }
    if (line1q != NULL) line1q[i] = 0;
    if (line3q != NULL) line3q[i] = 0;
    if (baq1 != NULL) baq1[i] = 0;
    if (baq3 != NULL) baq3[i] = 0;
    return min(65535,max(0,score));
  }
  
  if (line1q != NULL) {
    Forward(&pFWTable, pQuals1, iSeq1, iSeq2, iGapExtend, iGapOpen, iInsNucPrior, iLen1, iLen2, iStartMean, iStartSD);
    Backward(&pBWTable, pQuals1, iSeq1, iSeq2, iGapExtend, iGapOpen, iInsNucPrior, iLen1, iLen2, iStartMean, iStartSD);
  }

  Path* path = &Viterbi_trace(pViterbiTable, pQuals1, iSeq1, iSeq2, iGapExtend, iGapOpen, iInsNucPrior, iLen1, iLen2, iStartMean, iStartSD, 1.001);
  int pad1id = pViterbiTable->getId("pad1");

  int seq2locus, y1;
  *amapq = get_alignment( path, iSeq1, iSeq2, pad1id, &seq2locus, &y1, iLen1, line1, line2, line1q, baq1, pQuals1, pFWTable, pBWTable );
  delete path;

  if (line3 != NULL) {
    Path* altpath = &Viterbi_trace(pViterbiTable, pQuals1, iSeq1, iSeq2, iGapExtend, iGapOpen, iInsNucPrior, iLen1, iLen2, iStartMean, iStartSD, 0.999);
    int y0;
    get_alignment( altpath, iSeq1, iSeq2, pad1id, &y0, &y1, iLen1, line3, line4, line3q, baq3, pQuals1, pFWTable, pBWTable );
    delete altpath;
  }

  delete pViterbiTable;
  if (line1q) {
    delete pFWTable;
    delete pBWTable;
  }

  // return both score and seq2locus; don't clutter parameter list with another pointer
  // negative scores should be impossible, but be safe
  return max(0,score) + (seq2locus << 16);

}




unsigned int sensitiveAlign_0( char* iSeq1, char* iSeq2, int iLen1, int iLen2, char* pQuals1, int iInsNucPrior, int iGapOpen, int iGapExtend, 
			       char* line1, char* line2, char* line1q, char* baq1, int maxlikelihood, 
			       char* line3, char* line4, char* line3q, char* baq3, int iStartMean, int iStartSD, int band, double* amapq ) {

  AlignBandingDPTable* pViterbiTable;
  AlignBandingDPTable* pFWTable;
  AlignBandingDPTable* pBWTable;
  const double db = 10.0/log(0.1);
  int i;

  int score = int (0.5 + log(ViterbiBanding_recurse(&pViterbiTable, pQuals1, iSeq1, iSeq2, iGapExtend, iGapOpen,iInsNucPrior, iLen1, iLen2, iStartMean, iStartSD, band)) * db);

  if (score > maxlikelihood) {
    // bail out, with correct score but bogus alignment
    delete pViterbiTable;
    for (i=0;i<iLen1;i++) {
      if (line1 != NULL) {
	line1[i] = iSeq1[i];
	line2[i] = iSeq2[i];
      }
      if (line3 != NULL) {
	line3[i] = iSeq1[i];
	line4[i] = iSeq2[i];
      }
      if (line1q != NULL) line1q[i] = 'z';
      if (line3q != NULL) line3q[i] = 'z';
      if (baq1 != NULL) baq1[i] = pQuals1[i];
      if (baq3 != NULL) baq3[i] = pQuals1[i];
    }
    if (line1 != NULL) {
      line1[iLen1] = 0;
      line2[iLen1] = 0;
    }
    if (line3 != NULL) {
      line3[iLen1] = 0;
      line4[iLen1] = 0;
    }
    if (line1q != NULL) line1q[i] = 0;
    if (line3q != NULL) line3q[i] = 0;
    if (baq1 != NULL) baq1[i] = 0;
    if (baq3 != NULL) baq3[i] = 0;
    return min(65535,max(0,score));
  }
  
  if (line1q != NULL) {
    ForwardBanding(&pFWTable, pQuals1, iSeq1, iSeq2, iGapExtend, iGapOpen, iInsNucPrior, iLen1, iLen2, iStartMean, iStartSD, band);
    BackwardBanding(&pBWTable, pQuals1, iSeq1, iSeq2, iGapExtend, iGapOpen, iInsNucPrior, iLen1, iLen2, iStartMean, iStartSD, band);
  }

  Path* path = &ViterbiBanding_trace(pViterbiTable, pQuals1, iSeq1, iSeq2, iGapExtend, iGapOpen, iInsNucPrior, iLen1, iLen2, iStartMean, iStartSD, 1.001);
  int pad1id = pViterbiTable->getId("pad1");

  int seq2locus, y1;
  *amapq = get_alignment( path, iSeq1, iSeq2, pad1id, &seq2locus, &y1, iLen1, line1, line2, line1q, baq1, pQuals1, pFWTable, pBWTable );
  delete path;

  if (line3 != NULL) {
    Path* altpath = &ViterbiBanding_trace(pViterbiTable, pQuals1, iSeq1, iSeq2, iGapExtend, iGapOpen, iInsNucPrior, iLen1, iLen2, iStartMean, iStartSD, 0.999);
    int y0;
    get_alignment( altpath, iSeq1, iSeq2, pad1id, &y0, &y1, iLen1, line3, line4, line3q, baq3, pQuals1, pFWTable, pBWTable );
    delete altpath;
  }

  delete pViterbiTable;
  if (line1q) {
    delete pFWTable;
    delete pBWTable;
  }

  // return both score and seq2locus; don't clutter parameter list with another pointer
  // negative scores should be impossible, but be safe
  return max(0,score) + (seq2locus << 16);

}






// Returns an alignment phred score, computed in the obvious way, except that a penalty of 6 is levied for every unaligned
// nucleotide within the read footprint, representing the equilibrium distribution.  Effectively this increases the gap extension penalty by 6.

extern "C" unsigned int sensitiveAlign( char* iSeq1, char* iSeq2, int iLen1, int iLen2, char* pQuals1, int iInsNucPrior, int iGapOpen, int iGapExtend, 
					char* line1, char* line2, char* line1q, char* baq1, int maxlikelihood, 
					char* line3, char* line4, char* line3q, char* baq3, int iStartMean, int iStartSD, int band, double* amapq ) {

  if (band < 0) 
    return sensitiveAlign(iSeq1, iSeq2, iLen1, iLen2, pQuals1, iInsNucPrior, iGapOpen, iGapExtend, 
			  line1, line2, line1q, baq1, maxlikelihood, line3, line4, line3q, baq3, iStartMean, iStartSD, amapq);

  // hardcoded maximum
  const int global_maxband = 150;
  const int global_minband = 25;

  if (iLen2-iLen1>band) band = iLen2-iLen1;

  // First, check if band is small anyway; if so align as usual
  if (band <= 15) {
    return sensitiveAlign_0( iSeq1, iSeq2, iLen1, iLen2, pQuals1, iInsNucPrior, iGapOpen, iGapExtend, line1, line2, line1q, baq1, maxlikelihood,
			     line3, line4, line3q, baq3, iStartMean, iStartSD, band, amapq );
  }

  int yposmin;
  int iNewlen2;
  int newWidth;

  findband( iSeq1, iSeq2, iLen1, iLen2, pQuals1, iInsNucPrior, iGapOpen, iGapExtend, band, global_maxband, global_minband, &yposmin, &iNewlen2, &newWidth );

  unsigned int r = sensitiveAlign_0( iSeq1, iSeq2+yposmin, iLen1, iNewlen2, pQuals1, iInsNucPrior, iGapOpen, iGapExtend, 
				     line1, line2, line1q, baq1, 
				     maxlikelihood,
				     line3, line4, line3q, baq3, 
				     iStartMean-yposmin, iStartSD, newWidth, amapq );

  r += yposmin << 16;
  
  // Done
  return r;
  
}
