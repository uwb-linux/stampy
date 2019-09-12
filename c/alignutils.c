#include "../ext/readalign/frontend.h"
#include "alignutils.h"

/*****************************************************************************************************************
 This code is copyright (c) Gerton Lunter, Jan 2009
 It may not be distributed, made public, or used in other software without the permission of the copyright holder
******************************************************************************************************************/

//
// Todo: move some code from the indel caller to the SSE2 aligner:
//  - deletion followed by insertion
//  - gapq instead of gap opening
//  - dealing with Ns
// Move some code to the indel caller:
//  - Longer reads through -0x8000 initialization and higher ceiling


#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef __SSE2__
#include <emmintrin.h>
#endif

#include "Python.h"




// KYT->0, C->1, WMRA->2, G->3, N->4
static int dnaconvert[256] = {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
			      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
			      4, 2, 4, 1, 4, 4, 4, 3, 4, 4, 4, 0, 4, 2, 4, 4, 4, 4, 2, 4, 0, 4, 4, 2, 4, 0, 4, 4, 4, 4, 4, 4,
			      4, 2, 4, 1, 4, 4, 4, 3, 4, 4, 4, 0, 4, 2, 4, 4, 4, 4, 2, 4, 0, 4, 4, 2, 4, 0, 4, 4, 4, 4, 4, 4,
			      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
			      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
			      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
			      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

void inline convert( char* seq, char* output ) {
  int c;
  while ((c = *(seq++))) *(output++) = dnaconvert[c];
}

void inline revseq( char* seq, int len ) {
  int i;
  for (i=0; i<len/2; i++) {
    char c = seq[i];
    seq[i] = seq[ len - i - 1 ];
    seq[ len - i - 1] = c;
  }
}

inline int min(int i, int j) { return (i<j ? i : j); }
inline int max(int i, int j) { return (i>j ? i : j); }

void printseq( char* seq, int dim ) {
  int i;
  for (i=0; i<dim-1; i++) {
    printf("%c", (unsigned int)seq[i] < 4 ? "TGAC"[(int)seq[i]] : 'N');
  }
  printf("\n");
}

#define entry(i,j) dp[ (i) + dim1*(j) ]

void maketable( int* dp, char* seq1, char* seq2, int dim1, int dim2, int* scoretable, int gapstart ) {

  // scoretable is a 5x5 table of scores for matching nucleotides
  // code: 0=T 1=C 2=A 3=G 4=N

  // it is essential that the score table only contains multples of 4,
  // if not traceback will not work properly.  for speed this is not asserted.

  int i,j;

  // Banded version
  // band=1 means straight diagonal; 2 means 1 bp indels only

  // row j=0
  int maxi = min(BAND, dim1);
  int score;
  int belowentry;
  int temp;
  gapstart *= SCOREFACTOR;
  const int mininf = -0x3FF0;
  const int horiz = 1;
  const int vert = 2;
  const int diag = 3;
  for (i=0; i<maxi; i++) {
    dp[i] = gapstart*i + horiz;
  }

  if (BAND<dim1) dp[BAND] = mininf;

  // other rows
  for (j=1; j<dim2; j++) {

    // for speed
    int nuc2times5 = 5 * seq2[j-1];
    // the coordinate of the first entry to be filled
    i = max(j-BAND+1, 0);
    // points to (not beyond) the last entry to be filled
    // this assumes len(seq2) < len(seq1)
    maxi = min(j+BAND, dim1) - 1;
    // treat first element outside the main loop
    int* thisrowptr = &entry(i,j);
    int* prevrowptr = thisrowptr -1 -dim1;  /* &entry(i-1,j-1);  */     // can be out of bounds
    if (i>0) {
      // diagonal entry
      score = (*prevrowptr + scoretable[ seq1[i-1] + nuc2times5 ]) | diag;
      // drop a sentinel, so as not to complicate traceback
      *(thisrowptr-1) = mininf;
    } else {
      score = mininf;
    }
    // move to entry[i,j-1], below current
    prevrowptr++;
    // get entry below - the next round's diagonal entry
    belowentry = *prevrowptr;
    // add contribution - if i==maxi, band==1 and no valid left coordinate exists
    if ((i < maxi) && ((temp = belowentry + gapstart) > score)) {
      score = (temp & -4) + vert;
    }
    *thisrowptr = score;
    // now enter the main loop
    while (++i < maxi) {
      // score holds the value of the left entry, which is valid
      // the diagonal entry is valid too, and is held in belowentry
      if ((score += gapstart) > (temp = belowentry+scoretable[ seq1[i-1] + nuc2times5 ] )) {
	score = (score & -4) + horiz;
      } else {
	score = temp | diag;
      }
      belowentry = *(++prevrowptr);
      if ((temp = belowentry + gapstart) > score) {
	score = (temp & -4) + vert;
      }
      *(++thisrowptr) = score;
    }
    // if i==maxi (and not maxi+1 which happens when bottom == top; band==1), handle the top entry
    if (i==maxi) {
      if ((score += gapstart) > (temp = belowentry+scoretable[ seq1[i-1] + nuc2times5 ] )) {
	*(++thisrowptr) = (score & -4) + horiz;
      } else {
	*(++thisrowptr) = temp | diag;
      }
      ++i;
    }
    // now i == maxi+1; if still in range, drop a sentinel
    if (i<dim1) *(++thisrowptr) = mininf;
  }
}


void maketable_varband( int* dp, char* seq1, char* seq2, int dim1, int dim2, int* scoretable, int gapstart, int band ) {

  // scoretable is a 5x5 table of scores for matching nucleotides
  // code: 0=T 1=C 2=A 3=G 4=N

  // it is essential that the score table only contains multples of 4,
  // if not traceback will not work properly.  for speed this is not asserted.

  int i,j;

  // Banded version
  // band=1 means straight diagonal; 2 means 1 bp indels only

  // row j=0
  int maxi = min(band, dim1);
  int score;
  int belowentry;
  int temp;
  gapstart *= SCOREFACTOR;
  const int mininf = -0x3FF0;
  const int horiz = 1;
  const int vert = 2;
  const int diag = 3;
  for (i=0; i<maxi; i++) {
    dp[i] = gapstart*i + horiz;
  }

  if (band<dim1) dp[band] = mininf;

  // other rows
  for (j=1; j<dim2; j++) {

    // for speed
    int nuc2times5 = 5 * seq2[j-1];
    // the coordinate of the first entry to be filled
    i = max(j-band+1, 0);
    // points to (not beyond) the last entry to be filled
    // this assumes len(seq2) < len(seq1)
    maxi = min(j+band, dim1) - 1;
    // treat first element outside the main loop
    int* thisrowptr = &entry(i,j);
    int* prevrowptr = thisrowptr -1 -dim1;  /* &entry(i-1,j-1);  */     // can be out of bounds
    if (i>0) {
      // diagonal entry
      score = (*prevrowptr + scoretable[ seq1[i-1] + nuc2times5 ]) | diag;
      // drop a sentinel, so as not to complicate traceback
      *(thisrowptr-1) = mininf;
    } else {
      score = mininf;
    }
    // move to entry[i,j-1], below current
    prevrowptr++;
    // get entry below - the next round's diagonal entry
    belowentry = *prevrowptr;
    // add contribution - if i==maxi, band==1 and no valid left coordinate exists
    if ((i < maxi) && ((temp = belowentry + gapstart) > score)) {
      score = (temp & -4) + vert;
    }
    *thisrowptr = score;
    // now enter the main loop
    while (++i < maxi) {
      // score holds the value of the left entry, which is valid
      // the diagonal entry is valid too, and is held in belowentry
      if ((score += gapstart) > (temp = belowentry+scoretable[ seq1[i-1] + nuc2times5 ] )) {
	score = (score & -4) + horiz;
      } else {
	score = temp | diag;
      }
      belowentry = *(++prevrowptr);
      if ((temp = belowentry + gapstart) > score) {
	score = (temp & -4) + vert;
      }
      *(++thisrowptr) = score;
    }
    // if i==maxi (and not maxi+1 which happens when bottom == top; band==1), handle the top entry
    if (i==maxi) {
      if ((score += gapstart) > (temp = belowentry+scoretable[ seq1[i-1] + nuc2times5 ] )) {
	*(++thisrowptr) = (score & -4) + horiz;
      } else {
	*(++thisrowptr) = temp | diag;
      }
      ++i;
    }
    // now i == maxi+1; if still in range, drop a sentinel
    if (i<dim1) *(++thisrowptr) = mininf;
  }
}




#ifdef __SSE2__

int traceback_sse2( char* seq1, char* seq2, int* qual2, char* aln1, char* aln2, int len1, int len2, 
		    int gapopen, int gapextend, int nucprior, int* firstpos, int ab_cutoff ) {

  // seq2 is the read; the shorter of the sequences
  // no checks for overflow are done

  // the bottom-left and top-right corners of the DP table are just
  // included at the extreme ends of the diagonal, which measures
  // n=8 entries diagonally across.  This fixes the length of the
  // longer (horizontal) sequence to 14 (2*8-2) more than the shorter

  assert( len1 == len2 + 2*8 - 1);

  // make sure that special cases at beginning and end (initialization!)
  // do not mix.  (Todo: Check this assertion.)
  assert( len1 > 8);

  short gap_open = gapopen*4;
  short gap_extend = gapextend*4;
  short nuc_prior = nucprior*4;
  char traceback = (aln1 != NULL);

  register __m128i _m1;
  register __m128i _i1;
  register __m128i _d1;
  register __m128i _m2;
  register __m128i _i2;
  register __m128i _d2;
  __m128i _seq1win;
  __m128i _seq2win;
  __m128i _qual2win;

  __m128i _gap_extend_minus_gap_open = _mm_set1_epi16( gap_extend - gap_open );
  __m128i _gap_open = _mm_set1_epi16( gap_open );
  __m128i _gap_open_nuc_prior = _mm_set1_epi16( gap_open + nuc_prior );
  __m128i _three = _mm_set1_epi16( 3 );
  __m128i _abcutoff = _mm_set1_epi16( (ab_cutoff - 0x2000) << 2 );
  __m128i _initmask = _mm_set_epi16( 0,0,0,0,0,0,0,-1 );
  __m128i _initmask2 = _mm_set_epi16( 0,0,0,0,0,0,0,-0x8000 );

  __m128i _backpointers[ 2*(len1+8) ];
  const short pos_inf = 0x7800;

#define match_label 0
#define insert_label 1
#define delete_label 3

  // initialization
  _m1 = _mm_set_epi16( pos_inf, pos_inf, pos_inf, pos_inf, pos_inf, pos_inf, pos_inf, pos_inf );
  _i1 = _m1;
  _d1 = _m1;
  _m2 = _m1;
  _i2 = _m1;
  _d2 = _m1;

  // sequence 1 is initialized with the n-long prefix, in forward direction
  // sequence 2 is initialized as empty; reverse direction
  _seq1win = _mm_set_epi16( seq1[7], seq1[6], seq1[5], seq1[4], seq1[3], seq1[2], seq1[1], seq1[0] );
  _seq2win = _m1;
  _qual2win = _mm_set1_epi16(64*4);

  short _score = 0;
  short minscore = pos_inf;
  short minscoreidx = -1;
  
  // main loop.  Do one extra iteration, with nucs from sequence 2 just moved out
  // of the seq2win/qual arrays, to simplify getting back pointers
  int s = 0;
  for (; s<2*(len2+8-1 + 1); s+=2) {

    // seq1 is current; seq2 needs updating
    _seq2win = _mm_slli_si128( _seq2win, 2 );     
    _qual2win = _mm_slli_si128( _qual2win, 2 );
    if (s/2 < len2) {
      _seq2win = _mm_insert_epi16( _seq2win, seq2[ s/2 ], 0 );
      _qual2win = _mm_insert_epi16( _qual2win, 4*qual2[ s/2 ], 0 );
    } else {
      _seq2win = _mm_insert_epi16( _seq2win, '0', 0 );
      _qual2win = _mm_insert_epi16( _qual2win, 64*4, 0 );
    }

    //
    // S even
    //
    
    // initialize to -0x8000
    _m1 = _mm_or_si128( _initmask2, _mm_andnot_si128( _initmask, _m1 ) );
    // moved here, from below, to allow alignments to start with insertions or deletions
    _m2 = _mm_or_si128( _initmask2, _mm_andnot_si128( _initmask, _m2 ) );
    _m1 = _mm_min_epi16( _m1, _mm_min_epi16( _i1, _d1 ) );
    // at this point, extract minimum score.  Referred-to position must
    // be y==len2-1, so that current position has y==len2; i==0 so d=0 and y=s/2

    if (s/2 >= len2) {
      switch (s/2 - len2) {
      case 0: _score = _mm_extract_epi16( _m1, 0 ); break;
      case 1: _score = _mm_extract_epi16( _m1, 1 ); break;
      case 2: _score = _mm_extract_epi16( _m1, 2 ); break;
      case 3: _score = _mm_extract_epi16( _m1, 3 ); break;
      case 4: _score = _mm_extract_epi16( _m1, 4 ); break;
      case 5: _score = _mm_extract_epi16( _m1, 5 ); break;
      case 6: _score = _mm_extract_epi16( _m1, 6 ); break;
      case 7: _score = _mm_extract_epi16( _m1, 7 ); break;
      }
      if (_score < minscore) {
	minscore = _score;
	minscoreidx = s;     // point back to the match state at this entry, so as not to
      }                      // have to store the state at s-2
    }
    
    _m1 = _mm_add_epi16( _m1, _mm_andnot_si128( _mm_cmpeq_epi16( _seq2win, _seq1win ), _qual2win ) );
    _d1 = _mm_add_epi16( _mm_min_epi16( _mm_add_epi16( _d2, _gap_extend_minus_gap_open ), _m2 ), _gap_open );
    _d1 = _mm_insert_epi16( _mm_slli_si128( _d1, 2 ), pos_inf, 0 );
    _i1 = _mm_add_epi16( _mm_min_epi16( _mm_add_epi16( _i2, _gap_extend_minus_gap_open ), _m2 ), _gap_open_nuc_prior );

    // get back-pointers and store
    if (traceback) {
      _backpointers[ s ] = _mm_or_si128( _mm_or_si128( _mm_and_si128( _three, _m1 ),
						       _mm_slli_epi16( _mm_and_si128( _three, _i1 ), 2*insert_label ) ),
					 _mm_slli_epi16( _mm_and_si128( _three, _d1 ), 2*delete_label ) );
    
      // set state labels
      _m1 = _mm_andnot_si128( _three, _m1 );
      _i1 = _mm_or_si128( _mm_andnot_si128( _three, _i1 ), _mm_srli_epi16( _three, 1 ) );           // _mm_set1_epi16( insert_label ) );
      _d1 = _mm_or_si128( _mm_andnot_si128( _three, _d1 ), _three );                                // _mm_set1_epi16( delete_label ) );
    } else {
      // alpha-beta cutoff
      if ( (s/2 < len2) &&
	   (_mm_movemask_epi8(_mm_cmpgt_epi16(_m1,_abcutoff)) == 0xFFFF) &&
	   (_mm_movemask_epi8(_mm_cmpgt_epi16(_m2,_abcutoff)) == 0xFFFF)) {
	return ab_cutoff;
      }
    }

    //
    // S odd
    //

    // seq1 needs updating; seq2 is current
    char c = (8 + s/2 < len1) ? seq1[ 8+(s/2) ] : 'G';
    _seq1win = _mm_insert_epi16( _mm_srli_si128( _seq1win, 2 ), c, 8-1 );

    // Moved up:
    //_m2 = _mm_andnot_si128( _initmask, _m2 );
    _initmask = _mm_slli_si128( _initmask, 2 );
    _initmask2 = _mm_slli_si128( _initmask2, 2 );
    _m2 = _mm_min_epi16( _m2, _mm_min_epi16( _i2, _d2 ) );

    // at this point, extract minimum score.  Referred-to position must
    // be y==len2-1, so that current position has y==len2; i==0 so d=0 and y=s/2
    if (s/2 >= len2) {
      switch (s/2 - len2) {
      case 0: _score = _mm_extract_epi16( _m2, 0 ); break;
      case 1: _score = _mm_extract_epi16( _m2, 1 ); break;
      case 2: _score = _mm_extract_epi16( _m2, 2 ); break;
      case 3: _score = _mm_extract_epi16( _m2, 3 ); break;
      case 4: _score = _mm_extract_epi16( _m2, 4 ); break;
      case 5: _score = _mm_extract_epi16( _m2, 5 ); break;
      case 6: _score = _mm_extract_epi16( _m2, 6 ); break;
      case 7: _score = _mm_extract_epi16( _m2, 7 ); break;
      }
      if (_score < minscore) {
	minscore = _score;
	minscoreidx = s+1;
      }
    }

    _m2 = _mm_add_epi16( _m2, _mm_andnot_si128( _mm_cmpeq_epi16( _seq2win, _seq1win ), _qual2win ) );
    _d2 = _mm_add_epi16( _mm_min_epi16( _mm_add_epi16( _d1, _gap_extend_minus_gap_open ), _m1 ), _gap_open );
    _i2 = _mm_insert_epi16( _mm_srli_si128( _mm_add_epi16( _mm_min_epi16( _mm_add_epi16( _i1, _gap_extend_minus_gap_open ), _m1 ), _gap_open_nuc_prior ), 2 ), pos_inf, 8-1 );
    
    // get back-pointers and store
    if (traceback) {
      _backpointers[ s+1 ] = _mm_or_si128( _mm_or_si128( _mm_and_si128( _three, _m2 ),
							 _mm_slli_epi16( _mm_and_si128( _three, _i2 ), 2*insert_label ) ),
					   _mm_slli_epi16( _mm_and_si128( _three, _d2 ), 2*delete_label ) );
      
      // set state labels
      _m2 = _mm_andnot_si128( _three, _m2 );
      _i2 = _mm_or_si128( _mm_andnot_si128( _three, _i2 ), _mm_srli_epi16( _three, 1 ) );          //_mm_set1_epi16( insert_label ) );
      _d2 = _mm_or_si128( _mm_andnot_si128( _three, _d2 ), _three );                               // mm_set1_epi16( delete_label ) );
    }

  }

  // Backtrace.
  if (!traceback) {
    return (minscore + 0x8000) >> 2;
  }

  //DEBUG
  //printf("Backtracing starting from s=%u length=%u\n",minscoreidx,2*(len1+8));

  s = minscoreidx;    // point to the dummy match transition
  short i = s/2 - len2;
  short y = len2;
  short x = s - y;
  short alnidx = 0;
  short state = ((((short*)( _backpointers + s))[i]) >> (2*match_label)) & 3;
  s -= 2;

  // this is 2*y (s even) or 2*y+1 (s odd)
  while (y > 0) {
    short newstate = ((((short*)( _backpointers + s))[i]) >> (2*state)) & 3;
    if (state == match_label) {
      s -= 2;
      aln1[alnidx] = seq1[--x];
      aln2[alnidx] = seq2[--y];
    } else if (state == insert_label) {
      i += s&1;
      s -= 1;
      aln1[alnidx] = '-';
      aln2[alnidx] = seq2[--y];
    } else {
      s -= 1;
      i -= s&1;
      aln1[alnidx] = seq1[--x];
      aln2[alnidx] = '-';
    }
    state = newstate;
    alnidx++;
  }
  aln1[alnidx] = 0;
  aln2[alnidx] = 0;

  if (firstpos != NULL) *firstpos = x;

  // reverse them
  int j;
  for (i=0, j=alnidx-1; i<j; i++, j--) {
    x = aln1[i]; 
    y = aln2[i]; 
    aln1[i]=aln1[j]; 
    aln2[i]=aln2[j]; 
    aln1[j] = x;
    aln2[j] = y;
  }

  return (minscore + 0x8000) >> 2;
}

#endif  // __SSE2__



PyObject* traceback( int* dp, char* seq1, char* seq2, int dim1, int dim2, int *scoretable, int gapstart, char reverse, int band ) {

  int i;

  // this is different from the python version
  if (!reverse) {
    revseq( seq1, dim1-1 );
    revseq( seq2, dim2-1 );
  }

  maketable_varband( dp, seq1, seq2, dim1, dim2, scoretable, gapstart, band );

  int j = dim2-1;
  int mini = min(dim2-1-band, dim1);
  int maxi = min(dim2-1+band, dim1);
  int startpos = mini;
  int startmax = entry(startpos,j) & -4;
  int temp;
  for (i=mini+1; i<maxi; i++) {
    if ((temp = entry(i,j) & -4) > startmax) {
      startmax = temp;
      startpos = i;
    }
  }

  i = startpos;

  char aln1[ 2*MAXSEQLEN+1 ];
  char aln2[ 2*MAXSEQLEN+1 ];
  int alignlen = 0;
  int direction;
  char* convert = "TCAGN";

  while ( (i>0) || (j>0) ) {

    direction = entry(i,j) & 3;
    if (direction & 1) {
      i -= 1;
      aln1[ alignlen ] = convert[ (int)seq1[i] ];
    } else {
      aln1[ alignlen ] = '-';
    }

    if (direction & 2) {
      j -= 1;
      aln2[ alignlen ] = convert[ (int)seq2[j] ];
    } else {
      aln2[ alignlen ] = '-';
    }

    alignlen += 1;
    
  }

  // for debugging
  aln1[ alignlen ] = 0;
  aln2[ alignlen ] = 0;

  if (reverse) {
    revseq( aln1, alignlen );
    revseq( aln2, alignlen );
  }

  PyObject* align1 = PyString_FromStringAndSize( aln1, alignlen );
  PyObject* align2 = PyString_FromStringAndSize( aln2, alignlen );
  PyObject* firstcoord = PyInt_FromLong( startpos-1 );
  PyObject* retvalue = PyTuple_New(3);
  PyTuple_SET_ITEM(retvalue, 0, align2);     // reversed.  it would be more proper to do this
  PyTuple_SET_ITEM(retvalue, 1, align1);     // in the calling routine, as is done in Python,
  PyTuple_SET_ITEM(retvalue, 2, firstcoord); // but that invovles more overhead.

  // in contrast to python version, this implementation only returns the
  // first 'local' coordinate rather than a list of all coordinates

  return retvalue;

}



PyObject* alignprefix_c_rawscore( int* dp, char* seq1, int seq1len, char* seq2, int seq2len, int* scoretable, int gapstart, int band) {

  char numseq1[MAXSEQLEN+1];
  char numseq2[MAXSEQLEN+1];

  if (seq2len < seq1len || seq1len > MAXSEQLEN || seq2len > MAXSEQLEN) {
    printf("Error (alignprefix_c_rawscore): genome sequence shorter than read sequence, or either too long\n");
    Py_INCREF(Py_None);
    return Py_None;
  }

  convert( seq1, numseq1 );
  convert( seq2, numseq2 );

  return traceback( dp, numseq2, numseq1, seq2len+1, seq1len+1, scoretable, gapstart, 1, band );

}


PyObject* alignsuffix_c_rawscore( int* dp, char* seq1, int seq1len, char* seq2, int seq2len, int* scoretable, int gapstart, int band) {

  char numseq1[MAXSEQLEN+1];
  char numseq2[MAXSEQLEN+1];

  if (seq2len < seq1len || seq1len > MAXSEQLEN || seq2len > MAXSEQLEN) {
    printf("Error (alignsuffix_c_rawscore): genome sequence shorter than read sequence, or either too long\n");
    Py_INCREF(Py_None);
    return Py_None;
  }

  convert( seq1, numseq1 );
  convert( seq2, numseq2 );

  return traceback( dp, numseq2, numseq1, seq2len+1, seq1len+1, scoretable, gapstart, 0, band );

}


void test_sensitiveAlign(void)
{
  int n,m,k,band;
  for (n=1; n<10; n++) {
    for (m=1; m<10; m++) {
      char seq1[12+n+1];
      char seq2[30+12+m+1];
      char quals1[12+n+1];
      char line1[1000];
      char line2[1000];
      char line1q[1000];
      char baq[1000];
      double amapq[1000];
      for (k=0; k<6; k++) {
	seq1[k] = 'A';
	seq1[k+6+n] = 'T';
	quals1[k] = quals1[k+6+n] = '@';
	seq2[k+15] = 'A';
	seq2[k+15+6+m] = 'T';
      }
      for (k=0; k<n; k++) {
	seq1[k+6] = 'C';
	quals1[k+6] = '@';
      }
      for (k=0; k<m; k++) {
	seq2[k+15+6] = 'C';
      }
      for (k=0; k<15; k++) {
	seq2[k] = 'G';
	seq2[k+15+12+m] = 'G';
      }
      seq1[12+n] = 0;
      seq2[30+12+m] = 0;
      // print
      printf("Sequences:\n");
      for (k=0; seq1[k]; k++) printf("%c",seq1[k]); printf("\n");
      for (k=0; seq2[k]; k++) printf("%c",seq2[k]); printf("\n");
      for (band=1; band<10; band++) {
	printf ("Band=%u\n",band);
	sensitiveAlign( seq1, seq2, strlen(seq1), strlen(seq2), quals1, 3, 40, 3,
			line1, line2, line1q, baq, 99999, NULL, NULL, NULL, NULL, 10, 100, band, amapq);
      }
    }
  }
}

// mode == 0: return None instead of the alignment in entries 1,2
// mode == 1: return (score,alignment1,alignment2,seq2_locus)
// mode == 2: return (score,alignment1,alignment2,seq2_locus,alt_alignment1,alt_alignment2)
// mode == 4: return (score,alignment1,alignment2,seq2_locus,alt_alignment1,alt_alignment2,qual1,qual2)
PyObject* align_sensitive( char* seq1, int seq1len, char* seq2, int seq2len, char* pQuals1, 
			   int iGapOne, int iGapOpen, int iGapExtend, int maxlikelihood, int mode, int iStartMean, int iStartSD, int band ) {

  int maxalnlen = seq1len + seq2len + 1;

  char line1[maxalnlen];
  char line2[maxalnlen];
  char line3[maxalnlen];
  char line4[maxalnlen];
  char line1q[maxalnlen];
  char line3q[maxalnlen];
  char baq1[seq1len+1];
  char baq3[seq1len+1];
  double amapq;

  if ((mode<0)||(mode>5)) mode=0;

  //DEBUG
  //test_sensitiveAlign();

  unsigned int result = sensitiveAlign( seq1, seq2, seq1len, seq2len, pQuals1, iGapOne, iGapOpen, iGapExtend, 
					line1, line2, mode == 4 ? line1q : NULL, mode == 4 ? baq1 : NULL, maxlikelihood,
					mode >= 2 ? line3 : NULL, mode >= 2 ? line4 : NULL, mode == 4 ? line3q : NULL, mode == 4 ? baq3 : NULL,
					iStartMean, iStartSD, band, &amapq );

  PyObject* retval = PyTuple_New(4 + 2*(mode>=2) + 5*(mode==4) );
  PyTuple_SET_ITEM(retval, 0, PyInt_FromLong( result & 0xFFFF ) );
  if (mode >= 1) {
    PyTuple_SET_ITEM(retval, 1, PyString_FromString( line1 ) );
    PyTuple_SET_ITEM(retval, 2, PyString_FromString( line2 ) );
  } else {
    Py_INCREF(Py_None);    
    Py_INCREF(Py_None);    
    PyTuple_SET_ITEM(retval, 1, Py_None);
    PyTuple_SET_ITEM(retval, 2, Py_None);
  }
  PyTuple_SET_ITEM(retval, 3, PyInt_FromLong( result >> 16 ) );

  if (mode >= 2) {
    PyTuple_SET_ITEM(retval, 4, PyString_FromString( line3 ) );
    PyTuple_SET_ITEM(retval, 5, PyString_FromString( line4 ) );
  }

  if (mode == 4) {
    PyTuple_SET_ITEM(retval, 6, PyString_FromString( line1q ) );
    PyTuple_SET_ITEM(retval, 7, PyString_FromString( line3q ) );
    PyTuple_SET_ITEM(retval, 8, PyString_FromString( baq1 ) );
    PyTuple_SET_ITEM(retval, 9, PyString_FromString( baq3 ) );
    PyTuple_SET_ITEM(retval, 10, PyFloat_FromDouble( amapq ) );
  }

  return retval;
  
}



void findband( char* iSeq1, char* iSeq2, int iLen1, int iLen2, char* pQuals1, int iInsNucPrior, int iGapOpen, int iGapExtend, int band, int global_maxband, int global_minband,
	       int* yposmin, int* iNewlen2, int* newWidth )
{

  char genomealignment[ 2*iLen1+16 ];
  char readalignment[ 2*iLen1+16 ];
  int intquals1[ iLen1 ];
  char genomeseq[ iLen2 + 32 ];
  int i,k;
  int minband = band;   // bottom-most extent of band
  int maxband = -band;  // top-most extent of band
  int maxbandlen=0, spine=0;

  // copy genome sequence from iSeq2, with padding; make int quals
  for (i=0;i<16;i++) genomeseq[i] = genomeseq[iLen2+16+i] = 'G';
  for (i=0; i<iLen2;i++) genomeseq[i+16] = iSeq2[i];
  for (i=0; i<iLen1; i++) intquals1[i] = pQuals1[i] - 33;

  // Code below copied from maputils.c!
  int sumquality = 0;
  for (i=0; i<iLen1; i++) sumquality += intquals1[i];
  // we have 14 bits of score to play with; algorithm uses 0x200 (0x800/4) headroom; allow 16*64=0x400 for flanks;
  // altogether 0x3a00.  Scale down quality scores as required, and scale them up later
  int divisor = 1 + (sumquality + 256) / 0x3a00;
  if (divisor > 1) {
    for (i=0; i<iLen1; i++) intquals1[i] /= divisor;
    iGapOpen = iGapOpen==0 ? 0 : 1 + (iGapOpen-1)/divisor;
    iGapExtend = iGapExtend==0 ? 0 : 1 + (iGapExtend-1)/divisor;
    iInsNucPrior = iInsNucPrior==0 ? 0 : 1 + (iInsNucPrior-1)/divisor;
  }

  // Pre-align on 15 bp wide bands

  // SSE2 aligner needs vertical (seq2, genome) boundaries.  So calculate the maximum extent of the top
  // and bottom diagonals; calculate the maximum seq1 length, calculate seq2 length as the length
  // of seq1 plus 15, and introduce dummy nucs where needed.

  int maxlocallen[iLen1];
  int localband[iLen1];
  int bandlen[band+32];
  for (i=0; i<iLen1; i++) {
    maxlocallen[i] = 0;
    localband[i] = 0;
  }
  for (i=0; i<band+32; i++) bandlen[i] = 0;
  // Loop over width index: -width/2 ... width/2; bail out when a perfect match was found
  // Traverse it from middle to either end, to optimize bail-out
  //for (i=-band/2;(i<band/2 + 14) && (global_minband > 1);i+=15) {
  int half_idx = 0;
  int max_index = ((band/2 + 14) - (-band/2) + 14) / 15;
  int max_half_index = max_index - (max_index / 2);   // right half of interval
  int leftright_index;
  // loop over all indices in [0, max_index), but starting in the middle and alternating to each end
  for (half_idx = 0; (half_idx < max_half_index) && (global_minband > 1); half_idx++) {
   // first do right index, then left, but only if left index is in the interval [0, max_index/2 - 1 ]
   for (leftright_index = 0; leftright_index < 1 + (max_index/2 - 1 - half_idx >= 0 && global_minband > 1); leftright_index++) {
    if (leftright_index == 0) { 
      i = (-band/2) + (max_index/2 + half_idx)*15; 
    } else { 
      i = (-band/2) + (max_index/2 - 1 - half_idx)*15; 
    }

    //DEBUG
    //printf("Testing i=%d\n",i);

    // Last band may overlap by at most 14bp, to ensure full coverage of seq2
    if (i >= band/2) i = band/2-1;

    //printf("Testing band i=%d\n",i);
    
    // Calculate Y index into iSeq2 at botton of current band
    int ypos = (iLen2 - iLen1)/2 + i;

    //printf("Initial ypos=%d\n",ypos);
    
    // Calculate minimum X index.  Require longer hits at extreme ends of window
    int matchsize = 8;
    int xposmin = 0;
    if (ypos < -15) {
      xposmin += (-15-ypos);
      ypos = -15;
      matchsize = 10;
    }
    //printf("Truncated ypos=%d\n",ypos);

    // Calculate maximum X index (+1)
    int xposmax = iLen1;
    if (ypos + 15 + iLen1 >= iLen2+15) {
      xposmax = iLen2 - ypos;
      matchsize = 10;
    }

    // Calculate lengths of subsequences
    int seq1len = xposmax - xposmin;
    int seq2len = seq1len + 15;
    char* seq1 = iSeq1 + xposmin;
    int* pquals1 = intquals1 + xposmin;
    char* seq2 = genomeseq + 16 + ypos;
    int firstpos;

    // Align
    //printf("Aligning seq1len=%d seq2len=%d\n",seq1len,seq2len);
    //printf("xposmin/max=%d %d\n",xposmin, xposmax);
    //int jj;
    //for (jj=0; jj<seq1len; jj++) printf("%c",seq1[jj]); printf("\n");
    //for (jj=0; jj<seq2len; jj++) printf("%c",seq2[jj]); printf("\n");
    int score = traceback_sse2( seq2, seq1, pquals1, genomealignment, readalignment, seq2len, seq1len,
				iGapOpen, iGapExtend, iInsNucPrior, &firstpos, 0x3c00/divisor );

    if (score < iGapOpen) {
      //printf("*** perfect score found -- bailing out with minimal band\n");
      global_minband = 1;
    }

    //DEBUG
    //for (jj=0;genomealignment[jj];jj++) printf("%c",genomealignment[jj]); printf("\n");
    //for (jj=0;readalignment[jj];jj++) printf("%c",readalignment[jj]); printf("\n");

    // Find perfect matches to seed band
    int x=0, y=0, nummatches = 0, j=0;
    for (j=0; genomealignment[j]; j++) {
      // update length of unbroken match
      if (genomealignment[j] == readalignment[j]) {
	nummatches += 1;
      }
      // if a maximum match found (or the end of alignment is reached), process the match
      if ((genomealignment[j] != readalignment[j]) || (!genomealignment[j+1])) {
	if (nummatches > matchsize) {
	  // update x/y positions if end-of-alignment was the reason for entering this code,
	  // as normally this code is only entered at the first mis-match
	  if (!genomealignment[j+1]) {
	    x += 1;
	    y += 1;
	  }
	  int bandy = firstpos + y - x + i;
	  bandlen[ bandy + band/2 ] += nummatches;
	  //printf ("bandlen[]=%d maxbandlen=%d\n",bandlen[bandy+band/2],maxbandlen);
	  if (bandlen[ bandy + band/2 ] >maxbandlen) {
	    maxbandlen = bandlen[ bandy + band/2 ];
	    spine = bandy;
	    //printf ("spine/maxbandlen now %d %d\n",spine,maxbandlen);
	  }
	  //printf("Found match seed bandy=%d size=%d \n",bandy,nummatches);
	  //int z;
	  //for (z=0; z<nummatches;z++) printf("%c",genomealignment[j-nummatches+z]); printf(" ");
	  //for (z=0; z<nummatches;z++) printf("%c",readalignment[j-nummatches+z]); printf("\n");
	  for (k=1; k<=nummatches; k++) {
	    //DEBUG
	    if ( x-k<0 || x-k >= iLen1 ) {
	      printf("Out of range (1): k=%i x=%i iLen1=%i nummatches=%i\n",k,i,iLen1,nummatches);
	      exit(1);
	    }
	    if (maxlocallen[ x - k ] < nummatches) {
	      maxlocallen[ x - k ] = nummatches;
	      localband[ x - k ] = bandy;
	    }
	  }
	}
	nummatches = 0;
      }
      // update x and y coordinates
      if (genomealignment[j] == readalignment[j]) {
	x += 1;
	y += 1;
      } else {
	if (genomealignment[j] != '-') y += 1;
	if (readalignment[j] != '-') x += 1;
      }
    }
   }
  }

  // Now find band range
  for (k=0; k<iLen1; k++) {
    //printf("Array index=%d locallen=%d band=%d (min/maxband = %d %d)\n",k,maxlocallen[k],localband[k],minband,maxband);
    //printf("spine - globalmaxband = %d  spine + globalmaxband = %d\n",spine-global_maxband/2,spine+global_maxband/2);
    if ((maxlocallen[k] > 0) && (spine-global_maxband/2 < localband[k]) && (localband[k] < spine+global_maxband/2)) {
      if (minband > localband[k]) minband = localband[k];
      if (maxband < localband[k]) maxband = localband[k];
    }
  }

  //printf ("Band range %d - %d  band=%d\n",minband,maxband,band);

  // Calculate minimum and maximum (+1) index into iSeq2 - the segment to be aligned against
  if (minband == band) {
    // No seeds found.  Select small band
    minband = max(-band/2,-global_minband/2);
    maxband = min(band/2,global_minband/2);
  }
  // implement minimum alignment band
  int bandsize = maxband - minband;
  if (bandsize < global_minband) {
    int extra_band = global_minband - bandsize;
    minband = max(minband - extra_band/2, -band/2);
    maxband = min(band/2, maxband + extra_band/2);
  }
  //printf ("Band range %d - %d  band=%d\n",minband,maxband,band);
  int yposdiagonal = (iLen2 - iLen1)/2;
  *yposmin = yposdiagonal + minband;
  int yposmax = yposdiagonal + maxband + iLen1;
  //printf ("yposmin-max = %d - %d\n",*yposmin,yposmax);
  if (*yposmin<0) *yposmin = 0;
  if (yposmax>iLen2) yposmax = iLen2;
  //printf ("after truncating: yposmin-max = %d - %d\n",*yposmin,yposmax);
  
  // Calculate smallest symmetric band that includes the calculated band
  *iNewlen2 = yposmax - *yposmin;
  //printf ("newlen=%d\n",*iNewlen2);
  int abs_newyposdiagonal = (*iNewlen2 - iLen1)/2 + *yposmin;
  //printf ("abs_newyposidagonald=%d\n",abs_newyposdiagonal);
  int newminband = yposdiagonal + minband - abs_newyposdiagonal;
  int newmaxband = yposdiagonal + maxband - abs_newyposdiagonal;
  //printf ("new min/max band (rel to diag) = %d %d\n",newminband,newmaxband);
  *newWidth = max(-2*newminband+1,2*newmaxband+1);

  //printf ("newwidth %d oldband %d\n",*newWidth, band);
  //printf("Aligning\n");

  // Align
  if (*iNewlen2-iLen1>*newWidth) {
    fprintf(stderr,"Problem -- new width too small\n");
  }
  //int jj;
  //for (jj=0;jj<iLen1;jj++) printf("(%c)",iSeq1[jj]); printf("\n");
  //for (jj=0;jj<iLen1;jj++) printf("(%c)",pQuals1[jj]); printf("\n");
  //for (jj=0;jj<*iNewlen2;jj++) printf("(%c)",(iSeq2+*yposmin)[jj]); printf("\n");
  //printf("Exit findband.\n");
}
