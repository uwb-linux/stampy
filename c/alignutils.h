#include "Python.h"

/*****************************************************************************************************************
 This code is copyright (c) Gerton Lunter, Jan 2009
 It may not be distributed, made public, or used in other software without the permission of the copyright holder
******************************************************************************************************************/

// radius of the alignment band; minimum is 1
#define BAND 3

// maximum length of a single sequence (checked at runtime)
#define MAXSEQLEN 255

// to make room for back pointers
#define SCOREFACTOR 4

PyObject* alignprefix_c_rawscore( int* dp, char* seq1, int seq1len, char* seq2, int seq2len, int* scoretable, int gapstart, int band);
PyObject* alignsuffix_c_rawscore( int* dp, char* seq1, int seq1len, char* seq2, int seq2len, int* scoretable, int gapstart, int band);

int traceback_sse2( char* seq1, char* seq2, int* qual2, char* aln1, char* aln2, int len1, int len2, 
		    int gapopen, int gapextend, int nucprior, int* firstpos, int ab_cutoff );

PyObject* align_sensitive( char* seq1, int seq1len, char* seq2, int seq2len, char* pQuals1, int iGapOne, int iGapOpen, int iGapExtend, int maxlikelihood, int mode, int iStartMean, int iStartSD, int band );

