#include "maputils.h"
#include "alignutils.h"
#include "../ext/qsort/qsort.h"

#include <time.h>

/*****************************************************************************************************************
 This code is copyright (c) Gerton Lunter, Jan 2009
 It may not be distributed, made public, or used in other software without the permission of the copyright holder
******************************************************************************************************************/

#ifdef __SSE__

#include "xmmintrin.h"
#define PREFETCH_NTA(address) _mm_prefetch(address,_MM_HINT_NTA)
#define PREFETCH_T2(address) _mm_prefetch(address,_MM_HINT_T2)

#else

// non-SSE2 code is not supported; but keep it to avoid compile errors for ppc platform on Macs

#define PREFETCH_NTA(address)
#define PREFETCH_T2(address)

#endif


#define HIGH_COUNT_FLAG  0xFFFFFFFE
#define LINEAR_HASH_FLAG 0xFFFFFFFF

#define winsize 15     // must be odd
#define UNDUP_CANDIDATES 10
#define MAX_LIKELIHOOD_DIFF 200
#define ALPHA_BETA_CUTOFF 180     // broken pairs are highly unlikely

//static clock_t start, end;
//static long locate_time=0, qsort1_time=0, align_time=0, qsort2_time=0;

unsigned long revcomptable[256] = {0xaa, 0xea, 0x2a, 0x6a, 0xba, 0xfa, 0x3a, 0x7a, 0x8a, 0xca, 0xa, 0x4a, 0x9a, 0xda, 0x1a, 0x5a, 0xae, 0xee, 0x2e, 0x6e, 0xbe, 0xfe, 0x3e, 0x7e, 0x8e, 0xce, 0xe, 0x4e, 0x9e, 0xde, 0x1e, 0x5e, 0xa2, 0xe2, 0x22, 0x62, 0xb2, 0xf2, 0x32, 0x72, 0x82, 0xc2, 0x2, 0x42, 0x92, 0xd2, 0x12, 0x52, 0xa6, 0xe6, 0x26, 0x66, 0xb6, 0xf6, 0x36, 0x76, 0x86, 0xc6, 0x6, 0x46, 0x96, 0xd6, 0x16, 0x56, 0xab, 0xeb, 0x2b, 0x6b, 0xbb, 0xfb, 0x3b, 0x7b, 0x8b, 0xcb, 0xb, 0x4b, 0x9b, 0xdb, 0x1b, 0x5b, 0xaf, 0xef, 0x2f, 0x6f, 0xbf, 0xff, 0x3f, 0x7f, 0x8f, 0xcf, 0xf, 0x4f, 0x9f, 0xdf, 0x1f, 0x5f, 0xa3, 0xe3, 0x23, 0x63, 0xb3, 0xf3, 0x33, 0x73, 0x83, 0xc3, 0x3, 0x43, 0x93, 0xd3, 0x13, 0x53, 0xa7, 0xe7, 0x27, 0x67, 0xb7, 0xf7, 0x37, 0x77, 0x87, 0xc7, 0x7, 0x47, 0x97, 0xd7, 0x17, 0x57, 0xa8, 0xe8, 0x28, 0x68, 0xb8, 0xf8, 0x38, 0x78, 0x88, 0xc8, 0x8, 0x48, 0x98, 0xd8, 0x18, 0x58, 0xac, 0xec, 0x2c, 0x6c, 0xbc, 0xfc, 0x3c, 0x7c, 0x8c, 0xcc, 0xc, 0x4c, 0x9c, 0xdc, 0x1c, 0x5c, 0xa0, 0xe0, 0x20, 0x60, 0xb0, 0xf0, 0x30, 0x70, 0x80, 0xc0, 0x0, 0x40, 0x90, 0xd0, 0x10, 0x50, 0xa4, 0xe4, 0x24, 0x64, 0xb4, 0xf4, 0x34, 0x74, 0x84, 0xc4, 0x4, 0x44, 0x94, 0xd4, 0x14, 0x54, 0xa9, 0xe9, 0x29, 0x69, 0xb9, 0xf9, 0x39, 0x79, 0x89, 0xc9, 0x9, 0x49, 0x99, 0xd9, 0x19, 0x59, 0xad, 0xed, 0x2d, 0x6d, 0xbd, 0xfd, 0x3d, 0x7d, 0x8d, 0xcd, 0xd, 0x4d, 0x9d, 0xdd, 0x1d, 0x5d, 0xa1, 0xe1, 0x21, 0x61, 0xb1, 0xf1, 0x31, 0x71, 0x81, 0xc1, 0x1, 0x41, 0x91, 0xd1, 0x11, 0x51, 0xa5, 0xe5, 0x25, 0x65, 0xb5, 0xf5, 0x35, 0x75, 0x85, 0xc5, 0x5, 0x45, 0x95, 0xd5, 0x15, 0x55};

// 
// Bit patterns to compute a fingerprint with good behaviour under substitutions and indels
// For a 2bit pattern X of 4 nucleotides, fingerprints[x] gives three bytes, counting the
// occurrences of T, C and A (going from LSB to MSB)
//
long fingerprints[256] = {0x4, 0x103, 0x10003, 0x3, 0x103, 0x202, 0x10102, 0x102, 0x10003, 0x10102, 0x20002, 0x10002, 0x3, 0x102, 0x10002, 0x2, 0x103, 0x202, 0x10102, 0x102, 0x202, 0x301, 0x10201, 0x201, 0x10102, 0x10201, 0x20101, 0x10101, 0x102, 0x201, 0x10101, 0x101, 0x10003, 0x10102, 0x20002, 0x10002, 0x10102, 0x10201, 0x20101, 0x10101, 0x20002, 0x20101, 0x30001, 0x20001, 0x10002, 0x10101, 0x20001, 0x10001, 0x3, 0x102, 0x10002, 0x2, 0x102, 0x201, 0x10101, 0x101, 0x10002, 0x10101, 0x20001, 0x10001, 0x2, 0x101, 0x10001, 0x1, 0x103, 0x202, 0x10102, 0x102, 0x202, 0x301, 0x10201, 0x201, 0x10102, 0x10201, 0x20101, 0x10101, 0x102, 0x201, 0x10101, 0x101, 0x202, 0x301, 0x10201, 0x201, 0x301, 0x400, 0x10300, 0x300, 0x10201, 0x10300, 0x20200, 0x10200, 0x201, 0x300, 0x10200, 0x200, 0x10102, 0x10201, 0x20101, 0x10101, 0x10201, 0x10300, 0x20200, 0x10200, 0x20101, 0x20200, 0x30100, 0x20100, 0x10101, 0x10200, 0x20100, 0x10100, 0x102, 0x201, 0x10101, 0x101, 0x201, 0x300, 0x10200, 0x200, 0x10101, 0x10200, 0x20100, 0x10100, 0x101, 0x200, 0x10100, 0x100, 0x10003, 0x10102, 0x20002, 0x10002, 0x10102, 0x10201, 0x20101, 0x10101, 0x20002, 0x20101, 0x30001, 0x20001, 0x10002, 0x10101, 0x20001, 0x10001, 0x10102, 0x10201, 0x20101, 0x10101, 0x10201, 0x10300, 0x20200, 0x10200, 0x20101, 0x20200, 0x30100, 0x20100, 0x10101, 0x10200, 0x20100, 0x10100, 0x20002, 0x20101, 0x30001, 0x20001, 0x20101, 0x20200, 0x30100, 0x20100, 0x30001, 0x30100, 0x40000, 0x30000, 0x20001, 0x20100, 0x30000, 0x20000, 0x10002, 0x10101, 0x20001, 0x10001, 0x10101, 0x10200, 0x20100, 0x10100, 0x20001, 0x20100, 0x30000, 0x20000, 0x10001, 0x10100, 0x20000, 0x10000, 0x3, 0x102, 0x10002, 0x2, 0x102, 0x201, 0x10101, 0x101, 0x10002, 0x10101, 0x20001, 0x10001, 0x2, 0x101, 0x10001, 0x1, 0x102, 0x201, 0x10101, 0x101, 0x201, 0x300, 0x10200, 0x200, 0x10101, 0x10200, 0x20100, 0x10100, 0x101, 0x200, 0x10100, 0x100, 0x10002, 0x10101, 0x20001, 0x10001, 0x10101, 0x10200, 0x20100, 0x10100, 0x20001, 0x20100, 0x30000, 0x20000, 0x10001, 0x10100, 0x20000, 0x10000, 0x2, 0x101, 0x10001, 0x1, 0x101, 0x200, 0x10100, 0x100, 0x10001, 0x10100, 0x20000, 0x10000, 0x1, 0x100, 0x10000, 0x0};

double phredscores[80] = {1.0, 0.79432823472428149, 0.63095734448019325, 0.50118723362727235, 0.39810717055349726, 0.316227766016838, 0.25118864315095807, 0.199526231496888, 0.15848931924611137, 0.12589254117941673, 0.10000000000000002, 0.079432823472428138, 0.063095734448019344, 0.050118723362727227, 0.039810717055349748, 0.031622776601683798, 0.025118864315095805, 0.019952623149688809, 0.015848931924611141, 0.012589254117941677, 0.010000000000000004, 0.007943282347242819, 0.0063095734448019303, 0.0050118723362727264, 0.0039810717055349751, 0.003162277660168382, 0.0025118864315095799, 0.0019952623149688798, 0.0015848931924611149, 0.0012589254117941686, 0.0010000000000000002, 0.00079432823472428175, 0.00063095734448019353, 0.00050118723362727296, 0.00039810717055349779, 0.00031622776601683837, 0.00025118864315095817, 0.00019952623149688793, 0.00015848931924611147, 0.00012589254117941693, 0.00010000000000000009, 7.9432823472428302e-05, 6.3095734448019388e-05, 5.0118723362727326e-05, 3.9810717055349701e-05, 3.1622776601683802e-05, 2.5118864315095835e-05, 1.9952623149688807e-05, 1.5848931924611158e-05, 1.258925411794168e-05, 1.0000000000000016e-05, 7.9432823472428353e-06, 6.3095734448019322e-06, 5.0118723362727275e-06, 3.9810717055349725e-06, 3.1622776601683826e-06, 2.5118864315095852e-06, 1.9952623149688821e-06, 1.5848931924611168e-06, 1.2589254117941689e-06, 1.0000000000000004e-06, 7.9432823472428262e-07, 6.309573444801936e-07, 5.0118723362727303e-07, 3.9810717055349756e-07, 3.1622776601683845e-07, 2.511886431509587e-07, 1.9952623149688832e-07, 1.5848931924611178e-07, 1.2589254117941675e-07, 1.0000000000000029e-07, 7.9432823472428458e-08, 6.3095734448019402e-08, 5.0118723362727342e-08, 3.9810717055349709e-08, 3.1622776601683812e-08, 2.5118864315095841e-08, 1.9952623149688846e-08, 1.5848931924611189e-08, 1.2589254117941683e-08};

// hard coded: 15 nt hash; 5 nt stride; max 1 mismatch; no repeats
// first 5: length=15; etc.
//static int mapstatscacheI[190] = {1,1,0,0,0,2,2,0,0,0,4,4,1,0,0,7,7,1,0,0,100,100,2,0,0,100,100,3,1,0,100,100,4,1,0,100,100,5,2,1,100,100,6,2,1,100,100,8,3,1,100,100,9,4,2,100,100,10,5,2,100,100,12,6,3,100,100,13,7,3,100,100,16,8,4,100,100,17,9,5,100,100,19,11,5,100,100,21,13,6,100,100,24,16,7,100,100,100,100,9,100,100,100,100,10,100,100,100,100,11,100,100,100,100,13,100,100,100,100,15,100,100,100,100,18,100,100,100,100,19,100,100,100,100,21,100,100,100,100,22,100,100,100,100,25,100,100,100,100,31,100,100,100,100,33,100,100,100,100,34,100,100,100,100,36,100,100,100,100,40,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100};

double mapstatscacheD[190] = {0.7943, 0.7943, 1, 1, 1, 0.631, 0.631, 1, 1, 1, 0.3981, 0.3981, 0.7943, 1, 1, 0.1995, 0.1995, 0.7943, 1, 1, 1e-10, 1e-10, 0.631, 1, 1, 1e-10, 1e-10, 0.5012, 0.7943,1, 1e-10, 1e-10, 0.3981, 0.7943, 1, 1e-10, 1e-10, 0.3162, 0.631, 0.7943, 1e-10, 1e-10, 0.2512, 0.631, 0.7943, 1e-10, 1e-10, 0.1585, 0.5012, 0.7943, 1e-10, 1e-10, 0.1259, 0.3981, 0.631, 1e-10, 1e-10, 0.1, 0.3162, 0.631, 1e-10, 1e-10, 0.0631, 0.2512, 0.5012, 1e-10, 1e-10, 0.05012, 0.1995, 0.5012, 1e-10, 1e-10, 0.02512, 0.1585, 0.3981, 1e-10, 1e-10, 0.01995, 0.1259, 0.3162, 1e-10, 1e-10, 0.01259, 0.07943, 0.3162, 1e-10, 1e-10, 0.007943, 0.05012, 0.2512, 1e-10, 1e-10, 0.003981, 0.02512, 0.1995, 1e-10, 1e-10, 1e-10, 1e-10, 0.1259, 1e-10, 1e-10, 1e-10, 1e-10, 0.1, 1e-10, 1e-10, 1e-10, 1e-10, 0.07943, 1e-10, 1e-10, 1e-10, 1e-10, 0.05012, 1e-10, 1e-10, 1e-10, 1e-10, 0.03162, 1e-10, 1e-10, 1e-10, 1e-10, 0.01585, 1e-10, 1e-10, 1e-10, 1e-10, 0.01259, 1e-10, 1e-10, 1e-10, 1e-10, 0.007943, 1e-10, 1e-10, 1e-10, 1e-10, 0.00631, 1e-10, 1e-10, 1e-10, 1e-10, 0.003162, 1e-10, 1e-10, 1e-10, 1e-10, 0.0007943, 1e-10, 1e-10, 1e-10, 1e-10, 0.0005012, 1e-10, 1e-10, 1e-10, 1e-10, 0.0003981, 1e-10, 1e-10, 1e-10, 1e-10, 0.0002512, 1e-10, 1e-10, 1e-10, 1e-10, 0.0001, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10};


//static double db = -0.23025850929940456;   // log(0.1)/10.0;

inline int gu_min(int i, int j) { return (i<j ? i : j); }

inline int gu_max(int i, int j) { return (i>j ? i : j); }

inline double gd_max( double i, double j) { return (i>j ? i : j); }

inline double gd_min( double i, double j) { return (i<j ? i : j); }



u_int32_t inline get_rotated( const char* const data, u_int64_t locus, int rotation ) {

    return ((((u_int64_t)(((u_int32_t*)data)[locus>>4]) << 32) + ((u_int32_t*)data)[1 + (locus>>4)])
    >> (rotation - (locus & 15)*2)) & 0x3FFFFFFF;

    /*
      u_int64_t fifteen0;
      asm (" movq (%1,%2,4), %0; \n\t rorq %%cl, %0;"
      : "=r"(fifteen0)
      : "r"(data), "r"((u_int64_t)(locus >> 4)), "c"( 32 + rotation - 2*(locus & 15) ) );
      return fifteen0;
    */
}


u_int32_t inline get_fifteen( const char* const data, u_int64_t locus ) {
  return get_rotated( data, locus, 34 ) & 0x3FFFFFFF;
}


// computes a fingerprint around locus uloc.  
// For numnucnibbles>0, the fingerprint is computed on the 4*numnucnibbles nucleotides at uloc extending rightward
// For numnucnibbles<0, the fingerprint is computed on the -4*numnucnibbles nucs left of, and not including, uloc
// |numnucnibbles| should be <= 3
inline long fingerprint_fromstring( const char* data, u_int64_t uloc, int numnucnibbles ) {

  if (numnucnibbles == 0) return 0;
  if (numnucnibbles < 0) {
    numnucnibbles = -numnucnibbles;
    uloc -= 4*numnucnibbles;
    // this can cause roll-over.  To avoid this, pad the genome with 16 N bytes and make sure 16*N is always a high count
    // when this is done, remove this fix:
    if (uloc > 0xF0000000ul) uloc = 0;
  }
  unsigned long longword = get_rotated( data, uloc, 32 + 8 );
  long fingerprint = fingerprints[ (longword>>16) & 0xFF ];
  if (numnucnibbles > 1) fingerprint += fingerprints[ (longword>>8) & 0xFF ];
  if (numnucnibbles > 2) fingerprint += fingerprints[ longword & 0xFF ];
  return fingerprint;

}


// computes the fingerprint of a 15-nucleotide word
long inline fingerprint15_fromstring( const char* data, u_int64_t uloc ) {

  //u_int64_t twowords = ((u_int64_t)(((u_int32_t*)data)[uloc>>4]) << 32) + ((u_int32_t*)data)[1 + (uloc>>4)];
  //twowords >>= (32 - (uloc & 15)*2);
  //unsigned long longword = twowords;
  unsigned long longword = get_rotated( data, uloc, 32 );
  return fingerprints[ (longword>>24) & 0xFF ] +
    fingerprints[ (longword>>16) & 0xFF ] +
    fingerprints[ (longword>>8) & 0xFF ] +
    fingerprints[ (longword | 0x03) & 0xFF ];    // make the 16th nucleotide a G, which is not counted
}


long fingerprint( PyObject* mmap, int datastart, u_int64_t locus, int numnucnibbles ) {

  u_int32_t * data = (u_int32_t*) ((mmap_object_head*)mmap)->data;
  u_int64_t size = (u_int64_t)((mmap_object_head*)mmap)->size;
  u_int64_t uloc = locus + 4*datastart;

  if (uloc+16 >= 4*size) {
    return 0;
  }

  return fingerprint_fromstring( (const char*)data, uloc, numnucnibbles );

}
  

char* genome_slice_charptr( PyObject* mmap, int datastart, u_int64_t ustart, u_int64_t uend ) {

  static char convert[4] = "TCAG";

  u_int32_t * data = (u_int32_t*) ((mmap_object_head*)mmap)->data;
  u_int64_t size = (u_int64_t)((mmap_object_head*)mmap)->size;
  u_int64_t remaining = 0;

  ustart += datastart*4;
  uend += datastart*4;

  if (ustart >= uend)
    return NULL;

  if (4*size < uend) {
    remaining = uend - 4*size;
  }

  char* string = malloc( uend - ustart + 1 );
  u_int64_t idx = ustart;
    
  for (; idx < uend-remaining; idx++)
    string[ idx-ustart ] = convert[ ( data[idx>>4] >> (30 - (idx&15)*2)) & 3 ];

  for (; idx < uend; idx++)
    string[ idx-ustart ] = 'N';

  string[ idx-ustart ] = 0;

  return string;
}


void genome_slice_buffer( char* string, PyObject* mmap, int datastart, u_int64_t ustart, u_int64_t uend ) {

  static char convert[4] = "TCAG";

  u_int32_t * data = (u_int32_t*) ((mmap_object_head*)mmap)->data;

  ustart += datastart*4;
  uend += datastart*4;
  u_int64_t idx = ustart;
    
  for (; idx < uend; idx++)
    string[ idx-ustart ] = convert[ ( data[idx>>4] >> (30 - (idx&15)*2)) & 3 ];
  string[ idx-ustart ] = 0;

}


char* sequencetobitstring( char* string, int length, int* bitstringlen ) {

  const int padding = 2;
  *bitstringlen = (padding + (length + 15)/16)*4;
  u_int32_t * bs = (u_int32_t*)malloc( *bitstringlen );
  int i;
  char c;
  u_int32_t cc;
  for (i=0; i<*bitstringlen/4; i++) bs[i] = 0;
  for (i=0; i<*bitstringlen*4; i++) {
    if (i < length) c = string[i]; else c = 'N';
    switch (c & 0xDF) {
    case 'T': cc = 0; break;
    case 'C': cc = 1; break;
    case 'A': cc = 2; break;
    default: cc =3; break;
    }
    bs[ i>>4 ] |= (cc << (30 - (i & 15)*2));
  }
  return (char*)bs;
}


inline unsigned long revcomp( unsigned long x ) {
  return ( ( revcomptable[x & 0xFF] << 24 ) + 
	   ( revcomptable[(x >> 8) & 0xFF] << 16 ) + 
	   (revcomptable[(x >> 16) & 0xFF] << 8) + 
	   revcomptable[ x >> 24 ] ) >> 2;
}


// make a shorter index word, using rev-comp symmetry
inline unsigned long fold( unsigned long x ) {

  const unsigned long halfmask = (1 << winsize) - 1;      // lower half, including one bit of middle 2-bit char 

  if (x & (1<<winsize)) {
    x = revcomp(x);
  }
      
  // remove zero bit
  unsigned long y = x & halfmask;
  return y + ((x-y)>>1);

}


inline long genome_hash_29_from15( u_int32_t hash ) {

  return ( (u_int64_t)( fold(hash) & 0x1FFFFFFF ) * 0x5a36261UL) % 0x1ffffffdUL;

}

inline long genome_hash_29_fromstring( const char* data, u_int64_t uloc ) {

  /* using the inline asm version here doesn't work -- as GCC optimizes memory writes to data[] away */
  u_int64_t twowords = ((u_int64_t)(((u_int32_t*)data)[uloc>>4]) << 32) + ((u_int32_t*)data)[1 + (uloc>>4)];
  twowords >>= (34 - (uloc & 15)*2);
  u_int32_t hash = ((u_int32_t)twowords) & 0x3FFFFFFF;
  //u_int32_t hash = get_fifteen( data, uloc );  /* Doesn't work */

  // remove reverse-complement symmetry to project the hash to 29 bits, and
  // pseudo-randomize to get even loading of hash table -- p = 0x1fff fffd is a prime,
  // and 0x5a36261 is relative prime to p-1 and thus generates Z / pZ, and is
  // a nice 'random' bit pattern 
  return genome_hash_29_from15( hash );

  /*
  return ( (u_int64_t)( fold(hash) & 0x1FFFFFFF ) * 0x5a36261UL) % 0x1ffffffdUL;
  */

}


inline int strand_bit_29_fromstring( const char* data, u_int64_t uloc ) {

  u_int64_t twowords = ((u_int64_t)(((u_int32_t*)data)[uloc>>4]) << 32) + ((u_int32_t*)data)[1 + (uloc>>4)];
  twowords >>= (34 + winsize - (uloc & 15)*2);

  return ((u_int32_t)twowords) & 1;

}


int read_nuc_fromstring( char* data, u_int64_t uloc ) {

  u_int64_t twowords = ((u_int64_t)(((u_int32_t*)data)[uloc>>4]) << 32) + ((u_int32_t*)data)[1 + (uloc>>4)];
  twowords >>= (34 - (uloc & 15)*2);
  return (int) (((unsigned int)twowords) >> 28) & 3;

}

void set_nuc_fromstring_inplace( char* data, u_int64_t uloc, int nuc ) {

  u_int32_t mask = ~ (3 << (30 - (uloc & 15)*2));
  u_int32_t bitnuc = nuc << (30 - (uloc & 15)*2);
  ((u_int32_t*)data)[uloc>>4] &= mask;
  ((u_int32_t*)data)[uloc>>4] |= bitnuc;

}

void set_bit_fromstring_inplace( char* data, u_int64_t uloc, int state ) {

  u_int32_t mask = 1 << (31-(uloc & 31));
  if (state == 0) {
    ((u_int32_t*)data)[uloc>>5] &= ~mask;
  } else {
    ((u_int32_t*)data)[uloc>>5] |= mask;
  }

}

int read_bit_fromstring( char* data, u_int64_t uloc ) {

  return (int) (((u_int32_t*)data)[uloc>>5] >> (31-(uloc & 31))) & 1;

}

  

// accessing the genome sequence

inline char revcompchar( char in ) {
  switch (in) {
  case 'A': return 'T';
  case 'T': return 'A';
  case 'C': return 'G';
  case 'G': return 'C';
  };
  return 'N';
}


void revcompinplace( char * string ) {

  int len = strlen( string );
  int i=0, j=len-1;
  while (i <= j) {
    char n1 = revcompchar( string[i] );
    char n2 = revcompchar( string[j] );  // may be the same
    string[j] = n1;
    string[i] = n2;
    i += 1;
    j -= 1;
  }
}


char* getgenome( PyObject* mmap, int datastart, int64_t start, int64_t end, int64_t maxlen ) {

  if (start >= 0) {
    if (end > maxlen) {
      end = maxlen;
    }
    if (end < 0)
      end = 0;
    return genome_slice_charptr( mmap, datastart, start, end );
  }

  if (end <= 0) {
    if (start < -maxlen) {
      start = -maxlen;
    }
    if (start > -1) {
      start = -1;
    }
    
    char * string = genome_slice_charptr( mmap, datastart, -end+1, -start+1 );
    revcompinplace( string );
    return string;

  }

  // out of range
  char * string = malloc(1);
  string[0] = 0;
  return string;

}



void getgenome_buffer( char* string, PyObject* mmap, int datastart, int64_t start, int64_t end, int64_t maxlen ) {

  if (start >= 0) {
    if (end > maxlen) {
      //end = maxlen;
      string[0]=0; return;
    }
    if (end < 0) {
      //end = 0;
      string[0]=0; return;
    }
    genome_slice_buffer( string, mmap, datastart, start, end );
    return;
  }

  if (end <= 0) {
    if (start < -maxlen) {
      //start = -maxlen;
      string[0]=0; return;
    }
    if (start > -1) {
      //start = -1;
      string[0]=0; return;
    }
    
    genome_slice_buffer( string, mmap, datastart, -end+1, -start+1 );
    // revcompinplace
    int i=0, j=end-start-1;
    while (i <= j) {
      char n1 = revcompchar( string[i] );
      char n2 = revcompchar( string[j] );  // may be the same
      string[j] = n1;
      string[i] = n2;
      i += 1;
      j -= 1;
    }
    return;

  }
  // out of range
  string[0] = 0;
}


//DEBUG
//#define DEBUG0(s) fprintf(stdout,s);
//#define DEBUG1(s,v) fprintf(stdout,s,v);
#define DEBUG0(s) ;
#define DEBUG1(s,v) ;



void inline addmapping( struct mapping** mappings, int* nummappings, int* maxnummappings, 
			int64_t locus, long start, long end ) {

  if (*nummappings == *maxnummappings) {
    // double the space
    *maxnummappings *= 2;
    *mappings = (struct mapping*)realloc( (void*)*mappings, sizeof(struct mapping)*(*maxnummappings) );
  }

  //DEBUG1("Adding mapping for locus %li\n",locus);
  //DEBUG1("  (start = %li)\n", start)

  (*mappings)[*nummappings].readalignment = NULL;     // mark as unaligned
  (*mappings)[*nummappings].genomealignment = NULL;   // mark as unaligned
  (*mappings)[*nummappings].locus = locus;
  (*mappings)[*nummappings].start = start;
  (*mappings)[*nummappings].end = end;

  (*nummappings)++;

}



//
// this one uses the 3rd hash table format,
//
void addcandidateloci_fingerprint_hash3( long hashhit, char* seqbitstring, int offset, struct mapping** mappings, 
					 int* nummappings, int* maxnummappings, 
					 u_int32_t* hashmap, PyObject* genomemmap, int datastart, long mask,
					 unsigned long fingerprint, int leftnibs, int rightnibs, int maxdiff ) {

  const unsigned int hashfootprint = winsize;
  const unsigned long stride = 5;
  long hashptr = hashhit;

  // so that we can call genome_hash_29_fromstring, rather than genome_hash_29
  const char* const data = ((mmap_object_head*)genomemmap)->data + datastart;
  unsigned long matchfp;
  int t,c,a,dt,dc,da, mindiff;

  // get the bitstring for the 15-nuc sequence and its reverse instead
  const u_int64_t twowords = (((u_int64_t)(((u_int32_t*)seqbitstring)[offset>>4]) << 32) + ((u_int32_t*)seqbitstring)[1 + (offset>>4)]) >> (34 - (offset & 15)*2);
  const u_int32_t fwword = ((u_int32_t)twowords) & 0x3FFFFFFF;
  const u_int32_t bwword = revcomp( fwword );

  //if (maxdiff==3) printf("ENTERING offset=%i hashptr=%li\n",offset,hashptr);

// hash3
//#define ADVANCEHASH( hashptr, hashhit, bucketentry ) hashptr = ((hashptr+hashhit) % 0x1FFFFFFD) & mask
// linear
//#define ADVANCEHASH( hashptr, hashhit, bucketentry ) hashptr = ((hashptr + 1) % 0x1FFFFFFD) & mask
// quadratic
//#define ADVANCEHASH( hashptr, hashhit, bucketentry ) hashptr = ((hashptr + (++bucketentry)) % 0x1FFFFFFD) & mask
// subquadratic
//#define ADVANCEHASH( hashptr, hashhit, bucketentry ) hashptr = ((hashptr + 1 + (++bucketentry)/2) % 0x1FFFFFFD) & mask
// linquad -- at cutoff 50 it makes things 1% faster, probably a little more when running many jobs together
#define ADVANCEHASH( hashptr, hashhit, bucketentry ) hashptr = ((hashptr + (bucketentry==-1 ? 1 : ++bucketentry)) % 0x1FFFFFFD) & mask

  long bucketentry = 0;

  //DEBUG1("Entering with hashptr %lu\n", hashptr)

    // Use for LINQUAD hash
    if (hashmap[hashptr] == LINEAR_HASH_FLAG) {
      bucketentry = -1;
      hashptr = ((hashptr + 1) % 0x1FFFFFFD) & mask;   // advance to next entry
    }

  u_int64_t previous_locus = 0;

  do {

    // read from hash
    u_int32_t entry = hashmap[ hashptr ];
    //DEBUG1(" Hashptr now %lu\n",hashptr);
    //DEBUG1(" Entry now %u\n",entry);
    //DEBUG1(" Low 2bits %u\n",entry&3);

    // calculate next hash ptr
    ADVANCEHASH( hashptr, hashhit, bucketentry );

    // compute actual locus
    u_int64_t locus = (entry >> 2) * stride;

    if (entry & 1) {
      PREFETCH_NTA( &hashmap[ hashptr ] );
    }

    // this should not happen
    if (entry == 0) {
      printf("Internal error: Zero entry found in hash table\n");
      break;
    }
    
    // skip high count flags (belonging to other hashes)
    // also, skip entries with locus left of previous locus, which must belong to other hashes since loci are entered in-order.
    if (entry == HIGH_COUNT_FLAG || entry == LINEAR_HASH_FLAG || locus < previous_locus) {
      continue;
    }

    u_int32_t fifteen = get_fifteen( data, locus );

    // get strand, compute locus and strand of the read, and compute fingerprints;
    // continue if the hash entry belongs to another chain
    int64_t signedlocus;
    //DEBUG1(" locus %lli\n",locus);
    //DEBUG1(" fifteen %u\n",fifteen);
    //DEBUG1(" fwword %u\n",fwword);
    //DEBUG1(" bwword %u\n",bwword);
      
    if (fifteen == bwword) {
      signedlocus = -(int64_t)(locus + hashfootprint - 1 + offset);
      // get counts in the reverse-complement direction
      matchfp = (unsigned long)(fingerprint_fromstring( data, locus, -rightnibs ) + 
				fingerprint_fromstring( data, locus + hashfootprint, leftnibs ));
      // reverse complement the counts; from LSB: T, C, A
      t = matchfp & 0xFF;
      c = (matchfp>>8) & 0xFF;
      a = (matchfp>>16) & 0xFF;
      matchfp = a + (((leftnibs+rightnibs)*4 - (a+c+t))<<8) + (t<<16);
    } else if (fifteen == fwword) {
      // when mapping long reads, occassionally loci walk left off the genome
      signedlocus = (int64_t)locus - (int64_t)offset;
      if (signedlocus<0) signedlocus=0;
      matchfp = (unsigned long)(fingerprint_fromstring( data, locus, -leftnibs ) + 
				fingerprint_fromstring( data, locus + hashfootprint, rightnibs ));
    } else {
      // it may happen that the sequence doesn't match, but the hash does.
      // this is particularly true if the mask isn't full i.e. 0x1fff ffff, but also happens for full mask
      // Check whether this entry ends the chain
      if ((entry & 1)==0 && (genome_hash_29_from15( fifteen ) & mask) == hashhit) break;
      DEBUG0(" Skipping entry\n");
      continue;
    }

    // valid locus
    previous_locus = locus;
    
    // compare fingerprints
    matchfp += 0x101010 - fingerprint;
    dt = (matchfp & 0xFF) - 0x10;
    dc = ((matchfp>>8) & 0xFF) - 0x10;
    da = ((matchfp>>16) & 0xFF) - 0x10;
    
    // the minimum amount of substitutions and 1nt indels compatible with these fingerprints
    mindiff = gu_max(dt,0) + gu_max(dc,0) + gu_max(da,0) + gu_max( -dt-dc-da, 0 );

    if (mindiff <= maxdiff) {
      addmapping( mappings, nummappings, maxnummappings, 
		  signedlocus, offset, offset+hashfootprint );
    } else {
      DEBUG1(" not added; mindiff=%u\n",mindiff);
    }

    // is this the last entry of the chain
    if ((entry & 1) == 0) {
      DEBUG0(" End of chain\n");
      break;
    }
  
    // it is not - advance and check the next one
  } while (1);
  DEBUG1("Done with chain %i\n",0);
}



void extern inline addcandidateloci_nofingerprint_hash3( long hashhit, char* seqbitstring, int offset, struct mapping** mappings, 
						  int* nummappings, int* maxnummappings, 
						  u_int32_t* hashmap, PyObject* genomemmap, int datastart, long mask ) {

  const unsigned int hashfootprint = winsize;
  const unsigned long stride = 5;
  long hashptr = hashhit;

  // so that we can call genome_hash_29_fromstring, rather than genome_hash_29
  const char* const data = ((mmap_object_head*)genomemmap)->data + datastart;

  // get the bitstring for the 15-nuc sequence and its reverse instead
  const u_int64_t twowords = (((u_int64_t)(((u_int32_t*)seqbitstring)[offset>>4]) << 32) + ((u_int32_t*)seqbitstring)[1 + (offset>>4)]) >> (34 - (offset & 15)*2);
  const u_int32_t fwword = ((u_int32_t)twowords) & 0x3FFFFFFF;
  const u_int32_t bwword = revcomp( fwword );

// quadratic
//#define ADVANCEHASH( hashptr, hashhit, bucketentry ) hashptr = ((hashptr + (++bucketentry)) % 0x1FFFFFFD) & mask
// linquad
#define ADVANCEHASH( hashptr, hashhit, bucketentry ) hashptr = ((hashptr + (bucketentry==-1 ? 1 : ++bucketentry)) % 0x1FFFFFFD) & mask

  long bucketentry = 0;

  DEBUG1("Entering (nofingerprint) with hashptr %lu\n", hashptr)

  // Use for LINQUAD hash
  if (hashmap[hashptr] == LINEAR_HASH_FLAG) {
    bucketentry = -1;
    hashptr = ((hashptr + 1) % 0x1FFFFFFD) & mask;   // advance to next entry
  }

  do {

    // read from hash
    u_int32_t entry = hashmap[ hashptr ];

    // calculate next hash ptr
    ADVANCEHASH( hashptr, hashhit, bucketentry );

    // compute actual locus
    u_int64_t locus = (entry >> 2) * stride;

    // this should not happen
    if (entry == 0) {
      printf("Internal error: Zero entry found in hash table\n");
      break;
    }
    
    // skip high count flags (belonging to other hashes)
    if (entry == HIGH_COUNT_FLAG || entry == LINEAR_HASH_FLAG) continue;

    u_int32_t fifteen = get_fifteen( data, locus );

    // get strand, compute locus and strand of the read;
    // continue if the hash entry belongs to another chain
    int64_t signedlocus;
    if (fifteen == bwword) {
      signedlocus = -(int64_t)(locus + hashfootprint - 1 + offset);
    } else if (fifteen == fwword) {
      signedlocus = locus - offset;
    } else {
      // it may happen that the sequence doesn't match, but the hash does.
      // this is particularly true if the mask isn't full i.e. 0x1fff ffff, but also happens for full mask
      // Check whether this entry ends the chain
      if ((entry & 1)==0 && (genome_hash_29_from15( fifteen ) & mask) == hashhit) break;
      continue;
    }
    
    addmapping( mappings, nummappings, maxnummappings, 
		signedlocus, offset, offset+hashfootprint );

    // is this the last entry of the chain
    if ((entry & 1) == 0) break;
  
    // it is not - advance and check the next one
  } while (1);
}





void locate_candidates_fingerprint_hash3( struct mapping** mappings, int* nummappings, int* maxnummappings, int* directhighcounts, int* indirecthighcounts,
					  char* sequence, char* seqbitstring, int seqbitstringlen, PyObject* hashtablemmap, PyObject* genomemmap, 
					  int datastart, long mask, char tryvariants, int maxfingerprintvariants, int maxdirecthighcounts, 
					  int leftreliable, int rightreliable, int* intquality, int* repeatmask) {

  assert( maxfingerprintvariants >= 1 );

  const int hashfootprint = winsize;
  unsigned const int stride = 5;
  u_int32_t * hashmap = (u_int32_t*) ((mmap_object_head*)hashtablemmap)->data;
  u_int32_t entry;
  int i,j,k;
  char foundhighcount;
  int maxidx = strlen(sequence) - hashfootprint + 1;
  unsigned long fingerprint[maxidx];
  int leftnibs, rightnibs, ncount[maxidx];

  // all of the hash footprint must fall in reliable sequence, for repeat counters to be updated
  rightreliable -= hashfootprint;

  // if autoselect, choose tryvariant setting depending on sequence length
  // maxidx = 5  seqlen = 19: 1 x 5nt
  // maxidx = 10 seqlen = 24: 2 x 5nt
  // maxidx = 15 seqlen = 29: 3 x 5nt
  // maxidx = 20 seqlen = 34: 4 x 5nt
  // maxidx = 25 seqlen = 39: 3 x 5nt
  // maxidx = 35 seqlen = 49: 4 x 5nt
  // maxidx = 40 seqlen = 54: 3 x 5nt
  // maxidx = 50 seqlen = 64: 4 x 5nt
  // maxidx = 65 seqlen = 79: 5 x 5nt
  //DEBUG1("tryvariants set: %i\n",tryvariants);

  if (tryvariants == -1) tryvariants = (maxidx <= 20 ? 1 : (maxidx <= 35 ? 2 : 3));

  //DEBUG1("tryvariants now %i\n",tryvariants);

  // first do not consider variants
  for (i=0; i<maxidx; i++) {
    
    //DEBUG1("Now trying index=%i\n",i);

    // count number of Ns
    ncount[i] = 0;
    for (j=0; j<hashfootprint; j++) ncount[i] += (sequence[i+j] == 'N');

    // compute fingerprint
    leftnibs = gu_min(i/4,3);
    rightnibs = gu_min( 3-leftnibs, (maxidx + hashfootprint - 1 - 4*leftnibs)/4 );
    fingerprint[i] = (unsigned long)(fingerprint_fromstring( seqbitstring, i, -leftnibs ) + 
				     fingerprint_fromstring( seqbitstring, i+hashfootprint, rightnibs ));

    if (ncount[i] == 0) {

      // standard case
      long h = genome_hash_29_fromstring( seqbitstring, i ) & mask;
      entry = hashmap[h];

      if (entry == HIGH_COUNT_FLAG) {
	repeatmask[i] = 0;
	if ((leftreliable<=i) && (i<=rightreliable)) {
	  (*directhighcounts)++;
	}
      } else if (entry & 2) {
	if (maxfingerprintvariants == -1) {
	  addcandidateloci_nofingerprint_hash3(h, seqbitstring, i, mappings, nummappings, maxnummappings, 
					       hashmap, genomemmap, datastart, mask );
	} else {
	  addcandidateloci_fingerprint_hash3( h, seqbitstring, i, mappings, nummappings, maxnummappings, 
					      hashmap, genomemmap, datastart, mask, 
					      fingerprint[i], leftnibs, rightnibs, maxfingerprintvariants );
	}
      }

    } else if (ncount[i] == 1) {

      for (j=0; sequence[i+j] != 'N'; j++) {}
      
      foundhighcount = 0;

      for (k=0; k<4; k++) {

	set_nuc_fromstring_inplace( seqbitstring, i+j, k);
	long h = genome_hash_29_fromstring( seqbitstring, i ) & mask;
	entry = hashmap[h];
	if (entry == HIGH_COUNT_FLAG) {
	  repeatmask[i] = 0;
	  if (!foundhighcount) {
	    if ((leftreliable<=i) && (i<=rightreliable))
	      (*directhighcounts)++;
	    foundhighcount = 1;
	  }
	} else if (entry & 2) {
	  if (maxfingerprintvariants == -1) {
	    addcandidateloci_nofingerprint_hash3(h, seqbitstring, i, mappings, nummappings, maxnummappings, 
						 hashmap, genomemmap, datastart, mask );
	  } else {
	    addcandidateloci_fingerprint_hash3( h, seqbitstring, i, mappings, nummappings, maxnummappings, 
						hashmap, genomemmap, datastart, mask, 
						fingerprint[i], leftnibs, rightnibs, maxfingerprintvariants );
	  }
	}
      }
      // leave it in state 3 = G, which codes for N
    }
  }

  // do we need to consider variants?

  // Change in v1.0.20: trying variants is no longer dependent on being in a repeat.
  // This improves the specificity on human data: fewer (0.5%) wrong high-mapq maps in repetitive regions, 
  //  takes 5% longer to run
  
  //if (tryvariants && (*directhighcounts < maxdirecthighcounts)) {
  if (tryvariants) {

    for (i=0; i<maxidx; i++) {

      // skip positions 5-9, 15-19 etc. for speed, if tryvariants == 2
      // skip positions 5-14, 20-29 etc. for speed, if tryvariants == 3
      if (ncount[i]==0 && ((i/stride)%tryvariants) == 0) {

	leftnibs = gu_min(i/4,3);
	rightnibs = gu_min( 3-leftnibs, (maxidx + hashfootprint - 1 - 4*leftnibs)/4 );
	foundhighcount = 0;
	
	//DEBUG1("Now considering VARIANTS at idx %i\n",i);

	for (j=0; j<hashfootprint; j++) {

	  int original = read_nuc_fromstring( seqbitstring, i+j );
	  
	  for (k=0; k<4; k++) {

	    if (k != original) {
	      
	      set_nuc_fromstring_inplace( seqbitstring, i+j, k );
	      long h = genome_hash_29_fromstring( seqbitstring, i ) & mask;
	      entry = hashmap[h];

	      if (entry == HIGH_COUNT_FLAG) {
		repeatmask[i] = gu_min( repeatmask[i], intquality[i+j] );
		if ((leftreliable<=i) && (i<=rightreliable) && (!foundhighcount)) {
		  (*indirecthighcounts)++;
		  foundhighcount = 1;
		}
	      } else if (entry & 2) {
		if (maxfingerprintvariants == -1) {
		  addcandidateloci_nofingerprint_hash3(h, seqbitstring, i, mappings, nummappings, maxnummappings, 
						       hashmap, genomemmap, datastart, mask );
		} else {
		  addcandidateloci_fingerprint_hash3( h, seqbitstring, i, mappings, nummappings, maxnummappings, 
						      hashmap, genomemmap, datastart, mask, 
						      fingerprint[i], leftnibs, rightnibs, maxfingerprintvariants );
		}
	      }
	    }
	  }

	  // revert to original sequence
	  set_nuc_fromstring_inplace( seqbitstring, i+j, original );
	    
	}

      }

    }

  }

}



void align_candidate_c( int64_t locus, int start, int end, char* sequence, int overhang, PyObject* mmap, 
			int datastart, int64_t maxlen, int gapstartpenalty, int* scoretable,
			char** alignment0, char** alignment1, int64_t* alignlocus, int band ) {


  const int hashoverlap = 7;             // amount of overlap of alignment with hash window; must be smaller than winsize/2
  int readlen = strlen( sequence );
  char genomesequence[ 2*MAXSEQLEN ];    // must also fit overhang

  // If no alignment required, just return a dummy alignment pair
  if (overhang == 0) {
    *alignment0 = NULL;
    *alignment1 = getgenome( mmap, datastart, locus, locus + readlen, maxlen );
    if (*alignment1 == NULL)  return;
    *alignment0 = malloc( readlen );
    strncpy( *alignment0, sequence, readlen+1 );
    *alignlocus = locus;
    return;
  }
  if (overhang<0) {
    printf("Error: align_candidate_c: negative overhang.  Program will probably crash now.\n");
    *alignment0 = NULL;
    *alignment1 = NULL;
    return; // error
  }

  if (readlen + 2*overhang + 1 > 2*MAXSEQLEN) {
    printf ("Read too long; maximum is %u -- aborting\n",2*MAXSEQLEN - 1 - 2*overhang);
    exit(1);
  }

  // Prepare to align the left and right overhang from the hash footprint
  getgenome_buffer( genomesequence, mmap, datastart, locus - overhang, locus + readlen + overhang, maxlen );

  if (genomesequence[0] == 0) {
    *alignment0 = NULL;
    *alignment1 = NULL;
    return;
  }

  // prefix of read sequence
  int dp[ (2*MAXSEQLEN+1)*(2*MAXSEQLEN+1) ];

  PyObject* alnl = alignsuffix_c_rawscore( dp, sequence, start + hashoverlap, genomesequence, start + hashoverlap + overhang, scoretable, gapstartpenalty, band );
  PyObject* alnr = alignprefix_c_rawscore( dp, sequence + end - hashoverlap, readlen - end + hashoverlap, genomesequence + overhang + end - hashoverlap, readlen + overhang + hashoverlap - end, scoretable, gapstartpenalty, band );

  // retrieve the partial alignments
  char* alnl0 = PyString_AsString( PyTuple_GetItem( alnl, 0 ) );
  char* alnl1 = PyString_AsString( PyTuple_GetItem( alnl, 1 ) );
  long alnlpos = PyInt_AsLong( PyTuple_GetItem( alnl, 2 ) );

  char* alnr0 = PyString_AsString( PyTuple_GetItem( alnr, 0 ) );
  char* alnr1 = PyString_AsString( PyTuple_GetItem( alnr, 1 ) );

  // compute the alignment's left locus
  int64_t genomepos = (start + hashoverlap + overhang - 1 - alnlpos) + locus - overhang;

  // glue the partials together
  int alnl0len = strlen( alnl0 );
  int alnr0len = strlen( alnr0 );
  int midlen = end - start - 2*hashoverlap;
  char* aln0 = malloc( alnl0len + alnr0len + midlen + 1 );
  char* aln1 = malloc( alnl0len + alnr0len + midlen + 1 );
  strncpy( aln0, alnl0, alnl0len );
  strncpy( aln1, alnl1, alnl0len );
  strncpy( aln0 + alnl0len, sequence + start + hashoverlap, midlen );
  strncpy( aln1 + alnl0len, genomesequence + start + overhang + hashoverlap, midlen );
  strncpy( aln0 + alnl0len + midlen, alnr0, alnr0len+1 );
  strncpy( aln1 + alnl0len + midlen, alnr1, alnr0len+1 );

  // deallocate
  Py_DECREF( alnl );
  Py_DECREF( alnr );
  // build return value
  *alignment0 = aln0;
  *alignment1 = aln1;
  *alignlocus = genomepos;

}



// Nucleotides aligned to gaps are scored with 6 additional units (1/4; equilibrium distribution)
// Ns are not penalized this way
int score_alignment( char* readaln, char* genomealn, PyObject* quality, int gapopen, int gapextend ) {

  int idx = 0;
  int readpos = 0;
  int score = 0;
  char gap = 0;
  int alnlen = strlen( readaln );
  PyObject* qualseq = PySequence_Fast( quality, "Expected sequence of quality values" );
  while (idx < alnlen) {
    char readnuc = readaln[idx];
    char genomenuc = genomealn[idx];
    idx += 1;
    if (gap) {
      if (readnuc == '-') {
	score += gapextend;
	if (genomenuc != 'N') score += 6;
	continue;
      }
      if (genomenuc == '-') {
	score += gapextend;
	if (readnuc != 'N') score += 6;
	readpos += 1;
	continue;
      }
      gap = 0;
    }

    if (readnuc != genomenuc) {
      if (readnuc == '-' || genomenuc == '-') {
	// gap
	score += gapopen;
	// do not penalize double gaps, because the aligner uses linear gap penalties
	gapopen = gapextend;
	if (readnuc != 'N' && genomenuc != 'N')
	  score += 6;
	gap = 1;
      } else {
	// two nucleotides; not identical
	if (readnuc != 'N' && genomenuc != 'N') 
	  score += PyInt_AsLong( PySequence_Fast_GET_ITEM( qualseq, readpos ) );
	else
	  score += 6;
      }
    }

    if (readnuc != '-') readpos += 1;

  }
  Py_DECREF( qualseq );
  return score;

}



int score_alignment_c( char* readaln, char* genomealn, int* quality, int gapopen, int gapextend, int nucprior ) {

  int idx = 0;
  int readpos = 0;
  int score = 0;
  char gap = 0;
  char readnuc;
  while ( (readnuc = readaln[idx]) ) {
    char genomenuc = genomealn[idx];
    idx += 1;
    if (gap) {
      if (readnuc == '-') {
	score += gapextend;
	continue;
      }
      if (genomenuc == '-') {
	score += gapextend + nucprior;
	readpos += 1;
	continue;
      }
      gap = 0;
    }

    if (readnuc != genomenuc) {
      if (readnuc == '-' || genomenuc == '-') {
	score += gapopen;
	gap = 1;
      } else {
	if (readnuc != 'N' && genomenuc != 'N') score += quality[readpos];
      }
    }

    if (readnuc != '-') readpos += 1;

  }
  return score;

}





///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Python / Cython glue code
//
///////////////////////////////////////////////////////////////////////////////////////////////////
PyObject* rev_comp( char* string, int length ) {

  char rc[length+1];
  int j = length;
  int i=0;
  char c;
  rc[j] = 0;
  for (; i<length; i++) {
    switch (string[i]) {
    case 'a': c = 't'; break;
    case 'c': c = 'g'; break;
    case 't': c = 'a'; break;
    case 'g': c = 'c'; break;
    case 'A': c = 'T'; break;
    case 'C': c = 'G'; break;
    case 'T': c = 'A'; break;
    case 'G': c = 'C'; break;
    case 'n': c = 'n'; break;
    case 'N': c = 'N'; break;
    case '-': c = '-'; break;
    default: c = 'N';
    }
    rc[--j] = c;
  }
  PyObject* ret = PyString_FromStringAndSize( rc, length );
  return ret;
}


PyObject* genome_slice( PyObject* mmap, int datastart, u_int64_t ustart, u_int64_t uend ) {

  char* string = genome_slice_charptr( mmap, datastart, ustart, uend );
  if (string == NULL) {
    Py_INCREF(Py_None);
    return Py_None;
  }
  PyObject* ret = PyString_FromStringAndSize( string, uend-ustart );
  free( string );
  return ret;

}
  

// computes a hash for 15 nucleotides, into a 29 bit long
long genome_hash_29( PyObject* mmap, int datastart, u_int64_t locus ) {
  
  u_int32_t * data = (u_int32_t*) ((mmap_object_head*)mmap)->data;
  u_int64_t size = (u_int64_t)((mmap_object_head*)mmap)->size;
  u_int64_t uloc = locus + 4*datastart;

  if (uloc+16 >= 4*size) {
    return 0;
  }

  return genome_hash_29_fromstring( (char*)data, uloc );

}

// computes a hash for 15 nucleotides, into a 29 bit long
int strand_bit_29( PyObject* mmap, int datastart, u_int64_t locus ) {
  
  u_int32_t * data = (u_int32_t*) ((mmap_object_head*)mmap)->data;
  u_int64_t size = (u_int64_t)((mmap_object_head*)mmap)->size;
  u_int64_t uloc = locus + 4*datastart;

  if (uloc+16 >= 4*size) {
    return 0;
  }

  return strand_bit_29_fromstring( (char*)data, uloc );

}


PyObject* set_nuc_fromstring( char* data, int size, u_int64_t uloc, int nuc ) {

  char* newstring = malloc( size );
  memcpy( newstring, data, size );
  set_nuc_fromstring_inplace( newstring, uloc, nuc );
  PyObject* ret = PyString_FromStringAndSize( newstring, size );
  free( newstring );
  return ret;

}



void set_bit_mmap( PyObject* mmap, unsigned long offset, u_int64_t uloc, int bitstate ) {

  char* data = ((mmap_object_head*)mmap)->data;
  size_t size = ((mmap_object_head*)mmap)->size;

  if (offset + (uloc >> 3) > size) return;

  set_bit_fromstring_inplace( data+offset, uloc, bitstate );

}



PyObject* read_unsigned_long_frommmap( PyObject* mmap, int index ) {

  u_int32_t * data = (u_int32_t*) ((mmap_object_head*)mmap)->data;
  return PyLong_FromUnsignedLong( data[ index ] );

}




PyObject* align_candidate( int64_t locus, int start, int end, char* sequence, int overhang, PyObject* mmap, 
			   int datastart, int64_t maxlen, int gapstartpenalty, PyObject* pyscoretable, int band ) {

  char* aln1;
  char* aln2;
  int64_t alignlocus;

  // convert python list of integers to integer array
  int i;
  int scoretable[ 5*5 ];
  for (i=0; i< 25; i++) {
    scoretable[i] = SCOREFACTOR * PyInt_AsLong( PyList_GetItem( pyscoretable, i ) );
  }

  align_candidate_c( locus, start, end, sequence, overhang, mmap, datastart, maxlen, gapstartpenalty, scoretable,
		     &aln1, &aln2, &alignlocus, band );

  if (aln1 == NULL) {
    Py_INCREF( Py_None );
    return Py_None;
  }

  PyObject* align0 = PyString_FromString( aln1 );
  PyObject* align1 = PyString_FromString( aln2 );
  PyObject* gp = PyLong_FromLongLong( alignlocus );
  PyObject* st = PyInt_FromLong( start );
  PyObject* nd = PyInt_FromLong( end );
  PyObject* ret = PyTuple_New(5);
  PyTuple_SET_ITEM(ret, 0, align0);      
  PyTuple_SET_ITEM(ret, 1, align1);      
  PyTuple_SET_ITEM(ret, 2, gp);      
  PyTuple_SET_ITEM(ret, 3, st);      
  PyTuple_SET_ITEM(ret, 4, nd);      
  free (aln1 );
  free (aln2 );
  // return
  return ret;

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Implementation of mapper.mapper
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void score_alignments( struct mapping* mappings, int nummappings, PyObject* quality, int gapopen, int gapextend ) {

  int i=0;
  for (; i<nummappings; i++) {
    mappings[i].score = score_alignment( mappings[i].readalignment, mappings[i].genomealignment, quality, gapopen, gapextend );
  }
}



int compare_mapping_locus_then_start_ptr( const void* mapping1, const void* mapping2 ) {

  int64_t dlocus = (*(const struct mapping**)mapping1)->locus - (*(const struct mapping**)mapping2)->locus;
  if (dlocus<0) return -1;
  if (dlocus>0) return 1;
  return (*(const struct mapping**)mapping1)->start - (*(const struct mapping**)mapping2)->start;
}

int compare_mapping_score_then_locus_ptr( const void* mapping1, const void* mapping2 ) {

  int dscore = (*(const struct mapping**)mapping1)->score - (*(const struct mapping**)mapping2)->score;
  if (dscore!=0) return dscore;
  int64_t dlocus = (*(const struct mapping**)mapping1)->locus - (*(const struct mapping**)mapping2)->locus;
  if (dlocus<0) return -1;
  if (dlocus>0) return 1;
  return 0;
}



void unique_candidates_ptr( struct mapping** mappings, int* nummappings, int maxdistance, int minhashmatches ) {

  struct mapping** writer = mappings;
  int written = 0;
  int i=0;
  int maxhashmatches = 0;

  //DEBUG1("** Uniqueify %u mappings\n",*nummappings);

  while (i<*nummappings) {

    int have_align = i;
    int j = i+1;
    int hashmatches = winsize;     // nucleotides under matching hashes
    char keep = 0;                 // make sure to keep manually entered mappings (from externally mapped reads)

    while (j < *nummappings  &&  mappings[j]->locus - mappings[i]->locus < maxdistance) {
      hashmatches += mappings[j]->end - gu_max( mappings[j-1]->end, mappings[j]->end - winsize );
      if (mappings[j]->readalignment) have_align=j;  // prefer to keep a candidate with alignment
      keep |= (mappings[j]->end < 0);                // end<0 signifies manually entered mapping
      j++;
    }
    maxhashmatches = gu_max( maxhashmatches, hashmatches );
    if (keep)
      hashmatches = 9999;  // does not affect maxhashmatches; but ensures mapping is kept
    mappings[ have_align ]->hashmatches = hashmatches;
    *writer = mappings[ have_align ];
    //DEBUG1("  adding locus=%lli\n", (*writer)->locus);
    writer++;
    written++;
    i = j;

  }

  *nummappings = written;

  //DEBUG1("** Uniqueify pass 2; now %u mappings\n",*nummappings);

  // now remove candidates with few hash matches
  
  if (minhashmatches > -1) {

    writer = mappings;
    written = 0;
    i = 0;
    while (i<*nummappings) {
      if (mappings[i]->hashmatches > maxhashmatches - minhashmatches) {
	*writer = mappings[i];
	//DEBUG1("  adding locus=%lli\n", (*writer)->locus);
	writer++;
	written++;
      }
      i++;
    }
    *nummappings = written;

  }

}




// this also implements pseudo-random placement of repetitive reads
void unique_alignments_ptr( struct mapping** mappings, int* nummappings ) {

  if (*nummappings == 0) return;

  struct mapping** writer = mappings + 1;
  int written = 1;
  int i=1;
  static int pseudo_random = 0;
  struct mapping* temp;
  int bestscore = (*nummappings ? mappings[0]->score : 0);
  int placed_random = 0;

  pseudo_random += 1;

  for (; i<*nummappings; i++) {

    if ((!placed_random) && (mappings[i]->score > bestscore)) {
      // first non-optimal hit found; swap random element among 0..written (excludes i) with first one
      int idx = pseudo_random % written;
      temp = mappings[idx];
      mappings[idx] = mappings[0];
      mappings[0] = temp;
      placed_random = 1;
    }

    // ideally, if mappings have different scores, then their locus is different too.
    // however because of the way alignments are pegged to their potentially imperfect
    // hash anchors, small discrepancies may occur.
    if (mappings[i]->locus != mappings[i-1]->locus) {
      *writer = mappings[i];
      writer++;
      written++;

    }
  }

  if (!placed_random) {
      // swap random element with first one
      int idx = pseudo_random % written;
      temp = mappings[idx];
      mappings[idx] = mappings[0];
      mappings[0] = temp;
  }

  *nummappings = written;

}


void  align_score_candidates_sse2_ptr( struct mapping** mappings, int *nummappings, int readlength, char* sequence, int* intquality, 
				       int overhang, PyObject* genomemmap, int datastart, int64_t genomesize, 
				       int gapopen, int gapextend, int nucprior, char* alignmentbuffer, char adjustLocus, int start_record ) {

  // alignmentbuffer should be at least 8*readlength + 64 bytes long

  struct mapping** writer = mappings + start_record;
  int written = start_record;
  char write;
  char genomesequence[ readlength+32 ];    // must also fit overhang
  int sumquality = 0;
  int lowestscore = 0x3c00;
  int abcutoff,score;
  int i;
  for (i=0; i<readlength; i++) sumquality += intquality[i];
  // we have 14 bits of score to play with; algorithm uses 0x200 (0x800/4) headroom; allow 16*64=0x400 for flanks;
  // altogether 0x3a00.  Scale down quality scores as required, and scale them up later
  int divisor = 1 + (sumquality + 256) / 0x3a00;
  if (divisor > 1) {
    for (i=0; i<readlength; i++) intquality[i] /= divisor;
    gapopen = gapopen==0 ? 0 : 1 + (gapopen-1)/divisor;
    gapextend = gapextend==0 ? 0 : 1 + (gapextend-1)/divisor;
    nucprior = nucprior==0 ? 0 : 1 + (nucprior-1)/divisor;
  }

  for (i=written; i<*nummappings; i++) {

    write = 1;
    if (mappings[i]->genomealignment) {
      free(mappings[i]->genomealignment);
      free(mappings[i]->readalignment);
      mappings[i]->genomealignment = mappings[i]->readalignment = NULL;
    }
    if (overhang < 0) {
      printf("Error: align_candidate_c: negative overhang.  Program will probably crash now.\n");
    } else if (overhang == 0) {
      // if overhang==0, use the read + the genome, and compute the score in the old-fashioned way
      //mappings[i]->readalignment = NULL;
      mappings[i]->genomealignment = getgenome( genomemmap, datastart, 
						mappings[i]->locus, mappings[i]->locus + readlength, genomesize );
      if (mappings[i]->genomealignment != NULL) {
	mappings[i]->readalignment = malloc( readlength );
	strncpy( mappings[i]->readalignment, sequence, readlength+1 );
	mappings[i]->score = divisor*score_alignment_c( mappings[i]->readalignment, mappings[i]->genomealignment, 
							intquality, gapopen, gapextend, nucprior );
	if (!alignmentbuffer) {
	  free( mappings[i]->readalignment );
	  free( mappings[i]->genomealignment );
	  //mappings[i]->readalignment = NULL;
	  //mappings[i]->genomealignment = NULL;
	}
      }
    } else {
      // main case -- 15bp band
      if (overhang == 1) {
	getgenome_buffer( genomesequence, genomemmap, datastart, 
			  mappings[i]->locus - 8, mappings[i]->locus + readlength + 7, genomesize ); 
	if (genomesequence[0] == 0) {
	  //mappings[i]->readalignment = NULL;
	  //mappings[i]->genomealignment = NULL;
	  mappings[i]->score = 9999;
	} else {
	  int newstart;
	  if (alignmentbuffer) {
	    mappings[i]->genomealignment = (char*)malloc( (2*readlength+16)*sizeof(char) );
	    mappings[i]->readalignment = (char*)malloc( (2*readlength+16)*sizeof(char) );
	  
	    mappings[i]->score = divisor*traceback_sse2( genomesequence, sequence, intquality, 
							 mappings[i]->genomealignment, mappings[i]->readalignment,
							 readlength+15, readlength, gapopen, gapextend, nucprior, &newstart,
							 (lowestscore + ALPHA_BETA_CUTOFF)/divisor ); 
	    if (mappings[i]->score < lowestscore) lowestscore = mappings[i]->score;

	    if (adjustLocus) mappings[i]->locus += newstart - 8;
	  } else {
	    //mappings[i]->genomealignment = NULL;
	    //mappings[i]->readalignment = NULL;
	    abcutoff = (lowestscore + ALPHA_BETA_CUTOFF)/divisor;
	    score = traceback_sse2( genomesequence, sequence, intquality, 
				    NULL, NULL, 
				    readlength+15, readlength, gapopen, gapextend, nucprior, &newstart, abcutoff );
	    write = (score < abcutoff);
	    score *= divisor;
	    //mappings[i]->score = score;
	    if (score < lowestscore) lowestscore = score;
	    if (adjustLocus) mappings[i]->locus += newstart - 8;
	  }
	}
      } else {
	// main case -- 30 bp band
	getgenome_buffer( genomesequence, genomemmap, datastart, 
			  mappings[i]->locus - 15, mappings[i]->locus + readlength + 15, genomesize ); 
	if (genomesequence[0] == 0) {
	  //mappings[i]->readalignment = NULL;
	  //mappings[i]->genomealignment = NULL;
	  mappings[i]->score = 9999;
	} else {
	  int newstart1, newstart2;
	  int score1, score2;
	  if (alignmentbuffer) {
	    char* ga1 = alignmentbuffer;
	    char* ga2 = ga1 + (2*readlength+16);
	    char* ra1 = ga1 + (4*readlength+32);
	    char* ra2 = ga1 + (6*readlength+48);
	    char* gresult = (char*)malloc( (2*readlength+16)*sizeof(char) );
	    char* rresult = (char*)malloc( (2*readlength+16)*sizeof(char) );

	    score1 = divisor*traceback_sse2( genomesequence, sequence, intquality, 
					     ga1, ra1,
					     readlength+15, readlength, gapopen, gapextend, nucprior, &newstart1,
					     (lowestscore + ALPHA_BETA_CUTOFF)/divisor ); 

	    score2 = divisor*traceback_sse2( genomesequence+15, sequence, intquality, 
					     ga2, ra2,
					     readlength+15, readlength, gapopen, gapextend, nucprior, &newstart2,
					     (lowestscore + ALPHA_BETA_CUTOFF)/divisor ); 
	    if (score1 < score2) {
	      mappings[i]->score = score1;
	      mappings[i]->genomealignment = strcpy(gresult,ga1);
	      mappings[i]->readalignment = strcpy(rresult,ra1);
	      if (adjustLocus) mappings[i]->locus += newstart1 - 15;
	      //free( ga2 );
	      //free( ra2 );
	    } else {
	      mappings[i]->score = score2;
	      mappings[i]->genomealignment = strcpy(gresult,ga2);
	      mappings[i]->readalignment = strcpy(rresult,ra2);
	      if (adjustLocus) mappings[i]->locus += newstart2;
	      //free( ga1 );
	      //free( ra1 );
	    }
	    if (mappings[i]->score < lowestscore) lowestscore = mappings[i]->score;

	  } else {
	    abcutoff = (lowestscore + ALPHA_BETA_CUTOFF)/divisor;
	    score1 = traceback_sse2( genomesequence, sequence, intquality, 
				     NULL, NULL,
				     readlength+15, readlength, gapopen, gapextend, nucprior, &newstart1, abcutoff);
	    score2 = traceback_sse2( genomesequence+15, sequence, intquality, 
				     NULL, NULL,
				     readlength+15, readlength, gapopen, gapextend, nucprior, &newstart2, abcutoff);
	    if (score2 < score1) score1 = score2;
	    write = (score1 < abcutoff);
	    score1 *= divisor;
	    mappings[i]->score = score1;
	    // do not change loci if we're not making alignments
	    if (score1 < lowestscore) lowestscore = score1;
	    mappings[i]->genomealignment = NULL;
	    mappings[i]->readalignment = NULL;
	  }
	}
      }
    }
 
    if ((write) && (mappings[i]->readalignment != NULL || (!alignmentbuffer))) {
      *writer = mappings[i];
      writer++;
      written++;
    }
  }
  *nummappings = written;
}



void align_candidates_ptr( struct mapping** mappings, int *nummappings, char* sequence, int overhang, 
			   PyObject* genomemmap, int datastart, int64_t genomesize, 
			   int gapstartpenalty, PyObject* pyscoretable, int band ) {

  struct mapping** writer = mappings;
  int written = 0;
  int scoretable[ 25 ];
  int i;

  for (i=0; i< 25; i++) {
    scoretable[i] = SCOREFACTOR * PyInt_AsLong( PyList_GetItem( pyscoretable, i ) );
  }

  for (i=0; i<*nummappings; i++) {

    if (mappings[i]->readalignment == NULL) {

      align_candidate_c( mappings[i]->locus, mappings[i]->start, mappings[i]->end,
			 sequence, overhang, genomemmap, datastart, genomesize,
			 gapstartpenalty, scoretable, 
			 &(mappings[i]->readalignment), &(mappings[i]->genomealignment), &(mappings[i]->locus), band );
    }

    if (mappings[i]->readalignment != NULL) {
      *writer = mappings[i];
      writer++;
      written++;
    }
  }

  *nummappings = written;

}


void score_alignments_ptr( struct mapping** mappings, int nummappings, int readlength, PyObject* quality, 
			   int gapopen, int gapextend, int nucprior ) {

  int intquality[ readlength ];
  int i;
  PyObject* qualseq = PySequence_Fast( quality, "Expected sequence of quality values" );
  for (i=0; i<readlength; i++) {
	intquality[i] = PyInt_AsLong( PySequence_Fast_GET_ITEM( qualseq, i ) );
  }
  Py_DECREF( qualseq );
  for (i=0; i<nummappings; i++) {
    mappings[i]->score = score_alignment_c( mappings[i]->readalignment, mappings[i]->genomealignment, intquality, 
					    gapopen, gapextend, nucprior );
  }
}



// See if an indel exists; if so build the implied reference sequence
char get_implied_read( char* readalignment, char* genomealignment, int* intquality, char *impliedseq, int *impliedq) {

  int seqlen = 0;
  int indel = 0;
  int alnidx = 0;
  int readidx = 0;
  while (readalignment[alnidx]) {
    if (readalignment[alnidx] != '-' && genomealignment[alnidx] != '-') {
      impliedseq[seqlen] = readalignment[alnidx];
      impliedq[seqlen] = intquality[readidx];
      seqlen++;
      readidx++;
    } else {
      indel = 1;
      if (readalignment[alnidx] == '-') {
	impliedseq[seqlen] = genomealignment[alnidx];
	impliedq[seqlen] = 30;
	seqlen++;
      } else {
	readidx++;
      }
    }
    alnidx++;
  }
  impliedseq[seqlen] = 0;
  //printf ("seqlen=%u alnidx=%u readidx=%u lastq=%u lastreadq=%u\n",seqlen,alnidx,readidx,impliedq[seqlen-1],intquality[readidx-1]);
  return (indel && seqlen > winsize+1);
}

    

PyObject* locate_align_score_candidates_c( char* sequence, 
					   PyObject* hashtablemmap, PyObject* genomemmap, int datastart, 
					   long mask, char tryvariants, int maxfingerprintvariants, int maxdistance, int minhashmatches,
					   int64_t genomesize, int gapstartpenalty, PyObject* scoretable,
					   PyObject* quality, int gapopen_score, int gapextend_score, int nucprior_score,
					   struct mapping** savedmappings, int* numsavedmappings, struct mapping*** savedmappingptrs, 
					   int* numsavedmappingptrs, int lowqthreshold, int band, int max_candidates,
					   int64_t locus_to_consider ) {

  int nummappings = 0;
  int mappingsize = 1000;
  int directhighcounts = 0;
  int indirecthighcounts = 0;
  int i;

  int seqbitstringlen;
  int readlength = strlen(sequence);
  char* seqbitstring = sequencetobitstring( sequence, readlength, &seqbitstringlen );

  struct mapping* mappings = malloc( sizeof(struct mapping)*mappingsize );

  // add locus from pre-mapper for consideration, if any
  // we should really enter the position within the read where the match occurred, rather than -1;
  // this will work in the majority of cases though.
  if (locus_to_consider != 0) {
    addmapping( &mappings, &nummappings, &mappingsize, locus_to_consider, -1, -1 );
  }

  // get coordinates of reliable portion of read
  int qlen = PyObject_Length( quality );
  int intquality[qlen];
  int repeatmask[qlen];
  PyObject* qualseq = PySequence_Fast( quality, "Expected sequence of quality values" );
  for (i=0; i<qlen; i++) {
    intquality[i] = PyInt_AsLong( PySequence_Fast_GET_ITEM( qualseq, i ) );
    repeatmask[i] = 999;
  }
  Py_DECREF( qualseq );
  int leftreliable = 0;
  while (leftreliable < qlen && intquality[leftreliable] <= lowqthreshold) leftreliable++;
  int rightreliable = qlen-1;
  while (rightreliable > leftreliable && intquality[rightreliable] <= lowqthreshold) rightreliable--;
  rightreliable++;

  locate_candidates_fingerprint_hash3( &mappings, &nummappings, &mappingsize, &directhighcounts, &indirecthighcounts,
				       sequence, seqbitstring, seqbitstringlen, hashtablemmap, genomemmap, datastart, mask, 
				       tryvariants, maxfingerprintvariants, 21, leftreliable, rightreliable,
				       intquality, repeatmask) ;
  free (seqbitstring);
  int numpointers = nummappings;
  //printf("Candidates: %u\n", numpointers);
  struct mapping** mappingptrs = malloc( sizeof( void* )*numpointers);
  for (i=0; i<numpointers; i++) mappingptrs[i] = &mappings[i];

  int64_t templocus;
#define compare_ls(a,b) ((templocus=(*(a))->locus - (*(b))->locus ) < 0 || (templocus ==0 && (*(a))->start < (*(b))->start) )
  QSORT( struct mapping*, mappingptrs, numpointers, compare_ls )

  // uniquify, and remove candidates with few hash matches
  unique_candidates_ptr( mappingptrs, &numpointers, maxdistance, minhashmatches );
  //printf("Unique; Candidates: %u\n", numpointers);

  // sort by number of hash matches, to optimize AB shortcut
#define compare_hm(a,b) ( (*(a))->hashmatches > (*(b))->hashmatches )
  QSORT( struct mapping*, mappingptrs, numpointers, compare_hm );

#ifdef __SSE2__
  align_score_candidates_sse2_ptr( mappingptrs, &numpointers, strlen(sequence), sequence, intquality,
  				   maxdistance, genomemmap, datastart, genomesize,
  				   gapopen_score, gapextend_score, nucprior_score, NULL, 0, 0 );
#else
  align_candidates_ptr( mappingptrs, &numpointers, sequence, maxdistance, genomemmap, datastart, genomesize,
			gapstartpenalty, scoretable, band );
  score_alignments_ptr( mappingptrs, numpointers, strlen(sequence), quality, gapopen_score, gapextend_score, nucprior_score);
#endif

  //printf("Aligned; Candidates: %u\n", numpointers);
  // sort by (score then) locus
  int tempscore;
#define compare_sl(a,b) ((tempscore=(*(a))->score - (*(b))->score) < 0 || (tempscore ==0 && (*(a))->locus < (*(b))->locus) )
  QSORT( struct mapping*, mappingptrs, numpointers, compare_sl );

  // remove duplicates again
  unique_alignments_ptr( mappingptrs, &numpointers );
  //printf("Unique again; Candidates: %u\n", numpointers);

  // Addition: if top candidate contains indel, build new candidates from part of read that aligns to the genome
  // This finds additional repetitive candidates if these exist; the entropy filter won't catch these cases
  struct mapping* impliedmappings = NULL;
  int impliednummappings = 0;

#ifdef __SSE2__
  if (numpointers>0) {
    int testmaps = gu_min(UNDUP_CANDIDATES, numpointers);

    char* align_buffer = (char*)malloc( (8*readlength+64)*sizeof(char) );
    align_score_candidates_sse2_ptr( mappingptrs, &testmaps, readlength, sequence, intquality,
				     maxdistance, genomemmap, datastart, genomesize,
				     gapopen_score, gapextend_score, nucprior_score, align_buffer, 0, 0 );
    free(align_buffer);
    //unique_alignments_ptr( mappingptrs, &testmaps );
    //    printf("Unique again!; Candidates: %u\n", numpointers);

    if (testmaps>0) {
      int alen = strlen( mappingptrs[0]->readalignment );
      char impliedseq[alen + 1];
      int impliedq[alen + 1];
      if (get_implied_read( mappingptrs[0]->readalignment, mappingptrs[0]->genomealignment, intquality, impliedseq, impliedq )) {
	//printf("Found alignment (%u candidates)\n", numpointers);
	//for (i=0; i<gu_min(5,numpointers); i++) printf(" cand %u loc=%lu hm=%u\n",i,mappingptrs[i]->locus, mappingptrs[i]->hashmatches);
	//printf("Alen=%u ",alen);
	alen = strlen(impliedseq);
	//printf("implseqlen=%u lastq=%u\n",alen,impliedq[alen-1]);
	int impliedrepeatmask[alen + 1];
	for (i=0; i<alen; i++) impliedrepeatmask[i]=999;
	int impliedseqbitstringlen;
	char* impliedseqbitstring = sequencetobitstring( impliedseq, alen, &impliedseqbitstringlen );
	int impliedleftreliable = 0, impliedrightreliable = alen-1;
	//printf ("Starting left=%u right=%u\n",impliedleftreliable,impliedrightreliable);
	while (impliedleftreliable < alen && impliedq[impliedleftreliable] <= lowqthreshold) impliedleftreliable++;
	while (impliedrightreliable > impliedleftreliable && impliedq[impliedrightreliable] <= lowqthreshold) impliedrightreliable--;
	//printf ("End left=%u right=%u\n",impliedleftreliable,impliedrightreliable);
	impliedrightreliable++;
	// make new candidates (in new buffer, since realloc may move old buffer, and we keep pointers)
	int impliedmappingsize = mappingsize;
	impliedmappings = malloc( sizeof(struct mapping)*impliedmappingsize );
	locate_candidates_fingerprint_hash3( &impliedmappings, &impliednummappings, &impliedmappingsize, &directhighcounts, &indirecthighcounts,
					     impliedseq, impliedseqbitstring, impliedseqbitstringlen, hashtablemmap, genomemmap, datastart, mask, 
					     tryvariants, maxfingerprintvariants, 21, impliedleftreliable, impliedrightreliable,
					     impliedq, impliedrepeatmask) ;
	free (impliedseqbitstring);

	//printf("Found %u new candidates; total %u\n",impliednummappings,impliednummappings+numpointers);
	mappingptrs = realloc( mappingptrs, sizeof( void* )*(numpointers + impliednummappings) );
	for (i=0; i<impliednummappings; i++) mappingptrs[i + numpointers] = &impliedmappings[i];
	numpointers += impliednummappings;
	
	QSORT( struct mapping*, mappingptrs, numpointers, compare_ls ) ;

	unique_candidates_ptr( mappingptrs, &numpointers, maxdistance, minhashmatches );
	//printf("Uniquified by locus: total %u\n",numpointers);
	//for (i=0; i<gu_min(5,numpointers); i++) printf(" cand %u loc=%lu hm=%u\n",i,mappingptrs[i]->locus, mappingptrs[i]->hashmatches);

	QSORT( struct mapping*, mappingptrs, numpointers, compare_hm );

	align_score_candidates_sse2_ptr( mappingptrs, &numpointers, strlen(sequence), sequence, intquality,
					 maxdistance, genomemmap, datastart, genomesize,
					 gapopen_score, gapextend_score, nucprior_score, NULL, 0, 0 );
	//printf("Aligned: total %u\n",numpointers);

	QSORT( struct mapping*, mappingptrs, numpointers, compare_sl );

	unique_alignments_ptr( mappingptrs, &numpointers );

	//printf("Uniqueified again: total %u\n",numpointers);

      }
    }
  }
#endif

  // First select a limited number of candidates which are not too unlikely
  int likelihood = (numpointers>0) ? mappingptrs[0]->score : 0;
  while (numpointers > 0 && mappingptrs[numpointers-1]->score > likelihood + MAX_LIKELIHOOD_DIFF) numpointers--;

  // collect records up to a maximum number; realign; sort; uniqueify; and repeat until buffer full
  char* align_buffer = (char*)malloc( (8*readlength+64)*sizeof(char) );
  int unprocessed_records = 0;
  int total_records = numpointers;
  numpointers = 0;
  while (unprocessed_records < total_records && numpointers < max_candidates) {  // records available, and room for more
    int start_of_buffer = numpointers;
    while (unprocessed_records < total_records && numpointers < max_candidates) {
      mappingptrs[numpointers] = mappingptrs[unprocessed_records];
      ++numpointers;
      ++unprocessed_records;
    }
    align_score_candidates_sse2_ptr( mappingptrs, &numpointers, readlength, sequence, intquality,
				     maxdistance, genomemmap, datastart, genomesize,
				     gapopen_score, gapextend_score, nucprior_score, align_buffer, 1, start_of_buffer );
    QSORT( struct mapping*, mappingptrs, numpointers, compare_sl );
    unique_alignments_ptr( mappingptrs, &numpointers );
  }
  free (align_buffer);

  PyObject* pylist = PyList_New( numpointers );

  if (numsavedmappingptrs != NULL) *numsavedmappingptrs = numpointers;

  double posterior = 0.0;

  while (numpointers > 0) {

    --numpointers;
    //printf ("Pointer %u score %u locus %lld likelihood excess %u\n",numpointers, mappingptrs[numpointers]->score, (long long int)mappingptrs[numpointers]->locus, mappingptrs[numpointers]->score-likelihood);
    posterior += exp( (mappingptrs[numpointers]->score - likelihood) * MATH_DB );
    //printf ("Posterior now %f\n",posterior);
    // cap scores at 49000, so that paired-end scores are guaranteed to be below sentinel 99999
    PyObject* pyscore = PyInt_FromLong( gu_min( 49000, mappingptrs[numpointers]->score) );
    PyObject* pyaln1 = PyString_FromString( mappingptrs[numpointers]->readalignment );
    PyObject* pyaln2 = PyString_FromString( mappingptrs[numpointers]->genomealignment );
    PyObject* pylocus = PyLong_FromLongLong( mappingptrs[numpointers]->locus );
    PyObject* pyoffset = PyInt_FromLong( mappingptrs[numpointers]->start );
    PyObject* pyoffset2 = PyInt_FromLong( mappingptrs[numpointers]->end );
    PyObject* pyhashmatches = PyInt_FromLong( mappingptrs[numpointers]->hashmatches );
    PyObject* tuple = PyTuple_New(7);
    PyTuple_SET_ITEM(tuple, 0, pyscore);
    PyTuple_SET_ITEM(tuple, 1, pyaln1);
    PyTuple_SET_ITEM(tuple, 2, pyaln2); 
    PyTuple_SET_ITEM(tuple, 3, pylocus);
    PyTuple_SET_ITEM(tuple, 4, pyoffset);
    PyTuple_SET_ITEM(tuple, 5, pyoffset2);
    PyTuple_SET_ITEM(tuple, 6, pyhashmatches);
    PyList_SET_ITEM( pylist, numpointers, tuple );   // stolen

  }

  //printf("## TIMINGS: %lu %lu %lu %lu (mappings=%u)\n",locate_time,qsort1_time,align_time,qsort2_time,nummappings);

  if (savedmappings == NULL) {

    // free mappings array and malloc'ed memory for alignment strings

    for (i=0; i<nummappings; i++) {
      free( mappings[i].readalignment );
      free( mappings[i].genomealignment );
    }
    for (i=0; i<impliednummappings; i++) {
      free( impliedmappings[i].readalignment );
      free( impliedmappings[i].genomealignment );
    }

    free( mappingptrs );
    free( mappings );
    free( impliedmappings );

  } else {

    // return the raw results to the calling routine
    assert(0);

    *savedmappings = mappings;
    *numsavedmappings = nummappings;
    *savedmappingptrs = mappingptrs;

  }

  // build repeat mask
  PyObject* pyrepeatlist = PyList_New( qlen );
  for (i=0; i<qlen; i++) PyList_SET_ITEM( pyrepeatlist, i, PyInt_FromLong( repeatmask[i] ) );

  // build tuple of high hash counters, and repeat mask
  PyObject* tuple = PyTuple_New(3);
  PyObject* pydirect = PyInt_FromLong( directhighcounts );
  PyObject* pyindirect = PyInt_FromLong( indirecthighcounts );
  PyTuple_SET_ITEM(tuple, 0, pydirect);
  PyTuple_SET_ITEM(tuple, 1, pyindirect);
  PyTuple_SET_ITEM(tuple, 2, pyrepeatlist);

  // build tuple of return values
  PyObject* retval = PyTuple_New(3);
  PyObject* pyposterior = PyFloat_FromDouble( posterior );
  PyTuple_SET_ITEM(retval, 0, tuple);
  PyTuple_SET_ITEM(retval, 1, pylist);
  PyTuple_SET_ITEM(retval, 2, pyposterior);

  return retval;

}




inline PyObject* locate_align_score_candidates( char* sequence, /* char* seqbitstring, int seqbitstringlen,  */
						PyObject* hashtablemmap, PyObject* genomemmap, int datastart, 
						long mask, char tryvariants, int maxfingerprintvariants, int maxdistance, int minhashmatches,
						int64_t genomesize, int gapstartpenalty, PyObject* scoretable,
						PyObject* quality, int gapopen_score, int gapextend_score, int nucprior_score, 
						int lowqthreshold, int band, int max_candidates, int64_t locus_to_consider ) {

  PyObject* pyobj = locate_align_score_candidates_c(  sequence, 
						      hashtablemmap, genomemmap, datastart, 
						      mask, tryvariants, maxfingerprintvariants, maxdistance, minhashmatches,
						      genomesize, gapstartpenalty, scoretable,
						      quality, gapopen_score, gapextend_score, nucprior_score,
						      NULL, NULL, NULL, NULL, lowqthreshold, band, max_candidates, locus_to_consider );
  return pyobj;
}



double sympol2(double sp,double spp) { return (sp*sp-spp)/2.0; }
double sympol3(double sp,double spp,double sppp) { return (sp*sp*sp - 3*sp*spp + 2*sppp)/6.0; }
double sympol4(double sp,double spp,double sppp,double spppp) { return (sp*sp*sp*sp + 8*sppp*sp + 3*spp*spp - 6*spp*sp*sp - 6*spppp)/24.0; }
double sympol5(double sp,double spp,double sppp,double spppp,double sppppp) { return (sp*sp*sp*sp*sp - 10*spp*sp*sp*sp + 15*spp*spp*sp + 20*sppp*sp*sp - 20*sppp*spp - 30*spppp*sp + 24*sppppp)/120.0; }


PyObject* calc_mutprobs( PyObject* quals, PyObject* repmask, int tryvariants ) {

  int qlen = PyObject_Length( quals );
  int qvals[ qlen ];
  int rmask[ qlen ];
  int i,j,k;
  double sp,spp,sppp,spppp,sppppp;
  //double p1=1.0,p2=1.0,p3=1.0,p4=1.0;
  double e2,e3,e4,e5;
  double q,minp;
  double expectedQ;
  double oneminusp;
  int beststart,bestend;

  for (i=0; i<qlen; i++) {
    qvals[i] = PyInt_AsLong( PyList_GetItem( quals, i ) );
    rmask[i] = PyInt_AsLong( PyList_GetItem( repmask, i ) );
  }

  int maxidx = qlen - 15 + 1;
  if (tryvariants == -1) tryvariants = (maxidx <= 20 ? 1 : (maxidx <= 35 ? 2 : 3));

  //double lowsum = phredscores[gu_min(5*lowqthreshold,79)];

  double norepeatpa1[qlen], truepa1[qlen], pa1[qlen], pa2[qlen], pa3[qlen], pa4[qlen], pa5[qlen];
  expectedQ = 0.0;
  for (i=0; i<qlen; i++) {
    norepeatpa1[i] = truepa1[i] = qvals[i]<=0 ? 0.80 : (qvals[i]<80 ? phredscores[qvals[i]] : 1.0e-8);
    expectedQ += gu_max(0,qvals[i]) * truepa1[i];
    //truepa1[i] += (1.0-truepa1[i]) * (rmask[ gu_max(i-14,0) ]<80 ? phredscores[rmask[gu_max(i-14,0)]] : 1.0e-8 );
    truepa1[i] += (1.0-truepa1[i]) * (rmask[ i ]<80 ? phredscores[rmask[ i ]] : 1.0e-8 );
    truepa1[i] = gd_min(0.9, truepa1[i]); // repeat masking can make this 1.0, breaking p/1-p
    // use p/1-p as the variable for the Newton polynomials
    pa1[i] = truepa1[i] / (1.0-truepa1[i]);
    //printf ("i=%u p=%f  p/1-p=%f\n",i,truepa1[i],pa1[i]);
    pa2[i] = pa1[i]*pa1[i];
    pa3[i] = pa2[i]*pa1[i];
    pa4[i] = pa3[i]*pa1[i];
    pa5[i] = pa4[i]*pa1[i];
  }

  if (tryvariants == 3) {
    const int stride=5;
    // non-overlapping windows; use different algorithm to compute probability of not finding candidate
    minp = 0.0;
    for (i=0; i<stride; i++) {
      e2 = 1.0;
      for (j=i; j<qlen-winsize+1; j+=winsize) {
	// use first 2 probabilities including repeat mask
	q = (1.0 - truepa1[j]) * (1.0 - truepa1[j+1]);
	sp = pa1[j] + pa1[j+1];
	// DEBUG - for printing below
	// spp = q;
	// for the rest of window, use probabilities without repeat mask,
	// as the repeat probability pertains to the entire 15bp window
	for (k=j+2; k<j+winsize; k++) {
	  q *= (1.0 - norepeatpa1[k]);
	  sp += norepeatpa1[k] / (1.0-norepeatpa1[k]);
	}
	q *= (1.0 + sp); // probability of 0 or 1 mutation
	e2 *= (1.0 - q); // probability of >= 2 mutations in all windows
	//printf (" Calcmutprobs Phase=%u Window=%u (repeat=%f) Prob01=%f Prob >=2 muts in all %f\n",i,j,spp,q,e2);
      }
      minp += e2 / stride;
      //printf(" minp now %f\n",minp);
    }
  } else {
    // algorithm for overlapping windows
    beststart = bestend = 0;
    minp = 1.0;
    for (i=15; i<=qlen; i++) {
      sp=spp=sppp=spppp=sppppp=0;
      oneminusp = 1.0;
      //printf("i=%u\n",i);
      for (j=i-1; j>=0 && j>=i-52; j--) {
	oneminusp *= 1.0 - truepa1[j];
	//printf("i=%u j=%u pa1[j]=%f oneminusp=%f\t\t(sympol5=%f)\n",i,j,pa1[j],oneminusp, sympol5(sp,spp,sppp,spppp,sppppp));
	sp += pa1[j];
	spp += pa2[j];
	sppp += pa3[j];
	spppp += pa4[j];
	sppppp += pa5[j];
	e2 = sympol2(sp,spp); e3 = sympol3(sp,spp,sppp); e4 = sympol4(sp,spp,sppp,spppp);
	e5 = 1.0 - oneminusp*(1.0 + sp + e2 + e3 + e4);  // probability of 5 or more errors
	if (e5 > 0.01) {
	  //printf ("i=%u j=%u 5ormore = %f\n",i,j,e5);
	  break;
	}
	if (j > i-15) continue;
	// i-52 < j <= i-15
	k = 5*(i-15-j);
	q = oneminusp*mapstatscacheD[k];
	q += gd_max(0.0,1-q)*oneminusp*mapstatscacheD[k+1]*sp;
	q += gd_max(0.0,1-q)*oneminusp*mapstatscacheD[k+2]*e2;
	q += gd_max(0.0,1-q)*oneminusp*mapstatscacheD[k+3]*e3;
	q += gd_max(0.0,1-q)*oneminusp*mapstatscacheD[k+4]*e4;
	q += gd_max(0.0,1-q)*e5;  // safety - 5 mutations always bad
	//printf("i=%u j=%u q=%f oneminusp=%f\n",i,j,q,oneminusp);
	if (q < minp || (q < minp + 1e-8 && i-j > bestend-beststart)) {
	  minp = q;
	  beststart = j;
	  bestend = i;
	  /*
	    p1 = gd_min( gd_max( oneminusp*sp, 1.0e-10 ), 1.0 );
	    p2 = gd_min( gd_max( oneminusp*sympol2(sp,spp), 1.0e-10 ), 1.0 );
	    p3 = gd_min( gd_max( oneminusp*sympol3(sp,spp,sppp), 1.0e-10 ), 1.0 );
	    p4 = gd_min( gd_max( oneminusp*sympol4(sp,spp,sppp,spppp), 1.0e-10 ), 1.0 );
	  */
	  //printf("Best minp=%f start/end=%u/%u\n",minp,j,i);
	  //printf("p1-4: %f %f %f %f\n",oneminusp*sp,oneminusp*sympol2(sp,spp),oneminusp*sympol3(sp,spp,sppp),oneminusp*sympol4(sp,spp,sppp,spppp));
	}
      }
    }
    //printf("Longest stretch: %u  Phred %u    Minp %f\n",bestend-beststart,(int)(10*log(minp)/log(0.1)),minp);
    //int efflen = bestend-beststart;
  }

  // underflow occurs for very long reads
  if (minp < 1e-10) minp = 1e-10;

  PyObject* retval = PyTuple_New(2);  
  PyTuple_SET_ITEM(retval, 0, PyInt_FromLong( (int)(10*log(minp)/log(0.1)) ) );
  PyTuple_SET_ITEM(retval, 1, PyInt_FromLong( (int)(expectedQ + 0.5) ) );
  return retval;

}


PyObject* calc_mutprobs2( PyObject* quals, PyObject* repmask, int lowqthreshold ) {

  int qlen = PyObject_Length( quals );
  int qvals[ qlen ];
  int rmask[ qlen ];
  int i,j,k;
  double sp,spp,sppp,spppp,sppppp;
  //double p1=1.0,p2=1.0,p3=1.0,p4=1.0;
  double e2,e3,e4,e5;
  double q,minp;
  double oneminusp;
  int beststart,bestend;

  for (i=0; i<qlen; i++) {
    qvals[i] = PyInt_AsLong( PyList_GetItem( quals, i ) );
    rmask[i] = PyInt_AsLong( PyList_GetItem( repmask, i ) );
  }

  //double lowsum = phredscores[gu_min(5*lowqthreshold,79)];

  double truepa1[qlen], pa1[qlen], pa2[qlen], pa3[qlen], pa4[qlen], pa5[qlen];
  for (i=0; i<qlen; i++) {
    truepa1[i] = qvals[i]<=0 ? 0.80 : (qvals[i]<80 ? phredscores[qvals[i]] : 1.0e-8);
    truepa1[i] += (1.0-truepa1[i]) * (rmask[ gu_max(i-14,0) ]<80 ? phredscores[rmask[gu_max(i-14,0)]] : 1.0e-8 );
    // use p/1-p as the variable for the Newton polynomials
    pa1[i] = truepa1[i] / (1.0-truepa1[i]);
    pa2[i] = pa1[i]*pa1[i];
    pa3[i] = pa2[i]*pa1[i];
    pa4[i] = pa3[i]*pa1[i];
    pa5[i] = pa4[i]*pa1[i];
  }

  beststart = bestend = 0;
  minp = 1.0;
  for (i=15; i<=qlen; i++) {
    sp=spp=sppp=spppp=sppppp=0;
    oneminusp = 1.0;
    //printf("i=%u\n",i);
    for (j=i-1; j>=0 && j>=i-52; j--) {
      oneminusp *= 1.0 - truepa1[j];
      //printf("i=%u j=%u pa1[j]=%f oneminusp=%f\n",i,j,pa1[j],oneminusp);
      sp += pa1[j];
      spp += pa2[j];
      sppp += pa3[j];
      spppp += pa4[j];
      sppppp += pa5[j];
      e2 = sympol2(sp,spp); e3 = sympol3(sp,spp,sppp); e4 = sympol4(sp,spp,sppp,spppp);
      e5 = 1.0 - oneminusp*(1.0 + sp + e2 + e3 + e4);
      if (e5 > 0.01) {
	//printf ("i=%u j=%u sympol5=%f\n",i,j,sympol5(sp,spp,sppp,spppp,sppppp));
	break;
      }
      if (j > i-15) continue;
      // i-52 <= j <= i-15
      k = 5*(i-15-j);
      q = oneminusp*mapstatscacheD[k];
      q += gd_max(0.0,1-q)*oneminusp*mapstatscacheD[k+1]*sp;
      q += gd_max(0.0,1-q)*oneminusp*mapstatscacheD[k+2]*e2;
      q += gd_max(0.0,1-q)*oneminusp*mapstatscacheD[k+3]*e3;
      q += gd_max(0.0,1-q)*oneminusp*mapstatscacheD[k+4]*e4;
      q += gd_max(0.0,1-q)*e5;  // safety - 5 mutations always bad
      //printf("i=%u j=%u q=%f oneminusp=%f\n",i,j,q,oneminusp);
      if (q < minp || (q < minp + 1e-8 && i-j > bestend-beststart)) {
	minp = q;
	beststart = j;
	bestend = i;
	/*
	p1 = gd_min( gd_max( oneminusp*sp, 1.0e-10 ), 1.0 );
	p2 = gd_min( gd_max( oneminusp*sympol2(sp,spp), 1.0e-10 ), 1.0 );
	p3 = gd_min( gd_max( oneminusp*sympol3(sp,spp,sppp), 1.0e-10 ), 1.0 );
	p4 = gd_min( gd_max( oneminusp*sympol4(sp,spp,sppp,spppp), 1.0e-10 ), 1.0 );
	*/
	//printf("Best minp=%f start/end=%u/%u\n",minp,j,i);
	//printf("p1-4: %f %f %f %f\n",oneminusp*sp,oneminusp*sympol2(sp,spp),oneminusp*sympol3(sp,spp,sppp),oneminusp*sympol4(sp,spp,sppp,spppp));
      }
    }
  }
  //int efflen = bestend-beststart;

  //printf("Longest stretch: %u  Phred %u    Minp %f\n",bestend-beststart,(int)(10*log(minp)/log(0.1)),minp);

  PyObject* retval = PyTuple_New(2);  
  PyTuple_SET_ITEM(retval, 0, PyInt_FromLong( (int)(10*log(minp)/log(0.1)) ) );

  // find 2nd interval
  int beststart2 = 0;
  int bestend2 = 0;
  //p1=p2=p3=p4=1.0;
  minp = 1.0;
  for (i=15; i<=qlen; i++) {
    sp=spp=sppp=spppp=sppppp=0;
    oneminusp = 1.0;
    //printf("i=%u\n",i);
    for (j=i-1; j>=0 && j>=i-52; j--) {
      // check if inside 1st interval; use intervals of starting positions [j,i-14) and [beststart,bestend-14)
      if ((i-14 > beststart) && (j < bestend-14)) continue;
      oneminusp *= 1.0 - truepa1[j];
      //printf("i=%u j=%u pa1[j]=%f oneminusp=%f\n",i,j,pa1[j],oneminusp);
      sp += pa1[j];
      spp += pa2[j];
      sppp += pa3[j];
      spppp += pa4[j];
      sppppp += pa5[j];
      e2 = sympol2(sp,spp); e3 = sympol3(sp,spp,sppp); e4 = sympol4(sp,spp,sppp,spppp);
      e5 = 1.0 - oneminusp*(1.0 + sp + e2 + e3 + e4);
      if (e5 > 0.01) {
	//printf ("i=%u j=%u sympol5=%f\n",i,j,sympol5(sp,spp,sppp,spppp,sppppp));
	break;
      }
      if (j > i-15) continue;
      // i-52 < j <= i-15
      k = 5*(i-15-j);
      q = oneminusp*mapstatscacheD[k];
      q += gd_max(0.0,1-q)*oneminusp*mapstatscacheD[k+1]*sp;
      q += gd_max(0.0,1-q)*oneminusp*mapstatscacheD[k+2]*e2;
      q += gd_max(0.0,1-q)*oneminusp*mapstatscacheD[k+3]*e3;
      q += gd_max(0.0,1-q)*oneminusp*mapstatscacheD[k+4]*e4;
      q += gd_max(0.0,1-q)*e5; // safety - 5 mutations always bad
      //printf("i=%u j=%u q=%f oneminusp=%f\n",i,j,q,oneminusp);
      if (q < minp || (q < minp + 1e-8 && i-j > bestend2-beststart2)) {
	minp = q;
	beststart2 = j;
	bestend2 = i;
	/*
	p1 = gd_min( gd_max( oneminusp*sp, 1.0e-10 ), 1.0 );
	p2 = gd_min( gd_max( oneminusp*sympol2(sp,spp), 1.0e-10 ), 1.0 );
	p3 = gd_min( gd_max( oneminusp*sympol3(sp,spp,sppp), 1.0e-10 ), 1.0 );
	p4 = gd_min( gd_max( oneminusp*sympol4(sp,spp,sppp,spppp), 1.0e-10 ), 1.0 );
	*/
	//printf("Best minp=%f start/end=%u/%u\n",minp,j,i);
	//printf("p1-4: %f %f %f %f\n",oneminusp*sp,oneminusp*sympol2(sp,spp),oneminusp*sympol3(sp,spp,sppp),oneminusp*sympol4(sp,spp,sppp,spppp));
      }
    }
  }

  //efflen = bestend2-beststart2;
  //printf("longest stretch %u\n",efflen);

  PyTuple_SET_ITEM(retval, 1, PyInt_FromLong( (int)(10*log(minp)/log(0.1)) ) );
  return retval;

}

  
		   
int add_maq_record( char* address, long long length, long index,
		    char* seq, int seqlen, char* quals,
		    int se_qual, int pe_qual, int alt_qual,
		    int nummismatches, int likelihood, 
		    int c0, int c1, int flag, int chromid, long chrompos, int strand, int distance, char* readid, int readidlen,
		    char storevariants, int64_t bitmap, char maqlong) {

  // Ah now, this is -so- ugly.  But what do I do.
  if (maqlong) {

    // check that the buffer is sufficient
    if ((index+1)*sizeof(maq_map_long_record_t) > length) return -1;
    // compute address
    maq_map_long_record_t* p = &((maq_map_long_record_t*)address)[ index ];
    // fill the record
    int i;
    unsigned char b;
    char c;
    for (i=0; i<MAQ_MAX_READLEN; i++) {
      if (i<seqlen) {
	c = seq[i] & 0xDF;
	switch (c) {
	case 'A': b=0; break;
	case 'C': b=1<<6; break;
	case 'T': b=3<<6; break;
	default: b=2<<6;   // this also catches 'N'
	}
	b |= (quals[i] - 33) & 0x3F;
      } else b = 0;
      p->seq[i] = b;
    }
    p->seq[MAQ_MAX_READLEN-1] = se_qual;
    p->length = seqlen;
    p->pe_map_qual = pe_qual;
    p->info1 = nummismatches & 0xF;  // what's in the upper 4 bits?
    p->info2 = likelihood;
    p->c0 = c0;
    p->c1 = c1;
    p->flag = flag;
    p->alt_qual = alt_qual;
    p->seqid = chromid;
    p->pos = (chrompos << 1) | strand;
    p->dist = distance;
    if (readidlen > MAQ_MAX_NAMELEN-1) {
      readid += MAQ_MAX_NAMELEN-1-readidlen;
      readidlen = MAQ_MAX_NAMELEN-1;
    }
    for (i=0; i<readidlen; i++) p->name[i] = readid[i];
    for (; i<MAQ_MAX_NAMELEN; i++) p->name[i] = 0;
    if (storevariants) {
      *((int64_t*)&p->seq[MAQ_MAX_READLEN-9]) = bitmap;
    }
    
  } else {
    
    // check that the buffer is sufficient
    if ((index+1)*sizeof(maq_map_record_t) > length) return -1;
    // compute address
    maq_map_record_t* p = &((maq_map_record_t*)address)[ index ];
    // fill the record
    int i;
    unsigned char b;
    char c;
    for (i=0; i<MAQ_MAX_READLEN; i++) {
      if (i<seqlen) {
	c = seq[i] & 0xDF;
	switch (c) {
	case 'A': b=0; break;
	case 'C': b=1<<6; break;
	case 'T': b=3<<6; break;
	default: b=2<<6;   // this also catches 'N'
	}
	b |= (quals[i] - 33) & 0x3F;
      } else b = 0;
      p->seq[i] = b;
    }
    p->seq[MAQ_MAX_READLEN-1] = se_qual;
    p->length = seqlen;
    p->pe_map_qual = pe_qual;
    p->info1 = nummismatches & 0xF;  // what's in the upper 4 bits?
    p->info2 = likelihood;
    p->c0 = c0;
    p->c1 = c1;
    p->flag = flag;
    p->alt_qual = alt_qual;
    p->seqid = chromid;
    p->pos = (chrompos << 1) | strand;
    p->dist = distance;
    if (readidlen > MAQ_MAX_NAMELEN-1) {
      readid += MAQ_MAX_NAMELEN-1-readidlen;
      readidlen = MAQ_MAX_NAMELEN-1;
    }
    for (i=0; i<readidlen; i++) p->name[i] = readid[i];
    for (; i<MAQ_MAX_NAMELEN; i++) p->name[i] = 0;
    if (storevariants) {
      *((int64_t*)&p->seq[MAQ_MAX_READLEN-9]) = bitmap;
    }
    
  }
  return 0;
}

// returns size of struct
int sort_maq_records( char * address, long long length, long index, char maqlong ) {

#define compare_maq(a,b) ( (templocus = (a)->seqid - (b)->seqid ) < 0 || (templocus == 0 && ((a)->pos < (b)->pos)) )

  if (maqlong) {
    // sanity check
    if (index*sizeof(maq_map_long_record_t) > length) return 0;
    int32_t templocus;
    QSORT( maq_map_long_record_t, ((maq_map_long_record_t*)address), index, compare_maq )
      
      return sizeof(maq_map_long_record_t);

  } else {
    // sanity check
    if (index*sizeof(maq_map_record_t) > length) return 0;
    int32_t templocus;
    QSORT( maq_map_record_t, ((maq_map_record_t*)address), index, compare_maq )
      
      return sizeof(maq_map_record_t);
  }
}


PyObject* get_maq_records( char * address, long long length, long idxstart, long numrecords, char maqlong ) {

  if (maqlong) {
    // sanity check
    if ((idxstart+numrecords)*sizeof(maq_map_long_record_t) > length) {
      Py_INCREF(Py_None);
      return Py_None;
    }
    PyObject* ret = PyString_FromStringAndSize( (char*)&(((maq_map_long_record_t*)address)[idxstart]),
						numrecords * sizeof(maq_map_long_record_t) );
    return ret;
  } else {
    // sanity check
    if ((idxstart+numrecords)*sizeof(maq_map_record_t) > length) {
      Py_INCREF(Py_None);
      return Py_None;
    }
    PyObject* ret = PyString_FromStringAndSize( (char*)&(((maq_map_record_t*)address)[idxstart]),
						numrecords * sizeof(maq_map_record_t) );
    return ret;

  }

}

  
  
