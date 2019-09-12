#include "Python.h"

/*****************************************************************************************************************
 This code is copyright (c) Gerton Lunter, Jan 2009
 It may not be distributed, made public, or used in other software without the permission of the copyright holder
******************************************************************************************************************/

// from mmapmodule.c

typedef struct {
  PyObject_HEAD
  char *  data;
  size_t  size;
} mmap_object_head;


struct mapping {
  char* readalignment;
  char* genomealignment;
  int64_t locus;
  long start;
  long end;
  int score;
  int hashmatches;
};

// struct and defs for Maq's .map file format

#define MAQ_MAGIC -1               
#define MAQ_MAX_READLEN 64        
#define MAQ_MAX_LONG_READLEN 128
#define MAQ_MAX_NAMELEN 36

typedef struct {
  char seq[MAQ_MAX_READLEN];       /* bit 6,7: ACGT; bit 0-5: quality score 0-63 */
  // int64_t mismatch_bitmap;      /* at seq[55]: ones at mismatch positions */
  // unsigned char se_map_qual;    /* at seq[63]: single-end map quality, or indel length (0 for no indel) */
  unsigned char length;            /* sequence length */
  unsigned char pe_map_qual;       /* paired-end map quality, or indel position (0 for no indel) */
  unsigned char info1;             /* bits 0-3: # mismatches in best hit; bits 4-7 --? ; 0xff for no match */
  unsigned char info2;             /* sum of qualities at mismatch positions; 0xff for no match */
  unsigned char c0;                /* number of 0-mismatch hits to first 24bp in reference */
  unsigned char c1;                /* number of 1-mismatch hits to first 24bp in reference */
  unsigned char flag;              /* bits 0-7: FF; FR (normal); RF; RR; paired; different chromosomes; no mate; smith-waterman aligned */
  unsigned char alt_qual;          /* single-end map quality of mate */
  int32_t seqid;                   /* chromosome id */
  int32_t pos;                     /* bit 0 contains strand (0=fw), bits 1-31 the position */
  int32_t dist;                    /* physical length of read; 1+distance of 5' ends on chromosome */
  char name[MAQ_MAX_NAMELEN];      /* 0-terminated read identifier; truncated at beginning if necessary */
} maq_map_record_t;


typedef struct {
  char seq[MAQ_MAX_LONG_READLEN];       /* bit 6,7: ACGT; bit 0-5: quality score 0-63 */
  // int64_t mismatch_bitmap;      /* at seq[55]: ones at mismatch positions */
  // unsigned char se_map_qual;    /* at seq[63]: single-end map quality, or indel length (0 for no indel) */
  unsigned char length;            /* sequence length */
  unsigned char pe_map_qual;       /* paired-end map quality, or indel position (0 for no indel) */
  unsigned char info1;             /* bits 0-3: # mismatches in best hit; bits 4-7 --? ; 0xff for no match */
  unsigned char info2;             /* sum of qualities at mismatch positions; 0xff for no match */
  unsigned char c0;                /* number of 0-mismatch hits to first 24bp in reference */
  unsigned char c1;                /* number of 1-mismatch hits to first 24bp in reference */
  unsigned char flag;              /* bits 0-7: FF; FR (normal); RF; RR; paired; different chromosomes; no mate; smith-waterman aligned */
  unsigned char alt_qual;          /* single-end map quality of mate */
  int32_t seqid;                   /* chromosome id */
  int32_t pos;                     /* bit 0 contains strand (0=fw), bits 1-31 the position */
  int32_t dist;                    /* physical length of read; 1+distance of 5' ends on chromosome */
  char name[MAQ_MAX_NAMELEN];      /* 0-terminated read identifier; truncated at beginning if necessary */
} maq_map_long_record_t;



#define MATH_DB -0.23025850929940456    /* log(0.1) / 10.0 */

typedef int32_t int32_type;
typedef int64_t int64_type;
typedef uint32_t uint32_type;
typedef uint64_t uint64_type;

PyObject* genome_slice( PyObject* mmap, int datastart, u_int64_t start, u_int64_t end );

long genome_hash_29( PyObject* mmap, int datastart, u_int64_t locus );

long extern genome_hash_29_fromstring( const char* data, u_int64_t locus );

int strand_bit_29( PyObject* mmap, int datastart, u_int64_t locus );

int extern strand_bit_29_fromstring( const char* data, u_int64_t locus );

int read_nuc_fromstring( char* data, u_int64_t uloc );

PyObject* set_nuc_fromstring( char* data, int size, u_int64_t uloc, int nuc );

void set_bit_mmap( PyObject* mmap, unsigned long offset, u_int64_t uloc, int bitstate );

PyObject* read_unsigned_long_frommmap( PyObject* mmap, int index );

long fingerprint_fromstring( const char* data, u_int64_t uloc, int numnucnibbles );

long fingerprint15_fromstring( const char* data, u_int64_t uloc );

long fingerprint( PyObject* mmap, int datastart, u_int64_t locuc, int numnucnibbles );

PyObject* align_candidate( int64_t locus, int start, int end, char* sequence, int overhang, PyObject* mmap, int datastart, int64_t maxlen, int gapstartpenalty, PyObject* scoretable, int band );

int score_alignment( char* readaln, char* genomealn, PyObject* quality, int gapopen, int gapextend );

PyObject* locate_align_score_candidates( char* sequence, /* char* seqbitstring, int seqbitstringlen,  */
					 PyObject* hashtablemmap, PyObject* genomemmap, int datastart, 
					 long mask, char tryvariants, int maxfingerprintvariants, int maxdistance, int minhashmatches,
					 int64_t genomesize, int gapstartpenalty, PyObject* scoretable,
					 PyObject* quality, int gapopen_score, int gapextend_score, int nucprior_score, 
					 int lowqthreshold, int band, int max_candidates, int64_t locus_to_consider);

PyObject* map_region( char* sequence, int64_t start, int size, PyObject* genomemmap, int datastart, int64_t genomesize, PyObject* scoretable, int gapstart_map, int gapopen_score, int gapextend_score, PyObject* quality );

PyObject* calc_mutprobs( PyObject* quals, PyObject* rmask, int lowqthreshold );

PyObject* calc_mutprobs2( PyObject* quals, PyObject* rmask, int lowqthreshold );

int add_maq_record( char* address, long long length, long index,
		    char* seq, int seqlen, char* quals,
		    int se_qual, int pe_qual, int alt_qual,
		    int nummismatches, int likelihood, 
		    int c0, int c1, int flag, int chromid, long chrompos, int strand, int distance, char* readid, int readidlen,
		    char storevariants, int64_t bitmap, char maqlong);

int sort_maq_records( char * address, long long length, long index, char maqlong );

PyObject* get_maq_records( char * address, long long length, long idxstart, long numindices, char maqlong );

PyObject* rev_comp( char*, int );

