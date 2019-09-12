Notes for Stampy v1.0.31   --   Gerton Lunter, December 2016
------------------------------------------------------------



1. Summary 
==========

Stampy has the following features:

- Maps single, paired-end, mate pair Illumina reads to a reference
- Fast: about 10 (with BWA) or 15 hours (without) per Gbase
- Low memory footprint: 2.7 Gb shared memory for a 3Gbase genome 
- High sensitivity for indels and divergent reads, up to 10-15% 
- Low mapping bias for reads with SNPs or indels
- Well calibrated mapping quality scores 
- Strips adapters to improve miRNA read mapping
- Input: Fastq and Fasta; gzipped or plain; SAM and BAM
- Output: SAM, Maq's map file 
- Optionally calculates per-base alignment posteriors 
- Optionally processes part of the input 
- Handles reads up to 4500 bases

Bug reports are highly appreciated.  If you can, please design a 
small test file that reproduces the bug, and email the precise 
command line and details of the system you ran the program on.  
However, please do read this documentation before submitting a bug 
report.



2. Building 
===========

Just type "make"

Currently only the linux x86_64 platform is supported.

Stampy needs Python version 2.7.  Both 2-byte and 4-byte Unicode 
encodings are supported.

Most errors that occur at this stage are related to the Python 
installation.  Here is a checklist in case something goes wrong:

- Check that python version 2.7 is installed on your system (type 
  "python" in a shell).
- Check that the executables python2.7 and the related python2.7-config
  can be found in your path, and that they live in the same directory 
  as the Python executable you use.
- If you get linking errors, check that the -L<path> bit of the output 
  of "python-config--ldflags" points to a directory that contains 
  libpython.  If not, your python installation is faulty; setting 
  paths properly or re-installing Python might resolve the problem.



3. Quick-start guide
====================


Build a genome (.stidx) file:

     ./stampy.py --species=human --assembly=hg18_ncbi36 \
                 -G hg18 /data/genomes/hg18/*.fa.gz

Build a hash (.sthash) file:

     ./stampy.py -g hg18 -H hg18

Single-end mapping:

     ./stampy.py -g hg18 -h hg18 -M reads.fastq.gz

Paired-end mapping:

    ./stampy.py -g hg18 -h hg18 -M reads_1.fastq reads_2.fastq

To speed up mapping, use BWA (recommended) and multithreading:

    bwa aln -q10 -t8 hg18 reads_1.fastq > 1.sai
    bwa aln -q10 -t8 hg18 reads_2.fastq > 2.sai
    bwa sampe hg18 1.sai 2.sai reads_1.fastq reads_2.fastq | \
        samtools view -Sb - > bwa.bam

    ./stampy.py -g hg18 -h hg18 -t8 --bamkeepgoodreads -M bwa.bam

Note that this requires the fastq files to be in Sanger (base-64)
format; use the -A fastq conversion command if necessary.  Note 
also that it is not necessary to sort the BWA BAM file, but it is
allowed; it is possible to re-map from (one or more) sorted BAM 
file(s) directly, and Stampy will pair up reads by label name 
before mapping:

   ./stampy.py -g hg18 -h hg18 -M sample.BAM sample_unmapped.BAM

Set the divergence for mapping to a foreign reference:

     ./stampy.py -g hg18 -h hg18 --substitutionrate=0.05 \
     		 -M reads.fastq.gz

Set the initial insert size distribution -- only inserts within 
3sd of mean are considered for training the actual size distribution:

    ./stampy.py -g hg18 -h hg18 --insertsize=400 --insertsd=75 \
    		-M reads_1.fastq reads_2.fastq

Use (post v1.3) Solexa quality scores; default is Sanger qualities:

     ./stampy.py -g hg18 -h hg18 --solexa -M reads.fastq.gz

Process the 3rd of each set of 8 reads or read pairs:

     ./stampy.py -g hg18 -h hg18 --processpart=3/8 \
     		 -M reads.fastq.gz

To process miRNA libraries, first strip the adapters from the reads,
then map as usual.  Stampy automatically identifies the adapter 
sequences, and "strips" them by setting base qualities to 1.

     ./stampy.py --solexa --strip-adapter -A reads_1.fq.gz | gzip > r1_strip.fq.gz
     ./stampy.py --solexa --strip-adapter -A reads_2.fq.gz | gzip > r2_strip.fq.gz
     ./stampy.py -g hg18 -h hg18 -M r1_strip.fq.gz r2_strip.fq.gz

After stripping you may want to first map with BWA to speed things up.
The -A command outputs Sanger (base-33) fastq files, ready for 
processing with BWA.



4. Running Stampy 
=================

First you have to build a genome file:

      ./stampy.py -G hg18 /data/genomes/hg18/*.fa.gz

You may provide .fa or gzipped .fa files, and the .fa files may
contain more than one chromosome or contig.  As chromosome or contig
identifer, the first word on the Fasta header line is taken.  If the
header follows NCBI formatting, the "ref" field is taken.  The genome
file has extension .stidx; building this takes a few minutes for large
genomes.

After building the genome file, you need to build a hash table:

      ./stampy.py -g hg18 -H hg18

This takes about 10 minutes for large genomes, and produces a file
hg18.sthash.  Finally, you're ready to map some reads:

     ./stampy.py -g hg18 -h hg18 -M illuminareads.fastq.gz

Note that the hash file may be large (up to 2 Gb); to improve startup
speed you may want to keep it on a local drive.  Note also that the
file system must support memory mapped files for Stampy to work; NFS
does, but e.g. GlusterFS does not.  Some file systems support memory
mapped files but become exceedingly slow; if this happens try moving
both files to a local drive.

You need not unzip the input fastq files.  You can use .fasta files
too, using the option --inputformat=fasta.  By default Stampy assumes
that the quality scores use Sanger encoding (score 0 is '!', ascii
33), and use the phred (=log p) scale, not the logit (=log p/(1-p))
scale.  

Paired-end mapping is done by supplying two fastq files, one for each
mate:

    ./stampy.py -g hg18 -h hg18 -M solexareads_1.fastq \
                solexareads_2.fastq

(If -M is not the last option on the command line, the input files
need to be separated by a comma rather than a space.)  It is not
required to explicitly set parameters for the insert size
distribution; Stampy automatically determines these within the first
few hundred mappings, and adapts the scoring accordingly.  The default
settings are an average separation of 250 and standard deviation 60,
which is broad enough to capture all currently used libraries within
the default 4 standard deviations.  If a library with an insert size
distribution outside this range is used, you need to change these
defaults for autocalibration to work.



5.  Faster mapping with BWA 
============================

For large reference genomes and not very divergent data sets (<1% say), 
it is highly recommended to use BWA to speed up mapping.  First, create 
a BWA index for the same reference genome as used by Stampy, then map
the reads using BWA in the usual way, and use samtools to convert the
resulting SAM file into a BAM file.  Then re-map the BAM file using 
Stampy, but keep the well-mapped reads:

    ./stampy.py -g hg18 -h hg18 --bamkeepgoodreads -M bwa.bam

To further speed up mapping, tell both BWA and Stampy to use multiple
threads using the -t option.

Previous version of Stampy launched BWA itself, and processed reads in
batches.  For backward compatibility this is still supported:

    ./stampy.py --bwaoptions="-q10 BWAindex/hg18.fa" -g hg18 
             	-h hg18 -M solexareads_1.fastq solexareads_2.fastq

However this mode is not recommended: it takes longer because BWA has to
be launched many times, and the BWA and Stampy tables are needed
simultaneously so that the memory requirement is about twice as high.
Finally, this mode only works with BWA versions 0.5.6, 0.5.7 and 0.5.9;
the 0.6.x versions don't work as they occasionally segfault when BWA is
given a small number of reads to process.

The achieved speed-up depends on the quality of the input reads; if
BWA fails to map a large proportion of reads, the speed-up will be
comparatively low.  

On some filesystems, Stampy is particularly slow.  In order to conserve 
main memory, Stampy shares its two large tables across multiple 
instances running on the same shared-memory node.  This is achieved by 
memory-mapping a shared file.  On standard filesystems (NFS, Linux 
etx3/4) this works fine, but certain higher-end distributed file 
systems this causes excessive slow-down.  The solution is to copy the 
.stidx and .sthash files to a local drive, e.g. /tmp, and share those 
copies locally rather than through the globally shared file system.



6.  Mapping divergent reads 
===========================

To calculate correct mapping qualities, Stampy needs to know the
expected divergence from the reference.  This is set with the
--substitutionrate= option.  The default is 0.001 substitutions per
site.

Increasing the read length, and using paired-end reads, helps mapping
divergent reads.  The following table gives an indication of the
divergence at which a reasonable proportion of reads can be correctly
mapped.  These numbers were obtained by simulation, using the human
genome as reference, and should be taken as an indication only; they
are dependent on error rates, the repetitiveness of the genome, the
insert size distribution, and local variations in divergence; in
addition no indel mutations were included.

		   36bp   36bp   72bp   72bp
      divergence | single paired single paired
      -------------------------------------------------------
      0%         | 82%    95%    87%    96%
      3%         | 73%    91%    80%    94%
      6%         | 60%    83%    72%    92%
      9%         | 41%    56%    56%    88%
      12%        | 28%    51%    48%    80%



7.  Mapping mate pair libraries
===============================

Mate pair libraries contain a mixture of ordinary (--> <--) paired-end
reads, and mate pairs (<-- -->).  To map these, Stampy trains separate
insert size distributions for the two types of pairs, and attempts 
realignment within both regions.

To allow Stampy to accurately train the two distributions, both
distributions need to be seeded properly.  This is best done by mapping 
a small fraction of the data, starting with a very wide distribution:

  ./stampy.py -g ref -h ref --numrecords=5000 \
              --insertsize2=-2000 --insertsd2=750 -M read1.fq read2.fq

Stampy trains a distribution with reads that fall within 3 standard 
deviations of the mean, using the current estimates of each.  Stampy 
reports the estimated distributions when its done:

  stampy: # Paired-end insert size: 266.4 +/- 91.0  (903 pairs)
  stampy: # Mate pair insert size: -1519.1 +/- 164.5  (3374 pairs)

The values -1519 and -165 can be used as new initial parameters, and
the procedure repeated to obtain more accurate approximations.

Note that the mean mate-pair insert size is negative, corresponding to
the fact that the read mapping to the reverse strand maps to the left
of the forward-mapping read.  Protocols other than the currently
standard Illumina mate-pair protocol may need positive insert sizes.

Note also that with the current protocol, the standard deviation of
the mate pair insert size distribution is always larger than that of the
paired-end distribution, and that wider distributions slow down the
realignment process, so that mate pair mapping can be slow.

It is not recommended to use BWA for pre-mapping when mapping mate
pair reads, as BWA does not consider the alternative positions, causing
incorrect mappings in certain cases.



8.  Reporting alternative mapping positions
===========================================

Stampy can output alternative mapping locations for reads and read pairs,
using the XA tag.  By default this is switched off; it can be enabled with
the options

  --xa-max=3 --xa-max-discordant=10

to report at most 3 alternative placements for single reads and concordant
read pairs, and at most 10 for discordant read pairs (these are the BWA
defaults).

The XA tag reports, for each additional hit, the chromosome name, position 
and strand, CIGAR string, and the number of base mismatches plus the total 
length of any insertions or deletions.

Note that when using BWA for pre-mapping, there is no need to set BWA's -n 
or -N options; Stampy automatically sets these to the correct values.


9.  Testing Stampy
==================


Stampy includes test code that generates reads under an empirical read
error model and introduce SNPs van indel variants.  From the output,
stampy estimates the sensitivity for mapping reads back to their
correct location, and it tabulates results by mapping posterior score
to see if these are well calibrated.

To run these tests, you need to create a genome file and corresponding
hash file, and provide one or two .fastq files with short reads and
quality scores.  Stampy will use the length and qualities of these
reads, but generate sequences from random locations in the genome
provided.

You run a test by using the -S command to simulate reads, followed by a 
-P command to parse the .sam output; or run the two together using -T:

    ./stampy.py -g hg18 -h hg18 -T solexareads_1.fastq,solexareads_2.fastq

This uses paired-end reads; single-end reads are used if just one
.fastq file is provided.  The following options are useful in test
mode:

   --substitutionrate=S		  Introduce an expected fraction S of 
                                  Poisson-distributed substitutions 
				  (default: 0.001)
   --insertsize=N		  Set the mean insert size for 
				  paired-end reads (default: 250)
   --insertsd=N                   Set the standard deviation of the 
				  insert size distribution (default: 60)
   --numrecords=N		  Only map the first N reads, or read 
				  pairs (default: all)
   --simulate-minindellen=N       Set the lower bound for simulated indel 
				  lengths (default: 0)
   --simulate-maxindellen=N       Set the upper bound for simulated indel 
				  lengths (default: 0)
   --simulate-duplications        Introduce duplications rather than 
				  insertions of random sequence (default)
   --simulate-numsubstitutions=N  Introduce N substitutions in each read, 
				  rather than a Poisson-distributed number

The settings of the first three options are used in simulations, and
again when computing mapping statistics from the mapped reads.  These
options are therefore also meaningful outside of testing.

The default for librarysd is generous, to train the model on a wide
range of input data.  For simulations this needs to be adjusted.

Simulated indels are generated by drawing the indel length from a
uniform distribution.  This is not meant to be close to a real
distribution, but is useful for testing the behaviour of stampy
conditional on indels being present.  Negative indel length are
deletions, positive ones are insertions.  When single-end reads are
used, each read will contain a simulated indel; for paired-end reads,
one of the mate pairs will contain an indel.



10.  Mapping output 
===================


Output is written to stdout by default.  An output file can be chosen
with the -o or --output= option, however the SAM format is now the
accepted standard and use of any other format is discouraged.  

A number of formats can be chosen using th -f (or --format=) option:

 -f sam            :  SAM format (default)
 -f maqtxt         :  Maq's text output format (produced by 'maq mapview')
 -f maqmap         :  Maq's binary .map format, new version (long reads)
 -f maqmapShort    :  Maq's binary .map format, old version (short reads)
 -f maqmapShortN   :  Maq's binary .map format, old version, including 
		      variant positions (produced with 'maq map -N')

By default all reads are represented in the output.  The default
output format is the SAM format.  See below for details on the various
formats.

The SAM format is the most comprehensive format, and is recommended.
The Samtools program (samtools.sourceforge.net) is recommended for
dealing with .SAM files.

Several tools have been developed to use Maq's .map format as well,
including those that are included with the Maq package, and therefore
the .maq format was included as a convenience.  However, Maq's .map
format cannot represent all useful information.  In particular indels
are better represented in the SAM format.



10.1.  A note on likelihoods and posteriors

The single most often used statistic to judge the trustworthiness of a
read map location is its "mapping quality".  This is an approximation
of the probability that a read is mapped to the wrong location
(represented as a Phred score).

The probabilistic model used by Stampy is a hybrid of three models:

  (1) a Bayesian model, which considers all candidate locations
      weighted by their likelihood.
  (2) an error model, which considers the possibility that read errors
      cause reads to be incorrectly mapped
  (3) a random model, which predicts how well a random sequence would
      match to the genome

The first model deals with errors due to repetitive and
nearly-repetitive sequence, and assumes that the correct mapping
location was considered among the candidates.  The second model
estimates the probability that the correct candidate was missed
because of (single-nucleotide) read errors.  The third model acts as a
post-hoc filter, and assesses whether a candidate locations looks
better than a random best match.

Together these models capture most of the error modes of read mapping.
The most obvious exception is that the error model does not consider
indel errors or mutations; these often lead to the correct candidate
location being missed, particularly for short single-end reads.  The
"best" map that results is often caught by the third model; however if
the sequence is mildly repetitive, it will also pass this filter.

As a result, mapping quality is well calibrated for almost all cases,
except short single-end maps of reads that contain indel mutations, in
which case the mapping quality is overly optimistic by about an order
of magnitude.  Consequently, indels that are supported by single-end
maps only should be treated with caution.

Stampy computes posteriors and likelihoods both for pairs of reads,
and for reads considered by themselves.  The following table
summarizes the SAM tags for these statistics for paired reads:


	  Read                     | Posterior   | Likelihood |
	  -------------------------+-------------+------------+
	  Pair                     | MAPQ column | PQ:i:      |
	  This read as single read | SM:i:       | UQ:i:      |
	  Mate as single read      | MQ:i:       | XQ:i:      |

For single reads, only the MAPQ column and SM:i: tag are present.



10.2.  SAM output format

This tab-delimited output format is described in detail on
samtools.sourceforge.net.  A few things to note:

- The MAPQ field is the phred-scaled estimated probability that the
  read was mapped to the wrong location.  See 8.1 for
  details.

- In addition to the mapping quality (the probability that a read
  was mapped incorrectly), Stampy also reports the read likelihood: 
  the likelihood that a read was produced from the reference, 
  conditional on the mapping location being correct.  This score is 
  the sum of phred qualities on mismatching sites, and includes 
  probabilities for indels and read separation as well.  The single-
  read, paired-read and mate likelihoods are reported in the
  optional "UQ", "PQ" and "XQ" fields.  See 8.1 for details.

- The SAM format requires that paired reads share identifiers; if a
  trailer like "/1" is present, it will be removed from the identifier
  to conform to the standard.

- If an identifier contains spaces, only the first word will be used
  (and paired-end /1, /2 trailers silently added if necessary).  The
  option --keeplabel changes this behaviour, and instead replaces spaces
  by underscores.  It is not possible to keep spaces, since BWA also 
  only uses the first word.

- The "proper pair" flag bit (value 2) is set if two reads are
  correctly oriented, and their separation is within 5 standard
  deviations from the mean

- The "mate unmapped" bit (value 8) is never set by itself; pairs are
  either both mapped or both unmapped

- The "UQ" optional field (single read likelihood) is always present.

- The "PQ" and "XQ" fields are always present for paired-end reads,
  and represent the paired likelihood (which includes terms for
  mismatches, indels and read separation), and the single read
  posterior, respectively.

- The optional "SM" and "MQ" fields (mapping quality of the read or
  its mate, considered as a single read) are always present for
  paired-end reads



11. Hints and features


11.1 Ordering of chromosomes in SAM output

  The SAM/BAM format does not specify the order in which the reference
  sequences appear in the SAM header.  By default Stampy up to version
  1.0.13 ordered them alphabetically; from 1.0.14 the default order
  is as the original .fa input file.

  Some downstream tools require the references to appear in the same
  order as in the original .fa input file.  To enable this behaviour,
  you had to use --keepreforder option while mapping before version 
  1.0.14; this is now the default.  To re-enable alphabetic mapping,
  use --alphabetical.

  Note that .stidx files created before Stampy v1.0.4 will not work
  correctly; you need to re-create the .stidx file.

11.2 Parsing NCBI fasta files

  NCBI fasta files use a >gi|nnn|ref|xxx identifier.  By default
  Stampy parses this and uses only the "xxx" part.  Use --noparseNCBI
  to switch off this behaviour and use the full NCBI identifier.

11.3 Base-level alignment qualities

  The placement of an indel into an alignment is uncertain, because of
  actual ambiguity, polymorphisms, and read errors.  A probabilistic
  alignment considers all possibilities (under a suitable model) and
  this information can be used to compute a posterior score.  To 
  have Stampy compute this, use the --alignquals option.

  Note that this increases the output file size, and also increases 
  the runtime, since the required forward and backward iterations 
  take time.

11.4 Base quality recalibration

  For very old Solexa/Illumina reads it may be advantageous to 
  recalibrate the quality scores:

     ./stampy.py -g hg18 -h hg18 -R solexareads.fastq.gz

  By default this maps 1 percent of the fastq file onto the reference,
  collects statistics, and writes a *.recaldata file.  If you now start 
  the mapper again,

    ./stampy.py -g hg18 -h hg18 -M solexareads.fastq.gz,

  Stampy will use the .recaldata file and apply the recalibration before
  mapping, which should improve the mapping quality statistics.  (Note
  that if the original file is not ascii-33 based, or uses logit scores,
  you need to provide the relevant options both times.)

  For current Illumina data, the quality scores are sufficiently well 
  calibrated, and will degrade by this recalibration procedure because
  polymorphisms are not taken into account.

11.5 Multiple mapping locations

  Normally, when Stampy identifies several equally good mapping 
  locations for a read or read pair, it reports one of these at random
  (and assigns the choice a low mapping quality, of 3 or less).

  Alternatively Stampy can report a limited number of alternative
  mapping locations if you set the --xa-max option to a nonzero
  value.  The alternative locations are reported in a single XA:Z:...
  tag, using the format that is used by BWA.

  Normally, no alternative mappings are reported for discordant read
  pairs. However, if you set --xa-max-discordant to a nonzero value, 
  these will be reported as well.

11.6 Marking BWA-mapped reads

  The option --bwamark tells Stampy to mark reads that were mapped by
  BWA, and directly copied by Stampy, with a XP:Z:BWA tag.  This 
  feature causes the output file to increase in size, however it is
  more efficient than keeping separate BAM files for Stampy and BWA.

  In cases where Stampy re-considers reads that were also mapped by
  BWA (for instance, if these overlapped a large indel), such reads
  are reported as unmapped (flag value 4) and as secondary (flag
  value 256); the original flag value is reported as XP:Z:BWA:<value>.
  In this case, the read mapped by Stampy is reported without XP:Z:BWA
  tag.

  Downstream tools should ignore reads that are marked as either
  unmapped or secondary.  Since BWA-specific reads have both flags
  set, there is a good chance that such tools will indeed do so.
  However since many mappers don't include unmapped reads or secondary
  mappings, it is possible that downstream tools ignore these flags.
  Please be aware of this when using this feature.

  To reconstruct the BWA output, keep all reads with an XP:Z:BWA tag,
  reconstructing the flag value is necessary.

  Note that this feature is experimental, is not very well tested, and
  will for the moment not be supported.

11.7 Mapping from BAM files

  BAM files are useful containers for storing and archiving both
  mapping results and original read data.  Stampy has a few options
  to simplify re-mapping reads contained in BAM files.

  Sometimes, reads from a sequencing experiment are split between 
  multiple BAM files, to store them by chromosome for example, or
  unmapped reads may be stored in a separate BAM file.  Stampy can
  map multiple BAM files simultaneously, and will merge the contents
  before mapping.

  BAMs often contain data from multiple sequencing runs, and/or multiple
  libraries.  It is important to separate the data from different 
  libraries, since these often have different insert size distributions.
  To only re-map data from (a) particular read group(s), use

    --readgroup=ID:xxx,ID=yyy,...

  where xxx, yyy are the read group identifiers.  To map all reads from 
  a single library across readgroups, use the same option but specify
  the library identifier.  In the output, the reads will be assigned
  to a readgroup with the library identifier.

  If reads are mapped with BWA, a considerable speed-up can be achieved
  with the --bamkeepgoodreads option, which copies the well-mapped reads 
  to the output without remapping them.

  Stampy pairs up the reads from a single fragment before re-mapping them
  using a generous read-ahead window, allowing most reads to be paired
  on the fly.  Those mapping to very different locations are spilled to
  a temporary file, sorted, and mapped when all input has been read.
  Occasionally, this temporary file may become large.  To avoid problems,
  you can specify a location for this temporary file using

    --bamsortprefix=S

  You can limit the amount of memory used when sorting by specifying

    --bamsortmemory=N

  This value is passed to samtools via its -m option.


12. License
===========

This is a release version.  Permission is granted for the normal 
use of the program and its output in an academic setting, including
in publications.  If the program is used to generate data for a 
publication, please cite this paper:

  G. Lunter and M. Goodson.  Stampy: A statistical algorithm for 
  sensitive and fast mapping of Illumina sequence reads.  Genome
  Res. 2011 21:936-939.

The program itself may not be modified in any way, and may not be
reverse-engineered.

This license does not allow the use of this program for any 
commercial purpose.  If you wish to use this program for commercial 
purposes, please contact the author.

No guarantees are given as to the program's correctness, or the 
accuracy or completeness of its output.  The author accepts no 
liability for damage or otherwise following from using and 
interpreting the output of this program.  The software is supplied 
"as is", without obligation by the author to provide any services 
or support.



13. Revision history
====================


1.0.5, rev 786 (1 October 2010)  
       - first release

1.0.6, rev 803 (20 October 2010)
       - added NM (number of mismatches) tag
       - added support for sorted BAM/SAM input
       - bugfix: unmapped reads were sometimes reversed

1.0.7, rev 827 (3 November 2010)
       - added FAQ section to README
       - bugfix: BWA support was broken in 1.0.6

1.0.8, rev 828 (4 November 2010)
       - bugfix: occasional undefined variable reference

1.0.9, rev 852 (22 November 2010)
       - added --ignore-improper-pairs option, useful when
         remapping broken BAM files
       - bugfix: memory leak when mapping with BWA
       - bugfix: improper read pairing when remapping part 
         of a paired BAM file 
       - bugfix: v1.0.7/8 occasionally printed debug output

1.0.10, rev 854 (23 November 2010)
       - added --bwatmpdir option
       - better adherence to SAM v1.3 format specification;
         added option --referenceuri
       - bugfix: BWA paired-end mapping broken in 1.0.9

1.0.11, rev 880 (22 December 2010)
       - Faster, particularly for longer reads
       - Improved reporting of errors
       - Base Alignment Quality scores (BAQ, --baq)
       - Improved scoring of low-quality inserted bases
       - Improved indexing for references with many contigs
       - Fixed a bug requiring index and hash to be writeable
       - Fixed a bug throwing occasional segmentation faults
         for references with many contigs

1.0.12, rev 1010 (20 March 2011)
       - Improved specificity in low complexity regions
       - More informative error messages
       - Bugfix: underflow could cause negative BAQ values
       - Bugfix: all-N reads sometimes caused segfaults
       - Bugfix: -P sometimes failed to parse read label
       - Bugfix: Cigars sometimes started with I..D..

1.0.13, rev 1160 (16 June 2011)
       - Added support for mate pairs
       - Added XA (alternative mappings) output tag
       - Added optional tagging of BWA-mapped reads
       - Bugfix: simulation of duplications did not work
       - Bugfix: estimation of insert sizes improved

1.0.14, rev 1307 (9 December 2011)
       - BWA tags allow reconstruction of BWA-mapped reads
       - bz2 input supported
       - Bugfix: highly repetitive reads occasionally were
                 assigned high mapping qualities
       - Read pairs consisting of two identical sequences 
         are now mapped as a single read (for GenomeSTRiP)
       - Option --keepreforder now default; added --alphabetical
       - Added option --maxpairseeds to allow increasing the
         sensitivity for mapping to highly repetitive genomes
       - The 'properly paired' flag is now set for insert 
         sizes within 3 sd of the mean; this used to be 4.
       - Added option --overwrite
       - Fixed bug when read pairs were labeled as /1 and /3
       - Fixed bug in auto-sensing code for FASTA files

1.0.15, rev 1360 (2 February 2012)
       - Added --casava8 option to keep (but not map) reads 
               not passing QC filters
       - Bugfix: --overwrite now also works for .sam output
       - Bugfix: occassional math domain errors
       - Bugfix: the PQ tag value was float rather than int

1.0.16, rev 1430 (16 February 2012)
       - Bugfix: memory leak (present from revision 1307)
       - Bugfix: now does not require Python bz2 module

1.0.17, rev 1481 (4 April 2012)
       - Bugfix: --bwamark did not properly copy all flag bits
       - Bugfix: mapping positions of reads starting with adjacent 
               insertions and deletions were sometimes slightly off
       - Changed computation of ISIZE field to match BWA's sign and
               value for nonstandard directions and overlap
       - Changed computation of "NM" tag to agree with Samtools:
               count now includes bases involved in indels
       - Added option --labelfilter to map a subset of reads
       - Added option --maxbasequal to accept high base qualities
       - Added option --gatkcigarworkaround removing adjacent I/D
               events from CIGAR strings, which trips up GATK
       - Option --casava8 with --keeplabel keeps multiplex tag

1.0.18, rev 1526 (17 May 2012)
        - Bugfix: BWA-mapped reads lost their fragment size sign
	- Regression: simulation stopped working in 1.0.17
	- Added informative error message when mem-mapping fails

1.0.19, rev 1617 (10 October 2012)
        - Added capability to map from arbitrary BAM input
        - Added multithreading support
        - Added option --separator[1,2] to deal with alternatively
                formatted Casava v1.8 tags
        - Bugfix: issue with numeric-only read labels

1.0.20, rev 1642 (18 October 2012)
        - Increased specificity for reads mapping to highly 
                repetitive regions
        - Added --xcigar option for extended CIGAR output (XC tag)
        - Added --untrimq option to tweak soft trimming removal
                for BWA/BAM input
        - Bugfix: (--bamkeepgoodreads) Reads never trimmed, causing
                slowdown
        - Bugfix: (--bamkeepgoodreads) Sometimes accepted large 
                insert sizes
        - Bugfix: (-M bam) Meaningful error message when re-mapping 
                from incorrect reference
        - Bugfix: (-M bam) --readgroup can be used to annotate reads
                without readgroup annotation

1.0.21, rev 1715 (3 January 2013)
        - Speed improvements for unzipping, parsing and formatting;
                greatest effect on runs with high threading
	- Added -A command to convert fastq qualities, and 
                optionally strip adapters from miRNA reads
	- Bugfix: (-M bam) Reads mapped to chromosome boundary could
	        occasionally be mapped outside chromosome
        - Bugfix: Isize sometimes computed incorrectly
	- Bugfix: (-M bam) Lost MAPQ values from single-end BAM input
        - Bugfix: (-M bam) Lost 2 characters from single-end read label
	- Bugfix: (-M bam) Fixed crash when mapping 0 reads
	- Bugfix: (--bwaoptions) Lower-case fastq caused problems

1.0.22, rev 1848 (5 June 2013)
	- Added check for well-mated pairs in certain run modes
	- Added saner error message for out-of-memory
	- Bugfix: Stampy required input BAM reference to match
                mapping reference, even without --bamkeepgoodreads
        - Bugfix: Did not recognize .BAM files without .BAM extension
	- Bugfix: Unique @PG IDs when re-mapping from Stampy BAMs

1.0.23, rev 2049 (13 Dec 2013)
	- Bugfix: Ignore non-primary alignments in BAM from BWA mem
	- Bugfix: Preserves @PG records in input BAM
	- Bugfix: Warn about and ignore duplicate read labels in BAM

1.0.24, rev 3183 (1 Dec 2014)
	- Bugfix: Crash on BWA pre-mapping of SE reads from BAM input
	- Bugfix: "Internal error" messages for large genomes.
	- Bugfix: SAM parsing issues re-mapping Stampy-mapped BAMs
	- Improvement: can now handle genomes up to approximately 5.3 Gb
	- Improvement: Now uses BWA mapping as "hint" amongst other
	   candidates when mapping from premapped BWA BAMs, improving 
	   performance in small number of highly repetitive regions.

1.0.25, rev 3363 (27 March 2015)
        - Bugfix: threading and --bamkeepgoodreads would fail

1.0.26, rev 3393 (22 April 2015)
        - Bugfix: threading and remapping unmapped reads could give
           'too many values to unpack' error

1.0.27, rev 3417 (18 May 2015)
        - Bugfix: non-uniquely mapping reads were occasionally
          mapped to the wrong location and with wrong MAPQ values

1.0.28, rev 3475 (3 August 2015)
	- Bugfix: now deals correctly with '-' in fasta files which
	  e.g. exist in human genome build hg8.

1.0.29, rev 3755 (16 July 2016)
	- Bugfix: now handles mapping hints while re-mapping BAMs
	  that were mapped to different reference genomes

1.0.30, rev 3757 (26 October 2016)
        - Stampy now works with samtools version 1.3, as well as
          with older versions

1.0.31, rev 3759 (1 December 2016)
        - Fixed linking problem in newer Ubuntu version
        - Fixed clang compiler issues; now Mac-compatible again

1.0.32, rev 3761 (6 November 2017)
        - Fixed issue with old version of samtools



12. FAQ
=======

1. Q: Stampy complains that there are no input files, if I don't 
      use a comma between the file names
   A: Use -M as the last option on the command line to avoid this.


2. Q: Stampy starts but doesn't produce data for hours - is this
      normal?
   A: No. Stampy loads its index files in memory, and uses memory
      mapping to allow the memory to be shared between instances
      of Stampy.  However, this does not work on all file systems
      To avoid the problem, copy the .index and .hash files to a 
      local drive, e.g. /tmp


3. Q: Using Stampy with BWA seems to be slow, and it's supposed
      to be fast?
   A1: Running BWA and Stampy together (using --bwaoptions) is
      not recommended; running BWA first and re-mapping the
      resulting BAM is faster and uses less memory.  But if you
      must:
   A2: Make sure you actually run BWA (use --bwaoptions="...")
   A3: You might be running out of memory.  BWA and Stampy each 
      need 3 Gb of memory.  To run n copies side-by-side, and 
      mapping to a human-size genome, 4+3n Gb of free memory is 
      recommended.


4. Q: Stampy exits with "Fatal Python error: PyImport_GetModuleDict:
      no module dictionary".  Should I worry?
   A: This is an obscure Python bug that occurs rarely at shutdown.  
      As long as Stampy says "Stampy: Done", it finished OK.


5. Q: Stampy cannot find "python*-config.py" during installation;
      how can I solve this?
   A: This happens on Ubuntu installations, where python*-config is
      not always installed.  Try installing the "python-dev" 
      package.  If that fails, try installing python from source 
      (be sure to do "make install", not just "make")


6. Q: I'm getting "Suspiciously high Q scores" warnings, but I'm
      sure the fastq data is fine.
   A: Your fastq data is in Solexa format (base-64), while Stampy
      defaults to Sanger format (base-33).  Use the option --solexa.


7. Q: Stampy requires python 2.6 -- how do I choose this version?
   A: Ensure Python 2.6 is installed; then use "python2.6 stampy.py ..."


8. Q: Stampy complains about "mixed single-and paired-end in input",
      but I'm using only paired-end reads
   A: Stampy likely got confused by Casava 1.8 fastq headers.  Try
      using the --casava8 option.


9. Q: Can I use "-" (stdin) as input?
   A: Not quite, but you can build a pipe in this way:
         mkfifo --mode=0666 pipe.sam
	 bwa sampe [arguments] > pipe.sam &
	 ./stampy.py [arguments] --inputformat=sam -M pipe.sam
	 rm pipe.sam
      Note: (i) it is necessary to specify the input format, as
      format sniffing relies on seeks which don't work on fifos; 
      (ii) if you use sam or bam files, the fifos must have a .sam/
      .bam extension; (iii) this also works with two FASTQ files, 
      or multiple sam/bam files; just create a fifo for each.
