#!/usr/bin/env python

import sys, getopt, re, os, traceback, math, random

if not sys.version.startswith("2.") or int(sys.version[2]) not in [7]:
   sys.stderr.write("Stampy requires Python version 2.7\n")
   sys.exit(1)

if sys.maxint == 2147483647:
   print >> sys.stderr,"Stampy requires a 64-bit Python install to run; 32-bit installations are not supported"
   sys.exit(1)

try:
   import maptools
except ImportError:
   print >> sys.stderr,"Could not import the maptools module -- does the maptools.so library exist?"
   raise

try:
    from Stampy import hashtable, genome, reader, formatter, mapper, singleendedmapper, \
                       recalibrate, mapstats, fastqformatter, utils, multithread, filez
except ImportError:
    print >> sys.stderr,"Could not import main python modules -- please check that the Python interpreter matches the installed Stampy version"

import plugins

from ext import makeconstants


short_opts =  "?!g:G:d:v:H:h:f:i:s:R:M:T:S:P:A:o:f:t:"
string_opts = ["genome=","logfile=","hash=","outputformat=","output=","inputformat=","stats=","keep-recalibrated=",
               "map=","comment=","recalibrate=","qualitybase=","recaldatasuffix=","testmap=","recaldatadir=","remote=",
               "readgroup=","processpart=","parse=","species=","assembly=","referenceuri=","pg=","labelfilter=","separator1=",
               "separator2=","bamsortprefix=","bamsortmemory=","adapter-strip="]
multi_opts = []
float_opts =  ["recalfraction=","substitutionrate="]
int_opts =    ["verbosity=","bits=","maxcount=","maxscore=","minposterior=","numrecords=","lowqthreshold=","seed=",
               "insertsize=","insertsd=","maxfingerprintvariants=","linearalignmentband=","simulate-minindellen=",
               "simulate-maxindellen=","tryvariants=","fastaqual=","simulate-numsubstitutions=","gapopen=","gapextend=",
               "recalscoreprefix=","svprior=","longindelprior=","baseentropy=","banding=","xa-max=","xa-max-discordant=",
               "insertsize2=","insertsd2=","padding=","maxpairseeds=","paircandlikethres=","maxbasequal=","threads=",
               "bamsortbuflen=","unclipq="]
bool_opts =   ["logit","sensitive","fast","simulate-duplications","solexa","solexaold","sanger","noautosense","norefoutput",
               "stripns","keeplabel","rightmost","alignquals","keepreforder","alphabetical","noparseNCBI","mdtag","strip-adapter",
               "ignore-improper-pairs","baq","overwrite","casava8","gatkcigarworkaround","bamkeepgoodreads","xcigar","cigarx"]
long_opts =   ["help","build-genome=","build-hash=","recalibrate=","map=","simulate="] + string_opts + int_opts + float_opts + bool_opts

oldhelp = """
 -R FILE, --recalibrate=FILE[,FILE] Compute recalibration data from fastq file(s) (not usually recommended)
 --recalfraction=F                  Fraction of reads to use for recalibration [0.01]
 --recaldatadir=DIR                 Directory for recalibration data [same as fastq file]
 --recalscoreprefix=N               Read prefix length over which likelihood is calculated for recalibration [20]
"""


help = """%s\n\nUsage: %s [options] [files]

Command options (one required)
 -G PREFIX, --build-genome=PREFIX   Build genome index PREFIX.stidx from fasta file(s) on command line
 -H PREFIX, --build-hash=PREFIX     Build hash PREFIX.sthash
 -M FILE, --map=FILE[,FILE]         Map fastq/fasta/BAM file(s)
 -S FILE, --simulate=FILE           Simulate reads following empirical qualities from read FILE
 -P FILE, --parse=FILE              Parse simulated SAM file, and report statistics
 -A FILE                            Convert qualities; strip adapters

Mapping options:
 -g PREFIX, --genome=PREFIX         Use genome index file PREFIX.stidx
 -h PREFIX, --hash=PREFIX           Use hash file PREFIX.sthash
 -t N, --threads=N                  Number of threads to use [1]
 --numrecords=N                     Number of records to process [all]
 --processpart=A/B                  Process the Ath read (pair) in every block of B read (pairs) [1/1]
 --minposterior=N                   Minimum read mapping phred posterior [no filtering]
 --maxscore=N                       Maximum read likelihood phred score [no filtering]
 --qualitybase=C                    Character or ASCII for quality score 0 in input [33 or !]
 --logit                            Qualities are logit (log p/1-p) rather than phred (log p) scores [phred]
 --solexa, --solexaold, --sanger    Short for: --qualitybase=@; --qualitybase=@ --logit; default
 --substitutionrate=F               Set substitution rate for mapping and simulation [0.001]
 --gapopen=N                        Gap open penalty (phred score) [40]
 --gapextend=N                      Gap extension penalty (phred score) [3]
 --insertsize=N                     (Initial) mean insert size for paired-end reads [250]
 --insertsd=N                       (Initial) standard deviation for insert size [60]
 --insertsize2=N                    (Initial) mean insert size for mate pairs [-2000]
 --insertsd2=N                      (Initial) standard deviation for mate pairs [-1 = deactivated]
 --maxpairseeds=N                   Number of ambiguous single mappings to take forward for paired realignment [25]
 --noautosense                      Do not auto-sense insert size distribution(s) [False]
 --sensitive                        More sensitive; about 25-50%%%% slower
 --fast                             Less sensitive; about 50%%%% (paired-end) to 100%%%% (single-ended) faster
%s
Output options:
 -o FILE, --output=FILE             Write mapping output to FILE [stdout]
 -f FMT, --outputformat=FMT         Select mapping output format (sam,maqtxt,maqmap,maqmapN) [sam]
 --rightmost                        Flush gaps to rightmost position in reference coordinates [leftmost]
 --xa-max                           Number of alternative hits to output for single/concordant paired reads [0, max 20]
 --xa-max-discordant                Number of alternative hits to output for discordant paired reads [0]
 --overwrite                        Allow overwriting of output files
 --baq                              (SAM format) Compute base-alignment quality (BAQ; BQ tag)
 --alignquals                       (SAM format) Compute posterior alignment probabilities (YQ tag)
 --xcigar                           (SAM format) Compute the extended CIGAR string (XC tag)
 --readgroup=ID:id,tag:value,...    (SAM format) Set readgroup tags (ID,SM,LB,DS,PU,PI,CN,DT,PL); ",," quotes ","
 --norefoutput                      (SAM format) Use = signs at read positions that match the reference
 --comment=comment1%%%%comment2%%%%...    (SAM format) Add @COmments.  %%%% separate lines; commas are converted into tabs
 --comment=@F                       (SAM format) Add @COmments from file F

Simulation options (use with -S):
 --seed=N                           Set random number seed [1]
 --simulate-minindellen=N           Minimum indel length to simulate; negative is deletion [0]
 --simulate-maxindellen=N           Maximum indel length to simulate; negative is deletion [0]
 --simulate-duplications            For insertions, simulate duplications rather than (default) random sequence insertion
 --simulate-numsubstitutions=N      Insert N substitutions in read (pairs) (default = Poisson distribution)
%s
Advanced options:
 --maxcount=N                       (-H) Maximum multiplicity of repeat words [200]
 --assembly=S                       (-G) Set assembly identifier (appears in @SQ field (AS tag) in SAM output)
 --species=S                        (-G) Set species identifier (appears in @SQ field (SP tag) in SAM output)
 --referenceuri=S                   (-G) Set URI for the reference (appears in @SQ field (UR tag) in SAM output)
 --noparseNCBI                      (-G) Don't parse but just copy gi|nnn|ref|xxx| labels in NCBI fasta files
 --strip-adapter                    (-A) Strip adapters from fastq files
 --inputformat=FMT                  Read input format (fasta, fastq, sam, bam) [fastq]
 --banding=N                        Band size for banded alignment; -1=none [60]
 --fastaqual=N                      Base quality for fasta input [30]
 --maxbasequal=N                    Maximum accepted base quality [50]
 --svprior=N                        Prior probability of read pair bridging a SV (phred score) [55]
 --longindelprior=N                 Prior probability of read pair bridging a long indel (phred score) [40]
 --baseentropy=N                    Entropy of inserted nucleotides (phred score) [5]
 --remote=U                         Command to get input reads from URL %%%%s; shortcuts are 'wget' and 'scp'
 --keepreforder                     (-M) (default) Use reference order in .fa/.stidx file, rather than alphabetical
 --alphabetical                     (-M) Order reference identifiers alphabetically
 --labelfilter=R                    (-M) Only map reads / read pairs whose label match the regular expression R
 --casava8                          (-M) Parse Casava v1.8 FastQ headers; include but do not map filtered reads
 --separator1=c                     (-M) Use 'c' as separator for Casava data columns [:]
 --separator2=c                     (-M) use 'c' as separator between first 7 columns, and the auxiliary columns [ ]
 --gatkcigarworkaround              (-M) Remove adjacent I/D CIGAR operators (valid in SAM spec, but trips up GATK)
 --unclipq                          (-M X.bam / -M --bamoptions) Undo soft clipping over low Q bases (-1=disable) [10]
 --bamsortprefix=S                  (-M X.bam) Temporary storage for sorting [--output argument, or /tmp/bamsort]
 --bamsortbuflen=N                  (-M X.bam) Entries in memory sort buffer [10000]
 --bamsortmemory=N                  (-M X.bam) Value for option -m to samtools sort [200000000] (samtools 1.3: 200M)
 --bamkeepgoodreads                 (-M X.bam) If set, do not re-map already well-mapped reads from BAM [False]
%s
General options:
 -?, --help                         Display this help page
 -v N, --verbosity=N                Set verbosity level to N [2]
 --logfile=FILE                     Send logging information to FILE  [stderr]
%s
"""


shorthelp = """%s\n\nUsage: %s [options] [.fa files]


Option summary (--help for all):

Command options
 -G PREFIX file1.fa [...]           Build genome index PREFIX.stidx from fasta file(s) on command line
 -H PREFIX                          Build hash PREFIX.sthash
 -M FILE[,FILE]                     Map fastq/fasta/BAM file(s)
 -A FILE                            Convert qualities; strip adapters

Mapping/output options
 -g PREFIX                          Use genome index file PREFIX.stidx
 -h PREFIX                          Use hash file PREFIX.sthash
 -o FILE                            Write mapping output to FILE [stdout]
 --readgroup=ID:id,tag:value,...    Set read-group tags (ID,SM,LB,DS,PU,PI,CN,DT,PL)  (SAM format)
 --solexa, --solexaold, --sanger    Solexa read qualities (@-based); pre-v1.3 Solexa; and Sanger (!-based, default)
 --substitutionrate=F               Set substitution rate for mapping and simulation [0.001]
 --gapopen=N                        Gap open penalty (phred score) [40]
 --gapextend=N                      Gap extension penalty (phred score) [3]
%s
General options
 --help                             Full help
 -v N                               Set verbosity level (0-3) [2]
%s
"""


default_options = """
--svprior=55
--longindelprior=40
--baseentropy=5
--logfile=stderr
--verbosity=2
--bits=-1
--maxcount=200
--outputformat=sam
--inputformat=none
--output=-
--stats=mapstats.cache
--qualitybase=!
--numrecords=-1
--recalscoreprefix=20
--recalfraction=0.01
--recaldatasuffix=.recaldata
--maxscore=99999
--minposterior=-99999
--lowqthreshold=10
--seed=1
--insertsize=250
--insertsd=80
--insertsize2=-2000
--insertsd2=-1
--maxfingerprintvariants=3
--linearalignmentband=3
--simulate-minindellen=0
--simulate-maxindellen=0
--simulate-numsubstitutions=0
--substitutionrate=0.001
--tryvariants=-1
--gapopen=40
--gapextend=3
--fastaqual=30
--maxbasequal=50
--banding=60
--xa-max=0
--xa-max-discordant=0
--padding=160
--maxpairseeds=25
--paircandlikethres=100
--separator1=:
--threads=1
--bamsortbuflen=10000
--bamsortmemory=200000000
--unclipq=10
%s
"""

outputformats = ["sam","maqtxt","maqmap","maqmapShort","maqmapShortN"]
inputformats = ["fastq","fasta","sam","bam","none"]


def install_plugins():
   global help, shorthelp
   global default_options
   global short_opts, string_opts, float_opts, int_opts, bool_opts, long_opts, multi_opts
   
   mapping_opts=plugins.hook_add("get_mapping_opts","")
   simulation_opts = plugins.hook_add("get_simulation_opts","")
   advanced_opts = plugins.hook_add("get_advanced_opts","")
   general_opts = plugins.hook_add("get_general_opts","")
   default_opts = plugins.hook_add("get_default_opts","")
   new_short_opts = plugins.hook_add("get_new_short_opts","")
   new_string_opts = plugins.hook_add("get_new_string_opts",[])
   new_multi_opts =  plugins.hook_add("get_new_multi_opts",[])
   new_float_opts = plugins.hook_add("get_new_float_opts",[])
   new_int_opts = plugins.hook_add("get_new_int_opts",[])
   new_bool_opts = plugins.hook_add("get_new_bool_opts",[])
   new_long_opts = plugins.hook_add("get_newlong_opts",[])
   
   help = help % ("%s", "%s", mapping_opts, simulation_opts, advanced_opts, general_opts)
   shorthelp = shorthelp % ("%s", "%s", mapping_opts, general_opts)
   default_options += default_opts

   short_opts += new_short_opts
   string_opts += new_string_opts
   multi_opts += new_multi_opts
   float_opts += new_float_opts
   int_opts += new_int_opts
   bool_opts += new_bool_opts
   long_opts += new_long_opts + new_string_opts + new_int_opts + new_float_opts + new_bool_opts


def validate( settings ):
        if settings.verbosity < 0 or settings.verbosity > 3:
            raise ValueError("--verbosity out of range (must be between 0 and 3)")
        if settings.bits not in range(16,30) + [-1]:
            raise ValueError("--bits out of range (must be between 16 and 29)")
        if settings.maxcount < 1:
            raise ValueError("--maxcount must be >= 1")
        if settings.outputformat not in outputformats:
            raise ValueError("--outputformat not one of " + str(outputformats))
        if settings.inputformat not in inputformats:
            raise ValueError("--inputformat not one of " + str(inputformats))
        if len(settings.qualitybase)>1:
            raise ValueError("--qualitybase must be a single character or integer")
        if settings.insertsize < 1:
            raise ValueError("--insertsize must be a positive integer")
        if settings.insertsd < 1:
            raise ValueError("--insertsd must be a positive integer")
        if not -1 <= settings.maxfingerprintvariants <= 4:
            raise ValueError("--maxfingerprintvariants must be an integer between 1 and 4, or -1 for no filtering")
        if not 0 <= settings.substitutionrate <= 0.5:
            raise ValueError("--substitutionrate must be between 0 and 0.5")
        if not -1 <= settings.tryvariants <= 3:
            raise ValueError("--tryvariants must be -1 (autoselect), 0 (none), 1 (all) or 2, 3 (half/third of locations)")
        if settings.keepreforder and settings.alphabetical:
            raise ValueError("Must give at most one of (--keepreforder, --alphabetical)")
        if len(settings.separator1) != 1:
            raise ValueError("Separator1 must be a single character")
        if len(settings.separator2) != 1:
            raise ValueError("Separator2 must be a single character")
        if settings.threads < 1:
           raise ValueError("--threads must be a positive integer")
        if not settings.bamsortprefix:
           settings.bamsortprefix = [settings.output, "/tmp/bamsort"][settings.output == "-"]
        settings.keepreforder = not settings.alphabetical
        if settings.readgroup:
            tags = ("ID","SM","LB","DS","PU","PI","CN","DT","PL")
            tagvalues = [ v.replace('#$&*',',').split(':') for v in settings.readgroup.replace(',,','#$&*').split(',') ]
            for tv in tagvalues:
                if len(tv) < 2 or len(tv[0]) != 2: 
                    raise ValueError("--readgroup: I don't understand '%s' (expect two-letter-TAG:VALUE pair)" % ':'.join(tv))
                tv = [tv[0], ":".join(tv[1:])]
                if tv[0] not in tags: 
                    raise ValueError("--readgroup: tag '%s' unknown, must be one of %s" % (tv[0],','.join(tags)))
                try: tv[0] == "PI" and int(tv[1])
                except ValueError: raise ValueError("--readgroup PI tag must be numeric (found: '%s')" % tv[1])
            if len( set(["ID"]).intersection( set( tv[0] for tv in tagvalues ) ) ) != 1:
                # removed the check for SM tag, since --readgroup can be used to select a readgroup from a BAM file
                raise ValueError("--readgroup must set at least ID tag")
            settings.readgroupid = [ tv[1] for tv in tagvalues if tv[0] == "ID" ][0]
            try:    settings.readgroupsm = [ tv[1] for tv in tagvalues if tv[0] == "SM" ][0]
            except: settings.readgroupsm = "(none)"
        if settings.comment:
            if settings.comment[0] == "@":
               try:
                  f = open(settings.comment[1:],'r')
                  lines = [ line.strip() for line in f if len(line)>1 ]
               except:
                  raise ValueError("--comment: Could not open file '%s' for reading" % settings.comment[1:])
               settings.comment = lines
            else:
               settings.comment = [ c.replace(",","\t") for c in settings.comment.split('%') ]
        if settings.processpart:
            try:
                pq = map(int, settings.processpart.split('/'))
                if len(pq) != 2 or not (1 <= pq[0] <= pq[1]): raise ValueError("")
            except:
                raise ValueError("--processpart: argument should be of the form 'p/q' with 1 <= p <= q, q > 0")
            settings.processpartpq = pq


def parse_options( options, settings, fail_duplicates = False ):

    actions = []

    for option, value in options:

        if option == "-v": option = "--verbosity"
        if option == "-d": option = "--defaults"
        if option == "-g": option = "--genome"
        if option == "-h": option = "--hash"
        if option == "-f": option = "--outputformat"
        if option == "-i": option = "--inputformat"
        if option == "-s": option = "--stats"
        if option == "-o": option = "--output"
        if option == "-t": option = "--threads"

        if (len(value)>0 and value[0] == "-" and option[min(2,len(option)):]+"=" not in int_opts
            and not option.startswith('--')):
            sys.stderr.write("Warning -- argument to option '%s' starts with '-'\n" % option)

        optkey = option[2:].replace('-','_')

        if option == "-?" or option == "-!" or option == "--help":
            if option == "-!": helpstr = shorthelp
            else:              helpstr = help
            print helpstr % (settings.program + " v" + settings.version + ", <%s>" % settings.author, sys.argv[0])
            sys.exit(0)

        elif option == "-G" or option == "--build-genome":
            actions.append( ("-G", value) )

        elif option == "-H" or option == "--build-hash":
            actions.append( ("-H", value) )

        elif option == "-R" or option == "--recalibrate":
            actions.append( ("-R", value) )

        elif option == "-M" or option == "--map":
            actions.append( ("-M", value) )

        elif option == "-T" or option == "--testmap":
            actions.append( ("-T", value) )

	elif option == "-S" or option == "--simulate":
            actions.append( ("-S", value) )

        elif option == "-P" or option == "--parse":
            actions.append( ("-P", value) )

        elif option == "-A":
            actions.append( ("-A", value) )

        elif option in ["--solexa","--solexaold","--sanger"]:
            if option[:8] == "--solexa": settings.qualitybase = "@"
            if option == "--solexaold": settings.logit = True

        elif option[2:] in bool_opts:
            settings.__dict__[optkey] = True

        elif option[2:]+"=" in string_opts:
            if settings.__dict__[optkey] == None:
                settings.__dict__[optkey] = value
            elif value and option[2:]+"=" in multi_opts:
               settings.__dict__[optkey] = settings.__dict__[optkey] + " " + value

        elif option[2:]+"=" in float_opts:
            if settings.__dict__[optkey] == None:
                try:
                    settings.__dict__[optkey] = float(value)
                except:
                    raise ValueError("Option %s requires a float argment, got '%s'" % (option, value))

        elif option[2:]+"=" in int_opts:
            if settings.__dict__[optkey] == None:
                try:
                    settings.__dict__[optkey] = int(value)
                except:
                    raise ValueError("Option %s requires an integer argument, got '%s'" % (option, value))
    
        else:
            raise ValueError("Internal error: unrecognized option " + option)

    return actions


def build_genome( settings, logger, genomefile, arguments ):
    if len(arguments) == 0:
        raise ValueError("Please provide .fa[.gz] files to process")
    genomefile += ".stidx"
    if not settings.overwrite:
       utils.test_write( genomefile )
    logger.log("Building genome...",1)
    logger.log("Input files: " + str( [a.split(':')[-1] for a in arguments] ), 2)
    # chop off .gz/.dz extension, if it exists
    for i in range(len(arguments)):
        if arguments[i].lower().endswith(".gz"): arguments[i] = arguments[i][:-3]
    g = genome.genome( genomefile, parseNCBI = not settings.noparseNCBI, write=True, padding = settings.padding )
    g.openFromFa( arguments, assembly=settings.assembly, species=settings.species, uri=settings.referenceuri )


def build_hash( settings, logger, hashfile, arguments ):
    if len(arguments)>0: raise ValueError("Unexpected arguments on command line")
    hashfile += ".sthash"
    utils.open_genome( settings, logger )
    if not settings.overwrite:
       utils.test_write( hashfile )
    hashtable.buildhash3( hashfile, settings, logger )


def stripadapter( settings, logger, readfiles, arguments ):
    input_filenames = readfiles.split(',') + arguments
    if len(input_filenames) != 1:
       raise ValueError("Please provide one .fa[.gz] file to process")
    input_format = reader.sniff( input_filenames[0], logger )
    if input_format != "fastq":
       raise ValueError("Input file %s does not appear to be a fastq file" % input_filenames[0] )
    if settings.output == "-": 
       output = sys.stdout
    else: 
       if not settings.overwrite:
          utils.test_write(settings.output)
       output = open(settings.output,'w')
    logger.log("Converting fastq file...",1)
    maptools.stripadapter( input_filenames[0], output, settings, logger, filez )


def recalibratereads( settings, logger, readfiles, arguments ):
    input_filenames = readfiles.split(',') + arguments
    if len(input_filenames) == 1: input_format = reader.sniff( input_filenames[0], logger )
    else:                         input_format = "fastq"
    recaldata_filenames = [ utils.make_recaldatafilename(settings,name) for name in input_filenames ]
    if not settings.overwrite:
       for f in recaldata_filenames: 
          if not settings.overwrite:
             utils.test_write( f, filetype="Recalibration data file" )
    utils.open_genome( settings, logger )

    if input_format == "fastq":
        utils.open_hash( settings, logger )
        datareaders = [ reader.make( f, logger, settings, process_fraction=settings.recalfraction ) for f in input_filenames ]
        mapinstance = mapper.mapper( settings.hashmap, settings.genomemap, settings,
                                     lowqthreshold=settings.lowqthreshold,
                                     maxfingerprintvariants = settings.maxfingerprintvariants,
                                     phred_gapopen = settings.gapopen, phred_gapextend = settings.gapextend,
                                     tryvariants = settings.tryvariants,
                                     baseentropy = settings.baseentropy)
        if len(input_filenames) == 1:
	    semapper = singleendedmapper.SingleEndMapper( settings )
	    outputgenerator = semapper.generator( datareaders[0], mapinstance, genome )
        else:
	    pairedreader = reader.PairedReader( datareaders[0], datareaders[1], logger )
	    pemapper = singleendedmapper.PairedEndMapper( settings )
	    outputgenerator = pemapper.generator( pairedreader, mapinstance, genome )
    else:
        outputgenerator = reader.SamReader( input_filenames[0], settings, process_fraction=settings.recalfraction,
                                            logger=logger ).generator()
        
    logger.log("Recalibrating...",1)
    recalibrate.recalibrator( recaldata_filenames, outputgenerator, settings, maxreadlength=1000, maxq=50, pseudocount = 1 )


def mapreads_background( settings, logger, queue, pipe, process_idx ):
    """ maps reads from the queue and puts results in the pipe """

    logger.setprocess( process_idx )
    mapinstance = mapper.mapper( settings.hashmap, settings.genomemap, settings,
                                 lowqthreshold=settings.lowqthreshold,
                                 maxfingerprintvariants = settings.maxfingerprintvariants,
                                 phred_gapopen = settings.gapopen, phred_gapextend = settings.gapextend,
                                 tryvariants = settings.tryvariants,
                                 baseentropy = settings.baseentropy)
    readmapper = singleendedmapper.PairedEndMapper( settings )
    datareader = reader.QueueReader( queue, readmapper, settings )
    outputgenerator = readmapper.generator( datareader, mapinstance, settings.genomemap )
    for data in outputgenerator:
       update = readmapper.get_pair_stats_updates()
       pipe.send( (data, update) )
    pipe.send( None )
    logger.log("Child exiting",3)


def mapreads( settings, logger, readfiles, arguments ):
    input_filenames = readfiles.split(',') + arguments
    bams = sum( ("1234"+name.upper())[-4:] in [".BAM",".SAM"] for name in input_filenames ) == len(input_filenames)
    if not bams and len(input_filenames) not in [1,2]: raise ValueError("Map: please provide one or two input files")
    utils.open_genome( settings, logger )
    utils.open_hash( settings, logger )
    recaldata_filenames = [ utils.make_recaldatafilename(settings,name) for name in input_filenames ]
    mapinstance = mapper.mapper( settings.hashmap, settings.genomemap, settings,
                                 lowqthreshold=settings.lowqthreshold,
                                 maxfingerprintvariants = settings.maxfingerprintvariants,
                                 phred_gapopen = settings.gapopen, phred_gapextend = settings.gapextend,
                                 tryvariants = settings.tryvariants,
                                 baseentropy = settings.baseentropy)
    formatinstance = formatter.make( settings, logger )
    if len(input_filenames) == 1 or bams:
        # single-ended mapping
        if os.path.exists(recaldata_filenames[0]):
            recalhandle = open(recaldata_filenames[0],'r')
            recaldata = recalibrate.get_single_read_calibration( recalibrate.read_calibration_file( recalhandle ), 0 )
            recalhandle.close()
        else:
            recaldata = [None]
        if bams:
           datareader = reader.make( input_filenames, logger, settings, recal_data = None, dep_data = [None] )
        else:
           datareader = reader.make( input_filenames[0], logger, settings, recal_data = recaldata[0], dep_data = recaldata )
    else:
        # paired-end mapping
        recaldata = [[None], [None]]
        for idx in [0,1]:
            if os.path.exists(recaldata_filenames[idx]):
                recalhandle = open(recaldata_filenames[idx],'r')
                recaldata[idx] = recalibrate.get_single_read_calibration( recalibrate.read_calibration_file( recalhandle ), idx )
                recalhandle.close()
        datareaders = ( reader.make( input_filenames[0], logger, settings, 
				     recal_data = recaldata[0][0], dep_data=recaldata[0] ),
                        reader.make( input_filenames[1], logger, settings, 
				     recal_data = recaldata[1][0], dep_data=recaldata[1] ) )
        datareader = reader.PairedReader( datareaders[0], datareaders[1], logger )
    readmapper = singleendedmapper.PairedEndMapper( settings )
    readmapper = multithread.MultiMapper( settings, readmapper, logger )
    outputgenerator = plugins.hook( 'bwa', settings, logger, readmapper, datareader, mapinstance, singleendedmapper )
    if outputgenerator == None: outputgenerator = readmapper.generator( datareader, mapinstance, settings.genomemap )
    logger.log("Mapping...",1)
    formatgenerator = formatinstance.formatter( outputgenerator )
    for output in formatgenerator: pass
    formatinstance.dumpbuffer()
    if len(input_filenames) == 2:
       readmapper.dump_pair_stats( logger )


def testmap( settings, logger, readfiles, arguments ):
    input_filenames = readfiles.split(',') + arguments
    if len(input_filenames) not in [1,2]: raise ValueError("TestMap: please provide one or two input files")
    utils.open_genome( settings, logger )
    utils.open_hash( settings, logger )
    recaldata_filenames = [ utils.make_recaldatafilename(settings,name) for name in input_filenames ]
    mapinstance = mapper.mapper( settings.hashmap, settings.genomemap, settings,
                                 lowqthreshold=settings.lowqthreshold,
                                 maxfingerprintvariants = settings.maxfingerprintvariants,
                                 phred_gapopen = settings.gapopen, phred_gapextend = settings.gapextend,
                                 tryvariants = settings.tryvariants,
                                 baseentropy =settings.baseentropy)
    mapinstance.make_scoretable()
    formatinstance_output = formatter.make( settings, logger )
    formatinstance = formatter.TestDataCollector( formatinstance_output )
    recaldata = [[None] for f in recaldata_filenames]
    for idx in range(len(recaldata)):
        if os.path.exists(recaldata_filenames[idx]):
            recalhandle = open(recaldata_filenames[idx],'r')
            recaldata[idx] = recalibrate.get_single_read_calibration( recalibrate.read_calibration_file( recalhandle ), idx )
            recalhandle.close()
    logger.log("Testing...",1)
    if len(input_filenames) == 1:
        # single-ended mapping
        datareader = reader.make( input_filenames[0], logger, settings, 
				      recal_data = recaldata[0][0], dep_data = recaldata[0] )
        randomreader = reader.RandomReader( datareader, settings )
        themapper = singleendedmapper.SingleEndMapper( settings )
    else:
        # paired-end mapping
        datareaders = ( reader.make( input_filenames[0], logger, settings, 
				     recal_data = recaldata[0][0], dep_data = recaldata[0] ),
                        reader.make( input_filenames[1], logger, settings, 
				     recal_data = recaldata[1][0], dep_data = recaldata[1] ) )
        pairedreader = reader.PairedReader( datareaders[0], datareaders[1], logger )
        randomreader = reader.RandomReader( pairedreader, settings )
        themapper = singleendedmapper.PairedEndMapper( settings )

    outputgenerator = themapper.generator( randomreader, mapinstance, settings.genomemap )
    formatgenerator = formatinstance.formatter( outputgenerator )
    for output in formatgenerator: pass
    formatinstance.dumpbuffer()
    formatinstance.report( sys.stderr )


def simulatedata( settings, logger, readfiles, arguments ):
    input_filenames = readfiles.split(',') + arguments
    if len(input_filenames) not in [1,2]: raise ValueError("SimulateData: please provide one or two input files")
    recaldata_filenames = [ utils.make_recaldatafilename(settings,name) for name in input_filenames ]
    recaldata = [[None] for f in recaldata_filenames]
    for idx in range(len(recaldata)):
        if os.path.exists(recaldata_filenames[idx]):
            recalhandle = open(recaldata_filenames[idx],'r')
	    logger.log("Using recalibration data from %s" % recaldata_filenames[idx], 3)
	else:
            #This is a hack.  A better solution is to read all files in the recaldatadir, read the readgroup
            #from the SAM file, select the right one (allow multiple, or enforce a single readgroup per sam?)
            if not settings.recaldatadir: settings.recaldatadir = ''
            recalfile = os.path.join( settings.recaldatadir, os.path.basename( recaldata_filenames[idx] ) )
	    if os.path.exists( recalfile ):
	        recalhandle = open(recalfile,'r')
		logger.log("Using recalibration data from %s" % recalfile, 3)
	    else:
                recalhandle = None
		logger.log("Not using recalibration data", 3)
	if recalhandle != None:
            recaldata[idx] = recalibrate.get_single_read_calibration( recalibrate.read_calibration_file( recalhandle ), idx )
            recalhandle.close()
    utils.open_genome( settings, logger )
    if len(input_filenames) == 1:
        # single-ended mapping
        datareader = reader.make( input_filenames[0], logger, settings, 
				  recal_data = None, dep_data = recaldata[0] )
        randomreader = reader.RandomReader( datareader, settings )
    else:
        # paired-end mapping
        datareaders = ( reader.make( input_filenames[0], logger, settings, 
				     recal_data = None, dep_data = recaldata[0] ),
                        reader.make( input_filenames[1], logger, settings, 
				     recal_data = None, dep_data = recaldata[1] ) )
        pairedreader = reader.PairedReader( datareaders[0], datareaders[1], logger )
        randomreader = reader.RandomReader( pairedreader, settings )

    theformatter = fastqformatter.FastqFormatter( input_filenames, settings )
    logger.log("Simulating...",1)
    theformatter.format( randomreader )


def parsetest( settings, logger, readfiles, arguments ):
    input_filenames = readfiles.split(',') + arguments
    if len(input_filenames) != 1: raise ValueError("Parse: please provide one input file")
    utils.open_genome( settings, logger )
    formatinstance_output = formatter.make( settings, logger )
    formatinstance = formatter.TestDataCollector( formatinstance_output )
    logger.log("Parsing...",1)
    samreader = reader.SamReader( input_filenames[0], settings, logger=logger )
    formatgenerator = formatinstance.formatter( samreader.generator() )
    for output in formatgenerator: pass
    formatinstance.dumpbuffer()
    formatinstance.report( sys.stderr )






##################################################################################
#
# top-level code
#
##################################################################################



def main():

    global logger, settings

    options, arguments = getopt.getopt( sys.argv[1:], short_opts, long_opts )
    actions = parse_options( options, settings, fail_duplicates = True )
    defaultopts, defaultargs = [], []

    for line in default_options.split("\n"):
        o, a = getopt.gnu_getopt( re.split("\s*",line), short_opts, long_opts )
        defaultopts += o

    defaultactions = parse_options( defaultopts, settings )
    if not settings.separator2:
       settings.separator2 = " "

    if len(settings.qualitybase)>1:
        try:    settings.qualitybase = chr( int( settings.qualitybase ) )
        except: pass     # caught by validate()

    if settings.sensitive and settings.fast:
        raise ValueError("Can only set one of --sensitive and --fast")

    if settings.sensitive: 
        settings.linearalignmentband = 4       # no effect when SSE available
        settings.maxfingerprintvariants = -1
        settings.tryvariants = 1
        settings.maxpairseeds = max(settings.maxpairseeds, 45)
    if settings.fast: 
        settings.paircandlikethres = 50
        settings.maxfingerprintvariants = 2
        settings.linearalignmentband = 2       # no effect when SSE available
        settings.tryvariants = 3

    if len(actions) == 0:
        if len(options) == 0: parse_options([("-!","")], settings)
        else:
           if len(arguments)>0: 
              print "Nothing to do, with unused arguments on command line"
              print "First unused argument: '%s'" % arguments[0]
           else: print "Nothing to do...\n"
    elif len(actions) > 1:
        raise ValueError("Please provide exactly one of -G, -H, -M, ...")
    else:
        settings.cmdline = ' '.join( sys.argv[1:] )
        p1, p2 = settings.stats, os.path.join( os.path.dirname(os.path.abspath(sys.argv[0])), settings.stats )
        if os.path.exists(p1): statspath = p1
        elif os.path.exists(p2): statspath = p2
        else: statspath = p1
        settings.mapstats = mapstats.MapStats( statspath )

    logger.setfilename( settings.logfile )
    logger.log( "Starting Stampy with the following options:\n" + str( settings ), 3 )

    validate( settings )
    random.seed( settings.seed )
    actiondict = dict( actions )

    while len( actiondict ) > 0:

        if "-G" in actiondict:
            build_genome( settings, logger, actiondict['-G'], arguments )
            del actiondict['-G']

        elif "-H" in actiondict:
            logger.log("Building hash table...",1)
            build_hash( settings, logger, actiondict['-H'], arguments )
            del actiondict['-H']

        elif "-A" in actiondict:
            stripadapter( settings, logger, actiondict['-A'], arguments )
            del actiondict['-A']

        elif "-R" in actiondict:
            recalibratereads( settings, logger, actiondict['-R'], arguments )
            del actiondict['-R']

        elif "-M" in actiondict:
            mapreads( settings, logger, actiondict['-M'], arguments )
            del actiondict['-M']

        elif "-T" in actiondict:
            testmap( settings, logger, actiondict['-T'], arguments )
            del actiondict['-T']
	    
	elif "-S" in actiondict:
	     simulatedata( settings, logger, actiondict['-S'], arguments )
	     del actiondict['-S']

        elif "-P" in actiondict:
            parsetest( settings, logger, actiondict['-P'], arguments )
            del actiondict['-P']

        else:
            raise ValueError("Internal error: unrecognized actions: " + str(actiondict))

    if len(actions)>0: logger.log("Done",1)




if __name__ == '__main__':
   try:
       install_plugins()
       settings = utils.Settings( string_opts, int_opts, float_opts, bool_opts )
       logger = utils.Logger( settings )
       makeconstants.bind_all( sys.modules['Stampy.formatter'], verbose=False )
       makeconstants.bind_all( sys.modules['Stampy.singleendedmapper'], verbose=False )
       makeconstants.bind_all( sys.modules['Stampy.genome'], verbose=False )
       makeconstants.bind_all( sys.modules['Stampy.mapstats'], verbose=False )
       main()
   except:
       errtype, value, tb = sys.exc_info()
       if errtype == SystemExit:
           raise
       if (errtype != KeyboardInterrupt) or True:
           tbstr = "Traceback:\n" + ''.join( traceback.format_tb( tb ))
           logger.log( tbstr[:-1], 3 )
           if settings.verbosity < 3:
              tbstr =  "Internal problem: " + (["(no info)"]+traceback.format_tb( tb ))[-1][:-1].replace('\n',' : ')
           else: tbstr = None
       if errtype in (ValueError,getopt.GetoptError,IOError):
           logger.log("\nError: "+str(value),0)
           if errtype == getopt.GetoptError: logger.log("Use -? for help\n",0)
       elif errtype == KeyboardInterrupt:
           logger.log("\nInterrupted by user",0)
       else:
           if tbstr: logger.log(tbstr,0)
           logger.log("Internal problem: " + repr(value),0)
           logger.log("Please report this problem to <%s> mentioning the version number (%s) and command line - thanks!\n" % 
                      (settings.author,settings.version),0)
       sys.exit(1)

