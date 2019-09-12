#
# - Python 2.7 with 2- and 4-byte Unicode representations (UCS2 and UCS4) are supported.
#   Compiles with gcc 5.4.0 and Apple LLVM version 6.1.0
#   

SHELL=/bin/bash
python?=python
platform=`uname`
bdir=build/

#############################
# User targets
#############################

all:
	@if [[ `$(python) --version 2>&1 |awk -F "." '{print $$2}'` -eq 7 ]]; then \
	  $(MAKE) ucs=`$(python) -c "import sys;print sys.maxunicode > 65536 and 'ucs4' or 'ucs2'"` build-stampy ; \
	else \
	  echo "Python version 2.7 is required.  (If it is installed as python2.7, try 'make python=python2.7')" ; \
	fi

clean:
	-rm maptools.so
	-rm build/*/*.pyc

#############################
# End of user targets
#############################

objs=$(bdir)/pyx/maptools.o $(bdir)/c/maputils.o $(bdir)/c/alignutils.o $(bdir)/readalign.o $(bdir)/algebras.o $(bdir)/frontend.o

$(bdir)/c/%.o:
	mkdir -p $(bdir)/c
	gcc -pthread -fno-strict-aliasing -DNDEBUG -g -fwrapv -O2 -Wall -fPIC `$(python)-config --includes` -c c/$*.c -o $@

$(bdir)/pyx/%.o:
	mkdir -p $(bdir)/pyx
	ln -s -f `pwd`/pyx/$(ucs)/$*.c pyx/$*.c
	gcc -pthread -fno-strict-aliasing -DNDEBUG -g -fwrapv -O2 -Wall -fPIC `$(python)-config --includes` -c pyx/$*.c -o $@

$(bdir)/%.o:
	gcc -pthread -fno-strict-aliasing -DNDEBUG -g -fwrapv -O2 -Wall -fPIC `$(python)-config --includes` -c ext/readalign/$*.cc -o $@

build-stampy: $(objs)
	cp -r build/python2.7/* .
	@if [ $(platform) = Linux ]; then \
	 g++ `$(python)-config --ldflags` -pthread -shared $(objs) -o maptools.so ; \
	else \
	 g++ `$(python)-config --ldflags` -pthread -dynamiclib $(objs) -o maptools.so ; \
	fi

