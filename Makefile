#***************************************************************
#wilcoxon-command-line (c) A. Favorov 2015	
#$Id$
#***************************************************************

name=wilcoxon-command-line
exename=wilcoxon-command-line

cd=.
#current dir
wd=./probability
#wilcoxon and other distributions
od=./obj

srcdirlist=$(cd):$(wd)

empty=
space=$(empty) $(empty)
includeflags = $(foreach dir,$(subst :,$(space),$(srcdirlist)),$(INCLUDEKEY)$(dir)) $(INCLUDECLOSETERM)
#this strange invocation is just preparing -I flag from srcdirlist.

include ccvars

.PHONY: all clean

vpath %.c $(srcdirlist)
vpath %.cc $(srcdirlist)

all: $(od) $(exename)$(EXEEXT) 

$(od):
	mkdir -p $(od)

OBJS=$(od)/$(name).o \
$(od)/gauss.o \
$(od)/wilcoxon.o

$(od)/%.o: %.c
	$(CC) $(CCFLAGS) $< -o $@

$(od)/%.o: %.cc
	$(CPP) $(CPPFLAGS) $< -o $@

$(exename)$(EXEEXT): $(OBJS)
	$(CPP) -o $(exename)$(EXEEXT) $(OBJS) $(LINKFLAGS) 
	chmod 755 $(exename)$(EXEEXT) 

$(od)/wilcoxon.o $(od)/gauss.o : $(wd)/gauss.h
$(od)/$(name).o $(od)/wicoxon.o : $(wd)/wilcoxon.h

$(od)/$(name).o : $(cd)/$(name).cc
$(od)/gauss.o : $(wd)/gauss.h
$(od)/wicoxon.o : $(wd)/wilcoxon.c

clean:
	rm -f $(OBJS)
	rm -f `find . -name '*~'`
