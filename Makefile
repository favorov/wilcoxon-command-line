#*****************************************************************
#  APSampler. Looking for complex genetic interaction patterns 
# by a Metropolis-Hastings MCMC project. (c) A. Favorov 1999-2010
#    $Id$
#*****************************************************************

#this file can be started only from parent directory
name=wilcoxon_command_line
exename=wilcoxon_command_line

cd=.
wd=./probability
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

all: $(exename)$(EXEEXT) 

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
