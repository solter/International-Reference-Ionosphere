FFLAGS=
FC=gfortran
LINKLIBS=

SRCDIR=src
OBJDIR=obj
EXEDIR=bin
DATADIR=data
INDEXDIR=indices

SRCFILES=iritest.for irisub.for irifun.for iritec.for\
		 iridreg.for igrf.for cira.for iriflip.for
OBJECTS=$(patsubst %.for, $(OBJDIR)/%.o, $(SRCFILES))

all: $(OBJDIR) $(EXEDIR) $(EXEDIR)/iri

#actually build the executable
$(EXEDIR)/iri: $(OBJECTS)
	$(FC) $^ -o $@
	cp $(DATADIR)/* $(EXEDIR)
	cp $(INDEXDIR)/* $(EXEDIR)

#Instructions to build arbitrary objects into the OBJDIR
$(OBJDIR)/%.o: $(SRCDIR)/%.for
	$(FC) -c $^ -o $@

$(OBJDIR):
	mkdir -p $@

$(EXEDIR):
	mkdir -p $@

clean:
	rm -rf $(OBJDIR)
	rm -rf $(EXEDIR)
