FFLAGS=
FC=gfortran
LINKLIBS=

SRCDIR=src
OBJDIR=obj
EXEDIR=bin
DATADIR=data
INDEXDIR=indices
SRCDIR90=src_f90
OBJDIR90=obj90
EXEDIR90=bin90

SRCFILES=iritest.for irisub.for irifun.for iritec.for\
		 iridreg.for igrf.for cira.for iriflip.for
OBJECTS=$(patsubst %.for, $(OBJDIR)/%.o, $(SRCFILES))

SRCFILES90= cira.1.f90 elec_dens.1.f90 epstein_lay.1.f90\
			igrf.1.f90 ion_dens.0.f90 iridreg.1.f90\
			iriflip.1.f90 irifun.0.f90 irisub.1.f90 iritec.1.f90\
			iritest.1.f90 peaks.0.f90 profile.0.f90\
			temperature.0.f90 time.0.f90 vdrift.1.f90
OBJECTS90=$(patsubst %.f90, $(OBJDIR90)/%.o, $(SRCFILES90))

all: $(OBJDIR)   $(EXEDIR)   $(EXEDIR)/iri \
	 $(OBJDIR90) $(EXEDIR90) $(EXEDIR90)/iri90

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

#actually build the executable
$(EXEDIR90)/iri90: $(OBJECTS90)
	$(FC) $^ -o $@
	cp $(DATADIR)/* $(EXEDIR90)
	cp $(INDEXDIR)/* $(EXEDIR90)

#Instructions to build arbitrary objects into the OBJDIR90
$(OBJDIR90)/%.o: $(SRCDIR90)/%.f90
	$(FC) -c $^ -o $@

$(OBJDIR90):
	mkdir -p $@

$(EXEDIR90):
	mkdir -p $@

clean:
	rm -rf $(OBJDIR)
	rm -rf $(EXEDIR)
	rm -rf $(OBJDIR90)
	rm -rf $(EXEDIR90)
