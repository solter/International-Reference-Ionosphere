OBJDIR=obj
FC=f77
LIBS=m
FFLAGS=
FLIBS=$(patsubst %,-l%, $(LIBS))

# Define that objects are in the OBJDIR
VPATH=%.o $(OBJDIR)

SOURCES=iritest.for irisub.for irifun.for iritec.for\
		iridreg.for igrf.for cira.for iriflip.for

OBJECTS=$(patsubst %.for, $(OBJDIR)/%.o, $(SOURCES))

all: iri.x

iri.x: $(OBJECTS)
	$(FC) $(FFLAGS) $(FLIBS) -o $@ $^

$(OBJDIR)/%.o: %.for | $(OBJDIR)
	$(FC) -c $< -o $@

$(OBJDIR):
	mkdir -p $@

clean:
	rm -f *.o
	rm -f iri.x
