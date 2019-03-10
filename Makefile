LATEX=lualatex
BIBER=biber
REPORT_NAME=cosc6326-pa1-michael-yantosca
INCS=-I/usr/include/openmpi
LIBS=-lmpi
MPICPP=mpic++
CFLAGS=-g -std=c++11

all: $(REPORT_NAME).pdf nkmax

nkmax: nkmax.o
	$(MPICPP) $(CFLAGS) $^ -o $@ $(LIBS)

%.o: %.cpp
	$(MPICPP) $(CFLAGS) $(INCS) -c $*.cpp

%.pdf: %.tex
	@$(LATEX) -shell-escape $*.tex
	@$(BIBER) $*
	@$(LATEX) -shell-escape $*.tex
#	@$(LATEX) -shell-escape $*.tex

superclean: clean superclean-doc-$(REPORT_NAME)

clean: clean-doc-$(REPORT_NAME)
	@rm -f *~
	@rm -f *.o
	@rm -f nkmax

superclean-doc-%:
	rm -f $*.pdf

clean-doc-%:
	@rm -f $*.aux
	@rm -f $*.bbl
	@rm -f $*.bcf
	@rm -f $*.log
	@rm -f $*.run.xml
	@rm -f $*.dvi
	@rm -f $*.blg
	@rm -f $*.auxlock
	@rm -f $*.pyg
	@rm -f $*-figure*
	@rm -f $*.toc
	@rm -f $*.out
	@rm -f $*.snm
	@rm -f $*.nav
