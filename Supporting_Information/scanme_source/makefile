FC     = gfortran

FFLAGS = -ffixed-form -ffixed-line-length-0 -fcray-pointer -fbounds-check -ff2c -O3 -Wall

# if code becomes too large, may need to add this to FFLAGS: -mcmodel=medium

SRC1 = searchme.f srchbonds.f readpcm.f nummol.f gettoken.f\
       writesrchcon.f scanme.f delatom.f writepcm.f

# Main scanme
OBJ1 = searchme.o srchbonds.o readpcm.o nummol.o gettoken.o\
       writesrchcon.o scanme.o delatom.o writepcm.o
       
OBJS = $(OBJ1) 

scanme: $(OBJS)
	$(FC) $(FFLAGS) $(FOPT) $(LIBS) $(OBJS) -o scanme 

.f.o:
	$(FC) $(FFLAGS) $(FOPT) -c $*.f

clean:
	rm -f *.o scanme
