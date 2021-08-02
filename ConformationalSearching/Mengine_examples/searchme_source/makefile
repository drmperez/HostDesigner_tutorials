FC     = gfortran

FFLAGS = -ffixed-form -ffixed-line-length-0 -fcray-pointer -fbounds-check -ff2c -O3 -Wall

# if code becomes too large, may need to add this to FFLAGS: -mcmodel=medium

SRC1 = searchme.f srchbonds.f readpcm.f nummol.f gettoken.f\
       writesrchcon.f

# Main searchme
OBJ1 = searchme.o srchbonds.o readpcm.o nummol.o gettoken.o\
       writesrchcon.o
       
OBJS = $(OBJ1) 

hd4: $(OBJS)
	$(FC) $(FFLAGS) $(FOPT) $(LIBS) $(OBJS) -o searchme 

.f.o:
	$(FC) $(FFLAGS) $(FOPT) -c $*.f

clean:
	rm -f *.o searchme
