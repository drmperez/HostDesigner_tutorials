FC     = gfortran
FFLAGS =  -ffixed-form -ffixed-line-length-0 -fcray-pointer -fbounds-check -ff2c  -Wall

SRC1 = makedata.f gettoken.f readpcm.f writepcm.f delatom.f\
       confanal.f srchbonds.f writesrchcon.f nummol.f getlow.f\
       makeaquo.f superimpose.f center.f jacobi.f quatfit.f\
       davies.f upper.f loadmet.f atomtyper.f

OBJ1 = makedata.o gettoken.o readpcm.o writepcm.o delatom.o\
       confanal.o srchbonds.o writesrchcon.o nummol.o getlow.o\
       makeaquo.o superimpose.o center.o jacobi.o quatfit.o\
       davies.o upper.o loadmet.o atomtyper.o

OBJS = $(OBJ1) 

makedata: $(OBJS)
	$(FC) $(FFLAGS) $(FOPT) $(LIBS) $(OBJS) -o makedata

.f.o:
	$(FC) $(FFLAGS) $(FOPT) -c $*.f

clean:
	rm -f *.o makedata
