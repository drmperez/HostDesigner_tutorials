c***********************************************************************
c***********************************************************************
      subroutine davies(zlig,zmet,oldi,logk,newi,newk)
c***********************************************************************
c***********************************************************************
c  computes log K at a ionic strength of newi given the
c  charge of the ligand, charge of the metal, and old ionic 
c  strength of oldi, and old log K using Davies Equation. 
c***********************************************************************

      implicit none
      
      integer zlig     ! ligand charge
      integer zmet     ! metal species charge
      integer zcom     ! complex charge

      real    oldi     ! old ionic strength
      real    newi     ! new ionic strength
      real    logk     ! old log K
      real    sslogk   ! standard state log K
      real    A        ! davies parameter
      real    B        ! davies parameter
      real    gammet
      real    gamlig
      real    gamcom
      real    sqoldi
      real    sqnewi
      real    newk
      
      zcom = zmet - zlig
      
      A=0.51e0            !original Davies used a value of 0.500
      B=0.20e0            !original Davies used a value of 0.300
      sqoldi=sqrt(oldi)
      sqnewi=sqrt(newi)
      zcom=zlig+zmet

      gammet = -A*zmet*zmet*(sqoldi/(sqoldi+1)-B*oldi)
      gamlig = -A*zlig*zlig*(sqoldi/(sqoldi+1)-B*oldi)
      gamcom = -A*zcom*zcom*(sqoldi/(sqoldi+1)-B*oldi)

      sslogk = logk + gamcom - gammet - gamlig

      gammet = -A*zmet*zmet*(sqnewi/(sqnewi+1)-B*newi)
      gamlig = -A*zlig*zlig*(sqnewi/(sqnewi+1)-B*newi)
      gamcom = -A*zcom*zcom*(sqnewi/(sqnewi+1)-B*newi)

      newk = sslogk - gamcom + gammet + gamlig
      
      return
      end

