c***********************************************************************
c     'metals.i' - metal ion properties
c***********************************************************************

      integer metcharge,mettype,metcn
      real sir,rdhe,rdhc,metsolv
      character*2 metatlab


      common  /metals/ metcharge(MAXMTYPE),mettype(MAXMTYPE),
     &        metcn(MAXMTYPE),sir(MAXMTYPE),metsolv(MAXMTYPE),
     &        rdhe(MAXMTYPE),rdhc(MAXMTYPE),metatlab(MAXMTYPE)
         
