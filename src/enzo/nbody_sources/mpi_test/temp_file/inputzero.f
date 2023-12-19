      SUBROUTINE INPUTZERO
*
*       version1 (2023/01/08) : version with  
*       putting input variables insed code
*       ---------------------
*
      INCLUDE 'common6.h'

*       parameters needed to run the code

* original input from file

*1 100000.0 1.E8 40 40 640
*10001 1 20 10000 100 1 1
*0.02 0.02 0.19839410 1.00 1.00 1000.00 1.0E10 2.34663365 10.09899010
*2 2 1 0 1 0 4 2 0 0
*0 0 0 2 2 1 1 0 0 0
*1 2 1 0 0 2 0 0 0 2
*0 0 0 0 1 0 1 1 0 1
*0 0 0 0 0 0 0 0 0 0
*1.0E-05 1.0E-04 0.1 1.0 5.0E-05 0.01 1.0
*0.0 1.0 1.0 0 0 0.02 0.0 100.0
*0.50 0.0 0.0 0.0 0.125
*9.56543900e+10 8.50000000


*       input needed for nbody6.F

      KSTART = 1
      TCOMP = 100000.0
      TCRITp = 1.E8
      isernb = 40
      iserreg = 40
      iserks = 640


*       input needed for input.F

      N = 10001
      NFIX = 1
      NCRIT = 20
      NRAND = 10000
      NNBOPT = 100
      NRUN = 1 
*       1 is left out... why? 

      ETAI = 0.02
      ETAR = 0.02
      RS0 = 0.19839410
      DTADJ = 1.00
      DELTAT = 1.00
      TCRIT = 1000.00
      QE = 1.0E10
      RMBAR = 2.34663365
      ZMBAR = 10.09899010

      KZ(1) = 2
      KZ(2) = 2
      KZ(3) = 1
      KZ(4) = 0
      KZ(5) = 1
      KZ(6) = 0
      KZ(7) = 4
      KZ(8) = 2
      KZ(9) = 0
      KZ(10) = 0

      KZ(11) = 0
      KZ(12) = 0
      KZ(13) = 0
      KZ(14) = 2
      KZ(15) = 2
      KZ(16) = 1
      KZ(17) = 1
      KZ(18) = 0
      KZ(19) = 0
      KZ(20) = 0

      KZ(21) = 1
      KZ(22) = 2
      KZ(23) = 1
      KZ(24) = 0
      KZ(25) = 0
      KZ(26) = 2
      KZ(27) = 0
      KZ(28) = 0
      KZ(29) = 0
      KZ(30) = 2

      KZ(31) = 0
      KZ(32) = 0
      KZ(33) = 0
      KZ(34) = 0
      KZ(35) = 1
      KZ(36) = 0
      KZ(37) = 1
      KZ(38) = 1
      KZ(39) = 0
      KZ(40) = 1

      KZ(41) = 0
      KZ(42) = 0
      KZ(43) = 0
      KZ(44) = 0
      KZ(45) = 0
      KZ(46) = 0
      KZ(47) = 0
      KZ(48) = 0
      KZ(49) = 0
      KZ(50) = 0

      DTMIN = 1.0E-05
      RMIN = 1.0E-04
      ETAU = 0.1
      ECLOSE = 1.0
      GMIN = 5.0E-05
      GMAX = 0.01
      SMAX = 1.0


*       input needed for data.F

      ALPHA = 0.0
      BODY1 = 1.0
      BODYN = 1.0
      NBIN0 = 0
      NHI0 = 0
      ZMET = 0.02
      EPOCH0 = 0.0
      DTPLOT = 100.0


*       input needed for setup.F -> not used currently

*      AP0 = 
*      ECC = 
*      N2 = 
*      SCALE = 


*       input needed for scale.F

      Q = 0.5
      VXROT = 0.0 
      VZROT = 0.0
      RTIDE = 0.0


*       input needed for xtrnl10.F (KZ14 = 2, standard solar input)

      GMG = 9.56543900e+10 
      RG0 = 8.50000000


*       input needed for binpop.F

*      SEMI = 
*      ECC = 
*      RATIO = 
*      RANGE = 
*      NSKIP = 
*      IDORM = 


*       input needed for hipop.F


*      SEMI =  
*      ECC = 
*      RATIO = 
*      RANGE = 
*      NSKIP = 
*      IDORM = 


*       input needed for imbhinit.F

*      MMBH = 
*      XBH(1) = 
*      XBH(2) = 
*      XBH(3) = 
*      VBH(1) = 
*      VBH(2) = 
*      VBH(3) = 
*      DTBH = 


*       input needed for cloud0.F


*      NCL = 
*      RB2 = 
*      VCL = 
*      SIGMA = 
*      CLM = 
*      RCL2 = 


