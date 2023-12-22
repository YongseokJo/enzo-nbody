      MODULE POINTERS

      REAL*8, POINTER :: X(:,:),XDOT(:,:),X0(:,:), 
     & X0DOT(:,:), F(:,:),FDOT(:,:), 
     & BODY(:),RS(:),FI(:,:),D1(:,:), 
     & D2(:,:),D3(:,:),FR(:,:),D1R(:,:), 
     & D2R(:,:),D3R(:,:),STEP(:),T0(:), 
     & STEPR(:),T0R(:),TIMENW(:), 
     & RADIUS(:),TEV(:),TEV0(:), 
     & BODY0(:),EPOCH(:),SPIN(:),XSTAR(:), 
     & ZLMSTY(:),FIDOT(:,:),D0(:,:), 
     & FRDOT(:,:),D0R(:,:),KSTAR(:),FENZO(:,:)
 
      INTEGER, POINTER :: NAME(:),IMINR(:),LIST(:,:)
 
      REAL*8, POINTER :: U(:,:),U0(:,:),UDOT(:,:),FU(:,:), 
     & FUDOT(:,:),FUDOT2(:,:),FUDOT3(:,:), 
     & H(:),HDOT(:),HDOT2(:),HDOT3(:), 
     & HDOT4(:),DTAU(:),TDOT2(:),TDOT3(:), 
     & R(:),R0(:),GAMMA(:),SF(:),H0(:), 
     & FP0(:,:),FD0(:,:),KBLIST(:), 
     & KSLOW(:),TBLIST
 
      END MODULE POINTERS
