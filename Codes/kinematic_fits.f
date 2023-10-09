C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE fit_interf_signal(P,DP,M,CHISQR,C,D,loopcount)
C bind(C, name="carray")
C Description:
C ------------
C
C Kinematic Fit for Phi --> KlKs --> pi+pi-pi0pi0 events
C
C Number of input parameters: 38
C * Curvature, cotg(theta), phi for each track   (3x2)
C * E, X, Y, Z, T for each photon                (5x4)
C * Neutral and phi vertex coordinates           (3x2)
C * Sqrt(s)                                      (1x1)
C * Phi momentum                                 ( 3 )
C
C Total Constraints: 10
C * Time of flight of photons (4)
C * 4-momentum conservation   (4)
C * M(pi+pi-) = M(ko)         (1)
C * M( gggg ) = M(ko)         (1)
C
C Input parameters:
C -----------------
C M      = Number of constraints
C N      = Number of input parameters
C P (N)  = Measured Parameters
C DP(N)  = Parameter Errors
C
C Output parameters:
C ------------------
C CHISQR  = Chi2 value
C P (N)   = Adjusted Parameters
C loopcount = number of iterations
C
C-----------------------------------------------------------------------
        IMPLICIT NONE
C========================= Include klspm00.inc =========================
      INTEGER NCLMIN
      PARAMETER( NCLMIN = 4 )
      REAL ECLMIN
      PARAMETER( ECLMIN = 20. )
      REAL vtxdist
      PARAMETER( vtxdist = 50. )
      INTEGER testcounter
      PARAMETER( testcounter = 100 )
      INTEGER MaxNumFitPar
      PARAMETER( MaxNumFitPar = 38 )
      INTEGER MaxNumConstr
      PARAMETER( MaxNumConstr = 10 )
      INTEGER MaxNumComega
      PARAMETER( MaxNumComega = 8 )
      INTEGER cutN
      PARAMETER( cutN = 6 )
      INTEGER errN
      PARAMETER( errN = 17 )
      CHARACTER*40 messageCut(testcounter),messageErr(testcounter)
      LOGICAL    McFlag,EventOK,ErrFlag,truthmc,truthregg,truthomegaa,tr
     1  uthsemii,truththreee,truthelsee                                 
      INTEGER    nmctruth0,ntruth,nmctruth(7),ntruthregg,ntruthsemii,ntr
     1  uththreee,ntruthomegaa,ntruthelsee                              
      INTEGER    selected,analised
      INTEGER    counter(testcounter),counterBcg(testcounter)
      INTEGER    ErrFlagCount(testcounter),ErrFlagCountBcg(testcounter)
      LOGICAL    cutchar,cutneu,onlytruemc,makecuts
      COMMON / flagcommon / McFlag,EventOK,ErrFlag,truthmc
      COMMON / truthcommon / nmctruth0,ntruth,nmctruth
      COMMON / counterscommon / selected,analised,
     &                          counter,counterBcg
      COMMON / errcommon / ErrFlagCount,ErrFlagCountBcg,
     &                     messageCut,messageErr
      COMMON / talk / cutchar,cutneu,onlytruemc,makecuts
      SAVE /flagcommon/, /truthcommon/, /counterscommon/
      SAVE /errcommon/, /talk/
C-----------------------------------------------------------------------
C
C External functions
C
        REAL VDOTN, VMOD
C
C Local declarations
C
        REAL XI2
        COMMON/PHIFIT/XI2(50)
        INTEGER I, J, K, L, N
        PARAMETER( N = 36 )                  
C Number of constraints in fit
        INTEGER M                        
C Parameters & errors
        REAL*8    P(N), DP(N), P0(N), P1(N) 
C Parameter Corrections
        REAL*8  A(N)                      
C Derivative of constraint M w.r.t. par. N
        REAL*8  D(N,M)                    
C Balance of constraints
        REAL*8  C(M)                      
C Covariance Matrix (?)
        REAL*8  MTX(N+M,N+M)              
C (0,-C)
        REAL*8  Z(N+M)                    
        REAL*8  C1(M)
        INTEGER ITERATIONMAX, IFAIL, loopcount
        REAL*8    CHISQR, DX2, DX2MAX, R(100), C2
        REAL*8    PTEMP(N)
        LOGICAL FailFlag
C LAPACK
        INTEGER  IPIV(N+M), INFO
        REAL*8   WORK(N+M)

C------------------------------------------------------------------------------
C in:
C out:
        FailFlag = .FALSE.
C       ITERATIONMAX = 50
        ITERATIONMAX = 100
        CHISQR = 0.
        DX2 = 1.0e6
        DX2MAX = 0.01
        DO I = 1,N
          A(I)  = 0.
          P0(I) = P(I)
          Z(I)  = 0.
          DO J = 1,M
            D(I,J) = 0.d0
          ENDDO
        ENDDO

        WRITE(*,*) CHISQR
C
        CALL CONSTRAINT(M,N,P,C)
C
        loopcount = 0
        DO WHILE( loopcount.LT.ITERATIONMAX .AND. DX2.GT.DX2MAX )
          loopcount = loopcount+1
          DO I = 1,N
            P1(I) = P(I)
          ENDDO
         CALL DERIVATIVE(N,M,D,C,P,DP,PTEMP,C1)
          DO J = 1,M
            Z(N+J) = -C(J)
          ENDDO
          DO I = 1,N+M
            DO J = I,N+M
              IF( I.GT.N )THEN
                MTX(I,J) = 0.
              ELSEIF( J.GT.N )THEN
                MTX(I,J) = D(I,J-N)
              ELSEIF( I.EQ.J .AND. DP(I).GT.0 )THEN
                MTX(I,J) = 1/(DP(I)*DP(I))
              ELSE
                MTX(I,J) = 0.
              ENDIF
              IF( J.GT.I )THEN
                MTX(J,I) = MTX(I,J)
              ENDIF
            ENDDO
          ENDDO

          CALL DGETRF(N+M,N+M,MTX,N+M,IPIV,INFO)
          CALL DGETRI(N+M,MTX,N+M,IPIV,WORK,N+M,IFAIL)

          IF(IFAIL.NE.0) FailFlag = .TRUE.
          DX2 = CHISQR
          CHISQR = 0.
          DO I = 1,N
            A(I) = 0.
            DO J = 1,N+M
              A(I) = A(I) + MTX(I,J)*Z(J)
            ENDDO
C            P(I)= P1(I) + A(I)
            IF( ABS(A(I)).LT.1.0E3 )THEN
              P(I) = P1(I) + A(I)
            ENDIF
CC--- To avoid negative energy photons... ------------------
C!            IF( MOD(I+4,5).EQ.0 )THEN
C!               P(I) = ABS(P(I))
C!            ENDIF
            IF( I.EQ.7.or.I.EQ.12.or.I.EQ.17.or.I.EQ.22 )THEN
              IF( P(I).LT. ECLMIN ) THEN
                P(I) = ECLMIN
              ENDIF
            ENDIF
CC----------------------------------------------------------
            IF( DP(I).GT.0 )THEN
              CHISQR = CHISQR + ((P(I)-P0(I))/DP(I))**2
            ENDIF
          ENDDO
          DX2 = ABS(CHISQR-DX2)
          CALL CONSTRAINT(M,N,P,C)
          C2 = 0.
          DO J = 1,M
            C2 = C2 + C(J)*C(J)
          ENDDO
          XI2(loopcount) = CHISQR
        ENDDO
C
        IF( FailFlag ) loopcount = 100 + loopcount
C
        RETURN
        END
C
C
C=======================================================================
C=======================================================================
C=======================================================================
        SUBROUTINE DERIVATIVE(N,M,D,C,P,DP,P1,C1)
C
        IMPLICIT NONE
C
C Local declarations
C
        INTEGER I,J,K,L,M,N
        REAL*8    P(N), DP(N), P1(N)
        REAL*8  D(N,M), C(M), C1(M), D1, STEP
C------------------------------------------------------------------------------
        DO I=1,N
          P1(I)=P(I)
        ENDDO
        DO J=1,M
          DO I=1,N
            IF(DP(I).EQ.0)THEN
              IF(P1(I).EQ.0)THEN
                STEP=0.1
              ELSE
                STEP=0.001*P1(I)
              ENDIF
            ELSE
              STEP=0.01*DP(I)
            ENDIF
            L=0
            DO WHILE(ABS(D1-D(I,J)).GT.0.001*ABS(D(I,J)).OR.L.LE.1)
              D1=D(I,J)
              P1(I)=P(I)+STEP
              CALL CONSTRAINT(M,N,P1,C1)
              D(I,J)=(C1(J)-C(J))/STEP
              STEP=STEP/2
              L=L+1
            ENDDO
            P1(I)=P(I)
          ENDDO
        ENDDO
C
        RETURN
        END
C
C==============================================================================
C==============================================================================
C==============================================================================
C==============================================================================
        SUBROUTINE CONSTRAINT(M,N,P,C)
C
        IMPLICIT NONE
C-----------------------------------------------------------------------
C========================= Include klspm00.inc =========================
      INTEGER NCLMIN
      PARAMETER( NCLMIN = 4 )
      REAL ECLMIN
      PARAMETER( ECLMIN = 20. )
      REAL vtxdist
      PARAMETER( vtxdist = 50. )
      INTEGER testcounter
      PARAMETER( testcounter = 100 )
      INTEGER MaxNumFitPar
      PARAMETER( MaxNumFitPar = 38 )
      INTEGER MaxNumConstr
      PARAMETER( MaxNumConstr = 10 )
      INTEGER MaxNumComega
      PARAMETER( MaxNumComega = 8 )
      INTEGER cutN
      PARAMETER( cutN = 6 )
      INTEGER errN
      PARAMETER( errN = 17 )
      CHARACTER*40 messageCut(testcounter),messageErr(testcounter)
      LOGICAL    McFlag,EventOK,ErrFlag,truthmc,truthregg,truthomegaa,tr
     1  uthsemii,truththreee,truthelsee                                 
      INTEGER    nmctruth0,ntruth,nmctruth(7),ntruthregg,ntruthsemii,ntr
     1  uththreee,ntruthomegaa,ntruthelsee                              
      INTEGER    selected,analised
      INTEGER    counter(testcounter),counterBcg(testcounter)
      INTEGER    ErrFlagCount(testcounter),ErrFlagCountBcg(testcounter)
      LOGICAL    cutchar,cutneu,onlytruemc,makecuts
      COMMON / flagcommon / McFlag,EventOK,ErrFlag,truthmc
      COMMON / truthcommon / nmctruth0,ntruth,nmctruth
      COMMON / counterscommon / selected,analised,
     &                          counter,counterBcg
      COMMON / errcommon / ErrFlagCount,ErrFlagCountBcg,
     &                     messageCut,messageErr
      COMMON / talk / cutchar,cutneu,onlytruemc,makecuts
      SAVE /flagcommon/, /truthcommon/, /counterscommon/
      SAVE /errcommon/, /talk/
C======================== Include constans.inc =========================
      REAL       Mpip
      PARAMETER( Mpip = 139.57 )             
      REAL       Mpio
      PARAMETER( Mpio = 134.98 )             
      REAL       Mko
      PARAMETER( Mko  = 497.61 )             
      REAL       Momega
      PARAMETER( Momega  = 782.65 )          
      REAL       Cvel
      PARAMETER( Cvel = 29.9792458 )         
      REAL       EMCvel
      PARAMETER( EMCvel = 28.17 )            
      REAL       TauKs
      PARAMETER( TauKs = 0.08953 )   
      REAL       Mphi
      PARAMETER( Mphi  = 1019.460 )          
C= Include /kloe/soft/off/offline/inc/development/tls/maxstructdim.cin =
      INTEGER     NeleCluMax
      PARAMETER ( NeleCluMax = 2000 )
      INTEGER     MaxNumClu
      PARAMETER ( MaxNumClu = 100 )
      INTEGER    MaxNumVtx
      PARAMETER (MaxNumVtx = 20)
      INTEGER    MaxNumTrkV
      PARAMETER (MaxNumTrkV = 30)
      INTEGER    MaxNumTrk
      PARAMETER (MaxNumTrk = 100)
      INTEGER    MaxNumDHSP
      PARAMETER (MaxNumDHSP = 1000)
      Integer    nMaxDC
      Parameter (nMaxDC=1500)
      INTEGER    NQihiMax
      PARAMETER (NQihiMax=1000)
      INTEGER    NQcalMax
      PARAMETER (NQcalMax= 32 )
      INTEGER    MaxNumFirstHit
      PARAMETER (MaxNumFirstHit = 300 )
      INTEGER    MaxNtrkGen
      PARAMETER (MaxNtrkGen =50)
      INTEGER    MaxNvtxGen
      PARAMETER (MaxNvTxGen =30)
      integer TriggerElements
      parameter (TriggerElements = 300 )
C-----------------------------------------------------------------------
C
C External functions
C
C
C Local declarations
C
        INTEGER    M, N, ILoop, OffSet
        REAL*8       P(N), Ptrk(2,4), Pgam(4,4), Tgam(4)
        REAL*8       Ltot, PxTot, PyTot, PzTot, EnTot
        REAL*8       Mgggg, Mpipi, Vkne, Tkne
        REAL*8       Mkne,Vgam
             REAL*8       SinThetaBoost, CosThetaBoost
        REAL*8       PxPhi, PyPhi, PzPhi, EnPhi
        REAL*8     C(M)
C------------------------------------------------------------------------------
        Vgam = Cvel
        Mkne = Mko
C
C------------------------------------------------------------------------------
C Track related parameters
C------------------------------------------------------------------------------
C
C Px, Py, Py and E of the two pion from track parameters
C
        DO ILoop = 1,2
           Ptrk(ILoop,1) = (1000./ABS(P(3*(ILoop-1)+1)))*COS(P(3*(ILoop-
     &1)+2));
           Ptrk(ILoop,2) = (1000./ABS(P(3*(ILoop-1)+1)))*SIN(P(3*(ILoop-
     &1)+2));
           Ptrk(ILoop,3) = (1000./ABS(P(3*(ILoop-1)+1)))*P(3*(ILoop-1)+3
     &);
           Ptrk(ILoop,4) = SQRT( Ptrk(ILoop,1)**2 + Ptrk(ILoop,2)**2 +
     &          Ptrk(ILoop,3)**2 + Mpip**2 )
        ENDDO
C
C Invariant mass of pi+pi-
C
        PxTot = Ptrk(1,1) + Ptrk(2,1)
        PyTot = Ptrk(1,2) + Ptrk(2,2)
        PzTot = Ptrk(1,3) + Ptrk(2,3)
        EnTot = Ptrk(1,4) + Ptrk(2,4)
C
        Mpipi = SQRT( EnTot**2 - (PxTot**2+PyTot**2+PzTot**2) )
C
C------------------------------------------------------------------------------
C Photon related parameters
C------------------------------------------------------------------------------
C
C Px, Py, Pz, E and Time-of-Flight from cluster variables
C
        DO ILoop = 1,4
           OffSet = 5*(ILoop-1) + 6
C Cluster-NeuVtx distance
           Ltot = SQRT( (P(2+OffSet)-P(27))**2 +
     &          (P(3+OffSet)-P(28))**2 + (P(4+OffSet)-P(29))**2 )
           Pgam(ILoop,1) = P(1+OffSet) * (P(2+OffSet)-P(27))/Ltot
           Pgam(ILoop,2) = P(1+OffSet) * (P(3+OffSet)-P(28))/Ltot
           Pgam(ILoop,3) = P(1+OffSet) * (P(4+OffSet)-P(29))/Ltot
           Pgam(ILoop,4) = P(1+OffSet)
           Tgam(ILoop)   = Ltot / Vgam
        ENDDO
C
C Invariant mass of 4 photons
C
        PxTot = Pgam(1,1) + Pgam(2,1) + Pgam(3,1) + Pgam(4,1)
        PyTot = Pgam(1,2) + Pgam(2,2) + Pgam(3,2) + Pgam(4,2)
        PzTot = Pgam(1,3) + Pgam(2,3) + Pgam(3,3) + Pgam(4,3)
        EnTot = Pgam(1,4) + Pgam(2,4) + Pgam(3,4) + Pgam(4,4)
C
        Mgggg = SQRT( EnTot**2 - (PxTot**2+PyTot**2+PzTot**2) )
C
C Time-of-flight of the kaon decaying into photons
C
C NeuVtx-PhiVtx distance
        Ltot = SQRT( (P(27)-P(30))**2 +
     &       (P(28)-P(31))**2 + (P(29)-P(32))**2 )
        Vkne = Vgam * SQRT(PxTot**2+PyTot**2+PzTot**2) / EnTot
        Tkne = Ltot / Vkne
C Centroid-Apex distance - PROBABLY NEEDED, TO CHECK
C        Lcorr = SQRT( (P(27)-P(30))**2 +
C     &       (P(28)-P(31))**2 + (P(29)-P(32))**2 )
C        Vcorr = EMCvel
C        Tcorr = Ltot / Vcorr
C
C------------------------------------------------------------------------------
C Constraints
C------------------------------------------------------------------------------
C
C Total ToF for photons
C
        DO ILoop = 1,4
           OffSet = 5*(ILoop-1) + 6
           C(ILoop) = Tkne + Tgam(ILoop) - P(OffSet+5)
        ENDDO
C
C 4-momentum conservation
C
        PxTot = Ptrk(1,1) + Ptrk(2,1) + Pgam(1,1) + Pgam(2,1) +
     &       Pgam(3,1) + Pgam(4,1)
        PyTot = Ptrk(1,2) + Ptrk(2,2) + Pgam(1,2) + Pgam(2,2) +
     &       Pgam(3,2) + Pgam(4,2)
        PzTot = Ptrk(1,3) + Ptrk(2,3) + Pgam(1,3) + Pgam(2,3) +
     &       Pgam(3,3) + Pgam(4,3)
        EnTot = Ptrk(1,4) + Ptrk(2,4) + Pgam(1,4) + Pgam(2,4) +
     &       Pgam(3,4) + Pgam(4,4)
C
C        SinThetaBoost = P(38) / (P(36)+P(37))
C        CosThetaBoost = SQRT(1-SinThetaBoost**2)
C
C == SinThetaBoost * (P(36)+P(37))
        PxPhi = P(34)      
        PyPhi = P(35)
C eryk ERYK: in MC  PzPhi = 0
        PzPhi = P(36)
        EnPhi = P(33)
C
        C(5) = PxTot - PxPhi
        C(6) = PyTot - PyPhi
        C(7) = PzTot - PzPhi
        C(8) = EnTot - EnPhi
C
C Invariant masses
C
        C(9)  = Mgggg-Mkne
        C(10) = Mpipi-Mkne
C
        RETURN
        END
