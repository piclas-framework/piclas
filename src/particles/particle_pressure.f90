#include "boltzplatz.h"

MODULE MOD_part_pressure
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE
   
   PRIVATE
   PUBLIC       :: ParticlePressureIni, ParticlePressure, ParticleInsideCheck,ParticlePressureCellIni,ParticlePressureRem

CONTAINS   
   
SUBROUTINE ParticlePressureIni()

  USE MOD_Particle_Vars
  USE MOD_Globals
  USE MOD_Mesh_Vars,      ONLY : nElems,ElemToSide,SideToElem, nSides, nInnerSides, nBCSides

  IMPLICIT NONE
 
  REAL                 :: BN(3), BV1(3), BV2(3), BV3(3), SV1(3), SV2(3), OV(3), epsi
  REAL                 :: det1, det2, det3, RandVal(2), dist1, dist2, dete(6,2)
  INTEGER, ALLOCATABLE :: TempElemTotalInside(:), TempElemPartlyInside(:)
  INTEGER              :: nNodesInside, nBoundNodes, nInterest, nInterOld, Element, Side, ExamElem
  INTEGER              :: iShot, iElem, iNode, iSpec, iInit, iSide, i, n
  LOGICAL              :: ElementFound, InElementCheck
  
  
  epsi = 100.*epsilon(epsi)
  DO iSpec = 1,nSpecies
    DO iInit = Species(iSpec)%StartnumberOfInits, Species(iSpec)%NumberOfInits
      Species(iSpec)%Init(iInit)%ConstPress%InitialTemp = Species(iSpec)%Init(iInit)%MWTemperatureIC
      !Species(iSpec)%Init(iInit)%ConstPress%OrthoVector(1) = 0
      IF ((Species(iSpec)%Init(iInit)%ParticleEmissionType .EQ.3).OR.(Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.5)) THEN
        ALLOCATE (TempElemTotalInside(nElems))
        ALLOCATE (TempElemPartlyInside(nElems))
        ALLOCATE (Species(iSpec)%Init(iInit)%ConstPress%ElemStat(nElems))
        
        Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside = 0
        Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside = 0
        Species(iSpec)%Init(iInit)%ConstPress%ElemStat(1:nElems) = 3
        SELECT CASE (TRIM(Species(iSpec)%Init(iInit)%SpaceIC))
        CASE ('cuboid')
          
          BV1 = Species(iSpec)%Init(iInit)%BaseVector1IC
          BV2 = Species(iSpec)%Init(iInit)%BaseVector2IC
          
          Species(iSpec)%Init(iInit)%ConstPress%OrthoVector(1) = BV1(2) * BV2(3) - BV1(3) * BV2(2)
          Species(iSpec)%Init(iInit)%ConstPress%OrthoVector(2) = BV1(3) * BV2(1) - BV1(1) * BV2(3)
          Species(iSpec)%Init(iInit)%ConstPress%OrthoVector(3) = BV1(1) * BV2(2) - BV1(2) * BV2(1)
          OV = Species(iSpec)%Init(iInit)%ConstPress%OrthoVector
          
          IF ((OV(1) .EQ. 0) .AND. (OV(2) .EQ. 0) .AND. (OV(3) .EQ. 0)) THEN
            WRITE(*,*) 'Error in InitializeVariables: Cannot calculate Normal Vector of InitVolume in EmissionCase(3 or 4)!'
            STOP
          END IF
          
          Species(iSpec)%Init(iInit)%ConstPress%OrthoVector &
               = OV*(Species(iSpec)%Init(iInit)%CuboidHeightIC/SQRT(OV(1)**2 + OV(2)**2 + OV(3)**2))
          OV   = Species(iSpec)%Init(iInit)%ConstPress%OrthoVector
          Species(iSpec)%Init(iInit)%ConstPress%Determinant = ABS(BV1(1)*BV2(2)*OV(3) + BV1(2)*BV2(3)*OV(1) + BV1(3)*BV2(1)*OV(2) &
               - BV1(3)*BV2(2)*OV(1) - BV1(1)*BV2(3)*OV(2) - BV1(2)*BV2(1)*OV(3))
          !WRITE (*,*) 'ConstPress', Species(iSpec)%ConstantPressure , 'Deter', Species(iSpec)%ConstPress%Determinant, 'ispec', &
          !iSpec &
          !,'MPF',Species(iSpec)%MacroParticleFactor
          Species(iSpec)%Init(iInit)%ParticleEmission &
               = Species(iSpec)%Init(iInit)%ConstantPressure * Species(iSpec)%Init(iInit)%ConstPress%Determinant / &
               (BoltzmannConst * Species(iSpec)%Init(iInit)%MWTemperatureIC * Species(iSpec)%MacroParticleFactor)
          Species(iSpec)%Init(iInit)%ConstPress%EkinInside &
               = 1.5*Species(iSpec)%Init(iInit)%ConstantPressure*Species(iSpec)%Init(iInit)%ConstPress%Determinant + &
               INT(Species(iSpec)%Init(iInit)%ParticleEmission)*0.5*Species(iSpec)%MassIC*Species(iSpec)%Init(iInit)%VeloIC**2* &
               Species(iSpec)%MacroParticleFactor
          
          WRITE (*,'(A49,I3.3,A2,g0)') 'Number of Particles inside ConstPressArea ',iInit, ': ', &
               Species(iSpec)%Init(iInit)%ParticleEmission
          
          DO iElem = 1,nElems
            nNodesInside = 0
            nBoundNodes = 0
            DO iNode = 1,8
              BN(1:3) = GEO%NodeCoords(:,GEO%ElemToNodeID(iNode,iElem)) - Species(iSpec)%Init(iInit)%BasePointIC
              det1 = BN(1)*BV2(2)*OV(3) + BN(2)*BV2(3)*OV(1) + BN(3)*BV2(1)*OV(2) - &
                   BN(3)*BV2(2)*OV(1) - BN(1)*BV2(3)*OV(2) - BN(2)*BV2(1)*OV(3)
              det2 = BV1(1)*BN(2)*OV(3) + BV1(2)*BN(3)*OV(1) + BV1(3)*BN(1)*OV(2) - &
                   BV1(3)*BN(2)*OV(1) - BV1(1)*BN(3)*OV(2) - BV1(2)*BN(1)*OV(3)
              det3 = BV1(1)*BV2(2)*BN(3) + BV1(2)*BV2(3)*BN(1) + BV1(3)*BV2(1)*BN(2) - &
                   BV1(3)*BV2(2)*BN(1) - BV1(1)*BV2(3)*BN(2) - BV1(2)*BV2(1)*BN(3)
              
              det1 = det1/Species(iSpec)%Init(iInit)%ConstPress%Determinant
              det2 = det2/Species(iSpec)%Init(iInit)%ConstPress%Determinant
              det3 = det3/Species(iSpec)%Init(iInit)%ConstPress%Determinant
              
              IF (((det1 .LE. 1. + epsi) .AND. (det1 .GE. 1.-epsi)) .OR. ((det1 .LE. epsi) .AND. (det1 .GE. -epsi))) THEN
                IF (((det2 .LE. 1.+epsi) .AND. (det2 .GE. -epsi)) .AND. ((det3 .LE. 1.+epsi) .AND. (det3 .GE. -epsi))) THEN
                  nBoundNodes = nBoundNodes + 1
                END IF
              ELSE IF (((det2 .LE. 1.+epsi).AND.(det2 .GE.1-epsi)) .OR. ((det2 .GE. -epsi).AND.(det2 .LE. epsi))) THEN
                IF ((det1 .LE. 1.+epsi) .AND. (det1 .GE. -epsi) .AND. (det3 .LE. 1.+epsi) .AND. (det3 .GE. -epsi)) THEN
                  nBoundNodes = nBoundNodes + 1
                END IF
              ELSE IF (((det3 .GE. 1.-epsi).AND.(det3 .LE.1+epsi)) .OR. ((det3 .GE.-epsi).AND.(det3.LE.epsi))) THEN
                IF ((det2 .LE. 1.+epsi) .AND. (det2 .GE. -epsi) .AND. (det1 .LE. 1.+epsi) .AND. (det1 .GE. -epsi)) THEN
                  nBoundNodes = nBoundNodes + 1
                END IF
              ELSE IF ((det1 .LT. 1.) .AND. (det1 .GT. 0.) .AND. (det2 .LT. 1.) .AND. (det2 .GT. 0.) &
                   .AND. (det3 .LT. 1.) .AND. (det3 .GT. 0.)) THEN
                nNodesInside = nNodesInside + 1
              END IF
            END DO
            
            IF (nNodesInside .EQ. 8) THEN
              Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside = Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside + 1
              TempElemTotalInside(Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside)=iElem
              Species(iSpec)%Init(iInit)%ConstPress%ElemStat(iELem) = 1
            ELSE IF ((nNodesInside .GE. 1) .AND. (nNodesInside .LE. 7)) THEN
              IF (nNodesInside + nBoundNodes .EQ. 8) THEN
                Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside = Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside + 1
                TempElemTotalInside(Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside)=iElem
                Species(iSpec)%Init(iInit)%ConstPress%ElemStat(iELem) = 1
              ELSE
                Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside +1
                TempElemPartlyInside(Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside) = iElem
                Species(iSpec)%Init(iInit)%ConstPress%ElemStat(iELem) = 2
              END IF
            ELSE
              IF (nBoundNodes .EQ. 8) THEN
                Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside = Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside + 1
                TempElemTotalInside(Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside)=iElem
                Species(iSpec)%Init(iInit)%ConstPress%ElemStat(iELem) = 1
              ELSE
                Species(iSpec)%Init(iInit)%ConstPress%ElemStat(iELem) = 3
              END IF
            END IF
          END DO
          
          ! If no Element has been found, search for Basepoint ========================================
          IF ((Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside .EQ. 0) .AND. &
               (Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside .EQ. 0)) THEN
            CALL PointInsideQuad3D(iSpec,iInit,Element,InElementCheck,dete)
            IF (InElementCheck) THEN
              Species(iSpec)%Init(iInit)%ConstPress%ElemStat(Element) = 2
              Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside = 1
              TempElemPartlyInside(1) = Element
            END IF
          END IF
          
          ! Shoot on Sides of Neighbour-Elements of totally inside Elements ===========================
          
          DO iElem = 1,Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside
            Element = TempElemTotalInside(iElem)
            DO iSide = 1,6
              Side = ElemToSide(E2S_SIDE_ID, iSide, Element)
              IF (SideToElem(S2E_ELEM_ID,Side) .EQ. Element) THEN
                ExamElem = SideToElem(S2E_NB_ELEM_ID, Side)
              ELSE IF (SideToElem(S2E_NB_ELEM_ID,Side) .EQ. Element) THEN
                ExamElem = SideToElem(S2E_ELEM_ID,Side)
              END IF
              IF (Species(iSpec)%Init(iInit)%ConstPress%ElemStat(ExamElem) .EQ. 3) THEN
                DO iShot = 1,200
                  IF (iShot .LE. 100) THEN
                    !Shooting (123)
                    SV1 = GEO%NodeCoords(:,GEO%ElemSideNodeID(2, iSide, Element)) - &
                         GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                  ELSE
                    !Shooting (134)
                    SV1 = GEO%NodeCoords(:,GEO%ElemSideNodeID(4, iSide, Element)) - &
                         GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                  END IF
                  SV2 = GEO%NodeCoords(:,GEO%ElemSideNodeID(3, iSide, Element)) - &
                       GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                  DO
                    CALL RANDOM_NUMBER(RandVal)
                    IF (RandVal(1) + RandVal(2) .LE. 1) EXIT
                  END DO
                  BN = GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element)) + RandVal(1)*SV1 + RandVal(2)*SV2
                  BN = BN - Species(iSpec)%Init(iInit)%BasePointIC
                  !Lokalisierung
                  det1 = BN(1)*BV2(2)*OV(3) + BN(2)*BV2(3)*OV(1) + &
                       BN(3)*BV2(1)*OV(2) - BN(3)*BV2(2)*OV(1) - &
                       BN(1)*BV2(3)*OV(2) - BN(2)*BV2(1)*OV(3)
                  det2 = BV1(1)*BN(2)*OV(3) + BV1(2)*BN(3)*OV(1) + &
                       BV1(3)*BN(1)*OV(2) - BV1(3)*BN(2)*OV(1) - &
                       BV1(1)*BN(3)*OV(2) - BV1(2)*BN(1)*OV(3)
                  det3 = BV1(1)*BV2(2)*BN(3) + BV1(2)*BV2(3)*BN(1) + &
                       BV1(3)*BV2(1)*BN(2) - BV1(3)*BV2(2)*BN(1) - &
                       BV1(1)*BV2(3)*BN(2) - BV1(2)*BV2(1)*BN(3)
                  
                  det1 = det1/Species(iSpec)%Init(iInit)%ConstPress%Determinant
                  det2 = det2/Species(iSpec)%Init(iInit)%ConstPress%Determinant
                  det3 = det3/Species(iSpec)%Init(iInit)%ConstPress%Determinant
                  
                  IF ((det1 .LT. 1.-epsi) .AND. (det1 .GT. epsi) .AND. (det2 .LT. 1-epsi) .AND. &
                       (det2 .GT. epsi) .AND. (det3 .LT. 1.-epsi) .AND. (det3 .GT. epsi)) THEN
                    Species(iSpec)%Init(iInit)%ConstPress%ElemStat(ExamElem) = 2
                    Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside &
                         = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside + 1
                    TempElemPartlyInside(Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside) = ExamElem
                    EXIT
                  END IF
                END DO
              END IF
            END DO
          END DO

#ifdef MPI

          DO Side=1, nSides
            IF (Side.GT.nInnerSides+nBCSides) THEN
              IF (SideToElem(S2E_ELEM_ID,Side).NE.-1) THEN
                Element = SideToElem(S2E_ELEM_ID,Side)
                iSide = SideToElem(S2E_LOC_SIDE_ID,Side)
              ELSE
                Element = SideToElem(S2E_NB_ELEM_ID,Side)
                iSide = SideToElem(S2E_NB_LOC_SIDE_ID,Side)
              END IF
              IF (Species(iSpec)%Init(iInit)%ConstPress%ElemStat(Element) .EQ. 3) THEN
                DO iShot = 1,200
                  IF (iShot .LE. 100) THEN
                    !Shooting (123)
                    SV1 = GEO%NodeCoords(:,GEO%ElemSideNodeID(2, iSide, Element)) - &
                         GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                  ELSE
                    !Shooting (134)
                    SV1 = GEO%NodeCoords(:,GEO%ElemSideNodeID(4, iSide, Element)) - &
                         GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                  END IF
                  SV2 = GEO%NodeCoords(:,GEO%ElemSideNodeID(3, iSide, Element)) - &
                       GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                  DO
                    CALL RANDOM_NUMBER(RandVal)
                    IF (RandVal(1) + RandVal(2) .LE. 1) EXIT
                  END DO
                  BN = GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element)) + RandVal(1)*SV1 + RandVal(2)*SV2
                  BN = BN - Species(iSpec)%Init(iInit)%BasePointIC
                  !Lokalisierung
                  det1 = BN(1)*BV2(2)*OV(3) + BN(2)*BV2(3)*OV(1) + &
                       BN(3)*BV2(1)*OV(2) - BN(3)*BV2(2)*OV(1) - &
                       BN(1)*BV2(3)*OV(2) - BN(2)*BV2(1)*OV(3)
                  det2 = BV1(1)*BN(2)*OV(3) + BV1(2)*BN(3)*OV(1) + &
                       BV1(3)*BN(1)*OV(2) - BV1(3)*BN(2)*OV(1) - &
                       BV1(1)*BN(3)*OV(2) - BV1(2)*BN(1)*OV(3)
                  det3 = BV1(1)*BV2(2)*BN(3) + BV1(2)*BV2(3)*BN(1) + &
                       BV1(3)*BV2(1)*BN(2) - BV1(3)*BV2(2)*BN(1) - &
                       BV1(1)*BV2(3)*BN(2) - BV1(2)*BV2(1)*BN(3)
                  
                  det1 = det1/Species(iSpec)%Init(iInit)%ConstPress%Determinant
                  det2 = det2/Species(iSpec)%Init(iInit)%ConstPress%Determinant
                  det3 = det3/Species(iSpec)%Init(iInit)%ConstPress%Determinant
                  
                  IF ((det1 .LT. 1.-epsi) .AND. (det1 .GT. epsi) .AND. (det2 .LT. 1-epsi) .AND. &
                       (det2 .GT. epsi) .AND. (det3 .LT. 1.-epsi) .AND. (det3 .GT. epsi)) THEN
                    Species(iSpec)%Init(iInit)%ConstPress%ElemStat(Element) = 2
                    Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside &
                         = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside + 1
                    TempElemPartlyInside(Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside) = Element
                    EXIT
                  END IF
                END DO
              END IF
            END IF
          END DO
#endif

          ! Shoot on Sides of Neighbour-Elements of partly inside Elements ==========================================
          nInterest = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside
          nInterOld = 1
          DO WHILE (nInterest .GE. nInterOld)
            DO iElem = nInterOld, nInterest
              
              Element = TempElemPartlyInside(iElem)
              DO iSide = 1, 6
                Side = ElemToSide(E2S_SIDE_ID, iSide, Element)
                IF (SideToElem(S2E_ELEM_ID,Side) .EQ. Element) THEN
                  ExamElem = SideToElem(S2E_NB_ELEM_ID, Side)
                ELSE IF (SideToElem(S2E_NB_ELEM_ID,Side) .EQ. Element) THEN
                  ExamElem = SideToElem(S2E_ELEM_ID,Side)
                END IF
                
                IF (Species(iSpec)%Init(iInit)%ConstPress%ElemStat(ExamElem) .EQ. 3) THEN
                  DO iShot = 1,200
                    IF (iShot .LE. 100) THEN
                      !Shooting (123)
                      SV1 = GEO%NodeCoords(:,GEO%ElemSideNodeID(2, iSide, Element)) - &
                           GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                    ELSE
                      !Shooting (134)
                      SV1 = GEO%NodeCoords(:,GEO%ElemSideNodeID(4, iSide, Element)) - &
                           GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                    END IF
                    SV2 = GEO%NodeCoords(:,GEO%ElemSideNodeID(3, iSide, Element)) - &
                         GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                    
                    DO
                      CALL RANDOM_NUMBER(RandVal)
                      IF (RandVal(1) + RandVal(2) .LE. 1.-epsi) EXIT
                    END DO
                    BN = GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element)) + RandVal(1)*SV1 + RandVal(2)*SV2
                    BN = BN - Species(iSpec)%Init(iInit)%BasePointIC
                    !Lokalisierung
                    det1 = BN(1)*BV2(2)*OV(3) + BN(2)*BV2(3)*OV(1) + &
                         BN(3)*BV2(1)*OV(2) - BN(3)*BV2(2)*OV(1) - &
                         BN(1)*BV2(3)*OV(2) - BN(2)*BV2(1)*OV(3)
                    det2 = BV1(1)*BN(2)*OV(3) + BV1(2)*BN(3)*OV(1) + &
                         BV1(3)*BN(1)*OV(2) - BV1(3)*BN(2)*OV(1) - &
                         BV1(1)*BN(3)*OV(2) - BV1(2)*BN(1)*OV(3)
                    det3 = BV1(1)*BV2(2)*BN(3) + BV1(2)*BV2(3)*BN(1) + &
                         BV1(3)*BV2(1)*BN(2) - BV1(3)*BV2(2)*BN(1) - &
                         BV1(1)*BV2(3)*BN(2) - BV1(2)*BV2(1)*BN(3)
                    
                    det1 = det1/Species(iSpec)%Init(iInit)%ConstPress%Determinant
                    det2 = det2/Species(iSpec)%Init(iInit)%ConstPress%Determinant
                    det3 = det3/Species(iSpec)%Init(iInit)%ConstPress%Determinant
                    
                    IF (((det1-0.5)**2 .LT. 0.25-epsi).AND.((det2-0.5)**2 .LT. 0.25-epsi).AND.((det3-0.5)**2 .LT. 0.25-epsi)) THEN
                      Species(iSpec)%Init(iInit)%ConstPress%ElemStat(ExamElem) = 2
                      Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside &
                           = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside + 1
                      TempElemPartlyInside(Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside) = ExamElem
                      EXIT
                    END IF
                  END DO
                END IF
              END DO
            END DO
            nInterOld = nInterest + 1
            nInterest = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside
          END DO
          
          
        CASE ('cylinder')
          !       Species(iSpec)%ConstPress%BV1 = Species(iSpec)%CylinderHeightIC * Species(iSpec)%NormalIC
          
          Species(iSpec)%Init(iInit)%ConstPress%OrthoVector(1) = &
               Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(3) - &
               Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(2)
          Species(iSpec)%Init(iInit)%ConstPress%OrthoVector(2) = &
               Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(1) - &
               Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(3)
          Species(iSpec)%Init(iInit)%ConstPress%OrthoVector(3) = &
               Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(2) - &
               Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(1)
          OV(1:3) = Species(iSpec)%Init(iInit)%ConstPress%OrthoVector(1:3)
          !       WRITE(*,*) 'iSpec', iSpec
          !       WRITE(*,*) 'BaseVec1', Species(iSpec)%BaseVector1IC
          !       WRITE(*,*) 'BaseVec2', Species(iSpec)%BaseVector2IC
          !       WRITE(*,*) 'BV1', Species(iSpec)%ConstPress%OrthoVector
          !
          IF ((OV(1) .EQ. 0) .AND. (OV(2) .EQ. 0) .AND. (OV(3) .EQ. 0)) THEN
            WRITE(*,*) 'Error in InitializeVariables: Cannot calculate NormalVector(Cyl) of InitVolume in EmissionCase(3 or 4)!'
            STOP
          END IF
          Species(iSpec)%Init(iInit)%ConstPress%OrthoVector &
                  = OV*(Species(iSpec)%Init(iInit)%CylinderHeightIC/SQRT(OV(1)**2 + OV(2)**2 + OV(3)**2))
          OV(1:3) = Species(iSpec)%Init(iInit)%ConstPress%OrthoVector(1:3)
          
          Species(iSpec)%Init(iInit)%ParticleEmission = &
               Species(iSpec)%Init(iInit)%ConstantPressure * Species(iSpec)%Init(iInit)%RadiusIC**2 * &
               3.1415926535 * Species(iSpec)%Init(iInit)%CylinderHeightIC / (BoltzmannConst * &
               Species(iSpec)%Init(iInit)%MWTemperatureIC * Species(iSpec)%MacroParticleFactor)
          Species(iSpec)%Init(iInit)%ConstPress%EkinInside = &
               1.5*Species(iSpec)%Init(iInit)%ConstantPressure*Species(iSpec)%Init(iInit)%RadiusIC**2 * &
               3.1415926535 * Species(iSpec)%Init(iInit)%CylinderHeightIC + &
               INT(Species(iSpec)%Init(iInit)%ParticleEmission)*0.5*Species(iSpec)%MassIC*Species(iSpec)%Init(iInit)%VeloIC**2 &
               * Species(iSpec)%MacroParticleFactor
          WRITE (*,'(A49,I3.3,A2,g0)') 'Number of Particles inside ConstPressArea ',iInit, ': ', &
               Species(iSpec)%Init(iInit)%ParticleEmission
          
          DO iElem = 1,nElems
            nNodesInside = 0
            nBoundNodes = 0
            DO iNode = 1,8
              BN  = GEO%NodeCoords(:, GEO%ElemToNodeID(iNode,iElem)) - Species(iSpec)%Init(iInit)%BasePointIC
              BV2(1) = BN(2) * OV(3) - BN(3) * OV(2)                   !Vector orthogonal on BN and NormalIC
              BV2(2) = BN(3) * OV(1) - BN(1) * OV(3)
              BV2(3) = BN(1) * OV(2) - BN(2) * OV(1)
              dist1 = SQRT((BV2(1)**2 + BV2(2)**2 + BV2(3)**2)/(OV(1)**2 + OV(2)**2 + OV(3)**2))
              
              IF (dist1 .LE. Species(iSpec)%Init(iInit)%RadiusIC + epsi) THEN
                BV3(1) = OV(2) * BV2(3) - OV(3) * BV2(2)
                BV3(2) = OV(3) * BV2(1) - OV(1) * BV2(3)
                BV3(3) = OV(1) * BV2(2) - OV(2) * BV2(1)
                IF (BV3(1)**2 + BV3(2)**2 + BV3(3)**2 .NE. 0.) THEN
                  BV3    = dist1 * BV3/SQRT(BV3(1)**2 + BV3(2)**2 + BV3(3)**2)   !Shortest Vector from Node to Cylinder-Axis
                ELSE
                  BV3(:) = 0.
                END IF
                IF (OV(1) .NE. 0.) THEN
                  dist2 = (BN(1) - BV3(1))/OV(1)
                ELSE IF (OV(2) .NE. 0.) THEN
                  dist2 = (BN(2) - BV3(2))/OV(2)
                ELSE IF (OV(3) .NE. 0.) THEN
                  dist2 = (BN(3) - BV3(3))/OV(3)
                ELSE
                  dist2 = 0.
                END IF
                
                IF (((dist2 .GE. 1.-epsi).AND.(dist2 .LE. 1+epsi)) .OR. ((dist2 .GE. -epsi).AND.(dist2 .LE. epsi))) THEN
                  nBoundNodes = nBoundNodes + 1
                ELSE IF ((dist1 .LE. Species(iSpec)%Init(iInit)%RadiusIC+epsi) &
                     .AND. (dist1 .GE. Species(iSpec)%Init(iInit)%RadiusIC-epsi) &
                     .AND. (dist2 .GT. 0.) .AND. (dist2 .LT. 1.)) THEN
                  nBoundNodes = nBoundNodes + 1
                ELSE IF ((dist2 .LT. 1.) .AND. (dist2 .GT. 0.)) THEN
                  nNodesInside = nNodesInside + 1
                END IF
              END IF
            END DO
            
            IF (nNodesInside .EQ. 8) THEN
              Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside = Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside + 1
              TempElemTotalInside(Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside)=iElem
              Species(iSpec)%Init(iInit)%ConstPress%ElemStat(iELem) = 1
            ELSE IF ((nNodesInside .GE. 1) .AND. (nNodesInside .LE. 7)) THEN
              IF (nNodesInside + nBoundNodes .EQ. 8) THEN
                Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside = Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside + 1
                TempElemTotalInside(Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside)=iElem
                Species(iSpec)%Init(iInit)%ConstPress%ElemStat(iELem) = 1
              ELSE
                Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside +1
                TempElemPartlyInside(Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside) = iElem
                Species(iSpec)%Init(iInit)%ConstPress%ElemStat(iELem) = 2
              END IF
            ELSE
              IF (nBoundNodes .EQ. 8) THEN
                Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside = Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside + 1
                TempElemTotalInside(Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside)=iElem
                Species(iSpec)%Init(iInit)%ConstPress%ElemStat(iELem) = 1
              ELSE
                Species(iSpec)%Init(iInit)%ConstPress%ElemStat(iELem) = 3
              END IF
            END IF
          END DO
          
          ! If no Element has been found, search for Basepoint ========================================
          IF ((Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside .EQ. 0) .AND. &
               (Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside .EQ. 0)) THEN
            CALL PointInsideQuad3D(iSpec,iInit,Element,InElementCheck,dete)
            IF (InElementCheck) THEN
              Species(iSpec)%Init(iInit)%ConstPress%ElemStat(Element) = 2
              Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside = 1
              TempElemPartlyInside(1) = Element
            END IF
          END IF
          
          !Schießen auf Nachbarzellen der totalen ===================================================
          DO iElem = 1,Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside
            Element = TempElemTotalInside(iElem)
            DO iSide = 1,6
              Side = ElemToSide(E2S_SIDE_ID, iSide, Element)
              IF (SideToElem(S2E_ELEM_ID,Side) .EQ. Element) THEN
                ExamElem = SideToElem(S2E_NB_ELEM_ID, Side)
              ELSE IF (SideToElem(S2E_NB_ELEM_ID,Side) .EQ. Element) THEN
                ExamElem = SideToElem(S2E_ELEM_ID,Side)
              END IF
              IF (Species(iSpec)%Init(iInit)%ConstPress%ElemStat(ExamElem) .EQ. 3) THEN
                DO iShot = 1,200
                  IF (iShot .LE. 100) THEN
                    !Shooting (123)
                    SV1 = GEO%NodeCoords(:,GEO%ElemSideNodeID(2, iSide, Element)) - &
                         GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                  ELSE
                    !Shooting (134)
                    SV1 = GEO%NodeCoords(:,GEO%ElemSideNodeID(4, iSide, Element)) - &
                         GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                  END IF
                  SV2 = GEO%NodeCoords(:,GEO%ElemSideNodeID(3, iSide, Element)) - &
                       GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                  DO
                    CALL RANDOM_NUMBER(RandVal)
                    IF (RandVal(1) + RandVal(2) .LE. 1) EXIT
                  END DO
                  BN = GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element)) + RandVal(1)*SV1 + RandVal(2)*SV2
                  BN = BN - Species(iSpec)%Init(iInit)%BasePointIC
                  !Lokalisierung
                  BV2(1) = BN(2) * OV(3) - BN(3) * OV(2)                   !Vector orthogonal on BN and NormalIC
                  BV2(2) = BN(3) * OV(1) - BN(1) * OV(3)
                  BV2(3) = BN(1) * OV(2) - BN(2) * OV(1)
                  dist1  = SQRT((BV2(1)**2 + BV2(2)**2 + BV2(3)**2)/ (OV(1)**2 + OV(2)**2 + OV(3)**2))
                  
                  IF (dist1 .LE. Species(iSpec)%Init(iInit)%RadiusIC + epsi) THEN
                    BV3(1) = OV(2) * BV2(3) - OV(3) * BV2(2)
                    BV3(2) = OV(3) * BV2(1) - OV(1) * BV2(3)
                    BV3(3) = OV(1) * BV2(2) - OV(2) * BV2(1)
                    IF (BV3(1)**2 + BV3(2)**2 + BV3(3)**2 .NE. 0.) THEN
                      BV3    = dist1 * BV3/SQRT(BV3(1)**2 + BV3(2)**2 + BV3(3)**2)   !Shortest Vector from Node to Cylinder-Axis
                    ELSE
                      BV3(:) = 0.
                    END IF
                    IF (OV(1) .NE. 0.) THEN
                      dist2 = (BN(1) - BV3(1))/OV(1)
                    ELSE IF (OV(2) .NE. 0.) THEN
                      dist2 = (BN(2) - BV3(2))/OV(2)
                    ELSE IF (OV(3) .NE. 0.) THEN
                      dist2 = (BN(3) - BV3(3))/OV(3)
                    ELSE
                      dist2 = 0.
                    END IF
                    
                    IF ((dist2 .LT. 1.-epsi) .AND. (dist2 .GT. epsi)) THEN
                      Species(iSpec)%Init(iInit)%ConstPress%ElemStat(ExamElem) = 2
                      Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside &
                           = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside + 1
                      TempElemPartlyInside(Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside) = ExamElem
                      EXIT
                    END IF
                  END IF
                END DO
              END IF
            END DO
          END DO
          
#ifdef MPI
          ! Schiessen auf MPI Seiten
          DO Side=1, nSides
            IF (Side.GT.nInnerSides+nBCSides) THEN
              IF (SideToElem(S2E_ELEM_ID,Side).NE.-1) THEN
                Element = SideToElem(S2E_ELEM_ID,Side)
                iSide = SideToElem(S2E_LOC_SIDE_ID,Side)
              ELSE
                Element = SideToElem(S2E_NB_ELEM_ID,Side)
                iSide = SideToElem(S2E_NB_LOC_SIDE_ID,Side)
              END IF
              IF (Species(iSpec)%Init(iInit)%ConstPress%ElemStat(Element) .EQ. 3) THEN
                DO iShot = 1,200
                  IF (iShot .LE. 100) THEN
                    !Shooting (123)
                    SV1 = GEO%NodeCoords(:,GEO%ElemSideNodeID(2, iSide, Element)) - &
                         GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                  ELSE
                    !Shooting (134)
                    SV1 = GEO%NodeCoords(:,GEO%ElemSideNodeID(4, iSide, Element)) - &
                         GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                  END IF
                  SV2 = GEO%NodeCoords(:,GEO%ElemSideNodeID(3, iSide, Element)) - &
                       GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                  DO
                    CALL RANDOM_NUMBER(RandVal)
                    IF (RandVal(1) + RandVal(2) .LE. 1) EXIT
                  END DO
                  BN = GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element)) + RandVal(1)*SV1 + RandVal(2)*SV2
                  BN = BN - Species(iSpec)%Init(iInit)%BasePointIC
                  !Lokalisierung
                  BV2(1) = BN(2) * OV(3) - BN(3) * OV(2)                   !Vector orthogonal on BN and NormalIC
                  BV2(2) = BN(3) * OV(1) - BN(1) * OV(3)
                  BV2(3) = BN(1) * OV(2) - BN(2) * OV(1)
                  dist1 = SQRT((BV2(1)**2 + BV2(2)**2 + BV2(3)**2)/(OV(1)**2 + OV(2)**2 + OV(3)**2))
                  IF (dist1 .LE. Species(iSpec)%Init(iInit)%RadiusIC + epsi) THEN
                    BV3(1) = OV(2) * BV2(3) - OV(3) * BV2(2)
                    BV3(2) = OV(3) * BV2(1) - OV(1) * BV2(3)
                    BV3(3) = OV(1) * BV2(2) - OV(2) * BV2(1)
                    IF (BV3(1)**2 + BV3(2)**2 + BV3(3)**2 .NE. 0.) THEN
                      BV3    = dist1 * BV3/SQRT(BV3(1)**2 + BV3(2)**2 + BV3(3)**2)   !Shortest Vector from Node to Cylinder-Axis
                    ELSE
                      BV3(:) = 0.
                    END IF
                    IF (OV(1) .NE. 0.) THEN
                      dist2 = (BN(1) - BV3(1))/OV(1)
                    ELSE IF (OV(2) .NE. 0.) THEN
                      dist2 = (BN(2) - BV3(2))/OV(2)
                    ELSE IF (OV(3) .NE. 0.) THEN
                      dist2 = (BN(3) - BV3(3))/OV(3)
                    ELSE
                      dist2 = 0.
                    END IF
                    IF ((dist2 .LT. 1.-epsi) .AND. (dist2 .GT. epsi)) THEN
                      Species(iSpec)%Init(iInit)%ConstPress%ElemStat(Element) = 2
                      Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside &
                           = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside + 1
                      TempElemPartlyInside(Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside) = Element
                      EXIT
                    END IF
                  END IF
                END DO
              END IF
            END IF
          END DO
#endif

          !!!Schießen auf Nachbarelemente der teilweisen ===============================================
          nInterest = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside
          nInterOld = 1
          DO WHILE (nInterest .GE. nInterOld)
            DO iElem = nInterOld,nInterest
              Element = TempElemPartlyInside(iElem)
              DO iSide = 1,6
                Side = ElemToSide(E2S_SIDE_ID, iSide, Element)
                IF (SideToElem(S2E_ELEM_ID,Side) .EQ. Element) THEN
                  ExamElem = SideToElem(S2E_NB_ELEM_ID, Side)
                ELSE IF (SideToElem(S2E_NB_ELEM_ID,Side) .EQ. Element) THEN
                  ExamElem = SideToElem(S2E_ELEM_ID,Side)
                END IF
                IF (Species(iSpec)%Init(iInit)%ConstPress%ElemStat(ExamElem) .EQ. 3) THEN
                  DO iShot = 1,200
                    IF (iShot .LE. 100) THEN
                      !Shooting (123)
                      SV1 = GEO%NodeCoords(:,GEO%ElemSideNodeID(2, iSide, Element)) - &
                           GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                    ELSE
                      !Shooting (134)
                      SV1 = GEO%NodeCoords(:,GEO%ElemSideNodeID(4, iSide, Element)) - &
                           GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                    END IF
                    SV2 = GEO%NodeCoords(:,GEO%ElemSideNodeID(3, iSide, Element)) - &
                         GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                    DO
                      CALL RANDOM_NUMBER(RandVal)
                      IF (RandVal(1) + RandVal(2) .LE. 1-epsi) EXIT
                    END DO
                    BN = GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element)) + RandVal(1)*SV1 + RandVal(2)*SV2
                    BN = BN - Species(iSpec)%Init(iInit)%BasePointIC
                    !Lokalisierung
                    BV2(1) = BN(2) * OV(3) - BN(3) * OV(2)                   !Vector orthogonal on BN and NormalIC
                    BV2(2) = BN(3) * OV(1) - BN(1) * OV(3)
                    BV2(3) = BN(1) * OV(2) - BN(2) * OV(1)
                    dist1 = SQRT((BV2(1)**2 + BV2(2)**2 + BV2(3)**2)/(OV(1)**2 + OV(2)**2 + OV(3)**2))
                    IF (dist1 .LE. Species(iSpec)%Init(iInit)%RadiusIC + epsi) THEN
                      BV3(1) = OV(2) * BV2(3) - OV(3) * BV2(2)
                      BV3(2) = OV(3) * BV2(1) - OV(1) * BV2(3)
                      BV3(3) = OV(1) * BV2(2) - OV(2) * BV2(1)
                      IF (BV3(1)**2 + BV3(2)**2 + BV3(3)**2 .NE. 0.) THEN
                        BV3    = dist1 * BV3/SQRT(BV3(1)**2 + BV3(2)**2 + BV3(3)**2)   !Shortest Vector from Node to Cylinder-Axis
                      ELSE
                        BV3(:) = 0.
                      END IF
                      IF (OV(1) .NE. 0.) THEN
                        dist2 = (BN(1) - BV3(1))/OV(1)
                      ELSE IF (OV(2) .NE. 0.) THEN
                        dist2 = (BN(2) - BV3(2))/OV(2)
                      ELSE IF (OV(3) .NE. 0.) THEN
                        dist2 = (BN(3) - BV3(3))/OV(3)
                      ELSE
                        dist2 = 0.
                      END IF
                      IF ((dist2 .LT. 1.-epsi) .AND. (dist2 .GT. epsi)) THEN
                        Species(iSpec)%Init(iInit)%ConstPress%ElemStat(ExamElem) = 2
                        Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside &
                             = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside + 1
                        TempElemPartlyInside(Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside) = ExamElem
                        EXIT
                      END IF
                    END IF
                  END DO
                END IF
              END DO
            END DO
            nInterOld = nInterest + 1
            nInterest = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside
          END DO
          
          
          
        END SELECT
        
        ALLOCATE (Species(iSpec)%Init(iInit)%ConstPress%ElemTotalInside(Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside))
        ALLOCATE (Species(iSpec)%Init(iInit)%ConstPress%ElemPartlyInside(Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside))
        Species(iSpec)%Init(iInit)%ConstPress%ElemTotalInside(1:Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside) = &
             TempElemTotalInside(1:Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside)
        Species(iSpec)%Init(iInit)%ConstPress%ElemPartlyInside(1:Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside) = &
             TempElemPartlyInside(1:Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside)
        DEALLOCATE (TempElemTotalInside)
        DEALLOCATE (TempElemPartlyInside)
        
        WRITE (*,'(A49,I3.3,A2,I0)') 'Number of Elements inside ConstPressArea ',iInit, ': ', &
             Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside
        WRITE (*,'(A49,I3.3,A2,I0)') 'Number of Elements partly inside ConstPressArea ',iInit, ': ', &
             Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside
        
      END IF
    END DO
  END DO
  
END SUBROUTINE ParticlePressureIni




SUBROUTINE ParticlePressureCellIni()

  USE MOD_Particle_Vars
  USE MOD_Globals
  USE MOD_Mesh_Vars,      ONLY : nElems,ElemToSide,SideToElem, nSides, nInnerSides, nBCSides
#ifdef MPI
  USE MOD_part_MPI_Vars, ONLY : PMPIVAR
#endif

  IMPLICIT NONE
 
  REAL                 :: BN(3), BV1(3), BV2(3), BV3(3), SV1(3), SV2(3), OV(3), epsi
  REAL                 :: det1, det2, det3, RandVal(2), dist1, dist2, dete(6,2)
  INTEGER, ALLOCATABLE :: TempElemTotalInside(:), TempElemPartlyInside(:)
  INTEGER              :: nNodesInside, nBoundNodes, nInterest, nInterOld, Element, Side, ExamElem
  INTEGER              :: iShot, iElem, iNode, iSpec, iInit, iSide, i, n
  LOGICAL              :: ElementFound, InElementCheck
  
  
  epsi = 100.*epsilon(epsi)
  DO iSpec = 1,nSpecies
    DO iInit = Species(iSpec)%StartnumberOfInits, Species(iSpec)%NumberOfInits
      Species(iSpec)%Init(iInit)%ConstPress%InitialTemp = Species(iSpec)%Init(iInit)%MWTemperatureIC
      IF ((Species(iSpec)%Init(iInit)%ParticleEmissionType .EQ. 4).OR.(Species(iSpec)%Init(iInit)%ParticleEmissionType .EQ. 6)) THEN
        ALLOCATE (TempElemTotalInside(nElems))
        ALLOCATE (TempElemPartlyInside(nElems))
        ALLOCATE (Species(iSpec)%Init(iInit)%ConstPress%ElemStat(nElems))
        
        Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside = 0
        Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside = 0
        Species(iSpec)%Init(iInit)%ConstPress%ElemStat(1:nElems) = 3
        SELECT CASE (TRIM(Species(iSpec)%Init(iInit)%SpaceIC))
        CASE ('cuboid')
          BV1 = Species(iSpec)%Init(iInit)%BaseVector1IC
          BV2 = Species(iSpec)%Init(iInit)%BaseVector2IC
          
          Species(iSpec)%Init(iInit)%ConstPress%OrthoVector(1) = BV1(2) * BV2(3) - BV1(3) * BV2(2)
          Species(iSpec)%Init(iInit)%ConstPress%OrthoVector(2) = BV1(3) * BV2(1) - BV1(1) * BV2(3)
          Species(iSpec)%Init(iInit)%ConstPress%OrthoVector(3) = BV1(1) * BV2(2) - BV1(2) * BV2(1)
          OV = Species(iSpec)%Init(iInit)%ConstPress%OrthoVector
          
          IF ((OV(1) .EQ. 0) .AND. (OV(2) .EQ. 0) .AND. (OV(3) .EQ. 0)) THEN
            WRITE(*,*) 'Error in InitializeVariables: Cannot calculate Normal Vector of InitVolume in EmissionCase(3 or 4)!'
            STOP
          END IF
          
          Species(iSpec)%Init(iInit)%ConstPress%OrthoVector &
               = OV*(Species(iSpec)%Init(iInit)%CuboidHeightIC/SQRT(OV(1)**2 + OV(2)**2 + OV(3)**2))
          OV = Species(iSpec)%Init(iInit)%ConstPress%OrthoVector
          Species(iSpec)%Init(iInit)%ConstPress%Determinant = ABS(BV1(1)*BV2(2)*OV(3) + BV1(2)*BV2(3)*OV(1) + BV1(3)*BV2(1)*OV(2) &
               - BV1(3)*BV2(2)*OV(1) - BV1(1)*BV2(3)*OV(2) - BV1(2)*BV2(1)*OV(3))
          !WRITE (*,*) 'ConstPress', Species(iSpec)%ConstantPressure , 'Deter', Species(iSpec)%ConstPress%Determinant, 'ispec', &
          !iSpec &
          !,'MPF',Species(iSpec)%MacroParticleFactor
          
          ! ParticleEmission multiplied by cell volume equals particles to be in cells
          Species(iSpec)%Init(iInit)%ParticleEmission = Species(iSpec)%Init(iInit)%ConstantPressure/ &
               (BoltzmannConst * Species(iSpec)%Init(iInit)%MWTemperatureIC * Species(iSpec)%MacroParticleFactor)
          !        Species(iSpec)%ConstPress%EkinInside =1.5*Species(iSpec)%ConstantPressure*Species(iSpec)%ConstPress%Determinant+&
          !                            INT(Species(iSpec)%ParticleEmission)*0.5*Species(iSpec)%MassIC*Species(iSpec)%VeloIC**2* &
          !        Species(iSpec)%MacroParticleFactor
          
          !        WRITE(*,*) 'Number of Particles in Constant-Pressure-Area:', Species(iSpec)%ParticleEmission
          
          DO iElem = 1,nElems
            nNodesInside = 0
            nBoundNodes = 0
            DO iNode = 1,8
              BN(1:3) = GEO%NodeCoords(:,GEO%ElemToNodeID(iNode,iElem)) - Species(iSpec)%Init(iInit)%BasePointIC
              det1 = BN(1)*BV2(2)*OV(3) + BN(2)*BV2(3)*OV(1) + BN(3)*BV2(1)*OV(2) - &
                   BN(3)*BV2(2)*OV(1) - BN(1)*BV2(3)*OV(2) - BN(2)*BV2(1)*OV(3)
              det2 = BV1(1)*BN(2)*OV(3) + BV1(2)*BN(3)*OV(1) + BV1(3)*BN(1)*OV(2) - &
                   BV1(3)*BN(2)*OV(1) - BV1(1)*BN(3)*OV(2) - BV1(2)*BN(1)*OV(3)
              det3 = BV1(1)*BV2(2)*BN(3) + BV1(2)*BV2(3)*BN(1) + BV1(3)*BV2(1)*BN(2) - &
                   BV1(3)*BV2(2)*BN(1) - BV1(1)*BV2(3)*BN(2) - BV1(2)*BV2(1)*BN(3)
              
              det1 = det1/Species(iSpec)%Init(iInit)%ConstPress%Determinant
              det2 = det2/Species(iSpec)%Init(iInit)%ConstPress%Determinant
              det3 = det3/Species(iSpec)%Init(iInit)%ConstPress%Determinant
              
              IF (((det1 .LE. 1. + epsi) .AND. (det1 .GE. 1.-epsi)) .OR. ((det1 .LE. epsi) .AND. (det1 .GE. -epsi))) THEN
                IF (((det2 .LE. 1.+epsi) .AND. (det2 .GE. -epsi)) .AND. ((det3 .LE. 1.+epsi) .AND. (det3 .GE. -epsi))) THEN
                  nBoundNodes = nBoundNodes + 1
                END IF
              ELSE IF (((det2 .LE. 1.+epsi).AND.(det2 .GE.1-epsi)) .OR. ((det2 .GE. -epsi).AND.(det2 .LE. epsi))) THEN
                IF ((det1 .LE. 1.+epsi) .AND. (det1 .GE. -epsi) .AND. (det3 .LE. 1.+epsi) .AND. (det3 .GE. -epsi)) THEN
                  nBoundNodes = nBoundNodes + 1
                END IF
              ELSE IF (((det3 .GE. 1.-epsi).AND.(det3 .LE.1+epsi)) .OR. ((det3 .GE.-epsi).AND.(det3.LE.epsi))) THEN
                IF ((det2 .LE. 1.+epsi) .AND. (det2 .GE. -epsi) .AND. (det1 .LE. 1.+epsi) .AND. (det1 .GE. -epsi)) THEN
                  nBoundNodes = nBoundNodes + 1
                END IF
              ELSE IF ((det1 .LT. 1.) .AND. (det1 .GT. 0.) .AND. (det2 .LT. 1.) .AND. (det2 .GT. 0.) &
                   .AND. (det3 .LT. 1.) .AND. (det3 .GT. 0.)) THEN
                nNodesInside = nNodesInside + 1
              END IF
            END DO
            
            IF (nNodesInside .EQ. 8) THEN
              Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside = Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside + 1
              TempElemTotalInside(Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside)=iElem
              Species(iSpec)%Init(iInit)%ConstPress%ElemStat(iELem) = 1
            ELSE IF ((nNodesInside .GE. 1) .AND. (nNodesInside .LE. 7)) THEN
              IF (nNodesInside + nBoundNodes .EQ. 8) THEN
                Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside = Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside + 1
                TempElemTotalInside(Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside)=iElem
                Species(iSpec)%Init(iInit)%ConstPress%ElemStat(iELem) = 1
              ELSE
                Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside +1
                TempElemPartlyInside(Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside) = iElem
                Species(iSpec)%Init(iInit)%ConstPress%ElemStat(iELem) = 2
              END IF
            ELSE
              IF (nBoundNodes .EQ. 8) THEN
                Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside = Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside + 1
                TempElemTotalInside(Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside)=iElem
                Species(iSpec)%Init(iInit)%ConstPress%ElemStat(iELem) = 1
              ELSE
                Species(iSpec)%Init(iInit)%ConstPress%ElemStat(iELem) = 3
              END IF
            END IF
          END DO
          
          ! If no Element has been found, search for Basepoint ========================================
          IF ((Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside .EQ. 0) .AND. &
               (Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside .EQ. 0)) THEN
            CALL PointInsideQuad3D(iSpec,iInit,Element,InElementCheck,dete)
            IF (InElementCheck) THEN
              Species(iSpec)%Init(iInit)%ConstPress%ElemStat(Element) = 2
              Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside = 1
              TempElemPartlyInside(1) = Element
            END IF
          END IF
          
          ! Shoot on Sides of Neighbour-Elements of totally inside Elements ===========================
          
          DO iElem = 1,Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside
            Element = TempElemTotalInside(iElem)
            DO iSide = 1,6
              Side = ElemToSide(E2S_SIDE_ID, iSide, Element)
              IF (SideToElem(S2E_ELEM_ID,Side) .EQ. Element) THEN
                ExamElem = SideToElem(S2E_NB_ELEM_ID, Side)
              ELSE IF (SideToElem(S2E_NB_ELEM_ID,Side) .EQ. Element) THEN
                ExamElem = SideToElem(S2E_ELEM_ID,Side)
              END IF
              IF (Species(iSpec)%Init(iInit)%ConstPress%ElemStat(ExamElem) .EQ. 3) THEN
                DO iShot = 1,200
                  IF (iShot .LE. 100) THEN
                    !Shooting (123)
                    SV1 = GEO%NodeCoords(:,GEO%ElemSideNodeID(2, iSide, Element)) - &
                         GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                  ELSE
                    !Shooting (134)
                    SV1 = GEO%NodeCoords(:,GEO%ElemSideNodeID(4, iSide, Element)) - &
                         GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                  END IF
                  SV2 = GEO%NodeCoords(:,GEO%ElemSideNodeID(3, iSide, Element)) - &
                       GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                  DO
                    CALL RANDOM_NUMBER(RandVal)
                    IF (RandVal(1) + RandVal(2) .LE. 1) EXIT
                  END DO
                  BN = GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element)) + RandVal(1)*SV1 + RandVal(2)*SV2
                  BN = BN - Species(iSpec)%Init(iInit)%BasePointIC
                  !Lokalisierung
                  det1 = BN(1)*BV2(2)*OV(3) + BN(2)*BV2(3)*OV(1) + &
                       BN(3)*BV2(1)*OV(2) - BN(3)*BV2(2)*OV(1) - &
                       BN(1)*BV2(3)*OV(2) - BN(2)*BV2(1)*OV(3)
                  det2 = BV1(1)*BN(2)*OV(3) + BV1(2)*BN(3)*OV(1) + &
                       BV1(3)*BN(1)*OV(2) - BV1(3)*BN(2)*OV(1) - &
                       BV1(1)*BN(3)*OV(2) - BV1(2)*BN(1)*OV(3)
                  det3 = BV1(1)*BV2(2)*BN(3) + BV1(2)*BV2(3)*BN(1) + &
                       BV1(3)*BV2(1)*BN(2) - BV1(3)*BV2(2)*BN(1) - &
                       BV1(1)*BV2(3)*BN(2) - BV1(2)*BV2(1)*BN(3)
                  
                  det1 = det1/Species(iSpec)%Init(iInit)%ConstPress%Determinant
                  det2 = det2/Species(iSpec)%Init(iInit)%ConstPress%Determinant
                  det3 = det3/Species(iSpec)%Init(iInit)%ConstPress%Determinant
                  
                  IF ((det1 .LT. 1.-epsi) .AND. (det1 .GT. epsi) .AND. (det2 .LT. 1-epsi) .AND. &
                       (det2 .GT. epsi) .AND. (det3 .LT. 1.-epsi) .AND. (det3 .GT. epsi)) THEN
                    Species(iSpec)%Init(iInit)%ConstPress%ElemStat(ExamElem) = 2
                    Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside &
                         = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside + 1
                    TempElemPartlyInside(Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside) = ExamElem
                    EXIT
                  END IF
                END DO
              END IF
            END DO
          END DO
          
#ifdef MPI

          DO Side=1, nSides
            IF (Side.GT.nInnerSides+nBCSides) THEN
              IF (SideToElem(S2E_ELEM_ID,Side).NE.-1) THEN
                Element = SideToElem(S2E_ELEM_ID,Side)
                iSide = SideToElem(S2E_LOC_SIDE_ID,Side)
              ELSE
                Element = SideToElem(S2E_NB_ELEM_ID,Side)
                iSide = SideToElem(S2E_NB_LOC_SIDE_ID,Side)
              END IF
              IF (Species(iSpec)%Init(iInit)%ConstPress%ElemStat(Element) .EQ. 3) THEN
                DO iShot = 1,200
                  IF (iShot .LE. 100) THEN
                    !Shooting (123)
                    SV1 = GEO%NodeCoords(:,GEO%ElemSideNodeID(2, iSide, Element)) - &
                         GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                  ELSE
                    !Shooting (134)
                    SV1 = GEO%NodeCoords(:,GEO%ElemSideNodeID(4, iSide, Element)) - &
                         GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                  END IF
                  SV2 = GEO%NodeCoords(:,GEO%ElemSideNodeID(3, iSide, Element)) - &
                       GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                  DO
                    CALL RANDOM_NUMBER(RandVal)
                    IF (RandVal(1) + RandVal(2) .LE. 1) EXIT
                  END DO
                  BN = GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element)) + RandVal(1)*SV1 + RandVal(2)*SV2
                  BN = BN - Species(iSpec)%Init(iInit)%BasePointIC
                  !Lokalisierung
                  det1 = BN(1)*BV2(2)*OV(3) + BN(2)*BV2(3)*OV(1) + &
                       BN(3)*BV2(1)*OV(2) - BN(3)*BV2(2)*OV(1) - &
                       BN(1)*BV2(3)*OV(2) - BN(2)*BV2(1)*OV(3)
                  det2 = BV1(1)*BN(2)*OV(3) + BV1(2)*BN(3)*OV(1) + &
                       BV1(3)*BN(1)*OV(2) - BV1(3)*BN(2)*OV(1) - &
                       BV1(1)*BN(3)*OV(2) - BV1(2)*BN(1)*OV(3)
                  det3 = BV1(1)*BV2(2)*BN(3) + BV1(2)*BV2(3)*BN(1) + &
                       BV1(3)*BV2(1)*BN(2) - BV1(3)*BV2(2)*BN(1) - &
                       BV1(1)*BV2(3)*BN(2) - BV1(2)*BV2(1)*BN(3)
                  
                  det1 = det1/Species(iSpec)%Init(iInit)%ConstPress%Determinant
                  det2 = det2/Species(iSpec)%Init(iInit)%ConstPress%Determinant
                  det3 = det3/Species(iSpec)%Init(iInit)%ConstPress%Determinant
                  
                  IF ((det1 .LT. 1.-epsi) .AND. (det1 .GT. epsi) .AND. (det2 .LT. 1-epsi) .AND. &
                       (det2 .GT. epsi) .AND. (det3 .LT. 1.-epsi) .AND. (det3 .GT. epsi)) THEN
                    Species(iSpec)%Init(iInit)%ConstPress%ElemStat(Element) = 2
                    Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside &
                         = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside + 1
                    TempElemPartlyInside(Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside) = Element
                    EXIT
                  END IF
                END DO
              END IF
            END IF
          END DO
#endif
          
          ! Shoot on Sides of Neighbour-Elements of partly inside Elements ==========================================
          nInterest = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside
          nInterOld = 1
          DO WHILE (nInterest .GE. nInterOld)
            DO iElem = nInterOld, nInterest
              
              Element = TempElemPartlyInside(iElem)
              DO iSide = 1, 6
                Side = ElemToSide(E2S_SIDE_ID, iSide, Element)
                IF (SideToElem(S2E_ELEM_ID,Side) .EQ. Element) THEN
                  ExamElem = SideToElem(S2E_NB_ELEM_ID, Side)
                ELSE IF (SideToElem(S2E_NB_ELEM_ID,Side) .EQ. Element) THEN
                  ExamElem = SideToElem(S2E_ELEM_ID,Side)
                END IF
                
                IF (Species(iSpec)%Init(iInit)%ConstPress%ElemStat(ExamElem) .EQ. 3) THEN
                  DO iShot = 1,200
                    IF (iShot .LE. 100) THEN
                      !Shooting (123)
                      SV1 = GEO%NodeCoords(:,GEO%ElemSideNodeID(2, iSide, Element)) - &
                           GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                    ELSE
                      !Shooting (134)
                      SV1 = GEO%NodeCoords(:,GEO%ElemSideNodeID(4, iSide, Element)) - &
                           GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                    END IF
                    SV2 = GEO%NodeCoords(:,GEO%ElemSideNodeID(3, iSide, Element)) - &
                         GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                    
                    DO
                      CALL RANDOM_NUMBER(RandVal)
                      IF (RandVal(1) + RandVal(2) .LE. 1.-epsi) EXIT
                    END DO
                    BN = GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element)) + RandVal(1)*SV1 + RandVal(2)*SV2
                    BN = BN - Species(iSpec)%Init(iInit)%BasePointIC
                    !Lokalisierung
                    det1 = BN(1)*BV2(2)*OV(3) + BN(2)*BV2(3)*OV(1) + &
                         BN(3)*BV2(1)*OV(2) - BN(3)*BV2(2)*OV(1) - &
                         BN(1)*BV2(3)*OV(2) - BN(2)*BV2(1)*OV(3)
                    det2 = BV1(1)*BN(2)*OV(3) + BV1(2)*BN(3)*OV(1) + &
                         BV1(3)*BN(1)*OV(2) - BV1(3)*BN(2)*OV(1) - &
                         BV1(1)*BN(3)*OV(2) - BV1(2)*BN(1)*OV(3)
                    det3 = BV1(1)*BV2(2)*BN(3) + BV1(2)*BV2(3)*BN(1) + &
                         BV1(3)*BV2(1)*BN(2) - BV1(3)*BV2(2)*BN(1) - &
                         BV1(1)*BV2(3)*BN(2) - BV1(2)*BV2(1)*BN(3)
                    
                    det1 = det1/Species(iSpec)%Init(iInit)%ConstPress%Determinant
                    det2 = det2/Species(iSpec)%Init(iInit)%ConstPress%Determinant
                    det3 = det3/Species(iSpec)%Init(iInit)%ConstPress%Determinant
                    
                    IF (((det1-0.5)**2 .LT. 0.25-epsi).AND.((det2-0.5)**2 .LT. 0.25-epsi).AND.((det3-0.5)**2 .LT. 0.25-epsi)) THEN
                      Species(iSpec)%Init(iInit)%ConstPress%ElemStat(ExamElem) = 2
                      Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside &
                           = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside + 1
                      TempElemPartlyInside(Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside) = ExamElem
                      EXIT
                    END IF
                  END DO
                END IF
              END DO
            END DO
            nInterOld = nInterest + 1
            nInterest = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside
          END DO
          
          
        CASE ('cylinder')
          !       Species(iSpec)%ConstPress%BV1 = Species(iSpec)%CylinderHeightIC * Species(iSpec)%NormalIC
          
          Species(iSpec)%Init(iInit)%ConstPress%OrthoVector(1) = &
               Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(3) - &
               Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(2)
          Species(iSpec)%Init(iInit)%ConstPress%OrthoVector(2) = &
               Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(1) - &
               Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(3)
          Species(iSpec)%Init(iInit)%ConstPress%OrthoVector(3) = &
               Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(2) - &
               Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(1)
          OV(1:3) = Species(iSpec)%Init(iInit)%ConstPress%OrthoVector(1:3)
          !       WRITE(*,*) 'iSpec', iSpec
          !       WRITE(*,*) 'BaseVec1', Species(iSpec)%BaseVector1IC
          !       WRITE(*,*) 'BaseVec2', Species(iSpec)%BaseVector2IC
          !       WRITE(*,*) 'BV1', Species(iSpec)%ConstPress%OrthoVector
          !
          IF ((OV(1) .EQ. 0) .AND. (OV(2) .EQ. 0) .AND. (OV(3) .EQ. 0)) THEN
            WRITE(*,*) 'Error in InitializeVariables: Cannot calculate NormalVector(Cyl) of InitVolume in EmissionCase(3 or 4)!'
            STOP
          END IF
          Species(iSpec)%Init(iInit)%ConstPress%OrthoVector = &
               OV*(Species(iSpec)%Init(iInit)%CylinderHeightIC/SQRT(OV(1)**2 + OV(2)**2 + OV(3)**2))
          OV(1:3) = Species(iSpec)%Init(iInit)%ConstPress%OrthoVector(1:3)
          
          Species(iSpec)%Init(iInit)%ParticleEmission = &
               Species(iSpec)%Init(iInit)%ConstantPressure * Species(iSpec)%Init(iInit)%RadiusIC**2 * &
               3.1415926535 * Species(iSpec)%Init(iInit)%CylinderHeightIC / (BoltzmannConst * &
               Species(iSpec)%Init(iInit)%MWTemperatureIC * Species(iSpec)%MacroParticleFactor)
          Species(iSpec)%Init(iInit)%ConstPress%EkinInside = &
               1.5*Species(iSpec)%Init(iInit)%ConstantPressure*Species(iSpec)%Init(iInit)%RadiusIC**2 * &
               3.1415926535 * Species(iSpec)%Init(iInit)%CylinderHeightIC + &
               INT(Species(iSpec)%Init(iInit)%ParticleEmission)*0.5*Species(iSpec)%MassIC*Species(iSpec)%Init(iInit)%VeloIC**2 &
               * Species(iSpec)%MacroParticleFactor
          !        WRITE(*,*) 'Number of Particles in Constant-Pressure-Area:', Species(iSpec)%ParticleEmission
          
          DO iElem = 1,nElems
            nNodesInside = 0
            nBoundNodes = 0
            DO iNode = 1,8
              BN  = GEO%NodeCoords(:, GEO%ElemToNodeID(iNode,iElem)) - Species(iSpec)%Init(iInit)%BasePointIC
              BV2(1) = BN(2) * OV(3) - BN(3) * OV(2)                   !Vector orthogonal on BN and NormalIC
              BV2(2) = BN(3) * OV(1) - BN(1) * OV(3)
              BV2(3) = BN(1) * OV(2) - BN(2) * OV(1)
              dist1 = SQRT((BV2(1)**2 + BV2(2)**2 + BV2(3)**2)/(OV(1)**2 + OV(2)**2 + OV(3)**2))
              
              IF (dist1 .LE. Species(iSpec)%Init(iInit)%RadiusIC + epsi) THEN
                BV3(1) = OV(2) * BV2(3) - OV(3) * BV2(2)
                BV3(2) = OV(3) * BV2(1) - OV(1) * BV2(3)
                BV3(3) = OV(1) * BV2(2) - OV(2) * BV2(1)
                IF (BV3(1)**2 + BV3(2)**2 + BV3(3)**2 .NE. 0.) THEN
                  BV3    = dist1 * BV3/SQRT(BV3(1)**2 + BV3(2)**2 + BV3(3)**2)   !Shortest Vector from Node to Cylinder-Axis
                ELSE
                  BV3(:) = 0.
                END IF
                IF (OV(1) .NE. 0.) THEN
                  dist2 = (BN(1) - BV3(1))/OV(1)
                ELSE IF (OV(2) .NE. 0.) THEN
                  dist2 = (BN(2) - BV3(2))/OV(2)
                ELSE IF (OV(3) .NE. 0.) THEN
                  dist2 = (BN(3) - BV3(3))/OV(3)
                ELSE
                  dist2 = 0.
                END IF
                
                IF (((dist2 .GE. 1.-epsi).AND.(dist2 .LE. 1+epsi)) .OR. ((dist2 .GE. -epsi).AND.(dist2 .LE. epsi))) THEN
                  nBoundNodes = nBoundNodes + 1
                ELSE IF ((dist1 .LE. Species(iSpec)%Init(iInit)%RadiusIC+epsi) &
                     .AND. (dist1 .GE. Species(iSpec)%Init(iInit)%RadiusIC-epsi) &
                     .AND. (dist2 .GT. 0.) .AND. (dist2 .LT. 1.)) THEN
                  nBoundNodes = nBoundNodes + 1
                ELSE IF ((dist2 .LT. 1.) .AND. (dist2 .GT. 0.)) THEN
                  nNodesInside = nNodesInside + 1
                END IF
              END IF
            END DO
            
            IF (nNodesInside .EQ. 8) THEN
              Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside = Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside + 1
              TempElemTotalInside(Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside)=iElem
              Species(iSpec)%Init(iInit)%ConstPress%ElemStat(iELem) = 1
            ELSE IF ((nNodesInside .GE. 1) .AND. (nNodesInside .LE. 7)) THEN
              IF (nNodesInside + nBoundNodes .EQ. 8) THEN
                Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside = Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside + 1
                TempElemTotalInside(Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside)=iElem
                Species(iSpec)%Init(iInit)%ConstPress%ElemStat(iELem) = 1
              ELSE
                Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside +1
                TempElemPartlyInside(Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside) = iElem
                Species(iSpec)%Init(iInit)%ConstPress%ElemStat(iELem) = 2
              END IF
            ELSE
              IF (nBoundNodes .EQ. 8) THEN
                Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside = Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside + 1
                TempElemTotalInside(Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside)=iElem
                Species(iSpec)%Init(iInit)%ConstPress%ElemStat(iELem) = 1
              ELSE
                Species(iSpec)%Init(iInit)%ConstPress%ElemStat(iELem) = 3
              END IF
            END IF
          END DO
          
          ! If no Element has been found, search for Basepoint ========================================
          IF ((Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside .EQ. 0) .AND. &
               (Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside .EQ. 0)) THEN
            CALL PointInsideQuad3D(iSpec,iInit,Element,InElementCheck,dete)
            IF (InElementCheck) THEN
              Species(iSpec)%Init(iInit)%ConstPress%ElemStat(Element) = 2
              Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside = 1
              TempElemPartlyInside(1) = Element
            END IF
          END IF
          
          !Schießen auf Nachbarzellen der totalen ===================================================
          DO iElem = 1,Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside
            Element = TempElemTotalInside(iElem)
            DO iSide = 1,6
              Side = ElemToSide(E2S_SIDE_ID, iSide, Element)
              IF (SideToElem(S2E_ELEM_ID,Side) .EQ. Element) THEN
                ExamElem = SideToElem(S2E_NB_ELEM_ID, Side)
              ELSE IF (SideToElem(S2E_NB_ELEM_ID,Side) .EQ. Element) THEN
                ExamElem = SideToElem(S2E_ELEM_ID,Side)
              END IF
              IF (Species(iSpec)%Init(iInit)%ConstPress%ElemStat(ExamElem) .EQ. 3) THEN
                DO iShot = 1,200
                  IF (iShot .LE. 100) THEN
                    !Shooting (123)
                    SV1 = GEO%NodeCoords(:,GEO%ElemSideNodeID(2, iSide, Element)) - &
                         GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                  ELSE
                    !Shooting (134)
                    SV1 = GEO%NodeCoords(:,GEO%ElemSideNodeID(4, iSide, Element)) - &
                         GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                  END IF
                  SV2 = GEO%NodeCoords(:,GEO%ElemSideNodeID(3, iSide, Element)) - &
                       GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                  DO
                    CALL RANDOM_NUMBER(RandVal)
                    IF (RandVal(1) + RandVal(2) .LE. 1) EXIT
                  END DO
                  BN = GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element)) + RandVal(1)*SV1 + RandVal(2)*SV2
                  BN = BN - Species(iSpec)%Init(iInit)%BasePointIC
                  !Lokalisierung
                  BV2(1) = BN(2) * OV(3) - BN(3) * OV(2)                   !Vector orthogonal on BN and NormalIC
                  BV2(2) = BN(3) * OV(1) - BN(1) * OV(3)
                  BV2(3) = BN(1) * OV(2) - BN(2) * OV(1)
                  dist1  = SQRT((BV2(1)**2 + BV2(2)**2 + BV2(3)**2)/ (OV(1)**2 + OV(2)**2 + OV(3)**2))
                  
                  IF (dist1 .LE. Species(iSpec)%Init(iInit)%RadiusIC + epsi) THEN
                    BV3(1) = OV(2) * BV2(3) - OV(3) * BV2(2)
                    BV3(2) = OV(3) * BV2(1) - OV(1) * BV2(3)
                    BV3(3) = OV(1) * BV2(2) - OV(2) * BV2(1)
                    IF (BV3(1)**2 + BV3(2)**2 + BV3(3)**2 .NE. 0.) THEN
                      BV3    = dist1 * BV3/SQRT(BV3(1)**2 + BV3(2)**2 + BV3(3)**2)   !Shortest Vector from Node to Cylinder-Axis
                    ELSE
                      BV3(:) = 0.
                    END IF
                    IF (OV(1) .NE. 0.) THEN
                      dist2 = (BN(1) - BV3(1))/OV(1)
                    ELSE IF (OV(2) .NE. 0.) THEN
                      dist2 = (BN(2) - BV3(2))/OV(2)
                    ELSE IF (OV(3) .NE. 0.) THEN
                      dist2 = (BN(3) - BV3(3))/OV(3)
                    ELSE
                      dist2 = 0.
                    END IF
                    
                    IF ((dist2 .LT. 1.-epsi) .AND. (dist2 .GT. epsi)) THEN
                      Species(iSpec)%Init(iInit)%ConstPress%ElemStat(ExamElem) = 2
                      Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside &
                           = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside + 1
                      TempElemPartlyInside(Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside) = ExamElem
                      EXIT
                    END IF
                  END IF
                END DO
              END IF
            END DO
          END DO
          
#ifdef MPI
          ! Schiessen auf MPI Seiten
          DO Side=1, nSides
            IF (Side.GT.nInnerSides+nBCSides) THEN
              IF (SideToElem(S2E_ELEM_ID,Side).NE.-1) THEN
                Element = SideToElem(S2E_ELEM_ID,Side)
                iSide = SideToElem(S2E_LOC_SIDE_ID,Side)
              ELSE
                Element = SideToElem(S2E_NB_ELEM_ID,Side)
                iSide = SideToElem(S2E_NB_LOC_SIDE_ID,Side)
              END IF
              IF (Species(iSpec)%Init(iInit)%ConstPress%ElemStat(Element) .EQ. 3) THEN
                DO iShot = 1,200
                  IF (iShot .LE. 100) THEN
                    !Shooting (123)
                    SV1 = GEO%NodeCoords(:,GEO%ElemSideNodeID(2, iSide, Element)) - &
                         GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                  ELSE
                    !Shooting (134)
                    SV1 = GEO%NodeCoords(:,GEO%ElemSideNodeID(4, iSide, Element)) - &
                         GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                  END IF
                  SV2 = GEO%NodeCoords(:,GEO%ElemSideNodeID(3, iSide, Element)) - &
                       GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                  DO
                    CALL RANDOM_NUMBER(RandVal)
                    IF (RandVal(1) + RandVal(2) .LE. 1) EXIT
                  END DO
                  BN = GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element)) + RandVal(1)*SV1 + RandVal(2)*SV2
                  BN = BN - Species(iSpec)%Init(iInit)%BasePointIC
                  !Lokalisierung
                  BV2(1) = BN(2) * OV(3) - BN(3) * OV(2)                   !Vector orthogonal on BN and NormalIC
                  BV2(2) = BN(3) * OV(1) - BN(1) * OV(3)
                  BV2(3) = BN(1) * OV(2) - BN(2) * OV(1)
                  dist1 = SQRT((BV2(1)**2 + BV2(2)**2 + BV2(3)**2)/(OV(1)**2 + OV(2)**2 + OV(3)**2))
                  IF (dist1 .LE. Species(iSpec)%Init(iInit)%RadiusIC + epsi) THEN
                    BV3(1) = OV(2) * BV2(3) - OV(3) * BV2(2)
                    BV3(2) = OV(3) * BV2(1) - OV(1) * BV2(3)
                    BV3(3) = OV(1) * BV2(2) - OV(2) * BV2(1)
                    IF (BV3(1)**2 + BV3(2)**2 + BV3(3)**2 .NE. 0.) THEN
                      BV3    = dist1 * BV3/SQRT(BV3(1)**2 + BV3(2)**2 + BV3(3)**2)   !Shortest Vector from Node to Cylinder-Axis
                    ELSE
                      BV3(:) = 0.
                    END IF
                    IF (OV(1) .NE. 0.) THEN
                      dist2 = (BN(1) - BV3(1))/OV(1)
                    ELSE IF (OV(2) .NE. 0.) THEN
                      dist2 = (BN(2) - BV3(2))/OV(2)
                    ELSE IF (OV(3) .NE. 0.) THEN
                      dist2 = (BN(3) - BV3(3))/OV(3)
                    ELSE
                      dist2 = 0.
                    END IF
                    IF ((dist2 .LT. 1.-epsi) .AND. (dist2 .GT. epsi)) THEN
                      Species(iSpec)%Init(iInit)%ConstPress%ElemStat(Element) = 2
                      Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside &
                           = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside + 1
                      TempElemPartlyInside(Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside) = Element
                      EXIT
                    END IF
                  END IF
                END DO
              END IF
            END IF
          END DO
#endif

          !!!Schießen auf Nachbarelemente der teilweisen ===============================================
          nInterest = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside
          nInterOld = 1
          DO WHILE (nInterest .GE. nInterOld)
            DO iElem = nInterOld,nInterest
              Element = TempElemPartlyInside(iElem)
              DO iSide = 1,6
                Side = ElemToSide(E2S_SIDE_ID, iSide, Element)
                IF (SideToElem(S2E_ELEM_ID,Side) .EQ. Element) THEN
                  ExamElem = SideToElem(S2E_NB_ELEM_ID, Side)
                ELSE IF (SideToElem(S2E_NB_ELEM_ID,Side) .EQ. Element) THEN
                  ExamElem = SideToElem(S2E_ELEM_ID,Side)
                END IF
                IF (Species(iSpec)%Init(iInit)%ConstPress%ElemStat(ExamElem) .EQ. 3) THEN
                  DO iShot = 1,200
                    IF (iShot .LE. 100) THEN
                      !Shooting (123)
                      SV1 = GEO%NodeCoords(:,GEO%ElemSideNodeID(2, iSide, Element)) - &
                           GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                    ELSE
                      !Shooting (134)
                      SV1 = GEO%NodeCoords(:,GEO%ElemSideNodeID(4, iSide, Element)) - &
                           GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                    END IF
                    SV2 = GEO%NodeCoords(:,GEO%ElemSideNodeID(3, iSide, Element)) - &
                         GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element))
                    DO
                      CALL RANDOM_NUMBER(RandVal)
                      IF (RandVal(1) + RandVal(2) .LE. 1-epsi) EXIT
                    END DO
                    BN = GEO%NodeCoords(:,GEO%ElemSideNodeID(1, iSide, Element)) + RandVal(1)*SV1 + RandVal(2)*SV2
                    BN = BN - Species(iSpec)%Init(iInit)%BasePointIC
                    !Lokalisierung
                    BV2(1) = BN(2) * OV(3) - BN(3) * OV(2)                   !Vector orthogonal on BN and NormalIC
                    BV2(2) = BN(3) * OV(1) - BN(1) * OV(3)
                    BV2(3) = BN(1) * OV(2) - BN(2) * OV(1)
                    dist1 = SQRT((BV2(1)**2 + BV2(2)**2 + BV2(3)**2)/(OV(1)**2 + OV(2)**2 + OV(3)**2))
                    IF (dist1 .LE. Species(iSpec)%Init(iInit)%RadiusIC + epsi) THEN
                      BV3(1) = OV(2) * BV2(3) - OV(3) * BV2(2)
                      BV3(2) = OV(3) * BV2(1) - OV(1) * BV2(3)
                      BV3(3) = OV(1) * BV2(2) - OV(2) * BV2(1)
                      IF (BV3(1)**2 + BV3(2)**2 + BV3(3)**2 .NE. 0.) THEN
                        BV3    = dist1 * BV3/SQRT(BV3(1)**2 + BV3(2)**2 + BV3(3)**2)   !Shortest Vector from Node to Cylinder-Axis
                      ELSE
                        BV3(:) = 0.
                      END IF
                      IF (OV(1) .NE. 0.) THEN
                        dist2 = (BN(1) - BV3(1))/OV(1)
                      ELSE IF (OV(2) .NE. 0.) THEN
                        dist2 = (BN(2) - BV3(2))/OV(2)
                      ELSE IF (OV(3) .NE. 0.) THEN
                        dist2 = (BN(3) - BV3(3))/OV(3)
                      ELSE
                        dist2 = 0.
                      END IF
                      IF ((dist2 .LT. 1.-epsi) .AND. (dist2 .GT. epsi)) THEN
                        Species(iSpec)%Init(iInit)%ConstPress%ElemStat(ExamElem) = 2
                        Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside &
                             = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside + 1
                        TempElemPartlyInside(Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside) = ExamElem
                        EXIT
                      END IF
                    END IF
                  END DO
                END IF
              END DO
            END DO
            nInterOld = nInterest + 1
            nInterest = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside
          END DO
          
          
          
        END SELECT
        
        ALLOCATE (Species(iSpec)%Init(iInit)%ConstPress%ElemTotalInside(1:Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside+ &
             Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside))
        Species(iSpec)%Init(iInit)%ConstPress%ElemTotalInside(1:Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside) = &
             TempElemTotalInside(1:Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside)
        Species(iSpec)%Init(iInit)%ConstPress%ElemTotalInside(Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside + 1: &
             Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside + Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside) = &
             TempElemPartlyInside(1:Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside)
        Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside = Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside + &
             Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside
        IF (Species(iSpec)%Init(iInit)%ParticleEmissionType .EQ. 6) THEN
        ALLOCATE (Species(iSpec)%Init(iInit)%ConstPress%ConstPressureSamp(1:Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside &
               + Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside,6))
        END IF
        DEALLOCATE (TempElemTotalInside)
        DEALLOCATE (TempElemPartlyInside)
#ifdef MPI
        IF(Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside.NE.0)THEN
          WRITE (*,'(A49,I3.3,A2,I0,A3,I0)') 'Proc | Number of Elements inside ConstPressArea ',iInit,': ', &
               PMPIVAR%iProc,' | ',Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside
        END IF
#else
        WRITE (*,'(A49,I3.3,A2,I0)') 'Number of Elements inside ConstPressArea ',iInit,': ', &
               Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside
#endif
      END IF
    END DO
  END DO
END SUBROUTINE ParticlePressureCellIni

   
! SUBROUTINE ParticlePressureiniInsert(i)
!   
!   USE MOD_Particle_Vars
!   USE MOD_part_emission
! !--------------------------------------------------------------------------------------------------!
!   IMPLICIT NONE
! !--------------------------------------------------------------------------------------------------!
!   INTEGER, INTENT(IN)       :: i
! !--------------------------------------------------------------------------------------------------!
!   INTEGER          :: nPartInside, NbrOfParticle
!   REAL             :: TempInside, EkinInside, EInside
! !--------------------------------------------------------------------------------------------------!  
!        CALL ParticleInsideCheck(i, nPartInside, TempInside, EInside)
!        IF (Species(i)%ParticleEmission .GT. nPartInside) THEN
!         NbrOfParticle = Species(i)%ParticleEmission - nPartInside
!         WRITE(*,*) 'Emission PartNum (Spec ',i,')', NbrOfParticle
!         CALL SetParticlePosition(i,NbrOfParticle)
!         CALL SetParticleVelocity(i,NbrOfParticle)
!         CALL SetParticleChargeAndMass(i,NbrOfParticle)
!         IF (usevMPF) CALL SetParticleMPF(i,NbrOfParticle)
!         !IF (useDSMC) CALL SetParticleIntEnergy(i,NbrOfParticle)
!         PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfParticle
!         CALL UpdateNextFreePosition()
!        END IF
!    
! END SUBROUTINE ParticlePressureiniInsert

  
SUBROUTINE ParticlePressure (i, iInit, NbrOfParticle)
  USE MOD_Particle_Vars
#ifdef MPI
  USE MOD_part_MPI_Vars,      ONLY : PMPIVAR
#endif
!--------------------------------------------------------------------------------------------------!
  IMPLICIT NONE

#ifdef MPI
    INCLUDE 'mpif.h'                                                                               
#endif
  
  INTEGER, INTENT(IN)    :: i, iInit
  INTEGER, INTENT(OUT)   :: NbrOfParticle
!--------------------------------------------------------------------------------------------------!
  INTEGER          :: nPartInside
  REAL             :: TempInside, EkinInside, EInside
  REAL             :: TempComSend(2)
  REAL             :: TempComRec(2)
  INTEGER          :: IERROR

    CALL ParticleInsideCheck(i, iInit, nPartInside, TempInside, EInside)
#ifdef MPI
    TempComSend(1) = REAL(nPartInside)
    TempComSend(2) = EInside
    CALL MPI_ALLREDUCE(TempComSend, TempComRec, 2, MPI_DOUBLE_PRECISION, MPI_SUM, PMPIVAR%COMM, IERROR)
    nPartInside = INT(TempComRec(1))
    EInside = TempComRec(2)
#endif    
    IF (nPartInside .LT. Species(i)%Init(iInit)%ParticleEmission) THEN
      NbrOfParticle = INT(Species(i)%Init(iInit)%ParticleEmission) - nPartInside
      IF ((Species(i)%Init(iInit)%ConstPress%EkinInside - EInside).GT. 0.) THEN   !Both Factors should be bigger than 0
        Species(i)%Init(iInit)%MWTemperatureIC=2./3.*(Species(i)%Init(iInit)%ConstPress%EkinInside-EInside - &
                      NbrOfParticle*0.5*Species(i)%MassIC*Species(i)%Init(iInit)%VeloIC**2*Species(i)%MacroParticleFactor)/ &
                      (NbrOfParticle*Species(i)%MacroParticleFactor*BoltzmannConst)
      Species(i)%Init(iInit)%MWTemperatureIC=MAX(0.0, Species(i)%Init(iInit)%MWTemperatureIC)
      ELSE
        Species(i)%Init(iInit)%MWTemperatureIC = Species(i)%Init(iInit)%ConstPress%InitialTemp
      END IF
    ELSE
      NbrOfParticle = 0
    END IF
    
!    WRITE(*,*) 'Temperature  inside ConstPress:', TempInside
!    WRITE(*,*) 'EmissionTemperature:', Species(i)%MWTemperatureIC, 'nPartInside:', nPartInside
END SUBROUTINE ParticlePressure


SUBROUTINE ParticlePressureRem (i, iInit, NbrOfParticle)
! Removes all Particles in Pressure Area and sets NbrOfParticles to be inserted
  USE MOD_Particle_Vars
!--------------------------------------------------------------------------------------------------!
  IMPLICIT NONE
!--------------------------------------------------------------------------------------------------! 
  INTEGER, INTENT(IN)    :: i, iInit
  INTEGER, INTENT(OUT)   :: NbrOfParticle
!--------------------------------------------------------------------------------------------------!
  INTEGER          :: Particle
  REAL             :: BV1(3), BV2(3), BV3(3), OV(3), BN(3), vau(3), vauquad(3), TempVec(3)
  REAL             :: det1, det2, det3, dist1, dist2
!--------------------------------------------------------------------------------------------------

  SELECT CASE (TRIM(Species(i)%Init(iInit)%SpaceIC))
  CASE ('cuboid')
    BV1 = Species(i)%Init(iInit)%BaseVector1IC
    BV2 = Species(i)%Init(iInit)%BaseVector2IC
    OV  = Species(i)%Init(iInit)%ConstPress%OrthoVector
    DO Particle = 1,PDM%ParticleVecLength
      IF ((PartSpecies(Particle) .EQ. i) .AND. (PDM%ParticleInside(Particle))) THEN
        SELECT CASE (Species(i)%Init(iInit)%ConstPress%ElemStat(PEM%Element(Particle)))
        CASE (1)
          PDM%ParticleInside(Particle)=.false.
        CASE (2)
          BN(1:3) = PartState(Particle,1:3) - Species(i)%Init(iInit)%BasePointIC
          det1 = BN(1)*BV2(2)*OV(3) + BN(2)*BV2(3)*OV(1) + BN(3)*BV2(1)*OV(2) - &
                 BN(3)*BV2(2)*OV(1) - BN(1)*BV2(3)*OV(2) - BN(2)*BV2(1)*OV(3)
          det2 = BV1(1)*BN(2)*OV(3) + BV1(2)*BN(3)*OV(1) + BV1(3)*BN(1)*OV(2) - &
                 BV1(3)*BN(2)*OV(1) - BV1(1)*BN(3)*OV(2) - BV1(2)*BN(1)*OV(3)
          det3 = BV1(1)*BV2(2)*BN(3) + BV1(2)*BV2(3)*BN(1) + BV1(3)*BV2(1)*BN(2) - &
                 BV1(3)*BV2(2)*BN(1) - BV1(1)*BV2(3)*BN(2) - BV1(2)*BV2(1)*BN(3)
      
          det1 = det1/Species(i)%Init(iInit)%ConstPress%Determinant
          det2 = det2/Species(i)%Init(iInit)%ConstPress%Determinant
          det3 = det3/Species(i)%Init(iInit)%ConstPress%Determinant
          
          IF (((det1-0.5)**2 .LE. 0.25).AND.((det2-0.5)**2 .LE. 0.25).AND.((det3-0.5)**2 .LE. 0.25)) THEN
            PDM%ParticleInside(Particle)=.false.
          END IF
        END SELECT
      END IF
    END DO
  CASE ('cylinder')
    OV  = Species(i)%Init(iInit)%ConstPress%OrthoVector
    DO Particle = 1,PDM%ParticleVecLength
      IF ((PartSpecies(Particle) .EQ. i) .AND. (PDM%ParticleInside(Particle))) THEN
        SELECT CASE (Species(i)%Init(iInit)%ConstPress%ElemStat(PEM%Element(Particle)))
        CASE (1)
          PDM%ParticleInside(Particle)=.false.
        CASE (2)
          BN(1:3) = PartState(Particle,1:3) - Species(i)%Init(iInit)%BasePointIC
          BV2(1) = BN(2) * OV(3) - BN(3) * OV(2)                   !Vector orthogonal on BN and OrthoVector
          BV2(2) = BN(3) * OV(1) - BN(1) * OV(3)
          BV2(3) = BN(1) * OV(2) - BN(2) * OV(1)
          dist1  = SQRT((BV2(1)**2 + BV2(2)**2 + BV2(3)**2)/(OV(1)**2 + OV(2)**2 + OV(3)**2))
      
          IF (dist1 .LE. Species(i)%Init(iInit)%RadiusIC) THEN
            BV3(1) = OV(2) * BV2(3) - OV(3) * BV2(2)
            BV3(2) = OV(3) * BV2(1) - OV(1) * BV2(3)
            BV3(3) = OV(1) * BV2(2) - OV(2) * BV2(1)
            IF (BV3(1)**2 + BV3(2)**2 + BV3(3)**2 .NE. 0.) THEN
              BV3    = dist1 * BV3/SQRT(BV3(1)**2 + BV3(2)**2 + BV3(3)**2)   !Shortest Vector from Node to Cylinder-Axis
            ELSE
              BV3(:) = 0.
            END IF

            IF (OV(1) .NE. 0.) THEN
              dist2 = (BN(1) - BV3(1))/OV(1)
            ELSE IF (OV(2) .NE. 0.) THEN
              dist2 = (BN(2) - BV3(2))/OV(2)
            ELSE IF (OV(3) .NE. 0.) THEN
              dist2 = (BN(3) - BV3(3))/OV(3)
            ELSE
              dist2 = 0.
            END IF
            IF ((dist2 .LT. 1.) .AND. (dist2 .GT. 0.)) THEN
              PDM%ParticleInside(Particle)=.false.
            END IF
          END IF
        END SELECT
      END IF
    END DO
  END SELECT
  NbrOfParticle = INT(Species(i)%Init(iInit)%ParticleEmission)
END SUBROUTINE ParticlePressureRem


SUBROUTINE ParticleInsideCheck(i, iInit, nPartInside, TempInside, EkinInside)
  USE MOD_Particle_Vars
!--------------------------------------------------------------------------------------------------!
  IMPLICIT NONE
!--------------------------------------------------------------------------------------------------!
  INTEGER          :: nPartInside, Particle, i, iInit
  REAL             :: BV1(3), BV2(3), BV3(3), OV(3), BN(3), vau(3), vauquad(3), TempVec(3)
  REAL             :: det1, det2, det3, dist1, dist2, TempInside, EkinInside
!--------------------------------------------------------------------------------------------------!  
  
  nPartInside = 0
  vau(:) = 0
  vauquad(:) = 0
  SELECT CASE (TRIM(Species(i)%Init(iInit)%SpaceIC))
  CASE ('cuboid')
    BV1 = Species(i)%Init(iInit)%BaseVector1IC
    BV2 = Species(i)%Init(iInit)%BaseVector2IC
    OV  = Species(i)%Init(iInit)%ConstPress%OrthoVector
    DO Particle = 1,PDM%ParticleVecLength
      IF ((PartSpecies(Particle) .EQ. i) .AND. (PDM%ParticleInside(Particle))) THEN
        SELECT CASE (Species(i)%Init(iInit)%ConstPress%ElemStat(PEM%Element(Particle)))
        CASE (1)
          nPartInside = nPartInside + 1
          vau(1) = vau(1) + PartState(Particle,4)
          vau(2) = vau(2) + PartState(Particle,5)
          vau(3) = vau(3) + PartState(Particle,6)
          vauquad(1) = vauquad(1) + PartState(Particle,4)**2
          vauquad(2) = vauquad(2) + PartState(Particle,5)**2
          vauquad(3) = vauquad(3) + PartState(Particle,6)**2
        CASE (2)
          BN(1:3) = PartState(Particle,1:3) - Species(i)%Init(iInit)%BasePointIC
          det1 = BN(1)*BV2(2)*OV(3) + BN(2)*BV2(3)*OV(1) + BN(3)*BV2(1)*OV(2) - &
                 BN(3)*BV2(2)*OV(1) - BN(1)*BV2(3)*OV(2) - BN(2)*BV2(1)*OV(3)
          det2 = BV1(1)*BN(2)*OV(3) + BV1(2)*BN(3)*OV(1) + BV1(3)*BN(1)*OV(2) - &
                 BV1(3)*BN(2)*OV(1) - BV1(1)*BN(3)*OV(2) - BV1(2)*BN(1)*OV(3)
          det3 = BV1(1)*BV2(2)*BN(3) + BV1(2)*BV2(3)*BN(1) + BV1(3)*BV2(1)*BN(2) - &
                 BV1(3)*BV2(2)*BN(1) - BV1(1)*BV2(3)*BN(2) - BV1(2)*BV2(1)*BN(3)
      
          det1 = det1/Species(i)%Init(iInit)%ConstPress%Determinant
          det2 = det2/Species(i)%Init(iInit)%ConstPress%Determinant
          det3 = det3/Species(i)%Init(iInit)%ConstPress%Determinant
          
          IF (((det1-0.5)**2 .LE. 0.25).AND.((det2-0.5)**2 .LE. 0.25).AND.((det3-0.5)**2 .LE. 0.25)) THEN
            nPartInside = nPartInside + 1
            vau(1) = vau(1) + PartState(Particle,4)
            vau(2) = vau(2) + PartState(Particle,5)
            vau(3) = vau(3) + PartState(Particle,6)
            vauquad(1) = vauquad(1) + PartState(Particle,4)**2
            vauquad(2) = vauquad(2) + PartState(Particle,5)**2
            vauquad(3) = vauquad(3) + PartState(Particle,6)**2
          END IF
        CASE (3)
        END SELECT
      END IF
    END DO
  
  CASE ('cylinder')
    OV  = Species(i)%Init(iInit)%ConstPress%OrthoVector
    DO Particle = 1,PDM%ParticleVecLength
      IF ((PartSpecies(Particle) .EQ. i) .AND. (PDM%ParticleInside(Particle))) THEN
        SELECT CASE (Species(i)%Init(iInit)%ConstPress%ElemStat(PEM%Element(Particle)))
        CASE (1)
          nPartInside = nPartInside + 1
          vau = vau + PartState(Particle,4:6)
          vauquad = vauquad + PartState(Particle,4:6)**2
        CASE (2)
          BN(1:3) = PartState(Particle,1:3) - Species(i)%Init(iInit)%BasePointIC
          BV2(1) = BN(2) * OV(3) - BN(3) * OV(2)                   !Vector orthogonal on BN and OrthoVector
          BV2(2) = BN(3) * OV(1) - BN(1) * OV(3)
          BV2(3) = BN(1) * OV(2) - BN(2) * OV(1)
          dist1  = SQRT((BV2(1)**2 + BV2(2)**2 + BV2(3)**2)/(OV(1)**2 + OV(2)**2 + OV(3)**2))
      
          IF (dist1 .LE. Species(i)%Init(iInit)%RadiusIC) THEN
            BV3(1) = OV(2) * BV2(3) - OV(3) * BV2(2)
            BV3(2) = OV(3) * BV2(1) - OV(1) * BV2(3)
            BV3(3) = OV(1) * BV2(2) - OV(2) * BV2(1)
            IF (BV3(1)**2 + BV3(2)**2 + BV3(3)**2 .NE. 0.) THEN
              BV3    = dist1 * BV3/SQRT(BV3(1)**2 + BV3(2)**2 + BV3(3)**2)   !Shortest Vector from Node to Cylinder-Axis
            ELSE
              BV3(:) = 0.
            END IF

            IF (OV(1) .NE. 0.) THEN
              dist2 = (BN(1) - BV3(1))/OV(1)
            ELSE IF (OV(2) .NE. 0.) THEN
              dist2 = (BN(2) - BV3(2))/OV(2)
            ELSE IF (OV(3) .NE. 0.) THEN
              dist2 = (BN(3) - BV3(3))/OV(3)
            ELSE
              dist2 = 0.
            END IF
            IF ((dist2 .LT. 1.) .AND. (dist2 .GT. 0.)) THEN
              nPartInside = nPartInside + 1
              vau = vau + PartState(Particle,4:6)
              vauquad = vauquad + PartState(Particle,4:6)**2
            END IF
          END IF
        END SELECT
      END IF
    END DO
  END SELECT
  vau = vau/nPartInside
  EkinInside = 0.5*Species(i)%MassIC*(vauquad(1) + vauquad(2) +vauquad(3))*Species(i)%MacroParticleFactor
  vauquad = vauquad/nPartInside
  TempVec(1) = Species(i)%MassIC/BoltzmannConst*(vauquad(1) - vau(1)**2)
  TempVec(2) = Species(i)%MassIC/BoltzmannConst*(vauquad(2) - vau(2)**2)
  TempVec(3) = Species(i)%MassIC/BoltzmannConst*(vauquad(3) - vau(3)**2)
  TempInside = (TempVec(1) + TempVec(2) + TempVec(3))/3

END SUBROUTINE ParticleInsideCheck

SUBROUTINE PointInsideQuad3D(iSpec,iInit,Element,InElementCheck,dete)                                      !
  !DEC$ ATTRIBUTES FORCEINLINE :: ParticleInsideQuad3D
    USE MOD_Particle_Vars
    USE MOD_Mesh_Vars,     ONLY : ElemToSide
  !--------------------------------------------------------------------------------------------------!
    IMPLICIT NONE                                                                                    !
  !--------------------------------------------------------------------------------------------------!
  ! argument list declaration                                                                        !
    INTEGER                          :: iSpec, iInit, Element, kk, CellX, CellY, CellZ                      !
    LOGICAL                          :: InElementCheck                                               !
    REAL                             :: dete(6,2)                                                       !
  ! Local variable declaration                                                                       !
    INTEGER                          :: iLocSide, NodeNum
    LOGICAL                          :: PosCheck, NegCheck                                           !
    REAL                             :: A(1:3,1:4), cross(3)
  !--------------------------------------------------------------------------------------------------!
                                        !
  !--------------------------------------------------------------------------------------------------!


  InElementCheck = .FALSE.

  IF ( (Species(iSpec)%Init(iInit)%BasePointIC(1).LT.GEO%xmin).OR.(Species(iSpec)%Init(iInit)%BasePointIC(1).GT.GEO%xmax).OR. &
       (Species(iSpec)%Init(iInit)%BasePointIC(2).LT.GEO%ymin).OR.(Species(iSpec)%Init(iInit)%BasePointIC(2).GT.GEO%ymax).OR. &
       (Species(iSpec)%Init(iInit)%BasePointIC(3).LT.GEO%zmin).OR.(Species(iSpec)%Init(iInit)%BasePointIC(3).GT.GEO%zmax)) THEN
    RETURN
  END IF

  !--- get background mesh cell of Basepoint
  CellX = INT((Species(iSpec)%Init(iInit)%BasePointIC(1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1 
  CellX = MIN(GEO%FIBGMimax,CellX)                             
  CellY = INT((Species(iSpec)%Init(iInit)%BasePointIC(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
  CellY = MIN(GEO%FIBGMkmax,CellY) 
  CellZ = INT((Species(iSpec)%Init(iInit)%BasePointIC(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
  CellZ = MIN(GEO%FIBGMlmax,CellZ)

  !--- check all cells associated with this beckground mesh cell
  DO kk = 1, GEO%FIBGM(CellX,CellY,CellZ)%nElem
    Element = GEO%FIBGM(CellX,CellY,CellZ)%Element(kk)
    InElementCheck = .TRUE.
    DO iLocSide = 1,6                 ! for all 6 sides of the element
      !--- initialize flags for side checks
      PosCheck = .FALSE.
      NegCheck = .FALSE.
      !--- A = vector from particle to node coords
      DO NodeNum = 1,4
        A(:,NodeNum) = GEO%NodeCoords(:,GEO%ElemSideNodeID(NodeNum,iLocSide,Element)) - Species(iSpec)%Init(iInit)%BasePointIC(1:3)
      END DO

      !--- compute cross product for vector 1 and 3
      cross(1) = A(2,1) * A(3,3) - A(3,1) * A(2,3)
      cross(2) = A(3,1) * A(1,3) - A(1,1) * A(3,3)
      cross(3) = A(1,1) * A(2,3) - A(2,1) * A(1,3)

      !--- negative determinant of triangle 1 (points 1,3,2):
      dete(iLocSide,1) = cross(1) * A(1,2) + &
                        cross(2) * A(2,2) + &
                        cross(3) * A(3,2)
      dete(iLocSide,1) = -dete(iLocSide,1)
      !--- determinant of triangle 2 (points 1,3,4):
      dete(iLocSide,2) = cross(1) * A(1,4) + &
                        cross(2) * A(2,4) + &
                        cross(3) * A(3,4)
      IF (dete(iLocSide,1).lt.0) THEN
        NegCheck = .TRUE.
      ELSE
        PosCheck = .TRUE.
      END IF
      IF (dete(iLocSide,2).lt.0) THEN
        NegCheck = .TRUE.
      ELSE
        PosCheck = .TRUE.
      END IF

      !--- final determination whether point is in element
      IF (GEO%ConcaveElemSide(iLocSide,Element)) THEN
        IF (.NOT.PosCheck) InElementCheck = .FALSE.
      ELSE
        IF (NegCheck) InElementCheck = .FALSE.
      END IF
    END DO
    IF (InElementCheck) EXIT
  END DO
RETURN
END SUBROUTINE PointInsideQuad3D

END MODULE MOD_part_pressure
