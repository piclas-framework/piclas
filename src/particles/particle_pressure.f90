#include "boltzplatz.h"

MODULE MOD_part_pressure
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE
   
   PRIVATE
   PUBLIC       :: ParticlePressureIni, ParticlePressure, ParticleInsideCheck,ParticlePressureCellIni

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
  INTEGER              :: iShot, iElem, iNode, iSpec, iSide, i, n
  LOGICAL              :: ElementFound, InElementCheck
  
  
  epsi = 100.*epsilon(epsi)
  DO iSpec = 1,nSpecies
    Species(iSpec)%ConstPress%InitialTemp = Species(iSpec)%MWTemperatureIC
    !Species(iSpec)%ConstPress%OrthoVector(1) = 0
    IF (Species(iSpec)%ParticleEmissionType .EQ. 3) THEN
      ALLOCATE (TempElemTotalInside(nElems))
      ALLOCATE (TempElemPartlyInside(nElems))
      ALLOCATE (Species(iSpec)%ConstPress%ElemStat(nElems))
      
      Species(iSpec)%ConstPress%nElemTotalInside = 0
      Species(iSpec)%ConstPress%nElemPartlyInside = 0
      Species(iSpec)%ConstPress%ElemStat(1:nElems) = 3
      SELECT CASE (TRIM(Species(iSpec)%SpaceIC))
      CASE ('cuboid')
        
        BV1 = Species(iSpec)%BaseVector1IC
        BV2 = Species(iSpec)%BaseVector2IC

        Species(iSpec)%ConstPress%OrthoVector(1) = BV1(2) * BV2(3) - BV1(3) * BV2(2)
        Species(iSpec)%ConstPress%OrthoVector(2) = BV1(3) * BV2(1) - BV1(1) * BV2(3)
        Species(iSpec)%ConstPress%OrthoVector(3) = BV1(1) * BV2(2) - BV1(2) * BV2(1)
        OV = Species(iSpec)%ConstPress%OrthoVector
        
        IF ((OV(1) .EQ. 0) .AND. (OV(2) .EQ. 0) .AND. (OV(3) .EQ. 0)) THEN
          WRITE(*,*) 'Error in InitializeVariables: Cannot calculate Normal Vector of InitVolume in EmissionCase(3 or 4)!'
          STOP
        END IF
        
        Species(iSpec)%ConstPress%OrthoVector = OV*(Species(iSpec)%CuboidHeightIC/SQRT(OV(1)**2 + OV(2)**2 + OV(3)**2))
        OV = Species(iSpec)%ConstPress%OrthoVector
        Species(iSpec)%ConstPress%Determinant = ABS(BV1(1)*BV2(2)*OV(3) + BV1(2)*BV2(3)*OV(1) + BV1(3)*BV2(1)*OV(2) - &
                                                    BV1(3)*BV2(2)*OV(1) - BV1(1)*BV2(3)*OV(2) - BV1(2)*BV2(1)*OV(3))
        !WRITE (*,*) 'ConstPress', Species(iSpec)%ConstantPressure , 'Deter', Species(iSpec)%ConstPress%Determinant, 'ispec', &
        !iSpec &
                      !,'MPF',Species(iSpec)%MacroParticleFactor
        Species(iSpec)%ParticleEmission = Species(iSpec)%ConstantPressure * Species(iSpec)%ConstPress%Determinant / &
                            (BoltzmannConst * Species(iSpec)%MWTemperatureIC * Species(iSpec)%MacroParticleFactor)
        Species(iSpec)%ConstPress%EkinInside = 1.5*Species(iSpec)%ConstantPressure*Species(iSpec)%ConstPress%Determinant + &
                            INT(Species(iSpec)%ParticleEmission)*0.5*Species(iSpec)%MassIC*Species(iSpec)%VeloIC**2* &
        Species(iSpec)%MacroParticleFactor

        WRITE(*,*) 'Number of Particles in Constant-Pressure-Area:', Species(iSpec)%ParticleEmission
                        
        DO iElem = 1,nElems
          nNodesInside = 0
          nBoundNodes = 0
          DO iNode = 1,8
            BN(1:3) = GEO%NodeCoords(:,GEO%ElemToNodeID(iNode,iElem)) - Species(iSpec)%BasePointIC
            det1 = BN(1)*BV2(2)*OV(3) + BN(2)*BV2(3)*OV(1) + BN(3)*BV2(1)*OV(2) - &
                  BN(3)*BV2(2)*OV(1) - BN(1)*BV2(3)*OV(2) - BN(2)*BV2(1)*OV(3)
            det2 = BV1(1)*BN(2)*OV(3) + BV1(2)*BN(3)*OV(1) + BV1(3)*BN(1)*OV(2) - &
                  BV1(3)*BN(2)*OV(1) - BV1(1)*BN(3)*OV(2) - BV1(2)*BN(1)*OV(3)
            det3 = BV1(1)*BV2(2)*BN(3) + BV1(2)*BV2(3)*BN(1) + BV1(3)*BV2(1)*BN(2) - &
                  BV1(3)*BV2(2)*BN(1) - BV1(1)*BV2(3)*BN(2) - BV1(2)*BV2(1)*BN(3)
            
            det1 = det1/Species(iSpec)%ConstPress%Determinant
            det2 = det2/Species(iSpec)%ConstPress%Determinant
            det3 = det3/Species(iSpec)%ConstPress%Determinant

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
            Species(iSpec)%ConstPress%nElemTotalInside = Species(iSpec)%ConstPress%nElemTotalInside + 1
            TempElemTotalInside(Species(iSpec)%ConstPress%nElemTotalInside)=iElem
            Species(iSpec)%ConstPress%ElemStat(iELem) = 1
          ELSE IF ((nNodesInside .GE. 1) .AND. (nNodesInside .LE. 7)) THEN
            IF (nNodesInside + nBoundNodes .EQ. 8) THEN
              Species(iSpec)%ConstPress%nElemTotalInside = Species(iSpec)%ConstPress%nElemTotalInside + 1
              TempElemTotalInside(Species(iSpec)%ConstPress%nElemTotalInside)=iElem
              Species(iSpec)%ConstPress%ElemStat(iELem) = 1
            ELSE
              Species(iSpec)%ConstPress%nElemPartlyInside = Species(iSpec)%ConstPress%nElemPartlyInside + 1
              TempElemPartlyInside(Species(iSpec)%ConstPress%nElemPartlyInside) = iElem
              Species(iSpec)%ConstPress%ElemStat(iELem) = 2
            END IF
          ELSE
            IF (nBoundNodes .EQ. 8) THEN
              Species(iSpec)%ConstPress%nElemTotalInside = Species(iSpec)%ConstPress%nElemTotalInside + 1
              TempElemTotalInside(Species(iSpec)%ConstPress%nElemTotalInside)=iElem
              Species(iSpec)%ConstPress%ElemStat(iELem) = 1
            ELSE
              Species(iSpec)%ConstPress%ElemStat(iELem) = 3
            END IF
          END IF
        END DO
      
  ! If no Element has been found, search for Basepoint ========================================   
        IF ((Species(iSpec)%ConstPress%nElemTotalInside .EQ. 0) .AND. &
            (Species(iSpec)%ConstPress%nElemPartlyInside .EQ. 0)) THEN
          CALL PointInsideQuad3D(iSpec,Element,InElementCheck,dete)
          IF (InElementCheck) THEN
            Species(iSpec)%ConstPress%ElemStat(Element) = 2
            Species(iSpec)%ConstPress%nElemPartlyInside = 1
            TempElemPartlyInside(1) = Element
          END IF
        END IF

  ! Shoot on Sides of Neighbour-Elements of totally inside Elements ===========================
        
        DO iElem = 1,Species(iSpec)%ConstPress%nElemTotalInside
          Element = TempElemTotalInside(iElem)
          DO iSide = 1,6
            Side = ElemToSide(E2S_SIDE_ID, iSide, Element)
            IF (SideToElem(S2E_ELEM_ID,Side) .EQ. Element) THEN
              ExamElem = SideToElem(S2E_NB_ELEM_ID, Side)
            ELSE IF (SideToElem(S2E_NB_ELEM_ID,Side) .EQ. Element) THEN
              ExamElem = SideToElem(S2E_ELEM_ID,Side)
            END IF
            IF (Species(iSpec)%ConstPress%ElemStat(ExamElem) .EQ. 3) THEN
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
                BN = BN - Species(iSpec)%BasePointIC
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
                
                det1 = det1/Species(iSpec)%ConstPress%Determinant
                det2 = det2/Species(iSpec)%ConstPress%Determinant
                det3 = det3/Species(iSpec)%ConstPress%Determinant
                
                IF ((det1 .LT. 1.-epsi) .AND. (det1 .GT. epsi) .AND. (det2 .LT. 1-epsi) .AND. &
                    (det2 .GT. epsi) .AND. (det3 .LT. 1.-epsi) .AND. (det3 .GT. epsi)) THEN
                  Species(iSpec)%ConstPress%ElemStat(ExamElem) = 2
                  Species(iSpec)%ConstPress%nElemPartlyInside = Species(iSpec)%ConstPress%nElemPartlyInside + 1
                  TempElemPartlyInside(Species(iSpec)%ConstPress%nElemPartlyInside) = ExamElem
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
          IF (Species(iSpec)%ConstPress%ElemStat(Element) .EQ. 3) THEN
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
              BN = BN - Species(iSpec)%BasePointIC
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
              
              det1 = det1/Species(iSpec)%ConstPress%Determinant
              det2 = det2/Species(iSpec)%ConstPress%Determinant
              det3 = det3/Species(iSpec)%ConstPress%Determinant
                
              IF ((det1 .LT. 1.-epsi) .AND. (det1 .GT. epsi) .AND. (det2 .LT. 1-epsi) .AND. &
                  (det2 .GT. epsi) .AND. (det3 .LT. 1.-epsi) .AND. (det3 .GT. epsi)) THEN
                Species(iSpec)%ConstPress%ElemStat(Element) = 2
                Species(iSpec)%ConstPress%nElemPartlyInside = Species(iSpec)%ConstPress%nElemPartlyInside + 1
                TempElemPartlyInside(Species(iSpec)%ConstPress%nElemPartlyInside) = Element
                EXIT
              END IF
            END DO
          END IF
        END IF
      END DO
#endif

  ! Shoot on Sides of Neighbour-Elements of partly inside Elements ==========================================
      nInterest = Species(iSpec)%ConstPress%nElemPartlyInside
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

            IF (Species(iSpec)%ConstPress%ElemStat(ExamElem) .EQ. 3) THEN
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
                BN = BN - Species(iSpec)%BasePointIC
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
                
                det1 = det1/Species(iSpec)%ConstPress%Determinant
                det2 = det2/Species(iSpec)%ConstPress%Determinant
                det3 = det3/Species(iSpec)%ConstPress%Determinant
              
                IF (((det1-0.5)**2 .LT. 0.25-epsi).AND.((det2-0.5)**2 .LT. 0.25-epsi).AND.((det3-0.5)**2 .LT. 0.25-epsi)) THEN
                  Species(iSpec)%ConstPress%ElemStat(ExamElem) = 2
                  Species(iSpec)%ConstPress%nElemPartlyInside = Species(iSpec)%ConstPress%nElemPartlyInside + 1
                  TempElemPartlyInside(Species(iSpec)%ConstPress%nElemPartlyInside) = ExamElem
                  EXIT
                END IF
              END DO
            END IF
          END DO
        END DO
        nInterOld = nInterest + 1
        nInterest = Species(iSpec)%ConstPress%nElemPartlyInside
      END DO

        
      CASE ('cylinder')
  !       Species(iSpec)%ConstPress%BV1 = Species(iSpec)%CylinderHeightIC * Species(iSpec)%NormalIC
        
        Species(iSpec)%ConstPress%OrthoVector(1) = Species(iSpec)%BaseVector1IC(2) * Species(iSpec)%BaseVector2IC(3) - &
                            Species(iSpec)%BaseVector1IC(3) * Species(iSpec)%BaseVector2IC(2)
        Species(iSpec)%ConstPress%OrthoVector(2) = Species(iSpec)%BaseVector1IC(3) * Species(iSpec)%BaseVector2IC(1) - &
                            Species(iSpec)%BaseVector1IC(1) * Species(iSpec)%BaseVector2IC(3)
        Species(iSpec)%ConstPress%OrthoVector(3) = Species(iSpec)%BaseVector1IC(1) * Species(iSpec)%BaseVector2IC(2) - &
                            Species(iSpec)%BaseVector1IC(2) * Species(iSpec)%BaseVector2IC(1)
        OV(1:3) = Species(iSpec)%ConstPress%OrthoVector(1:3)
  !       WRITE(*,*) 'iSpec', iSpec                  
  !       WRITE(*,*) 'BaseVec1', Species(iSpec)%BaseVector1IC
  !       WRITE(*,*) 'BaseVec2', Species(iSpec)%BaseVector2IC
  !       WRITE(*,*) 'BV1', Species(iSpec)%ConstPress%OrthoVector
  !       
        IF ((OV(1) .EQ. 0) .AND. (OV(2) .EQ. 0) .AND. (OV(3) .EQ. 0)) THEN
          WRITE(*,*) 'Error in InitializeVariables: Cannot calculate NormalVector(Cyl) of InitVolume in EmissionCase(3 or 4)!'
          STOP
        END IF
        Species(iSpec)%ConstPress%OrthoVector = OV*(Species(iSpec)%CylinderHeightIC/SQRT(OV(1)**2 + OV(2)**2 + OV(3)**2))
        OV(1:3) = Species(iSpec)%ConstPress%OrthoVector(1:3)
        
        Species(iSpec)%ParticleEmission = Species(iSpec)%ConstantPressure * Species(iSpec)%RadiusIC**2 * &
                                          3.1415926535 * Species(iSpec)%CylinderHeightIC / (BoltzmannConst * &
                                          Species(iSpec)%MWTemperatureIC * Species(iSpec)%MacroParticleFactor)
        Species(iSpec)%ConstPress%EkinInside = 1.5*Species(iSpec)%ConstantPressure*Species(iSpec)%RadiusIC**2 * &
                                              3.1415926535 * Species(iSpec)%CylinderHeightIC + &
                                  INT(Species(iSpec)%ParticleEmission)*0.5*Species(iSpec)%MassIC*Species(iSpec)%VeloIC**2 &
* Species(iSpec)%MacroParticleFactor
        WRITE(*,*) 'Number of Particles in Constant-Pressure-Area:', Species(iSpec)%ParticleEmission
        
        DO iElem = 1,nElems
          nNodesInside = 0
          nBoundNodes = 0
          DO iNode = 1,8
            BN  = GEO%NodeCoords(:, GEO%ElemToNodeID(iNode,iElem)) - Species(iSpec)%BasePointIC
            BV2(1) = BN(2) * OV(3) - BN(3) * OV(2)                   !Vector orthogonal on BN and NormalIC
            BV2(2) = BN(3) * OV(1) - BN(1) * OV(3)
            BV2(3) = BN(1) * OV(2) - BN(2) * OV(1)
            dist1 = SQRT((BV2(1)**2 + BV2(2)**2 + BV2(3)**2)/(OV(1)**2 + OV(2)**2 + OV(3)**2))

            IF (dist1 .LE. Species(iSpec)%RadiusIC + epsi) THEN
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
              ELSE IF ((dist1 .LE. Species(iSpec)%RadiusIC+epsi) .AND. (dist1 .GE. Species(iSpec)%RadiusIC-epsi) &
                      .AND. (dist2 .GT. 0.) .AND. (dist2 .LT. 1.)) THEN
                nBoundNodes = nBoundNodes + 1
              ELSE IF ((dist2 .LT. 1.) .AND. (dist2 .GT. 0.)) THEN
                nNodesInside = nNodesInside + 1
              END IF
            END IF
          END DO
          
          IF (nNodesInside .EQ. 8) THEN
            Species(iSpec)%ConstPress%nElemTotalInside = Species(iSpec)%ConstPress%nElemTotalInside + 1
            TempElemTotalInside(Species(iSpec)%ConstPress%nElemTotalInside)=iElem
            Species(iSpec)%ConstPress%ElemStat(iELem) = 1
          ELSE IF ((nNodesInside .GE. 1) .AND. (nNodesInside .LE. 7)) THEN
            IF (nNodesInside + nBoundNodes .EQ. 8) THEN
              Species(iSpec)%ConstPress%nElemTotalInside = Species(iSpec)%ConstPress%nElemTotalInside + 1
              TempElemTotalInside(Species(iSpec)%ConstPress%nElemTotalInside)=iElem
              Species(iSpec)%ConstPress%ElemStat(iELem) = 1
            ELSE
              Species(iSpec)%ConstPress%nElemPartlyInside = Species(iSpec)%ConstPress%nElemPartlyInside + 1
              TempElemPartlyInside(Species(iSpec)%ConstPress%nElemPartlyInside) = iElem
              Species(iSpec)%ConstPress%ElemStat(iELem) = 2
            END IF
          ELSE
            IF (nBoundNodes .EQ. 8) THEN
              Species(iSpec)%ConstPress%nElemTotalInside = Species(iSpec)%ConstPress%nElemTotalInside + 1
              TempElemTotalInside(Species(iSpec)%ConstPress%nElemTotalInside)=iElem
              Species(iSpec)%ConstPress%ElemStat(iELem) = 1
            ELSE
              Species(iSpec)%ConstPress%ElemStat(iELem) = 3
            END IF
          END IF 
        END DO 
        
  ! If no Element has been found, search for Basepoint ========================================   
        IF ((Species(iSpec)%ConstPress%nElemTotalInside .EQ. 0) .AND. &
            (Species(iSpec)%ConstPress%nElemPartlyInside .EQ. 0)) THEN
          CALL PointInsideQuad3D(iSpec,Element,InElementCheck,dete)
          IF (InElementCheck) THEN
            Species(iSpec)%ConstPress%ElemStat(Element) = 2
            Species(iSpec)%ConstPress%nElemPartlyInside = 1
            TempElemPartlyInside(1) = Element
          END IF
        END IF
        
    !Schießen auf Nachbarzellen der totalen ===================================================   
        DO iElem = 1,Species(iSpec)%ConstPress%nElemTotalInside
          Element = TempElemTotalInside(iElem)
          DO iSide = 1,6
            Side = ElemToSide(E2S_SIDE_ID, iSide, Element)
            IF (SideToElem(S2E_ELEM_ID,Side) .EQ. Element) THEN
              ExamElem = SideToElem(S2E_NB_ELEM_ID, Side)
            ELSE IF (SideToElem(S2E_NB_ELEM_ID,Side) .EQ. Element) THEN
              ExamElem = SideToElem(S2E_ELEM_ID,Side)
            END IF
            IF (Species(iSpec)%ConstPress%ElemStat(ExamElem) .EQ. 3) THEN
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
                BN = BN - Species(iSpec)%BasePointIC
                !Lokalisierung
                BV2(1) = BN(2) * OV(3) - BN(3) * OV(2)                   !Vector orthogonal on BN and NormalIC
                BV2(2) = BN(3) * OV(1) - BN(1) * OV(3)
                BV2(3) = BN(1) * OV(2) - BN(2) * OV(1)
                dist1  = SQRT((BV2(1)**2 + BV2(2)**2 + BV2(3)**2)/ (OV(1)**2 + OV(2)**2 + OV(3)**2))

                IF (dist1 .LE. Species(iSpec)%RadiusIC + epsi) THEN
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
                    Species(iSpec)%ConstPress%ElemStat(ExamElem) = 2
                    Species(iSpec)%ConstPress%nElemPartlyInside = Species(iSpec)%ConstPress%nElemPartlyInside + 1
                    TempElemPartlyInside(Species(iSpec)%ConstPress%nElemPartlyInside) = ExamElem
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
          IF (Species(iSpec)%ConstPress%ElemStat(Element) .EQ. 3) THEN
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
              BN = BN - Species(iSpec)%BasePointIC
              !Lokalisierung
              BV2(1) = BN(2) * OV(3) - BN(3) * OV(2)                   !Vector orthogonal on BN and NormalIC
              BV2(2) = BN(3) * OV(1) - BN(1) * OV(3)
              BV2(3) = BN(1) * OV(2) - BN(2) * OV(1)
              dist1 = SQRT((BV2(1)**2 + BV2(2)**2 + BV2(3)**2)/(OV(1)**2 + OV(2)**2 + OV(3)**2))
              IF (dist1 .LE. Species(iSpec)%RadiusIC + epsi) THEN
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
                  Species(iSpec)%ConstPress%ElemStat(Element) = 2
                  Species(iSpec)%ConstPress%nElemPartlyInside = Species(iSpec)%ConstPress%nElemPartlyInside + 1
                  TempElemPartlyInside(Species(iSpec)%ConstPress%nElemPartlyInside) = Element
                  EXIT
                END IF
              END IF
            END DO
          END IF
        END IF
      END DO
#endif

    !!!Schießen auf Nachbarelemente der teilweisen ===============================================
        nInterest = Species(iSpec)%ConstPress%nElemPartlyInside
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
              IF (Species(iSpec)%ConstPress%ElemStat(ExamElem) .EQ. 3) THEN
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
                  BN = BN - Species(iSpec)%BasePointIC
                  !Lokalisierung
                  BV2(1) = BN(2) * OV(3) - BN(3) * OV(2)                   !Vector orthogonal on BN and NormalIC
                  BV2(2) = BN(3) * OV(1) - BN(1) * OV(3)
                  BV2(3) = BN(1) * OV(2) - BN(2) * OV(1)
                  dist1 = SQRT((BV2(1)**2 + BV2(2)**2 + BV2(3)**2)/(OV(1)**2 + OV(2)**2 + OV(3)**2))
                  IF (dist1 .LE. Species(iSpec)%RadiusIC + epsi) THEN
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
                      Species(iSpec)%ConstPress%ElemStat(ExamElem) = 2
                      Species(iSpec)%ConstPress%nElemPartlyInside = Species(iSpec)%ConstPress%nElemPartlyInside + 1
                      TempElemPartlyInside(Species(iSpec)%ConstPress%nElemPartlyInside) = ExamElem
                      EXIT
                    END IF
                  END IF
                END DO
              END IF
            END DO
          END DO
          nInterOld = nInterest + 1
          nInterest = Species(iSpec)%ConstPress%nElemPartlyInside
        END DO



      END SELECT
      
      ALLOCATE (Species(iSpec)%ConstPress%ElemTotalInside(Species(iSpec)%ConstPress%nElemTotalInside))
      ALLOCATE (Species(iSpec)%ConstPress%ElemPartlyInside(Species(iSpec)%ConstPress%nElemPartlyInside))
      Species(iSpec)%ConstPress%ElemTotalInside(1:Species(iSpec)%ConstPress%nElemTotalInside) = &
                                                    TempElemTotalInside(1:Species(iSpec)%ConstPress%nElemTotalInside)
      Species(iSpec)%ConstPress%ElemPartlyInside(1:Species(iSpec)%ConstPress%nElemPartlyInside) = &
                                                    TempElemPartlyInside(1:Species(iSpec)%ConstPress%nElemPartlyInside)
      DEALLOCATE (TempElemTotalInside)
      DEALLOCATE (TempElemPartlyInside)

      WRITE (*,*) 'Number of Elements inside ConstPressArea:', Species(iSpec)%ConstPress%nElemTotalInside
      WRITE (*,*) 'Number of Elements partly inside ConstPressArea:', Species(iSpec)%ConstPress%nElemPartlyInside
    
    END IF
  END DO  
  
END SUBROUTINE ParticlePressureIni




SUBROUTINE ParticlePressureCellIni()

  USE MOD_Particle_Vars
  USE MOD_Globals
  USE MOD_Mesh_Vars,      ONLY : nElems,ElemToSide,SideToElem, nSides, nInnerSides, nBCSides

  IMPLICIT NONE
 
  REAL                 :: BN(3), BV1(3), BV2(3), BV3(3), SV1(3), SV2(3), OV(3), epsi
  REAL                 :: det1, det2, det3, RandVal(2), dist1, dist2, dete(6,2)
  INTEGER, ALLOCATABLE :: TempElemTotalInside(:), TempElemPartlyInside(:)
  INTEGER              :: nNodesInside, nBoundNodes, nInterest, nInterOld, Element, Side, ExamElem
  INTEGER              :: iShot, iElem, iNode, iSpec, iSide, i, n
  LOGICAL              :: ElementFound, InElementCheck
  
  
  epsi = 100.*epsilon(epsi)
  DO iSpec = 1,nSpecies
    Species(iSpec)%ConstPress%InitialTemp = Species(iSpec)%MWTemperatureIC
    IF (Species(iSpec)%ParticleEmissionType .EQ. 4) THEN
      ALLOCATE (TempElemTotalInside(nElems))
      ALLOCATE (TempElemPartlyInside(nElems))
      ALLOCATE (Species(iSpec)%ConstPress%ElemStat(nElems))
      
      Species(iSpec)%ConstPress%nElemTotalInside = 0
      Species(iSpec)%ConstPress%nElemPartlyInside = 0
      Species(iSpec)%ConstPress%ElemStat(1:nElems) = 3
      SELECT CASE (TRIM(Species(iSpec)%SpaceIC))
      CASE ('cuboid')
        BV1 = Species(iSpec)%BaseVector1IC
        BV2 = Species(iSpec)%BaseVector2IC

        Species(iSpec)%ConstPress%OrthoVector(1) = BV1(2) * BV2(3) - BV1(3) * BV2(2)
        Species(iSpec)%ConstPress%OrthoVector(2) = BV1(3) * BV2(1) - BV1(1) * BV2(3)
        Species(iSpec)%ConstPress%OrthoVector(3) = BV1(1) * BV2(2) - BV1(2) * BV2(1)
        OV = Species(iSpec)%ConstPress%OrthoVector
        
        IF ((OV(1) .EQ. 0) .AND. (OV(2) .EQ. 0) .AND. (OV(3) .EQ. 0)) THEN
          WRITE(*,*) 'Error in InitializeVariables: Cannot calculate Normal Vector of InitVolume in EmissionCase(3 or 4)!'
          STOP
        END IF
        
        Species(iSpec)%ConstPress%OrthoVector = OV*(Species(iSpec)%CuboidHeightIC/SQRT(OV(1)**2 + OV(2)**2 + OV(3)**2))
        OV = Species(iSpec)%ConstPress%OrthoVector
        Species(iSpec)%ConstPress%Determinant = ABS(BV1(1)*BV2(2)*OV(3) + BV1(2)*BV2(3)*OV(1) + BV1(3)*BV2(1)*OV(2) - &
                                                    BV1(3)*BV2(2)*OV(1) - BV1(1)*BV2(3)*OV(2) - BV1(2)*BV2(1)*OV(3))
        !WRITE (*,*) 'ConstPress', Species(iSpec)%ConstantPressure , 'Deter', Species(iSpec)%ConstPress%Determinant, 'ispec', &
        !iSpec &
                      !,'MPF',Species(iSpec)%MacroParticleFactor

        ! ParticleEmission multiplied by cell volume equals particles to be in cells
        Species(iSpec)%ParticleEmission = Species(iSpec)%ConstantPressure/ &
                            (BoltzmannConst * Species(iSpec)%MWTemperatureIC * Species(iSpec)%MacroParticleFactor)
!        Species(iSpec)%ConstPress%EkinInside = 1.5*Species(iSpec)%ConstantPressure*Species(iSpec)%ConstPress%Determinant + &
!                            INT(Species(iSpec)%ParticleEmission)*0.5*Species(iSpec)%MassIC*Species(iSpec)%VeloIC**2* &
!        Species(iSpec)%MacroParticleFactor

!        WRITE(*,*) 'Number of Particles in Constant-Pressure-Area:', Species(iSpec)%ParticleEmission
                        
        DO iElem = 1,nElems
          nNodesInside = 0
          nBoundNodes = 0
          DO iNode = 1,8
            BN(1:3) = GEO%NodeCoords(:,GEO%ElemToNodeID(iNode,iElem)) - Species(iSpec)%BasePointIC
            det1 = BN(1)*BV2(2)*OV(3) + BN(2)*BV2(3)*OV(1) + BN(3)*BV2(1)*OV(2) - &
                  BN(3)*BV2(2)*OV(1) - BN(1)*BV2(3)*OV(2) - BN(2)*BV2(1)*OV(3)
            det2 = BV1(1)*BN(2)*OV(3) + BV1(2)*BN(3)*OV(1) + BV1(3)*BN(1)*OV(2) - &
                  BV1(3)*BN(2)*OV(1) - BV1(1)*BN(3)*OV(2) - BV1(2)*BN(1)*OV(3)
            det3 = BV1(1)*BV2(2)*BN(3) + BV1(2)*BV2(3)*BN(1) + BV1(3)*BV2(1)*BN(2) - &
                  BV1(3)*BV2(2)*BN(1) - BV1(1)*BV2(3)*BN(2) - BV1(2)*BV2(1)*BN(3)
            
            det1 = det1/Species(iSpec)%ConstPress%Determinant
            det2 = det2/Species(iSpec)%ConstPress%Determinant
            det3 = det3/Species(iSpec)%ConstPress%Determinant

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
            Species(iSpec)%ConstPress%nElemTotalInside = Species(iSpec)%ConstPress%nElemTotalInside + 1
            TempElemTotalInside(Species(iSpec)%ConstPress%nElemTotalInside)=iElem
            Species(iSpec)%ConstPress%ElemStat(iELem) = 1
          ELSE IF ((nNodesInside .GE. 1) .AND. (nNodesInside .LE. 7)) THEN
            IF (nNodesInside + nBoundNodes .EQ. 8) THEN
              Species(iSpec)%ConstPress%nElemTotalInside = Species(iSpec)%ConstPress%nElemTotalInside + 1
              TempElemTotalInside(Species(iSpec)%ConstPress%nElemTotalInside)=iElem
              Species(iSpec)%ConstPress%ElemStat(iELem) = 1
            ELSE
              Species(iSpec)%ConstPress%nElemPartlyInside = Species(iSpec)%ConstPress%nElemPartlyInside + 1
              TempElemPartlyInside(Species(iSpec)%ConstPress%nElemPartlyInside) = iElem
              Species(iSpec)%ConstPress%ElemStat(iELem) = 2
            END IF
          ELSE
            IF (nBoundNodes .EQ. 8) THEN
              Species(iSpec)%ConstPress%nElemTotalInside = Species(iSpec)%ConstPress%nElemTotalInside + 1
              TempElemTotalInside(Species(iSpec)%ConstPress%nElemTotalInside)=iElem
              Species(iSpec)%ConstPress%ElemStat(iELem) = 1
            ELSE
              Species(iSpec)%ConstPress%ElemStat(iELem) = 3
            END IF
          END IF
        END DO
      
  ! If no Element has been found, search for Basepoint ========================================   
        IF ((Species(iSpec)%ConstPress%nElemTotalInside .EQ. 0) .AND. &
            (Species(iSpec)%ConstPress%nElemPartlyInside .EQ. 0)) THEN
          CALL PointInsideQuad3D(iSpec,Element,InElementCheck,dete)
          IF (InElementCheck) THEN
            Species(iSpec)%ConstPress%ElemStat(Element) = 2
            Species(iSpec)%ConstPress%nElemPartlyInside = 1
            TempElemPartlyInside(1) = Element
          END IF
        END IF

  ! Shoot on Sides of Neighbour-Elements of totally inside Elements ===========================
        
        DO iElem = 1,Species(iSpec)%ConstPress%nElemTotalInside
          Element = TempElemTotalInside(iElem)
          DO iSide = 1,6
            Side = ElemToSide(E2S_SIDE_ID, iSide, Element)
            IF (SideToElem(S2E_ELEM_ID,Side) .EQ. Element) THEN
              ExamElem = SideToElem(S2E_NB_ELEM_ID, Side)
            ELSE IF (SideToElem(S2E_NB_ELEM_ID,Side) .EQ. Element) THEN
              ExamElem = SideToElem(S2E_ELEM_ID,Side)
            END IF
            IF (Species(iSpec)%ConstPress%ElemStat(ExamElem) .EQ. 3) THEN
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
                BN = BN - Species(iSpec)%BasePointIC
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
                
                det1 = det1/Species(iSpec)%ConstPress%Determinant
                det2 = det2/Species(iSpec)%ConstPress%Determinant
                det3 = det3/Species(iSpec)%ConstPress%Determinant
                
                IF ((det1 .LT. 1.-epsi) .AND. (det1 .GT. epsi) .AND. (det2 .LT. 1-epsi) .AND. &
                    (det2 .GT. epsi) .AND. (det3 .LT. 1.-epsi) .AND. (det3 .GT. epsi)) THEN
                  Species(iSpec)%ConstPress%ElemStat(ExamElem) = 2
                  Species(iSpec)%ConstPress%nElemPartlyInside = Species(iSpec)%ConstPress%nElemPartlyInside + 1
                  TempElemPartlyInside(Species(iSpec)%ConstPress%nElemPartlyInside) = ExamElem
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
          IF (Species(iSpec)%ConstPress%ElemStat(Element) .EQ. 3) THEN
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
              BN = BN - Species(iSpec)%BasePointIC
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
              
              det1 = det1/Species(iSpec)%ConstPress%Determinant
              det2 = det2/Species(iSpec)%ConstPress%Determinant
              det3 = det3/Species(iSpec)%ConstPress%Determinant
                
              IF ((det1 .LT. 1.-epsi) .AND. (det1 .GT. epsi) .AND. (det2 .LT. 1-epsi) .AND. &
                  (det2 .GT. epsi) .AND. (det3 .LT. 1.-epsi) .AND. (det3 .GT. epsi)) THEN
                Species(iSpec)%ConstPress%ElemStat(Element) = 2
                Species(iSpec)%ConstPress%nElemPartlyInside = Species(iSpec)%ConstPress%nElemPartlyInside + 1
                TempElemPartlyInside(Species(iSpec)%ConstPress%nElemPartlyInside) = Element
                EXIT
              END IF
            END DO
          END IF
        END IF
      END DO
#endif

  ! Shoot on Sides of Neighbour-Elements of partly inside Elements ==========================================
      nInterest = Species(iSpec)%ConstPress%nElemPartlyInside
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

            IF (Species(iSpec)%ConstPress%ElemStat(ExamElem) .EQ. 3) THEN
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
                BN = BN - Species(iSpec)%BasePointIC
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
                
                det1 = det1/Species(iSpec)%ConstPress%Determinant
                det2 = det2/Species(iSpec)%ConstPress%Determinant
                det3 = det3/Species(iSpec)%ConstPress%Determinant
              
                IF (((det1-0.5)**2 .LT. 0.25-epsi).AND.((det2-0.5)**2 .LT. 0.25-epsi).AND.((det3-0.5)**2 .LT. 0.25-epsi)) THEN
                  Species(iSpec)%ConstPress%ElemStat(ExamElem) = 2
                  Species(iSpec)%ConstPress%nElemPartlyInside = Species(iSpec)%ConstPress%nElemPartlyInside + 1
                  TempElemPartlyInside(Species(iSpec)%ConstPress%nElemPartlyInside) = ExamElem
                  EXIT
                END IF
              END DO
            END IF
          END DO
        END DO
        nInterOld = nInterest + 1
        nInterest = Species(iSpec)%ConstPress%nElemPartlyInside
      END DO

        
      CASE ('cylinder')
  !       Species(iSpec)%ConstPress%BV1 = Species(iSpec)%CylinderHeightIC * Species(iSpec)%NormalIC
        
        Species(iSpec)%ConstPress%OrthoVector(1) = Species(iSpec)%BaseVector1IC(2) * Species(iSpec)%BaseVector2IC(3) - &
                            Species(iSpec)%BaseVector1IC(3) * Species(iSpec)%BaseVector2IC(2)
        Species(iSpec)%ConstPress%OrthoVector(2) = Species(iSpec)%BaseVector1IC(3) * Species(iSpec)%BaseVector2IC(1) - &
                            Species(iSpec)%BaseVector1IC(1) * Species(iSpec)%BaseVector2IC(3)
        Species(iSpec)%ConstPress%OrthoVector(3) = Species(iSpec)%BaseVector1IC(1) * Species(iSpec)%BaseVector2IC(2) - &
                            Species(iSpec)%BaseVector1IC(2) * Species(iSpec)%BaseVector2IC(1)
        OV(1:3) = Species(iSpec)%ConstPress%OrthoVector(1:3)
  !       WRITE(*,*) 'iSpec', iSpec                  
  !       WRITE(*,*) 'BaseVec1', Species(iSpec)%BaseVector1IC
  !       WRITE(*,*) 'BaseVec2', Species(iSpec)%BaseVector2IC
  !       WRITE(*,*) 'BV1', Species(iSpec)%ConstPress%OrthoVector
  !       
        IF ((OV(1) .EQ. 0) .AND. (OV(2) .EQ. 0) .AND. (OV(3) .EQ. 0)) THEN
          WRITE(*,*) 'Error in InitializeVariables: Cannot calculate NormalVector(Cyl) of InitVolume in EmissionCase(3 or 4)!'
          STOP
        END IF
        Species(iSpec)%ConstPress%OrthoVector = OV*(Species(iSpec)%CylinderHeightIC/SQRT(OV(1)**2 + OV(2)**2 + OV(3)**2))
        OV(1:3) = Species(iSpec)%ConstPress%OrthoVector(1:3)
        
        Species(iSpec)%ParticleEmission = Species(iSpec)%ConstantPressure * Species(iSpec)%RadiusIC**2 * &
                                          3.1415926535 * Species(iSpec)%CylinderHeightIC / (BoltzmannConst * &
                                          Species(iSpec)%MWTemperatureIC * Species(iSpec)%MacroParticleFactor)
        Species(iSpec)%ConstPress%EkinInside = 1.5*Species(iSpec)%ConstantPressure*Species(iSpec)%RadiusIC**2 * &
                                              3.1415926535 * Species(iSpec)%CylinderHeightIC + &
                                  INT(Species(iSpec)%ParticleEmission)*0.5*Species(iSpec)%MassIC*Species(iSpec)%VeloIC**2 &
* Species(iSpec)%MacroParticleFactor
        WRITE(*,*) 'Number of Particles in Constant-Pressure-Area:', Species(iSpec)%ParticleEmission
        
        DO iElem = 1,nElems
          nNodesInside = 0
          nBoundNodes = 0
          DO iNode = 1,8
            BN  = GEO%NodeCoords(:, GEO%ElemToNodeID(iNode,iElem)) - Species(iSpec)%BasePointIC
            BV2(1) = BN(2) * OV(3) - BN(3) * OV(2)                   !Vector orthogonal on BN and NormalIC
            BV2(2) = BN(3) * OV(1) - BN(1) * OV(3)
            BV2(3) = BN(1) * OV(2) - BN(2) * OV(1)
            dist1 = SQRT((BV2(1)**2 + BV2(2)**2 + BV2(3)**2)/(OV(1)**2 + OV(2)**2 + OV(3)**2))

            IF (dist1 .LE. Species(iSpec)%RadiusIC + epsi) THEN
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
              ELSE IF ((dist1 .LE. Species(iSpec)%RadiusIC+epsi) .AND. (dist1 .GE. Species(iSpec)%RadiusIC-epsi) &
                      .AND. (dist2 .GT. 0.) .AND. (dist2 .LT. 1.)) THEN
                nBoundNodes = nBoundNodes + 1
              ELSE IF ((dist2 .LT. 1.) .AND. (dist2 .GT. 0.)) THEN
                nNodesInside = nNodesInside + 1
              END IF
            END IF
          END DO
          
          IF (nNodesInside .EQ. 8) THEN
            Species(iSpec)%ConstPress%nElemTotalInside = Species(iSpec)%ConstPress%nElemTotalInside + 1
            TempElemTotalInside(Species(iSpec)%ConstPress%nElemTotalInside)=iElem
            Species(iSpec)%ConstPress%ElemStat(iELem) = 1
          ELSE IF ((nNodesInside .GE. 1) .AND. (nNodesInside .LE. 7)) THEN
            IF (nNodesInside + nBoundNodes .EQ. 8) THEN
              Species(iSpec)%ConstPress%nElemTotalInside = Species(iSpec)%ConstPress%nElemTotalInside + 1
              TempElemTotalInside(Species(iSpec)%ConstPress%nElemTotalInside)=iElem
              Species(iSpec)%ConstPress%ElemStat(iELem) = 1
            ELSE
              Species(iSpec)%ConstPress%nElemPartlyInside = Species(iSpec)%ConstPress%nElemPartlyInside + 1
              TempElemPartlyInside(Species(iSpec)%ConstPress%nElemPartlyInside) = iElem
              Species(iSpec)%ConstPress%ElemStat(iELem) = 2
            END IF
          ELSE
            IF (nBoundNodes .EQ. 8) THEN
              Species(iSpec)%ConstPress%nElemTotalInside = Species(iSpec)%ConstPress%nElemTotalInside + 1
              TempElemTotalInside(Species(iSpec)%ConstPress%nElemTotalInside)=iElem
              Species(iSpec)%ConstPress%ElemStat(iELem) = 1
            ELSE
              Species(iSpec)%ConstPress%ElemStat(iELem) = 3
            END IF
          END IF 
        END DO 
        
  ! If no Element has been found, search for Basepoint ========================================   
        IF ((Species(iSpec)%ConstPress%nElemTotalInside .EQ. 0) .AND. &
            (Species(iSpec)%ConstPress%nElemPartlyInside .EQ. 0)) THEN
          CALL PointInsideQuad3D(iSpec,Element,InElementCheck,dete)
          IF (InElementCheck) THEN
            Species(iSpec)%ConstPress%ElemStat(Element) = 2
            Species(iSpec)%ConstPress%nElemPartlyInside = 1
            TempElemPartlyInside(1) = Element
          END IF
        END IF
        
    !Schießen auf Nachbarzellen der totalen ===================================================   
        DO iElem = 1,Species(iSpec)%ConstPress%nElemTotalInside
          Element = TempElemTotalInside(iElem)
          DO iSide = 1,6
            Side = ElemToSide(E2S_SIDE_ID, iSide, Element)
            IF (SideToElem(S2E_ELEM_ID,Side) .EQ. Element) THEN
              ExamElem = SideToElem(S2E_NB_ELEM_ID, Side)
            ELSE IF (SideToElem(S2E_NB_ELEM_ID,Side) .EQ. Element) THEN
              ExamElem = SideToElem(S2E_ELEM_ID,Side)
            END IF
            IF (Species(iSpec)%ConstPress%ElemStat(ExamElem) .EQ. 3) THEN
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
                BN = BN - Species(iSpec)%BasePointIC
                !Lokalisierung
                BV2(1) = BN(2) * OV(3) - BN(3) * OV(2)                   !Vector orthogonal on BN and NormalIC
                BV2(2) = BN(3) * OV(1) - BN(1) * OV(3)
                BV2(3) = BN(1) * OV(2) - BN(2) * OV(1)
                dist1  = SQRT((BV2(1)**2 + BV2(2)**2 + BV2(3)**2)/ (OV(1)**2 + OV(2)**2 + OV(3)**2))

                IF (dist1 .LE. Species(iSpec)%RadiusIC + epsi) THEN
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
                    Species(iSpec)%ConstPress%ElemStat(ExamElem) = 2
                    Species(iSpec)%ConstPress%nElemPartlyInside = Species(iSpec)%ConstPress%nElemPartlyInside + 1
                    TempElemPartlyInside(Species(iSpec)%ConstPress%nElemPartlyInside) = ExamElem
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
          IF (Species(iSpec)%ConstPress%ElemStat(Element) .EQ. 3) THEN
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
              BN = BN - Species(iSpec)%BasePointIC
              !Lokalisierung
              BV2(1) = BN(2) * OV(3) - BN(3) * OV(2)                   !Vector orthogonal on BN and NormalIC
              BV2(2) = BN(3) * OV(1) - BN(1) * OV(3)
              BV2(3) = BN(1) * OV(2) - BN(2) * OV(1)
              dist1 = SQRT((BV2(1)**2 + BV2(2)**2 + BV2(3)**2)/(OV(1)**2 + OV(2)**2 + OV(3)**2))
              IF (dist1 .LE. Species(iSpec)%RadiusIC + epsi) THEN
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
                  Species(iSpec)%ConstPress%ElemStat(Element) = 2
                  Species(iSpec)%ConstPress%nElemPartlyInside = Species(iSpec)%ConstPress%nElemPartlyInside + 1
                  TempElemPartlyInside(Species(iSpec)%ConstPress%nElemPartlyInside) = Element
                  EXIT
                END IF
              END IF
            END DO
          END IF
        END IF
      END DO
#endif

    !!!Schießen auf Nachbarelemente der teilweisen ===============================================
        nInterest = Species(iSpec)%ConstPress%nElemPartlyInside
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
              IF (Species(iSpec)%ConstPress%ElemStat(ExamElem) .EQ. 3) THEN
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
                  BN = BN - Species(iSpec)%BasePointIC
                  !Lokalisierung
                  BV2(1) = BN(2) * OV(3) - BN(3) * OV(2)                   !Vector orthogonal on BN and NormalIC
                  BV2(2) = BN(3) * OV(1) - BN(1) * OV(3)
                  BV2(3) = BN(1) * OV(2) - BN(2) * OV(1)
                  dist1 = SQRT((BV2(1)**2 + BV2(2)**2 + BV2(3)**2)/(OV(1)**2 + OV(2)**2 + OV(3)**2))
                  IF (dist1 .LE. Species(iSpec)%RadiusIC + epsi) THEN
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
                      Species(iSpec)%ConstPress%ElemStat(ExamElem) = 2
                      Species(iSpec)%ConstPress%nElemPartlyInside = Species(iSpec)%ConstPress%nElemPartlyInside + 1
                      TempElemPartlyInside(Species(iSpec)%ConstPress%nElemPartlyInside) = ExamElem
                      EXIT
                    END IF
                  END IF
                END DO
              END IF
            END DO
          END DO
          nInterOld = nInterest + 1
          nInterest = Species(iSpec)%ConstPress%nElemPartlyInside
        END DO



      END SELECT
      
      ALLOCATE (Species(iSpec)%ConstPress%ElemTotalInside(1:Species(iSpec)%ConstPress%nElemTotalInside+ &
           Species(iSpec)%ConstPress%nElemPartlyInside))
      Species(iSpec)%ConstPress%ElemTotalInside(1:Species(iSpec)%ConstPress%nElemTotalInside) = &
           TempElemTotalInside(1:Species(iSpec)%ConstPress%nElemTotalInside)
      Species(iSpec)%ConstPress%ElemTotalInside(Species(iSpec)%ConstPress%nElemTotalInside + 1: &
           Species(iSpec)%ConstPress%nElemPartlyInside + Species(iSpec)%ConstPress%nElemTotalInside) = &
           TempElemPartlyInside(1:Species(iSpec)%ConstPress%nElemPartlyInside)
      Species(iSpec)%ConstPress%nElemTotalInside = Species(iSpec)%ConstPress%nElemTotalInside + &
                                                   Species(iSpec)%ConstPress%nElemPartlyInside
      DEALLOCATE (TempElemTotalInside)
      DEALLOCATE (TempElemPartlyInside)

      WRITE (*,*) 'Number of Elements inside ConstPressArea:', Species(iSpec)%ConstPress%nElemTotalInside
    END IF
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

  
SUBROUTINE ParticlePressure (i, NbrOfParticle)
  USE MOD_Particle_Vars
#ifdef MPI
  USE MOD_part_MPI_Vars,      ONLY : PMPIVAR
#endif
!--------------------------------------------------------------------------------------------------!
  IMPLICIT NONE

#ifdef MPI
    INCLUDE 'mpif.h'                                                                               
#endif
  
  INTEGER, INTENT(IN)    :: i
  INTEGER, INTENT(OUT)   :: NbrOfParticle
!--------------------------------------------------------------------------------------------------!
  INTEGER          :: nPartInside
  REAL             :: TempInside, EkinInside, EInside
  REAL             :: TempComSend(2)
  REAL             :: TempComRec(2)
  INTEGER          :: IERROR

    CALL ParticleInsideCheck(i, nPartInside, TempInside, EInside)
#ifdef MPI
    TempComSend(1) = REAL(nPartInside)
    TempComSend(2) = EInside
    CALL MPI_ALLREDUCE(TempComSend, TempComRec, 2, MPI_DOUBLE_PRECISION, MPI_SUM, PMPIVAR%COMM, IERROR)
    nPartInside = INT(TempComRec(1))
    EInside = TempComRec(2)
#endif    
    IF (nPartInside .LT. Species(i)%ParticleEmission) THEN
      NbrOfParticle = INT(Species(i)%ParticleEmission) - nPartInside
      IF ((Species(i)%ConstPress%EkinInside - EInside).GT. 0.) THEN   !Both Factors should be bigger than 0
        Species(i)%MWTemperatureIC=2./3.*(Species(i)%ConstPress%EkinInside-EInside - &
                      NbrOfParticle*0.5*Species(i)%MassIC*Species(i)%VeloIC**2*Species(i)%MacroParticleFactor)/ &
                      (NbrOfParticle*Species(i)%MacroParticleFactor*BoltzmannConst)
      Species(i)%MWTemperatureIC=MAX(0.0, Species(i)%MWTemperatureIC)
      ELSE
        Species(i)%MWTemperatureIC = Species(i)%ConstPress%InitialTemp
      END IF
    ELSE
      NbrOfParticle = 0
    END IF
    
!    WRITE(*,*) 'Temperature  inside ConstPress:', TempInside
!    WRITE(*,*) 'EmissionTemperature:', Species(i)%MWTemperatureIC, 'nPartInside:', nPartInside
END SUBROUTINE ParticlePressure
  
SUBROUTINE ParticleInsideCheck(i, nPartInside, TempInside, EkinInside)
  USE MOD_Particle_Vars
!--------------------------------------------------------------------------------------------------!
  IMPLICIT NONE
!--------------------------------------------------------------------------------------------------!
  INTEGER          :: nPartInside, Particle, i 
  REAL             :: BV1(3), BV2(3), BV3(3), OV(3), BN(3), vau(3), vauquad(3), TempVec(3)
  REAL             :: det1, det2, det3, dist1, dist2, TempInside, EkinInside
!--------------------------------------------------------------------------------------------------!  
  
  nPartInside = 0
  vau(:) = 0
  vauquad(:) = 0
  SELECT CASE (TRIM(Species(i)%SpaceIC))
  CASE ('cuboid')
    BV1 = Species(i)%BaseVector1IC
    BV2 = Species(i)%BaseVector2IC
    OV  = Species(i)%ConstPress%OrthoVector
    DO Particle = 1,PDM%ParticleVecLength
      IF ((PartSpecies(Particle) .EQ. i) .AND. (PDM%ParticleInside(Particle))) THEN
        SELECT CASE (Species(i)%ConstPress%ElemStat(PEM%Element(Particle)))
        CASE (1)
          nPartInside = nPartInside + 1
          vau(1) = vau(1) + PartState(Particle,4)
          vau(2) = vau(2) + PartState(Particle,5)
          vau(3) = vau(3) + PartState(Particle,6)
          vauquad(1) = vauquad(1) + PartState(Particle,4)**2
          vauquad(2) = vauquad(2) + PartState(Particle,5)**2
          vauquad(3) = vauquad(3) + PartState(Particle,6)**2
        CASE (2)
          BN(1:3) = PartState(Particle,1:3) - Species(i)%BasePointIC
          det1 = BN(1)*BV2(2)*OV(3) + BN(2)*BV2(3)*OV(1) + BN(3)*BV2(1)*OV(2) - &
                 BN(3)*BV2(2)*OV(1) - BN(1)*BV2(3)*OV(2) - BN(2)*BV2(1)*OV(3)
          det2 = BV1(1)*BN(2)*OV(3) + BV1(2)*BN(3)*OV(1) + BV1(3)*BN(1)*OV(2) - &
                 BV1(3)*BN(2)*OV(1) - BV1(1)*BN(3)*OV(2) - BV1(2)*BN(1)*OV(3)
          det3 = BV1(1)*BV2(2)*BN(3) + BV1(2)*BV2(3)*BN(1) + BV1(3)*BV2(1)*BN(2) - &
                 BV1(3)*BV2(2)*BN(1) - BV1(1)*BV2(3)*BN(2) - BV1(2)*BV2(1)*BN(3)
      
          det1 = det1/Species(i)%ConstPress%Determinant
          det2 = det2/Species(i)%ConstPress%Determinant
          det3 = det3/Species(i)%ConstPress%Determinant
          
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
    OV  = Species(i)%ConstPress%OrthoVector
    DO Particle = 1,PDM%ParticleVecLength
      IF ((PartSpecies(Particle) .EQ. i) .AND. (PDM%ParticleInside(Particle))) THEN
        SELECT CASE (Species(i)%ConstPress%ElemStat(PEM%Element(Particle)))
        CASE (1)
          nPartInside = nPartInside + 1
          vau = vau + PartState(Particle,4:6)
          vauquad = vauquad + PartState(Particle,4:6)**2
        CASE (2)
          BN(1:3) = PartState(Particle,1:3) - Species(i)%BasePointIC
          BV2(1) = BN(2) * OV(3) - BN(3) * OV(2)                   !Vector orthogonal on BN and OrthoVector
          BV2(2) = BN(3) * OV(1) - BN(1) * OV(3)
          BV2(3) = BN(1) * OV(2) - BN(2) * OV(1)
          dist1  = SQRT((BV2(1)**2 + BV2(2)**2 + BV2(3)**2)/(OV(1)**2 + OV(2)**2 + OV(3)**2))
      
          IF (dist1 .LE. Species(i)%RadiusIC) THEN
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

SUBROUTINE PointInsideQuad3D(iSpec,Element,InElementCheck,dete)                                      !
  !DEC$ ATTRIBUTES FORCEINLINE :: ParticleInsideQuad3D
    USE MOD_Particle_Vars
    USE MOD_Mesh_Vars,     ONLY : ElemToSide
  !--------------------------------------------------------------------------------------------------!
    IMPLICIT NONE                                                                                    !
  !--------------------------------------------------------------------------------------------------!
  ! argument list declaration                                                                        !
    INTEGER                          :: iSpec, Element, kk, CellX, CellY, CellZ                      !
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

  IF ( (Species(iSpec)%BasePointIC(1).LT.GEO%xmin).OR.(Species(iSpec)%BasePointIC(1).GT.GEO%xmax).OR. &
       (Species(iSpec)%BasePointIC(2).LT.GEO%ymin).OR.(Species(iSpec)%BasePointIC(2).GT.GEO%ymax).OR. &
       (Species(iSpec)%BasePointIC(3).LT.GEO%zmin).OR.(Species(iSpec)%BasePointIC(3).GT.GEO%zmax)) THEN
    RETURN
  END IF

  !--- get background mesh cell of Basepoint
  CellX = INT((Species(iSpec)%BasePointIC(1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1 
  CellX = MIN(GEO%FIBGMimax,CellX)                             
  CellY = INT((Species(iSpec)%BasePointIC(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
  CellY = MIN(GEO%FIBGMkmax,CellY) 
  CellZ = INT((Species(iSpec)%BasePointIC(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
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
        A(:,NodeNum) = GEO%NodeCoords(:,GEO%ElemSideNodeID(NodeNum,iLocSide,Element)) - Species(iSpec)%BasePointIC(1:3)
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
