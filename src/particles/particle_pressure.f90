!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (gitlab.com/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"

MODULE MOD_part_pressure
!===================================================================================================================================
! Module for constant pressure emission types
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE ParticlePressureIni
  MODULE PROCEDURE ParticlePressureIni
END INTERFACE

INTERFACE ParticlePressure
  MODULE PROCEDURE ParticlePressure
END INTERFACE

INTERFACE ParticleInsideCheck
  MODULE PROCEDURE ParticleInsideCheck
END INTERFACE

INTERFACE ParticlePressureCellIni
  MODULE PROCEDURE ParticlePressureCellIni
END INTERFACE

INTERFACE ParticlePressureRem
  MODULE PROCEDURE ParticlePressureRem
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC       :: ParticlePressureIni, ParticlePressure, ParticleInsideCheck,ParticlePressureCellIni,ParticlePressureRem
!===================================================================================================================================

CONTAINS

SUBROUTINE ParticlePressureIni()
!===================================================================================================================================
! Initialization of constant pressure emission types
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_Preproc
  USE MOD_Globals_Vars,            ONLY:BoltzmannConst, Pi
  USE MOD_Particle_Vars
  USE MOD_Particle_Tracking_Vars,  ONLY:DoRefMapping
  USE MOD_Particle_Mesh_Vars,      ONLY:epsInCell
  USE MOD_Particle_Mesh,           ONLY:PointToExactElement
  USE MOD_Mesh_Vars,               ONLY:nElems,ElemToSide
  USE MOD_Mesh_Vars,               ONLY:NGeo,XCL_NGeo,XiCL_NGeo,wBaryCL_NGeo
  USE MOD_Eval_XYZ,                ONLY:TensorProductInterpolation
  USE MOD_Particle_Mesh_Vars,      ONLY:PartElemToElemAndSide
#if USE_MPI
  USE MOD_Mesh_Vars,               ONLY : nInnerSides, nBCSides
#endif
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                 :: BN(3), BV1(3), BV2(3), BV3(3), OV(3), epsi
  REAL                 :: det1, det2, det3, RandVal(2), dist1, dist2, xi(3)
  INTEGER, ALLOCATABLE :: TempElemTotalInside(:), TempElemPartlyInside(:)
  INTEGER              :: nNodesInside, nBoundNodes, nInterest, nInterOld, Element, ExamElem
  INTEGER              :: iShot, iElem, iSpec, iInit, iElem2
  INTEGER              :: iLocSide,SideID,locSideID
  INTEGER              :: i,j,k
  LOGICAL              :: InElementCheck,Marked
!===================================================================================================================================

  IF(.NOT.DoRefMapping) CALL abort(&
__STAMP__&
,' Particle pressure is only possible with tracking via reference element mapping!')


   CALL abort(&
__STAMP__&
,' GEO has to be exchanged by XCL_NGeo!' )
  !epsi = 100.*epsilon(epsi)
  ! accuracy of mapping into ref element
  epsi = epsInCell
  DO iSpec = 1,nSpecies
    DO iInit = Species(iSpec)%StartnumberOfInits, Species(iSpec)%NumberOfInits
      Species(iSpec)%Init(iInit)%ConstPress%InitialTemp = Species(iSpec)%Init(iInit)%MWTemperatureIC
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
            CALL Abort(&
__STAMP__&
,'Error in InitializeVariables: Cannot calculate Normal Vector of InitVolume in EmissionCase(3 or 4)!')
          END IF

          Species(iSpec)%Init(iInit)%ConstPress%OrthoVector &
               = OV*(Species(iSpec)%Init(iInit)%CuboidHeightIC/SQRT(OV(1)**2 + OV(2)**2 + OV(3)**2))
          OV   = Species(iSpec)%Init(iInit)%ConstPress%OrthoVector
          Species(iSpec)%Init(iInit)%ConstPress%Determinant = ABS(BV1(1)*BV2(2)*OV(3) + BV1(2)*BV2(3)*OV(1) + BV1(3)*BV2(1)*OV(2) &
               - BV1(3)*BV2(2)*OV(1) - BV1(1)*BV2(3)*OV(2) - BV1(2)*BV2(1)*OV(3))
          Species(iSpec)%Init(iInit)%ParticleEmission &
               = Species(iSpec)%Init(iInit)%ConstantPressure * Species(iSpec)%Init(iInit)%ConstPress%Determinant / &
               (BoltzmannConst * Species(iSpec)%Init(iInit)%MWTemperatureIC * Species(iSpec)%MacroParticleFactor)
          Species(iSpec)%Init(iInit)%ConstPress%EkinInside &
               = 1.5*Species(iSpec)%Init(iInit)%ConstantPressure*Species(iSpec)%Init(iInit)%ConstPress%Determinant + &
               INT(Species(iSpec)%Init(iInit)%ParticleEmission)*0.5*Species(iSpec)%MassIC*Species(iSpec)%Init(iInit)%VeloIC**2* &
               Species(iSpec)%MacroParticleFactor

          IPWRITE (UNIT_StdOut,'(I4,A49,I3.3,A2,g0)') 'Number of Particles inside ConstPressArea ',iInit, ': ', &
               Species(iSpec)%Init(iInit)%ParticleEmission

          DO iElem = 1,nElems
            nNodesInside = 0
            nBoundNodes = 0
            ! loop over all 8 nodes
            DO k=0,NGeo,NGeo
              DO j=0,NGeo,NGeo
                DO i=0,NGeo,NGeo
                  BN(1:3) = XCL_NGeo(1:3,i,j,k,iElem) - Species(iSpec)%Init(iInit)%BasePointIC
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
                 END DO ! i=0,NGeo
              END DO ! j=0,NGeo
            END DO ! k=0,NGeo

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
          END DO ! iElem

          ! If no Element has been found, search for Basepoint ========================================
          IF ((Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside .EQ. 0) .AND. &
               (Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside .EQ. 0)) THEN
            CALL PointToExactElement(Species(iSpec)%Init(iInit)%BasePointIC,Element,InElementCheck,doHalo=.FALSE.)
            IF (InElementCheck) THEN
              Species(iSpec)%Init(iInit)%ConstPress%ElemStat(Element) = 2
              Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside = 1
              TempElemPartlyInside(1) = Element
            END IF
          END IF

          ! Shoot on Sides of Neighbour-Elements of totally inside Elements ===========================
          DO iElem = 1,PP_nElems
            marked=.FALSE.
            DO iElem2=1,Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside
              IF(iElem.EQ.TempElemTotalInside(iElem2))THEN
                marked=.TRUE.
                EXIT
              END IF
            END DO ! iElem
#if USE_MPI
            IF(.NOT.Marked)THEN
              DO ilocSide=1,6
                SideID=ElemToSide(E2S_SIDE_ID,ilocSide,iElem)
                IF (SideID.GT.nInnerSides+nBCSides) marked=.TRUE.
              END DO ! ilocSide
            END IF
#endif /*USE_MPI*/
            IF(.NOT.marked) CYCLE
            DO ilocSide=1,6
              DO iShot = 1,200
                CALL RANDOM_NUMBER(RandVal)
                RandVal=2.0*RandVal-1.0
                SELECT CASE(ilocSide)
                CASE(XI_MINUS)
                  XI=(/-1.0,RandVal(1),RandVal(2)/)
                CASE(XI_PLUS)
                  XI=(/ 1.0,RandVal(1),RandVal(2)/)
                CASE(ETA_MINUS)
                  XI=(/RandVal(1),-1.0,RandVal(2)/)
                CASE(ETA_PLUS)
                  XI=(/RandVal(1), 1.0,RandVal(2)/)
                CASE(ZETA_MINUS)
                  XI=(/RandVal(1),RandVal(2),-1.0/)
                CASE(ZETA_PLUS)
                  XI=(/RandVal(1),RandVal(2), 1.0/)
                END SELECT
                CALL TensorProductInterpolation(Xi,3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,iElem),BN)
                ! check if point is in element
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
                  Species(iSpec)%Init(iInit)%ConstPress%ElemStat(iElem) = 2
                  Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside &
                       = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside + 1
                  TempElemPartlyInside(Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside) = iElem
                  EXIT
                END IF
              END DO ! iShot=1,200
            END DO ! ilocSide
          END DO ! iElem=1,Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside

          ! Shoot on Sides of Neighbour-Elements of partly inside Elements ==========================================
          nInterest = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside
          nInterOld = 1
          DO WHILE (nInterest .GE. nInterOld)
            DO iElem = nInterOld, nInterest
              Element = TempElemPartlyInside(iElem)
              DO ilocSide = 1, 6
                ExamElem =PartElemToElemAndSide(1,ilocSide,Element)
                locSideID=PartElemToElemAndSide(5,ilocSide,Element)
                IF(ExamElem.EQ.-1) CYCLE
                IF(ExamElem.GT.PP_nElems) CYCLE
                IF (Species(iSpec)%Init(iInit)%ConstPress%ElemStat(ExamElem) .EQ. 3) THEN
                  DO iShot = 1,200
                    CALL RANDOM_NUMBER(RandVal)
                    RandVal=2.0*RandVal-1.0
                    SELECT CASE(locSideID)
                    CASE(XI_MINUS)
                      XI=(/-1.0,RandVal(1),RandVal(2)/)
                    CASE(XI_PLUS)
                      XI=(/ 1.0,RandVal(1),RandVal(2)/)
                    CASE(ETA_MINUS)
                      XI=(/RandVal(1),-1.0,RandVal(2)/)
                    CASE(ETA_PLUS)
                      XI=(/RandVal(1), 1.0,RandVal(2)/)
                    CASE(ZETA_MINUS)
                      XI=(/RandVal(1),RandVal(2),-1.0/)
                    CASE(ZETA_PLUS)
                      XI=(/RandVal(1),RandVal(2), 1.0/)
                    END SELECT
                    CALL TensorProductInterpolation(Xi,3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,ExamElem),BN)
                    ! check if point is in element
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
          IF ((OV(1) .EQ. 0) .AND. (OV(2) .EQ. 0) .AND. (OV(3) .EQ. 0)) THEN
            CALL Abort(&
__STAMP__&
,'Error in InitializeVariables: Cannot calculate NormalVector(Cyl) of InitVolume in EmissionCase(3 or 4)!')
          END IF
          Species(iSpec)%Init(iInit)%ConstPress%OrthoVector &
                  = OV*(Species(iSpec)%Init(iInit)%CylinderHeightIC/SQRT(OV(1)**2 + OV(2)**2 + OV(3)**2))
          OV(1:3) = Species(iSpec)%Init(iInit)%ConstPress%OrthoVector(1:3)

          Species(iSpec)%Init(iInit)%ParticleEmission = &
               Species(iSpec)%Init(iInit)%ConstantPressure * Species(iSpec)%Init(iInit)%RadiusIC**2 * &
               Pi * Species(iSpec)%Init(iInit)%CylinderHeightIC / (BoltzmannConst * &
               Species(iSpec)%Init(iInit)%MWTemperatureIC * Species(iSpec)%MacroParticleFactor)
          Species(iSpec)%Init(iInit)%ConstPress%EkinInside = &
               1.5*Species(iSpec)%Init(iInit)%ConstantPressure*Species(iSpec)%Init(iInit)%RadiusIC**2 * &
               Pi * Species(iSpec)%Init(iInit)%CylinderHeightIC + &
               INT(Species(iSpec)%Init(iInit)%ParticleEmission)*0.5*Species(iSpec)%MassIC*Species(iSpec)%Init(iInit)%VeloIC**2 &
               * Species(iSpec)%MacroParticleFactor
          IPWRITE (UNIT_StdOut,'(I4,A49,I3.3,A2,g0)') 'Number of Particles inside ConstPressArea ',iInit, ': ', &
               Species(iSpec)%Init(iInit)%ParticleEmission

          DO iElem = 1,nElems
            nNodesInside = 0
            nBoundNodes = 0
            ! loop over all 8 nodes
            DO k=0,NGeo,NGeo
              DO j=0,NGeo,NGeo
                DO i=0,NGeo,NGeo
                  BN(1:3) = XCL_NGeo(1:3,i,j,k,iElem) - Species(iSpec)%Init(iInit)%BasePointIC
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
                 END DO ! i=0,NGeo
              END DO ! j=0,NGeo
            END DO ! k=0,NGeo

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
          END DO ! iElem

          ! If no Element has been found, search for Basepoint ========================================
          IF ((Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside .EQ. 0) .AND. &
               (Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside .EQ. 0)) THEN
            CALL PointToExactElement(Species(iSpec)%Init(iInit)%BasePointIC,Element,InElementCheck,doHalo=.FALSE.)
            IF (InElementCheck) THEN
              Species(iSpec)%Init(iInit)%ConstPress%ElemStat(Element) = 2
              Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside = 1
              TempElemPartlyInside(1) = Element
            END IF
          END IF

          !Schießen auf Nachbarzellen der totalen ===================================================
          ! Shoot on Sides of Neighbour-Elements of totally inside Elements ===========================
          DO iElem = 1,PP_nElems
            marked=.FALSE.
            DO iElem2=1,Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside
              IF(iElem.EQ.TempElemTotalInside(iElem2))THEN
                marked=.TRUE.
                EXIT
              END IF
            END DO ! iElem
#if USE_MPI
            IF(.NOT.Marked)THEN
              DO ilocSide=1,6
                SideID=ElemToSide(E2S_SIDE_ID,ilocSide,iElem)
                IF (SideID.GT.nInnerSides+nBCSides) marked=.TRUE.
              END DO ! ilocSide
            END IF
#endif /*USE_MPI*/
            IF(.NOT.marked) CYCLE
            DO ilocSide=1,6
              DO iShot = 1,200
                CALL RANDOM_NUMBER(RandVal)
                RandVal=2.0*RandVal-1.0
                SELECT CASE(ilocSide)
                CASE(XI_MINUS)
                  XI=(/-1.0,RandVal(1),RandVal(2)/)
                CASE(XI_PLUS)
                  XI=(/ 1.0,RandVal(1),RandVal(2)/)
                CASE(ETA_MINUS)
                  XI=(/RandVal(1),-1.0,RandVal(2)/)
                CASE(ETA_PLUS)
                  XI=(/RandVal(1), 1.0,RandVal(2)/)
                CASE(ZETA_MINUS)
                  XI=(/RandVal(1),RandVal(2),-1.0/)
                CASE(ZETA_PLUS)
                  XI=(/RandVal(1),RandVal(2), 1.0/)
                END SELECT
                CALL TensorProductInterpolation(Xi,3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,iElem),BN)
                ! check if point is in element
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
              END DO ! iShot=1,200
            END DO ! ilocSide
          END DO ! iElem=1,Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside


          !!!Schießen auf Nachbarelemente der teilweisen ===============================================
          nInterest = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside
          nInterOld = 1
          DO WHILE (nInterest .GE. nInterOld)
            DO iElem = nInterOld,nInterest
              Element = TempElemPartlyInside(iElem)
              DO ilocSide = 1, 6
                ExamElem =PartElemToElemAndSide(1    ,ilocSide,Element)
                locSideID=PartElemToElemAndSide(5,ilocSide,Element)
                IF(ExamElem.EQ.-1) CYCLE
                IF(ExamElem.GT.PP_nElems) CYCLE
                IF (Species(iSpec)%Init(iInit)%ConstPress%ElemStat(ExamElem) .EQ. 3) THEN
                  DO iShot = 1,200
                    CALL RANDOM_NUMBER(RandVal)
                    RandVal=2.0*RandVal-1.0
                    SELECT CASE(locSideID)
                    CASE(XI_MINUS)
                      XI=(/-1.0,RandVal(1),RandVal(2)/)
                    CASE(XI_PLUS)
                      XI=(/ 1.0,RandVal(1),RandVal(2)/)
                    CASE(ETA_MINUS)
                      XI=(/RandVal(1),-1.0,RandVal(2)/)
                    CASE(ETA_PLUS)
                      XI=(/RandVal(1), 1.0,RandVal(2)/)
                    CASE(ZETA_MINUS)
                      XI=(/RandVal(1),RandVal(2),-1.0/)
                    CASE(ZETA_PLUS)
                      XI=(/RandVal(1),RandVal(2), 1.0/)
                    END SELECT
                    CALL TensorProductInterpolation(Xi,3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,ExamElem),BN)
                    ! check if point is in element
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

        IPWRITE (UNIT_StdOut,'(I4,A49,I3.3,A2,I0)') 'Number of Elements inside ConstPressArea ',iInit, ': ', &
             Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside
        IPWRITE (UNIT_StdOut,'(I4,A49,I3.3,A2,I0)') 'Number of Elements partly inside ConstPressArea ',iInit, ': ', &
             Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside

      END IF
    END DO
  END DO

END SUBROUTINE ParticlePressureIni


SUBROUTINE ParticlePressureCellIni()
!===================================================================================================================================
! Initialization of constant pressure in a cell
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars
  USE MOD_Globals
  USE MOD_Preproc
  USE MOD_Globals_Vars
  USE MOD_Particle_Tracking_Vars,  ONLY:DoRefMapping
  USE MOD_Particle_Mesh_Vars,      ONLY:epsInCell
  USE MOD_Particle_Mesh,           ONLY:PointToExactElement
  USE MOD_Mesh_Vars,               ONLY:nElems,ElemToSide
  USE MOD_Mesh_Vars,               ONLY:NGeo,XCL_NGeo,XiCL_NGeo,wBaryCL_NGeo
  USE MOD_Eval_XYZ,                ONLY:TensorProductInterpolation
  USE MOD_Particle_Mesh_Vars,      ONLY:PartElemToElemAndSide
#if USE_MPI
  USE MOD_Mesh_Vars,               ONLY : nInnerSides, nBCSides
#endif
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                 :: BN(3), BV1(3), BV2(3), BV3(3), OV(3), epsi,XI(3)
  REAL                 :: det1, det2, det3, RandVal(2), dist1, dist2
  INTEGER, ALLOCATABLE :: TempElemTotalInside(:), TempElemPartlyInside(:)
  INTEGER              :: nNodesInside, nBoundNodes, nInterest, nInterOld, Element, ExamElem
  INTEGER              :: iShot, iElem, iSpec, iInit, iElem2
  INTEGER              :: iLocSide,SideID,locSideID
  INTEGER              :: i,j,k
  LOGICAL              :: InElementCheck,Marked
!===================================================================================================================================
  IF(.NOT.DoRefMapping) CALL abort(&
__STAMP__&
,' Particle pressure is only possible with tracking via reference element mapping!')


   CALL abort(&
__STAMP__&
,' GEO has to be exchanged by XCL_NGeo!' )
  !epsi = 100.*epsilon(epsi)
  ! accuracy of mapping into ref element
  epsi = epsInCell

  epsi = 100.*epsilon(epsi)
  DO iSpec = 1,nSpecies
    DO iInit = Species(iSpec)%StartnumberOfInits, Species(iSpec)%NumberOfInits
      Species(iSpec)%Init(iInit)%ConstPress%InitialTemp = Species(iSpec)%Init(iInit)%MWTemperatureIC
      IF ((Species(iSpec)%Init(iInit)%ParticleEmissionType .EQ. 4) &
        .OR.(Species(iSpec)%Init(iInit)%ParticleEmissionType .EQ. 6)) THEN
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
            CALL Abort(&
__STAMP__&
,'Error in InitializeVariables: Cannot calculate Normal Vector of InitVolume in EmissionCase(3 or 4)!')
          END IF

          Species(iSpec)%Init(iInit)%ConstPress%OrthoVector &
               = OV*(Species(iSpec)%Init(iInit)%CuboidHeightIC/SQRT(OV(1)**2 + OV(2)**2 + OV(3)**2))
          OV = Species(iSpec)%Init(iInit)%ConstPress%OrthoVector
          Species(iSpec)%Init(iInit)%ConstPress%Determinant = ABS(BV1(1)*BV2(2)*OV(3) + BV1(2)*BV2(3)*OV(1) + BV1(3)*BV2(1)*OV(2) &
               - BV1(3)*BV2(2)*OV(1) - BV1(1)*BV2(3)*OV(2) - BV1(2)*BV2(1)*OV(3))

          ! ParticleEmission multiplied by cell volume equals particles to be in cells
          Species(iSpec)%Init(iInit)%ParticleEmission = Species(iSpec)%Init(iInit)%ConstantPressure/ &
               (BoltzmannConst * Species(iSpec)%Init(iInit)%MWTemperatureIC * Species(iSpec)%MacroParticleFactor)

          DO iElem = 1,nElems
            nNodesInside = 0
            nBoundNodes = 0
            ! loop over all 8 nodes
            DO k=0,NGeo,NGeo
              DO j=0,NGeo,NGeo
                DO i=0,NGeo,NGeo
                  BN(1:3) = XCL_NGeo(1:3,i,j,k,iElem) - Species(iSpec)%Init(iInit)%BasePointIC
               !   BN(1:3) = GEO%NodeCoords(:,GEO%ElemToNodeID(iNode,iElem)) - Species(iSpec)%Init(iInit)%BasePointIC
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
                 END DO ! i=0,NGeo
              END DO ! j=0,NGeo
            END DO ! k=0,NGeo

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
            CALL PointToExactElement(Species(iSpec)%Init(iInit)%BasePointIC,Element,InElementCheck,doHalo=.FALSE.)
            IF (InElementCheck) THEN
              Species(iSpec)%Init(iInit)%ConstPress%ElemStat(Element) = 2
              Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside = 1
              TempElemPartlyInside(1) = Element
            END IF
          END IF

          ! Shoot on Sides of Neighbour-Elements of totally inside Elements ===========================
          DO iElem = 1,PP_nElems
            marked=.FALSE.
            DO iElem2=1,Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside
              IF(iElem.EQ.TempElemTotalInside(iElem2))THEN
                marked=.TRUE.
                EXIT
              END IF
            END DO ! iElem
#if USE_MPI
            IF(.NOT.Marked)THEN
              DO ilocSide=1,6
                SideID=ElemToSide(E2S_SIDE_ID,ilocSide,iElem)
                IF (SideID.GT.nInnerSides+nBCSides) marked=.TRUE.
              END DO ! ilocSide
            END IF
#endif /*USE_MPI*/
            IF(.NOT.marked) CYCLE
            DO ilocSide=1,6
              DO iShot = 1,200
                CALL RANDOM_NUMBER(RandVal)
                RandVal=2.0*RandVal-1.0
                SELECT CASE(ilocSide)
                CASE(XI_MINUS)
                  XI=(/-1.0,RandVal(1),RandVal(2)/)
                CASE(XI_PLUS)
                  XI=(/ 1.0,RandVal(1),RandVal(2)/)
                CASE(ETA_MINUS)
                  XI=(/RandVal(1),-1.0,RandVal(2)/)
                CASE(ETA_PLUS)
                  XI=(/RandVal(1), 1.0,RandVal(2)/)
                CASE(ZETA_MINUS)
                  XI=(/RandVal(1),RandVal(2),-1.0/)
                CASE(ZETA_PLUS)
                  XI=(/RandVal(1),RandVal(2), 1.0/)
                END SELECT
                CALL TensorProductInterpolation(Xi,3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,iElem),BN)
                ! check if point is in element
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
                  Species(iSpec)%Init(iInit)%ConstPress%ElemStat(iElem) = 2
                  Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside &
                       = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside + 1
                  TempElemPartlyInside(Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside) = iElem
                  EXIT
                END IF
              END DO ! iShot=1,200
            END DO ! ilocSide
          END DO ! iElem=1,Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside

          ! Shoot on Sides of Neighbour-Elements of partly inside Elements ==========================================
          nInterest = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside
          nInterOld = 1
          DO WHILE (nInterest .GE. nInterOld)
            DO iElem = nInterOld, nInterest
              Element = TempElemPartlyInside(iElem)
              DO ilocSide = 1, 6
                ExamElem =PartElemToElemAndSide(1,ilocSide,Element)
                locSideID=PartElemToElemAndSide(5,ilocSide,Element)
                IF(ExamElem.EQ.-1) CYCLE
                IF(ExamElem.GT.PP_nElems) CYCLE
                IF (Species(iSpec)%Init(iInit)%ConstPress%ElemStat(ExamElem) .EQ. 3) THEN
                  DO iShot = 1,200
                    CALL RANDOM_NUMBER(RandVal)
                    RandVal=2.0*RandVal-1.0
                    SELECT CASE(locSideID)
                    CASE(XI_MINUS)
                      XI=(/-1.0,RandVal(1),RandVal(2)/)
                    CASE(XI_PLUS)
                      XI=(/ 1.0,RandVal(1),RandVal(2)/)
                    CASE(ETA_MINUS)
                      XI=(/RandVal(1),-1.0,RandVal(2)/)
                    CASE(ETA_PLUS)
                      XI=(/RandVal(1), 1.0,RandVal(2)/)
                    CASE(ZETA_MINUS)
                      XI=(/RandVal(1),RandVal(2),-1.0/)
                    CASE(ZETA_PLUS)
                      XI=(/RandVal(1),RandVal(2), 1.0/)
                    END SELECT
                    CALL TensorProductInterpolation(Xi,3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,ExamElem),BN)
                    ! check if point is in element
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

          IF ((OV(1) .EQ. 0) .AND. (OV(2) .EQ. 0) .AND. (OV(3) .EQ. 0)) THEN
            CALL Abort(&
__STAMP__&
,'Error in InitializeVariables: Cannot calculate NormalVector(Cyl) of InitVolume in EmissionCase(3 or 4)!')
          END IF
          Species(iSpec)%Init(iInit)%ConstPress%OrthoVector = &
               OV*(Species(iSpec)%Init(iInit)%CylinderHeightIC/SQRT(OV(1)**2 + OV(2)**2 + OV(3)**2))
          OV(1:3) = Species(iSpec)%Init(iInit)%ConstPress%OrthoVector(1:3)

          Species(iSpec)%Init(iInit)%ParticleEmission = &
               Species(iSpec)%Init(iInit)%ConstantPressure * Species(iSpec)%Init(iInit)%RadiusIC**2 * &
               Pi * Species(iSpec)%Init(iInit)%CylinderHeightIC / (BoltzmannConst * &
               Species(iSpec)%Init(iInit)%MWTemperatureIC * Species(iSpec)%MacroParticleFactor)
          Species(iSpec)%Init(iInit)%ConstPress%EkinInside = &
               1.5*Species(iSpec)%Init(iInit)%ConstantPressure*Species(iSpec)%Init(iInit)%RadiusIC**2 * &
               Pi * Species(iSpec)%Init(iInit)%CylinderHeightIC + &
               INT(Species(iSpec)%Init(iInit)%ParticleEmission)*0.5*Species(iSpec)%MassIC*Species(iSpec)%Init(iInit)%VeloIC**2 &
               * Species(iSpec)%MacroParticleFactor

          DO iElem = 1,nElems
            nNodesInside = 0
            nBoundNodes = 0
            ! loop over all 8 nodes
            DO k=0,NGeo,NGeo
              DO j=0,NGeo,NGeo
                DO i=0,NGeo,NGeo
                  BN(1:3) = XCL_NGeo(1:3,i,j,k,iElem) - Species(iSpec)%Init(iInit)%BasePointIC
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
                 END DO ! i=0,NGeo
              END DO ! j=0,NGeo
            END DO ! k=0,NGeo

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
            CALL PointToExactElement(Species(iSpec)%Init(iInit)%BasePointIC,Element,InElementCheck,doHalo=.FALSE.)
            IF (InElementCheck) THEN
              Species(iSpec)%Init(iInit)%ConstPress%ElemStat(Element) = 2
              Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside = 1
              TempElemPartlyInside(1) = Element
            END IF
          END IF

          !Schießen auf Nachbarzellen der totalen ===================================================
          ! Shoot on Sides of Neighbour-Elements of totally inside Elements ===========================
          DO iElem = 1,PP_nElems
            marked=.FALSE.
            DO iElem2=1,Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside
              IF(iElem.EQ.TempElemTotalInside(iElem2))THEN
                marked=.TRUE.
                EXIT
              END IF
            END DO ! iElem
#if USE_MPI
            IF(.NOT.Marked)THEN
              DO ilocSide=1,6
                SideID=ElemToSide(E2S_SIDE_ID,ilocSide,iElem)
                IF (SideID.GT.nInnerSides+nBCSides) marked=.TRUE.
              END DO ! ilocSide
            END IF
#endif /*USE_MPI*/
            IF(.NOT.marked) CYCLE
            DO ilocSide=1,6
              DO iShot = 1,200
                CALL RANDOM_NUMBER(RandVal)
                RandVal=2.0*RandVal-1.0
                SELECT CASE(ilocSide)
                CASE(XI_MINUS)
                  XI=(/-1.0,RandVal(1),RandVal(2)/)
                CASE(XI_PLUS)
                  XI=(/ 1.0,RandVal(1),RandVal(2)/)
                CASE(ETA_MINUS)
                  XI=(/RandVal(1),-1.0,RandVal(2)/)
                CASE(ETA_PLUS)
                  XI=(/RandVal(1), 1.0,RandVal(2)/)
                CASE(ZETA_MINUS)
                  XI=(/RandVal(1),RandVal(2),-1.0/)
                CASE(ZETA_PLUS)
                  XI=(/RandVal(1),RandVal(2), 1.0/)
                END SELECT
                CALL TensorProductInterpolation(Xi,3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,iElem),BN)
                ! check if point is in element
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
              END DO ! iShot=1,200
            END DO ! ilocSide
          END DO ! iElem=1,Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside


          !!!Schießen auf Nachbarelemente der teilweisen ===============================================
          nInterest = Species(iSpec)%Init(iInit)%ConstPress%nElemPartlyInside
          nInterOld = 1
          DO WHILE (nInterest .GE. nInterOld)
            DO iElem = nInterOld,nInterest
              Element = TempElemPartlyInside(iElem)
              DO ilocSide = 1, 6
                ExamElem =PartElemToElemAndSide(1,ilocSide,Element)
                locSideID=PartElemToElemAndSide(5,ilocSide,Element)
                IF(ExamElem.EQ.-1) CYCLE
                IF(ExamElem.GT.PP_nElems) CYCLE
                IF (Species(iSpec)%Init(iInit)%ConstPress%ElemStat(ExamElem) .EQ. 3) THEN
                  DO iShot = 1,200
                    CALL RANDOM_NUMBER(RandVal)
                    RandVal=2.0*RandVal-1.0
                    SELECT CASE(locSideID)
                    CASE(XI_MINUS)
                      XI=(/-1.0,RandVal(1),RandVal(2)/)
                    CASE(XI_PLUS)
                      XI=(/ 1.0,RandVal(1),RandVal(2)/)
                    CASE(ETA_MINUS)
                      XI=(/RandVal(1),-1.0,RandVal(2)/)
                    CASE(ETA_PLUS)
                      XI=(/RandVal(1), 1.0,RandVal(2)/)
                    CASE(ZETA_MINUS)
                      XI=(/RandVal(1),RandVal(2),-1.0/)
                    CASE(ZETA_PLUS)
                      XI=(/RandVal(1),RandVal(2), 1.0/)
                    END SELECT
                    CALL TensorProductInterpolation(Xi,3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,ExamElem),BN)
                    ! check if point is in element
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
#if USE_MPI
        IF(Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside.NE.0)THEN
          IPWRITE (UNIT_StdOut,'(I4,A49,I3.3,A2,I0)') 'Number of Elements inside ConstPressArea ',iInit,': ', &
                              Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside
        END IF
#else
        WRITE (UNIT_StdOut,'(A49,I3.3,A2,I0)') 'Number of Elements inside ConstPressArea ',iInit,': ', &
               Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside
#endif
      END IF
    END DO
  END DO
END SUBROUTINE ParticlePressureCellIni


SUBROUTINE ParticlePressure (iSpec, iInit, NbrOfParticle)
!===================================================================================================================================
! Performs constant pressure calculations
!===================================================================================================================================
! MODULES
  USE MOD_Globals_Vars,       ONLY: BoltzmannConst
  USE MOD_Particle_Vars
#if USE_MPI
  USE MOD_Particle_MPI_Vars,  ONLY: PartMPI
  USE mpi
#endif
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)    :: iSpec, iInit
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  INTEGER, INTENT(OUT)   :: NbrOfParticle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER          :: nPartInside
  REAL             :: TempInside, EInside
#if USE_MPI
  REAL             :: TempComSend(2)
  REAL             :: TempComRec(2)
  INTEGER          :: IERROR
#endif
!===================================================================================================================================

    CALL ParticleInsideCheck(iSpec, iInit, nPartInside, TempInside, EInside)
#if USE_MPI
    TempComSend(1) = REAL(nPartInside)
    TempComSend(2) = EInside
    CALL MPI_ALLREDUCE(TempComSend, TempComRec, 2, MPI_DOUBLE_PRECISION, MPI_SUM, PartMPI%COMM, IERROR)
    nPartInside = INT(TempComRec(1))
    EInside = TempComRec(2)
#endif
    IF (nPartInside .LT. Species(iSpec)%Init(iInit)%ParticleEmission) THEN
      NbrOfParticle = INT(Species(iSpec)%Init(iInit)%ParticleEmission) - nPartInside
      IF ((Species(iSpec)%Init(iInit)%ConstPress%EkinInside - EInside).GT. 0.) THEN   !Both Factors should be bigger than 0
        Species(iSpec)%Init(iInit)%MWTemperatureIC=2./3.*(Species(iSpec)%Init(iInit)%ConstPress%EkinInside-EInside &
          - NbrOfParticle*0.5*Species(iSpec)%MassIC*Species(iSpec)%Init(iInit)%VeloIC**2*Species(iSpec)%MacroParticleFactor) &
          / (NbrOfParticle*Species(iSpec)%MacroParticleFactor*BoltzmannConst)
      Species(iSpec)%Init(iInit)%MWTemperatureIC=MAX(0.0, Species(iSpec)%Init(iInit)%MWTemperatureIC)
      ELSE
        Species(iSpec)%Init(iInit)%MWTemperatureIC = Species(iSpec)%Init(iInit)%ConstPress%InitialTemp
      END IF
    ELSE
      NbrOfParticle = 0
    END IF
END SUBROUTINE ParticlePressure


SUBROUTINE ParticlePressureRem (iSpec, iInit, NbrOfParticle)
!===================================================================================================================================
! Removes all particles in Pressure Area and sets NbrOfParticles to be inserted
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)    :: iSpec, iInit
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  INTEGER, INTENT(OUT)   :: NbrOfParticle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER          :: Particle
  REAL             :: BV1(3), BV2(3), BV3(3), OV(3), BN(3)
  REAL             :: det1, det2, det3, dist1, dist2
!===================================================================================================================================

  SELECT CASE (TRIM(Species(iSpec)%Init(iInit)%SpaceIC))
  CASE ('cuboid')
    BV1 = Species(iSpec)%Init(iInit)%BaseVector1IC
    BV2 = Species(iSpec)%Init(iInit)%BaseVector2IC
    OV  = Species(iSpec)%Init(iInit)%ConstPress%OrthoVector
    DO Particle = 1,PDM%ParticleVecLength
      IF ((PartSpecies(Particle) .EQ. iSpec) .AND. (PDM%ParticleInside(Particle))) THEN
        SELECT CASE (Species(iSpec)%Init(iInit)%ConstPress%ElemStat(PEM%Element(Particle)))
        CASE (1)
          PDM%ParticleInside(Particle)=.false.
        CASE (2)
          BN(1:3) = PartState(Particle,1:3) - Species(iSpec)%Init(iInit)%BasePointIC
          det1 = BN(1)*BV2(2)*OV(3) + BN(2)*BV2(3)*OV(1) + BN(3)*BV2(1)*OV(2) - &
                 BN(3)*BV2(2)*OV(1) - BN(1)*BV2(3)*OV(2) - BN(2)*BV2(1)*OV(3)
          det2 = BV1(1)*BN(2)*OV(3) + BV1(2)*BN(3)*OV(1) + BV1(3)*BN(1)*OV(2) - &
                 BV1(3)*BN(2)*OV(1) - BV1(1)*BN(3)*OV(2) - BV1(2)*BN(1)*OV(3)
          det3 = BV1(1)*BV2(2)*BN(3) + BV1(2)*BV2(3)*BN(1) + BV1(3)*BV2(1)*BN(2) - &
                 BV1(3)*BV2(2)*BN(1) - BV1(1)*BV2(3)*BN(2) - BV1(2)*BV2(1)*BN(3)

          det1 = det1/Species(iSpec)%Init(iInit)%ConstPress%Determinant
          det2 = det2/Species(iSpec)%Init(iInit)%ConstPress%Determinant
          det3 = det3/Species(iSpec)%Init(iInit)%ConstPress%Determinant

          IF (((det1-0.5)**2 .LE. 0.25).AND.((det2-0.5)**2 .LE. 0.25).AND.((det3-0.5)**2 .LE. 0.25)) THEN
            PDM%ParticleInside(Particle)=.false.
          END IF
        END SELECT
      END IF
    END DO
  CASE ('cylinder')
    OV  = Species(iSpec)%Init(iInit)%ConstPress%OrthoVector
    DO Particle = 1,PDM%ParticleVecLength
      IF ((PartSpecies(Particle) .EQ. iSpec) .AND. (PDM%ParticleInside(Particle))) THEN
        SELECT CASE (Species(iSpec)%Init(iInit)%ConstPress%ElemStat(PEM%Element(Particle)))
        CASE (1)
          PDM%ParticleInside(Particle)=.false.
        CASE (2)
          BN(1:3) = PartState(Particle,1:3) - Species(iSpec)%Init(iInit)%BasePointIC
          BV2(1) = BN(2) * OV(3) - BN(3) * OV(2)                   !Vector orthogonal on BN and OrthoVector
          BV2(2) = BN(3) * OV(1) - BN(1) * OV(3)
          BV2(3) = BN(1) * OV(2) - BN(2) * OV(1)
          dist1  = SQRT((BV2(1)**2 + BV2(2)**2 + BV2(3)**2)/(OV(1)**2 + OV(2)**2 + OV(3)**2))

          IF (dist1 .LE. Species(iSpec)%Init(iInit)%RadiusIC) THEN
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
  NbrOfParticle = INT(Species(iSpec)%Init(iInit)%ParticleEmission)
END SUBROUTINE ParticlePressureRem


SUBROUTINE ParticleInsideCheck(iSpec, iInit, nPartInside, TempInside, EkinInside)
!===================================================================================================================================
! Checks how many particles are inside including their temperature/energy
!===================================================================================================================================
! MODULES
  USE MOD_Globals_Vars,         ONLY: BoltzmannConst
  USE MOD_Particle_Vars
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER          :: nPartInside, Particle, iSpec, iInit
  REAL             :: BV1(3), BV2(3), BV3(3), OV(3), BN(3), vau(3), vauquad(3), TempVec(3)
  REAL             :: det1, det2, det3, dist1, dist2, TempInside, EkinInside
!===================================================================================================================================

  nPartInside = 0
  vau(:) = 0
  vauquad(:) = 0
  SELECT CASE (TRIM(Species(iSpec)%Init(iInit)%SpaceIC))
  CASE ('cuboid')
    BV1 = Species(iSpec)%Init(iInit)%BaseVector1IC
    BV2 = Species(iSpec)%Init(iInit)%BaseVector2IC
    OV  = Species(iSpec)%Init(iInit)%ConstPress%OrthoVector
    DO Particle = 1,PDM%ParticleVecLength
      IF ((PartSpecies(Particle) .EQ. iSpec) .AND. (PDM%ParticleInside(Particle))) THEN
        SELECT CASE (Species(iSpec)%Init(iInit)%ConstPress%ElemStat(PEM%Element(Particle)))
        CASE (1)
          nPartInside = nPartInside + 1
          vau(1) = vau(1) + PartState(Particle,4)
          vau(2) = vau(2) + PartState(Particle,5)
          vau(3) = vau(3) + PartState(Particle,6)
          vauquad(1) = vauquad(1) + PartState(Particle,4)**2
          vauquad(2) = vauquad(2) + PartState(Particle,5)**2
          vauquad(3) = vauquad(3) + PartState(Particle,6)**2
        CASE (2)
          BN(1:3) = PartState(Particle,1:3) - Species(iSpec)%Init(iInit)%BasePointIC
          det1 = BN(1)*BV2(2)*OV(3) + BN(2)*BV2(3)*OV(1) + BN(3)*BV2(1)*OV(2) - &
                 BN(3)*BV2(2)*OV(1) - BN(1)*BV2(3)*OV(2) - BN(2)*BV2(1)*OV(3)
          det2 = BV1(1)*BN(2)*OV(3) + BV1(2)*BN(3)*OV(1) + BV1(3)*BN(1)*OV(2) - &
                 BV1(3)*BN(2)*OV(1) - BV1(1)*BN(3)*OV(2) - BV1(2)*BN(1)*OV(3)
          det3 = BV1(1)*BV2(2)*BN(3) + BV1(2)*BV2(3)*BN(1) + BV1(3)*BV2(1)*BN(2) - &
                 BV1(3)*BV2(2)*BN(1) - BV1(1)*BV2(3)*BN(2) - BV1(2)*BV2(1)*BN(3)

          det1 = det1/Species(iSpec)%Init(iInit)%ConstPress%Determinant
          det2 = det2/Species(iSpec)%Init(iInit)%ConstPress%Determinant
          det3 = det3/Species(iSpec)%Init(iInit)%ConstPress%Determinant

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
    OV  = Species(iSpec)%Init(iInit)%ConstPress%OrthoVector
    DO Particle = 1,PDM%ParticleVecLength
      IF ((PartSpecies(Particle) .EQ. iSpec) .AND. (PDM%ParticleInside(Particle))) THEN
        SELECT CASE (Species(iSpec)%Init(iInit)%ConstPress%ElemStat(PEM%Element(Particle)))
        CASE (1)
          nPartInside = nPartInside + 1
          vau = vau + PartState(Particle,4:6)
          vauquad = vauquad + PartState(Particle,4:6)**2
        CASE (2)
          BN(1:3) = PartState(Particle,1:3) - Species(iSpec)%Init(iInit)%BasePointIC
          BV2(1) = BN(2) * OV(3) - BN(3) * OV(2)                   !Vector orthogonal on BN and OrthoVector
          BV2(2) = BN(3) * OV(1) - BN(1) * OV(3)
          BV2(3) = BN(1) * OV(2) - BN(2) * OV(1)
          dist1  = SQRT((BV2(1)**2 + BV2(2)**2 + BV2(3)**2)/(OV(1)**2 + OV(2)**2 + OV(3)**2))

          IF (dist1 .LE. Species(iSpec)%Init(iInit)%RadiusIC) THEN
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
  IF (nPartInside.eq.0) THEN
    EkinInside = 0.
    TempInside = 0.
  ELSE
    vau = vau/nPartInside
    EkinInside = 0.5*Species(iSpec)%MassIC*(vauquad(1) + vauquad(2) +vauquad(3))*Species(iSpec)%MacroParticleFactor
    vauquad = vauquad/nPartInside
    TempVec(1) = Species(iSpec)%MassIC/BoltzmannConst*(vauquad(1) - vau(1)**2)
    TempVec(2) = Species(iSpec)%MassIC/BoltzmannConst*(vauquad(2) - vau(2)**2)
    TempVec(3) = Species(iSpec)%MassIC/BoltzmannConst*(vauquad(3) - vau(3)**2)
    TempInside = (TempVec(1) + TempVec(2) + TempVec(3))/3
  END IF

END SUBROUTINE ParticleInsideCheck

END MODULE MOD_part_pressure
