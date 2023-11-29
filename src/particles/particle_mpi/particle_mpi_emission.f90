!==================================================================================================================================
! Copyright (c) 2010 - 2019 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (piclas.boltzplatz.eu/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"

!===================================================================================================================================
! module for MPI communication during particle emission
!===================================================================================================================================
MODULE MOD_Particle_MPI_Emission
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

#if USE_MPI
INTERFACE InitEmissionComm
  MODULE PROCEDURE InitEmissionComm
END INTERFACE

INTERFACE SendEmissionParticlesToProcs
  MODULE PROCEDURE SendEmissionParticlesToProcs
END INTERFACE

PUBLIC :: InitEmissionComm
PUBLIC :: SendEmissionParticlesToProcs
!===================================================================================================================================
CONTAINS


SUBROUTINE InitEmissionComm()
!===================================================================================================================================
! build emission communicators for particle emission regions
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_MPI_Vars  ,ONLY: PartMPIInitGroup,MPI_halo_eps
USE MOD_Particle_Vars      ,ONLY: Species,nSpecies
USE MOD_Particle_Mesh_Vars ,ONLY: GEO,SideInfo_Shared
USE MOD_Mesh_Vars          ,ONLY: nElems,BoundaryName
USE MOD_Particle_Vars      ,ONLY: NeutralizationSource,nNeutralizationElems,isNeutralizationElem,NeutralizationBalanceElem
#if ! (USE_HDG)
USE MOD_CalcTimeStep       ,ONLY: CalcTimeStep
#endif /*USE_HDG*/
USE MOD_Mesh_Vars          ,ONLY: offsetElem
USE MOD_Particle_Mesh_Vars ,ONLY: ElemInfo_Shared
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: GlobalElemID,iElem,BCID,iSide,iSpec,iInit,iNode,iRank
INTEGER                         :: nInitRegions
LOGICAL                         :: RegionOnProc,RegionExists
REAL                            :: xCoords(3,8),lineVector(3),radius,height
REAL                            :: xlen,ylen,zlen
REAL                            :: dt
INTEGER                         :: color
CHARACTER(50)                   :: hilf
!===================================================================================================================================

! get number of total init regions
nInitRegions=0
DO iSpec=1,nSpecies
  nInitRegions=nInitRegions+Species(iSpec)%NumberOfInits
END DO ! iSpec
IF(nInitRegions.EQ.0) RETURN

! allocate communicators
ALLOCATE( PartMPIInitGroup(1:nInitRegions))

! Default value for neutralization regions (Landmark and Liu2010)
nNeutralizationElems = -1

nInitRegions=0
DO iSpec=1,nSpecies
  RegionOnProc=.FALSE.
  DO iInit=1, Species(iSpec)%NumberOfInits
    nInitRegions=nInitRegions+1
    SELECT CASE(TRIM(Species(iSpec)%Init(iInit)%SpaceIC))
    CASE ('point')
       xCoords(1:3,1)=Species(iSpec)%Init(iInit)%BasePointIC
       RegionOnProc=PointInProc(xCoords(1:3,1))
    CASE ('line_with_equidistant_distribution')
      xCoords(1:3,1)=Species(iSpec)%Init(iInit)%BasePointIC
      xCoords(1:3,2)=Species(iSpec)%Init(iInit)%BasePointIC+Species(iSpec)%Init(iInit)%BaseVector1IC
      RegionOnProc=BoxInProc(xCoords(1:3,1:2),2)
    CASE ('line')
      xCoords(1:3,1)=Species(iSpec)%Init(iInit)%BasePointIC
      xCoords(1:3,2)=Species(iSpec)%Init(iInit)%BasePointIC+Species(iSpec)%Init(iInit)%BaseVector1IC
      RegionOnProc=BoxInProc(xCoords(1:3,1:2),2)
    CASE('disc')
      xlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(1)*Species(iSpec)%Init(iInit)%NormalIC(1))
      ylen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(2)*Species(iSpec)%Init(iInit)%NormalIC(2))
      zlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(3)*Species(iSpec)%Init(iInit)%NormalIC(3))
      ! all 8 edges
      xCoords(1:3,1) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,-zlen/)
      xCoords(1:3,2) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,-zlen/)
      xCoords(1:3,3) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,-zlen/)
      xCoords(1:3,4) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,-zlen/)
      xCoords(1:3,5) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,+zlen/)
      xCoords(1:3,6) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,+zlen/)
      xCoords(1:3,7) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,+zlen/)
      xCoords(1:3,8) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,+zlen/)
      RegionOnProc=BoxInProc(xCoords(1:3,1:8),8)
    CASE('photon_SEE_disc','photon_SEE_honeycomb')
      ASSOCIATE( O  => Species(iSpec)%Init(iInit)%BasePointIC ,&
                 v1 => Species(iSpec)%Init(iInit)%NormalIC     )
        ! 1. Check if inside outer radius
      radius=Species(iSpec)%Init(iInit)%RadiusIC
        xlen=radius * SQRT(1.0 - v1(1)*v1(1))
        ylen=radius * SQRT(1.0 - v1(2)*v1(2))
        zlen=radius * SQRT(1.0 - v1(3)*v1(3)) + 0.1*radius ! 10 percent of radius as height
      ! all 8 edges
        xCoords(1:3,1) = O+(/-xlen,-ylen,-zlen/)
        xCoords(1:3,2) = O+(/+xlen,-ylen,-zlen/)
        xCoords(1:3,3) = O+(/-xlen,+ylen,-zlen/)
        xCoords(1:3,4) = O+(/+xlen,+ylen,-zlen/)
        xCoords(1:3,5) = O+(/-xlen,-ylen,+zlen/)
        xCoords(1:3,6) = O+(/+xlen,-ylen,+zlen/)
        xCoords(1:3,7) = O+(/-xlen,+ylen,+zlen/)
        xCoords(1:3,8) = O+(/+xlen,+ylen,+zlen/)
      ! Check if inside box
      RegionOnProc=BoxInProc(xCoords(1:3,1:8),8)
      END ASSOCIATE
    CASE('2D_landmark','2D_landmark_copy')
       ! Ionization profile from T. Charoy, 2D axial-azimuthal particle-in-cell benchmark
       ! for low-temperature partially magnetized plasmas (2019)
       ASSOCIATE( x2 => 1.0e-2       ,& ! m
                  x1 => 0.25e-2      ,& ! m
                  y2 => GEO%ymaxglob ,& ! m
                  y1 => GEO%yminglob ,& ! m
                  z2 => GEO%zmaxglob ,& ! m
                  z1 => GEO%zminglob )
        ! Check all 8 edges
        xCoords(1:3,1) = (/x1,y1,z1/)
        xCoords(1:3,2) = (/x2,y1,z1/)
        xCoords(1:3,3) = (/x1,y2,z1/)
        xCoords(1:3,4) = (/x2,y2,z1/)
        xCoords(1:3,5) = (/x1,y1,z2/)
        xCoords(1:3,6) = (/x2,y1,z2/)
        xCoords(1:3,7) = (/x1,y2,z2/)
        xCoords(1:3,8) = (/x2,y2,z2/)
        RegionOnProc=BoxInProc(xCoords(1:3,1:8),8)
      END ASSOCIATE
    CASE('2D_landmark_neutralization')
      ! Neutralization at const. x-position from T. Charoy, 2D axial-azimuthal particle-in-cell benchmark
      ! for low-temperature partially magnetized plasmas (2019)
      ! Check 1st region (emission at fixed x-position x=2.4cm)
      ASSOCIATE( &
                 x2 => 2.4001e-2    ,& ! m
                 x1 => 2.3999e-2-MPI_halo_eps ,& ! m
                 y2 => GEO%ymaxglob ,& ! m
                 y1 => GEO%yminglob ,& ! m
                 z2 => GEO%zmaxglob ,& ! m
                 z1 => GEO%zminglob )
       ! Check all 8 edges
       xCoords(1:3,1) = (/x1,y1,z1/)
       xCoords(1:3,2) = (/x2,y1,z1/)
       xCoords(1:3,3) = (/x1,y2,z1/)
       xCoords(1:3,4) = (/x2,y2,z1/)
       xCoords(1:3,5) = (/x1,y1,z2/)
       xCoords(1:3,6) = (/x2,y1,z2/)
       xCoords(1:3,7) = (/x1,y2,z2/)
       xCoords(1:3,8) = (/x2,y2,z2/)
       RegionOnProc=BoxInProc(xCoords(1:3,1:8),8)
      END ASSOCIATE

      ! Check 2nd region (left boundary where the exiting particles are counted)
      IF(.NOT.RegionOnProc)THEN
        ASSOCIATE(&
                   x2 => 0.0001e-2    ,& ! m
                   x1 => -0.001e-2    ,& ! m
                   y2 => GEO%ymaxglob ,& ! m
                   y1 => GEO%yminglob ,& ! m
                   z2 => GEO%zmaxglob ,& ! m
                   z1 => GEO%zminglob )
         ! Check all 8 edges
         xCoords(1:3,1) = (/x1,y1,z1/)
         xCoords(1:3,2) = (/x2,y1,z1/)
         xCoords(1:3,3) = (/x1,y2,z1/)
         xCoords(1:3,4) = (/x2,y2,z1/)
         xCoords(1:3,5) = (/x1,y1,z2/)
         xCoords(1:3,6) = (/x2,y1,z2/)
         xCoords(1:3,7) = (/x1,y2,z2/)
         xCoords(1:3,8) = (/x2,y2,z2/)
         RegionOnProc=BoxInProc(xCoords(1:3,1:8),8)
      END ASSOCIATE
      END IF ! .NOT.RegionOnProc
    CASE('2D_Liu2010_neutralization')
      ! Neutralization at right BC (max. x-position) H. Liu "Particle-in-cell simulation of a Hall thruster" (2010)
      ! Check one region (emission at fixed x-position x=30 mm)
      ASSOCIATE( &
                 x2 => 30.01e-3    ,& ! m
                 x1 => 29.99e-3-MPI_halo_eps ,& ! m
                 y2 => GEO%ymaxglob ,& ! m
                 y1 => GEO%yminglob ,& ! m
                 z2 => GEO%zmaxglob ,& ! m
                 z1 => GEO%zminglob )
       ! Check all 8 edges
       xCoords(1:3,1) = (/x1,y1,z1/)
       xCoords(1:3,2) = (/x2,y1,z1/)
       xCoords(1:3,3) = (/x1,y2,z1/)
       xCoords(1:3,4) = (/x2,y2,z1/)
       xCoords(1:3,5) = (/x1,y1,z2/)
       xCoords(1:3,6) = (/x2,y1,z2/)
       xCoords(1:3,7) = (/x1,y2,z2/)
       xCoords(1:3,8) = (/x2,y2,z2/)
       RegionOnProc=BoxInProc(xCoords(1:3,1:8),8)
      END ASSOCIATE
    CASE('2D_Liu2010_neutralization_Szabo','3D_Liu2010_neutralization_Szabo')

      ! Find all elements that have a neutralization BC and add up the number
      nNeutralizationElems = 0
      ALLOCATE(isNeutralizationElem(1:nElems))
      isNeutralizationElem = .FALSE.
      ELEMLOOP: DO iElem=1,nElems ! loop over all local elems
        GlobalElemID = iElem + offsetElem
        ! Check 6 local sides
        DO iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID)+1,ElemInfo_Shared(ELEM_LASTSIDEIND,GlobalElemID)
          ! Get BC index of the global side index
          BCID   = SideInfo_Shared(SIDE_BCID,iSide)
          ! Only check BC sides with BC index > 0
          IF(BCID.GT.0)THEN
            ! Check if neutralization BC is found
            IF(TRIM(BoundaryName(BCID)).EQ.TRIM(NeutralizationSource))THEN
              ! Add up the number of neutralization elems
              nNeutralizationElems = nNeutralizationElems + 1
              ! Flag element
              isNeutralizationElem(iElem) = .TRUE.
              ! Go to the next element
              CYCLE ELEMLOOP
            END IF
          END IF ! BCID.GT.0
        END DO
      END DO ELEMLOOP

      ! Only processors with neutralization elements are part of the communicator (+root which outputs global information to .csv)
      RegionOnProc = nNeutralizationElems.GT.0
      ! If no neutralization elements are present, deallocate the logical array
      IF(nNeutralizationElems.EQ.0)THEN
        DEALLOCATE(isNeutralizationElem)
      ELSE
        ALLOCATE(NeutralizationBalanceElem(1:nElems))
      END IF

    CASE('3D_Liu2010_neutralization')
      ! Neutralization at right BC (max. z-position) H. Liu "Particle-in-cell simulation of a Hall thruster" (2010)
      ! Check one region (emission at fixed z-position x=30 mm)
      ASSOCIATE( &
                 x2 => GEO%xmaxglob  ,& ! m
                 x1 => GEO%xminglob  ,& ! m
                 y2 => GEO%ymaxglob ,& ! m
                 y1 => GEO%yminglob ,& ! m
                 z2 => 30.01e-3 ,& ! m
                 z1 => 29.99e-3-MPI_halo_eps)
       ! Check all 8 edges
       xCoords(1:3,1) = (/x1,y1,z1/)
       xCoords(1:3,2) = (/x2,y1,z1/)
       xCoords(1:3,3) = (/x1,y2,z1/)
       xCoords(1:3,4) = (/x2,y2,z1/)
       xCoords(1:3,5) = (/x1,y1,z2/)
       xCoords(1:3,6) = (/x2,y1,z2/)
       xCoords(1:3,7) = (/x1,y2,z2/)
       xCoords(1:3,8) = (/x2,y2,z2/)
       RegionOnProc=BoxInProc(xCoords(1:3,1:8),8)
      END ASSOCIATE
    CASE('circle')
      xlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(1)*Species(iSpec)%Init(iInit)%NormalIC(1))
      ylen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(2)*Species(iSpec)%Init(iInit)%NormalIC(2))
      zlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(3)*Species(iSpec)%Init(iInit)%NormalIC(3))
      ! all 8 edges
      xCoords(1:3,1) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,-zlen/)
      xCoords(1:3,2) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,-zlen/)
      xCoords(1:3,3) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,-zlen/)
      xCoords(1:3,4) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,-zlen/)
      xCoords(1:3,5) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,+zlen/)
      xCoords(1:3,6) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,+zlen/)
      xCoords(1:3,7) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,+zlen/)
      xCoords(1:3,8) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,+zlen/)
      RegionOnProc=BoxInProc(xCoords(1:3,1:8),8)
    CASE('gyrotron_circle')
      Radius=Species(iSpec)%Init(iInit)%RadiusIC+Species(iSpec)%Init(iInit)%RadiusICGyro
      !xlen=Species(iSpec)%Init(iInit)%RadiusIC * &
      !     SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(1)*Species(iSpec)%Init(iInit)%NormalIC(1))
      !ylen=Species(iSpec)%Init(iInit)%RadiusIC * &
      !     SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(2)*Species(iSpec)%Init(iInit)%NormalIC(2))
      xlen=Radius
      ylen=Radius
      zlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(3)*Species(iSpec)%Init(iInit)%NormalIC(3))
      IF(Species(iSpec)%Init(iInit)%ParticleNumber.NE.0)THEN
        lineVector(1:3)=(/0.,0.,Species(iSpec)%Init(iInit)%CylinderHeightIC/)
      ELSE
#if !(USE_HDG)
        dt = CALCTIMESTEP()
#endif /*USE_HDG*/
        lineVector(1:3)= dt* Species(iSpec)%Init(iInit)%VeloIC/Species(iSpec)%Init(iInit)%alpha
        zlen=0.
      END IF
      xCoords(1:3,1) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,-zlen/)
      xCoords(1:3,2) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,-zlen/)
      xCoords(1:3,3) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,-zlen/)
      xCoords(1:3,4) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,-zlen/)
      xCoords(1:3,5) = Species(iSpec)%Init(iInit)%BasePointIC+lineVector+(/-xlen,-ylen,+zlen/)
      xCoords(1:3,6) = Species(iSpec)%Init(iInit)%BasePointIC+lineVector+(/+xlen,-ylen,+zlen/)
      xCoords(1:3,7) = Species(iSpec)%Init(iInit)%BasePointIC+lineVector+(/-xlen,+ylen,+zlen/)
      xCoords(1:3,8) = Species(iSpec)%Init(iInit)%BasePointIC+lineVector+(/+xlen,+ylen,+zlen/)
      RegionOnProc=BoxInProc(xCoords(1:3,1:8),8)
    CASE('circle_equidistant')
      xlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(1)*Species(iSpec)%Init(iInit)%NormalIC(1))
      ylen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(2)*Species(iSpec)%Init(iInit)%NormalIC(2))
      zlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(3)*Species(iSpec)%Init(iInit)%NormalIC(3))
      ! all 8 edges
      xCoords(1:3,1) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,-zlen/)
      xCoords(1:3,2) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,-zlen/)
      xCoords(1:3,3) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,-zlen/)
      xCoords(1:3,4) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,-zlen/)
      xCoords(1:3,5) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,+zlen/)
      xCoords(1:3,6) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,+zlen/)
      xCoords(1:3,7) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,+zlen/)
      xCoords(1:3,8) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,+zlen/)
      RegionOnProc=BoxInProc(xCoords(1:3,1:8),8)
    CASE('cuboid','photon_rectangle','photon_SEE_rectangle')
      ASSOCIATE( O => Species(iSpec)%Init(iInit)%BasePointIC   ,&
                v2 => Species(iSpec)%Init(iInit)%BaseVector1IC ,&
                v3 => Species(iSpec)%Init(iInit)%BaseVector2IC ,&
                normal => Species(iSpec)%Init(iInit)%NormalIC)

        ! Use NormalIC if it is non-zero. Else, use the cross-product calculated 3rd vector to span the coordinate space
        ! The first choice results in a non-rectangular 3rd coordinate vector
        IF(VECNORM(normal).LE.0.)THEN
          lineVector(1) = v2(2) * v3(3) - v2(3) * v3(2)
          lineVector(2) = v2(3) * v3(1) - v2(1) * v3(3)
          lineVector(3) = v2(1) * v3(2) - v2(2) * v3(1)
          lineVector = UNITVECTOR(lineVector)
          IF(VECNORM(lineVector).LE.0.) CALL ABORT(__STAMP__,'BaseVectors are parallel!')
        ELSE
          ! Normalize the vector even though it is probably already normalized for safety reasons
          lineVector = UNITVECTOR(normal)
        END IF ! VECNORM(lineVector).LE.0.

        xCoords(1:3,1)=O(1:3)
        xCoords(1:3,2)=O(1:3)+v2(1:3)
        xCoords(1:3,3)=O(1:3)+v3(1:3)
        xCoords(1:3,4)=O(1:3)+v2(1:3)+v3(1:3)

        height= Species(iSpec)%Init(iInit)%CuboidHeightIC
        DO iNode=1,4
          xCoords(1:3,iNode+4)=xCoords(1:3,iNode)+lineVector*height
        END DO ! iNode
        RegionOnProc=BoxInProc(xCoords,8)
      END ASSOCIATE
    CASE('sphere')
      ASSOCIATE ( radius => Species(iSpec)%Init(iInit)%RadiusIC        ,&
                  origin => Species(iSpec)%Init(iInit)%BasePointIC(1:3) )
        ! Set the 8 bounding box coordinates depending on the origin and radius
        xCoords(1:3,1)=origin + (/ radius  , -radius , -radius/)
        xCoords(1:3,2)=origin + (/ radius  , radius  , -radius/)
        xCoords(1:3,3)=origin + (/ -radius , radius  , -radius/)
        xCoords(1:3,4)=origin + (/ -radius , -radius , -radius/)
        xCoords(1:3,5)=origin + (/ radius  , -radius , radius /)
        xCoords(1:3,6)=origin + (/ radius  , radius  , radius /)
        xCoords(1:3,7)=origin + (/ -radius , radius  , radius /)
        xCoords(1:3,8)=origin + (/ -radius , -radius , radius /)
      END ASSOCIATE
      RegionOnProc=BoxInProc(xCoords,8)
    CASE('cylinder','photon_cylinder','photon_honeycomb')
      ASSOCIATE( v1 => Species(iSpec)%Init(iInit)%BaseVector1IC,&
                 v2 => Species(iSpec)%Init(iInit)%BaseVector2IC)
      ! Cross-product of 1IC and 2IC
        lineVector(1) = v1(2) * v2(3) - v1(3) * v2(2)
        lineVector(2) = v1(3) * v2(1) - v1(1) * v2(3)
        lineVector(3) = v1(1) * v2(2) - v1(2) * v2(1)
      ! Sanity check
      IF ((lineVector(1).eq.0).AND.(lineVector(2).eq.0).AND.(lineVector(3).eq.0)) THEN
         CALL ABORT(__STAMP__,'BaseVectors are parallel!')
      ELSE
        lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * lineVector(2) + &
          lineVector(3) * lineVector(3))
      END IF
      ! 1. Check if inside outer radius
      radius = Species(iSpec)%Init(iInit)%RadiusIC
      ! here no radius, already included
      xCoords(1:3,1)=Species(iSpec)%Init(iInit)%BasePointIC -v1(1:3) -v2(1:3)

      xCoords(1:3,2)=xCoords(1:3,1)+2.0*v1(1:3)
      xCoords(1:3,3)=xCoords(1:3,1)+2.0*v2(1:3)
      xCoords(1:3,4)=xCoords(1:3,1)+2.0*v1(1:3)+2.0*v2(1:3)

      height= Species(iSpec)%Init(iInit)%CylinderHeightIC
      DO iNode=1,4
        xCoords(1:3,iNode+4)=xCoords(1:3,iNode)+lineVector*height
      END DO ! iNode
      ! Check if inside box
      RegionOnProc=BoxInProc(xCoords,8)
      END ASSOCIATE
    CASE('cell_local')
      RegionOnProc=.TRUE.
    CASE('cuboid_equal')
       xlen = SQRT(Species(iSpec)%Init(iInit)%BaseVector1IC(1)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector1IC(2)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector1IC(3)**2 )
       ylen = SQRT(Species(iSpec)%Init(iInit)%BaseVector2IC(1)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector2IC(2)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector2IC(3)**2 )
       zlen = ABS(Species(iSpec)%Init(iInit)%CuboidHeightIC)

       ! make sure the vectors correspond to x,y,z-dir
       IF ((xlen.NE.Species(iSpec)%Init(iInit)%BaseVector1IC(1)).OR. &
           (ylen.NE.Species(iSpec)%Init(iInit)%BaseVector2IC(2)).OR. &
           (zlen.NE.Species(iSpec)%Init(iInit)%CuboidHeightIC)) THEN
          CALL ABORT(__STAMP__&
          ,'Basevectors1IC,-2IC and CuboidHeightIC have to be in x,y,z-direction, respectively for emission condition')
       END IF
       DO iNode=1,8
        xCoords(1:3,iNode) = Species(iSpec)%Init(iInit)%BasePointIC(1:3)
       END DO
       xCoords(1:3,2) = xCoords(1:3,1) + (/xlen,0.,0./)
       xCoords(1:3,3) = xCoords(1:3,1) + (/0.,ylen,0./)
       xCoords(1:3,4) = xCoords(1:3,1) + (/xlen,ylen,0./)
       xCoords(1:3,5) = xCoords(1:3,1) + (/0.,0.,zlen/)
       xCoords(1:3,6) = xCoords(1:3,5) + (/xlen,0.,0./)
       xCoords(1:3,7) = xCoords(1:3,5) + (/0.,ylen,0./)
       xCoords(1:3,8) = xCoords(1:3,5) + (/xlen,ylen,0./)
       RegionOnProc=BoxInProc(xCoords,8)

     !~j CALL ABORT(&
     !~j __STAMP__&
     !~j ,'ERROR in ParticleEmission_parallel: cannot deallocate particle_positions!')
    CASE ('cuboid_with_equidistant_distribution')
       xlen = SQRT(Species(iSpec)%Init(iInit)%BaseVector1IC(1)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector1IC(2)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector1IC(3)**2 )
       ylen = SQRT(Species(iSpec)%Init(iInit)%BaseVector2IC(1)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector2IC(2)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector2IC(3)**2 )
       zlen = ABS(Species(iSpec)%Init(iInit)%CuboidHeightIC)

       ! make sure the vectors correspond to x,y,z-dir
       IF ((xlen.NE.Species(iSpec)%Init(iInit)%BaseVector1IC(1)).OR. &
           (ylen.NE.Species(iSpec)%Init(iInit)%BaseVector2IC(2)).OR. &
           (zlen.NE.Species(iSpec)%Init(iInit)%CuboidHeightIC)) THEN
          CALL ABORT(__STAMP__&
          ,'Basevectors1IC,-2IC and CuboidHeightIC have to be in x,y,z-direction, respectively for emission condition')
       END IF
       DO iNode=1,8
        xCoords(1:3,iNode) = Species(iSpec)%Init(iInit)%BasePointIC(1:3)
       END DO
       xCoords(1:3,2) = xCoords(1:3,1) + (/xlen,0.,0./)
       xCoords(1:3,3) = xCoords(1:3,1) + (/0.,ylen,0./)
       xCoords(1:3,4) = xCoords(1:3,1) + (/xlen,ylen,0./)
       xCoords(1:3,5) = xCoords(1:3,1) + (/0.,0.,zlen/)
       xCoords(1:3,6) = xCoords(1:3,5) + (/xlen,0.,0./)
       xCoords(1:3,7) = xCoords(1:3,5) + (/0.,ylen,0./)
       xCoords(1:3,8) = xCoords(1:3,5) + (/xlen,ylen,0./)
       RegionOnProc=BoxInProc(xCoords,8)
    CASE('sin_deviation')
       IF(Species(iSpec)%Init(iInit)%ParticleNumber.NE. &
            (Species(iSpec)%Init(iInit)%maxParticleNumberX * Species(iSpec)%Init(iInit)%maxParticleNumberY &
            * Species(iSpec)%Init(iInit)%maxParticleNumberZ)) THEN
         SWRITE(*,*) 'for species ',iSpec,' does not match number of particles in each direction!'
         CALL ABORT(__STAMP__,'ERROR: Number of particles in init / emission region',iInit)
       END IF
       xlen = abs(GEO%xmaxglob  - GEO%xminglob)
       ylen = abs(GEO%ymaxglob  - GEO%yminglob)
       zlen = abs(GEO%zmaxglob  - GEO%zminglob)
       xCoords(1:3,1) = (/GEO%xminglob,GEO%yminglob,GEO%zminglob/)
       xCoords(1:3,2) = xCoords(1:3,1) + (/xlen,0.,0./)
       xCoords(1:3,3) = xCoords(1:3,1) + (/0.,ylen,0./)
       xCoords(1:3,4) = xCoords(1:3,1) + (/xlen,ylen,0./)
       xCoords(1:3,5) = xCoords(1:3,1) + (/0.,0.,zlen/)
       xCoords(1:3,6) = xCoords(1:3,5) + (/xlen,0.,0./)
       xCoords(1:3,7) = xCoords(1:3,5) + (/0.,ylen,0./)
       xCoords(1:3,8) = xCoords(1:3,5) + (/xlen,ylen,0./)
       RegionOnProc=BoxInProc(xCoords,8)
    CASE ('IMD')
       RegionOnProc=.TRUE.
    CASE ('background')
       RegionOnProc=.TRUE.
    CASE ('EmissionDistribution')
       RegionOnProc=.TRUE.
    CASE DEFAULT
      IPWRITE(*,*) 'ERROR: Species ', iSpec, 'of', iInit, 'is using an unknown SpaceIC!'
      CALL ABORT(__STAMP__,'ERROR: Given SpaceIC is not implemented: '//TRIM(Species(iSpec)%Init(iInit)%SpaceIC))
    END SELECT

    ! Sanity check if at least one proc will be on the new emission communicator
    CALL MPI_ALLREDUCE(RegionOnProc,RegionExists,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_PICLAS,iError)
    IF (.NOT. RegionExists) THEN
      WRITE(hilf,'(A,I0,A,I0)') 'Species',iSpec,'-Init',iInit
      CALL CollectiveStop(__STAMP__,'The emission region was not found on any processor.  No processor in range for '//TRIM(hilf))
    END IF

    ! Add MPIRoot to specific inits automatically for output of analysis data to disk
    ! The root sometimes also reads data during restart and broadcasts it to the other processors in the communicator
    SELECT CASE(TRIM(Species(iSpec)%Init(iInit)%SpaceIC))
    CASE('2D_landmark_neutralization','2D_Liu2010_neutralization','3D_Liu2010_neutralization','2D_Liu2010_neutralization_Szabo',&
         '3D_Liu2010_neutralization_Szabo')
      IF(MPIRoot) RegionOnProc=.TRUE.
    END SELECT

    ! create new communicator
    color = MERGE(nInitRegions,MPI_UNDEFINED,RegionOnProc)

    ! set communicator id
    Species(iSpec)%Init(iInit)%InitCOMM=nInitRegions
    ! create new emission communicator for emission communication. Pass MPI_INFO_NULL as rank to follow the original ordering
    CALL MPI_COMM_SPLIT(MPI_COMM_PICLAS,color,MPI_INFO_NULL,PartMPIInitGroup(nInitRegions)%COMM,iError)

    ! Find my rank on the shared communicator, comm size and proc name
    IF (RegionOnProc) THEN
      CALL MPI_COMM_RANK(PartMPIInitGroup(nInitRegions)%COMM,PartMPIInitGroup(nInitRegions)%MyRank,iError)
      CALL MPI_COMM_SIZE(PartMPIInitGroup(nInitRegions)%COMM,PartMPIInitGroup(nInitRegions)%nProcs,iError)

      ! inform about size of emission communicator
      IF (PartMPIInitGroup(nInitRegions)%MyRank.EQ.0) THEN
#if USE_LOADBALANCE
        IF(.NOT.PerformLoadBalance)&
#endif /*USE_LOADBALANCE*/
            WRITE(UNIT_StdOut,'(A,I0,A,I0,A,I0,A)') ' Emission-Region,Emission-Communicator: ',nInitRegions,' on ',&
      PartMPIInitGroup(nInitRegions)%nProcs,' procs ('//TRIM(Species(iSpec)%Init(iInit)%SpaceIC)//', iSpec=',iSpec,')'
      END IF
    END IF

    ! build mapping for procs on emission communicator
    IF(PartMPIInitGroup(nInitRegions)%COMM.NE.MPI_COMM_NULL) THEN
      PartMPIInitGroup(nInitRegions)%MPIRoot=MERGE(.TRUE.,.FALSE.,PartMPIInitGroup(nInitRegions)%MyRank.EQ.0)

      ALLOCATE(PartMPIInitGroup(nInitRegions)%GroupToComm(0:PartMPIInitGroup(nInitRegions)%nProcs-1))
      PartMPIInitGroup(nInitRegions)%GroupToComm(PartMPIInitGroup(nInitRegions)%MyRank) = myRank
      CALL MPI_ALLGATHER(myRank,1,MPI_INTEGER&
                        ,PartMPIInitGroup(nInitRegions)%GroupToComm(0:PartMPIInitGroup(nInitRegions)%nProcs-1)&
                       ,1,MPI_INTEGER,PartMPIInitGroup(nInitRegions)%COMM,iERROR)

      ALLOCATE(PartMPIInitGroup(nInitRegions)%CommToGroup(0:nProcessors-1))
      PartMPIInitGroup(nInitRegions)%CommToGroup(0:nProcessors-1) = -1
      DO iRank = 0,PartMPIInitGroup(nInitRegions)%nProcs-1
        PartMPIInitGroup(nInitRegions)%CommToGroup(PartMPIInitGroup(nInitRegions)%GroupToComm(iRank))=iRank
      END DO ! iRank

    END IF
  END DO ! iniT
END DO ! iSpec

END SUBROUTINE InitEmissionComm


SUBROUTINE SendEmissionParticlesToProcs(chunkSize,DimSend,FractNbr,iInit,mySumOfMatchedParticles,particle_positions)
!----------------------------------------------------------------------------------------------------------------------------------!
! A particle's host cell in the FIBGM is found and the corresponding procs are notified.
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
USE MOD_Particle_Localization  ,ONLY: LocateParticleInElement,SinglePointToElement
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemInfo_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: FIBGMToProc,FIBGMProcs
!USE MOD_Mesh_Tools             ,ONLY: GetCNElemID
!USE MOD_Particle_Mesh_Vars     ,ONLY: FIBGM_nElems, FIBGM_offsetElem, FIBGM_Element
USE MOD_Particle_Mesh_Vars     ,ONLY: FIBGM_nElems,FIBGM_nTotalElems
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPIInitGroup,PartMPIInsert,PartMPILocate
USE MOD_Particle_MPI_Vars      ,ONLY: EmissionSendBuf,EmissionRecvBuf
USE MOD_Particle_Vars          ,ONLY: PDM,PEM,PartState,PartPosRef,Species
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
USE MOD_Part_Tools             ,ONLY: GetNextFreePosition
#if defined(MEASURE_MPI_WAIT)
USE MOD_Particle_MPI_Vars      ,ONLY: MPIW8TimePart,MPIW8CountPart
#endif /*defined(MEASURE_MPI_WAIT)*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: chunkSize
INTEGER,INTENT(IN)            :: DimSend
INTEGER,INTENT(IN)            :: FractNbr
INTEGER,INTENT(IN)            :: iInit
REAL,INTENT(IN),OPTIONAL      :: particle_positions(1:chunkSize*DimSend)
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)           :: mySumOfMatchedParticles
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Counters
INTEGER                       :: i,iPos,iProc,iDir,ElemID,ProcID
! BGM
INTEGER                       :: ijkBGM(3,chunkSize)
INTEGER                       :: TotalNbrOfRecvParts!,iBGMElem,nBGMElems
LOGICAL                       :: InsideMyBGM(2,chunkSize)
! Temporary state arrays
REAL,ALLOCATABLE              :: chunkState(:,:)
! MPI Communication
INTEGER                       :: ALLOCSTAT,PartCommSize,ParticleIndexNbr
INTEGER                       :: InitGroup,tProc
INTEGER                       :: msg_status(1:MPI_STATUS_SIZE),messageSize
INTEGER                       :: nRecvParticles,nSendParticles
REAL,ALLOCATABLE              :: recvPartPos(:)
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=8)               :: CounterStart,CounterEnd
REAL(KIND=8)                  :: Rate
#endif /*defined(MEASURE_MPI_WAIT)*/
!===================================================================================================================================
InitGroup = Species(FractNbr)%Init(iInit)%InitCOMM

! Arrays for communication of particles not located in final element
ALLOCATE( PartMPIInsert%nPartsSend  (2,0:PartMPIInitGroup(InitGroup)%nProcs-1) &
        , PartMPIInsert%nPartsRecv  (1,0:PartMPIInitGroup(InitGroup)%nProcs-1) &
        , PartMPIInsert%SendRequest (2,0:PartMPIInitGroup(InitGroup)%nProcs-1) &
        , PartMPIInsert%RecvRequest (2,0:PartMPIInitGroup(InitGroup)%nProcs-1) &
        , PartMPIInsert%send_message(  0:PartMPIInitGroup(InitGroup)%nProcs-1) &
        , STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) &
  CALL ABORT(__STAMP__,' Cannot allocate particle emission MPI arrays! ALLOCSTAT',ALLOCSTAT)

PartMPIInsert%nPartsSend=0
PartMPIInsert%nPartsRecv=0

! Inter-CN communication
ALLOCATE( PartMPILocate%nPartsSend (2,0:PartMPIInitGroup(InitGroup)%nProcs-1) &
        , PartMPILocate%nPartsRecv (1,0:PartMPIInitGroup(InitGroup)%nProcs-1) &
        , PartMPILocate%SendRequest(2,0:PartMPIInitGroup(InitGroup)%nProcs-1) &
        , PartMPILocate%RecvRequest(2,0:PartMPIInitGroup(InitGroup)%nProcs-1) &
        , EmissionRecvBuf          (  0:PartMPIInitGroup(InitGroup)%nProcs-1) &
        , EmissionSendBuf          (  0:PartMPIInitGroup(InitGroup)%nProcs-1) &
        , STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) &
  CALL ABORT(__STAMP__,' Cannot allocate particle emission MPI arrays! ALLOCSTAT',ALLOCSTAT)

PartMPILocate%nPartsSend=0
PartMPILocate%nPartsRecv=0

! Arrays for communication of particles located in final element. Reuse particle_mpi infrastructure wherever possible
PartCommSize   = 0
PartCommSize   = PartCommSize + 3                              ! Emission position (physical space)
IF(TrackingMethod.EQ.REFMAPPING) PartCommSize = PartCommSize+3 ! Emission position (reference space)
!PartCommSize   = PartCommSize + 1                             ! Species-ID
PartCommSize   = PartCommSize + 1                              ! ID of element

! Temporary array to hold ElemID of located particles
ALLOCATE( chunkState(PartCommSize,chunkSize)                                                &
        , STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) &
  CALL ABORT(__STAMP__,' Cannot allocate particle emission MPI arrays! ALLOCSTAT',ALLOCSTAT)

chunkState = -1

!--- 1/10 Open receive buffer (located and non-located particles)
DO iProc=0,PartMPIInitGroup(InitGroup)%nProcs-1
  IF (iProc.EQ.PartMPIInitGroup(InitGroup)%myRank) CYCLE

  !--- MPI_IRECV lengths of lists of particles entering local mesh
  CALL MPI_IRECV( PartMPIInsert%nPartsRecv(:,iProc)                           &
                , 1                                                           &
                , MPI_INTEGER                                                 &
                , iProc                                                       &
                , 1011                                                        &
                , PartMPIInitGroup(InitGroup)%COMM                           &
                , PartMPIInsert%RecvRequest(1,iProc)                          &
                , IERROR)

  ! Inter-CN communication
  CALL MPI_IRECV( PartMPILocate%nPartsRecv(:,iProc)                           &
                , 1                                                           &
                , MPI_INTEGER                                                 &
                , iProc                                                       &
                , 1111                                                        &
                , PartMPIInitGroup(InitGroup)%COMM                           &
                , PartMPILocate%RecvRequest(1,iProc)                          &
                , IERROR)
END DO

! Identify particles that are on the node (or in the halo region of the node) or on other nodes
DO i=1,chunkSize
  ! Set BGM cell index
  ASSOCIATE( xMin   => (/GEO%xminglob  , GEO%yminglob  , GEO%zminglob/)  ,    &
             BGMMin => (/GEO%FIBGMimin , GEO%FIBGMjmin , GEO%FIBGMkmin/) ,    &
             BGMMax => (/GEO%FIBGMimax , GEO%FIBGMjmax , GEO%FIBGMkmax/) )
    DO iDir = 1, 3
      ijkBGM(iDir,i) = INT((particle_positions(DimSend*(i-1)+iDir)-xMin(iDir))/GEO%FIBGMdeltas(iDir))+1
    END DO ! iDir = 1, 3

    ! Check BGM cell index
    InsideMyBGM(:,i)=.TRUE.
    DO iDir = 1, 3
      IF(ijkBGM(iDir,i).LT.BGMMin(iDir)) THEN
        InsideMyBGM(1,i)=.FALSE.
        EXIT
      END IF
      IF(ijkBGM(iDir,i).GT.BGMMax(iDir)) THEN
        InsideMyBGM(1,i)=.FALSE.
        EXIT
      END IF
    END DO ! iDir = 1, 3
  END ASSOCIATE

  IF (InsideMyBGM(1,i)) THEN
    ASSOCIATE(iBGM => ijkBGM(1,i), &
              jBGM => ijkBGM(2,i), &
              kBGM => ijkBGM(3,i))

    !--- check if BGM cell contains elements not within the halo region
    IF (FIBGM_nElems(iBGM,jBGM,kBGM).NE.FIBGM_nTotalElems(iBGM,jBGM,kBGM)) THEN
      !--- if any elements are found, communicate particle to all procs
      InsideMyBGM(2,i) = .FALSE.
    END IF ! (FIBGM_nElems(iBGM,jBGM,kBGM).NE.FIBGM_nTotalElems(iBGM,jBGM,kBGM))
    END ASSOCIATE
  END IF ! InsideMyBGM(i)
END DO ! i = 1, chunkSize

!--- Find non-local particles for sending to other nodes
DO i = 1, chunkSize
  IF(ANY(.NOT.InsideMyBGM(:,i))) THEN
    ! Inter-CN communication
    ASSOCIATE(iBGM => ijkBGM(1,i), &
              jBGM => ijkBGM(2,i), &
              kBGM => ijkBGM(3,i))

    ! Sanity check if the emission is within the global FIBGM region
    IF (iBGM.LT.GEO%FIBGMiminglob .OR. iBGM.GT.GEO%FIBGMimaxglob .OR. &
        jBGM.LT.GEO%FIBGMjminglob .OR. jBGM.GT.GEO%FIBGMjmaxglob .OR. &
        kBGM.LT.GEO%FIBGMkminglob .OR. kBGM.GT.GEO%FIBGMkmaxglob) THEN
      CYCLE
    END IF

    !-- Find all procs associated with the background mesh cell. Then loop over all procs and count number of particles per proc for
    !-- sending
    DO iProc = FIBGMToProc(FIBGM_FIRSTPROCIND,iBGM,jBGM,kBGM)+1, &
               FIBGMToProc(FIBGM_FIRSTPROCIND,iBGM,jBGM,kBGM)+FIBGMToProc(FIBGM_NPROCS,iBGM,jBGM,kBGM)
      ProcID = FIBGMProcs(iProc)
      IF (ProcID.EQ.myRank) CYCLE

      tProc=PartMPIInitGroup(InitGroup)%CommToGroup(ProcID)
      ! Processor is not on emission communicator
      IF(tProc.EQ.-1) CYCLE

      PartMPIInsert%nPartsSend(1,tProc) = PartMPIInsert%nPartsSend(1,tProc)+1
    END DO

    END ASSOCIATE
  END IF ! .NOT.InsideMyBGM(i)
END DO ! i = 1, chunkSize

!--- 2/10 Send number of non-located particles
DO iProc=0,PartMPIInitGroup(InitGroup)%nProcs-1
  IF (iProc.EQ.PartMPIInitGroup(InitGroup)%myRank) CYCLE

  ! send particles
  !--- MPI_ISEND lengths of lists of particles leaving local mesh
  CALL MPI_ISEND( PartMPIInsert%nPartsSend( 1,iProc)                          &
                , 1                                                           &
                , MPI_INTEGER                                                 &
                , iProc                                                       &
                , 1011                                                        &
                , PartMPIInitGroup(InitGroup)%COMM                           &
                , PartMPIInsert%SendRequest(1,iProc)                          &
                , IERROR)
  IF (PartMPIInsert%nPartsSend(1,iProc).GT.0) THEN
    MessageSize = DimSend*PartMPIInsert%nPartsSend(1,iProc)
    ALLOCATE( PartMPIInsert%send_message(iProc)%content(MessageSize), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) &
      CALL ABORT(__STAMP__,'  Cannot allocate emission PartSendBuf, local ProcId, ALLOCSTAT',iProc,REAL(ALLOCSTAT))
  END IF
END DO


!--- 3/10 Send actual non-located particles
PartMPIInsert%nPartsSend(2,:)=0
DO i = 1, chunkSize
  IF(ANY(.NOT.InsideMyBGM(:,i))) THEN
    ! Inter-CN communication
    ASSOCIATE(iBGM => ijkBGM(1,i), &
          jBGM => ijkBGM(2,i), &
          kBGM => ijkBGM(3,i))

    ! Sanity check if the emission is within the global FIBGM region
    IF (iBGM.LT.GEO%FIBGMiminglob .OR. iBGM.GT.GEO%FIBGMimaxglob .OR. &
        jBGM.LT.GEO%FIBGMjminglob .OR. jBGM.GT.GEO%FIBGMjmaxglob .OR. &
        kBGM.LT.GEO%FIBGMkminglob .OR. kBGM.GT.GEO%FIBGMkmaxglob) THEN
      CYCLE
    END IF

    !-- Find all procs associated with the background mesh cell. Then loop over all procs and count number of particles per proc for
    !-- sending
    DO iProc = FIBGMToProc(FIBGM_FIRSTPROCIND,iBGM,jBGM,kBGM)+1, &
               FIBGMToProc(FIBGM_FIRSTPROCIND,iBGM,jBGM,kBGM)+FIBGMToProc(FIBGM_NPROCS,iBGM,jBGM,kBGM)
      ProcID = FIBGMProcs(iProc)
      IF (ProcID.EQ.myRank) CYCLE

      tProc=PartMPIInitGroup(InitGroup)%CommToGroup(ProcID)
      ! Processor is not on emission communicator
      IF (tProc.EQ.-1) CYCLE

      ! Assemble message
      iPos = PartMPIInsert%nPartsSend(2,tProc) * DimSend
      PartMPIInsert%send_message(tProc)%content(iPos+1:iPos+3) = particle_positions(DimSend*(i-1)+1:DimSend*i)

      ! Counter of previous particles on proc
      PartMPIInsert%nPartsSend(2,tProc)=PartMPIInsert%nPartsSend(2,tProc) + 1
    END DO

    END ASSOCIATE
  END IF ! .NOT.InsideMyBGM(i)
END DO ! i = 1, chunkSize


!--- 4/10 Receive actual non-located particles
#if defined(MEASURE_MPI_WAIT)
CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/
DO iProc=0,PartMPIInitGroup(InitGroup)%nProcs-1
  IF (iProc.EQ.PartMPIInitGroup(InitGroup)%myRank) CYCLE

  CALL MPI_WAIT(PartMPIInsert%SendRequest(1,iProc),msg_status(:),IERROR)
  IF (IERROR.NE.MPI_SUCCESS) CALL abort(__STAMP__,' MPI Communication error', IERROR)
  CALL MPI_WAIT(PartMPIInsert%RecvRequest(1,iProc),msg_status(:),IERROR)
  IF (IERROR.NE.MPI_SUCCESS) CALL abort(__STAMP__,' MPI Communication error', IERROR)
END DO
#if defined(MEASURE_MPI_WAIT)
CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
MPIW8TimePart(5)  = MPIW8TimePart(5) + REAL(CounterEnd-CounterStart,8)/Rate
MPIW8CountPart(5) = MPIW8CountPart(5) + 1_8
#endif /*defined(MEASURE_MPI_WAIT)*/

! recvPartPos holds particles from ALL procs
! Inter-CN communication
ALLOCATE(recvPartPos(1:SUM(PartMPIInsert%nPartsRecv(1,:)*DimSend)), STAT=ALLOCSTAT)
TotalNbrOfRecvParts = 0
DO iProc=0,PartMPIInitGroup(InitGroup)%nProcs-1
  IF (iProc.EQ.PartMPIInitGroup(InitGroup)%myRank) CYCLE

  IF (PartMPIInsert%nPartsRecv(1,iProc).GT.0) THEN
  !--- MPI_IRECV lengths of lists of particles entering local mesh
    CALL MPI_IRECV( recvPartPos(TotalNbrOfRecvParts*DimSend+1)                &
                  , DimSend*PartMPIInsert%nPartsRecv(1,iProc)                 &
                  , MPI_DOUBLE_PRECISION                                      &
                  , iProc                                                     &
                  , 1022                                                      &
                  , PartMPIInitGroup(InitGroup)%COMM                         &
                  , PartMPIInsert%RecvRequest(2,iProc)                        &
                  , IERROR)
    TotalNbrOfRecvParts = TotalNbrOfRecvParts + PartMPIInsert%nPartsRecv(1,iProc)
  END IF
  !--- (non-blocking:) send messages to all procs receiving particles from myself
  IF (PartMPIInsert%nPartsSend(2,iProc).GT.0) THEN
    CALL MPI_ISEND( PartMPIInsert%send_message(iProc)%content                 &
                  , DimSend*PartMPIInsert%nPartsSend(2,iProc)                 &
                  , MPI_DOUBLE_PRECISION                                      &
                  , iProc                                                     &
                  , 1022                                                      &
                  , PartMPIInitGroup(InitGroup)%COMM                         &
                  , PartMPIInsert%SendRequest(2,iProc)                        &
                  , IERROR)
  END IF
END DO

mySumOfMatchedParticles = 0
ParticleIndexNbr        = 1

!--- Locate local (node or halo of node) particles
DO i = 1, chunkSize
  IF(InsideMyBGM(1,i))THEN
    ! We cannot call LocateParticleInElement because we do not know the final PartID yet. Locate the position and fill PartState
    ! manually if we got a hit
    ElemID = SinglePointToElement(particle_positions(DimSend*(i-1)+1:DimSend*(i-1)+3),doHALO=.TRUE.,doEmission_opt=.TRUE.)
    ! Checked every possible cell and didn't find it. Apparently, we emitted the particle outside the domain
    IF(ElemID.EQ.-1) CYCLE

    ! Only keep the particle if it belongs on the current proc. Otherwise prepare to send it to the correct proc
    ! TODO: Implement U_Shared, so we can finish emission on this proc and send the fully initialized particle (i.e. including
    ! velocity)
    ProcID = ElemInfo_Shared(ELEM_RANK,ElemID)
    IF (ProcID.NE.myRank) THEN
      ! Particle was sent to every potential proc, so trust the other proc to find it and do not send it again
      IF (.NOT.InsideMyBGM(2,i)) CYCLE

      ! ProcID on emission communicator
      tProc=PartMPIInitGroup(InitGroup)%CommToGroup(ProcID)
      ! Processor is not on emission communicator
      IF(tProc.EQ.-1) CALL ABORT(__STAMP__,'Error in particle_mpi_emission: proc not on emission communicator')

      PartMPILocate%nPartsSend(1,tProc)= PartMPILocate%nPartsSend(1,tProc)+1

      ! Assemble temporary PartState to send the final particle position
      chunkState(1:3,i) = particle_positions(DimSend*(i-1)+1:DimSend*(i-1)+3)
      IF (TrackingMethod.EQ.REFMAPPING) THEN
        CALL GetPositionInRefElem(chunkState(1:3,i),chunkState(4:6,i),ElemID)
        chunkState(7,i) = REAL(ElemID,KIND=8)
      ELSE
        chunkState(4,i) = REAL(ElemID,KIND=8)
      END IF ! TrackingMethod.EQ.REFMAPPING
    ! Located particle on local proc.
    ELSE
      ! Get the next free position in the PDM array
      ParticleIndexNbr = GetNextFreePosition(mySumOfMatchedParticles+1)
      ! Fill the PartState manually to avoid a second localization
      PartState(1:DimSend,ParticleIndexNbr) = particle_positions(DimSend*(i-1)+1:DimSend*(i-1)+DimSend)
      PDM%ParticleInside( ParticleIndexNbr) = .TRUE.
      IF (TrackingMethod.EQ.REFMAPPING) THEN
        CALL GetPositionInRefElem(PartState(1:3,ParticleIndexNbr),PartPosRef(1:3,ParticleIndexNbr),ElemID)
      END IF ! TrackingMethod.EQ.REFMAPPING
      PEM%GlobalElemID(ParticleIndexNbr)         = ElemID
      mySumOfMatchedParticles = mySumOfMatchedParticles + 1
    END IF ! ElemID.EQ.-1
  END IF ! InsideMyBGM(i)
END DO ! i = 1, chunkSize

!---  /  Send number of located particles
! Inter-CN communication
DO iProc=0,PartMPIInitGroup(InitGroup)%nProcs-1
  IF (iProc.EQ.PartMPIInitGroup(InitGroup)%myRank) CYCLE

  ! send particles
  !--- MPI_ISEND lengths of lists of particles leaving local mesh
  CALL MPI_ISEND( PartMPILocate%nPartsSend( 1,iProc)                          &
                , 1                                                           &
                , MPI_INTEGER                                                 &
                , iProc                                                       &
                , 1111                                                        &
                , PartMPIInitGroup(InitGroup)%COMM                           &
                , PartMPILocate%SendRequest(1,iProc)                          &
                , IERROR)
  IF (PartMPILocate%nPartsSend(1,iProc).GT.0) THEN
    MessageSize = PartMPILocate%nPartsSend(1,iProc)*PartCommSize
    ALLOCATE(EmissionSendBuf(iProc)%content(MessageSize),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) &
      CALL ABORT(__STAMP__,'  Cannot allocate emission EmissionSendBuf, local ProcId, ALLOCSTAT',iProc,REAL(ALLOCSTAT))
  END IF
END DO

!--- 5/10 Send actual located particles. PartState is filled in LocateParticleInElement
PartMPILocate%nPartsSend(2,:) = 0
DO i = 1, chunkSize
  ElemID = INT(chunkState(PartCommSize,i))
  ! Skip non-located particles
  IF(ElemID.EQ.-1) CYCLE
  ProcID = ElemInfo_Shared(ELEM_RANK,ElemID)
  IF (ProcID.NE.myRank) THEN
    ! ProcID on emission communicator
    tProc=PartMPIInitGroup(InitGroup)%CommToGroup(ProcID)
    ! Processor is not on emission communicator
    IF(tProc.EQ.-1) CYCLE

    ! Assemble message
    iPos = PartMPILocate%nPartsSend(2,tProc) * PartCommSize
    EmissionSendBuf(tProc)%content(1+iPos:PartCommSize+iPos) = chunkState(1:PartCommSize,i)

    ! Counter of previous particles on proc
    PartMPILocate%nPartsSend(2,tProc) = PartMPILocate%nPartsSend(2,tProc) + 1
  END IF ! ProcID.NE.myRank
END DO ! i = 1, chunkSize

!--- 6/10 Receive actual non-located particles
#if defined(MEASURE_MPI_WAIT)
CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/
DO iProc=0,PartMPIInitGroup(InitGroup)%nProcs-1
  IF (iProc.EQ.PartMPIInitGroup(InitGroup)%myRank) CYCLE

  CALL MPI_WAIT(PartMPILocate%SendRequest(1,iProc),msg_status(:),IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL abort(__STAMP__,' MPI Communication error', IERROR)
  CALL MPI_WAIT(PartMPILocate%RecvRequest(1,iProc),msg_status(:),IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL abort(__STAMP__,' MPI Communication error', IERROR)
END DO
#if defined(MEASURE_MPI_WAIT)
CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
MPIW8TimePart(5)  = MPIW8TimePart(5) + REAL(CounterEnd-CounterStart,8)/Rate
MPIW8CountPart(5) = MPIW8CountPart(5) + 1_8
#endif /*defined(MEASURE_MPI_WAIT)*/

DO iProc=0,PartMPIInitGroup(InitGroup)%nProcs-1
  IF (iProc.EQ.PartMPIInitGroup(InitGroup)%myRank) CYCLE

  ! Allocate receive array and open receive buffer if expecting particles from iProc
  IF (PartMPILocate%nPartsRecv(1,iProc).GT.0) THEN
    nRecvParticles = PartMPILocate%nPartsRecv(1,iProc)
    MessageSize    = nRecvParticles * PartCommSize
    ALLOCATE(EmissionRecvBuf(iProc)%content(MessageSize),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) &
      CALL ABORT(__STAMP__,'  Cannot allocate emission EmissionRecvBuf, local ProcId, ALLOCSTAT',iProc,REAL(ALLOCSTAT))

    !--- MPI_IRECV lengths of lists of particles entering local mesh
    CALL MPI_IRECV( EmissionRecvBuf(iProc)%content                             &
                  , MessageSize                                                &
                  , MPI_DOUBLE_PRECISION                                       &
                  , iProc                                                      &
                  , 1122                                                       &
                  , PartMPIInitGroup(InitGroup)%COMM                          &
                  , PartMPILocate%RecvRequest(2,iProc)                         &
                  , IERROR )
    IF(IERROR.NE.MPI_SUCCESS) &
      CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  END IF
  !--- (non-blocking:) send messages to all procs receiving particles from myself
  IF (PartMPILocate%nPartsSend(2,iProc).GT.0) THEN
    nSendParticles = PartMPILocate%nPartsSend(1,iProc)
    MessageSize    = nSendParticles * PartCommSize
    CALL MPI_ISEND( EmissionSendBuf(iProc)%content                             &
                  , MessageSize                                                &
                  , MPI_DOUBLE_PRECISION                                       &
                  , iProc                                                      &
                  , 1122                                                       &
                  , PartMPIInitGroup(InitGroup)%COMM                          &
                  , PartMPILocate%SendRequest(2,iProc)                         &
                  , IERROR )
    IF(IERROR.NE.MPI_SUCCESS) &
      CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  END IF
END DO

!--- 7/10 Finish communication of actual non-located particles
#if defined(MEASURE_MPI_WAIT)
CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/
DO iProc=0,PartMPIInitGroup(InitGroup)%nProcs-1
  IF (iProc.EQ.PartMPIInitGroup(InitGroup)%myRank) CYCLE

  IF (PartMPIInsert%nPartsSend(1,iProc).GT.0) THEN
    CALL MPI_WAIT(PartMPIInsert%SendRequest(2,iProc),msg_status(:),IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(__STAMP__,' MPI Communication error', IERROR)
  END IF
  IF (PartMPIInsert%nPartsRecv(1,iProc).GT.0) THEN
    CALL MPI_WAIT(PartMPIInsert%RecvRequest(2,iProc),msg_status(:),IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(__STAMP__,' MPI Communication error', IERROR)
  END IF
END DO
#if defined(MEASURE_MPI_WAIT)
CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
MPIW8TimePart(5)  = MPIW8TimePart(5) + REAL(CounterEnd-CounterStart,8)/Rate
MPIW8CountPart(5) = MPIW8CountPart(5) + 1_8
#endif /*defined(MEASURE_MPI_WAIT)*/

!--- 8/10 Try to locate received non-located particles
TotalNbrOfRecvParts = SUM(PartMPIInsert%nPartsRecv(1,:))
DO i = 1,TotalNbrOfRecvParts
  ! We cannot call LocateParticleInElement because we do not know the final PartID yet. Locate the position and fill PartState
  ! manually if we got a hit
  ElemID = SinglePointToElement(recvPartPos(DimSend*(i-1)+1:DimSend*(i-1)+3),doHALO=.FALSE.,doEmission_opt=.TRUE.)
  ! Checked every possible cell and didn't find it. Apparently, we emitted the particle outside the domain
  IF(ElemID.EQ.-1) CYCLE

  ! Only keep the particle if it belongs on the current proc. Trust the other procs to do their jobs and locate it if needed
  IF (ElemInfo_Shared(ELEM_RANK,ElemID).NE.myRank) CYCLE

  ! Find a free position in the PDM array
  ParticleIndexNbr = GetNextFreePosition(mySumOfMatchedParticles+1)
  ! Fill the PartState manually to avoid a second localization
  PartState(1:3,ParticleIndexNbr) = recvPartPos(DimSend*(i-1)+1:DimSend*(i-1)+3)
  PDM%ParticleInside( ParticleIndexNbr) = .TRUE.
  IF (TrackingMethod.EQ.REFMAPPING) THEN
    CALL GetPositionInRefElem(PartState(1:3,ParticleIndexNbr),PartPosRef(1:3,ParticleIndexNbr),ElemID)
  END IF ! TrackingMethod.EQ.REFMAPPING
  PEM%GlobalElemID(ParticleIndexNbr)    = ElemID
  mySumOfMatchedParticles = mySumOfMatchedParticles + 1
END DO

!--- 9/10 Finish communication of actual non-located particles
#if defined(MEASURE_MPI_WAIT)
CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/
DO iProc=0,PartMPIInitGroup(InitGroup)%nProcs-1
  IF (iProc.EQ.PartMPIInitGroup(InitGroup)%myRank) CYCLE

  IF (PartMPILocate%nPartsSend(1,iProc).GT.0) THEN
    CALL MPI_WAIT(PartMPILocate%SendRequest(2,iProc),msg_status(:),IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(__STAMP__,' MPI Communication error', IERROR)
  END IF
  IF (PartMPILocate%nPartsRecv(1,iProc).GT.0) THEN
    CALL MPI_WAIT(PartMPILocate%RecvRequest(2,iProc),msg_status(:),IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(__STAMP__,' MPI Communication error', IERROR)
  END IF
END DO
#if defined(MEASURE_MPI_WAIT)
CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
MPIW8TimePart(5)  = MPIW8TimePart(5) + REAL(CounterEnd-CounterStart,8)/Rate
MPIW8CountPart(5) = MPIW8CountPart(5) + 1_8
#endif /*defined(MEASURE_MPI_WAIT)*/

!--- 10/10 Write located particles
DO iProc=0,PartMPIInitGroup(InitGroup)%nProcs-1
  IF (iProc.EQ.PartMPIInitGroup(InitGroup)%myRank) CYCLE
  IF (PartMPILocate%nPartsRecv(1,iProc).EQ.0) CYCLE

  DO i = 1,PartMPILocate%nPartsRecv(1,iProc)
    ! Find a free position in the PDM array
    ParticleIndexNbr = GetNextFreePosition(mySumOfMatchedParticles+1)
    ! Fill the PartState manually to avoid a second localization
    PartState(1:3,ParticleIndexNbr) = EmissionRecvBuf(iProc)%content(PartCommSize*(i-1)+1:PartCommSize*(i-1)+3)
    IF (TrackingMethod.EQ.REFMAPPING) THEN
      PartPosRef(1:3,ParticleIndexNbr) = EmissionRecvBuf(iProc)%content(PartCommSize*(i-1)+4:PartCommSize*(i-1)+6)
    END IF ! TrackingMethod.EQ.REFMAPPING
    PEM%GlobalElemID(ParticleIndexNbr)    = INT(EmissionRecvBuf(iProc)%content(PartCommSize*(i)),KIND=4)
    PDM%ParticleInside( ParticleIndexNbr) = .TRUE.
    mySumOfMatchedParticles = mySumOfMatchedParticles + 1
  END DO
END DO

!--- Clean up
SDEALLOCATE(recvPartPos)
SDEALLOCATE(chunkState)
DO iProc=0,PartMPIInitGroup(InitGroup)%nProcs-1
  SDEALLOCATE(EmissionRecvBuf(iProc)%content)
  SDEALLOCATE(EmissionSendBuf(iProc)%content)
END DO
SDEALLOCATE(PartMPIInsert%nPartsSend)
SDEALLOCATE(PartMPIInsert%nPartsRecv)
SDEALLOCATE(PartMPIInsert%SendRequest)
SDEALLOCATE(PartMPIInsert%RecvRequest)
SDEALLOCATE(PartMPIInsert%send_message)

SDEALLOCATE(PartMPILocate%nPartsSend)
SDEALLOCATE(PartMPILocate%nPartsRecv)
SDEALLOCATE(PartMPILocate%SendRequest)
SDEALLOCATE(PartMPILocate%RecvRequest)
SDEALLOCATE(EmissionRecvBuf)
SDEALLOCATE(EmissionSendBuf)

END SUBROUTINE SendEmissionParticlesToProcs


!===================================================================================================================================
!> Check if bounding box is on processor. The bounding box is built from the min/max extents of the input nodes.
!> The number of input nodes is nNodes can be any integer number > 1.
!> The bounding box is compared with the GEO%xmin, GEO%xmax, ... etc. of each processor.
!===================================================================================================================================
PURE FUNCTION BoxInProc(CartNodes,nNodes)
! MODULES
USE MOD_Particle_Mesh_Vars,       ONLY:GEO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN):: nNodes
REAL,INTENT(IN)   :: CartNodes(1:3,1:nNodes)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL           :: BoxInProc
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!INTEGER           :: xmin,xmax,ymin,ymax,zmin,zmax,testval
REAL,DIMENSION(6)  :: xCoords
!===================================================================================================================================

BoxInProc=.FALSE.
!! get background of nodes
!xmin = HUGE(1)
!xmax =-HUGE(1)
!ymin = HUGE(1)
!ymax =-HUGE(1)
!zmin = HUGE(1)
!zmax =-HUGE(1)
!
!testval = FLOOR((MINVAL(CartNodes(1,:)) - GEO%xminglob)/GEO%FIBGMdeltas(1)) + 1
!xmin    = MIN(xmin,testval)
!testval = FLOOR((MAXVAL(CartNodes(1,:)) - GEO%xminglob)/GEO%FIBGMdeltas(1)) + 1
!xmax    = MAX(xmax,testval)
!testval = FLOOR((MINVAL(CartNodes(2,:)) - GEO%yminglob)/GEO%FIBGMdeltas(2)) + 1
!ymin    = MIN(ymin,testval)
!testval = FLOOR((MAXVAL(CartNodes(2,:)) - GEO%yminglob)/GEO%FIBGMdeltas(2)) + 1
!ymax    = MAX(ymax,testval)
!testval = FLOOR((MINVAL(CartNodes(3,:)) - GEO%zminglob)/GEO%FIBGMdeltas(3)) + 1
!zmin    = MIN(zmin,testval)
!testval = FLOOR((MAXVAL(CartNodes(3,:)) - GEO%zminglob)/GEO%FIBGMdeltas(3)) + 1
!zmax    = MAX(zmax,testval)
!
!IF(    ((xmin.LE.GEO%FIBGMimax).AND.(xmax.GE.GEO%FIBGMimin)) &
!  .AND.((ymin.LE.GEO%FIBGMjmax).AND.(ymax.GE.GEO%FIBGMjmin)) &
!  .AND.((zmin.LE.GEO%FIBGMkmax).AND.(zmax.GE.GEO%FIBGMkmin)) ) BoxInProc = .TRUE.

! Calculate directly with global coordinates
xCoords(1) = MINVAL(CartNodes(1,:))
xCoords(2) = MAXVAL(CartNodes(1,:))
xCoords(3) = MINVAL(CartNodes(2,:))
xCoords(4) = MAXVAL(CartNodes(2,:))
xCoords(5) = MINVAL(CartNodes(3,:))
xCoords(6) = MAXVAL(CartNodes(3,:))

IF(    ((xCoords(1).LE.GEO%xmax).AND.(xCoords(2).GE.GEO%xmin)) &
  .AND.((xCoords(3).LE.GEO%ymax).AND.(xCoords(4).GE.GEO%ymin)) &
  .AND.((xCoords(5).LE.GEO%zmax).AND.(xCoords(6).GE.GEO%zmin)) ) BoxInProc = .TRUE.

END FUNCTION BoxInProc


PURE FUNCTION PointInProc(CartNode)
!===================================================================================================================================
! check if point is on proc
!===================================================================================================================================
! MODULES
USE MOD_Particle_Mesh_Vars,       ONLY:GEO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)   :: CartNode(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL           :: PointInProc
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!INTEGER           :: xmin,xmax,ymin,ymax,zmin,zmax,testval
!===================================================================================================================================

PointInProc=.FALSE.
!! get background of nodes
!xmin = HUGE(1)
!xmax =-HUGE(1)
!ymin = HUGE(1)
!ymax =-HUGE(1)
!zmin = HUGE(1)
!zmax =-HUGE(1)
!
!testval = FLOOR((CartNode(1)-GEO%xminglob)/GEO%FIBGMdeltas(1)) + 1
!xmin    = MIN(xmin,testval)
!testval = FLOOR((CartNode(1)-GEO%xminglob)/GEO%FIBGMdeltas(1)) + 1
!xmax    = MAX(xmax,testval)
!testval = FLOOR((CartNode(2)-GEO%yminglob)/GEO%FIBGMdeltas(2)) + 1
!ymin    = MIN(ymin,testval)
!testval = FLOOR((CartNode(2)-GEO%yminglob)/GEO%FIBGMdeltas(2)) + 1
!ymax    = MAX(ymax,testval)
!testval = FLOOR((CartNode(3)-GEO%zminglob)/GEO%FIBGMdeltas(3)) + 1
!zmin    = MIN(zmin,testval)
!testval = FLOOR((CartNode(3)-GEO%zminglob)/GEO%FIBGMdeltas(3)) + 1
!zmax    = MAX(zmax,testval)
!
!IF(    ((xmin.LE.GEO%FIBGMimax).AND.(xmax.GE.GEO%FIBGMimin)) &
!  .AND.((ymin.LE.GEO%FIBGMjmax).AND.(ymax.GE.GEO%FIBGMjmin)) &
!  .AND.((zmin.LE.GEO%FIBGMkmax).AND.(zmax.GE.GEO%FIBGMkmin)) ) PointInProc = .TRUE.

IF(    ((CartNode(1).LE.GEO%xmax).AND.(CartNode(1).GE.GEO%xmin)) &
  .AND.((CartNode(2).LE.GEO%ymax).AND.(CartNode(2).GE.GEO%ymin)) &
  .AND.((CartNode(3).LE.GEO%zmax).AND.(CartNode(3).GE.GEO%zmin)) ) PointInProc = .TRUE.

END FUNCTION PointInProc
#endif /*USE_MPI*/

END MODULE MOD_Particle_MPI_Emission
