#include "boltzplatz.h"


MODULE                              MOD_BoundaryTools                                              !
!===================================================================================================================================
!
!===================================================================================================================================
   IMPLICIT NONE                                                                                   !
   PRIVATE                                                                                         !
!----------------------------------------------------------------------------------------------------------------------------------
! PUBLIC
   PUBLIC                        :: ParticleThroughSideCheck3DFast                                 !
   PUBLIC                        :: ParticleThroughSideCheck3DFastold                              !
   PUBLIC                        :: ParticleThroughSideLastPosCheck
   PUBLIC                        :: ParticleThroughSideCheck3D                                     !
   PUBLIC                        :: PerfectReflection3D                                            !
   PUBLIC                        :: PeriodicWallBnd3D                                              !
   PUBLIC                        :: ParticleInsideQuad3D                                           !
   PUBLIC                        :: ParticleInsideQuad3Dold                                           !
   PUBLIC                        :: SingleParticleToExactElement                                   !
   PUBLIC                        :: DiffuseReflection3D                                            !
#ifdef MPI
   PUBLIC                        :: ParticleThroughSideCheck3DFast_halocells                       !
   PUBLIC                        :: ParticleThroughSideLastPosCheck_halocells
   PUBLIC                        :: PerfectReflection3D_halocells                                  !
   PUBLIC                        :: ParticleInsideQuad3D_halocells                                 !
   PUBLIC                        :: DiffuseReflection3D_halocells
#endif 
!===================================================================================================================================

                                                                                                   !
CONTAINS                                                                                           !

SUBROUTINE ParticleThroughSideCheck3DFast(i,iLocSide,Element,ThroughSide,TriNum)                   !
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars                                                                            !
  USE MOD_Particle_Surfaces_Vars, ONLY: nPartCurved, SuperSampledNodes,nTriangles
  USE MOD_Mesh_Vars,     ONLY : ElemToSide
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
   IMPLICIT NONE                                                                                   !
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION                                                                        !
   INTEGER                          :: i                                                           !
   INTEGER                          :: iLocSide                                                    !
   INTEGER                          :: Element                                                     !
   INTEGER                          :: TriNum                                                      !
   LOGICAL                          :: ThroughSide                                                 !
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION                                                                       !
   INTEGER                          :: n, j, m                                                     !
   REAL                             :: Px, Py, Pz                                                  !
   REAL                             :: Vx, Vy, Vz, Vall                                            !
   REAL                             :: xNode(3), yNode(3), zNode(3), Ax(3), Ay(3), Az(3)           !
   REAL                             :: det(3)                                                      !
   REAL                             :: eps                                                         !
   INTEGER                          :: QuadID,p,q,SideID,flip,h
!----------------------------------------------------------------------------------------------------------------------------------
   INTENT(IN)                       :: i                                                           !
   INTENT(OUT)                      :: ThroughSide
!===================================================================================================================================


   eps = 0.
   SideID=ElemToSide(E2S_SIDE_ID,ilocSide,Element)
   flip  =ElemToSide(E2S_FLIP,ilocSide,Element)

   ThroughSide = .FALSE.

   Px = lastPartPos(i,1)
   Py = lastPartPos(i,2)
   Pz = lastPartPos(i,3)

   Vx = PartState(i,1)-lastPartPos(i,1)
   Vy = PartState(i,2)-lastPartPos(i,2)
   Vz = PartState(i,3)-lastPartPos(i,3)

   Vall = SQRT(Vx*Vx + Vy*Vy + Vz*Vz)

   Vx = Vx/Vall
   Vy = Vy/Vall
   Vz = Vz/Vall

!   xNode(1) = GEO%NodeCoords(1,GEO%ElemSideNodeID(1,iLocSide,Element))
!   yNode(1) = GEO%NodeCoords(2,GEO%ElemSideNodeID(1,iLocSide,Element))
!   zNode(1) = GEO%NodeCoords(3,GEO%ElemSideNodeID(1,iLocSide,Element))
!   Ax(1) = xNode(1) - Px
!   Ay(1) = yNode(1) - Py
!   Az(1) = zNode(1) - Pz

   ! we have to compute the p,q coordinates

   ! quadrilateral id
   QuadID=triNum/2+MOD(triNum,2)
   q=(QuadID-1)/NPartCurved ! fortran takes floor of integer devision
   p=MOD(QuadID-1,NPartCurved)
   
   !print*,'p,q',p,q
   !xNode=0.
   !yNode=0.
   !zNode=0.
   xNode(1) = SuperSampledNodes(1,p  ,q  ,SideID)
   yNode(1) = SuperSampledNodes(2,p  ,q  ,SideID)
   zNode(1) = SuperSampledNodes(3,p  ,q  ,SideID)
   Ax(1) = xNode(1) - Px
   Ay(1) = yNode(1) - Py
   Az(1) = zNode(1) - Pz

   !print*,'Node Coord in x,y,z'
   !print*,'x1',SuperSampledNodes(:,0  ,0  ,SideID)
   !print*,'x2',SuperSampledNodes(:,1  ,0  ,SideID)
   !print*,'x3',SuperSampledNodes(:,1  ,1  ,SideID)
   !print*,'x4',SuperSampledNodes(:,0  ,1  ,SideID)

   IF(MOD(triNum,2).EQ.0)THEN
     !print*,'tri2,p2',p+1,q+1
     xNode(2) = SuperSampledNodes(1,p+1,q+1,SideID)
     yNode(2) = SuperSampledNodes(2,p+1,q+1,SideID)
     zNode(2) = SuperSampledNodes(3,p+1,q+1,SideID)
     Ax(2) = xNode(2) - Px
     Ay(2) = yNode(2) - Py
     Az(2) = zNode(2) - Pz

     m=0!MOD(triNum,2)
     n=1!MOD(triNum+1,2)
     !print*,'tri2,p3',p+m,q+n
     xNode(3) = SuperSampledNodes(1,p+m,q+n,SideID)
     yNode(3) = SuperSampledNodes(2,p+m,q+n,SideID)
     zNode(3) = SuperSampledNodes(3,p+m,q+n,SideID)
     Ax(3) = xNode(3) - Px
     Ay(3) = yNode(3) - Py
       Az(3) = zNode(3) - Pz
   ELSE
     !print*,'not here'
     m=1!MOD(triNum,2)
     n=0!MOD(triNum+1,2)
     xNode(2) = SuperSampledNodes(1,p+m,q+n,SideID)
     yNode(2) = SuperSampledNodes(2,p+m,q+n,SideID)
     zNode(2) = SuperSampledNodes(3,p+m,q+n,SideID)
     Ax(2) = xNode(2) - Px
     Ay(2) = yNode(2) - Py
     Az(2) = zNode(2) - Pz

     xNode(3) = SuperSampledNodes(1,p+1,q+1,SideID)
     yNode(3) = SuperSampledNodes(2,p+1,q+1,SideID)
     zNode(3) = SuperSampledNodes(3,p+1,q+1,SideID)
     Ax(3) = xNode(3) - Px
     Ay(3) = yNode(3) - Py
     Az(3) = zNode(3) - Pz

   END IF


   !print*,'xNode',xNode(1),xNode(2),xNode(3)
   !print*,'yNode',yNode(1),yNode(2),yNode(3)
   !print*,'zNode',zNode(1),zNode(2),zNode(3)
   !print*,'dad',p+m,q+n
   !xNode(2) = SuperSampledNodes(1,p+m,q+n,SideID)
   !yNode(2) = SuperSampledNodes(2,p+m,q+n,SideID)
   !zNode(2) = SuperSampledNodes(3,p+m,q+n,SideID)
   !Ax(2) = xNode(2) - Px
   !Ay(2) = yNode(2) - Py
   !Az(2) = zNode(2) - Pz
   !print*,'a1',ax(1),ay(1),az(1)
   !print*,'a2',ax(2),ay(2),az(2)
   !print*,'a3',ax(3),ay(3),az(3)

   !xNode(3) = SuperSampledNodes(1,p+1,q+1,SideID)
   !yNode(3) = SuperSampledNodes(2,p+1,q+1,SideID)
   !zNode(3) = SuperSampledNodes(3,p+1,q+1,SideID)
   !Ax(3) = xNode(3) - Px
   !Ay(3) = yNode(3) - Py
   !Az(3) = zNode(3) - Pz

!!$IF((i.EQ.3389).AND.(iLocSide.EQ.1).AND.(Element.EQ.10913))THEN
!!$WRITE(*,*) 'Triangle',j            
!!$WRITE(*,*) xNode(1),yNode(1),zNode(1)
!!$WRITE(*,*) xNode(2),yNode(2),zNode(2)
!!$WRITE(*,*) xNode(3),yNode(3),zNode(3)
!!$!WRITE(*,*) Vx,Vy,Vz
!!$!WRITE(*,*) Ax(1),Ay(1),Az(1)
!!$!WRITE(*,*) Ax(2),Ay(2),Az(2)
!!$!WRITE(*,*) Ax(3),Ay(3),Az(3)
!!$END IF         
   !--- check whether v and the vectors from the particle to the two edge nodes build
   !--- a right-hand-system. If yes for all edges: vector goes potentially through side
   det(1) = ((Ay(1) * Vz - Az(1) * Vy) * Ax(3)  + &
             (Az(1) * Vx - Ax(1) * Vz) * Ay(3)  + &
             (Ax(1) * Vy - Ay(1) * Vx) * Az(3))
   
   det(2) = ((Ay(2) * Vz - Az(2) * Vy) * Ax(1)  + &
             (Az(2) * Vx - Ax(2) * Vz) * Ay(1)  + &
             (Ax(2) * Vy - Ay(2) * Vx) * Az(1))
   
   det(3) = ((Ay(3) * Vz - Az(3) * Vy) * Ax(2)  + &
             (Az(3) * Vx - Ax(3) * Vz) * Ay(2)  + &
             (Ax(3) * Vy - Ay(3) * Vx) * Az(2))
!!$IF((i.EQ.3389).AND.(iLocSide.EQ.1).AND.(Element.EQ.10913))THEN
!!$WRITE(*,*) 'det:'
!!$WRITE(*,*) det(:)
!!$END IF
!IF((i.EQ.359616).AND.(Element.EQ.13718))THEN
!WRITE(*,*) '======='
!WRITE(*,*) iLocSide,TriNum
!WRITE(*,*) det(:)
!WRITE(*,*) XNODE(1),YNODE(1),ZNODE(1)
!WRITE(*,*) XNODE(2),YNODE(2),ZNODE(2)
!WRITE(*,*) XNODE(3),YNODE(3),ZNODE(3)
!WRITE(*,*) VX,VY,VZ
!WRITE(*,*) AX(1),AY(1),AZ(1)
!WRITE(*,*) AX(2),AY(2),AZ(2)
!WRITE(*,*) AX(3),AY(3),AZ(3)
!WRITE(*,*) '=-=-=-='
!END IF
   IF(flip.NE.0) det=-det
   IF ((det(1).ge.-eps).AND.(det(2).ge.-eps).AND.(det(3).ge.-eps)) THEN
   !  print*,'here',ilocSide
     ThroughSide = .TRUE.
   END IF
 RETURN
END SUBROUTINE ParticleThroughSideCheck3DFast

SUBROUTINE ParticleThroughSideCheck3DFastold(i,iLocSide,Element,ThroughSide,TriNum)                   !
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars                                                                            !
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
   IMPLICIT NONE                                                                                   !
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION                                                                        !
   INTEGER                          :: i                                                           !
   INTEGER                          :: iLocSide                                                    !
   INTEGER                          :: Element                                                     !
   INTEGER                          :: TriNum                                                      !
   LOGICAL                          :: ThroughSide                                                 !
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION                                                                       !
   INTEGER                          :: n, j, m                                                     !
   REAL                             :: Px, Py, Pz                                                  !
   REAL                             :: Vx, Vy, Vz, Vall                                            !
   REAL                             :: xNode(3), yNode(3), zNode(3), Ax(3), Ay(3), Az(3)           !
   REAL                             :: det(3)                                                      !
   REAL                             :: eps                                                         !
!----------------------------------------------------------------------------------------------------------------------------------
   INTENT(IN)                       :: i                                                           !
   INTENT(OUT)                      :: ThroughSide
!===================================================================================================================================


   eps = 0.

   ThroughSide = .FALSE.

   Px = lastPartPos(i,1)
   Py = lastPartPos(i,2)
   Pz = lastPartPos(i,3)

   Vx = PartState(i,1)-lastPartPos(i,1)
   Vy = PartState(i,2)-lastPartPos(i,2)
   Vz = PartState(i,3)-lastPartPos(i,3)

   Vall = SQRT(Vx*Vx + Vy*Vy + Vz*Vz)

   Vx = Vx/Vall
   Vy = Vy/Vall
   Vz = Vz/Vall

   xNode(1) = GEO%NodeCoords(1,GEO%ElemSideNodeID(1,iLocSide,Element))
   yNode(1) = GEO%NodeCoords(2,GEO%ElemSideNodeID(1,iLocSide,Element))
   zNode(1) = GEO%NodeCoords(3,GEO%ElemSideNodeID(1,iLocSide,Element))
   Ax(1) = xNode(1) - Px
   Ay(1) = yNode(1) - Py
   Az(1) = zNode(1) - Pz

   DO n = 2,3
     m = n+TriNum-1       ! m = true node number of the sides
     xNode(n) = GEO%NodeCoords(1,GEO%ElemSideNodeID(m,iLocSide,Element))
     yNode(n) = GEO%NodeCoords(2,GEO%ElemSideNodeID(m,iLocSide,Element))
     zNode(n) = GEO%NodeCoords(3,GEO%ElemSideNodeID(m,iLocSide,Element))

     Ax(n) = xNode(n) - Px
     Ay(n) = yNode(n) - Py
     Az(n) = zNode(n) - Pz
   END DO
   !print*,'xNode',xNode(1),xNode(2),xNode(3)
   !print*,'yNode',yNode(1),yNode(2),yNode(3)
   !print*,'zNode',zNode(1),zNode(2),zNode(3)

   !print*,'a1',ax(1),ay(1),az(1)
   !print*,'a2',ax(2),ay(2),az(2)
   !print*,'a3',ax(3),ay(3),az(3)
!!$IF((i.EQ.3389).AND.(iLocSide.EQ.1).AND.(Element.EQ.10913))THEN
!!$WRITE(*,*) 'Triangle',j            
!!$WRITE(*,*) xNode(1),yNode(1),zNode(1)
!!$WRITE(*,*) xNode(2),yNode(2),zNode(2)
!!$WRITE(*,*) xNode(3),yNode(3),zNode(3)
!!$!WRITE(*,*) Vx,Vy,Vz
!!$!WRITE(*,*) Ax(1),Ay(1),Az(1)
!!$!WRITE(*,*) Ax(2),Ay(2),Az(2)
!!$!WRITE(*,*) Ax(3),Ay(3),Az(3)
!!$END IF         
   !--- check whether v and the vectors from the particle to the two edge nodes build
   !--- a right-hand-system. If yes for all edges: vector goes potentially through side
   det(1) = ((Ay(1) * Vz - Az(1) * Vy) * Ax(3)  + &
             (Az(1) * Vx - Ax(1) * Vz) * Ay(3)  + &
             (Ax(1) * Vy - Ay(1) * Vx) * Az(3))
   
   det(2) = ((Ay(2) * Vz - Az(2) * Vy) * Ax(1)  + &
             (Az(2) * Vx - Ax(2) * Vz) * Ay(1)  + &
             (Ax(2) * Vy - Ay(2) * Vx) * Az(1))
   
   det(3) = ((Ay(3) * Vz - Az(3) * Vy) * Ax(2)  + &
             (Az(3) * Vx - Ax(3) * Vz) * Ay(2)  + &
             (Ax(3) * Vy - Ay(3) * Vx) * Az(2))
!!$IF((i.EQ.3389).AND.(iLocSide.EQ.1).AND.(Element.EQ.10913))THEN
!!$WRITE(*,*) 'det:'
!!$WRITE(*,*) det(:)
!!$END IF
!IF((i.EQ.359616).AND.(Element.EQ.13718))THEN
!WRITE(*,*) '======='
!WRITE(*,*) iLocSide,TriNum
!WRITE(*,*) det(:)
!WRITE(*,*) XNODE(1),YNODE(1),ZNODE(1)
!WRITE(*,*) XNODE(2),YNODE(2),ZNODE(2)
!WRITE(*,*) XNODE(3),YNODE(3),ZNODE(3)
!WRITE(*,*) VX,VY,VZ
!WRITE(*,*) AX(1),AY(1),AZ(1)
!WRITE(*,*) AX(2),AY(2),AZ(2)
!WRITE(*,*) AX(3),AY(3),AZ(3)
!WRITE(*,*) '=-=-=-='
!END IF
   IF ((det(1).ge.-eps).AND.(det(2).ge.-eps).AND.(det(3).ge.-eps)) THEN
   !  print*,'here',ilocSide
     ThroughSide = .TRUE.
   END IF
 RETURN
END SUBROUTINE ParticleThroughSideCheck3DFastold

!!$SUBROUTINE ParticleThroughSideScalarCheck(i,iLocSide,Element,ThroughSide,TriNum,det)            !
!!$  USE MOD_Particle_Vars                                                                            !
!!$!--------------------------------------------------------------------------------------------------!
!!$   IMPLICIT NONE                                                                                   !
!!$!--------------------------------------------------------------------------------------------------!
!!$! argument list declaration
!!$   INTEGER                          :: i
!!$   INTEGER                          :: iLocSide
!!$   INTEGER                          :: Element
!!$   INTEGER                          :: TriNum
!!$   LOGICAL                          :: ThroughSide
!!$   REAL                             :: det(1:3)
!!$! Local variable declaration
!!$   INTEGER                          :: n, j, m, counter,negcount
!!$   REAL                             :: Px, Py, Pz
!!$   REAL                             :: Vx, Vy, Vz, Vall
!!$   REAL                             :: Ax, Ay, Az, tempdet(1:3)
!!$!--------------------------------------------------------------------------------------------------!
!!$   INTENT(IN)                       :: i,TriNum
!!$   INTENT(OUT)                      :: ThroughSide
!!$!--------------------------------------------------------------------------------------------------!
!!$
!!$ !--- if the particle ends up in this subroutine, it means that it apparently has crossed more than one
!!$ !    sides. This can happen if the hexagon has concave sides
!!$!IF(i.EQ.99981)THEN
!!$!WRITE(*,*) '---',iLocSide, TriNum
!!$!WRITE(*,*)  GEO%NodeCoords(:,GEO%ElemSideNodeID(1,iLocSide,Element))
!!$!WRITE(*,*)  GEO%NodeCoords(:,GEO%ElemSideNodeID(2,iLocSide,Element))
!!$!WRITE(*,*)  GEO%NodeCoords(:,GEO%ElemSideNodeID(3,iLocSide,Element))
!!$!WRITE(*,*)  GEO%NodeCoords(:,GEO%ElemSideNodeID(4,iLocSide,Element))
!!$!WRITE(*,*) PartState(i,1:3)
!!$!WRITE(*,*) lastPartPos(i,1:3)
!!$!END IF
!!$
!!$ ! Case 1: Particle has already been on outside of the side at LastPos -> not really through side
!!$ ! Check by using the scalar product of the velocity vector and the three vectors to the nodes
!!$
!!$ ThroughSide = .FALSE.
!!$
!!$ Px = lastPartPos(i,1)
!!$ Py = lastPartPos(i,2)
!!$ Pz = lastPartPos(i,3)
!!$
!!$ Vx = PartState(i,1)-lastPartPos(i,1)
!!$ Vy = PartState(i,2)-lastPartPos(i,2)
!!$ Vz = PartState(i,3)-lastPartPos(i,3)
!!$
!!$ Vall = SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
!!$
!!$ Vx = Vx/Vall
!!$ Vy = Vy/Vall
!!$ Vz = Vz/Vall
!!$
!!$ ! Check Node 1:
!!$ Ax = GEO%NodeCoords(1,GEO%ElemSideNodeID(1,iLocSide,Element)) - Px
!!$ Ay = GEO%NodeCoords(2,GEO%ElemSideNodeID(1,iLocSide,Element)) - Py
!!$ Az = GEO%NodeCoords(3,GEO%ElemSideNodeID(1,iLocSide,Element)) - Pz
!!$ det(1) = Ax*Vx + Ay*Vy + Az*Vz
!!$ IF (det(1).GE.0) ThroughSide = .TRUE.
!!$
!!$ ! Check Nodes 2+3 or 3+4 (depending on Triangle TriNum):
!!$ DO n = 2,3
!!$   m = n+TriNum-1       ! m = true node number of the sides
!!$   Ax = GEO%NodeCoords(1,GEO%ElemSideNodeID(m,iLocSide,Element)) - Px
!!$   Ay = GEO%NodeCoords(2,GEO%ElemSideNodeID(m,iLocSide,Element)) - Py
!!$   Az = GEO%NodeCoords(3,GEO%ElemSideNodeID(m,iLocSide,Element)) - Pz
!!$   det(n) = Ax*Vx + Ay*Vy + Az*Vz
!!$   IF (det(n).GE.0) ThroughSide = .TRUE.
!!$ END DO
!!$ ! Sort det by size, positives first
!!$ negcount = 0
!!$ DO n = 1,3
!!$   IF (det(n).LT.0) negcount = negcount + 1
!!$ END DO
!!$ counter = 1
!!$ DO WHILE (counter.LT.4)
!!$   DO n = 1,3
!!$     IF(MIN(det(1),det(2),det(3)).EQ.det(n))THEN
!!$       tempdet(counter) = det(n)
!!$       counter = counter + 1
!!$       det(n) = HUGE(det(n))
!!$     END IF
!!$   END DO
!!$ END DO
!!$ DO n = negcount+1 , 3   ! positives first
!!$   det(n-negcount) = tempdet(n)
!!$ END DO
!!$ counter = 1
!!$ DO n = negcount,1,-1       ! then negatives
!!$   det(3-negcount+counter) = tempdet(n)
!!$   counter = counter + 1
!!$ END DO
!!$!IF(i.EQ.99981)THEN
!!$!WRITE(*,*) 'det',det(:)
!!$!END IF
!!$ RETURN
!!$END SUBROUTINE ParticleThroughSideScalarCheck

#ifdef MPI
SUBROUTINE ParticleThroughSideCheck3DFast_halocells(i,iLocSide,Element,ThroughSide,TriNum)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars, ONLY : PartState, lastPartPos
USE MOD_part_MPI_Vars, ONLY : MPIGEO 
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
   IMPLICIT NONE                                                                                   !
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION                                                                        !
   INTEGER                          :: i                                                           !
   INTEGER                          :: iLocSide                                                    !
   INTEGER                          :: Element                                                     !
   INTEGER                          :: TriNum                                                      !
   LOGICAL                          :: ThroughSide                                                 !
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION                                                                       !
   INTEGER                          :: n, j, m                                                     !
   REAL                             :: Px, Py, Pz                                                  !
   REAL                             :: Vx, Vy, Vz, Vall                                            !
   REAL                             :: xNode(3), yNode(3), zNode(3), Ax(3), Ay(3), Az(3)           !
   REAL                             :: det(3)                                                      !
   REAL                             :: eps                                                         !
!--------------------------------------------------------------------------------------------------!
   INTENT(IN)                       :: i                                                           !
   INTENT(OUT)                      :: ThroughSide                                         !
!===================================================================================================================================


   eps = 0.

   ThroughSide = .FALSE.

   Px = lastPartPos(i,1)
   Py = lastPartPos(i,2)
   Pz = lastPartPos(i,3)

   Vx = PartState(i,1)-lastPartPos(i,1)
   Vy = PartState(i,2)-lastPartPos(i,2)
   Vz = PartState(i,3)-lastPartPos(i,3)

   Vall = SQRT(Vx*Vx + Vy*Vy + Vz*Vz)

   Vx = Vx/Vall
   Vy = Vy/Vall
   Vz = Vz/Vall

   xNode(1) = MPIGEO%NodeCoords(1,MPIGEO%ElemSideNodeID(1,iLocSide,Element))
   yNode(1) = MPIGEO%NodeCoords(2,MPIGEO%ElemSideNodeID(1,iLocSide,Element))
   zNode(1) = MPIGEO%NodeCoords(3,MPIGEO%ElemSideNodeID(1,iLocSide,Element))
   Ax(1) = xNode(1) - Px
   Ay(1) = yNode(1) - Py
   Az(1) = zNode(1) - Pz

   DO n = 2,3
     m = n+TriNum-1       ! m = true node number of the sides
     xNode(n) = MPIGEO%NodeCoords(1,MPIGEO%ElemSideNodeID(m,iLocSide,Element))
     yNode(n) = MPIGEO%NodeCoords(2,MPIGEO%ElemSideNodeID(m,iLocSide,Element))
     zNode(n) = MPIGEO%NodeCoords(3,MPIGEO%ElemSideNodeID(m,iLocSide,Element))
     
     Ax(n) = xNode(n) - Px
     Ay(n) = yNode(n) - Py
     Az(n) = zNode(n) - Pz
   END DO
   
   !--- check whether v and the vectors from the particle to the two edge nodes build
   !--- a right-hand-system. If yes for all edges: vector goes through side
   det(1) = ((Ay(1) * Vz - Az(1) * Vy) * Ax(3)  + &
             (Az(1) * Vx - Ax(1) * Vz) * Ay(3)  + &
             (Ax(1) * Vy - Ay(1) * Vx) * Az(3))
         
   det(2) = ((Ay(2) * Vz - Az(2) * Vy) * Ax(1)  + &
             (Az(2) * Vx - Ax(2) * Vz) * Ay(1)  + &
             (Ax(2) * Vy - Ay(2) * Vx) * Az(1))
         
   det(3) = ((Ay(3) * Vz - Az(3) * Vy) * Ax(2)  + &
             (Az(3) * Vx - Ax(3) * Vz) * Ay(2)  + &
             (Ax(3) * Vy - Ay(3) * Vx) * Az(2))
   IF ((det(1).ge.-eps).AND.(det(2).ge.-eps).AND.(det(3).ge.-eps)) THEN
     ThroughSide = .TRUE.
   END IF
   RETURN
 END SUBROUTINE ParticleThroughSideCheck3DFast_halocells

!!$SUBROUTINE ParticleThroughSideScalarCheck_halocells(i,iLocSide,Element,ThroughSide,TriNum,det)
!!$  USE MOD_Particle_Vars , ONLY : PartState, lastPartPos
!!$  USE MOD_part_MPI_Vars,  ONLY : MPIGEO                                                                            !
!!$!--------------------------------------------------------------------------------------------------!
!!$   IMPLICIT NONE                                                                                   !
!!$!--------------------------------------------------------------------------------------------------!
!!$! argument list declaration
!!$   INTEGER                          :: i
!!$   INTEGER                          :: iLocSide
!!$   INTEGER                          :: Element
!!$   INTEGER                          :: TriNum
!!$   LOGICAL                          :: ThroughSide
!!$   REAL                             :: det(1:3)
!!$! Local variable declaration
!!$   INTEGER                          :: n, j, m, negcount,counter
!!$   REAL                             :: Px, Py, Pz, tempdet(3)
!!$   REAL                             :: Vx, Vy, Vz, Vall
!!$   REAL                             :: Ax, Ay, Az
!!$!--------------------------------------------------------------------------------------------------!
!!$   INTENT(IN)                       :: i,TriNum                                                    !
!!$   INTENT(OUT)                      :: ThroughSide
!!$!--------------------------------------------------------------------------------------------------!
!!$
!!$ !--- if the particle ends up in this subroutine, it means that it apparently has crossed more than one
!!$ !    sides. This can happen if the hexagon has concave sides
!!$
!!$ ! Case 1: Particle has already been on outside of the side at LastPos -> not really through side
!!$ ! Check by using the scalar product of the velocity vector and the three vectors to the nodes
!!$
!!$ Px = lastPartPos(i,1)
!!$ Py = lastPartPos(i,2)
!!$ Pz = lastPartPos(i,3)
!!$
!!$ Vx = PartState(i,1)-lastPartPos(i,1)
!!$ Vy = PartState(i,2)-lastPartPos(i,2)
!!$ Vz = PartState(i,3)-lastPartPos(i,3)
!!$
!!$ Vall = SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
!!$
!!$ Vx = Vx/Vall
!!$ Vy = Vy/Vall
!!$ Vz = Vz/Vall
!!$
!!$ ThroughSide = .FALSE.
!!$ ! Check Node 1:
!!$ Ax = MPIGEO%NodeCoords(1,MPIGEO%ElemSideNodeID(1,iLocSide,Element)) - Px
!!$ Ay = MPIGEO%NodeCoords(2,MPIGEO%ElemSideNodeID(1,iLocSide,Element)) - Py
!!$ Az = MPIGEO%NodeCoords(3,MPIGEO%ElemSideNodeID(1,iLocSide,Element)) - Pz
!!$ det(1) = Ax*Vx + Ay*Vy + Az*Vz
!!$ IF (det(1).GE.0) ThroughSide = .TRUE.
!!$
!!$ ! Check Nodes 2+3 or 3+4 (depending on Triangle TriNum):
!!$ DO n = 2,3
!!$   m = n+TriNum-1       ! m = true node number of the sides
!!$   Ax = MPIGEO%NodeCoords(1,MPIGEO%ElemSideNodeID(m,iLocSide,Element)) - Px
!!$   Ay = MPIGEO%NodeCoords(2,MPIGEO%ElemSideNodeID(m,iLocSide,Element)) - Py
!!$   Az = MPIGEO%NodeCoords(3,MPIGEO%ElemSideNodeID(m,iLocSide,Element)) - Pz
!!$   det(n) = Ax*Vx + Ay*Vy + Az*Vz
!!$   IF (det(n).GE.0) ThroughSide = .TRUE.
!!$ END DO
!!$ ! Sort det by size, put to huge if negative
!!$ ! Sort det by size, positives first
!!$ negcount = 0
!!$ DO n = 1,3
!!$   IF (det(n).LT.0) negcount = negcount + 1
!!$ END DO
!!$ counter = 1
!!$ DO WHILE (counter.LT.4)
!!$   DO n = 1,3
!!$     IF(MIN(det(1),det(2),det(3)).EQ.det(n))THEN
!!$       tempdet(counter) = det(n)
!!$       counter = counter + 1
!!$       det(n) = HUGE(det(n))
!!$     END IF
!!$   END DO
!!$ END DO
!!$ DO n = negcount+1 , 3   ! positives first
!!$   det(n-negcount) = tempdet(n)
!!$ END DO
!!$ counter = 1
!!$ DO n = negcount,1,-1       ! then negatives
!!$   det(3-negcount+counter) = tempdet(n)
!!$   counter = counter + 1
!!$ END DO
!!$ RETURN
!!$END SUBROUTINE ParticleThroughSideScalarCheck_halocells


#endif /*MPI*/
                                                                                                   !
SUBROUTINE ParticleThroughSideCheck3DFast___OLD(i,iLocSide,Element,ThroughSide)                          !
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars                                                                            !
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
   IMPLICIT NONE                                                                                   !
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION                                                                        !
   INTEGER                          :: i                                                           !
   INTEGER                          :: iLocSide                                                    !
   INTEGER                          :: Element                                                     !
   LOGICAL                          :: ThroughSide                                                 !
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION                                                                       !
   INTEGER                          :: n, j, m                                                     !
   REAL                             :: Px, Py, Pz                                                  !
   REAL                             :: Vx, Vy, Vz, Vall                                            !
   REAL                             :: xNode(3), yNode(3), zNode(3), Ax(3), Ay(3), Az(3)           !
   REAL                             :: det(3)                                                      !
   REAL                             :: eps, Flux                                                   !
!----------------------------------------------------------------------------------------------------------------------------------
   INTENT(IN)                       :: i                                                           !
   INTENT(OUT)                      :: ThroughSide                                                 !
!===================================================================================================================================


   ThroughSide = .FALSE.

   eps = 0

   Px = lastPartPos(i,1)
   Py = lastPartPos(i,2)
   Pz = lastPartPos(i,3)

   Vx = PartState(i,1)-lastPartPos(i,1)
   Vy = PartState(i,2)-lastPartPos(i,2)
   Vz = PartState(i,3)-lastPartPos(i,3)

   Vall = SQRT(Vx*Vx + Vy*Vy + Vz*Vz)

   Vx = Vx/Vall
   Vy = Vy/Vall
   Vz = Vz/Vall

   xNode(1) = GEO%NodeCoords(1,GEO%ElemSideNodeID(1,iLocSide,Element))
   yNode(1) = GEO%NodeCoords(2,GEO%ElemSideNodeID(1,iLocSide,Element))
   zNode(1) = GEO%NodeCoords(3,GEO%ElemSideNodeID(1,iLocSide,Element))
   Ax(1) = xNode(1) - Px
   Ay(1) = yNode(1) - Py
   Az(1) = zNode(1) - Pz

   DO j = 1, 2

      DO n = 2,3
         m = n+j-1
         xNode(n) = GEO%NodeCoords(1,GEO%ElemSideNodeID(m,iLocSide,Element))
         yNode(n) = GEO%NodeCoords(2,GEO%ElemSideNodeID(m,iLocSide,Element))
         zNode(n) = GEO%NodeCoords(3,GEO%ElemSideNodeID(m,iLocSide,Element))

         Ax(n) = xNode(n) - Px
         Ay(n) = yNode(n) - Py
         Az(n) = zNode(n) - Pz
      END DO

      Flux = 1 ! REAL (MESH%ELEM%FluxSign(GlobElemNbr,LocSideNbr))
      det(1) = ((Ay(1) * Vz - Az(1) * Vy) * Ax(3)  + &
                (Az(1) * Vx - Ax(1) * Vz) * Ay(3)  + &
                (Ax(1) * Vy - Ay(1) * Vx) * Az(3)) * &
                Flux

      det(2) = ((Ay(2) * Vz - Az(2) * Vy) * Ax(1)  + &
                (Az(2) * Vx - Ax(2) * Vz) * Ay(1)  + &
                (Ax(2) * Vy - Ay(2) * Vx) * Az(1)) * &
                Flux

      det(3) = ((Ay(3) * Vz - Az(3) * Vy) * Ax(2)  + &
                (Az(3) * Vx - Ax(3) * Vz) * Ay(2)  + &
                (Ax(3) * Vy - Ay(3) * Vx) * Az(2)) * &
                Flux
      IF ((det(1).ge.-eps).AND.(det(2).ge.-eps).AND.(det(3).ge.-eps)) THEN
         ThroughSide = .TRUE.
      END IF
   END DO
   RETURN
  END SUBROUTINE ParticleThroughSideCheck3DFast___OLD

FUNCTION ParticleThroughSideCheck3D(i,iLocSide,Element) RESULT(ThroughSide)                        !
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
   IMPLICIT NONE                                                                                   !
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION                                                                        !
   INTEGER                          :: i                                                           !
   INTEGER                          :: iLocSide                                                    !
   INTEGER                          :: Element                                                     !
   LOGICAL                          :: ThroughSide                                                 !
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION                                                                       !
   INTEGER                          :: n, j, m                                                     !
   REAL                             :: Px, Py, Pz                                                  !
   REAL                             :: Vx, Vy, Vz, Vall                                            !
   REAL                             :: xNode(3), yNode(3), zNode(3), Ax(3), Ay(3), Az(3)           !
   REAL                             :: Point(3,3)                                                  !
   REAL                             :: det(5)                                                      !
   REAL                             :: eps, Flux                                                   !
!--------------------------------------------------------------------------------------------------!
   INTENT(IN)                       :: i                                                           !
!===================================================================================================================================


   ThroughSide = .FALSE.

   eps = 0

   Px = lastPartPos(i,1)
   Py = lastPartPos(i,2)
   Pz = lastPartPos(i,3)

   Vx = PartState(i,1)-lastPartPos(i,1)
   Vy = PartState(i,2)-lastPartPos(i,2)
   Vz = PartState(i,3)-lastPartPos(i,3)

   Vall = SQRT(Vx*Vx + Vy*Vy + Vz*Vz)

   Vx = Vx/Vall
   Vy = Vy/Vall
   Vz = Vz/Vall

   xNode(1) = GEO%NodeCoords(1,GEO%ElemSideNodeID(1,iLocSide,Element))
   yNode(1) = GEO%NodeCoords(2,GEO%ElemSideNodeID(1,iLocSide,Element))
   zNode(1) = GEO%NodeCoords(3,GEO%ElemSideNodeID(1,iLocSide,Element))
   Ax(1) = xNode(1) - Px
   Ay(1) = yNode(1) - Py
   Az(1) = zNode(1) - Pz

   DO j = 1, 2

      DO n = 2,3
         m = n+j-1
         xNode(n) = GEO%NodeCoords(1,GEO%ElemSideNodeID(m,iLocSide,Element))
         yNode(n) = GEO%NodeCoords(2,GEO%ElemSideNodeID(m,iLocSide,Element))
         zNode(n) = GEO%NodeCoords(3,GEO%ElemSideNodeID(m,iLocSide,Element))

         Ax(n) = xNode(n) - Px
         Ay(n) = yNode(n) - Py
         Az(n) = zNode(n) - Pz
      END DO

      Flux = 1 ! REAL (MESH%ELEM%FluxSign(GlobElemNbr,LocSideNbr))
      det(1) = ((Ay(1) * Vz - Az(1) * Vy) * Ax(3)  + &
                (Az(1) * Vx - Ax(1) * Vz) * Ay(3)  + &
                (Ax(1) * Vy - Ay(1) * Vx) * Az(3)) * &
                Flux

      det(2) = ((Ay(2) * Vz - Az(2) * Vy) * Ax(1)  + &
                (Az(2) * Vx - Ax(2) * Vz) * Ay(1)  + &
                (Ax(2) * Vy - Ay(2) * Vx) * Az(1)) * &
                Flux

      det(3) = ((Ay(3) * Vz - Az(3) * Vy) * Ax(2)  + &
                (Az(3) * Vx - Ax(3) * Vz) * Ay(2)  + &
                (Ax(3) * Vy - Ay(3) * Vx) * Az(2)) * &
                Flux
      IF ((det(1).ge.-eps).AND.(det(2).ge.-eps).AND.(det(3).ge.-eps)) THEN
         ThroughSide = .TRUE.
      END IF
   END DO
   !--- Now double-check that the particle has actually crossed the side within this timestep
   IF (ThroughSide) THEN
     Point(1,:)=0.25*(GEO%NodeCoords(:,GEO%ElemSideNodeID(1,iLocSide,Element))+ & ! Barycenter
                      GEO%NodeCoords(:,GEO%ElemSideNodeID(2,iLocSide,Element))+ &
                      GEO%NodeCoords(:,GEO%ElemSideNodeID(3,iLocSide,Element))+ &
                      GEO%NodeCoords(:,GEO%ElemSideNodeID(4,iLocSide,Element)))
     Point(2,:)=Point(1,:) + GEO%NodeCoords(:,GEO%ElemSideNodeID(3,iLocSide,Element)) &
                           - GEO%NodeCoords(:,GEO%ElemSideNodeID(1,iLocSide,Element))
     Point(3,:)=Point(1,:) + GEO%NodeCoords(:,GEO%ElemSideNodeID(4,iLocSide,Element)) &
                           - GEO%NodeCoords(:,GEO%ElemSideNodeID(2,iLocSide,Element))
     Ax(:) = Point(:,1) - Px
     Ay(:) = Point(:,2) - Py
     Az(:) = Point(:,3) - Pz
     det(4) = ((Ay(1) * Az(2) - Az(1) * Ay(2)) * Ax(3) +     &
               (Az(1) * Ax(2) - Ax(1) * Az(2)) * Ay(3) +     &
               (Ax(1) * Ay(2) - Ay(1) * Ax(2)) * Az(3))

     Ax(:) = Point(:,1) - PartState(i,1)
     Ay(:) = Point(:,2) - PartState(i,2)
     Az(:) = Point(:,3) - PartState(i,3)
     det(5) = ((Ay(1) * Az(2) - Az(1) * Ay(2)) * Ax(3) +     &
               (Az(1) * Ax(2) - Ax(1) * Az(2)) * Ay(3) +     &
               (Ax(1) * Ay(2) - Ay(1) * Ax(2)) * Az(3))

     IF ((det(4).lt.-eps).OR.(det(4).NE.det(4)).OR.(det(5).ge.-eps).OR.(det(5).NE.det(5))) THEN
        ThroughSide = .FALSE.
     END IF
   END IF
   RETURN
  END FUNCTION ParticleThroughSideCheck3D

SUBROUTINE ParticleInElementCheck3D(i,Element,InElementCheck,det)                                  !
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
!DEC$ ATTRIBUTES FORCEINLINE :: ParticleInElementCheck3D
   USE MOD_Particle_Vars
   USE MOD_Mesh_Vars,     ONLY : ElemToSide
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
   IMPLICIT NONE                                                                                   !
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION                                                                        !
   INTEGER                          :: i, Element                                                  !
   LOGICAL                          :: InElementCheck                                              !
   REAL                             :: det(6)                                                      !
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION                                                                       !
   INTEGER                          :: iLocSide,n                                                  !
   REAL                             :: epsilon, xNode(3), yNode(3), zNode(3)                       !
   REAL                             :: Ax(3),Ay(3),Az(3)                                           !
!--------------------------------------------------------------------------------------------------!
   INTENT(IN)                       :: Element                                                     !
   INTENT(OUT)                      :: InElementCheck, det                                         !
!===================================================================================================================================


   epsilon = 0

   InElementCheck = .TRUE.
   DO iLocSide = 1,6
      DO n = 1,3
         xNode(n) = GEO%NodeCoords(1,GEO%ElemSideNodeID(n,iLocSide,Element))
         yNode(n) = GEO%NodeCoords(2,GEO%ElemSideNodeID(n,iLocSide,Element))
         zNode(n) = GEO%NodeCoords(3,GEO%ElemSideNodeID(n,iLocSide,Element))
         Ax(n) = xNode(n) - PartState(i,1)
         Ay(n) = yNode(n) - PartState(i,2)
         Az(n) = zNode(n) - PartState(i,3)
      END DO

      det(iLocSide) = ((Ay(1) * Az(2) - Az(1) * Ay(2)) * Ax(3) +     &
                       (Az(1) * Ax(2) - Ax(1) * Az(2)) * Ay(3) +     &
                       (Ax(1) * Ay(2) - Ay(1) * Ax(2)) * Az(3))
      IF ((det(iLocSide).lt.-epsilon).OR.(det(iLocSide).NE.det(iLocSide))) THEN
         InElementCheck = .FALSE.
      END IF
   END DO
 RETURN
END SUBROUTINE ParticleInElementCheck3D


SUBROUTINE ParticleInsideQuad3Dold(i,Element,InElementCheck,det)                                      !
!DEC$ ATTRIBUTES FORCEINLINE :: ParticleInsideQuad3D
  USE MOD_Particle_Vars
  USE MOD_Mesh_Vars,     ONLY : ElemToSide
!--------------------------------------------------------------------------------------------------!
  IMPLICIT NONE                                                                                    !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
  INTEGER                          :: i, Element                                                   !
  LOGICAL                          :: InElementCheck                                               !
  REAL                             :: det(6,2)                                                       !
! Local variable declaration                                                                       !
  INTEGER                          :: iLocSide, NodeNum
  LOGICAL                          :: PosCheck, NegCheck                                           !
  REAL                             :: A(1:3,1:4), cross(3), PartStateLoc(3)
!--------------------------------------------------------------------------------------------------!
  INTENT(IN)                       :: Element                                                      !
  INTENT(OUT)                      :: InElementCheck, det                                          !
!--------------------------------------------------------------------------------------------------!

  InElementCheck = .TRUE.
  PartStateLoc(1:3) = PartState(i,1:3)
  DO iLocSide = 1,6                 ! for all 6 sides of the element
     !--- initialize flags for side checks
     PosCheck = .FALSE.
     NegCheck = .FALSE.
     !--- A = vector from particle to node coords
     DO NodeNum = 1,4
       A(:,NodeNum) = GEO%NodeCoords(:,GEO%ElemSideNodeID(NodeNum,iLocSide,Element)) - PartStateLoc(1:3)
     END DO

     !--- compute cross product for vector 1 and 3
     cross(1) = A(2,1) * A(3,3) - A(3,1) * A(2,3)
     cross(2) = A(3,1) * A(1,3) - A(1,1) * A(3,3)
     cross(3) = A(1,1) * A(2,3) - A(2,1) * A(1,3)

     !--- negative determinant of triangle 1 (points 1,3,2):
     det(iLocSide,1) = cross(1) * A(1,2) + &
                       cross(2) * A(2,2) + &
                       cross(3) * A(3,2)
     det(iLocSide,1) = -det(iLocSide,1)
     !--- determinant of triangle 2 (points 1,3,4):
     det(iLocSide,2) = cross(1) * A(1,4) + &
                       cross(2) * A(2,4) + &
                       cross(3) * A(3,4)
     !print*,iLocSide,det(iLocSide,:)
     IF (det(iLocSide,1).lt.0) THEN
       NegCheck = .TRUE.
     ELSE
       PosCheck = .TRUE.
     END IF
     IF (det(iLocSide,2).lt.0) THEN
       NegCheck = .TRUE.
     ELSE
       PosCheck = .TRUE.
     END IF

     !--- final determination whether particle is in element
     IF (GEO%ConcaveElemSide(iLocSide,Element)) THEN
       IF (.NOT.PosCheck) InElementCheck = .FALSE.
     ELSE
       IF (NegCheck) InElementCheck = .FALSE.
     END IF
  END DO
    !print*,InElementCheck
    !read*
 RETURN
END SUBROUTINE ParticleInsideQuad3Dold


SUBROUTINE ParticleInsideQuad3D(i,Element,InElementCheck,det)  
!DEC$ ATTRIBUTES FORCEINLINE :: ParticleInsideQuad3D
  USE MOD_Particle_Vars
  USE MOD_Particle_Surfaces_Vars, ONLY: nPartCurved, SuperSampledNodes,nTriangles
  USE MOD_Mesh_Vars,     ONLY : ElemToSide
!--------------------------------------------------------------------------------------------------!
  IMPLICIT NONE                                                                                    !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
  INTEGER                          :: i, Element!
  LOGICAL                          :: InElementCheck                                               !
    REAL                             :: det(6,ntriangles)                                                       !
! Local variable declaration                                                                       !
  INTEGER                          :: iLocSide, NodeNum,p,q,nPosCheck,nNegCheck
  LOGICAL                          :: PosCheck, NegCheck                                           !
    REAL                             :: A(1:3,1:4), cross(3), PartStateLoc(3)
  INTEGER                          :: SideID,flip
  INTEGER                          :: nlocTriangles
!--------------------------------------------------------------------------------------------------!
  INTENT(IN)                       :: Element!
  INTENT(OUT)                      :: InElementCheck, det                                          !
!--------------------------------------------------------------------------------------------------!

  InElementCheck = .TRUE.
  PartStateLoc(1:3) = PartState(i,1:3)
  det=0.

  DO iLocSide = 1,6                 ! for all 6 sides of the element
     !--- initialize flags for side checks
     PosCheck = .FALSE.
     NegCheck = .FALSE.
     nPosCheck=0
     nNegCheck=0
     nloctriangles=0
     SideID=ElemToSide(E2S_SIDE_ID,ilocSide,Element)
     flip  =ElemToSide(E2S_FLIP,ilocSide,Element)
!     print*,'nodes'
!     print*,SuperSampledNodes(1:3,0  ,0  ,SideID)
!     print*,SuperSampledNodes(1:3,0  ,nPartCurved  ,SideID)
!     print*,SuperSampledNodes(1:3,nPartCurved  ,0  ,SideID)
!     print*,SuperSampledNodes(1:3,nPartCurved  ,nPartCurved  ,SideID)
!     read*

     DO q=0,NPartCurved-1
       DO p=0,NPartCurved-1
!print*,p,q
         A(:,1)=SuperSampledNodes(1:3,p  ,q  ,SideID)-PartStateLoc(1:3)
         A(:,2)=SuperSampledNodes(1:3,p+1,q  ,SideID)-PartStateLoc(1:3)
         A(:,3)=SuperSampledNodes(1:3,p+1,q+1,SideID)-PartStateLoc(1:3)
         A(:,4)=SuperSampledNodes(1:3,p  ,q+1,SideID)-PartStateLoc(1:3)

         cross(1) = A(2,1) * A(3,3) - A(3,1) * A(2,3)
         cross(2) = A(3,1) * A(1,3) - A(1,1) * A(3,3)
         cross(3) = A(1,1) * A(2,3) - A(2,1) * A(1,3)
    
         det(iLocSide,nlocTriangles+1) = cross(1) * A(1,2) + &
                                         cross(2) * A(2,2) + &
                                         cross(3) * A(3,2)
         det(iLocSide,nlocTriangles+1) = -det(iLocSide,nlocTriangles+1)
         !--- determinant of triangle 2 (points 1,3,4):
         det(iLocSide,nlocTriangles+2) = cross(1) * A(1,4) + &
                                         cross(2) * A(2,4) + &
                                         cross(3) * A(3,4)

        !print*,iLocSide,det(iLocSide,:)
 !       IF(flip.EQ.0)THEN ! master side
        IF(flip.NE.0)THEN ! master side
          det(iLocSide,nlocTriangles+1:nlocTriangles+2)= -det(iLocSide,nlocTriangles+1:nlocTriangles+2)
        END IF
      !j  print*,det(iLocSide,nlocTriangles+1:nlocTriangles+2)
      !j  read*
        IF (det(iLocSide,nlocTriangles+1).lt.0) THEN
           nNegCheck = nNegCheck+1
         ELSE
           nPosCheck = nPosCheck+1
         END IF
         IF (det(iLocSide,nlocTriangles+2).lt.0) THEN
           nNegCheck = nNegCheck+1
         ELSE
           nPosCheck = nPosCheck+1
         END IF
 !        ELSE ! slave sides
 !         IF (det(iLocSide,nlocTriangles+1).GE.0) THEN
 !            nNegCheck = nNegCheck+1
 !          ELSE
 !            nPosCheck = nPosCheck+1
 !          END IF
 !          IF (det(iLocSide,nlocTriangles+2).GE.0) THEN
 !            nNegCheck = nNegCheck+1
 !          ELSE
 !            nPosCheck = nPosCheck+1
 !          END IF
 !        END IF ! flip
         nloctriangles=nlocTriangles+2
         END DO !p
     END DO !q

      ! print*,nNegCheck
      ! print*,nPosCheck
     !print*,'ilocside,flip,det',iLocSide,flip,det(iLocSide,:)
     IF(nNegCheck.EQ.nTriangles)THEN
       NegCheck=.TRUE.
     ELSE IF(nPosCheck.EQ.nTriangles)THEN
       PosCheck=.TRUE.
     ELSE
       NegCheck=.FALSE.
       PosCheck=.FALSE.
       !print*,'wtf'
      ! print*,'false'
       !stop
     END IF

!     IF(nPosCheck.EQ.nTriangles)THEN
!       NegCheck=.TRUE.
!     ELSE
!       PosCheck=.TRUE.
!     END IF
!     IF(nPosCheck.NE.nTriangles.AND.nNegCheck.NE.nTriangles)THEN
!       print*,'wft'
!       stop
!     END IF
!     print*,"ici"

     !--- final determination whether particle is in element
     IF (GEO%ConcaveElemSide(iLocSide,Element)) THEN
       IF (.NOT.PosCheck) InElementCheck = .FALSE.
     ELSE
       IF (NegCheck) InElementCheck = .FALSE.
     END IF
  END DO
  !stop
!print*,InElementCheck
 RETURN
  END SUBROUTINE ParticleInsideQuad3D

#ifdef MPI
SUBROUTINE ParticleInsideQuad3D_halocells(i,Element,InElementCheck,det)                                      !
!DEC$ ATTRIBUTES FORCEINLINE :: ParticleInsideQuad3D
  USE MOD_part_MPI_Vars, ONLY : MPIGEO , PMPIVAR
  USE MOD_Particle_Vars
!--------------------------------------------------------------------------------------------------!
  IMPLICIT NONE                                                                                    !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
  INTEGER                          :: i, Element                                                   !
  LOGICAL                          :: InElementCheck                                               !
  REAL                             :: det(6,2)                                                       !
! Local variable declaration                                                                       !
  INTEGER                          :: iLocSide, NodeNum
  LOGICAL                          :: PosCheck, NegCheck                                           !
  REAL                             :: A(1:3,1:4), cross(1:3)
!--------------------------------------------------------------------------------------------------!
  INTENT(IN)                       :: Element                                                      !
  INTENT(OUT)                      :: InElementCheck, det                                          !
!--------------------------------------------------------------------------------------------------!

  InElementCheck = .TRUE.
  DO iLocSide = 1,6                 ! for all 6 sides of the element
     !--- initialize flags for side checks
     PosCheck = .FALSE.
     NegCheck = .FALSE.
     !--- A = vector from particle to node coords
     DO NodeNum = 1,4
       A(:,NodeNum) = MPIGEO%NodeCoords(:,MPIGEO%ElemSideNodeID(NodeNum,iLocSide,Element)) - PartState(i,1:3)
     END DO

     !--- compute cross product for vector 1 and 3
     cross(1) = A(2,1) * A(3,3) - A(3,1) * A(2,3)
     cross(2) = A(3,1) * A(1,3) - A(1,1) * A(3,3)
     cross(3) = A(1,1) * A(2,3) - A(2,1) * A(1,3)

     !--- negative determinant of triangle 1 (points 1,3,2):
     det(iLocSide,1) = cross(1) * A(1,2) + &
                       cross(2) * A(2,2) + &
                       cross(3) * A(3,2)
     det(iLocSide,1) = -det(iLocSide,1)
     !--- determinant of triangle 2 (points 1,3,4):
     det(iLocSide,2) = cross(1) * A(1,4) + &
                       cross(2) * A(2,4) + &
                       cross(3) * A(3,4)
     IF (det(iLocSide,1).lt.0) THEN
       NegCheck = .TRUE.
     ELSE
       PosCheck = .TRUE.
     END IF
     IF (det(iLocSide,2).lt.0) THEN
       NegCheck = .TRUE.
     ELSE
       PosCheck = .TRUE.
     END IF

     !--- final determination whether particle is in element
     IF (MPIGEO%ConcaveElemSide(iLocSide,Element)) THEN
       IF (.NOT.PosCheck) InElementCheck = .FALSE.
     ELSE
       IF (NegCheck) InElementCheck = .FALSE.
     END IF
  END DO
 RETURN
END SUBROUTINE ParticleInsideQuad3D_halocells
#endif /*MPI*/

SUBROUTINE ParticleThroughSideLastPosCheck(i,iLocSide,Element,InElementCheck,TriNum,det)
!DEC$ ATTRIBUTES FORCEINLINE :: ParticleThroughSideLastPosCheck
  USE MOD_Particle_Vars
  USE MOD_Particle_Surfaces_Vars, ONLY: nPartCurved, SuperSampledNodes,nTriangles
  USE MOD_Mesh_Vars,     ONLY : ElemToSide
!--------------------------------------------------------------------------------------------------!
  IMPLICIT NONE                                                                                    !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
  INTEGER                          :: i, Element, iLocSide, TriNum
  LOGICAL                          :: InElementCheck
  REAL                             :: det
! Local variable declaration                                                                       !
  INTEGER                          :: NodeNum, ind, iNode
  REAL                             :: Ax(3),Ay(3),Az(3)
  REAL                             :: NodeCoord(1:3,1:3)
  INTEGER                          :: QuadID,p,q,SideID,flip
!--------------------------------------------------------------------------------------------------!
  INTENT(IN)                       :: Element, iLocSide, i, TriNum
  INTENT(OUT)                      :: InElementCheck, det
!--------------------------------------------------------------------------------------------------!

  InElementCheck = .TRUE.

  SideID=ElemToSide(E2S_SIDE_ID,ilocSide,Element)
  flip  =ElemToSide(E2S_FLIP,ilocSide,Element)

  !--- coords of first node:
   QuadID=triNum/2+MOD(triNum,2)
   q=(QuadID-1)/NPartCurved ! fortran takes floor of integer devision
   p=MOD(QuadID-1,NPartCurved)

!  QuadID=triNum/2+MOD(triNum,2)
!  IF(nPartCurved.EQ.1)THEN
!    q=0
!    p=0
!  ELSE
!    q=QuadID/NPartCurved-MOD(QuadID+1,2) 
!    p=QuadID-NPartCurved*q-1
!  END IF

   NodeCoord(:,1) = SuperSampledNodes(:,p  ,q  ,SideID)

  !!--- coords of other two nodes (depending on triangle):
   IF(MOD(triNum,2).EQ.0)THEN
     !print*,'tri2,p2',p+1,q+1
     NodeCoord(:,2) = SuperSampledNodes(:,p+1,q+1,SideID)

     NodeCoord(:,3) = SuperSampledNodes(:,p,q+1,SideID)
   ELSE
     NodeCoord(:,2) = SuperSampledNodes(:,p+1,q,SideID)

     NodeCoord(:,3) = SuperSampledNodes(:,p+1,q+1,SideID)

   END IF


  !DO ind = 1,3
  !  NodeCoord(ind,1) = GEO%NodeCoords(ind,GEO%ElemSideNodeID(1,iLocSide,Element))
  !END DO

  !!--- coords of other two nodes (depending on triangle):
  !DO iNode = 2,3
  !  NodeNum = iNode + TriNum - 1
  !  DO ind = 1,3
  !    NodeCoord(ind,iNode) = GEO%NodeCoords(ind,GEO%ElemSideNodeID(NodeNum,iLocSide,Element))
  !  END DO
  !END DO

  !--- vector from lastPos(!) to triangle nodes
  DO ind = 1,3
    Ax(ind) = NodeCoord(1,ind) - lastPartPos(i,1)
    Ay(ind) = NodeCoord(2,ind) - lastPartPos(i,2)
    Az(ind) = NodeCoord(3,ind) - lastPartPos(i,3)
  END DO

  !--- determine whether particle is on inner side (rel. to element) of triangle
  !--- set corresponding "flag" (see below)
  det = ((Ay(1) * Az(2) - Az(1) * Ay(2)) * Ax(3) +     &
         (Az(1) * Ax(2) - Ax(1) * Az(2)) * Ay(3) +     &
         (Ax(1) * Ay(2) - Ay(1) * Ax(2)) * Az(3))

  IF(flip.NE.0) det=-det
  IF ((det.lt.0).OR.(det.NE.det)) THEN
    InElementCheck = .FALSE.
  END IF
 RETURN
END SUBROUTINE ParticleThroughSideLastPosCheck

SUBROUTINE ParticleThroughSideLastPosCheckold(i,iLocSide,Element,InElementCheck,TriNum,det)
!DEC$ ATTRIBUTES FORCEINLINE :: ParticleThroughSideLastPosCheck
  USE MOD_Particle_Vars
!--------------------------------------------------------------------------------------------------!
  IMPLICIT NONE                                                                                    !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
  INTEGER                          :: i, Element, iLocSide, TriNum
  LOGICAL                          :: InElementCheck
  REAL                             :: det
! Local variable declaration                                                                       !
  INTEGER                          :: NodeNum, ind, iNode
  REAL                             :: Ax(3),Ay(3),Az(3)
  REAL                             :: NodeCoord(1:3,1:3)
!--------------------------------------------------------------------------------------------------!
  INTENT(IN)                       :: Element, iLocSide, i, TriNum
  INTENT(OUT)                      :: InElementCheck, det
!--------------------------------------------------------------------------------------------------!

  InElementCheck = .TRUE.

  !--- coords of first node:

  DO ind = 1,3
    NodeCoord(ind,1) = GEO%NodeCoords(ind,GEO%ElemSideNodeID(1,iLocSide,Element))
  END DO

  !--- coords of other two nodes (depending on triangle):
  DO iNode = 2,3
    NodeNum = iNode + TriNum - 1
    DO ind = 1,3
      NodeCoord(ind,iNode) = GEO%NodeCoords(ind,GEO%ElemSideNodeID(NodeNum,iLocSide,Element))
    END DO
  END DO

  !--- vector from lastPos(!) to triangle nodes
  DO ind = 1,3
    Ax(ind) = NodeCoord(1,ind) - lastPartPos(i,1)
    Ay(ind) = NodeCoord(2,ind) - lastPartPos(i,2)
    Az(ind) = NodeCoord(3,ind) - lastPartPos(i,3)
  END DO

  !--- determine whether particle is on inner side (rel. to element) of triangle
  !--- set corresponding "flag" (see below)
  det = ((Ay(1) * Az(2) - Az(1) * Ay(2)) * Ax(3) +     &
         (Az(1) * Ax(2) - Ax(1) * Az(2)) * Ay(3) +     &
         (Ax(1) * Ay(2) - Ay(1) * Ax(2)) * Az(3))

  IF ((det.lt.0).OR.(det.NE.det)) THEN
    InElementCheck = .FALSE.
  END IF
 RETURN
END SUBROUTINE ParticleThroughSideLastPosCheckold

#ifdef MPI
SUBROUTINE ParticleThroughSideLastPosCheck_halocells(i,iLocSide,Element,InElementCheck,TriNum,det)
!DEC$ ATTRIBUTES FORCEINLINE :: ParticleThroughSideLastPosCheck
  USE MOD_Particle_Vars
  USE MOD_part_MPI_Vars, ONLY : MPIGEO , PMPIVAR
!--------------------------------------------------------------------------------------------------!
  IMPLICIT NONE                                                                                    !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
  INTEGER                          :: i, Element, iLocSide, TriNum
  LOGICAL                          :: InElementCheck
  REAL                             :: det
! Local variable declaration                                                                       !
  INTEGER                          :: NodeNum, ind, iNode
  REAL                             :: Ax(3),Ay(3),Az(3)
  REAL                             :: NodeCoord(1:3,1:3)
!--------------------------------------------------------------------------------------------------!
  INTENT(IN)                       :: Element, iLocSide, i, TriNum
  INTENT(OUT)                      :: InElementCheck, det
!--------------------------------------------------------------------------------------------------!

  InElementCheck = .TRUE.

  !--- coords of first node:

  DO ind = 1,3
    NodeCoord(ind,1) = MPIGEO%NodeCoords(ind,MPIGEO%ElemSideNodeID(1,iLocSide,Element))
  END DO

  !--- coords of other two nodes (depending on triangle):
  DO iNode = 2,3
    NodeNum = iNode + TriNum - 1
    DO ind = 1,3
      NodeCoord(ind,iNode) = MPIGEO%NodeCoords(ind,MPIGEO%ElemSideNodeID(NodeNum,iLocSide,Element))
    END DO
  END DO

  !--- vector from lastPos(!) to triangle nodes
  DO ind = 1,3
    Ax(ind) = NodeCoord(1,ind) - lastPartPos(i,1)
    Ay(ind) = NodeCoord(2,ind) - lastPartPos(i,2)
    Az(ind) = NodeCoord(3,ind) - lastPartPos(i,3)
  END DO

  !--- determine whether particle is on inner side (rel. to element) of triangle
  !--- set corresponding "flag" (see below)
  det = ((Ay(1) * Az(2) - Az(1) * Ay(2)) * Ax(3) +     &
         (Az(1) * Ax(2) - Ax(1) * Az(2)) * Ay(3) +     &
         (Ax(1) * Ay(2) - Ay(1) * Ax(2)) * Az(3))

  IF ((det.lt.0).OR.(det.NE.det)) THEN
    InElementCheck = .FALSE.
  END IF
 RETURN
END SUBROUTINE ParticleThroughSideLastPosCheck_halocells
#endif

SUBROUTINE SingleParticleToExactElement(iPart)                                                         
!===================================================================================================================================
! this subroutine maps each particle to an element
! currently, a background mesh is used to find possible elements. if multiple elements are possible, the element with the smallest
! distance is picked as an initial guess
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars
USE MOD_TimeDisc_Vars,          ONLY:dt
USE MOD_Equation_Vars,          ONLY:c_inv
USE MOD_Particle_Surfaces_Vars, ONLY:epsilontol,OneMepsilon,epsilonOne,SuperSampledNodes,NPartCurved
USE MOD_Mesh_Vars,              ONLY:ElemToSide,XCL_NGeo,xBaryCL_NGeo
USE MOD_Eval_xyz,               ONLY:eval_xyz_elemcheck
USE MOD_Utils,                  ONLY:BubbleSortID
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE                                                                                   
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)                :: iPart
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: iBGMElem,nBGMElems, ElemID, CellX,CellY,CellZ,iDist                       
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                           :: ilocSide,SideID
LOGICAL                           :: InElementCheck,ParticleFound                                
REAL                              :: xi(1:3),vBary(1:3)
REAL,ALLOCATABLE                  :: Distance(:)
INTEGER,ALLOCATABLE               :: ListDistance(:)
REAL,PARAMETER                    :: eps=1e-8 ! same value as in eval_xyz_elem
REAL                              :: epsOne,OneMeps
!===================================================================================================================================

epsOne=1.0+eps
OneMeps=1.0-eps
ParticleFound = .FALSE.
IF ( (PartState(iPart,1).LT.GEO%xmin).OR.(PartState(iPart,1).GT.GEO%xmax).OR. &
     (PartState(iPart,2).LT.GEO%ymin).OR.(PartState(iPart,2).GT.GEO%ymax).OR. &
     (PartState(iPart,3).LT.GEO%zmin).OR.(PartState(iPart,3).GT.GEO%zmax)) THEN
   PDM%ParticleInside(iPart) = .FALSE.
   RETURN
END IF

! --- get background mesh cell of particle
CellX = CEILING((PartState(iPart,1)-GEO%xminglob)/GEO%FIBGMdeltas(1)) 
CellX = MIN(GEO%FIBGMimax,CellX)                             
CellY = CEILING((PartState(iPart,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))
CellY = MIN(GEO%FIBGMkmax,CellY) 
CellZ = CEILING((PartState(iPart,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))
CellZ = MIN(GEO%FIBGMlmax,CellZ)
!   print*,'cell indices',CellX,CellY,CellZ
!   print*,'number of cells in bgm',GEO%FIBGM(CellX,CellY,CellZ)%nElem
!   read*

!--- check all cells associated with this beckground mesh cell
nBGMElems=GEO%FIBGM(CellX,CellY,CellZ)%nElem
ALLOCATE( Distance(1:nBGMElems) &
        , ListDistance(1:nBGMElems) )

! get closest element barycenter
Distance=0.
ListDistance=0.
DO iBGMElem = 1, nBGMElems
  ElemID = GEO%FIBGM(CellX,CellY,CellZ)%Element(iBGMElem)
  Distance(iBGMElem)=(PartState(iPart,1)-xBaryCL_NGeo(1,ElemID))*(PartState(iPart,1)-xBaryCL_NGeo(1,ElemID)) &
                    +(PartState(iPart,2)-xBaryCL_NGeo(2,ElemID))*(PartState(iPart,2)-xBaryCL_NGeo(2,ElemID)) &
                    +(PartState(iPart,3)-xBaryCL_NGeo(3,ElemID))*(PartState(iPart,3)-xBaryCL_NGeo(3,ElemID)) 
  Distance(iBGMElem)=SQRT(Distance(iBGMElem))
  ListDistance(iBGMElem)=ElemID
END DO ! nBGMElems

!print*,'earlier',Distance,ListDistance
CALL BubbleSortID(Distance,ListDistance,nBGMElems)
!print*,'after',Distance,ListDistance
!read*

! loop through sorted list and start by closest element  
DO iBGMElem=1,nBGMElems
  ElemID=ListDistance(iBGMElem)
  CALL Eval_xyz_elemcheck(PartState(iPart,1:3),xi,ElemID)
  !print*,'xi',xi
  IF(ALL(ABS(Xi).LE.OneMEps)) THEN ! particle inside
    InElementCheck=.TRUE.
  ELSE IF(ANY(ABS(Xi).GT.epsOne))THEN ! particle outside
  !  print*,'ici'
    InElementCheck=.FALSE.
  ELSE ! particle at face,edge or node, check most possible point
    ! alter particle position
    ! 1) compute vector to cell centre
    vBary=xBaryCL_NGeo(1:3,ElemID)-PartState(iPart,1:3)
    ! 2) move particle pos along vector
    PartState(iPart,1:3) = PartState(iPart,1:3)+eps*VBary(1:3)
    CALL Eval_xyz_elemcheck(PartState(iPart,1:3),xi,ElemID)
    !print*,xi
    IF(ALL(ABS(Xi).LT.1.0)) THEN ! particle inside
      InElementCheck=.TRUE.
    ELSE
      SWRITE(*,*) ' Particle not located!'
      SWRITE(*,*) ' PartPos', PartState(iPart,1:3)
      InElementCheck=.FALSE.
    END IF
  END IF
  IF (InElementCheck) THEN !  !     print*,Element
 ! read*
    PEM%Element(iPart) = ElemID
    ParticleFound = .TRUE.
    EXIT
  END IF
END DO ! iBGMElem



! particle not found
IF (.NOT.ParticleFound) THEN
  PDM%ParticleInside(iPart) = .FALSE.
END IF
! deallocate lists
DEALLOCATE( Distance,ListDistance)
!read*
END SUBROUTINE SingleParticleToExactElement


SUBROUTINE PeriodicWallBnd3D (i,SideID)                                                            !
  USE MOD_Particle_Vars
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
   INTEGER                          :: i                                                           !
   INTEGER                          :: SideID                                                      !
! Local variable declaration                                                                       !
   REAL                             :: deltax, deltay, deltaz                                      !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

!   IF (Side%BC%BCalphaind.GT.0) THEN
!      deltax = MESH%vv(1,Side%BC%BCalphaind)
!      deltay = MESH%vv(2,Side%BC%BCalphaind)
!      deltaz = MESH%vv(3,Side%BC%BCalphaind)
!   ELSE
!      deltax = -MESH%vv(1,-Side%BC%BCalphaind)
!      deltay = -MESH%vv(2,-Side%BC%BCalphaind)
!      deltaz = -MESH%vv(3,-Side%BC%BCalphaind)
!   END IF
WRITE (*,*)'ERROR: Periodic Boundaries not implemented yet (for particles)!'
STOP

   !---- Assign new values to "old" variables to continue loop

   lastPartPos(i,1) = lastPartPos(i,1) + deltax
   lastPartPos(i,2) = lastPartPos(i,2) + deltay
   lastPartPos(i,3) = lastPartPos(i,3) + deltaz
   PartState(i,1)   = PartState(i,1) + deltax
   PartState(i,2)   = PartState(i,2) + deltay
   PartState(i,3)   = PartState(i,3) + deltaz

 RETURN
END SUBROUTINE PeriodicWallBnd3D

SUBROUTINE PerfectReflection3D (i,iLocSide,Element,TriNum, WallVelo)                                         !
  USE MOD_Particle_Vars
  USE MOD_Particle_Surfaces_Vars, ONLY: nPartCurved, SuperSampledNodes,nTriangles
  USE MOD_Mesh_Vars,     ONLY : ElemToSide
!  USE MOD_LD_Vars,    ONLY : UseLD
!  USE MOD_LD,         ONLY : LD_PerfectReflection
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
   INTEGER                          :: i                                                           !
   INTEGER                          :: iLocSide                                                    !
   INTEGER                          :: Element                                                     !
   INTEGER                          :: TriNum                                                      !
   REAL                             :: WallVelo(1:3)
! Local variable declaration                                                                       !
   INTEGER                          :: Node1, Node2                                                !
   REAL                             :: PoldX, PoldY, PoldZ, PnewX, PnewY, PnewZ, nx, ny, nz, nVal  !
   REAL                             :: bx,by,bz, ax,ay,az, dist, PoldStarX, PoldStarY, PoldStarZ   !
   REAL                             :: xNod, yNod, zNod, PnewStarX, PnewStarY, PnewStarZ, Velo     !
   REAL                             :: VelX, VelY, VelZ, NewVelocity                               !
   REAL                             :: Vector1(1:3), Vector2(1:3)                                  !
   INTEGER                          :: QuadID,p,q,SideID
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

   SideID=ElemToSide(E2S_SIDE_ID,ilocSide,Element)
   PoldX = lastPartPos(i,1)
   PoldY = lastPartPos(i,2)
   PoldZ = lastPartPos(i,3)
   PnewX = PartState(i,1)
   PnewY = PartState(i,2)
   PnewZ = PartState(i,3)

   !print*,trinum
   QuadID=triNum/2+MOD(triNum,2)
   q=(QuadID-1)/NPartCurved ! fortran takes floor of integer devision
   p=MOD(QuadID-1,NPartCurved)

!   QuadID=triNum/2+MOD(triNum,2)
!   IF(nPartCurved.EQ.1)THEN
!      q=0
!      p=0
!    ELSE
!      q=QuadID/NPartCurved-MOD(QuadID+1,2) 
!      p=QuadID-NPartCurved*q-1
!    END IF

    !print*,p,q,SideID
   xNod= SuperSampledNodes(1,p  ,q  ,SideID)
   yNod= SuperSampledNodes(2,p  ,q  ,SideID)
   zNod= SuperSampledNodes(3,p  ,q  ,SideID)

   !xNod = GEO%NodeCoords(1,GEO%ElemSideNodeID(1,iLocSide,Element))
   !yNod = GEO%NodeCoords(2,GEO%ElemSideNodeID(1,iLocSide,Element))
   !zNod = GEO%NodeCoords(3,GEO%ElemSideNodeID(1,iLocSide,Element))

   !---- Calculate normal vector:

   IF(MOD(triNum,2).EQ.0)THEN
     !print*,'tri2,p2',p+1,q+1
     Vector1(1) = SuperSampledNodes(1,p+1,q+1,SideID)-xNod
     Vector1(2) = SuperSampledNodes(2,p+1,q+1,SideID)-yNod
     Vector1(3) = SuperSampledNodes(3,p+1,q+1,SideID)-zNod

     !print*,'tri2,p3',p+m,q+n
     Vector2(1) = SuperSampledNodes(1,p,q+1,SideID)-xNod
     Vector2(2) = SuperSampledNodes(2,p,q+1,SideID)-yNod
     Vector2(3) = SuperSampledNodes(3,p,q+1,SideID)-zNod
   ELSE
     !print*,'not here'
     Vector1(1) = SuperSampledNodes(1,p+1,q+0,SideID)-xNod
     Vector1(2) = SuperSampledNodes(2,p+1,q+0,SideID)-yNod
     Vector1(3) = SuperSampledNodes(3,p+1,q+0,SideID)-zNod

     Vector2(1) = SuperSampledNodes(1,p+1,q+1,SideID)-xNod
     Vector2(2) = SuperSampledNodes(2,p+1,q+1,SideID)-yNod
     Vector2(3) = SuperSampledNodes(3,p+1,q+1,SideID)-zNod
   END IF

   !Node1 = TriNum+1     ! normal = cross product of 1-2 and 1-3 for first triangle
   !Node2 = TriNum+2     !          and 1-3 and 1-4 for second triangle

   !Vector1(1) = GEO%NodeCoords(1,GEO%ElemSideNodeID(Node1,iLocSide,Element)) - xNod
   !Vector1(2) = GEO%NodeCoords(2,GEO%ElemSideNodeID(Node1,iLocSide,Element)) - yNod
   !Vector1(3) = GEO%NodeCoords(3,GEO%ElemSideNodeID(Node1,iLocSide,Element)) - zNod

   !Vector2(1) = GEO%NodeCoords(1,GEO%ElemSideNodeID(Node2,iLocSide,Element)) - xNod
   !Vector2(2) = GEO%NodeCoords(2,GEO%ElemSideNodeID(Node2,iLocSide,Element)) - yNod
   !Vector2(3) = GEO%NodeCoords(3,GEO%ElemSideNodeID(Node2,iLocSide,Element)) - zNod

   nx = Vector1(2) * Vector2(3) - Vector1(3) * Vector2(2)
   ny = Vector1(3) * Vector2(1) - Vector1(1) * Vector2(3)
   nz = Vector1(1) * Vector2(2) - Vector1(2) * Vector2(1)

   nVal = SQRT(nx*nx + ny*ny + nz*nz)

   nx = nx/nVal
   ny = ny/nVal
   nz = nz/nVal

   !---- Calculate part. positions PoldStar and PnewStar (symmetric to side)

!--- the following has been used for impulse computations, not implemented yet?
!   IF (nx.NE.0) PIC%InverseImpulseX(i) = .NOT.(PIC%InverseImpulseX(i))
!   IF (ny.NE.0) PIC%InverseImpulseY(i) = .NOT.(PIC%InverseImpulseY(i))
!   IF (nz.NE.0) PIC%InverseImpulseZ(i) = .NOT.(PIC%InverseImpulseZ(i))

   bx = PoldX - xNod
   by = PoldY - yNod
   bz = PoldZ - zNod

   ax = bx - nx * (bx * nx + by * ny + bz * nz)
   ay = by - ny * (bx * nx + by * ny + bz * nz)
   az = bz - nz * (bx * nx + by * ny + bz * nz)

   dist = SQRT(((ay * bz - az * by) * (ay * bz - az * by) +   &
          (az * bx - ax * bz) * (az * bx - ax * bz) +   &
          (ax * by - ay * bx) * (ax * by - ay * bx))/   &
          (ax * ax + ay * ay + az * az))

   ! If vector from old point to new point goes through the node, a will be zero
   ! dist is then simply length of vector b instead of |axb|/|a|
   IF (dist.NE.dist) dist = SQRT(bx*bx+by*by+bz*bz)

   PoldStarX = PoldX + 2 * dist * nx
   PoldStarY = PoldY + 2 * dist * ny
   PoldStarZ = PoldZ + 2 * dist * nz

   bx = PnewX - xNod
   by = PnewY - yNod
   bz = PnewZ - zNod

   ax = bx - nx * (bx * nx + by * ny + bz * nz)
   ay = by - ny * (bx * nx + by * ny + bz * nz)
   az = bz - nz * (bx * nx + by * ny + bz * nz)

   dist = SQRT(((ay * bz - az * by) * (ay * bz - az * by) +   &
        (az * bx - ax * bz) * (az * bx - ax * bz) +   &
        (ax * by - ay * bx) * (ax * by - ay * bx))/   &
        (ax * ax + ay * ay + az * az))

   ! If vector from old point to new point goes through the node, a will be zero
   ! dist is then simply length of vector b instead of |axb|/|a|
   IF (dist.NE.dist) dist = SQRT(bx*bx+by*by+bz*bz)

   PnewStarX = PnewX - 2 * dist * nx
   PnewStarY = PnewY - 2 * dist * ny
   PnewStarZ = PnewZ - 2 * dist * nz

   !---- Calculate new velocity vector

   Velo = SQRT(PartState(i,4) * PartState(i,4) + &
               PartState(i,5) * PartState(i,5) + &
               PartState(i,6) * PartState(i,6))

   VelX = PnewStarX - PoldStarX
   VelY = PnewStarY - PoldStarY
   VelZ = PnewStarZ - PoldStarZ

   NewVelocity = SQRT(VelX * VelX + VelY * VelY + VelZ * VelZ)

   VelX = VelX/NewVelocity * Velo
   VelY = VelY/NewVelocity * Velo
   VelZ = VelZ/NewVelocity * Velo

   !---- Assign new values to "old" variables to continue loop

   lastPartPos(i,1) = PoldStarX
   lastPartPos(i,2) = PoldStarY
   lastPartPos(i,3) = PoldStarZ

   PartState(i,1)   = PnewStarX
   PartState(i,2)   = PnewStarY
   PartState(i,3)   = PnewStarZ

   PartState(i,4)   = VelX + WallVelo(1) 
   PartState(i,5)   = VelY + WallVelo(2)
   PartState(i,6)   = VelZ + WallVelo(3)

   !IF(UseLD) CALL LD_PerfectReflection(nx,ny,nz,xNod,yNod,zNod,PoldStarX,PoldStarY,PoldStarZ,i)

 RETURN
END SUBROUTINE PerfectReflection3D

SUBROUTINE PerfectReflection3Dold (i,iLocSide,Element,TriNum, WallVelo)                                         !
  USE MOD_Particle_Vars
  USE MOD_Particle_Surfaces_Vars, ONLY: nPartCurved, SuperSampledNodes,nTriangles
!  USE MOD_LD_Vars,    ONLY : UseLD
!  USE MOD_LD,         ONLY : LD_PerfectReflection
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
   INTEGER                          :: i                                                           !
   INTEGER                          :: iLocSide                                                    !
   INTEGER                          :: Element                                                     !
   INTEGER                          :: TriNum                                                      !
   REAL                             :: WallVelo(1:3)
! Local variable declaration                                                                       !
   INTEGER                          :: Node1, Node2                                                !
   REAL                             :: PoldX, PoldY, PoldZ, PnewX, PnewY, PnewZ, nx, ny, nz, nVal  !
   REAL                             :: bx,by,bz, ax,ay,az, dist, PoldStarX, PoldStarY, PoldStarZ   !
   REAL                             :: xNod, yNod, zNod, PnewStarX, PnewStarY, PnewStarZ, Velo     !
   REAL                             :: VelX, VelY, VelZ, NewVelocity                               !
   REAL                             :: Vector1(1:3), Vector2(1:3)                                  !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

   PoldX = lastPartPos(i,1)
   PoldY = lastPartPos(i,2)
   PoldZ = lastPartPos(i,3)
   PnewX = PartState(i,1)
   PnewY = PartState(i,2)
   PnewZ = PartState(i,3)

   print*,trinum
   xNod = GEO%NodeCoords(1,GEO%ElemSideNodeID(1,iLocSide,Element))
   yNod = GEO%NodeCoords(2,GEO%ElemSideNodeID(1,iLocSide,Element))
   zNod = GEO%NodeCoords(3,GEO%ElemSideNodeID(1,iLocSide,Element))

   !---- Calculate normal vector:

   Node1 = TriNum+1     ! normal = cross product of 1-2 and 1-3 for first triangle
   Node2 = TriNum+2     !          and 1-3 and 1-4 for second triangle

   Vector1(1) = GEO%NodeCoords(1,GEO%ElemSideNodeID(Node1,iLocSide,Element)) - xNod
   Vector1(2) = GEO%NodeCoords(2,GEO%ElemSideNodeID(Node1,iLocSide,Element)) - yNod
   Vector1(3) = GEO%NodeCoords(3,GEO%ElemSideNodeID(Node1,iLocSide,Element)) - zNod

   Vector2(1) = GEO%NodeCoords(1,GEO%ElemSideNodeID(Node2,iLocSide,Element)) - xNod
   Vector2(2) = GEO%NodeCoords(2,GEO%ElemSideNodeID(Node2,iLocSide,Element)) - yNod
   Vector2(3) = GEO%NodeCoords(3,GEO%ElemSideNodeID(Node2,iLocSide,Element)) - zNod

   nx = Vector1(2) * Vector2(3) - Vector1(3) * Vector2(2)
   ny = Vector1(3) * Vector2(1) - Vector1(1) * Vector2(3)
   nz = Vector1(1) * Vector2(2) - Vector1(2) * Vector2(1)

   nVal = SQRT(nx*nx + ny*ny + nz*nz)

   nx = nx/nVal
   ny = ny/nVal
   nz = nz/nVal

   !---- Calculate part. positions PoldStar and PnewStar (symmetric to side)

!--- the following has been used for impulse computations, not implemented yet?
!   IF (nx.NE.0) PIC%InverseImpulseX(i) = .NOT.(PIC%InverseImpulseX(i))
!   IF (ny.NE.0) PIC%InverseImpulseY(i) = .NOT.(PIC%InverseImpulseY(i))
!   IF (nz.NE.0) PIC%InverseImpulseZ(i) = .NOT.(PIC%InverseImpulseZ(i))

   bx = PoldX - xNod
   by = PoldY - yNod
   bz = PoldZ - zNod

   ax = bx - nx * (bx * nx + by * ny + bz * nz)
   ay = by - ny * (bx * nx + by * ny + bz * nz)
   az = bz - nz * (bx * nx + by * ny + bz * nz)

   dist = SQRT(((ay * bz - az * by) * (ay * bz - az * by) +   &
          (az * bx - ax * bz) * (az * bx - ax * bz) +   &
          (ax * by - ay * bx) * (ax * by - ay * bx))/   &
          (ax * ax + ay * ay + az * az))

   ! If vector from old point to new point goes through the node, a will be zero
   ! dist is then simply length of vector b instead of |axb|/|a|
   IF (dist.NE.dist) dist = SQRT(bx*bx+by*by+bz*bz)

   PoldStarX = PoldX + 2 * dist * nx
   PoldStarY = PoldY + 2 * dist * ny
   PoldStarZ = PoldZ + 2 * dist * nz

   bx = PnewX - xNod
   by = PnewY - yNod
   bz = PnewZ - zNod

   ax = bx - nx * (bx * nx + by * ny + bz * nz)
   ay = by - ny * (bx * nx + by * ny + bz * nz)
   az = bz - nz * (bx * nx + by * ny + bz * nz)

   dist = SQRT(((ay * bz - az * by) * (ay * bz - az * by) +   &
        (az * bx - ax * bz) * (az * bx - ax * bz) +   &
        (ax * by - ay * bx) * (ax * by - ay * bx))/   &
        (ax * ax + ay * ay + az * az))

   ! If vector from old point to new point goes through the node, a will be zero
   ! dist is then simply length of vector b instead of |axb|/|a|
   IF (dist.NE.dist) dist = SQRT(bx*bx+by*by+bz*bz)

   PnewStarX = PnewX - 2 * dist * nx
   PnewStarY = PnewY - 2 * dist * ny
   PnewStarZ = PnewZ - 2 * dist * nz

   !---- Calculate new velocity vector

   Velo = SQRT(PartState(i,4) * PartState(i,4) + &
               PartState(i,5) * PartState(i,5) + &
               PartState(i,6) * PartState(i,6))

   VelX = PnewStarX - PoldStarX
   VelY = PnewStarY - PoldStarY
   VelZ = PnewStarZ - PoldStarZ

   NewVelocity = SQRT(VelX * VelX + VelY * VelY + VelZ * VelZ)

   VelX = VelX/NewVelocity * Velo
   VelY = VelY/NewVelocity * Velo
   VelZ = VelZ/NewVelocity * Velo

   !---- Assign new values to "old" variables to continue loop

   lastPartPos(i,1) = PoldStarX
   lastPartPos(i,2) = PoldStarY
   lastPartPos(i,3) = PoldStarZ

   PartState(i,1)   = PnewStarX
   PartState(i,2)   = PnewStarY
   PartState(i,3)   = PnewStarZ

   PartState(i,4)   = VelX + WallVelo(1) 
   PartState(i,5)   = VelY + WallVelo(2)
   PartState(i,6)   = VelZ + WallVelo(3)

   !IF(UseLD) CALL LD_PerfectReflection(nx,ny,nz,xNod,yNod,zNod,PoldStarX,PoldStarY,PoldStarZ,i)

 RETURN
END SUBROUTINE PerfectReflection3Dold

#ifdef MPI
SUBROUTINE PerfectReflection3D_halocells(i,iLocSide,Element,TriNum, WallVelo)                                         !
  USE MOD_part_MPI_Vars, ONLY : MPIGEO 
  USE MOD_Particle_Vars, ONLY : PartState, LastPartPos
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
   INTEGER                          :: i                                                           !
   INTEGER                          :: iLocSide                                                    !
   INTEGER                          :: Element                                                     !
   INTEGER                          :: TriNum                                                      !
   REAL                             :: WallVelo(1:3)
! Local variable declaration                                                                       !
   INTEGER                          :: Node1, Node2                                                !
   REAL                             :: PoldX, PoldY, PoldZ, PnewX, PnewY, PnewZ, nx, ny, nz, nVal  !
   REAL                             :: bx,by,bz, ax,ay,az, dist, PoldStarX, PoldStarY, PoldStarZ   !
   REAL                             :: xNod, yNod, zNod, PnewStarX, PnewStarY, PnewStarZ, Velo     !
   REAL                             :: VelX, VelY, VelZ, NewVelocity                               !
   REAL                             :: Vector1(1:3), Vector2(1:3)                                  !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

   PoldX = lastPartPos(i,1)
   PoldY = lastPartPos(i,2)
   PoldZ = lastPartPos(i,3)
   PnewX = PartState(i,1)
   PnewY = PartState(i,2)
   PnewZ = PartState(i,3)

   xNod = MPIGEO%NodeCoords(1,MPIGEO%ElemSideNodeID(1,iLocSide,Element))
   yNod = MPIGEO%NodeCoords(2,MPIGEO%ElemSideNodeID(1,iLocSide,Element))
   zNod = MPIGEO%NodeCoords(3,MPIGEO%ElemSideNodeID(1,iLocSide,Element))

   !---- Calculate normal vector:

   Node1 = TriNum+1     ! normal = cross product of 1-2 and 1-3 for first triangle
   Node2 = TriNum+2     !          and 1-3 and 1-4 for second triangle

   Vector1(1) = MPIGEO%NodeCoords(1,MPIGEO%ElemSideNodeID(Node1,iLocSide,Element)) - xNod
   Vector1(2) = MPIGEO%NodeCoords(2,MPIGEO%ElemSideNodeID(Node1,iLocSide,Element)) - yNod
   Vector1(3) = MPIGEO%NodeCoords(3,MPIGEO%ElemSideNodeID(Node1,iLocSide,Element)) - zNod

   Vector2(1) = MPIGEO%NodeCoords(1,MPIGEO%ElemSideNodeID(Node2,iLocSide,Element)) - xNod
   Vector2(2) = MPIGEO%NodeCoords(2,MPIGEO%ElemSideNodeID(Node2,iLocSide,Element)) - yNod
   Vector2(3) = MPIGEO%NodeCoords(3,MPIGEO%ElemSideNodeID(Node2,iLocSide,Element)) - zNod

   nx = Vector1(2) * Vector2(3) - Vector1(3) * Vector2(2)
   ny = Vector1(3) * Vector2(1) - Vector1(1) * Vector2(3)
   nz = Vector1(1) * Vector2(2) - Vector1(2) * Vector2(1)

   nVal = SQRT(nx*nx + ny*ny + nz*nz)

   nx = nx/nVal
   ny = ny/nVal
   nz = nz/nVal

   !---- Calculate part. positions PoldStar and PnewStar (symmetric to side)

!--- the following has been used for impulse computations, not implemented yet?
!   IF (nx.NE.0) PIC%InverseImpulseX(i) = .NOT.(PIC%InverseImpulseX(i))
!   IF (ny.NE.0) PIC%InverseImpulseY(i) = .NOT.(PIC%InverseImpulseY(i))
!   IF (nz.NE.0) PIC%InverseImpulseZ(i) = .NOT.(PIC%InverseImpulseZ(i))

   bx = PoldX - xNod
   by = PoldY - yNod
   bz = PoldZ - zNod

   ax = bx - nx * (bx * nx + by * ny + bz * nz)
   ay = by - ny * (bx * nx + by * ny + bz * nz)
   az = bz - nz * (bx * nx + by * ny + bz * nz)

   dist = SQRT(((ay * bz - az * by) * (ay * bz - az * by) +   &
          (az * bx - ax * bz) * (az * bx - ax * bz) +   &
          (ax * by - ay * bx) * (ax * by - ay * bx))/   &
          (ax * ax + ay * ay + az * az))

   ! If vector from old point to new point goes through the node, a will be zero
   ! dist is then simply length of vector b instead of |axb|/|a|
   IF (dist.NE.dist) dist = SQRT(bx*bx+by*by+bz*bz)

   PoldStarX = PoldX + 2 * dist * nx
   PoldStarY = PoldY + 2 * dist * ny
   PoldStarZ = PoldZ + 2 * dist * nz

   bx = PnewX - xNod
   by = PnewY - yNod
   bz = PnewZ - zNod

   ax = bx - nx * (bx * nx + by * ny + bz * nz)
   ay = by - ny * (bx * nx + by * ny + bz * nz)
   az = bz - nz * (bx * nx + by * ny + bz * nz)

   dist = SQRT(((ay * bz - az * by) * (ay * bz - az * by) +   &
        (az * bx - ax * bz) * (az * bx - ax * bz) +   &
        (ax * by - ay * bx) * (ax * by - ay * bx))/   &
        (ax * ax + ay * ay + az * az))

   ! If vector from old point to new point goes through the node, a will be zero
   ! dist is then simply length of vector b instead of |axb|/|a|
   IF (dist.NE.dist) dist = SQRT(bx*bx+by*by+bz*bz)

   PnewStarX = PnewX - 2 * dist * nx
   PnewStarY = PnewY - 2 * dist * ny
   PnewStarZ = PnewZ - 2 * dist * nz

   !---- Calculate new velocity vector

   Velo = SQRT(PartState(i,4) * PartState(i,4) + &
               PartState(i,5) * PartState(i,5) + &
               PartState(i,6) * PartState(i,6))

   VelX = PnewStarX - PoldStarX
   VelY = PnewStarY - PoldStarY
   VelZ = PnewStarZ - PoldStarZ

   NewVelocity = SQRT(VelX * VelX + VelY * VelY + VelZ * VelZ)

   VelX = VelX/NewVelocity * Velo
   VelY = VelY/NewVelocity * Velo
   VelZ = VelZ/NewVelocity * Velo

   !---- Assign new values to "old" variables to continue loop

   lastPartPos(i,1) = PoldStarX
   lastPartPos(i,2) = PoldStarY
   lastPartPos(i,3) = PoldStarZ

   PartState(i,1)   = PnewStarX
   PartState(i,2)   = PnewStarY
   PartState(i,3)   = PnewStarZ

   PartState(i,4)   = VelX + WallVelo(1)
   PartState(i,5)   = VelY + WallVelo(2)
   PartState(i,6)   = VelZ + WallVelo(3)

 RETURN
END SUBROUTINE PerfectReflection3D_halocells
#endif /*MPI*/

SUBROUTINE DiffuseReflection3D (i,iLocSide,Element,TriNum, &
                                WallTemp,                  &
                                TransACC,                  &
                                VibACC,                    &
                                RotACC,                    &
                                WallVelo                   )                                        !
  USE MOD_Particle_Vars
  USE MOD_TimeDisc_Vars,          ONLY : dt, TEnd
  USE MOD_DSMC_Vars,              ONLY : PartStateIntEn,SpecDSMC, DSMC, SampWall, SurfMesh
  USE MOD_Mesh_Vars,              ONLY : ElemToSide
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
   INTEGER                          :: i, iTest                                                    !
   INTEGER                          :: iLocSide                                                    !
   INTEGER                          :: Element                                                     !
   INTEGER                          :: TriNum                                                      !
   REAL                             :: WallTemp                                                    !
   REAL                             :: TransACC                                                    !
   REAL                             :: VibACC                                                      !
   REAL                             :: RotACC                                                      !
   REAL                             :: WallVelo(1:3)
! Local variable declaration                                                                       !
   INTEGER                          :: Node1, Node2                                                !
   INTEGER                          :: VibQuant,VibQuantWall,VibQuantNew                           !
   REAL                             :: VibQuantNewR                                                !
   REAL                             :: PoldX, PoldY, PoldZ, PnewX, PnewY, PnewZ, nx, ny, nz, nVal  !
   REAL                             :: xNod, yNod, zNod, VeloReal, VeloCrad, VeloCx, VeloCy ,VeloCz!
   REAL                             :: EtraOld, EtraNew, RanNum, Cmr, Phi, Fak_D                   !
   REAL                             :: EtraWall, ErotWall, EvibNew, ErotNew                        !
   REAL                             :: VelX, VelY, VelZ,VecX, VecY, VecZ                           !
   REAL                             :: Vector1(1:3), Vector2(1:3)                                  !
   REAL                             :: POI_X, POI_Y, POI_Z, POI_fak                                !
   REAL, PARAMETER                  :: PI=3.14159265358979323846	                                 !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
   PoldX = lastPartPos(i,1)
   PoldY = lastPartPos(i,2)
   PoldZ = lastPartPos(i,3)
   PnewX = PartState(i,1)
   PnewY = PartState(i,2)
   PnewZ = PartState(i,3)

   VelX = PnewX - PoldX
   VelY = PnewY - PoldY
   VelZ = PnewZ - PoldZ

   xNod = GEO%NodeCoords(1,GEO%ElemSideNodeID(1,iLocSide,Element))
   yNod = GEO%NodeCoords(2,GEO%ElemSideNodeID(1,iLocSide,Element))
   zNod = GEO%NodeCoords(3,GEO%ElemSideNodeID(1,iLocSide,Element))

   !---- Calculate normal vector:

   Node1 = TriNum+1     ! normal = cross product of 1-2 and 1-3 for first triangle
   Node2 = TriNum+2     !          and 1-3 and 1-4 for second triangle

   Vector1(1) = GEO%NodeCoords(1,GEO%ElemSideNodeID(Node1,iLocSide,Element)) - xNod
   Vector1(2) = GEO%NodeCoords(2,GEO%ElemSideNodeID(Node1,iLocSide,Element)) - yNod
   Vector1(3) = GEO%NodeCoords(3,GEO%ElemSideNodeID(Node1,iLocSide,Element)) - zNod

   Vector2(1) = GEO%NodeCoords(1,GEO%ElemSideNodeID(Node2,iLocSide,Element)) - xNod
   Vector2(2) = GEO%NodeCoords(2,GEO%ElemSideNodeID(Node2,iLocSide,Element)) - yNod
   Vector2(3) = GEO%NodeCoords(3,GEO%ElemSideNodeID(Node2,iLocSide,Element)) - zNod

   nx = Vector1(2) * Vector2(3) - Vector1(3) * Vector2(2)
   ny = Vector1(3) * Vector2(1) - Vector1(1) * Vector2(3)
   nz = Vector1(1) * Vector2(2) - Vector1(2) * Vector2(1)

   nVal = SQRT(nx*nx + ny*ny + nz*nz)

   nx = nx/nVal
   ny = ny/nVal
   nz = nz/nVal

   !---- Calculate Point of Intersection (POI)

   POI_fak = (Vector2(2)*(Vector1(1)*(zNod-PoldZ)+Vector1(3)*(PoldX-xNod)) &
           +Vector1(2)*(Vector2(1)*(PoldZ-zNod)+Vector2(3)*(xNod-PoldX)) &
           +yNod*(Vector1(3)*Vector2(1)-Vector1(1)*Vector2(3)) &
           +PoldY*(Vector1(1)*Vector2(3)-Vector1(3)*Vector2(1))) &
           /(Vector1(2)*(Vector2(3)*VelX-Vector2(1)*VelZ) &
           + Vector2(2)*(Vector1(1)*VelZ-Vector1(3)*VelX) &
           + VelY*(Vector1(3)*Vector2(1)-Vector1(1)*Vector2(3)))

   POI_X = PoldX + POI_fak * VelX
   POI_Y = PoldY + POI_fak * VelY
   POI_Z = PoldZ + POI_fak * VelZ

   !---- Calculate new velocity vector (Extended Maxwellian Model)

   VeloReal = SQRT(PartState(i,4) * PartState(i,4) + &
             PartState(i,5) * PartState(i,5) + &
             PartState(i,6) * PartState(i,6))
   EtraOld = 0.5 * Species(PartSpecies(i))%MassIC * VeloReal**2
   CALL RANDOM_NUMBER(RanNum)
   VeloCrad    = SQRT(-LOG(RanNum))
   CALL RANDOM_NUMBER(RanNum)
   VeloCz      = SQRT(-LOG(RanNum))
   Fak_D       = VeloCrad**2 + VeloCz**2
   EtraWall    = BoltzmannConst * WallTemp * Fak_D
   EtraNew = EtraOld + TransACC * (EtraWall - EtraOld)
   Cmr     = SQRT(2.0 * EtraNew / (Species(PartSpecies(i))%MassIC * Fak_D))
   CALL RANDOM_NUMBER(RanNum)
   Phi     = 2 * PI * RanNum
   VeloCx  = Cmr * VeloCrad * COS(Phi)
   VeloCy  = Cmr * VeloCrad * SIN(Phi)
   VeloCz  = Cmr * VeloCz
 
   IF (DSMC%CalcSurfaceVal.AND.(Time.ge.(1-DSMC%TimeFracSamp)*TEnd)) THEN
   !----  Sampling for energy (translation) accommodation at walls

   SampWall(SurfMesh%GlobSideToSurfSideMap(ElemToSide(1,iLocSide,Element)))%Energy(1) = &
     SampWall(SurfMesh%GlobSideToSurfSideMap(ElemToSide(1,iLocSide,Element)))%Energy(1) + &
     EtraOld * Species(PartSpecies(i))%MacroParticleFactor  
   SampWall(SurfMesh%GlobSideToSurfSideMap(ElemToSide(1,iLocSide,Element)))%Energy(2) = &
     SampWall(SurfMesh%GlobSideToSurfSideMap(ElemToSide(1,iLocSide,Element)))%Energy(2) + &
     EtraWall * Species(PartSpecies(i))%MacroParticleFactor 
   SampWall(SurfMesh%GlobSideToSurfSideMap(ElemToSide(1,iLocSide,Element)))%Energy(3) = &
     SampWall(SurfMesh%GlobSideToSurfSideMap(ElemToSide(1,iLocSide,Element)))%Energy(3) + & 
     EtraNew * Species(PartSpecies(i))%MacroParticleFactor  
   END IF

   !---- Transformation local distribution -> global coordinates
    
   VecX = Vector1(1) / SQRT( Vector1(1)**2 + Vector1(2)**2 + Vector1(3)**2 )
   VecY = Vector1(2) / SQRT( Vector1(1)**2 + Vector1(2)**2 + Vector1(3)**2 )
   VecZ = Vector1(3) / SQRT( Vector1(1)**2 + Vector1(2)**2 + Vector1(3)**2 )
   
   VelX = VecX*VeloCx + (nz*VecY-ny*VecZ)*VeloCy - nx*VeloCz
   VelY = VecY*VeloCx + (nx*VecZ-nz*VecX)*VeloCy - ny*VeloCz
   VelZ = VecZ*VeloCx + (ny*VecX-nx*VecY)*VeloCy - nz*VeloCz

   lastPartPos(i,1) = POI_X
   lastPartPos(i,2) = POI_Y
   lastPartPos(i,3) = POI_Z

   PartState(i,1)   = POI_X + (1 - POI_fak) * dt * VelX
   PartState(i,2)   = POI_Y + (1 - POI_fak) * dt * VelY
   PartState(i,3)   = POI_Z + (1 - POI_fak) * dt * VelZ

   !---- Internal energy accommodation

   IF (SpecDSMC(PartSpecies(i))%InterID.EQ.2) THEN

   !---- Rotational energy accommodation

     CALL RANDOM_NUMBER(RanNum)
     ErotWall = - BoltzmannConst * WallTemp * LOG(RanNum)
     ErotNew  = PartStateIntEn(i,2) + RotACC *(ErotWall - PartStateIntEn(i,2))

     IF (DSMC%CalcSurfaceVal.AND.(Time.ge.(1-DSMC%TimeFracSamp)*TEnd)) THEN
   !----  Sampling for internal energy accommodation at walls

       SampWall(SurfMesh%GlobSideToSurfSideMap(ElemToSide(1,iLocSide,Element)))%Energy(4) = &
         SampWall(SurfMesh%GlobSideToSurfSideMap(ElemToSide(1,iLocSide,Element)))%Energy(4) &
         + PartStateIntEn(i,2) * Species(PartSpecies(i))%MacroParticleFactor   
       SampWall(SurfMesh%GlobSideToSurfSideMap(ElemToSide(1,iLocSide,Element)))%Energy(5) = &
         SampWall(SurfMesh%GlobSideToSurfSideMap(ElemToSide(1,iLocSide,Element)))%Energy(5) &
         + ErotWall * Species(PartSpecies(i))%MacroParticleFactor  
       SampWall(SurfMesh%GlobSideToSurfSideMap(ElemToSide(1,iLocSide,Element)))%Energy(6) = &
         SampWall(SurfMesh%GlobSideToSurfSideMap(ElemToSide(1,iLocSide,Element)))%Energy(6) &
         + ErotNew * Species(PartSpecies(i))%MacroParticleFactor 
     END IF  

     PartStateIntEn(i,2) = ErotNew

   !---- Vibrational energy accommodation

     VibQuant     = NINT(PartStateIntEn(i,1)/(BoltzmannConst*SpecDSMC(PartSpecies(i))%CharaTVib) &
                  - DSMC%GammaQuant)
     CALL RANDOM_NUMBER(RanNum)
     VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(PartSpecies(i))%CharaTVib)
     DO WHILE (VibQuantWall.GE.SpecDSMC(PartSpecies(i))%MaxVibQuant)
       CALL RANDOM_NUMBER(RanNum)
       VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(PartSpecies(i))%CharaTVib)
     END DO
     VibQuantNewR = VibQuant + VibACC*(VibQuantWall - VibQuant)
     VibQuantNew = INT(VibQuantNewR)
     CALL RANDOM_NUMBER(RanNum)
     IF (RanNum.LT.(VibQuantNewR - VibQuantNew)) THEN
       EvibNew = (VibQuantNew + DSMC%GammaQuant + 1.0d0)*BoltzmannConst*SpecDSMC(PartSpecies(i))%CharaTVib
     ELSE
       EvibNew = (VibQuantNew + DSMC%GammaQuant)*BoltzmannConst*SpecDSMC(PartSpecies(i))%CharaTVib
     END IF

     IF (DSMC%CalcSurfaceVal.AND.(Time.ge.(1-DSMC%TimeFracSamp)*TEnd)) THEN
   !----  Sampling for internal energy accommodation at walls
 
       SampWall(SurfMesh%GlobSideToSurfSideMap(ElemToSide(1,iLocSide,Element)))%Energy(7) = &
         SampWall(SurfMesh%GlobSideToSurfSideMap(ElemToSide(1,iLocSide,Element)))%Energy(7) + (VibQuant + DSMC%GammaQuant) &
                           * BoltzmannConst * SpecDSMC(PartSpecies(i))%CharaTVib * Species(PartSpecies(i))%MacroParticleFactor
       SampWall(SurfMesh%GlobSideToSurfSideMap(ElemToSide(1,iLocSide,Element)))%Energy(8) = &
         SampWall(SurfMesh%GlobSideToSurfSideMap(ElemToSide(1,iLocSide,Element)))%Energy(8) + VibQuantWall &
                           * BoltzmannConst * SpecDSMC(PartSpecies(i))%CharaTVib * Species(PartSpecies(i))%MacroParticleFactor
       SampWall(SurfMesh%GlobSideToSurfSideMap(ElemToSide(1,iLocSide,Element)))%Energy(9) = &
         SampWall(SurfMesh%GlobSideToSurfSideMap(ElemToSide(1,iLocSide,Element)))%Energy(9) &
                           + EvibNew * Species(PartSpecies(i))%MacroParticleFactor
     END IF

     PartStateIntEn(i,1) = EvibNew
   END IF

   IF (DSMC%CalcSurfaceVal.AND.(Time.ge.(1-DSMC%TimeFracSamp)*TEnd)) THEN
   !----  Sampling force at walls
     SampWall(SurfMesh%GlobSideToSurfSideMap(ElemToSide(1,iLocSide,Element)))%Force(1) = &
       SampWall(SurfMesh%GlobSideToSurfSideMap(ElemToSide(1,iLocSide,Element)))%Force(1) + &
       Species(PartSpecies(i))%MassIC * (PartState(i,4) - VelX) * Species(PartSpecies(i))%MacroParticleFactor
     SampWall(SurfMesh%GlobSideToSurfSideMap(ElemToSide(1,iLocSide,Element)))%Force(2) = &
       SampWall(SurfMesh%GlobSideToSurfSideMap(ElemToSide(1,iLocSide,Element)))%Force(2) + &
       Species(PartSpecies(i))%MassIC * (PartState(i,5) - VelY) * Species(PartSpecies(i))%MacroParticleFactor
     SampWall(SurfMesh%GlobSideToSurfSideMap(ElemToSide(1,iLocSide,Element)))%Force(3) = &
       SampWall(SurfMesh%GlobSideToSurfSideMap(ElemToSide(1,iLocSide,Element)))%Force(3) + &
       Species(PartSpecies(i))%MassIC * (PartState(i,6) - VelZ) * Species(PartSpecies(i))%MacroParticleFactor
   !---- Counter for collisions
     SampWall(SurfMesh%GlobSideToSurfSideMap(ElemToSide(1,iLocSide,Element)))%Counter(1) = &
       SampWall(SurfMesh%GlobSideToSurfSideMap(ElemToSide(1,iLocSide,Element)))%Counter(1) + 1
   END IF

   !----  saving new particle velocity

   PartState(i,4)   = VelX + WallVelo(1)
   PartState(i,5)   = VelY + WallVelo(2)
   PartState(i,6)   = VelZ + WallVelo(3)
  
 RETURN
END SUBROUTINE DiffuseReflection3D

#ifdef MPI
SUBROUTINE DiffuseReflection3D_halocells(i,iLocSide,Element,TriNum, &
                                         WallTemp,                  &
                                         TransACC,                  &
                                         VibACC,                    &
                                         RotACC,                    &
                                         WallVelo                   )
  USE MOD_part_MPI_Vars, ONLY : MPIGEO                                          !
  USE MOD_Particle_Vars
  USE MOD_TimeDisc_Vars,          ONLY : dt, TEnd
  USE MOD_DSMC_Vars,              ONLY : PartStateIntEn,SpecDSMC, DSMC, SampWallHaloCell, SurfMesh
  USE MOD_Mesh_Vars,              ONLY : ElemToSide
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
   INTEGER                          :: i, iTest                                                    !
   INTEGER                          :: iLocSide                                                    !
   INTEGER                          :: Element                                                     !
   INTEGER                          :: TriNum                                                      !
   REAL                             :: WallTemp                                                    !
   REAL                             :: TransACC                                                    !
   REAL                             :: VibACC                                                      !
   REAL                             :: RotACC                                                      !
   REAL                             :: WallVelo(1:3)
! Local variable declaration                                                                       !
   INTEGER                          :: Node1, Node2                                                !
   INTEGER                          :: VibQuant,VibQuantWall,VibQuantNew                           !
   REAL                             :: VibQuantNewR                                                !
   REAL                             :: PoldX, PoldY, PoldZ, PnewX, PnewY, PnewZ, nx, ny, nz, nVal  !
   REAL                             :: xNod, yNod, zNod, VeloReal, VeloCrad, VeloCx, VeloCy ,VeloCz!
   REAL                             :: EtraOld, EtraNew, RanNum, Cmr, Phi, Fak_D                   !
   REAL                             :: EtraWall, ErotWall, EvibNew, ErotNew                        !
   REAL                             :: VelX, VelY, VelZ,VecX, VecY, VecZ                           !
   REAL                             :: Vector1(1:3), Vector2(1:3)                                  !
   REAL                             :: POI_X, POI_Y, POI_Z, POI_fak                                !
   REAL, PARAMETER                  :: PI=3.14159265358979323846	                                 !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
   PoldX = lastPartPos(i,1)
   PoldY = lastPartPos(i,2)
   PoldZ = lastPartPos(i,3)
   PnewX = PartState(i,1)
   PnewY = PartState(i,2)
   PnewZ = PartState(i,3)

   VelX = PnewX - PoldX
   VelY = PnewY - PoldY
   VelZ = PnewZ - PoldZ

   xNod = MPIGEO%NodeCoords(1,MPIGEO%ElemSideNodeID(1,iLocSide,Element))
   yNod = MPIGEO%NodeCoords(2,MPIGEO%ElemSideNodeID(1,iLocSide,Element))
   zNod = MPIGEO%NodeCoords(3,MPIGEO%ElemSideNodeID(1,iLocSide,Element))

   !---- Calculate normal vector:

   Node1 = TriNum+1     ! normal = cross product of 1-2 and 1-3 for first triangle
   Node2 = TriNum+2     !          and 1-3 and 1-4 for second triangle

   Vector1(1) = MPIGEO%NodeCoords(1,MPIGEO%ElemSideNodeID(Node1,iLocSide,Element)) - xNod
   Vector1(2) = MPIGEO%NodeCoords(2,MPIGEO%ElemSideNodeID(Node1,iLocSide,Element)) - yNod
   Vector1(3) = MPIGEO%NodeCoords(3,MPIGEO%ElemSideNodeID(Node1,iLocSide,Element)) - zNod

   Vector2(1) = MPIGEO%NodeCoords(1,MPIGEO%ElemSideNodeID(Node2,iLocSide,Element)) - xNod
   Vector2(2) = MPIGEO%NodeCoords(2,MPIGEO%ElemSideNodeID(Node2,iLocSide,Element)) - yNod
   Vector2(3) = MPIGEO%NodeCoords(3,MPIGEO%ElemSideNodeID(Node2,iLocSide,Element)) - zNod

   nx = Vector1(2) * Vector2(3) - Vector1(3) * Vector2(2)
   ny = Vector1(3) * Vector2(1) - Vector1(1) * Vector2(3)
   nz = Vector1(1) * Vector2(2) - Vector1(2) * Vector2(1)

   nVal = SQRT(nx*nx + ny*ny + nz*nz)

   nx = nx/nVal
   ny = ny/nVal
   nz = nz/nVal

   !---- Calculate Point of Intersection (POI)

   POI_fak = (Vector2(2)*(Vector1(1)*(zNod-PoldZ)+Vector1(3)*(PoldX-xNod)) &
           +Vector1(2)*(Vector2(1)*(PoldZ-zNod)+Vector2(3)*(xNod-PoldX)) &
           +yNod*(Vector1(3)*Vector2(1)-Vector1(1)*Vector2(3)) &
           +PoldY*(Vector1(1)*Vector2(3)-Vector1(3)*Vector2(1))) &
           /(Vector1(2)*(Vector2(3)*VelX-Vector2(1)*VelZ) &
           + Vector2(2)*(Vector1(1)*VelZ-Vector1(3)*VelX) &
           + VelY*(Vector1(3)*Vector2(1)-Vector1(1)*Vector2(3)))

   POI_X = PoldX + POI_fak * VelX
   POI_Y = PoldY + POI_fak * VelY
   POI_Z = PoldZ + POI_fak * VelZ

   !---- Calculate new velocity vector (Extended Maxwellian Model)

   VeloReal = SQRT(PartState(i,4) * PartState(i,4) + &
             PartState(i,5) * PartState(i,5) + &
             PartState(i,6) * PartState(i,6))
   EtraOld = 0.5 * Species(PartSpecies(i))%MassIC * VeloReal**2
   CALL RANDOM_NUMBER(RanNum)
   VeloCrad    = SQRT(-LOG(RanNum))
   CALL RANDOM_NUMBER(RanNum)
   VeloCz      = SQRT(-LOG(RanNum))
   Fak_D       = VeloCrad**2 + VeloCz**2
   EtraWall    = BoltzmannConst * WallTemp * Fak_D
   EtraNew = EtraOld + TransACC * (EtraWall - EtraOld)
   Cmr     = SQRT(2.0 * EtraNew / (Species(PartSpecies(i))%MassIC * Fak_D))
   CALL RANDOM_NUMBER(RanNum)
   Phi     = 2 * PI * RanNum
   VeloCx  = Cmr * VeloCrad * COS(Phi)
   VeloCy  = Cmr * VeloCrad * SIN(Phi)
   VeloCz  = Cmr * VeloCz

   IF (DSMC%CalcSurfaceVal.AND.(Time.ge.(1-DSMC%TimeFracSamp)*TEnd)) THEN
   !----  Sampling for energy (translation) accommodation at walls

   SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(1) = &
     SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(1) + &
     EtraOld * Species(PartSpecies(i))%MacroParticleFactor  
   SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(2) = &
     SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(2) + &
     EtraWall * Species(PartSpecies(i))%MacroParticleFactor 
   SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(3) = &
     SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(3) + & 
     EtraNew * Species(PartSpecies(i))%MacroParticleFactor  
   END IF
    
!     Transformation local distribution -> global coordinates
    
   VecX = Vector1(1) / SQRT( Vector1(1)**2 + Vector1(2)**2 + Vector1(3)**2 )
   VecY = Vector1(2) / SQRT( Vector1(1)**2 + Vector1(2)**2 + Vector1(3)**2 )
   VecZ = Vector1(3) / SQRT( Vector1(1)**2 + Vector1(2)**2 + Vector1(3)**2 )
   
   VelX = VecX*VeloCx + (nz*VecY-ny*VecZ)*VeloCy - nx*VeloCz
   VelY = VecY*VeloCx + (nx*VecZ-nz*VecX)*VeloCy - ny*VeloCz
   VelZ = VecZ*VeloCx + (ny*VecX-nx*VecY)*VeloCy - nz*VeloCz

   lastPartPos(i,1) = POI_X
   lastPartPos(i,2) = POI_Y
   lastPartPos(i,3) = POI_Z

   PartState(i,1)   = POI_X + (1 - POI_fak) * dt * VelX
   PartState(i,2)   = POI_Y + (1 - POI_fak) * dt * VelY
   PartState(i,3)   = POI_Z + (1 - POI_fak) * dt * VelZ

   !---- Internal energy accommodation

   IF (SpecDSMC(PartSpecies(i))%InterID.EQ.2) THEN

   !---- Rotational energy accommodation

     CALL RANDOM_NUMBER(RanNum)
     ErotWall = - BoltzmannConst * WallTemp * LOG(RanNum)
     ErotNew  = PartStateIntEn(i,2) + RotACC *(ErotWall - PartStateIntEn(i,2))

     IF (DSMC%CalcSurfaceVal.AND.(Time.ge.(1-DSMC%TimeFracSamp)*TEnd)) THEN
   !----  Sampling for internal energy accommodation at walls

       SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(4) = &
         SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(4) &
         + PartStateIntEn(i,2) * Species(PartSpecies(i))%MacroParticleFactor   
       SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(5) = &
         SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(5) &
         + ErotWall * Species(PartSpecies(i))%MacroParticleFactor  
       SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(6) = &
         SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(6) &
         + ErotNew * Species(PartSpecies(i))%MacroParticleFactor 
     END IF 

     PartStateIntEn(i,2) = ErotNew

   !---- Vibrational energy accommodation

     VibQuant     = NINT(PartStateIntEn(i,1)/(BoltzmannConst*SpecDSMC(PartSpecies(i))%CharaTVib) &
                  - DSMC%GammaQuant)
     CALL RANDOM_NUMBER(RanNum)
     VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(PartSpecies(i))%CharaTVib)
     DO WHILE (VibQuantWall.GE.SpecDSMC(PartSpecies(i))%MaxVibQuant)
       CALL RANDOM_NUMBER(RanNum)
       VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(PartSpecies(i))%CharaTVib)
     END DO
     VibQuantNewR = VibQuant + VibACC*(VibQuantWall - VibQuant)
     VibQuantNew = INT(VibQuantNewR)
     CALL RANDOM_NUMBER(RanNum)
     IF (RanNum.LT.(VibQuantNewR - VibQuantNew)) THEN
       EvibNew = (VibQuantNew + DSMC%GammaQuant + 1.0d0)*BoltzmannConst*SpecDSMC(PartSpecies(i))%CharaTVib
     ELSE
       EvibNew = (VibQuantNew + DSMC%GammaQuant)*BoltzmannConst*SpecDSMC(PartSpecies(i))%CharaTVib
     END IF

     IF (DSMC%CalcSurfaceVal.AND.(Time.ge.(1-DSMC%TimeFracSamp)*TEnd)) THEN
   !----  Sampling for internal energy accommodation at walls
 
       SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(7) = &
         SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(7) & 
                           + (VibQuant + DSMC%GammaQuant) &
                           * BoltzmannConst * SpecDSMC(PartSpecies(i))%CharaTVib * Species(PartSpecies(i))%MacroParticleFactor
       SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(8) = &
         SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(8) + VibQuantWall &
                           * BoltzmannConst * SpecDSMC(PartSpecies(i))%CharaTVib * Species(PartSpecies(i))%MacroParticleFactor
       SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(9) = &
         SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(9) &
                           + EvibNew * Species(PartSpecies(i))%MacroParticleFactor
     END IF

     PartStateIntEn(i,1) = EvibNew
   END IF

   IF (DSMC%CalcSurfaceVal.AND.(Time.ge.(1-DSMC%TimeFracSamp)*TEnd)) THEN
   !----  Sampling force at walls
     SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Force(1) = &
       SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Force(1) + &
       Species(PartSpecies(i))%MassIC * (PartState(i,4) - VelX) * Species(PartSpecies(i))%MacroParticleFactor
     SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Force(2) = &
       SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Force(2) + &
       Species(PartSpecies(i))%MassIC * (PartState(i,5) - VelY) * Species(PartSpecies(i))%MacroParticleFactor
     SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Force(3) = &
       SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Force(3) + &
       Species(PartSpecies(i))%MassIC * (PartState(i,6) - VelZ) * Species(PartSpecies(i))%MacroParticleFactor
   !---- Counter for collisions
     SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Counter(1) = &
       SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Counter(1) + 1
   END IF

   !----  saving new particle velocity

   PartState(i,4)   = VelX + WallVelo(1)
   PartState(i,5)   = VelY + WallVelo(2)
   PartState(i,6)   = VelZ + WallVelo(3)

 RETURN
END SUBROUTINE DiffuseReflection3D_halocells
#endif /*MPI*/

END MODULE MOD_BoundaryTools
