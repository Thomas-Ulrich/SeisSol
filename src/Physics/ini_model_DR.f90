!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Alice Gabriel (gabriel AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/gabriel)
!!
!! @section LICENSE
!! Copyright (c) 2007-2016, SeisSol Group
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions are met:
!!
!! 1. Redistributions of source code must retain the above copyright notice,
!!    this list of conditions and the following disclaimer.
!!
!! 2. Redistributions in binary form must reproduce the above copyright notice,
!!    this list of conditions and the following disclaimer in the documentation
!!    and/or other materials provided with the distribution.
!!
!! 3. Neither the name of the copyright holder nor the names of its
!!    contributors may be used to endorse or promote products derived from this
!!    software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
!! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
!! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.
!!
!! @section DESCRIPTION
!! Module containing Dynamic Rupture initial model setups
!! includes background stress and nucleation types
!!
!! Can be edited by users: Please add your models as a new subroutines

#ifdef BG
#include "../Initializer/preProcessorMacros.fpp"
#else
#include "Initializer/preProcessorMacros.fpp"
#endif

MODULE ini_model_DR_mod
  !---------------------------------------------------------------------------!
  USE TypesDef
  USE DGBasis_mod
  USE read_backgroundstress_mod
  !use StressReader
  !---------------------------------------------------------------------------!
  IMPLICIT NONE
  PRIVATE
  !---------------------------------------------------------------------------!
  REAL, PARAMETER :: ZERO = 0.0D0
  !---------------------------------------------------------------------------!
  INTERFACE DR_setup
     MODULE PROCEDURE DR_setup
  END INTERFACE
  !---------------------------------------------------------------------------!
  PUBLIC  :: DR_setup
  PRIVATE :: DR_basic_ini
  private :: rotateStressToFaultCS
  !---------------------------------------------------------------------------!
  PRIVATE :: friction_RSF34
  PRIVATE :: friction_RSF7
  PRIVATE :: friction_RSF101
  PRIVATE :: friction_RSF103
  PRIVATE :: friction_LSW
  PRIVATE :: friction_LSW6
  !---------------------------------------------------------------------------!

  CONTAINS

  !> Interface to dynamic rupture initial models
  !<
  SUBROUTINE DR_setup(EQN,DISC,MESH,IO,BND)
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    TYPE(tEquations), target       :: EQN
    TYPE(tDiscretization), target  :: DISC
    TYPE(tUnstructMesh)            :: MESH
    TYPE(tInputOutput)             :: IO
    TYPE (tBoundary)               :: BND
    !-------------------------------------------------------------------------!
    INTENT(IN)                      :: MESH, BND
    INTENT(INOUT)                   :: IO, EQN, DISC
    ! -------------------------------------------------------------------------

    ! Basic DR setup, valid for all models
    CALL DR_basic_ini(DISC,EQN,MESH,BND)
    !-------------------------------------------------------------------------!

    ! Initialize model dependent (space dependent) friction law parameters
    SELECT CASE(EQN%FL)
    CASE(2,16)
      ! Initialization of friction for linear slip weakening
      CALL friction_LSW(DISC,EQN,MESH,BND)
    CASE(3,4)
      ! Initialization of initial slip rate and friction for rate and state friction
      CALL friction_RSF34(DISC,EQN,MESH,BND)
    CASE(6)
      ! Initialization of friction and fault strength for bi-material linear slip weakening
      CALL friction_LSW6(DISC,EQN,MESH,BND)
    CASE(7)
      ! Initialization of initial slip rate and friction for fast velocity weakening friction
      CALL friction_RSF7(DISC,EQN,MESH,BND)
    CASE(101)
     ! Initialization of initial slip rate and friction for SCEC TPV103
     CALL friction_RSF101(DISC,EQN,MESH,BND)
    CASE(103)
     ! Initialization of initial slip rate and friction for SCEC TPV103
     CALL friction_RSF103(DISC,EQN,MESH,BND)
    END SELECT  ! Initialize model dependent rate-and-state friction law parameters type

  END SUBROUTINE DR_setup


  !> Initialization of basic dynamic rupture setup
  !<
  SUBROUTINE DR_basic_ini(DISC,EQN,MESH,BND)                 ! global variables
    use JacobiNormal_mod, only: RotationMatrix3D
    use f_ftoc_bind_interoperability
    use iso_c_binding, only: c_null_char, c_bool
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    TYPE(tEquations)               :: EQN
    TYPE(tDiscretization), target  :: DISC
    TYPE(tUnstructMesh)            :: MESH
    TYPE (tBoundary)               :: BND
    !-------------------------------------------------------------------------!
    ! Local variable declaration
    integer                             :: i
    real, allocatable, dimension(:,:)   :: nuc_xx,nuc_yy,nuc_zz,nuc_xy,nuc_yz,nuc_xz
    logical(kind=c_bool)                :: faultParameterizedByTraction
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: MESH, BND
    INTENT(INOUT) :: EQN,DISC
    !-------------------------------------------------------------------------!

    ! Allocation of DR fields
    ALLOCATE(  EQN%IniMu(DISC%Galerkin%nBndGP,MESH%Fault%nSide),            &
               EQN%IniBulk_xx(DISC%Galerkin%nBndGP,MESH%Fault%nSide),       &
               EQN%IniBulk_yy(DISC%Galerkin%nBndGP,MESH%Fault%nSide),       &
               EQN%IniBulk_zz(DISC%Galerkin%nBndGP,MESH%Fault%nSide),       &
               EQN%IniStateVar(DISC%Galerkin%nBndGP,MESH%Fault%nSide),      &
               EQN%IniShearXY(DISC%Galerkin%nBndGP,MESH%Fault%nSide),       &
               EQN%IniShearYZ(DISC%Galerkin%nBndGP,MESH%Fault%nSide),       &
               EQN%IniShearXZ(DISC%Galerkin%nBndGP,MESH%Fault%nSide)        )
    ALLOCATE(  DISC%DynRup%Strength(DISC%Galerkin%nBndGP,MESH%Fault%nSide)  )
    ALLOCATE(  DISC%DynRup%RF(DISC%Galerkin%nBndGP,MESH%Fault%nSide)        )
    ALLOCATE(  DISC%DynRup%DS(DISC%Galerkin%nBndGP,MESH%Fault%nSide)        )
    ALLOCATE(  DISC%DynRup%cohesion(DISC%Galerkin%nBndGP,MESH%Fault%nSide)  )

    ALLOCATE(  DISC%DynRup%filter(DISC%Galerkin%nBndGP, DISC%Galerkin%nBndGP))

    IF (DISC%Galerkin%nBndGP.EQ.16) THEN
      DISC%DynRup%filter(1,1) = 0.15642717369410010645D0
      DISC%DynRup%filter(1,2) = 0.27464400686739390330D0
      DISC%DynRup%filter(1,3) = 0.25043771334395920740D0
      DISC%DynRup%filter(1,4) = 0.12372329374147723503D0
      DISC%DynRup%filter(1,5) = 0.22191620992279913356D0
      DISC%DynRup%filter(1,6) = 0.27426709126275440381D0
      DISC%DynRup%filter(1,7) = 0.92520038728438841932D-1
      DISC%DynRup%filter(1,8) = -0.23632898290357492543D-1
      DISC%DynRup%filter(1,9) = -0.15218840265827830345D-1
      DISC%DynRup%filter(1,10) = -0.14503490290801995555D0
      DISC%DynRup%filter(1,11) = -0.28178784003296133999D0
      DISC%DynRup%filter(1,12) = -0.19997871003163378446D0
      DISC%DynRup%filter(1,13) = -0.18280149283407890828D-1
      DISC%DynRup%filter(1,14) = 0.32078300425808074406D-1
      DISC%DynRup%filter(1,15) = 0.13591339660113770065D0
      DISC%DynRup%filter(1,16) = 0.12200611622433968643D0
      DISC%DynRup%filter(2,1) = 0.14649537417313707426D0
      DISC%DynRup%filter(2,2) = 0.26735493420410152682D0
      DISC%DynRup%filter(2,3) = 0.25779816505778056539D0
      DISC%DynRup%filter(2,4) = 0.13358371421191128566D0
      DISC%DynRup%filter(2,5) = 0.14629432703155492513D0
      DISC%DynRup%filter(2,6) = 0.22059025444176653702D0
      DISC%DynRup%filter(2,7) = 0.14883558359911262259D0
      DISC%DynRup%filter(2,8) = 0.49350276551200802096D-1
      DISC%DynRup%filter(2,9) = -0.77361754993378137306D-1
      DISC%DynRup%filter(2,10) = -0.18018093697265854643D0
      DISC%DynRup%filter(2,11) = -0.23417169970215916505D0
      DISC%DynRup%filter(2,12) = -0.15030590157024706139D0
      DISC%DynRup%filter(2,13) = 0.17110595921308508773D-1
      DISC%DynRup%filter(2,14) = 0.70558051658065245595D-1
      DISC%DynRup%filter(2,15) = 0.11155268152458996386D0
      DISC%DynRup%filter(2,16) = 0.72496334863913852336D-1
      DISC%DynRup%filter(3,1) = 0.13358371421191128567D0
      DISC%DynRup%filter(3,2) = 0.25779816505778056538D0
      DISC%DynRup%filter(3,3) = 0.26735493420410152684D0
      DISC%DynRup%filter(3,4) = 0.14649537417313707425D0
      DISC%DynRup%filter(3,5) = 0.49350276551200802111D-1
      DISC%DynRup%filter(3,6) = 0.14883558359911262262D0
      DISC%DynRup%filter(3,7) = 0.22059025444176653698D0
      DISC%DynRup%filter(3,8) = 0.14629432703155492511D0
      DISC%DynRup%filter(3,9) = -0.15030590157024706139D0
      DISC%DynRup%filter(3,10) = -0.23417169970215916499D0
      DISC%DynRup%filter(3,11) = -0.18018093697265854646D0
      DISC%DynRup%filter(3,12) = -0.77361754993378137338D-1
      DISC%DynRup%filter(3,13) = 0.72496334863913852338D-1
      DISC%DynRup%filter(3,14) = 0.11155268152458996385D0
      DISC%DynRup%filter(3,15) = 0.70558051658065245653D-1
      DISC%DynRup%filter(3,16) = 0.17110595921308508782D-1
      DISC%DynRup%filter(4,1) = 0.12372329374147723503D0
      DISC%DynRup%filter(4,2) = 0.25043771334395920738D0
      DISC%DynRup%filter(4,3) = 0.27464400686739390329D0
      DISC%DynRup%filter(4,4) = 0.15642717369410010644D0
      DISC%DynRup%filter(4,5) = -0.23632898290357492489D-1
      DISC%DynRup%filter(4,6) = 0.92520038728438842004D-1
      DISC%DynRup%filter(4,7) = 0.27426709126275440381D0
      DISC%DynRup%filter(4,8) = 0.22191620992279913355D0
      DISC%DynRup%filter(4,9) = -0.19997871003163378440D0
      DISC%DynRup%filter(4,10) = -0.28178784003296133985D0
      DISC%DynRup%filter(4,11) = -0.14503490290801995553D0
      DISC%DynRup%filter(4,12) = -0.15218840265827830338D-1
      DISC%DynRup%filter(4,13) = 0.12200611622433968641D0
      DISC%DynRup%filter(4,14) = 0.13591339660113770057D0
      DISC%DynRup%filter(4,15) = 0.32078300425808074402D-1
      DISC%DynRup%filter(4,16) = -0.18280149283407890807D-1
      DISC%DynRup%filter(5,1) = 0.53289900521001361894D-1
      DISC%DynRup%filter(5,2) = 0.65861191549103713634D-1
      DISC%DynRup%filter(5,3) = 0.22217320950789920454D-1
      DISC%DynRup%filter(5,4) = -0.56750915102335963452D-2
      DISC%DynRup%filter(5,5) = 0.34022128530768793354D0
      DISC%DynRup%filter(5,6) = 0.35817824692139797986D0
      DISC%DynRup%filter(5,7) = 0.21954007630124436518D-1
      DISC%DynRup%filter(5,8) = -0.11403401603774197004D0
      DISC%DynRup%filter(5,9) = 0.31353865609575951337D0
      DISC%DynRup%filter(5,10) = 0.25752454437518961863D0
      DISC%DynRup%filter(5,11) = -0.38018265592858047193D-1
      DISC%DynRup%filter(5,12) = -0.85754055475244188938D-1
      DISC%DynRup%filter(5,13) = -0.78971837571295276373D-1
      DISC%DynRup%filter(5,14) = -0.16692809384499648324D0
      DISC%DynRup%filter(5,15) = -0.38263746835464955940D-1
      DISC%DynRup%filter(5,16) = 0.94859953516780039538D-1
      DISC%DynRup%filter(6,1) = 0.35130422140007264244D-1
      DISC%DynRup%filter(6,2) = 0.52971491894141318636D-1
      DISC%DynRup%filter(6,3) = 0.35740667375047113459D-1
      DISC%DynRup%filter(6,4) = 0.11850740101465703301D-1
      DISC%DynRup%filter(6,5) = 0.19105261717494247462D0
      DISC%DynRup%filter(6,6) = 0.26814983829069693035D0
      DISC%DynRup%filter(6,7) = 0.13540677982975222973D0
      DISC%DynRup%filter(6,8) = 0.11710288526076745226D-1
      DISC%DynRup%filter(6,9) = 0.13736383661641423648D0
      DISC%DynRup%filter(6,10) = 0.22344392138011551347D0
      DISC%DynRup%filter(6,11) = 0.10676210039065944372D0
      DISC%DynRup%filter(6,12) = -0.20278978984342297582D-1
      DISC%DynRup%filter(6,13) = -0.89039603837554314230D-1
      DISC%DynRup%filter(6,14) = -0.65325773873654813997D-1
      DISC%DynRup%filter(6,15) = -0.14528428111413097370D-1
      DISC%DynRup%filter(6,16) = -0.20409918912354450312D-1
      DISC%DynRup%filter(7,1) = 0.11850740101465703292D-1
      DISC%DynRup%filter(7,2) = 0.35740667375047113450D-1
      DISC%DynRup%filter(7,3) = 0.52971491894141318630D-1
      DISC%DynRup%filter(7,4) = 0.35130422140007264243D-1
      DISC%DynRup%filter(7,5) = 0.11710288526076745218D-1
      DISC%DynRup%filter(7,6) = 0.13540677982975222973D0
      DISC%DynRup%filter(7,7) = 0.26814983829069693040D0
      DISC%DynRup%filter(7,8) = 0.19105261717494247465D0
      DISC%DynRup%filter(7,9) = -0.20278978984342297593D-1
      DISC%DynRup%filter(7,10) = 0.10676210039065944379D0
      DISC%DynRup%filter(7,11) = 0.22344392138011551355D0
      DISC%DynRup%filter(7,12) = 0.13736383661641423649D0
      DISC%DynRup%filter(7,13) = -0.20409918912354450315D-1
      DISC%DynRup%filter(7,14) = -0.14528428111413097363D-1
      DISC%DynRup%filter(7,15) = -0.65325773873654814011D-1
      DISC%DynRup%filter(7,16) = -0.89039603837554314263D-1
      DISC%DynRup%filter(8,1) = -0.56750915102335963582D-2
      DISC%DynRup%filter(8,2) = 0.22217320950789920446D-1
      DISC%DynRup%filter(8,3) = 0.65861191549103713624D-1
      DISC%DynRup%filter(8,4) = 0.53289900521001361893D-1
      DISC%DynRup%filter(8,5) = -0.11403401603774197004D0
      DISC%DynRup%filter(8,6) = 0.21954007630124436534D-1
      DISC%DynRup%filter(8,7) = 0.35817824692139797993D0
      DISC%DynRup%filter(8,8) = 0.34022128530768793360D0
      DISC%DynRup%filter(8,9) = -0.85754055475244188938D-1
      DISC%DynRup%filter(8,10) = -0.38018265592858047067D-1
      DISC%DynRup%filter(8,11) = 0.25752454437518961879D0
      DISC%DynRup%filter(8,12) = 0.31353865609575951342D0
      DISC%DynRup%filter(8,13) = 0.94859953516780039547D-1
      DISC%DynRup%filter(8,14) = -0.38263746835464955931D-1
      DISC%DynRup%filter(8,15) = -0.16692809384499648328D0
      DISC%DynRup%filter(8,16) = -0.78971837571295276408D-1
      DISC%DynRup%filter(9,1) = -0.23322892082007864531D-2
      DISC%DynRup%filter(9,2) = -0.22226617334591226444D-1
      DISC%DynRup%filter(9,3) = -0.43184022358574655443D-1
      DISC%DynRup%filter(9,4) = -0.30646762771008261828D-1
      DISC%DynRup%filter(9,5) = 0.20009491565312789465D0
      DISC%DynRup%filter(9,6) = 0.16434768403684776945D0
      DISC%DynRup%filter(9,7) = -0.24262595693329025957D-1
      DISC%DynRup%filter(9,8) = -0.54726746331374303345D-1
      DISC%DynRup%filter(9,9) = 0.42753068173757869379D0
      DISC%DynRup%filter(9,10) = 0.26033788677614015865D0
      DISC%DynRup%filter(9,11) = -0.37428258638082778209D-1
      DISC%DynRup%filter(9,12) = 0.25234135834989557741D-1
      DISC%DynRup%filter(9,13) = 0.23125258610935440010D0
      DISC%DynRup%filter(9,14) = -0.16239363945634718633D-1
      DISC%DynRup%filter(9,15) = -0.14076441672720108739D0
      DISC%DynRup%filter(9,16) = 0.63013182859958368461D-1
      DISC%DynRup%filter(10,1) = -0.11855698801416825724D-1
      DISC%DynRup%filter(10,2) = -0.27612751529328107880D-1
      DISC%DynRup%filter(10,3) = -0.35886842791018210509D-1
      DISC%DynRup%filter(10,4) = -0.23034398550611786050D-1
      DISC%DynRup%filter(10,5) = 0.87663210794514795800D-1
      DISC%DynRup%filter(10,6) = 0.14259802334581429403D0
      DISC%DynRup%filter(10,7) = 0.68133715117077374507D-1
      DISC%DynRup%filter(10,8) = -0.12941691592134129442D-1
      DISC%DynRup%filter(10,9) = 0.13886447612573873859D0
      DISC%DynRup%filter(10,10) = 0.33716692177727082467D0
      DISC%DynRup%filter(10,11) = 0.21960731652177836357D0
      DISC%DynRup%filter(10,12) = -0.19964268714162294378D-1
      DISC%DynRup%filter(10,13) = -0.86620921558945794733D-2
      DISC%DynRup%filter(10,14) = 0.13508554258884769330D0
      DISC%DynRup%filter(10,15) = 0.85922411471118344831D-1
      DISC%DynRup%filter(10,16) = -0.75083873607594495806D-1
      DISC%DynRup%filter(11,1) = -0.23034398550611786062D-1
      DISC%DynRup%filter(11,2) = -0.35886842791018210518D-1
      DISC%DynRup%filter(11,3) = -0.27612751529328107885D-1
      DISC%DynRup%filter(11,4) = -0.11855698801416825722D-1
      DISC%DynRup%filter(11,5) = -0.12941691592134129485D-1
      DISC%DynRup%filter(11,6) = 0.68133715117077374465D-1
      DISC%DynRup%filter(11,7) = 0.14259802334581429408D0
      DISC%DynRup%filter(11,8) = 0.87663210794514795854D-1
      DISC%DynRup%filter(11,9) = -0.19964268714162294389D-1
      DISC%DynRup%filter(11,10) = 0.21960731652177836357D0
      DISC%DynRup%filter(11,11) = 0.33716692177727082473D0
      DISC%DynRup%filter(11,12) = 0.13886447612573873864D0
      DISC%DynRup%filter(11,13) = -0.75083873607594495780D-1
      DISC%DynRup%filter(11,14) = 0.85922411471118344831D-1
      DISC%DynRup%filter(11,15) = 0.13508554258884769326D0
      DISC%DynRup%filter(11,16) = -0.86620921558945794850D-2
      DISC%DynRup%filter(12,1) = -0.30646762771008261838D-1
      DISC%DynRup%filter(12,2) = -0.43184022358574655446D-1
      DISC%DynRup%filter(12,3) = -0.22226617334591226454D-1
      DISC%DynRup%filter(12,4) = -0.23322892082007864531D-2
      DISC%DynRup%filter(12,5) = -0.54726746331374303343D-1
      DISC%DynRup%filter(12,6) = -0.24262595693329025948D-1
      DISC%DynRup%filter(12,7) = 0.16434768403684776946D0
      DISC%DynRup%filter(12,8) = 0.20009491565312789468D0
      DISC%DynRup%filter(12,9) = 0.25234135834989557745D-1
      DISC%DynRup%filter(12,10) = -0.37428258638082778188D-1
      DISC%DynRup%filter(12,11) = 0.26033788677614015875D0
      DISC%DynRup%filter(12,12) = 0.42753068173757869386D0
      DISC%DynRup%filter(12,13) = 0.63013182859958368475D-1
      DISC%DynRup%filter(12,14) = -0.14076441672720108742D0
      DISC%DynRup%filter(12,15) = -0.16239363945634718683D-1
      DISC%DynRup%filter(12,16) = 0.23125258610935440011D0
      DISC%DynRup%filter(13,1) = -0.42063743416009481622D-2
      DISC%DynRup%filter(13,2) = 0.73814134524469705475D-2
      DISC%DynRup%filter(13,3) = 0.31274505217622546784D-1
      DISC%DynRup%filter(13,4) = 0.28074354801372348083D-1
      DISC%DynRup%filter(13,5) = -0.75673625640753048507D-1
      DISC%DynRup%filter(13,6) = -0.15995644106858496663D0
      DISC%DynRup%filter(13,7) = -0.36665684156398591098D-1
      DISC%DynRup%filter(13,8) = 0.90898183852533976783D-1
      DISC%DynRup%filter(13,9) = 0.34722735794426513092D0
      DISC%DynRup%filter(13,10) = -0.24383517315009205045D-1
      DISC%DynRup%filter(13,11) = -0.21135874558236743419D0
      DISC%DynRup%filter(13,12) = 0.94614730015497685021D-1
      DISC%DynRup%filter(13,13) = 0.69753698315067017878D0
      DISC%DynRup%filter(13,14) = 0.33177066808066622373D0
      DISC%DynRup%filter(13,15) = -0.15560396125650684607D0
      DISC%DynRup%filter(13,16) = 0.39070152846145977810D-1
      DISC%DynRup%filter(14,1) = 0.39372529478313001736D-2
      DISC%DynRup%filter(14,2) = 0.16235839953299873504D-1
      DISC%DynRup%filter(14,3) = 0.25668955435047758305D-1
      DISC%DynRup%filter(14,4) = 0.16681850793661985247D-1
      DISC%DynRup%filter(14,5) = -0.85320917623590375310D-1
      DISC%DynRup%filter(14,6) = -0.62597481695225032235D-1
      DISC%DynRup%filter(14,7) = -0.13921656933196206479D-1
      DISC%DynRup%filter(14,8) = -0.19557510761191015326D-1
      DISC%DynRup%filter(14,9) = -0.13006191300013112603D-1
      DISC%DynRup%filter(14,10) = 0.20283187677482893780D0
      DISC%DynRup%filter(14,11) = 0.12901309527067642041D0
      DISC%DynRup%filter(14,12) = -0.11273895568310606814D0
      DISC%DynRup%filter(14,13) = 0.17696678953426330769D0
      DISC%DynRup%filter(14,14) = 0.50561218183789380393D0
      DISC%DynRup%filter(14,15) = 0.31319417220581804065D0
      DISC%DynRup%filter(14,16) = -0.82999300756999617497D-1
      DISC%DynRup%filter(15,1) = 0.16681850793661985256D-1
      DISC%DynRup%filter(15,2) = 0.25668955435047758308D-1
      DISC%DynRup%filter(15,3) = 0.16235839953299873517D-1
      DISC%DynRup%filter(15,4) = 0.39372529478313001736D-2
      DISC%DynRup%filter(15,5) = -0.19557510761191015330D-1
      DISC%DynRup%filter(15,6) = -0.13921656933196206485D-1
      DISC%DynRup%filter(15,7) = -0.62597481695225032243D-1
      DISC%DynRup%filter(15,8) = -0.85320917623590375330D-1
      DISC%DynRup%filter(15,9) = -0.11273895568310606812D0
      DISC%DynRup%filter(15,10) = 0.12901309527067642041D0
      DISC%DynRup%filter(15,11) = 0.20283187677482893774D0
      DISC%DynRup%filter(15,12) = -0.13006191300013112643D-1
      DISC%DynRup%filter(15,13) = -0.82999300756999617490D-1
      DISC%DynRup%filter(15,14) = 0.31319417220581804065D0
      DISC%DynRup%filter(15,15) = 0.50561218183789380397D0
      DISC%DynRup%filter(15,16) = 0.17696678953426330770D0
      DISC%DynRup%filter(16,1) = 0.28074354801372348087D-1
      DISC%DynRup%filter(16,2) = 0.31274505217622546782D-1
      DISC%DynRup%filter(16,3) = 0.73814134524469705536D-2
      DISC%DynRup%filter(16,4) = -0.42063743416009481579D-2
      DISC%DynRup%filter(16,5) = 0.90898183852533976769D-1
      DISC%DynRup%filter(16,6) = -0.36665684156398591090D-1
      DISC%DynRup%filter(16,7) = -0.15995644106858496669D0
      DISC%DynRup%filter(16,8) = -0.75673625640753048539D-1
      DISC%DynRup%filter(16,9) = 0.94614730015497684996D-1
      DISC%DynRup%filter(16,10) = -0.21135874558236743427D0
      DISC%DynRup%filter(16,11) = -0.24383517315009205072D-1
      DISC%DynRup%filter(16,12) = 0.34722735794426513094D0
      DISC%DynRup%filter(16,13) = 0.39070152846145977807D-1
      DISC%DynRup%filter(16,14) = -0.15560396125650684608D0
      DISC%DynRup%filter(16,15) = 0.33177066808066622375D0
      DISC%DynRup%filter(16,16) = 0.69753698315067017878D0
    ELSEIF (DISC%Galerkin%nBndGP.EQ.25) THEN
      DISC%DynRup%filter(1,1) = 0.126649237400676823D0
      DISC%DynRup%filter(1,2) = 0.237393198874224859D0
      DISC%DynRup%filter(1,3) = 0.250296762019001029D0
      DISC%DynRup%filter(1,4) = 0.184037732673873800D0
      DISC%DynRup%filter(1,5) = 0.822015547979396549D-1
      DISC%DynRup%filter(1,6) = 0.192765467851087169D0
      DISC%DynRup%filter(1,7) = 0.265221272830132526D0
      DISC%DynRup%filter(1,8) = 0.107734191062446513D0
      DISC%DynRup%filter(1,9) = -0.753056662980977065D-1
      DISC%DynRup%filter(1,10) = -0.909251707016350552D-1
      DISC%DynRup%filter(1,11) = -0.209072054007676221D-1
      DISC%DynRup%filter(1,12) = -0.102500193808265555D0
      DISC%DynRup%filter(1,13) = -0.194309374123471851D0
      DISC%DynRup%filter(1,14) = -0.192400677049793745D0
      DISC%DynRup%filter(1,15) = -0.959973575161842396D-1
      DISC%DynRup%filter(1,16) = -0.397654123591174302D-1
      DISC%DynRup%filter(1,17) = -0.784581135967018714D-2
      DISC%DynRup%filter(1,18) = 0.135814554816169469D0
      DISC%DynRup%filter(1,19) = 0.254190180307654945D0
      DISC%DynRup%filter(1,20) = 0.177862684481334965D0
      DISC%DynRup%filter(1,21) = 0.220874455520700293D-1
      DISC%DynRup%filter(1,22) = 0.166056465974928485D-1
      DISC%DynRup%filter(1,23) = -0.441325793207340208D-1
      DISC%DynRup%filter(1,24) = -0.107127102450584252D0
      DISC%DynRup%filter(1,25) = -0.816433803739799319D-1
      DISC%DynRup%filter(2,1) = 0.117512456212712182D0
      DISC%DynRup%filter(2,2) = 0.226589272103784295D0
      DISC%DynRup%filter(2,3) = 0.250491568897271366D0
      DISC%DynRup%filter(2,4) = 0.194884323077293392D0
      DISC%DynRup%filter(2,5) = 0.911008655413345914D-1
      DISC%DynRup%filter(2,6) = 0.131287683365774621D0
      DISC%DynRup%filter(2,7) = 0.196842932248576463D0
      DISC%DynRup%filter(2,8) = 0.114157311352126090D0
      DISC%DynRup%filter(2,9) = -0.552063347576991874D-2
      DISC%DynRup%filter(2,10) = -0.372772006135837011D-1
      DISC%DynRup%filter(2,11) = -0.507388138016544657D-1
      DISC%DynRup%filter(2,12) = -0.118028921472924805D0
      DISC%DynRup%filter(2,13) = -0.170458053715789426D0
      DISC%DynRup%filter(2,14) = -0.171648397585313944D0
      DISC%DynRup%filter(2,15) = -0.952406210524022190D-1
      DISC%DynRup%filter(2,16) = -0.388377047893724570D-2
      DISC%DynRup%filter(2,17) = 0.470998822092249109D-1
      DISC%DynRup%filter(2,18) = 0.149066251482438378D0
      DISC%DynRup%filter(2,19) = 0.202146672341430467D0
      DISC%DynRup%filter(2,20) = 0.125827164261824109D0
      DISC%DynRup%filter(2,21) = 0.821999185053042154D-2
      DISC%DynRup%filter(2,22) = -0.949537322978761915D-2
      DISC%DynRup%filter(2,23) = -0.562299378230998528D-1
      DISC%DynRup%filter(2,24) = -0.836754621142230676D-1
      DISC%DynRup%filter(2,25) = -0.530291913431370845D-1
      DISC%DynRup%filter(3,1) = 0.104241852520110284D0
      DISC%DynRup%filter(3,2) = 0.210748441185749730D0
      DISC%DynRup%filter(3,3) = 0.250597899157088744D0
      DISC%DynRup%filter(3,4) = 0.210748440589761582D0
      DISC%DynRup%filter(3,5) = 0.104241851877704683D0
      DISC%DynRup%filter(3,6) = 0.448683880269845695D-1
      DISC%DynRup%filter(3,7) = 0.960450550243799467D-1
      DISC%DynRup%filter(3,8) = 0.117663211698128198D0
      DISC%DynRup%filter(3,9) = 0.960450513534612527D-1
      DISC%DynRup%filter(3,10) = 0.448683860426770040D-1
      DISC%DynRup%filter(3,11) = -0.809246143634732157D-1
      DISC%DynRup%filter(3,12) = -0.143413087355267144D0
      DISC%DynRup%filter(3,13) = -0.157439406278221533D0
      DISC%DynRup%filter(3,14) = -0.143413086370521609D0
      DISC%DynRup%filter(3,15) = -0.809246134848602000D-1
      DISC%DynRup%filter(3,16) = 0.565630990215231569D-1
      DISC%DynRup%filter(3,17) = 0.125415317706323054D0
      DISC%DynRup%filter(3,18) = 0.156299359165836316D0
      DISC%DynRup%filter(3,19) = 0.125415322241949317D0
      DISC%DynRup%filter(3,20) = 0.565631034497972585D-1
      DISC%DynRup%filter(3,21) = -0.183800298803868216D-1
      DISC%DynRup%filter(3,22) = -0.473084648881170189D-1
      DISC%DynRup%filter(3,23) = -0.628329788611468876D-1
      DISC%DynRup%filter(3,24) = -0.473084678308004145D-1
      DISC%DynRup%filter(3,25) = -0.183800320790184941D-1
      DISC%DynRup%filter(4,1) = 0.911008660099973944D-1
      DISC%DynRup%filter(4,2) = 0.194884323077293392D0
      DISC%DynRup%filter(4,3) = 0.250491568188891234D0
      DISC%DynRup%filter(4,4) = 0.226589271139397641D0
      DISC%DynRup%filter(4,5) = 0.117512455581199049D0
      DISC%DynRup%filter(4,6) = -0.372771992277989850D-1
      DISC%DynRup%filter(4,7) = -0.552063104683602011D-2
      DISC%DynRup%filter(4,8) = 0.114157312633576757D0
      DISC%DynRup%filter(4,9) = 0.196842928619177043D0
      DISC%DynRup%filter(4,10) = 0.131287681890711327D0
      DISC%DynRup%filter(4,11) = -0.952406216797501565D-1
      DISC%DynRup%filter(4,12) = -0.171648398134837676D0
      DISC%DynRup%filter(4,13) = -0.170458054558396438D0
      DISC%DynRup%filter(4,14) = -0.118028920448862026D0
      DISC%DynRup%filter(4,15) = -0.507388129240293334D-1
      DISC%DynRup%filter(4,16) = 0.125827160288683765D0
      DISC%DynRup%filter(4,17) = 0.202146668787253042D0
      DISC%DynRup%filter(4,18) = 0.149066251804333055D0
      DISC%DynRup%filter(4,19) = 0.470998860903535585D-1
      DISC%DynRup%filter(4,20) = -0.388376669244747076D-2
      DISC%DynRup%filter(4,21) = -0.530291896317336489D-1
      DISC%DynRup%filter(4,22) = -0.836754595075475116D-1
      DISC%DynRup%filter(4,23) = -0.562299377216796750D-1
      DISC%DynRup%filter(4,24) = -0.949537568893815860D-2
      DISC%DynRup%filter(4,25) = 0.821999003711472616D-2
      DISC%DynRup%filter(5,1) = 0.822015547979396688D-1
      DISC%DynRup%filter(5,2) = 0.184037731727102999D0
      DISC%DynRup%filter(5,3) = 0.250296760476510782D0
      DISC%DynRup%filter(5,4) = 0.237393197598471384D0
      DISC%DynRup%filter(5,5) = 0.126649236820880584D0
      DISC%DynRup%filter(5,6) = -0.909251708509276196D-1
      DISC%DynRup%filter(5,7) = -0.753056657239956073D-1
      DISC%DynRup%filter(5,8) = 0.107734193782441326D0
      DISC%DynRup%filter(5,9) = 0.265221270845967805D0
      DISC%DynRup%filter(5,10) = 0.192765467846866073D0
      DISC%DynRup%filter(5,11) = -0.959973579469476923D-1
      DISC%DynRup%filter(5,12) = -0.192400676168223872D0
      DISC%DynRup%filter(5,13) = -0.194309374185924005D0
      DISC%DynRup%filter(5,14) = -0.102500192899394174D0
      DISC%DynRup%filter(5,15) = -0.209072052288645713D-1
      DISC%DynRup%filter(5,16) = 0.177862681278273899D0
      DISC%DynRup%filter(5,17) = 0.254190178912143627D0
      DISC%DynRup%filter(5,18) = 0.135814555279988591D0
      DISC%DynRup%filter(5,19) = -0.784581041632487902D-2
      DISC%DynRup%filter(5,20) = -0.397654107522127792D-1
      DISC%DynRup%filter(5,21) = -0.816433799164070634D-1
      DISC%DynRup%filter(5,22) = -0.107127100743489925D0
      DISC%DynRup%filter(5,23) = -0.441325792317234025D-1
      DISC%DynRup%filter(5,24) = 0.166056454468267271D-1
      DISC%DynRup%filter(5,25) = 0.220874454753936383D-1
      DISC%DynRup%filter(6,1) = 0.410729335866528780D-1
      DISC%DynRup%filter(6,2) = 0.565112404057486473D-1
      DISC%DynRup%filter(6,3) = 0.229551462405463563D-1
      DISC%DynRup%filter(6,4) = -0.160455323242020172D-1
      DISC%DynRup%filter(6,5) = -0.193736126358492571D-1
      DISC%DynRup%filter(6,6) = 0.318330067485885715D0
      DISC%DynRup%filter(6,7) = 0.416325595974134488D0
      DISC%DynRup%filter(6,8) = 0.161868426782236258D0
      DISC%DynRup%filter(6,9) = -0.830906258294618200D-1
      DISC%DynRup%filter(6,10) = -0.981783510346773275D-1
      DISC%DynRup%filter(6,11) = 0.268672706961027352D0
      DISC%DynRup%filter(6,12) = 0.229690751699140772D0
      DISC%DynRup%filter(6,13) = -0.331568189947115191D-1
      DISC%DynRup%filter(6,14) = -0.570872422524531450D-1
      DISC%DynRup%filter(6,15) = 0.239015517736891996D-1
      DISC%DynRup%filter(6,16) = -0.853498027248205393D-1
      DISC%DynRup%filter(6,17) = -0.228497707117464310D0
      DISC%DynRup%filter(6,18) = -0.198609745651033398D0
      DISC%DynRup%filter(6,19) = 0.306148318075481658D-1
      DISC%DynRup%filter(6,20) = 0.111018982826401144D0
      DISC%DynRup%filter(6,21) = 0.155020881098146181D-1
      DISC%DynRup%filter(6,22) = 0.795529909043319328D-1
      DISC%DynRup%filter(6,23) = 0.109659629591701974D0
      DISC%DynRup%filter(6,24) = 0.209493443660965480D-2
      DISC%DynRup%filter(6,25) = -0.683824402961213385D-1
      DISC%DynRup%filter(7,1) = 0.279737360410679432D-1
      DISC%DynRup%filter(7,2) = 0.419417273000458404D-1
      DISC%DynRup%filter(7,3) = 0.243237335796138605D-1
      DISC%DynRup%filter(7,4) = -0.117629217999131960D-2
      DISC%DynRup%filter(7,5) = -0.794272945333899041D-2
      DISC%DynRup%filter(7,6) = 0.206086120464411188D0
      DISC%DynRup%filter(7,7) = 0.320192878019005422D0
      DISC%DynRup%filter(7,8) = 0.207148354129096518D0
      DISC%DynRup%filter(7,9) = 0.229586075943931864D-1
      DISC%DynRup%filter(7,10) = -0.411308486625149156D-1
      DISC%DynRup%filter(7,11) = 0.113699653258892044D0
      DISC%DynRup%filter(7,12) = 0.193048804850437933D0
      DISC%DynRup%filter(7,13) = 0.136763336126311030D0
      DISC%DynRup%filter(7,14) = 0.167680191209812861D-1
      DISC%DynRup%filter(7,15) = -0.282588632982551831D-1
      DISC%DynRup%filter(7,16) = -0.113109081653725782D0
      DISC%DynRup%filter(7,17) = -0.155079256835467277D0
      DISC%DynRup%filter(7,18) = -0.971962260149470125D-1
      DISC%DynRup%filter(7,19) = -0.205935780143297054D-1
      DISC%DynRup%filter(7,20) = 0.151547071305784906D-1
      DISC%DynRup%filter(7,21) = 0.393796770124243689D-1
      DISC%DynRup%filter(7,22) = 0.651358733509134480D-1
      DISC%DynRup%filter(7,23) = 0.331242005405212384D-1
      DISC%DynRup%filter(7,24) = -0.249568306305498360D-3
      DISC%DynRup%filter(7,25) = 0.103701663237614036D-2
      DISC%DynRup%filter(8,1) = 0.956019896221255334D-2
      DISC%DynRup%filter(8,2) = 0.204645159208947484D-1
      DISC%DynRup%filter(8,3) = 0.250707418370476227D-1
      DISC%DynRup%filter(8,4) = 0.204645161506151821D-1
      DISC%DynRup%filter(8,5) = 0.956019920358153197D-2
      DISC%DynRup%filter(8,6) = 0.674138352247070360D-1
      DISC%DynRup%filter(8,7) = 0.174282084248812780D0
      DISC%DynRup%filter(8,8) = 0.231863274121814872D0
      DISC%DynRup%filter(8,9) = 0.174282085494603944D0
      DISC%DynRup%filter(8,10) = 0.674138358368934903D-1
      DISC%DynRup%filter(8,11) = -0.138089217177095845D-1
      DISC%DynRup%filter(8,12) = 0.115064392209549851D0
      DISC%DynRup%filter(8,13) = 0.229510012036522743D0
      DISC%DynRup%filter(8,14) = 0.115064391624459708D0
      DISC%DynRup%filter(8,15) = -0.138089221908465429D-1
      DISC%DynRup%filter(8,16) = -0.827156032280315368D-1
      DISC%DynRup%filter(8,17) = -0.817750202395710474D-1
      DISC%DynRup%filter(8,18) = -0.418421899820246163D-1
      DISC%DynRup%filter(8,19) = -0.817750223458549863D-1
      DISC%DynRup%filter(8,20) = -0.827156050371562340D-1
      DISC%DynRup%filter(8,21) = 0.456702789708136811D-1
      DISC%DynRup%filter(8,22) = 0.278686970414667654D-1
      DISC%DynRup%filter(8,23) = -0.865075201928662031D-2
      DISC%DynRup%filter(8,24) = 0.278686983746651712D-1
      DISC%DynRup%filter(8,25) = 0.456702802384039430D-1
      DISC%DynRup%filter(9,1) = -0.794272951389137297D-2
      DISC%DynRup%filter(9,2) = -0.117629269752925333D-2
      DISC%DynRup%filter(9,3) = 0.243237326499413410D-1
      DISC%DynRup%filter(9,4) = 0.419417265267222755D-1
      DISC%DynRup%filter(9,5) = 0.279737358317917462D-1
      DISC%DynRup%filter(9,6) = -0.411308478021554444D-1
      DISC%DynRup%filter(9,7) = 0.229586075943931552D-1
      DISC%DynRup%filter(9,8) = 0.207148355609820034D0
      DISC%DynRup%filter(9,9) = 0.320192876418365469D0
      DISC%DynRup%filter(9,10) = 0.206086119340495144D0
      DISC%DynRup%filter(9,11) = -0.282588630967663991D-1
      DISC%DynRup%filter(9,12) = 0.167680195972687643D-1
      DISC%DynRup%filter(9,13) = 0.136763335404885777D0
      DISC%DynRup%filter(9,14) = 0.193048805236201410D0
      DISC%DynRup%filter(9,15) = 0.113699653159450978D0
      DISC%DynRup%filter(9,16) = 0.151547068189641556D-1
      DISC%DynRup%filter(9,17) = -0.205935783165272467D-1
      DISC%DynRup%filter(9,18) = -0.971962254194515402D-1
      DISC%DynRup%filter(9,19) = -0.155079255914575997D0
      DISC%DynRup%filter(9,20) = -0.113109081831030092D0
      DISC%DynRup%filter(9,21) = 0.103701694336937717D-2
      DISC%DynRup%filter(9,22) = -0.249568283606738605D-3
      DISC%DynRup%filter(9,23) = 0.331242007781653719D-1
      DISC%DynRup%filter(9,24) = 0.651358734885949114D-1
      DISC%DynRup%filter(9,25) = 0.393796764387531473D-1
      DISC%DynRup%filter(10,1) = -0.193736126040391816D-1
      DISC%DynRup%filter(10,2) = -0.160455329206967597D-1
      DISC%DynRup%filter(10,3) = 0.229551452253534335D-1
      DISC%DynRup%filter(10,4) = 0.565112397708249919D-1
      DISC%DynRup%filter(10,5) = 0.410729335857534794D-1
      DISC%DynRup%filter(10,6) = -0.981783510346773275D-1
      DISC%DynRup%filter(10,7) = -0.830906275675200173D-1
      DISC%DynRup%filter(10,8) = 0.161868428252166824D0
      DISC%DynRup%filter(10,9) = 0.416325593703651575D0
      DISC%DynRup%filter(10,10) = 0.318330066793356181D0
      DISC%DynRup%filter(10,11) = 0.239015527468505844D-1
      DISC%DynRup%filter(10,12) = -0.570872416366516242D-1
      DISC%DynRup%filter(10,13) = -0.331568218111994231D-1
      DISC%DynRup%filter(10,14) = 0.229690751322169984D0
      DISC%DynRup%filter(10,15) = 0.268672706849347243D0
      DISC%DynRup%filter(10,16) = 0.111018983539645413D0
      DISC%DynRup%filter(10,17) = 0.306148334861245514D-1
      DISC%DynRup%filter(10,18) = -0.198609743314420523D0
      DISC%DynRup%filter(10,19) = -0.228497705364670539D0
      DISC%DynRup%filter(10,20) = -0.853498028179669882D-1
      DISC%DynRup%filter(10,21) = -0.683824411383817804D-1
      DISC%DynRup%filter(10,22) = 0.209493339616607759D-2
      DISC%DynRup%filter(10,23) = 0.109659630370997824D0
      DISC%DynRup%filter(10,24) = 0.795529913202112221D-1
      DISC%DynRup%filter(10,25) = 0.155020865568915294D-1
      DISC%DynRup%filter(11,1) = -0.224914038019239813D-2
      DISC%DynRup%filter(11,2) = -0.110266922918833866D-1
      DISC%DynRup%filter(11,3) = -0.209032746068422634D-1
      DISC%DynRup%filter(11,4) = -0.206979420736091192D-1
      DISC%DynRup%filter(11,5) = -0.103271350719277167D-1
      DISC%DynRup%filter(11,6) = 0.135649326428708611D0
      DISC%DynRup%filter(11,7) = 0.115967848370969892D0
      DISC%DynRup%filter(11,8) = -0.167404438608757186D-1
      DISC%DynRup%filter(11,9) = -0.288225993379228780D-1
      DISC%DynRup%filter(11,10) = 0.120675805420787929D-1
      DISC%DynRup%filter(11,11) = 0.447079842257498195D0
      DISC%DynRup%filter(11,12) = 0.322881930886912305D0
      DISC%DynRup%filter(11,13) = -0.114766623540389578D0
      DISC%DynRup%filter(11,14) = -0.694634658962941520D-1
      DISC%DynRup%filter(11,15) = 0.833294461220658172D-1
      DISC%DynRup%filter(11,16) = 0.301026109944167009D0
      DISC%DynRup%filter(11,17) = 0.131382303513677384D0
      DISC%DynRup%filter(11,18) = -0.138007787659341158D0
      DISC%DynRup%filter(11,19) = -0.433428086041094951D-1
      DISC%DynRup%filter(11,20) = 0.330022168844885042D-1
      DISC%DynRup%filter(11,21) = -0.817567013272564452D-1
      DISC%DynRup%filter(11,22) = -0.108426639158082774D0
      DISC%DynRup%filter(11,23) = 0.566779571413645877D-1
      DISC%DynRup%filter(11,24) = 0.766918619987364131D-1
      DISC%DynRup%filter(11,25) = -0.492251785992543350D-1
      DISC%DynRup%filter(12,1) = -0.545834389870728052D-2
      DISC%DynRup%filter(12,2) = -0.126972308434024385D-1
      DISC%DynRup%filter(12,3) = -0.183374145828355400D-1
      DISC%DynRup%filter(12,4) = -0.184654685294081827D-1
      DISC%DynRup%filter(12,5) = -0.102457275235443793D-1
      DISC%DynRup%filter(12,6) = 0.574054638940168019D-1
      DISC%DynRup%filter(12,7) = 0.974678099228973860D-1
      DISC%DynRup%filter(12,8) = 0.690500154611690953D-1
      DISC%DynRup%filter(12,9) = 0.846595319849919280D-2
      DISC%DynRup%filter(12,10) = -0.142675295558888429D-1
      DISC%DynRup%filter(12,11) = 0.159830395048467605D0
      DISC%DynRup%filter(12,12) = 0.298110715235178014D0
      DISC%DynRup%filter(12,13) = 0.217442213896502506D0
      DISC%DynRup%filter(12,14) = 0.280630503918998277D-1
      DISC%DynRup%filter(12,15) = -0.343852414259764802D-1
      DISC%DynRup%filter(12,16) = 0.650358030676682353D-1
      DISC%DynRup%filter(12,17) = 0.183257864394674730D0
      DISC%DynRup%filter(12,18) = 0.100258504320620959D0
      DISC%DynRup%filter(12,19) = -0.430369290660866446D-1
      DISC%DynRup%filter(12,20) = -0.214552060282828146D-1
      DISC%DynRup%filter(12,21) = -0.536724750915595350D-1
      DISC%DynRup%filter(12,22) = -0.203823642488661534D-1
      DISC%DynRup%filter(12,23) = -0.378850818109888535D-1
      DISC%DynRup%filter(12,24) = -0.320621634596619243D-1
      DISC%DynRup%filter(12,25) = 0.379633832353791223D-1
      DISC%DynRup%filter(13,1) = -0.870565032588176656D-2
      DISC%DynRup%filter(13,2) = -0.154279904146435148D-1
      DISC%DynRup%filter(13,3) = -0.169369037754391501D-1
      DISC%DynRup%filter(13,4) = -0.154279904909070467D-1
      DISC%DynRup%filter(13,5) = -0.870565032867980960D-2
      DISC%DynRup%filter(13,6) = -0.697194317073457012D-2
      DISC%DynRup%filter(13,7) = 0.580945022185465135D-1
      DISC%DynRup%filter(13,8) = 0.115876595331782170D0
      DISC%DynRup%filter(13,9) = 0.580945019120985978D-1
      DISC%DynRup%filter(13,10) = -0.697194376296248700D-2
      DISC%DynRup%filter(13,11) = -0.477972042009845582D-1
      DISC%DynRup%filter(13,12) = 0.182942714743913409D0
      DISC%DynRup%filter(13,13) = 0.398770113117102420D0
      DISC%DynRup%filter(13,14) = 0.182942714965008052D0
      DISC%DynRup%filter(13,15) = -0.477972040325519223D-1
      DISC%DynRup%filter(13,16) = -0.574765220732079835D-1
      DISC%DynRup%filter(13,17) = 0.843514360953923070D-1
      DISC%DynRup%filter(13,18) = 0.230310211005658239D0
      DISC%DynRup%filter(13,19) = 0.843514362916361060D-1
      DISC%DynRup%filter(13,20) = -0.574765221685975894D-1
      DISC%DynRup%filter(13,21) = 0.236048405921999022D-1
      DISC%DynRup%filter(13,22) = -0.318742139362529764D-1
      DISC%DynRup%filter(13,23) = -0.894999555825313070D-1
      DISC%DynRup%filter(13,24) = -0.318742143099685626D-1
      DISC%DynRup%filter(13,25) = 0.236048407259937219D-1
      DISC%DynRup%filter(14,1) = -0.102457275704897731D-1
      DISC%DynRup%filter(14,2) = -0.184654684702919118D-1
      DISC%DynRup%filter(14,3) = -0.183374144569217433D-1
      DISC%DynRup%filter(14,4) = -0.126972307332365499D-1
      DISC%DynRup%filter(14,5) = -0.545834385030802890D-2
      DISC%DynRup%filter(14,6) = -0.142675297097930500D-1
      DISC%DynRup%filter(14,7) = 0.846595295802789541D-2
      DISC%DynRup%filter(14,8) = 0.690500151100571502D-1
      DISC%DynRup%filter(14,9) = 0.974678101176642969D-1
      DISC%DynRup%filter(14,10) = 0.574054637998023740D-1
      DISC%DynRup%filter(14,11) = -0.343852415808582404D-1
      DISC%DynRup%filter(14,12) = 0.280630503918998277D-1
      DISC%DynRup%filter(14,13) = 0.217442214159291353D0
      DISC%DynRup%filter(14,14) = 0.298110715176488239D0
      DISC%DynRup%filter(14,15) = 0.159830394732174030D0
      DISC%DynRup%filter(14,16) = -0.214552057294015268D-1
      DISC%DynRup%filter(14,17) = -0.430369287877111345D-1
      DISC%DynRup%filter(14,18) = 0.100258504280521812D0
      DISC%DynRup%filter(14,19) = 0.183257864441884744D0
      DISC%DynRup%filter(14,20) = 0.650358025559096276D-1
      DISC%DynRup%filter(14,21) = 0.379633833113437250D-1
      DISC%DynRup%filter(14,22) = -0.320621635911863753D-1
      DISC%DynRup%filter(14,23) = -0.378850816484505859D-1
      DISC%DynRup%filter(14,24) = -0.203823639962634855D-1
      DISC%DynRup%filter(14,25) = -0.536724748911700522D-1
      DISC%DynRup%filter(15,1) = -0.103271350255873442D-1
      DISC%DynRup%filter(15,2) = -0.206979419372722191D-1
      DISC%DynRup%filter(15,3) = -0.209032743798916697D-1
      DISC%DynRup%filter(15,4) = -0.110266921011555902D-1
      DISC%DynRup%filter(15,5) = -0.224914036169953984D-2
      DISC%DynRup%filter(15,6) = 0.120675800507423665D-1
      DISC%DynRup%filter(15,7) = -0.288225995434311759D-1
      DISC%DynRup%filter(15,8) = -0.167404444344558048D-1
      DISC%DynRup%filter(15,9) = 0.115967848269545079D0
      DISC%DynRup%filter(15,10) = 0.135649326372322798D0
      DISC%DynRup%filter(15,11) = 0.833294461220658172D-1
      DISC%DynRup%filter(15,12) = -0.694634655834092257D-1
      DISC%DynRup%filter(15,13) = -0.114766623135963322D0
      DISC%DynRup%filter(15,14) = 0.322881930247950699D0
      DISC%DynRup%filter(15,15) = 0.447079840847797616D0
      DISC%DynRup%filter(15,16) = 0.330022172126200192D-1
      DISC%DynRup%filter(15,17) = -0.433428089153108792D-1
      DISC%DynRup%filter(15,18) = -0.138007787081388583D0
      DISC%DynRup%filter(15,19) = 0.131382304696690155D0
      DISC%DynRup%filter(15,20) = 0.301026109534524133D0
      DISC%DynRup%filter(15,21) = -0.492251778710173296D-1
      DISC%DynRup%filter(15,22) = 0.766918610898476572D-1
      DISC%DynRup%filter(15,23) = 0.566779574285141002D-1
      DISC%DynRup%filter(15,24) = -0.108426637880082957D0
      DISC%DynRup%filter(15,25) = -0.817567021351359030D-1
      DISC%DynRup%filter(16,1) = -0.374591696706182186D-2
      DISC%DynRup%filter(16,2) = -0.739078521134447987D-3
      DISC%DynRup%filter(16,3) = 0.127937824806232496D-1
      DISC%DynRup%filter(16,4) = 0.239448113757112779D-1
      DISC%DynRup%filter(16,5) = 0.167547322177997032D-1
      DISC%DynRup%filter(16,6) = -0.377336391122209333D-1
      DISC%DynRup%filter(16,7) = -0.101020151339744485D0
      DISC%DynRup%filter(16,8) = -0.878065112039025375D-1
      DISC%DynRup%filter(16,9) = 0.135349943079551994D-1
      DISC%DynRup%filter(16,10) = 0.490821317185344830D-1
      DISC%DynRup%filter(16,11) = 0.263594440134837515D0
      DISC%DynRup%filter(16,12) = 0.115045319729165776D0
      DISC%DynRup%filter(16,13) = -0.120846943561141687D0
      DISC%DynRup%filter(16,14) = -0.379532639956148068D-1
      DISC%DynRup%filter(16,15) = 0.288984931273315836D-1
      DISC%DynRup%filter(16,16) = 0.527242469825473292D0
      DISC%DynRup%filter(16,17) = 0.277985272993021559D0
      DISC%DynRup%filter(16,18) = -0.351171284614251247D-1
      DISC%DynRup%filter(16,19) = 0.417692271172955126D-1
      DISC%DynRup%filter(16,20) = -0.253831418766758792D-1
      DISC%DynRup%filter(16,21) = 0.186299453489571887D0
      DISC%DynRup%filter(16,22) = -0.755807621439533450D-1
      DISC%DynRup%filter(16,23) = -0.899736571387409756D-1
      DISC%DynRup%filter(16,24) = 0.890579670575030874D-1
      DISC%DynRup%filter(16,25) = -0.301029060822143452D-1
      DISC%DynRup%filter(17,1) = -0.365852602792934552D-3
      DISC%DynRup%filter(17,2) = 0.443682681358834475D-2
      DISC%DynRup%filter(17,3) = 0.140420973758204129D-1
      DISC%DynRup%filter(17,4) = 0.190422930649534194D-1
      DISC%DynRup%filter(17,5) = 0.118529664168880444D-1
      DISC%DynRup%filter(17,6) = -0.500061762584364730D-1
      DISC%DynRup%filter(17,7) = -0.685614321214681466D-1
      DISC%DynRup%filter(17,8) = -0.429710110079540852D-1
      DISC%DynRup%filter(17,9) = -0.910453951546024570D-2
      DISC%DynRup%filter(17,10) = 0.669998302715054533D-2
      DISC%DynRup%filter(17,11) = 0.569488006627246121D-1
      DISC%DynRup%filter(17,12) = 0.160470313205551607D0
      DISC%DynRup%filter(17,13) = 0.877916682582307317D-1
      DISC%DynRup%filter(17,14) = -0.376854192030500884D-1
      DISC%DynRup%filter(17,15) = -0.187873170059287380D-1
      DISC%DynRup%filter(17,16) = 0.137606015486930994D0
      DISC%DynRup%filter(17,17) = 0.438845946797728370D0
      DISC%DynRup%filter(17,18) = 0.246733343578125458D0
      DISC%DynRup%filter(17,19) = -0.573648649435567898D-1
      DISC%DynRup%filter(17,20) = 0.206762641254219295D-1
      DISC%DynRup%filter(17,21) = -0.374133766985503616D-1
      DISC%DynRup%filter(17,22) = 0.148031626103087149D0
      DISC%DynRup%filter(17,23) = 0.326880279849608851D-1
      DISC%DynRup%filter(17,24) = -0.107690934053554080D0
      DISC%DynRup%filter(17,25) = 0.440847526207119311D-1
      DISC%DynRup%filter(18,1) = 0.532826553350770633D-2
      DISC%DynRup%filter(18,2) = 0.118141707822741644D-1
      DISC%DynRup%filter(18,3) = 0.147234590671309826D-1
      DISC%DynRup%filter(18,4) = 0.118141708077857607D-1
      DISC%DynRup%filter(18,5) = 0.532826555170422007D-2
      DISC%DynRup%filter(18,6) = -0.365690446133524263D-1
      DISC%DynRup%filter(18,7) = -0.361532071788102094D-1
      DISC%DynRup%filter(18,8) = -0.184986730525766238D-1
      DISC%DynRup%filter(18,9) = -0.361532069573091142D-1
      DISC%DynRup%filter(18,10) = -0.365690441831233032D-1
      DISC%DynRup%filter(18,11) = -0.503294938014961155D-1
      DISC%DynRup%filter(18,12) = 0.738625945148532609D-1
      DISC%DynRup%filter(18,13) = 0.201671845377789150D0
      DISC%DynRup%filter(18,14) = 0.738625944853113503D-1
      DISC%DynRup%filter(18,15) = -0.503294935907249649D-1
      DISC%DynRup%filter(18,16) = -0.146253371186536221D-1
      DISC%DynRup%filter(18,17) = 0.207586497951167931D0
      DISC%DynRup%filter(18,18) = 0.400574386041718855D0
      DISC%DynRup%filter(18,19) = 0.207586497313406948D0
      DISC%DynRup%filter(18,20) = -0.146253369517749688D-1
      DISC%DynRup%filter(18,21) = -0.374716030125250982D-1
      DISC%DynRup%filter(18,22) = 0.275017279016362623D-1
      DISC%DynRup%filter(18,23) = 0.996398462830784537D-1
      DISC%DynRup%filter(18,24) = 0.275017275510161356D-1
      DISC%DynRup%filter(18,25) = -0.374716022647492847D-1
      DISC%DynRup%filter(19,1) = 0.118529664819611684D-1
      DISC%DynRup%filter(19,2) = 0.190422933997582873D-1
      DISC%DynRup%filter(19,3) = 0.140420978836507740D-1
      DISC%DynRup%filter(19,4) = 0.443682717919213678D-2
      DISC%DynRup%filter(19,5) = -0.365852558804451178D-3
      DISC%DynRup%filter(19,6) = 0.669998265979812949D-2
      DISC%DynRup%filter(19,7) = -0.910453938185697836D-2
      DISC%DynRup%filter(19,8) = -0.429710121147609050D-1
      DISC%DynRup%filter(19,9) = -0.685614317143368041D-1
      DISC%DynRup%filter(19,10) = -0.500061758748417476D-1
      DISC%DynRup%filter(19,11) = -0.187873168710357974D-1
      DISC%DynRup%filter(19,12) = -0.376854194468104753D-1
      DISC%DynRup%filter(19,13) = 0.877916684624782395D-1
      DISC%DynRup%filter(19,14) = 0.160470313246891177D0
      DISC%DynRup%filter(19,15) = 0.569488011755117721D-1
      DISC%DynRup%filter(19,16) = 0.206762640757735157D-1
      DISC%DynRup%filter(19,17) = -0.573648649435567898D-1
      DISC%DynRup%filter(19,18) = 0.246733342820094931D0
      DISC%DynRup%filter(19,19) = 0.438845947047239837D0
      DISC%DynRup%filter(19,20) = 0.137606016172368090D0
      DISC%DynRup%filter(19,21) = 0.440847527709886944D-1
      DISC%DynRup%filter(19,22) = -0.107690934118601270D0
      DISC%DynRup%filter(19,23) = 0.326880276348220922D-1
      DISC%DynRup%filter(19,24) = 0.148031625525316513D0
      DISC%DynRup%filter(19,25) = -0.374133762530053579D-1
      DISC%DynRup%filter(20,1) = 0.167547325195292716D-1
      DISC%DynRup%filter(20,2) = 0.239448121317968089D-1
      DISC%DynRup%filter(20,3) = 0.127937834822369582D-1
      DISC%DynRup%filter(20,4) = -0.739077800568372972D-3
      DISC%DynRup%filter(20,5) = -0.374591681569079243D-2
      DISC%DynRup%filter(20,6) = 0.490821314032050884D-1
      DISC%DynRup%filter(20,7) = 0.135349945862646426D-1
      DISC%DynRup%filter(20,8) = -0.878065131243736946D-1
      DISC%DynRup%filter(20,9) = -0.101020151498098801D0
      DISC%DynRup%filter(20,10) = -0.377336391534014853D-1
      DISC%DynRup%filter(20,11) = 0.288984928400022253D-1
      DISC%DynRup%filter(20,12) = -0.379532645243219266D-1
      DISC%DynRup%filter(20,13) = -0.120846943761702588D0
      DISC%DynRup%filter(20,14) = 0.115045318823888548D0
      DISC%DynRup%filter(20,15) = 0.263594439776132505D0
      DISC%DynRup%filter(20,16) = -0.253831418766758514D-1
      DISC%DynRup%filter(20,17) = 0.417692272175929219D-1
      DISC%DynRup%filter(20,18) = -0.351171280607301614D-1
      DISC%DynRup%filter(20,19) = 0.277985274377709646D0
      DISC%DynRup%filter(20,20) = 0.527242472177020360D0
      DISC%DynRup%filter(20,21) = -0.301029060199849013D-1
      DISC%DynRup%filter(20,22) = 0.890579674018394385D-1
      DISC%DynRup%filter(20,23) = -0.899736577821243083D-1
      DISC%DynRup%filter(20,24) = -0.755807629947486681D-1
      DISC%DynRup%filter(20,25) = 0.186299453123348308D0
      DISC%DynRup%filter(21,1) = 0.359398100570748396D-2
      DISC%DynRup%filter(21,2) = 0.270200433404546407D-2
      DISC%DynRup%filter(21,3) = -0.718107742536976904D-2
      DISC%DynRup%filter(21,4) = -0.174312946802516813D-1
      DISC%DynRup%filter(21,5) = -0.132846850021471653D-1
      DISC%DynRup%filter(21,6) = 0.118384269957140745D-1
      DISC%DynRup%filter(21,7) = 0.607519633432323097D-1
      DISC%DynRup%filter(21,8) = 0.837433971113602516D-1
      DISC%DynRup%filter(21,9) = 0.159983067687452908D-2
      DISC%DynRup%filter(21,10) = -0.522213866590600997D-1
      DISC%DynRup%filter(21,11) = -0.123661089711676331D0
      DISC%DynRup%filter(21,12) = -0.164000699602106198D0
      DISC%DynRup%filter(21,13) = 0.857282367401541734D-1
      DISC%DynRup%filter(21,14) = 0.116000266648824429D0
      DISC%DynRup%filter(21,15) = -0.744555374416959875D-1
      DISC%DynRup%filter(21,16) = 0.321802349339813776D0
      DISC%DynRup%filter(21,17) = -0.130553615998884903D0
      DISC%DynRup%filter(21,18) = -0.155415028578829939D0
      DISC%DynRup%filter(21,19) = 0.153833318255190699D0
      DISC%DynRup%filter(21,20) = -0.519979296650433312D-1
      DISC%DynRup%filter(21,21) = 0.772527849787711274D0
      DISC%DynRup%filter(21,22) = 0.289352773781101635D0
      DISC%DynRup%filter(21,23) = -0.165204631783187367D0
      DISC%DynRup%filter(21,24) = 0.694046397267487347D-1
      DISC%DynRup%filter(21,25) = -0.174720478217044832D-1
      DISC%DynRup%filter(22,1) = 0.133752438264239876D-2
      DISC%DynRup%filter(22,2) = -0.154504924288674950D-2
      DISC%DynRup%filter(22,3) = -0.914951099262413178D-2
      DISC%DynRup%filter(22,4) = -0.136153368837328825D-1
      DISC%DynRup%filter(22,5) = -0.862869798203701535D-2
      DISC%DynRup%filter(22,6) = 0.300729438053620353D-1
      DISC%DynRup%filter(22,7) = 0.497420912419659009D-1
      DISC%DynRup%filter(22,8) = 0.252958458220055374D-1
      DISC%DynRup%filter(22,9) = -0.190586963767069346D-3
      DISC%DynRup%filter(22,10) = 0.791935206743413635D-3
      DISC%DynRup%filter(22,11) = -0.811822972103533741D-1
      DISC%DynRup%filter(22,12) = -0.308293428343421674D-1
      DISC%DynRup%filter(22,13) = -0.573030758037888430D-1
      DISC%DynRup%filter(22,14) = -0.484956220630015575D-1
      DISC%DynRup%filter(22,15) = 0.574215110691921257D-1
      DISC%DynRup%filter(22,16) = -0.646255909890447455D-1
      DISC%DynRup%filter(22,17) = 0.255700831051602684D0
      DISC%DynRup%filter(22,18) = 0.564633124153447463D-1
      DISC%DynRup%filter(22,19) = -0.186018772310679009D0
      DISC%DynRup%filter(22,20) = 0.761493217634536118D-1
      DISC%DynRup%filter(22,21) = 0.143233063616833328D0
      DISC%DynRup%filter(22,22) = 0.599987294359594858D0
      DISC%DynRup%filter(22,23) = 0.308515129718571812D0
      DISC%DynRup%filter(22,24) = -0.137483022124214871D0
      DISC%DynRup%filter(22,25) = 0.343561216704810907D-1
      DISC%DynRup%filter(23,1) = -0.299072503025367463D-2
      DISC%DynRup%filter(23,2) = -0.769784473298415137D-2
      DISC%DynRup%filter(23,3) = -0.102239316011195173D-1
      DISC%DynRup%filter(23,4) = -0.769784471909979329D-2
      DISC%DynRup%filter(23,5) = -0.299072502422171411D-2
      DISC%DynRup%filter(23,6) = 0.348768668525689315D-1
      DISC%DynRup%filter(23,7) = 0.212823932001052728D-1
      DISC%DynRup%filter(23,8) = -0.660629042478183635D-2
      DISC%DynRup%filter(23,9) = 0.212823933527922990D-1
      DISC%DynRup%filter(23,10) = 0.348768671004212813D-1
      DISC%DynRup%filter(23,11) = 0.357034995661505913D-1
      DISC%DynRup%filter(23,12) = -0.482113393187932476D-1
      DISC%DynRup%filter(23,13) = -0.135373148102624796D0
      DISC%DynRup%filter(23,14) = -0.482113391119522935D-1
      DISC%DynRup%filter(23,15) = 0.357034997470364462D-1
      DISC%DynRup%filter(23,16) = -0.647261674463895931D-1
      DISC%DynRup%filter(23,17) = 0.475048125792903800D-1
      DISC%DynRup%filter(23,18) = 0.172111812664594271D0
      DISC%DynRup%filter(23,19) = 0.475048120704411533D-1
      DISC%DynRup%filter(23,20) = -0.647261679092332454D-1
      DISC%DynRup%filter(23,21) = -0.688032746558202707D-1
      DISC%DynRup%filter(23,22) = 0.259565952514154707D0
      DISC%DynRup%filter(23,23) = 0.567083235538852648D0
      DISC%DynRup%filter(23,24) = 0.259565952931338939D0
      DISC%DynRup%filter(23,25) = -0.688032749355810130D-1
      DISC%DynRup%filter(24,1) = -0.862869811953724530D-2
      DISC%DynRup%filter(24,2) = -0.136153373078807278D-1
      DISC%DynRup%filter(24,3) = -0.914951156174244502D-2
      DISC%DynRup%filter(24,4) = -0.154504964302991962D-2
      DISC%DynRup%filter(24,5) = 0.133752428996043097D-2
      DISC%DynRup%filter(24,6) = 0.791935600056112342D-3
      DISC%DynRup%filter(24,7) = -0.190586981101358054D-3
      DISC%DynRup%filter(24,8) = 0.252958470321224434D-1
      DISC%DynRup%filter(24,9) = 0.497420913471086495D-1
      DISC%DynRup%filter(24,10) = 0.300729439625743965D-1
      DISC%DynRup%filter(24,11) = 0.574215117497045707D-1
      DISC%DynRup%filter(24,12) = -0.484956218640642811D-1
      DISC%DynRup%filter(24,13) = -0.573030764756502151D-1
      DISC%DynRup%filter(24,14) = -0.308293424522680590D-1
      DISC%DynRup%filter(24,15) = -0.811822962534763382D-1
      DISC%DynRup%filter(24,16) = 0.761493214690276005D-1
      DISC%DynRup%filter(24,17) = -0.186018772198320470D0
      DISC%DynRup%filter(24,18) = 0.564633116954927697D-1
      DISC%DynRup%filter(24,19) = 0.255700830053596784D0
      DISC%DynRup%filter(24,20) = -0.646255917165202143D-1
      DISC%DynRup%filter(24,21) = 0.343561219316504476D-1
      DISC%DynRup%filter(24,22) = -0.137483022124214871D0
      DISC%DynRup%filter(24,23) = 0.308515130214429056D0
      DISC%DynRup%filter(24,24) = 0.599987294366183588D0
      DISC%DynRup%filter(24,25) = 0.143233063918103615D0
      DISC%DynRup%filter(25,1) = -0.132846850766015975D-1
      DISC%DynRup%filter(25,2) = -0.174312952428093623D-1
      DISC%DynRup%filter(25,3) = -0.718107828437513068D-2
      DISC%DynRup%filter(25,4) = 0.270200373795519231D-2
      DISC%DynRup%filter(25,5) = 0.359398099323100751D-2
      DISC%DynRup%filter(25,6) = -0.522213860158539239D-1
      DISC%DynRup%filter(25,7) = 0.159983019709787405D-2
      DISC%DynRup%filter(25,8) = 0.837433994356796135D-1
      DISC%DynRup%filter(25,9) = 0.607519624582160828D-1
      DISC%DynRup%filter(25,10) = 0.118384258097985803D-1
      DISC%DynRup%filter(25,11) = -0.744555385431907718D-1
      DISC%DynRup%filter(25,12) = 0.116000266416708298D0
      DISC%DynRup%filter(25,13) = 0.857282372260675662D-1
      DISC%DynRup%filter(25,14) = -0.164000698989799493D0
      DISC%DynRup%filter(25,15) = -0.123661090933634352D0
      DISC%DynRup%filter(25,16) = -0.519979297725346595D-1
      DISC%DynRup%filter(25,17) = 0.153833317730801389D0
      DISC%DynRup%filter(25,18) = -0.155415025477398716D0
      DISC%DynRup%filter(25,19) = -0.130553614444159916D0
      DISC%DynRup%filter(25,20) = 0.321802348707221408D0
      DISC%DynRup%filter(25,21) = -0.174720478217045040D-1
      DISC%DynRup%filter(25,22) = 0.694046391991465367D-1
      DISC%DynRup%filter(25,23) = -0.165204632454925310D0
      DISC%DynRup%filter(25,24) = 0.289352774389713852D0
      DISC%DynRup%filter(25,25) = 0.772527849362193320D0
    ELSEIF (DISC%Galerkin%nBndGP.EQ.36) THEN
      DISC%DynRup%filter(1,1) = 0.10344794618160320508D0
      DISC%DynRup%filter(1,2) = 0.20261469682791441808D0
      DISC%DynRup%filter(1,3) = 0.23237333056443609384D0
      DISC%DynRup%filter(1,4) = 0.19849681271608108924D0
      DISC%DynRup%filter(1,5) = 0.13023928376406040949D0
      DISC%DynRup%filter(1,6) = 0.54975960270491760471D-1
      DISC%DynRup%filter(1,7) = 0.16627029334871033120D0
      DISC%DynRup%filter(1,8) = 0.25588401720935710860D0
      DISC%DynRup%filter(1,9) = 0.15146543604223195391D0
      DISC%DynRup%filter(1,10) = -0.37776490834449138549D-1
      DISC%DynRup%filter(1,11) = -0.14851045968043987590D0
      DISC%DynRup%filter(1,12) = -0.10463256411369722162D0
      DISC%DynRup%filter(1,13) = -0.32385173695976836512D-1
      DISC%DynRup%filter(1,14) = -0.96631204372499433219D-1
      DISC%DynRup%filter(1,15) = -0.15179509630842771613D0
      DISC%DynRup%filter(1,16) = -0.13758525949993794141D0
      DISC%DynRup%filter(1,17) = -0.67274566826509101101D-1
      DISC%DynRup%filter(1,18) = -0.13487001776236168076D-1
      DISC%DynRup%filter(1,19) = -0.34115111515859130473D-1
      DISC%DynRup%filter(1,20) = -0.23948494500240581391D-1
      DISC%DynRup%filter(1,21) = 0.71069834211476229585D-1
      DISC%DynRup%filter(1,22) = 0.19099218186869933284D0
      DISC%DynRup%filter(1,23) = 0.22974764380404569350D0
      DISC%DynRup%filter(1,24) = 0.13388046100906394094D0
      DISC%DynRup%filter(1,25) = 0.44571992496248307149D-1
      DISC%DynRup%filter(1,26) = 0.58434604602357794117D-1
      DISC%DynRup%filter(1,27) = -0.14173314461728794553D-1
      DISC%DynRup%filter(1,28) = -0.14414066714251249842D0
      DISC%DynRup%filter(1,29) = -0.22028779273274315184D0
      DISC%DynRup%filter(1,30) = -0.14289919456203149383D0
      DISC%DynRup%filter(1,31) = -0.18139581719029019985D-1
      DISC%DynRup%filter(1,32) = -0.26148532006393934710D-1
      DISC%DynRup%filter(1,33) = -0.10706534884184134139D-2
      DISC%DynRup%filter(1,34) = 0.49948460794472411537D-1
      DISC%DynRup%filter(1,35) = 0.84034263868231536166D-1
      DISC%DynRup%filter(1,36) = 0.56553939657648840853D-1
      DISC%DynRup%filter(2,1) = 0.96221057551409406467D-1
      DISC%DynRup%filter(2,2) = 0.19201576221703826374D0
      DISC%DynRup%filter(2,3) = 0.22769442271473357606D0
      DISC%DynRup%filter(2,4) = 0.20367220525450506805D0
      DISC%DynRup%filter(2,5) = 0.14069437223834175731D0
      DISC%DynRup%filter(2,6) = 0.61850210348558904462D-1
      DISC%DynRup%filter(2,7) = 0.12151848376181202864D0
      DISC%DynRup%filter(2,8) = 0.19491508262141035110D0
      DISC%DynRup%filter(2,9) = 0.13146964824992001195D0
      DISC%DynRup%filter(2,10) = -0.27915679713924568763D-2
      DISC%DynRup%filter(2,11) = -0.91884282593422018384D-1
      DISC%DynRup%filter(2,12) = -0.70527132096614758577D-1
      DISC%DynRup%filter(2,13) = -0.45889843248070265827D-1
      DISC%DynRup%filter(2,14) = -0.98636286282456634094D-1
      DISC%DynRup%filter(2,15) = -0.12684540339885629196D0
      DISC%DynRup%filter(2,16) = -0.11753087452817692194D0
      DISC%DynRup%filter(2,17) = -0.78307423032072959443D-1
      DISC%DynRup%filter(2,18) = -0.31948471989954123408D-1
      DISC%DynRup%filter(2,19) = -0.11373061794892404058D-1
      DISC%DynRup%filter(2,20) = 0.13154221494323543761D-1
      DISC%DynRup%filter(2,21) = 0.90891747804624335086D-1
      DISC%DynRup%filter(2,22) = 0.17402024774773752471D0
      DISC%DynRup%filter(2,23) = 0.19182695420236718155D0
      DISC%DynRup%filter(2,24) = 0.10910640542302530325D0
      DISC%DynRup%filter(2,25) = 0.27750402811168068894D-1
      DISC%DynRup%filter(2,26) = 0.22543001068960017992D-1
      DISC%DynRup%filter(2,27) = -0.47799465518217907118D-1
      DISC%DynRup%filter(2,28) = -0.14076212485409768055D0
      DISC%DynRup%filter(2,29) = -0.17561223225222670677D0
      DISC%DynRup%filter(2,30) = -0.10461395305599563000D0
      DISC%DynRup%filter(2,31) = -0.12417852418716866442D-1
      DISC%DynRup%filter(2,32) = -0.12625789637249118815D-1
      DISC%DynRup%filter(2,33) = 0.13577584202894867827D-1
      DISC%DynRup%filter(2,34) = 0.50655396038971459605D-1
      DISC%DynRup%filter(2,35) = 0.66080961561950799267D-1
      DISC%DynRup%filter(2,36) = 0.39907597358660279311D-1
      DISC%DynRup%filter(3,1) = 0.85082405031954511660D-1
      DISC%DynRup%filter(3,2) = 0.17555236560302527657D0
      DISC%DynRup%filter(3,3) = 0.22023636327875971760D0
      DISC%DynRup%filter(3,4) = 0.21156697992925871093D0
      DISC%DynRup%filter(3,5) = 0.15703123956096640827D0
      DISC%DynRup%filter(3,6) = 0.72678676920622351347D-1
      DISC%DynRup%filter(3,7) = 0.55458358953603155809D-1
      DISC%DynRup%filter(3,8) = 0.10136307899024199321D0
      DISC%DynRup%filter(3,9) = 0.95168905854787693825D-1
      DISC%DynRup%filter(3,10) = 0.46693871585625664893D-1
      DISC%DynRup%filter(3,11) = -0.21522984852973878068D-2
      DISC%DynRup%filter(3,12) = -0.13831684927247961862D-1
      DISC%DynRup%filter(3,13) = -0.55579062513789173414D-1
      DISC%DynRup%filter(3,14) = -0.97797786906874173670D-1
      DISC%DynRup%filter(3,15) = -0.10395433240734094307D0
      DISC%DynRup%filter(3,16) = -0.10083463668241205035D0
      DISC%DynRup%filter(3,17) = -0.90616286551136080713D-1
      DISC%DynRup%filter(3,18) = -0.50376197418034776753D-1
      DISC%DynRup%filter(3,19) = 0.26021886441301099776D-1
      DISC%DynRup%filter(3,20) = 0.70077523861382173110D-1
      DISC%DynRup%filter(3,21) = 0.11901629366961804516D0
      DISC%DynRup%filter(3,22) = 0.14841033457698928250D0
      DISC%DynRup%filter(3,23) = 0.13416958479134085593D0
      DISC%DynRup%filter(3,24) = 0.69930891536554027588D-1
      DISC%DynRup%filter(3,25) = -0.51894926097970481778D-2
      DISC%DynRup%filter(3,26) = -0.36853380711904490645D-1
      DISC%DynRup%filter(3,27) = -0.90671901228624376934D-1
      DISC%DynRup%filter(3,28) = -0.12447559102023029462D0
      DISC%DynRup%filter(3,29) = -0.10852757705183020996D0
      DISC%DynRup%filter(3,30) = -0.52776429178023416707D-1
      DISC%DynRup%filter(3,31) = -0.39201475285147048003D-3
      DISC%DynRup%filter(3,32) = 0.10468315374499617203D-1
      DISC%DynRup%filter(3,33) = 0.32045661403442071543D-1
      DISC%DynRup%filter(3,34) = 0.45712237733247235354D-1
      DISC%DynRup%filter(3,35) = 0.39055302713060840371D-1
      DISC%DynRup%filter(3,36) = 0.18288394635113126338D-1
      DISC%DynRup%filter(4,1) = 0.72678676920622351354D-1
      DISC%DynRup%filter(4,2) = 0.15703123956096640828D0
      DISC%DynRup%filter(4,3) = 0.21156697992925871093D0
      DISC%DynRup%filter(4,4) = 0.22023636327875971765D0
      DISC%DynRup%filter(4,5) = 0.17555236560302527661D0
      DISC%DynRup%filter(4,6) = 0.85082405031954511675D-1
      DISC%DynRup%filter(4,7) = -0.13831684927247961891D-1
      DISC%DynRup%filter(4,8) = -0.21522984852973878081D-2
      DISC%DynRup%filter(4,9) = 0.46693871585625664666D-1
      DISC%DynRup%filter(4,10) = 0.95168905854787694137D-1
      DISC%DynRup%filter(4,11) = 0.10136307899024199326D0
      DISC%DynRup%filter(4,12) = 0.55458358953603155779D-1
      DISC%DynRup%filter(4,13) = -0.50376197418034776808D-1
      DISC%DynRup%filter(4,14) = -0.90616286551136080888D-1
      DISC%DynRup%filter(4,15) = -0.10083463668241205041D0
      DISC%DynRup%filter(4,16) = -0.10395433240734094333D0
      DISC%DynRup%filter(4,17) = -0.97797786906874173969D-1
      DISC%DynRup%filter(4,18) = -0.55579062513789173330D-1
      DISC%DynRup%filter(4,19) = 0.69930891536554027357D-1
      DISC%DynRup%filter(4,20) = 0.13416958479134085550D0
      DISC%DynRup%filter(4,21) = 0.14841033457698928236D0
      DISC%DynRup%filter(4,22) = 0.11901629366961804510D0
      DISC%DynRup%filter(4,23) = 0.70077523861382173338D-1
      DISC%DynRup%filter(4,24) = 0.26021886441301099883D-1
      DISC%DynRup%filter(4,25) = -0.52776429178023417009D-1
      DISC%DynRup%filter(4,26) = -0.10852757705183021045D0
      DISC%DynRup%filter(4,27) = -0.12447559102023029464D0
      DISC%DynRup%filter(4,28) = -0.90671901228624376874D-1
      DISC%DynRup%filter(4,29) = -0.36853380711904490416D-1
      DISC%DynRup%filter(4,30) = -0.51894926097970480669D-2
      DISC%DynRup%filter(4,31) = 0.18288394635113126501D-1
      DISC%DynRup%filter(4,32) = 0.39055302713060840616D-1
      DISC%DynRup%filter(4,33) = 0.45712237733247235496D-1
      DISC%DynRup%filter(4,34) = 0.32045661403442071481D-1
      DISC%DynRup%filter(4,35) = 0.10468315374499617098D-1
      DISC%DynRup%filter(4,36) = -0.39201475285147055436D-3
      DISC%DynRup%filter(5,1) = 0.61850210348558904450D-1
      DISC%DynRup%filter(5,2) = 0.14069437223834175731D0
      DISC%DynRup%filter(5,3) = 0.20367220525450506802D0
      DISC%DynRup%filter(5,4) = 0.22769442271473357614D0
      DISC%DynRup%filter(5,5) = 0.19201576221703826385D0
      DISC%DynRup%filter(5,6) = 0.96221057551409406512D-1
      DISC%DynRup%filter(5,7) = -0.70527132096614758807D-1
      DISC%DynRup%filter(5,8) = -0.91884282593422018526D-1
      DISC%DynRup%filter(5,9) = -0.27915679713924569718D-2
      DISC%DynRup%filter(5,10) = 0.13146964824992001259D0
      DISC%DynRup%filter(5,11) = 0.19491508262141035151D0
      DISC%DynRup%filter(5,12) = 0.12151848376181202875D0
      DISC%DynRup%filter(5,13) = -0.31948471989954123679D-1
      DISC%DynRup%filter(5,14) = -0.78307423032072960249D-1
      DISC%DynRup%filter(5,15) = -0.11753087452817692225D0
      DISC%DynRup%filter(5,16) = -0.12684540339885629233D0
      DISC%DynRup%filter(5,17) = -0.98636286282456634770D-1
      DISC%DynRup%filter(5,18) = -0.45889843248070265729D-1
      DISC%DynRup%filter(5,19) = 0.10910640542302530310D0
      DISC%DynRup%filter(5,20) = 0.19182695420236718118D0
      DISC%DynRup%filter(5,21) = 0.17402024774773752438D0
      DISC%DynRup%filter(5,22) = 0.90891747804624334463D-1
      DISC%DynRup%filter(5,23) = 0.13154221494323543418D-1
      DISC%DynRup%filter(5,24) = -0.11373061794892404233D-1
      DISC%DynRup%filter(5,25) = -0.10461395305599563027D0
      DISC%DynRup%filter(5,26) = -0.17561223225222670713D0
      DISC%DynRup%filter(5,27) = -0.14076212485409768034D0
      DISC%DynRup%filter(5,28) = -0.47799465518217906969D-1
      DISC%DynRup%filter(5,29) = 0.22543001068960018186D-1
      DISC%DynRup%filter(5,30) = 0.27750402811168069013D-1
      DISC%DynRup%filter(5,31) = 0.39907597358660279306D-1
      DISC%DynRup%filter(5,32) = 0.66080961561950799387D-1
      DISC%DynRup%filter(5,33) = 0.50655396038971459705D-1
      DISC%DynRup%filter(5,34) = 0.13577584202894867854D-1
      DISC%DynRup%filter(5,35) = -0.12625789637249118755D-1
      DISC%DynRup%filter(5,36) = -0.12417852418716866468D-1
      DISC%DynRup%filter(6,1) = 0.54975960270491760469D-1
      DISC%DynRup%filter(6,2) = 0.13023928376406040951D0
      DISC%DynRup%filter(6,3) = 0.19849681271608108921D0
      DISC%DynRup%filter(6,4) = 0.23237333056443609388D0
      DISC%DynRup%filter(6,5) = 0.20261469682791441819D0
      DISC%DynRup%filter(6,6) = 0.10344794618160320512D0
      DISC%DynRup%filter(6,7) = -0.10463256411369722182D0
      DISC%DynRup%filter(6,8) = -0.14851045968043987591D0
      DISC%DynRup%filter(6,9) = -0.37776490834449138542D-1
      DISC%DynRup%filter(6,10) = 0.15146543604223195454D0
      DISC%DynRup%filter(6,11) = 0.25588401720935710894D0
      DISC%DynRup%filter(6,12) = 0.16627029334871033120D0
      DISC%DynRup%filter(6,13) = -0.13487001776236168331D-1
      DISC%DynRup%filter(6,14) = -0.67274566826509102163D-1
      DISC%DynRup%filter(6,15) = -0.13758525949993794180D0
      DISC%DynRup%filter(6,16) = -0.15179509630842771650D0
      DISC%DynRup%filter(6,17) = -0.96631204372499434045D-1
      DISC%DynRup%filter(6,18) = -0.32385173695976836372D-1
      DISC%DynRup%filter(6,19) = 0.13388046100906394061D0
      DISC%DynRup%filter(6,20) = 0.22974764380404569284D0
      DISC%DynRup%filter(6,21) = 0.19099218186869933230D0
      DISC%DynRup%filter(6,22) = 0.71069834211476228896D-1
      DISC%DynRup%filter(6,23) = -0.23948494500240581782D-1
      DISC%DynRup%filter(6,24) = -0.34115111515859130715D-1
      DISC%DynRup%filter(6,25) = -0.14289919456203149415D0
      DISC%DynRup%filter(6,26) = -0.22028779273274315231D0
      DISC%DynRup%filter(6,27) = -0.14414066714251249823D0
      DISC%DynRup%filter(6,28) = -0.14173314461728794500D-1
      DISC%DynRup%filter(6,29) = 0.58434604602357794071D-1
      DISC%DynRup%filter(6,30) = 0.44571992496248307092D-1
      DISC%DynRup%filter(6,31) = 0.56553939657648840931D-1
      DISC%DynRup%filter(6,32) = 0.84034263868231536427D-1
      DISC%DynRup%filter(6,33) = 0.49948460794472411703D-1
      DISC%DynRup%filter(6,34) = -0.10706534884184133147D-2
      DISC%DynRup%filter(6,35) = -0.26148532006393934523D-1
      DISC%DynRup%filter(6,36) = -0.18139581719029019911D-1
      DISC%DynRup%filter(7,1) = 0.33054590688524898006D-1
      DISC%DynRup%filter(7,2) = 0.50869829373861301843D-1
      DISC%DynRup%filter(7,3) = 0.30111380036689864635D-1
      DISC%DynRup%filter(7,4) = -0.75099791852939196535D-2
      DISC%DynRup%filter(7,5) = -0.29523929734136698171D-1
      DISC%DynRup%filter(7,6) = -0.20800989219496839078D-1
      DISC%DynRup%filter(7,7) = 0.28089536413802033096D0
      DISC%DynRup%filter(7,8) = 0.41641097276451850671D0
      DISC%DynRup%filter(7,9) = 0.24939310387237093665D0
      DISC%DynRup%filter(7,10) = 0.10363862498204974387D-1
      DISC%DynRup%filter(7,11) = -0.96632847278931262231D-1
      DISC%DynRup%filter(7,12) = -0.64510174046214008338D-1
      DISC%DynRup%filter(7,13) = 0.23350747898705560797D0
      DISC%DynRup%filter(7,14) = 0.21830190817035490730D0
      DISC%DynRup%filter(7,15) = -0.48278671672797815304D-1
      DISC%DynRup%filter(7,16) = -0.13244683565230673016D0
      DISC%DynRup%filter(7,17) = 0.12854248701525850873D-1
      DISC%DynRup%filter(7,18) = 0.76401461295571615913D-1
      DISC%DynRup%filter(7,19) = -0.73234897927734001422D-1
      DISC%DynRup%filter(7,20) = -0.19105206876581538567D0
      DISC%DynRup%filter(7,21) = -0.19962690645621961660D0
      DISC%DynRup%filter(7,22) = -0.51823883897052149640D-1
      DISC%DynRup%filter(7,23) = 0.60142452390046288376D-1
      DISC%DynRup%filter(7,24) = 0.45828894820281634737D-1
      DISC%DynRup%filter(7,25) = 0.12101878943131175453D-1
      DISC%DynRup%filter(7,26) = 0.96547273571849273958D-1
      DISC%DynRup%filter(7,27) = 0.20187831218367065602D0
      DISC%DynRup%filter(7,28) = 0.14719203851871083685D0
      DISC%DynRup%filter(7,29) = -0.48742503266255410696D-1
      DISC%DynRup%filter(7,30) = -0.10686825102041950498D0
      DISC%DynRup%filter(7,31) = 0.95096704866001113056D-3
      DISC%DynRup%filter(7,32) = -0.25974425535746353074D-1
      DISC%DynRup%filter(7,33) = -0.75627227904733683067D-1
      DISC%DynRup%filter(7,34) = -0.69859031553306201934D-1
      DISC%DynRup%filter(7,35) = 0.15205556008914496048D-1
      DISC%DynRup%filter(7,36) = 0.50501049104496412903D-1
      DISC%DynRup%filter(8,1) = 0.24157915770394060035D-1
      DISC%DynRup%filter(8,2) = 0.38749184507412221831D-1
      DISC%DynRup%filter(8,3) = 0.26136210644384219737D-1
      DISC%DynRup%filter(8,4) = -0.55496465914119449192D-3
      DISC%DynRup%filter(8,5) = -0.18266626531202157877D-1
      DISC%DynRup%filter(8,6) = -0.14020817771698541871D-1
      DISC%DynRup%filter(8,7) = 0.19775221048966361194D0
      DISC%DynRup%filter(8,8) = 0.32920828594743790461D0
      DISC%DynRup%filter(8,9) = 0.26094813056674798238D0
      DISC%DynRup%filter(8,10) = 0.89647549170853697479D-1
      DISC%DynRup%filter(8,11) = -0.35745270767766595012D-1
      DISC%DynRup%filter(8,12) = -0.45890623458967122475D-1
      DISC%DynRup%filter(8,13) = 0.10367086296549589197D0
      DISC%DynRup%filter(8,14) = 0.15391542168066790112D0
      DISC%DynRup%filter(8,15) = 0.90032987978306844644D-1
      DISC%DynRup%filter(8,16) = 0.10863035751443319831D-1
      DISC%DynRup%filter(8,17) = -0.42471592581343909218D-2
      DISC%DynRup%filter(8,18) = 0.61044407116238695482D-2
      DISC%DynRup%filter(8,19) = -0.90730003261533762894D-1
      DISC%DynRup%filter(8,20) = -0.14533844980439390425D0
      DISC%DynRup%filter(8,21) = -0.12965021920634913386D0
      DISC%DynRup%filter(8,22) = -0.73944164122999387974D-1
      DISC%DynRup%filter(8,23) = 0.13349724130328170733D-2
      DISC%DynRup%filter(8,24) = 0.28561454145750140391D-1
      DISC%DynRup%filter(8,25) = 0.45849984785055837925D-1
      DISC%DynRup%filter(8,26) = 0.10772732340819268300D0
      DISC%DynRup%filter(8,27) = 0.12062433672280619522D0
      DISC%DynRup%filter(8,28) = 0.60209561695383524159D-1
      DISC%DynRup%filter(8,29) = -0.91548021285703139402D-2
      DISC%DynRup%filter(8,30) = -0.23147655552180899199D-1
      DISC%DynRup%filter(8,31) = -0.12335169824637907770D-1
      DISC%DynRup%filter(8,32) = -0.37266529091756304699D-1
      DISC%DynRup%filter(8,33) = -0.47220300283489775644D-1
      DISC%DynRup%filter(8,34) = -0.21188884541088466275D-1
      DISC%DynRup%filter(8,35) = 0.59867020394229667707D-2
      DISC%DynRup%filter(8,36) = 0.72210688698341693921D-2
      DISC%DynRup%filter(9,1) = 0.11025140562085025338D-1
      DISC%DynRup%filter(9,2) = 0.20151014468496785364D-1
      DISC%DynRup%filter(9,3) = 0.18919610749151107968D-1
      DISC%DynRup%filter(9,4) = 0.92827575018972629529D-2
      DISC%DynRup%filter(9,5) = -0.42787766857325447964D-3
      DISC%DynRup%filter(9,6) = -0.27497436529083198626D-2
      DISC%DynRup%filter(9,7) = 0.91314115196887632047D-1
      DISC%DynRup%filter(9,8) = 0.20119096934610798801D0
      DISC%DynRup%filter(9,9) = 0.24644753626435536806D0
      DISC%DynRup%filter(9,10) = 0.18405473459668821134D0
      DISC%DynRup%filter(9,11) = 0.69118246902226537876D-1
      DISC%DynRup%filter(9,12) = 0.37946796417037422225D-2
      DISC%DynRup%filter(9,13) = -0.17677009180408916015D-1
      DISC%DynRup%filter(9,14) = 0.69415420164693236099D-1
      DISC%DynRup%filter(9,15) = 0.19173874038189951229D0
      DISC%DynRup%filter(9,16) = 0.15698182716390740230D0
      DISC%DynRup%filter(9,17) = 0.83753989274710172054D-2
      DISC%DynRup%filter(9,18) = -0.48494787628158816256D-1
      DISC%DynRup%filter(9,19) = -0.73092455442834823102D-1
      DISC%DynRup%filter(9,20) = -0.99960299471808714163D-1
      DISC%DynRup%filter(9,21) = -0.82622353785270520782D-1
      DISC%DynRup%filter(9,22) = -0.78105289732597335960D-1
      DISC%DynRup%filter(9,23) = -0.57010939396588560158D-1
      DISC%DynRup%filter(9,24) = -0.18975072007393278549D-1
      DISC%DynRup%filter(9,25) = 0.73916797089653601146D-1
      DISC%DynRup%filter(9,26) = 0.93001345437058255497D-1
      DISC%DynRup%filter(9,27) = 0.31773219700346755558D-1
      DISC%DynRup%filter(9,28) = 0.31021492142714602556D-2
      DISC%DynRup%filter(9,29) = 0.46421562994489277240D-1
      DISC%DynRup%filter(9,30) = 0.53893674494867677615D-1
      DISC%DynRup%filter(9,31) = -0.27690554765492817997D-1
      DISC%DynRup%filter(9,32) = -0.36406844403200360945D-1
      DISC%DynRup%filter(9,33) = -0.75748988164467163188D-2
      DISC%DynRup%filter(9,34) = 0.87843686690850453692D-2
      DISC%DynRup%filter(9,35) = -0.16336626788341425456D-1
      DISC%DynRup%filter(9,36) = -0.25578556727319042878D-1
      DISC%DynRup%filter(10,1) = -0.27497436529083198631D-2
      DISC%DynRup%filter(10,2) = -0.42787766857325446467D-3
      DISC%DynRup%filter(10,3) = 0.92827575018972629983D-2
      DISC%DynRup%filter(10,4) = 0.18919610749151108025D-1
      DISC%DynRup%filter(10,5) = 0.20151014468496785465D-1
      DISC%DynRup%filter(10,6) = 0.11025140562085025386D-1
      DISC%DynRup%filter(10,7) = 0.37946796417037420534D-2
      DISC%DynRup%filter(10,8) = 0.69118246902226537696D-1
      DISC%DynRup%filter(10,9) = 0.18405473459668821134D0
      DISC%DynRup%filter(10,10) = 0.24644753626435536811D0
      DISC%DynRup%filter(10,11) = 0.20119096934610798821D0
      DISC%DynRup%filter(10,12) = 0.91314115196887632145D-1
      DISC%DynRup%filter(10,13) = -0.48494787628158816030D-1
      DISC%DynRup%filter(10,14) = 0.83753989274710171571D-2
      DISC%DynRup%filter(10,15) = 0.15698182716390740216D0
      DISC%DynRup%filter(10,16) = 0.19173874038189951190D0
      DISC%DynRup%filter(10,17) = 0.69415420164693235620D-1
      DISC%DynRup%filter(10,18) = -0.17677009180408916224D-1
      DISC%DynRup%filter(10,19) = -0.18975072007393278412D-1
      DISC%DynRup%filter(10,20) = -0.57010939396588559833D-1
      DISC%DynRup%filter(10,21) = -0.78105289732597335943D-1
      DISC%DynRup%filter(10,22) = -0.82622353785270521016D-1
      DISC%DynRup%filter(10,23) = -0.99960299471808714493D-1
      DISC%DynRup%filter(10,24) = -0.73092455442834823204D-1
      DISC%DynRup%filter(10,25) = 0.53893674494867677439D-1
      DISC%DynRup%filter(10,26) = 0.46421562994489277088D-1
      DISC%DynRup%filter(10,27) = 0.31021492142714601998D-2
      DISC%DynRup%filter(10,28) = 0.31773219700346755777D-1
      DISC%DynRup%filter(10,29) = 0.93001345437058255632D-1
      DISC%DynRup%filter(10,30) = 0.73916797089653601210D-1
      DISC%DynRup%filter(10,31) = -0.25578556727319042865D-1
      DISC%DynRup%filter(10,32) = -0.16336626788341425461D-1
      DISC%DynRup%filter(10,33) = 0.87843686690850453490D-2
      DISC%DynRup%filter(10,34) = -0.75748988164467162860D-2
      DISC%DynRup%filter(10,35) = -0.36406844403200360927D-1
      DISC%DynRup%filter(10,36) = -0.27690554765492818111D-1
      DISC%DynRup%filter(11,1) = -0.14020817771698541871D-1
      DISC%DynRup%filter(11,2) = -0.18266626531202157850D-1
      DISC%DynRup%filter(11,3) = -0.55496465914119448933D-3
      DISC%DynRup%filter(11,4) = 0.26136210644384219753D-1
      DISC%DynRup%filter(11,5) = 0.38749184507412221911D-1
      DISC%DynRup%filter(11,6) = 0.24157915770394060069D-1
      DISC%DynRup%filter(11,7) = -0.45890623458967122644D-1
      DISC%DynRup%filter(11,8) = -0.35745270767766595013D-1
      DISC%DynRup%filter(11,9) = 0.89647549170853697711D-1
      DISC%DynRup%filter(11,10) = 0.26094813056674798266D0
      DISC%DynRup%filter(11,11) = 0.32920828594743790482D0
      DISC%DynRup%filter(11,12) = 0.19775221048966361195D0
      DISC%DynRup%filter(11,13) = 0.61044407116238696387D-2
      DISC%DynRup%filter(11,14) = -0.42471592581343914217D-2
      DISC%DynRup%filter(11,15) = 0.10863035751443319610D-1
      DISC%DynRup%filter(11,16) = 0.90032987978306844250D-1
      DISC%DynRup%filter(11,17) = 0.15391542168066790033D0
      DISC%DynRup%filter(11,18) = 0.10367086296549589177D0
      DISC%DynRup%filter(11,19) = 0.28561454145750140293D-1
      DISC%DynRup%filter(11,20) = 0.13349724130328169562D-2
      DISC%DynRup%filter(11,21) = -0.73944164122999388340D-1
      DISC%DynRup%filter(11,22) = -0.12965021920634913440D0
      DISC%DynRup%filter(11,23) = -0.14533844980439390472D0
      DISC%DynRup%filter(11,24) = -0.90730003261533763083D-1
      DISC%DynRup%filter(11,25) = -0.23147655552180899401D-1
      DISC%DynRup%filter(11,26) = -0.91548021285703141970D-2
      DISC%DynRup%filter(11,27) = 0.60209561695383524213D-1
      DISC%DynRup%filter(11,28) = 0.12062433672280619545D0
      DISC%DynRup%filter(11,29) = 0.10772732340819268310D0
      DISC%DynRup%filter(11,30) = 0.45849984785055837991D-1
      DISC%DynRup%filter(11,31) = 0.72210688698341694435D-2
      DISC%DynRup%filter(11,32) = 0.59867020394229669037D-2
      DISC%DynRup%filter(11,33) = -0.21188884541088466214D-1
      DISC%DynRup%filter(11,34) = -0.47220300283489775571D-1
      DISC%DynRup%filter(11,35) = -0.37266529091756304686D-1
      DISC%DynRup%filter(11,36) = -0.12335169824637907881D-1
      DISC%DynRup%filter(12,1) = -0.20800989219496839039D-1
      DISC%DynRup%filter(12,2) = -0.29523929734136698075D-1
      DISC%DynRup%filter(12,3) = -0.75099791852939196343D-2
      DISC%DynRup%filter(12,4) = 0.30111380036689864621D-1
      DISC%DynRup%filter(12,5) = 0.50869829373861301894D-1
      DISC%DynRup%filter(12,6) = 0.33054590688524898008D-1
      DISC%DynRup%filter(12,7) = -0.64510174046214008334D-1
      DISC%DynRup%filter(12,8) = -0.96632847278931261874D-1
      DISC%DynRup%filter(12,9) = 0.10363862498204974846D-1
      DISC%DynRup%filter(12,10) = 0.24939310387237093693D0
      DISC%DynRup%filter(12,11) = 0.41641097276451850672D0
      DISC%DynRup%filter(12,12) = 0.28089536413802033072D0
      DISC%DynRup%filter(12,13) = 0.76401461295571615861D-1
      DISC%DynRup%filter(12,14) = 0.12854248701525849953D-1
      DISC%DynRup%filter(12,15) = -0.13244683565230673060D0
      DISC%DynRup%filter(12,16) = -0.48278671672797815750D-1
      DISC%DynRup%filter(12,17) = 0.21830190817035490634D0
      DISC%DynRup%filter(12,18) = 0.23350747898705560778D0
      DISC%DynRup%filter(12,19) = 0.45828894820281634407D-1
      DISC%DynRup%filter(12,20) = 0.60142452390046287706D-1
      DISC%DynRup%filter(12,21) = -0.51823883897052150410D-1
      DISC%DynRup%filter(12,22) = -0.19962690645621961717D0
      DISC%DynRup%filter(12,23) = -0.19105206876581538589D0
      DISC%DynRup%filter(12,24) = -0.73234897927734001490D-1
      DISC%DynRup%filter(12,25) = -0.10686825102041950508D0
      DISC%DynRup%filter(12,26) = -0.48742503266255410907D-1
      DISC%DynRup%filter(12,27) = 0.14719203851871083711D0
      DISC%DynRup%filter(12,28) = 0.20187831218367065629D0
      DISC%DynRup%filter(12,29) = 0.96547273571849274000D-1
      DISC%DynRup%filter(12,30) = 0.12101878943131175466D-1
      DISC%DynRup%filter(12,31) = 0.50501049104496412946D-1
      DISC%DynRup%filter(12,32) = 0.15205556008914496235D-1
      DISC%DynRup%filter(12,33) = -0.69859031553306201851D-1
      DISC%DynRup%filter(12,34) = -0.75627227904733683007D-1
      DISC%DynRup%filter(12,35) = -0.25974425535746353152D-1
      DISC%DynRup%filter(12,36) = 0.95096704866001099305D-3
      DISC%DynRup%filter(13,1) = -0.28683166523695777717D-2
      DISC%DynRup%filter(13,2) = -0.85585118437885783570D-2
      DISC%DynRup%filter(13,3) = -0.13444312714728335990D-1
      DISC%DynRup%filter(13,4) = -0.12185764221894136923D-1
      DISC%DynRup%filter(13,5) = -0.59584290676013194340D-2
      DISC%DynRup%filter(13,6) = -0.11945278462447171106D-2
      DISC%DynRup%filter(13,7) = 0.10403142273835132519D0
      DISC%DynRup%filter(13,8) = 0.97257090830559120190D-1
      DISC%DynRup%filter(13,9) = -0.21508942342344918395D-1
      DISC%DynRup%filter(13,10) = -0.59007243835928168418D-1
      DISC%DynRup%filter(13,11) = 0.57267792297413695823D-2
      DISC%DynRup%filter(13,12) = 0.34038107697218555919D-1
      DISC%DynRup%filter(13,13) = 0.43524383591115798629D0
      DISC%DynRup%filter(13,14) = 0.41105676050149142962D0
      DISC%DynRup%filter(13,15) = -0.31564449457014120674D-1
      DISC%DynRup%filter(13,16) = -0.16556894090670261634D0
      DISC%DynRup%filter(13,17) = 0.99542907191781193916D-3
      DISC%DynRup%filter(13,18) = 0.66379927208593731291D-1
      DISC%DynRup%filter(13,19) = 0.27885643693439252321D0
      DISC%DynRup%filter(13,20) = 0.14828203744988478017D0
      DISC%DynRup%filter(13,21) = -0.10089593060253854107D0
      DISC%DynRup%filter(13,22) = -0.48832598774410679931D-2
      DISC%DynRup%filter(13,23) = 0.37228966235896802952D-1
      DISC%DynRup%filter(13,24) = -0.36249712121918346030D-1
      DISC%DynRup%filter(13,25) = -0.10765859366938658523D0
      DISC%DynRup%filter(13,26) = -0.20343742162514419869D0
      DISC%DynRup%filter(13,27) = -0.54665075825707132844D-1
      DISC%DynRup%filter(13,28) = 0.13530827225983809852D0
      DISC%DynRup%filter(13,29) = 0.48676656652699504924D-1
      DISC%DynRup%filter(13,30) = -0.55874582497571185011D-1
      DISC%DynRup%filter(13,31) = 0.25275352189884356682D-1
      DISC%DynRup%filter(13,32) = 0.75687257086546390722D-1
      DISC%DynRup%filter(13,33) = 0.48243253046719999712D-1
      DISC%DynRup%filter(13,34) = -0.64384533332292553240D-1
      DISC%DynRup%filter(13,35) = -0.42179374373919118810D-1
      DISC%DynRup%filter(13,36) = 0.39800337769641430697D-1
      DISC%DynRup%filter(14,1) = -0.40644093126917938304D-2
      DISC%DynRup%filter(14,2) = -0.87360995845765704421D-2
      DISC%DynRup%filter(14,3) = -0.11234547829233182816D-1
      DISC%DynRup%filter(14,4) = -0.10409570988918571502D-1
      DISC%DynRup%filter(14,5) = -0.69355961340712809559D-2
      DISC%DynRup%filter(14,6) = -0.28296384971352657176D-2
      DISC%DynRup%filter(14,7) = 0.46187074682139227115D-1
      DISC%DynRup%filter(14,8) = 0.68571852037770439738D-1
      DISC%DynRup%filter(14,9) = 0.40111177052651735197D-1
      DISC%DynRup%filter(14,10) = 0.48396611079976173685D-2
      DISC%DynRup%filter(14,11) = -0.18921793089315763802D-2
      DISC%DynRup%filter(14,12) = 0.27196287459698406118D-2
      DISC%DynRup%filter(14,13) = 0.19520951257895503175D0
      DISC%DynRup%filter(14,14) = 0.32014999902461881942D0
      DISC%DynRup%filter(14,15) = 0.21036905067096118808D0
      DISC%DynRup%filter(14,16) = 0.21237920202203348685D-1
      DISC%DynRup%filter(14,17) = -0.30896646161334001832D-1
      DISC%DynRup%filter(14,18) = 0.47272601403983637432D-3
      DISC%DynRup%filter(14,19) = 0.70418655125613338786D-1
      DISC%DynRup%filter(14,20) = 0.18486811561989458950D0
      DISC%DynRup%filter(14,21) = 0.13477387092957948092D0
      DISC%DynRup%filter(14,22) = -0.37549931591339888912D-1
      DISC%DynRup%filter(14,23) = -0.47852086498362494163D-1
      DISC%DynRup%filter(14,24) = 0.17679914432891125667D-1
      DISC%DynRup%filter(14,25) = -0.96611766869650293034D-1
      DISC%DynRup%filter(14,26) = -0.91899732666443250163D-1
      DISC%DynRup%filter(14,27) = -0.51029992031605233439D-1
      DISC%DynRup%filter(14,28) = -0.33224993178662543793D-1
      DISC%DynRup%filter(14,29) = 0.11999354830556000781D-1
      DISC%DynRup%filter(14,30) = 0.23116385210533819981D-1
      DISC%DynRup%filter(14,31) = 0.35943631108943100860D-1
      DISC%DynRup%filter(14,32) = 0.35770179977729635932D-1
      DISC%DynRup%filter(14,33) = 0.26578475489506295606D-2
      DISC%DynRup%filter(14,34) = 0.16130581481627257875D-1
      DISC%DynRup%filter(14,35) = 0.11970899353252597617D-1
      DISC%DynRup%filter(14,36) = -0.20030847083922715487D-1
      DISC%DynRup%filter(15,1) = -0.49225720395378198376D-2
      DISC%DynRup%filter(15,2) = -0.86618346834658673229D-2
      DISC%DynRup%filter(15,3) = -0.92071126599200610125D-2
      DISC%DynRup%filter(15,4) = -0.89308048876615674343D-2
      DISC%DynRup%filter(15,5) = -0.80257776638945756738D-2
      DISC%DynRup%filter(15,6) = -0.44617604121467741134D-2
      DISC%DynRup%filter(15,7) = -0.78753983502977198810D-2
      DISC%DynRup%filter(15,8) = 0.30925711463459350581D-1
      DISC%DynRup%filter(15,9) = 0.85422762656326430151D-1
      DISC%DynRup%filter(15,10) = 0.69937986118348521724D-1
      DISC%DynRup%filter(15,11) = 0.37313779850039497991D-2
      DISC%DynRup%filter(15,12) = -0.21605225555243248168D-1
      DISC%DynRup%filter(15,13) = -0.11557175114670196001D-1
      DISC%DynRup%filter(15,14) = 0.16219450636794285983D0
      DISC%DynRup%filter(15,15) = 0.35414056088734394884D0
      DISC%DynRup%filter(15,16) = 0.25601252801996371683D0
      DISC%DynRup%filter(15,17) = 0.16374433275672082585D-1
      DISC%DynRup%filter(15,18) = -0.60622291106808189574D-1
      DISC%DynRup%filter(15,19) = -0.36942571734671261526D-1
      DISC%DynRup%filter(15,20) = 0.10391063417836409684D0
      DISC%DynRup%filter(15,21) = 0.21230605583637152748D0
      DISC%DynRup%filter(15,22) = 0.73803394035985101408D-1
      DISC%DynRup%filter(15,23) = -0.28950991598728067158D-1
      DISC%DynRup%filter(15,24) = -0.17879826990452455564D-2
      DISC%DynRup%filter(15,25) = -0.20015361105372665830D-1
      DISC%DynRup%filter(15,26) = -0.39344116166935647061D-1
      DISC%DynRup%filter(15,27) = -0.85766504974335971340D-1
      DISC%DynRup%filter(15,28) = -0.11645078694472243773D0
      DISC%DynRup%filter(15,29) = -0.25616464734255283323D-1
      DISC%DynRup%filter(15,30) = 0.49542489220350506558D-1
      DISC%DynRup%filter(15,31) = 0.17664040816602672609D-1
      DISC%DynRup%filter(15,32) = 0.20492000597443421138D-2
      DISC%DynRup%filter(15,33) = 0.19687920643503144107D-1
      DISC%DynRup%filter(15,34) = 0.54178549857021562623D-1
      DISC%DynRup%filter(15,35) = 0.12436675891705013668D-1
      DISC%DynRup%filter(15,36) = -0.23574094881996229436D-1
      DISC%DynRup%filter(16,1) = -0.44617604121467741006D-2
      DISC%DynRup%filter(16,2) = -0.80257776638945756520D-2
      DISC%DynRup%filter(16,3) = -0.89308048876615674280D-2
      DISC%DynRup%filter(16,4) = -0.92071126599200610350D-2
      DISC%DynRup%filter(16,5) = -0.86618346834658673497D-2
      DISC%DynRup%filter(16,6) = -0.49225720395378198497D-2
      DISC%DynRup%filter(16,7) = -0.21605225555243248097D-1
      DISC%DynRup%filter(16,8) = 0.37313779850039498760D-2
      DISC%DynRup%filter(16,9) = 0.69937986118348521784D-1
      DISC%DynRup%filter(16,10) = 0.85422762656326429973D-1
      DISC%DynRup%filter(16,11) = 0.30925711463459350441D-1
      DISC%DynRup%filter(16,12) = -0.78753983502977199533D-2
      DISC%DynRup%filter(16,13) = -0.60622291106808189558D-1
      DISC%DynRup%filter(16,14) = 0.16374433275672082673D-1
      DISC%DynRup%filter(16,15) = 0.25601252801996371683D0
      DISC%DynRup%filter(16,16) = 0.35414056088734394891D0
      DISC%DynRup%filter(16,17) = 0.16219450636794286003D0
      DISC%DynRup%filter(16,18) = -0.11557175114670195934D-1
      DISC%DynRup%filter(16,19) = -0.17879826990452455048D-2
      DISC%DynRup%filter(16,20) = -0.28950991598728067087D-1
      DISC%DynRup%filter(16,21) = 0.73803394035985101424D-1
      DISC%DynRup%filter(16,22) = 0.21230605583637152763D0
      DISC%DynRup%filter(16,23) = 0.10391063417836409694D0
      DISC%DynRup%filter(16,24) = -0.36942571734671261468D-1
      DISC%DynRup%filter(16,25) = 0.49542489220350506753D-1
      DISC%DynRup%filter(16,26) = -0.25616464734255283044D-1
      DISC%DynRup%filter(16,27) = -0.11645078694472243772D0
      DISC%DynRup%filter(16,28) = -0.85766504974335971438D-1
      DISC%DynRup%filter(16,29) = -0.39344116166935647303D-1
      DISC%DynRup%filter(16,30) = -0.20015361105372665954D-1
      DISC%DynRup%filter(16,31) = -0.23574094881996229505D-1
      DISC%DynRup%filter(16,32) = 0.12436675891705013525D-1
      DISC%DynRup%filter(16,33) = 0.54178549857021562603D-1
      DISC%DynRup%filter(16,34) = 0.19687920643503144115D-1
      DISC%DynRup%filter(16,35) = 0.20492000597443421390D-2
      DISC%DynRup%filter(16,36) = 0.17664040816602672650D-1
      DISC%DynRup%filter(17,1) = -0.28296384971352656729D-2
      DISC%DynRup%filter(17,2) = -0.69355961340712808847D-2
      DISC%DynRup%filter(17,3) = -0.10409570988918571482D-1
      DISC%DynRup%filter(17,4) = -0.11234547829233182851D-1
      DISC%DynRup%filter(17,5) = -0.87360995845765705025D-2
      DISC%DynRup%filter(17,6) = -0.40644093126917938653D-2
      DISC%DynRup%filter(17,7) = 0.27196287459698408062D-2
      DISC%DynRup%filter(17,8) = -0.18921793089315761590D-2
      DISC%DynRup%filter(17,9) = 0.48396611079976173942D-2
      DISC%DynRup%filter(17,10) = 0.40111177052651734920D-1
      DISC%DynRup%filter(17,11) = 0.68571852037770439386D-1
      DISC%DynRup%filter(17,12) = 0.46187074682139226912D-1
      DISC%DynRup%filter(17,13) = 0.47272601403983624560D-3
      DISC%DynRup%filter(17,14) = -0.30896646161334001835D-1
      DISC%DynRup%filter(17,15) = 0.21237920202203348577D-1
      DISC%DynRup%filter(17,16) = 0.21036905067096118834D0
      DISC%DynRup%filter(17,17) = 0.32014999902461882000D0
      DISC%DynRup%filter(17,18) = 0.19520951257895503196D0
      DISC%DynRup%filter(17,19) = 0.17679914432891125724D-1
      DISC%DynRup%filter(17,20) = -0.47852086498362494229D-1
      DISC%DynRup%filter(17,21) = -0.37549931591339888912D-1
      DISC%DynRup%filter(17,22) = 0.13477387092957948127D0
      DISC%DynRup%filter(17,23) = 0.18486811561989458989D0
      DISC%DynRup%filter(17,24) = 0.70418655125613338933D-1
      DISC%DynRup%filter(17,25) = 0.23116385210533820360D-1
      DISC%DynRup%filter(17,26) = 0.11999354830556001254D-1
      DISC%DynRup%filter(17,27) = -0.33224993178662543603D-1
      DISC%DynRup%filter(17,28) = -0.51029992031605233664D-1
      DISC%DynRup%filter(17,29) = -0.91899732666443250556D-1
      DISC%DynRup%filter(17,30) = -0.96611766869650293276D-1
      DISC%DynRup%filter(17,31) = -0.20030847083922715695D-1
      DISC%DynRup%filter(17,32) = 0.11970899353252597350D-1
      DISC%DynRup%filter(17,33) = 0.16130581481627257810D-1
      DISC%DynRup%filter(17,34) = 0.26578475489506295529D-2
      DISC%DynRup%filter(17,35) = 0.35770179977729635973D-1
      DISC%DynRup%filter(17,36) = 0.35943631108943101026D-1
      DISC%DynRup%filter(18,1) = -0.11945278462447170893D-2
      DISC%DynRup%filter(18,2) = -0.59584290676013193831D-2
      DISC%DynRup%filter(18,3) = -0.12185764221894136911D-1
      DISC%DynRup%filter(18,4) = -0.13444312714728335970D-1
      DISC%DynRup%filter(18,5) = -0.85585118437885783406D-2
      DISC%DynRup%filter(18,6) = -0.28683166523695777579D-2
      DISC%DynRup%filter(18,7) = 0.34038107697218555941D-1
      DISC%DynRup%filter(18,8) = 0.57267792297413694959D-2
      DISC%DynRup%filter(18,9) = -0.59007243835928168696D-1
      DISC%DynRup%filter(18,10) = -0.21508942342344918649D-1
      DISC%DynRup%filter(18,11) = 0.97257090830559119992D-1
      DISC%DynRup%filter(18,12) = 0.10403142273835132511D0
      DISC%DynRup%filter(18,13) = 0.66379927208593731289D-1
      DISC%DynRup%filter(18,14) = 0.99542907191781220824D-3
      DISC%DynRup%filter(18,15) = -0.16556894090670261638D0
      DISC%DynRup%filter(18,16) = -0.31564449457014120485D-1
      DISC%DynRup%filter(18,17) = 0.41105676050149143005D0
      DISC%DynRup%filter(18,18) = 0.43524383591115798653D0
      DISC%DynRup%filter(18,19) = -0.36249712121918345984D-1
      DISC%DynRup%filter(18,20) = 0.37228966235896802968D-1
      DISC%DynRup%filter(18,21) = -0.48832598774410677886D-2
      DISC%DynRup%filter(18,22) = -0.10089593060253854070D0
      DISC%DynRup%filter(18,23) = 0.14828203744988478058D0
      DISC%DynRup%filter(18,24) = 0.27885643693439252328D0
      DISC%DynRup%filter(18,25) = -0.55874582497571185201D-1
      DISC%DynRup%filter(18,26) = 0.48676656652699504557D-1
      DISC%DynRup%filter(18,27) = 0.13530827225983809850D0
      DISC%DynRup%filter(18,28) = -0.54665075825707132955D-1
      DISC%DynRup%filter(18,29) = -0.20343742162514419879D0
      DISC%DynRup%filter(18,30) = -0.10765859366938658532D0
      DISC%DynRup%filter(18,31) = 0.39800337769641430811D-1
      DISC%DynRup%filter(18,32) = -0.42179374373919118577D-1
      DISC%DynRup%filter(18,33) = -0.64384533332292553160D-1
      DISC%DynRup%filter(18,34) = 0.48243253046719999646D-1
      DISC%DynRup%filter(18,35) = 0.75687257086546390659D-1
      DISC%DynRup%filter(18,36) = 0.25275352189884356787D-1
      DISC%DynRup%filter(19,1) = -0.21173572948963040120D-2
      DISC%DynRup%filter(19,2) = -0.14863653459932543269D-2
      DISC%DynRup%filter(19,3) = 0.44109552989380144413D-2
      DISC%DynRup%filter(19,4) = 0.11853945995746162642D-1
      DISC%DynRup%filter(19,5) = 0.14259307033705388251D-1
      DISC%DynRup%filter(19,6) = 0.83093021879715510094D-2
      DISC%DynRup%filter(19,7) = -0.22863794823487263319D-1
      DISC%DynRup%filter(19,8) = -0.59646089835132517317D-1
      DISC%DynRup%filter(19,9) = -0.62323137733684479290D-1
      DISC%DynRup%filter(19,10) = -0.16179317264122364339D-1
      DISC%DynRup%filter(19,11) = 0.18776358410224887463D-1
      DISC%DynRup%filter(19,12) = 0.14307693160056797255D-1
      DISC%DynRup%filter(19,13) = 0.19541017092851854964D0
      DISC%DynRup%filter(19,14) = 0.10390944746428145168D0
      DISC%DynRup%filter(19,15) = -0.70703374330471962204D-1
      DISC%DynRup%filter(19,16) = -0.34219710250534489552D-2
      DISC%DynRup%filter(19,17) = 0.26088401385406191331D-1
      DISC%DynRup%filter(19,18) = -0.25402183717639008698D-1
      DISC%DynRup%filter(19,19) = 0.52480866411765823214D0
      DISC%DynRup%filter(19,20) = 0.28454758783745889304D0
      DISC%DynRup%filter(19,21) = -0.13939269254802411990D0
      DISC%DynRup%filter(19,22) = 0.43087364640303242242D-1
      DISC%DynRup%filter(19,23) = 0.86101744442090686655D-1
      DISC%DynRup%filter(19,24) = -0.56016583365258560429D-1
      DISC%DynRup%filter(19,25) = 0.27405939406932164936D0
      DISC%DynRup%filter(19,26) = 0.39904791356480491370D-1
      DISC%DynRup%filter(19,27) = -0.15295085716424642227D0
      DISC%DynRup%filter(19,28) = 0.94278977283770631369D-2
      DISC%DynRup%filter(19,29) = 0.27196128916243868072D-1
      DISC%DynRup%filter(19,30) = -0.82591284197194606640D-2
      DISC%DynRup%filter(19,31) = -0.76129198764368349368D-1
      DISC%DynRup%filter(19,32) = -0.63705854233972002600D-1
      DISC%DynRup%filter(19,33) = 0.85750384064317256150D-1
      DISC%DynRup%filter(19,34) = 0.23784044912797163836D-1
      DISC%DynRup%filter(19,35) = -0.62207479053885976242D-1
      DISC%DynRup%filter(19,36) = 0.26811800970057952442D-1
      DISC%DynRup%filter(20,1) = -0.70587004663688111868D-3
      DISC%DynRup%filter(20,2) = 0.81641787472217521162D-3
      DISC%DynRup%filter(20,3) = 0.56412040503086694949D-2
      DISC%DynRup%filter(20,4) = 0.10800581462471442747D-1
      DISC%DynRup%filter(20,5) = 0.11905756211563501136D-1
      DISC%DynRup%filter(20,6) = 0.67716983230426505504D-2
      DISC%DynRup%filter(20,7) = -0.28325733190111435131D-1
      DISC%DynRup%filter(20,8) = -0.45374385577356714419D-1
      DISC%DynRup%filter(20,9) = -0.40476550041473295987D-1
      DISC%DynRup%filter(20,10) = -0.23085226370777535044D-1
      DISC%DynRup%filter(20,11) = 0.41677582969688498319D-3
      DISC%DynRup%filter(20,12) = 0.89168312638771560763D-2
      DISC%DynRup%filter(20,13) = 0.49346257113261419294D-1
      DISC%DynRup%filter(20,14) = 0.12954734152690908086D0
      DISC%DynRup%filter(20,15) = 0.94443526011355531205D-1
      DISC%DynRup%filter(20,16) = -0.26313319610923153178D-1
      DISC%DynRup%filter(20,17) = -0.33532610918826520275D-1
      DISC%DynRup%filter(20,18) = 0.12389296583265415012D-1
      DISC%DynRup%filter(20,19) = 0.13513071980497497650D0
      DISC%DynRup%filter(20,20) = 0.37622629793840107229D0
      DISC%DynRup%filter(20,21) = 0.27989672989626038351D0
      DISC%DynRup%filter(20,22) = -0.38023524107024085203D-1
      DISC%DynRup%filter(20,23) = -0.50983575536805291638D-1
      DISC%DynRup%filter(20,24) = 0.40889437128421318422D-1
      DISC%DynRup%filter(20,25) = 0.18950655040340081968D-1
      DISC%DynRup%filter(20,26) = 0.17502918052298709194D0
      DISC%DynRup%filter(20,27) = 0.74944527189351450008D-1
      DISC%DynRup%filter(20,28) = -0.83350522096298500113D-1
      DISC%DynRup%filter(20,29) = -0.91109669384735631507D-2
      DISC%DynRup%filter(20,30) = 0.12915352768550628602D-1
      DISC%DynRup%filter(20,31) = -0.30253702039269859028D-1
      DISC%DynRup%filter(20,32) = -0.10897347256154196903D-1
      DISC%DynRup%filter(20,33) = -0.54113828344145422162D-1
      DISC%DynRup%filter(20,34) = 0.37458178976918125988D-2
      DISC%DynRup%filter(20,35) = 0.55364886198483933706D-1
      DISC%DynRup%filter(20,36) = -0.29542128561660223893D-1
      DISC%DynRup%filter(21,1) = 0.16150505929854226008D-2
      DISC%DynRup%filter(21,2) = 0.43493674727456238632D-2
      DISC%DynRup%filter(21,3) = 0.73867563790825963558D-2
      DISC%DynRup%filter(21,4) = 0.92111000255270865140D-2
      DISC%DynRup%filter(21,5) = 0.83272466800845975211D-2
      DISC%DynRup%filter(21,6) = 0.43402667250462311917D-2
      DISC%DynRup%filter(21,7) = -0.22819324552604982758D-1
      DISC%DynRup%filter(21,8) = -0.31207413982784705099D-1
      DISC%DynRup%filter(21,9) = -0.25794540556935979625D-1
      DISC%DynRup%filter(21,10) = -0.24384321813861053405D-1
      DISC%DynRup%filter(21,11) = -0.17798706053282259896D-1
      DISC%DynRup%filter(21,12) = -0.59239811266759588716D-2
      DISC%DynRup%filter(21,13) = -0.25887708874762666799D-1
      DISC%DynRup%filter(21,14) = 0.72815944323573782980D-1
      DISC%DynRup%filter(21,15) = 0.14877462796351260619D0
      DISC%DynRup%filter(21,16) = 0.51718131387692338913D-1
      DISC%DynRup%filter(21,17) = -0.20287565455010715998D-1
      DISC%DynRup%filter(21,18) = -0.12529386399635719979D-2
      DISC%DynRup%filter(21,19) = -0.51037980550771541460D-1
      DISC%DynRup%filter(21,20) = 0.21580033657390066261D0
      DISC%DynRup%filter(21,21) = 0.42432032316334500895D0
      DISC%DynRup%filter(21,22) = 0.16759330083426631452D0
      DISC%DynRup%filter(21,23) = -0.29316131356957501581D-1
      DISC%DynRup%filter(21,24) = 0.15776236460445431446D-1
      DISC%DynRup%filter(21,25) = -0.56002238930015551191D-1
      DISC%DynRup%filter(21,26) = 0.57782219170006729138D-1
      DISC%DynRup%filter(21,27) = 0.19616737615976376869D0
      DISC%DynRup%filter(21,28) = 0.52242129531908858621D-1
      DISC%DynRup%filter(21,29) = -0.64263239976609246014D-1
      DISC%DynRup%filter(21,30) = 0.34519805314026300846D-2
      DISC%DynRup%filter(21,31) = 0.31397100910350776926D-1
      DISC%DynRup%filter(21,32) = -0.41721753499217609172D-1
      DISC%DynRup%filter(21,33) = -0.60081860120312681494D-1
      DISC%DynRup%filter(21,34) = -0.68862309155885235802D-2
      DISC%DynRup%filter(21,35) = 0.28880250346834590499D-2
      DISC%DynRup%filter(21,36) = 0.87084164850306226862D-2
      DISC%DynRup%filter(22,1) = 0.43402667250462312040D-2
      DISC%DynRup%filter(22,2) = 0.83272466800845975369D-2
      DISC%DynRup%filter(22,3) = 0.92111000255270865232D-2
      DISC%DynRup%filter(22,4) = 0.73867563790825963519D-2
      DISC%DynRup%filter(22,5) = 0.43493674727456238334D-2
      DISC%DynRup%filter(22,6) = 0.16150505929854225851D-2
      DISC%DynRup%filter(22,7) = -0.59239811266759587827D-2
      DISC%DynRup%filter(22,8) = -0.17798706053282259806D-1
      DISC%DynRup%filter(22,9) = -0.24384321813861053409D-1
      DISC%DynRup%filter(22,10) = -0.25794540556935979698D-1
      DISC%DynRup%filter(22,11) = -0.31207413982784705231D-1
      DISC%DynRup%filter(22,12) = -0.22819324552604982824D-1
      DISC%DynRup%filter(22,13) = -0.12529386399635720513D-2
      DISC%DynRup%filter(22,14) = -0.20287565455010716000D-1
      DISC%DynRup%filter(22,15) = 0.51718131387692338897D-1
      DISC%DynRup%filter(22,16) = 0.14877462796351260629D0
      DISC%DynRup%filter(22,17) = 0.72815944323573783165D-1
      DISC%DynRup%filter(22,18) = -0.25887708874762666707D-1
      DISC%DynRup%filter(22,19) = 0.15776236460445431458D-1
      DISC%DynRup%filter(22,20) = -0.29316131356957501614D-1
      DISC%DynRup%filter(22,21) = 0.16759330083426631452D0
      DISC%DynRup%filter(22,22) = 0.42432032316334500905D0
      DISC%DynRup%filter(22,23) = 0.21580033657390066284D0
      DISC%DynRup%filter(22,24) = -0.51037980550771541349D-1
      DISC%DynRup%filter(22,25) = 0.34519805314026302425D-2
      DISC%DynRup%filter(22,26) = -0.64263239976609245770D-1
      DISC%DynRup%filter(22,27) = 0.52242129531908858583D-1
      DISC%DynRup%filter(22,28) = 0.19616737615976376857D0
      DISC%DynRup%filter(22,29) = 0.57782219170006729169D-1
      DISC%DynRup%filter(22,30) = -0.56002238930015551227D-1
      DISC%DynRup%filter(22,31) = 0.87084164850306226038D-2
      DISC%DynRup%filter(22,32) = 0.28880250346834589425D-2
      DISC%DynRup%filter(22,33) = -0.68862309155885236599D-2
      DISC%DynRup%filter(22,34) = -0.60081860120312681487D-1
      DISC%DynRup%filter(22,35) = -0.41721753499217609180D-1
      DISC%DynRup%filter(22,36) = 0.31397100910350776902D-1
      DISC%DynRup%filter(23,1) = 0.67716983230426505696D-2
      DISC%DynRup%filter(23,2) = 0.11905756211563501158D-1
      DISC%DynRup%filter(23,3) = 0.10800581462471442781D-1
      DISC%DynRup%filter(23,4) = 0.56412040503086695131D-2
      DISC%DynRup%filter(23,5) = 0.81641787472217519048D-3
      DISC%DynRup%filter(23,6) = -0.70587004663688113012D-3
      DISC%DynRup%filter(23,7) = 0.89168312638771561750D-2
      DISC%DynRup%filter(23,8) = 0.41677582969688502132D-3
      DISC%DynRup%filter(23,9) = -0.23085226370777535176D-1
      DISC%DynRup%filter(23,10) = -0.40476550041473296122D-1
      DISC%DynRup%filter(23,11) = -0.45374385577356714567D-1
      DISC%DynRup%filter(23,12) = -0.28325733190111435166D-1
      DISC%DynRup%filter(23,13) = 0.12389296583265415004D-1
      DISC%DynRup%filter(23,14) = -0.33532610918826520234D-1
      DISC%DynRup%filter(23,15) = -0.26313319610923153242D-1
      DISC%DynRup%filter(23,16) = 0.94443526011355531291D-1
      DISC%DynRup%filter(23,17) = 0.12954734152690908114D0
      DISC%DynRup%filter(23,18) = 0.49346257113261419429D-1
      DISC%DynRup%filter(23,19) = 0.40889437128421318521D-1
      DISC%DynRup%filter(23,20) = -0.50983575536805291638D-1
      DISC%DynRup%filter(23,21) = -0.38023524107024085163D-1
      DISC%DynRup%filter(23,22) = 0.27989672989626038382D0
      DISC%DynRup%filter(23,23) = 0.37622629793840107265D0
      DISC%DynRup%filter(23,24) = 0.13513071980497497656D0
      DISC%DynRup%filter(23,25) = 0.12915352768550628789D-1
      DISC%DynRup%filter(23,26) = -0.91109669384735629143D-2
      DISC%DynRup%filter(23,27) = -0.83350522096298500148D-1
      DISC%DynRup%filter(23,28) = 0.74944527189351449935D-1
      DISC%DynRup%filter(23,29) = 0.17502918052298709190D0
      DISC%DynRup%filter(23,30) = 0.18950655040340081833D-1
      DISC%DynRup%filter(23,31) = -0.29542128561660224008D-1
      DISC%DynRup%filter(23,32) = 0.55364886198483933560D-1
      DISC%DynRup%filter(23,33) = 0.37458178976918125154D-2
      DISC%DynRup%filter(23,34) = -0.54113828344145422240D-1
      DISC%DynRup%filter(23,35) = -0.10897347256154196879D-1
      DISC%DynRup%filter(23,36) = -0.30253702039269858931D-1
      DISC%DynRup%filter(24,1) = 0.83093021879715510304D-2
      DISC%DynRup%filter(24,2) = 0.14259307033705388272D-1
      DISC%DynRup%filter(24,3) = 0.11853945995746162682D-1
      DISC%DynRup%filter(24,4) = 0.44109552989380144605D-2
      DISC%DynRup%filter(24,5) = -0.14863653459932543498D-2
      DISC%DynRup%filter(24,6) = -0.21173572948963040265D-2
      DISC%DynRup%filter(24,7) = 0.14307693160056797359D-1
      DISC%DynRup%filter(24,8) = 0.18776358410224887529D-1
      DISC%DynRup%filter(24,9) = -0.16179317264122364458D-1
      DISC%DynRup%filter(24,10) = -0.62323137733684479380D-1
      DISC%DynRup%filter(24,11) = -0.59646089835132517445D-1
      DISC%DynRup%filter(24,12) = -0.22863794823487263342D-1
      DISC%DynRup%filter(24,13) = -0.25402183717639008733D-1
      DISC%DynRup%filter(24,14) = 0.26088401385406191245D-1
      DISC%DynRup%filter(24,15) = -0.34219710250534490540D-2
      DISC%DynRup%filter(24,16) = -0.70703374330471962090D-1
      DISC%DynRup%filter(24,17) = 0.10390944746428145190D0
      DISC%DynRup%filter(24,18) = 0.19541017092851854969D0
      DISC%DynRup%filter(24,19) = -0.56016583365258560433D-1
      DISC%DynRup%filter(24,20) = 0.86101744442090686444D-1
      DISC%DynRup%filter(24,21) = 0.43087364640303242202D-1
      DISC%DynRup%filter(24,22) = -0.13939269254802411960D0
      DISC%DynRup%filter(24,23) = 0.28454758783745889319D0
      DISC%DynRup%filter(24,24) = 0.52480866411765823197D0
      DISC%DynRup%filter(24,25) = -0.82591284197194606822D-2
      DISC%DynRup%filter(24,26) = 0.27196128916243867886D-1
      DISC%DynRup%filter(24,27) = 0.94278977283770633420D-2
      DISC%DynRup%filter(24,28) = -0.15295085716424642203D0
      DISC%DynRup%filter(24,29) = 0.39904791356480491121D-1
      DISC%DynRup%filter(24,30) = 0.27405939406932164912D0
      DISC%DynRup%filter(24,31) = 0.26811800970057952455D-1
      DISC%DynRup%filter(24,32) = -0.62207479053885976177D-1
      DISC%DynRup%filter(24,33) = 0.23784044912797163914D-1
      DISC%DynRup%filter(24,34) = 0.85750384064317256025D-1
      DISC%DynRup%filter(24,35) = -0.63705854233972002534D-1
      DISC%DynRup%filter(24,36) = -0.76129198764368349089D-1
      DISC%DynRup%filter(25,1) = 0.28735159143547768185D-2
      DISC%DynRup%filter(25,2) = 0.37672259387559907928D-2
      DISC%DynRup%filter(25,3) = -0.91374072335582649341D-3
      DISC%DynRup%filter(25,4) = -0.92926180263219641094D-2
      DISC%DynRup%filter(25,5) = -0.14201754121916461751D-1
      DISC%DynRup%filter(25,6) = -0.92125814154940427232D-2
      DISC%DynRup%filter(25,7) = 0.39245262798201844907D-2
      DISC%DynRup%filter(25,8) = 0.31309378829373423632D-1
      DISC%DynRup%filter(25,9) = 0.65467250599151124760D-1
      DISC%DynRup%filter(25,10) = 0.47733002954458998524D-1
      DISC%DynRup%filter(25,11) = -0.15806738433887099108D-1
      DISC%DynRup%filter(25,12) = -0.34656375392525733145D-1
      DISC%DynRup%filter(25,13) = -0.78364505388654980916D-1
      DISC%DynRup%filter(25,14) = -0.14808174972223305140D0
      DISC%DynRup%filter(25,15) = -0.39790614786128181935D-1
      DISC%DynRup%filter(25,16) = 0.98490658936136682558D-1
      DISC%DynRup%filter(25,17) = 0.35431654757411740389D-1
      DISC%DynRup%filter(25,18) = -0.40671012614804785464D-1
      DISC%DynRup%filter(25,19) = 0.28467471655593505883D0
      DISC%DynRup%filter(25,20) = 0.41450449845760067559D-1
      DISC%DynRup%filter(25,21) = -0.15887520315835503990D0
      DISC%DynRup%filter(25,22) = 0.97930746824362183648D-2
      DISC%DynRup%filter(25,23) = 0.28249534437385817950D-1
      DISC%DynRup%filter(25,24) = -0.85790346646099366343D-2
      DISC%DynRup%filter(25,25) = 0.60059998121258159445D0
      DISC%DynRup%filter(25,26) = 0.28551291624184362033D0
      DISC%DynRup%filter(25,27) = -0.34733860275297734582D-1
      DISC%DynRup%filter(25,28) = 0.30139986668067570991D-1
      DISC%DynRup%filter(25,29) = -0.46491353612005898879D-1
      DISC%DynRup%filter(25,30) = 0.19941211192762176947D-1
      DISC%DynRup%filter(25,31) = 0.15468541431126348333D0
      DISC%DynRup%filter(25,32) = -0.10089914519510413440D0
      DISC%DynRup%filter(25,33) = -0.48425309420414399792D-1
      DISC%DynRup%filter(25,34) = 0.80730970476171568057D-1
      DISC%DynRup%filter(25,35) = -0.51323276482609233046D-1
      DISC%DynRup%filter(25,36) = 0.15543403600048403319D-1
      DISC%DynRup%filter(26,1) = 0.17890432902311675544D-2
      DISC%DynRup%filter(26,2) = 0.14533268247863452185D-2
      DISC%DynRup%filter(26,3) = -0.30815881716710943038D-2
      DISC%DynRup%filter(26,4) = -0.90748064704688631315D-2
      DISC%DynRup%filter(26,5) = -0.11321561273587122550D-1
      DISC%DynRup%filter(26,6) = -0.67443666332679602363D-2
      DISC%DynRup%filter(26,7) = 0.14868721713700324957D-1
      DISC%DynRup%filter(26,8) = 0.34934964542022827799D-1
      DISC%DynRup%filter(26,9) = 0.39117345469995966681D-1
      DISC%DynRup%filter(26,10) = 0.19525398351807503252D-1
      DISC%DynRup%filter(26,11) = -0.29688168018339204432D-2
      DISC%DynRup%filter(26,12) = -0.75065684393018032516D-2
      DISC%DynRup%filter(26,13) = -0.70323539138121147624D-1
      DISC%DynRup%filter(26,14) = -0.66893657536261199919D-1
      DISC%DynRup%filter(26,15) = -0.37144643536996837280D-1
      DISC%DynRup%filter(26,16) = -0.24184415458583994130D-1
      DISC%DynRup%filter(26,17) = 0.87343097678497353753D-2
      DISC%DynRup%filter(26,18) = 0.16826377083840866842D-1
      DISC%DynRup%filter(26,19) = 0.19684683207004173017D-1
      DISC%DynRup%filter(26,20) = 0.18180870071469128563D0
      DISC%DynRup%filter(26,21) = 0.77847402777409199027D-1
      DISC%DynRup%filter(26,22) = -0.86578992605344442112D-1
      DISC%DynRup%filter(26,23) = -0.94638680041174067456D-2
      DISC%DynRup%filter(26,24) = 0.13415611608909377503D-1
      DISC%DynRup%filter(26,25) = 0.13558915111034668147D0
      DISC%DynRup%filter(26,26) = 0.51313637418937058775D0
      DISC%DynRup%filter(26,27) = 0.24580281143772826799D0
      DISC%DynRup%filter(26,28) = -0.66019040237055266576D-1
      DISC%DynRup%filter(26,29) = 0.48538176863777574238D-1
      DISC%DynRup%filter(26,30) = -0.22078591936216515258D-1
      DISC%DynRup%filter(26,31) = -0.47916674400731639084D-1
      DISC%DynRup%filter(26,32) = 0.15000741484863905285D0
      DISC%DynRup%filter(26,33) = -0.10903861003931910787D-1
      DISC%DynRup%filter(26,34) = -0.91430832035197967452D-1
      DISC%DynRup%filter(26,35) = 0.74929266201437867281D-1
      DISC%DynRup%filter(26,36) = -0.24373256320859715155D-1
      DISC%DynRup%filter(27,1) = -0.33456187992792772872D-3
      DISC%DynRup%filter(27,2) = -0.23759040160106504161D-2
      DISC%DynRup%filter(27,3) = -0.58455352021157061965D-2
      DISC%DynRup%filter(27,4) = -0.80248284115962505635D-2
      DISC%DynRup%filter(27,5) = -0.69966744212982570538D-2
      DISC%DynRup%filter(27,6) = -0.34024485030287354266D-2
      DISC%DynRup%filter(27,7) = 0.23970526730738042043D-1
      DISC%DynRup%filter(27,8) = 0.30159467462986829401D-1
      DISC%DynRup%filter(27,9) = 0.10303758308480350285D-1
      DISC%DynRup%filter(27,10) = 0.10059980084532167531D-2
      DISC%DynRup%filter(27,11) = 0.15054079187067490430D-1
      DISC%DynRup%filter(27,12) = 0.17477215138664969961D-1
      DISC%DynRup%filter(27,13) = -0.14569146964845249749D-1
      DISC%DynRup%filter(27,14) = -0.28638514569900227102D-1
      DISC%DynRup%filter(27,15) = -0.62429291635254095380D-1
      DISC%DynRup%filter(27,16) = -0.84764327769941324212D-1
      DISC%DynRup%filter(27,17) = -0.18646180674350213910D-1
      DISC%DynRup%filter(27,18) = 0.36061892796018533402D-1
      DISC%DynRup%filter(27,19) = -0.58171410427432349887D-1
      DISC%DynRup%filter(27,20) = 0.60020335810980717981D-1
      DISC%DynRup%filter(27,21) = 0.20376565596462229618D0
      DISC%DynRup%filter(27,22) = 0.54265658242726807151D-1
      DISC%DynRup%filter(27,23) = -0.66752390252603003376D-1
      DISC%DynRup%filter(27,24) = 0.35856883602577184779D-2
      DISC%DynRup%filter(27,25) = -0.12717640019566366918D-1
      DISC%DynRup%filter(27,26) = 0.18951393057979956445D0
      DISC%DynRup%filter(27,27) = 0.47207779536413002275D0
      DISC%DynRup%filter(27,28) = 0.24595985083232387817D0
      DISC%DynRup%filter(27,29) = -0.50900670074719462095D-1
      DISC%DynRup%filter(27,30) = 0.11035614745983693279D-1
      DISC%DynRup%filter(27,31) = -0.17730699903889888432D-1
      DISC%DynRup%filter(27,32) = -0.84068751909879858012D-2
      DISC%DynRup%filter(27,33) = 0.10053944845327242728D0
      DISC%DynRup%filter(27,34) = 0.16844075484119729164D-1
      DISC%DynRup%filter(27,35) = -0.70493157721922772854D-1
      DISC%DynRup%filter(27,36) = 0.29559266168764178822D-1
      DISC%DynRup%filter(28,1) = -0.34024485030287354306D-2
      DISC%DynRup%filter(28,2) = -0.69966744212982570638D-2
      DISC%DynRup%filter(28,3) = -0.80248284115962505617D-2
      DISC%DynRup%filter(28,4) = -0.58455352021157061918D-2
      DISC%DynRup%filter(28,5) = -0.23759040160106504090D-2
      DISC%DynRup%filter(28,6) = -0.33456187992792772730D-3
      DISC%DynRup%filter(28,7) = 0.17477215138664969929D-1
      DISC%DynRup%filter(28,8) = 0.15054079187067490414D-1
      DISC%DynRup%filter(28,9) = 0.10059980084532167709D-2
      DISC%DynRup%filter(28,10) = 0.10303758308480350355D-1
      DISC%DynRup%filter(28,11) = 0.30159467462986829459D-1
      DISC%DynRup%filter(28,12) = 0.23970526730738042075D-1
      DISC%DynRup%filter(28,13) = 0.36061892796018533405D-1
      DISC%DynRup%filter(28,14) = -0.18646180674350214024D-1
      DISC%DynRup%filter(28,15) = -0.84764327769941324216D-1
      DISC%DynRup%filter(28,16) = -0.62429291635254095451D-1
      DISC%DynRup%filter(28,17) = -0.28638514569900227225D-1
      DISC%DynRup%filter(28,18) = -0.14569146964845249780D-1
      DISC%DynRup%filter(28,19) = 0.35856883602577184023D-2
      DISC%DynRup%filter(28,20) = -0.66752390252603003353D-1
      DISC%DynRup%filter(28,21) = 0.54265658242726807200D-1
      DISC%DynRup%filter(28,22) = 0.20376565596462229605D0
      DISC%DynRup%filter(28,23) = 0.60020335810980717925D-1
      DISC%DynRup%filter(28,24) = -0.58171410427432349794D-1
      DISC%DynRup%filter(28,25) = 0.11035614745983693210D-1
      DISC%DynRup%filter(28,26) = -0.50900670074719462095D-1
      DISC%DynRup%filter(28,27) = 0.24595985083232387816D0
      DISC%DynRup%filter(28,28) = 0.47207779536413002269D0
      DISC%DynRup%filter(28,29) = 0.18951393057979956463D0
      DISC%DynRup%filter(28,30) = -0.12717640019566366742D-1
      DISC%DynRup%filter(28,31) = 0.29559266168764178889D-1
      DISC%DynRup%filter(28,32) = -0.70493157721922772784D-1
      DISC%DynRup%filter(28,33) = 0.16844075484119729137D-1
      DISC%DynRup%filter(28,34) = 0.10053944845327242725D0
      DISC%DynRup%filter(28,35) = -0.84068751909879858970D-2
      DISC%DynRup%filter(28,36) = -0.17730699903889888501D-1
      DISC%DynRup%filter(29,1) = -0.67443666332679602218D-2
      DISC%DynRup%filter(29,2) = -0.11321561273587122527D-1
      DISC%DynRup%filter(29,3) = -0.90748064704688630912D-2
      DISC%DynRup%filter(29,4) = -0.30815881716710942848D-2
      DISC%DynRup%filter(29,5) = 0.14533268247863452311D-2
      DISC%DynRup%filter(29,6) = 0.17890432902311675529D-2
      DISC%DynRup%filter(29,7) = -0.75065684393018032188D-2
      DISC%DynRup%filter(29,8) = -0.29688168018339203584D-2
      DISC%DynRup%filter(29,9) = 0.19525398351807503316D-1
      DISC%DynRup%filter(29,10) = 0.39117345469995966737D-1
      DISC%DynRup%filter(29,11) = 0.34934964542022827834D-1
      DISC%DynRup%filter(29,12) = 0.14868721713700324963D-1
      DISC%DynRup%filter(29,13) = 0.16826377083840866968D-1
      DISC%DynRup%filter(29,14) = 0.87343097678497350324D-2
      DISC%DynRup%filter(29,15) = -0.24184415458583994388D-1
      DISC%DynRup%filter(29,16) = -0.37144643536996837509D-1
      DISC%DynRup%filter(29,17) = -0.66893657536261200206D-1
      DISC%DynRup%filter(29,18) = -0.70323539138121147659D-1
      DISC%DynRup%filter(29,19) = 0.13415611608909377596D-1
      DISC%DynRup%filter(29,20) = -0.94638680041174069912D-2
      DISC%DynRup%filter(29,21) = -0.86578992605344442439D-1
      DISC%DynRup%filter(29,22) = 0.77847402777409199067D-1
      DISC%DynRup%filter(29,23) = 0.18180870071469128559D0
      DISC%DynRup%filter(29,24) = 0.19684683207004172891D-1
      DISC%DynRup%filter(29,25) = -0.22078591936216515045D-1
      DISC%DynRup%filter(29,26) = 0.48538176863777574238D-1
      DISC%DynRup%filter(29,27) = -0.66019040237055266579D-1
      DISC%DynRup%filter(29,28) = 0.24580281143772826822D0
      DISC%DynRup%filter(29,29) = 0.51313637418937058758D0
      DISC%DynRup%filter(29,30) = 0.13558915111034668125D0
      DISC%DynRup%filter(29,31) = -0.24373256320859715163D-1
      DISC%DynRup%filter(29,32) = 0.74929266201437867110D-1
      DISC%DynRup%filter(29,33) = -0.91430832035197967442D-1
      DISC%DynRup%filter(29,34) = -0.10903861003931910710D-1
      DISC%DynRup%filter(29,35) = 0.15000741484863905289D0
      DISC%DynRup%filter(29,36) = -0.47916674400731638926D-1
      DISC%DynRup%filter(30,1) = -0.92125814154940427023D-2
      DISC%DynRup%filter(30,2) = -0.14201754121916461714D-1
      DISC%DynRup%filter(30,3) = -0.92926180263219640556D-2
      DISC%DynRup%filter(30,4) = -0.91374072335582647307D-3
      DISC%DynRup%filter(30,5) = 0.37672259387559908081D-2
      DISC%DynRup%filter(30,6) = 0.28735159143547768150D-2
      DISC%DynRup%filter(30,7) = -0.34656375392525733108D-1
      DISC%DynRup%filter(30,8) = -0.15806738433887098968D-1
      DISC%DynRup%filter(30,9) = 0.47733002954458998677D-1
      DISC%DynRup%filter(30,10) = 0.65467250599151124817D-1
      DISC%DynRup%filter(30,11) = 0.31309378829373423676D-1
      DISC%DynRup%filter(30,12) = 0.39245262798201844918D-2
      DISC%DynRup%filter(30,13) = -0.40671012614804785327D-1
      DISC%DynRup%filter(30,14) = 0.35431654757411739807D-1
      DISC%DynRup%filter(30,15) = 0.98490658936136682171D-1
      DISC%DynRup%filter(30,16) = -0.39790614786128182184D-1
      DISC%DynRup%filter(30,17) = -0.14808174972223305176D0
      DISC%DynRup%filter(30,18) = -0.78364505388654980980D-1
      DISC%DynRup%filter(30,19) = -0.85790346646099366202D-2
      DISC%DynRup%filter(30,20) = 0.28249534437385817531D-1
      DISC%DynRup%filter(30,21) = 0.97930746824362179159D-2
      DISC%DynRup%filter(30,22) = -0.15887520315835503999D0
      DISC%DynRup%filter(30,23) = 0.41450449845760067269D-1
      DISC%DynRup%filter(30,24) = 0.28467471655593505858D0
      DISC%DynRup%filter(30,25) = 0.19941211192762176943D-1
      DISC%DynRup%filter(30,26) = -0.46491353612005899324D-1
      DISC%DynRup%filter(30,27) = 0.30139986668067571178D-1
      DISC%DynRup%filter(30,28) = -0.34733860275297734104D-1
      DISC%DynRup%filter(30,29) = 0.28551291624184361987D0
      DISC%DynRup%filter(30,30) = 0.60059998121258159416D0
      DISC%DynRup%filter(30,31) = 0.15543403600048403433D-1
      DISC%DynRup%filter(30,32) = -0.51323276482609233061D-1
      DISC%DynRup%filter(30,33) = 0.80730970476171568248D-1
      DISC%DynRup%filter(30,34) = -0.48425309420414399704D-1
      DISC%DynRup%filter(30,35) = -0.10089914519510413439D0
      DISC%DynRup%filter(30,36) = 0.15468541431126348354D0
      DISC%DynRup%filter(31,1) = -0.21920676927230929096D-2
      DISC%DynRup%filter(31,2) = -0.31599048484796100412D-2
      DISC%DynRup%filter(31,3) = -0.12938252702934486222D-3
      DISC%DynRup%filter(31,4) = 0.60359940435644087828D-2
      DISC%DynRup%filter(31,5) = 0.10155074012212573228D-1
      DISC%DynRup%filter(31,6) = 0.68342294734224565078D-2
      DISC%DynRup%filter(31,7) = 0.57806282501462950790D-3
      DISC%DynRup%filter(31,8) = -0.15789032674142386355D-1
      DISC%DynRup%filter(31,9) = -0.45971402555153452068D-1
      DISC%DynRup%filter(31,10) = -0.42465098227528581015D-1
      DISC%DynRup%filter(31,11) = 0.92429771092665881613D-2
      DISC%DynRup%filter(31,12) = 0.30697992272900191591D-1
      DISC%DynRup%filter(31,13) = 0.34486017574539201535D-1
      DISC%DynRup%filter(31,14) = 0.10326867291289158034D0
      DISC%DynRup%filter(31,15) = 0.65823718693343353136D-1
      DISC%DynRup%filter(31,16) = -0.87847090372674444892D-1
      DISC%DynRup%filter(31,17) = -0.57550084169517498836D-1
      DISC%DynRup%filter(31,18) = 0.54304095843450664710D-1
      DISC%DynRup%filter(31,19) = -0.14822811905773394560D0
      DISC%DynRup%filter(31,20) = -0.12403912164234631101D0
      DISC%DynRup%filter(31,21) = 0.16696114427360991870D0
      DISC%DynRup%filter(31,22) = 0.46308962897671453687D-1
      DISC%DynRup%filter(31,23) = -0.12112169523839343332D0
      DISC%DynRup%filter(31,24) = 0.52204185658684982482D-1
      DISC%DynRup%filter(31,25) = 0.28995092583146593806D0
      DISC%DynRup%filter(31,26) = -0.18913095779059291923D0
      DISC%DynRup%filter(31,27) = -0.90771087646768294662D-1
      DISC%DynRup%filter(31,28) = 0.15132661173687805181D0
      DISC%DynRup%filter(31,29) = -0.96203197949176431952D-1
      DISC%DynRup%filter(31,30) = 0.29135418387524136306D-1
      DISC%DynRup%filter(31,31) = 0.81888008130588519358D0
      DISC%DynRup%filter(31,32) = 0.25501534696734011669D0
      DISC%DynRup%filter(31,33) = -0.16230024600529912255D0
      DISC%DynRup%filter(31,34) = 0.81174053365296697594D-1
      DISC%DynRup%filter(31,35) = -0.34292093188078363764D-1
      DISC%DynRup%filter(31,36) = 0.88070164006750962802D-2
      DISC%DynRup%filter(32,1) = -0.15006284886666967974D-2
      DISC%DynRup%filter(32,2) = -0.15257565465193808885D-2
      DISC%DynRup%filter(32,3) = 0.16407756329447726557D-2
      DISC%DynRup%filter(32,4) = 0.61214232411234699129D-2
      DISC%DynRup%filter(32,5) = 0.79855171518135014700D-2
      DISC%DynRup%filter(32,6) = 0.48226114702717243728D-2
      DISC%DynRup%filter(32,7) = -0.74981600318463736816D-2
      DISC%DynRup%filter(32,8) = -0.22653145674842769538D-1
      DISC%DynRup%filter(32,9) = -0.28703728713182951400D-1
      DISC%DynRup%filter(32,10) = -0.12880053492904447830D-1
      DISC%DynRup%filter(32,11) = 0.36391270321154385338D-2
      DISC%DynRup%filter(32,12) = 0.43894596310180937000D-2
      DISC%DynRup%filter(32,13) = 0.49041955372311554775D-1
      DISC%DynRup%filter(32,14) = 0.48805296404538784615D-1
      DISC%DynRup%filter(32,15) = 0.36264015866113523466D-2
      DISC%DynRup%filter(32,16) = 0.22008774092792547019D-1
      DISC%DynRup%filter(32,17) = 0.16333249973249843161D-1
      DISC%DynRup%filter(32,18) = -0.27330346947471225626D-1
      DISC%DynRup%filter(32,19) = -0.58905773613803103260D-1
      DISC%DynRup%filter(32,20) = -0.21217789136311250450D-1
      DISC%DynRup%filter(32,21) = -0.10536287154804548430D0
      DISC%DynRup%filter(32,22) = 0.72933322604142015577D-2
      DISC%DynRup%filter(32,23) = 0.10779875627546752893D0
      DISC%DynRup%filter(32,24) = -0.57520297346229227368D-1
      DISC%DynRup%filter(32,25) = -0.89817673935954167950D-1
      DISC%DynRup%filter(32,26) = 0.28118222400352450666D0
      DISC%DynRup%filter(32,27) = -0.20438802244572473751D-1
      DISC%DynRup%filter(32,28) = -0.17138302609967884901D0
      DISC%DynRup%filter(32,29) = 0.14045157524200589435D0
      DISC%DynRup%filter(32,30) = -0.45686584395994429070D-1
      DISC%DynRup%filter(32,31) = 0.12110595510195225558D0
      DISC%DynRup%filter(32,32) = 0.66380750893729051744D0
      DISC%DynRup%filter(32,33) = 0.28851607821684572701D0
      DISC%DynRup%filter(32,34) = -0.15432153254200960260D0
      DISC%DynRup%filter(32,35) = 0.64461352149421274974D-1
      DISC%DynRup%filter(32,36) = -0.16285203017680555203D-1
      DISC%DynRup%filter(33,1) = -0.47372805399094641562D-4
      DISC%DynRup%filter(33,2) = 0.12650377657608624588D-2
      DISC%DynRup%filter(33,3) = 0.38725401799501301999D-2
      DISC%DynRup%filter(33,4) = 0.55240700171168299558D-2
      DISC%DynRup%filter(33,5) = 0.47196163965022166328D-2
      DISC%DynRup%filter(33,6) = 0.22100509070364461071D-2
      DISC%DynRup%filter(33,7) = -0.16832213415299792679D-1
      DISC%DynRup%filter(33,8) = -0.22130570512655272360D-1
      DISC%DynRup%filter(33,9) = -0.46045416770278180923D-2
      DISC%DynRup%filter(33,10) = 0.53397401897116150273D-2
      DISC%DynRup%filter(33,11) = -0.99305193022041127255D-2
      DISC%DynRup%filter(33,12) = -0.15548396532167629315D-1
      DISC%DynRup%filter(33,13) = 0.24101045851402028776D-1
      DISC%DynRup%filter(33,14) = 0.27959550768346187303D-2
      DISC%DynRup%filter(33,15) = 0.26862453674917082805D-1
      DISC%DynRup%filter(33,16) = 0.73921914460209249730D-1
      DISC%DynRup%filter(33,17) = 0.16968761509160574066D-1
      DISC%DynRup%filter(33,18) = -0.32164800090490697775D-1
      DISC%DynRup%filter(33,19) = 0.61132039839428827488D-1
      DISC%DynRup%filter(33,20) = -0.81234758086984939549D-1
      DISC%DynRup%filter(33,21) = -0.11698298760096923886D0
      DISC%DynRup%filter(33,22) = -0.13407904885144406681D-1
      DISC%DynRup%filter(33,23) = 0.56231580737867737387D-2
      DISC%DynRup%filter(33,24) = 0.16955809551375648781D-1
      DISC%DynRup%filter(33,25) = -0.33235407977301521444D-1
      DISC%DynRup%filter(33,26) = -0.15758313450753417622D-1
      DISC%DynRup%filter(33,27) = 0.18845672225405551529D0
      DISC%DynRup%filter(33,28) = 0.31573469956048635961D-1
      DISC%DynRup%filter(33,29) = -0.13213628729807620532D0
      DISC%DynRup%filter(33,30) = 0.55407529085357474170D-1
      DISC%DynRup%filter(33,31) = -0.59425473800572702934D-1
      DISC%DynRup%filter(33,32) = 0.22244585283026923780D0
      DISC%DynRup%filter(33,33) = 0.62435731992100222211D0
      DISC%DynRup%filter(33,34) = 0.26916683934081425512D0
      DISC%DynRup%filter(33,35) = -0.11898187833601689667D0
      DISC%DynRup%filter(33,36) = 0.29721498890323502492D-1
      DISC%DynRup%filter(34,1) = 0.22100509070364460994D-2
      DISC%DynRup%filter(34,2) = 0.47196163965022166237D-2
      DISC%DynRup%filter(34,3) = 0.55240700171168299374D-2
      DISC%DynRup%filter(34,4) = 0.38725401799501301930D-2
      DISC%DynRup%filter(34,5) = 0.12650377657608624609D-2
      DISC%DynRup%filter(34,6) = -0.47372805399094637587D-4
      DISC%DynRup%filter(34,7) = -0.15548396532167629333D-1
      DISC%DynRup%filter(34,8) = -0.99305193022041127524D-2
      DISC%DynRup%filter(34,9) = 0.53397401897116150400D-2
      DISC%DynRup%filter(34,10) = -0.46045416770278180714D-2
      DISC%DynRup%filter(34,11) = -0.22130570512655272323D-1
      DISC%DynRup%filter(34,12) = -0.16832213415299792665D-1
      DISC%DynRup%filter(34,13) = -0.32164800090490697816D-1
      DISC%DynRup%filter(34,14) = 0.16968761509160574135D-1
      DISC%DynRup%filter(34,15) = 0.73921914460209249760D-1
      DISC%DynRup%filter(34,16) = 0.26862453674917082814D-1
      DISC%DynRup%filter(34,17) = 0.27959550768346187203D-2
      DISC%DynRup%filter(34,18) = 0.24101045851402028742D-1
      DISC%DynRup%filter(34,19) = 0.16955809551375648727D-1
      DISC%DynRup%filter(34,20) = 0.56231580737867738613D-2
      DISC%DynRup%filter(34,21) = -0.13407904885144406543D-1
      DISC%DynRup%filter(34,22) = -0.11698298760096923884D0
      DISC%DynRup%filter(34,23) = -0.81234758086984939676D-1
      DISC%DynRup%filter(34,24) = 0.61132039839428827400D-1
      DISC%DynRup%filter(34,25) = 0.55407529085357474044D-1
      DISC%DynRup%filter(34,26) = -0.13213628729807620533D0
      DISC%DynRup%filter(34,27) = 0.31573469956048636015D-1
      DISC%DynRup%filter(34,28) = 0.18845672225405551523D0
      DISC%DynRup%filter(34,29) = -0.15758313450753417508D-1
      DISC%DynRup%filter(34,30) = -0.33235407977301521387D-1
      DISC%DynRup%filter(34,31) = 0.29721498890323502563D-1
      DISC%DynRup%filter(34,32) = -0.11898187833601689655D0
      DISC%DynRup%filter(34,33) = 0.26916683934081425512D0
      DISC%DynRup%filter(34,34) = 0.62435731992100222202D0
      DISC%DynRup%filter(34,35) = 0.22244585283026923789D0
      DISC%DynRup%filter(34,36) = -0.59425473800572702903D-1
      DISC%DynRup%filter(35,1) = 0.48226114702717243588D-2
      DISC%DynRup%filter(35,2) = 0.79855171518135014572D-2
      DISC%DynRup%filter(35,3) = 0.61214232411234698733D-2
      DISC%DynRup%filter(35,4) = 0.16407756329447726401D-2
      DISC%DynRup%filter(35,5) = -0.15257565465193808825D-2
      DISC%DynRup%filter(35,6) = -0.15006284886666967867D-2
      DISC%DynRup%filter(35,7) = 0.43894596310180936465D-2
      DISC%DynRup%filter(35,8) = 0.36391270321154384521D-2
      DISC%DynRup%filter(35,9) = -0.12880053492904447829D-1
      DISC%DynRup%filter(35,10) = -0.28703728713182951388D-1
      DISC%DynRup%filter(35,11) = -0.22653145674842769532D-1
      DISC%DynRup%filter(35,12) = -0.74981600318463737034D-2
      DISC%DynRup%filter(35,13) = -0.27330346947471225777D-1
      DISC%DynRup%filter(35,14) = 0.16333249973249843526D-1
      DISC%DynRup%filter(35,15) = 0.22008774092792547274D-1
      DISC%DynRup%filter(35,16) = 0.36264015866113523935D-2
      DISC%DynRup%filter(35,17) = 0.48805296404538784668D-1
      DISC%DynRup%filter(35,18) = 0.49041955372311554737D-1
      DISC%DynRup%filter(35,19) = -0.57520297346229227432D-1
      DISC%DynRup%filter(35,20) = 0.10779875627546752922D0
      DISC%DynRup%filter(35,21) = 0.72933322604142018456D-2
      DISC%DynRup%filter(35,22) = -0.10536287154804548432D0
      DISC%DynRup%filter(35,23) = -0.21217789136311250399D-1
      DISC%DynRup%filter(35,24) = -0.58905773613803103189D-1
      DISC%DynRup%filter(35,25) = -0.45686584395994429058D-1
      DISC%DynRup%filter(35,26) = 0.14045157524200589467D0
      DISC%DynRup%filter(35,27) = -0.17138302609967884917D0
      DISC%DynRup%filter(35,28) = -0.20438802244572473995D-1
      DISC%DynRup%filter(35,29) = 0.28118222400352450673D0
      DISC%DynRup%filter(35,30) = -0.89817673935954167968D-1
      DISC%DynRup%filter(35,31) = -0.16285203017680555271D-1
      DISC%DynRup%filter(35,32) = 0.64461352149421274981D-1
      DISC%DynRup%filter(35,33) = -0.15432153254200960276D0
      DISC%DynRup%filter(35,34) = 0.28851607821684572712D0
      DISC%DynRup%filter(35,35) = 0.66380750893729051773D0
      DISC%DynRup%filter(35,36) = 0.12110595510195225543D0
      DISC%DynRup%filter(36,1) = 0.68342294734224564976D-2
      DISC%DynRup%filter(36,2) = 0.10155074012212573232D-1
      DISC%DynRup%filter(36,3) = 0.60359940435644087289D-2
      DISC%DynRup%filter(36,4) = -0.12938252702934489145D-3
      DISC%DynRup%filter(36,5) = -0.31599048484796100486D-2
      DISC%DynRup%filter(36,6) = -0.21920676927230929003D-2
      DISC%DynRup%filter(36,7) = 0.30697992272900191570D-1
      DISC%DynRup%filter(36,8) = 0.92429771092665880931D-2
      DISC%DynRup%filter(36,9) = -0.42465098227528581036D-1
      DISC%DynRup%filter(36,10) = -0.45971402555153452253D-1
      DISC%DynRup%filter(36,11) = -0.15789032674142386496D-1
      DISC%DynRup%filter(36,12) = 0.57806282501462942431D-3
      DISC%DynRup%filter(36,13) = 0.54304095843450664568D-1
      DISC%DynRup%filter(36,14) = -0.57550084169517498236D-1
      DISC%DynRup%filter(36,15) = -0.87847090372674444629D-1
      DISC%DynRup%filter(36,16) = 0.65823718693343353293D-1
      DISC%DynRup%filter(36,17) = 0.10326867291289158081D0
      DISC%DynRup%filter(36,18) = 0.34486017574539201679D-1
      DISC%DynRup%filter(36,19) = 0.52204185658684982449D-1
      DISC%DynRup%filter(36,20) = -0.12112169523839343284D0
      DISC%DynRup%filter(36,21) = 0.46308962897671454144D-1
      DISC%DynRup%filter(36,22) = 0.16696114427360991857D0
      DISC%DynRup%filter(36,23) = -0.12403912164234631063D0
      DISC%DynRup%filter(36,24) = -0.14822811905773394505D0
      DISC%DynRup%filter(36,25) = 0.29135418387524136098D-1
      DISC%DynRup%filter(36,26) = -0.96203197949176431915D-1
      DISC%DynRup%filter(36,27) = 0.15132661173687805147D0
      DISC%DynRup%filter(36,28) = -0.90771087646768295021D-1
      DISC%DynRup%filter(36,29) = -0.18913095779059291863D0
      DISC%DynRup%filter(36,30) = 0.28995092583146593845D0
      DISC%DynRup%filter(36,31) = 0.88070164006750962802D-2
      DISC%DynRup%filter(36,32) = -0.34292093188078363632D-1
      DISC%DynRup%filter(36,33) = 0.81174053365296697392D-1
      DISC%DynRup%filter(36,34) = -0.16230024600529912247D0
      DISC%DynRup%filter(36,35) = 0.25501534696734011638D0
      DISC%DynRup%filter(36,36) = 0.81888008130588519277D0
    ELSE
      logError(*) 'DISC%DynRup%filter not implemented for this order'
      STOP
    ENDIF

    ! Allocate and initialize magnitude output
    ALLOCATE(  DISC%DynRup%magnitude_out(MESH%Fault%nSide)                  )
    DISC%DynRup%magnitude_out(:) = .FALSE.

    IF (DISC%DynRup%magnitude_output_on.EQ.1) THEN
       ALLOCATE(  DISC%DynRup%averaged_Slip(MESH%Fault%nSide)        )
       !ini magnitude output
       DISC%DynRup%magnitude_out(:) = .TRUE.
       DISC%DynRup%averaged_Slip(:) = 0.0D0
    ENDIF

    ! ini rupture front output
    DISC%DynRup%RF = .FALSE.
    !ini dyn.stress ouput
    DISC%DynRup%DS = .FALSE.

    ! Initialize '+'side elements for RF und DS output
    IF ((DISC%DynRup%RFtime_on .EQ. 1) .AND. (DISC%DynRup%DS_output_on .EQ. 1) ) THEN
       ! Loop over every mesh element
           DO i = 1, MESH%Fault%nSide
              IF (MESH%FAULT%Face(i,1,1) .NE. 0) THEN
                 DISC%DynRup%RF(:,i) = .TRUE.
                 DISC%DynRup%DS(:,i) = .TRUE.
              ENDIF
           ENDDO
     ELSEIF ((DISC%DynRup%RFtime_on .EQ. 1) .AND. (DISC%DynRup%DS_output_on .EQ. 0 )) THEN
           DO i = 1, MESH%Fault%nSide
              IF (MESH%FAULT%Face(i,1,1) .NE. 0) THEN
                 DISC%DynRup%RF(:,i) = .TRUE.
              ENDIF
          ENDDO
    ENDIF

    faultParameterizedByTraction = c_interoperability_faultParameterizedByTraction(trim(DISC%DynRup%ModelFileName) // c_null_char)    
    
    if (faultParameterizedByTraction) then
      call c_interoperability_addFaultParameter("T_n" // c_null_char, EQN%IniBulk_xx)
      call c_interoperability_addFaultParameter("T_s" // c_null_char, EQN%IniShearXY)
      call c_interoperability_addFaultParameter("T_d" // c_null_char, EQN%IniShearXZ)
      EQN%IniBulk_yy(:,:) = 0.0d0
      EQN%IniBulk_zz(:,:) = 0.0d0
      EQN%IniShearYZ(:,:) = 0.0d0
    else
      call c_interoperability_addFaultParameter("s_xx" // c_null_char, EQN%IniBulk_xx)
      call c_interoperability_addFaultParameter("s_yy" // c_null_char, EQN%IniBulk_yy)
      call c_interoperability_addFaultParameter("s_zz" // c_null_char, EQN%IniBulk_zz)
      call c_interoperability_addFaultParameter("s_xy" // c_null_char, EQN%IniShearXY)
      call c_interoperability_addFaultParameter("s_yz" // c_null_char, EQN%IniShearYZ)
      call c_interoperability_addFaultParameter("s_xz" // c_null_char, EQN%IniShearXZ)
    endif

    !frictional parameter initialization
    SELECT CASE(EQN%FL)
    CASE(0)
       CONTINUE
    CASE(2,6,16)
       ALLOCATE(  DISC%DynRup%D_C(DISC%Galerkin%nBndGP,MESH%Fault%nSide)       )
       ALLOCATE(  DISC%DynRup%Mu_S(DISC%Galerkin%nBndGP,MESH%Fault%nSide)      )
       ALLOCATE(  DISC%DynRup%Mu_D(DISC%Galerkin%nBndGP,MESH%Fault%nSide)      )
       call c_interoperability_addFaultParameter("cohesion" // c_null_char, DISC%DynRup%cohesion)
       call c_interoperability_addFaultParameter("d_c" // c_null_char, DISC%DynRup%D_C)
       call c_interoperability_addFaultParameter("mu_s" // c_null_char, DISC%DynRup%Mu_S)
       call c_interoperability_addFaultParameter("mu_d" // c_null_char, DISC%DynRup%Mu_D)
       if (EQN%FL == 16) then
         ALLOCATE(  DISC%DynRup%forced_rupture_time(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
         call c_interoperability_addFaultParameter("forced_rupture_time" // c_null_char, DISC%DynRup%forced_rupture_time)
       end if

    CASE(3,4,7,101,103)
      ALLOCATE(  DISC%DynRup%RS_a_array(DISC%Galerkin%nBndGP, MESH%Fault%nSide)        )
      call c_interoperability_addFaultParameter("rs_a" // c_null_char, DISC%DynRup%RS_a_array)
      if (EQN%FL == 103) then
        allocate( DISC%DynRup%RS_srW_array(DISC%Galerkin%nBndGP, MESH%Fault%nSide), &
                  DISC%DynRup%RS_sl0_array(DISC%Galerkin%nBndGP,MESH%Fault%nSide),  &
                  nuc_xx(DISC%Galerkin%nBndGP,MESH%Fault%nSide),                    &
                  nuc_yy(DISC%Galerkin%nBndGP,MESH%Fault%nSide),                    &
                  nuc_zz(DISC%Galerkin%nBndGP,MESH%Fault%nSide),                    &
                  nuc_xy(DISC%Galerkin%nBndGP,MESH%Fault%nSide),                    &
                  nuc_yz(DISC%Galerkin%nBndGP,MESH%Fault%nSide),                    &
                  nuc_xz(DISC%Galerkin%nBndGP,MESH%Fault%nSide)                     )
        call c_interoperability_addFaultParameter("rs_srW" // c_null_char, DISC%DynRup%RS_srW_array)
        call c_interoperability_addFaultParameter("RS_sl0" // c_null_char, DISC%DynRup%RS_sl0_array)
        if (faultParameterizedByTraction) then
          call c_interoperability_addFaultParameter("Tnuc_n" // c_null_char, nuc_xx)
          call c_interoperability_addFaultParameter("Tnuc_s" // c_null_char, nuc_xy)
          call c_interoperability_addFaultParameter("Tnuc_d" // c_null_char, nuc_xz)
          nuc_yy(:,:) = 0.0d0
          nuc_zz(:,:) = 0.0d0
          nuc_yz(:,:) = 0.0d0
        else
          call c_interoperability_addFaultParameter("nuc_xx" // c_null_char, nuc_xx)
          call c_interoperability_addFaultParameter("nuc_yy" // c_null_char, nuc_yy)
          call c_interoperability_addFaultParameter("nuc_zz" // c_null_char, nuc_zz)
          call c_interoperability_addFaultParameter("nuc_xy" // c_null_char, nuc_xy)
          call c_interoperability_addFaultParameter("nuc_yz" // c_null_char, nuc_yz)
          call c_interoperability_addFaultParameter("nuc_xz" // c_null_char, nuc_xz)
        endif
      end if
    END SELECT

    call c_interoperability_initializeFault(  trim(DISC%DynRup%ModelFileName) // c_null_char, &
                                              EQN%GPwise,                                     &
                                              MESH%ELEM%BndGP_Tri,                            &
                                              DISC%Galerkin%nBndGP                            )

    ! Rotate initial stresses to fault coordinate system
    allocate(EQN%InitialStressInFaultCS(DISC%Galerkin%nBndGP,6,MESH%Fault%nSide))
    call rotateStressToFaultCS(EQN,MESH,DISC%Galerkin%nBndGP,EQN%IniBulk_xx,EQN%IniBulk_yy,EQN%IniBulk_zz,EQN%IniShearXY,EQN%IniShearYZ,EQN%IniShearXZ,EQN%InitialStressInFaultCS,faultParameterizedByTraction)
    
    if (EQN%FL == 103) then
      allocate(EQN%NucleationStressInFaultCS(DISC%Galerkin%nBndGP,6,MESH%Fault%nSide))
      call rotateStressToFaultCS(EQN,MESH,DISC%Galerkin%nBndGP,nuc_xx,nuc_yy,nuc_zz,nuc_xy,nuc_yz,nuc_xz,EQN%NucleationStressInFaultCS,faultParameterizedByTraction)
      deallocate(nuc_xx,nuc_yy,nuc_zz,nuc_xy,nuc_yz,nuc_xz)
    end if
  END SUBROUTINE DR_basic_ini
  
  SUBROUTINE rotateStressToFaultCS(EQN,MESH,nBndGP,s_xx,s_yy,s_zz,s_xy,s_yz,s_xz,stressInFaultCS,faultParameterizedByTraction)
    use JacobiNormal_mod, only: RotationMatrix3D
    use iso_c_binding, only: c_bool
    use create_fault_rotationmatrix_mod, only: create_fault_rotationmatrix
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    TYPE(tEquations)                      :: EQN
    TYPE(tUnstructMesh)                   :: MESH
    integer                               :: nBndGP
    real, allocatable, dimension(:,:)     :: s_xx,s_yy,s_zz,s_xy,s_yz,s_xz
    real, allocatable, dimension(:,:,:)   :: stressInFaultCS
    logical(kind=c_bool)                  :: faultParameterizedByTraction
    !-------------------------------------------------------------------------!
    ! Local variable declaration
    integer                             :: iFace, iBndGP
    real                                :: normal(3)
    real                                :: tangent1(3)
    real                                :: tangent2(3)
    real                                :: T(9,9)
    real                                :: iT(9,9)
    real                                :: rotmat(6,6)
    real                                :: iRotmat(6,6)
    real                                :: Stress(1:6,1:nBndGP)
    real                                :: StressinFaultCSTmp(6)
    !-------------------------------------------------------------------------!
    intent(in)                          :: EQN,MESH,nBndGP
    intent(inout)                       :: s_xx,s_yy,s_zz,s_xy,s_yz,s_xz
    intent(inout)                       :: stressInFaultCS
    
    do iFace = 1, MESH%Fault%nSide
      normal   = MESH%Fault%geoNormals( 1:3, iFace)
      tangent1 = MESH%Fault%geoTangent1(1:3, iFace)
      tangent2 = MESH%Fault%geoTangent2(1:3, iFace)
      CALL RotationMatrix3D(normal, tangent1, tangent2, T(:,:), iT(:,:), EQN)

      Stress(1,:) = s_xx(:,iFace)
      Stress(2,:) = s_yy(:,iFace)
      Stress(3,:) = s_zz(:,iFace)
      Stress(4,:) = s_xy(:,iFace)
      Stress(5,:) = s_yz(:,iFace)
      Stress(6,:) = s_xz(:,iFace)
      
      if (faultParameterizedByTraction) then
        call create_fault_rotationmatrix(rotmat,iFace,EQN,MESH,iRotmat)
        do iBndGP=1,nBndGP
          StressinFaultCSTmp = MATMUL(iRotmat, Stress(:,iBndGP))
          Stress(:,iBndGP) = StressinFaultCSTmp
        enddo
        s_xx(:,iFace) = Stress(1,:)
        s_yy(:,iFace) = Stress(2,:)
        s_zz(:,iFace) = Stress(3,:)
        s_xy(:,iFace) = Stress(4,:)
        s_yz(:,iFace) = Stress(5,:)
        s_xz(:,iFace) = Stress(6,:)
      endif

      do iBndGP=1,nBndGP
        StressinFaultCSTmp = MATMUL(iT(1:6,1:6), Stress(:,iBndGP))
        stressInFaultCS(iBndGP,:,iFace) = StressinFaultCSTmp
      enddo
    enddo
  END SUBROUTINE rotateStressToFaultCS


  !> Initialization of initial slip rate and friction for rate and state friction
  !<
  SUBROUTINE friction_RSF34(DISC,EQN,MESH,BND)
  !-------------------------------------------------------------------------!
  USE DGBasis_mod
  !-------------------------------------------------------------------------!
  IMPLICIT NONE
  !-------------------------------------------------------------------------!
  TYPE(tDiscretization), target  :: DISC
  TYPE(tEquations)               :: EQN
  TYPE(tUnstructMesh)            :: MESH
  TYPE (tBoundary)               :: BND
  !-------------------------------------------------------------------------!
  ! Local variable declaration
  INTEGER                        :: i
  INTEGER                        :: iSide,iElem,iBndGP
  INTEGER                        :: iLocalNeighborSide,iNeighbor
  INTEGER                        :: MPIIndex, iObject
  REAL                           :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
  REAL                           :: chi,tau
  REAL                           :: xi, eta, zeta, XGp, YGp, ZGp
  REAL                           :: iniSlipRate, X2
  !-------------------------------------------------------------------------!
  INTENT(IN)    :: MESH,BND
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------!

  EQN%IniSlipRate1 = DISC%DynRup%RS_iniSlipRate1
  EQN%IniSlipRate2 = DISC%DynRup%RS_iniSlipRate2
  iniSlipRate = SQRT(EQN%IniSlipRate1**2 + EQN%IniSlipRate2**2)

  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide

      ! element ID
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)

      ! get vertices of complete tet
      IF (MESH%Fault%Face(i,1,1) == 0) THEN
          ! iElem is in the neighbor domain
          ! The neighbor element belongs to a different MPI domain
          iNeighbor           = MESH%Fault%Face(i,1,2)          ! iNeighbor denotes "-" side
          iLocalNeighborSide  = MESH%Fault%Face(i,2,2)
          iObject  = MESH%ELEM%BoundaryToObject(iLocalNeighborSide,iNeighbor)
          MPIIndex = MESH%ELEM%MPINumber(iLocalNeighborSide,iNeighbor)
          !
          xV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(1,1:4,MPIIndex)
          yV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(2,1:4,MPIIndex)
          zV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(3,1:4,MPIIndex)
      ELSE
          !
          ! get vertices
          xV(1:4) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(1:4,iElem))
          yV(1:4) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(1:4,iElem))
          zV(1:4) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(1:4,iElem))
      ENDIF
      !
      DO iBndGP = 1,DISC%Galerkin%nBndGP ! Loop over all Gauss integration points
          !
          ! Transformation of boundary GP's into XYZ coordinate system
          chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
          tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
          CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
          CALL TetraTrafoXiEtaZeta2XYZ(xGp,yGp,zGp,xi,eta,zeta,xV,yV,zV)
          !
          EQN%IniStateVar(i,iBndGP) = DISC%DynRup%NucRS_sv0
          !EQN%IniStateVar(i,iBndGP) = DISC%DynRup%RS_sl0/DISC%DynRup%RS_sr0*EXP((sstress/(nstress*DISC%DynRup%RS_b))-DISC%DynRup%RS_f0/DISC%DynRup%RS_b-DISC%DynRup%RS_a_array(iBndGP,i)/DISC%DynRup%RS_b*LOG(iniSlipRate/DISC%DynRup%RS_sr0))
          X2  = iniSlipRate*0.5/DISC%DynRup%RS_sr0 * EXP((DISC%DynRup%RS_f0 + DISC%DynRup%RS_b*LOG(DISC%DynRup%RS_sr0*EQN%IniStateVar(i,iBndGP)/DISC%DynRup%RS_sl0)) / DISC%DynRup%RS_a)
          EQN%IniMu(iBndGP,i)=DISC%DynRup%RS_a * LOG(X2 + SQRT(X2**2 + 1.0))

      ENDDO ! iBndGP

  ENDDO !    MESH%Fault%nSide

  END SUBROUTINE friction_RSF34      ! Initialization of initial slip rate and friction for rate and state friction

  !> Initialization of initial slip rate and friction for fast velocity weakening friction
  !<
  SUBROUTINE friction_RSF7(DISC,EQN,MESH,BND)
  !-------------------------------------------------------------------------!
  USE DGBasis_mod
  !-------------------------------------------------------------------------!
  IMPLICIT NONE
  !-------------------------------------------------------------------------!
  TYPE(tDiscretization), target  :: DISC
  TYPE(tEquations)               :: EQN
  TYPE(tUnstructMesh)            :: MESH
  TYPE (tBoundary)               :: BND
  !-------------------------------------------------------------------------!
  ! Local variable declaration
  INTEGER                        :: i
  INTEGER                        :: iSide,iElem,iBndGP
  INTEGER                        :: iLocalNeighborSide,iNeighbor
  INTEGER                        :: MPIIndex, iObject
  REAL                           :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
  REAL                           :: chi,tau
  REAL                           :: xi, eta, zeta, XGp, YGp, ZGp
  REAL                           :: iniSlipRate
  !-------------------------------------------------------------------------!
  INTENT(IN)    :: MESH,BND
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------!

  EQN%IniSlipRate1 = DISC%DynRup%RS_iniSlipRate1
  EQN%IniSlipRate2 = DISC%DynRup%RS_iniSlipRate2
  iniSlipRate = SQRT(EQN%IniSlipRate1**2 + EQN%IniSlipRate2**2)

  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide

      ! element ID
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)

      ! get vertices of complete tet
      IF (MESH%Fault%Face(i,1,1) == 0) THEN
          ! iElem is in the neighbor domain
          ! The neighbor element belongs to a different MPI domain
          iNeighbor           = MESH%Fault%Face(i,1,2)          ! iNeighbor denotes "-" side
          iLocalNeighborSide  = MESH%Fault%Face(i,2,2)
          iObject  = MESH%ELEM%BoundaryToObject(iLocalNeighborSide,iNeighbor)
          MPIIndex = MESH%ELEM%MPINumber(iLocalNeighborSide,iNeighbor)
          !
          xV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(1,1:4,MPIIndex)
          yV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(2,1:4,MPIIndex)
          zV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(3,1:4,MPIIndex)
      ELSE
          !
          ! get vertices
          xV(1:4) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(1:4,iElem))
          yV(1:4) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(1:4,iElem))
          zV(1:4) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(1:4,iElem))
      ENDIF
      !
      DO iBndGP = 1,DISC%Galerkin%nBndGP ! Loop over all Gauss integration points
          !
          ! Transformation of boundary GP's into XYZ coordinate system
          chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
          tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
          CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
          CALL TetraTrafoXiEtaZeta2XYZ(xGp,yGp,zGp,xi,eta,zeta,xV,yV,zV)
                   !
          EQN%IniStateVar(i,iBndGP)= (DISC%DynRup%RS_sl0*(DISC%DynRup%RS_a*EQN%IniBulk_yy(i,iBndGP)*iniSlipRate &
              + (iniSlipRate + DISC%DynRup%RS_sr0)*(DISC%DynRup%RS_f0*EQN%IniBulk_yy(i,iBndGP)-EQN%IniShearXY(i,iBndGP)))) &
              / (DISC%DynRup%RS_a*EQN%IniBulk_yy(i,iBndGP)*iniSlipRate-(iniSlipRate + DISC%DynRup%RS_sr0) &
              * (DISC%DynRup%RS_b*EQN%IniBulk_yy(i,iBndGP)-DISC%DynRup%RS_f0*EQN%IniBulk_yy(i,iBndGP)+EQN%IniShearXY(i,iBndGP)))
          EQN%IniMu(iBndGP,i)=DISC%DynRup%RS_f0+DISC%DynRup%RS_a*iniSlipRate/(iniSlipRate+DISC%DynRup%RS_sr0)-DISC%DynRup%RS_b*EQN%IniStateVar(i,iBndGP)/(EQN%IniStateVar(i,iBndGP)+DISC%DynRup%RS_sl0)
      ENDDO ! iBndGP

  ENDDO !    MESH%Fault%nSide

  END SUBROUTINE friction_RSF7      ! Initialization of initial slip rate and friction for fast velocity weakening friction

  !> Initialization of initial slip rate and friction for SCEC TPV101
  !<
  SUBROUTINE friction_RSF101(DISC,EQN,MESH,BND)
  !-------------------------------------------------------------------------!
  USE DGBasis_mod
  !-------------------------------------------------------------------------!
  IMPLICIT NONE
  !-------------------------------------------------------------------------!
  TYPE(tDiscretization), target  :: DISC
  TYPE(tEquations)               :: EQN
  TYPE(tUnstructMesh)            :: MESH
  TYPE (tBoundary)               :: BND
  !-------------------------------------------------------------------------!
  ! Local variable declaration
  INTEGER                        :: i
  INTEGER                        :: iSide,iElem,iBndGP
  INTEGER                        :: iLocalNeighborSide,iNeighbor
  INTEGER                        :: MPIIndex, iObject
  REAL                           :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
  REAL                           :: chi,tau
  REAL                           :: xi, eta, zeta, XGp, YGp, ZGp
  REAL                           :: iniSlipRate, tmp
  !-------------------------------------------------------------------------!
  INTENT(IN)    :: MESH,BND
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------!

  EQN%IniSlipRate1 = DISC%DynRup%RS_iniSlipRate1
  EQN%IniSlipRate2 = DISC%DynRup%RS_iniSlipRate2
  iniSlipRate = SQRT(EQN%IniSlipRate1**2 + EQN%IniSlipRate2**2)

  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide

      ! element ID
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)

      ! get vertices of complete tet
      IF (MESH%Fault%Face(i,1,1) == 0) THEN
          ! iElem is in the neighbor domain
          ! The neighbor element belongs to a different MPI domain
          iNeighbor           = MESH%Fault%Face(i,1,2)          ! iNeighbor denotes "-" side
          iLocalNeighborSide  = MESH%Fault%Face(i,2,2)
          iObject  = MESH%ELEM%BoundaryToObject(iLocalNeighborSide,iNeighbor)
          MPIIndex = MESH%ELEM%MPINumber(iLocalNeighborSide,iNeighbor)
          !
          xV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(1,1:4,MPIIndex)
          yV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(2,1:4,MPIIndex)
          zV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(3,1:4,MPIIndex)
      ELSE
          !
          ! get vertices
          xV(1:4) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(1:4,iElem))
          yV(1:4) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(1:4,iElem))
          zV(1:4) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(1:4,iElem))
      ENDIF
      !
      DO iBndGP = 1,DISC%Galerkin%nBndGP ! Loop over all Gauss integration points
          !
          ! Transformation of boundary GP's into XYZ coordinate system
          chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
          tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
          CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
          CALL TetraTrafoXiEtaZeta2XYZ(xGp,yGp,zGp,xi,eta,zeta,xV,yV,zV)
               !
          EQN%IniStateVar(i,iBndGP) = (DISC%DynRup%RS_sl0/DISC%DynRup%RS_sr0) * EXP((-EQN%IniShearXY(i,iBndGP)/EQN%IniBulk_yy(i,iBndGP)-DISC%DynRup%RS_f0-DISC%DynRup%RS_a_array(iBndGP,i)*LOG(iniSlipRate/DISC%DynRup%RS_sr0))/DISC%DynRup%RS_b)
          ! ASINH(X)=LOG(X+SQRT(X^2+1))
          tmp  = iniSlipRate*0.5/DISC%DynRup%RS_sr0 * EXP((DISC%DynRup%RS_f0 + DISC%DynRup%RS_b*LOG(DISC%DynRup%RS_sr0*EQN%IniStateVar(i,iBndGP)/DISC%DynRup%RS_sl0)) / DISC%DynRup%RS_a_array(iBndGP,i))
          EQN%IniMu(iBndGP,i)=DISC%DynRup%RS_a_array(iBndGP,i) * LOG(tmp + SQRT(tmp**2 + 1.0))

      ENDDO ! iBndGP

  ENDDO !    MESH%Fault%nSide

  END SUBROUTINE friction_RSF101      ! Initialization of initial slip rate and friction for SCEC TPV101

  !> Initialization of initial slip rate and friction for SCEC TPV103
  !<
  SUBROUTINE friction_RSF103(DISC,EQN,MESH,BND)
  !-------------------------------------------------------------------------!
  USE DGBasis_mod
  USE JacobiNormal_mod
  !-------------------------------------------------------------------------!
  IMPLICIT NONE
  !-------------------------------------------------------------------------!
  TYPE(tDiscretization), target  :: DISC
  TYPE(tEquations)               :: EQN
  TYPE(tUnstructMesh)            :: MESH
  TYPE (tBoundary)               :: BND
  !-------------------------------------------------------------------------!
  ! Local variable declaration
  INTEGER                        :: iFace, iBndGP
  REAL                           :: iniSlipRate, tmp

  !-------------------------------------------------------------------------!
  INTENT(IN)    :: MESH,BND
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------!

  EQN%IniSlipRate1 = DISC%DynRup%RS_iniSlipRate1
  EQN%IniSlipRate2 = DISC%DynRup%RS_iniSlipRate2
  iniSlipRate = SQRT(EQN%IniSlipRate1**2 + EQN%IniSlipRate2**2)


  ! Loop over every mesh element
  DO iFace = 1, MESH%Fault%nSide
      DO iBndGP = 1,DISC%Galerkin%nBndGP ! Loop over all Gauss integration points
          tmp = ABS(SQRT(EQN%InitialStressInFaultCS(iBndGP,4,iFace)**2+EQN%InitialStressInFaultCS(iBndGP,6,iFace)**2)/(DISC%DynRup%RS_a_array(iBndGP,iFace)*EQN%InitialStressInFaultCS(iBndGP,1,iFace)))
          EQN%IniStateVar(iBndGP,iFace)=DISC%DynRup%RS_a_array(iBndGP,iFace)*LOG(2.0D0*DISC%DynRup%RS_sr0/iniSlipRate * (EXP(tmp)-EXP(-tmp))/2.0D0)
          ! ASINH(X)=LOG(X+SQRT(X^2+1))
          tmp  = iniSlipRate*0.5/DISC%DynRup%RS_sr0 * EXP(EQN%IniStateVar(iBndGP,iFace)/ DISC%DynRup%RS_a_array(iBndGP,iFace))
          EQN%IniMu(iBndGP,iFace)=DISC%DynRup%RS_a_array(iBndGP,iFace) * LOG(tmp + SQRT(tmp**2 + 1.0D0))
      ENDDO ! iBndGP
  ENDDO !    MESH%Fault%nSide

  END SUBROUTINE friction_RSF103      ! Initialization of initial slip rate and friction for SCEC TPV103

  !> Initialization of friction for linear slip weakening
  !<
  SUBROUTINE friction_LSW(DISC,EQN,MESH,BND)
  !-------------------------------------------------------------------------!
  USE DGBasis_mod
  !-------------------------------------------------------------------------!
  IMPLICIT NONE
  !-------------------------------------------------------------------------!
  TYPE(tDiscretization), target  :: DISC
  TYPE(tEquations)               :: EQN
  TYPE(tUnstructMesh)            :: MESH
  TYPE (tBoundary)               :: BND
  !-------------------------------------------------------------------------!
  INTENT(IN)    :: MESH,BND
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------!

  EQN%IniMu(:,:) = DISC%DynRup%Mu_S(:,:)

  END SUBROUTINE friction_LSW

    SUBROUTINE friction_LSW6(DISC,EQN,MESH,BND)
  !-------------------------------------------------------------------------!
  USE DGBasis_mod
  !-------------------------------------------------------------------------!
  IMPLICIT NONE
  !-------------------------------------------------------------------------!
  TYPE(tDiscretization), target  :: DISC
  TYPE(tEquations)               :: EQN
  TYPE(tUnstructMesh)            :: MESH
  TYPE (tBoundary)               :: BND
  !-------------------------------------------------------------------------!
  ! Local variable declaration
  INTEGER                        :: iFace,iBndGP
  !-------------------------------------------------------------------------!
  INTENT(IN)    :: MESH,BND
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------!
  
  EQN%IniMu(:,:) = DISC%DynRup%Mu_S(:,:)
  do iFace = 1, MESH%Fault%nSide
    do iBndGP = 1,DISC%Galerkin%nBndGP
      DISC%DynRup%Strength(iBndGP,iFace) = EQN%IniMu(iBndGP,iFace) * EQN%InitialStressInFaultCS(iBndGP,1,iFace)
    enddo
  enddo

  END SUBROUTINE friction_LSW6    ! Initialization of friction for bimaterial linear slip weakening



  END MODULE ini_model_DR_mod
