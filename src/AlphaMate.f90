
#ifdef BINARY
#define BINFILE ,form="unformatted"
#else
#define BINFILE
#endif

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#ifdef OS_UNIX
#define DASH "/"
#define COPY "cp"
#define MKDIR "mkdir -p"
#define RMDIR "rm -r"
#define RM "rm"
#define RENAME "mv"
#else
#define DASH "\"
#define COPY "copy"
#define MKDIR "md"
#define RMDIR "rmdir /S"
#define RM "del"
#define RENAME "move /Y"
#endif

!###############################################################################

module AlphaMateMod

  use ISO_Fortran_Env, STDIN=>input_unit,STDOUT=>output_unit,STDERR=>error_unit
  use IFPort,only : SystemQQ
  use OrderPackMod,only : MrgRnk
  use AlphaHouseMod,only : CountLines,Int2Char,Real2Char,RandomOrder,SetSeed,ToLower
  use AlphaStatMod,only : CalcDescStat,DescStatD,CalcDescStatSymMatrix,DescStatMatrixD
  use AlphaEvolveMod,only : EvolAlgDE,EvolveCrit

  implicit none

  integer(int32) :: nInd,nMat,nPar,nPotPar1,nPotPar2,nMal,nFem,nPar1,nPar2,nFrontierSteps
  integer(int32) :: EvolAlgNSol,EvolAlgNGen,EvolAlgNGenBurnIn,EvolAlgNGenStop,EvolAlgNGenPrint
  integer(int32) :: PAGEPar1Max,PAGEPar2Max
  integer(int32),allocatable :: Gender(:),IdPotPar1(:),IdPotPar2(:)
  integer(int32),allocatable :: nVecPar1(:),nVecPar2(:),nVec(:),MateAlloc(:,:)

  real(real64) :: LimitPar1Min,LimitPar1Max,LimitPar2Min,LimitPar2Max
  real(real64) :: EvolAlgStopTol,EvolAlgCRBurnIn,EvolAlgCRLate,EvolAlgFBase,EvolAlgFHigh1,EvolAlgFHigh2
  real(real64) :: PopInbOld,PopInbTarget,RatePopInbTarget,GainMinStand
  real(real64) :: PopInbPenalty,PrgInbPenalty,SelfingPenalty,LimitPar1Penalty,LimitPar2Penalty
  real(real64) :: PAGEPar1Cost,PAGEPar2Cost
  real(real64),allocatable :: Bv(:),BvStand(:),BvPAGE(:),BvPAGEStand(:)
  real(real64),allocatable :: RelMtx(:,:),RatePopInbFrontier(:),xVec(:),GeneEdit(:)

  logical :: ModeMin,ModeOpt,BvAvailable,GenderMatters,EqualizePar1,EqualizePar2
  logical :: SelfingAllowed,PopInbPenaltyBellow,InferPopInbOld,EvaluateFrontier
  logical :: PAGE,PAGEPar1,PAGEPar2

  character(len=100),allocatable :: IdC(:)
  CHARACTER(len=100),PARAMETER :: FMTREAL2CHAR="(f11.5)"
  CHARACTER(len=100),PARAMETER :: FMTLOGHEADERSTDOUT="(11a12)"
  CHARACTER(len=100),PARAMETER :: FMTLOGSTDOUT="(i12,10(1x,f11.5))"
  CHARACTER(len=100),PARAMETER :: FMTLOGHEADERUNIT="(12a,11a22)"
  CHARACTER(len=100),PARAMETER :: FMTLOGUNIT="(i12,11(1x,es21.14))"
  CHARACTER(len=100),PARAMETER :: FMTINDHEAD="(6a12)"
  CHARACTER(len=100),PARAMETER :: FMTINDHEADEDIT="(8a12)"
  CHARACTER(len=100),PARAMETER :: FMTIND="(a12,1x,i11,3(1x,f11.5),1x,i11)"
  CHARACTER(len=100),PARAMETER :: FMTINDEDIT="(a12,1x,i11,3(1x,f11.5),2(1x,i11),1x,f11.5)"
  CHARACTER(len=100),PARAMETER :: FMTMATHEAD="(3a12)"
  CHARACTER(len=100),PARAMETER :: FMTMAT="(i12,2(1x,a11))"
  CHARACTER(len=100),PARAMETER :: FMTFROHEAD="(a12,7a22)"
  CHARACTER(len=100),PARAMETER :: FMTFRO="(a12,7(1x,es21.14))"

  private
  public :: AlphaMateTitles,ReadSpecAndDataForAlphaMate,SetInbreedingParameters
  public :: AlphaMateSearch,EvolAlgLogHeaderForAlphaMate,EvolAlgLogForAlphaMate
  public :: FixSolEtcMateAndCalcCrit

  contains

    !###########################################################################

    ! This note is just to make clear what the inbreeding calculations are doing.
    !
    ! The total inbreeding that we see today is:
    !
    ! F_total = F_recent + (1 - F_recent) * F_ancient,
    !
    ! where F_recent is recent/new inbreeding after the previous time point
    !       F_ancient is ancient/old inbreeding before to the previous time point
    !
    ! Equivalent formulas with some other notation are:
    !
    ! F_T = F_IS + (1 - F_IS) * F_ST
    !
    ! where F_T is (total) inbreeding relative to ancient/old base
    !       F_IS is (individual-to-subtotal=recent/new) inbreeding relative to recent base
    !       F-ST is (subtotal-to-total=ancient/old) inbreeding at recent base relative to ancient/old base
    !
    ! F_t = DeltaF + (1 - DeltaF) * F_t-1
    !
    ! where F_t is (total) inbreeding in generation t
    !       DeltaF is recent/new inbreeding after generation t-1 (=rate of inbreeding)
    !       F-t-1 is ancient/old inbreeding before generation t
    !
    ! Clearly:
    !
    ! F_total   = F_T  = F_t
    ! F_recent  = F_IS = DeltaF
    ! F_ancient = F_ST = F_t-1

    !###########################################################################

    subroutine AlphaMateTitles
      write(STDOUT,"(a)") ""
      write(STDOUT,"(a30,a,a30)") " ","**********************"," "
      write(STDOUT,"(a30,a,a30)") " ","*                    *"," "
      write(STDOUT,"(a30,a,a30)") " ","*     AlphaMate      *"," "
      write(STDOUT,"(a30,a,a30)") " ","*                    *"," "
      write(STDOUT,"(a30,a,a30)") " ","**********************"
      write(STDOUT,"(a30,a,a30)") " ","VERSION:"//TOSTRING(VERS)," "
      write(STDOUT,"(a15,a)")     " ","Software for optimizing contributions to the next generation"
      write(STDOUT,"(a)") ""
      write(STDOUT,"(a35,a)")     " ","No Liability"
      write(STDOUT,"(a25,a)")     " ","Bugs to Gregor.Gorjanc@roslin.ed.ac.uk"
      write(STDOUT,"(a)") ""
    end subroutine

    !###########################################################################

    subroutine ReadSpecAndDataForAlphaMate

      implicit none

      integer(int32) :: i,j,DumI,jMal,jFem,nIndTmp,GenderTmp,Seed
      integer(int32) :: UnitSpec,UnitRelMtx,UnitBv,UnitGender
      integer(int32),allocatable :: Order(:)

      real(real64) :: BvTmp,BvTmp2

      logical :: Success

      character(len=1000) :: RelMtxFile,BvFile,GenderFile,SeedFile
      character(len=100) :: DumC,DumC2,DumC3,IdCTmp

      type(DescStatD) :: BvDescStat
      type(DescStatMatrixD) :: RelDescStat

      ! --- Spec file ---

      write(STDOUT,"(a)") "--- Specifications ---"
      write(STDOUT,"(a)") " "

      Success=SystemQQ(COPY//" AlphaMateSpec.txt AlphaMateResults"//DASH//"AlphaMateSpec.txt")
      if (.not.Success) then
        write(STDERR,"(a)") "ERROR: Failed to copy the AlphaMateSpec.txt file in the output folder!"
        write(STDERR,"(a)") " "
        stop 1
      end if

      open(newunit=UnitSpec,file="AlphaMateSpec.txt",status="old")
      write(STDOUT,"(a)") "SpecFile: AlphaMateSpec.txt"

      ! Mode
      read(UnitSpec,*) DumC,DumC
      if (ToLower(trim(DumC)) == "minthenopt") then
        ModeMin=.true.
        ModeOpt=.true.
        write(STDOUT,"(a)") "Mode: MinThenOpt"
      else if (ToLower(trim(DumC)) == "min") then
        ModeMin=.true.
        ModeOpt=.false.
        write(STDOUT,"(a)") "Mode: Min"
      else if (ToLower(trim(DumC)) == "opt") then
        ModeMin=.false.
        ModeOpt=.true.
        write(STDOUT,"(a)") "Mode: Opt"
      else
        write(STDERR,"(a)") "ERROR: Mode must be: Min, Opt, or MinThenOpt!"
        write(STDERR,"(a)") " "
        stop 1
      end if

      ! RelationshipMatrixFile
      read(UnitSpec,*) DumC,RelMtxFile
      write(STDOUT,"(2a)") "RelationshipMatrixFile: ",trim(RelMtxFile)
      nInd=CountLines(RelMtxFile)

      ! BreedingValueFile
      read(UnitSpec,*) DumC,BvFile
      if (ToLower(trim(BvFile)) /= "none") then
        BvAvailable=.true.
        write(STDOUT,"(2a)") "BreedingValueFile: ",trim(BvFile)
        nIndTmp=CountLines(BvFile)
        if (nIndTmp /= nInd) then
          write(STDERR,"(a)") "ERROR: Number of individuals in Ebv file and Relationship Matrix file is not the same!"
          write(STDERR,"(a)") " "
          stop 1
        end if
      else
        BvAvailable=.false.
      end if

      ! GenderFile
      read(UnitSpec,*) DumC,GenderFile
      if (ToLower(trim(GenderFile)) /= "none") then
        GenderMatters=.true.
        write(STDOUT,"(2a)") "GenderFile: ",trim(GenderFile)
        nIndTmp=CountLines(GenderFile)
        if (nIndTmp /= nInd) then
          write(STDERR,"(a)") "ERROR: Number of individuals in Gender file and Relationship Matrix file is not the same!"
          write(STDERR,"(a)") " "
          stop 1
        end if
      else
        GenderMatters=.false.
      end if

      ! NumberOfIndividuals
      read(UnitSpec,*) DumC,nInd
      DumC=Int2Char(nInd)
      write(STDOUT,"(2a)") "NumberOfIndividuals: ",trim(adjustl(DumC))

      ! NumberOfMatings
      read(UnitSpec,*) DumC,nMat
      DumC=Int2Char(nMat)
      write(STDOUT,"(2a)") "NumberOfMatings: ",trim(adjustl(DumC))
      ! TODO: In animals one would not be able to generate more matings than there is individuals, i.e.,
      !       10 males and 10 females can give 10 matings only (females are a bottleneck). But if we do
      !       plants or collect ova from females, then we could technically do 10*10=100 matings (each male
      !       with each female).
      ! if (nMat > nInd) then
      !   write(STDOUT,"(a)") "NOTE: The number of matings is larger than the number of all individuals! Was this really the intention?"
      !   DumC=Int2Char(nMat)
      !   write(STDOUT,"(2a)") "NOTE: Number of     matings: ",trim(adjustl(DumC))
      !   DumC=Int2Char(nInd)
      !   write(STDOUT,"(2a)") "NOTE: Number of individuals: ",trim(adjustl(DumC))
      !   write(STDOUT,"(a)") " "
      ! end if

      ! NumberOfParents
      read(UnitSpec,*) DumC,nPar
      DumC=Int2Char(nPar)
      write(STDOUT,"(2a)") "NumberOfParents: ",trim(adjustl(DumC))
      if (nPar > nInd) then
        write(STDERR,"(a)") "ERROR: The number of parents can not be larger than the number of individuals!"
        DumC=Int2Char(nPar)
        write(STDERR,"(2a)") "ERROR: Number of     parents: ",trim(adjustl(DumC))
        DumC=Int2Char(nInd)
        write(STDERR,"(2a)") "ERROR: Number of individuals: ",trim(adjustl(DumC))
        write(STDERR,"(a)") " "
        stop 1
      end if
      if (nMat > nPar) then
        write(STDOUT,"(a)") "NOTE: The number of matings is larger than the number of parents! Was this really the intention?"
        DumC=Int2Char(nMat)
        write(STDOUT,"(2a)") "NOTE: Number of matings: ",trim(adjustl(DumC))
        DumC=Int2Char(nPar)
        write(STDOUT,"(2a)") "NOTE: Number of parents: ",trim(adjustl(DumC))
        write(STDOUT,"(a)") " "
      end if

      ! NumberOfMaleParents
      read(UnitSpec,*) DumC,nPar1

      ! NumberOfFemaleParents
      read(UnitSpec,*) DumC,nPar2

      if (.not.GenderMatters) then
        nPar1=nPar
      else
        DumC=Int2Char(nPar1)
        write(STDOUT,"(2a)") "NumberOfMaleParents: ",trim(adjustl(DumC))
        DumC=Int2Char(nPar2)
        write(STDOUT,"(2a)") "NumberOfFemaleParents: ",trim(adjustl(DumC))
      end if
      if (GenderMatters) then
        if ((nPar1+nPar2) > nInd) then
          write(STDERR,"(a)") "ERROR: The number of parents can not be larger than the number of individuals!"
          DumC=Int2Char(nPar1+nPar2)
          write(STDERR,"(2a)") "ERROR: Number of        parents: ",trim(adjustl(DumC))
          DumC=Int2Char(nPar1)
          write(STDERR,"(2a)") "ERROR: Number of   male parents: ",trim(adjustl(DumC))
          DumC=Int2Char(nPar2)
          write(STDERR,"(2a)") "ERROR: Number of female parents: ",trim(adjustl(DumC))
          DumC=Int2Char(nInd)
          write(STDERR,"(2a)") "ERROR: Number of    individuals: ",trim(adjustl(DumC))
          write(STDERR,"(a)") " "
          stop 1
        end if
        if ((nPar1+nPar2) /= nPar) then
          write(STDOUT,"(a)") "NOTE: The number of male and female parents does not match with the total number of parents - redefined!"
          DumC=Int2Char(nPar1)
          write(STDOUT,"(2a)") "NOTE: Number of   male parents: ",trim(adjustl(DumC))
          DumC=Int2Char(nPar2)
          write(STDOUT,"(2a)") "NOTE: Number of female parents: ",trim(adjustl(DumC))
          DumC=Int2Char(nPar)
          write(STDOUT,"(3a)") "NOTE: Number of        parents: ",trim(adjustl(DumC))," (defined)"
          nPar=nPar1+nPar2
          DumC=Int2Char(nPar)
          write(STDOUT,"(3a)") "NOTE: Number of        parents: ",trim(adjustl(DumC))," (redefined)"
          write(STDOUT,"(a)") " "
        end if
        if ((nMat > nPar1) .and. (nMat > nPar2)) then
          write(STDOUT,"(a)") "NOTE: The number of matings is larger than the number of male and female parents! Was this really the intention?"
          DumC=Int2Char(nMat)
          write(STDOUT,"(2a)") "NOTE: Number of        matings: ",trim(adjustl(DumC))
          DumC=Int2Char(nPar1)
          write(STDOUT,"(2a)") "NOTE: Number of   male parents: ",trim(adjustl(DumC))
          DumC=Int2Char(nPar2)
          write(STDOUT,"(2a)") "NOTE: Number of female parents: ",trim(adjustl(DumC))
          write(STDOUT,"(a)") " "
        end if
      end if

      ! EqualizeParentContributions
      read(UnitSpec,*) DumC,DumC
      if (.not.GenderMatters) then
        if      (ToLower(trim(DumC)) == "no") then
          EqualizePar1=.true.
          write(STDOUT,"(a)") "EqualizeParentContributions: yes"
        else if (ToLower(trim(DumC)) == "yes") then
          EqualizePar1=.false.
          write(STDOUT,"(a)") "EqualizeParentContributions: no"
        else
          write(STDERR,"(a)") "ERROR: EqualizeParentContributions must be: Yes or no!"
          write(STDERR,"(a)") " "
          stop 1
        end if
      end if

      ! EqualizeMaleParentContributions
      read(UnitSpec,*) DumC,DumC
      if (GenderMatters) then
        if      (ToLower(trim(DumC)) == "yes") then
          if (mod(nMat,nPar1) /= 0) then
            ! TODO: might consider handling this better at some point
            write(STDERR,"(a)") "ERROR: When contributions are equalized the number of matings needs to divide into the number of parents"
            DumC=Int2Char(nMat)
            write(STDERR,"(2a)") "ERROR: Number of       matings: ",trim(adjustl(DumC))
            DumC=Int2Char(nPar1)
            write(STDERR,"(2a)") "ERROR: Number of  male parents: ",trim(adjustl(DumC))
            DumC=Int2Char(mod(nMat,nPar1))
            write(STDERR,"(2a)") "ERROR: Modulo (should be zero): ",trim(adjustl(DumC))
            write(STDERR,"(a)") " "
            stop 1
          end if
          EqualizePar1=.true.
          write(STDOUT,"(a)") "EqualizeMaleParentContributions: yes"
        else if (ToLower(trim(DumC)) == "no") then
          EqualizePar1=.false.
          write(STDOUT,"(a)") "EqualizeMaleParentContributions: no"
        else
          write(STDERR,"(a)") "ERROR: EqualizeMaleParentContributions must be: Yes or No!"
          write(STDERR,"(a)") " "
          stop 1
        end if
      end if

      ! EqualizeFemaleParentContributions
      read(UnitSpec,*) DumC,DumC
      if (GenderMatters) then
        if      (ToLower(trim(DumC)) == "yes") then
          if (mod(nMat,nPar2) /= 0) then
            ! TODO: might consider handling this better at some point
            write(STDERR,"(a)") "ERROR: When contributions are equalized the number of matings needs to divide into the number of parents"
            DumC=Int2Char(nMat)
            write(STDERR,"(2a)") "ERROR: Number of        matings: ",trim(adjustl(DumC))
            DumC=Int2Char(nPar2)
            write(STDERR,"(2a)") "ERROR: Number of female parents: ",trim(adjustl(DumC))
            DumC=Int2Char(mod(nMat,nPar2))
            write(STDERR,"(2a)") "ERROR: Modulo  (should be zero): ",trim(adjustl(DumC))
            write(STDERR,"(a)") " "
            stop 1
          end if
          EqualizePar2=.true.
          write(STDOUT,"(a)") "EqualizeFemaleParentContributions: yes"
        else if (ToLower(trim(DumC)) == "no") then
          EqualizePar2=.false.
          write(STDOUT,"(a)") "EqualizeFemaleParentContributions: no"
        else
          write(STDERR,"(a)") "ERROR: EqualizeFemaleParentContributions must be: Yes or No!"
          write(STDERR,"(a)") " "
          stop 1
        end if
      end if

      ! LimitParentContributions
      read(UnitSpec,*) DumC,DumC
      LimitPar1Min=1.0d0
      LimitPar1Max=huge(LimitPar1Max)-1.0d0
      LimitPar1Penalty=0.0d0
      if (.not.GenderMatters) then
        if      (ToLower(trim(DumC)) == "yes") then
          if (EqualizePar1) then
            write(STDOUT,"(a)") "LimitParentContributions: no"
            write(STDOUT,"(a)") "NOTE: Limit parent contributions option ignored when contributions are to be equalized!"
            write(STDOUT,"(a)") " "
          else
            backspace(UnitSpec)
            read(UnitSpec,*) DumC,DumC,LimitPar1Min,LimitPar1Max,LimitPar1Penalty
            DumC=Int2Char(nint(LimitPar1Min))
            DumC2=Int2Char(nint(LimitPar1Max))
            DumC3=Real2Char(LimitPar1Penalty,fmt=FMTREAL2CHAR)
            write(STDOUT,"(6a)") "LimitParentContributions: yes, min ",&
              trim(adjustl(DumC)),", max ",trim(adjustl(DumC2)),", penalty ",trim(adjustl(DumC3))
          end if
        else if (ToLower(trim(DumC)) == "no") then
          write(STDOUT,"(a)") "LimitParentContributions: no"
        else
          write(STDERR,"(a)") "ERROR: LimitParentContributions must be: Yes or No!"
          write(STDERR,"(a)") " "
          stop 1
        end if
      end if

      ! LimitMaleParentContributions
      read(UnitSpec,*) DumC,DumC
      if (GenderMatters) then
        if      (ToLower(trim(DumC)) == "yes") then
          if (EqualizePar1) then
            write(STDOUT,"(a)") "LimitMaleParentContributions: no"
            write(STDOUT,"(a)") "NOTE: Limit male parent contributions option ignored when contributions are to be equalized!"
            write(STDOUT,"(a)") " "
          else
            backspace(UnitSpec)
            read(UnitSpec,*) DumC,DumC,LimitPar1Min,LimitPar1Max,LimitPar1Penalty
            DumC=Int2Char(nint(LimitPar1Min))
            DumC2=Int2Char(nint(LimitPar1Max))
            DumC3=Real2Char(LimitPar1Penalty,fmt=FMTREAL2CHAR)
            write(STDOUT,"(6a)") "LimitMaleParentContributions: yes, min ",&
              trim(adjustl(DumC)),", max ",trim(adjustl(DumC)),", penalty ",trim(adjustl(DumC3))
          end if
        else if (ToLower(trim(DumC)) == "no") then
          write(STDOUT,"(a)") "LimitMaleParentContributions: no"
        else
          write(STDERR,"(a)") "ERROR: LimitMaleParentContributions must be: Yes or No!"
          write(STDERR,"(a)") " "
          stop 1
        end if
      end if

      ! LimitFemaleParentContributions
      read(UnitSpec,*) DumC,DumC
      LimitPar2Min=1.0d0
      LimitPar2Max=huge(LimitPar2Max)-1.0d0
      LimitPar2Penalty=0.0d0
      if (GenderMatters) then
        if      (ToLower(trim(DumC)) == "yes") then
          if (EqualizePar2) then
            write(STDOUT,"(a)") "LimitFemaleParentContributions: no"
            write(STDOUT,"(a)") "NOTE: Limit female parent contributions option ignored when contributions are to be equalized!"
            write(STDOUT,"(a)") " "
          else
            backspace(UnitSpec)
            read(UnitSpec,*) DumC,DumC,LimitPar2Min,LimitPar2Max,LimitPar2Penalty
            DumC=Int2Char(nint(LimitPar2Min))
            DumC2=Int2Char(nint(LimitPar2Max))
            DumC3=Real2Char(LimitPar2Penalty,fmt=FMTREAL2CHAR)
            write(STDOUT,"(6a)") "LimitFemaleParentContributions: yes, min ",&
              trim(adjustl(DumC)),", max ",trim(adjustl(DumC)),", penalty ",trim(adjustl(DumC3))
          end if
        else if (ToLower(trim(DumC)) == "no") then
          write(STDOUT,"(a)") "LimitFemaleParentContributions: no"
        else
          write(STDERR,"(a)") "ERROR: LimitFemaleParentContributions must be: Yes or No!"
          write(STDERR,"(a)") " "
          stop 1
        end if
      end if

      ! AllowSelfing
      read(UnitSpec,*) DumC,DumC
      if      (ToLower(trim(DumC)) == "yes") then
        SelfingAllowed=.true.
        if (.not.GenderMatters) then
          write(STDOUT,"(a)") "AllowSelfing: Yes"
        else
          write(STDOUT,"(a)") "NOTE: When gender matters, AlphaMate can not perform selfing! See the manual for a solution."
          write(STDOUT,"(a)") " "
        end if
      else if (ToLower(trim(DumC)) == "no") then
        SelfingAllowed=.false.
        backspace(UnitSpec)
        read(UnitSpec,*) DumC,DumC,SelfingPenalty
        DumC=Real2Char(SelfingPenalty,fmt=FMTREAL2CHAR)
        write(STDOUT,"(2a)") "AllowSelfing: no, penalty ",trim(adjustl(DumC))
      else
        write(STDERR,"(a)") "ERROR: AllowSelfing must be: Yes or No!"
        write(STDERR,"(a)") " "
        stop 1
      end if

      ! PAGE
      read(UnitSpec,*) DumC,DumC
      if (.not.GenderMatters) then
        if      (ToLower(trim(DumC)) == "yes") then
          PAGEPar1=.true.
          if (.not.BvAvailable) then
            write(STDERR,"(a)") "ERROR: No breeding values file available!"
            write(STDERR,"(a)") " "
            stop 1
          end if
          backspace(UnitSpec)
          read(UnitSpec,*) DumC,DumC,PAGEPar1Max,PAGEPar1Cost
          DumC=Int2Char(PAGEPar1Max)
          DumC2=Real2Char(PAGEPar1Cost,fmt=FMTREAL2CHAR)
          if (PAGEPar1Max <= nPar) then
            write(STDOUT,"(4a)") "PAGE: yes, no. of individuals ",trim(adjustl(DumC)),", cost ",trim(adjustl(DumC2))
          else
            write(STDERR,"(a)") "ERROR: The max number of individuals to edit must not be greater than the total number of parents!"
            DumC=Int2Char(nPar)
            write(STDERR,"(2a)") "ERROR: Number of             parents: ",trim(adjustl(DumC))
            DumC=Int2Char(PAGEPar1Max)
            write(STDERR,"(2a)") "ERROR: Number of individuals to edit: ",trim(adjustl(DumC))
            write(STDERR,"(a)") " "
          end if
        else if (ToLower(trim(DumC)) == "no") then
          PAGEPar1=.false.
          write(STDOUT,"(a)") "PAGE: no"
        else
          write(STDERR,"(a)") "ERROR: PAGE must be: Yes or No!"
          write(STDERR,"(a)") " "
          stop 1
        end if
      end if

      ! PAGEMales
      read(UnitSpec,*) DumC,DumC
      if (GenderMatters) then
        if      (ToLower(trim(DumC)) == "yes") then
          PAGEPar1=.true.
          if (.not.BvAvailable) then
            write(STDERR,"(a)") "ERROR: No breeding values file available!"
            write(STDERR,"(a)") " "
            stop 1
          end if
          backspace(UnitSpec)
          read(UnitSpec,*) DumC,DumC,PAGEPar1Max,PAGEPar1Cost
          DumC=Int2Char(PAGEPar1Max)
          DumC2=Real2Char(PAGEPar1Cost,fmt=FMTREAL2CHAR)
          if (PAGEPar1Max <= nPar1) then
            write(STDOUT,"(4a)") "PAGEMales: yes, no. of individuals ",trim(adjustl(DumC)),", cost ",trim(adjustl(DumC2))
          else
            write(STDERR,"(a)") "ERROR: The max number of male individuals to edit must not be greater than the total number of male parents!"
            DumC=Int2Char(nPar1)
            write(STDERR,"(2a)") "ERROR: Number of male             parents: ",trim(adjustl(DumC))
            DumC=Int2Char(PAGEPar1Max)
            write(STDERR,"(2a)") "ERROR: Number of male individuals to edit: ",trim(adjustl(DumC))
            write(STDERR,"(a)") " "
          end if
        else if (ToLower(trim(DumC)) == "no") then
          PAGEPar1=.false.
          write(STDOUT,"(a)") "PAGEMales: no"
        else
          write(STDERR,"(a)") "ERROR: PAGEMales must be: Yes or No!"
          write(STDERR,"(a)") " "
          stop 1
        end if
      end if

      ! PAGEFemales
      read(UnitSpec,*) DumC,DumC
      if (GenderMatters) then
        if      (ToLower(trim(DumC)) == "yes") then
          PAGEPar2=.false.
          if (.not.BvAvailable) then
            write(STDERR,"(a)") "ERROR: No breeding values file available!"
            write(STDERR,"(a)") " "
            stop 1
          end if
          backspace(UnitSpec)
          read(UnitSpec,*) DumC,DumC,PAGEPar2Max,PAGEPar2Cost
          DumC=Int2Char(PAGEPar2Max)
          DumC2=Real2Char(PAGEPar2Cost,fmt=FMTREAL2CHAR)
          if (PAGEPar2Max <= nPar2) then
            write(STDOUT,"(4a)") "PAGEFemales: yes, no. of individuals ",trim(adjustl(DumC)),", cost ",trim(adjustl(DumC2))
          else
            write(STDERR,"(a)") "ERROR: The max number of female individuals to edit must not be greater than the total number of female parents!"
            DumC=Int2Char(nPar2)
            write(STDERR,"(2a)") "ERROR: Number of female             parents: ",trim(adjustl(DumC))
            DumC=Int2Char(PAGEPar2Max)
            write(STDERR,"(2a)") "ERROR: Number of female individuals to edit: ",trim(adjustl(DumC))
            write(STDERR,"(a)") " "
          end if
        else if (ToLower(trim(DumC)) == "no") then
          PAGEPar2=.false.
          write(STDOUT,"(a)") "PAGEFemales: no"
        else
          write(STDERR,"(a)") "ERROR: PAGEFemales must be: Yes or No!"
          write(STDERR,"(a)") " "
          stop 1
        end if
      end if

      if (PAGEPar1 .or. PAGEPar2) then
        PAGE=.true.
      else
        PAGE=.false.
      end if

      ! OldCoancestry
      read(UnitSpec,*) DumC,DumC
      if (ToLower(trim(DumC)) == "unknown") then
        InferPopInbOld=.true.
        write(STDOUT,"(a)") "OldCoancestry: unknown"
      else
        InferPopInbOld=.false.
        backspace(UnitSpec)
        read(UnitSpec,*) DumC,PopInbOld
        DumC=Real2Char(PopInbOld,fmt=FMTREAL2CHAR)
        write(STDOUT,"(2a)") "OldCoancestry: ",trim(adjustl(DumC))
      end if

      ! TargetedRateOfPopulationInbreeding
      read(UnitSpec,*) DumC,RatePopInbTarget,PopInbPenalty,DumC
      if      (ToLower(trim(DumC)) == "above") then
        PopInbPenaltyBellow=.false.
      else if (ToLower(trim(DumC)) == "aboveandbellow") then
        PopInbPenaltyBellow=.true.
      else
        write(STDERR,"(a)") "ERROR: PopInbPenaltyMode must be: Above or AboveAndBellow!"
        write(STDERR,"(a)") " "
        stop 1
      end if
      DumC2=Real2Char(RatePopInbTarget,fmt=FMTREAL2CHAR)
      DumC3=Real2Char(PopInbPenalty,fmt=FMTREAL2CHAR)
      write(STDOUT,"(6a)") "TargetedRateOfPopulationInbreeding: ",trim(adjustl(DumC2)),", penalty ",trim(adjustl(DumC3)), ", mode "//trim(adjustl(DumC))

      ! ProgenyInbreedingPenalty
      read(UnitSpec,*) DumC,PrgInbPenalty
      DumC=Real2Char(PrgInbPenalty,fmt=FMTREAL2CHAR)
      write(STDOUT,"(2a)") "ProgenyInbreedingPenalty: ",trim(adjustl(DumC))

      ! EvaluateFrontier
      read(UnitSpec,*) DumC,DumC
      if      (ToLower(trim(DumC)) == "no") then
        EvaluateFrontier=.false.
        write(STDOUT,"(a)") "EvaluateFrontier: no"
      else if (ToLower(trim(DumC)) == "yes") then
        EvaluateFrontier=.true.
        backspace(UnitSpec)
        read(UnitSpec,*) DumC,DumC,nFrontierSteps
        allocate(RatePopInbFrontier(nFrontierSteps))
        backspace(UnitSpec)
        read(UnitSpec,*) DumC,DumC,nFrontierSteps,RatePopInbFrontier(:)
        DumC=Int2Char(nFrontierSteps)
        write(STDOUT,"(2a)") "EvaluateFrontier: yes, #steps: ",trim(adjustl(DumC))
      else
        write(STDERR,"(a)") "ERROR: EvaluateFrontier must be: Yes or No!"
        write(STDERR,"(a)") " "
        stop 1
      end if

      ! EvolutionaryAlgorithmIterations
      read(UnitSpec,*) DumC,EvolAlgNSol,EvolAlgNGen,EvolAlgNGenBurnIn,EvolAlgNGenStop,EvolAlgStopTol,EvolAlgNGenPrint

      ! EvolutionaryAlgorithmParameters
      read(UnitSpec,*) DumC,EvolAlgCRBurnIn,EvolAlgCRLate,EvolAlgFBase,EvolAlgFHigh1,EvolAlgFHigh2

      ! Seed
      read(UnitSpec,*) DumC,DumC
      SeedFile="AlphaMateResults"//DASH//"Seed.txt"
      if ((ToLower(trim(DumC)) == "unknown") .or. (ToLower(trim(DumC)) == "none")) then
        call SetSeed(SeedFile=SeedFile,Out=Seed)
      else
        backspace(UnitSpec)
        read(UnitSpec,*) DumC,DumI
        call SetSeed(Seed=DumI,SeedFile=SeedFile,Out=Seed)
      end if
      DumC=Int2Char(Seed)
      write(STDOUT,"(2a)") "Seed: ",trim(adjustl(DumC))
      write(STDOUT,"(a)") " "

      close(UnitSpec)

      allocate(Bv(nInd))
      allocate(BvStand(nInd))
      if (PAGE) then
        allocate(BvPAGE(nInd))
        allocate(BvPAGEStand(nInd))
        allocate(GeneEdit(nInd))
      end if
      allocate(RelMtx(nInd,nInd))
      allocate(IdC(nInd))
      allocate(Gender(nInd))
      allocate(xVec(nInd))
      allocate(nVec(nInd))
      allocate(MateAlloc(2,nMat))

      write(STDOUT,"(a)") "--- Data ---"
      write(STDOUT,"(a)") " "

      ! --- Relationships ---

      open(newunit=UnitRelMtx,file=trim(RelMtxFile),status="old")
      do i=1,nInd
        read(UnitRelMtx,*) IdC(i),RelMtx(:,i)
      end do
      close(UnitRelMtx)

      RelDescStat=CalcDescStatSymMatrix(RelMtx)
      write(STDOUT,"(a)") "Relationships"
      write(STDOUT,"(a)") "  - self-relationships (diagonal)"
      DumC=Real2Char(RelDescStat%Diag%Mean,fmt=FMTREAL2CHAR)
      write(STDOUT,"(2a)") "    - average: ",trim(adjustl(DumC))
      DumC=Real2Char(RelDescStat%Diag%SD,fmt=FMTREAL2CHAR)
      write(STDOUT,"(2a)") "    - st.dev.: ",trim(adjustl(DumC))
      DumC=Real2Char(RelDescStat%Diag%Min,fmt=FMTREAL2CHAR)
      write(STDOUT,"(2a)") "    - minimum: ",trim(adjustl(DumC))
      DumC=Real2Char(RelDescStat%Diag%Max,fmt=FMTREAL2CHAR)
      write(STDOUT,"(2a)") "    - maximum: ",trim(adjustl(DumC))
      write(STDOUT,"(a)") "  - co-relationships (off-diagonal)"
      DumC=Real2Char(RelDescStat%OffDiag%Mean,fmt=FMTREAL2CHAR)
      write(STDOUT,"(2a)") "    - average: ",trim(adjustl(DumC))
      DumC=Real2Char(RelDescStat%OffDiag%SD,fmt=FMTREAL2CHAR)
      write(STDOUT,"(2a)") "    - st.dev.: ",trim(adjustl(DumC))
      DumC=Real2Char(RelDescStat%OffDiag%Min,fmt=FMTREAL2CHAR)
      write(STDOUT,"(2a)") "    - minimum: ",trim(adjustl(DumC))
      DumC=Real2Char(RelDescStat%OffDiag%Max,fmt=FMTREAL2CHAR)
      write(STDOUT,"(2a)") "    - maximum: ",trim(adjustl(DumC))
      write(STDOUT,"(a)") " "

      ! --- Breeding values ---

      if (BvAvailable) then
        open(newunit=UnitBv,file=trim(BvFile),status="old")
        do i=1,nInd
          if (.not.PAGE) then
            read(UnitBv,*) IdCTmp,BvTmp
          else
            read(UnitBv,*) IdCTmp,BvTmp,BvTmp2
          end if
          do j=1,nInd
            if (trim(IdCTmp) == trim(IdC(j))) then
              Bv(j)=BvTmp
              if (PAGE) then
                BvPAGE(j)=BvTmp2
              end if
              exit
            end if
          end do
        end do
        close(UnitBv)

        BvDescStat=CalcDescStat(Bv)
        BvStand(:)=(Bv(:)-BvDescStat%Mean)/BvDescStat%SD
        write(STDOUT,"(a)") "Breeding values"
        DumC=Real2Char(BvDescStat%Mean,fmt=FMTREAL2CHAR)
        write(STDOUT,"(2a)") "  - average: ",trim(adjustl(DumC))
        DumC=Real2Char(BvDescStat%SD,fmt=FMTREAL2CHAR)
        write(STDOUT,"(2a)") "  - st.dev.: ",trim(adjustl(DumC))
        DumC=Real2Char(BvDescStat%Min,fmt=FMTREAL2CHAR)
        write(STDOUT,"(2a)") "  - minimum: ",trim(adjustl(DumC))
        DumC=Real2Char(BvDescStat%Max,fmt=FMTREAL2CHAR)
        write(STDOUT,"(2a)") "  - maximum: ",trim(adjustl(DumC))
        write(STDOUT,"(a)") " "

        if (PAGE) then
          ! must have the same scaling!!!!
          BvPAGEStand(:)=(BvPAGE(:)-BvDescStat%Mean)/BvDescStat%SD
          ! only the PAGE bit of Bv
          BvPAGE(:)=BvPAGE(:)-Bv(:)
          BvPAGEStand(:)=BvPAGEStand(:)-BvStand(:)
          BvDescStat=CalcDescStat(BvPAGE)
          write(STDOUT,"(a)") "Gene edit increments"
          DumC=Real2Char(BvDescStat%Mean,fmt=FMTREAL2CHAR)
          write(STDOUT,"(2a)") "  - average: ",trim(adjustl(DumC))
          DumC=Real2Char(BvDescStat%SD,fmt=FMTREAL2CHAR)
          write(STDOUT,"(2a)") "  - st.dev.: ",trim(adjustl(DumC))
          DumC=Real2Char(BvDescStat%Min,fmt=FMTREAL2CHAR)
          write(STDOUT,"(2a)") "  - minimum: ",trim(adjustl(DumC))
          DumC=Real2Char(BvDescStat%Max,fmt=FMTREAL2CHAR)
          write(STDOUT,"(2a)") "  - maximum: ",trim(adjustl(DumC))
          write(STDOUT,"(a)") " "
        end if
      else
        Bv(:)=0.0d0
        BvStand(:)=0.0d0
      end if

      ! --- Gender ---

      if (.not.GenderMatters) then
        Gender(:)=0
      else
        nMal=0
        nFem=0
        open(newunit=UnitGender,file=trim(GenderFile),status="old")

        do i=1,nInd
          read(UnitGender,*) IdCTmp,GenderTmp
          if (GenderTmp == 1) then
            nMal=nMal+1
          elseif (GenderTmp == 2) then
            nFem=nFem+1
          else
            write(STDERR,"(a)") "ERROR: Gender code must be either 1 for male individuals or 2 for female individuals!"
            write(STDERR,"(a,i6,a,i3)") "ERROR: ",i,IdCTmp,GenderTmp
            write(STDERR,"(a)") " "
            stop 1
          end if
          do j=1,nInd
            if (trim(IdCTmp) == trim(IdC(j))) then
              Gender(j)=GenderTmp
              exit
            end if
          end do
        end do
        close(UnitGender)
        DumC=Int2Char(nMal)
        write(STDOUT,"(2a)") "Number of   males in data: ",trim(adjustl(DumC))
        DumC=Int2Char(nFem)
        write(STDOUT,"(2a)") "Number of females in data: ",trim(adjustl(DumC))
        write(STDOUT,"(a)") " "
        if (nPar1 > nMat) then
          write(STDERR,"(a)") "ERROR: The number of male parents can not be larger than the number of males"
          DumC=Int2Char(nPar1)
          write(STDERR,"(2a)") "ERROR: Number of male parents: ",trim(adjustl(DumC))
          DumC=Int2Char(nMal)
          write(STDERR,"(2a)") "ERROR: Number of        males: ",trim(adjustl(DumC))
          write(STDERR,"(a)") " "
          stop 1
        end if
        if (nPar2 > nFem) then
          write(STDERR,"(a)") "ERROR: The number of female parents can not be larger than the number of females"
          DumC=Int2Char(nPar2)
          write(STDERR,"(2a)") "ERROR: Number of female parents: ",trim(adjustl(DumC))
          DumC=Int2Char(nFem)
          write(STDERR,"(2a)") "ERROR: Number of        females: ",trim(adjustl(DumC))
          write(STDERR,"(a)") " "
          stop 1
        end if
      end if

      ! --- Shuffle the data ---

      ! To avoid having good animals together - better for cross-overs in Evolutionary algorithms

      allocate(Order(nInd))
      Order=RandomOrder(nInd)
      IdC(:)=IdC(Order)
      Bv(:)=Bv(Order)
      BvStand(:)=BvStand(Order)
      if (PAGE) then
        BvPAGE(:)=BvPAGE(Order)
        BvPAGEStand(:)=BvPAGEStand(Order)
      end if
      Gender(:)=Gender(Order)
      RelMtx(:,:)=RelMtx(Order,Order)
      deallocate(Order)

      ! --- Define ---

      if (.not.GenderMatters) then
        nPotPar1=nInd
        allocate(IdPotPar1(nPotPar1))
        allocate(nVecPar1(nPotPar1))
        do i=1,nInd
          IdPotPar1(i)=i
        end do
      else
        nPotPar1=nMal
        nPotPar2=nFem
        allocate(IdPotPar1(nPotPar1))
        allocate(IdPotPar2(nPotPar2))
        allocate(nVecPar1(nPotPar1))
        allocate(nVecPar2(nPotPar2))
        jMal=0
        jFem=0
        do i=1,nInd
          if (Gender(i) == 1) then
            jMal=jMal+1
            IdPotPar1(jMal)=i
          else
            jFem=jFem+1
            IdPotPar2(jFem)=i
          end if
        end do
        if (nPar1 > nPotPar1) then
          nPar1=nPotPar1
          write(STDOUT,"(a)") "NOTE: The number of male parents reduced to the number of male individuals!"
          write(STDOUT,"(a)") " "
        end if
        if (nPar2 > nPotPar2) then
          nPar2=nPotPar2
          write(STDOUT,"(a)") "NOTE: The number of female parents reduced to the number of female individuals!"
          write(STDOUT,"(a)") " "
        end if
      end if
    end subroutine

    !###########################################################################

    subroutine SetInbreedingParameters

      implicit none

      integer(int32) :: i,UnitInbree

      real(real64) :: Tmp

      character(len=100) :: DumC

      ! Old inbreeding
      if (InferPopInbOld) then
        PopInbOld=0.0d0
        do i=1,nInd
          Tmp=RelMtx(i,i)-1.0d0
          if (Tmp < 0.0d0) then
            write(STDERR,"(a)") "ERROR: Relationship matrix must have diagonals equal or more than 1.0!"
            write(STDERR,"(a)") " "
            stop 1
          end if
          PopInbOld=PopInbOld+Tmp
        end do
        PopInbOld=PopInbOld/dble(nInd)
      end if

      ! Targeted level of inbreeding
      ! F_t = DeltaF + (1 - DeltaF) * F_t-1
      PopInbTarget=RatePopInbTarget+(1.0d0-RatePopInbTarget)*PopInbOld

      ! Report
      DumC=Real2Char(PopInbOld,fmt=FMTREAL2CHAR)
      write(STDOUT,"(2a)") "Old coancestry: ",trim(adjustl(DumC))
      DumC=Real2Char(RatePopInbTarget,fmt=FMTREAL2CHAR)
      write(STDOUT,"(2a)") "Targeted rate of population inbreeding: ",trim(adjustl(DumC))
      DumC=Real2Char(PopInbTarget,fmt=FMTREAL2CHAR)
      write(STDOUT,"(2a)") "Targeted population inbreeding: ",trim(adjustl(DumC))
      write(STDOUT,"(a)") " "

      open(newunit=UnitInbree,file="AlphaMateResults"//DASH//"ConstraintPopulationInbreeding.txt",status="unknown")
      write(UnitInbree,"(a,f)") "Old_coancestry_defined, ",PopInbOld
      write(UnitInbree,"(a,f)") "Targeted_rate_of_population_inbreeding_defined, ",RatePopInbTarget
      write(UnitInbree,"(a,f)") "Targeted_population_inbreeding_defined, ",PopInbTarget
      close(UnitInbree)
    end subroutine

    !###########################################################################

    subroutine AlphaMateSearch

      implicit none

      integer(int32) :: i,j,k,nTmp,DumI,Rank(nInd)
      integer(int32) :: UnitInbree,UnitMating,UnitContri,UnitLog,UnitLog2,UnitFrontier

      real(real64) :: PopInbTargetHold,RatePopInbTargetHold,DumR(8)
      real(real64),allocatable :: InitEqual(:,:)

      character(len=1000) :: EvolAlgLogFile,EvolAlgLogFile2
      character(len=100) :: DumC,DumC2,DumC3,DumC4

      logical :: Success

      type(EvolveCrit) :: CritMin,CritOpt,Crit

      ! --- Optimise for minimum inbreeding ---

      if (ModeMin) then
        write(STDOUT,"(a)") "--- Optimise for minimum inbreeding --- "
        write(STDOUT,"(a)") " "

        EvolAlgLogFile="AlphaMateResults"//DASH//"OptimisationLogMinimumInbreeding.txt"
        EvolAlgLogFile2="AlphaMateResults"//DASH//"OptimisationLogMinimumInbreedingInitial.txt"
        if (GenderMatters) then
          nTmp=nPotPar1+nPotPar2+nMat
        else
          nTmp=nPotPar1+nMat
        end if
        if (PAGE) then
          nTmp=nTmp+nInd
        end if
        allocate(InitEqual(nTmp,nint(real(EvolAlgNSol*0.1))))
        InitEqual(:,:)=1.0d0 ! A couple of solutions that would give equal contribution and the rest for everybody
        call EvolAlgDE(nParam=nTmp,nSol=EvolAlgNSol,Init=InitEqual,nGen=EvolAlgNGen,nGenBurnIn=EvolAlgNGenBurnIn,&
                       nGenStop=EvolAlgNGenStop,StopTolerance=EvolAlgStopTol,&
                       nGenPrint=EvolAlgNGenPrint,File=EvolAlgLogFile,CritType="Min",&
                       CRBurnIn=EvolAlgCRBurnIn,CRLate=EvolAlgCRLate,&
                       FBase=EvolAlgFBase,FHigh1=EvolAlgFHigh1,FHigh2=EvolAlgFHigh2,&
                       CalcCriterion=FixSolEtcMateAndCalcCrit,&
                       LogHeader=EvolAlgLogHeaderForAlphaMate,Log=EvolAlgLogForAlphaMate,&
                       BestCriterion=CritMin)
        deallocate(InitEqual)
        GainMinStand=CritMin%GainStand

        open(newunit=UnitContri,file="AlphaMateResults"//DASH//"IndividualResultsMinimumInbreeding.txt",status="unknown")
        Rank=MrgRnk(nVec)
        !                             1234567890123456789012
        if (.not.PAGE) then
          write(UnitContri,FMTINDHEAD) "          Id",&
                                       "      Gender",&
                                       "       Merit",&
                                       " AvgCoancest",&
                                       "  Contribute",&
                                       "    nMatings"
          do i=nInd,1,-1 ! MrgRnk ranks small to large
            j=Rank(i)
            write(UnitContri,FMTIND) IdC(j),Gender(j),Bv(j),0.5d0*sum(RelMtx(:,j))/dble(nInd),xVec(j),nVec(j)
          end do
        else
          !                                 1234567890123456789012
          write(UnitContri,FMTINDHEADEDIT) "          Id",&
                                           "      Gender",&
                                           "       Merit",&
                                           " AvgCoancest",&
                                           "  Contribute",&
                                           "    nMatings",&
                                           "    GeneEdit",&
                                           " EditedMerit"
          do i=nInd,1,-1 ! MrgRnk ranks small to large
            j=Rank(i)
            write(UnitContri,FMTINDEDIT) IdC(j),Gender(j),Bv(j),0.5d0*sum(RelMtx(:,j))/dble(nInd),xVec(j),nVec(j),0,Bv(j)
          end do
        end if
        close(UnitContri)

        open(newunit=UnitMating,file="AlphaMateResults"//DASH//"MatingResultsMinimumInbreeding.txt",status="unknown")
        !                             1234567890123456789012
        write(UnitMating,FMTMATHEAD) "      Mating",&
                                     "     Parent1",&
                                     "     Parent2"
        do i=1,nMat
          write(UnitMating,FMTMAT) i,IdC(MateAlloc(1,i)),IdC(MateAlloc(2,i))
        end do
        close(UnitMating)

        if (PopInbOld > CritMin%PopInb) then

          write(STDOUT,"(a)") "NOTE: Old coancestry is higher than the minimum group coancestry (x'Ax/2) under no selection."
          write(STDOUT,"(a)") "NOTE: Resetting the old coancestry to the minimum group coancestry under no selection and"
          write(STDOUT,"(a)") "NOTE:   recomputing the log file values and the targeted population inbreeding."
          write(STDOUT,"(a)") " "
          PopInbOld=CritMin%PopInb
          ! F_t = DeltaF + (1 - DeltaF) * F_t-1
          PopInbTarget=RatePopInbTarget+(1.0d0-RatePopInbTarget)*PopInbOld
          CritMin%RatePopInb=0.0d0

          Success=SystemQQ(COPY//" "//EvolAlgLogFile//" "//EvolAlgLogFile2)
          if (.not.Success) then
            write(STDERR,"(a)") "ERROR: Failed to make a backup of the "//EvolAlgLogFile//" file in the output folder!"
            write(STDERR,"(a)") " "
            stop 1
          end if

          open(newunit=UnitLog,file=trim(EvolAlgLogFile),status="unknown")
          open(newunit=UnitLog2,file=trim(EvolAlgLogFile2),status="unknown")
          nTmp=CountLines(EvolAlgLogFile2)
          read(UnitLog2,*) DumC
          call EvolAlgLogHeaderForAlphaMate(UnitLog)
          do i=2,nTmp
            read(UnitLog2,*) DumI,DumR(:)
            ! F_t = DeltaF + (1 - DeltaF) * F_t-1
            ! DeltaF = (F_t - F_t-1) / (1 - F_t-1)
            DumR(7)=(DumR(6)-PopInbOld)/(1.0d0-PopInbOld)
            write(STDOUT, FMTLOGSTDOUT) DumI,DumR(:)
            write(UnitLog,FMTLOGUNIT)   DumI,DumR(:)
          end do
          write(STDOUT,"(a)") " "
          close(UnitLog)
          close(UnitLog2)

          DumC=Real2Char(PopInbOld,fmt=FMTREAL2CHAR)
          write(STDOUT,"(2a)") "Old coancestry: ",trim(adjustl(DumC))
          DumC=Real2Char(RatePopInbTarget,fmt=FMTREAL2CHAR)
          write(STDOUT,"(2a)") "Targeted rate of population inbreeding: ",trim(adjustl(DumC))
          DumC=Real2Char(PopInbTarget,fmt=FMTREAL2CHAR)
          write(STDOUT,"(2a)") "Targeted population inbreeding: ",trim(adjustl(DumC))
          write(STDOUT,"(a)") " "

        end if

        if (PopInbTarget < CritMin%PopInb) then
          write(STDERR,"(a)") "ERROR: Targeted population inbreeding is lower than the group coancestry (x'Ax/2) under no selection."
          write(STDERR,"(a)") "ERROR: Can not optimise!"
          write(STDERR,"(a)") " "
          stop 1
        end if

        open(newunit=UnitInbree,file="AlphaMateResults"//DASH//"ConstraintPopulationInbreeding.txt",status="old")
        write(UnitInbree,"(a,f)") "Old_coancestry_redefined, ",PopInbOld
        write(UnitInbree,"(a,f)") "Targeted_rate_of_population_inbreeding_redefined, ",RatePopInbTarget
        write(UnitInbree,"(a,f)") "Targeted_population_inbreeding_redefined, ",PopInbTarget
        close(UnitInbree)
      else
        GainMinStand=0.0d0
      end if

      ! --- Optimise for maximum gain with constraint on inbreeding ---

      if (ModeOpt) then
        write(STDOUT,"(a)") "--- Optimise for maximum gain with constraint on inbreeding ---"
        write(STDOUT,"(a)") " "

        EvolAlgLogFile="AlphaMateResults"//DASH//"OptimisationLogOptimumGain.txt"
        if (GenderMatters) then
          nTmp=nPotPar1+nPotPar2+nMat
        else
          nTmp=nPotPar1+nMat
        end if
        if (PAGE) then
          nTmp=nTmp+nInd
        end if
        call EvolAlgDE(nParam=nTmp,nSol=EvolAlgNSol,nGen=EvolAlgNGen,nGenBurnIn=EvolAlgNGenBurnIn,&
                       nGenStop=EvolAlgNGenStop,StopTolerance=EvolAlgStopTol,&
                       nGenPrint=EvolAlgNGenPrint,File=EvolAlgLogFile,CritType="Opt",&
                       CRBurnIn=EvolAlgCRBurnIn,CRLate=EvolAlgCRLate,&
                       FBase=EvolAlgFBase,FHigh1=EvolAlgFHigh1,FHigh2=EvolAlgFHigh2,&
                       CalcCriterion=FixSolEtcMateAndCalcCrit,&
                       LogHeader=EvolAlgLogHeaderForAlphaMate,Log=EvolAlgLogForAlphaMate,&
                       BestCriterion=CritOpt)

        open(newunit=UnitContri,file="AlphaMateResults"//DASH//"IndividualResultsOptimumGain.txt",status="unknown")
        Rank=MrgRnk(nVec)
        !                             1234567890123456789012
        if (.not.PAGE) then
          write(UnitContri,FMTINDHEAD) "          Id",&
                                       "      Gender",&
                                       "       Merit",&
                                       " AvgCoancest",&
                                       "  Contribute",&
                                       "    nMatings"
          do i=nInd,1,-1 ! MrgRnk ranks small to large
            j=Rank(i)
            write(UnitContri,FMTIND) IdC(j),Gender(j),Bv(j),0.5d0*sum(RelMtx(:,j))/dble(nInd),xVec(j),nVec(j)
          end do
        else
          !                                 1234567890123456789012
          write(UnitContri,FMTINDHEADEDIT) "          Id",&
                                           "      Gender",&
                                           "       Merit",&
                                           " AvgCoancest",&
                                           "  Contribute",&
                                           "    nMatings",&
                                           "    GeneEdit",&
                                           " EditedMerit"
          do i=nInd,1,-1 ! MrgRnk ranks small to large
            j=Rank(i)
            write(UnitContri,FMTINDEDIT) IdC(j),Gender(j),Bv(j),0.5d0*sum(RelMtx(:,j))/dble(nInd),xVec(j),nVec(j),nint(GeneEdit(j)),Bv(j)+GeneEdit(j)*BvPAGE(j)
          end do
        end if
        close(UnitContri)

        open(newunit=UnitMating,file="AlphaMateResults"//DASH//"MatingResultsOptimumGain.txt",status="unknown")
        !                             1234567890123456789012
        write(UnitMating,FMTMATHEAD) "      Mating",&
                                     "     Parent1",&
                                     "     Parent2"
        do i=1,nMat
          write(UnitMating,FMTMAT) i,IdC(MateAlloc(1,i)),IdC(MateAlloc(2,i))
        end do
        close(UnitMating)
      end if

      ! --- Evaluate the full frontier ---

      if (EvaluateFrontier) then
        write(STDOUT,"(a)") "--- Evaluate the full frontier (this might take some time!) ---"
        write(STDOUT,"(a)") " "

        PopInbPenaltyBellow=.true. ! we want to target certain rates of inbreeding

        open(newunit=UnitFrontier,file="AlphaMateResults"//DASH//"Frontier.txt",status="unknown")
        !                               1234567890123456789012
        write(UnitFrontier,FMTFROHEAD) "        Step",&
                                       "             Objective",&
                                       "             Penalties",&
                                       "                  Gain",&
                                       "             GainStand",&
                                       "            PopInbreed",&
                                       "            RatePopInb",&
                                       "            PrgInbreed"
        if (ModeMin) then
          DumC="Min"
          write(UnitFrontier,FMTFRO) adjustl(DumC),CritMin%Value,CritMin%Penalty,CritMin%Gain,CritMin%GainStand,CritMin%PopInb,CritMin%RatePopInb,CritMin%PrgInb
        end if
        if (ModeOpt) then
          DumC="Opt"
          write(UnitFrontier,FMTFRO) adjustl(DumC),CritOpt%Value,CritOpt%Penalty,CritOpt%Gain,CritOpt%GainStand,CritOpt%PopInb,CritOpt%RatePopInb,CritOpt%PrgInb
        end if

        ! Hold old results
        PopInbTargetHold=PopInbTarget
        RatePopInbTargetHold=RatePopInbTarget

        ! Evaluate
        do k=1,nFrontierSteps
          RatePopInbTarget=RatePopInbFrontier(k)
          ! F_t = DeltaF + (1 - DeltaF) * F_t-1
          PopInbTarget=RatePopInbTarget+(1.0d0-RatePopInbTarget)*PopInbOld
          DumC=Int2Char(k)
          DumC2=Int2Char(nFrontierSteps)
          DumC3=Real2Char(RatePopInbTarget,fmt=FMTREAL2CHAR)
          DumC4=Real2Char(PopInbTarget,fmt=FMTREAL2CHAR)
          write(STDOUT,"(9a)") "Step ",trim(adjustl(DumC))," out of ",trim(adjustl(DumC2)),&
                               " for the rate of population inbreeding of ",trim(adjustl(DumC3)),&
                               " (=pop. inbreed. of ",trim(adjustl(DumC4)),")"
          write(STDOUT,"(a)") ""
          EvolAlgLogFile="AlphaMateResults"//DASH//"OptimisationLogFrontier"//Int2Char(k)//".txt"
          if (GenderMatters) then
            nTmp=nPotPar1+nPotPar2+nMat
          else
            nTmp=nPotPar1+nMat
          end if
          if (PAGE) then
            nTmp=nTmp+nInd
          end if
          call EvolAlgDE(nParam=nTmp,nSol=EvolAlgNSol,nGen=EvolAlgNGen,nGenBurnIn=EvolAlgNGenBurnIn,&
                         nGenStop=EvolAlgNGenStop,StopTolerance=EvolAlgStopTol,&
                         nGenPrint=EvolAlgNGenPrint,File=EvolAlgLogFile,CritType="Opt",&
                         CRBurnIn=EvolAlgCRBurnIn,CRLate=EvolAlgCRLate,&
                         FBase=EvolAlgFBase,FHigh1=EvolAlgFHigh1,FHigh2=EvolAlgFHigh2,&
                         CalcCriterion=FixSolEtcMateAndCalcCrit,&
                         LogHeader=EvolAlgLogHeaderForAlphaMate,Log=EvolAlgLogForAlphaMate,&
                         BestCriterion=Crit)
          DumC="Frontier"//Int2Char(k)
          write(UnitFrontier,FMTFRO) adjustl(DumC),Crit%Value,Crit%Penalty,Crit%Gain,Crit%GainStand,Crit%PopInb,Crit%RatePopInb,Crit%PrgInb

          open(newunit=UnitContri,file="AlphaMateResults"//DASH//"IndividualResultsFrontier"//Int2Char(k)//".txt",status="unknown")
          Rank=MrgRnk(nVec)
          !                             1234567890123456789012
          if (.not.PAGE) then
            write(UnitContri,FMTINDHEAD) "          Id",&
                                         "      Gender",&
                                         "       Merit",&
                                         " AvgCoancest",&
                                         "  Contribute",&
                                         "    nMatings"
            do i=nInd,1,-1 ! MrgRnk ranks small to large
              j=Rank(i)
              write(UnitContri,FMTIND) IdC(j),Gender(j),Bv(j),0.5d0*sum(RelMtx(:,j))/dble(nInd),xVec(j),nVec(j)
            end do
          else
            !                                 1234567890123456789012
            write(UnitContri,FMTINDHEADEDIT) "          Id",&
                                             "      Gender",&
                                             "       Merit",&
                                             " AvgCoancest",&
                                             "  Contribute",&
                                             "    nMatings",&
                                             "    GeneEdit",&
                                             " EditedMerit"
            do i=nInd,1,-1 ! MrgRnk ranks small to large
              j=Rank(i)
              write(UnitContri,FMTINDEDIT) IdC(j),Gender(j),Bv(j),0.5d0*sum(RelMtx(:,j))/dble(nInd),xVec(j),nVec(j),nint(GeneEdit(j)),Bv(j)+GeneEdit(j)*BvPAGE(j)
            end do
          end if
          close(UnitContri)

          open(newunit=UnitMating,file="AlphaMateResults"//DASH//"MatingResultsFrontier"//Int2Char(k)//".txt",status="unknown")
          !                             1234567890123456789012
          write(UnitMating,FMTMATHEAD) "      Mating",&
                                       "     Parent1",&
                                       "     Parent2"
          do i=1,nMat
            write(UnitMating,FMTMAT) i,IdC(MateAlloc(1,i)),IdC(MateAlloc(2,i))
          end do
          close(UnitMating)

          if ((RatePopInbTarget-Crit%RatePopInb) > 0.01d0) then
            DumC=Real2Char(RatePopInbTarget,fmt=FMTREAL2CHAR)
            write(STDOUT,"(2a)") "NOTE: Could not achieve the rate of population inbreeding of ",trim(adjustl(DumC))
            write(STDOUT,"(a)") "NOTE: Stopping the frontier evaluation."
            write(STDOUT,"(a)") ""
            exit
          end if
        end do

        ! Put back old results
        PopInbTarget=PopInbTargetHold
        RatePopInbTarget=RatePopInbTargetHold

        close(UnitFrontier)

      end if
    end subroutine

    !###########################################################################

    subroutine EvolAlgLogHeaderForAlphaMate(LogUnit)
      implicit none
      integer(int32),intent(in) :: LogUnit
      ! With ifort this way of formating the code was requied to get desired layout
      ! in STDOUT and in the file. Odd.
      write(STDOUT,FMTLOGHEADERSTDOUT) "Step","AcceptRate","Criterion","Penalties","Gain","GainStand","PopInbreed","RatePopInb","PrgInbreed"
      !                                 1234567890123456789012
      write(LogUnit,FMTLOGHEADERUNIT)  "        Step",&
                                       "            AcceptRate",&
                                       "             Criterion",&
                                       "             Penalties",&
                                       "                  Gain",&
                                       "             GainStand",&
                                       "            PopInbreed",&
                                       "            RatePopInb",&
                                       "            PrgInbreed"
    end subroutine

    !###########################################################################

    subroutine EvolAlgLogForAlphaMate(LogUnit,Gen,AcceptRate,Criterion)
      implicit none
      integer(int32),intent(in)   :: LogUnit
      integer(int32),intent(in)   :: Gen
      real(real64),intent(in)     :: AcceptRate
      type(EvolveCrit),intent(in) :: Criterion
      write(STDOUT, FMTLOGSTDOUT) Gen,AcceptRate,Criterion%Value,Criterion%Penalty,Criterion%Gain,Criterion%GainStand,Criterion%PopInb,Criterion%RatePopInb,Criterion%PrgInb
      write(LogUnit,FMTLOGUNIT)   Gen,AcceptRate,Criterion%Value,Criterion%Penalty,Criterion%Gain,Criterion%GainStand,Criterion%PopInb,Criterion%RatePopInb,Criterion%PrgInb
    end subroutine

    !###########################################################################

    subroutine InitialiseAlphaMateCrit(This,Value)
      implicit none
      type(EvolveCrit),intent(out)     :: This
      real(real64),intent(in),optional :: Value
      if (present(Value)) then
        This%Value=Value
      else
        This%Value=0.0d0
      end if
      This%Penalty=0.0d0
      This%Gain=0.0d0
      This%GainStand=0.0d0
      This%PopInb=0.0d0
      This%RatePopInb=0.0d0
      This%PrgInb=0.0d0
    end subroutine

    !###########################################################################

    subroutine FixSolEtcMateAndCalcCrit(Sol,CritType,Criterion)

      implicit none

      ! Arguments
      real(real64),intent(inout)   :: Sol(:)    ! Solution
      character(len=*),intent(in)  :: CritType  ! Type of criterion (Min,Opt)
      type(EvolveCrit),intent(out) :: Criterion ! Criterion of the solution

      ! Other
      integer(int32) :: i,j,k,l,g,nCumMat,RankSol(nInd),SolInt(nInd),MatPar2(nMat)
      integer(int32) :: TmpMin,TmpMax,TmpI

      real(real64) :: TmpVec(nInd,1),TmpR!,RanNum

      ! Criterion base
      call InitialiseAlphaMateCrit(Criterion)

      ! The solution (based on the mate selection driver) has:
      ! - nInd individual contributions
      !   - nPotPar1 individual contributions for "parent1" (males   when GenderMatters, all ind when .not. GenderMatters)
      !   - nPotPar2 individual contributions for "parent2" (females when GenderMatters, present only when GenderMatters)
      ! - nMat     rankings of parent1 1:nMat matings to pair with 1:nPotPar2 "parent2" (see bellow)
      ! - nInd edit indicators
      !   - nPotPar1 edit indicators for "parent1" (males   when GenderMatters, all ind when .not. GenderMatters)
      !   - nPotPar2 edit indicators for "parent2" (females when GenderMatters, present only when GenderMatters)

      ! Say we have Sol=(|0,2,0,1|...|2.5,1.5,1.0|0,1,0,0|...) then we:
      ! - mate male 2 with the first  available female (rank 2.5)
      ! - mate male 2 with the second available female (rank 1.5)
      ! - mate male 4 with the third  available female (rank 1.0)
      ! - edit male 2

      ! TODO: consider spliting the Sol() vector internally into a type with
      !       separate vectors to simplify the code, e.g.,
      ! Sol2%ContPar1
      ! Sol2%ContPar2
      ! Sol2%MateRank
      ! Sol2%EditPar1
      ! Sol2%EditPar2

      ! --- Parse the mate selection driver (=Is the solution valid?) ---

      ! The approach bellow assures that we have nMat contributions for each of
      ! the two parent sets. It does this by ranking internal solution values and
      ! traverses from top to the defined number of parents checking when the sum of
      ! interegrised values gives nMat. If values bellow 0.5 are found, they are
      ! changed to 1 contribution. If this still does not give nMat, then we start
      ! adding on contribution to each parent (starting from least contributing
      ! parents to avoid local minima) until we reach nMat. How to treat values
      ! for the individuals that do not contribute is unlcear. None of the tested
      ! methods seemed to be very different. Intuitively, using properly ordered
      ! negative values should inform optim. alg. which individuals should less
      ! likely contribute, but this did not seem to be the case - better final
      ! solution was found when this strategy was not implemented - either zeroing
      ! values for those individuals (was the fastest) or giving random value (was
      ! marginally better, but slower). Potential advantage of not preserving the
      ! order is that this gives more randomness and more solutions being explored.

      ! "Parent1"
      if (GenderMatters) then
        g=1
      else
        g=2
      end if
      ! ... find ranks to find the top values
      if (.not.(EqualizePar1 .and. (nPar1 == nPotPar1))) then
        RankSol(1:nPotPar1)=MrgRnk(Sol(1:nPotPar1))
        RankSol(1:nPotPar1)=RankSol(nPotPar1:1:-1) ! MrgRnk ranks small to large
      end if
      ! ... handle cases with equalized contributions
      if (EqualizePar1) then
        if (nPar1 == nPotPar1) then
          ! ... set integers to all the values (no need for sorting here)
          Sol(1:nPotPar1)=dble(nMat*g)/dble(nPar1)
        else
          ! ... set integers to the top values
          Sol(RankSol(1:nPar1))=dble(nMat*g)/dble(nPar1)
          ! TODO: anything better to preserve the order of non contributing individuals? See below!
          Sol(RankSol((nPar1+1):nPotPar1))=0.0d0
          ! Sol(RankSol((nPar1+1):nPotPar1))=-1.0d0
        end if
      else
        ! ... handle cases with unequal contributions
        ! ... work for the defined number or parents
        nCumMat=0
        do i=1,nPar1
          j=RankSol(i)
          ! ... these would have "zero" contributions and if we hit them then they need at least minimum contrib
          if (Sol(j) < LimitPar1Min) then
            Sol(j)=LimitPar1Min
          end if
          ! ... but not above max allowed
          if (Sol(j) > LimitPar1Max) then
            Sol(j)=LimitPar1Max
          end if
          ! ... accumulate and check if we reached nMat
          nCumMat=nCumMat+nint(Sol(j)) ! internally real, externally integer
          if (nCumMat >= nMat*g) then
            ! ... there should be exactly nMat contributions
            if (nCumMat > nMat*g) then
              Sol(j)=Sol(j)-dble(nCumMat-nMat*g)
              if (nint(Sol(j)) < LimitPar1Min) then
                TmpR=LimitPar1Penalty*(LimitPar1Min-nint(Sol(j)))
                Criterion%Value=Criterion%Value-TmpR
                Criterion%Penalty=Criterion%Penalty+TmpR
              end if
              nCumMat=nMat*g
            end if
            exit
          end if
        end do
        ! ... increment i if we have hit the exit, do loop would have ended with i=nPar1+1
        if (i <= nPar1) then
          i=i+1
        end if
        ! ... the other individuals do not contribute
        if (i <= nPotPar1) then ! "=" to capture incremented i+1 on the do loop exit
          ! ... zero (the same for all ind so no order)
          Sol(RankSol(i:nPotPar1))=0.0d0
          ! ... negative (the same for all ind so no order)
          ! Sol(RankSol(i:nPotPar1))=-1.0d0
          ! ... negative (variable with partially preserving order)
          !     Found faster convergence than with properly decreasing negative values?
          !     I guess it adds some more randomness, i.e., it causes suffling of individuals that just did not make it onto the mating list.
          ! Sol(RankSol(i:nPotPar1))=sign(Sol(RankSol(i:nPotPar1)),-1.0d0)
          ! ... negative and properly decreasing
          ! TmpR=maxval(Sol(RankSol(i:nPotPar1)))
          ! if (TmpR > 0.0d0) then
          !     Sol(RankSol(i:nPotPar1))=Sol(RankSol(i:nPotPar1))-abs(TmpR)
          ! end if
          ! ... negative (random so no order)
          ! do j=i,nPotPar1 ! TODO: really need this loop?
          !   call random_number(RanNum)
          !   Sol(RankSol(j))=-1.0d0*RanNum
          ! end do
        end if
        ! ... nMat still not reached?
        do while (nCumMat < nMat*g)
          ! ... add more contributions
          do i=nPar1,1,-1 ! to bottom ranked selected individuals (to avoid local optima)
            j=RankSol(i)
            Sol(j)=Sol(j)+1.0d0
            ! ... accumulate and check if we reached nMat
            nCumMat=nCumMat+1
            if (nCumMat >= nMat*g) then
              ! To cater for real vs. integer issues
              TmpI=sum(nint(Sol(RankSol(1:nPar1))))
              if (TmpI /= nMat*g) then
                if (TmpI > nMat*g) then
                  Sol(j)=dble(nint(Sol(j))-1)
                else
                  Sol(j)=dble(nint(Sol(j))+1)
                end if
              end if
              exit
            end if
          end do
        end do
      end if

      ! "Parent2"
      if (GenderMatters) then
        ! ... find ranks to find the top values
        if (.not.(EqualizePar2 .and. (nPar2 == nPotPar2))) then
          RankSol(1:nPotPar2)=MrgRnk(Sol((nPotPar1+1):(nPotPar1+nPotPar2)))
          RankSol(1:nPotPar2)=RankSol(nPotPar2:1:-1) ! MrgRnk ranks small to large
        end if
        ! ... handle cases with equalized contributions
        if (EqualizePar2) then
          if (nPar2 == nPotPar2) then
            ! ... set integers to all the values (no need for sorting here)
            Sol((nPotPar1+1):(nPotPar1+nPotPar2))=dble(nMat)/dble(nPar2)
          else
            ! ... set integers to the top values
            Sol(nPotPar1+RankSol(1:nPar2))=dble(nMat)/dble(nPar2)
            ! TODO: anything better to preserve the order of non contributing individuals? See below!
            Sol(nPotPar1+RankSol((nPar2+1):nPotPar2))=0.0d0
            ! Sol(nPotPar1+RankSol((nPar2+1):nPotPar2))=-1.0d0
          end if
        else
          ! ... handle cases with unequal contributions
          ! ... work for the defined number or parents
          nCumMat=0
          do i=1,nPar2
            j=nPotPar1+RankSol(i)
            ! ... these would have "zero" contributions and if we hit them then they need at least minimum contrib
            if (Sol(j) < LimitPar2Min) then
              Sol(j)=LimitPar2Min
            end if
            ! ... but not above max allowed
            if (Sol(j) > LimitPar2Max) then
              Sol(j)=LimitPar2Max
            end if
            ! ... accumulate and check if we reached nMat
            nCumMat=nCumMat+nint(Sol(j)) ! internally real, externally integer
            if (nCumMat >= nMat) then
              ! ... there should be exactly nMat contributions
              if (nCumMat > nMat) then
                Sol(j)=Sol(j)-dble(nCumMat-nMat)
                if (nint(Sol(j)) < LimitPar2Min) then
                  TmpR=LimitPar2Penalty*(LimitPar2Min-nint(Sol(j)))
                  Criterion%Value=Criterion%Value-TmpR
                  Criterion%Penalty=Criterion%Penalty+TmpR
                end if
                nCumMat=nMat
              end if
              exit
            end if
          end do
          ! ... increment i if we have hit the exit, do loop would have ended with i=nPar2+1
          if (i <= nPar2) then
            i=i+1
          end if
          ! ... the other individuals do not contribute
          if (i <= nPotPar2) then ! "="" to capture incremented i+1 on the do loop exit
            ! ... zero (the same for all ind so no order)
            Sol(nPotPar1+(RankSol(i:nPotPar2)))=0.0d0
            ! ... negative (the same for all ind so no order)
            ! Sol(nPotPar1+(RankSol(i:nPotPar2)))=-1.0d0
            ! ... negative (variable with partially preserving order, i.e., ~large positives become ~large negatives)
            !     Found faster convergence than with properly decreasing negative values?
            !     I guess it adds some more randomness, i.e., it causes suffling of individuals that just did not make it onto the mating list.
            ! Sol(nPotPar1+(RankSol(i:nPotPar2)))=sign(Sol(nPotPar1+(RankSol(i:nPotPar2))),-1.0d0)
            ! ... negative and properly decreasing
            ! TmpR=maxval(Sol(nPotPar1+(RankSol(i:nPotPar2))))
            ! if (TmpR > 0.0d0) then
            !     Sol(nPotPar1+(RankSol(i:nPotPar2)))=Sol(nPotPar1+(RankSol(i:nPotPar2)))-abs(TmpR)
            ! end if
            ! ... negative (random so no order)
            ! do j=i,nPotPar2 ! TODO: really need this loop?
            !   call random_number(RanNum)
            !   Sol(nPotPar1+RankSol(j))=-1.0d0*RanNum
            ! end do
          end if
          ! ... nMat still not reached?
          do while (nCumMat < nMat)
            ! ... add more contributions
            do i=nPar2,1,-1 ! to bottom ranked selected individuals (to avoid local optima)
              j=nPotPar1+RankSol(i)
              Sol(j)=Sol(j)+1.0d0
              ! ... accumulate and check if we reached nMat
              nCumMat=nCumMat+1
              if (nCumMat == nMat) then
                ! To cater for real vs. integer issues
                TmpI=sum(nint(Sol(nPotPar1+RankSol(1:nPar2))))
                if (TmpI /= nMat) then
                  if (TmpI > nMat) then
                    Sol(j)=dble(nint(Sol(j))-1)
                  else
                    Sol(j)=dble(nint(Sol(j))+1)
                  end if
                end if
                exit
              end if
            end do
          end do
        end if
      end if

      ! --- Genetic contributions (nVec) ---

      nVecPar1(:)=0

      ! "Parent1"
      ! ... get integer values
      SolInt(1:nPotPar1)=nint(Sol(1:nPotPar1))
      ! ... remove negatives
      do i=1,nPotPar1
        if (SolInt(i) < 0) then
          SolInt(i)=0
        end if
      end do
      ! ... map internal to external order
      nVecPar1(:)=SolInt(1:nPotPar1)
      if (.not.GenderMatters) then
        nVec(:)=nVecPar1(:)
      else
        nVec(IdPotPar1)=nVecPar1(:)
      end if

      ! "Parent2"
      if (GenderMatters) then
        nVecPar2(:)=0
        ! ... get integer values
        SolInt(1:nPotPar2)=nint(Sol((nPotPar1+1):(nPotPar1+nPotPar2)))
        ! ... remove negatives
        do i=1,nPotPar2
          if (SolInt(i) < 0) then
            SolInt(i)=0
          end if
        end do
        ! ... map internal to external order
        nVecPar2(:)=SolInt(1:nPotPar2)
        nVec(IdPotPar2)=nVecPar2(:)
      end if

      xVec(:)=dble(nVec(:))/(dble(2*nMat))

      ! --- PAGE ---

      if (PAGE) then
        GeneEdit(:)=0.0d0
        if (.not.GenderMatters) then
          RankSol(1:nInd)=MrgRnk(Sol((nPotPar1+nMat+1):(nPotPar1+nMat+nInd)))
          GeneEdit(RankSol(nInd:(nInd-PAGEPar1Max+1):-1))=1.0d0 ! MrgRnk ranks small to large
        else
          if (PAGEPar1) then
            RankSol(1:nPotPar1)=MrgRnk(Sol((nPotPar1+nPotPar2+nMat+1):(nPotPar1+nPotPar2+nMat+nPotPar1)))
            GeneEdit(IdPotPar1(RankSol(nPotPar1:(nPotPar1-PAGEPar1Max+1):-1)))=1.0d0 ! MrgRnk ranks small to large
          end if
          if (PAGEPar2) then
            RankSol(1:nPotPar2)=MrgRnk(Sol((nPotPar1+nPotPar2+nMat+nPotPar1+1):(nPotPar1+nPotPar2+nMat+nPotPar1+nPotPar2)))
            GeneEdit(IdPotPar2(RankSol(nPotPar2:(nPotPar2-PAGEPar2Max+1):-1)))=1.0d0 ! MrgRnk ranks small to large
          end if
        end if
      end if

      ! --- Genetic gain ---

      Criterion%Gain=dot_product(xVec,Bv)
      Criterion%GainStand=dot_product(xVec,BvStand)

      if (PAGE) then
        Criterion%Gain=Criterion%Gain+dot_product(xVec,BvPAGE(:)*GeneEdit(:))
        Criterion%GainStand=Criterion%GainStand+dot_product(xVec,BvPAGEStand(:)*GeneEdit(:))
        ! TODO: how do we handle costs?
      end if

      TmpR=Criterion%GainStand!-GainMinStand ! TODO: do we need this difference to the gain under minimum inbreeding?
      if (ToLower(trim(CritType)) == "min") then
        Criterion%Value=Criterion%Value
      else
        Criterion%Value=Criterion%Value+TmpR
      end if

      ! --- Future population inbreeding (=selected group coancestry) ---

      ! xA
      do i=1,nInd
        TmpVec(i,1)=dot_product(xVec,RelMtx(:,i))
      end do
      ! xAx
      Criterion%PopInb=0.5d0*dot_product(TmpVec(:,1),xVec)
      if (Criterion%PopInb < 0.0d0) then
        write(STDERR,"(a)") "ERROR: negative inbreeding examples have not been tested yet! Stopping."
        write(STDERR,"(a)") " "
        stop 1
      endif

      ! Matrix multiplication with symmetric matrix using BLAS routine
      ! (it was ~5x slower than the above with 1.000 individuals, might be
      !  benefical with larger cases so kept in commented.)
      ! http://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3.html#ga253c8edb8b21d1b5b1783725c2a6b692
      ! Ax
      ! call dsymm(side="l",uplo="l",m=nInd,n=1,alpha=1.0d0,A=RelMtx,lda=nInd,b=xVec,ldb=nInd,beta=0,c=TmpVec,ldc=nInd)
      ! call dsymm(     "l",     "l",  nInd,  1,      1.0d0,  RelMtx,    nInd,  xVec,    nInd,     0,  TmpVec,    nInd)
      ! xAx
      ! PopInb=0.5d0*dot_product(xVec,TmpVec(:,1))
      ! print*,xVec,TmpVec,PopInb
      ! stop 1

      ! F_t = DeltaF + (1 - DeltaF) * F_t-1
      ! DeltaF = (F_t - F_t-1) / (1 - F_t-1)
      Criterion%RatePopInb=(Criterion%PopInb-PopInbOld)/(1.0d0-PopInbOld)

      if (ToLower(trim(CritType)) == "min") then
        Criterion%Value=Criterion%Value-Criterion%PopInb
      else
        TmpR=Criterion%RatePopInb/RatePopInbTarget
        if (TmpR > 1.0d0) then
          TmpR=PopInbPenalty*abs(1.0d0-      TmpR)
        else
          if (PopInbPenaltyBellow) then
            TmpR=PopInbPenalty*abs(1.0d0-1.0d0/TmpR)
          else
            TmpR=0.0d0
          end if
        end if
        Criterion%Value=Criterion%Value-TmpR
        Criterion%Penalty=Criterion%Penalty+TmpR
      endif

      ! --- Mate allocation ---

      if (GenderMatters) then
        ! Distribute parent2 (=female) contributions into matings
        k=0
        do i=1,nPotPar2 ! need to loop whole nVecPar2 as some entries are zero
          do j=1,nVecPar2(i)
            k=k+1
            MatPar2(k)=IdPotPar2(i)
          end do
        end do
        ! Reorder parent2 contributions according to the rank of matings
        RankSol(1:nMat)=MrgRnk(Sol((nPotPar1+nPotPar2+1):(nPotPar1+nPotPar2+nMat)))
        MatPar2(:)=MatPar2(RankSol(1:nMat))
      else
        ! Distribute one half of contributions into matings
        k=0
        do while (k < nMat)
          do i=1,nPotPar1 ! need to loop whole nVecPar1 as some entries are zero
            l=nVecPar1(i)/2
            do j=1,l
              if (k == nMat) then
                exit
              end if
              k=k+1
              MatPar2(k)=IdPotPar1(i)
              nVecPar1(i)=nVecPar1(i)-1
            end do
          end do
        end do
        ! Reorder one half of contributions according to the rank of matings
        RankSol(1:nMat)=MrgRnk(Sol((nPotPar1+1):(nPotPar1+nMat)))
        MatPar2(:)=MatPar2(RankSol(1:nMat))
      end if

      ! Pair the contributions
      k=nMat ! MrgRnk ranks small to large
      MateAlloc(:,:)=0
      if (GenderMatters .or. SelfingAllowed) then
        ! When GenderMatters selfing can not happen (we have two distinct sets of parents,
        ! unless the user adds individuals of one sex in both sets) and when SelfingAllowed
        ! we do not need to care about it - faster code
        do i=1,nPotPar1
          do j=1,nVecPar1(i)
            !if (k<2) print*,k,i,j,nVecPar1(i),nMat,sum(nVecPar1(:))
            MateAlloc(1,k)=IdPotPar1(i)
            MateAlloc(2,k)=MatPar2(k)
            k=k-1
          end do
        end do
      else
        ! When .not. GenderMatters selfing can happen (we have one set of parents)
        ! and when .not. SelfingAllowed we do need to care about it - slower code
        do i=1,nPotPar1
          do j=1,nVecPar1(i)
            MateAlloc(1,k)=IdPotPar1(i)
            if (MatPar2(k) == IdPotPar1(i)) then
              ! Try to avoid selfing by swapping the MatPar2 and Rank elements
              do l=k,1,-1
                if (MatPar2(l) /= IdPotPar1(i)) then
                  MatPar2([k,l])=MatPar2([l,k])
                  Sol(nPotPar1+RankSol([k,l]))=Sol(nPotPar1+RankSol([l,k]))
                  exit
                end if
              end do
              if (l < 1) then ! Above loop ran out without finding a swap
                Criterion%Value=Criterion%Value-SelfingPenalty
                Criterion%Penalty=Criterion%Penalty+SelfingPenalty
              end if
            end if
            MateAlloc(2,k)=MatPar2(k)
            k=k-1
          end do
        end do
      end if

      ! --- Progeny inbreeding ---

      TmpR=0.0d0
      do i=1,nMat
        ! Try to speed-up retrieval by targeting lower triangle
        TmpMin=minval([MateAlloc(1,i),MateAlloc(2,i)])
        TmpMax=maxval([MateAlloc(1,i),MateAlloc(2,i)])
        TmpR=TmpR+0.5d0*RelMtx(TmpMax,TmpMin)
      end do
      Criterion%PrgInb=TmpR/dble(nMat)
      TmpR=PrgInbPenalty*Criterion%PrgInb
      Criterion%Value=Criterion%Value-TmpR
      Criterion%Penalty=Criterion%Penalty+TmpR
    end subroutine

    !###########################################################################
end module

!###############################################################################

program AlphaMate

  use ISO_Fortran_Env, STDIN=>input_unit,STDOUT=>output_unit,STDERR=>error_unit
  use IFPort,only : SystemQQ
  use AlphaMateMod
  use AlphaHouseMod, only : Int2Char

  implicit none

  real(real32) :: Start,Finish

  logical :: Success

  character(len=100) :: DumC

  call cpu_time(Start)
  call AlphaMateTitles

  Success=SystemQQ(RMDIR//" AlphaMateResults")
  if (.not.Success) then
    write(STDERR,"(a)") "ERROR: Failed to remove the old output folder (AlphaMateResults)!"
    write(STDERR,"(a)") " "
    stop 1
  end if
  Success=SystemQQ(MKDIR//" AlphaMateResults")
  if (.not.Success) then
    write(STDERR,"(a)") "ERROR: Failed to make the output folder (AlphaMateResults)!"
    write(STDERR,"(a)") " "
    stop 1
  end if

  call ReadSpecAndDataForAlphaMate
  call SetInbreedingParameters
  call AlphaMateSearch
  call cpu_time(Finish)

  DumC=Int2Char(nint(Finish-Start))
  write(STDOUT,"(3a)") "Time duration of AlphaMate: ",trim(adjustl(DumC))," seconds"
  write(STDOUT,"(a)") " "
end program

!###############################################################################
