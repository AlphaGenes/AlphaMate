
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

  use ISO_Fortran_Env, STDIN => input_unit, STDOUT => output_unit, STDERR => error_unit
  use IFPort, only : SystemQQ
  use OrderPackMod, only : MrgRnk
  use AlphaHouseMod, only : CountLines, Int2Char, Real2Char, RandomOrder, SetSeed, ToLower, FindLoc
  use AlphaStatMod, only : CalcDescStat, DescStatD, CalcDescStatMatrix, CalcDescStatSymMatrix, CalcDescStatLowTriMatrix, DescStatMatrixD
  use AlphaEvolveMod, only : AlphaEvolveSol, DifferentialEvolution, RandomSearch

  implicit none

  type, extends(AlphaEvolveSol) :: AlphaMateSol
    real(real64)                :: Penalty
    real(real64)                :: Gain
    real(real64)                :: GainStand
    real(real64)                :: PopInb
    real(real64)                :: RatePopInb
    real(real64)                :: PrgInb
    real(real64), allocatable   :: GenericIndVal(:)
    real(real64), allocatable   :: GenericMatVal(:)
    real(real64)                :: Cost
    integer(int32), allocatable :: nVec(:)
    real(real64), allocatable   :: xVec(:)
    integer(int32), allocatable :: MatingPlan(:,:)
    real(real64), allocatable   :: GenomeEdit(:)
    contains
      procedure         :: Initialise    => InitialiseAlphaMateSol
      procedure         :: Assign        => AssignAlphaMateSol
      procedure         :: UpdateMean    => UpdateMeanAlphaMateSol
      procedure         :: CalcCriterion => FixSolEtcMateAndCalcCrit
      procedure, nopass :: LogHead       => LogHeadAlphaMateSol
      procedure         :: Log           => LogAlphaMateSol
  end type

  integer(int32) :: nInd, nMat, nPotMat, nPar, nPotPar1, nPotPar2, nMal, nFem, nPar1, nPar2, nFrontierSteps
  integer(int32) :: EvolAlgNSol, EvolAlgNGen, EvolAlgNGenBurnIn, EvolAlgNGenStop, EvolAlgNGenPrint, RanAlgStricter
  integer(int32) :: PAGEPar1Max, PAGEPar2Max, nGenericIndVal, nGenericMatVal
  integer(int32), allocatable :: Gender(:), IdPotPar1(:), IdPotPar2(:), IdPotParSeq(:)

  real(real64) :: LimitPar1Min, LimitPar1Max, LimitPar2Min, LimitPar2Max
  real(real64) :: EvolAlgStopTol, EvolAlgCRBurnIn, EvolAlgCRLate, EvolAlgFBase, EvolAlgFHigh1, EvolAlgFHigh2
  real(real64) :: PopInbOld, PopInbTarget, RatePopInbTarget
  real(real64) :: PopInbWeight, PrgInbWeight, SelfingWeight, LimitPar1Weight, LimitPar2Weight
  real(real64), allocatable :: GenericIndValWeight(:), GenericMatValWeight(:)
  real(real64) :: PAGEPar1Cost, PAGEPar2Cost
  real(real64), allocatable :: BreedVal(:), BreedValStand(:), BreedValPAGE(:), BreedValPAGEStand(:)
  real(real64), allocatable :: RatePopInbFrontier(:), GenericIndValTmp(:), GenericMatValTmp(:)
  real(real64), allocatable :: RelMtx(:,:), GenericIndVal(:,:), GenericMatVal(:,:,:)

  logical :: ModeMin, ModeRan, ModeOpt, BreedValAvailable, GenderMatters, EqualizePar1, EqualizePar2
  logical :: SelfingAllowed, PopInbWeightBellow, EvaluateFrontier
  logical :: PAGE, PAGEPar1, PAGEPar2, GenericIndValAvailable, GenericMatValAvailable

  character(len=100), allocatable :: IdC(:)
  CHARACTER(len=100), PARAMETER :: FMTREAL2CHAR = "(f11.5)"

  CHARACTER(len=100), PARAMETER  :: FMTLOGSTDOUTHEADA = "("
  CHARACTER(len=100), PARAMETER  :: FMTLOGSTDOUTHEADB = "a12)"
  CHARACTER(len=100)             :: FMTLOGSTDOUTHEAD
  CHARACTER(len=100), PARAMETER  :: FMTLOGSTDOUTA = "(i12, "
  CHARACTER(len=100), PARAMETER  :: FMTLOGSTDOUTB = "(1x, f11.5))"
  CHARACTER(len=100)             :: FMTLOGSTDOUT
  CHARACTER(len=12), ALLOCATABLE :: COLNAMELOGSTDOUT(:)

  CHARACTER(len=100), PARAMETER  :: FMTLOGUNITHEADA = "("
  CHARACTER(len=100), PARAMETER  :: FMTLOGUNITHEADB = "a22)"
  CHARACTER(len=100)             :: FMTLOGUNITHEAD
  CHARACTER(len=100), PARAMETER  :: FMTLOGUNITA = "(i22, "
  CHARACTER(len=100), PARAMETER  :: FMTLOGUNITB = "(1x, es21.13e3))"
  CHARACTER(len=100)             :: FMTLOGUNIT
  CHARACTER(len=22), ALLOCATABLE :: COLNAMELOGUNIT(:)

  CHARACTER(len=100), PARAMETER  :: FMTINDHEAD = "(6a12)"
  CHARACTER(len=100), PARAMETER  :: FMTINDHEADEDIT = "(8a12)"
  CHARACTER(len=100), PARAMETER  :: FMTIND = "(a12, 1x, i11, 3(1x, f11.5), 1x, i11)"
  CHARACTER(len=100), PARAMETER  :: FMTINDEDIT = "(a12, 1x, i11, 3(1x, f11.5), 2(1x, i11), 1x, f11.5)"
  CHARACTER(len=100), PARAMETER  :: FMTMATHEAD = "(3a12)"
  CHARACTER(len=100), PARAMETER  :: FMTMAT = "(i12, 2(1x, a11))"
  CHARACTER(len=100), PARAMETER  :: FMTFROHEAD = "(a12, 7a22)"
  CHARACTER(len=100), PARAMETER  :: FMTFRO = "(a12, 7(1x, es21.13e3))"

  private
  public :: AlphaMateTitles, ReadSpecAndDataForAlphaMate, SetupColNamesAndFormats, AlphaMateSearch

  contains

    !###########################################################################

    ! This note is just to make clear what the inbreeding calculations are doing.
    !
    ! The total inbreeding that we see today is:
    !
    ! F_total = F_new + (1 - F_new) * F_old,
    !
    ! where F_new is new/recent  inbreeding built up after  some time point
    !       F_old is old/ancient inbreeding built up before some time point
    !
    ! Equivalent formulas with some other notation are:
    !
    ! F_T = F_IS + (1 - F_IS) * F_ST
    !
    ! where F_T is (total) inbreeding relative to old/ancient base time point
    !       F_IS is (individual-to-subtotal=new/recent) inbreeding relative to recent base time point
    !       F-ST is (subtotal-to-total=old/ancient)     inbreeding at recent base time point relative to old/ancient base time point
    !
    ! F_t = DeltaF + (1 - DeltaF) * F_t-1
    !     = DeltaF * (1 - F_t-1)  + F_t-1
    !
    ! where F_t is (total) inbreeding in generation t
    !       DeltaF is new/recent inbreeding built up between the generations t and t-1 (=rate of inbreeding)
    !       F-t-1 is old/ancient inbreeding built up before generation t
    !
    ! Clearly:
    !
    ! F_total   = F_T  = F_t
    ! F_new     = F_IS = DeltaF
    ! F_old     = F_ST = F_t-1

    !###########################################################################

    subroutine AlphaMateTitles
      write(STDOUT, "(a)") ""
      write(STDOUT, "(a30, a, a30)") " ", "**********************", " "
      write(STDOUT, "(a30, a, a30)") " ", "*                    *", " "
      write(STDOUT, "(a30, a, a30)") " ", "*     AlphaMate      *", " "
      write(STDOUT, "(a30, a, a30)") " ", "*                    *", " "
      write(STDOUT, "(a30, a, a30)") " ", "**********************"
      write(STDOUT, "(a30, a, a30)") " ", "Version:"//TOSTRING(VERS), " "
      write(STDOUT, "(a15, a)")      " ", "Software for optimizing contributions to the next generation"
      write(STDOUT, "(a)") " "
      write(STDOUT, "(a35, a)")      " ", "No Liability"
      write(STDOUT, "(a25, a)")      " ", "Bugs to Gregor.Gorjanc@roslin.ed.ac.uk"
      write(STDOUT, "(a)") ""
    end subroutine

    !###########################################################################

    subroutine ReadSpecAndDataForAlphaMate

      implicit none

      integer(int32) :: i, j, k, l, m, DumI, jMal, jFem, nIndTmp, GenderTmp, Seed
      integer(int32) :: SpecUnit, RelMtxUnit, BreedValUnit, GenderUnit, InbreedUnit
      integer(int32) :: GenericIndValUnit, GenericMatValUnit
      integer(int32), allocatable :: Order(:)

      real(real64) :: BreedValTmp, BreedValTmp2

      logical :: Success

      character(len=1000) :: RelMtxFile, BreedValFile, GenderFile, SeedFile
      character(len=1000) :: GenericIndValFile, GenericMatValFile
      character(len=100) :: DumC, IdCTmp, IdCTmp2

      type(DescStatD) :: VecDescStat
      type(DescStatMatrixD) :: MtxDescStat

      write(STDOUT, "(a)") "--- Specifications ---"
      write(STDOUT, "(a)") " "

      Success=SystemQQ(COPY//" AlphaMateSpec.txt AlphaMateResults"//DASH//"AlphaMateSpec.txt")
      if (.not.Success) then
        write(STDERR, "(a)") "ERROR: Failed to copy the AlphaMateSpec.txt file in the output folder!"
        write(STDERR, "(a)") " "
        stop 1
      end if

      open(newunit=SpecUnit, file="AlphaMateSpec.txt", status="old")
      write(STDOUT, "(a)") "SpecFile: AlphaMateSpec.txt"

      ! --- Mode ---

      ModeMin = .false.
      ModeRan = .false.
      ModeOpt = .false.

      read(SpecUnit, *) DumC, DumC
      if (index(ToLower(trim(DumC)), "min") > 0) then
        ModeMin = .true.
        write(STDOUT, "(a)") "Mode: Min"
      end if
      if (index(ToLower(trim(DumC)), "ran") > 0) then
        ModeRan = .true.
        write(STDOUT, "(a)") "Mode: Ran"
      end if
      if (index(ToLower(trim(DumC)), "opt") > 0) then
        ModeOpt = .true.
        write(STDOUT, "(a)") "Mode: Opt"
      end if
      if (.not.ModeMin .and. .not.ModeRan .and. .not.ModeOpt) then
        write(STDERR, "(a)") "ERROR: Mode must be: Min, Ran, Opt, or a combination of the three, e.g., MinOpt, RanOpt, ...!"
        write(STDERR, "(a)") " "
        stop 1
      end if

      ! --- RelationshipMatrixFile ---

      read(SpecUnit, *) DumC, RelMtxFile
      write(STDOUT, "(a)") "RelationshipMatrixFile: "//trim(RelMtxFile)

      ! --- BreedingValueFile ---

      read(SpecUnit, *) DumC, BreedValFile
      if (ToLower(trim(BreedValFile)) /= "none") then
        BreedValAvailable = .true.
        write(STDOUT, "(a)") "BreedingValueFile: "//trim(BreedValFile)
      else
        BreedValAvailable = .false.
      end if

      ! --- GenderFile ---

      read(SpecUnit, *) DumC, GenderFile
      if (ToLower(trim(GenderFile)) /= "none") then
        GenderMatters = .true.
        write(STDOUT, "(a)") "GenderFile: "//trim(GenderFile)
      else
        GenderMatters = .false.
      end if

      ! --- NumberOfIndividuals ---

      read(SpecUnit, *) DumC, nInd
      write(STDOUT, "(a)") "NumberOfIndividuals: "//trim(Int2Char(nInd))

      ! --- NumberOfMatings ---

      read(SpecUnit, *) DumC, nMat
      write(STDOUT, "(a)") "NumberOfMatings: "//trim(Int2Char(nMat))

      ! --- NumberOfParents ---

      read(SpecUnit, *) DumC, nPar
      write(STDOUT, "(a)") "NumberOfParents: "//trim(Int2Char(nPar))
      if (nPar > nInd) then
        write(STDERR, "(a)") "ERROR: The number of parents can not be larger than the number of individuals!"
        write(STDERR, "(a)") "ERROR: Number of     parents: "//trim(Int2Char(nPar))
        write(STDERR, "(a)") "ERROR: Number of individuals: "//trim(Int2Char(nInd))
        write(STDERR, "(a)") " "
        stop 1
      end if
      if (nMat > nPar) then
        write(STDOUT, "(a)") "NOTE: The number of matings is larger than the number of parents! Was this really the intention?"
        write(STDOUT, "(a)") "NOTE: Number of matings: "//trim(Int2Char(nMat))
        write(STDOUT, "(a)") "NOTE: Number of parents: "//trim(Int2Char(nPar))
        write(STDOUT, "(a)") " "
      end if

      ! --- NumberOfMaleParents ---

      read(SpecUnit, *) DumC, nPar1

      ! --- NumberOfFemaleParents ---

      read(SpecUnit, *) DumC, nPar2

      if (.not.GenderMatters) then
        nPar1 = nPar
      else
        write(STDOUT, "(a)") "NumberOfMaleParents:   "//trim(Int2Char(nPar1))
        write(STDOUT, "(a)") "NumberOfFemaleParents: "//trim(Int2Char(nPar2))
      end if
      if (GenderMatters) then
        if ((nPar1 + nPar2) > nInd) then
          write(STDERR, "(a)") "ERROR: The number of parents can not be larger than the number of individuals!"
          write(STDERR, "(a)") "ERROR: Number of        parents: "//trim(Int2Char(nPar1 + nPar2))
          write(STDERR, "(a)") "ERROR: Number of   male parents: "//trim(Int2Char(nPar1))
          write(STDERR, "(a)") "ERROR: Number of female parents: "//trim(Int2Char(nPar2))
          write(STDERR, "(a)") "ERROR: Number of    individuals: "//trim(Int2Char(nInd))
          write(STDERR, "(a)") " "
          stop 1
        end if
        if ((nPar1 + nPar2) /= nPar) then
          write(STDOUT, "(a)") "NOTE: The number of male and female parents does not match with the total number of parents - redefined!"
          write(STDOUT, "(a)") "NOTE: Number of   male parents: "//trim(Int2Char(nPar1))
          write(STDOUT, "(a)") "NOTE: Number of female parents: "//trim(Int2Char(nPar2))
          write(STDOUT, "(a)") "NOTE: Number of        parents: "//trim(Int2Char(nPar))//" (defined)"
          nPar = nPar1 + nPar2
          write(STDOUT, "(a)") "NOTE: Number of        parents: "//trim(Int2Char(nPar))//" (redefined)"
          write(STDOUT, "(a)") " "
        end if
        if ((nMat > nPar1) .and. (nMat > nPar2)) then
          write(STDOUT, "(a)") "NOTE: The number of matings is larger than the number of male and female parents! Was this really the intention?"
          write(STDOUT, "(a)") "NOTE: Number of        matings: "//trim(Int2Char(nMat))
          write(STDOUT, "(a)") "NOTE: Number of   male parents: "//trim(Int2Char(nPar1))
          write(STDOUT, "(a)") "NOTE: Number of female parents: "//trim(Int2Char(nPar2))
          write(STDOUT, "(a)") " "
        end if
      end if

      ! --- EqualizeParentContributions ---

      read(SpecUnit, *) DumC, DumC
      if (.not.GenderMatters) then
        if      (ToLower(trim(DumC)) == "yes") then
          EqualizePar1 = .true.
          write(STDOUT, "(a)") "EqualizeParentContributions: yes"
        else if (ToLower(trim(DumC)) == "no") then
          EqualizePar1 = .false.
          !write(STDOUT, "(a)") "EqualizeParentContributions: no"
        else
          write(STDERR, "(a)") "ERROR: EqualizeParentContributions must be: Yes or no!"
          write(STDERR, "(a)") " "
          stop 1
        end if
      end if

      ! --- EqualizeMaleParentContributions ---

      read(SpecUnit, *) DumC, DumC
      if (GenderMatters) then
        if      (ToLower(trim(DumC)) == "yes") then
          EqualizePar1 = .true.
          write(STDOUT, "(a)") "EqualizeMaleParentContributions: yes"
          if (mod(nMat, nPar1) /= 0) then
            ! TODO: might consider handling this better at some point
            write(STDERR, "(a)") "ERROR: When contributions are equalized the number of matings needs to divide into the number of parents"
            write(STDERR, "(a)") "ERROR: Number of       matings: "//trim(Int2Char(nMat))
            write(STDERR, "(a)") "ERROR: Number of  male parents: "//trim(Int2Char(nPar1))
            write(STDERR, "(a)") "ERROR: Modulo (should be zero): "//trim(Int2Char(mod(nMat, nPar1)))
            write(STDERR, "(a)") " "
            stop 1
          end if
        else if (ToLower(trim(DumC)) == "no") then
          EqualizePar1 = .false.
          !write(STDOUT, "(a)") "EqualizeMaleParentContributions: no"
        else
          write(STDERR, "(a)") "ERROR: EqualizeMaleParentContributions must be: Yes or No!"
          write(STDERR, "(a)") " "
          stop 1
        end if
      end if

      ! --- EqualizeFemaleParentContributions ---

      read(SpecUnit, *) DumC, DumC
      if (GenderMatters) then
        if      (ToLower(trim(DumC)) == "yes") then
          EqualizePar2 = .true.
          write(STDOUT, "(a)") "EqualizeFemaleParentContributions: yes"
          if (mod(nMat, nPar2) /= 0) then
            ! TODO: might consider handling this better at some point
            write(STDERR, "(a)") "ERROR: When contributions are equalized the number of matings needs to divide into the number of parents"
            write(STDERR, "(a)") "ERROR: Number of        matings: "//trim(Int2Char(nMat))
            write(STDERR, "(a)") "ERROR: Number of female parents: "//trim(Int2Char(nPar2))
            write(STDERR, "(a)") "ERROR: Modulo  (should be zero): "//trim(Int2Char(mod(nMat, nPar2)))
            write(STDERR, "(a)") " "
            stop 1
          end if
        else if (ToLower(trim(DumC)) == "no") then
          EqualizePar2 = .false.
          !write(STDOUT, "(a)") "EqualizeFemaleParentContributions: no"
        else
          write(STDERR, "(a)") "ERROR: EqualizeFemaleParentContributions must be: Yes or No!"
          write(STDERR, "(a)") " "
          stop 1
        end if
      end if

      ! --- LimitParentContributions ---

      read(SpecUnit, *) DumC, DumC
      LimitPar1Min = 1.0d0
      LimitPar1Max = huge(LimitPar1Max) - 1.0d0
      LimitPar1Weight = 0.0d0
      if (.not.GenderMatters) then
        if      (ToLower(trim(DumC)) == "yes") then
          if (EqualizePar1) then
            write(STDOUT, "(a)") "LimitParentContributions: no"
            write(STDOUT, "(a)") "NOTE: Limit parent contributions option ignored when contributions are to be equalized!"
            write(STDOUT, "(a)") " "
          else
            backspace(SpecUnit)
            read(SpecUnit, *) DumC, DumC, LimitPar1Min, LimitPar1Max, LimitPar1Weight
            write(STDOUT, "(a)") "LimitParentContributions: yes, min "//trim(Int2Char(nint(LimitPar1Min)))//", max "//&
              trim(Int2Char(nint(LimitPar1Max)))//", penalty weight "//trim(Real2Char(LimitPar1Weight, fmt=FMTREAL2CHAR))
            if (LimitPar1Weight > 0.0) then
              write(STDERR, "(a)") "ERROR: Penalty weight for limiting parent contributions should be zero or negative!"
              write(STDERR, "(a)") " "
              stop 1
            end if
          end if
        else if (ToLower(trim(DumC)) == "no") then
          !write(STDOUT, "(a)") "LimitParentContributions: no"
        else
          write(STDERR, "(a)") "ERROR: LimitParentContributions must be: Yes or No!"
          write(STDERR, "(a)") " "
          stop 1
        end if
      end if

      ! --- LimitMaleParentContributions ---

      read(SpecUnit, *) DumC, DumC
      if (GenderMatters) then
        if      (ToLower(trim(DumC)) == "yes") then
          if (EqualizePar1) then
            write(STDOUT, "(a)") "LimitMaleParentContributions: no"
            write(STDOUT, "(a)") "NOTE: Limit male parent contributions option ignored when contributions are to be equalized!"
            write(STDOUT, "(a)") " "
          else
            backspace(SpecUnit)
            read(SpecUnit, *) DumC, DumC, LimitPar1Min, LimitPar1Max, LimitPar1Weight
            write(STDOUT, "(a)") "LimitMaleParentContributions: yes, min "//trim(Int2Char(nint(LimitPar1Min)))//", max "//&
              trim(Int2Char(nint(LimitPar1Max)))//", penalty weight "//trim(Real2Char(LimitPar1Weight, fmt=FMTREAL2CHAR))
            if (LimitPar1Weight > 0.0) then
              write(STDERR, "(a)") "ERROR: Penalty weight for limiting parent contributions should be zero or negative!"
              write(STDERR, "(a)") " "
              stop 1
            end if
          end if
        else if (ToLower(trim(DumC)) == "no") then
          !write(STDOUT, "(a)") "LimitMaleParentContributions: no"
        else
          write(STDERR, "(a)") "ERROR: LimitMaleParentContributions must be: Yes or No!"
          write(STDERR, "(a)") " "
          stop 1
        end if
      end if

      ! --- LimitFemaleParentContributions ---

      read(SpecUnit, *) DumC, DumC
      LimitPar2Min = 1.0d0
      LimitPar2Max = huge(LimitPar2Max) - 1.0d0
      LimitPar2Weight = 0.0d0
      if (GenderMatters) then
        if      (ToLower(trim(DumC)) == "yes") then
          if (EqualizePar2) then
            write(STDOUT, "(a)") "LimitFemaleParentContributions: no"
            write(STDOUT, "(a)") "NOTE: Limit female parent contributions option ignored when contributions are to be equalized!"
            write(STDOUT, "(a)") " "
          else
            backspace(SpecUnit)
            read(SpecUnit, *) DumC, DumC, LimitPar2Min, LimitPar2Max, LimitPar2Weight
            write(STDOUT, "(a)") "LimitFemaleParentContributions: yes, min "//trim(Int2Char(nint(LimitPar2Min)))//", max "//&
              trim(Int2Char(nint(LimitPar2Max)))//", penalty weight "//trim(Real2Char(LimitPar2Weight, fmt=FMTREAL2CHAR))
            if (LimitPar2Weight > 0.0) then
              write(STDERR, "(a)") "ERROR: Penalty weight for limiting parent contributions should be zero or negative!"
              write(STDERR, "(a)") " "
              stop 1
            end if
          end if
        else if (ToLower(trim(DumC)) == "no") then
          !write(STDOUT, "(a)") "LimitFemaleParentContributions: no"
        else
          write(STDERR, "(a)") "ERROR: LimitFemaleParentContributions must be: Yes or No!"
          write(STDERR, "(a)") " "
          stop 1
        end if
      end if

      ! --- AllowSelfing ---

      read(SpecUnit, *) DumC, DumC
      if      (ToLower(trim(DumC)) == "yes") then
        SelfingAllowed = .true.
        write(STDOUT, "(a)") "AllowSelfing: Yes"
        if (GenderMatters) then
          write(STDOUT, "(a)") "ERROR: When gender matters, AlphaMate can not perform selfing! See the manual for a solution."
          write(STDOUT, "(a)") " "
          stop 1
        end if
      else if (ToLower(trim(DumC)) == "no") then
        SelfingAllowed = .false.
        if (.not.GenderMatters) then
          backspace(SpecUnit)
          read(SpecUnit, *) DumC, DumC, SelfingWeight
          write(STDOUT, "(a)") "AllowSelfing: no, penalty weight "//trim(Real2Char(SelfingWeight, fmt=FMTREAL2CHAR))
          if (SelfingWeight > 0.0) then
            write(STDERR, "(a)") "ERROR: Penalty weight for selfing should be zero or negative!"
            write(STDERR, "(a)") " "
            stop 1
          end if
        end if
      else
        write(STDERR, "(a)") "ERROR: AllowSelfing must be: Yes or No!"
        write(STDERR, "(a)") " "
        stop 1
      end if

      ! --- PAGE ---

      read(SpecUnit, *) DumC, DumC
      if (.not.GenderMatters) then
        if      (ToLower(trim(DumC)) == "yes") then
          PAGEPar1 = .true.
          backspace(SpecUnit)
          read(SpecUnit, *) DumC, DumC, PAGEPar1Max!, PAGEPar1Cost
          write(STDOUT, "(a)") "PAGE: yes, no. of individuals "//trim(Int2Char(PAGEPar1Max))//", cost "//trim(Real2Char(PAGEPar1Cost, fmt=FMTREAL2CHAR))
          if (.not.BreedValAvailable) then
            write(STDERR, "(a)") "ERROR: Can not use PAGE when breeding value file is None!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          if (PAGEPar1Max > nPar) then
            write(STDERR, "(a)") "ERROR: The max number of individuals to edit must not be greater than the total number of parents!"
            write(STDERR, "(a)") "ERROR: Number of             parents: "//trim(Int2Char(nPar))
            write(STDERR, "(a)") "ERROR: Number of individuals to edit: "//trim(Int2Char(PAGEPar1Max))
            write(STDERR, "(a)") " "
          end if
        else if (ToLower(trim(DumC)) == "no") then
          PAGEPar1 = .false.
          !write(STDOUT, "(a)") "PAGE: no"
        else
          write(STDERR, "(a)") "ERROR: PAGE must be: Yes or No!"
          write(STDERR, "(a)") " "
          stop 1
        end if
      end if

      ! --- PAGEMales ---

      read(SpecUnit, *) DumC, DumC
      if (GenderMatters) then
        if      (ToLower(trim(DumC)) == "yes") then
          PAGEPar1 = .true.
          backspace(SpecUnit)
          read(SpecUnit, *) DumC, DumC, PAGEPar1Max!, PAGEPar1Cost
          write(STDOUT, "(a)") "PAGEMales: yes, no. of individuals "//trim(Int2Char(PAGEPar1Max))//", cost "//trim(Real2Char(PAGEPar1Cost, fmt=FMTREAL2CHAR))
          if (.not.BreedValAvailable) then
            write(STDERR, "(a)") "ERROR: Can not use PAGE when breeding value file is None!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          if (PAGEPar1Max > nPar1) then
            write(STDERR, "(a)") "ERROR: The max number of male individuals to edit must not be greater than the total number of male parents!"
            write(STDERR, "(a)") "ERROR: Number of male             parents: "//trim(Int2Char(nPar1))
            write(STDERR, "(a)") "ERROR: Number of male individuals to edit: "//trim(Int2Char(PAGEPar1Max))
            write(STDERR, "(a)") " "
          end if
        else if (ToLower(trim(DumC)) == "no") then
          PAGEPar1 = .false.
          !write(STDOUT, "(a)") "PAGEMales: no"
        else
          write(STDERR, "(a)") "ERROR: PAGEMales must be: Yes or No!"
          write(STDERR, "(a)") " "
          stop 1
        end if
      end if

      ! --- PAGEFemales ---

      read(SpecUnit, *) DumC, DumC
      if (GenderMatters) then
        if      (ToLower(trim(DumC)) == "yes") then
          PAGEPar2 = .false.
          backspace(SpecUnit)
          read(SpecUnit, *) DumC, DumC, PAGEPar2Max!, PAGEPar2Cost
          write(STDOUT, "(a)") "PAGEFemales: yes, no. of individuals "//trim(Int2Char(PAGEPar2Max))//", cost "//trim(Real2Char(PAGEPar2Cost, fmt=FMTREAL2CHAR))
          if (.not.BreedValAvailable) then
            write(STDERR, "(a)") "ERROR: Can not use PAGE when breeding value file is None!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          if (PAGEPar2Max > nPar2) then
            write(STDERR, "(a)") "ERROR: The max number of female individuals to edit must not be greater than the total number of female parents!"
            write(STDERR, "(a)") "ERROR: Number of female             parents: "//trim(Int2Char(nPar2))
            write(STDERR, "(a)") "ERROR: Number of female individuals to edit: "//trim(Int2Char(PAGEPar2Max))
            write(STDERR, "(a)") " "
          end if
        else if (ToLower(trim(DumC)) == "no") then
          PAGEPar2 = .false.
          !write(STDOUT, "(a)") "PAGEFemales: no"
        else
          write(STDERR, "(a)") "ERROR: PAGEFemales must be: Yes or No!"
          write(STDERR, "(a)") " "
          stop 1
        end if
      end if

      if (PAGEPar1 .or. PAGEPar2) then
        PAGE = .true.
      else
        PAGE = .false.
      end if

      ! --- TargetedRateOfPopulationInbreeding ---

      read(SpecUnit, *) DumC, RatePopInbTarget, PopInbWeight, DumC
      write(STDOUT, "(a)") "TargetedRateOfPopulationInbreeding: "//trim(Real2Char(RatePopInbTarget, fmt=FMTREAL2CHAR))//&
        ", penalty weight "//trim(Real2Char(PopInbWeight, fmt=FMTREAL2CHAR))//", mode "//trim(DumC)
      if      (ToLower(trim(DumC)) == "above") then
        PopInbWeightBellow = .false.
      else if (ToLower(trim(DumC)) == "aboveandbellow") then
        PopInbWeightBellow = .true.
      else
        write(STDERR, "(a)") "ERROR: PopInbWeightMode must be: Above or AboveAndBellow!"
        write(STDERR, "(a)") " "
        stop 1
      end if
      if (RatePopInbTarget == 0.0) then
        write(STDERR, "(a)") "ERROR: Can not work with the targeted rate of population inbreeding equal to zero - it is numerically unstable!"
        write(STDERR, "(a)") " "
        stop 1
      end if
      if (PopInbWeight > 0.0) then
        write(STDERR, "(a)") "ERROR: Penalty weight for the targeted rate of population inbreeding should be zero or negative!"
        write(STDERR, "(a)") " "
        stop 1
      end if

      ! --- ProgenyInbreedingWeight (=inbreeding of a mating) ---

      read(SpecUnit, *) DumC, PrgInbWeight
      write(STDOUT, "(a)") "ProgenyInbreedingWeight: "//trim(Real2Char(PrgInbWeight, fmt=FMTREAL2CHAR))
      if (PrgInbWeight > 0.0) then
        write(STDERR, "(a)") "ERROR: Penalty weight for the progeny inbreeding should be zero or negative!"
        write(STDERR, "(a)") " "
        stop 1
      end if

      ! --- EvaluateFrontier ---

      read(SpecUnit, *) DumC, DumC
      if      (ToLower(trim(DumC)) == "no") then
        EvaluateFrontier = .false.
        !write(STDOUT, "(a)") "EvaluateFrontier: no"
      else if (ToLower(trim(DumC)) == "yes") then
        EvaluateFrontier = .true.
        backspace(SpecUnit)
        read(SpecUnit, *) DumC, DumC, nFrontierSteps
        allocate(RatePopInbFrontier(nFrontierSteps))
        backspace(SpecUnit)
        read(SpecUnit, *) DumC, DumC, nFrontierSteps, RatePopInbFrontier(:)
        write(STDOUT, "("//Int2Char(1+nFrontierSteps)//"a)") "EvaluateFrontier: yes, #steps: "//trim(Int2Char(nFrontierSteps))//&
          ", rates of pop. inbreeding: ", (trim(Real2Char(RatePopInbFrontier(i), fmt=FMTREAL2CHAR)), i = 1, nFrontierSteps)
        if (any(RatePopInbFrontier(:) == 0.0)) then
          write(STDERR, "(a)") "ERROR: Can not work with RatePopInbFrontier equal to zero - it is numerically unstable!"
          write(STDERR, "(a)") " "
          stop 1
        end if
      else
        write(STDERR, "(a)") "ERROR: EvaluateFrontier must be: Yes or No!"
        write(STDERR, "(a)") " "
        stop 1
      end if

      ! --- EvolutionaryAlgorithmIterations ---

      read(SpecUnit, *) DumC, EvolAlgNSol, EvolAlgNGen, EvolAlgNGenBurnIn, EvolAlgNGenStop, EvolAlgStopTol, EvolAlgNGenPrint

      ! --- EvolutionaryAlgorithmParameters ---

      read(SpecUnit, *) DumC, EvolAlgCRBurnIn, EvolAlgCRLate, EvolAlgFBase, EvolAlgFHigh1, EvolAlgFHigh2

      ! --- RandomSearchIterationsStricter ---

      read(SpecUnit, *) DumC, RanAlgStricter

      ! --- Seed ---

      read(SpecUnit, *) DumC, DumC
      SeedFile = "AlphaMateResults"//DASH//"Seed.txt"
      if ((ToLower(trim(DumC)) == "unknown") .or. (ToLower(trim(DumC)) == "none")) then
        call SetSeed(SeedFile=SeedFile, Out=Seed)
      else
        backspace(SpecUnit)
        read(SpecUnit, *) DumC, DumI
        call SetSeed(Seed=DumI, SeedFile=SeedFile, Out=Seed)
      end if
      write(STDOUT, "(a)") "Seed: "//trim(Int2Char(Seed))

      ! --- GenericIndividualValuesFile ---
      ! --- GenericIndividualValuesWeight ---

      nGenericIndVal = 0
      read(SpecUnit, *) DumC, GenericIndValFile
      if (ToLower(trim(GenericIndValFile)) == "none") then
        GenericIndValAvailable = .false.
        read(SpecUnit, *) DumC
      else
        GenericIndValAvailable = .true.
        backspace(SpecUnit)
        read(SpecUnit, *) DumC, GenericIndValFile, nGenericIndVal
        allocate(GenericIndValWeight(nGenericIndVal))
        read(SpecUnit, *) DumC, GenericIndValWeight(:)
        write(STDOUT, "("//Int2Char(1 + nGenericIndVal)//"a)") "GenericIndividualValuesFile: "//trim(GenericIndValFile)//&
          ", weight(s): ", (trim(Real2Char(GenericIndValWeight(i), fmt=FMTREAL2CHAR)), i = 1, nGenericIndVal)
      end if

      ! --- GenericMatingValuesFile ---
      ! --- GenericMatingValuesWeight ---

      nGenericMatVal = 0
      read(SpecUnit, *) DumC, GenericMatValFile
      if (ToLower(trim(GenericMatValFile)) == "none") then
        GenericMatValAvailable = .false.
        read(SpecUnit, *) DumC
      else
        GenericMatValAvailable = .true.
        backspace(SpecUnit)
        read(SpecUnit, *) DumC, GenericMatValFile, nGenericMatVal
        allocate(GenericMatValWeight(nGenericMatVal))
        read(SpecUnit, *) DumC, GenericMatValWeight(:)
        write(STDOUT, "("//Int2Char(1 + nGenericMatVal)//"a)") "GenericMatingValuesFile: "//trim(GenericMatValFile)//&
          ", weight(s): ", (trim(Real2Char(GenericMatValWeight(i), fmt=FMTREAL2CHAR)), i = 1, nGenericMatVal)
      end if

      write(STDOUT, "(a)") " "
      close(SpecUnit)

      write(STDOUT, "(a)") "--- Data ---"
      write(STDOUT, "(a)") " "

      ! --- Relationships ---

      write(STDOUT, "(a)") "Relationships"
      allocate(IdC(nInd))
      allocate(RelMtx(nInd, nInd))

      nIndTmp = CountLines(RelMtxFile)
      if (nIndTmp < nInd) then
        write(STDERR, "(a)") "ERROR: The Relationship Matrix file has less rows than there are defined number of individuals!"
        write(STDERR, "(a)") "ERROR: Number of defined individuals:                  "//trim(Int2Char(nInd))
        write(STDERR, "(a)") "ERROR: Number of rows in the relationship matrix file: "//trim(Int2Char(nIndTmp))
        write(STDERR, "(a)") " "
      end if
      open(newunit=RelMtxUnit, file=trim(RelMtxFile), status="old")
      do i = 1, nInd
        read(RelMtxUnit, *) IdC(i), RelMtx(:,i)
      end do
      close(RelMtxUnit)

      MtxDescStat = CalcDescStatSymMatrix(RelMtx)
      write(STDOUT, "(a)") "  - all-relationships"
      write(STDOUT, "(a)") "    - average: "//trim(Real2Char(MtxDescStat%All%Mean, fmt=FMTREAL2CHAR))
      write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(MtxDescStat%All%SD,   fmt=FMTREAL2CHAR))
      write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(MtxDescStat%All%Min,  fmt=FMTREAL2CHAR))
      write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(MtxDescStat%All%Max,  fmt=FMTREAL2CHAR))
      write(STDOUT, "(a)") "  - self-relationships (diagonal)"
      write(STDOUT, "(a)") "    - average: "//trim(Real2Char(MtxDescStat%Diag%Mean, fmt=FMTREAL2CHAR))
      write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(MtxDescStat%Diag%SD,   fmt=FMTREAL2CHAR))
      write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(MtxDescStat%Diag%Min,  fmt=FMTREAL2CHAR))
      write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(MtxDescStat%Diag%Max,  fmt=FMTREAL2CHAR))
      write(STDOUT, "(a)") "  - co-relationships (off-diagonal)"
      write(STDOUT, "(a)") "    - average: "//trim(Real2Char(MtxDescStat%OffDiag%Mean, fmt=FMTREAL2CHAR))
      write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(MtxDescStat%OffDiag%SD,   fmt=FMTREAL2CHAR))
      write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(MtxDescStat%OffDiag%Min,  fmt=FMTREAL2CHAR))
      write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(MtxDescStat%OffDiag%Max,  fmt=FMTREAL2CHAR))
      write(STDOUT, "(a)") " "

      ! --- Group coancestry (=current population inbreeding) ---

      PopInbOld = MtxDescStat%All%Mean / 2.0d0

      ! Targeted future population inbreeding
      ! F_t = DeltaF + (1 - DeltaF) * F_t-1
      PopInbTarget = RatePopInbTarget + (1.0d0 - RatePopInbTarget) * PopInbOld

      ! Report
      write(STDOUT, "(a)") "Coancestry/Inbreeding"
      write(STDOUT, "(a)") "  - current group coancestry/inbreeding:    "//trim(Real2Char(PopInbOld,        fmt=FMTREAL2CHAR))
      write(STDOUT, "(a)") "  - targeted rate of population inbreeding: "//trim(Real2Char(RatePopInbTarget, fmt=FMTREAL2CHAR))
      write(STDOUT, "(a)") "  - targeted future population inbreeding:  "//trim(Real2Char(PopInbTarget,     fmt=FMTREAL2CHAR))
      write(STDOUT, "(a)") " "

      open(newunit=InbreedUnit, file="AlphaMateResults"//DASH//"ConstraintPopulationInbreeding.txt", status="unknown")
      write(InbreedUnit, "(a, f)") "Old_coancestry, ", PopInbOld
      write(InbreedUnit, "(a, f)") "Targeted_rate_of_population_inbreeding, ", RatePopInbTarget
      write(InbreedUnit, "(a, f)") "Targeted_future_population_inbreeding, ", PopInbTarget
      close(InbreedUnit)

      ! --- Shuffle the data ---

      ! To avoid having good animals together - better for Evolutionary algorithm
      ! This is done upfront for Id and RelMtx only. Then all the data is read and mapped to this order
      ! TODO: is it worth it?
      ! TODO: does this scale well with large data sets?

      allocate(Order(nInd))
      Order = RandomOrder(nInd)
      IdC(:) = IdC(Order)
      RelMtx(:,:) = RelMtx(Order, Order)
      deallocate(Order)

      ! --- Breeding values ---

      allocate(BreedVal(nInd))
      allocate(BreedValStand(nInd))
      if (PAGE) then
        allocate(BreedValPAGE(nInd))
        allocate(BreedValPAGEStand(nInd))
      end if

      if (.not.BreedValAvailable) then
        BreedVal(:) = 0.0d0
        BreedValStand(:) = 0.0d0
      else
        write(STDOUT, "(a)") "Breeding values"
        nIndTmp = CountLines(BreedValFile)
        if (nIndTmp /= nInd) then
          write(STDERR, "(a)") "ERROR: Number of individuals in the Breeding Value file and the Relationship Matrix file is not the same!"
          write(STDERR, "(a)") "ERROR: Number of individuals in the Relationship Matrix file: "//trim(Int2Char(nInd))
          write(STDERR, "(a)") "ERROR: Number of individuals in the Breeding Value file:      "//trim(Int2Char(nIndTmp))
          write(STDERR, "(a)") " "
          stop 1
        end if
        open(newunit=BreedValUnit, file=trim(BreedValFile), status="old")
        do i = 1, nInd
          if (PAGE) then
            read(BreedValUnit, *) IdCTmp, BreedValTmp, BreedValTmp2
          else
            read(BreedValUnit, *) IdCTmp, BreedValTmp
          end if
          j = FindLoc(IdCTmp, IdC)
          if (j == 0) then
            write(STDERR, "(a)") "ERROR: Individual "//trim(IdCTmp)//" from the Breeding Value file not present in the Relationship Matrix File!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          BreedVal(j) = BreedValTmp
          if (PAGE) then
            BreedValPAGE(j) = BreedValTmp2
          end if
        end do
        close(BreedValUnit)

        VecDescStat = CalcDescStat(BreedVal)
        BreedValStand(:) = (BreedVal(:) - VecDescStat%Mean) / VecDescStat%SD
        write(STDOUT, "(a)") "  - average: "//trim(Real2Char(VecDescStat%Mean, fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "  - st.dev.: "//trim(Real2Char(VecDescStat%SD,   fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "  - minimum: "//trim(Real2Char(VecDescStat%Min,  fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "  - maximum: "//trim(Real2Char(VecDescStat%Max,  fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") " "

        if (PAGE) then
          ! must have the same scaling as breeding values!!!!
          BreedValPAGEStand(:) = (BreedValPAGE(:) - VecDescStat%Mean) / VecDescStat%SD
          ! only the PAGE bit of BreedVal
          BreedValPAGE(:) = BreedValPAGE(:) - BreedVal(:)
          BreedValPAGEStand(:) = BreedValPAGEStand(:) - BreedValStand(:)
          VecDescStat = CalcDescStat(BreedValPAGE)
          write(STDOUT, "(a)") "Genome editing increments"
          write(STDOUT, "(a)") "  - average: "//trim(Real2Char(VecDescStat%Mean, fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "  - st.dev.: "//trim(Real2Char(VecDescStat%SD,   fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "  - minimum: "//trim(Real2Char(VecDescStat%Min,  fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "  - maximum: "//trim(Real2Char(VecDescStat%Max,  fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") " "
        end if
      end if

      ! --- Gender ---

      allocate(Gender(nInd))

      Gender(:) = 0
      if (GenderMatters) then
        write(STDOUT, "(a)") "Gender"
        nIndTmp = CountLines(GenderFile)
        if (nIndTmp /= nInd) then
          write(STDERR, "(a)") "ERROR: Number of individuals in the Gender file and the Relationship Matrix file is not the same!"
          write(STDERR, "(a)") "ERROR: Number of individuals in the Relationship Matrix file: "//trim(Int2Char(nInd))
          write(STDERR, "(a)") "ERROR: Number of individuals in the Gender file:              "//trim(Int2Char(nIndTmp))
          write(STDERR, "(a)") " "
          stop 1
        end if

        nMal = 0
        nFem = 0

        open(newunit=GenderUnit, file=trim(GenderFile), status="old")
        do i = 1, nInd
          read(GenderUnit, *) IdCTmp, GenderTmp
          if      (GenderTmp == 1) then
            nMal = nMal + 1
          else if (GenderTmp == 2) then
            nFem = nFem + 1
          else
            write(STDERR, "(a)") "ERROR: Gender code must be either 1 for male individuals or 2 for female individuals!"
            write(STDERR, "(a)") "ERROR: "//trim(Int2Char(i))//" "//trim(IdCTmp)//" "//trim(Int2Char(GenderTmp))
            write(STDERR, "(a)") " "
            stop 1
          end if
          j = FindLoc(IdCTmp, IdC)
          if (j == 0) then
            write(STDERR, "(a)") "ERROR: Individual "//trim(IdCTmp)//" from the Gender file not present in the Relationship Matrix File!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          Gender(j) = GenderTmp
        end do
        close(GenderUnit)

        write(STDOUT, "(a)") "  - number of   males in data: "//trim(Int2Char(nMal))
        write(STDOUT, "(a)") "  - number of females in data: "//trim(Int2Char(nFem))
        write(STDOUT, "(a)") " "

        if (nPar1 > nMat) then
          write(STDERR, "(a)") "ERROR: The number of male parents can not be larger than the number of males"
          write(STDERR, "(a)") "ERROR: Number of male parents: "//trim(Int2Char(nPar1))
          write(STDERR, "(a)") "ERROR: Number of        males: "//trim(Int2Char(nMal))
          write(STDERR, "(a)") " "
          stop 1
        end if
        if (nPar2 > nFem) then
          write(STDERR, "(a)") "ERROR: The number of female parents can not be larger than the number of females"
          write(STDERR, "(a)") "ERROR: Number of female parents: "//trim(Int2Char(nPar2))
          write(STDERR, "(a)") "ERROR: Number of        females: "//trim(Int2Char(nFem))
          write(STDERR, "(a)") " "
          stop 1
        end if
      end if

      ! --- Define potential parents ---

      if (.not.GenderMatters) then
        nPotPar1 = nInd
        nPotPar2 = nInd
        allocate(IdPotPar1(nPotPar1))
        do i = 1, nInd
          IdPotPar1(i) = i
        end do
      else
        nPotPar1 = nMal
        nPotPar2 = nFem
        allocate(IdPotPar1(nPotPar1))
        allocate(IdPotPar2(nPotPar2))
        allocate(IdPotParSeq(nInd))
        jMal = 0
        jFem = 0
        do i = 1, nInd
          if (Gender(i) == 1) then
            jMal = jMal + 1
            IdPotPar1(jMal) = i
            IdPotParSeq(i) = jMal
          else
            jFem = jFem + 1
            IdPotPar2(jFem) = i
            IdPotParSeq(i) = jFem
          end if
        end do
        if (nPar1 > nPotPar1) then
          nPar1 = nPotPar1
          write(STDOUT, "(a)") "NOTE: The number of male parents reduced to the number of male individuals!"
          write(STDOUT, "(a)") " "
        end if
        if (nPar2 > nPotPar2) then
          nPar2 = nPotPar2
          write(STDOUT, "(a)") "NOTE: The number of female parents reduced to the number of female individuals!"
          write(STDOUT, "(a)") " "
        end if
      end if

      ! --- Number of all potential matings ---

      if (GenderMatters) then
        nPotMat = nPotPar1 * nPotPar2
      else
        nPotMat = real(nPotPar1 * nPotPar1) / 2.0
        if (SelfingAllowed) then
          nPotMat = nint(nPotMat + real(nPotPar1) / 2.0)
        else
          nPotMat = nint(nPotMat - real(nPotPar1) / 2.0)
        end if
      end if
      if (nMat > nPotMat) then
        write(STDOUT, "(a)") "NOTE: The number of matings is larger than the number of all potential matings!"
        write(STDOUT, "(a)") "NOTE: Number of all potential matings: "//trim(Int2Char(nPotMat))
        if (GenderMatters) then
          write(STDOUT, "(a)") "NOTE: = all males with all females"
          write(STDOUT, "(a)") "NOTE: = no. of males ("//trim(Int2Char(nPotPar1))//") * no. of females ("//trim(Int2Char(nPotPar2))//")"
        else
          if (SelfingAllowed) then
            write(STDOUT, "(a)") "NOTE: = half-diallel including selfing"
            write(STDOUT, "(a)") "NOTE: = no. of individuals * no. of individuals / 2 + individuals / 2"
          else
            write(STDOUT, "(a)") "NOTE: = half-diallel excluding selfing"
            write(STDOUT, "(a)") "NOTE: = no. of individuals * no. of individuals / 2 - individuals / 2"
          end if
          write(STDOUT, "(a)") "NOTE:   (no. of individuals = "//trim(Int2Char(nPotPar1))//")"
        end if
        write(STDOUT, "(a)") "NOTE: Number of              matings: "//trim(Int2Char(nMat))
        write(STDOUT, "(a)") " "
      end if

      ! --- Generic individual values ---

      if (GenericIndValAvailable) then
        write(STDOUT, "(a)") "Generic individual values"
        nIndTmp = CountLines(GenericIndValFile)
        if (nIndTmp /= nInd) then
          write(STDERR, "(a)") "ERROR: Number of individuals in the Generic Individual Values file and the Relationship Matrix file is not the same!"
          write(STDERR, "(a)") "ERROR: Number of individuals in the Relationship Matrix file:       "//trim(Int2Char(nInd))
          write(STDERR, "(a)") "ERROR: Number of individuals in the Generic Individual Values file: "//trim(Int2Char(nIndTmp))
          write(STDERR, "(a)") " "
          stop 1
        end if
        allocate(GenericIndVal(nInd, nGenericIndVal))
        allocate(GenericIndValTmp(nGenericIndVal))
        GenericIndVal(:,:) = 0.0d0
        open(newunit=GenericIndValUnit, file=GenericIndValFile, status="unknown")
        do i = 1, nInd
          read(GenericIndValUnit, *) IdCTmp, GenericIndValTmp(:)
          j = FindLoc(IdCTmp, IdC)
          if (j == 0) then
            write(STDERR, "(a)") "ERROR: Individual "//trim(IdCTmp)//" from the Generic Individual Values file not present in the Relationship Matrix File!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          GenericIndVal(j,:) = GenericIndValTmp(:)
        end do
        close(GenericIndValUnit)

        do j = 1, nGenericIndVal
          write(STDOUT, "(a)") "  - column "//trim(Int2Char(j))
          VecDescStat = CalcDescStat(GenericIndVal(:,j))
          write(STDOUT, "(a)") "    - average: "//trim(Real2Char(VecDescStat%Mean, fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(VecDescStat%SD,   fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(VecDescStat%Min,  fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(VecDescStat%Max,  fmt=FMTREAL2CHAR))
        end do
        write(STDOUT, "(a)") " "
      end if

      ! --- Generic mating values ---

      if (GenericMatValAvailable) then
        write(STDOUT, "(a)") "Generic mating values"
        DumI = CountLines(GenericMatValFile)
        if (DumI /= nPotMat) then
          write(STDERR, "(a)") "ERROR: Number of matings in the Generic Mating Values file and the number of all potential matings is not the same!"
          write(STDERR, "(a)") "ERROR: Number of all potential matings:                         "//trim(Int2Char(nPotMat))
          write(STDERR, "(a)") "ERROR: Number of individuals in the Generic Mating Values file: "//trim(Int2Char(DumI))
          write(STDERR, "(a)") " "
          stop 1
        end if
        allocate(GenericMatVal(nPotPar1, nPotPar2, nGenericMatVal))
        allocate(GenericMatValTmp(nGenericMatVal))
        GenericMatVal(:,:,:) = 0.0d0
        open(newunit=GenericMatValUnit, file=GenericMatValFile, status="unknown")
        do i = 1, nPotMat
          read(GenericMatValUnit, *) IdCTmp, IdCTmp2, GenericMatValTmp(:)
          j = FindLoc(IdCTmp, IdC)
          if (j == 0) then
            write(STDERR, "(a)") "ERROR: Individual "//trim(IdCTmp)//" from the Generic Mating Values file not present in the Relationship Matrix File!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          k = FindLoc(IdCTmp2, IdC)
          if (k == 0) then
            write(STDERR, "(a)") "ERROR: Individual "//trim(IdCTmp2)//" from the Generic Mating Values file not present in the Relationship Matrix File!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          if (GenderMatters) then
            l = FindLoc(j, IdPotPar1)
            if (l == 0) then
              write(STDERR, "(a)") "ERROR: Individual "//trim(IdCTmp)//" from the first column in the Generic Mating Values file should be a male!"
              write(STDERR, "(a)") "ERROR: Generic Mating Values file (line "//trim(Int2Char(i))//"): "//trim(IdCTmp)//" "//trim(IdCTmp2)
              write(STDERR, "(a)") " "
              stop 1
            end if
            m = FindLoc(k, IdPotPar2)
            if (l == 0) then
              write(STDERR, "(a)") "ERROR: Individual "//trim(IdCTmp2)//" from the second column in the Generic Mating Values file should be a female!"
              write(STDERR, "(a)") "ERROR: Generic Mating Values file (line "//trim(Int2Char(i))//"): "//trim(IdCTmp)//" "//trim(IdCTmp2)
              write(STDERR, "(a)") " "
              stop 1
            end if
            ! fill full-matrix
            ! - l and m tell respectively the position within IdPotPar1 and IdPotPar2
            ! - values in IdPotPar1 and IdPotPar2 are "joint" Id of males and females
            ! - values in IdPotParSeq are separate Id of males and females (need these to find matching row and column)
            GenericMatVal(IdPotParSeq(IdPotPar1(l)), IdPotParSeq(IdPotPar2(m)), :) = GenericMatValTmp(:)
          else
            ! fill lower-triangle (half-diallel)
            GenericMatVal(maxval([j,k]), minval([j,k]),:) = GenericMatValTmp(:)
          end if
        end do
        close(GenericMatValUnit)

        do k = 1, nGenericMatVal
          write(STDOUT, "(a)") "  - column "//trim(Int2Char(k))
          if (GenderMatters) then
            MtxDescStat = CalcDescStatMatrix(GenericMatVal(:,:,k))
            write(STDOUT, "(a)") "    - average: "//trim(Real2Char(MtxDescStat%All%Mean, fmt=FMTREAL2CHAR))
            write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(MtxDescStat%All%SD,   fmt=FMTREAL2CHAR))
            write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(MtxDescStat%All%Min,  fmt=FMTREAL2CHAR))
            write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(MtxDescStat%All%Max,  fmt=FMTREAL2CHAR))
          else
            if (SelfingAllowed) then
              MtxDescStat = CalcDescStatLowTriMatrix(GenericMatVal(:,:,k))
              write(STDOUT, "(a)") "    - average: "//trim(Real2Char(MtxDescStat%All%Mean, fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(MtxDescStat%All%SD,   fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(MtxDescStat%All%Min,  fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(MtxDescStat%All%Max,  fmt=FMTREAL2CHAR))
            end if
              MtxDescStat = CalcDescStatLowTriMatrix(GenericMatVal(:,:,k), Diag=.false.)
              write(STDOUT, "(a)") "    - average: "//trim(Real2Char(MtxDescStat%OffDiag%Mean, fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(MtxDescStat%OffDiag%SD,   fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(MtxDescStat%OffDiag%Min,  fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(MtxDescStat%OffDiag%Max,  fmt=FMTREAL2CHAR))
          end if
        end do
        write(STDOUT, "(a)") " "
      end if
    end subroutine

    !###########################################################################

    subroutine SetupColNamesAndFormats()
      implicit none
      integer(int32) :: nCol, nColTmp, i

      ! --- Optimisation log ---

      nCol = 9
      if (GenericIndValAvailable) then
        nCol = nCol + nGenericIndVal
      end if
      if (GenericMatValAvailable) then
        nCol = nCol + nGenericMatVal
      end if
      allocate(COLNAMELOGUNIT(nCol))
      allocate(COLNAMELOGSTDOUT(nCol))
      !                    1234567890123456789012
      COLNAMELOGUNIT(1) = "                  Step"
      COLNAMELOGUNIT(2) = "            AcceptRate"
      COLNAMELOGUNIT(3) = "             Criterion"
      COLNAMELOGUNIT(4) = "             Penalties"
      COLNAMELOGUNIT(5) = "                  Gain"
      COLNAMELOGUNIT(6) = "             GainStand"
      COLNAMELOGUNIT(7) = "            PopInbreed"
      COLNAMELOGUNIT(8) = "            RatePopInb"
      COLNAMELOGUNIT(9) = "            PrgInbreed"
      nColTmp = 9
      if (GenericIndValAvailable) then
        do i = 1, nGenericIndVal
          nColTmp = nColTmp + 1
          COLNAMELOGUNIT(nColTmp) = "GenIndVal"//trim(Int2Char(i))
          COLNAMELOGUNIT(nColTmp) = adjustr(COLNAMELOGUNIT(nColTmp))
        end do
      end if
      if (GenericMatValAvailable) then
        do i = 1, nGenericMatVal
          nColTmp = nColTmp + 1
          COLNAMELOGUNIT(nColTmp) = "GenMatVal"//trim(Int2Char(i))
          COLNAMELOGUNIT(nColTmp) = adjustr(COLNAMELOGUNIT(nColTmp))
        end do
      end if
      do i = 1, nCol
        COLNAMELOGSTDOUT(i) = COLNAMELOGUNIT(i)(13:22)
        COLNAMELOGSTDOUT(i) = adjustr(COLNAMELOGSTDOUT(i))
      end do
      FMTLOGSTDOUTHEAD = trim(FMTLOGSTDOUTHEADA)//trim(Int2Char(nCol))  //trim(FMTLOGSTDOUTHEADB)
      FMTLOGSTDOUT     = trim(FMTLOGSTDOUTA)    //trim(Int2Char(nCol-1))//trim(FMTLOGSTDOUTB)
      FMTLOGUNITHEAD   = trim(FMTLOGUNITHEADA)  //trim(Int2Char(nCol)  )//trim(FMTLOGUNITHEADB)
      FMTLOGUNIT       = trim(FMTLOGUNITA)      //trim(Int2Char(nCol-1))//trim(FMTLOGUNITB)
    end subroutine

    !###########################################################################

    subroutine AlphaMateSearch

      implicit none

      integer(int32) :: nParam, k, FrontierUnit

      real(real64) :: PopInbTargetHold, RatePopInbTargetHold
      real(real64), allocatable :: InitEqual(:,:)

      character(len=1000) :: LogFile, ContribFile, MatingFile
      character(len=100) :: DumC

      type(AlphaMateSol) :: SolMin, SolRan, SolOpt, Sol

      ! --- Number of parameters to optimise ---

      if (GenderMatters) then
        nParam = nPotPar1 + nPotPar2 + nMat
      else
        nParam = nPotPar1 + nMat
      end if

      if (PAGE) then
        nParam = nParam + nInd
      end if

      ! --- Optimise for minimum inbreeding ---

      if (ModeMin) then
        write(STDOUT, "(a)") "--- Optimise for minimum inbreeding --- "
        write(STDOUT, "(a)") " "

        LogFile     = "AlphaMateResults"//DASH//"OptimisationLogMinimumInbreeding.txt"
        ContribFile = "AlphaMateResults"//DASH//"IndividualResultsMinimumInbreeding.txt"
        MatingFile  = "AlphaMateResults"//DASH//"MatingResultsMinimumInbreeding.txt"

        allocate(InitEqual(nParam, nint(real(EvolAlgNSol * 0.1))))
        InitEqual(:,:) = 1.0d0 ! A couple of solutions that would give equal contributions for everybody

        call DifferentialEvolution(nParam=nParam, nSol=EvolAlgNSol, Init=InitEqual, nGen=EvolAlgNGen, nGenBurnIn=EvolAlgNGenBurnIn, &
           nGenStop=EvolAlgNGenStop, StopTolerance=EvolAlgStopTol, nGenPrint=EvolAlgNGenPrint, File=LogFile, CritType="min", &
           CRBurnIn=EvolAlgCRBurnIn, CRLate=EvolAlgCRLate, FBase=EvolAlgFBase, FHigh1=EvolAlgFHigh1, FHigh2=EvolAlgFHigh2, &
           BestSol=SolMin)

        deallocate(InitEqual)

        call SaveSolution(SolMin, ContribFile, MatingFile)
      end if

      ! --- Evaluate 'average' inbreeding under random mating ---

      if (ModeRan) then
        write(STDOUT, "(a)") "--- Evaluate 'average' inbreeding under random mating --- "
        write(STDOUT, "(a)") " "

        LogFile = "AlphaMateResults"//DASH//"OptimisationLogRandomMating.txt"

        allocate(InitEqual(nParam, nint(real(EvolAlgNSol * 0.1))))
        InitEqual(:,:) = 1.0d0 ! A couple of solutions that would give equal contributions for everybody

        call RandomSearch(Mode="avg", nParam=nParam, Init=InitEqual, nSamp=EvolAlgNSol*EvolAlgNGen*RanAlgStricter, nSampStop=EvolAlgNGenStop*RanAlgStricter, &
          StopTolerance=EvolAlgStopTol/real(RanAlgStricter), nSampPrint=EvolAlgNGenPrint, File=LogFile, CritType="ran", BestSol=SolRan)

        deallocate(InitEqual)
      end if

      ! --- Optimise for maximum gain with constraint on inbreeding ---

      if (ModeOpt) then
        write(STDOUT, "(a)") "--- Optimise for maximum gain with constraint on inbreeding ---"
        write(STDOUT, "(a)") " "

        LogFile     = "AlphaMateResults"//DASH//"OptimisationLogOptimumGain.txt"
        ContribFile = "AlphaMateResults"//DASH//"IndividualResultsOptimumGain.txt"
        MatingFile  = "AlphaMateResults"//DASH//"MatingResultsOptimumGain.txt"

        call DifferentialEvolution(nParam=nParam, nSol=EvolAlgNSol, nGen=EvolAlgNGen, nGenBurnIn=EvolAlgNGenBurnIn, &
          nGenStop=EvolAlgNGenStop, StopTolerance=EvolAlgStopTol, nGenPrint=EvolAlgNGenPrint, File=LogFile, CritType="opt", &
          CRBurnIn=EvolAlgCRBurnIn, CRLate=EvolAlgCRLate, FBase=EvolAlgFBase, FHigh1=EvolAlgFHigh1, FHigh2=EvolAlgFHigh2, &
          BestSol=SolOpt)

        call SaveSolution(SolOpt, ContribFile, MatingFile)
      end if

      ! --- Evaluate the full frontier ---

      if (EvaluateFrontier) then
        write(STDOUT, "(a)") "--- Evaluate the full frontier (this might take some time!) ---"
        write(STDOUT, "(a)") " "

        PopInbWeightBellow = .true. ! we want to target certain rates of inbreeding

        open(newunit=FrontierUnit, file="AlphaMateResults"//DASH//"Frontier.txt", status="unknown")
        !                                1234567890123456789012
        write(FrontierUnit, FMTFROHEAD) "        Step", &
                                        "             Criterion", &
                                        "             Penalties", &
                                        "                  Gain", &
                                        "             GainStand", &
                                        "            PopInbreed", &
                                        "            RatePopInb", &
                                        "            PrgInbreed"
        ! TODO: add the generic stuff from the log? Just call This%Log?
        if (ModeMin) then
          DumC = "Min"
          write(FrontierUnit, FMTFRO) adjustl(DumC), SolMin%Criterion, SolMin%Penalty, SolMin%Gain, SolMin%GainStand, SolMin%PopInb, SolMin%RatePopInb, SolMin%PrgInb
        end if
        if (ModeRan) then
          DumC = "Ran"
          write(FrontierUnit, FMTFRO) adjustl(DumC), SolRan%Criterion, SolRan%Penalty, SolRan%Gain, SolRan%GainStand, SolRan%PopInb, SolRan%RatePopInb, SolRan%PrgInb
        end if
        if (ModeOpt) then
          DumC = "Opt"
          write(FrontierUnit, FMTFRO) adjustl(DumC), SolOpt%Criterion, SolOpt%Penalty, SolOpt%Gain, SolOpt%GainStand, SolOpt%PopInb, SolOpt%RatePopInb, SolOpt%PrgInb
        end if

        ! Hold old results
        PopInbTargetHold = PopInbTarget
        RatePopInbTargetHold = RatePopInbTarget

        ! Evaluate
        do k = 1, nFrontierSteps
          RatePopInbTarget = RatePopInbFrontier(k)
          ! F_t = DeltaF + (1 - DeltaF) * F_t-1
          PopInbTarget = RatePopInbTarget + (1.0d0 - RatePopInbTarget) * PopInbOld
          write(STDOUT, "(a)") "Step "//trim(Int2Char(k))//" out of "//trim(Int2Char(nFrontierSteps))//&
                               " for the rate of population inbreeding "//trim(Real2Char(RatePopInbTarget, fmt=FMTREAL2CHAR))//&
                               " (=future pop. inbreed. "//trim(Real2Char(PopInbTarget, fmt=FMTREAL2CHAR))//")"
          write(STDOUT, "(a)") ""
          LogFile = "AlphaMateResults"//DASH//"OptimisationLogFrontier"//trim(Int2Char(k))//".txt"

          call DifferentialEvolution(nParam=nParam, nSol=EvolAlgNSol, nGen=EvolAlgNGen, nGenBurnIn=EvolAlgNGenBurnIn, &
            nGenStop=EvolAlgNGenStop, StopTolerance=EvolAlgStopTol, nGenPrint=EvolAlgNGenPrint, File=LogFile, CritType="opt", &
            CRBurnIn=EvolAlgCRBurnIn, CRLate=EvolAlgCRLate, FBase=EvolAlgFBase, FHigh1=EvolAlgFHigh1, FHigh2=EvolAlgFHigh2, &
            BestSol=Sol)

          DumC = "Frontier"//trim(Int2Char(k))
          write(FrontierUnit, FMTFRO) adjustl(DumC), Sol%Criterion, Sol%Penalty, Sol%Gain, Sol%GainStand, Sol%PopInb, Sol%RatePopInb, Sol%PrgInb

          ContribFile = "AlphaMateResults"//DASH//"IndividualResultsFrontier"//trim(Int2Char(k))//".txt"
          MatingFile  = "AlphaMateResults"//DASH//"MatingResultsFrontier"//trim(Int2Char(k))//".txt"
          call SaveSolution(Sol, ContribFile, MatingFile)

          if ((RatePopInbTarget - Sol%RatePopInb) > 0.01d0) then
            write(STDOUT, "(a)") "NOTE: Could not achieve the rate of population inbreeding "//trim(Real2Char(RatePopInbTarget, fmt=FMTREAL2CHAR))
            write(STDOUT, "(a)") "NOTE: Stopping the frontier evaluation."
            write(STDOUT, "(a)") ""
            exit
          end if
        end do

        ! Put back old results
        PopInbTarget = PopInbTargetHold
        RatePopInbTarget = RatePopInbTargetHold

        close(FrontierUnit)

      end if
    end subroutine

    !###########################################################################

    subroutine SaveSolution(Sol, ContribFile, MatingFile)

      implicit none

      ! Arguments
      type(AlphaMateSol) :: Sol
      character(len=*)   :: ContribFile
      character(len=*)   :: MatingFile

      ! Other
      integer(int32) :: i, j, ContribUnit, MatingUnit, Rank(nInd)

      ! TODO: should we have constant output no matter which options are switched on?
      open(newunit=ContribUnit, file=ContribFile, status="unknown")
      Rank = MrgRnk(Sol%nVec + BreedValStand / 100.0d0)
      if (.not.PAGE) then
        !                               1234567890123456789012
        write(ContribUnit, FMTINDHEAD) "          Id", &
                                       "      Gender", &
                                       "       Merit", &
                                       " AvgCoancest", &
                                       "  Contribute", &
                                       "    nMatings"
        do i = nInd, 1, -1 ! MrgRnk ranks small to large
          j = Rank(i)
          write(ContribUnit, FMTIND) IdC(j), Gender(j), BreedVal(j), &
                                     0.5d0 * sum(RelMtx(:,j)) / dble(nInd), &
                                     Sol%xVec(j), Sol%nVec(j)
        end do
      else
        !                                   1234567890123456789012
        write(ContribUnit, FMTINDHEADEDIT) "          Id", &
                                           "      Gender", &
                                           "       Merit", &
                                           " AvgCoancest", &
                                           "  Contribute", &
                                           "    nMatings", &
                                           "  GenomeEdit", &
                                           " EditedMerit"
        do i = nInd, 1, -1 ! MrgRnk ranks small to large
          j = Rank(i)
          write(ContribUnit, FMTINDEDIT) IdC(j), Gender(j), BreedVal(j), &
                                         0.5d0 * sum(RelMtx(:,j)) / dble(nInd), &
                                         Sol%xVec(j), Sol%nVec(j), &
                                         nint(Sol%GenomeEdit(j)), BreedVal(j) + Sol%GenomeEdit(j) * BreedValPAGE(j)
        end do
      end if
      close(ContribUnit)

      open(newunit=MatingUnit, file=MatingFile, status="unknown")
      !                              1234567890123456789012
      write(MatingUnit, FMTMATHEAD) "      Mating", &
                                    "     Parent1", &
                                    "     Parent2"
      do i = 1, nMat
        write(MatingUnit, FMTMAT) i, IdC(Sol%MatingPlan(1,i)), IdC(Sol%MatingPlan(2,i))
      end do
      close(MatingUnit)
    end subroutine

    !###########################################################################

    subroutine InitialiseAlphaMateSol(This)
      implicit none

      ! Argument
      class(AlphaMateSol), intent(out) :: This

      ! Initialisation
      This%Criterion = 0.0d0
      This%Penalty = 0.0d0
      This%Gain = 0.0d0
      This%GainStand = 0.0d0
      This%PopInb = 0.0d0
      This%RatePopInb = 0.0d0
      This%PrgInb = 0.0d0
      if (GenericIndValAvailable) then
        allocate(This%GenericIndVal(nGenericIndVal))
        This%GenericIndVal(:) = 0.0d0
      else
        allocate(This%GenericIndVal(0))
      end if
      if (GenericMatValAvailable) then
        allocate(This%GenericMatVal(nGenericMatVal))
        This%GenericMatVal(:) = 0.0d0
      else
        allocate(This%GenericMatVal(0))
      end if
      This%Cost = 0.0d0
      allocate(This%nVec(nInd))
      This%nVec(:) = 0
      allocate(This%xVec(nInd))
      This%xVec(:) = 0.0d0
      allocate(This%MatingPlan(2, nMat))
      This%MatingPlan(:,:) = 0
      if (PAGE) then
        allocate(This%GenomeEdit(nInd))
        This%GenomeEdit(:) = 0.0d0
      else
        allocate(This%GenomeEdit(0))
      end if
    end subroutine

    !###########################################################################

    subroutine AssignAlphaMateSol(Out, In)
      implicit none

      ! Arguments
      class(AlphaMateSol), intent(out)   :: Out
      class(AlphaEvolveSol), intent(in)  :: In

      ! Assignments
      ! (Need to go via the select type stuff as all but the first arguments must
      !  be the same as in the base class/type)
      select type (In)
        class is (AlphaMateSol)
          Out%Criterion       = In%Criterion
          Out%Penalty         = In%Penalty
          Out%Gain            = In%Gain
          Out%GainStand       = In%GainStand
          Out%PopInb          = In%PopInb
          Out%RatePopInb      = In%RatePopInb
          Out%PrgInb          = In%PrgInb
          if (allocated(In%GenericIndVal)) then
            allocate(Out%GenericIndVal(size(In%GenericIndVal)))
            Out%GenericIndVal = In%GenericIndVal
          end if
          if (allocated(In%GenericMatVal)) then
            allocate(Out%GenericMatVal(size(In%GenericMatVal)))
            Out%GenericMatVal = In%GenericMatVal
          end if
          Out%Cost            = In%Cost
          if (allocated(In%nVec)) then
            allocate(Out%nVec(size(In%nVec)))
            Out%nVec          = In%nVec
          end if
          if (allocated(In%xVec)) then
            allocate(Out%xVec(size(In%xVec)))
            Out%xVec          = In%xVec
          end if
          if (allocated(In%MatingPlan)) then
            allocate(Out%MatingPlan(size(In%MatingPlan, dim=1), size(In%MatingPlan, dim=2)))
            Out%MatingPlan    = In%MatingPlan
          end if
          if (allocated(In%GenomeEdit)) then
            allocate(Out%GenomeEdit(size(In%GenomeEdit)))
            Out%GenomeEdit    = In%GenomeEdit
          end if
        class default
          write(STDERR, "(a)") "ERROR: Both Out and In must be of the AlphaMateSol class for AlphaMate!"
          write(STDERR, "(a)") " "
          stop 1
      end select
    end subroutine

    !###########################################################################

    subroutine UpdateMeanAlphaMateSol(This, Add, n)
      implicit none

      ! Arguments
      class(AlphaMateSol), intent(inout) :: This
      class(AlphaEvolveSol), intent(in)  :: Add
      integer(int32), intent(in)         :: n

      ! Other
      real(real64) :: nR, kR

      ! Updates
      nR = dble(n)
      kR = (nR - 1.0d0) / nR

      ! (Need to go via the select type stuff as all but the first arguments must
      !  be the same as in the base class/type)
      select type (Add)
        class is (AlphaMateSol)
          This%Criterion       = This%Criterion     * kR + Add%Criterion     / nR
          This%Penalty         = This%Penalty       * kR + Add%Penalty       / nR
          This%Gain            = This%Gain          * kR + Add%Gain          / nR
          This%GainStand       = This%GainStand     * kR + Add%GainStand     / nR
          This%PopInb          = This%PopInb        * kR + Add%PopInb        / nR
          This%RatePopInb      = This%RatePopInb    * kR + Add%RatePopInb    / nR
          This%PrgInb          = This%PrgInb        * kR + Add%PrgInb        / nR
          if (allocated(This%GenericIndVal)) then
            This%GenericIndVal = This%GenericIndVal * kR + Add%GenericIndVal / nR
          end if
          if (allocated(This%GenericMatVal)) then
            This%GenericMatVal = This%GenericMatVal * kR + Add%GenericMatVal / nR
          end if
          This%Cost            = This%Cost          * kR + Add%Cost          / nR
          if (allocated(This%nVec)) then
            This%nVec          = This%nVec          * kR + Add%nVec          / nR
          end if
          if (allocated(This%xVec)) then
            This%xVec          = This%xVec          * kR + Add%xVec          / nR
          end if
          if (allocated(This%MatingPlan)) then
            This%MatingPlan    = This%MatingPlan    * kR + Add%MatingPlan    / nR
          end if
          if (allocated(This%GenomeEdit)) then
            This%GenomeEdit    = This%GenomeEdit    * kR + Add%GenomeEdit    / nR
          end if
        class default
          write(STDERR, "(a)") "ERROR: Both This and Add must be of the AlphaMateSol class for AlphaMate!"
          write(STDERR, "(a)") " "
          stop 1
      end select
    end subroutine

    !###########################################################################

    subroutine FixSolEtcMateAndCalcCrit(This, Chrom, CritType)

      implicit none

      ! Arguments
      class(AlphaMateSol)          :: This      ! Solution
      real(real64), intent(inout)  :: Chrom(:)  ! Internal representation of the solution
      character(len=*), intent(in) :: CritType  ! Type of criterion (Min, Ran, Opt)

      ! Other
      integer(int32) :: i, j, k, l, g, nCumMat, Rank(nInd), ChromInt(nInd), MatPar2(nMat)
      integer(int32) :: nVecPar1(nPotPar1), nVecPar2(nPotPar2), TmpMin, TmpMax, TmpI

      real(real64) :: TmpVec(nInd,1), TmpR, RanNum

      ! Initialize the solution
      call This%Initialise()

      ! The solution (based on the mate selection driver) has:
      ! - nInd individual contributions
      !   - nPotPar1 individual contributions for "parent1" (males   when GenderMatters, all ind when .not. GenderMatters)
      !   - nPotPar2 individual contributions for "parent2" (females when GenderMatters, meaningful only when GenderMatters)
      ! - nMat     rankings of parent1 1:nMat matings to pair with 1:nPotPar2 "parent2" (see bellow)
      ! - nInd edit indicators
      !   - nPotPar1 edit indicators for "parent1" (males   when GenderMatters, all ind when .not. GenderMatters)
      !   - nPotPar2 edit indicators for "parent2" (females when GenderMatters, present only when GenderMatters)

      ! Say we have Chrom=(|0,2,0,1|...|2.5,1.5,1.0|0,1,0,0|...) then we:
      ! - mate male 2 with the first  available female (rank 2.5)
      ! - mate male 2 with the second available female (rank 1.5)
      ! - mate male 4 with the third  available female (rank 1.0)
      ! - edit male 2

      ! TODO: consider spliting the Chrom() vector internally into a type with
      !       separate vectors to simplify the code, e.g.,
      ! Chrom2%ContPar1
      ! Chrom2%ContPar2
      ! Chrom2%MateRank
      ! Chrom2%EditPar1
      ! Chrom2%EditPar2
      !       and then at the end combine it back. Since I modify some elements
      !       it would have to be put back.

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
        g = 1
      else
        g = 2
      end if
      ! ... find ranks to find the top values
      if (.not.(EqualizePar1 .and. (nPar1 == nPotPar1))) then
        Rank(1:nPotPar1) = MrgRnk(Chrom(1:nPotPar1))
        Rank(1:nPotPar1) = Rank(nPotPar1:1:-1) ! MrgRnk ranks small to large
      end if
      ! ... handle cases with equalized contributions
      if (EqualizePar1) then
        if (nPar1 == nPotPar1) then
          ! ... set integers to all the values (no need for sorting here)
          Chrom(1:nPotPar1) = dble(nMat * g) / dble(nPar1)
        else
          ! ... set integers to the top values
          Chrom(Rank(1:nPar1)) = dble(nMat * g) / dble(nPar1)
          ! TODO: anything better to preserve the order of non contributing individuals? See below!
          Chrom(Rank((nPar1+1):nPotPar1)) = 0.0d0
          ! Chrom(Rank((nPar1+1):nPotPar1)) = -1.0d0
        end if
      else
        ! ... handle cases with unequal contributions
        ! ... work for the defined number or parents
        nCumMat = 0
        do i = 1, nPar1
          j = Rank(i)
          ! ... these would have "zero" contributions and if we hit them then they need at least minimum contrib
          if (Chrom(j) < LimitPar1Min) then
            Chrom(j) = LimitPar1Min
          end if
          ! ... but not above max allowed
          if (Chrom(j) > LimitPar1Max) then
            Chrom(j) = LimitPar1Max
          end if
          ! ... accumulate and check if we reached nMat
          nCumMat = nCumMat + nint(Chrom(j)) ! internally real, externally integer
          if (nCumMat >= nMat * g) then
            ! ... there should be exactly nMat contributions
            if (nCumMat > nMat * g) then
              Chrom(j) = Chrom(j) - dble(nCumMat - nMat * g)
              if (nint(Chrom(j)) < LimitPar1Min) then
                TmpR = LimitPar1Weight * (LimitPar1Min - nint(Chrom(j)))
                This%Criterion = This%Criterion + TmpR
                if (LimitPar1Weight < 0.0d0) then
                  This%Penalty = This%Penalty + abs(TmpR)
                end if
              end if
              nCumMat = nMat * g
            end if
            exit
          end if
        end do
        ! ... increment i if we have hit the exit, do loop would have ended with i=nPar1 + 1
        if (i <= nPar1) then
          i = i + 1
        end if
        ! ... the other individuals do not contribute
        if (i <= nPotPar1) then ! "=" to capture incremented i+1 on the do loop exit
          ! ... zero (the same for all ind so no order)
          Chrom(Rank(i:nPotPar1)) = 0.0d0
          ! ... negative (the same for all ind so no order)
          ! Chrom(Rank(i:nPotPar1)) = -1.0d0
          ! ... negative (variable with partially preserving order)
          !     Found faster convergence than with properly decreasing negative values?
          !     I guess it adds some more randomness, i.e., it causes suffling of individuals that just did not make it onto the mating list.
          ! Chrom(Rank(i:nPotPar1)) = sign(Chrom(Rank(i:nPotPar1)), -1.0d0)
          ! ... negative and properly decreasing
          ! TmpR = maxval(Chrom(Rank(i:nPotPar1)))
          ! if (TmpR > 0.0d0) then
          !     Chrom(Rank(i:nPotPar1)) = Chrom(Rank(i:nPotPar1)) - abs(TmpR)
          ! end if
          ! ... negative (random so no order)
          ! do j = i, nPotPar1 ! TODO: really need this loop?
          !   call random_number(RanNum)
          !   Chrom(Rank(j)) = -1.0d0 * RanNum
          ! end do
        end if
        ! ... nMat still not reached?
        do while (nCumMat < nMat * g)
          ! ... add more contributions
          do i = nPar1, 1, -1 ! to bottom ranked selected individuals (to avoid local optima)
            j = Rank(i)
            Chrom(j) = Chrom(j) + 1.0d0
            ! ... accumulate and check if we reached nMat
            nCumMat = nCumMat + 1
            if (nCumMat >= nMat * g) then
              ! To cater for real vs. integer issues
              TmpI = sum(nint(Chrom(Rank(1:nPar1))))
              if (TmpI /= nMat * g) then
                if (TmpI > nMat * g) then
                  Chrom(j) = dble(nint(Chrom(j)) - 1)
                else
                  Chrom(j) = dble(nint(Chrom(j)) + 1)
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
          Rank(1:nPotPar2) = MrgRnk(Chrom((nPotPar1+1):(nPotPar1+nPotPar2)))
          Rank(1:nPotPar2) = Rank(nPotPar2:1:-1) ! MrgRnk ranks small to large
        end if
        ! ... handle cases with equalized contributions
        if (EqualizePar2) then
          if (nPar2 == nPotPar2) then
            ! ... set integers to all the values (no need for sorting here)
            Chrom((nPotPar1+1):(nPotPar1+nPotPar2)) = dble(nMat) / dble(nPar2)
          else
            ! ... set integers to the top values
            Chrom(nPotPar1+Rank(1:nPar2)) = dble(nMat) / dble(nPar2)
            ! TODO: anything better to preserve the order of non contributing individuals? See below!
            Chrom(nPotPar1+Rank((nPar2+1):nPotPar2)) = 0.0d0
            ! Chrom(nPotPar1+Rank((nPar2+1):nPotPar2)) = -1.0d0
          end if
        else
          ! ... handle cases with unequal contributions
          ! ... work for the defined number or parents
          nCumMat = 0
          do i = 1, nPar2
            j = nPotPar1 + Rank(i)
            ! ... these would have "zero" contributions and if we hit them then they need at least minimum contrib
            if (Chrom(j) < LimitPar2Min) then
              Chrom(j) = LimitPar2Min
            end if
            ! ... but not above max allowed
            if (Chrom(j) > LimitPar2Max) then
              Chrom(j) = LimitPar2Max
            end if
            ! ... accumulate and check if we reached nMat
            nCumMat = nCumMat + nint(Chrom(j)) ! internally real, externally integer
            if (nCumMat >= nMat) then
              ! ... there should be exactly nMat contributions
              if (nCumMat > nMat) then
                Chrom(j) = Chrom(j) - dble(nCumMat - nMat)
                if (nint(Chrom(j)) < LimitPar2Min) then
                  TmpR = LimitPar2Weight * (LimitPar2Min - nint(Chrom(j)))
                  This%Criterion = This%Criterion + TmpR
                  if (LimitPar2Weight < 0.0d0) then
                    This%Penalty = This%Penalty + abs(TmpR)
                  end if
                end if
                nCumMat = nMat
              end if
              exit
            end if
          end do
          ! ... increment i if we have hit the exit, do loop would have ended with i=nPar2+1
          if (i <= nPar2) then
            i = i + 1
          end if
          ! ... the other individuals do not contribute
          if (i <= nPotPar2) then ! "="" to capture incremented i+1 on the do loop exit
            ! ... zero (the same for all ind so no order)
            Chrom(nPotPar1+(Rank(i:nPotPar2))) = 0.0d0
            ! ... negative (the same for all ind so no order)
            ! Chrom(nPotPar1+(Rank(i:nPotPar2))) = -1.0d0
            ! ... negative (variable with partially preserving order, i.e., ~large positives become ~large negatives)
            !     Found faster convergence than with properly decreasing negative values?
            !     I guess it adds some more randomness, i.e., it causes suffling of individuals that just did not make it onto the mating list.
            ! Chrom(nPotPar1+(Rank(i:nPotPar2))) = sign(Chrom(nPotPar1+(Rank(i:nPotPar2))), -1.0d0)
            ! ... negative and properly decreasing
            ! TmpR = maxval(Chrom(nPotPar1+(Rank(i:nPotPar2))))
            ! if (TmpR > 0.0d0) then
            !     Chrom(nPotPar1+(Rank(i:nPotPar2))) = Chrom(nPotPar1+(Rank(i:nPotPar2))) - abs(TmpR)
            ! end if
            ! ... negative (random so no order)
            ! do j = i, nPotPar2 ! TODO: really need this loop?
            !   call random_number(RanNum)
            !   Chrom(nPotPar1+Rank(j)) = -1.0d0 * RanNum
            ! end do
          end if
          ! ... nMat still not reached?
          do while (nCumMat < nMat)
            ! ... add more contributions
            do i = nPar2, 1, -1 ! to bottom ranked selected individuals (to avoid local optima)
              j = nPotPar1 + Rank(i)
              Chrom(j) = Chrom(j) + 1.0d0
              ! ... accumulate and check if we reached nMat
              nCumMat = nCumMat + 1
              if (nCumMat == nMat) then
                ! To cater for real vs. integer issues
                TmpI = sum(nint(Chrom(nPotPar1+Rank(1:nPar2))))
                if (TmpI /= nMat) then
                  if (TmpI > nMat) then
                    Chrom(j) = dble(nint(Chrom(j)) - 1)
                  else
                    Chrom(j) = dble(nint(Chrom(j)) + 1)
                  end if
                end if
                exit
              end if
            end do
          end do
        end if
      end if

      ! --- Genetic contributions (nVec & xVec) ---

      nVecPar1(:) = 0

      ! "Parent1"
      ! ... get integer values
      ChromInt(1:nPotPar1) = nint(Chrom(1:nPotPar1))
      ! ... remove negatives
      do i = 1, nPotPar1
        if (ChromInt(i) < 0) then
          ChromInt(i) = 0
        end if
      end do
      ! ... map internal to external order
      nVecPar1(:) = ChromInt(1:nPotPar1)
      if (.not.GenderMatters) then
        This%nVec(:) = nVecPar1(:)
      else
        This%nVec(IdPotPar1) = nVecPar1(:)
      end if

      ! "Parent2"
      if (GenderMatters) then
        nVecPar2(:) = 0
        ! ... get integer values
        ChromInt(1:nPotPar2) = nint(Chrom((nPotPar1+1):(nPotPar1+nPotPar2)))
        ! ... remove negatives
        do i = 1, nPotPar2
          if (ChromInt(i) < 0) then
            ChromInt(i) = 0
          end if
        end do
        ! ... map internal to external order
        nVecPar2(:) = ChromInt(1:nPotPar2)
        This%nVec(IdPotPar2) = nVecPar2(:)
      end if

      This%xVec(:) = dble(This%nVec(:)) / (dble(2 * nMat))

      ! --- PAGE ---

      if (PAGE) then
        if (.not.GenderMatters) then
          Rank(1:nInd) = MrgRnk(Chrom((nPotPar1+nMat+1):(nPotPar1+nMat+nInd)))
          This%GenomeEdit(Rank(nInd:(nInd-PAGEPar1Max+1):-1)) = 1.0d0 ! MrgRnk ranks small to large
        else
          if (PAGEPar1) then
            Rank(1:nPotPar1) = MrgRnk(Chrom((nPotPar1+nPotPar2+nMat+1):(nPotPar1+nPotPar2+nMat+nPotPar1)))
            This%GenomeEdit(IdPotPar1(Rank(nPotPar1:(nPotPar1-PAGEPar1Max+1):-1))) = 1.0d0 ! MrgRnk ranks small to large
          end if
          if (PAGEPar2) then
            Rank(1:nPotPar2) = MrgRnk(Chrom((nPotPar1+nPotPar2+nMat+nPotPar1+1):(nPotPar1+nPotPar2+nMat+nPotPar1+nPotPar2)))
            This%GenomeEdit(IdPotPar2(Rank(nPotPar2:(nPotPar2-PAGEPar2Max+1):-1))) = 1.0d0 ! MrgRnk ranks small to large
          end if
        end if
      end if

      ! --- Mate allocation ---

      MatPar2(:) = 0
      if (GenderMatters) then
        ! Distribute parent2 (=female) contributions into matings
        k = 0
        do i = 1, nPotPar2 ! need to loop whole nVecPar2 as some entries are zero
          do j = 1, nVecPar2(i)
            k = k + 1
            MatPar2(k) = IdPotPar2(i)
          end do
        end do
        ! Reorder parent2 contributions according to the rank of matings
        Rank(1:nMat) = MrgRnk(Chrom((nPotPar1+nPotPar2+1):(nPotPar1+nPotPar2+nMat)))
        MatPar2(:) = MatPar2(Rank(1:nMat))
      else
        ! Distribute one half of contributions into matings
        k = 0
        do while (k < nMat)
          do i = 1, nPotPar1 ! need to loop whole nVecPar1 as some entries are zero
            l = nVecPar1(i) / 2
            if (mod(nVecPar1(i), 2) == 1) then
              call random_number(RanNum)
              if (RanNum > 0.5) then
                l = l + 1
              end if
            end if
            do j = 1, l
              if (k == nMat) then
                exit
              end if
              k = k + 1
              MatPar2(k) = IdPotPar1(i)
              nVecPar1(i) = nVecPar1(i) - 1
            end do
          end do
        end do
        ! Reorder one half of contributions according to the rank of matings
        Rank(1:nMat) = MrgRnk(Chrom((nPotPar1+1):(nPotPar1+nMat)))
        MatPar2(:) = MatPar2(Rank(1:nMat))
      end if

      ! Pair the contributions (=Mating plan)
      k = nMat ! MrgRnk ranks small to large
      if (GenderMatters .or. SelfingAllowed) then
        ! When gender matters selfing can not happen (we have two distinct sets of parents,
        ! unless the user adds individuals of one sex in both sets) and when SelfingAllowed
        ! we do not need to care about it - faster code
        do i = 1, nPotPar1
          do j = 1, nVecPar1(i)
            !if (k<2) print*, k, i, j, nVecPar1(i), nMat, sum(nVecPar1(:))
            This%MatingPlan(1,k) = IdPotPar1(i)
            This%MatingPlan(2,k) = MatPar2(k)
            k = k - 1
          end do
        end do
      else
        ! When gender does not matter, selfing can happen (we have one set of parents)
        ! and when selfing is not allowed we need to avoid it - slower code
        do i = 1, nPotPar1
          do j = 1, nVecPar1(i)
            This%MatingPlan(1,k) = IdPotPar1(i)
            if (MatPar2(k) == IdPotPar1(i)) then
              ! Try to avoid selfing by swapping the MatPar2 and Rank elements
              do l = k, 1, -1
                if (MatPar2(l) /= IdPotPar1(i)) then
                  MatPar2([k,l]) = MatPar2([l,k])
                  Chrom(nPotPar1+Rank([k,l])) = Chrom(nPotPar1+Rank([l,k]))
                  exit
                end if
              end do
              if (l < 1) then ! Above loop ran out without finding a swap
                This%Criterion = This%Criterion + SelfingWeight
                if (SelfingWeight < 0.0d0) then
                  This%Penalty = This%Penalty + abs(SelfingWeight)
                end if
              end if
            end if
            This%MatingPlan(2,k) = MatPar2(k)
            k = k - 1
          end do
        end do
      end if

      ! --- Genetic gain ---

      if (BreedValAvailable) then
        This%Gain      = dot_product(This%xVec, BreedVal)
        This%GainStand = dot_product(This%xVec, BreedValStand)
        if (PAGE) then
          This%Gain      = This%Gain      + dot_product(This%xVec, BreedValPAGE(:)      * This%GenomeEdit(:))
          This%GainStand = This%GainStand + dot_product(This%xVec, BreedValPAGEStand(:) * This%GenomeEdit(:))
        end if
        if (CritType == "opt") then
          This%Criterion = This%Criterion + This%GainStand
        end if
      end if

      ! --- Generic individual values ---

      if (GenericIndValAvailable) then
        do j = 1, nGenericIndVal
          TmpR = dot_product(This%xVec, GenericIndVal(:,j))
          This%GenericIndVal(j) = TmpR
          TmpR = GenericIndValWeight(j) * This%GenericIndVal(j)
          This%Criterion = This%Criterion + TmpR
          if (GenericIndValWeight(j) < 0.0) then
            This%Penalty = This%Penalty + abs(TmpR)
          end if
        end do
      end if

      ! --- Selected group coancestry (=future population inbreeding) ---

      ! xA
      do i = 1, nInd
        TmpVec(i,1) = dot_product(This%xVec, RelMtx(:,i))
      end do
      ! TODO: consider using matmul instead of repeated dot_product?
      ! TODO: consider using BLAS/LAPACK - perhaps non-symmetric is more optimised?
      ! Matrix multiplication with a symmetric matrix using BLAS routine
      ! (it was ~5x slower than the above with 1.000 individuals, might be benefical with larger cases)
      ! http://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3.html#ga253c8edb8b21d1b5b1783725c2a6b692
      ! call dsymm(side="l", uplo="l", m=nInd, n=1, alpha=1.0d0, A=RelMtx, lda=nInd, b=This%xVec, ldb=nInd, beta=0, c=TmpVec, ldc=nInd)
      ! call dsymm(     "l",      "l",   nInd,   1,       1.0d0,   RelMtx,     nInd,   This%xVec,     nInd,      0,   TmpVec,     nInd)

      ! xAx
      This%PopInb = 0.5d0 * dot_product(TmpVec(:,1), This%xVec)

      if (This%PopInb < 0.0d0) then
        write(STDERR, "(a)") "ERROR: Negative inbreeding cases have not been well tested! Stopping."
        write(STDERR, "(a)") "ERROR: Contact the authors."
        write(STDERR, "(a)") " "
        stop 1
      end if

      ! F_t = DeltaF + (1 - DeltaF) * F_t-1
      ! DeltaF = (F_t - F_t-1) / (1 - F_t-1)
      This%RatePopInb = (This%PopInb - PopInbOld) / (1.0d0 - PopInbOld)

      if      (CritType == "min" .or. CritType == "ran") then
        This%Criterion = This%Criterion - This%PopInb
      else if (CritType == "opt") then
        ! We know the targeted inbreeding so we can work with relative values,
        ! which makes the PopInbWeight generic for ~any scenario.
        TmpR = This%RatePopInb / RatePopInbTarget
        if (TmpR > 1.0d0) then
          ! Rate of inbreeding for this solution is higher than the target
          TmpR = PopInbWeight * abs(1.0d0 - TmpR)
        else
          ! Rate of inbreeding for this solution is lower than the target
          if (PopInbWeightBellow) then
            TmpR = PopInbWeight * abs(1.0d0 - TmpR)
          else
            TmpR = 0.0d0
          end if
        end if
        This%Criterion = This%Criterion + TmpR
        if (PopInbWeight < 0.0d0) then
          This%Penalty = This%Penalty + abs(TmpR)
        end if
      else
        write(STDERR, "(a)") "ERROR: Wrong CritType!!!"
        write(STDERR, "(a)") " "
        stop 1
      end if

      ! --- Progeny inbreeding (=inbreeding of a mating) ---

      TmpR = 0.0d0
      do j = 1, nMat
        ! Lower triangle
        TmpMax = maxval(This%MatingPlan(:,j))
        TmpMin = minval(This%MatingPlan(:,j))
        TmpR = TmpR + 0.5d0 * RelMtx(TmpMax, TmpMin)
      end do
      This%PrgInb = TmpR / dble(nMat)
      TmpR = PrgInbWeight * This%PrgInb
      This%Criterion = This%Criterion + TmpR
      if (PrgInbWeight < 0.0d0) then
        This%Penalty = This%Penalty + abs(TmpR)
      end if

      ! --- Generic mating values ---

      if (GenericMatValAvailable) then
        do k = 1, nGenericMatVal
          TmpR = 0.0d0
          if (GenderMatters) then
            do j = 1, nMat
              TmpR = TmpR + GenericMatVal(IdPotParSeq(This%MatingPlan(1,j)), &
                                          IdPotParSeq(This%MatingPlan(2,j)), k)
            end do
          else
            do j = 1, nMat
              TmpMax = maxval(This%MatingPlan(:,j))
              TmpMin = minval(This%MatingPlan(:,j))
              TmpR = TmpR + GenericMatVal(TmpMax, TmpMin, k)
            end do
          end if
          This%GenericMatVal(k) = TmpR / dble(nMat)
          TmpR = GenericMatValWeight(k) * This%GenericMatVal(k)
          This%Criterion = This%Criterion + TmpR
          if (GenericMatValWeight(k) < 0.0) then
            This%Penalty = This%Penalty + abs(TmpR)
          end if
        end do
      end if

      ! TODO: how should we handle costs?
    end subroutine

    !###########################################################################

    subroutine LogHeadAlphaMateSol(LogUnit)
      implicit none
      integer(int32), intent(in), optional :: LogUnit
      write(STDOUT, FMTLOGSTDOUTHEAD)   COLNAMELOGSTDOUT(:)
      if (present(LogUnit)) then
        write(LogUnit,  FMTLOGUNITHEAD) COLNAMELOGUNIT(:)
      end if
    end subroutine

    !###########################################################################

    subroutine LogAlphaMateSol(This, LogUnit, Gen, AcceptRate)
      implicit none
      class(AlphaMateSol), intent(in)      :: This
      integer(int32), intent(in), optional :: LogUnit
      integer(int32), intent(in)           :: Gen
      real(real64), intent(in)             :: AcceptRate
      if (GenericIndValAvailable) then
        if (GenericMatValAvailable) then
          write(STDOUT,  FMTLOGSTDOUT)  Gen, AcceptRate, This%Criterion, This%Penalty, This%Gain, This%GainStand, This%PopInb, This%RatePopInb, This%PrgInb, This%GenericIndVal, This%GenericMatVal
        else
          write(STDOUT,  FMTLOGSTDOUT)  Gen, AcceptRate, This%Criterion, This%Penalty, This%Gain, This%GainStand, This%PopInb, This%RatePopInb, This%PrgInb, This%GenericIndVal
        end if
      else
        if (GenericMatValAvailable) then
          write(STDOUT,  FMTLOGSTDOUT)  Gen, AcceptRate, This%Criterion, This%Penalty, This%Gain, This%GainStand, This%PopInb, This%RatePopInb, This%PrgInb,                     This%GenericMatVal
        else
          write(STDOUT,  FMTLOGSTDOUT)  Gen, AcceptRate, This%Criterion, This%Penalty, This%Gain, This%GainStand, This%PopInb, This%RatePopInb, This%PrgInb
        end if
      end if
      if (present(LogUnit)) then
        if (GenericIndValAvailable) then
          if (GenericMatValAvailable) then
            write(LogUnit,  FMTLOGUNIT) Gen, AcceptRate, This%Criterion, This%Penalty, This%Gain, This%GainStand, This%PopInb, This%RatePopInb, This%PrgInb, This%GenericIndVal, This%GenericMatVal
          else
            write(LogUnit,  FMTLOGUNIT) Gen, AcceptRate, This%Criterion, This%Penalty, This%Gain, This%GainStand, This%PopInb, This%RatePopInb, This%PrgInb, This%GenericIndVal
          end if
        else
          if (GenericMatValAvailable) then
            write(LogUnit,  FMTLOGUNIT) Gen, AcceptRate, This%Criterion, This%Penalty, This%Gain, This%GainStand, This%PopInb, This%RatePopInb, This%PrgInb,                     This%GenericMatVal
          else
            write(LogUnit,  FMTLOGUNIT) Gen, AcceptRate, This%Criterion, This%Penalty, This%Gain, This%GainStand, This%PopInb, This%RatePopInb, This%PrgInb
          end if
        end if
      end if
    end subroutine

    !###########################################################################
end module

!###############################################################################

program AlphaMate

  use ISO_Fortran_Env, STDIN => input_unit, STDOUT => output_unit, STDERR => error_unit
  use IFPort, only : SystemQQ
  use AlphaMateMod
  use AlphaHouseMod, only : Int2Char

  implicit none

  real(real32) :: Start,Finish

  logical :: Success

  call cpu_time(Start)
  call AlphaMateTitles

  Success=SystemQQ(RMDIR//" AlphaMateResults")
  if (.not.Success) then
    write(STDERR, "(a)") "ERROR: Failed to remove the old output folder (AlphaMateResults)!"
    write(STDERR, "(a)") " "
    stop 1
  end if
  Success=SystemQQ(MKDIR//" AlphaMateResults")
  if (.not.Success) then
    write(STDERR, "(a)") "ERROR: Failed to make the output folder (AlphaMateResults)!"
    write(STDERR, "(a)") " "
    stop 1
  end if

  call ReadSpecAndDataForAlphaMate
  call SetupColNamesAndFormats
  call AlphaMateSearch
  call cpu_time(Finish)

  write(STDOUT, "(a)") "Time duration of AlphaMate: "//trim(Int2Char(nint(Finish - Start)))//" seconds"
  write(STDOUT, "(a)") " "
end program

!###############################################################################
