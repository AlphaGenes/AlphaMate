#ifdef SINGLEPREC
#define FLOATTYPE real32
#define FLOATFUN real
#define DESCSTATFLOAT DescStatReal32
#define DESCSTATMATRIXFLOAT DescStatMatrixReal32
#else
#define FLOATTYPE real64
#define FLOATFUN dble
#define DESCSTATFLOAT DescStatReal64
#define DESCSTATMATRIXFLOAT DescStatMatrixReal64
#endif

! ADD in alphahouse preprocessor to get ah types
#ifdef SINGLEPRECAH
#define FLOATTYPEAH real32
#define FLOATFUNAH real
#define SAMPLEINTELUNIFORMAH SampleIntelUniformS
#else
#define FLOATTYPEAH real64
#define FLOATFUNAH dble
#define SAMPLEINTELUNIFORMAH SampleIntelUniformD
#endif

#ifdef _WIN32

#define STRINGIFY(x)#x
#define TOSTRING(x) STRINGIFY(x)

#define DASH "\"
#define COPY "copy"
#define MD "md"
#define RMDIR "RMDIR /S /Q"
#define RM "del"
#define RENAME "MOVE /Y"
#define SH "BAT"
#define EXE ".exe"
#define NULL " >NUL"

#else

#define STRINGIFY(x)#x
#define TOSTRING(x) STRINGIFY(x)

#define DASH "/"
#define COPY "cp"
#define MD "mkdir"
#define RMDIR "rm -r"
#define RM "rm"
#define RENAME "mv"
#define SH "sh"
#define EXE ""
#define NULL ""

#endif

!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     AlphaMateModule.f90
!
! DESCRIPTION:
!> @brief    Mate selection / Optimum contribution selection
!
!> @details  Optimise selection, maintenance of diversity, and mate allocation
!!           in breeding programs
!
!> @author   Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
!
!> @date      2018-11-02
!
!> @version  0.2.0 (alpha)
!
!-------------------------------------------------------------------------------
module AlphaMateModule
  use Iso_Fortran_Env, STDOUT => output_unit, STDERR => error_unit
  use, intrinsic :: Ieee_Arithmetic
  use ConstantModule, only : FILELENGTH, SPECOPTIONLENGTH, IDLENGTH, RAD2DEG, DEG2RAD
  use OrderPackModule, only : MrgRnk, RapKnr, MulCnt
  use AlphaHouseMod, only : Append, CountLines, GeneratePairing, &
                            Char2Int,Char2Real, Char2Double, Int2Char, Real2Char, &
                            SetSeed, ToLower, ParseToFirstWhitespace, SplitLineIntoTwoParts
  use HashModule
  use IntelRngMod, only : IntitialiseIntelRng, SampleIntelUniformD, RandomOrderIntel
  use Mkl_Vsl_Type
  use Mkl_Vsl
  use AlphaStatMod, only : Mean, StdDev, DescStat, DescStatReal32, DescStatReal64, &
                           DescStatMatrix, DescStatLowTriMatrix,DescStatMatrixReal32, DescStatMatrixReal64
  use AlphaEvolveModule, only : AlphaEvolveSol, AlphaEvolveSpec, AlphaEvolveData, &
                                DifferentialEvolution, RandomSearch
  use AlphaRelateModule
  use Blas95, only : dot , symv
  use Omp_Lib

  implicit none

  private

  ! --- Export ---

  ! Types
  public :: AlphaMateSpec, AlphaMateData, AlphaMateSol
  ! Functions
  public :: AlphaMateTitle, AlphaMateSystem, AlphaMateSearch
  public :: MinCoancestryPct2Degree, Degree2MinCoancestryPct, MaxCriterionPct2Degree, Degree2MaxCriterionPct
  public :: Degree2SelCriterionStd, SelCriterionStd2Degree
  public :: Degree2CoancestryRate, CoancestryRate2Degree
  public :: MinCoancestryPct2CoancestryRate, CoancestryRate2MinCoancestryPct
  public :: MaxCriterionPct2SelCriterionStd, SelCriterionStd2MaxCriterionPct
  public :: Coancestry2CoancestryRate, CoancestryRate2Coancestry
  public :: SelCriterion2SelCriterionStd, SelCriterionStd2SelCriterion

  ! --- Module global objects ---

  type(vsl_stream_state) :: AlphaMateStream

  ! --- Module parameters ---

  REAL(FLOATTYPE), PARAMETER :: TARGETDEGREEFRONTIER(8) = [80, 70, 60, 50, 40, 30, 20, 10]

  INTEGER,                   PARAMETER :: CHARLENGTH = 100
  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTREAL2CHAR = "(f11.5)"
  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTINT2CHAR  = "(i11)"

  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTLOGSTDOUTHEADA = "("
  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTLOGSTDOUTHEADB = "a11)"

  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTLOGSTDOUTA = "(i11, f11.1, 2(f11.5), f11.1, 3(2(f11.5), f11.1), "
  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTLOGSTDOUTB = "(f11.5))"

  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTLOGUNITHEADA = "("
  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTLOGUNITHEADB = "a17)"

  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTLOGUNITA = "(i17, "
  ! CHARACTER(len=CHARLENGTH), PARAMETER :: FMTLOGUNITB = "(1x, es21.13e3))"
  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTLOGUNITB = "(f17.5))"

  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTLOGPOPUNITA = "(2i17, "

  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTCONTRIBUTIONHEADA = "(a"
  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTCONTRIBUTIONHEADB = ", 5a15)"

  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTCONTRIBUTIONHEADEDITA = "(a"
  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTCONTRIBUTIONHEADEDITB = ", 7a15)"

  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTCONTRIBUTIONA = "(a"
  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTCONTRIBUTIONB = ", 4x, i11, 3(4x, f11.5), 4x, i11)"

  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTCONTRIBUTIONEDITA = "(a"
  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTCONTRIBUTIONEDITB = ", 4x, i11, 3(4x, f11.5), 2(4x, i11), 4x, f11.5)"

  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTMATINGHEADA = "(a15, 2a"
  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTMATINGHEADB = ")"

  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTMATINGA = "(i15, 2(1x, a"
  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTMATINGB = "))"

  ! --- Module types ---

  !> @brief Optimisation mode specifications
  type :: AlphaMateModeSpec
    character(len=SPECOPTIONLENGTH) :: Name
    logical                         :: ObjectiveCriterion
    logical                         :: ObjectiveCoancestry
    logical                         :: ObjectiveInbreeding
    real(FLOATTYPE)                 :: TargetDegree
    real(FLOATTYPE)                 :: TargetSelCriterion
    real(FLOATTYPE)                 :: TargetSelCriterionStd
    real(FLOATTYPE)                 :: TargetMaxCriterionPct
    real(FLOATTYPE)                 :: TargetCoancestry
    real(FLOATTYPE)                 :: TargetCoancestryRate
    logical                         :: CoancestryWeightBelow
    real(FLOATTYPE)                 :: TargetMinCoancestryPct
    real(FLOATTYPE)                 :: TargetInbreeding
    real(FLOATTYPE)                 :: TargetInbreedingRate
    logical                         :: InbreedingWeightBelow
    real(FLOATTYPE)                 :: TargetMinInbreedingPct
    real(FLOATTYPE)                 :: Degree
    real(FLOATTYPE)                 :: SelCriterion
    real(FLOATTYPE)                 :: SelCriterionStd
    real(FLOATTYPE)                 :: MaxCriterionPct
    real(FLOATTYPE)                 :: Coancestry
    real(FLOATTYPE)                 :: CoancestryRate
    real(FLOATTYPE)                 :: MinCoancestryPct
    real(FLOATTYPE)                 :: Inbreeding
    real(FLOATTYPE)                 :: InbreedingRate
    real(FLOATTYPE)                 :: MinInbreedingPct
    contains
      procedure :: Initialise       => InitialiseAlphaMateModeSpec
      procedure :: Assign           => AssignAlphaMateModeSpec
      procedure :: SetTargets       => SetTargetsAlphaMateModeSpec
      procedure :: SaveSol2ModeSpec => SaveSol2ModeSpecAlphaMateModeSpec
      procedure :: Write            => WriteAlphaMateModeSpec
      procedure :: LogTargets       => LogTargetsAlphaMateModeSpec
  end type

  !> @brief AlphaMate specifications
  type, extends(AlphaEvolveSpec) :: AlphaMateSpec
    ! Inputs
    character(len=FILELENGTH) :: SpecFile, RelMtxFile, SelCriterionFile, GenderFile, SeedFile, &
                                 GenericIndCritFile, GenericMatCritFile, OutputBasename
    logical :: RelMtxGiven, NrmInsteadOfCoancestry, SelCriterionGiven, GenderGiven, SeedFileGiven, &
               GenericIndCritGiven, GenericMatCritGiven, SeedGiven, nThreadsGiven
    integer(int32) :: Seed, nThreads
    integer(int32) :: nGenericIndCrit, nGenericMatCrit

    ! Search mode specifications
    logical :: ModeMinCoancestry, ModeMinInbreeding, ModeMaxCriterion, ModeOpt, ModeRan, EvaluateFrontier
    integer(int32) :: nTargets
    character(len=SPECOPTIONLENGTH), allocatable :: AllTargets(:)
    real(FLOATTYPE), allocatable :: AllTargetValues(:)

    ! Targets
    logical :: TargetDegreeGiven,                                                             &
               TargetSelCriterionGiven, TargetSelCriterionStdGiven, TargetMaxCriterionPctGiven,  &
               TargetCoancestryGiven, TargetCoancestryRateGiven, TargetMinCoancestryPctGiven, &
               TargetInbreedingGiven, TargetInbreedingRateGiven, TargetMinInbreedingPctGiven
    real(FLOATTYPE), allocatable :: TargetDegree(:),                                                        &
                                    TargetSelCriterion(:), TargetSelCriterionStd(:), TargetMaxCriterionPct(:), &
                                    TargetCoancestry(:), TargetCoancestryRate(:), TargetMinCoancestryPct(:)
    real(FLOATTYPE) :: TargetInbreeding, TargetInbreedingRate, TargetMinInbreedingPct
    real(FLOATTYPE) :: CoancestryWeight, InbreedingWeight, MultipleMatingsWeight, SelfingWeight
    logical :: CoancestryWeightBelow, InbreedingWeightBelow
    real(FLOATTYPE), allocatable :: GenericIndCritWeight(:), GenericMatCritWeight(:)

    ! Biological specifications
    integer(int32) :: nInd, nMat, nPar, nPar1, nPar2 ! NOTE: nInd is here just for OO-flexibility (do not use it; the main one is in Data!!!)
    logical :: PreselectPar, PreselectPar1, PreselectPar2, &
               EqualizePar, EqualizePar1, EqualizePar2,    &
               LimitPar, LimitPar1, LimitPar2,             &
               MateAllocation, RandomMateAllocation,       &
               ReciprocalMatingsAllowed, MultipleMatingsAllowed, SelfingAllowed
    real(FLOATTYPE) :: LimitParMin,       LimitPar1Min,       LimitPar2Min, &
                       LimitParMax,       LimitPar1Max,       LimitPar2Max, &
                       LimitParMinWeight, LimitPar1MinWeight, LimitPar2MinWeight, &
                       PreselectParPct,   PreselectPar1Pct,   PreselectPar2Pct
    logical :: PAGEPar, PAGEPar1, PAGEPar2
    integer(int32) :: PreselectParN,   PreselectPar1N,     PreselectPar2N, &
                      PAGEParMax,      PAGEPar1Max,        PAGEPar2Max

    ! Algorithm specifications
    ! ... generic evolutionary parameters
    integer(int32) :: EvolAlgNSol, EvolAlgNIter, EvolAlgNIterStop, EvolAlgNIterPrint
    real(FLOATTYPEAH) :: EvolAlgStopTolCoancestry, EvolAlgStopTol
    logical :: EvolAlgLogPop
    character(len=SPECOPTIONLENGTH) :: EvolAlg
    ! ... differential evolution
    integer(int32) :: DiffEvolNIterBurnIn
    real(FLOATTYPEAH) :: DiffEvolParamCrBurnIn, DiffEvolParamCr1, DiffEvolParamCr2, &
                         DiffEvolParamFBase, DiffEvolParamFHigh1, DiffEvolParamFHigh2
    ! ... random search
    integer(int32) :: RanAlgStricter

    ! Modes
    ! @todo do we need ModeRanSpec?
    type(AlphaMateModeSpec) :: ModeSpec, ModeMinCoancestrySpec, ModeMinInbreedingSpec, ModeMaxCriterionSpec, ModeRanSpec, ModeOptSpec

    ! Column headers and formats for logging
    character(len=CHARLENGTH)      :: FmtLogStdoutHead,        &
                                      FmtLogStdout,            &
                                      FmtLogUnitHead,          &
                                      FmtLogUnit,              &
                                      FmtLogPopUnit,           &
                                      FmtContributionHead,     &
                                      FmtContributionHeadEdit, &
                                      FmtContribution,         &
                                      FmtContributionEdit,     &
                                      FmtMatingHead,           &
                                      FmtMating
    character(len=11), allocatable :: ColnameLogStdout(:)
    character(len=17), allocatable :: ColnameLogUnit(:), ColnameLogPopUnit(:)

    contains
      procedure :: Initialise => InitialiseAlphaMateSpec
      procedure :: Read       => ReadAlphaMateSpec
      procedure :: Write      => WriteAlphaMateSpec
      procedure :: SetupMode  => SetupModeAlphaMateSpec
      procedure :: SetupColNamesAndFormats
      procedure :: LogHead    => LogHeadAlphaMateSpec
      procedure :: LogPopHead => LogPopHeadAlphaMateSpec
  end type

  !> @brief AlphaMate data
  type, extends(AlphaEvolveData) :: AlphaMateData
    ! Raw data
    type(RelMat) :: Coancestry
    real(FLOATTYPE), allocatable :: SelCriterion(:),     SelCriterionStd(:),     &
                                    SelCriterionPAGE(:), SelCriterionStdPAGE(:), &
                                    AvgCoancestry(:),                         &
                                    GenericIndCrit(:, :), GenericMatCrit(:, :, :)
    integer(int32), allocatable :: Gender(:)
    ! Data summaries
    type(DESCSTATFLOAT) :: InbreedingStat, SelCriterionStat, SelCriterionStdStat, SelCriterionPAGEStat
    type(DESCSTATFLOAT), allocatable :: GenericIndCritStat(:)
    type(DESCSTATMATRIXFLOAT) :: CoancestryStat, CoancestryStatGender1, CoancestryStatGender2, CoancestryStatGenderDiff
    type(DESCSTATMATRIXFLOAT), allocatable :: GenericMatCritStat(:)
    ! Derived data
    integer(int32) :: nInd, nPotMat, nPotPar1, nPotPar2, nPotPar, nMal, nFem
    integer(int32), allocatable :: IdPotPar1(:), IdPotPar2(:), IdPotParSeq(:)
    real(FLOATTYPE) :: CoancestryRanMate, CoancestryRanMateNoSelf, CoancestryGenderMate, Inbreeding
    contains
      procedure :: Read  => ReadAlphaMateData
      procedure :: Write => WriteAlphaMateData
  end type

  !> @brief AlphaMate solution
  type, extends(AlphaEvolveSol) :: AlphaMateSol
    ! Solution results
    real(FLOATTYPE)              :: Penalty
    real(FLOATTYPE)              :: PenaltyCoancestry
    real(FLOATTYPE)              :: PenaltyInbreeding
    real(FLOATTYPE)              :: PenaltyLimitPar1
    real(FLOATTYPE)              :: PenaltyLimitPar2
    real(FLOATTYPE)              :: PenaltyMultipleMatings
    real(FLOATTYPE)              :: PenaltySelfing
    real(FLOATTYPE)              :: PenaltyGenericIndCrit
    real(FLOATTYPE)              :: PenaltyGenericMatCrit
    real(FLOATTYPE)              :: Degree
    real(FLOATTYPE)              :: SelCriterion
    real(FLOATTYPE)              :: SelCriterionStd
    real(FLOATTYPE)              :: MaxCriterionPct
    real(FLOATTYPE)              :: CoancestryRanMate
    real(FLOATTYPE)              :: CoancestryRateRanMate
    real(FLOATTYPE)              :: MinCoancestryPct
    real(FLOATTYPE)              :: Inbreeding
    real(FLOATTYPE)              :: InbreedingRate
    real(FLOATTYPE)              :: MinInbreedingPct
    real(FLOATTYPE), allocatable :: GenericIndCrit(:)
    real(FLOATTYPE), allocatable :: GenericMatCrit(:)
    integer(int32), allocatable  :: nVec(:)
    integer(int32), allocatable  :: MatingPlan(:, :)
    real(FLOATTYPE), allocatable :: GenomeEdit(:)

    contains
      procedure :: Initialise   => InitialiseAlphaMateSol
      procedure :: Assign       => AssignAlphaMateSol
      procedure :: UpdateMean   => UpdateMeanAlphaMateSol
      procedure :: Evaluate     => FixSolEtcMateAndEvaluateAlphaMateSol
      procedure :: Write        => WriteAlphaMateSol
      procedure :: WriteContributions
      procedure :: WriteMatingPlan
      procedure :: Log          => LogAlphaMateSol
      procedure :: LogPop       => LogPopAlphaMateSol
  end type

  !> @brief AlphaMate chromosome
  type :: AlphaMateChrom
    real(FLOATTYPE), allocatable :: ContPar1(:)
    real(FLOATTYPE), allocatable :: ContPar2(:)
    real(FLOATTYPE), allocatable :: MateRank(:)
    real(FLOATTYPE), allocatable :: EditPar1(:)
    real(FLOATTYPE), allocatable :: EditPar2(:)

    contains
      procedure :: Write        => WriteAlphaMateChrom
  end type

  contains

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Print AlphaMate title
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !> @return Printout
    !---------------------------------------------------------------------------
    subroutine AlphaMateTitle ! not pure due to IO
      implicit none
      write(STDOUT, "(a)") " "
      write(STDOUT, "(a)") "                            ***********************                           "
      write(STDOUT, "(a)") "                            *                     *                           "
      write(STDOUT, "(a)") "                            *      AlphaMate      *                           "
      write(STDOUT, "(a)") "                            *                     *                           "
      write(STDOUT, "(a)") "                            ***********************                           "
      write(STDOUT, "(a)") "                                                                              "
      write(STDOUT, "(a)") "       Software for optimising selection, maintenance of diversity, and       "
      write(STDOUT, "(a)") "                     mate allocation in breeding programs                     "
      write(STDOUT, "(a)") "                       http://AlphaGenes.Roslin.ed.ac.uk                      "
      write(STDOUT, "(a)") "                                 No liability                                 "
      write(STDOUT, "(a)") " "
      write(STDOUT, "(a)") "                       Commit:   "//TOSTRING(COMMIT)//"                       "
      write(STDOUT, "(a)") "                       Compiled: "//__DATE__//", "//__TIME__
      write(STDOUT, "(a)") " "
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Initialise AlphaMate specifications
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !---------------------------------------------------------------------------
    pure subroutine InitialiseAlphaMateSpec(This)
      implicit none
      class(AlphaMateSpec), intent(out)  :: This !< @return AlphaMateSpec holder

      ! Inputs

      This%OutputBasename = ""
      This%SpecFile = ""
      This%RelMtxFile = ""
      This%SelCriterionFile = ""
      This%GenderFile = ""
      This%GenericIndCritFile = ""
      This%nGenericIndCrit = 0
      This%GenericMatCritFile = ""
      This%nGenericMatCrit = 0
      This%SeedFile = ""

      This%RelMtxGiven = .false.
      This%NrmInsteadOfCoancestry = .false.
      This%SelCriterionGiven = .false.
      This%GenderGiven = .false.
      This%GenericIndCritGiven = .false.
      This%GenericMatCritGiven = .false.
      This%SeedFileGiven = .false.
      This%SeedGiven = .false.
      This%nThreads = 0
      This%nThreadsGiven = .false.

      ! Search mode specifications

      This%ModeMinCoancestry = .false.
      This%ModeMinInbreeding = .false.
      This%ModeMaxCriterion = .false.
      This%ModeOpt = .false.
      This%ModeRan = .false.
      This%EvaluateFrontier = .false.

      ! Targets

      This%nTargets = 0
      ! This%AllTargets ! allocatable so skip here
      ! This%AllTargetValues ! allocatable so skip here

      This%TargetDegreeGiven = .false.
      ! This%TargetDegree ! allocatable so skip here

      This%TargetSelCriterionGiven = .false.
      ! This%TargetSelCriterion ! allocatable so skip here
      This%TargetSelCriterionStdGiven = .false.
      ! This%TargetSelCriterionStd ! allocatable so skip here
      This%TargetMaxCriterionPctGiven = .false.
      ! This%TargetMaxCriterionPct ! allocatable so skip here

      This%TargetCoancestryGiven = .false.
      ! This%TargetCoancestry ! allocatable so skip here
      This%TargetCoancestryRateGiven = .false.
      ! This%TargetCoancestryRate ! allocatable so skip here
      This%CoancestryWeight = -1.0
      This%CoancestryWeightBelow = .false.
      This%TargetMinCoancestryPctGiven = .false.
      ! This%TargetMinCoancestryPct ! allocatable so skip here

      This%TargetInbreedingGiven = .false.
      This%TargetInbreeding = 1.0 ! set it high so that it does not have any effect
      This%TargetInbreedingRateGiven = .false.
      This%TargetInbreedingRate = 1.0 ! set it high so that it does not have any effect
      This%InbreedingWeight = 0.0 ! set it to zero so that it does not have any effect
      This%InbreedingWeightBelow = .false.
      This%TargetMinInbreedingPctGiven = .false.
      This%TargetMinInbreedingPct = 100.0 ! not used anyhow as TargetMinInbreedingPctGiven = .false.

      ! This%GenericIndCritWeight ! allocatable so skip here
      ! This%GenericMatCritWeight ! allocatable so skip here

      ! Biological specifications

      This%nInd = 0
      This%nMat = 0
      This%nPar = 0
      This%nPar1 = 0
      This%nPar2 = 0

      This%PreselectPar      = .false.
      This%PreselectPar1     = .false.
      This%PreselectPar2     = .false.
      This%PreselectParPct   = 100.0
      This%PreselectPar1Pct  = 100.0
      This%PreselectPar2Pct  = 100.0
      This%PreselectParN     = 0
      This%PreselectPar1N    = 0
      This%PreselectPar2N    = 0

      This%EqualizePar  = .false.
      This%EqualizePar1 = .false.
      This%EqualizePar2 = .false.

      This%LimitPar  = .false.
      This%LimitPar1 = .false.
      This%LimitPar2 = .false.
      This%LimitParMin  = 1.0
      This%LimitPar1Min = 1.0
      This%LimitPar2Min = 1.0
      This%LimitParMax  = huge(This%LimitParMax)  - 1.0
      This%LimitPar1Max = huge(This%LimitPar1Max) - 1.0
      This%LimitPar2Max = huge(This%LimitPar2Max) - 1.0
      This%LimitParMinWeight  = -1.0
      This%LimitPar1MinWeight = -1.0
      This%LimitPar2MinWeight = -1.0

      This%MateAllocation = .true.
      This%RandomMateAllocation = .false.

      This%ReciprocalMatingsAllowed = .false.

      This%MultipleMatingsAllowed = .false.
      This%MultipleMatingsWeight = -1.0

      This%SelfingAllowed = .false.
      This%SelfingWeight = -1.0

      This%PAGEPar  = .false.
      This%PAGEPar1 = .false.
      This%PAGEPar2 = .false.
      This%PAGEParMax  = 0
      This%PAGEPar1Max = 0
      This%PAGEPar2Max = 0

      ! Search algorithm specifications

      ! Ideally n >> p, say n = 2-20 (or even 50) * p
      ! However, when p is very large (say 1000), we need to consider computation resources/time.
      ! For high-dimensional settings n ~ 0.5*p
      ! Chen et al. (2014) Measuring the curse of dimensionality and its effects on particle swarm optimization
      ! and differential evolution https://link.springer.com/article/10.1007/s10489-014-0613-2; see also
      ! https://www.researchgate.net/post/What_is_the_optimal_recommended_population_size_for_differential_evolution2
      This%EvolAlgNSol = 0 ! this is set in AlphaMateSearch afer we know how many unknowns we have

      This%EvolAlgNIter = 100000
      This%EvolAlgNIterStop = 100
      This%EvolAlgNIterPrint = 100
      This%EvolAlgStopTolCoancestry = 0.001
      This%EvolAlgStopTol = 0.0010
      This%EvolAlgLogPop = .false.
      This%EvolAlg = "DE"

      This%DiffEvolNIterBurnIn = 100     ! 1000

      ! Cr [0, 1] should be low (<=0.1) for separable problems and large (0.9) for non-separable (epistatic) problems.
      ! Also, large values encourage exploration (large moves away from the target vector),
      ! while low values encourage exploitation (small moves away from the target vector, more closer to the base (a) vector).
      ! Montgomery and Chen (2010) An analysis of the operation of differential evolution at high and low
      ! crossover rates http://ieeexplore.ieee.org/document/5586128/
      This%DiffEvolParamCrBurnIn = 0.9 ! 0.4
      This%DiffEvolParamCr1 = 0.9      ! 0.2
      This%DiffEvolParamCr2 = 0.9      ! 0.2

      ! F should be [0 or 2/n, 1.2]
      ! Large values mean more exploration away from the base (a) vector. Should be large for large Cr.
      This%DiffEvolParamFBase  = 0.2   ! 0.1
      This%DiffEvolParamFHigh1 = 0.2   ! 1.0
      This%DiffEvolParamFHigh2 = 0.9   ! 4.0

      This%RanAlgStricter = 10
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Write AlphaMate specifications to a file or standard output
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !> @return Output to a file or standard output
    !---------------------------------------------------------------------------
    subroutine WriteAlphaMateSpec(This, File) ! not pure due to IO
      implicit none
      class(AlphaMateSpec), intent(in)       :: This !< AlphaMateSpec holder
      character(len=*), intent(in), optional :: File !< File (if missing use standard output)

      integer(int32) :: Unit
      if (present(File)) then
        open(newunit=Unit, file=File, action="write", status="unknown")
      else
        Unit = STDOUT
      end if

      ! Inputs

      write(Unit, *) "SpecFile: ",           trim(This%SpecFile)
      write(Unit, *) "RelMtxFile: ",         trim(This%RelMtxFile)
      write(Unit, *) "SelCriterionFile: ",   trim(This%SelCriterionFile)
      write(Unit, *) "GenderFile: ",         trim(This%GenderFile)
      write(Unit, *) "GenericIndCritFile: ", trim(This%GenericIndCritFile)
      write(Unit, *) "nGenericIndCrit: ",         This%nGenericIndCrit
      write(Unit, *) "GenericMatCritFile: ", trim(This%GenericMatCritFile)
      write(Unit, *) "nGenericMatCrit: ",         This%nGenericMatCrit
      write(Unit, *) "SeedFile: ",           trim(This%SeedFile)
      write(Unit, *) "Seed: ",                    This%Seed
      write(Unit, *) "nThreads: ",                This%nThreads
      write(Unit, *) "OutputBasename: ",     trim(This%OutputBasename)

      write(Unit, *) "RelMtxGiven: ",            This%RelMtxGiven
      write(Unit, *) "NrmInsteadOfCoancestry: ", This%NrmInsteadOfCoancestry
      write(Unit, *) "SelCriterionGiven: ",      This%SelCriterionGiven
      write(Unit, *) "GenderGiven: ",            This%GenderGiven
      write(Unit, *) "GenericIndCritGiven: ",    This%GenericIndCritGiven
      write(Unit, *) "GenericMatCritGiven: ",    This%GenericMatCritGiven
      write(Unit, *) "SeedFileGiven: ",          This%SeedFileGiven
      write(Unit, *) "SeedGiven: ",              This%SeedGiven
      write(Unit, *) "nThreadsGiven: ",          This%nThreadsGiven

      ! Search mode specifications

      write(Unit, *) "ModeMinCoancestry: ", This%ModeMinCoancestry
      write(Unit, *) "ModeMinInbreeding: ", This%ModeMinInbreeding
      write(Unit, *) "ModeMaxCriterion: ",  This%ModeMaxCriterion
      write(Unit, *) "ModeRan: ",           This%ModeRan
      write(Unit, *) "ModeOpt: ",           This%ModeOpt
      write(Unit, *) "EvaluateFrontier: ",  This%EvaluateFrontier

      ! Targets

      write(Unit, *) "nTargets: ", This%nTargets

      if (allocated(This%AllTargets)) then
        write(Unit, *) "AllTargets: ", This%AllTargets
      else
        write(Unit, *) "AllTargets: not allocated"
      end if
      if (allocated(This%AllTargetValues)) then
        write(Unit, *) "AllTargetValues: ", This%AllTargetValues
      else
        write(Unit, *) "AllTargetValues: not allocated"
      end if

      write(Unit, *) "TargetDegreeGiven: ",               This%TargetDegreeGiven
      if (allocated(This%TargetDegree)) then
        write(Unit, *) "TargetDegree: ",                  This%TargetDegree
      else
        write(Unit, *) "TargetDegree: not allocated"
      end if

      write(Unit, *) "TargetSelCriterionGiven: ",         This%TargetSelCriterionGiven
      if (allocated(This%TargetSelCriterion)) then
        write(Unit, *) "TargetSelCriterion: ",            This%TargetSelCriterion
      else
        write(Unit, *) "TargetSelCriterion: not allocated"
      end if
      write(Unit, *) "TargetSelCriterionStdGiven: ",      This%TargetSelCriterionStdGiven
      if (allocated(This%TargetSelCriterionStd)) then
        write(Unit, *) "TargetSelCriterionStd: ",         This%TargetSelCriterionStd
      else
        write(Unit, *) "TargetSelCriterionStd: not allocated"
      end if
      write(Unit, *) "TargetMaxCriterionPctGiven: ",      This%TargetMaxCriterionPctGiven
      if (allocated(This%TargetMaxCriterionPct)) then
        write(Unit, *) "TargetMaxCriterionPct: ",         This%TargetMaxCriterionPct
      else
        write(Unit, *) "TargetMaxCriterionPct: not allocated"
      end if

      write(Unit, *) "TargetCoancestryGiven: ",           This%TargetCoancestryGiven
      if (allocated(This%TargetCoancestry)) then
        write(Unit, *) "TargetCoancestry: ",              This%TargetCoancestry
      else
        write(Unit, *) "TargetCoancestry: not allocated"
      end if
      write(Unit, *) "TargetCoancestryRateGiven: ",       This%TargetCoancestryRateGiven
      if (allocated(This%TargetCoancestryRate)) then
        write(Unit, *) "TargetCoancestryRate: ",          This%TargetCoancestryRate
      else
        write(Unit, *) "TargetCoancestryRate: not allocated"
      end if
      write(Unit, *) "CoancestryWeight: ",                This%CoancestryWeight
      write(Unit, *) "CoancestryWeightBelow: ",           This%CoancestryWeightBelow
      write(Unit, *) "TargetMinCoancestryPctGiven: ",     This%TargetMinCoancestryPctGiven
      if (allocated(This%TargetMinCoancestryPct)) then
        write(Unit, *) "TargetMinCoancestryPct: ",        This%TargetMinCoancestryPct
      else
        write(Unit, *) "TargetMinCoancestryPct: not allocated"
      end if

      write(Unit, *) "TargetInbreedingGiven: ",           This%TargetInbreedingGiven
      write(Unit, *) "TargetInbreeding: ",                This%TargetInbreeding
      write(Unit, *) "TargetInbreedingRateGiven: ",       This%TargetInbreedingRateGiven
      write(Unit, *) "TargetInbreedingRate: ",            This%TargetInbreedingRate
      write(Unit, *) "InbreedingWeight: ",                This%InbreedingWeight
      write(Unit, *) "InbreedingWeightBelow: ",           This%InbreedingWeightBelow
      write(Unit, *) "TargetMinInbreedingPctGiven: ",     This%TargetMinInbreedingPctGiven
      write(Unit, *) "TargetMinInbreedingPct: ",          This%TargetMinInbreedingPct

      if (allocated(This%GenericIndCritWeight)) then
        write(Unit, *) "GenericIndCritWeight: ", This%GenericIndCritWeight
      else
        write(Unit, *) "GenericIndCritWeight: not allocated"
      end if
      if (allocated(This%GenericMatCritWeight)) then
        write(Unit, *) "GenericMatCritWeight: ", This%GenericMatCritWeight
      else
        write(Unit, *) "GenericMatCritWeight: not allocated"
      end if

      ! Biological specifications

      write(Unit, *) "nInd:  ", This%nInd
      write(Unit, *) "nMat:  ", This%nMat
      write(Unit, *) "nPar:  ", This%nPar
      write(Unit, *) "nPar1: ", This%nPar1
      write(Unit, *) "nPar2: ", This%nPar2

      write(Unit, *) "PreselectPar:     ", This%PreselectPar
      write(Unit, *) "PreselectPar1:    ", This%PreselectPar1
      write(Unit, *) "PreselectPar2:    ", This%PreselectPar2
      write(Unit, *) "PreselectParPct:  ", This%PreselectParPct
      write(Unit, *) "PreselectPar1Pct: ", This%PreselectPar1Pct
      write(Unit, *) "PreselectPar2Pct: ", This%PreselectPar2Pct
      write(Unit, *) "PreselectParN:    ", This%PreselectParN
      write(Unit, *) "PreselectPar1N:   ", This%PreselectPar1N
      write(Unit, *) "PreselectPar2N:   ", This%PreselectPar2N

      write(Unit, *) "EqualizePar:  ", This%EqualizePar
      write(Unit, *) "EqualizePar1: ", This%EqualizePar1
      write(Unit, *) "EqualizePar2: ", This%EqualizePar2

      write(Unit, *) "LimitPar:           ", This%LimitPar
      write(Unit, *) "LimitPar1:          ", This%LimitPar1
      write(Unit, *) "LimitPar2:          ", This%LimitPar2
      write(Unit, *) "LimitParMin:        ", This%LimitParMin
      write(Unit, *) "LimitPar1Min:       ", This%LimitPar1Min
      write(Unit, *) "LimitPar2Min:       ", This%LimitPar2Min
      write(Unit, *) "LimitParMax:        ", This%LimitParMax
      write(Unit, *) "LimitPar1Max:       ", This%LimitPar1Max
      write(Unit, *) "LimitPar2Max:       ", This%LimitPar2Max
      write(Unit, *) "LimitParMinWeight:  ", This%LimitParMinWeight
      write(Unit, *) "LimitPar1MinWeight: ", This%LimitPar1MinWeight
      write(Unit, *) "LimitPar2MinWeight: ", This%LimitPar2MinWeight

      write(Unit, *) "MateAllocation:       ", This%MateAllocation
      write(Unit, *) "RandomMateAllocation: ", This%RandomMateAllocation

      write(Unit, *) "ReciprocalMatingsAllowed: ", This%ReciprocalMatingsAllowed

      write(Unit, *) "MultipleMatingsAllowed: ", This%MultipleMatingsAllowed
      write(Unit, *) "MultipleMatingsWeight:  ", This%MultipleMatingsWeight

      write(Unit, *) "SelfingAllowed: ", This%SelfingAllowed
      write(Unit, *) "SelfingWeight:  ", This%SelfingWeight

      write(Unit, *) "PAGEPar:     ", This%PAGEPar
      write(Unit, *) "PAGEPar1:    ", This%PAGEPar1
      write(Unit, *) "PAGEPar2:    ", This%PAGEPar2
      write(Unit, *) "PAGEParMax:  ", This%PAGEParMax
      write(Unit, *) "PAGEPar1Max: ", This%PAGEPar1Max
      write(Unit, *) "PAGEPar2Max: ", This%PAGEPar2Max

      ! Algorithm specifications

      write(Unit, *) "EvolAlgNSol:              ", This%EvolAlgNSol
      write(Unit, *) "EvolAlgNIter:             ", This%EvolAlgNIter
      write(Unit, *) "EvolAlgNIterStop:         ", This%EvolAlgNIterStop
      write(Unit, *) "EvolAlgNIterPrint:        ", This%EvolAlgNIterPrint
      write(Unit, *) "EvolAlgStopTolCoancestry: ", This%EvolAlgStopTolCoancestry
      write(Unit, *) "EvolAlgStopTol:           ", This%EvolAlgStopTol
      write(Unit, *) "EvolAlgLogPop:            ", This%EvolAlgLogPop
      write(Unit, *) "EvolAlg:                  ", This%EvolAlg

      write(Unit, *) "DiffEvolNIterBurnIn:   ", This%DiffEvolNIterBurnIn
      write(Unit, *) "DiffEvolParamCrBurnIn: ", This%DiffEvolParamCrBurnIn
      write(Unit, *) "DiffEvolParamCr1:      ", This%DiffEvolParamCr1
      write(Unit, *) "DiffEvolParamCr2:      ", This%DiffEvolParamCr2
      write(Unit, *) "DiffEvolParamFBase:    ", This%DiffEvolParamFBase
      write(Unit, *) "DiffEvolParamFHigh1:   ", This%DiffEvolParamFHigh1
      write(Unit, *) "DiffEvolParamFHigh2:   ", This%DiffEvolParamFHigh2

      write(Unit, *) "RanAlgStricter: ", This%RanAlgStricter

      ! Data&Spec derived quantities

      write(Unit, *) "ModeSpec: "
      call This%ModeSpec%Write(Unit)
      write(Unit, *) "ModeMinCoancestrySpec: "
      call This%ModeMinCoancestrySpec%Write(Unit)
      write(Unit, *) "ModeMinInbreedingSpec: "
      call This%ModeMinInbreedingSpec%Write(Unit)
      write(Unit, *) "ModeMaxCriterionSpec: "
      call This%ModeMaxCriterionSpec%Write(Unit)
      ! @todo do we need ModeRanSpec?
      write(Unit, *) "ModeRanSpec: "
      call This%ModeRanSpec%Write(Unit)
      write(Unit, *) "ModeOptSpec: "
      call This%ModeOptSpec%Write(Unit)

      ! Formats for logging

      write(Unit, *) "FmtLogStdoutHead:        ", trim(This%FmtLogStdoutHead)
      write(Unit, *) "FmtLogStdout:            ", trim(This%FmtLogStdout)
      if (allocated(This%ColnameLogStdout)) then
        write(Unit, *) "ColnameLogStdout:        ", This%ColnameLogStdout
      else
        write(Unit, *) "ColnameLogStdout: not allocated"
      end if
      write(Unit, *) "FmtLogUnitHead:            ", trim(This%FmtLogUnitHead)
      if (allocated(This%ColnameLogUnit)) then
        write(Unit, *) "ColnameLogUnit:          ", This%ColnameLogUnit
      else
        write(Unit, *) "ColnameLogUnit: not allocated"
      end if
      write(Unit, *) "FmtLogPopUnit:           ", trim(This%FmtLogPopUnit)
      if (allocated(This%ColnameLogPopUnit)) then
        write(Unit, *) "ColnameLogPopUnit:       ", This%ColnameLogPopUnit
      else
        write(Unit, *) "ColnameLogPopUnit: not allocated"
      end if
      write(Unit, *) "FmtContributionHead:     ", trim(This%FmtContributionHead)
      write(Unit, *) "FmtContributionHeadEdit: ", trim(This%FmtContributionHeadEdit)
      write(Unit, *) "FmtContribution:         ", trim(This%FmtContribution)
      write(Unit, *) "FmtContributionEdit:     ", trim(This%FmtContributionEdit)
      write(Unit, *) "FmtMatingHead:           ", trim(This%FmtMatingHead)
      write(Unit, *) "FmtMating:               ", trim(This%FmtMating)

      if (present(File)) then
        close(Unit)
      end if
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Read AlphaMate specifications from a file
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !---------------------------------------------------------------------------
    subroutine ReadAlphaMateSpec(This, SpecFile, LogStdout) ! not pure due to IO
      implicit none
      class(AlphaMateSpec), intent(out) :: This      !< @return AlphaMateSpec holder
      character(len=*), intent(in)      :: SpecFile  !< Spec file; when missing, a stub with defaults is created
      logical, optional                 :: LogStdout !< Log process on stdout (default .false.)

      ! Other
      character(len=:), allocatable :: DumString
      character(len=SPECOPTIONLENGTH) :: Line
      character(len=SPECOPTIONLENGTH) :: First
      character(len=SPECOPTIONLENGTH), allocatable, dimension(:) :: Second

      integer(int32) :: SpecUnit, Stat, nGenericIndCrit, nGenericMatCrit

      logical :: LogStdoutInternal

      if (present(LogStdout)) then
        LogStdoutInternal = LogStdout
      else
        LogStdoutInternal = .false.
      end if

      if (LogStdoutInternal) then
        write(STDOUT, "(a)") " "
      end if

      ! Defaults
      call This%Initialise

      This%SpecFile = SpecFile
      open(newunit=SpecUnit, file=This%SpecFile, action="read", status="old")

      Stat = 0
      ReadSpec: do while (Stat .eq. 0)
        read(SpecUnit, "(a)", iostat=Stat) Line
        if (len_trim(Line) .eq. 0) then
          cycle
        end if
        call SplitLineIntoTwoParts(trim(adjustl(Line)), First, Second)
        DumString = ParseToFirstWhitespace(First)
        if (First(1:1) .eq. "=") then
          cycle
        else
          select case (ToLower(trim(DumString)))
            ! Inputs
            case ("outputbasename")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .ne. "none") then
                  write(This%OutputBasename, *) trim(adjustl(Second(1)))
                  This%OutputBasename = adjustl(This%OutputBasename)
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Output basename: "//trim(This%OutputBasename)
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a string for OutputBasename, for example, OutputBasename, AnalysisX"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("coancestrymatrixfile")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .ne. "none") then
                  write(This%RelMtxFile, *) trim(adjustl(Second(1)))
                  This%RelMtxFile = adjustl(This%RelMtxFile)
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Coancestry matrix file: "//trim(This%RelMtxFile)
                  end if
                  This%RelMtxGiven = .true.
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a file for CoancestryMatrixFile, for example, CoancestryMatrixFile, CoaMtx.txt"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("nrmmatrixfile")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .ne. "none") then
                  write(This%RelMtxFile, *) trim(adjustl(Second(1)))
                  This%RelMtxFile = adjustl(This%RelMtxFile)
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Numerator relationship matrix file: "//trim(This%RelMtxFile)
                  end if
                  This%RelMtxGiven = .true.
                  This%NrmInsteadOfCoancestry = .true.
                else
                  write(STDERR, "(a)") " ERROR: Must specify a file for NrmMatrixFile, for example, NrmMatrixFile, NrmMtx.txt"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a file for NrmMatrixFile, for example, NrmMatrixFile, NrmMtx.txt"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("selcriterionfile")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .ne. "none") then
                  write(This%SelCriterionFile, *) trim(adjustl(Second(1)))
                  This%SelCriterionFile = adjustl(This%SelCriterionFile)
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Selection criterion file: "//trim(This%SelCriterionFile)
                  end if
                  This%SelCriterionGiven = .true.
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a file for SelCriterionFile, for example, SelCriterionFile, SelCrit.txt"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("genderfile")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .ne. "none") then
                  write(This%GenderFile, *) trim(adjustl(Second(1)))
                  This%GenderFile = adjustl(This%GenderFile)
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Gender file: "//trim(This%GenderFile)
                  end if
                  This%GenderGiven = .true.
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a file for GenderFile, for example, GenderFile, Gender.txt"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("seedfile")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .ne. "none") then
                  write(This%SeedFile, *) trim(adjustl(Second(1)))
                  This%SeedFileGiven = .true.
                  This%SeedFile = adjustl(This%SeedFile)
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Seed file: "//trim(This%SeedFile)
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a file for SeedFile, for example, SeedFile, Seed.txt"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("seed")
              if (allocated(Second)) then
                This%SeedGiven = .true.
                This%Seed = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Seed: "//trim(Int2Char(This%Seed))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for Seed, for example, Seed, 19791123"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("numberofthreads")
              if (allocated(Second)) then
                This%nThreadsGiven = .true.
                This%nThreads = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " NumberOfThreads: "//trim(Int2Char(This%nThreads))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfThreads, for example, NumberOfThreads, 8"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("genericindividualcriterionfile")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .ne. "none") then
                  write(This%GenericIndCritFile, *) trim(adjustl(Second(1)))
                  This%GenericIndCritFile = adjustl(This%GenericIndCritFile)
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Generic individual criterion file: "//trim(This%GenericIndCritFile)
                  end if
                  This%GenericIndCritGiven = .true.
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a file for GenericIndividualCriterionFile, for example, GenericIndividualCriterionFile, IndividualCriterion.txt"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("genericindividualcriterioncolumns")
              if (This%GenericIndCritGiven) then
                if (allocated(Second)) then
                  This%nGenericIndCrit = Char2Int(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Generic individual criterion - number of columns: "//trim(Int2Char(This%nGenericIndCrit))
                  end if
                  allocate(This%GenericIndCritWeight(This%nGenericIndCrit))
                  nGenericIndCrit = 0
                else
                  write(STDERR, "(a)") " ERROR: Must specify a number for GenericIndividualCriterionColumns, for example, GenericIndividualCriterionColumns, 2"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              end if

            case ("genericindividualcriterionweight")
              if (This%GenericIndCritGiven) then
                if (allocated(Second)) then
                  nGenericIndCrit = nGenericIndCrit + 1
                  This%GenericIndCritWeight(nGenericIndCrit) = Char2Real(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Generic individual criterion - weight ("//trim(Int2Char(nGenericIndCrit))//"): "//trim(Real2Char(This%GenericIndCritWeight(nGenericIndCrit), fmt=FMTREAL2CHAR))
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify a value for GenericIndividualCriterionWeight, for example, GenericIndividualCriterionWeight, 2"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              end if

            case ("genericmatingcriterionfile")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .ne. "none") then
                  write(This%GenericMatCritFile, *) trim(adjustl(Second(1)))
                  This%GenericMatCritFile = adjustl(This%GenericMatCritFile)
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Generic mating/crossing criterion file: "//trim(This%GenericMatCritFile)
                  end if
                  This%GenericMatCritGiven = .true.
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a file for GenericMatingCriterionFile, for example, GenericMatingCriterionFile, MatingCriterion.txt"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("genericmatingcriterioncolumns")
              if (This%GenericMatCritGiven) then
                if (allocated(Second)) then
                  This%nGenericMatCrit = Char2Int(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Generic mating/crossing criterion - number of columns: "//trim(Int2Char(This%nGenericMatCrit))
                  end if
                  allocate(This%GenericMatCritWeight(This%nGenericMatCrit))
                  nGenericMatCrit = 0
                else
                  write(STDERR, "(a)") " ERROR: Must specify a number for GenericMatingCriterionColumns, for example, GenericMatingCriterionColumns, 2"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              end if

            case ("genericmatingcriterionweight")
              if (This%GenericMatCritGiven) then
                if (allocated(Second)) then
                  nGenericMatCrit = nGenericMatCrit + 1
                  This%GenericMatCritWeight(nGenericMatCrit) = Char2Real(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Generic mating/crossing criterion - weight ("//trim(Int2Char(nGenericMatCrit))//"): "//trim(Real2Char(This%GenericMatCritWeight(nGenericMatCrit), fmt=FMTREAL2CHAR))
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify a value for GenericMatingCriterionWeight, for example, GenericMatingCriterionWeight, 2"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              end if

            ! Search mode specifications
            case ("modemincoancestry")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%ModeMinCoancestry = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " ModeMinCoancestry"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for ModeMinCoancestry, for example, ModeMinCoancestry, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("modemininbreeding")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%ModeMinInbreeding = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " ModeMinInbreeding"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for ModeMinInbreeding, for example, ModeMinInbreeding, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("modemaxcriterion")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%ModeMaxCriterion = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " ModeMaxCriterion"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for ModeMaxCriterion, for example, ModeMaxCriterion, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            ! case ("moderan")
            !   if (allocated(Second)) then
            !     if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
            !       This%ModeRan = .true.
            !       if (LogStdoutInternal) then
            !         write(STDOUT, "(a)") " ModeRan"
            !       end if
            !     end if
            !   else
            !     write(STDERR, "(a)") " ERROR: Must specify Yes or No for ModeRan, for example, ModeRan, Yes"
            !     write(STDERR, "(a)") " "
            !     stop 1
            !   end if

            case ("modeopt")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%ModeOpt = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " ModeOpt"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for ModeOpt, for example, ModeOpt, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("evaluatefrontier")
              if (allocated(Second)) then
                This%ModeOpt = .true.
                This%ModeMaxCriterion = .true.
                This%ModeMinCoancestry = .true.
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%EvaluateFrontier = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " EvaluateFrontier"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for EvaluateFrontier, for example, EvaluateFrontier, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            ! Targets
            case ("targetdegree")
              if (allocated(Second)) then
                This%ModeOpt = .true.
                This%ModeMaxCriterion = .true.
                This%ModeMinCoancestry = .true.
                This%TargetDegreeGiven = .true.
                This%nTargets = This%nTargets + 1
                block
                  integer(int32) :: n
                  real(FLOATTYPE) :: Tmp
                  Tmp = Char2Real(trim(adjustl(Second(1))))
                  call Append(x=This%TargetDegree, y=Tmp)
                  n = size(This%TargetDegree)
                  call Append(x=This%AllTargetValues, y=Tmp)
                  call Append(x=This%AllTargets,      y="Degree", Len=SPECOPTIONLENGTH)
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Targeted degree: "//trim(Real2Char(This%TargetDegree(n), fmt=FMTREAL2CHAR))
                  end if
                  if ((This%TargetDegree(n) .lt. 0.0) .or. (This%TargetDegree(n) .gt. 90.0)) then
                    write(STDERR, "(a)") " ERROR: TargetDegree must be between 0 and 90!"
                    write(STDERR, "(a)") " "
                    stop 1
                  end if
                end block
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetDegree, for example, TargetDegree, 45"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("targetselcriterion")
              if (allocated(Second)) then
                This%ModeOpt = .true.
                This%ModeMaxCriterion = .true.
                This%TargetSelCriterionGiven = .true.
                This%nTargets = This%nTargets + 1
                block
                  integer(int32) :: n
                  real(FLOATTYPE) :: Tmp
                  Tmp = Char2Real(trim(adjustl(Second(1))))
                  call Append(x=This%TargetSelCriterion, y=Tmp)
                  n = size(This%TargetSelCriterion)
                  call Append(x=This%AllTargetValues, y=Tmp)
                  call Append(x=This%AllTargets,      y="SelCriterion", Len=SPECOPTIONLENGTH)
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Targeted selection criterion: "//trim(Real2Char(This%TargetSelCriterion(n), fmt=FMTREAL2CHAR))
                  end if
                end block
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetSelCriterion, for example, TargetSelCriterion, 100"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("targetselcriterionstd")
              if (allocated(Second)) then
                This%ModeOpt = .true.
                This%ModeMaxCriterion = .true.
                This%TargetSelCriterionStdGiven = .true.
                This%nTargets = This%nTargets + 1
                block
                  integer(int32) :: n
                  real(FLOATTYPE) :: Tmp
                  Tmp = Char2Real(trim(adjustl(Second(1))))
                  call Append(x=This%TargetSelCriterionStd, y=Tmp)
                  n = size(This%TargetSelCriterionStd)
                  call Append(x=This%AllTargetValues, y=Tmp)
                  call Append(x=This%AllTargets,      y="SelCriterionStd", Len=SPECOPTIONLENGTH)
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Targeted stand. selection criterion: "//trim(Real2Char(This%TargetSelCriterionStd(n), fmt=FMTREAL2CHAR))
                  end if
                  if ((This%TargetSelCriterionStd(n) .lt. -5.0) .or. (This%TargetSelCriterionStd(n) .gt. 5.0)) then
                    write(STDERR, "(a)") " ERROR: TargetSelCriterionStd must be above -5 and below 5!"
                    write(STDERR, "(a)") " "
                    stop 1
                  end if
                end block
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetSelCriterionStd, for example, TargetSelCriterionStd, 2"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("targetmaxcriterionpct")
              if (allocated(Second)) then
                This%ModeOpt = .true.
                This%ModeMaxCriterion = .true.
                This%TargetMaxCriterionPctGiven = .true.
                This%nTargets = This%nTargets + 1
                block
                  integer(int32) :: n
                  real(FLOATTYPE) :: Tmp
                  Tmp = Char2Real(trim(adjustl(Second(1))))
                  call Append(x=This%TargetMaxCriterionPct, y=Tmp)
                  n = size(This%TargetMaxCriterionPct)
                  call Append(x=This%AllTargetValues, y=Tmp)
                  call Append(x=This%AllTargets,      y="MaxCriterionPct", Len=SPECOPTIONLENGTH)
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Targeted percentage of maximum criterion: "//trim(Real2Char(This%TargetMaxCriterionPct(n), fmt=FMTREAL2CHAR))
                  end if
                  if ((This%TargetMaxCriterionPct(n) .lt. 0.0) .or. (This%TargetMaxCriterionPct(n) .gt. 100.0)) then
                    write(STDERR, "(a)") " ERROR: TargetMaxCriterionPct must be between 0 and 100!"
                    write(STDERR, "(a)") " "
                    stop 1
                  end if
                end block
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetMaxCriterionPct, for example, TargetMaxCriterionPct, 90"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("targetcoancestry")
              if (allocated(Second)) then
                This%ModeOpt = .true.
                This%ModeMinCoancestry = .true.
                This%TargetCoancestryGiven = .true.
                This%nTargets = This%nTargets + 1
                block
                  integer(int32) :: n
                  real(FLOATTYPE) :: Tmp
                  Tmp = Char2Real(trim(adjustl(Second(1))))
                  call Append(x=This%TargetCoancestry, y=Tmp)
                  n = size(This%TargetCoancestry)
                  call Append(x=This%AllTargetValues, y=Tmp)
                  call Append(x=This%AllTargets,      y="Coancestry", Len=SPECOPTIONLENGTH)
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Targeted coancestry: "//trim(Real2Char(This%TargetCoancestry(n), fmt=FMTREAL2CHAR))
                  end if
                  if ((This%TargetCoancestry(n) .lt. -1.0) .or. (This%TargetCoancestry(n) .gt. 1.0)) then
                    write(STDERR, "(a)") " ERROR: TargetCoancestry must be between -1 and +1!"
                    write(STDERR, "(a)") " "
                    stop 1
                  end if
                end block
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetCoancestry, for example, TargetCoancestry, 0.31"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("targetcoancestryrate")
              if (allocated(Second)) then
                This%ModeOpt = .true.
                This%ModeMinCoancestry = .true.
                This%TargetCoancestryRateGiven = .true.
                This%nTargets = This%nTargets + 1
                block
                  integer(int32) :: n
                  real(FLOATTYPE) :: Tmp
                  Tmp = Char2Real(trim(adjustl(Second(1))))
                  call Append(x=This%TargetCoancestryRate, y=Tmp)
                  n = size(This%TargetCoancestryRate)
                  call Append(x=This%AllTargetValues, y=Tmp)
                  call Append(x=This%AllTargets,      y="CoancestryRate", Len=SPECOPTIONLENGTH)
                  This%TargetCoancestryRate(n) = Char2Real(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Targeted rate of coancestry: "//trim(Real2Char(This%TargetCoancestryRate(n), fmt=FMTREAL2CHAR))
                  end if
                  if ((This%TargetCoancestryRate(n) .lt. -1.0) .or. (This%TargetCoancestryRate(n) .gt. 1.0)) then
                    write(STDERR, "(a)") " ERROR: TargetCoancestryRate must be between -1 and +1!"
                    write(STDERR, "(a)") " "
                    stop 1
                  end if
                end block
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetCoancestryRate, for example, TargetCoancestryRate, 0.01"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("targetmincoancestrypct")
              if (allocated(Second)) then
                This%ModeOpt = .true.
                This%ModeMinCoancestry = .true.
                This%TargetMinCoancestryPctGiven = .true.
                This%nTargets = This%nTargets + 1
                block
                  integer(int32) :: n
                  real(FLOATTYPE) :: Tmp
                  Tmp = Char2Real(trim(adjustl(Second(1))))
                  call Append(x=This%TargetMinCoancestryPct, y=Tmp)
                  n = size(This%TargetMinCoancestryPct)
                  call Append(x=This%AllTargetValues, y=Tmp)
                  call Append(x=This%AllTargets,      y="MinCoancestryPct", Len=SPECOPTIONLENGTH)
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Targeted percentage of minimum coancestry: "//trim(Real2Char(This%TargetMinCoancestryPct(n), fmt=FMTREAL2CHAR))
                  end if
                  if ((This%TargetMinCoancestryPct(n) .lt. 0.0) .or. (This%TargetMinCoancestryPct(n) .gt. 100.0)) then
                    write(STDERR, "(a)") " ERROR: TargetMinCoancestryPct must be between 0 and 100!"
                    write(STDERR, "(a)") " "
                    stop 1
                  end if
                end block
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetMinCoancestryPct, for example, TargetMinCoancestryPct, 90"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("coancestryweight")
              if (allocated(Second)) then
                This%CoancestryWeight = Char2Real(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Targeted coancestry - weight: "//trim(Real2Char(This%CoancestryWeight, fmt=FMTREAL2CHAR))
                end if
                if (This%CoancestryWeight .gt. 0.0) then
                  write(STDOUT, "(a)") " NOTE: Positive weight for the targeted coancestry, that is, encourage higher coancestry. Was this intended?"
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for CoancestryWeight, for example, CoancestryWeight, -2"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("coancestryweightbelow")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%CoancestryWeightBelow = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Targeted coancestry - weight also values below the target"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for CoancestryWeightBelow, for example, CoancestryWeightBelow, No"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("targetinbreeding")
              if (allocated(Second)) then
                This%ModeOpt = .true.
                This%ModeMinInbreeding = .true.
                This%TargetInbreedingGiven = .true.
                ! This%nTargets = This%nTargets + 1
                ! block
                !   @todo
                ! end block
                This%TargetInbreeding = Char2Real(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Targeted inbreding: "//trim(Real2Char(This%TargetInbreeding, fmt=FMTREAL2CHAR))
                end if
                if ((This%TargetInbreeding .lt. -1.0) .or. (This%TargetInbreeding .gt. 1.0)) then
                  write(STDERR, "(a)") " ERROR: TargetInbreeding must be between -1 and +1!"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetInbreeding, for example, TargetInbreeding, 0.31"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("targetinbreedingrate")
              if (allocated(Second)) then
                This%ModeOpt = .true.
                This%ModeMinInbreeding = .true.
                This%TargetInbreedingRateGiven = .true.
                ! This%nTargets = This%nTargets + 1
                ! block
                !   @todo
                ! end block
                This%TargetInbreedingRate = Char2Real(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Targeted rate of inbreeding: "//trim(Real2Char(This%TargetInbreedingRate, fmt=FMTREAL2CHAR))
                end if
                if ((This%TargetInbreedingRate .lt. -1.0) .or. (This%TargetInbreedingRate .gt. 1.0)) then
                  write(STDERR, "(a)") " ERROR: TargetInbreedingRate must be between -1 and +1!"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetInbreedingRate, for example, TargetInbreedingRate, 0.01"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("targetmininbreedingpct")
              if (allocated(Second)) then
                This%ModeOpt = .true.
                This%ModeMinInbreeding = .true.
                This%TargetMinInbreedingPctGiven = .true.
                ! This%nTargets = This%nTargets + 1
                !   @todo
                ! end block
                This%TargetMinInbreedingPct = Char2Real(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Targeted percentage of minimum inbreeding: "//trim(Real2Char(This%TargetMinInbreedingPct, fmt=FMTREAL2CHAR))
                end if
                if ((This%TargetMinInbreedingPct .lt. 0.0) .or. (This%TargetMinInbreedingPct .gt. 100.0)) then
                  write(STDERR, "(a)") " ERROR: TargetMinInbreedingPct must be between 0 and 100!"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetMinInbreedingPct, for example, TargetMinInbreedingPct, 90"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("inbreedingweight")
              if (allocated(Second)) then
                This%InbreedingWeight = Char2Real(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Targeted inbreeding - weight: "//trim(Real2Char(This%InbreedingWeight, fmt=FMTREAL2CHAR))
                end if
                if (This%InbreedingWeight .gt. 0.0) then
                  write(STDOUT, "(a)") " NOTE: Positive weight for targeted inbreeding, that is, encourage higher inbreeding. Was this intended?"
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for InbreedingWeight, for example, InbreedingWeight, -2"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("inbreedingweightbelow")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%InbreedingWeightBelow = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Targeted inbreeding - weight also values below the target"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for InbreedingWeightBelow, for example, InbreedingWeightBelow, No"
                write(STDERR, "(a)") " "
                stop 1
              end if

            ! Biological specifications
            case ("numberofmatings")
              if (allocated(Second)) then
                This%nMat = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Number of matings/crosses: "//trim(Int2Char(This%nMat))
                end if
                if (This%nMat .lt. 1) then
                  write(STDERR, "(a)") " ERROR: NumberOfMatings must be at least 1"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfMatings, for example, NumberOfMatings, 10"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("numberofparents")
              if (allocated(Second)) then
                This%nPar = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Number of parents: "//trim(Int2Char(This%nPar))
                end if
                if (This%nPar .lt. 1) then
                  write(STDERR, "(a)") " ERROR: NumberOfParents must be at least 1"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfParents, for example, NumberOfParents, 20"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("numberofmaleparents")
              if (allocated(Second)) then
                This%nPar1 = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Number of male parents: "//trim(Int2Char(This%nPar1))
                end if
                if (This%nPar1 .lt. 1) then
                  write(STDERR, "(a)") " ERROR: NumberOfMaleParents must be at least 1"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfMaleParents, for example, NumberOfMaleParents, 10"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("numberoffemaleparents")
              if (allocated(Second)) then
                This%nPar2 = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Number of female parents: "//trim(Int2Char(This%nPar2))
                end if
                if (This%nPar2 .lt. 1) then
                  write(STDERR, "(a)") " ERROR: NumberOfFemaleParents must be at least 1"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfFemaleParents, for example, NumberOfFemaleParents, 10"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("preselect")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%PreselectPar = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Preselect"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for Preselect, for example, Preselect, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("preselectpercentage")
              if (allocated(Second)) then
                This%PreselectParPct = Char2Real(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Preselect percentage: "//trim(Real2Char(This%PreselectParPct, fmt=FMTREAL2CHAR))
                end if
                if (This%PreselectParPct .lt. 0.0 .or. This%PreselectParPct .gt. 100.) then
                  write(STDERR, "(a)") " ERROR: The percentage value must be between 0 and 100"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for PreselectPercentage, for example, PreselectPercentage, 10"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("preselectmales")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%PreselectPar  = .true.
                  This%PreselectPar1 = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Preselect males"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for PreselectMales, for example, PreselectMales, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("preselectmalespercentage")
              if (allocated(Second)) then
                This%PreselectPar1Pct = Char2Real(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Preselect males percentage: "//trim(Real2Char(This%PreselectPar1Pct, fmt=FMTREAL2CHAR))
                end if
                if (This%PreselectPar1Pct .lt. 0.0 .or. This%PreselectPar1Pct .gt. 100.) then
                  write(STDERR, "(a)") " ERROR: The percentage value must be between 0 and 100"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for PreselectMalesPercentage, for example, PreselectMalesPercentage, 10"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("preselectfemales")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%PreselectPar  = .true.
                  This%PreselectPar2 = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Preselect females"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for PreselectFemales, for example, PreselectFemales, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("preselectfemalespercentage")
              if (allocated(Second)) then
                This%PreselectPar2Pct = Char2Real(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Preselect females percentage: "//trim(Real2Char(This%PreselectPar2Pct, fmt=FMTREAL2CHAR))
                end if
                if (This%PreselectPar2Pct .lt. 0.0 .or. This%PreselectPar2Pct .gt. 100.) then
                  write(STDERR, "(a)") " ERROR: The percentage value must be between 0 and 100"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for PreselectFemalesPercentage, for example, PreselectFemalesPercentage, 10"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("equalizecontributions")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%EqualizePar = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Equalize contributions"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for EqualizeContributions, for example, EqualizeContributions, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("equalizemalecontributions")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%EqualizePar  = .true.
                  This%EqualizePar1 = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Equalize contributions of males"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for EqualizeMaleContributions, for example, EqualizeMaleContributions, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("equalizefemalecontributions")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%EqualizePar  = .true.
                  This%EqualizePar2 = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Equalize contributions of females"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for EqualizeFemaleContributions, for example, EqualizeFemaleContributions, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("limitcontributions")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%LimitPar = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Limit contributions"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for LimitContributions, for example, LimitContributions, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("limitcontributionsmin")
              if (allocated(Second)) then
                This%LimitParMin = Char2Real(trim(adjustl(Second(1)))) ! real because of continous solution representation
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Limit contributions - minimum: "//trim(Int2Char(nint(This%LimitParMin))) ! nint because of continous solution representation
                end if
                if (This%LimitParMin .lt. 1.0) then
                  write(STDERR, "(a)") " ERROR: LimitContributionsMin must be at least 1"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for LimitContributionsMin, for example, LimitContributionsMin, 1"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("limitcontributionsmax")
              if (allocated(Second)) then
                This%LimitParMax = Char2Real(trim(adjustl(Second(1)))) ! Real because of continous solution representation
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Limit contributions - maximum: "//trim(Int2Char(nint(This%LimitParMax))) ! nint because of continous solution representation
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for LimitContributionsMax, for example, LimitContributionsMax, 10"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("limitcontributionsminweight")
              if (allocated(Second)) then
                This%LimitParMinWeight = Char2Real(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Limit contributions - weight for contributions below minimum: "//trim(Real2Char(This%LimitParMinWeight, fmt=FMTREAL2CHAR))
                end if
                if (This%LimitParMinWeight .gt. 0.0) then
                  write(STDOUT, "(a)") " NOTE: Positive weight for limit on minimum contributions, that is, encourage smaller contributions than defined minimum. Was this intended?"
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for LimitContributionsMinWeight, for example, LimitContributionsMinWeight, -2"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("limitmalecontributions")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%LimitPar  = .true.
                  This%LimitPar1 = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Limit contributions of males"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for LimitMaleContributions, for example, LimitMaleContributions, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("limitmalecontributionsmin")
              if (allocated(Second)) then
                This%LimitPar1Min = Char2Real(trim(adjustl(Second(1)))) ! Real because of continous solution representation
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Limit contributions of males - minimum: "//trim(Int2Char(nint(This%LimitPar1Min))) ! nint because of continous solution representation
                end if
                if (This%LimitPar1Min .lt. 1.0) then
                  write(STDERR, "(a)") " ERROR: LimitMaleContributionsMin must be at least 1"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for LimitMaleContributionsMin, for example, LimitMaleContributionsMin, 1"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("limitmalecontributionsmax")
              if (allocated(Second)) then
                This%LimitPar1Max = Char2Real(trim(adjustl(Second(1)))) ! Real because of continous solution representation
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Limit contributions of males - maximum: "//trim(Int2Char(nint(This%LimitPar1Max))) ! nint because of continous solution representation
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for LimitMaleContributionsMax, for example, LimitMaleContributionsMax, 10"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("limitmalecontributionsminweight")
              if (allocated(Second)) then
                This%LimitPar1MinWeight = Char2Real(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Limit contributions of males - weight for contributions below minimum: "//trim(Real2Char(This%LimitPar1MinWeight, fmt=FMTREAL2CHAR))
                end if
                if (This%LimitPar1MinWeight .gt. 0.0) then
                  write(STDOUT, "(a)") " NOTE: Positive weight for limit on minimum contributions, that is, encourage smaller contributions than defined minimum. Was this intended?"
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for LimitMaleContributionsMinWeight, for example, LimitMaleContributionsMinWeight, -2"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("limitfemalecontributions")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%LimitPar2 = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Limit contributions of females"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for LimitFemaleContributions, for example, LimitFemaleContributions, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("limitfemalecontributionsmin")
              if (allocated(Second)) then
                This%LimitPar2Min = Char2Real(trim(adjustl(Second(1)))) ! Real because of continous solution representation
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Limit contributions of females - minimum: "//trim(Int2Char(nint(This%LimitPar2Min))) ! nint because of continous solution representation
                end if
                if (This%LimitPar2Min .lt. 1.0) then
                  write(STDERR, "(a)") " ERROR: LimitFemaleContributionsMin must be at least 1"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for LimitFemaleContributionsMin, for example, LimitFemaleContributionsMin, 1"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("limitfemalecontributionsmax")
              if (allocated(Second)) then
                This%LimitPar2Max = Char2Real(trim(adjustl(Second(1)))) ! Real because of continous solution representation
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Limit contributions of females - maximum: "//trim(Int2Char(nint(This%LimitPar2Max))) ! nint because of continous solution representation
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for LimitFemaleContributionsMax, for example, LimitFemaleContributionsMax, 10"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("limitfemalecontributionsminweight")
              if (allocated(Second)) then
                This%LimitPar2MinWeight = Char2Real(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Limit contributions of females - weight for contributions below minimum: "//trim(Real2Char(This%LimitPar2MinWeight, fmt=FMTREAL2CHAR))
                end if
                if (This%LimitPar2MinWeight .gt. 0.0) then
                  write(STDOUT, "(a)") " NOTE: Positive weight for limit on minimum contributions, that is, encourage smaller contributions than defined minimum. Was this intended?"
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for LimitFemaleContributionsMinWeight, for example, LimitFemaleContributionsMinWeight, -2"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("mateallocation")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .ne. "yes") then
                  This%MateAllocation = .false.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Do not allocate matings/crosses"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for MateAllocation, for example, MateAllocation, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("randommateallocation")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%RandomMateAllocation = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Randomly allocate matings/crosses"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for RandomMateAllocation, for example, RandomMateAllocation, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("allowreciprocalmatings")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%ReciprocalMatingsAllowed = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Reciprocal matings/crosses allowed"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for AllowReciprocalMatings, for example, AllowReciprocalMatings, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("allowmultiplematings")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%MultipleMatingsAllowed = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Multiple (repeated) matings/crosses allowed"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for AllowMultipleMatings, for example, AllowMultipleMatings, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("multiplematingsweight")
              if (allocated(Second)) then
                This%MultipleMatingsWeight = Char2Real(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") "Multiple (repeated) matings/crosses - weight: "//trim(Real2Char(This%MultipleMatingsWeight, fmt=FMTREAL2CHAR))
                end if
                if (This%MultipleMatingsWeight .gt. 0.0) then
                  write(STDOUT, "(a)") " NOTE: Positive weight for multiple (repeated) matings/crosses, that is, encourage multiple matings/crosses. Was this intended?"
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for MultipleMatingsWeight, for example, MultipleMatingsWeight, -2"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("allowselfing")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%SelfingAllowed = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Selfing allowed"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for AllowSelfing, for example, AllowSelfing, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("selfingweight")
              if (allocated(Second)) then
                This%SelfingWeight = Char2Real(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") "Selfing - weight: "//trim(Real2Char(This%SelfingWeight, fmt=FMTREAL2CHAR))
                end if
                if (This%SelfingWeight .gt. 0.0) then
                  write(STDOUT, "(a)") " NOTE: Positive weight for selfing, that is, encourage selfing. Was this intended?"
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for SelfingWeight, for example, SelfingWeight, -2"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("page")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%PAGEPar = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Promotion of Alleles by Genome Editing (PAGE)"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for PAGE, for example, PAGE, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("pagemax")
              if (allocated(Second)) then
                This%PAGEParMax = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Promotion of Alleles by Genome Editing (PAGE) - maximum number of individuals: "//trim(Int2Char(This%PAGEParMax))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify number for PAGEMax, for example, PAGEMax, 10"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("pagemales")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%PAGEPar  = .true.
                  This%PAGEPar1 = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Promotion of Alleles by Genome Editing (PAGE) in males"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for PAGEMales, for example, PAGEMales, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("pagemalesmax")
              if (allocated(Second)) then
                This%PAGEPar1Max = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Promotion of Alleles by Genome Editing (PAGE) in males - maxium number of individuals: "//trim(Int2Char(This%PAGEPar1Max))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify number for PAGEMalesMax, for example, PAGEMalesMax, 10"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("pagefemales")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%PAGEPar  = .true.
                  This%PAGEPar2 = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Promotion of Alleles by Genome Editing (PAGE) in females"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for PAGEFemales, for example, PAGEFemales, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("pagefemalesmax")
              if (allocated(Second)) then
                This%PAGEPar2Max = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Promotion of Alleles by Genome Editing (PAGE) in females - maximum number of individuals: "//trim(Int2Char(This%PAGEPar2Max))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify number for PAGEFemalesMax, for example, PAGEFemalesMax, 10"
                write(STDERR, "(a)") " "
                stop 1
              end if

            ! Algorithm specifications
            case ("evolalgnumberofsolutions")
              if (allocated(Second)) then
                This%EvolAlgNSol = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Evolutionary algorithm - number of solutions: "//trim(Int2Char(This%EvolAlgNSol))
                end if
                if (This%EvolAlgNSol .lt. 3) then
                  write(STDERR, "(a)") " ERROR: EvolAlgNumberOfSolutions must be at least 3 (preferably more)"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for EvolAlgNumberOfSolutions, for example, EvolAlgNumberOfSolutions, 100"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("evolalgnumberofiterations")
              if (allocated(Second)) then
                This%EvolAlgNIter = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Evolutionary algorithm - number of iterations: "//trim(Int2Char(This%EvolAlgNIter))
                end if
                if (This%EvolAlgNIter .lt. 1) then
                  write(STDERR, "(a)") " ERROR: EvolAlgNumberOfIterations must be at least 1 (preferably more)"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for EvolAlgNumberOfIterations, for example, EvolAlgNumberOfIterations, 10000"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("evolalgnumberofiterationsprint")
              if (allocated(Second)) then
                This%EvolAlgNIterPrint = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Evolutionary algorithm - number of iterations to print optimisation status: "//trim(Int2Char(This%EvolAlgNIterPrint))
                end if
                if (This%EvolAlgNIterPrint .lt. 1) then
                  write(STDERR, "(a)") " ERROR: EvolAlgNumberOfIterationsPrint must be at least 1 (preferably more)"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for EvolAlgNumberOfIterationsPrint, for example, EvolAlgNumberOfIterationsPrint, 100"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("evolalgnumberofiterationsstop")
              if (allocated(Second)) then
                This%EvolAlgNIterStop = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Evolutionary algorithm - number of iterations to stop upon no improvement of objective: "//trim(Int2Char(This%EvolAlgNIterStop))
                end if
                if (This%EvolAlgNIterStop .lt. 1) then
                  write(STDERR, "(a)") " ERROR: EvolAlgNumberOfIterationsStop must be at least 1 (preferably more)"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for EvolAlgNumberOfIterationsStop, for example, EvolAlgNumberOfIterationsStop, 100"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("evolalgstoptolerancecoancestry")
              if (allocated(Second)) then
                This%EvolAlgStopTolCoancestry = Char2Real(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Evolutionary algorithm - stopping tolerance for minimum coancestry or minimum inbreeding optimisations: "//trim(Real2Char(This%EvolAlgStopTolCoancestry, fmt=FMTREAL2CHAR))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for EvolAlgStopToleranceCoancestry, for example, EvolAlgStopToleranceCoancestry, 0.01"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("evolalgstoptolerance")
              if (allocated(Second)) then
                This%EvolAlgStopTol = Char2Real(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Evolutionary algorithm - stopping tolerance for maximum criterion or optimum/balance optimisations: "//trim(Real2Char(This%EvolAlgStopTol, fmt=FMTREAL2CHAR))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for EvolAlgStopTol, for example, EvolAlgStopTol, 0.01"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("evolalglogallsolutions")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%EvolAlgLogPop = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Evolutionary algorithm - log all evaluated solutions"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for EvolAlgLogAllSolutions, for example, EvolAlgLogAllSolutions, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("evolalg")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .ne. "none") then
                  write(This%EvolAlg, *) trim(adjustl(Second(1)))
                  This%EvolAlg = adjustl(This%EvolAlg)
                  if (.not. (This%EvolAlg .eq. "DE")) then
                    write(STDERR, "(a)") " ERROR: Must specify a valid algorithm [DE/???] for EvolAlg, for example, EvolAlg, DE"
                    write(STDERR, "(a)") " "
                    stop 1
                  end if
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Evolutionary algorithm: "//trim(This%EvolAlg)
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a algorithm for EvolAlg, for example, EvolAlg, DE"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("diffevolnumberofiterationsburnin")
              if (allocated(Second)) then
                This%DiffEvolNIterBurnIn = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Differential evolution algorithm - number of warming iterations (burn-in): "//trim(Int2Char(This%DiffEvolNIterBurnIn))
                end if
                if (This%DiffEvolNIterBurnIn .lt. 1) then
                  write(STDERR, "(a)") " ERROR: EvolAlgNumberOfIterationsBurnin must be at least 1 (preferably more)"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for EvolAlgNumberOfIterationsBurnin, for example, EvolAlgNumberOfIterationsBurnin, 1000"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("diffevolparametercrburnin")
              if (allocated(Second)) then
                This%DiffEvolParamCrBurnIn = Char2Real(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Differential evolution algorithm - cross-over parameter for warmup (burn-in): "//trim(Real2Char(This%DiffEvolParamCrBurnIn, fmt=FMTREAL2CHAR))
                end if
                if (This%DiffEvolParamCrBurnIn .lt. 0) then
                  write(STDERR, "(a)") " ERROR: DiffEvolParameterCrBurnin must be at least 0 (preferably more)"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for DiffEvolParameterCrBurnin, for example, DiffEvolParameterCrBurnin, 0.4"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("diffevolparametercr1")
              if (allocated(Second)) then
                This%DiffEvolParamCr1 = Char2Real(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Differential evolution algorithm - cross-over parameter 1 (common): "//trim(Real2Char(This%DiffEvolParamCr1, fmt=FMTREAL2CHAR))
                end if
                if (This%DiffEvolParamCr1 .lt. 0) then
                  write(STDERR, "(a)") " ERROR: DiffEvolParameterCr1 must be at least 0 (preferably more)"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for DiffEvolParameterCr1, for example, DiffEvolParameterCr1, 0.2"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("diffevolparametercr2")
              if (allocated(Second)) then
                This%DiffEvolParamCr2 = Char2Real(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Differential evolution algorithm - cross-over parameter 2 (rare): "//trim(Real2Char(This%DiffEvolParamCr2, fmt=FMTREAL2CHAR))
                end if
                if (This%DiffEvolParamCr2 .lt. 0) then
                  write(STDERR, "(a)") " ERROR: DiffEvolParameterCr2 must be at least 0 (preferably more)"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for DiffEvolParameterCr2, for example, DiffEvolParameterCr2, 0.2"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("diffevolparameterfbase")
              if (allocated(Second)) then
                This%DiffEvolParamFBase = Char2Real(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Differential evolution algorithm - parameter F (base value): "//trim(Real2Char(This%DiffEvolParamFBase, fmt=FMTREAL2CHAR))
                end if
                if (This%DiffEvolParamFBase .lt. 0) then
                  write(STDERR, "(a)") " ERROR: DiffEvolParameterFBase must be at least 0 (preferably more)"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for DiffEvolParameterFBase, for example, DiffEvolParameterFBase, 0.1"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("diffevolparameterfhigh1")
              if (allocated(Second)) then
                This%DiffEvolParamFHigh1 = Char2Real(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Differential evolution algorithm - parameter F (high value 1): "//trim(Real2Char(This%DiffEvolParamFHigh1, fmt=FMTREAL2CHAR))
                end if
                if (This%DiffEvolParamFHigh1 .lt. 0) then
                  write(STDERR, "(a)") " ERROR: DiffEvolParameterFHigh1 must be at least 0 (preferably more)"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for DiffEvolParameterFHigh1, for example, DiffEvolParameterFHigh1, 1.0"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("diffevolparameterfhigh2")
              if (allocated(Second)) then
                This%DiffEvolParamFHigh2 = Char2Real(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Differential evolution algorithm - parameter F (high value 2): "//trim(Real2Char(This%DiffEvolParamFHigh2, fmt=FMTREAL2CHAR))
                end if
                if (This%DiffEvolParamFHigh2 .lt. 0) then
                  write(STDERR, "(a)") " ERROR: DiffEvolParameterFHigh2 must be at least 0 (preferably more)"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for DiffEvolParameterFHigh2, for example, DiffEvolParameterFHigh2, 4.0"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("randomsearchstricter")
              if (allocated(Second)) then
                This%RanAlgStricter = Char2Real(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Random search algorithm - perform k times more iterations than with the evolutionary algorithm: k="//trim(Int2Char(This%RanAlgStricter))
                end if
                if (This%RanAlgStricter .lt. 1) then
                  write(STDERR, "(a)") " ERROR: RandomSearchStricter must be at least 1"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for RandomSearchStricter, for example, RandomSearchStricter, 10"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("stop")
              if (LogStdoutInternal) then
                write(STDOUT, "(3a)") " NOTE: Encountered Stop specification - the rest of specifications will be ignored"
              end if
              exit

            case default
              if (LogStdoutInternal) then
                write(STDOUT, "(a)") " NOTE: Specification '"//trim(Line)//"' was ignored"
                write(STDOUT, "(a)") " "
              end if
          end select
        end if
      end do ReadSpec
      close(SpecUnit)

      if ((This%ModeMinCoancestry .and. This%ModeMaxCriterion) .or. &
          (This%ModeMinInbreeding .and. This%ModeMaxCriterion) .or. &
          (This%ModeMinCoancestry .and. This%ModeMinInbreeding)) then
        This%ModeOpt = .true.
      end if

      if (.not. (This%ModeMinCoancestry .or. &
                 This%ModeMinInbreeding .or. &
                 This%ModeMaxCriterion  .or. &
                 This%ModeRan           .or. &
                 This%ModeOpt)) then
        write(STDERR, "(a)") " ERROR: One of the modes must be activated!"
        write(STDERR, "(a)") " ERROR: ModeMinCoancestry, ModeMinInbreeding, ModeMaxCriterion, or ModeOpt"
        write(STDERR, "(a)") " "
        stop 1
      end if

      if (This%ModeOpt .and. .not. (This%TargetDegreeGiven       .or. &
                                    This%TargetSelCriterionGiven .or. This%TargetSelCriterionStdGiven .or. This%TargetMaxCriterionPctGiven  .or. &
                                    This%TargetCoancestryGiven   .or. This%TargetCoancestryRateGiven  .or. This%TargetMinCoancestryPctGiven .or. &
                                    This%TargetInbreedingGiven   .or. This%TargetInbreedingRateGiven  .or. This%TargetMinInbreedingPctGiven .or. &
                                    This%EvaluateFrontier)) then
        write(STDERR, "(a)") " ERROR: One of targets must be provided when ModeOpt is activated!"
        write(STDERR, "(a)") " ERROR: TargetDegree,"
        write(STDERR, "(a)") " ERROR: TargetSelCriterion, TargetSelCriterionStd,    TargetMaxCriterionPct,"
        write(STDERR, "(a)") " ERROR: TargetCoancestry,   TargetCoancestryRate,     TargetMinCoancestryPct"
        write(STDERR, "(a)") " ERROR: TargetInbreeding,   TargetInbreedingRate,  or TargetMinInbreedingPct"
        write(STDERR, "(a)") " "
        stop 1
      end if

      if (.not. This%RelMtxGiven) then
        write(STDERR, "(a)") " ERROR: One of CoancestryMatrixFile or NrmMatrixFile must be specified!"
        write(STDERR, "(a)") " "
        stop 1
      end if

      if (.not. This%SelCriterionGiven .and. &
          (This%ModeMaxCriterion .or. This%ModeOpt)) then
        write(STDERR, "(a)") " ERROR: Selection criterion is needed for modes: ModeMaxCriterion or ModeOpt!"
        write(STDERR, "(a)") " "
        stop 1
      end if

      if (This%SeedFileGiven .and. This%SeedGiven) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " NOTE: The specification Seed has priority over SeedFile."
          write(STDOUT, "(a)") " "
        end if
        This%SeedFile = ""
        This%SeedFileGiven = .false.
      end if

      if (This%nMat .le. 0) then
        write(STDERR, "(a)") " ERROR: Number of matings/crosses must be larger than zero!"
        write(STDERR, "(a)") " "
        stop 1
      end if

      ! The nPar tests are in ReadAlphaMateData where we count number of individuals and males and females

      if (.not. This%GenderGiven .and. (This%LimitParMin .eq. This%LimitParMax)) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " NOTE: Since LimitContributionsMin equals LimitContributionsMax, option EqualizeContributions is activated."
          write(STDOUT, "(a)") " "
        end if
        This%EqualizePar = .true.
      end if

      if (This%GenderGiven .and. (This%LimitPar1Min .eq. This%LimitPar1Max)) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " NOTE: Since LimitMaleContributionsMin equals LimitMaleContributionsMax, option EqualizeMaleContributions is activated."
          write(STDOUT, "(a)") " "
        end if
        This%EqualizePar1 = .true.
      end if

      if (This%GenderGiven .and. (This%LimitPar2Min .eq. This%LimitPar2Max)) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " NOTE: Since LimitFemaleContributionsMin equals LimitFemaleContributionsMax, option EqualizeFemaleContributions is activated."
          write(STDOUT, "(a)") " "
        end if
        This%EqualizePar2 = .true.
      end if

      if (.not. This%GenderGiven .and. (This%LimitPar .and. This%EqualizePar)) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " NOTE: The specification EqualizeContributions has priority over LimitContributions."
          write(STDOUT, "(a)") " "
        end if
        ! ... therefore reset all limit specifications to default values
        This%LimitPar          = .false.
        This%LimitParMin       = 1.0
        This%LimitParMax       = huge(This%LimitParMax) - 1.0
        This%LimitParMinWeight = -1
      end if

      if (This%GenderGiven .and. (This%LimitPar1 .and. This%EqualizePar1)) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " NOTE: The specification EqualizeMaleContributions has priority over LimitMaleContributions."
          write(STDOUT, "(a)") " "
        end if
        ! ... therefore reset all limit specifications to default values
        This%LimitPar1          = .false.
        This%LimitPar1Min       = 1.0
        This%LimitPar1Max       = huge(This%LimitPar1Max) - 1.0
        This%LimitPar1MinWeight = -1.0
      end if

      if (This%GenderGiven .and. (This%LimitPar2 .and. This%EqualizePar2)) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " NOTE: The specification EqualizeFemaleContributions has priority over LimitFemaleContributions."
          write(STDOUT, "(a)") " "
        end if
        ! ... therefore reset all limit specifications to default values
        This%LimitPar2          = .false.
        This%LimitPar2Min       = 1.0
        This%LimitPar2Max       = huge(This%LimitPar2Max) - 1.0
        This%LimitPar2MinWeight = -1.0
      end if

      ! Impose upper limits (when it is not given by the user) to avoid explosion in optimisation
      if (.not. This%GenderGiven .and. (.not. This%LimitPar .and. .not. This%EqualizePar)) then
        This%LimitPar    = .true.
        This%LimitParMax = FLOATFUN(This%nMat) * 2.0
      end if

      if (This%GenderGiven .and. (.not. This%LimitPar1 .and. .not. This%EqualizePar1)) then
        This%LimitPar1    = .true.
        This%LimitPar1Max = FLOATFUN(This%nMat)
      end if

      if (This%GenderGiven .and. (.not. This%LimitPar2 .and. .not. This%EqualizePar2)) then
        This%LimitPar2    = .true.
        This%LimitPar2Max = FLOATFUN(This%nMat)
      end if

      if (.not. This%GenderGiven) then
        This%EqualizePar1       = This%EqualizePar

        This%LimitPar1          = This%LimitPar
        This%LimitPar1Min       = This%LimitParMin
        This%LimitPar1Max       = This%LimitParMax
        This%LimitPar1MinWeight = This%LimitParMinWeight

        This%PreselectPar1      = This%PreselectPar
        This%PreselectPar1Pct   = This%PreselectParPct
        This%PreselectPar1N     = This%PreselectParN

        This%PAGEPar1           = This%PAGEPar
        This%PAGEPar1Max        = This%PAGEParMax
      end if

      if ((.not. This%SelCriterionGiven) .and. This%PAGEPar) then
        write(STDERR, "(a)") " ERROR: Can not use PAGE when selection criterion file is not given!"
        write(STDERR, "(a)") " "
        stop 1
      end if

      if (This%GenderGiven .and. This%SelfingAllowed) then
        write(STDERR, "(a)") " ERROR: When gender matters, AlphaMate can not perform selfing!"
        write(STDERR, "(a)") " NOTE: If you really want to do this, you could present individuals in the data once as males and once as females"
        write(STDERR, "(a)") " "
        stop 1
      end if

      if ((.not. This%MateAllocation .or. This%RandomMateAllocation) .and. &
          (This%TargetInbreedingGiven .or. This%TargetInbreedingRateGiven .or. This%TargetMinInbreedingPctGiven)) then
        This%ModeMinInbreeding           = .false.
        This%TargetInbreedingGiven       = .false.
        This%TargetInbreedingRateGiven   = .false.
        This%TargetMinInbreedingPctGiven = .false.
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " NOTE: Inbreeding target is not active when mate allocation is inactive or random."
          write(STDOUT, "(a)") " "
        end if
      end if

      call This%SetupColNamesAndFormats
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Setup AlphaMate optimisation mode
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   April 16, 2017
    !---------------------------------------------------------------------------
    subroutine SetupModeAlphaMateSpec(This, Mode, Data, Degree, &
                                      SelCriterion, SelCriterionStd, MaxCriterionPct,                         ModeMaxCriterionSpec,&
                                      Coancestry,   CoancestryRate,  MinCoancestryPct, CoancestryWeightBelow, ModeMinCoancestrySpec, &
                                      Inbreeding,   InbreedingRate,  MinInbreedingPct, InbreedingWeightBelow, ModeMinInbreedingSpec)
      ! @todo not pure due to error stop in here and in SetTargets()
      ! @todo This will work with Fortran 2015 standard (or at least with new? ifort) that allows error stop in PURE subroutines
      implicit none
      class(AlphaMateSpec), intent(inout)           :: This                  !< @return AlphaMateSpec holder
      character(len=*), intent(in)                  :: Mode                  !< Mode definition/name
      type(AlphaMateData), intent(in), optional     :: Data                  !< AlphaMateData holder
      real(FLOATTYPE), intent(in), optional         :: Degree                !< Targeted degree
      real(FLOATTYPE), intent(in), optional         :: SelCriterion          !< Targeted selection criterion
      real(FLOATTYPE), intent(in), optional         :: SelCriterionStd       !< Targeted stand. selection criterion
      real(FLOATTYPE), intent(in), optional         :: MaxCriterionPct       !< Targeted maximum criterion percentage
      real(FLOATTYPE), intent(in), optional         :: Coancestry            !< Targeted coancestry
      real(FLOATTYPE), intent(in), optional         :: CoancestryRate        !< Targeted coancestry rate
      real(FLOATTYPE), intent(in), optional         :: MinCoancestryPct      !< Targeted minimum coancestry percentage
      type(AlphaMateModeSpec), intent(in), optional :: ModeMaxCriterionSpec  !< Maximum criterion  solution specs
      logical, intent(in), optional                 :: CoancestryWeightBelow !< Weight deviations below the targeted coancestry
      type(AlphaMateModeSpec), intent(in), optional :: ModeMinCoancestrySpec !< Minimum coancestry solution specs
      real(FLOATTYPE), intent(in), optional         :: Inbreeding            !< Targeted inbreeding
      real(FLOATTYPE), intent(in), optional         :: InbreedingRate        !< Targeted inbreeding rate
      real(FLOATTYPE), intent(in), optional         :: MinInbreedingPct      !< Targeted minimum inbreeding percentage
      logical, intent(in), optional                 :: InbreedingWeightBelow !< Weight deviations below the targeted inbreeding
      type(AlphaMateModeSpec), intent(in), optional :: ModeMinInbreedingSpec !< Minimum inbreeding solution specs

      select case (trim(Mode))
        case ("MinCoancestry") ! Only coancestry!!!
          call This%ModeMinCoancestrySpec%Initialise(Name="MinCoancestry")
          This%ModeMinCoancestrySpec%ObjectiveCoancestry = .true.
          ! @todo Do these lines still make sense when we go above two objectives?
          This%ModeMinCoancestrySpec%TargetDegree           =  90.0
          This%ModeMinCoancestrySpec%TargetMinCoancestryPct = 100.0
          This%ModeMinCoancestrySpec%TargetMinInbreedingPct =   0.0
          This%ModeMinCoancestrySpec%TargetMaxCriterionPct  =   0.0
          call This%ModeSpec%Assign(In=This%ModeMinCoancestrySpec)

        case ("MinInbreeding") ! Only inbreeding!!!
          call This%ModeMinInbreedingSpec%Initialise(Name="MinInbreeding")
          This%ModeMinInbreedingSpec%ObjectiveInbreeding = .true.
          ! @todo Do these lines still make sense when we go above two objectives?
          This%ModeMinInbreedingSpec%TargetDegree           =  45.0
          This%ModeMinInbreedingSpec%TargetMinCoancestryPct =   0.0
          This%ModeMinInbreedingSpec%TargetMinInbreedingPct = 100.0
          This%ModeMinInbreedingSpec%TargetMaxCriterionPct  =   0.0
          call This%ModeSpec%Assign(In=This%ModeMinInbreedingSpec)

        case ("MaxCriterion") ! Only criterion!!!
          call This%ModeMaxCriterionSpec%Initialise(Name="MaxCriterion")
          This%ModeMaxCriterionSpec%ObjectiveCriterion = .true.
          ! @todo Do these lines still make sense when we go above two objectives?
          This%ModeMaxCriterionSpec%TargetDegree           =   0.0
          This%ModeMaxCriterionSpec%TargetMinCoancestryPct =   0.0
          This%ModeMaxCriterionSpec%TargetMinInbreedingPct =   0.0
          This%ModeMaxCriterionSpec%TargetMaxCriterionPct  = 100.0
          call This%ModeSpec%Assign(In=This%ModeMaxCriterionSpec)

        case ("Opt") ! All objectives jointly!!!
          call This%ModeOptSpec%Initialise(Name="Opt")
          This%ModeOptSpec%ObjectiveCriterion  = .true.
          This%ModeOptSpec%ObjectiveCoancestry = .true.

          if      (present(Degree)) then
            if (Degree .gt. 90 .or. Degree .lt. 0) then
              write(STDERR, "(a)") " ERROR: TargetDegree must be between 0 and 90!"
              write(STDERR, "(a)") " ERROR: Target: "//trim(Real2Char(Degree, fmt=FMTREAL2CHAR))
              write(STDERR, "(a)") " "
              stop 1
            end if
            call This%ModeOptSpec%SetTargets(Degree=Degree, &
                                             Data=Data, ModeMinCoancestrySpec=ModeMinCoancestrySpec, ModeMaxCriterionSpec=ModeMaxCriterionSpec)
          else if (present(SelCriterion)) then
            if (SelCriterion .gt. ModeMaxCriterionSpec%SelCriterion) then
              write(STDERR, "(a)") " ERROR: TargetSelCriterion must be less than the the maximum achieved SelCriterion!"
              write(STDERR, "(a)") " ERROR: Target:           "//trim(Real2Char(SelCriterion,                      fmt=FMTREAL2CHAR))
              write(STDERR, "(a)") " ERROR: Maximum achieved: "//trim(Real2Char(ModeMaxCriterionSpec%SelCriterion, fmt=FMTREAL2CHAR))
              write(STDERR, "(a)") " "
              stop 1
            end if
            call This%ModeOptSpec%SetTargets(SelCriterion=SelCriterion, &
                                             Data=Data, ModeMinCoancestrySpec=ModeMinCoancestrySpec, ModeMaxCriterionSpec=ModeMaxCriterionSpec)
          else if (present(SelCriterionStd)) then
            if (SelCriterionStd .gt. ModeMaxCriterionSpec%SelCriterionStd) then
              write(STDERR, "(a)") " ERROR: TargetSelCriterionStd must be less than the maximum achieved SelCriterionStd!"
              write(STDERR, "(a)") " ERROR: Target:           "//trim(Real2Char(SelCriterionStd,                      fmt=FMTREAL2CHAR))
              write(STDERR, "(a)") " ERROR: Maximum achieved: "//trim(Real2Char(ModeMaxCriterionSpec%SelCriterionStd, fmt=FMTREAL2CHAR))
              write(STDERR, "(a)") " "
              stop 1
            end if
            call This%ModeOptSpec%SetTargets(SelCriterionStd=SelCriterionStd, &
                                             Data=Data, ModeMinCoancestrySpec=ModeMinCoancestrySpec, ModeMaxCriterionSpec=ModeMaxCriterionSpec)
          else if (present(MaxCriterionPct)) then
            if (MaxCriterionPct .gt. 100 .or. MaxCriterionPct .lt. 0) then
              write(STDERR, "(a)") " ERROR: TargetMaxCriterionPct must be between 0 and 100!"
              write(STDERR, "(a)") " ERROR: Target: "//trim(Real2Char(MaxCriterionPct, fmt=FMTREAL2CHAR))
              write(STDERR, "(a)") " "
              stop 1
            end if
            call This%ModeOptSpec%SetTargets(MaxCriterionPct=MaxCriterionPct, &
                                             Data=Data, ModeMinCoancestrySpec=ModeMinCoancestrySpec, ModeMaxCriterionSpec=ModeMaxCriterionSpec)
          else if (present(Coancestry)) then
            if (Coancestry .lt. ModeMinCoancestrySpec%Coancestry) then
              write(STDERR, "(a)") " ERROR: TargetCoancestry must be more than the minimum achieved Coancestry!"
              write(STDERR, "(a)") " ERROR: Target:           "//trim(Real2Char(Coancestry,                      fmt=FMTREAL2CHAR))
              write(STDERR, "(a)") " ERROR: Minimum achieved: "//trim(Real2Char(ModeMinCoancestrySpec%Coancestry, fmt=FMTREAL2CHAR))
              write(STDERR, "(a)") " "
              stop 1
            end if
            call This%ModeOptSpec%SetTargets(Coancestry=Coancestry, &
                                             Data=Data, ModeMinCoancestrySpec=ModeMinCoancestrySpec, ModeMaxCriterionSpec=ModeMaxCriterionSpec)
          else if (present(CoancestryRate)) then
            if (CoancestryRate .lt. ModeMinCoancestrySpec%CoancestryRate) then
              write(STDERR, "(a)") " ERROR: TargetCoancestryRate must be more than the minimum achieved CoancestryRate!"
              write(STDERR, "(a)") " ERROR: Target:           "//trim(Real2Char(CoancestryRate,                      fmt=FMTREAL2CHAR))
              write(STDERR, "(a)") " ERROR: Minimum achieved: "//trim(Real2Char(ModeMinCoancestrySpec%CoancestryRate, fmt=FMTREAL2CHAR))
              write(STDERR, "(a)") " "
              stop 1
            end if
            call This%ModeOptSpec%SetTargets(CoancestryRate=CoancestryRate, &
                                             Data=Data, ModeMinCoancestrySpec=ModeMinCoancestrySpec, ModeMaxCriterionSpec=ModeMaxCriterionSpec)
          else if (present(MinCoancestryPct)) then
            if (MinCoancestryPct .gt. 100 .or. MaxCriterionPct .lt. 0) then
              write(STDERR, "(a)") " ERROR: TargetMinCoancestryPct must be between 0 and 100!"
              write(STDERR, "(a)") " ERROR: Target: "//trim(Real2Char(MinCoancestryPct, fmt=FMTREAL2CHAR))
              write(STDERR, "(a)") " "
              stop 1
            end if
            call This%ModeOptSpec%SetTargets(MinCoancestryPct=MinCoancestryPct, &
                                             Data=Data, ModeMinCoancestrySpec=ModeMinCoancestrySpec, ModeMaxCriterionSpec=ModeMaxCriterionSpec)
          end if

          if (present(CoancestryWeightBelow)) then
            This%ModeOptSpec%CoancestryWeightBelow =      CoancestryWeightBelow
          else
            This%ModeOptSpec%CoancestryWeightBelow = This%CoancestryWeightBelow
          end if

          if      (This%TargetInbreedingGiven) then
            if (This%TargetInbreeding .lt. ModeMinInbreedingSpec%Inbreeding) then
              write(STDERR, "(a)") " ERROR: TargetInbreeding must be more than the minimum achieved Inbreeding!"
              write(STDERR, "(a)") " ERROR: Target:           "//trim(Real2Char(This%TargetInbreeding,            fmt=FMTREAL2CHAR))
              write(STDERR, "(a)") " ERROR: Minimum achieved: "//trim(Real2Char(ModeMinInbreedingSpec%Inbreeding, fmt=FMTREAL2CHAR))
              write(STDERR, "(a)") " "
              stop 1
            end if
            This%ModeOptSpec%ObjectiveInbreeding = .true.
            call This%ModeOptSpec%SetTargets(Inbreeding=This%TargetInbreeding, &
                                             Data=Data, ModeMinInbreedingSpec=ModeMinInbreedingSpec)
          else if (This%TargetInbreedingRateGiven) then
            if (This%TargetInbreedingRate .lt. ModeMinInbreedingSpec%InbreedingRate) then
              write(STDERR, "(a)") " ERROR: TargetInbreedingRate must be more than the minimum achieved InbreedingRate!"
              write(STDERR, "(a)") " ERROR: Target:           "//trim(Real2Char(This%TargetInbreedingRate,            fmt=FMTREAL2CHAR))
              write(STDERR, "(a)") " ERROR: Minimum achieved: "//trim(Real2Char(ModeMinInbreedingSpec%InbreedingRate, fmt=FMTREAL2CHAR))
              write(STDERR, "(a)") " "
              stop 1
            end if
            This%ModeOptSpec%ObjectiveInbreeding = .true.
            call This%ModeOptSpec%SetTargets(InbreedingRate=This%TargetInbreedingRate, &
                                             Data=Data, ModeMinInbreedingSpec=ModeMinInbreedingSpec)
          else if (This%TargetMinInbreedingPctGiven) then
            if (This%TargetMinInbreedingPct .gt. 100 .or. This%TargetMinInbreedingPct .lt. 0) then
              write(STDERR, "(a)") " ERROR: TargetMinInbreedingPct must be between 0 and 100!"
              write(STDERR, "(a)") " ERROR: Target: "//trim(Real2Char(This%TargetMinInbreedingPct, fmt=FMTREAL2CHAR))
              write(STDERR, "(a)") " "
              stop 1
            end if
            This%ModeOptSpec%ObjectiveInbreeding = .true.
            call This%ModeOptSpec%SetTargets(MinInbreedingPct=This%TargetMinInbreedingPct, &
                                             Data=Data, ModeMinInbreedingSpec=ModeMinInbreedingSpec, ModeMaxCriterionSpec=ModeMaxCriterionSpec)
          end if

          if (present(InbreedingWeightBelow)) then
            This%ModeOptSpec%InbreedingWeightBelow =      InbreedingWeightBelow
          else
            This%ModeOptSpec%InbreedingWeightBelow = This%InbreedingWeightBelow
          end if

          call This%ModeSpec%Assign(In=This%ModeOptSpec)

        case ("Ran")
          ! @todo do we need ModeRanSpec?
          call This%ModeRanSpec%Initialise(Name="Ran")
          ! @todo???
          This%ModeRanSpec%CoancestryWeightBelow = This%CoancestryWeightBelow
          ! @todo???
          This%ModeRanSpec%InbreedingWeightBelow = This%InbreedingWeightBelow
          call This%ModeSpec%Assign(In=This%ModeRanSpec)

      end select

    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Initialise AlphaMate optimisation mode specifications
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   April 16, 2017
    !---------------------------------------------------------------------------
    pure subroutine InitialiseAlphaMateModeSpec(This, Name)
      implicit none
      class(AlphaMateModeSpec), intent(out) :: This !< @return AlphaMateModeSpec holder
      character(len=*), intent(in)          :: Name !< Mode name
      real(FLOATTYPE) :: NANFLOATTYPE
      NANFLOATTYPE = IEEE_Value(x=NANFLOATTYPE, class=IEEE_Quiet_NaN)
      This%Name = Name
      This%ObjectiveCriterion = .false.
      This%ObjectiveCoancestry = .false.
      This%ObjectiveInbreeding = .false.
      This%TargetDegree = NANFLOATTYPE
      This%TargetSelCriterion = NANFLOATTYPE
      This%TargetSelCriterionStd = NANFLOATTYPE
      This%TargetMaxCriterionPct = NANFLOATTYPE
      This%TargetCoancestry = NANFLOATTYPE
      This%TargetCoancestryRate = NANFLOATTYPE
      This%CoancestryWeightBelow = .false.
      This%TargetMinCoancestryPct = NANFLOATTYPE
      This%TargetInbreeding = NANFLOATTYPE
      This%TargetInbreedingRate = NANFLOATTYPE
      This%InbreedingWeightBelow = .false.
      This%TargetMinInbreedingPct = NANFLOATTYPE
      This%Degree = NANFLOATTYPE
      This%SelCriterion = NANFLOATTYPE
      This%SelCriterionStd = NANFLOATTYPE
      This%MaxCriterionPct = NANFLOATTYPE
      This%Coancestry = NANFLOATTYPE
      This%CoancestryRate = NANFLOATTYPE
      This%MinCoancestryPct = NANFLOATTYPE
      This%Inbreeding = NANFLOATTYPE
      This%InbreedingRate = NANFLOATTYPE
      This%MinInbreedingPct = NANFLOATTYPE
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Assign AlphaMate optimisation mode specifications
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   April 16, 2017
    !---------------------------------------------------------------------------
    pure subroutine AssignAlphaMateModeSpec(Out, In)
      implicit none
      class(AlphaMateModeSpec), intent(out) :: Out !< @return AlphaMateModeSpec holder
      class(AlphaMateModeSpec), intent(in)  :: In  !< AlphaMateModeSpec holder
      Out%Name = In%Name
      Out%ObjectiveCriterion = In%ObjectiveCriterion
      Out%ObjectiveCoancestry = In%ObjectiveCoancestry
      Out%ObjectiveInbreeding = In%ObjectiveInbreeding
      Out%TargetDegree = In%TargetDegree
      Out%TargetSelCriterion = In%TargetSelCriterion
      Out%TargetSelCriterionStd = In%TargetSelCriterionStd
      Out%TargetMaxCriterionPct = In%TargetMaxCriterionPct
      Out%TargetCoancestry = In%TargetCoancestry
      Out%TargetCoancestryRate = In%TargetCoancestryRate
      Out%CoancestryWeightBelow = In%CoancestryWeightBelow
      Out%TargetMinCoancestryPct = In%TargetMinCoancestryPct
      Out%TargetInbreeding = In%TargetInbreeding
      Out%TargetInbreedingRate = In%TargetInbreedingRate
      Out%InbreedingWeightBelow = In%InbreedingWeightBelow
      Out%TargetMinInbreedingPct = In%TargetMinInbreedingPct
      Out%Degree = In%Degree
      Out%SelCriterion = In%SelCriterion
      Out%SelCriterionStd = In%SelCriterionStd
      Out%MaxCriterionPct = In%MaxCriterionPct
      Out%Coancestry = In%Coancestry
      Out%CoancestryRate = In%CoancestryRate
      Out%MinCoancestryPct = In%MinCoancestryPct
      Out%Inbreeding = In%Inbreeding
      Out%InbreedingRate = In%InbreedingRate
      Out%MinInbreedingPct = In%MinInbreedingPct
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Save AlphaMate optimisation mode specifications
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   April 16, 2017
    !---------------------------------------------------------------------------
    pure subroutine SaveSol2ModeSpecAlphaMateModeSpec(Out, In)
      implicit none
      class(AlphaMateModeSpec), intent(out) :: Out !< @return AlphaMateModeSpec holder
      class(AlphaMateSol), intent(in)       :: In  !< AlphaMateSol holder
      Out%Degree = In%Degree
      Out%SelCriterion = In%SelCriterion
      Out%SelCriterionStd = In%SelCriterionStd
      Out%MaxCriterionPct = In%MaxCriterionPct
      Out%Coancestry = In%CoancestryRanMate
      Out%CoancestryRate = In%CoancestryRateRanMate
      Out%MinCoancestryPct = In%MinCoancestryPct
      Out%Inbreeding = In%Inbreeding
      Out%InbreedingRate = In%InbreedingRate
      Out%MinInbreedingPct = In%MinInbreedingPct
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Set targets for an AlphaMate optimisation mode
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   April 13, 2017
    !---------------------------------------------------------------------------
    subroutine SetTargetsAlphaMateModeSpec(This, Data, Degree, &
                                           SelCriterion, SelCriterionStd, MaxCriterionPct, &
                                           Coancestry, CoancestryRate, MinCoancestryPct, &
                                           Inbreeding, InbreedingRate, MinInbreedingPct, &
                                           ModeMinCoancestrySpec, ModeMinInbreedingSpec, ModeMaxCriterionSpec)
      ! @todo not pure due to error stop
      ! @todo This will work with Fortran 2015 standard (or at least with new? ifort) that allows error stop in PURE subroutines
      implicit none
      class(AlphaMateModeSpec), intent(inout)       :: This                  !< @return AlphaMateModeSpec holder
      type(AlphaMateData), intent(in), optional     :: Data                  !< AlphaMateData holder
      real(FLOATTYPE), intent(in), optional         :: Degree                !< Targeted selection/Coancestry frontier degree
      real(FLOATTYPE), intent(in), optional         :: SelCriterion          !< Targeted selection criterion
      real(FLOATTYPE), intent(in), optional         :: SelCriterionStd       !< Targeted stand. selection criterion
      real(FLOATTYPE), intent(in), optional         :: MaxCriterionPct       !< Targeted maximum criterion percentage (100 means the maximum possible selection criterion)
      real(FLOATTYPE), intent(in), optional         :: Coancestry            !< Targeted coancestry
      real(FLOATTYPE), intent(in), optional         :: CoancestryRate        !< Targeted coancestry rate
      real(FLOATTYPE), intent(in), optional         :: MinCoancestryPct      !< Targeted minimum coancestry percentage (100 means the minimum possible coancestry)
      real(FLOATTYPE), intent(in), optional         :: Inbreeding            !< Targeted inbreeding
      real(FLOATTYPE), intent(in), optional         :: InbreedingRate        !< Targeted inbreeding rate
      real(FLOATTYPE), intent(in), optional         :: MinInbreedingPct      !< Targeted minimum inbreeding percentage (100 means the minimum possible inbreeding)
      type(AlphaMateModeSpec), intent(in), optional :: ModeMinCoancestrySpec !< Minimum coancestry solution specs
      type(AlphaMateModeSpec), intent(in), optional :: ModeMinInbreedingSpec !< Minimum inbreeding solution specs
      type(AlphaMateModeSpec), intent(in), optional :: ModeMaxCriterionSpec  !< Maximum criterion  solution specs

      if      (present(Degree)) then
        This%TargetDegree = Degree
        This%TargetMaxCriterionPct = Degree2MaxCriterionPct(Degree=This%TargetDegree)
        This%TargetMinCoancestryPct = Degree2MinCoancestryPct(Degree=This%TargetDegree)
        if (present(ModeMinCoancestrySpec) .and. present(ModeMaxCriterionSpec)) then
          This%TargetCoancestryRate = MinCoancestryPct2CoancestryRate(MinCoancestryPct=This%TargetMinCoancestryPct, &
                                                                      MinCoancestryRate=ModeMinCoancestrySpec%CoancestryRate, &
                                                                      MaxCoancestryRate=ModeMaxCriterionSpec%CoancestryRate)
          This%TargetSelCriterionStd = MaxCriterionPct2SelCriterionStd(MaxCriterionPct=This%TargetMaxCriterionPct, &
                                                                       MinSelCriterionStd=ModeMinCoancestrySpec%SelCriterionStd, &
                                                                       MaxSelCriterionStd=ModeMaxCriterionSpec%SelCriterionStd)
          if (present(Data)) then
            This%TargetCoancestry = CoancestryRate2Coancestry(CoancestryRate=This%TargetCoancestryRate, &
                                                              CurrentCoancestry=Data%CoancestryRanMate)
            This%TargetSelCriterion = SelCriterionStd2SelCriterion(SelCriterionStd=This%TargetSelCriterionStd, &
                                                                   Mean=Data%SelCriterionStat%Mean, &
                                                                   Sd=Data%SelCriterionStat%Sd)
          end if
        end if
      else if (present(SelCriterion)) then
        This%TargetSelCriterion = SelCriterion
        if (present(Data)) then
            This%TargetSelCriterionStd = SelCriterion2SelCriterionStd(SelCriterion=This%TargetSelCriterion, &
                                                                      Mean=Data%SelCriterionStat%Mean, &
                                                                      Sd=Data%SelCriterionStat%Sd)
        end if
        if (present(ModeMinCoancestrySpec) .and. present(ModeMaxCriterionSpec)) then
          This%TargetMaxCriterionPct = SelCriterionStd2MaxCriterionPct(SelCriterionStd=This%TargetSelCriterionStd, &
                                                                       MinSelCriterionStd=ModeMinCoancestrySpec%SelCriterionStd, &
                                                                       MaxSelCriterionStd=ModeMaxCriterionSpec%SelCriterionStd)
          This%TargetDegree = MaxCriterionPct2Degree(MaxCriterionPct=This%TargetMaxCriterionPct)
          This%TargetMinCoancestryPct = Degree2MinCoancestryPct(Degree=This%TargetDegree)
          This%TargetCoancestryRate = MinCoancestryPct2CoancestryRate(MinCoancestryPct=This%TargetMinCoancestryPct, &
                                                                      MinCoancestryRate=ModeMinCoancestrySpec%CoancestryRate, &
                                                                      MaxCoancestryRate=ModeMaxCriterionSpec%CoancestryRate)
          if (present(Data)) then
            This%TargetCoancestry = CoancestryRate2Coancestry(CoancestryRate=This%TargetCoancestryRate, &
                                                              CurrentCoancestry=Data%CoancestryRanMate)
          end if
        end if
      else if (present(SelCriterionStd)) then
        This%TargetSelCriterionStd = SelCriterionStd
        if (present(Data)) then
            This%TargetSelCriterion = SelCriterionStd2SelCriterion(SelCriterionStd=This%TargetSelCriterionStd, &
                                                                   Mean=Data%SelCriterionStat%Mean, &
                                                                   Sd=Data%SelCriterionStat%Sd)
        end if
        if (present(ModeMinCoancestrySpec) .and. present(ModeMaxCriterionSpec)) then
          This%TargetMaxCriterionPct = SelCriterionStd2MaxCriterionPct(SelCriterionStd=This%TargetSelCriterionStd, &
                                                                       MinSelCriterionStd=ModeMinCoancestrySpec%SelCriterionStd, &
                                                                       MaxSelCriterionStd=ModeMaxCriterionSpec%SelCriterionStd)
          This%TargetDegree = MaxCriterionPct2Degree(MaxCriterionPct=This%TargetMaxCriterionPct)
          This%TargetMinCoancestryPct = Degree2MinCoancestryPct(Degree=This%TargetDegree)
          This%TargetCoancestryRate = MinCoancestryPct2CoancestryRate(MinCoancestryPct=This%TargetMinCoancestryPct, &
                                                                      MinCoancestryRate=ModeMinCoancestrySpec%CoancestryRate, &
                                                                      MaxCoancestryRate=ModeMaxCriterionSpec%CoancestryRate)
          if (present(Data)) then
            This%TargetCoancestry = CoancestryRate2Coancestry(CoancestryRate=This%TargetCoancestryRate, &
                                                              CurrentCoancestry=Data%CoancestryRanMate)
          end if
        end if
      else if (present(MaxCriterionPct)) then
        This%TargetMaxCriterionPct = MaxCriterionPct
        This%TargetDegree = MaxCriterionPct2Degree(MaxCriterionPct=This%TargetMaxCriterionPct)
        This%TargetMinCoancestryPct = Degree2MinCoancestryPct(Degree=This%TargetDegree)
        if (present(ModeMinCoancestrySpec) .and. present(ModeMaxCriterionSpec)) then
          This%TargetCoancestryRate = MinCoancestryPct2CoancestryRate(MinCoancestryPct=This%TargetMinCoancestryPct, &
                                                                      MinCoancestryRate=ModeMinCoancestrySpec%CoancestryRate, &
                                                                      MaxCoancestryRate=ModeMaxCriterionSpec%CoancestryRate)
          This%TargetSelCriterionStd = MaxCriterionPct2SelCriterionStd(MaxCriterionPct=This%TargetMaxCriterionPct, &
                                                                       MinSelCriterionStd=ModeMinCoancestrySpec%SelCriterionStd, &
                                                                       MaxSelCriterionStd=ModeMaxCriterionSpec%SelCriterionStd)
          if (present(Data)) then
            This%TargetCoancestry = CoancestryRate2Coancestry(CoancestryRate=This%TargetCoancestryRate, &
                                                              CurrentCoancestry=Data%CoancestryRanMate)
            This%TargetSelCriterion = SelCriterionStd2SelCriterion(SelCriterionStd=This%TargetSelCriterionStd, &
                                                                   Mean=Data%SelCriterionStat%Mean, &
                                                                   Sd=Data%SelCriterionStat%Sd)
          end if
        end if
      else if (present(Coancestry)) then
        This%TargetCoancestry = Coancestry
        if (present(Data)) then
          This%TargetCoancestryRate = Coancestry2CoancestryRate(CurrentCoancestry=Data%CoancestryRanMate, &
                                                                FutureCoancestry=This%TargetCoancestry)
        end if
        if (present(ModeMinCoancestrySpec) .and. present(ModeMaxCriterionSpec)) then
          This%TargetMinCoancestryPct = CoancestryRate2MinCoancestryPct(CoancestryRate=This%TargetCoancestryRate, &
                                                                        MinCoancestryRate=ModeMinCoancestrySpec%CoancestryRate, &
                                                                        MaxCoancestryRate=ModeMaxCriterionSpec%CoancestryRate)
          This%TargetDegree = MinCoancestryPct2Degree(MinCoancestryPct=This%TargetMinCoancestryPct)
          This%TargetMaxCriterionPct = Degree2MaxCriterionPct(Degree=This%TargetDegree)
          This%TargetSelCriterionStd = MaxCriterionPct2SelCriterionStd(MaxCriterionPct=This%TargetMaxCriterionPct, &
                                                                       MinSelCriterionStd=ModeMinCoancestrySpec%SelCriterionStd, &
                                                                       MaxSelCriterionStd=ModeMaxCriterionSpec%SelCriterionStd)
          if (present(Data)) then
            This%TargetSelCriterion = SelCriterionStd2SelCriterion(SelCriterionStd=This%TargetSelCriterionStd, &
                                                                   Mean=Data%SelCriterionStat%Mean, &
                                                                   Sd=Data%SelCriterionStat%Sd)
          end if
        end if
      else if (present(CoancestryRate)) then
        This%TargetCoancestryRate = CoancestryRate
        if (present(Data)) then
          This%TargetCoancestry = CoancestryRate2Coancestry(CoancestryRate=This%TargetCoancestryRate, &
                                                            CurrentCoancestry=Data%CoancestryRanMate)
        end if
        if (present(ModeMinCoancestrySpec) .and. present(ModeMaxCriterionSpec)) then
          This%TargetMinCoancestryPct = CoancestryRate2MinCoancestryPct(CoancestryRate=This%TargetCoancestryRate, &
                                                                        MinCoancestryRate=ModeMinCoancestrySpec%CoancestryRate, &
                                                                        MaxCoancestryRate=ModeMaxCriterionSpec%CoancestryRate)
          This%TargetDegree = MinCoancestryPct2Degree(MinCoancestryPct=This%TargetMinCoancestryPct)
          This%TargetMaxCriterionPct = Degree2MaxCriterionPct(Degree=This%TargetDegree)
          This%TargetSelCriterionStd = MaxCriterionPct2SelCriterionStd(MaxCriterionPct=This%TargetMaxCriterionPct, &
                                                                       MinSelCriterionStd=ModeMinCoancestrySpec%SelCriterionStd, &
                                                                       MaxSelCriterionStd=ModeMaxCriterionSpec%SelCriterionStd)
          if (present(Data)) then
            This%TargetSelCriterion = SelCriterionStd2SelCriterion(SelCriterionStd=This%TargetSelCriterionStd, &
                                                                   Mean=Data%SelCriterionStat%Mean, &
                                                                   Sd=Data%SelCriterionStat%Sd)
          end if
        end if
      else if (present(MinCoancestryPct)) then
        This%TargetMinCoancestryPct = MinCoancestryPct
        This%TargetDegree = MinCoancestryPct2Degree(MinCoancestryPct=This%TargetMinCoancestryPct)
        This%TargetMaxCriterionPct = Degree2MaxCriterionPct(Degree=This%TargetDegree)
        if (present(ModeMinCoancestrySpec) .and. present(ModeMaxCriterionSpec)) then
          This%TargetCoancestryRate = MinCoancestryPct2CoancestryRate(MinCoancestryPct=This%TargetMinCoancestryPct, &
                                                                      MinCoancestryRate=ModeMinCoancestrySpec%CoancestryRate, &
                                                                      MaxCoancestryRate=ModeMaxCriterionSpec%CoancestryRate)
          This%TargetSelCriterionStd = MaxCriterionPct2SelCriterionStd(MaxCriterionPct=This%TargetMaxCriterionPct, &
                                                                       MinSelCriterionStd=ModeMinCoancestrySpec%SelCriterionStd, &
                                                                       MaxSelCriterionStd=ModeMaxCriterionSpec%SelCriterionStd)
          if (present(Data)) then
            This%TargetCoancestry = CoancestryRate2Coancestry(CoancestryRate=This%TargetCoancestryRate, &
                                                              CurrentCoancestry=Data%CoancestryRanMate)
            This%TargetSelCriterion = SelCriterionStd2SelCriterion(SelCriterionStd=This%TargetSelCriterionStd, &
                                                                   Mean=Data%SelCriterionStat%Mean, &
                                                                   Sd=Data%SelCriterionStat%Sd)
          end if
        end if
      end if

      if      (present(Inbreeding)) then
        This%TargetInbreeding = Inbreeding
        if (present(Data)) then
          This%TargetInbreedingRate = Coancestry2CoancestryRate(CurrentCoancestry=Data%Inbreeding, &
                                                                FutureCoancestry=This%TargetInbreeding)
        end if
        if (present(ModeMinInbreedingSpec)) then
          This%TargetMinInbreedingPct = CoancestryRate2MinCoancestryPct(CoancestryRate=This%TargetInbreedingRate, &
                                                                        MinCoancestryRate=ModeMinInbreedingSpec%InbreedingRate, &
                                                                        MaxCoancestryRate=ModeMaxCriterionSpec%InbreedingRate)
        end if
      else if (present(InbreedingRate)) then
        This%TargetInbreedingRate = InbreedingRate
        if (present(Data)) then
          This%TargetInbreeding = CoancestryRate2Coancestry(CoancestryRate=This%TargetInbreedingRate, &
                                                            CurrentCoancestry=Data%Inbreeding)
        end if
        if (present(ModeMinInbreedingSpec)) then
          This%TargetMinInbreedingPct = CoancestryRate2MinCoancestryPct(CoancestryRate=This%TargetInbreedingRate, &
                                                                        MinCoancestryRate=ModeMinInbreedingSpec%InbreedingRate, &
                                                                        MaxCoancestryRate=ModeMaxCriterionSpec%InbreedingRate)
        end if
      else if (present(MinInbreedingPct)) then
        This%TargetMinInbreedingPct = MinInbreedingPct
        if (present(ModeMinInbreedingSpec)) then
          This%TargetInbreedingRate = MinCoancestryPct2CoancestryRate(MinCoancestryPct=This%TargetMinInbreedingPct, &
                                                                      MinCoancestryRate=ModeMinInbreedingSpec%InbreedingRate, &
                                                                      MaxCoancestryRate=ModeMaxCriterionSpec%InbreedingRate)
          if (present(Data)) then
            This%TargetInbreeding = CoancestryRate2Coancestry(CoancestryRate=This%TargetInbreedingRate, &
                                                              CurrentCoancestry=Data%Inbreeding)
          end if
        end if
      end if

      if (This%TargetMinCoancestryPct .eq. 0) then
        error stop " ERROR: TargetMinCoancestryPct must be greater than zero!"
      end if
      if (This%TargetMinInbreedingPct .eq. 0) then
        error stop " ERROR: TargetMinInbreedingPct must be greater than zero!"
      end if
      if (This%TargetMaxCriterionPct .eq. 0) then
        error stop " ERROR: TargetMaxCriterionPct must be greater than zero!"
      end if
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Write an AlphaMate optimisation mode specifications
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   April 13, 2017
    !---------------------------------------------------------------------------
    subroutine WriteAlphaMateModeSpec(This, Unit) ! not pure due to IO
      implicit none
      class(AlphaMateModeSpec), intent(in) :: This !< AlphaMateModeSpec holder
      integer(int32), intent(in)           :: Unit !< Unit to write to
      write(Unit, *) "Name: ", trim(This%Name)
      write(Unit, *) "ObjectiveCriterion: ", This%ObjectiveCriterion
      write(Unit, *) "ObjectiveCoancestry: ", This%ObjectiveCoancestry
      write(Unit, *) "ObjectiveInbreeding: ", This%ObjectiveInbreeding
      write(Unit, *) "TargetDegree: ", This%TargetDegree
      write(Unit, *) "TargetSelCriterion: ", This%TargetSelCriterion
      write(Unit, *) "TargetSelCriterionStd: ", This%TargetSelCriterionStd
      write(Unit, *) "TargetMaxCriterionPct: ", This%TargetMaxCriterionPct
      write(Unit, *) "TargetCoancestry: ", This%TargetCoancestry
      write(Unit, *) "TargetCoancestryRate: ", This%TargetCoancestryRate
      write(Unit, *) "CoancestryWeightBelow: ", This%CoancestryWeightBelow
      write(Unit, *) "TargetMinCoancestryPct: ", This%TargetMinCoancestryPct
      write(Unit, *) "TargetInbreeding: ", This%TargetInbreeding
      write(Unit, *) "TargetInbreedingRate: ", This%TargetInbreedingRate
      write(Unit, *) "InbreedingWeightBelow: ", This%InbreedingWeightBelow
      write(Unit, *) "TargetMinInbreedingPct: ", This%TargetMinInbreedingPct
      write(Unit, *) "Degree: ", This%Degree
      write(Unit, *) "SelCriterion: ", This%SelCriterion
      write(Unit, *) "SelCriterionStd: ", This%SelCriterionStd
      write(Unit, *) "MaxCriterionPct: ", This%MaxCriterionPct
      write(Unit, *) "Coancestry: ", This%Coancestry
      write(Unit, *) "CoancestryRate: ", This%CoancestryRate
      write(Unit, *) "MinCoancestryPct: ", This%MinCoancestryPct
      write(Unit, *) "Inbreeding: ", This%Inbreeding
      write(Unit, *) "InbreedingRate: ", This%InbreedingRate
      write(Unit, *) "MinInbreedingPct: ", This%MinInbreedingPct
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Log an AlphaMate optimisation mode specifications
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   April 18, 2017
    !---------------------------------------------------------------------------
    subroutine LogTargetsAlphaMateModeSpec(This, Spec, Unit) ! not pure due to IO
      implicit none
      class(AlphaMateModeSpec), intent(in) :: This !< AlphaMateModeSpec holder
      type(AlphaMateSpec), intent(in)      :: Spec !< AlphaMateSpec holder
      integer(int32), intent(in)           :: Unit !< Unit to write to

      if (This%ObjectiveCriterion) then
        write(Unit, "(a)") "   Selection criterion / standardized"
        write(Unit, "(a)") "     @MinCoancestry: "//trim(Real2Char(Spec%ModeMinCoancestrySpec%SelCriterion,     fmt=FMTREAL2CHAR))//" /" &
                                                  //trim(Real2Char(Spec%ModeMinCoancestrySpec%SelCriterionStd, fmt=FMTREAL2CHAR))
        write(Unit, "(a)") "     @MaxCriterion:  "//trim(Real2Char(Spec%ModeMaxCriterionSpec%SelCriterion,      fmt=FMTREAL2CHAR))//" /" &
                                                  //trim(Real2Char(Spec%ModeMaxCriterionSpec%SelCriterionStd,  fmt=FMTREAL2CHAR))
        write(Unit, "(a)") "     Pct of max:     "//trim(Real2Char(This%TargetMaxCriterionPct,  fmt="(f7.1)"))
        write(Unit, "(a)") "     Target:         "//trim(Real2Char(This%TargetSelCriterion,     fmt=FMTREAL2CHAR))//" /" &
                                                  //trim(Real2Char(This%TargetSelCriterionStd, fmt=FMTREAL2CHAR))
        write(Unit, "(a)") " "
      end if

      if (This%ObjectiveCoancestry) then
        write(Unit, "(a)") "   Coancestry coefficient / rate"
        write(Unit, "(a)") "     @MinCoancestry: "//trim(Real2Char(Spec%ModeMinCoancestrySpec%Coancestry, fmt=FMTREAL2CHAR))//" /" &
                                                  //trim(Real2Char(Spec%ModeMinCoancestrySpec%CoancestryRate, fmt=FMTREAL2CHAR))
        write(Unit, "(a)") "     @MaxCriterion:  "//trim(Real2Char(Spec%ModeMaxCriterionSpec%Coancestry, fmt=FMTREAL2CHAR))//" /" &
                                                  //trim(Real2Char(Spec%ModeMaxCriterionSpec%CoancestryRate, fmt=FMTREAL2CHAR))
        write(Unit, "(a)") "     Pct of min:     "//trim(Real2Char(This%TargetMinCoancestryPct, fmt="(f7.1)"))
        write(Unit, "(a)") "     Target:         "//trim(Real2Char(This%TargetCoancestry,     fmt=FMTREAL2CHAR))//" /" &
                                                  //trim(Real2Char(This%TargetCoancestryRate, fmt=FMTREAL2CHAR))
        write(Unit, "(a)") " "
      end if

      if (This%ObjectiveCriterion .and. This%ObjectiveCoancestry) then
        write(Unit, "(a)") "   Degree"
        write(Unit, "(a)") "     Target:         "//trim(Real2Char(This%TargetDegree, fmt="(f7.1)"))
        write(Unit, "(a)") " "
      end if

      if (This%ObjectiveInbreeding) then
        write(Unit, "(a)") "   Inbreeding coefficient / rate"
        write(Unit, "(a)") "     @MinInbreeding: "//trim(Real2Char(Spec%ModeMinInbreedingSpec%Inbreeding, fmt=FMTREAL2CHAR))//" /" &
                                                  //trim(Real2Char(Spec%ModeMinInbreedingSpec%Inbreedingrate, fmt=FMTREAL2CHAR))
        write(Unit, "(a)") "     @MaxCriterion:  "//trim(Real2Char(Spec%ModeMaxCriterionSpec%Coancestry, fmt=FMTREAL2CHAR))//" /" &
                                                  //trim(Real2Char(Spec%ModeMaxCriterionSpec%CoancestryRate, fmt=FMTREAL2CHAR))
        write(Unit, "(a)") "     Pct of min:     "//trim(Real2Char(This%TargetMinInbreedingPct, fmt="(f7.1)"))
        write(Unit, "(a)") "     Target:         "//trim(Real2Char(This%TargetInbreeding,     fmt=FMTREAL2CHAR))//" /" &
                                                  //trim(Real2Char(This%TargetInbreedingRate, fmt=FMTREAL2CHAR))
        write(Unit, "(a)") " "
      end if
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Read AlphaMate data from a file, check specifications against the
    !!         data, and summarize data for further use
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !---------------------------------------------------------------------------
    subroutine ReadAlphaMateData(This, Spec, LogStdout) ! not pure due to IO
      implicit none
      class(AlphaMateData), intent(out)  :: This      !< @return AlphaMateData holder
      type(AlphaMateSpec), intent(inout) :: Spec      !< AlphaMateSpec holder (inout because we save some info in Spec that is based on Data)
      logical, optional                  :: LogStdout !< Log process on stdout (default .false.)

      integer(int32) :: Ind, IndLoc, IndLoc2, nIndTmp, Mat, nMatTmp, GenderTmp, jMal, jFem, IndPair(2), Crit
      integer(int32) :: SelCriterionUnit, GenderUnit, GenericIndCritUnit, GenericMatCritUnit
      integer(int32) :: CoancestrySummaryUnit, InbreedingSummaryUnit, CriterionSummaryUnit
      integer(int32) :: GenericIndCritSummaryUnit, GenericMatCritSummaryUnit

      real(FLOATTYPE) :: SelCriterionTmp, SelCriterionTmp2
      real(FLOATTYPE), allocatable :: GenericIndCritTmp(:), GenericMatCritTmp(:)

      character(len=IDLENGTH) :: IdCTmp, IdCTmp2

      logical :: LogStdoutInternal
      if (present(LogStdout)) then
        LogStdoutInternal = LogStdout
      else
        LogStdoutInternal = .false.
      end if

      if (LogStdoutInternal) then
        write(STDOUT, "(a)") " "
      end if

      ! --- Coancestry or Numerator Relationship Matrix ---

      call This%Coancestry%Read(File=Spec%RelMtxFile)
      if (Spec%NrmInsteadOfCoancestry) then
        call This%Coancestry%Nrm2Coancestry
      end if
      This%nInd = This%Coancestry%nInd
       ! Would not normally put data into spec, but need to do it for OO-flexibility (InitialiseAlphaMateSol)
      Spec%nInd = This%Coancestry%nInd
      allocate(This%AvgCoancestry(This%nInd))
      do Ind = 1, This%nInd
        This%AvgCoancestry(Ind) = Mean(This%Coancestry%Value(1:, Ind))
      end do

      if (LogStdoutInternal) then
        write(STDOUT, "(a)") " Number of individuals in the coancestry matrix file: "//trim(Int2Char(This%nInd))
      end if

      ! --- Selection criterion ---

      allocate(This%SelCriterion(This%nInd))
      allocate(This%SelCriterionStd(This%nInd))
      if (Spec%PAGEPar) then
        allocate(This%SelCriterionPAGE(This%nInd))
        allocate(This%SelCriterionStdPAGE(This%nInd))
      end if

      if (.not. Spec%SelCriterionGiven) then
        This%SelCriterion     = 0.0
        This%SelCriterionStd = 0.0
        if (Spec%PreselectPar) then
          write(STDERR, "(a)") " ERROR: Can not preselect when selection criterion information is not provided!"
          write(STDERR, "(a)") " "
          stop 1
        end if
      else
        nIndTmp = CountLines(Spec%SelCriterionFile)
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " Number of individuals in the selection criterion file: "//trim(Int2Char(nIndTmp))
        end if
        if (nIndTmp .ne. This%nInd) then
          write(STDERR, "(a)") " ERROR: Number of individuals in the selection criterion file and the coancestry matrix file is not the same!"
          write(STDERR, "(a)") " ERROR: Number of individuals in the coancestry matrix file:   "//trim(Int2Char(This%nInd))
          write(STDERR, "(a)") " ERROR: Number of individuals in the selection criterion file: "//trim(Int2Char(nIndTmp))
          write(STDERR, "(a)") " "
          stop 1
        end if
        open(newunit=SelCriterionUnit, file=Spec%SelCriterionFile, status="old")
        do Ind = 1, This%nInd
          if (Spec%PAGEPar) then
            read(SelCriterionUnit, *) IdCTmp, SelCriterionTmp, SelCriterionTmp2
          else
            read(SelCriterionUnit, *) IdCTmp, SelCriterionTmp
          end if
          IndLoc = This%Coancestry%OriginalIdDict%GetValue(key=IdCTmp);
          if (IndLoc .eq. DICT_NULL) then
            write(STDERR, "(a)") " ERROR: Individual "//trim(IdCTmp)//" from the selection criterion file not present in the coancestry matrix file!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          This%SelCriterion(IndLoc) = SelCriterionTmp
          if (Spec%PAGEPar) then
            This%SelCriterionPAGE(IndLoc) = SelCriterionTmp2
          end if
        end do
        close(SelCriterionUnit)
      end if

      ! --- Gender (and number of parents when not provided as a limit/constraint) ---

      allocate(This%Gender(This%nInd))
      if (.not. Spec%GenderGiven) then
        This%Gender = 0
        if (Spec%nPar .eq. 0) then
          write(STDOUT, "(a)") " NOTE: When number of parents is not provided all candidates could become parents (though other constraints apply)"
          Spec%nPar = This%nInd
        end if
        Spec%nPar1 = Spec%nPar
        if (Spec%EqualizePar) then
          if (Spec%nMat .lt. Spec%nPar) then
            write(STDERR, "(a)") " ERROR: Number of parents must be smaller or equal to the number of matings/crosses when EqualizeContributions is active!"
            write(STDERR, "(a)") " ERROR: The number of parents:         "//trim(Int2Char(Spec%nPar))
            write(STDERR, "(a)") " ERROR: The number of matings/crosses: "//trim(Int2Char(Spec%nMat))
            write(STDERR, "(a)") " "
            stop 1
          end if
          if (mod(Spec%nMat, Spec%nPar) .gt. 0) then
            write(STDERR, "(a)") " ERROR: The number of parents and the number of matings/crosses must divide without remainder when EqualizeContributions is active!"
            write(STDERR, "(a)") " ERROR: The number of parents:         "//trim(Int2Char(Spec%nPar))
            write(STDERR, "(a)") " ERROR: The number of matings/crosses: "//trim(Int2Char(Spec%nMat))
            write(STDERR, "(a)") " ERROR: The remainder:                 "//trim(Int2Char(mod(Spec%nMat, Spec%nPar)))
            write(STDERR, "(a)") " "
            stop 1
          end if
        end if
        if (Spec%LimitPar) then
          if ((Spec%nPar * Spec%LimitParMax) .lt. Spec%nMat) then
            write(STDERR, "(a)") " ERROR: The number of parents * LimitContributionsMax is to small to achieve specified number of matings/crosses!"
            write(STDERR, "(a)") " ERROR: The number of parents:         "//trim(Int2Char(Spec%nPar))
            write(STDERR, "(a)") " ERROR: LimitContributionsMax:         "//trim(Int2Char(nint(Spec%LimitParMax)))
            write(STDERR, "(a)") " ERROR:         their product:         "//trim(Int2Char(Spec%nPar * nint(Spec%LimitParMax)))
            write(STDERR, "(a)") " ERROR: The number of matings/crosses: "//trim(Int2Char(Spec%nMat))
            write(STDERR, "(a)") " "
            stop 1
          end if
        end if
      else
        nIndTmp = CountLines(Spec%GenderFile)
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " Number of individuals in the gender file: "//trim(Int2Char(nIndTmp))
        end if
        if (nIndTmp .ne. This%nInd) then
          write(STDERR, "(a)") " ERROR: Number of individuals in the gender file and the coancestry matrix file is not the same!"
          write(STDERR, "(a)") " ERROR: Number of individuals in the coancestry matrix file: "//trim(Int2Char(This%nInd))
          write(STDERR, "(a)") " ERROR: Number of individuals in the gender file:            "//trim(Int2Char(nIndTmp))
          write(STDERR, "(a)") " "
          stop 1
        end if

        This%nMal = 0
        This%nFem = 0

        open(newunit=GenderUnit, file=Spec%GenderFile, status="old")
        do Ind = 1, This%nInd
          read(GenderUnit, *) IdCTmp, GenderTmp
          if      (GenderTmp .eq. 1) then
            This%nMal = This%nMal + 1
          else if (GenderTmp .eq. 2) then
            This%nFem = This%nFem + 1
          else
            write(STDERR, "(a)") " ERROR: Gender code must be either 1 for male individuals or 2 for female individuals!"
            write(STDERR, "(a)") " ERROR: "//trim(Int2Char(Ind))//" "//trim(IdCTmp)//" "//trim(Int2Char(GenderTmp))
            write(STDERR, "(a)") " "
            stop 1
          end if
          IndLoc = This%Coancestry%OriginalIdDict%GetValue(key=IdCTmp);
          if (IndLoc .eq. DICT_NULL) then
            write(STDERR, "(a)") " ERROR: Individual "//trim(IdCTmp)//" from the gender file not present in the coancestry matrix file!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          This%Gender(IndLoc) = GenderTmp
        end do
        close(GenderUnit)

        write(STDOUT, "(a)") " Number of   males: "//trim(Int2Char(This%nMal))
        write(STDOUT, "(a)") " Number of females: "//trim(Int2Char(This%nFem))

        if (Spec%nPar1 .eq. 0) then
          write(STDOUT, "(a)") " NOTE: When number of male parents is not provided all male candidates could become parents (though other constraints apply)"
          Spec%nPar1 = This%nMal
        end if
        if (Spec%EqualizePar1) then
          if (Spec%nMat .lt. Spec%nPar1) then
            write(STDERR, "(a)") " ERROR: The number of male parents must be smaller or equal to the number of matings/crosses when EqualizeMaleContributions is active!"
            write(STDERR, "(a)") " ERROR: The number of male parents:    "//trim(Int2Char(Spec%nPar1))
            write(STDERR, "(a)") " ERROR: The number of matings/crosses: "//trim(Int2Char(Spec%nMat))
            write(STDERR, "(a)") " "
            stop 1
          end if
          if (mod(Spec%nMat, Spec%nPar1) .gt. 0) then
            write(STDERR, "(a)") " ERROR: The number of male parents and the number of matings/crosses must divide without remainder when EqualizeMaleContributions is active!"
            write(STDERR, "(a)") " ERROR: The number of male parents:    "//trim(Int2Char(Spec%nPar1))
            write(STDERR, "(a)") " ERROR: The number of matings/crosses: "//trim(Int2Char(Spec%nMat))
            write(STDERR, "(a)") " ERROR: The remainder:                 "//trim(Int2Char(mod(Spec%nMat, Spec%nPar1)))
            write(STDERR, "(a)") " "
            stop 1
          end if
        end if
        if (Spec%LimitPar1) then
          if ((Spec%nPar1 * Spec%LimitPar1Max) .lt. Spec%nMat) then
            write(STDERR, "(a)") " ERROR: The number of male parents * LimitMaleContributionsMax is to small to achieve specified number of matings/crosses!"
            write(STDERR, "(a)") " ERROR: The number of male parents:    "//trim(Int2Char(Spec%nPar1))
            write(STDERR, "(a)") " ERROR: LimitMaleContributionsMax:     "//trim(Int2Char(nint(Spec%LimitPar1Max)))
            write(STDERR, "(a)") " ERROR:             their product:     "//trim(Int2Char(Spec%nPar1 * nint(Spec%LimitPar1Max)))
            write(STDERR, "(a)") " ERROR: The number of matings/crosses: "//trim(Int2Char(Spec%nMat))
            write(STDERR, "(a)") " "
            stop 1
          end if
        end if

        if (Spec%nPar2 .eq. 0) then
          write(STDOUT, "(a)") " NOTE: When number of female parents is not provided all female candidates could become parents (though other constraints apply)"
          Spec%nPar2 = This%nFem
        end if
        if (Spec%EqualizePar2) then
          if (Spec%nMat .lt. Spec%nPar2) then
            write(STDERR, "(a)") " ERROR: Number of female parents must be smaller or equal to the number of matings/crosses when EqualizeFemaleContributions is active!"
            write(STDERR, "(a)") " ERROR: The number of female parents:  "//trim(Int2Char(Spec%nPar2))
            write(STDERR, "(a)") " ERROR: The number of matings/crosses: "//trim(Int2Char(Spec%nMat))
            write(STDERR, "(a)") " "
            stop 1
          end if
          if (mod(Spec%nMat, Spec%nPar2) .gt. 0) then
            write(STDERR, "(a)") " ERROR: The number of female parents and the number of matings/crosses must divide without remainder when EqualizeFemaleContributions is active!"
            write(STDERR, "(a)") " ERROR: The number of female parents:  "//trim(Int2Char(Spec%nPar2))
            write(STDERR, "(a)") " ERROR: The number of matings/crosses: "//trim(Int2Char(Spec%nMat))
            write(STDERR, "(a)") " ERROR: The remainder:                 "//trim(Int2Char(mod(Spec%nMat, Spec%nPar2)))
            write(STDERR, "(a)") " "
            stop 1
          end if
        end if
        if (Spec%LimitPar2) then
          if ((Spec%nPar2 * Spec%LimitPar2Max) .lt. Spec%nMat) then
            write(STDERR, "(a)") " ERROR: The number of female parents * LimitFemaleContributionsMax is to small to achieve specified number of matings/crosses!"
            write(STDERR, "(a)") " ERROR: The number of female parents:  "//trim(Int2Char(Spec%nPar2))
            write(STDERR, "(a)") " ERROR: LimitFemaleContributionsMax:   "//trim(Int2Char(nint(Spec%LimitPar2Max)))
            write(STDERR, "(a)") " ERROR:               their product:   "//trim(Int2Char(Spec%nPar2 * nint(Spec%LimitPar2Max)))
            write(STDERR, "(a)") " ERROR: The number of matings/crosses: "//trim(Int2Char(Spec%nMat))
            write(STDERR, "(a)") " "
            stop 1
          end if
        end if

        Spec%nPar = Spec%nPar1 + Spec%nPar2

        if (Spec%nPar1 .gt. This%nMal) then
          write(STDERR, "(a)") " ERROR: Number of male parents can not be larger than number of males"
          write(STDERR, "(a)") " ERROR: Number of male parents: "//trim(Int2Char(Spec%nPar1))
          write(STDERR, "(a)") " ERROR: Number of males:        "//trim(Int2Char(This%nMal))
          write(STDERR, "(a)") " "
          stop 1
        end if
        if (Spec%nPar2 .gt. This%nFem) then
          write(STDERR, "(a)") " ERROR: Number of female parents can not be larger than number of females"
          write(STDERR, "(a)") " ERROR: Number of female parents: "//trim(Int2Char(Spec%nPar2))
          write(STDERR, "(a)") " ERROR: Number of females:        "//trim(Int2Char(This%nFem))
          write(STDERR, "(a)") " "
          stop 1
        end if
      end if

      if (Spec%nPar .lt. 1) then
        write(STDERR, "(a)") " ERROR: Number of specified/derived parents must be larger than zero!"
        write(STDERR, "(a)") " "
        stop 1
      end if

      if (Spec%GenderGiven .and. ((Spec%nPar1 .lt. 1) .or. (Spec%nPar2 .lt. 1))) then
        write(STDERR, "(a)") " ERROR: Number of specified/derived parents must be larger than zero!"
        write(STDERR, "(a)") " ERROR: Number of   male parents: "//trim(Int2Char(Spec%nPar1))
        write(STDERR, "(a)") " ERROR: Number of female parents: "//trim(Int2Char(Spec%nPar2))
        write(STDERR, "(a)") " "
        stop 1
      end if

      if (This%nInd .lt. Spec%nPar) then
        write(STDERR, "(a)") " ERROR: Number of individuals can not be smaller than number of parents!"
        write(STDERR, "(a)") " ERROR: Number of individuals:               "//trim(Int2Char(This%nInd))
        write(STDERR, "(a)") " ERROR: Number of specified/derived parents: "//trim(Int2Char(Spec%nPar))
        write(STDERR, "(a)") " "
        stop 1
      end if

      if (Spec%PreselectPar) then
        if (Spec%GenderGiven) then
          Spec%PreselectPar1N = nint(Spec%PreselectPar1Pct / 100.0 * This%nMal)
          if (Spec%PreselectPar1N .lt. Spec%nPar1) then
            write(STDERR, "(a)") " ERROR: Can not preselect less males than there should be male parents!"
            write(STDERR, "(a)") " ERROR: Number of males:        "//trim(Int2Char(This%nMal))
            write(STDERR, "(a)") " ERROR: Number of male parents: "//trim(Int2Char(Spec%nPar1))
            write(STDERR, "(a)") " ERROR: Preselect percentage:   "//trim(Int2Char(nint(Spec%PreselectPar1Pct)))
            write(STDERR, "(a)") " ERROR: Preselect number:       "//trim(Int2Char(Spec%PreselectPar1N))
            write(STDERR, "(a)") " "
            stop 1
          end if
          Spec%PreselectPar2N = nint(Spec%PreselectPar2Pct / 100.0 * This%nFem)
          if (Spec%PreselectPar2N .lt. Spec%nPar2) then
            write(STDERR, "(a)") " ERROR: Can not preselect less females than there should be female parents!"
            write(STDERR, "(a)") " ERROR: Number of females:        "//trim(Int2Char(This%nFem))
            write(STDERR, "(a)") " ERROR: Number of female parents: "//trim(Int2Char(Spec%nPar2))
            write(STDERR, "(a)") " ERROR: Preselect percentage:     "//trim(Int2Char(nint(Spec%PreselectPar2Pct)))
            write(STDERR, "(a)") " ERROR: Preselect number:         "//trim(Int2Char(Spec%PreselectPar2N))
            write(STDERR, "(a)") " "
            stop 1
          end if
        else
          Spec%PreselectParN = nint(Spec%PreselectParPct / 100.0 * This%nInd)
          Spec%PreselectPar1N = Spec%PreselectParN
          if (Spec%PreselectParN .lt. Spec%nPar) then
            write(STDERR, "(a)") " ERROR: Can not preselect less individuals than there should be parents!"
            write(STDERR, "(a)") " ERROR: Number of individuals: "//trim(Int2Char(This%nInd))
            write(STDERR, "(a)") " ERROR: Number of parents:     "//trim(Int2Char(Spec%nPar))
            write(STDERR, "(a)") " ERROR: Preselect percentage:  "//trim(Int2Char(nint(Spec%PreselectParPct)))
            write(STDERR, "(a)") " ERROR: Preselect number:      "//trim(Int2Char(Spec%PreselectParN))
            write(STDERR, "(a)") " "
            stop 1
          end if
        end if
      end if

      if (Spec%PAGEPar) then
        if (Spec%GenderGiven) then
          if (Spec%PAGEPar1Max .gt. Spec%nPar1) then
            write(STDERR, "(a)") " ERROR: Can not PAGE more males than there are male parents!"
            write(STDERR, "(a)") " ERROR: Number of male parents:      "//trim(Int2Char(Spec%nPar1))
            write(STDERR, "(a)") " ERROR: Max number of male for PAGE: "//trim(Int2Char(Spec%PAGEPar1Max))
            write(STDERR, "(a)") " "
            stop 1
          end if
          if (Spec%PAGEPar2Max .gt. Spec%nPar2) then
            write(STDERR, "(a)") " ERROR: Can not PAGE more females than there are female parents!"
            write(STDERR, "(a)") " ERROR: Number of female parents:      "//trim(Int2Char(Spec%nPar2))
            write(STDERR, "(a)") " ERROR: Max number of female for PAGE: "//trim(Int2Char(Spec%PAGEPar2Max))
            write(STDERR, "(a)") " "
            stop 1
          end if
        else
          if (Spec%PAGEParMax .gt. Spec%nPar) then
            write(STDERR, "(a)") " ERROR: Can not PAGE more individuals than there are parents!"
            write(STDERR, "(a)") " ERROR: Number of parents:                  "//trim(Int2Char(Spec%nPar))
            write(STDERR, "(a)") " ERROR: Max number of individuals for PAGE: "//trim(Int2Char(Spec%PAGEParMax))
            write(STDERR, "(a)") " "
            stop 1
          end if
        end if
      end if

      ! --- Define potential parents ---

      ! IdPotPar1   holds   male ids in the form 1:nInd
      ! IdPotPar2   holds female ids in the form 1:nInd
      ! IdPotParSeq holds   both ids, but males are in the form 1:nMal and females in the form 1:nFem

      if (.not. Spec%GenderGiven) then
        This%nPotPar1 = This%nInd
        This%nPotPar2 = This%nInd
        This%nPotPar = This%nInd
        allocate(This%IdPotPar1(This%nPotPar1))
        do Ind = 1, This%nInd
          This%IdPotPar1(Ind) = Ind
        end do
      else
        This%nPotPar1 = This%nMal
        This%nPotPar2 = This%nFem
        This%nPotPar = This%nPotPar1 + This%nPotPar2
        allocate(This%IdPotPar1(This%nPotPar1))
        allocate(This%IdPotPar2(This%nPotPar2))
        allocate(This%IdPotParSeq(This%nInd))
        jMal = 0
        jFem = 0
        do Ind = 1, This%nInd
          if (This%Gender(Ind) .eq. 1) then
            jMal = jMal + 1
            This%IdPotPar1(jMal) = Ind
            This%IdPotParSeq(Ind) = jMal
          else
            jFem = jFem + 1
            This%IdPotPar2(jFem) = Ind
            This%IdPotParSeq(Ind) = jFem
          end if
        end do
      end if

      ! --- Number of all potential matings ---

      if (Spec%GenderGiven) then
        This%nPotMat = This%nPotPar1 * This%nPotPar2
      else
        This%nPotMat = FLOATFUN(This%nPotPar1 * This%nPotPar1) / 2.0
        if (Spec%SelfingAllowed) then
          This%nPotMat = nint(This%nPotMat + real(This%nPotPar1) / 2.0)
        else
          This%nPotMat = nint(This%nPotMat - real(This%nPotPar1) / 2.0)
        end if
      end if

      if (Spec%nMat .gt. This%nPotMat) then
        ! @todo Might need to change this if we would do in-vitro fertilisation
        write(STDERR, "(a)") " ERROR: Number of specified matings/crosses is larger than the number of all potential matings/crosses!"
        write(STDERR, "(a)") " ERROR: Number of     specified matings/crosses: "//trim(Int2Char(Spec%nMat))
        write(STDERR, "(a)") " ERROR: Number of all potential matings/crosses: "//trim(Int2Char(This%nPotMat))
        if (Spec%GenderGiven) then
          write(STDERR, "(a)") " ERROR: = no. of males * no. of females"
          write(STDERR, "(a)") " ERROR: = (no. of males = "//trim(Int2Char(This%nPotPar1))//", no. of females = "//trim(Int2Char(This%nPotPar2))
        else
          if (Spec%SelfingAllowed) then
            write(STDERR, "(a)") " ERROR: = half-diallel including selfing"
            write(STDERR, "(a)") " ERROR: = no. of individuals * no. of individuals / 2 + individuals / 2"
          else
            write(STDERR, "(a)") " ERROR: = half-diallel excluding selfing"
            write(STDERR, "(a)") " ERROR: = no. of individuals * no. of individuals / 2 - individuals / 2"
          end if
          write(STDERR, "(a)") " ERROR:   (no. of individuals = "//trim(Int2Char(This%nPotPar1))//")"
        end if
        write(STDERR, "(a)") " "
        stop 1
      end if

      ! --- Generic individual criterion ---

      ! @todo standardize these input values!!!

      if (Spec%GenericIndCritGiven) then
        nIndTmp = CountLines(Spec%GenericIndCritFile)
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " Number of individuals in the generic individual criterion file: "//trim(Int2Char(nIndTmp))
        end if
        if (nIndTmp .ne. This%nInd) then
          write(STDERR, "(a)") " ERROR: Number of individuals in the generic individual criterion file and the coancestry matrix file is not the same!"
          write(STDERR, "(a)") " ERROR: Number of individuals in the generic individual criterion file: "//trim(Int2Char(nIndTmp))
          write(STDERR, "(a)") " ERROR: Number of individuals in the coancestry matrix file:            "//trim(Int2Char(This%nInd))
          write(STDERR, "(a)") " "
          stop 1
        end if
        allocate(This%GenericIndCrit(This%nInd, Spec%nGenericIndCrit))
        allocate(GenericIndCritTmp(Spec%nGenericIndCrit))
        This%GenericIndCrit = 0.0
        open(newunit=GenericIndCritUnit, file=Spec%GenericIndCritFile, status="old")
        do Ind = 1, This%nInd
          read(GenericIndCritUnit, *) IdCTmp, GenericIndCritTmp
          IndLoc = This%Coancestry%OriginalIdDict%GetValue(key=IdCTmp)
          if (IndLoc .eq. DICT_NULL) then
            write(STDERR, "(a)") " ERROR: Individual "//trim(IdCTmp)//" from the generic individual criterion file not present in the coancestry matrix file!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          This%GenericIndCrit(IndLoc, :) = GenericIndCritTmp
        end do
        close(GenericIndCritUnit)
      end if

      ! --- Generic mating criterion ---

      ! @todo standardize these input values!!!

      if (Spec%GenericMatCritGiven) then
        nMatTmp = CountLines(Spec%GenericMatCritFile)
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " Number of matings/crosses in the generic mating/crossing criterion file: "//trim(Int2Char(nMatTmp))
        end if
        if (nMatTmp .ne. This%nPotMat) then
          write(STDERR, "(a)") " ERROR: Number of matings/crosses in the generic mating/crossing criterion file and the number of all potential matings/crosses is not the same!"
          write(STDERR, "(a)") " ERROR: Number of matings/crosses in the generic mating/crossing criterion file: "//trim(Int2Char(nMatTmp))
          write(STDERR, "(a)") " ERROR: Number of all potential matings/crosses:                                 "//trim(Int2Char(This%nPotMat))
          write(STDERR, "(a)") " "
          stop 1
        end if
        allocate(This%GenericMatCrit(This%nPotPar1, This%nPotPar2, Spec%nGenericMatCrit))
        allocate(GenericMatCritTmp(Spec%nGenericMatCrit))
        This%GenericMatCrit = 0.0
        open(newunit=GenericMatCritUnit, file=Spec%GenericMatCritFile, status="old")
        do Mat = 1, This%nPotMat
          read(GenericMatCritUnit, *) IdCTmp, IdCTmp2, GenericMatCritTmp
          IndLoc = This%Coancestry%OriginalIdDict%GetValue(key=IdCTmp)
          if (IndLoc .eq. DICT_NULL) then
            write(STDERR, "(a)") " ERROR: Individual "//trim(IdCTmp)//" from the generic mating/crossing criterion file not present in the coancestry matrix file!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          IndLoc2 = This%Coancestry%OriginalIdDict%GetValue(key=IdCTmp2)
          if (IndLoc2 .eq. DICT_NULL) then
            write(STDERR, "(a)") " ERROR: Individual "//trim(IdCTmp2)//" from the generic mating/crossing criterion file not present in the coancestry matrix file!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          if (Spec%GenderGiven) then
            jMal = This%IdPotParSeq(IndLoc)
            if (This%Gender(IndLoc) .ne. 1) then
              write(STDERR, "(a)") " ERROR: Individual "//trim(IdCTmp)//" from the first column in the generic mating/crossing criterion file should be a male!"
              write(STDERR, "(a)") " ERROR: Generic mating/crossing criterion file:"
              write(STDERR, "(a)") " ERROR:   - line:                  "//trim(Int2Char(Mat))
              write(STDERR, "(a)") " ERROR:   - individual 1   (male): "//trim(IdCTmp) //" gender "//trim(Int2Char(This%Gender(IndLoc)))
              write(STDERR, "(a)") " ERROR:   - individual 2 (female): "//trim(IdCTmp2)//" gender "//trim(Int2Char(This%Gender(IndLoc2)))
              write(STDERR, "(a)") " "
              stop 1
            end if
            jFem = This%IdPotParSeq(IndLoc2)
            if (This%Gender(IndLoc2) .ne. 2) then
              write(STDERR, "(a)") " ERROR: Individual "//trim(IdCTmp2)//" from the second column in the generic mating/crossing criterion file should be a female!"
              write(STDERR, "(a)") " ERROR: Generic mating/crossing criterion file:"
              write(STDERR, "(a)") " ERROR:   - line:                  "//trim(Int2Char(Mat))
              write(STDERR, "(a)") " ERROR:   - individual 1   (male): "//trim(IdCTmp) //" gender "//trim(Int2Char(This%Gender(IndLoc)))
              write(STDERR, "(a)") " ERROR:   - individual 2 (female): "//trim(IdCTmp2)//" gender "//trim(Int2Char(This%Gender(IndLoc2)))
              write(STDERR, "(a)") " "
              stop 1
            end if
            ! fill full-matrix
            This%GenericMatCrit(jMal, jFem, :) = GenericMatCritTmp
          else
            ! fill lower-triangle (half-diallel)
            IndPair = [IndLoc, IndLoc2]
            This%GenericMatCrit(maxval(IndPair), minval(IndPair), :) = GenericMatCritTmp
          end if
        end do
        close(GenericMatCritUnit)
      end if

      ! --- Current coancestry ---

      if (LogStdoutInternal) then
        write(STDOUT, "(a)") " "
        write(STDOUT, "(a)") " Current coancestry (average identity of the four genome combinations of two individuals)"
      end if

      This%CoancestryStat = DescStatMatrix(FLOATFUN(This%Coancestry%Value(1:, 1:)))
      if (Spec%GenderGiven) then
        This%CoancestryStatGender1    = DescStatMatrix(FLOATFUN(This%Coancestry%Value(This%IdPotPar1, This%IdPotPar1)))
        This%CoancestryStatGender2    = DescStatMatrix(FLOATFUN(This%Coancestry%Value(This%IdPotPar2, This%IdPotPar2)))
        This%CoancestryStatGenderDiff = DescStatMatrix(FLOATFUN(This%Coancestry%Value(This%IdPotPar1, This%IdPotPar2)))
      end if

      This%CoancestryRanMate       = This%CoancestryStat%All%Mean
      This%CoancestryRanMateNoSelf = This%CoancestryStat%OffDiag%Mean
      if (Spec%GenderGiven) then
        This%CoancestryGenderMate  = This%CoancestryStatGenderDiff%All%Mean
      end if

      if (LogStdoutInternal) then
        write(STDOUT, "(a)") " "
        write(STDOUT, "(a)") "  - coancestry among/between individuals"
        write(STDOUT, "(a)") "                     Among     Between"
        write(STDOUT, "(a)") "    - n:       "//trim( Int2Char(This%CoancestryStat%All%n,    fmt=FMTINT2CHAR)) //" "//trim( Int2Char(This%CoancestryStat%OffDiag%n,      fmt=FMTINT2CHAR))
        write(STDOUT, "(a)") "    - average: "//trim(Real2Char(This%CoancestryStat%All%Mean, fmt=FMTREAL2CHAR))//" "//trim(Real2Char(This%CoancestryStat%OffDiag%Mean,   fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(This%CoancestryStat%All%Sd,   fmt=FMTREAL2CHAR))//" "//trim(Real2Char(This%CoancestryStat%OffDiag%Sd,     fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(This%CoancestryStat%All%Min,  fmt=FMTREAL2CHAR))//" "//trim(Real2Char(This%CoancestryStat%OffDiag%Min,    fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(This%CoancestryStat%All%Max,  fmt=FMTREAL2CHAR))//" "//trim(Real2Char(This%CoancestryStat%OffDiag%Max,    fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "    Among   = coancestry among   individuals (including self-coancestry) = expected inbreeding in their progeny under random mating/crossing, including selfing"
        write(STDOUT, "(a)") "    Between = coancestry between individuals (excluding self-coancestry) = expected inbreeding in their progeny under random mating/crossing, excluding selfing"

        if (Spec%GenderGiven) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") "  - coancestry between males and females"
          write(STDOUT, "(a)") "    (=expected inbreeding in their progeny under random mating/crossing between genders)"
          write(STDOUT, "(a)") "    - n:       "//trim( Int2Char(This%CoancestryStatGenderDiff%All%n,    fmt=FMTINT2CHAR))
          write(STDOUT, "(a)") "    - average: "//trim(Real2Char(This%CoancestryStatGenderDiff%All%Mean, fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(This%CoancestryStatGenderDiff%All%Sd,   fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(This%CoancestryStatGenderDiff%All%Min,  fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(This%CoancestryStatGenderDiff%All%Max,  fmt=FMTREAL2CHAR))

          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") "  - coancestry among/between males"
          write(STDOUT, "(a)") "                     Among     Between"
          write(STDOUT, "(a)") "    - n:       "//trim( Int2Char(This%CoancestryStatGender1%All%n,    fmt=FMTINT2CHAR)) //" "//trim( Int2Char(This%CoancestryStatGender1%OffDiag%n,    fmt=FMTINT2CHAR))
          write(STDOUT, "(a)") "    - average: "//trim(Real2Char(This%CoancestryStatGender1%All%Mean, fmt=FMTREAL2CHAR))//" "//trim(Real2Char(This%CoancestryStatGender1%OffDiag%Mean, fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(This%CoancestryStatGender1%All%Sd,   fmt=FMTREAL2CHAR))//" "//trim(Real2Char(This%CoancestryStatGender1%OffDiag%Sd,   fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(This%CoancestryStatGender1%All%Min,  fmt=FMTREAL2CHAR))//" "//trim(Real2Char(This%CoancestryStatGender1%OffDiag%Min,  fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(This%CoancestryStatGender1%All%Max,  fmt=FMTREAL2CHAR))//" "//trim(Real2Char(This%CoancestryStatGender1%OffDiag%Max,  fmt=FMTREAL2CHAR))

          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") "  - coancestry among/between females"
          write(STDOUT, "(a)") "                     Among     Between"
          write(STDOUT, "(a)") "    - n:       "//trim( Int2Char(This%CoancestryStatGender2%All%n,    fmt=FMTINT2CHAR)) //" "//trim( Int2Char(This%CoancestryStatGender2%OffDiag%n,    fmt=FMTINT2CHAR))
          write(STDOUT, "(a)") "    - average: "//trim(Real2Char(This%CoancestryStatGender2%All%Mean, fmt=FMTREAL2CHAR))//" "//trim(Real2Char(This%CoancestryStatGender2%OffDiag%Mean, fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(This%CoancestryStatGender2%All%Sd,   fmt=FMTREAL2CHAR))//" "//trim(Real2Char(This%CoancestryStatGender2%OffDiag%Sd,   fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(This%CoancestryStatGender2%All%Min,  fmt=FMTREAL2CHAR))//" "//trim(Real2Char(This%CoancestryStatGender2%OffDiag%Min,  fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(This%CoancestryStatGender2%All%Max,  fmt=FMTREAL2CHAR))//" "//trim(Real2Char(This%CoancestryStatGender2%OffDiag%Max,  fmt=FMTREAL2CHAR))
        end if
      end if

      if (This%CoancestryStat%OffDiag%Sd .eq. 0.0) then
        write(STDERR, "(a)") " ERROR: There is no variation in coancestry between individuals!"
        write(STDERR, "(a)") " "
        stop 1
      end if

      ! Save means to a file
      open(newunit=CoancestrySummaryUnit, file=trim(Spec%OutputBasename)//"CoancestrySummary.txt", status="unknown")
      write(CoancestrySummaryUnit, "(a, f)") "Current (random mating/crossing),                 ",   This%CoancestryRanMate
      write(CoancestrySummaryUnit, "(a, f)") "Current (random mating/crossing, no selfing),     ",   This%CoancestryRanMateNoSelf
      if (Spec%GenderGiven) then
        write(CoancestrySummaryUnit, "(a, f)") "Current (random mating/crossing between genders), ", This%CoancestryGenderMate
      end if
      close(CoancestrySummaryUnit)

      ! --- Current inbreeding ---

      if (LogStdoutInternal) then
        write(STDOUT, "(a)") " "
        write(STDOUT, "(a)") " Current inbreeding (identity between the two genomes of an individual)"
      end if

      block
        type(InbVec) :: Inbreeding
        call This%Coancestry%Inbreeding(Out=Inbreeding, Nrm=.false.)
        This%InbreedingStat = DescStat(FLOATFUN(Inbreeding%Value(1:)))
        This%Inbreeding = This%InbreedingStat%Mean
      end block

      if (LogStdoutInternal) then
        write(STDOUT, "(a)") "  - n:       "//trim( Int2Char(This%InbreedingStat%n,    fmt=FMTINT2CHAR))
        write(STDOUT, "(a)") "  - average: "//trim(Real2Char(This%InbreedingStat%Mean, fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "  - st.dev.: "//trim(Real2Char(This%InbreedingStat%Sd,   fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "  - minimum: "//trim(Real2Char(This%InbreedingStat%Min,  fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "  - maximum: "//trim(Real2Char(This%InbreedingStat%Max,  fmt=FMTREAL2CHAR))
      end if

      open(newunit=InbreedingSummaryUnit, file=trim(Spec%OutputBasename)//"InbreedingSummary.txt", status="unknown")
      write(InbreedingSummaryUnit, "(a, f)") "Current, ", This%Inbreeding
      close(InbreedingSummaryUnit)

      ! --- Current selection criterion ---

      if (Spec%SelCriterionGiven) then

        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") " Current selection criterion / standardized"
        end if

        This%SelCriterionStat = DescStat(This%SelCriterion)
        This%SelCriterionStd = SelCriterion2SelCriterionStd(SelCriterion=This%SelCriterion, &
                                                            Mean=This%SelCriterionStat%Mean, &
                                                            Sd=This%SelCriterionStat%Sd)
        This%SelCriterionStdStat = DescStat(This%SelCriterionStd)

        if (LogStdoutInternal) then
          write(STDOUT, "(a)") "  - n:       "//trim( Int2Char(This%SelCriterionStat%n,       fmt=FMTINT2CHAR)) //" /"&
                                              //trim( Int2Char(This%SelCriterionStdStat%n,    fmt=FMTINT2CHAR))
          write(STDOUT, "(a)") "  - average: "//trim(Real2Char(This%SelCriterionStat%Mean,    fmt=FMTREAL2CHAR))//" /"&
                                              //trim(Real2Char(This%SelCriterionStdStat%Mean, fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "  - st.dev.: "//trim(Real2Char(This%SelCriterionStat%Sd,      fmt=FMTREAL2CHAR))//" /"&
                                              //trim(Real2Char(This%SelCriterionStdStat%Sd,   fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "  - minimum: "//trim(Real2Char(This%SelCriterionStat%Min,     fmt=FMTREAL2CHAR))//" /"&
                                              //trim(Real2Char(This%SelCriterionStdStat%Min,  fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "  - maximum: "//trim(Real2Char(This%SelCriterionStat%Max,     fmt=FMTREAL2CHAR))//" /"&
                                              //trim(Real2Char(This%SelCriterionStdStat%Max,  fmt=FMTREAL2CHAR))
        end if

        if (This%SelCriterionStat%Sd .eq. 0.0) then
          write(STDERR, "(a)") " ERROR: There is no variation in selection criterion!"
          write(STDERR, "(a)") " "
          stop 1
        end if

        open(newunit=CriterionSummaryUnit, file=trim(Spec%OutputBasename)//"SelCriterionSummary.txt", status="unknown")
        write(CriterionSummaryUnit, "(a, f)") "Mean, ", This%SelCriterionStat%Mean
        write(CriterionSummaryUnit, "(a, f)") "Sd, ", This%SelCriterionStat%Sd
        close(CriterionSummaryUnit)

        if (Spec%PAGEPar) then
          ! must have the same scale as selection criterion!!!!
          This%SelCriterionStdPAGE = SelCriterion2SelCriterionStd(SelCriterion=This%SelCriterionPAGE, &
                                                                  Mean=This%SelCriterionStat%Mean, &
                                                                  Sd=This%SelCriterionStat%Sd)
          ! only the PAGE bit of SelCriterion
          This%SelCriterionPAGE    = This%SelCriterionPAGE    - This%SelCriterion
          This%SelCriterionStdPAGE = This%SelCriterionStdPAGE - This%SelCriterionStd
          This%SelCriterionPAGEStat = DescStat(This%SelCriterionPAGE)
          if (LogStdoutInternal) then
            write(STDOUT, "(a)") " "
            write(STDOUT, "(a)") "Selection criterion increments with PAGE"
            write(STDOUT, "(a)") "  - n:       "//trim( Int2Char(This%SelCriterionPAGEStat%n,    fmt=FMTINT2CHAR))
            write(STDOUT, "(a)") "  - average: "//trim(Real2Char(This%SelCriterionPAGEStat%Mean, fmt=FMTREAL2CHAR))
            write(STDOUT, "(a)") "  - st.dev.: "//trim(Real2Char(This%SelCriterionPAGEStat%Sd,   fmt=FMTREAL2CHAR))
            write(STDOUT, "(a)") "  - minimum: "//trim(Real2Char(This%SelCriterionPAGEStat%Min,  fmt=FMTREAL2CHAR))
            write(STDOUT, "(a)") "  - maximum: "//trim(Real2Char(This%SelCriterionPAGEStat%Max,  fmt=FMTREAL2CHAR))
          end if

          if (This%SelCriterionPAGEStat%Sd .eq. 0.0) then
            write(STDERR, "(a)") " ERROR: There is no variation in selection criterion increments with PAGE!"
            write(STDERR, "(a)") " "
            stop 1
          end if

          open(newunit=CriterionSummaryUnit, file=trim(Spec%OutputBasename)//"PAGESummary.txt", status="unknown")
          write(CriterionSummaryUnit, "(a, f)") "Mean, ", This%SelCriterionPAGEStat%Mean
          write(CriterionSummaryUnit, "(a, f)") "Sd, ", This%SelCriterionPAGEStat%Sd
          close(CriterionSummaryUnit)
        end if

      end if

      ! --- Current generic individual values ---

      if (Spec%GenericIndCritGiven) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") " Current generic individual selection criterion"
        end if

        open(newunit=GenericIndCritSummaryUnit, file=trim(Spec%OutputBasename)//"GenericIndCritSummary.txt", status="unknown")

        allocate(This%GenericIndCritStat(Spec%nGenericIndCrit))
        do Crit = 1, Spec%nGenericIndCrit
          This%GenericIndCritStat(Crit) = DescStat(This%GenericIndCrit(:, Crit))
          if (LogStdoutInternal) then
            write(STDOUT, "(a)") " "
            write(STDOUT, "(a)") "  - criterion "//trim(Int2Char(Crit))
            write(STDOUT, "(a)") "    - n:       "//trim( Int2Char(This%GenericIndCritStat(Crit)%n,    fmt=FMTINT2CHAR))
            write(STDOUT, "(a)") "    - average: "//trim(Real2Char(This%GenericIndCritStat(Crit)%Mean, fmt=FMTREAL2CHAR))
            write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(This%GenericIndCritStat(Crit)%Sd,   fmt=FMTREAL2CHAR))
            write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(This%GenericIndCritStat(Crit)%Min,  fmt=FMTREAL2CHAR))
            write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(This%GenericIndCritStat(Crit)%Max,  fmt=FMTREAL2CHAR))
          end if
          if (This%GenericIndCritStat(Crit)%Sd .eq. 0.0) then
            write(STDERR, "(a)") " ERROR: There is no variation in generic individual selection criterion "//trim(Int2Char(Crit))//"!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          write(GenericIndCritSummaryUnit, "(a, f)") "Mean criterion "//trim(Int2Char(Crit))//",", This%GenericIndCritStat(Crit)%Mean
          write(GenericIndCritSummaryUnit, "(a, f)") "Sd criterion "//trim(Int2Char(Crit))//",", This%GenericIndCritStat(Crit)%Sd
        end do

        close(GenericIndCritSummaryUnit)
      end if

      ! --- Generic mating values ---

      if (Spec%GenericMatCritGiven) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") " Generic mating/crossing selection criterion"
        end if

        open(newunit=GenericMatCritSummaryUnit, file=trim(Spec%OutputBasename)//"GenericMatCritSummary.txt", status="unknown")

        allocate(This%GenericMatCritStat(Spec%nGenericMatCrit))
        do Crit = 1, Spec%nGenericMatCrit
          if (LogStdoutInternal) then
            write(STDOUT, "(a)") " "
            write(STDOUT, "(a)") "  - criterion "//trim(Int2Char(Crit))
          end if
          if (Spec%GenderGiven) then
            This%GenericMatCritStat(Crit) = DescStatMatrix(This%GenericMatCrit(:, :, Crit))
            if (LogStdoutInternal) then
              write(STDOUT, "(a)") "    - n:       "//trim( Int2Char(This%GenericMatCritStat(Crit)%All%n,    fmt=FMTINT2CHAR))
              write(STDOUT, "(a)") "    - average: "//trim(Real2Char(This%GenericMatCritStat(Crit)%All%Mean, fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(This%GenericMatCritStat(Crit)%All%Sd,   fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(This%GenericMatCritStat(Crit)%All%Min,  fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(This%GenericMatCritStat(Crit)%All%Max,  fmt=FMTREAL2CHAR))
            end if
            if (This%GenericMatCritStat(Crit)%All%Sd .eq. 0.0) then
              write(STDERR, "(a)") " ERROR: There is no variation in generic mating/crossing selection criterion "//trim(Int2Char(Crit))//"!"
              write(STDERR, "(a)") " "
              stop 1
            end if
            write(GenericMatCritSummaryUnit, "(a, f)") "Mean criterion "//trim(Int2Char(Crit))//",", This%GenericMatCritStat(Crit)%All%Mean
            write(GenericMatCritSummaryUnit, "(a, f)") "Sd criterion "//trim(Int2Char(Crit))//",", This%GenericMatCritStat(Crit)%All%Sd
          else
            if (Spec%SelfingAllowed) then
              This%GenericMatCritStat(Crit) = DescStatLowTriMatrix(This%GenericMatCrit(:, :, Crit))
              if (LogStdoutInternal) then
                write(STDOUT, "(a)") "    - n:       "//trim( Int2Char(This%GenericMatCritStat(Crit)%All%n,    fmt=FMTINT2CHAR))
                write(STDOUT, "(a)") "    - average: "//trim(Real2Char(This%GenericMatCritStat(Crit)%All%Mean, fmt=FMTREAL2CHAR))
                write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(This%GenericMatCritStat(Crit)%All%Sd,   fmt=FMTREAL2CHAR))
                write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(This%GenericMatCritStat(Crit)%All%Min,  fmt=FMTREAL2CHAR))
                write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(This%GenericMatCritStat(Crit)%All%Max,  fmt=FMTREAL2CHAR))
              end if
              if (This%GenericMatCritStat(Crit)%All%Sd .eq. 0.0) then
                write(STDERR, "(a)") " ERROR: There is no variation in generic mating/crossing selection criterion "//trim(Int2Char(Crit))//"!"
                write(STDERR, "(a)") " "
                stop 1
              end if
              write(GenericMatCritSummaryUnit, "(a, f)") "Mean criterion "//trim(Int2Char(Crit))//",", This%GenericMatCritStat(Crit)%All%Mean
              write(GenericMatCritSummaryUnit, "(a, f)") "Sd criterion "//trim(Int2Char(Crit))//",", This%GenericMatCritStat(Crit)%All%Sd
            end if
              This%GenericMatCritStat(Crit) = DescStatLowTriMatrix(This%GenericMatCrit(:, :, Crit), Diag=.false.)
              if (LogStdoutInternal) then
                write(STDOUT, "(a)") "    - n:       "//trim( Int2Char(This%GenericMatCritStat(Crit)%OffDiag%n,    fmt=FMTINT2CHAR))
                write(STDOUT, "(a)") "    - average: "//trim(Real2Char(This%GenericMatCritStat(Crit)%OffDiag%Mean, fmt=FMTREAL2CHAR))
                write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(This%GenericMatCritStat(Crit)%OffDiag%Sd,   fmt=FMTREAL2CHAR))
                write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(This%GenericMatCritStat(Crit)%OffDiag%Min,  fmt=FMTREAL2CHAR))
                write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(This%GenericMatCritStat(Crit)%OffDiag%Max,  fmt=FMTREAL2CHAR))
              end if
              if (This%GenericMatCritStat(Crit)%OffDiag%Sd .eq. 0.0) then
                write(STDERR, "(a)") " ERROR: There is no variation in generic mating/crossing selection criterion "//trim(Int2Char(Crit))//"!"
                write(STDERR, "(a)") " "
                stop 1
              end if
              write(GenericMatCritSummaryUnit, "(a, f)") "Mean criterion "//trim(Int2Char(Crit))//",", This%GenericMatCritStat(Crit)%OffDiag%Mean
              write(GenericMatCritSummaryUnit, "(a, f)") "Sd criterion "//trim(Int2Char(Crit))//",", This%GenericMatCritStat(Crit)%OffDiag%Sd
          end if
        end do

        close(GenericMatCritSummaryUnit)
      end if

    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Write AlphaMate data to a file or standard output
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   April 4, 2017
    !> @return Output to a file or standard output
    !---------------------------------------------------------------------------
    subroutine WriteAlphaMateData(This, File) ! not pure due to IO
      implicit none
      class(AlphaMateData), intent(in)       :: This !< AlphaMateData holder
      character(len=*), intent(in), optional :: File !< File (if missing use standard output)

      integer(int32) :: Unit, Ind
      if (present(File)) then
        open(newunit=Unit, file=File, action="write", status="unknown")
      else
        Unit = STDOUT
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " --- Raw data ---"

      write(Unit, "(a)") " "
      write(Unit, "(a)") " Coancestry:"
      call This%Coancestry%Write(File=File)

      write(Unit, "(a)") " "
      write(Unit, "(a)") " Gender:"
      if (allocated(This%Gender)) then
        do Ind = 1, This%nInd
          write(Unit, "(a, i)") This%Coancestry%OriginalId(Ind), This%Gender(Ind)
        end do
      else
        write(Unit, "(a)") " not allocated"
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " SelCriterion:"
      if (allocated(This%SelCriterion)) then
        do Ind = 1, This%nInd
          write(Unit, "(a, 2f)") This%Coancestry%OriginalId(Ind), This%SelCriterion(Ind), This%SelCriterionStd(Ind)
        end do
      else
        write(Unit, "(a)") " not allocated"
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " SelCriterionPAGE:"
      if (allocated(This%SelCriterionPAGE)) then
        do Ind = 1, This%nInd
          write(Unit, "(a, f)") This%Coancestry%OriginalId(Ind), This%SelCriterionPAGE(Ind), This%SelCriterionStdPAGE(Ind)
        end do
      else
        write(Unit, "(a)") " not allocated"
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " GenericIndCrit:"
      if (allocated(This%GenericIndCrit)) then
        do Ind = 1, This%nInd
          write(Unit, "(a, f)") This%Coancestry%OriginalId(Ind), This%GenericIndCrit(:, Ind)
        end do
      else
        write(Unit, "(a)") " not allocated"
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " GenericMatCrit:"
      if (allocated(This%GenericMatCrit)) then
        write(Unit, "(a)") " @todo"
        ! do Mat = 1, This%nPotMat
        !   write(Unit, "(i, f)") Mat, This%GenericMatCrit(:, :)
        ! end do
      else
        write(Unit, "(a)") " not allocated"
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " --- Data summaries ---"

      write(Unit, "(a)") " "
      write(Unit, "(a)") " CoancestryStat:"
      write(Unit, *) This%CoancestryStat

      write(Unit, "(a)") " "
      write(Unit, "(a)") " CoancestryStatGender1:"
      if (allocated(This%Gender)) then
        write(Unit, *) This%CoancestryStatGender1
      else
        write(Unit, "(a)") " Gender not allocated (so no CoancestryStatGender1)"
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " CoancestryStatGender2:"
      if (allocated(This%Gender)) then
        write(Unit, *) This%CoancestryStatGender2
      else
        write(Unit, "(a)") " Gender not allocated (so no CoancestryStatGender2)"
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " CoancestryStatGenderDiff:"
      if (allocated(This%Gender)) then
        write(Unit, *) This%CoancestryStatGenderDiff
      else
        write(Unit, "(a)") " Gender not allocated (so no CoancestryStatGenderDiff)"
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " InbreedingStat:"
      write(Unit, *) This%InbreedingStat

      write(Unit, "(a)") " "
      write(Unit, "(a)") " SelCriterionStat:"
      if (allocated(This%SelCriterion)) then
        write(Unit, *) This%SelCriterionStat
      else
        write(Unit, "(a)") " SelCriterion not allocated"
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " SelCriterionStdStat:"
      if (allocated(This%SelCriterion)) then
        write(Unit, *) This%SelCriterionStdStat
      else
        write(Unit, "(a)") " SelCriterion not allocated"
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " SelCriterionPAGEStat:"
      if (allocated(This%SelCriterionPAGE)) then
        write(Unit, *) This%SelCriterionPAGEStat
      else
        write(Unit, "(a)") " SelCriterionPAGE not allocated"
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " GenericIndCritStat:"
      if (allocated(This%GenericIndCrit)) then
        write(Unit, *) This%GenericIndCritStat
      else
        write(Unit, "(a)") " GenericIndCrit not allocated"
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " GenericMatCritStat:"
      if (allocated(This%GenericMatCrit)) then
        write(Unit, *) This%GenericMatCritStat
      else
        write(Unit, "(a)") " GenericMatCrit not allocated"
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " --- Derived data ---"

      write(Unit, "(a)") " "
      write(Unit, "(a)") " nInd:"
      write(Unit, "(i)") This%nInd

      write(Unit, "(a)") " "
      write(Unit, "(a)") " nPotMat:"
      write(Unit, "(i)") This%nPotMat

      write(Unit, "(a)") " "
      write(Unit, "(a)") " nPotPar1:"
      write(Unit, "(i)") This%nPotPar1

      write(Unit, "(a)") " "
      write(Unit, "(a)") " nPotPar2:"
      write(Unit, "(i)") This%nPotPar2

      write(Unit, "(a)") " "
      write(Unit, "(a)") " nPotPar:"
      write(Unit, "(i)") This%nPotPar

      write(Unit, "(a)") " "
      write(Unit, "(a)") " nMal:"
      if (allocated(This%Gender)) then
        write(Unit, "(i)") This%nMal
      else
        write(Unit, "(a)") " Gender not allocated (so no nMal)"
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " nFem:"
      if (allocated(This%Gender)) then
        write(Unit, "(i)") This%nFem
      else
        write(Unit, "(a)") " Gender not allocated (so no nFem)"
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " IdPotPar1:"
      write(Unit, "(i)") This%IdPotPar1

      write(Unit, "(a)") " "
      write(Unit, "(a)") " IdPotPar2:"
      if (allocated(This%Gender)) then
        write(Unit, "(i)") This%IdPotPar2
      else
        write(Unit, "(a)") " Gender not allocated (so no IdPotPar2)"
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " IdPotParSeq:"
      if (allocated(This%Gender)) then
        write(Unit, "(i)") This%IdPotParSeq
      else
        write(Unit, "(a)") " Gender not allocated (so no IdPotParSeq)"
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " CoancestryRanMate:"
      write(Unit, "(f)") This%CoancestryRanMate

      write(Unit, "(a)") " "
      write(Unit, "(a)") " CoancestryRanMateNoSelf:"
      write(Unit, "(f)") This%CoancestryRanMateNoSelf

      write(Unit, "(a)") " "
      write(Unit, "(a)") " CoancestryGenderMate:"
      if (allocated(This%Gender)) then
        write(Unit, "(f)") This%CoancestryGenderMate
      else
        write(Unit, "(a)") " Gender not allocated (so no CoancestryGenderMate)"
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " Inbreeding:"
      write(Unit, "(f)") This%Inbreeding

      if (present(File)) then
        close(Unit)
      end if
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Setup AlphaMate system (seeed, random number generator, and paralelisation)
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   July 14, 2018
    !> @return Seed, random number generator, and paralelisation setup
    !---------------------------------------------------------------------------
    subroutine AlphaMateSystem(Spec, LogStdout) ! not pure due to IO and RNG
      implicit none
      type(AlphaMateSpec), intent(inout) :: Spec !< AlphaMateSpec holder (inout because we save some info in Spec that
      logical, optional                  :: LogStdout !< Log process on stdout (default .false.)

      integer(int32) :: Unit, nProcs, EnvVarStatus

      character(len = 4) :: EnvVarOmp_Num_Threads

      logical :: LogStdoutInternal
      if (present(LogStdout)) then
        LogStdoutInternal = LogStdout
      else
        LogStdoutInternal = .false.
      end if

      if (LogStdoutInternal) then
        write(STDOUT, "(a)") " "
      end if

      ! --- System/Compiler Seed and Random Number Generator ---

      ! This part gets or sets seed for system/compiler random number generator. While AlphaMate
      ! uses Intel random number generator, it is handy to get a "random" seed, when seed is not
      ! specified by the user. Further, setting system/compiler seed enforces reproducibility in
      ! case random_number() would be used anywhere in the code.
      if (.not. (Spec%SeedGiven .or. Spec%SeedFileGiven)) then
        call SetSeed(Out=Spec%Seed) ! get a random seed from the system/compiler
        write(STDOUT, "(a)") " RNG seed (system): "//trim(Int2Char(Spec%Seed))
      else
        if ((.not. Spec%SeedGiven) .and. Spec%SeedFileGiven) then
          open(newunit=Unit, file=Spec%SeedFile, status="old")
          read(Unit, *) Spec%Seed
          close(Unit)
        end if
        call SetSeed(Seed=Spec%Seed) ! set user specified seed
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " Specified RNG seed: "//trim(Int2Char(Spec%Seed))
        end if
      end if
      open(newunit=Unit, file=trim(Spec%OutputBasename)//"SeedUsed.txt", status="unknown")
      write(Unit, *) Spec%Seed
      close(Unit)

      ! --- Intel Random Number Generator ---

      call IntitialiseIntelRng(Seed=Spec%Seed, Stream=AlphaMateStream)

      ! --- Processors & Threads ---

      nProcs = Omp_Get_Num_Procs()
      if (LogStdoutInternal) then
        write(STDOUT, "(a)") " Number of system processors: "//trim(Int2Char(nProcs))
      end if
      call Get_Environment_Variable(name="OMP_NUM_THREADS", value=EnvVarOmp_Num_Threads, status=EnvVarStatus)
      if (EnvVarStatus .eq. 0) then
        write(STDOUT, "(a)") " Number of available threads (via the OMP_NUM_THREADS environment variable): "//trim(EnvVarOmp_Num_Threads)
      end if
      if (Spec%nThreadsGiven) then
        if (LogStdoutInternal) then
          if (EnvVarStatus .eq. 0) then
            write(STDOUT, "(a)") " Number of specified threads (overrules the OMP_NUM_THREADS environment variable): "//trim(Int2Char(Spec%nThreads))
          else
            write(STDOUT, "(a)") " Number of specified threads: "//trim(Int2Char(Spec%nThreads))
          end if
        end if
        call Omp_Set_Num_Threads(Spec%nThreads)
      else
      end if

    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Initialise AlphaMate solution
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !---------------------------------------------------------------------------
    pure subroutine InitialiseAlphaMateSol(This, Chrom, Spec)
      implicit none

      ! Argument
      class(AlphaMateSol), intent(out)             :: This     !< @return AlphaMateSol holder
      real(FLOATTYPEAH), intent(in)                :: Chrom(:) !< Provided initial solution
      class(AlphaEvolveSpec), intent(in), optional :: Spec     !< AlphaEvolveSpec --> AlphaMateSpec holder

      real(FLOATTYPE) :: NANFLOATTYPE
      NANFLOATTYPE = IEEE_Value(x=NANFLOATTYPE, class=IEEE_Quiet_NaN)

      ! Initialisation
      select type (Spec)
        class default
          ! @todo This will work with Fortran 2015 standard (or at least with new? ifort) that allows error stop in PURE subroutines
          ! error stop " ERROR: InitialiseAlphaMateSol works only with argument Spec being of type AlphaMateSpec!"
        class is (AlphaMateSpec)
          This%Objective = -huge(This%Objective)
          This%nParam = size(Chrom)
          allocate(This%Chrom(This%nParam))
          This%Chrom = Chrom
          This%Penalty = 0.0
          This%PenaltyCoancestry = 0.0
          This%PenaltyInbreeding = 0.0
          This%PenaltyLimitPar1 = 0.0
          This%PenaltyLimitPar2 = 0.0
          This%PenaltyMultipleMatings = 0.0
          This%PenaltySelfing = 0.0
          This%PenaltyGenericIndCrit = 0.0
          This%PenaltyGenericMatCrit = 0.0
          This%Degree = NANFLOATTYPE
          This%SelCriterion = 0.0
          This%SelCriterionStd = 0.0
          This%MaxCriterionPct = NANFLOATTYPE
          This%CoancestryRanMate = 0.0
          This%CoancestryRateRanMate = 0.0
          This%MinCoancestryPct = NANFLOATTYPE
          This%Inbreeding = 0.0
          This%InbreedingRate = 0.0
          This%MinInbreedingPct = NANFLOATTYPE
          if (Spec%GenericIndCritGiven) then
            allocate(This%GenericIndCrit(Spec%nGenericIndCrit))
            This%GenericIndCrit = 0.0
          end if
          if (Spec%GenericMatCritGiven) then
            allocate(This%GenericMatCrit(Spec%nGenericMatCrit))
            This%GenericMatCrit = 0.0
          end if
          ! This%Cost = 0.0
          allocate(This%nVec(Spec%nInd))
          This%nVec = 0
          allocate(This%MatingPlan(2, Spec%nMat))
          This%MatingPlan = 0
          if (Spec%PAGEPar) then
            allocate(This%GenomeEdit(Spec%nInd))
            This%GenomeEdit = 0.0
          end if
      end select
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Assign one AlphaMate solution to another
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !---------------------------------------------------------------------------
    pure subroutine AssignAlphaMateSol(Out, In)
      implicit none

      ! Arguments
      class(AlphaMateSol), intent(out)   :: Out !< @return AlphaMateSol holder
      class(AlphaEvolveSol), intent(in)  :: In  !< AlphaEvolveSol --> AlphaMateSol holder

      ! Assignments
      select type (In)
        class default
          ! @todo This will work with Fortran 2015 standard (or at least with new? ifort) that allows error stop in PURE subroutines
          ! error stop " ERROR: AssignAlphaMateSol works only with argument In being of type AlphaMateSol!"
        class is (AlphaMateSol)
          Out%Objective = In%Objective
          Out%nParam = In%nParam
          Out%Chrom = In%Chrom
          Out%Penalty = In%Penalty
          Out%PenaltyCoancestry = In%PenaltyCoancestry
          Out%PenaltyInbreeding = In%PenaltyInbreeding
          Out%PenaltyLimitPar1 = In%PenaltyLimitPar1
          Out%PenaltyLimitPar2 = In%PenaltyLimitPar2
          Out%PenaltyMultipleMatings = In%PenaltyMultipleMatings
          Out%PenaltySelfing = In%PenaltySelfing
          Out%PenaltyGenericIndCrit = In%PenaltyGenericIndCrit
          Out%PenaltyGenericMatCrit = In%PenaltyGenericMatCrit
          Out%Degree = In%Degree
          Out%SelCriterion = In%SelCriterion
          Out%SelCriterionStd = In%SelCriterionStd
          Out%MaxCriterionPct = In%MaxCriterionPct
          Out%CoancestryRanMate = In%CoancestryRanMate
          Out%CoancestryRateRanMate = In%CoancestryRateRanMate
          Out%MinCoancestryPct = In%MinCoancestryPct
          Out%Inbreeding = In%Inbreeding
          Out%InbreedingRate = In%InbreedingRate
          Out%MinInbreedingPct = In%MinInbreedingPct
          if (allocated(In%GenericIndCrit)) then
            allocate(Out%GenericIndCrit(size(In%GenericIndCrit)))
            Out%GenericIndCrit = In%GenericIndCrit
          end if
          if (allocated(In%GenericMatCrit)) then
            allocate(Out%GenericMatCrit(size(In%GenericMatCrit)))
            Out%GenericMatCrit = In%GenericMatCrit
          end if
          ! Out%Cost = In%Cost
          allocate(Out%nVec(size(In%nVec)))
          Out%nVec = In%nVec
          allocate(Out%MatingPlan(size(In%MatingPlan, dim=1), size(In%MatingPlan, dim=2)))
          Out%MatingPlan = In%MatingPlan
          if (allocated(In%GenomeEdit)) then
            allocate(Out%GenomeEdit(size(In%GenomeEdit)))
            Out%GenomeEdit = In%GenomeEdit
          end if
      end select
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Update mean of AlphaMate solution (when performing random search)
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !---------------------------------------------------------------------------
    pure subroutine UpdateMeanAlphaMateSol(This, Add, n)
      implicit none

      ! Arguments
      class(AlphaMateSol), intent(inout) :: This !< @return AlphaMateSol holder
      class(AlphaEvolveSol), intent(in)  :: Add  !< AlphaEvolveSol --> AlphaMateSol holder
      integer(int32), intent(in)         :: n    !< Number of solutions averaged togehter

      ! Other
      real(FLOATTYPE) :: kR

      ! Updates
      kR = (FLOATFUN(n) - 1.0) / n

      ! (Need to go via the select type stuff as all but the first arguments must
      !  be the same as in the base class/type)
      select type (Add)
        class default
          ! @todo This will work with Fortran 2015 standard (or at least with new? ifort) that allows error stop in PURE subroutines
          ! error stop " ERROR: UpdateMeanAlphaMateSol works only with argument Add being of type AlphaMateSol!"
        class is (AlphaMateSol)
          This%Objective                 = This%Objective                 * kR + Add%Objective                 / n
          ! This%nParam                    = This%nParam                    * kR + Add%nParam                    / n ! the same all the time
          ! This%Chrom                     = This%Chrom                     * kR + Add%Chrom                     / n ! hmm, do we really want to average over chromosomes?
          This%Penalty                   = This%Penalty                   * kR + Add%Penalty                   / n
          This%PenaltyCoancestry         = This%PenaltyCoancestry         * kR + Add%PenaltyCoancestry         / n
          This%PenaltyInbreeding         = This%PenaltyInbreeding         * kR + Add%PenaltyInbreeding         / n
          This%PenaltyLimitPar1          = This%PenaltyLimitPar1          * kR + Add%PenaltyLimitPar1          / n
          This%PenaltyLimitPar2          = This%PenaltyLimitPar2          * kR + Add%PenaltyLimitPar2          / n
          This%PenaltyMultipleMatings    = This%PenaltyMultipleMatings    * kR + Add%PenaltyMultipleMatings    / n
          This%PenaltySelfing            = This%PenaltySelfing            * kR + Add%PenaltySelfing            / n
          This%PenaltyGenericIndCrit     = This%PenaltyGenericIndCrit     * kR + Add%PenaltyGenericIndCrit     / n
          This%PenaltyGenericMatCrit     = This%PenaltyGenericMatCrit     * kR + Add%PenaltyGenericMatCrit     / n
          This%Degree                    = This%Degree                    * kR + Add%Degree                    / n
          This%SelCriterion              = This%SelCriterion              * kR + Add%SelCriterion              / n
          This%SelCriterionStd           = This%SelCriterionStd           * kR + Add%SelCriterionStd           / n
          This%MaxCriterionPct           = This%MaxCriterionPct           * kR + Add%MaxCriterionPct           / n
          This%CoancestryRanMate         = This%CoancestryRanMate         * kR + Add%CoancestryRanMate         / n
          This%CoancestryRateRanMate     = This%CoancestryRateRanMate     * kR + Add%CoancestryRateRanMate     / n
          This%MinCoancestryPct          = This%MinCoancestryPct          * kR + Add%MinCoancestryPct          / n
          This%Inbreeding                = This%Inbreeding                * kR + Add%Inbreeding                / n
          This%InbreedingRate            = This%InbreedingRate            * kR + Add%InbreedingRate            / n
          This%MinInbreedingPct          = This%MinInbreedingPct          * kR + Add%MinInbreedingPct          / n
          if (allocated(This%GenericIndCrit)) then
            This%GenericIndCrit          = This%GenericIndCrit            * kR + Add%GenericIndCrit            / n
          end if
          if (allocated(This%GenericMatCrit)) then
            This%GenericMatCrit          = This%GenericMatCrit            * kR + Add%GenericMatCrit            / n
          end if
          ! This%Cost                      = This%Cost                      * kR + Add%Cost                      / n
          This%nVec                      = This%nVec                      * kR + Add%nVec                      / n
          ! @todo It does not make sense to average a mating plan, no?
          This%MatingPlan                = This%MatingPlan                * kR + Add%MatingPlan                / n
          if (allocated(This%GenomeEdit)) then
            This%GenomeEdit              = This%GenomeEdit                * kR + Add%GenomeEdit                / n
          end if
          ! Constant
          ! This%TargetDegree
          ! This%TargetSelCriterion
          ! This%TargetSelCriterionStd
          ! This%TargetMaxCriterionPct
          ! This%TargetCoancestry
          ! This%TargetCoancestryRate
          ! This%TargetMinCoancestryPct
          ! This%TargetInbreeding
          ! This%TargetInbreedingRate
      end select
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Write AlphaMate chromosome to a file or standard output
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   July 14, 2017
    !> @return Output to a file or standard output
    !---------------------------------------------------------------------------
    subroutine WriteAlphaMateChrom(This, File) ! not pure due to IO
      implicit none
      class(AlphaMateChrom), intent(in)      :: This !< AlphaMateChrom holder
      character(len=*), intent(in), optional :: File !< File (if missing use standard output)

      integer(int32) :: Unit, Param
      if (present(File)) then
        open(newunit=Unit, file=File, action="write", status="unknown")
      else
        Unit = STDOUT
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " ContPar1:"
      do Param = 1, size(This%ContPar1)
        write(Unit, "(i0, f)") Param, This%ContPar1(Param)
      end do

      write(Unit, "(a)") " "
      write(Unit, "(a)") " ContPar2:"
      if (allocated(This%ContPar2)) then
        do Param = 1, size(This%ContPar2)
          write(Unit, "(i0, f)") Param, This%ContPar2(Param)
        end do
      else
        write(Unit, "(a)") " not allocated"
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " MateRank:"
      if (allocated(This%MateRank)) then
        do Param = 1, size(This%MateRank)
          write(Unit, "(i0, f)") Param, This%MateRank(Param)
        end do
      else
        write(Unit, "(a)") " not allocated"
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " EditPar1:"
      if (allocated(This%EditPar1)) then
        do Param = 1, size(This%EditPar1)
          write(Unit, "(i0, f)") Param, This%EditPar1(Param)
        end do
      else
        write(Unit, "(a)") " not allocated"
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " EditPar2:"
      if (allocated(This%EditPar2)) then
        do Param = 1, size(This%EditPar2)
          write(Unit, "(i0, f)") Param, This%EditPar2(Param)
        end do
      else
        write(Unit, "(a)") " not allocated"
      end if

      if (present(File)) then
        close(Unit)
      end if
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  AlphaMate evaluate function plus much MORE (this is the core!!!!)
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !---------------------------------------------------------------------------
    subroutine FixSolEtcMateAndEvaluateAlphaMateSol(This, Chrom, Spec, Data, Stream) ! not pure due to RNG
      implicit none
      ! Arguments
      class(AlphaMateSol), intent(inout)           :: This     !< @return AlphaMateSol holder (out because we sometimes need to fix a solution)
      real(FLOATTYPEAH), intent(in)                :: Chrom(:) !< A solution
      class(AlphaEvolveSpec), intent(in)           :: Spec     !< AlphaEvolveSpec --> AlphaMateSpec holder
      class(AlphaEvolveData), intent(in), optional :: Data     !< AlphaEvolveData --> AlphaMateData holder
      type(vsl_stream_state), intent(inout)        :: Stream   !< Intel RNG stream

      ! Other
      integer(int32) :: i, j, k, l, GenderMode, Start, End, nCumMat, TmpMin, TmpMax, &
                        TmpI, nRanNum, RanNumLoc, Par1, Par2, TmpMate(2)
      integer(int32), allocatable :: Rank(:), MatPar2(:), nVecPar1(:), Pair(:, :)

      real(FLOATTYPE) :: TmpR, Diff, MaxDiff
      real(FLOATTYPE), allocatable :: TmpVec(:), RanNum(:) !, TmpVec2(:,:)

      type(AlphaMateChrom) :: SChrom

      select type (Spec)
        class default
          error stop " ERROR: FixSolEtcMateAndEvaluate works only with argument Spec being of type AlphaMateSpec!"
        class is (AlphaMateSpec)
          select type (Data)
            class default
              error stop " ERROR: FixSolEtcMateAndEvaluate works only with argument Data being of type AlphaMateData!"
            class is (AlphaMateData)

              ! --- Initialise ---

              ! Solution
              call This%Initialise(Chrom=Chrom, Spec=Spec)
              This%Objective = 0.0

              ! A solution (based on the mate selection driver) has:
              ! - Data%nInd individual contributions
              !   - Data%nPotPar1 individual contributions for "parent1" (males   when GenderGiven, all ind when .not. GenderGiven)
              !   - Data%nPotPar2 individual contributions for "parent2" (females when GenderGiven, meaningful only when GenderGiven)
              ! - Spec%nMat rankings of parent1 1:Spec%nMat matings to pair with 1:Data%nPotPar2 "parent2" (see below)
              ! - Data%nInd edit indicators
              !   - Data%nPotPar1 edit indicators for "parent1" (males   when GenderGiven, all ind when .not. GenderGiven)
              !   - Data%nPotPar2 edit indicators for "parent2" (females when GenderGiven, present only when GenderGiven)
              !
              ! Say we have Chrom=(| 0, 2, 0, 1 | 1, 2, 0 | 1.5, 2.5, 1.0 | 0, 1, 0, 0 | 0, 0, 0) then we:
              ! - sort female contributions based on ranks from (1st, 2nd, 2nd) female to (2nd, 1st, 2nd) female
              ! - then we
              !   - mate first  male 2 contribution with the first  female 2 contribution
              !   - mate second male 2 contribution with the first  female 1 contribution
              !   - mate first  male 4 contribution with the second female 2 contribution
              ! - edit male 2
              ! - do not edit any of the females

              ! Below we create a structured chromosome to simplify the evaluation code

              ! Allocate, Assign, and Setup some constants
              allocate(SChrom%ContPar1(Data%nPotPar1))
              Start = 1
              End = Data%nPotPar1
              SChrom%ContPar1 = Chrom(Start:End)
              if (Spec%GenderGiven) then
                allocate(SChrom%ContPar2(Data%nPotPar2))
                Start = End + 1
                End = Start - 1 + Data%nPotPar2
                SChrom%ContPar2 = Chrom(Start:End)
                GenderMode = 1
                nRanNum = GenderMode * Spec%nMat ! see nCumMat
              else
                GenderMode = 2
                nRanNum = GenderMode * Spec%nMat ! see nCumMat
              end if
              if (Spec%MateAllocation) then
                allocate(SChrom%MateRank(Spec%nMat))
                Start = End + 1
                End = Start - 1 + Spec%nMat
                SChrom%MateRank = Chrom(Start:End)
                if (.not. Spec%GenderGiven) then
                  nRanNum = nRanNum + Data%nPotPar1 ! see "Distribute one half of contributions into matings"
                end if
              end if
              if (Spec%PAGEPar) then
                if (Spec%PAGEPar1) then
                  allocate(SChrom%EditPar1(Data%nPotPar1))
                  Start = End + 1
                  End = Start - 1 + Data%nPotPar1
                  SChrom%EditPar1 = Chrom(Start:End)
                end if
                if (Spec%GenderGiven) then
                  if (Spec%PAGEPar2) then
                    allocate(SChrom%EditPar2(Data%nPotPar2))
                    Start = End + 1
                    End = Start - 1 + Data%nPotPar2
                    SChrom%EditPar2 = Chrom(Start:End)
                  end if
                end if
              end if

              ! do k = 1, size(Chrom)
              !   write(STDOUT, "(i3, f)") k, Chrom(k)
              ! end do
              ! call SChrom%Write

              ! Working vectors
              allocate(Rank(maxval([Data%nInd, Spec%nMat]))) ! for ranking many things
              allocate(nVecPar1(Data%nPotPar1))              ! for nVec
              allocate(MatPar2(Spec%nMat))                   ! for MatingPlan
              allocate(TmpVec(Data%nInd))                    ! for many things
              allocate(RanNum(nRanNum))                      ! for stochastic decisions
              RanNumLoc = 0
              RanNum = SampleIntelUniformD(n=nRanNum, Stream=Stream)

              ! --- Parse the mate selection driver (=Is the solution valid?) ---

              ! The approach below assures that we have Spec%nMat contributions for each of
              ! the two parent sets. It does this by ranking internal solution values and
              ! traverses from top to the defined number of parents checking when the sum of
              ! interegrised values gives Spec%nMat. If this still does not give Spec%nMat,
              ! then we start adding contributions to a random set of parent until Spec%nMat
              ! is reached. How to treat values for the individuals that do not contribute
              ! is unlcear. None of the tested methods seemed to be very different. Intuitively,
              ! using ordered negative values should inform optimisation which individuals
              ! should less likely contribute, but this did not seem to be the case - better final
              ! solution was found when this strategy was not implemented - either zeroing
              ! values for those individuals was the fastest or giving random value was
              ! marginally better, but slower. Potential advantage of not preserving the
              ! order is that this gives more randomness and more solutions being explored.

              ! "Parent1"
              ! ... preselect
              if (Spec%PreselectPar1) then
                ! if (Spec%ModeSpec%ObjectiveCoancestry) then
                !   Rank(1:Spec%PreselectPar1N) = RapKnr(-Data%AvgCoancestry(Data%IdPotPar1), Spec%PreselectPar1N) ! preselect contributors
                !   TmpVec(1:Spec%PreselectPar1N) = SChrom%ContPar1(Rank(1:Spec%PreselectPar1N))                   ! save contributors
                !   SChrom%ContPar1 = 0.0                                                                        ! set everyones contributions to zero
                !   SChrom%ContPar1(Rank(1:Spec%PreselectPar1N)) = TmpVec(1:Spec%PreselectPar1N)                   ! put contributors back
                ! end if
                if (Spec%ModeSpec%ObjectiveCriterion) then
                  Rank(1:Spec%PreselectPar1N) = RapKnr(Data%SelCriterionStd(Data%IdPotPar1),  Spec%PreselectPar1N) ! preselect contributors
                  TmpVec(1:Spec%PreselectPar1N) = SChrom%ContPar1(Rank(1:Spec%PreselectPar1N))                     ! save contributors
                  SChrom%ContPar1 = 0.0                                                                          ! set everyones contributions to zero
                  SChrom%ContPar1(Rank(1:Spec%PreselectPar1N)) = TmpVec(1:Spec%PreselectPar1N)                     ! put contributors back
                end if
              end if
              ! ... ranks to find contributors
              if (.not. (Spec%EqualizePar1 .and. (Spec%nPar1 .eq. Data%nPotPar1))) then
                Rank(1:Spec%nPar1) = RapKnr(SChrom%ContPar1, Spec%nPar1)
              end if
              if (Spec%EqualizePar1) then ! ... equal contributions
                if (Spec%nPar1 .eq. Data%nPotPar1) then
                  SChrom%ContPar1 = FLOATFUN(Spec%nMat * GenderMode) / Spec%nPar1 ! no need for indexing here, hence the above if (.not. ...)
                else
                  SChrom%ContPar1 = 0.0
                  SChrom%ContPar1(Rank(1:Spec%nPar1)) = FLOATFUN(Spec%nMat * GenderMode) / Spec%nPar1
                end if
              else                        ! ... unequal contributions
                TmpVec(1:Spec%nPar1) = SChrom%ContPar1(Rank(1:Spec%nPar1)) ! save top contributions
                SChrom%ContPar1 = 0.0                                    ! set everyones contributions to zero
                SChrom%ContPar1(Rank(1:Spec%nPar1)) = TmpVec(1:Spec%nPar1) ! put top contributions back
                nCumMat = 0
                do i = 1, Spec%nPar1
                  j = Rank(i)
                  ! .. cap minimum usage @todo could consider penalising solution instead?
                  if (SChrom%ContPar1(j) .lt. Spec%LimitPar1Min) then
                    SChrom%ContPar1(j) = Spec%LimitPar1Min
                  end if
                  ! .. cap maximum usage @todo could consider penalising solution instead?
                  if (SChrom%ContPar1(j) .gt. Spec%LimitPar1Max) then
                    SChrom%ContPar1(j) = Spec%LimitPar1Max
                  end if
                  ! ... accumulate
                  nCumMat = nCumMat + nint(SChrom%ContPar1(j)) ! internally real, externally integer
                  ! ... did we reach Spec%nMat
                  if (nCumMat .ge. Spec%nMat * GenderMode) then
                    ! ... there should be exactly Spec%nMat contributions
                    if (nCumMat .gt. Spec%nMat * GenderMode) then
                      SChrom%ContPar1(j) = SChrom%ContPar1(j) - FLOATFUN(nCumMat - Spec%nMat * GenderMode) ! internally real, externally integer
                      ! ... did we go below the minimum usage limit?
                      if (nint(SChrom%ContPar1(j)) .lt. Spec%LimitPar1Min) then
                        TmpR = Spec%LimitPar1MinWeight * (Spec%LimitPar1Min - nint(SChrom%ContPar1(j))) ! internally real, externally integer
                        This%Objective = This%Objective + TmpR
                        if (Spec%LimitPar1MinWeight .lt. 0.0) then
                          This%PenaltyLimitPar1 = This%PenaltyLimitPar1 + TmpR
                          This%Penalty          = This%Penalty          + TmpR
                        end if
                        ! ... make sure we do not have negative contributions
                        if (SChrom%ContPar1(j) .lt. 0) then
                          SChrom%ContPar1(j) = 0.0
                        end if
                      end if
                      nCumMat = Spec%nMat * GenderMode
                    end if
                    ! ... remove contributions for the remaining contributors
                    do k = i + 1, Spec%nPar1
                      SChrom%ContPar1(Rank(k)) = 0.0
                    end do
                    exit
                  end if
                end do

                ! call SChrom%Write
                ! print*, sum(SChrom%ContPar1(:))
                ! pause

                ! ... Spec%nMat still not reached?
                do while (nCumMat .lt. Spec%nMat * GenderMode)
                  ! ... add more contributions to randomly chosen individuals
                  RanNumLoc = RanNumLoc + 1
                  i = int(RanNum(RanNumLoc) * Spec%nPar1) + 1
                  j = Rank(i)
                  if (nint(SChrom%ContPar1(j) + 1) .le. Spec%LimitPar1Max) then ! make sure we do not go above max
                    SChrom%ContPar1(j) = SChrom%ContPar1(j) + 1.0
                    ! ... accumulate
                    nCumMat = nCumMat + 1
                    ! ... did we reach Spec%nMat
                    if (nCumMat .ge. Spec%nMat * GenderMode) then
                      ! Internally real, externally integer
                      TmpI = sum(nint(SChrom%ContPar1(Rank(1:Spec%nPar1))))
                      if (TmpI .ne. Spec%nMat * GenderMode) then
                        if (TmpI .gt. Spec%nMat * GenderMode) then
                          SChrom%ContPar1(j) = FLOATFUN(nint(SChrom%ContPar1(j)) - 1)
                        else
                          SChrom%ContPar1(j) = FLOATFUN(nint(SChrom%ContPar1(j)) + 1)
                        end if
                      end if
                      exit
                    end if
                  end if
                end do
              end if

              ! call SChrom%Write
              ! print*, sum(SChrom%ContPar1(:))
              ! pause

              ! Left here for debugging
              ! do i = 1, Data%nPotPar1
              !   if (nint(SChrom%ContPar1(i)) .gt. Spec%LimitPar1Max) then
              !     print*, "n", i, nint(SChrom%ContPar1(i)), SChrom%ContPar1(i)
              !     print*, i, nint(SChrom%ContPar1)
              !     print*, i, SChrom%ContPar1
              !     stop
              !   end if
              ! end do

              ! "Parent2"
              if (Spec%GenderGiven) then
                ! ... preselect
                if (Spec%PreselectPar2) then
                  ! if (Spec%ModeSpec%ObjectiveCoancestry) then
                  !   Rank(1:Spec%PreselectPar2N) = RapKnr(-Data%AvgCoancestry(Data%IdPotPar2), Spec%PreselectPar2N) ! preselect contributors
                  !   TmpVec(1:Spec%PreselectPar2N) = SChrom%ContPar2(Rank(1:Spec%PreselectPar2N))                   ! save contributors
                  !   SChrom%ContPar2 = 0.0                                                                        ! set everyones contributions to zero
                  !   SChrom%ContPar2(Rank(1:Spec%PreselectPar2N)) = TmpVec(1:Spec%PreselectPar2N)                   ! put contributors back
                  ! end if
                  if (Spec%ModeSpec%ObjectiveCriterion) then
                    Rank(1:Spec%PreselectPar2N) = RapKnr(Data%SelCriterionStd(Data%IdPotPar2),  Spec%PreselectPar2N) ! preselect contributors
                    TmpVec(1:Spec%PreselectPar2N) = SChrom%ContPar2(Rank(1:Spec%PreselectPar2N))                     ! save contributors
                    SChrom%ContPar2 = 0.0                                                                          ! set everyones contributions to zero
                    SChrom%ContPar2(Rank(1:Spec%PreselectPar2N)) = TmpVec(1:Spec%PreselectPar2N)                     ! put contributors back
                  end if
                end if
                ! ... ranks to find contributors
                if (.not. (Spec%EqualizePar2 .and. (Spec%nPar2 .eq. Data%nPotPar2))) then
                  Rank(1:Spec%nPar2) = RapKnr(SChrom%ContPar2, Spec%nPar2)
                end if
                if (Spec%EqualizePar2) then ! ... equal contributions
                  if (Spec%nPar2 .eq. Data%nPotPar2) then
                    SChrom%ContPar2 = FLOATFUN(Spec%nMat) / Spec%nPar2 ! no need for indexing here, hence the above if (.not. ...)
                  else
                    SChrom%ContPar2 = 0.0
                    SChrom%ContPar2(Rank(1:Spec%nPar2)) = FLOATFUN(Spec%nMat) / Spec%nPar2
                  end if
                else                        ! ... unequal contributions
                  TmpVec(1:Spec%nPar2) = SChrom%ContPar2(Rank(1:Spec%nPar2)) ! save top contributions
                  SChrom%ContPar2 = 0.0                                    ! set everyones contributions to zero
                  SChrom%ContPar2(Rank(1:Spec%nPar2)) = TmpVec(1:Spec%nPar2) ! put top contributions back
                  nCumMat = 0
                  do i = 1, Spec%nPar2
                    j = Rank(i)
                    ! .. cap minimum usage @todo could consider penalising solution instead?
                    if (SChrom%ContPar2(j) .lt. Spec%LimitPar2Min) then
                      SChrom%ContPar2(j) = Spec%LimitPar2Min
                    end if
                    ! .. cap maximum usage @todo could consider penalising solution instead?
                    if (SChrom%ContPar2(j) .gt. Spec%LimitPar2Max) then
                      SChrom%ContPar2(j) = Spec%LimitPar2Max
                    end if
                    ! ... accumulate
                    nCumMat = nCumMat + nint(SChrom%ContPar2(j)) ! internally real, externally integer
                    ! ... did we reach Spec%nMat
                    if (nCumMat .ge. Spec%nMat) then
                      ! ... there should be exactly Spec%nMat contributions
                      if (nCumMat .gt. Spec%nMat) then
                        SChrom%ContPar2(j) = SChrom%ContPar2(j) - FLOATFUN(nCumMat - Spec%nMat)! internally real, externally integer
                        ! ... did we go below the minimum usage limit?
                        if (nint(SChrom%ContPar2(j)) .lt. Spec%LimitPar2Min) then
                          TmpR = Spec%LimitPar2MinWeight * (Spec%LimitPar2Min - nint(SChrom%ContPar2(j)))! internally real, externally integer
                          This%Objective = This%Objective + TmpR
                          if (Spec%LimitPar2MinWeight .lt. 0.0) then
                            This%PenaltyLimitPar2 = This%PenaltyLimitPar2 + TmpR
                            This%Penalty          = This%Penalty          + TmpR
                          end if
                          ! ... make sure we do not have negative contributions
                          if (SChrom%ContPar2(j) .lt. 0) then
                            SChrom%ContPar2(j) = 0.0
                          end if
                        end if
                        nCumMat = Spec%nMat
                      end if
                      ! ... remove contributions for the remaining contributors
                      do k = i + 1, Spec%nPar2
                        SChrom%ContPar2(Rank(k)) = 0.0
                      end do
                      exit
                    end if
                  end do

                  ! call SChrom%Write

                  ! ... Spec%nMat still not reached?
                  do while (nCumMat .lt. Spec%nMat)
                    ! ... add more contributions to randomly chosen individuals
                    RanNumLoc = RanNumLoc + 1
                    i = int(RanNum(RanNumLoc) * Spec%nPar2) + 1
                    j = Rank(i)
                    if (nint(SChrom%ContPar2(j) + 1) .le. Spec%LimitPar2Max) then ! make sure we do not go above max
                      SChrom%ContPar2(j) = SChrom%ContPar2(j) + 1.0
                      ! ... accumulate
                      nCumMat = nCumMat + 1
                      ! ...did we reach Spec%nMat
                      if (nCumMat .eq. Spec%nMat) then
                        ! Internally real, externally integer
                        TmpI = sum(nint(SChrom%ContPar2(Rank(1:Spec%nPar2))))
                        if (TmpI .ne. Spec%nMat) then
                          if (TmpI .gt. Spec%nMat) then
                            SChrom%ContPar2(j) = FLOATFUN(nint(SChrom%ContPar2(j)) - 1)
                          else
                            SChrom%ContPar2(j) = FLOATFUN(nint(SChrom%ContPar2(j)) + 1)
                          end if
                        end if
                        exit
                      end if
                    end if
                  end do
                end if
              end if

              ! call SChrom%Write

              ! --- Contributions (nVec) ---

              ! "Parent1"
              ! ... get integer values
              nVecPar1 = nint(SChrom%ContPar1)
              ! ... map chromosome to data order
              if (.not. Spec%GenderGiven) then
                This%nVec = nVecPar1
              else
                This%nVec(Data%IdPotPar1) = nVecPar1
              end if

              ! "Parent2"
              if (Spec%GenderGiven) then
                ! ... get integer values and map chromosome to data order
                This%nVec(Data%IdPotPar2) = nint(SChrom%ContPar2)
              end if

              ! call This%Write

              ! --- PAGE ---

              if (Spec%PAGEPar) then
                if (.not. Spec%GenderGiven) then
                  This%GenomeEdit(RapKnr(SChrom%EditPar1, Spec%PAGEPar1Max)) = 1.0
                else
                  if (Spec%PAGEPar1) then
                    This%GenomeEdit(Data%IdPotPar1(RapKnr(SChrom%EditPar1, Spec%PAGEPar1Max))) = 1.0
                  end if
                  if (Spec%PAGEPar2) then
                    This%GenomeEdit(Data%IdPotPar2(RapKnr(SChrom%EditPar2, Spec%PAGEPar2Max))) = 1.0
                  end if
                end if
              end if

              ! call This%Write

              ! --- Mate allocation ---

              if (Spec%MateAllocation) then
                if (Spec%GenderGiven) then
                  ! Distribute parent2 (=female) contributions into matings
                  k = 0
                  do i = 1, Data%nPotPar2 ! need to loop all females as some do not contribute
                    do j = 1, This%nVec(Data%IdPotPar2(i))
                      k = k + 1
                      MatPar2(k) = Data%IdPotPar2(i)
                    end do
                  end do
                  ! Reorder parent2 contributions according to the rank of matings
                  if (Spec%RandomMateAllocation) then
                    Rank(1:Spec%nMat) = RandomOrderIntel(n=Spec%nMat, Stream=Stream)
                  else
                    Rank(1:Spec%nMat) = MrgRnk(SChrom%MateRank) ! MrgRnk ranks small to large
                  end if
                  MatPar2 = MatPar2(Rank(1:Spec%nMat))
                else
                  ! Distribute one half of contributions into matings
                  k = 0
                  do while (k .lt. Spec%nMat)
                    do i = 1, Data%nPotPar1 ! need to loop all individuals as some do not contribute
                      l = This%nVec(Data%IdPotPar1(i)) / 2
                      if (mod(This%nVec(Data%IdPotPar1(i)), 2) .eq. 1) then
                        RanNumLoc = RanNumLoc + 1
                        if (RanNum(RanNumLoc) .gt. 0.5) then
                          l = l + 1
                        end if
                      end if
                      do j = 1, l
                        if (k .eq. Spec%nMat) then
                          exit
                        end if
                        k = k + 1
                        MatPar2(k) = Data%IdPotPar1(i)
                        nVecPar1(i) = nVecPar1(i) - 1 ! we substract here to remove "female" contributions
                      end do
                    end do
                  end do
                  ! Reorder one half of contributions according to the rank of matings
                  if (Spec%RandomMateAllocation) then
                    Rank(1:Spec%nMat) = RandomOrderIntel(n=Spec%nMat, Stream=Stream)
                  else
                    Rank(1:Spec%nMat) = MrgRnk(SChrom%MateRank) ! MrgRnk ranks small to large
                  end if
                  MatPar2 = MatPar2(Rank(1:Spec%nMat))
                end if

                ! print*, MatPar2

                ! Pair the contributions (=Mating plan)
                k = Spec%nMat ! MrgRnk ranks small to large
                ! We use nVecPar1 below instead of This%nVec because when .not. GenderGiven we modify nVecPar1 above, while can not modify This%nVec
                if (Spec%GenderGiven .or. Spec%SelfingAllowed) then
                  ! When gender matters selfing can not happen (we have two distinct sets of parents;
                  ! unless the user adds individuals of one sex in both sets) and when gender does not matter,
                  ! but selfing is allowed we do not need to care about it - faster code
                  do i = 1, Data%nPotPar1 ! need to loop all individuals as some do not contribute
                    ! print*, i, "/", Data%nPotPar1
                    do j = 1, nVecPar1(i)
                      if (k .gt. 0) then
                        ! if (k .lt. 10) write(*, '(i6,a1,i6,i6,a1,i6,i6,i6,a1,i6)') i, "/", Data%nPotPar1, j, "/", nVecPar1(i), sum(nVecPar1), k, "/", Spec%nMat
                        !                write(*, '(i6,a1,i6,i6,a1,i6,i6,i6,a1,i6)') i, "/", Data%nPotPar1, j, "/", nVecPar1(i), sum(nVecPar1), k, "/", Spec%nMat
                        Par1 = Data%IdPotPar1(i)
                        Par2 = MatPar2(k)
                        if (Spec%SelfingAllowed) then ! .not. GenderGiven
                          if (.not. Spec%ReciprocalMatingsAllowed) then
                            TmpMate = [Par1, Par2]
                            Par1 = maxval(TmpMate)
                            Par2 = minval(TmpMate)
                          end if
                        end if
                        This%MatingPlan(1, k) = Par1
                        This%MatingPlan(2, k) = Par2
                        k = k - 1
                      end if
                    end do
                  end do
                else
                  ! When gender does not matter, selfing can happen (we have one set of parents)
                  ! and when selfing is not allowed we need to avoid it - slower code
                  do i = 1, Data%nPotPar1
                    ! print*, i, "/", Data%nPotPar1
                    do j = 1, nVecPar1(i)
                      if (k .gt. 0) then
                        if (MatPar2(k) .eq. Data%IdPotPar1(i)) then
                          ! Try to avoid selfing by swapping the MatPar2 and Rank elements
                          do l = k, 1, -1
                            if (MatPar2(l) .ne. Data%IdPotPar1(i)) then
                              MatPar2([k, l]) = MatPar2([l, k])
                              SChrom%MateRank(Rank([k, l])) = SChrom%MateRank(Rank([l, k]))
                              exit
                            end if
                          end do
                          if (l .lt. 1) then ! Above loop ran out without finding a swap
                            This%Objective = This%Objective + Spec%SelfingWeight
                            if (Spec%SelfingWeight .lt. 0.0) then
                              This%PenaltySelfing = This%PenaltySelfing + Spec%SelfingWeight
                              This%Penalty        = This%Penalty        + Spec%SelfingWeight
                            end if
                          end if
                        end if
                        Par1 = Data%IdPotPar1(i)
                        Par2 = MatPar2(k)
                        if (.not. Spec%ReciprocalMatingsAllowed) then
                          TmpMate = [Par1, Par2]
                          Par1 = maxval(TmpMate)
                          Par2 = minval(TmpMate)
                        end if
                        This%MatingPlan(1, k) = Par1
                        This%MatingPlan(2, k) = Par2
                        ! print*, k, This%MatingPlan(:, k)
                        k = k - 1
                      end if
                    end do
                  end do
                end if

                ! Avoid multiple matings between the same male and female
                if (.not. Spec%MultipleMatingsAllowed) then
                  allocate(Pair(Spec%nMat, 2))
                  do i = 1, Spec%nMat
                    Pair(i, 1) = GeneratePairing(xin=This%MatingPlan(1, i), &
                                                 yin=This%MatingPlan(2, i))
                    ! print*, "i", i, "Spec%nMat", Spec%nMat, "Pair(i)", Pair(i), This%MatingPlan(1:2, i)
                  end do

                  ! Count multiple matings
                  Pair(:, 2) = MulCnt(Pair(:, 1))
                  ! print*,"Start 0"
                  ! do i = 1, Spec%nMat
                  !   print*, This%MatingPlan(:, i), Pair(i, 2)
                  ! end do
                  print*,"Stop 0",sum(Pair(:, 2)),Spec%nMat,sum(Pair(:, 2)) - Spec%nMat

                  ! Fix multiples by swapping
                  TmpI = 0
                  if ((sum(Pair(:, 2)) - Spec%nMat) .gt. 0) then
                    do j = 1, 10 ! @todo
                      do i = 1, (Spec%nMat - 1)
                        if (Pair(i, 2) .gt. 1) then
                          This%MatingPlan(2, [i, i + 1]) = This%MatingPlan(2, [i + 1, i])
                          Pair(i,     1) = GeneratePairing(xin=This%MatingPlan(1, i), &
                                                          yin=This%MatingPlan(2, i))
                          Pair(i + 1, 1) = GeneratePairing(xin=This%MatingPlan(1, i + 1), &
                                                          yin=This%MatingPlan(2, i + 1))
                          SChrom%MateRank([i, i + 1]) = SChrom%MateRank([i + 1, i])
                          Pair(:, 2) = MulCnt(Pair(:, 1)) ! recalculate (this might get expensive!)
                        end if
                      end do
                      print*,"Start ", j
                      k = 0
                      TmpI = 0
                      do i = 1, Spec%nMat
                        ! print*, This%MatingPlan(:, i), Pair(i, 2)
                        if (Pair(i, 2) .gt. 1) then
                          k = k + 1
                          TmpI = TmpI + Pair(i, 2) - 1
                        end if
                      end do
                      print*,"Stop ", j,sum(Pair(:, 2)),Spec%nMat,sum(Pair(:, 2)) - Spec%nMat
                      if (k .gt. 0) then
                        TmpI = TmpI / k
                      end if
                      if (TmpI .eq. 0 .or. Pair(Spec%nMat, 2) .gt. 1) then
                        exit
                      end if
                    end do
                  end if

                  ! Penalise
                  if (TmpI .gt. 0) then ! Above attempts were not (fully) successful
                    TmpR = TmpI * Spec%MultipleMatingsWeight
                    This%Objective = This%Objective + TmpR
                    if (Spec%MultipleMatingsWeight .lt. 0.0) then
                      This%PenaltyMultipleMatings = This%PenaltyMultipleMatings + TmpR
                      This%Penalty                = This%Penalty                + TmpR
                    end if
                  end if
                end if
              end if

              ! call This%Write

              ! --- Selection criterion ---
              ! Note that dot_product() is faster than dot() when called MANY times!

              if (Spec%SelCriterionGiven) then
                This%SelCriterionStd = dot_product(FLOATFUN(This%nVec), Data%SelCriterionStd) / (2.0 * Spec%nMat)
                if (Spec%PAGEPar) then
                  This%SelCriterionStd = This%SelCriterionStd + dot_product(FLOATFUN(This%nVec), Data%SelCriterionStdPAGE * This%GenomeEdit) / (2.0 * Spec%nMat)
                end if

                ! Inlined SelCriterionStd2SelCriterion
                This%SelCriterion = This%SelCriterionStd * Data%SelCriterionStat%Sd + Data%SelCriterionStat%Mean

                ! Inlined SelCriterionStd2MaxCriterionPct START
                Diff    =                      This%SelCriterionStd - Spec%ModeMinCoancestrySpec%SelCriterionStd
                MaxDiff = Spec%ModeMaxCriterionSpec%SelCriterionStd - Spec%ModeMinCoancestrySpec%SelCriterionStd
                if (MaxDiff .eq. 0) then
                  if (Diff .ge. 0) then
                    This%MaxCriterionPct = 100.0
                  else
                    This%MaxCriterionPct = 0.0
                  end if
                  ! Not sure about the above fix, but the logic is that if MaxDiff is zero,
                  ! then whatever positive (or negative) Diff we get, we achieve 100% (or 0%).
                else
                  This%MaxCriterionPct = Diff / MaxDiff * 100.0
                end if
                ! Inlined SelCriterionStd2MaxCriterionPct STOP

                if (Spec%ModeSpec%ObjectiveCriterion) then
                  if      (trim(Spec%ModeSpec%Name) .eq. "MaxCriterion") then
                    ! This gives objective in the [0, 1] form of the max SelCriterionStd
                    This%Objective = This%Objective + This%SelCriterionStd / Data%SelCriterionStat%Max
                  else if (trim(Spec%ModeSpec%Name) .eq. "Opt") then
                    ! This gives objective in the [0, 1] form of the TargetMaxCriterionPct
                    This%Objective = This%Objective + This%MaxCriterionPct / Spec%ModeSpec%TargetMaxCriterionPct
                  end if
                end if
                ! Handle beyond the nadir point case so that degree calculation will be meaningful
                if (.not. ieee_is_nan(This%MaxCriterionPct)) then
                  if (This%MaxCriterionPct .lt. 0.0) then
                    This%MaxCriterionPct = 0.0
                  end if
                end if
                ! @todo Should we handle also cases above 100%?
                ! @todo Should we modify Spec%ModeMaxCriterionSpec and Spec%ModeMinCoancestrySpec on the fly?
              end if

              ! --- Generic individual criterion ---
              ! Note that dot_product() is faster than dot() when called MANY times!

              if (Spec%GenericIndCritGiven) then
                do j = 1, Spec%nGenericIndCrit
                  TmpR = dot_product(FLOATFUN(This%nVec), Data%GenericIndCrit(:, j)) / (2.0 * Spec%nMat)
                  This%GenericIndCrit(j) = TmpR
                  TmpR = Spec%GenericIndCritWeight(j) * This%GenericIndCrit(j)
                  This%Objective = This%Objective + TmpR
                  if (Spec%GenericIndCritWeight(j) .lt. 0.0) then
                    This%PenaltyGenericIndCrit = This%PenaltyGenericIndCrit + TmpR
                    This%Penalty               = This%Penalty               + TmpR
                  end if
                end do
              end if

              ! --- Coancestry ---

              ! @todo Enable running AlphaMate without coancestry matrix, that is,
              !       select most performant individuals, but consider some limits
              !       on the amount of use etc. without coancestry among them.
              !       --> A hack would be to setup identity coancetry matrix, when one is not given

              ! Group coancestry x'Cx

              ! Via repeated use of dot(), https://software.intel.com/en-us/mkl-developer-reference-fortran-dot
              ! Note that dot_product() is faster than dot() when called MANY times!
              ! ... w=x'C
              do i = 1, Data%nInd
                TmpVec(i) = dot_product(FLOATFUN(This%nVec), Data%Coancestry%Value(1:, i))
              end do
              ! ... wx
              This%CoancestryRanMate = dot_product(TmpVec, FLOATFUN(This%nVec)) / (4.0 * Spec%nMat * Spec%nMat)

              ! Via BLAS subroutine
              ! This is slower than the above dot() code on a test case with n=370 (~35 sec vs. ~130 sec).
              ! On a large case (n=5120, but with different parameters than the small case) symv fails? (~45 sec vs. 135 sec).
              ! ... w=Cx, symmetric matrix times a vector https://software.intel.com/en-us/mkl-developer-reference-fortran-symv
              ! TmpVec = 0.0
              ! TmpVec2 = Data%Coancestry%Value(1:, 1:) ! needed to avoid putting this largish matrix on the stack
              ! call symv(A=TmpVec2, x=FLOATFUN(This%nVec), y=TmpVec)
              ! ... x'w
              ! This%CoancestryRanMate = dot(x=FLOATFUN(This%nVec), y=TmpVec) / (4.0 * Spec%nMat * Spec%nMat)

              ! Inlined Coancestry2CoancestryRate START
              Diff    = This%CoancestryRanMate - Data%CoancestryRanMate
              MaxDiff =                    1.0 - Data%CoancestryRanMate
              if (MaxDiff .eq. 0) then
                if (Diff .ge. 0) then
                  This%CoancestryRateRanMate =  1.0
                else if (Diff .eq. 0) then
                  This%CoancestryRateRanMate =  0.0
                else
                  This%CoancestryRateRanMate = -1.0
                end if
                ! Not sure about the above fix, but the logic is that if MaxDiff is zero,
                ! then whatever positive, zero, or negative) Diff we get, we respetively
                ! achieve rate of 1, 0, -1
              else
                This%CoancestryRateRanMate = Diff / MaxDiff
              end if
              ! Inlined Coancestry2CoancestryRate STOP

              ! Inlined CoancestryRate2MinCoancestryPct START
              Diff    = Spec%ModeMaxCriterionSpec%CoancestryRate - This%CoancestryRateRanMate
              MaxDiff = Spec%ModeMaxCriterionSpec%CoancestryRate - Spec%ModeMinCoancestrySpec%CoancestryRate
              if (MaxDiff .eq. 0) then
                if (Diff .ge. 0) then
                  This%MinCoancestryPct = 100.0
                else
                  This%MinCoancestryPct = 0.0
                end if
                ! Not sure about the above fix, but the logic is that if MaxDiff is zero,
                ! then whatever positive (or negative) Diff we get, we achieve 100% (or 0%).
              else
                This%MinCoancestryPct = Diff / MaxDiff * 100.0
              end if
              ! Inlined CoancestryRate2MinCoancestryPct STOP

              ! Handle beyond the nadir point case so that degree calculation will be meaningful
              if (.not. ieee_is_nan(This%MinCoancestryPct)) then
                if (This%MinCoancestryPct .lt. 0.0) then
                  This%MinCoancestryPct = 0.0
                end if
              end if
              ! @todo Should we handle also cases above 100%?
              ! @todo Should we then modify Spec%ModeMaxCriterionSpec and Spec%ModeMinCoancestrySpec on the fly?

              ! Degree
              ! This calculation ASSUMES unit circular shape of selection/coancestry frontier
              ! - given results from ModeMinCoancestry and ModeMaxCriterion we can map the solution on the unit circle,
              !   that is, the circle center is at (MaxCoancestryRate, MinSelCriterionStd) and solution is
              !   at (SolCoancestryRate, SolSelCriterionStd), which is mapped to (MinCoancestryPct, 1 - MaxCriterionPct)
              !   on unit circle with center at (1, 0) (MinCoancestryMode is at (0, 0) and MaxCriterionMode is at (1, 1))
              ! - then we can calculate degrees of the angle between the max-line (MaxCoancestryRate, MinSelCriterionStd)-(MaxCoancestryRate, MaxSelCriterionStd)
              !   and the sol-line (MaxCoancestryRate, MinSelCriterionStd)-(SolCoancestryRate, SolSelCriterionStd)
              ! - the min-line is (MaxCoancestryRate, MinSelCriterionStd)-(MinCoancestryRate, MinSelCriterionStd)
              ! - the solution y-coordinate on the max-line is adjacent to the angle and
              !   the solution x-coordinate on the min-line is opposite to the angle (if we
              !   put the min-line parallely up) so we use the arctangent function to compute
              !   angle degrees
              if (This%MinCoancestryPct .eq. 0 .and. This%MaxCriterionPct .eq. 0) then
                This%Degree = 45.0 ! a hack as atan2(0, 0) is not defined; following the logic that atan2(1, 1) = 45
              else
                This%Degree = atan2(This%MinCoancestryPct, This%MaxCriterionPct) * RAD2DEG
              end if

              ! @todo ModeRan?

              if (Spec%ModeSpec%ObjectiveCoancestry) then
                if      (trim(Spec%ModeSpec%Name) .eq. "MinCoancestry") then
                  This%Objective = This%Objective - This%CoancestryRateRanMate
                else if (trim(Spec%ModeSpec%Name) .eq. "Opt") then

                  ! Weight * (CoancestryRate - TargetCoancestryRate) [negative = good, positive = bad]
                  !   * smaller CoancestryRate is better
                  !     c( 0.05,       0.03, 0.02) - 0.02 =  0.03,  0.01,        0.0
                  !     c(-0.05, 0.00, 0.01, 0.02) - 0.02 = -0.07, -0.02, -0.01, 0.0
                  !   * this gives coancestry objective in the [-2, 2] form
                  !   * this penalty formulation works ~well with large weights as differences tend to be small,
                  !     but optimisation for criterion is hampered at low coancestries (~80 degrees)
                  ! TmpR = This%CoancestryRateRanMate - Spec%ModeSpec%TargetCoancestryRate
                  ! if (This%CoancestryRateRanMate .gt. Spec%ModeSpec%TargetCoancestryRate) then
                  !   ! CANCEL SelCriterion from the objective (want to focus on coancestry for such solutions)
                  !   This%Objective = This%Objective - This%MaxCriterionPct / Spec%ModeSpec%TargetMaxCriterionPct
                  ! else ! if (This%CoancestryRateRanMate .le. Spec%ModeSpec%TargetCoancestryRate) then
                  !   if (Spec%ModeSpec%CoancestryWeightBelow) then
                  !     TmpR = abs(TmpR)
                  !   else
                  !     TmpR = 0.0
                  !   end if
                  ! end if

                  ! Weight * (1 - MinCoancestryPct / TargetMinCoancestryPct) [negative = good, positive = bad]
                  !   * bigger MinCoancestryPct is better
                  !     1 - c(50, 60, 70)/70 =  0.28,  0.14, 0
                  !     1 - c(90, 80, 70)/70 = -0.28, -0.14, 0
                  !     1 -            0 /70 =  1
                  !   * this gives coancestry objective in the [-Inf, Inf] form of the TargetMinCoancestryPct
                  !   * this penalty formulation works ~well with weights close to -1,
                  !     but optimisation for criterion is hampered at low coancestries (~80 degrees)
                  ! TmpR = 1.0 - This%MinCoancestryPct / Spec%ModeSpec%TargetMinCoancestryPct
                  ! if (This%MinCoancestryPct .lt. Spec%ModeSpec%TargetMinCoancestryPct) then
                  !   ! CANCEL SelCriterion from the objective (want to focus on coancestry for such solutions)
                  !   This%Objective = This%Objective - This%MaxCriterionPct / Spec%ModeSpec%TargetMaxCriterionPct
                  ! else ! if (This%MinCoancestryPct .ge. Spec%ModeSpec%TargetMinCoancestryPct) then
                  !   if (Spec%ModeSpec%CoancestryWeightBelow) then
                  !     TmpR = abs(TmpR)
                  !   else
                  !     TmpR = 0.0
                  !   end if
                  ! end if

                  ! Weight * (1 - Degree / TargetDegree) [negative = good for coancestry, positive = bad for criterion]
                  !   * larger Degree is better in terms of coancestry, but smaller degree is better in terms of criterion
                  !     1 - c(50, 60, 70)/70 =  0.28,  0.14, 0
                  !     1 - c(90, 80, 70)/70 = -0.28, -0.14, 0
                  !   * this gives coancestry objective in the [-Inf, Inf] form of the TargetDegree
                  !   * this penalty formulation works ~well with weights close to -1 and works
                  !     well on any part of the frontier (due to exploiting its geometry, albeit assuming circular shape)
                  TmpR = 1.0 - This%Degree / Spec%ModeSpec%TargetDegree
                  if (This%Degree .lt. Spec%ModeSpec%TargetDegree) then
                    ! CANCEL SelCriterion from the objective (want to focus on coancestry for such solutions)
                    This%Objective = This%Objective - This%MaxCriterionPct / Spec%ModeSpec%TargetMaxCriterionPct
                  else ! if (This%Degree .ge. Spec%ModeSpec%TargetDegree) then
                    if (Spec%ModeSpec%CoancestryWeightBelow) then
                      TmpR = abs(TmpR)
                    else
                      TmpR = 0.0
                    end if
                  end if

                  TmpR = Spec%CoancestryWeight * TmpR
                  This%Objective = This%Objective + TmpR
                  if (Spec%CoancestryWeight .lt. 0.0) then
                    This%PenaltyCoancestry = This%PenaltyCoancestry + TmpR
                    This%Penalty           = This%Penalty           + TmpR
                  end if
                end if
              end if

              ! --- Expected progeny inbreeding (=inbreeding of a mating) ---

              if (Spec%MateAllocation) then
                TmpR = 0.0
                do j = 1, Spec%nMat
                  ! Lower triangle to speedup lookup
                  TmpMax = maxval(This%MatingPlan(:, j))
                  TmpMin = minval(This%MatingPlan(:, j))
                  TmpR = TmpR + Data%Coancestry%Value(TmpMax, TmpMin)
                end do
                ! @todo different number of progeny per mating? Might be relevant for in-vitro fertilisation
                This%Inbreeding = TmpR / Spec%nMat

                ! Inlined Coancestry2CoancestryRate START
                Diff    = This%Inbreeding - Data%Inbreeding
                MaxDiff =             1.0 - Data%Inbreeding
                if (MaxDiff .eq. 0) then
                  if (Diff .ge. 0) then
                    This%InbreedingRate =  1.0
                  else if (Diff .eq. 0) then
                    This%InbreedingRate =  0.0
                  else
                    This%InbreedingRate = -1.0
                  end if
                  ! Not sure about the above fix, but the logic is that if MaxDiff is zero,
                  ! then whatever positive, zero, or negative) Diff we get, we respetively
                  ! achieve rate of 1, 0, -1
                else
                  This%InbreedingRate = Diff / MaxDiff
                end if
                ! Inlined Coancestry2CoancestryRate STOP

                ! Inlined CoancestryRate2MinCoancestryPct START
                Diff    = Spec%ModeMaxCriterionSpec%InbreedingRate - This%InbreedingRate
                MaxDiff = Spec%ModeMaxCriterionSpec%InbreedingRate - Spec%ModeMinInbreedingSpec%InbreedingRate
                if (MaxDiff .eq. 0) then
                  if (Diff .ge. 0) then
                    This%MinInbreedingPct = 100.0
                  else
                    This%MinInbreedingPct = 0.0
                  end if
                  ! Not sure about the above fix, but the logic is that if MaxDiff is zero,
                  ! then whatever positive (or negative) Diff we get, we achieve 100% (or 0%).
                else
                  This%MinInbreedingPct = Diff / MaxDiff * 100.0
                end if
                ! Inlined CoancestryRate2MinCoancestryPct STOP

                ! @todo Pareto front formulation for this objective too?
                if (Spec%ModeSpec%ObjectiveInbreeding) then
                  if      (trim(Spec%ModeSpec%Name) .eq. "MinInbreeding") then
                    This%Objective = This%Objective - This%InbreedingRate
                  else if (trim(Spec%ModeSpec%Name) .eq. "Opt") then

                    ! Weight * (InbreedingRate - TargetInbreedingRate) [negative = good, positive = bad]
                    !   * smaller InbreedingRate is better
                    !     c( 0.05,       0.03, 0.02) - 0.02 =  0.03,  0.01,        0.0
                    !     c(-0.05, 0.00, 0.01, 0.02) - 0.02 = -0.07, -0.02, -0.01, 0.0
                    !   * this gives inbreeding objective in the [-2, 2] form
                    !   * this penalty formulation works ~well with large weights as differences tend to be small
                    ! TmpR = This%InbreedingRate - Spec%ModeSpec%TargetInbreedingRate
                    ! if (This%InbreedingRate .gt. Spec%ModeSpec%TargetInbreedingRate) then
                    !   ! CANCEL SelCriterion from the objective (want to focus on inbreeding for such solutions)
                    !   This%Objective = This%Objective - This%MaxCriterionPct / Spec%ModeSpec%TargetMaxCriterionPct
                    ! else ! if (This%InbreedingRate .le. Spec%ModeSpec%TargetInbreedingRate) then
                    !   if (Spec%ModeSpec%InbreedingWeightBelow) then
                    !     TmpR = abs(TmpR)
                    !   else
                    !     TmpR = 0.0
                    !   end if
                    ! end if

                    ! Weight * (1 - MinInbreedingPct / TargetMinInbreedingPct) [negative = good, positive = bad]
                    !   * bigger MinInbreedingPct is better
                    !     1 - c(50, 60, 70)/70 =  0.28,  0.14, 0
                    !     1 - c(90, 80, 70)/70 = -0.28, -0.14, 0
                    !     1 -            0 /70 =  1
                    !   * this gives inbreeding objective in the [-Inf, Inf] form of the TargetMinInbreedingPct
                    !   * this penalty formulation works ~well with weights close to -1
                    TmpR = 1.0 - This%MinInbreedingPct / Spec%ModeSpec%TargetMinInbreedingPct
                    if (This%MinInbreedingPct .lt. Spec%ModeSpec%TargetMinInbreedingPct) then
                      ! CANCEL SelCriterion from the objective (want to focus on inbreeding for such solutions)
                      This%Objective = This%Objective - This%MaxCriterionPct / Spec%ModeSpec%TargetMaxCriterionPct
                    else ! if (This%MinInbreedingPct .ge. Spec%ModeSpec%TargetMinInbreedingPct) then
                      if (Spec%ModeSpec%InbreedingWeightBelow) then
                        TmpR = abs(TmpR)
                      else
                        TmpR = 0.0
                      end if
                    end if

                    TmpR = Spec%InbreedingWeight * TmpR
                    This%Objective = This%Objective + TmpR
                    if (Spec%InbreedingWeight .lt. 0.0) then
                      This%PenaltyInbreeding = This%PenaltyInbreeding + TmpR
                      This%Penalty           = This%Penalty           + TmpR
                    end if
                  end if
                end if
              end if

              ! --- Generic mating criterion ---

              if (Spec%GenericMatCritGiven .and. Spec%MateAllocation) then
                do k = 1, Spec%nGenericMatCrit
                  TmpR = 0.0
                  if (Spec%GenderGiven) then
                    do j = 1, Spec%nMat
                      TmpR = TmpR + Data%GenericMatCrit(Data%IdPotParSeq(This%MatingPlan(1, j)), &
                                                        Data%IdPotParSeq(This%MatingPlan(2, j)), k)
                    end do
                  else
                    do j = 1, Spec%nMat
                      ! Speedup lookup
                      TmpMax = maxval(This%MatingPlan(:, j))
                      TmpMin = minval(This%MatingPlan(:, j))
                      TmpR = TmpR + Data%GenericMatCrit(TmpMax, TmpMin, k)
                    end do
                  end if
                  This%GenericMatCrit(k) = TmpR / Spec%nMat
                  TmpR = Spec%GenericMatCritWeight(k) * This%GenericMatCrit(k)
                  This%Objective = This%Objective + TmpR
                  if (Spec%GenericMatCritWeight(k) .lt. 0.0) then
                    This%PenaltyGenericMatCrit = This%PenaltyGenericMatCrit + TmpR
                    This%Penalty               = This%Penalty               + TmpR
                  end if
                end do
              end if

              ! --- Assign structured chromosome into solutions' chromosome vector ---

              Start = 1
              End = Data%nPotPar1
              This%Chrom(Start:End) = SChrom%ContPar1
              if (Spec%GenderGiven) then
                Start = End + 1
                End = Start - 1 + Data%nPotPar2
                This%Chrom(Start:End) = SChrom%ContPar2
              end if
              if (Spec%MateAllocation) then
                Start = End + 1
                End = Start - 1 + Spec%nMat
                This%Chrom(Start:End) = SChrom%MateRank
              end if
              if (Spec%PAGEPar) then
                if (Spec%PAGEPar1) then
                  Start = End + 1
                  End = Start - 1 + Data%nPotPar1
                  This%Chrom(Start:End) = SChrom%EditPar1
                end if
                if (Spec%GenderGiven) then
                  if (Spec%PAGEPar2) then
                    Start = End + 1
                    End = Start - 1 + Data%nPotPar2
                    This%Chrom(Start:End) = SChrom%EditPar2
                  end if
                end if
              end if

              ! call SChrom%Write
              ! call This%Write

          end select
      end select
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Call various optimisations for AlphaMate
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !---------------------------------------------------------------------------
    ! @todo the best random mating?
    ! @todo return all solutions?
    subroutine AlphaMateSearch(Spec, Data, LogStdout) ! not pure due to IO
      implicit none
      type(AlphaMateSpec), intent(inout) :: Spec      !< AlphaMateSpec holder (out because we set and reset some parameters for different search modes)
      type(AlphaMateData), intent(in)    :: Data      !< AlphaMateData holder (out because we set and reset some parameters for different search modes)
      logical, intent(in), optional      :: LogStdout !< Log process on stdout (default .false.)

      type(AlphaMateSol) :: SolMinCoancestry, SolMinInbreeding, SolMaxCriterion, Sol !< For frontier modes and random mating (no optimisation) mode

      integer(int32) :: nParam, Point, iSol, Target, Unit, nRanNum, RanNumLoc

      real(FLOATTYPE) :: Tmp
      real(FLOATTYPE), allocatable :: RanNum(:), AvgCoancestryStd(:), SelCriterionStd(:)
      real(FLOATTYPEAH), allocatable :: InitChrom(:, :)
      logical :: LogStdoutInternal !, OptimOK

      character(len=FILELENGTH) :: LogFile, LogPopFile, ContribFile, MatingFile

      real(FLOATTYPE) :: NANFLOATTYPE
      NANFLOATTYPE = Ieee_Value(x=NANFLOATTYPE, class=Ieee_Quiet_NaN)

      if (present(LogStdout)) then
        LogStdoutInternal = LogStdout
      else
        LogStdoutInternal = .false.
      end if

      ! Seting these two modes here in case they are not activated
      ! (need to get those NANs in (to avoid propagating 0s) as other modes use these results)
      call Spec%SetupMode(Mode="MinCoancestry", Data=Data)
      call Spec%SetupMode(Mode="MinInbreeding", Data=Data)
      call Spec%SetupMode(Mode="MaxCriterion",  Data=Data)

      ! --- Number of parameters to optimise ---

      nParam = Data%nPotPar1

      if (Spec%GenderGiven) then
        nParam = nParam + Data%nPotPar2
      end if

      if (Spec%MateAllocation) then
        nParam = nParam + Spec%nMat
      end if

      if (Spec%PAGEPar) then
        if (Spec%PAGEPar1) then
          nParam = nParam + Data%nPotPar1
        end if
        if (Spec%PAGEPar2) then
          nParam = nParam + Data%nPotPar2
        end if
      end if

      if (Spec%EvolAlgNSol .eq. 0) then
        if      (nParam .lt.  100) then
          Spec%EvolAlgNSol = 100
        else if (nParam .lt. 1000) then
          Spec%EvolAlgNSol = maxval([100, nint(0.5 * nParam)])
        else
          Spec%EvolAlgNSol = 500
        end if
      end if

      ! Presample random numbers
      nRanNum = 10 * Spec%EvolAlgNSol
      allocate(RanNum(nRanNum))
      RanNum = SampleIntelUniformD(n=nRanNum, Accurate=.false., Stream=AlphaMateStream)
      RanNumLoc = 0

      ! Distribute contributions, ranks, ... at random (exact values not important as we use ranks to get top values in evaluate)
      allocate(InitChrom(nParam, Spec%EvolAlgNSol))
      do iSol = 1, Spec%EvolAlgNSol
        InitChrom(:, iSol) = SampleIntelUniformD(n=nParam, Accurate=.false., Stream=AlphaMateStream)
      end do

      ! --- Standardized average coancestry ---

      allocate(AvgCoancestryStd(Data%nInd))
      ! note the front minus as we want to boost less related individuals
      AvgCoancestryStd = - ((Data%AvgCoancestry - Mean(Data%AvgCoancestry)) / StdDev(Data%AvgCoancestry))
      ! reorder to chromosome structure
      if (Spec%GenderGiven) then
        AvgCoancestryStd = AvgCoancestryStd([Data%IdPotPar1, Data%IdPotPar2])
      else
        AvgCoancestryStd = AvgCoancestryStd(Data%IdPotPar1)
      end if

      ! --- Standardized selection criterion ---

      allocate(SelCriterionStd(Data%nInd))
      ! reorder to chromosome structure
      if (Spec%GenderGiven) then
        SelCriterionStd = Data%SelCriterionStd([Data%IdPotPar1, Data%IdPotPar2])
      else
        SelCriterionStd = Data%SelCriterionStd(Data%IdPotPar1)
      end if

      ! --- Minimum coancestry ---

      if (Spec%ModeMinCoancestry) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") " Optimise contributions for minimum future coancestry (ModeMinCoancestry) ..."
          write(STDOUT, "(a)") " "
        end if

        ! Setup
        call Spec%SetupMode(Mode="MinCoancestry")

        LogFile     = trim(Spec%OutputBasename)//"OptimisationLogModeMinCoancestry.txt"
        LogPopFile  = trim(Spec%OutputBasename)//"OptimisationLogPopModeMinCoancestry.txt"
        ContribFile = trim(Spec%OutputBasename)//"ContributorsModeMinCoancestry.txt"
        MatingFile  = trim(Spec%OutputBasename)//"MatingPlanModeMinCoancestry.txt"

        ! Initialise
        ! @todo initialise with an approximate programming/exact solution?
        ! ... approximate minimum coancestry solution with equal contributions
        iSol = 1
          InitChrom(1:Data%nPotPar, iSol) = 0.0
          InitChrom(                 RapKnr(AvgCoancestryStd(1:Data%nPotPar1),                  Spec%nPar1),  iSol) = FLOATFUN(Spec%nMat) / Spec%nPar1
        if (Spec%GenderGiven) then
          InitChrom((Data%nPotPar1 + RapKnr(AvgCoancestryStd((Data%nPotPar1 + 1):Data%nPotPar), Spec%nPar2)), iSol) = FLOATFUN(Spec%nMat) / Spec%nPar2
        end if
        ! ... another one
        iSol = iSol + 1
          InitChrom(1:Data%nPotPar, iSol) = 0.0
          InitChrom(                 RapKnr(AvgCoancestryStd(1:Data%nPotPar1),                  Spec%nPar1),  iSol) = FLOATFUN(Spec%nMat) / Spec%nPar1
        if (Spec%GenderGiven) then
          InitChrom((Data%nPotPar1 + RapKnr(AvgCoancestryStd((Data%nPotPar1 + 1):Data%nPotPar), Spec%nPar2)), iSol) = FLOATFUN(Spec%nMat) / Spec%nPar2
        end if
        ! ... approximate minimum coancestry solution - in a different way
        iSol = iSol + 1
        InitChrom(1:Data%nPotPar, iSol) = AvgCoancestryStd
        ! ... another one
        iSol = iSol + 1
        InitChrom(1:Data%nPotPar, iSol) = AvgCoancestryStd
        ! ... noiser solutions
        do iSol = iSol + 1, Spec%EvolAlgNSol
          RanNumLoc = RanNumLoc + 1
          if (RanNumLoc .gt. nRanNum) then
            RanNumLoc = 1
          end if
          if (RanNum(RanNumLoc) .lt. 0.75) then ! keep 25% of purely random solution
            ! Multiply by standardized average coancestry to boost less related individuals
            InitChrom(1:Data%nPotPar, iSol) = InitChrom(1:Data%nPotPar, iSol) * AvgCoancestryStd
          end if
        end do

        ! Search
        if (trim(Spec%EvolAlg) .eq. "DE") then
          call DifferentialEvolution(Spec=Spec, Data=Data, nParam=nParam, nSol=Spec%EvolAlgNSol, Init=InitChrom, &
            nIter=Spec%EvolAlgNIter, nIterBurnIn=Spec%DiffEvolNIterBurnIn, nIterStop=Spec%EvolAlgNIterStop, StopTolerance=Spec%EvolAlgStopTolCoancestry, nIterPrint=Spec%EvolAlgNIterPrint, &
            LogStdout=LogStdoutInternal, LogFile=LogFile, LogPop=Spec%EvolAlgLogPop, LogPopFile=LogPopFile, &
            CRBurnIn=Spec%DiffEvolParamCrBurnIn, CRLate1=Spec%DiffEvolParamCr1, CRLate2=Spec%DiffEvolParamCr2, &
            FBase=Spec%DiffEvolParamFBase, FHigh1=Spec%DiffEvolParamFHigh1, FHigh2=Spec%DiffEvolParamFHigh2, &
            Stream=AlphaMateStream, BestSol=SolMinCoancestry)!, Status=OptimOK)
          ! if (.not. OptimOK) then
          !   write(STDERR, "(a)") " ERROR: Optimisation failed!"
          !   write(STDERR, "(a)") " "
          !   stop 1
          ! end if
        end if

        ! Save
        ! @todo Do these lines still make sense when we go above two objectives?
        SolMinCoancestry%Degree           =  90.0
        SolMinCoancestry%MinCoancestryPct = 100.0
        SolMinCoancestry%MinInbreedingPct =   0.0
        SolMinCoancestry%MaxCriterionPct  =   0.0
        call Spec%ModeMinCoancestrySpec%SaveSol2ModeSpec(In=SolMinCoancestry)
        call SolMinCoancestry%WriteContributions(Data=Data, Spec=Spec, ContribFile=ContribFile)
        if (Spec%MateAllocation) then
          call SolMinCoancestry%WriteMatingPlan(Data=Data, Spec=Spec, MatingFile=MatingFile)
        end if
      end if

      ! --- Minimum inbreeding ---

      if (Spec%ModeMinInbreeding) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") " Optimise contributions for minimum future inbreeding (ModeMinInbreeding) ..."
          write(STDOUT, "(a)") " "
        end if

        ! Setup
        call Spec%SetupMode(Mode="MinInbreeding")

        LogFile     = trim(Spec%OutputBasename)//"OptimisationLogModeMinInbreeding.txt"
        LogPopFile  = trim(Spec%OutputBasename)//"OptimisationLogPopModeMinInbreeding.txt"
        ContribFile = trim(Spec%OutputBasename)//"ContributorsModeMinInbreeding.txt"
        MatingFile  = trim(Spec%OutputBasename)//"MatingPlanModeMinInbreeding.txt"

        ! Initialise
        ! @todo initialise with an approximate programming/exact solutions?
        ! ... approximate minimum coancestry solution with equal contributions
        iSol = 1
          InitChrom(1:Data%nPotPar, iSol) = 0.0
          InitChrom(                 RapKnr(AvgCoancestryStd(1:Data%nPotPar1),                  Spec%nPar1),  iSol) = FLOATFUN(Spec%nMat) / Spec%nPar1
        if (Spec%GenderGiven) then
          InitChrom((Data%nPotPar1 + RapKnr(AvgCoancestryStd((Data%nPotPar1 + 1):Data%nPotPar), Spec%nPar2)), iSol) = FLOATFUN(Spec%nMat) / Spec%nPar2
        end if
        ! ... another one
        iSol = iSol + 1
          InitChrom(1:Data%nPotPar, iSol) = 0.0
          InitChrom(                 RapKnr(AvgCoancestryStd(1:Data%nPotPar1),                  Spec%nPar1),  iSol) = FLOATFUN(Spec%nMat) / Spec%nPar1
        if (Spec%GenderGiven) then
          InitChrom((Data%nPotPar1 + RapKnr(AvgCoancestryStd((Data%nPotPar1 + 1):Data%nPotPar), Spec%nPar2)), iSol) = FLOATFUN(Spec%nMat) / Spec%nPar2
        end if
        ! ... approximate minimum coancestry solution - in a different way
        iSol = iSol + 1
        InitChrom(1:Data%nPotPar, iSol) = AvgCoancestryStd
        ! ... another one
        iSol = iSol + 1
        InitChrom(1:Data%nPotPar, iSol) = AvgCoancestryStd
        ! ... noiser solutions
        do iSol = iSol + 1, Spec%EvolAlgNSol
          RanNumLoc = RanNumLoc + 1
          if (RanNumLoc .gt. nRanNum) then
            RanNumLoc = 1
          end if
          if (RanNum(RanNumLoc) .lt. 0.75) then ! keep 25% of purely random solution
            ! Multiply by standardized average coancestry to boost less related individuals
            InitChrom(1:Data%nPotPar, iSol) = InitChrom(1:Data%nPotPar, iSol) * AvgCoancestryStd
          end if
        end do

        ! Search
        if (trim(Spec%EvolAlg) .eq. "DE") then
          call DifferentialEvolution(Spec=Spec, Data=Data, nParam=nParam, nSol=Spec%EvolAlgNSol, Init=InitChrom, &
            nIter=Spec%EvolAlgNIter, nIterBurnIn=Spec%DiffEvolNIterBurnIn, nIterStop=Spec%EvolAlgNIterStop, StopTolerance=Spec%EvolAlgStopTolCoancestry, nIterPrint=Spec%EvolAlgNIterPrint, &
            LogStdout=LogStdoutInternal, LogFile=LogFile, LogPop=Spec%EvolAlgLogPop, LogPopFile=LogPopFile, &
            CRBurnIn=Spec%DiffEvolParamCrBurnIn, CRLate1=Spec%DiffEvolParamCr1, CRLate2=Spec%DiffEvolParamCr2, &
            FBase=Spec%DiffEvolParamFBase, FHigh1=Spec%DiffEvolParamFHigh1, FHigh2=Spec%DiffEvolParamFHigh2, &
            Stream=AlphaMateStream, BestSol=SolMinInbreeding)!, Status=OptimOK)
          ! if (.not. OptimOK) then
          !   write(STDERR, "(a)") " ERROR: Optimisation failed!"
          !   write(STDERR, "(a)") " "
          !   stop 1
          ! end if
        end if

        ! Save
        ! @todo Do these lines still make sense when we go above two objectives?
        SolMinInbreeding%Degree           =  45.0
        SolMinInbreeding%MinCoancestryPct =   0.0
        SolMinInbreeding%MinInbreedingPct = 100.0
        SolMinInbreeding%MaxCriterionPct  =   0.0
        call Spec%ModeMinInbreedingSpec%SaveSol2ModeSpec(In=SolMinInbreeding)
        ! check the optimums
        if (Spec%ModeMinCoancestry) then
          if (Spec%ModeMinInbreedingSpec%Inbreeding .gt. Spec%ModeMinCoancestrySpec%Inbreeding) then
            Spec%ModeMinInbreedingSpec%Inbreeding     = Spec%ModeMinCoancestrySpec%Inbreeding
            Spec%ModeMinInbreedingSpec%InbreedingRate = Spec%ModeMinCoancestrySpec%InbreedingRate
          end if
          if (Spec%ModeMinCoancestrySpec%Coancestry .gt. Spec%ModeMinInbreedingSpec%Coancestry) then
            Spec%ModeMinCoancestrySpec%Coancestry     = Spec%ModeMinInbreedingSpec%Coancestry
            Spec%ModeMinCoancestrySpec%CoancestryRate = Spec%ModeMinInbreedingSpec%CoancestryRate
          end if
        end if
        call SolMinInbreeding%WriteContributions(Data=Data, Spec=Spec, ContribFile=ContribFile)
        if (Spec%MateAllocation) then
          call SolMinInbreeding%WriteMatingPlan(Data=Data, Spec=Spec, MatingFile=MatingFile)
        end if
      end if

      ! --- Maximum selection criterion ---

      if (Spec%ModeMaxCriterion) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") " Optimise contributions for maximum future selection criterion (ModeMaxCriterion) ..."
          write(STDOUT, "(a)") " "
        end if

        ! Setup
        call Spec%SetupMode(Mode="MaxCriterion")

        LogFile     = trim(Spec%OutputBasename)//"OptimisationLogModeMaxCriterion.txt"
        LogPopFile  = trim(Spec%OutputBasename)//"OptimisationLogPopModeMaxCriterion.txt"
        ContribFile = trim(Spec%OutputBasename)//"ContributorsModeMaxCriterion.txt"
        MatingFile  = trim(Spec%OutputBasename)//"MatingPlanModeMaxCriterion.txt"

        ! Initialise
        ! ... exact truncation selection solution with equal contributions
        iSol = 1
          InitChrom(1:Data%nPotPar, iSol) = 0.0
          InitChrom(                 RapKnr(SelCriterionStd(1:Data%nPotPar1),                  Spec%nPar1),  iSol) = FLOATFUN(Spec%nMat) / Spec%nPar1
        if (Spec%GenderGiven) then
          InitChrom((Data%nPotPar1 + RapKnr(SelCriterionStd((Data%nPotPar1 + 1):Data%nPotPar), Spec%nPar2)), iSol) = FLOATFUN(Spec%nMat) / Spec%nPar2
        end if
        ! ... another one
        iSol = iSol + 1
          InitChrom(1:Data%nPotPar, iSol) = 0.0
          InitChrom(                 RapKnr(SelCriterionStd(1:Data%nPotPar1),                  Spec%nPar1),  iSol) = FLOATFUN(Spec%nMat) / Spec%nPar1
        if (Spec%GenderGiven) then
          InitChrom((Data%nPotPar1 + RapKnr(SelCriterionStd((Data%nPotPar1 + 1):Data%nPotPar), Spec%nPar2)), iSol) = FLOATFUN(Spec%nMat) / Spec%nPar2
        end if
        ! ... approximate truncation selection solution
        iSol = iSol + 1
        InitChrom(1:Data%nPotPar, iSol) = SelCriterionStd
        ! ... another one
        iSol = iSol + 1
        InitChrom(1:Data%nPotPar, iSol) = SelCriterionStd
        ! ... noiser solutions
        do iSol = iSol + 1, Spec%EvolAlgNSol
          RanNumLoc = RanNumLoc + 1
          if (RanNumLoc .gt. nRanNum) then
            RanNumLoc = 1
          end if
          if (RanNum(RanNumLoc) .lt. 0.75) then ! keep 25% of purely random solution
            ! Multiply by standardized selection criterion to boost better individuals
            InitChrom(1:Data%nPotPar, iSol) = InitChrom(1:Data%nPotPar, iSol) * SelCriterionStd
          end if
        end do

        ! Search
        if (trim(Spec%EvolAlg) .eq. "DE") then
          call DifferentialEvolution(Spec=Spec, Data=Data, nParam=nParam, nSol=Spec%EvolAlgNSol, Init=InitChrom, &
            nIter=Spec%EvolAlgNIter, nIterBurnIn=Spec%DiffEvolNIterBurnIn, nIterStop=Spec%EvolAlgNIterStop, StopTolerance=Spec%EvolAlgStopTol, nIterPrint=Spec%EvolAlgNIterPrint, &
            LogStdout=LogStdoutInternal, LogFile=LogFile, LogPop=Spec%EvolAlgLogPop, LogPopFile=LogPopFile, &
            CRBurnIn=Spec%DiffEvolParamCrBurnIn, CRLate1=Spec%DiffEvolParamCr1, CRLate2=Spec%DiffEvolParamCr2, &
            FBase=Spec%DiffEvolParamFBase, FHigh1=Spec%DiffEvolParamFHigh1, FHigh2=Spec%DiffEvolParamFHigh2, &
            Stream=AlphaMateStream, BestSol=SolMaxCriterion)!, Status=OptimOK)
          ! if (.not. OptimOK) then
          !   write(STDERR, "(a)") " ERROR: Optimisation failed!"
          !   write(STDERR, "(a)") " "
          !   stop 1
          ! end if
        end if

        ! Save
        ! @todo Do these lines still make sense when we go above two objectives?
        SolMaxCriterion%Degree           =   0.0
        SolMaxCriterion%MinCoancestryPct =   0.0
        SolMinInbreeding%MinInbreedingPct =  0.0
        SolMaxCriterion%MaxCriterionPct  = 100.0
        call Spec%ModeMaxCriterionSpec%SaveSol2ModeSpec(In=SolMaxCriterion)
        call SolMaxCriterion%WriteContributions(Data=Data, Spec=Spec, ContribFile=ContribFile)
        if (Spec%MateAllocation) then
          call SolMaxCriterion%WriteMatingPlan(Data=Data, Spec=Spec, MatingFile=MatingFile)
        end if
      end if

      ! --- Evaluate frontier ---
      ! ModeMinCoancestry, ModeMinInbreeding, and ModeMaxCriterion must be run prior to this!
      ! @todo This code will have to change in light of more than two objectives using the NBI or AWS method

      if (Spec%EvaluateFrontier) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") " Evaluate frontier ..."
        end if

        open(newunit=Unit, file=trim(Spec%OutputBasename)//"Frontier.txt", status="unknown")

        ! Setup
        call Spec%LogHead(LogUnit=Unit, String="ModeOrPoint", StringNum=18)

        ! Add minimum coancestry solution to frontier output (90 degress with two objectives)
        call SolMinCoancestry%Log(Spec=Spec, LogUnit=Unit, Iteration=-1, AcceptPct=real(NANFLOATTYPE), String="ModeMinCoancestry", StringNum=18)

        ! Add minimum inbreeding solution to frontier output
        if (Spec%ModeMinInbreeding) then
          call SolMinInbreeding%Log(Spec=Spec, LogUnit=Unit, Iteration=-1, AcceptPct=real(NANFLOATTYPE), String="ModeMinInbreeding", StringNum=18)
        end if

        ! Add maximum criterion solution to frontier output (0 degress with two objectives)
        call SolMaxCriterion%Log(Spec=Spec, LogUnit=Unit, Iteration=-1, AcceptPct=real(NANFLOATTYPE), String="ModeMaxCriterion", StringNum=18)

        ! Frontier
        do Point = 1, size(TARGETDEGREEFRONTIER) ! 80, 70, ..., 10 degrees

          ! Setup
          call Spec%SetupMode(Mode="Opt", Data=Data, &
                              Degree=TARGETDEGREEFRONTIER(Point), &
                              ModeMinCoancestrySpec=Spec%ModeMinCoancestrySpec, &
                              ModeMinInbreedingSpec=Spec%ModeMinInbreedingSpec, &
                              ModeMaxCriterionSpec=Spec%ModeMaxCriterionSpec)
          if (LogStdoutInternal) then
            write(STDOUT, "(a)") " "
            write(STDOUT, "(a)") "  Point "//Int2Char(Point)//" of "//Int2Char(size(TARGETDEGREEFRONTIER))
            write(STDOUT, "(a)") " "
            call Spec%ModeSpec%LogTargets(Unit=STDOUT, Spec=Spec)
          end if

          LogFile     = trim(Spec%OutputBasename)//"OptimisationLogModeFrontier"//trim(Int2Char(Point))//".txt"
          LogPopFile  = trim(Spec%OutputBasename)//"OptimisationLogPopModeFrontier"//trim(Int2Char(Point))//".txt"
          ContribFile = trim(Spec%OutputBasename)//"ContributorsModeFrontier"//trim(Int2Char(Point))//".txt"
          MatingFile  = trim(Spec%OutputBasename)//"MatingPlanModeFrontier"//trim(Int2Char(Point))//".txt"

          ! Initialise
          ! @todo initialise with an approximate programming/exact solution?
          ! @todo implement evolutionary algorithm for multiple-objectives so that we evolve the whole frontier at once
          ! ... The MaxCriterion solution
          iSol = 1
          InitChrom(:, iSol) = SolMaxCriterion%Chrom
          ! ... exact truncation selection solution with equal contributions
          iSol = iSol + 1
            InitChrom(1:Data%nPotPar, iSol) = 0.0
            InitChrom(                 RapKnr(SelCriterionStd(1:Data%nPotPar1),                  Spec%nPar1),  iSol) = FLOATFUN(Spec%nMat) / Spec%nPar1
          if (Spec%GenderGiven) then
            InitChrom((Data%nPotPar1 + RapKnr(SelCriterionStd((Data%nPotPar1 + 1):Data%nPotPar), Spec%nPar2)), iSol) = FLOATFUN(Spec%nMat) / Spec%nPar2
          end if
          ! ... approximate truncation selection solution
          iSol = iSol + 1
          InitChrom(1:Data%nPotPar, iSol) = SelCriterionStd
          ! ... The previous target solution
          if (Point .gt. 1) then
            iSol = iSol + 1
            InitChrom(:, iSol) = Sol%Chrom
          end if
          ! ... approximate minimum coancestry solution with equal contributions
          iSol = iSol + 1
            InitChrom(1:Data%nPotPar, iSol) = 0.0
            InitChrom(                 RapKnr(AvgCoancestryStd(1:Data%nPotPar1),                  Spec%nPar1),  iSol) = FLOATFUN(Spec%nMat) / Spec%nPar1
          if (Spec%GenderGiven) then
            InitChrom((Data%nPotPar1 + RapKnr(AvgCoancestryStd((Data%nPotPar1 + 1):Data%nPotPar), Spec%nPar2)), iSol) = FLOATFUN(Spec%nMat) / Spec%nPar2
          end if
          ! ... approximate minimum coancestry solution
          iSol = iSol + 1
          InitChrom(1:Data%nPotPar, iSol) = AvgCoancestryStd
          ! ... The MinCoancestry solution
          iSol = iSol + 1
          InitChrom(:, iSol) = SolMinCoancestry%Chrom
          ! ... noiser solutions
          Tmp = (100.0 - 100.0 / 90.0 * TARGETDEGREEFRONTIER(Point)) / 100.0
          do iSol = iSol + 1, Spec%EvolAlgNSol
            RanNumLoc = RanNumLoc + 1
            if (RanNumLoc .gt. nRanNum) then
              RanNumLoc = 1
            end if
            if (RanNum(RanNumLoc) .lt. 0.75) then ! keep 25% of purely random solution
              RanNumLoc = RanNumLoc + 1
              if (RanNumLoc .gt. nRanNum) then
                RanNumLoc = 1
              end if
              if (RanNum(RanNumLoc) .lt. Tmp) then
                RanNumLoc = RanNumLoc + 1
                if (RanNumLoc .gt. nRanNum) then
                  RanNumLoc = 1
                end if
                if (RanNum(RanNumLoc) .lt. 0.5 .and. Point .gt. 1) then
                  ! Multiply by the previous target solution
                  InitChrom(:, iSol) =                  InitChrom(:, iSol) * Sol%Chrom
                else
                  RanNumLoc = RanNumLoc + 1
                  if (RanNumLoc .gt. nRanNum) then
                    RanNumLoc = 1
                  end if
                  if (RanNum(RanNumLoc) .lt. 0.5) then
                    ! Multiply by standardized selection criterion to boost better individuals
                    InitChrom(1:Data%nPotPar, iSol) =   InitChrom(1:Data%nPotPar, iSol) * SelCriterionStd
                  else
                    ! Multiply by product to boost better that are less individuals (note the - in front!)
                    InitChrom(1:Data%nPotPar, iSol) = - InitChrom(1:Data%nPotPar, iSol) * SelCriterionStd * AvgCoancestryStd
                  end if
                end if
              else
                RanNumLoc = RanNumLoc + 1
                if (RanNumLoc .gt. nRanNum) then
                  RanNumLoc = 1
                end if
                if (RanNum(RanNumLoc) .lt. 0.5 .and. Point .gt. 1) then
                  ! Multiply by the previous target solution
                  InitChrom(:, iSol) =                  InitChrom(:, iSol) * Sol%Chrom
                else
                  RanNumLoc = RanNumLoc + 1
                  if (RanNumLoc .gt. nRanNum) then
                    RanNumLoc = 1
                  end if
                  if (RanNum(RanNumLoc) .lt. 0.5) then
                    ! Multiply by product to boost better that are less individuals (note the - in front!)
                    InitChrom(1:Data%nPotPar, iSol) = - InitChrom(1:Data%nPotPar, iSol) * SelCriterionStd * AvgCoancestryStd
                  else
                    ! Multiply by standardized average coancestry to boost less related individuals
                    InitChrom(1:Data%nPotPar, iSol) =   InitChrom(1:Data%nPotPar, iSol) * AvgCoancestryStd
                  end if
                end if
              end if
            end if
          end do

          ! Search
          if (trim(Spec%EvolAlg) .eq. "DE") then
            call DifferentialEvolution(Spec=Spec, Data=Data, nParam=nParam, nSol=Spec%EvolAlgNSol, Init=InitChrom, &
              nIter=Spec%EvolAlgNIter, nIterBurnIn=Spec%DiffEvolNIterBurnIn, nIterStop=Spec%EvolAlgNIterStop, StopTolerance=Spec%EvolAlgStopTol, nIterPrint=Spec%EvolAlgNIterPrint, &
              LogStdout=LogStdoutInternal, LogFile=LogFile, LogPop=Spec%EvolAlgLogPop, LogPopFile=LogPopFile, &
              CRBurnIn=Spec%DiffEvolParamCrBurnIn, CRLate1=Spec%DiffEvolParamCr1, CRLate2=Spec%DiffEvolParamCr2, &
              FBase=Spec%DiffEvolParamFBase, FHigh1=Spec%DiffEvolParamFHigh1, FHigh2=Spec%DiffEvolParamFHigh2, &
              Stream=AlphaMateStream, BestSol=Sol)!, Status=OptimOK)
            ! if (.not. OptimOK) then
            !   write(STDERR, "(a)") " ERROR: Optimisation failed!"
            !   write(STDERR, "(a)") " "
            !   stop 1
            ! end if
          end if

          ! Save
          call Sol%Log(Spec=Spec, LogUnit=Unit, Iteration=-1, AcceptPct=real(NANFLOATTYPE), String=trim("ModeFrontier"//trim(Int2Char(Point))), StringNum=18)
          call Sol%WriteContributions(Data=Data, Spec=Spec, ContribFile=ContribFile)
          if (Spec%MateAllocation) then
            call Sol%WriteMatingPlan(Data=Data, Spec=Spec, MatingFile=MatingFile)
          end if

        end do

        close(Unit)

      end if

      ! --- Random mating ---

      ! @todo is this mode any good?
      if (Spec%ModeRan) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") " Evaluate random mating/crossing (ModeRan) ..."
          write(STDOUT, "(a)") " "
        end if

        ! Setup
        call Spec%SetupMode(Mode="Ran")

        LogFile = trim(Spec%OutputBasename)//"OptimisationLogModeRan.txt"
        ! @todo other reports from here - at least the best random solution?

        ! Search
        call RandomSearch(Mode="avg", Spec=Spec, Data=Data, nParam=nParam, Init=InitChrom, &
          nSamp=Spec%EvolAlgNSol*Spec%EvolAlgNIter*Spec%RanAlgStricter, nSampStop=Spec%EvolAlgNIterStop*Spec%RanAlgStricter, &
          StopTolerance=Spec%EvolAlgStopTol/Spec%RanAlgStricter, nSampPrint=Spec%EvolAlgNIterPrint, &
          LogStdout=LogStdoutInternal, LogFile=LogFile, BestSol=Sol)

      end if

      ! --- Maximum selection criterion with constraint on coancestry and inbreeding ---
      ! ModeMinCoancestry, ModeMinInbreeding, and ModeMaxCriterion must be run prior to this!
      ! @todo Will this code have to change in light of more than two objectives using the NBI or AWS method?

      if (Spec%ModeOpt .and. (Spec%nTargets .gt. 0)) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") " Optimise contributions for maximum future selection criterion with constraint on coancestry (and inbreeding) ..."
        end if

        open(newunit=Unit, file=trim(Spec%OutputBasename)//"Targets.txt", status="unknown")

        ! Setup
        call Spec%LogHead(LogUnit=Unit, String="Target", StringNum=18)

        ! Add minimum coancestry solution to target output (90 degress with two objectives)
        call SolMinCoancestry%Log(Spec=Spec, LogUnit=Unit, Iteration=-1, AcceptPct=real(NANFLOATTYPE), String="ModeMinCoancestry", StringNum=18)

        ! Add minimum inbreeding solution to target output
        if (Spec%ModeMinInbreeding) then
          call SolMinInbreeding%Log(Spec=Spec, LogUnit=Unit, Iteration=-1, AcceptPct=real(NANFLOATTYPE), String="ModeMinInbreeding", StringNum=18)
        end if

        ! Add maximum criterion solution to target output (0 degress with two objectives)
        call SolMaxCriterion%Log(Spec=Spec, LogUnit=Unit, Iteration=-1, AcceptPct=real(NANFLOATTYPE), String="ModeMaxCriterion", StringNum=18)

        ! Targets
        do Target = 1, Spec%nTargets

          ! Setup
          select case (trim(Spec%AllTargets(Target)))
            case ("Degree")
              call Spec%SetupMode(Mode="Opt", Data=Data, &
                                  Degree=Spec%AllTargetValues(Target), &
                                  ModeMinCoancestrySpec=Spec%ModeMinCoancestrySpec, &
                                  ModeMinInbreedingSpec=Spec%ModeMinInbreedingSpec, &
                                  ModeMaxCriterionSpec=Spec%ModeMaxCriterionSpec)
            case ("SelCriterion")
              call Spec%SetupMode(Mode="Opt", Data=Data, &
                                  SelCriterion=Spec%AllTargetValues(Target), &
                                  ModeMinCoancestrySpec=Spec%ModeMinCoancestrySpec, &
                                  ModeMinInbreedingSpec=Spec%ModeMinInbreedingSpec, &
                                  ModeMaxCriterionSpec=Spec%ModeMaxCriterionSpec)

            case ("SelCriterionStd")
              call Spec%SetupMode(Mode="Opt", Data=Data, &
                                  SelCriterionStd=Spec%AllTargetValues(Target), &
                                  ModeMinCoancestrySpec=Spec%ModeMinCoancestrySpec, &
                                  ModeMinInbreedingSpec=Spec%ModeMinInbreedingSpec, &
                                  ModeMaxCriterionSpec=Spec%ModeMaxCriterionSpec)

            case ("MaxCriterionPct")
              call Spec%SetupMode(Mode="Opt", Data=Data, &
                                  MaxCriterionPct=Spec%AllTargetValues(Target), &
                                  ModeMinCoancestrySpec=Spec%ModeMinCoancestrySpec, &
                                  ModeMinInbreedingSpec=Spec%ModeMinInbreedingSpec, &
                                  ModeMaxCriterionSpec=Spec%ModeMaxCriterionSpec)

            case ("Coancestry")
              call Spec%SetupMode(Mode="Opt", Data=Data, &
                                  Coancestry=Spec%AllTargetValues(Target), &
                                  ModeMinCoancestrySpec=Spec%ModeMinCoancestrySpec, &
                                  ModeMinInbreedingSpec=Spec%ModeMinInbreedingSpec, &
                                  ModeMaxCriterionSpec=Spec%ModeMaxCriterionSpec)

            case ("CoancestryRate")
              call Spec%SetupMode(Mode="Opt", Data=Data, &
                                  CoancestryRate=Spec%AllTargetValues(Target), &
                                  ModeMinCoancestrySpec=Spec%ModeMinCoancestrySpec, &
                                  ModeMinInbreedingSpec=Spec%ModeMinInbreedingSpec, &
                                  ModeMaxCriterionSpec=Spec%ModeMaxCriterionSpec)

            case ("MinCoancestryPct")
              call Spec%SetupMode(Mode="Opt", Data=Data, &
                                  MinCoancestryPct=Spec%AllTargetValues(Target), &
                                  ModeMinCoancestrySpec=Spec%ModeMinCoancestrySpec, &
                                  ModeMinInbreedingSpec=Spec%ModeMinInbreedingSpec, &
                                  ModeMaxCriterionSpec=Spec%ModeMaxCriterionSpec)
          end select

          if (LogStdoutInternal) then
            write(STDOUT, "(a)") " "
            write(STDOUT, "(a)") "  Target "//Int2Char(Target)//" of "//Int2Char(Spec%nTargets)
            write(STDOUT, "(a)") " "
            call Spec%ModeSpec%LogTargets(Unit=STDOUT, Spec=Spec)
          end if

          LogFile     = trim(Spec%OutputBasename)//"OptimisationLogModeOptTarget"//trim(Int2Char(Target))//".txt"
          LogPopFile  = trim(Spec%OutputBasename)//"OptimisationLogPopModeOptTarget"//trim(Int2Char(Target))//".txt"
          ContribFile = trim(Spec%OutputBasename)//"ContributorsModeOptTarget"//trim(Int2Char(Target))//".txt"
          MatingFile  = trim(Spec%OutputBasename)//"MatingPlanModeOptTarget"//trim(Int2Char(Target))//".txt"

          ! Initialise
          ! @todo initialise with an approximate programming/exact solution?
          ! ... The MaxCriterion solution
          iSol = 1
          InitChrom(:, iSol) = SolMaxCriterion%Chrom
          ! ... exact truncation selection solution with equal contributions
          iSol = iSol + 1
            InitChrom(1:Data%nPotPar, iSol) = 0.0
            InitChrom(                 RapKnr(SelCriterionStd(1:Data%nPotPar1),                  Spec%nPar1),  iSol) = FLOATFUN(Spec%nMat) / Spec%nPar1
          if (Spec%GenderGiven) then
            InitChrom((Data%nPotPar1 + RapKnr(SelCriterionStd((Data%nPotPar1 + 1):Data%nPotPar), Spec%nPar2)), iSol) = FLOATFUN(Spec%nMat) / Spec%nPar2
          end if
          ! ... approximate truncation selection solution
          iSol = iSol + 1
          InitChrom(1:Data%nPotPar, iSol) = SelCriterionStd
          ! ... approximate minimum coancestry solution with equal contributions
          iSol = iSol + 1
            InitChrom(1:Data%nPotPar, iSol) = 0.0
            InitChrom(                 RapKnr(AvgCoancestryStd(1:Data%nPotPar1),                  Spec%nPar1),  iSol) = FLOATFUN(Spec%nMat) / Spec%nPar1
          if (Spec%GenderGiven) then
            InitChrom((Data%nPotPar1 + RapKnr(AvgCoancestryStd((Data%nPotPar1 + 1):Data%nPotPar), Spec%nPar2)), iSol) = FLOATFUN(Spec%nMat) / Spec%nPar2
          end if
          ! ... approximate minimum coancestry solution
          iSol = iSol + 1
          InitChrom(1:Data%nPotPar, iSol) = AvgCoancestryStd
          ! ... The MinCoancestry solution
          iSol = iSol + 1
          InitChrom(:, iSol) = SolMinCoancestry%Chrom
          ! ... noiser solutions
          if (Spec%ModeSpec%ObjectiveCriterion .and. Spec%ModeSpec%ObjectiveCoancestry) then
            Tmp = (100.0 - 100.0 / 90.0 * Spec%ModeSpec%TargetDegree)
          else
            Tmp = 0.5
          end if
          do iSol = iSol + 1, Spec%EvolAlgNSol
            RanNumLoc = RanNumLoc + 1
            if (RanNumLoc .gt. nRanNum) then
              RanNumLoc = 1
            end if
            if (RanNum(RanNumLoc) .lt. 0.75) then ! keep 25% of purely random solution
              RanNumLoc = RanNumLoc + 1
              if (RanNumLoc .gt. nRanNum) then
                RanNumLoc = 1
              end if
              if (RanNum(RanNumLoc) .lt. Tmp) then
                RanNumLoc = RanNumLoc + 1
                if (RanNumLoc .gt. nRanNum) then
                  RanNumLoc = 1
                end if
                if (RanNum(RanNumLoc) .lt. 0.5) then
                  ! Multiply by standardized selection criterion to boost better individuals
                  InitChrom(1:Data%nPotPar, iSol) =   InitChrom(1:Data%nPotPar, iSol) * SelCriterionStd
                else
                  ! Multiply by product to boost better that are less individuals (note the - in front!)
                  InitChrom(1:Data%nPotPar, iSol) = - InitChrom(1:Data%nPotPar, iSol) * SelCriterionStd * AvgCoancestryStd
                end if
              else
                RanNumLoc = RanNumLoc + 1
                if (RanNumLoc .gt. nRanNum) then
                  RanNumLoc = 1
                end if
                if (RanNum(RanNumLoc) .lt. 0.5) then
                  ! Multiply by product to boost better that are less individuals (note the - in front!)
                  InitChrom(1:Data%nPotPar, iSol) = - InitChrom(1:Data%nPotPar, iSol) * SelCriterionStd * AvgCoancestryStd
                else
                  ! Multiply by standardized average coancestry to boost less related individuals
                  InitChrom(1:Data%nPotPar, iSol) =   InitChrom(1:Data%nPotPar, iSol) * AvgCoancestryStd
                end if
              end if
            end if
          end do

          ! Search
          if (trim(Spec%EvolAlg) .eq. "DE") then
            call DifferentialEvolution(Spec=Spec, Data=Data, nParam=nParam, nSol=Spec%EvolAlgNSol, Init=InitChrom, &
              nIter=Spec%EvolAlgNIter, nIterBurnIn=Spec%DiffEvolNIterBurnIn, nIterStop=Spec%EvolAlgNIterStop, StopTolerance=Spec%EvolAlgStopTol, nIterPrint=Spec%EvolAlgNIterPrint, &
              LogStdout=LogStdoutInternal, LogFile=LogFile, LogPop=Spec%EvolAlgLogPop, LogPopFile=LogPopFile, &
              CRBurnIn=Spec%DiffEvolParamCrBurnIn, CRLate1=Spec%DiffEvolParamCr1, CRLate2=Spec%DiffEvolParamCr2, &
              FBase=Spec%DiffEvolParamFBase, FHigh1=Spec%DiffEvolParamFHigh1, FHigh2=Spec%DiffEvolParamFHigh2, &
              Stream=AlphaMateStream, BestSol=Sol)!, Status=OptimOK)
            ! if (.not. OptimOK) then
            !   write(STDERR, "(a)") " ERROR: Optimisation failed!"
            !   write(STDERR, "(a)") " "
            !   stop 1
            ! end if
          end if

          ! Save
          call Sol%Log(Spec=Spec, LogUnit=Unit, Iteration=-1, AcceptPct=real(NANFLOATTYPE), String=trim("ModeOpt"//trim(Int2Char(Target))), StringNum=18)
          call Sol%WriteContributions(Data=Data, Spec=Spec, ContribFile=ContribFile)
          if (Spec%MateAllocation) then
            call Sol%WriteMatingPlan(Data=Data, Spec=Spec, MatingFile=MatingFile)
          end if

        end do

        close(Unit)

      end if

      deallocate(InitChrom)
      deallocate(AvgCoancestryStd)

    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Write solution to a file or standard output
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 25, 2017
    !> @return Output to a file or standard output
    !---------------------------------------------------------------------------
    subroutine WriteAlphaMateSol(This, File) ! not pure due to IO
      implicit none
      class(AlphaMateSol), intent(in)        :: This !< AlphaMateSol holder
      character(len=*), intent(in), optional :: File !< File, if missing use standard output

      integer(int32) :: Unit, Mat
      if (present(File)) then
        open(newunit=Unit, file=File, action="write", status="unknown")
      else
        Unit = STDOUT
      end if

      ! AlphaEvolveSol components
      write(Unit, *) "Objective: ", This%Objective
      write(Unit, *) "nParam: ", This%nParam
      write(Unit, *) "Chrom: ", This%Chrom

      ! AlphaMateSol components
      write(Unit, *) "Penalty: ", This%Penalty
      write(Unit, *) "PenaltyCoancestry: ", This%PenaltyCoancestry
      write(Unit, *) "PenaltyInbreeding: ", This%PenaltyInbreeding
      write(Unit, *) "PenaltyLimitPar1: ", This%PenaltyLimitPar1
      write(Unit, *) "PenaltyLimitPar2: ", This%PenaltyLimitPar2
      write(Unit, *) "PenaltyMultipleMatings: ", This%PenaltyMultipleMatings
      write(Unit, *) "PenaltySelfing: ", This%PenaltySelfing
      write(Unit, *) "PenaltyGenericIndCrit: ", This%PenaltyGenericIndCrit
      write(Unit, *) "PenaltyGenericMatCrit: ", This%PenaltyGenericMatCrit
      write(Unit, *) "Degree: ", This%Degree
      write(Unit, *) "SelCriterion: ", This%SelCriterion
      write(Unit, *) "SelCriterionStd: ", This%SelCriterionStd
      write(Unit, *) "MaxCriterionPct: ", This%MaxCriterionPct
      write(Unit, *) "CoancestryRanMate: ", This%CoancestryRanMate
      write(Unit, *) "CoancestryRateRanMate: ", This%CoancestryRateRanMate
      write(Unit, *) "MinCoancestryPct: ", This%MinCoancestryPct
      write(Unit, *) "Inbreeding: ", This%Inbreeding
      write(Unit, *) "InbreedingRate: ", This%InbreedingRate
      write(Unit, *) "MinInbreedingPct: ", This%MinInbreedingPct
      if (allocated(This%GenericIndCrit)) then
        write(Unit, *) "GenericIndCrit: ", This%GenericIndCrit
      else
        write(Unit, *) "GenericIndCrit: not allocated"
      end if
      if (allocated(This%GenericMatCrit)) then
        write(Unit, *) "GenericMatCrit: ", This%GenericMatCrit
      else
        write(Unit, *) "GenericMatCrit: not allocated"
      end if
      ! write(Unit, *) "Cost: ", This%Cost
      write(Unit, *) "nVec: ", This%nVec
      write(Unit, *) "Mating/crossing plan:"
      do Mat = 1, size(This%MatingPlan, dim=2)
        write(Unit, *) Mat, This%MatingPlan(:, Mat)
      end do
      if (allocated(This%GenomeEdit)) then
        write(Unit, *) "GenomeEdit: ", This%GenomeEdit
      else
        write(Unit, *) "GenomeEdit: not allocated"
      end if

      ! pause

      if (present(File)) then
        close(Unit)
      end if
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Write contributions to a file or standard output
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   April 7, 2017
    !> @return Output to file or standard output
    !---------------------------------------------------------------------------
    subroutine WriteContributions(This, Data, Spec, ContribFile) ! not pure due to IO
      implicit none

      ! Arguments
      class(AlphaMateSol), intent(in)        :: This        !< AlphaMateSol holder
      type(AlphaMateData), intent(in)        :: Data        !< AlphaMateData holder
      type(AlphaMateSpec), intent(in)        :: Spec        !< AlphaMateSpec holder
      character(len=*), intent(in), optional :: ContribFile !< File to write individual contributions to (default STDOUT)

      ! Other
      integer(int32) :: ContribUnit, Rank(Data%nInd), nCon, i, Ind

      if (.not. present(ContribFile)) then
        ContribUnit = STDOUT
      end if

      ! @todo should we have constant output no matter which options are switched on?
      if (present(ContribFile)) then
        open(newunit=ContribUnit, file=ContribFile, status="unknown")
      end if

      ! Find contributors
      nCon = 0
      do Ind = 1, Data%nInd
        if (This%nVec(Ind) .gt. 0) then
          nCon = nCon + 1
        end if
      end do

      ! Rank them by number of contributions
      Rank(1:nCon) = RapKnr(float(This%nVec), nCon) ! @todo float(This%nVec) due to RapKnr bug with integers
      ! call This%Write
      ! print*, This%nVec
      ! print*, Rank(1:nCon)
      ! print*, Rank
      ! print*, This%nVec(Rank(1:nCon))

      ! Write them out
      if (.not. allocated(This%GenomeEdit)) then
        !                                             12345678901234567890123456789012
        write(ContribUnit, Spec%FmtContributionHead) "                              Id", &
                                                     "         Gender", &
                                                     "   SelCriterion", &
                                                     "  AvgCoancestry", &
                                                     "   Contribution", &
                                                     "        nMating"
        if (.not. Spec%GenderGiven) then
          do i = 1, nCon
            Ind = Rank(i)
            write(ContribUnit, Spec%FmtContribution) Data%Coancestry%OriginalId(Ind),              &
                                                     Data%Gender(Ind),                             &
                                                     Data%SelCriterion(Ind),                       &
                                                     Data%AvgCoancestry(Ind),                      &
                                                     FLOATFUN(This%nVec(Ind)) / (2.0 * Spec%nMat), &
                                                     This%nVec(Ind)
          end do
        else
          do i = 1, nCon
            Ind = Rank(i)
            if (Data%Gender(Ind) .eq. 1) then
              write(ContribUnit, Spec%FmtContribution) Data%Coancestry%OriginalId(Ind),              &
                                                       Data%Gender(Ind),                             &
                                                       Data%SelCriterion(Ind),                       &
                                                       Data%AvgCoancestry(Ind),                      &
                                                       FLOATFUN(This%nVec(Ind)) / (2.0 * Spec%nMat), &
                                                       This%nVec(Ind)
            end if
          end do
          do i = 1, nCon
            Ind = Rank(i)
            if (Data%Gender(Ind) .eq. 2) then
              write(ContribUnit, Spec%FmtContribution) Data%Coancestry%OriginalId(Ind),              &
                                                       Data%Gender(Ind),                             &
                                                       Data%SelCriterion(Ind),                       &
                                                       Data%AvgCoancestry(Ind),                      &
                                                       FLOATFUN(This%nVec(Ind)) / (2.0 * Spec%nMat), &
                                                       This%nVec(Ind)
            end if
          end do
        end if
      else
        !                                                 12345678901234567890123456789012
        write(ContribUnit, Spec%FmtContributionHeadEdit) "                              Id", &
                                                         "         Gender", &
                                                         "   SelCriterion", &
                                                         "  AvgCoancestry", &
                                                         "   Contribution", &
                                                         "        nMating", &
                                                         "     GenomeEdit", &
                                                         "  EditedSelCrit"
        if (.not. Spec%GenderGiven) then
          do i = 1, nCon
            Ind = Rank(i)
            write(ContribUnit, Spec%FmtContributionEdit) Data%Coancestry%OriginalId(Ind),              &
                                                         Data%Gender(Ind),                             &
                                                         Data%SelCriterion(Ind),                       &
                                                         Data%AvgCoancestry(Ind),                      &
                                                         FLOATFUN(This%nVec(Ind)) / (2.0 * Spec%nMat), &
                                                         This%nVec(Ind),                               &
                                                         nint(This%GenomeEdit(Ind)),                   &
                                                         Data%SelCriterion(Ind) + This%GenomeEdit(Ind) * Data%SelCriterionPAGE(Ind)
          end do
        else
          do i = 1, nCon
            Ind = Rank(i)
            if (Data%Gender(Ind) .eq. 1) then
              write(ContribUnit, Spec%FmtContributionEdit) Data%Coancestry%OriginalId(Ind),              &
                                                           Data%Gender(Ind),                             &
                                                           Data%SelCriterion(Ind),                       &
                                                           Data%AvgCoancestry(Ind),                      &
                                                           FLOATFUN(This%nVec(Ind)) / (2.0 * Spec%nMat), &
                                                           This%nVec(Ind),                               &
                                                           nint(This%GenomeEdit(Ind)),                   &
                                                           Data%SelCriterion(Ind) + This%GenomeEdit(Ind) * Data%SelCriterionPAGE(Ind)
            end if
          end do
          do i = 1, nCon
            Ind = Rank(i)
            if (Data%Gender(Ind) .eq. 2) then
              write(ContribUnit, Spec%FmtContributionEdit) Data%Coancestry%OriginalId(Ind),              &
                                                           Data%Gender(Ind),                             &
                                                           Data%SelCriterion(Ind),                       &
                                                           Data%AvgCoancestry(Ind),                      &
                                                           FLOATFUN(This%nVec(Ind)) / (2.0 * Spec%nMat), &
                                                           This%nVec(Ind),                               &
                                                           nint(This%GenomeEdit(Ind)),                   &
                                                           Data%SelCriterion(Ind) + This%GenomeEdit(Ind) * Data%SelCriterionPAGE(Ind)
            end if
          end do
        end if
      end if

      if (present(ContribFile)) then
        close(ContribUnit)
      end if
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Write mating plan to a file or standard output
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   April 7, 2017
    !> @return Output to file or standard output
    !---------------------------------------------------------------------------
    subroutine WriteMatingPlan(This, Spec, Data, MatingFile) ! not pure due to IO
      implicit none

      ! Arguments
      class(AlphaMateSol), intent(in)        :: This       !< Solution holder
      type(AlphaMateSpec), intent(in)        :: Spec       !< Spec holder
      type(AlphaMateData), intent(in)        :: Data       !< Data holder
      character(len=*), intent(in), optional :: MatingFile !< File to write mating plan to (default STDOUT)

      ! Other
      integer(int32) :: nMat, MatingUnit, k, &
                        Pair(size(This%MatingPlan, dim=2))

      nMat = size(This%MatingPlan, dim=2)

      if (.not. present(MatingFile)) then
        MatingUnit = STDOUT
      end if

      if (present(MatingFile)) then
        open(newunit=MatingUnit, file=MatingFile, status="unknown")
      end if

      !                                      12345678901234567890123456789012
      write(MatingUnit, Spec%FmtMatingHead) "         Mating", &
                                            "                         Parent1", &
                                            "                         Parent2"

      ! Sort such that repeated matings would appear together, but otherwise in no particular order
      Pair = MrgRnk(This%MatingPlan(1, :) * Data%nInd + This%MatingPlan(2, :))

      ! Write out
      k = nMat ! MrgRnk ranks small to large
      do while (k .gt. 0)
        print*, "x", nMat - k + 1, This%MatingPlan(1:2, Pair(k))
        write(MatingUnit, Spec%FmtMating) nMat - k + 1, &
                                          Data%Coancestry%OriginalId(This%MatingPlan(1:2, Pair(k)))
        k = k - 1 ! MrgRnk ranks small to large
      end do

      if (present(MatingFile)) then
        close(MatingUnit)
      end if
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Setup colnames and formats for output
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !---------------------------------------------------------------------------
    pure subroutine SetupColNamesAndFormats(This)
      implicit none
      class(AlphaMateSpec), intent(inout) :: This !< @return Spec holder

      integer(int32) :: nCol, nColTmp, i
      character(len=CHARLENGTH) :: Tmp

      ! --- Optimisation log ---

      nCol = 14
      if (This%GenericIndCritGiven) then
        nCol = nCol + This%nGenericIndCrit
      end if
      if (This%GenericMatCritGiven) then
        nCol = nCol + This%nGenericMatCrit
      end if
      if (.not. allocated(This%ColnameLogUnit)) then
        allocate(This%ColnameLogUnit(nCol))
      end if
      if (.not. allocated(This%ColnameLogStdout)) then
        allocate(This%ColnameLogStdout(nCol))
      end if
      if (.not. allocated(This%ColnameLogPopUnit)) then
        allocate(This%ColnameLogPopUnit(nCol))
      end if

      !                          12345678901234567
      This%ColnameLogUnit(1)  = "        Iteration"
      This%ColnameLogUnit(2)  = "        AcceptPct"
      This%ColnameLogUnit(3)  = "        Objective"
      This%ColnameLogUnit(4)  = "        Penalties"
      This%ColnameLogUnit(5)  = "   FrontierDegree"
      This%ColnameLogUnit(6)  = "     SelCriterion"
      This%ColnameLogUnit(7)  = "  SelCriterionStd"
      This%ColnameLogUnit(8)  = "  MaxSelCriterPct"
      This%ColnameLogUnit(9)  = "       Coancestry"
      This%ColnameLogUnit(10) = "   CoancestryRate"
      This%ColnameLogUnit(11) = " MinCoancestryPct"
      This%ColnameLogUnit(12) = "       Inbreeding"
      This%ColnameLogUnit(13) = "   InbreedingRate"
      This%ColnameLogUnit(14) = " MinInbreedingPct"

      !                            12345678901
      This%ColnameLogStdout(1)  = "  Iteration"
      This%ColnameLogStdout(2)  = "  AcceptPct"
      This%ColnameLogStdout(3)  = "  Objective"
      This%ColnameLogStdout(4)  = "  Penalties"
      This%ColnameLogStdout(5)  = "     Degree"
      This%ColnameLogStdout(6)  = "  SelCriter"
      This%ColnameLogStdout(7)  = "     ...Std"
      This%ColnameLogStdout(8)  = "     ...Pct"
      This%ColnameLogStdout(9)  = " Coancestry"
      This%ColnameLogStdout(10) = "    ...Rate"
      This%ColnameLogStdout(11) = "     ...Pct"
      This%ColnameLogStdout(12) = " Inbreeding"
      This%ColnameLogStdout(13) = "    ...Rate"
      This%ColnameLogStdout(14) = "     ...Pct"

      nColTmp = nCol
      if (This%GenericIndCritGiven) then
        do i = 1, This%nGenericIndCrit
          nColTmp = nColTmp + 1
          !                               1234567890
          This%ColnameLogUnit(nColTmp) = "GenIndCrit"//trim(Int2Char(i))
          This%ColnameLogUnit(nColTmp) = adjustr(This%ColnameLogUnit(nColTmp))
          Tmp = This%ColnameLogUnit(nColTmp)
          This%ColnameLogStdout(nColTmp) = Tmp(7:17)
          This%ColnameLogStdout(nColTmp) = adjustr(This%ColnameLogStdout(nColTmp))
        end do
      end if
      if (This%GenericMatCritGiven) then
        do i = 1, This%nGenericMatCrit
          nColTmp = nColTmp + 1
          !                               1234567890
          This%ColnameLogUnit(nColTmp) = "GenMatCrit"//trim(Int2Char(i))
          This%ColnameLogUnit(nColTmp) = adjustr(This%ColnameLogUnit(nColTmp))
          Tmp = This%ColnameLogUnit(nColTmp)
          This%ColnameLogStdout(nColTmp) = Tmp(7:17)
          This%ColnameLogStdout(nColTmp) = adjustr(This%ColnameLogStdout(nColTmp))
        end do
      end if
      This%ColnameLogPopUnit = This%ColnameLogUnit
      !                            12345678901234567
      This%ColnameLogPopUnit(2) = "         Solution"
      This%FmtLogStdoutHead  = trim(FMTLOGSTDOUTHEADA) //trim(Int2Char(nCol))        //trim(FMTLOGSTDOUTHEADB)
      This%FmtLogStdout      = trim(FMTLOGSTDOUTA)     //trim(Int2Char(nColTmp-nCol))//trim(FMTLOGSTDOUTB)
      This%FmtLogUnitHead    = trim(FMTLOGUNITHEADA)   //trim(Int2Char(nCol)  )      //trim(FMTLOGUNITHEADB)
      This%FmtLogUnit        = trim(FMTLOGUNITA)       //trim(Int2Char(nCol-1))      //trim(FMTLOGUNITB)
      This%FmtLogPopUnit     = trim(FMTLOGPOPUNITA)    //trim(Int2Char(nCol-2))      //trim(FMTLOGUNITB)

      ! --- Contributions output ---

      This%FmtContributionHead     = trim(FMTCONTRIBUTIONHEADA)    //trim(Int2Char(IDLENGTH))//trim(FMTCONTRIBUTIONHEADB)
      This%FmtContributionHeadEdit = trim(FMTCONTRIBUTIONHEADEDITA)//trim(Int2Char(IDLENGTH))//trim(FMTCONTRIBUTIONHEADEDITB)

      This%FmtContribution         = trim(FMTCONTRIBUTIONA)        //trim(Int2Char(IDLENGTH))//trim(FMTCONTRIBUTIONB)
      This%FmtContributionEdit     = trim(FMTCONTRIBUTIONEDITA)    //trim(Int2Char(IDLENGTH))//trim(FMTCONTRIBUTIONEDITB)

      ! --- Mating output ---

      This%FmtMatingHead = trim(FMTMATINGHEADA)//trim(Int2Char(IDLENGTH))    //trim(FMTMATINGHEADB)
      This%FmtMating     = trim(FMTMATINGA)    //trim(Int2Char(IDLENGTH - 1))//trim(FMTMATINGB)
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Write head of the AlphaMate optimisation log - the best solution
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !> @return Output to file or standard output
    !---------------------------------------------------------------------------
    subroutine LogHeadAlphaMateSpec(This, LogUnit, String, StringNum) ! not pure due to IO
      implicit none
      class(AlphaMateSpec), intent(in)       :: This      !< AlphaMateSpec holder
      integer(int32), intent(in), optional   :: LogUnit   !< Unit to write to (default STDOUT)
      character(len=*), intent(in), optional :: String    !< Additional string that will be written before the head
      integer(int32), optional               :: StringNum !< How much space is needed for the String

      character(len=CHARLENGTH) :: StringFmt

      if (present(String)) then
        if (present(StringNum)) then
          StringFmt = "(a"//trim(Int2Char(StringNum))//")"
        else
          StringFmt = "(a)"
        end if
      end if
      if (present(LogUnit)) then
        if (present(String)) then
          write(LogUnit, StringFmt, Advance="No") adjustl(String)
        end if
        write(LogUnit, This%FmtLogUnitHead)  This%ColnameLogUnit
      else
        if (present(String)) then
          write(STDOUT, StringFmt, Advance="No") adjustl(String)
        end if
        write(STDOUT, This%FmtLogStdoutHead) This%ColnameLogStdout
      end if
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Write the AlphaMate optimisation log - the best solution
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !> @return Output to file or standard output
    !---------------------------------------------------------------------------
    subroutine LogAlphaMateSol(This, Spec, LogUnit, Iteration, AcceptPct, String, StringNum) ! not pure due to IO
      implicit none
      class(AlphaMateSol), intent(in)        :: This      !< Solution holder
      class(AlphaEvolveSpec), intent(in)     :: Spec      !< Spec holder
      integer(int32), intent(in), optional   :: LogUnit   !< Unit to write to (default STDOUT)
      integer(int32), intent(in)             :: Iteration !< Generation/Iteration
      real(real32), intent(in)               :: AcceptPct !< Acceptance rate
      character(len=*), intent(in), optional :: String    !< Additional string that will be written before the head
      integer(int32), optional               :: StringNum !< How much space is needed for the String

      integer(int32) :: Unit
      character(len=CHARLENGTH) :: Fmt, StringFmt

      select type (Spec)
        class default
          error stop " ERROR: LogAlphaMateSol works only with argument Spec being of type AlphaMateSpec!"
        class is (AlphaMateSpec)
          if (present(LogUnit)) then
            Unit = LogUnit
            Fmt = Spec%FmtLogUnit
          else
            Unit = STDOUT
            Fmt = Spec%FmtLogStdout
          end if
          if (present(String)) then
            if (present(StringNum)) then
              StringFmt = "(a"//trim(Int2Char(StringNum))//")"
            else
              StringFmt = "(a)"
            end if
          end if
          if (allocated(This%GenericIndCrit)) then
            if (allocated(This%GenericMatCrit)) then
              if (present(String)) then
                write(Unit, StringFmt, Advance="No") trim(adjustl(String))
              end if
              write(Unit, Fmt) Iteration, AcceptPct, This%Objective, This%Penalty, This%Degree, &
                               This%SelCriterion, This%SelCriterionStd, This%MaxCriterionPct, &
                               This%CoancestryRanMate, This%CoancestryRateRanMate, This%MinCoancestryPct, &
                               This%Inbreeding, This%InbreedingRate, This%MinInbreedingPct, &
                               This%GenericIndCrit, This%GenericMatCrit
            else
              if (present(String)) then
                write(Unit, StringFmt, Advance="No") trim(adjustl(String))
              end if
              write(Unit, Fmt) Iteration, AcceptPct, This%Objective, This%Penalty, This%Degree, &
                               This%SelCriterion, This%SelCriterionStd, This%MaxCriterionPct, &
                               This%CoancestryRanMate, This%CoancestryRateRanMate, This%MinCoancestryPct, &
                               This%Inbreeding, This%InbreedingRate, This%MinInbreedingPct, &
                               This%GenericIndCrit
            end if
          else
            if (allocated(This%GenericMatCrit)) then
              if (present(String)) then
                write(Unit, StringFmt, Advance="No") trim(adjustl(String))
              end if
              write(Unit, Fmt) Iteration, AcceptPct, This%Objective, This%Penalty, This%Degree, &
                               This%SelCriterion, This%SelCriterionStd, This%MaxCriterionPct, &
                               This%CoancestryRanMate, This%CoancestryRateRanMate, This%MinCoancestryPct, &
                               This%Inbreeding, This%InbreedingRate, This%MinInbreedingPct, &
                                                    This%GenericMatCrit
            else
              if (present(String)) then
                write(Unit, StringFmt, Advance="No") trim(adjustl(String))
              end if
              write(Unit, Fmt) Iteration, AcceptPct, This%Objective, This%Penalty, This%Degree, &
                               This%SelCriterion, This%SelCriterionStd, This%MaxCriterionPct, &
                               This%CoancestryRanMate, This%CoancestryRateRanMate, This%MinCoancestryPct, &
                               This%Inbreeding, This%InbreedingRate, This%MinInbreedingPct
            end if
          end if
      end select
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Write head of the AlphaMate optimisation log - all solutions
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !> @return Output to file or standard output
    !---------------------------------------------------------------------------
    subroutine LogPopHeadAlphaMateSpec(This, LogPopUnit) ! not pure due to IO
      implicit none
      class(AlphaMateSpec), intent(in)     :: This       !< Spec holder
      integer(int32), intent(in), optional :: LogPopUnit !< Log file unit (default STDOUT)
      integer(int32) :: Unit

      if (present(LogPopUnit)) then
        Unit = LogPopUnit
      else
        Unit = STDOUT
      end if
      write(Unit, This%FmtLogUnitHead) This%ColnameLogPopUnit
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Write the AlphaMate optimisation log - all solutions
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !> @return Output to file or standard output
    !---------------------------------------------------------------------------
    subroutine LogPopAlphaMateSol(This, Spec, LogPopUnit, Iteration, i) ! not pure due to IO
      implicit none
      class(AlphaMateSol), intent(in)      :: This       !< Solution holder
      class(AlphaEvolveSpec), intent(in)   :: Spec       !< Spec holder
      integer(int32), intent(in), optional :: LogPopUnit !< population log file unit (default STDOUT)
      integer(int32), intent(in)           :: Iteration  !< generation/iteration
      integer(int32), intent(in)           :: i          !< solution id
      integer(int32) :: Unit

      select type (Spec)
        class default
          error stop " ERROR: LogPopAlphaMateSol works only with argument Spec being of type AlphaMateSpec!"
        class is (AlphaMateSpec)
          if (present(LogPopUnit)) then
            Unit = LogPopUnit
          else
            Unit = STDOUT
          end if
          if (allocated(This%GenericIndCrit)) then
            if (allocated(This%GenericMatCrit)) then
              write(Unit, Spec%FmtLogPopUnit) Iteration, i, This%Objective, This%Penalty, This%Degree, &
                                              This%SelCriterion, This%SelCriterionStd, This%MaxCriterionPct, &
                                              This%CoancestryRanMate, This%CoancestryRateRanMate, This%MinCoancestryPct, &
                                              This%Inbreeding, This%InbreedingRate, This%MinInbreedingPct, &
                                              This%GenericIndCrit, This%GenericMatCrit
            else
              write(Unit, Spec%FmtLogPopUnit) Iteration, i, This%Objective, This%Penalty, This%Degree, &
                                              This%SelCriterion, This%SelCriterionStd, This%MaxCriterionPct, &
                                              This%CoancestryRanMate, This%CoancestryRateRanMate, This%MinCoancestryPct, &
                                              This%Inbreeding, This%InbreedingRate, This%MinInbreedingPct, &
                                              This%GenericIndCrit
            end if
          else
            if (allocated(This%GenericMatCrit)) then
              write(Unit, Spec%FmtLogPopUnit) Iteration, i, This%Objective, This%Penalty, This%Degree, &
                                              This%SelCriterion, This%SelCriterionStd, This%MaxCriterionPct, &
                                              This%CoancestryRanMate, This%CoancestryRateRanMate, This%MinCoancestryPct, &
                                              This%Inbreeding, This%InbreedingRate, This%MinInbreedingPct, &
                                                                  This%GenericMatCrit
            else
              write(Unit, Spec%FmtLogPopUnit) Iteration, i, This%Objective, This%Penalty, This%Degree, &
                                              This%SelCriterion, This%SelCriterionStd, This%MaxCriterionPct, &
                                              This%CoancestryRanMate, This%CoancestryRateRanMate, This%MinCoancestryPct, &
                                              This%Inbreeding, This%InbreedingRate, This%MinInbreedingPct
            end if
          end if
      end select
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert selection criterion to stand. selection criterion
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    April 18, 2017
    !---------------------------------------------------------------------------
    elemental function SelCriterion2SelCriterionStd(SelCriterion, Mean, Sd) result(SelCriterionStd)
      implicit none
      real(FLOATTYPE), intent(in) :: SelCriterion    !< selection criterion
      real(FLOATTYPE), intent(in) :: Mean            !< mean of selection criterion
      real(FLOATTYPE), intent(in) :: Sd              !< standard deviation of selection criterion
      real(FLOATTYPE)             :: SelCriterionStd !< @return stand. selection criterion
      SelCriterionStd = (SelCriterion - Mean) / Sd
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert stand selection criterion to selection criterion
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    April 18, 2017
    !---------------------------------------------------------------------------
    elemental function SelCriterionStd2SelCriterion(SelCriterionStd, Mean, Sd) result(SelCriterion)
      implicit none
      real(FLOATTYPE), intent(in) :: SelCriterionStd !< stand. selection criterion
      real(FLOATTYPE), intent(in) :: Mean            !< mean of selection criterion
      real(FLOATTYPE), intent(in) :: Sd              !< standard deviation of selection criterion
      real(FLOATTYPE)             :: SelCriterion    !< @return selection criterion
      SelCriterion = SelCriterionStd * Sd + Mean
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert current and future coancestry to coancestry rate
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    April 18, 2017
    !> @details dF = (F_t+1 - F_t) / (1 - F_t)
    !---------------------------------------------------------------------------
    pure function Coancestry2CoancestryRate(CurrentCoancestry, FutureCoancestry) result(CoancestryRate)
      implicit none
      real(FLOATTYPE), intent(in) :: CurrentCoancestry !< Current coancestry
      real(FLOATTYPE), intent(in) :: FutureCoancestry  !< Future coancestry
      real(FLOATTYPE)             :: CoancestryRate    !< @return CoancestryRate
      real(FLOATTYPE)             :: Diff, MaxDiff
      Diff    = FutureCoancestry - CurrentCoancestry
      MaxDiff =              1.0 - CurrentCoancestry
      ! @todo What should be done, when we have coancestry estimates that are above 1 or below 1?
      if (MaxDiff .eq. 0) then
        if (Diff .ge. 0) then
          CoancestryRate =  1.0
        else if (Diff .eq. 0) then
          CoancestryRate =  0.0
        else
          CoancestryRate = -1.0
        end if
        ! Not sure about the above fix, but the logic is that if MaxDiff is zero,
        ! then whatever positive, zero, or negative) Diff we get, we respetively
        ! achieve rate of 1, 0, -1
      else
        CoancestryRate = Diff / MaxDiff
      end if
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert coancestry rate and current coancestry to future coancestry
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    April 18, 2017
    !> @details F_t = DeltaF + (1 - DeltaF) * F_t-1
    !---------------------------------------------------------------------------
    pure function CoancestryRate2Coancestry(CoancestryRate, CurrentCoancestry) result(FutureCoancestry)
      implicit none
      real(FLOATTYPE), intent(in) :: CoancestryRate    !< CoancestryRate
      real(FLOATTYPE), intent(in) :: CurrentCoancestry !< Current coancestry
      real(FLOATTYPE)             :: FutureCoancestry  !< @return Future coancestry
      FutureCoancestry = CoancestryRate + (1.0 - CoancestryRate) * CurrentCoancestry
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert MinCoancestryPct to Frontier degree
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    April 13, 2017
    !> @details Assuming unit circular selection/coancestry frontier and given the
    !!          percentage of minimum inbreeding achieved (x-axis) we can evaluate
    !!          the angle between the maximum selection line and the solution line
    !!          recognising that sin(angle) = opposite/hypothenuse, where
    !!          opposite = MinCoancestryPct/100 and hypothenuse = 1 (unit circle). Then
    !!          angle = asin(MinCoancestryPct/100).
    !---------------------------------------------------------------------------
    pure function MinCoancestryPct2Degree(MinCoancestryPct) result(Degree)
      implicit none
      real(FLOATTYPE), intent(in) :: MinCoancestryPct !< percentage of minimum coancestry achieved (100 means we achieved the minimum possible coancestry)
      real(FLOATTYPE)             :: Degree           !< @return Frontier degree
      Degree = asin(MinCoancestryPct / 100.0) * RAD2DEG
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert Frontier degree to MinCoancestryPct
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    April 13, 2017
    !> @details Assuming unit circular selection/coancestry frontier and given the
    !!          the angle between the maximum selection line and the solution line
    !!          we can evaluate the percentage of maximum selection criterion achieved
    !!          (y-axis) by recognising that sin(angle) = opposite/hypothenuse, where
    !!          opposite = MinCoancestryPct/100 and hypothenuse = 1 (unit circle). Then
    !!          MinCoancestryPct = sin(angle) * 100.
    !---------------------------------------------------------------------------
    pure function Degree2MinCoancestryPct(Degree) result(MinCoancestryPct)
      implicit none
      real(FLOATTYPE), intent(in) :: Degree           !< Frontier degree
      real(FLOATTYPE)             :: MinCoancestryPct !< @return Percentage of minimum coancestry achieved (100 means we achieved the minimum possible coancestry)
      MinCoancestryPct = sin(Degree * DEG2RAD) * 100.0
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert MaxCriterionPct to Frontier degree
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    April 13, 2017
    !> @details Assuming unit circular selection/coancestry frontier and given the
    !!          percentage of maximum selection criterion achieved (y-axis) we can evaluate
    !!          the angle between the maximum selection line and the solution line
    !!          recognising that cos(angle) = adjacent/hypothenuse, where
    !!          adjacent = MaxCriterionPct/100 and hypothenuse = 1 (unit circle). Then
    !!          angle = acos(MaxCriterionPct/100).
    !---------------------------------------------------------------------------
    pure function MaxCriterionPct2Degree(MaxCriterionPct) result(Degree)
      implicit none
      real(FLOATTYPE), intent(in) :: MaxCriterionPct !< percentage of maximum criterion achieved (100 means we achieved the maximum possible selection criterion)
      real(FLOATTYPE)             :: Degree          !< @return Frontier degree
      Degree = acos(MaxCriterionPct / 100.0) * RAD2DEG
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert Frontier degree to MaxCriterionPct
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    April 13, 2017
    !> @details Assuming unit circular selection/coancestry frontier and given the
    !!          the angle between the maximum selection line and the solution line
    !!          we can evaluate the percentage of maximum selection criterion achieved
    !!          (y-axis) by recognising that cos(angle) = adjacent/hypothenuse, where
    !!          adjacent = MaxCriterionPct/100 and hypothenuse = 1 (unit circle). Then
    !!          MaxCriterionPct = cos(angle) * 100.
    !---------------------------------------------------------------------------
    pure function Degree2MaxCriterionPct(Degree) result(MaxCriterionPct)
      implicit none
      real(FLOATTYPE), intent(in) :: Degree          !< Frontier degree
      real(FLOATTYPE)             :: MaxCriterionPct !< @return Percentage of maximum criterion achieved (100 means we achieved the maximum possible selection criterion)
      MaxCriterionPct = cos(Degree * DEG2RAD) * 100.0
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert MinCoancestryPct to coancestry rate
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    April 13, 2017
    !> @details Compute difference between Max and Min and take (1 - MinCoancestryPct) of
    !!          the difference as we express MinCoancestryPct as 100% when solution achieves
    !!          Min coancestry
    !---------------------------------------------------------------------------
    pure function MinCoancestryPct2CoancestryRate(MinCoancestryPct, MinCoancestryRate, MaxCoancestryRate) result(CoancestryRate)
      implicit none
      real(FLOATTYPE), intent(in) :: MinCoancestryPct  !< MinCoancestryPct of a solution
      real(FLOATTYPE), intent(in) :: MinCoancestryRate !< Minimum possible coancestry rate
      real(FLOATTYPE), intent(in) :: MaxCoancestryRate !< Maximum possible coancestry rate
      real(FLOATTYPE)             :: CoancestryRate    !< @return Coancestry rate at a given MinCoancestryPct
      CoancestryRate = MinCoancestryRate + (100.0 - MinCoancestryPct) / 100.0 * (MaxCoancestryRate - MinCoancestryRate)
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert coancestry rate to MinCoancestryPct
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    April 18, 2017
    !> @details Compute difference between Max and Min and take (1 - MinCoancestryPct) of
    !!          the difference as we express MinCoancestryPct as 100% when solution achieves
    !!          Min coancestry
    !---------------------------------------------------------------------------
    pure function CoancestryRate2MinCoancestryPct(CoancestryRate, MinCoancestryRate, MaxCoancestryRate) result(MinCoancestryPct)
      implicit none
      real(FLOATTYPE), intent(in) :: CoancestryRate    !< Coancestry rate at a given MinCoancestryPct
      real(FLOATTYPE), intent(in) :: MinCoancestryRate !< Minimum possible coancestry rate
      real(FLOATTYPE), intent(in) :: MaxCoancestryRate !< Maximum possible coancestry rate
      real(FLOATTYPE)             :: MinCoancestryPct  !< @return MinCoancestryPct of a solution
      real(FLOATTYPE)             :: Diff, MaxDiff
      Diff    = MaxCoancestryRate - CoancestryRate
      MaxDiff = MaxCoancestryRate - MinCoancestryRate
      if (MaxDiff .eq. 0) then
        if (Diff .ge. 0) then
          MinCoancestryPct = 100.0
        else
          MinCoancestryPct = 0.0
        end if
        ! Not sure about the above fix, but the logic is that if MaxDiff is zero,
        ! then whatever positive (or negative) Diff we get, we achieve 100% (or 0%).
      else
        MinCoancestryPct = Diff / MaxDiff * 100.0
      end if
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert MaxCriterionPct to stand. selection criterion
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    April 13, 2017
    !> @details Compute difference between Max and Min and take MaxCriterionPct of
    !!          the difference as we express MaxCriterionPct as 100% when solution
    !!          achieves Max criterion
    !---------------------------------------------------------------------------
    pure function MaxCriterionPct2SelCriterionStd(MaxCriterionPct, MinSelCriterionStd, MaxSelCriterionStd) result(SelCriterionStd)
      implicit none
      real(FLOATTYPE), intent(in) :: MaxCriterionPct    !< MaxCriterionPct of a solution
      real(FLOATTYPE), intent(in) :: MinSelCriterionStd !< Minimum possible stand. selection criterion
      real(FLOATTYPE), intent(in) :: MaxSelCriterionStd !< Maximum possible stand. selection criterion
      real(FLOATTYPE)             :: SelCriterionStd    !< @return stand. selection criterion at a given MaxCriterionPct
      SelCriterionStd = MinSelCriterionStd + MaxCriterionPct / 100.0 * (MaxSelCriterionStd - MinSelCriterionStd)
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert stand. selection criterion to MaxCriterionPct
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    April 18, 2017
    !> @details See how close SelCriterionStd is to MaxSelCriterionStd
    !---------------------------------------------------------------------------
    pure function SelCriterionStd2MaxCriterionPct(SelCriterionStd, MinSelCriterionStd, MaxSelCriterionStd) result(MaxCriterionPct)
      implicit none
      real(FLOATTYPE), intent(in) :: SelCriterionStd    !< Stand. selection criterion
      real(FLOATTYPE), intent(in) :: MinSelCriterionStd !< Minimum possible stand. selection criterion
      real(FLOATTYPE), intent(in) :: MaxSelCriterionStd !< Maximum possible stand. selection criterion
      real(FLOATTYPE)             :: MaxCriterionPct    !< @return MaxCriterionPct
      real(FLOATTYPE)             :: Diff, MaxDiff
      Diff    =    SelCriterionStd - MinSelCriterionStd
      MaxDiff = MaxSelCriterionStd - MinSelCriterionStd
      if (MaxDiff .eq. 0) then
        if (Diff .ge. 0) then
          MaxCriterionPct = 100.0
        else
          MaxCriterionPct = 0.0
        end if
        ! Not sure about the above fix, but the logic is that if MaxDiff is zero,
        ! then whatever positive (or negative) Diff we get, we achieve 100% (or 0%).
      else
        MaxCriterionPct = Diff / MaxDiff * 100.0
      end if
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert degree to stand. selection criterion
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    April 18, 2017
    !---------------------------------------------------------------------------
    pure function Degree2SelCriterionStd(Degree, MinSelCriterionStd, MaxSelCriterionStd) result(SelCriterionStd)
      implicit none
      real(FLOATTYPE), intent(in) :: Degree             !< Degree
      real(FLOATTYPE), intent(in) :: MinSelCriterionStd !< Minimum possible stand. selection criterion
      real(FLOATTYPE), intent(in) :: MaxSelCriterionStd !< Maximum possible stand. selection criterion
      real(FLOATTYPE)             :: SelCriterionStd    !< @return stand. selection criterion
      SelCriterionStd = MaxCriterionPct2SelCriterionStd(MaxCriterionPct=Degree2MaxCriterionPct(Degree=Degree), &
                                                        MinSelCriterionStd=MinSelCriterionStd, &
                                                        MaxSelCriterionStd=MaxSelCriterionStd)
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert stand. selection criterion to degree
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    April 18, 2017
    !---------------------------------------------------------------------------
    pure function SelCriterionStd2Degree(SelCriterionStd, MinSelCriterionStd, MaxSelCriterionStd) result(Degree)
      implicit none
      real(FLOATTYPE), intent(in) :: SelCriterionStd    !< Stand. selection criterion
      real(FLOATTYPE), intent(in) :: MinSelCriterionStd !< Minimum possible stand. selection criterion
      real(FLOATTYPE), intent(in) :: MaxSelCriterionStd !< Maximum possible stand. selection criterion
      real(FLOATTYPE)             :: Degree          !< @return Degree
      Degree = MaxCriterionPct2Degree(MaxCriterionPct=SelCriterionStd2MaxCriterionPct(SelCriterionStd=SelCriterionStd, &
                                                                                      MinSelCriterionStd=MinSelCriterionStd, &
                                                                                      MaxSelCriterionStd=MaxSelCriterionStd))
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert degree to coancestry rate
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    April 18, 2017
    !---------------------------------------------------------------------------
    pure function Degree2CoancestryRate(Degree, MinCoancestryRate, MaxCoancestryRate) result(CoancestryRate)
      implicit none
      real(FLOATTYPE), intent(in) :: Degree            !< Degree
      real(FLOATTYPE), intent(in) :: MinCoancestryRate !< Minimum possible coancestry rate
      real(FLOATTYPE), intent(in) :: MaxCoancestryRate !< Maximum possible coancestry rate
      real(FLOATTYPE)             :: CoancestryRate    !< @return Coancestry rate
      CoancestryRate = MinCoancestryPct2CoancestryRate(MinCoancestryPct=Degree2MinCoancestryPct(Degree=Degree), &
                                                       MinCoancestryRate=MinCoancestryRate, &
                                                       MaxCoancestryRate=MaxCoancestryRate)
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert coancestry rate to degree
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    April 18, 2017
    !---------------------------------------------------------------------------
    pure function CoancestryRate2Degree(CoancestryRate, MinCoancestryRate, MaxCoancestryRate) result(Degree)
      implicit none
      real(FLOATTYPE), intent(in) :: CoancestryRate    !< Coancestry rate
      real(FLOATTYPE), intent(in) :: MinCoancestryRate !< Minimum possible coancestry rate
      real(FLOATTYPE), intent(in) :: MaxCoancestryRate !< Maximum possible coancestry rate
      real(FLOATTYPE)             :: Degree            !< @return Degree
      Degree = MinCoancestryPct2Degree(MinCoancestryPct=CoancestryRate2MinCoancestryPct(CoancestryRate=CoancestryRate, &
                                                                                        MinCoancestryRate=MinCoancestryRate, &
                                                                                        MaxCoancestryRate=MaxCoancestryRate))
    end function

    !###########################################################################

end module

!###############################################################################
