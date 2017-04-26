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
!> @details  Optimize contributions or individuals to the next generation and
!!           generate a mating plan
!
!> @author   Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
!
!> @date     March 9, 2017
!
!> @version  0.1.0 (alpha)
!
!-------------------------------------------------------------------------------
module AlphaMateModule
  use ISO_Fortran_Env, STDIN => input_unit, STDOUT => output_unit, STDERR => error_unit
  use, intrinsic :: IEEE_Arithmetic
  use ConstantModule, only : FILELENGTH, SPECOPTIONLENGTH, IDLENGTH, RAD2DEG, DEG2RAD
  use OrderPackModule, only : MrgRnk
  use AlphaHouseMod, only : CountLines, Char2Int, Char2Double, Int2Char, Real2Char, &
                            RandomOrder, SetSeed, ToLower, FindLoc, &
                            ParseToFirstWhitespace, SplitLineIntoTwoParts
  use AlphaStatMod, only : Mean, DescStat, DescStatReal64, &
                           DescStatMatrix, DescStatLowTriMatrix, DescStatMatrixReal64
  use AlphaEvolveModule, only : AlphaEvolveSol, AlphaEvolveSpec, AlphaEvolveData, &
                                DifferentialEvolution, RandomSearch
  use AlphaRelateModule

  implicit none

  private
  ! Types
  public :: AlphaMateSpec, AlphaMateData, AlphaMateSol
  ! Functions
  public :: AlphaMateTitle, AlphaMateSearch
  public :: MinCoancestryPct2Degree, Degree2MinCoancestryPct, MaxCriterionPct2Degree, Degree2MaxCriterionPct
  public :: Degree2SelIntensity, SelIntensity2Degree
  public :: Degree2CoancestryRate, CoancestryRate2Degree
  public :: MinCoancestryPct2CoancestryRate, CoancestryRate2MinCoancestryPct
  public :: MaxCriterionPct2SelIntensity, SelIntensity2MaxCriterionPct
  public :: Coancestry2CoancestryRate, CoancestryRate2Coancestry
  public :: SelCriterion2SelIntensity, SelIntensity2SelCriterion

  ! Module parameters
  REAL(real64), PARAMETER :: TARGETDEGREEFRONTIER(8) = [80, 70, 60, 50, 40, 30, 20, 10]

  INTEGER,                   PARAMETER :: CHARLENGTH = 100
  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTREAL2CHAR = "(f11.5)"
  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTINT2CHAR  = "(i11)"

  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTLOGSTDOUTHEADA = "("
  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTLOGSTDOUTHEADB = "a11)"

  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTLOGSTDOUTA = "(i11, f11.1, 2(f11.5), f11.1, 3(2(f11.5), f11.1), "
  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTLOGSTDOUTB = "(f11.5))"

  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTLOGUNITHEADA = "("
  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTLOGUNITHEADB = "a22)"

  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTLOGUNITA = "(i22, "
  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTLOGUNITB = "(1x, es21.13e3))"

  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTLOGPOPUNITA = "(2i22, "

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

  !> @brief Optimisation mode specifications
  type :: AlphaMateModeSpec
    character(len=SPECOPTIONLENGTH) :: Name
    logical                         :: ObjectiveCriterion
    logical                         :: ObjectiveCoancestry
    logical                         :: ObjectiveInbreeding
    real(real64)                    :: TargetDegree
    real(real64)                    :: TargetSelCriterion
    real(real64)                    :: TargetSelIntensity
    real(real64)                    :: TargetMaxCriterionPct
    real(real64)                    :: TargetCoancestry
    real(real64)                    :: TargetCoancestryRate
    logical                         :: CoancestryWeightBelow
    real(real64)                    :: TargetMinCoancestryPct
    real(real64)                    :: TargetInbreeding
    real(real64)                    :: TargetInbreedingRate
    logical                         :: InbreedingWeightBelow
    real(real64)                    :: TargetMinInbreedingPct
    real(real64)                    :: Degree
    real(real64)                    :: SelCriterion
    real(real64)                    :: SelIntensity
    real(real64)                    :: MaxCriterionPct
    real(real64)                    :: Coancestry
    real(real64)                    :: CoancestryRate
    real(real64)                    :: MinCoancestryPct
    real(real64)                    :: Inbreeding
    real(real64)                    :: InbreedingRate
    real(real64)                    :: MinInbreedingPct
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
    character(len=FILELENGTH) :: SpecFile, RelMtxFile, SelCriterionFile, GenderFile, SeedFile
    character(len=FILELENGTH) :: GenericIndCritFile, GenericMatCritFile
    character(len=FILELENGTH) :: OutputBasename
    logical :: RelMtxGiven, SelCriterionGiven, GenderGiven, SeedFileGiven, GenericIndCritGiven, GenericMatCritGiven
    integer(int32) :: Seed
    logical :: SeedGiven

    ! Search mode specifications
    logical :: ModeMinCoancestry, ModeMinInbreeding, ModeMaxCriterion, ModeOpt, ModeRan, EvaluateFrontier
    integer(int32) :: nTargets

    ! Biological specifications
    logical :: NrmInsteadOfCoancestry
    logical :: TargetDegreeGiven
    logical :: TargetSelCriterionGiven, TargetSelIntensityGiven, TargetMaxCriterionPctGiven
    logical :: TargetCoancestryGiven, TargetCoancestryRateGiven, TargetMinCoancestryPctGiven
    logical :: TargetInbreedingGiven, TargetInbreedingRateGiven, TargetMinInbreedingPctGiven
    real(real64), allocatable :: TargetDegree(:)
    real(real64), allocatable :: TargetSelCriterion(:), TargetSelIntensity(:), TargetMaxCriterionPct(:)
    real(real64), allocatable :: TargetCoancestry(:), TargetCoancestryRate(:), TargetMinCoancestryPct(:)
    real(real64) :: TargetInbreeding, TargetInbreedingRate, TargetMinInbreedingPct
    real(real64) :: CoancestryWeight, InbreedingWeight, SelfingWeight
    logical :: SelfingAllowed, CoancestryWeightBelow, InbreedingWeightBelow
    integer(int32) :: nInd, nMat, nPar, nPar1, nPar2 ! NOTE: nInd is here just for OO-flexibility (do not use it; the main one is in Data!!!)
    logical :: EqualizePar, EqualizePar1, EqualizePar2, LimitPar, LimitPar1, LimitPar2
    real(real64) :: LimitParMin, LimitPar1Min, LimitPar2Min, LimitParMax, LimitPar1Max, LimitPar2Max, LimitParMinWeight, LimitPar1MinWeight, LimitPar2MinWeight
    integer(int32) :: nGenericIndCrit, nGenericMatCrit
    real(real64), allocatable :: GenericIndCritWeight(:), GenericMatCritWeight(:)
    logical :: PAGEPar, PAGEPar1, PAGEPar2
    integer(int32) :: PAGEParMax, PAGEPar1Max, PAGEPar2Max
    real(real64) :: PAGEParCost, PAGEPar1Cost, PAGEPar2Cost

    ! Algorithm specifications
    ! ... generic evolutionary parameters
    integer(int32) :: EvolAlgNSol, EvolAlgNIter, EvolAlgNIterStop, EvolAlgNIterPrint
    real(real64) :: EvolAlgStopTol
    logical :: EvolAlgLogPop
    character(len=SPECOPTIONLENGTH) :: EvolAlg
    ! ... differential evolution
    integer(int32) :: DiffEvolNIterBurnIn
    real(real64) :: DiffEvolParamCrBurnIn, DiffEvolParamCr, DiffEvolParamFBase, DiffEvolParamFHigh1, DiffEvolParamFHigh2
    ! ... random search
    integer(int32) :: RanAlgStricter

    ! Modes
    !@todo do we need ModeRanSpec?
    type(AlphaMateModeSpec) :: ModeSpec, ModeMinCoancestrySpec, ModeMinInbreedingSpec, ModeMaxCriterionSpec, ModeRanSpec, ModeOptSpec

    contains
      procedure :: Initialise => InitialiseAlphaMateSpec
      procedure :: Read       => ReadAlphaMateSpec
      procedure :: Write      => WriteAlphaMateSpec
      procedure :: SetupMode  => SetupModeAlphaMateSpec
  end type

  !> @brief AlphaMate data
  type, extends(AlphaEvolveData) :: AlphaMateData
    ! Raw data
    type(RelMat) :: Coancestry
    real(real64), allocatable :: SelCriterion(:), SelIntensity(:), SelCriterionPAGE(:), SelIntensityPAGE(:)
    integer(int32), allocatable :: Gender(:)
    real(real64), allocatable :: GenericIndCrit(:, :), GenericMatCrit(:, :, :)
    ! Data summaries
    type(DescStatReal64) :: InbreedingStat, SelCriterionStat, SelCriterionPAGEStat
    type(DescStatReal64), allocatable :: GenericIndCritStat(:)
    type(DescStatMatrixReal64) :: CoancestryStat, CoancestryStatGender1, CoancestryStatGender2, CoancestryStatGenderDiff
    type(DescStatMatrixReal64), allocatable :: GenericMatCritStat(:)
    ! Derived data
    integer(int32) :: nInd, nPotMat, nPotPar1, nPotPar2, nMal, nFem
    integer(int32), allocatable :: IdPotPar1(:), IdPotPar2(:), IdPotParSeq(:)
    real(real64) :: CoancestryRanMate, CoancestryRanMateNoSelf, CoancestryGenderMate
    real(real64) :: Inbreeding
    contains
      procedure :: Read  => ReadAlphaMateData
      procedure :: Write => WriteAlphaMateData
  end type

  !> @brief AlphaMate solution
  type, extends(AlphaEvolveSol) :: AlphaMateSol
    ! Solution results
    real(real64)                :: Penalty
    real(real64)                :: PenaltyCoancestry
    real(real64)                :: PenaltyInbreeding
    real(real64)                :: PenaltySelfing
    real(real64)                :: PenaltyLimitPar1
    real(real64)                :: PenaltyLimitPar2
    real(real64)                :: PenaltyGenericIndCrit
    real(real64)                :: PenaltyGenericMatCrit
    real(real64)                :: Degree
    real(real64)                :: SelCriterion
    real(real64)                :: SelIntensity
    real(real64)                :: MaxCriterionPct
    real(real64)                :: CoancestryRanMate
    real(real64)                :: CoancestryRateRanMate
    real(real64)                :: MinCoancestryPct
    real(real64)                :: Inbreeding
    real(real64)                :: InbreedingRate
    real(real64)                :: MinInbreedingPct
    real(real64), allocatable   :: GenericIndCrit(:)
    real(real64), allocatable   :: GenericMatCrit(:)
    real(real64)                :: Cost
    integer(int32), allocatable :: nVec(:)
    real(real64), allocatable   :: xVec(:)
    integer(int32), allocatable :: MatingPlan(:, :)
    real(real64), allocatable   :: GenomeEdit(:)

    ! Column headers and formats for logging
    character(len=CHARLENGTH)      :: FmtLogStdoutHead
    character(len=CHARLENGTH)      :: FmtLogStdout
    character(len=15), allocatable :: ColnameLogStdout(:)
    character(len=CHARLENGTH)      :: FmtLogUnitHead
    character(len=CHARLENGTH)      :: FmtLogUnit
    character(len=22), allocatable :: ColnameLogUnit(:)
    character(len=CHARLENGTH)      :: FmtLogPopUnit
    character(len=22), allocatable :: ColnameLogPopUnit(:)
    character(len=CHARLENGTH)      :: FmtContributionHead     ! = "(6a15)"
    character(len=CHARLENGTH)      :: FmtContributionHeadEdit ! = "(8a15)"
    character(len=CHARLENGTH)      :: FmtContribution         ! = "(a??, 4x, i11, 3(4x, f11.5), 4x, i11)"
    character(len=CHARLENGTH)      :: FmtContributionEdit     != "(a??, 4x, i11, 3(4x, f11.5), 2(4x, i11), 4x, f11.5)"
    character(len=CHARLENGTH)      :: FmtMatingHead           != "(a15, 2a??)"
    character(len=CHARLENGTH)      :: FmtMating               != "(i15, 1x, 2a??)"

    contains
      procedure :: Initialise   => InitialiseAlphaMateSol
      procedure :: Assign       => AssignAlphaMateSol
      procedure :: UpdateMean   => UpdateMeanAlphaMateSol
      procedure :: Evaluate     => FixSolEtcMateAndEvaluateAlphaMateSol
      procedure :: Write        => WriteAlphaMateSol
      procedure :: WriteContributions
      procedure :: WriteMatingPlan
      procedure :: SetupColNamesAndFormats
      procedure :: LogHead      => LogHeadAlphaMateSol
      procedure :: Log          => LogAlphaMateSol
      procedure :: LogPopHead   => LogPopHeadAlphaMateSol
      procedure :: LogPop       => LogPopAlphaMateSol
  end type

  contains

    !###########################################################################

    ! TODO: make this clearer with objective being dG - l * dF etc.

    ! With two individuals there are four genome combinations, hence four coefficients
    ! that measure similarity between the two individuals. The COEFFICIENT OF COANCESTRY
    ! is an average of these four coefficients, f_i,j = 1/4 (f_i1,j1 + f_i1,j2 + f_i2,j1 + f_i2,j2).
    ! When the two individuals in comparison is just one individual, there are only
    ! two genomes to compare, hence f_i,i = 1/4 (f_i1,i1 + f_i1,i2 + f_i2,i1 + f_i2,i2) =
    ! 1/4 (1 + f_i1,i2 + f_i2,i1 + 1) = 1/4 (2 + 2f_i1,i2) = 1/2 (1 + f_i1,i2). The
    ! coefficient f_i1,i2 is the COEFFICIENT OF INBREEDING.

    !       Also Coancestry is expected future inbreeding = average coancestry between
    !       the parents (selection/mate candidates we have at hand). When working
    !       with dF we then need to evaluate average coancestry in the new generation.

    ! This note is just to make clear what the objective function and its components are.
    !
    ! AlphaMate works with the objective function:
    !
    ! Objective = x'a - l * x'Cx
    !
    !   x is a vector of contributions of individuals to the next generation (sum(x)=1)
    !   a is a vector of selection criterion values
    !   l is a penalty for the loss of genetic diversity
    !   C is the coancestry (kinship) matrix where
    !      the diagonal values are 1/2 * (1 + F_i,i), with F_i being inbreeding
    !        coefficient of individual i (note: F_i,i = C_f(i),m(i))
    !      the off-diagonal values are coancestry coefficients between individuals (C_i,j)
    !      (note: the numerator relationship matrix A = 2 * C)
    !
    !   x'a is average selection criterion  of the contributing gametes to the next generation,
    !     which is a proxy for the expected selection criterion value of the next generation (SelCriterion)
    !     (note: SelCriterion = x'a assumes random mating; @todo SelCriterion under
    !      designed mating will be added in future )
    !   x'Cx is average coancestry among the contributing gametes to the next generation,
    !     which is the expected population inbreeding in the next generation (CoancestryRandom)
    !     (note: CoancestryRandom = x'Cx assumes random mating (selfing included) with contributions x;
    !      CoancestryRandom under random mating without selfing is hard to calculate exactly;
    !      @todo CoancestryRandom under designed mating will be added in future)
    !
    ! We want to maximize selection criterion in the next generation (x'a), but control
    ! for loss of genetic diversity (x'Cx). For operational reasons optimisation
    ! does not work directly with x'a and x'Cx, but with x's and dF
    !
    !       s   is a standardized selection criterion, s = (a - mean(a)) / sd(a)
    !       x's is the genetic selection intensity (SelIntensity)
    !       dF  is the rate of population inbreeding (CoancestryRate)
    !
    ! Use of x's and dF essentially removes the need to specify the hard to define
    ! penalty value (l) in favour of easier to define rate of inbreeding (dF), which
    ! is calculated based on the following relations.
    !
    ! The total inbreeding that we see today is:
    !
    ! F_total = F_new + (1 - F_new) * F_old,
    !
    ! where F_new is new/recent  inbreeding built up after  some time point
    !       F_old is old/ancient inbreeding built up before some time point
    !
    ! Equivalent formulas with different notation are:
    !
    ! F_T = F_IS + (1 - F_IS) * F_ST
    !
    !   F_T is (total) inbreeding relative to old/ancient base time point
    !   F_IS is (individual-to-subtotal=new/recent) inbreeding relative to recent base time point
    !   F-ST is (subtotal-to-total=old/ancient)     inbreeding at recent base time point relative to old/ancient base time point
    !
    ! F_t+1 = dF + (1 - dF) * F_t
    !       = dF * (1 - F_t) + F_t
    !
    !   F_t+1 is (total) inbreeding in the next generation (t + 1)
    !   dF is new/recent inbreeding built up between the generations t and t + 1 (=rate of inbreeding)
    !   F-t is old/ancient inbreeding built up to the generation t
    !
    ! Then:
    !
    ! F_total   = F_T  = F_t+1
    ! F_new     = F_IS = dF
    ! F_old     = F_ST = F_t
    !
    ! Finally:
    !
    ! dF = (F_t+1 - F_t) / (1 - F_t)
    !
    ! In AlphaMate F_t is calculated as x'Cx with x set to 1/n (n number of ???), so we obtain average
    ! of the values in the matrix C, the average coancestry in current population, which is the expected
    !
    ! With available F_t and dF the F_t+1 is then calculated as shown above.

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
      write(STDOUT, "(a)") "         Software for optimizing contributions to the next generation         "
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
      This%SelCriterionGiven = .false.
      This%GenderGiven = .false.
      This%GenericIndCritGiven = .false.
      This%GenericMatCritGiven = .false.
      This%SeedFileGiven = .false.
      This%SeedGiven = .false.
      This%NrmInsteadOfCoancestry = .false.

      ! Search mode specifications

      This%ModeMinCoancestry = .false.
      This%ModeMinInbreeding = .false.
      This%ModeMaxCriterion = .false.
      This%ModeOpt = .false.
      This%ModeRan = .false.
      This%EvaluateFrontier = .false.

      ! Biological specifications

      This%nInd = 0
      This%nMat = 0
      This%nPar = 0
      This%nPar1 = 0
      This%nPar2 = 0

      This%nTargets = 0
      This%TargetDegreeGiven = .false.
      ! This%TargetDegree ! allocatable so skip here

      This%TargetSelCriterionGiven = .false.
      ! This%TargetSelCriterion ! allocatable so skip here
      This%TargetSelIntensityGiven = .false.
      ! This%TargetSelIntensity ! allocatable so skip here
      This%TargetMaxCriterionPctGiven = .false.
      ! This%TargetMaxCriterionPct ! allocatable so skip here

      This%TargetCoancestryGiven = .false.
      ! This%TargetCoancestry ! allocatable so skip here
      This%TargetCoancestryRateGiven = .false.
      ! This%TargetCoancestryRate ! allocatable so skip here
      This%CoancestryWeight = -1 ! -1000.0d0
      This%CoancestryWeightBelow = .false.
      This%TargetMinCoancestryPctGiven = .false.
      ! This%TargetMinCoancestryPct ! allocatable so skip here

      This%TargetInbreedingGiven = .false.
      This%TargetInbreeding = 1.0d0 ! set it high so that it does not have any effect
      This%TargetInbreedingRateGiven = .false.
      This%TargetInbreedingRate = 1.0d0 ! set it high so that it does not have any effect
      This%InbreedingWeight = 0.0d0 ! set it to zero so that it does not have any effect
      This%InbreedingWeightBelow = .false.
      This%TargetMinInbreedingPctGiven = .false.
      This%TargetMinInbreedingPct = 100.0d0 ! not used anyhow as TargetMinInbreedingPctGiven = .false.

      This%SelfingAllowed = .false.
      This%SelfingWeight = -1 ! -1000.0d0
      ! This%GenericIndCritWeight ! allocatable so skip here
      ! This%GenericMatCritWeight ! allocatable so skip here

      This%EqualizePar  = .false.
      This%EqualizePar1 = .false.
      This%EqualizePar2 = .false.

      This%LimitPar  = .false.
      This%LimitPar1 = .false.
      This%LimitPar2 = .false.
      This%LimitParMin  = 1.0d0
      This%LimitPar1Min = 1.0d0
      This%LimitPar2Min = 1.0d0
      This%LimitParMax  = huge(This%LimitParMax)  - 1.0d0
      This%LimitPar1Max = huge(This%LimitPar1Max) - 1.0d0
      This%LimitPar2Max = huge(This%LimitPar2Max) - 1.0d0
      This%LimitParMinWeight  = -1000.0d0
      This%LimitPar1MinWeight = -1000.0d0
      This%LimitPar2MinWeight = -1000.0d0

      This%PAGEPar  = .false.
      This%PAGEPar1 = .false.
      This%PAGEPar2 = .false.
      This%PAGEParMax  = 0
      This%PAGEPar1Max = 0
      This%PAGEPar2Max = 0
      This%PAGEParCost  = 0.0d0
      This%PAGEPar1Cost = 0.0d0
      This%PAGEPar2Cost = 0.0d0

      ! Search algorithm specifications

      This%EvolAlgNSol = 100
      This%EvolAlgNIter = 100000
      This%EvolAlgNIterStop = 1000
      This%EvolAlgNIterPrint = 100
      This%EvolAlgStopTol = 0.0001d0
      This%EvolAlgLogPop = .false.
      This%EvolAlg = "DE"

      This%DiffEvolNIterBurnIn = 1000
      This%DiffEvolParamCrBurnIn = 0.4d0
      This%DiffEvolParamCr = 0.2d0
      This%DiffEvolParamFBase = 0.1d0
      This%DiffEvolParamFHigh1 = 1.0d0
      This%DiffEvolParamFHigh2 = 4.0d0

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

      write(Unit, *) "OutputBasename: ",     trim(This%OutputBasename)
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

      write(Unit, *) "RelMtxGiven: ",            This%RelMtxGiven
      write(Unit, *) "SelCriterionGiven: ",      This%SelCriterionGiven
      write(Unit, *) "GenderGiven: ",            This%GenderGiven
      write(Unit, *) "GenericIndCritGiven: ",    This%GenericIndCritGiven
      write(Unit, *) "GenericMatCritGiven: ",    This%GenericMatCritGiven
      write(Unit, *) "SeedFileGiven: ",          This%SeedFileGiven
      write(Unit, *) "SeedGiven: ",              This%SeedGiven
      write(Unit, *) "NrmInsteadOfCoancestry: ", This%NrmInsteadOfCoancestry

      ! Search mode specifications

      write(Unit, *) "ModeMinCoancestry: ",  This%ModeMinCoancestry
      write(Unit, *) "ModeMinInbreeding: ",  This%ModeMinInbreeding
      write(Unit, *) "ModeMaxCriterion: ",  This%ModeMaxCriterion
      write(Unit, *) "ModeRan: ",  This%ModeRan
      write(Unit, *) "ModeOpt: ",  This%ModeOpt
      write(Unit, *) "nTargets: ", This%nTargets

      ! Biological specifications

      write(Unit, *) "nInd: ",  This%nInd
      write(Unit, *) "nMat: ",  This%nMat
      write(Unit, *) "nPar: ",  This%nPar
      write(Unit, *) "nPar1: ", This%nPar1
      write(Unit, *) "nPar2: ", This%nPar2

      write(Unit, *) "nTargets: ",                        This%nTargets

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
      write(Unit, *) "TargetSelIntensityGiven: ",         This%TargetSelIntensityGiven
      if (allocated(This%TargetSelIntensity)) then
        write(Unit, *) "TargetSelIntensity: ",            This%TargetSelIntensity
      else
        write(Unit, *) "TargetSelIntensity: not allocated"
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

      write(Unit, *) "SelfingAllowed: ",                  This%SelfingAllowed
      write(Unit, *) "SelfingWeight: ",                   This%SelfingWeight
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

      write(Unit, *) "EqualizePar:  ", This%EqualizePar
      write(Unit, *) "EqualizePar1: ", This%EqualizePar1
      write(Unit, *) "EqualizePar2: ", This%EqualizePar2

      write(Unit, *) "LimitPar:  ",          This%LimitPar
      write(Unit, *) "LimitPar1: ",          This%LimitPar1
      write(Unit, *) "LimitPar2: ",          This%LimitPar2
      write(Unit, *) "LimitParMin: ",        This%LimitParMin
      write(Unit, *) "LimitPar1Min: ",       This%LimitPar1Min
      write(Unit, *) "LimitPar2Min: ",       This%LimitPar2Min
      write(Unit, *) "LimitParMax: ",        This%LimitParMax
      write(Unit, *) "LimitPar1Max: ",       This%LimitPar1Max
      write(Unit, *) "LimitPar2Max: ",       This%LimitPar2Max
      write(Unit, *) "LimitParMinWeight: ",  This%LimitParMinWeight
      write(Unit, *) "LimitPar1MinWeight: ", This%LimitPar1MinWeight
      write(Unit, *) "LimitPar2MinWeight: ", This%LimitPar2MinWeight

      write(Unit, *) "PAGEPar: ",      This%PAGEPar
      write(Unit, *) "PAGEPar1: ",     This%PAGEPar1
      write(Unit, *) "PAGEPar2: ",     This%PAGEPar2
      write(Unit, *) "PAGEParMax: ",   This%PAGEParMax
      write(Unit, *) "PAGEPar1Max: ",  This%PAGEPar1Max
      write(Unit, *) "PAGEPar2Max: ",  This%PAGEPar2Max
      write(Unit, *) "PAGEParCost: ",  This%PAGEParCost
      write(Unit, *) "PAGEPar1Cost: ", This%PAGEPar1Cost
      write(Unit, *) "PAGEPar2Cost: ", This%PAGEPar2Cost

      ! Algorithm specifications

      write(Unit, *) "EvolAlgNSol: ",       This%EvolAlgNSol
      write(Unit, *) "EvolAlgNIter: ",      This%EvolAlgNIter
      write(Unit, *) "EvolAlgNIterStop: ",  This%EvolAlgNIterStop
      write(Unit, *) "EvolAlgNIterPrint: ", This%EvolAlgNIterPrint
      write(Unit, *) "EvolAlgStopTol: ",    This%EvolAlgStopTol
      write(Unit, *) "EvolAlgLogPop: ",     This%EvolAlgLogPop
      write(Unit, *) "EvolAlg: ",           This%EvolAlg

      write(Unit, *) "DiffEvolNIterBurnIn: ",   This%DiffEvolNIterBurnIn
      write(Unit, *) "DiffEvolParamCrBurnIn: ", This%DiffEvolParamCrBurnIn
      write(Unit, *) "DiffEvolParamCr: ",       This%DiffEvolParamCr
      write(Unit, *) "DiffEvolParamFBase: ",    This%DiffEvolParamFBase
      write(Unit, *) "DiffEvolParamFHigh1: ",   This%DiffEvolParamFHigh1
      write(Unit, *) "DiffEvolParamFHigh2: ",   This%DiffEvolParamFHigh2

      write(Unit, *) "RanAlgStricter: ",       This%RanAlgStricter

      ! Data&Spec derived quantities

      write(Unit, *) "ModeSpec: "
      call This%ModeSpec%Write(Unit)
      write(Unit, *) "ModeMinCoancestrySpec: "
      call This%ModeMinCoancestrySpec%Write(Unit)
      write(Unit, *) "ModeMinInbreedingSpec: "
      call This%ModeMinInbreedingSpec%Write(Unit)
      write(Unit, *) "ModeMaxCriterionSpec: "
      call This%ModeMaxCriterionSpec%Write(Unit)
      !@todo do we need ModeRanSpec?
      write(Unit, *) "ModeRanSpec: "
      call This%ModeRanSpec%Write(Unit)
      write(Unit, *) "ModeOptSpec: "
      call This%ModeOptSpec%Write(Unit)

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
        ! @todo why (len_trim(Line) .eq. 0)? if we use (len_trim(Line) .eq. 0) above
        if (First(1:1) .eq. "=" .or. len_trim(Line) .eq. 0) then
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
                write(STDERR, "(a)") " ERROR: Must specify a string for OutputBasename, i.e., OutputBasename, AnalysisX"
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
                write(STDERR, "(a)") " ERROR: Must specify a file for CoancestryMatrixFile, i.e., CoancestryMatrixFile, CoaMtx.txt"
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
                  write(STDERR, "(a)") " ERROR: Must specify a file for NrmMatrixFile, i.e., NrmMatrixFile, NrmMtx.txt"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a file for NrmMatrixFile, i.e., NrmMatrixFile, NrmMtx.txt"
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
                write(STDERR, "(a)") " ERROR: Must specify a file for SelCriterionFile, i.e., SelCriterionFile, SelCrit.txt"
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
                write(STDERR, "(a)") " ERROR: Must specify a file for GenderFile, i.e., GenderFile, Gender.txt"
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
                write(STDERR, "(a)") " ERROR: Must specify a file for SeedFile, i.e., SeedFile, Seed.txt"
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
                write(STDERR, "(a)") " ERROR: Must specify a number for Seed, i.e., Seed, 19791123"
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
                write(STDERR, "(a)") " ERROR: Must specify a file for GenericIndividualCriterionFile, i.e., GenericIndividualCriterionFile, IndividualCriterion.txt"
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
                  write(STDERR, "(a)") " ERROR: Must specify a number for GenericIndividualCriterionColumns, i.e., GenericIndividualCriterionColumns, 2"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              end if

            case ("genericindividualcriterionweight")
              if (This%GenericIndCritGiven) then
                if (allocated(Second)) then
                  nGenericIndCrit = nGenericIndCrit + 1
                  This%GenericIndCritWeight(nGenericIndCrit) = Char2Double(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Generic individual criterion - weight ("//trim(Int2Char(nGenericIndCrit))//"): "//trim(Real2Char(This%GenericIndCritWeight(nGenericIndCrit), fmt=FMTREAL2CHAR))
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify a value for GenericIndividualCriterionWeight, i.e., GenericIndividualCriterionWeight, 2"
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
                    write(STDOUT, "(a)") " Generic mating criterion file: "//trim(This%GenericMatCritFile)
                  end if
                  This%GenericMatCritGiven = .true.
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a file for GenericMatingCriterionFile, i.e., GenericMatingCriterionFile, MatingCriterion.txt"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("genericmatingcriterioncolumns")
              if (This%GenericMatCritGiven) then
                if (allocated(Second)) then
                  This%nGenericMatCrit = Char2Int(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Generic mating criterion - number of columns: "//trim(Int2Char(This%nGenericMatCrit))
                  end if
                  allocate(This%GenericMatCritWeight(This%nGenericMatCrit))
                  nGenericMatCrit = 0
                else
                  write(STDERR, "(a)") " ERROR: Must specify a number for GenericMatingCriterionColumns, i.e., GenericMatingCriterionColumns, 2"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              end if

            case ("genericmatingcriterionweight")
              if (This%GenericMatCritGiven) then
                if (allocated(Second)) then
                  nGenericMatCrit = nGenericMatCrit + 1
                  This%GenericMatCritWeight(nGenericMatCrit) = Char2Double(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Generic mating criterion - weight ("//trim(Int2Char(nGenericMatCrit))//"): "//trim(Real2Char(This%GenericMatCritWeight(nGenericMatCrit), fmt=FMTREAL2CHAR))
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify a value for GenericMatingCriterionWeight, i.e., GenericMatingCriterionWeight, 2"
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
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for ModeMinCoancestry, i.e., ModeMinCoancestry, Yes"
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
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for ModeMinInbreeding, i.e., ModeMinInbreeding, Yes"
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
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for ModeMaxCriterion, i.e., ModeMaxCriterion, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("moderan")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%ModeRan = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " ModeRan"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for ModeRan, i.e., ModeRan, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("modeopt")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%ModeOpt = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " ModeOpt"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for ModeOpt, i.e., ModeOpt, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("evaluatefrontier")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%EvaluateFrontier = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " EvaluateFrontier"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for EvaluateFrontier, i.e., EvaluateFrontier, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            ! Biological specifications
            case ("numberofmatings")
              if (allocated(Second)) then
                This%nMat = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Number of matings: "//trim(Int2Char(This%nMat))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfMatings, i.e., NumberOfMatings, 10"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("numberofparents")
              if (allocated(Second)) then
                This%nPar = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Number of parents: "//trim(Int2Char(This%nPar))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfParents, i.e., NumberOfParents, 20"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("numberofmaleparents")
              if (allocated(Second)) then
                This%nPar1 = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Number of male parents: "//trim(Int2Char(This%nPar1))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfMaleParents, i.e., NumberOfMaleParents, 10"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("numberoffemaleparents")
              if (allocated(Second)) then
                This%nPar2 = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Number of female parents: "//trim(Int2Char(This%nPar2))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfFemaleParents, i.e., NumberOfFemaleParents, 10"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("targetdegree")
              if (allocated(Second)) then
                This%ModeOpt = .true.
                This%TargetDegreeGiven = .true.
                This%nTargets = This%nTargets + 1
                block
                  integer(int32) :: n
                  real(real64), allocatable :: Tmp(:)
                  if (allocated(This%TargetDegree)) then
                    n = size(This%TargetDegree)
                    allocate(Tmp(n))
                    Tmp = This%TargetDegree
                    deallocate(This%TargetDegree)
                    allocate(This%TargetDegree(n + 1))
                    This%TargetDegree(1:n) = Tmp
                    n = n + 1
                  else
                    n = 1
                    allocate(This%TargetDegree(n))
                  end if
                  This%TargetDegree(n) = Char2Double(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Targeted degree: "//trim(Real2Char(This%TargetDegree(n), fmt=FMTREAL2CHAR))
                  end if
                  if ((This%TargetDegree(n) .lt. 0.0d0) .or. (This%TargetDegree(n) .gt. 90.0d0)) then
                    write(STDERR, "(a)") "ERROR: TargetDegree must be between 0 and 90!"
                    write(STDERR, "(a)") " "
                    stop 1
                  end if
                end block
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetDegree, i.e., TargetDegree, 45"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("targetselcriterion")
              if (allocated(Second)) then
                This%ModeOpt = .true.
                This%TargetSelCriterionGiven = .true.
                This%nTargets = This%nTargets + 1
                block
                  integer(int32) :: n
                  real(real64), allocatable :: Tmp(:)
                  if (allocated(This%TargetSelCriterion)) then
                    n = size(This%TargetSelCriterion)
                    allocate(Tmp(n))
                    Tmp = This%TargetSelCriterion
                    deallocate(This%TargetSelCriterion)
                    allocate(This%TargetSelCriterion(n + 1))
                    This%TargetSelCriterion(1:n) = Tmp
                    n = n + 1
                  else
                    n = 1
                    allocate(This%TargetSelCriterion(n))
                  end if
                  This%TargetSelCriterion(n) = Char2Double(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Targeted selection criterion: "//trim(Real2Char(This%TargetSelCriterion(n), fmt=FMTREAL2CHAR))
                  end if
                end block
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetSelCriterion, i.e., TargetSelCriterion, 100"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("targetselintensity")
              if (allocated(Second)) then
                This%ModeOpt = .true.
                This%TargetSelIntensityGiven = .true.
                This%nTargets = This%nTargets + 1
                block
                  integer(int32) :: n
                  real(real64), allocatable :: Tmp(:)
                  if (allocated(This%TargetSelIntensity)) then
                    n = size(This%TargetSelIntensity)
                    allocate(Tmp(n))
                    Tmp = This%TargetSelIntensity
                    deallocate(This%TargetSelIntensity)
                    allocate(This%TargetSelIntensity(n + 1))
                    This%TargetSelIntensity(1:n) = Tmp
                    n = n + 1
                  else
                    n = 1
                    allocate(This%TargetSelIntensity(n))
                  end if
                  This%TargetSelIntensity(n) = Char2Double(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Targeted selection intensity: "//trim(Real2Char(This%TargetSelIntensity(n), fmt=FMTREAL2CHAR))
                  end if
                  if ((This%TargetSelIntensity(n) .lt. 0.0d0) .or. (This%TargetSelIntensity(n) .gt. 5.0d0)) then
                    write(STDERR, "(a)") "ERROR: TargetSelIntensity must be above 0 and (probably) bellow 5!"
                    write(STDERR, "(a)") " "
                    stop 1
                  end if
                end block
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetSelIntensity, i.e., TargetSelIntensity, 2"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("targetmaxcriterionpct")
              if (allocated(Second)) then
                This%ModeOpt = .true.
                This%TargetMaxCriterionPctGiven = .true.
                This%nTargets = This%nTargets + 1
                block
                  integer(int32) :: n
                  real(real64), allocatable :: Tmp(:)
                  if (allocated(This%TargetMaxCriterionPct)) then
                    n = size(This%TargetMaxCriterionPct)
                    allocate(Tmp(n))
                    Tmp = This%TargetMaxCriterionPct
                    deallocate(This%TargetMaxCriterionPct)
                    allocate(This%TargetMaxCriterionPct(n + 1))
                    This%TargetMaxCriterionPct(1:n) = Tmp
                    n = n + 1
                  else
                    n = 1
                    allocate(This%TargetMaxCriterionPct(n))
                  end if
                  This%TargetMaxCriterionPct(n) = Char2Double(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Targeted percentage of maximum criterion: "//trim(Real2Char(This%TargetMaxCriterionPct(n), fmt=FMTREAL2CHAR))
                  end if
                  if ((This%TargetMaxCriterionPct(n) .lt. 0.0d0) .or. (This%TargetMaxCriterionPct(n) .gt. 100.0d0)) then
                    write(STDERR, "(a)") "ERROR: TargetMaxCriterionPct must be between 0 and 100!"
                    write(STDERR, "(a)") " "
                    stop 1
                  end if
                end block
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetMaxCriterionPct, i.e., TargetMaxCriterionPct, 90"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("targetcoancestry")
              if (allocated(Second)) then
                This%ModeOpt = .true.
                This%TargetCoancestryGiven = .true.
                This%nTargets = This%nTargets + 1
                block
                  integer(int32) :: n
                  real(real64), allocatable :: Tmp(:)
                  if (allocated(This%TargetCoancestry)) then
                    n = size(This%TargetCoancestry)
                    allocate(Tmp(n))
                    Tmp = This%TargetCoancestry
                    deallocate(This%TargetCoancestry)
                    allocate(This%TargetCoancestry(n + 1))
                    This%TargetCoancestry(1:n) = Tmp
                    n = n + 1
                  else
                    n = 1
                    allocate(This%TargetCoancestry(n))
                  end if
                  This%TargetCoancestry(n) = Char2Double(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Targeted coancestry: "//trim(Real2Char(This%TargetCoancestry(n), fmt=FMTREAL2CHAR))
                  end if
                  if ((This%TargetCoancestry(n) .lt. -1.0d0) .or. (This%TargetCoancestry(n) .gt. 1.0d0)) then
                    write(STDERR, "(a)") "ERROR: TargetCoancestry must be between -1 and +1!"
                    write(STDERR, "(a)") " "
                    stop 1
                  end if
                end block
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetCoancestry, i.e., TargetCoancestry, 0.31"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("targetcoancestryrate")
              if (allocated(Second)) then
                This%ModeOpt = .true.
                This%TargetCoancestryRateGiven = .true.
                This%nTargets = This%nTargets + 1
                block
                  integer(int32) :: n
                  real(real64), allocatable :: Tmp(:)
                  if (allocated(This%TargetCoancestryRate)) then
                    n = size(This%TargetCoancestryRate)
                    allocate(Tmp(n))
                    Tmp = This%TargetCoancestryRate
                    deallocate(This%TargetCoancestryRate)
                    allocate(This%TargetCoancestryRate(n + 1))
                    This%TargetCoancestryRate(1:n) = Tmp
                    n = n + 1
                  else
                    n = 1
                    allocate(This%TargetCoancestryRate(n))
                  end if
                  This%TargetCoancestryRate(n) = Char2Double(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Targeted rate of coancestry: "//trim(Real2Char(This%TargetCoancestryRate(n), fmt=FMTREAL2CHAR))
                  end if
                  if ((This%TargetCoancestryRate(n) .lt. -1.0d0) .or. (This%TargetCoancestryRate(n) .gt. 1.0d0)) then
                    write(STDERR, "(a)") "ERROR: TargetCoancestryRate must be between -1 and +1!"
                    write(STDERR, "(a)") " "
                    stop 1
                  end if
                end block
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetCoancestryRate, i.e., TargetCoancestryRate, 0.01"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("targetmincoancestrypct")
              if (allocated(Second)) then
                This%ModeOpt = .true.
                This%TargetMinCoancestryPctGiven = .true.
                This%nTargets = This%nTargets + 1
                block
                  integer(int32) :: n
                  real(real64), allocatable :: Tmp(:)
                  if (allocated(This%TargetMinCoancestryPct)) then
                    n = size(This%TargetMinCoancestryPct)
                    allocate(Tmp(n))
                    Tmp = This%TargetMinCoancestryPct
                    deallocate(This%TargetMinCoancestryPct)
                    allocate(This%TargetMinCoancestryPct(n + 1))
                    This%TargetMinCoancestryPct(1:n) = Tmp
                    n = n + 1
                  else
                    n = 1
                    allocate(This%TargetMinCoancestryPct(n))
                  end if
                  This%TargetMinCoancestryPct(n) = Char2Double(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Targeted percentage of minimum coancestry: "//trim(Real2Char(This%TargetMinCoancestryPct(n), fmt=FMTREAL2CHAR))
                  end if
                  if ((This%TargetMinCoancestryPct(n) .lt. 0.0d0) .or. (This%TargetMinCoancestryPct(n) .gt. 100.0d0)) then
                    write(STDERR, "(a)") "ERROR: TargetMinCoancestryPct must be between 0 and 100!"
                    write(STDERR, "(a)") " "
                    stop 1
                  end if
                end block
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetMinCoancestryPct, i.e., TargetMinCoancestryPct, 90"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("coancestryweight")
              if (allocated(Second)) then
                This%CoancestryWeight = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Targeted coancestry - weight: "//trim(Real2Char(This%CoancestryWeight, fmt=FMTREAL2CHAR))
                end if
                if (This%CoancestryWeight .gt. 0.0d0) then
                  write(STDOUT, "(a)") " NOTE: Positive weight for the targeted coancestry, i.e., encourage higher coancestry. Was this intended?"
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for CoancestryWeight, i.e., CoancestryWeight, -1000"
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
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for CoancestryWeightBelow, i.e., CoancestryWeightBelow, No"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("targetinbreeding")
              if (allocated(Second)) then
                This%TargetInbreedingGiven = .true.
                ! This%nTargets = This%nTargets + 1
                ! block
                !   integer(int32) :: n
                !   real(real64), allocatable :: Tmp(:)
                !   if (allocated(This%TargetInbreeding)) then
                !     n = size(This%TargetInbreeding)
                !     allocate(Tmp(n))
                !     Tmp = This%TargetInbreeding
                !     deallocate(This%TargetInbreeding)
                !     allocate(This%TargetInbreeding(n + 1))
                !     This%TargetInbreeding(1:n) = Tmp
                !     n = n + 1
                !   else
                !     n = 1
                !     allocate(This%TargetInbreeding(n))
                !   end if
                !   This%TargetInbreeding(n) = Char2Double(trim(adjustl(Second(1))))
                !   if (LogStdoutInternal) then
                !     write(STDOUT, "(a)") " Targeted inbreding: "//trim(Real2Char(This%TargetInbreeding(n), fmt=FMTREAL2CHAR))
                !   end if
                !   if ((This%TargetInbreeding(n) .lt. -1.0d0) .or. (This%TargetInbreeding(n) .gt. 1.0d0)) then
                !     write(STDERR, "(a)") "ERROR: TargetInbreeding must be between -1 and +1!"
                !     write(STDERR, "(a)") " "
                !     stop 1
                !   end if
                ! end block
                This%TargetInbreeding = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Targeted inbreding: "//trim(Real2Char(This%TargetInbreeding, fmt=FMTREAL2CHAR))
                end if
                if ((This%TargetInbreeding .lt. -1.0d0) .or. (This%TargetInbreeding .gt. 1.0d0)) then
                  write(STDERR, "(a)") "ERROR: TargetInbreeding must be between -1 and +1!"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetInbreeding, i.e., TargetInbreeding, 0.31"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("targetinbreedingrate")
              if (allocated(Second)) then
                This%ModeOpt = .true.
                This%TargetInbreedingRateGiven = .true.
                ! This%nTargets = This%nTargets + 1
                ! block
                !   integer(int32) :: n
                !   real(real64), allocatable :: Tmp(:)
                !   if (allocated(This%TargetInbreedingRate)) then
                !     n = size(This%TargetInbreedingRate)
                !     allocate(Tmp(n))
                !     Tmp = This%TargetInbreedingRate
                !     deallocate(This%TargetInbreedingRate)
                !     allocate(This%TargetInbreedingRate(n + 1))
                !     This%TargetInbreedingRate(1:n) = Tmp
                !     n = n + 1
                !   else
                !     n = 1
                !     allocate(This%TargetInbreedingRate(n))
                !   end if
                !   This%TargetInbreedingRate(n) = Char2Double(trim(adjustl(Second(1))))
                !   if (LogStdoutInternal) then
                !     write(STDOUT, "(a)") " Targeted rate of inbreeding: "//trim(Real2Char(This%TargetInbreedingRate(n), fmt=FMTREAL2CHAR))
                !   end if
                !   if ((This%TargetInbreedingRate(n) .lt. -1.0d0) .or. (This%TargetInbreedingRate(n) .gt. 1.0d0)) then
                !     write(STDERR, "(a)") "ERROR: TargetInbreedingRate must be between -1 and +1!"
                !     write(STDERR, "(a)") " "
                !     stop 1
                !   end if
                ! end block
                This%TargetInbreedingRate = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Targeted rate of inbreeding: "//trim(Real2Char(This%TargetInbreedingRate, fmt=FMTREAL2CHAR))
                end if
                if ((This%TargetInbreedingRate .lt. -1.0d0) .or. (This%TargetInbreedingRate .gt. 1.0d0)) then
                  write(STDERR, "(a)") "ERROR: TargetInbreedingRate must be between -1 and +1!"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetInbreedingRate, i.e., TargetInbreedingRate, 0.01"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("targetmininbreedingpct")
              if (allocated(Second)) then
                This%ModeOpt = .true.
                This%TargetMinInbreedingPctGiven = .true.
                ! This%nTargets = This%nTargets + 1
                ! block
                !   integer(int32) :: n
                !   real(real64), allocatable :: Tmp(:)
                !   if (allocated(This%TargetMinInbreedingPct)) then
                !     n = size(This%TargetMinInbreedingPct)
                !     allocate(Tmp(n))
                !     Tmp = This%TargetMinInbreedingPct
                !     deallocate(This%TargetMinInbreedingPct)
                !     allocate(This%TargetMinInbreedingPct(n + 1))
                !     This%TargetMinInbreedingPct(1:n) = Tmp
                !     n = n + 1
                !   else
                !     n = 1
                !     allocate(This%TargetMinInbreedingPct(n))
                !   end if
                !   This%TargetMinInbreedingPct(n) = Char2Double(trim(adjustl(Second(1))))
                !   if (LogStdoutInternal) then
                !     write(STDOUT, "(a)") " Targeted percentage of minimum inbreeding: "//trim(Real2Char(This%TargetMinInbreedingPct(n), fmt=FMTREAL2CHAR))
                !   end if
                !   if ((This%TargetMinInbreedingPct(n) .lt. 0.0d0) .or. (This%TargetMinInbreedingPct(n) .gt. 100.0d0)) then
                !     write(STDERR, "(a)") "ERROR: TargetMinInbreedingPct must be between 0 and 100!"
                !     write(STDERR, "(a)") " "
                !     stop 1
                !   end if
                ! end block
                This%TargetMinInbreedingPct = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Targeted percentage of minimum inbreeding: "//trim(Real2Char(This%TargetMinInbreedingPct, fmt=FMTREAL2CHAR))
                end if
                if ((This%TargetMinInbreedingPct .lt. 0.0d0) .or. (This%TargetMinInbreedingPct .gt. 100.0d0)) then
                  write(STDERR, "(a)") "ERROR: TargetMinInbreedingPct must be between 0 and 100!"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetMinInbreedingPct, i.e., TargetMinInbreedingPct, 90"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("inbreedingweight")
              if (allocated(Second)) then
                This%InbreedingWeight = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Targeted inbreeding - weight: "//trim(Real2Char(This%InbreedingWeight, fmt=FMTREAL2CHAR))
                end if
                if (This%InbreedingWeight .gt. 0.0d0) then
                  write(STDOUT, "(a)") " NOTE: Positive weight for targeted inbreeding, i.e., encourage higher inbreeding. Was this intended?"
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for InbreedingWeight, i.e., InbreedingWeight, -1000"
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
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for InbreedingWeightBelow, i.e., InbreedingWeightBelow, No"
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
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for EqualizeContributions, i.e., EqualizeContributions, Yes"
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
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for EqualizeMaleContributions, i.e., EqualizeMaleContributions, Yes"
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
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for EqualizeFemaleContributions, i.e., EqualizeFemaleContributions, Yes"
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
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for LimitContributions, i.e., LimitContributions, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("limitcontributionsmin")
              if (This%LimitPar) then
                if (allocated(Second)) then
                  This%LimitParMin = Char2Double(trim(adjustl(Second(1)))) ! real because of continous solution representation
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Limit contributions - minimum: "//trim(Int2Char(nint(This%LimitParMin))) ! nint because of continous solution representation
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify a number for LimitContributionsMin, i.e., LimitContributionsMin, 1"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              end if

            case ("limitcontributionsmax")
              if (This%LimitPar) then
                if (allocated(Second)) then
                  This%LimitParMax = Char2Double(trim(adjustl(Second(1)))) ! Real because of continous solution representation
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Limit contributions - maximum: "//trim(Int2Char(nint(This%LimitParMax))) ! nint because of continous solution representation
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify a number for LimitContributionsMax, i.e., LimitContributionsMax, 10"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              end if

            case ("limitcontributionsminweight")
              if (This%LimitPar) then
                if (allocated(Second)) then
                  This%LimitParMinWeight = Char2Double(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Limit contributions - weight for contributions bellow minimum: "//trim(Real2Char(This%LimitParMinWeight, fmt=FMTREAL2CHAR))
                  end if
                  if (This%LimitParMinWeight .gt. 0.0d0) then
                    write(STDOUT, "(a)") " NOTE: Positive weight for limit on minimum contributions, i.e., encourage smaller contributions than defined minimum. Was this intended?"
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify a value for LimitContributionsMinWeight, i.e., LimitContributionsMinWeight, -1000"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
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
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for LimitMaleContributions, i.e., LimitMaleContributions, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("limitmalecontributionsmin")
              if (This%LimitPar1) then
                if (allocated(Second)) then
                  This%LimitPar1Min = Char2Double(trim(adjustl(Second(1)))) ! Real because of continous solution representation
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Limit contributions of males - minimum: "//trim(Int2Char(nint(This%LimitPar1Min))) ! nint because of continous solution representation
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify a number for LimitMaleContributionsMin, i.e., LimitMaleContributionsMin, 1"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              end if

            case ("limitmalecontributionsmax")
              if (This%LimitPar1) then
                if (allocated(Second)) then
                  This%LimitPar1Max = Char2Double(trim(adjustl(Second(1)))) ! Real because of continous solution representation
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Limit contributions of males - maximum: "//trim(Int2Char(nint(This%LimitPar1Max))) ! nint because of continous solution representation
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify a number for LimitMaleContributionsMax, i.e., LimitMaleContributionsMax, 10"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              end if

            case ("limitmalecontributionsminweight")
              if (This%LimitPar1) then
                if (allocated(Second)) then
                  This%LimitPar1MinWeight = Char2Double(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Limit contributions of males - weight for contributions bellow minimum: "//trim(Real2Char(This%LimitPar1MinWeight, fmt=FMTREAL2CHAR))
                  end if
                  if (This%LimitPar1MinWeight .gt. 0.0d0) then
                    write(STDOUT, "(a)") " NOTE: Positive weight for limit on minimum contributions, i.e., encourage smaller contributions than defined minimum. Was this intended?"
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify a value for LimitMaleContributionsMinWeight, i.e., LimitMaleContributionsMinWeight, -1000"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              end if

            case ("limitfemalecontributions")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%LimitPar  = .true.
                  This%LimitPar2 = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Limit contributions of females"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for LimitFemaleContributions, i.e., LimitFemaleContributions, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("limitfemalecontributionsmin")
              if (This%LimitPar2) then
                if (allocated(Second)) then
                  This%LimitPar2Min = Char2Double(trim(adjustl(Second(1)))) ! Real because of continous solution representation
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Limit contributions of females - minimum: "//trim(Int2Char(nint(This%LimitPar2Min))) ! nint because of continous solution representation
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify a number for LimitFemaleContributionsMin, i.e., LimitFemaleContributionsMin, 1"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              end if

            case ("limitfemalecontributionsmax")
              if (This%LimitPar2) then
                if (allocated(Second)) then
                  This%LimitPar2Max = Char2Double(trim(adjustl(Second(1)))) ! Real because of continous solution representation
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Limit contributions of females - maximum: "//trim(Int2Char(nint(This%LimitPar2Max))) ! nint because of continous solution representation
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify a number for LimitFemaleContributionsMax, i.e., LimitFemaleContributionsMax, 10"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              end if

            case ("limitfemalecontributionsminweight")
              if (This%LimitPar2) then
                if (allocated(Second)) then
                  This%LimitPar2MinWeight = Char2Double(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Limit contributions of females - weight for contributions bellow minimum:: "//trim(Real2Char(This%LimitPar2MinWeight, fmt=FMTREAL2CHAR))
                  end if
                  if (This%LimitPar2MinWeight .gt. 0.0d0) then
                    write(STDOUT, "(a)") " NOTE: Positive weight for limit on minimum contributions, i.e., encourage smaller contributions than defined minimum. Was this intended?"
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify a value for LimitFemaleContributionsMinWeight, i.e., LimitFemaleContributionsMinWeight, -1000"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              end if

            case ("allowselfing")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%SelfingAllowed = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Allow selfing"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for AllowSelfing, i.e., AllowSelfing, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("selfingweight")
              if (This%SelfingAllowed) then
                if (allocated(Second)) then
                  This%SelfingWeight = Char2Double(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") "Selfing - weight: "//trim(Real2Char(This%SelfingWeight, fmt=FMTREAL2CHAR))
                  end if
                  if (This%SelfingWeight .gt. 0.0d0) then
                    write(STDOUT, "(a)") " NOTE: Positive weight for selfing, i.e., encourage selfing. Was this intended?"
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify a value for SelfingWeight, i.e., SelfingWeight, -1000"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
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
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for PAGE, i.e., PAGE, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("pagemax")
              if (This%PAGEPar) then
                if (allocated(Second)) then
                  This%PAGEParMax = Char2Int(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Promotion of Alleles by Genome Editing (PAGE) - maximum number of individuals: "//trim(Int2Char(This%PAGEParMax))
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify number for for PAGEMax, i.e., PAGEMax, 10"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              end if

            case ("pagecost")
              if (This%PAGEPar) then
                if (allocated(Second)) then
                  This%PAGEParCost = Char2Double(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Promotion of Alleles by Genome Editing (PAGE) - cost: "//trim(Real2Char(This%PAGEParCost, fmt=FMTREAL2CHAR))
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify a value for PAGECost, i.e., PAGECost, Value"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
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
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for PAGEMales, i.e., PAGEMales, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("pagemalesmax")
              if (This%PAGEPar1) then
                if (allocated(Second)) then
                  This%PAGEPar1Max = Char2Int(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Promotion of Alleles by Genome Editing (PAGE) in males - maxium number of individuals : "//trim(Int2Char(This%PAGEPar1Max))
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify number for for PAGEMalesMax, i.e., PAGEMalesMax, 10"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              end if

            case ("pagemalescost")
              if (This%PAGEPar1) then
                if (allocated(Second)) then
                  This%PAGEPar1Cost = Char2Double(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Promotion of Alleles by Genome Editing (PAGE) in males - cost: "//trim(Real2Char(This%PAGEPar1Cost, fmt=FMTREAL2CHAR))
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify a value for PAGEMalesCost, i.e., PAGEMalesCost, Value"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
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
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for PAGEFemales, i.e., PAGEFemales, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("pagefemalesmax")
              if (This%PAGEPar2) then
                if (allocated(Second)) then
                  This%PAGEPar2Max = Char2Int(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Promotion of Alleles by Genome Editing (PAGE) in females - maximum number of individuals: "//trim(Int2Char(This%PAGEPar2Max))
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify number for for PAGEFemalesMax, i.e., PAGEFemalesMax, 10"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              end if

            case ("pagefemalescost")
              if (This%PAGEPar2) then
                if (allocated(Second)) then
                  This%PAGEPar2Cost = Char2Double(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Promotion of Alleles by Genome Editing (PAGE) in females - cost: "//trim(Real2Char(This%PAGEPar2Cost, fmt=FMTREAL2CHAR))
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify a value for PAGEFemalesCost, i.e., PAGEFemalesCost, Value"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              end if

            ! Algorithm specifications
            case ("evolalgnumberofsolutions")
              if (allocated(Second)) then
                This%EvolAlgNSol = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Evolutionary algorithm - number of solutions: "//trim(Int2Char(This%EvolAlgNSol))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for EvolAlgNumberOfSolutions, i.e., EvolAlgNumberOfSolutions, 100"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("evolalgnumberofiterations")
              if (allocated(Second)) then
                This%EvolAlgNIter = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Evolutionary algorithm - number of iterations: "//trim(Int2Char(This%EvolAlgNIter))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for EvolAlgNumberOfIterations, i.e., EvolAlgNumberOfIterations, 10000"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("evolalgnumberofiterationsprint")
              if (allocated(Second)) then
                This%EvolAlgNIterPrint = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Evolutionary algorithm - number of iterations to print optimisation status: "//trim(Int2Char(This%EvolAlgNIterPrint))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for EvolAlgNumberOfIterationsPrint, i.e., EvolAlgNumberOfIterationsPrint, 100"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("evolalgnumberofiterationsstop")
              if (allocated(Second)) then
                This%EvolAlgNIterStop = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Evolutionary algorithm - number of iterations to stop upon no improvement of objective: "//trim(Int2Char(This%EvolAlgNIterStop))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for EvolAlgNumberOfIterationsStop, i.e., EvolAlgNumberOfIterationsStop, 100"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("evolalgstoptolerance")
              if (allocated(Second)) then
                This%EvolAlgStopTol = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Evolutionary algorithm - stopping tolerance: "//trim(Real2Char(This%EvolAlgStopTol, fmt=FMTREAL2CHAR))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for EvolAlgStopTolerance, i.e., EvolAlgStopTolerance, 0.01"
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
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for EvolAlgLogAllSolutions, i.e., EvolAlgLogAllSolutions, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("evolalg")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .ne. "none") then
                  write(This%EvolAlg, *) trim(adjustl(Second(1)))
                  This%EvolAlg = adjustl(This%EvolAlg)
                  if (.not. (This%EvolAlg .eq. "DE")) then
                    write(STDERR, "(a)") " ERROR: Must specify a valid algorithm [DE/???] for EvolAlg, i.e., EvolAlg, DE"
                    write(STDERR, "(a)") " "
                    stop 1
                  end if
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Evolutionary algorithm: "//trim(This%EvolAlg)
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a algorithm for EvolAlg, i.e., EvolAlg, DE"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("diffevolnumberofiterationsburnin")
              if (allocated(Second)) then
                This%DiffEvolNIterBurnIn = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Differential evolution algorithm - number of warming iterations (burn-in): "//trim(Int2Char(This%DiffEvolNIterBurnIn))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for EvolAlgNumberOfIterationsBurnin, i.e., EvolAlgNumberOfIterationsBurnin, 1000"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("diffevolparametercrburnin")
              if (allocated(Second)) then
                This%DiffEvolParamCrBurnIn = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Differential evolution algorithm - cross-over parameter for warmup (burn-in): "//trim(Real2Char(This%DiffEvolParamCrBurnIn, fmt=FMTREAL2CHAR))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for DiffEvolParameterCrBurnin, i.e., DiffEvolParameterCrBurnin, 0.4"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("diffevolparametercr")
              if (allocated(Second)) then
                This%DiffEvolParamCr = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Differential evolution algorithm - cross-over parameter: "//trim(Real2Char(This%DiffEvolParamCr, fmt=FMTREAL2CHAR))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for DiffEvolParameterCr, i.e., DiffEvolParameterCr, 0.2"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("diffevolparameterfbase")
              if (allocated(Second)) then
                This%DiffEvolParamFBase = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Differential evolution algorithm - parameter F (base value): "//trim(Real2Char(This%DiffEvolParamFBase, fmt=FMTREAL2CHAR))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for DiffEvolParameterFBase, i.e., DiffEvolParameterFBase, 0.1"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("diffevolparameterfhigh1")
              if (allocated(Second)) then
                This%DiffEvolParamFHigh1 = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Differential evolution algorithm - parameter F (high value 1): "//trim(Real2Char(This%DiffEvolParamFHigh1, fmt=FMTREAL2CHAR))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for DiffEvolParameterFHigh1, i.e., DiffEvolParameterFHigh1, 1.0"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("diffevolparameterfhigh2")
              if (allocated(Second)) then
                This%DiffEvolParamFHigh2 = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Differential evolution algorithm - parameter F (high value 2): "//trim(Real2Char(This%DiffEvolParamFHigh2, fmt=FMTREAL2CHAR))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for DiffEvolParameterFHigh2, i.e., DiffEvolParameterFHigh2, 4.0"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("randomsearchstricter")
              if (allocated(Second)) then
                This%RanAlgStricter = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Random search algorithm - perform k times more iterations than with the evolutionary algorithm: k="//trim(Int2Char(This%RanAlgStricter))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for RandomSearchStricter, i.e., RandomSearchStricter, 10"
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

      if (This%ModeOpt) then
        This%ModeMinCoancestry = .true.
        This%ModeMaxCriterion  = .true.
      end if

      if (This%EvaluateFrontier) then
        This%ModeMinCoancestry = .true.
        This%ModeMaxCriterion =  .true.
      end if

      if (.not. (This%ModeMinCoancestry .or. &
                 This%ModeMinInbreeding .or. &
                 This%ModeMaxCriterion  .or. &
                 This%ModeOpt           .or. &
                 This%ModeRan)) then
        write(STDERR, "(a)") " ERROR: One of the modes must be activated!"
        write(STDERR, "(a)") " ERROR: ModeMinCoancestry, ModeMinInbreeding, ModeMaxCriterion, ModeOpt, or ModeRan"
        write(STDERR, "(a)") " "
        stop 1
      end if

      if (.not. This%SelCriterionGiven .and. &
          (This%ModeMaxCriterion .or. This%ModeOpt)) then
        write(STDERR, "(a)") " ERROR: Selection criterion is needed for modes: ModeMaxCriterion or ModeOpt!"
        write(STDERR, "(a)") " "
        stop 1
      end if

      if (This%TargetDegreeGiven          .or. &
          This%TargetSelCriterionGiven    .or. &
          This%TargetSelIntensityGiven    .or. &
          This%TargetMaxCriterionPctGiven .or. &
          This%TargetMinCoancestryPctGiven) then
        This%ModeMinCoancestry = .true.
        This%ModeMaxCriterion  = .true.
      end if

      if (This%ModeOpt .and. .not. (This%TargetDegreeGiven          .or. This%TargetSelCriterionGiven .or. This%TargetSelIntensityGiven   .or. &
                                    This%TargetMaxCriterionPctGiven .or. This%TargetCoancestryGiven   .or. This%TargetCoancestryRateGiven .or. &
                                    This%TargetMinCoancestryPctGiven)) then
        write(STDERR, "(a)") " ERROR: One of targets must be provided when ModeOpt is activated!"
        write(STDERR, "(a)") " ERROR: TargetDegree, TargetSelCriterion, TargetSelIntensity, TargetMaxCriterionPct,"
        write(STDERR, "(a)") " ERROR: TargetCoancestry, TargetCoancestryRate, or TargetMinCoancestryPct"
        write(STDERR, "(a)") " "
        stop 1
      end if

      if (This%TargetInbreedingGiven .or. This%TargetInbreedingRateGiven) then
        This%ModeMinInbreeding = .true.
      end if

      if (.not. This%RelMtxGiven) then
        write(STDERR, "(a)") " ERROR: One of CoancestryMatrixFile or NrmMatrixFile must be specified!"
        write(STDERR, "(a)") " "
        stop 1
      end if

      if (This%LimitPar .and. This%EqualizePar) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " NOTE: The specification Equalize*Contributions has priority over Limit*Contributions."
          write(STDOUT, "(a)") " "
        end if
        ! ... therefore reset all limit specifications to default values
        This%LimitPar  = .false.
        This%LimitPar1 = .false.
        This%LimitPar2 = .false.
        This%LimitParMin  = 1.0d0
        This%LimitPar1Min = 1.0d0
        This%LimitPar2Min = 1.0d0
        This%LimitParMax  = huge(This%LimitParMax)  - 1.0d0
        This%LimitPar1Max = huge(This%LimitPar1Max) - 1.0d0
        This%LimitPar2Max = huge(This%LimitPar2Max) - 1.0d0
        This%LimitParMinWeight  = -1000.0d0
        This%LimitPar1MinWeight = -1000.0d0
        This%LimitPar2MinWeight = -1000.0d0
      end if

      if (.not. This%GenderGiven) then
        This%nPar1 = This%nPar
        This%EqualizePar1 = This%EqualizePar
        This%LimitPar1 = This%LimitPar
        This%PAGEPar1 = This%PAGEPar
      else
        This%nPar = This%nPar1 + This%nPar2
        if (This%EqualizePar) then
          if (.not. This%EqualizePar1) then
            This%EqualizePar1 = This%EqualizePar
          end if
          if (.not. This%EqualizePar2) then
            This%EqualizePar2 = This%EqualizePar
          end if
        end if
        if (This%LimitPar) then
          if (.not. This%LimitPar1) then
            This%LimitPar1 = This%LimitPar
            This%LimitPar1Min = This%LimitParMin
            This%LimitPar1Max = This%LimitParMax
            This%LimitPar1MinWeight = This%LimitParMinWeight
          end if
          if (.not. This%LimitPar2) then
            This%LimitPar2 = This%LimitPar
            This%LimitPar2Min = This%LimitParMin
            This%LimitPar2Max = This%LimitParMax
            This%LimitPar2MinWeight = This%LimitParMinWeight
          end if
        end if
        if (This%PAGEPar) then
          if (.not. This%PAGEPar1) then
            This%PAGEPar1 = This%PAGEPar
            This%PAGEPar1Max = This%PAGEParMax
            This%PAGEPar1Cost = This%PAGEParCost
          end if
          if (.not. This%PAGEPar2) then
            This%PAGEPar2 = This%PAGEPar
            This%PAGEPar2Max = This%PAGEParMax
            This%PAGEPar2Cost = This%PAGEParCost
          end if
        end if
      end if

      if (This%nMat .le. 0) then
        write(STDERR, "(a)") " ERROR: Number of matings must be larger than zero!"
        write(STDERR, "(a)") " "
        stop 1
      end if

      if (This%nPar .le. 0) then
        write(STDERR, "(a)") " ERROR: Number of parents must be larger than zero!"
        write(STDERR, "(a)") " "
        stop 1
      end if

      if (This%GenderGiven .and. ((This%nPar1 .le. 0) .or. (This%nPar2 .le. 0))) then
        write(STDERR, "(a)") " ERROR: Number of parents must be larger than zero!"
        write(STDERR, "(a)") " ERROR: Number of   male parents: "//trim(Int2Char(This%nPar1))
        write(STDERR, "(a)") " ERROR: Number of female parents: "//trim(Int2Char(This%nPar2))
        write(STDERR, "(a)") " "
        stop 1
      end if

      if (This%nMat > This%nPar) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " NOTE: Number of matings is larger than the number of parents! Was this intented?"
          write(STDOUT, "(a)") " NOTE: Number of matings: "//trim(Int2Char(This%nMat))
          write(STDOUT, "(a)") " NOTE: Number of parents: "//trim(Int2Char(This%nPar))
          write(STDOUT, "(a)") " "
        end if
      end if

      if (This%GenderGiven .and. This%SelfingAllowed) then
        write(STDERR, "(a)") " ERROR: When gender matters, AlphaMate can not perform selfing! See the manual for a solution."
        ! @todo: what is the solution? Provide the same individual both as male and a female?
        write(STDERR, "(a)") " "
        stop 1
      end if

      if ((.not. This%SelCriterionGiven) .and. This%PAGEPar) then
        write(STDERR, "(a)") " ERROR: Can not use PAGE when selection criterion file is not given!"
        ! @todo: what about using the GenericIndCrit values?
        write(STDERR, "(a)") " "
        stop 1
      end if

      if (This%PAGEPar) then
        if (This%GenderGiven) then
          if (This%PAGEPar1Max .gt. This%nPar1) then
            write(STDERR, "(a)") " ERROR: Can not PAGE more males than there are male parents!"
            write(STDERR, "(a)") " ERROR: Number of      male parents: "//trim(Int2Char(This%nPar1))
            write(STDERR, "(a)") " ERROR: Max number of male for PAGE: "//trim(Int2Char(This%PAGEPar1Max))
            write(STDERR, "(a)") " "
          end if
          if (This%PAGEPar2Max .gt. This%nPar2) then
            write(STDERR, "(a)") " ERROR: Can not PAGE more females than there are female parents!"
            write(STDERR, "(a)") " ERROR: Number of      female parents: "//trim(Int2Char(This%nPar2))
            write(STDERR, "(a)") " ERROR: Max number of female for PAGE: "//trim(Int2Char(This%PAGEPar2Max))
            write(STDERR, "(a)") " "
          end if
        else
          if (This%PAGEParMax .gt. This%nPar) then
            write(STDERR, "(a)") " ERROR: Can not PAGE more individuals than there are parents!"
            write(STDERR, "(a)") " ERROR: Number of                  parents: "//trim(Int2Char(This%nPar))
            write(STDERR, "(a)") " ERROR: Max number of individuals for PAGE: "//trim(Int2Char(This%PAGEParMax))
            write(STDERR, "(a)") " "
          end if
        end if
      end if

      if (This%SeedFileGiven .and. This%SeedGiven) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " NOTE: The specification Seed has priority over SeedFile."
          write(STDOUT, "(a)") " "
        end if
        This%SeedFile = ""
        This%SeedFileGiven = .false.
      end if
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Setup AlphaMate optimisation mode
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   April 16, 2017
    !---------------------------------------------------------------------------
    pure subroutine SetupModeAlphaMateSpec(This, Mode, Data, Degree, &
                                           SelCriterion, SelIntensity,   MaxCriterionPct,                         ModeMaxCriterionSpec,&
                                           Coancestry,   CoancestryRate, MinCoancestryPct, CoancestryWeightBelow, ModeMinCoancestrySpec, &
                                           Inbreeding,   InbreedingRate, MinInbreedingPct, InbreedingWeightBelow, ModeMinInbreedingSpec)
      implicit none
      class(AlphaMateSpec), intent(inout)           :: This                  !< @return AlphaMateSpec holder
      character(len=*), intent(in)                  :: Mode                  !< Mode definition/name
      type(AlphaMateData), intent(in), optional     :: Data                  !< AlphaMateData holder
      real(real64), intent(in), optional            :: Degree                !< Targeted degree
      real(real64), intent(in), optional            :: SelCriterion          !< Targeted selection criterion
      real(real64), intent(in), optional            :: SelIntensity          !< Targeted selection intensity
      real(real64), intent(in), optional            :: MaxCriterionPct       !< Targeted maximum criterion percentage
      real(real64), intent(in), optional            :: Coancestry            !< Targeted coancestry
      real(real64), intent(in), optional            :: CoancestryRate        !< Targeted coancestry rate
      real(real64), intent(in), optional            :: MinCoancestryPct      !< Targeted minimum coancestry percentage
      type(AlphaMateModeSpec), intent(in), optional :: ModeMaxCriterionSpec  !< Maximum criterion  solution specs
      logical, intent(in), optional                 :: CoancestryWeightBelow !< Weight deviations bellow the targeted coancestry
      type(AlphaMateModeSpec), intent(in), optional :: ModeMinCoancestrySpec !< Minimum coancestry solution specs
      real(real64), intent(in), optional            :: Inbreeding            !< Targeted inbreeding
      real(real64), intent(in), optional            :: InbreedingRate        !< Targeted inbreeding rate
      real(real64), intent(in), optional            :: MinInbreedingPct      !< Targeted minimum inbreeding percentage
      logical, intent(in), optional                 :: InbreedingWeightBelow !< Weight deviations bellow the targeted inbreeding
      type(AlphaMateModeSpec), intent(in), optional :: ModeMinInbreedingSpec !< Minimum inbreeding solution specs

      select case (trim(Mode))
        case ("MinCoancestry") ! Only coancestry!!!
          call This%ModeMinCoancestrySpec%Initialise(Name="MinCoancestry")
          This%ModeMinCoancestrySpec%ObjectiveCoancestry = .true.
          !@todo Do these lines still make sense when we go above two objectives?
          This%ModeMinCoancestrySpec%TargetDegree           =  90.0d0
          This%ModeMinCoancestrySpec%TargetMinCoancestryPct = 100.0d0
          This%ModeMinCoancestrySpec%TargetMinInbreedingPct =   0.0d0
          This%ModeMinCoancestrySpec%TargetMaxCriterionPct  =   0.0d0
          call This%ModeSpec%Assign(In=This%ModeMinCoancestrySpec)

        case ("MinInbreeding") ! Only inbreeding!!!
          call This%ModeMinInbreedingSpec%Initialise(Name="MinInbreeding")
          This%ModeMinInbreedingSpec%ObjectiveInbreeding = .true.
          !@todo Do these lines still make sense when we go above two objectives?
          This%ModeMinInbreedingSpec%TargetDegree           =  45.0d0
          This%ModeMinInbreedingSpec%TargetMinCoancestryPct =   0.0d0
          This%ModeMinInbreedingSpec%TargetMinInbreedingPct = 100.0d0
          This%ModeMinInbreedingSpec%TargetMaxCriterionPct  =   0.0d0
          call This%ModeSpec%Assign(In=This%ModeMinInbreedingSpec)

        case ("MaxCriterion") ! Only criterion!!!
          call This%ModeMaxCriterionSpec%Initialise(Name="MaxCriterion")
          This%ModeMaxCriterionSpec%ObjectiveCriterion = .true.
          !@todo Do these lines still make sense when we go above two objectives?
          This%ModeMaxCriterionSpec%TargetDegree           =   0.0d0
          This%ModeMaxCriterionSpec%TargetMinCoancestryPct =   0.0d0
          This%ModeMaxCriterionSpec%TargetMinInbreedingPct =   0.0d0
          This%ModeMaxCriterionSpec%TargetMaxCriterionPct  = 100.0d0
          call This%ModeSpec%Assign(In=This%ModeMaxCriterionSpec)

        case ("Opt") ! All objectives jointly!!!
          call This%ModeOptSpec%Initialise(Name="Opt")
          This%ModeOptSpec%ObjectiveCriterion  = .true.
          This%ModeOptSpec%ObjectiveCoancestry = .true.

          if      (present(Degree)) then
            call This%ModeOptSpec%SetTargets(Degree=Degree, &
                                             Data=Data, ModeMinCoancestrySpec=ModeMinCoancestrySpec, ModeMaxCriterionSpec=ModeMaxCriterionSpec)
          else if (present(SelCriterion)) then
            call This%ModeOptSpec%SetTargets(SelCriterion=SelCriterion, &
                                             Data=Data, ModeMinCoancestrySpec=ModeMinCoancestrySpec, ModeMaxCriterionSpec=ModeMaxCriterionSpec)
          else if (present(SelIntensity)) then
            call This%ModeOptSpec%SetTargets(SelIntensity=SelIntensity, &
                                             Data=Data, ModeMinCoancestrySpec=ModeMinCoancestrySpec, ModeMaxCriterionSpec=ModeMaxCriterionSpec)
          else if (present(MaxCriterionPct)) then
            call This%ModeOptSpec%SetTargets(MaxCriterionPct=MaxCriterionPct, &
                                             Data=Data, ModeMinCoancestrySpec=ModeMinCoancestrySpec, ModeMaxCriterionSpec=ModeMaxCriterionSpec)
          else if (present(Coancestry)) then
            call This%ModeOptSpec%SetTargets(Coancestry=Coancestry, &
                                             Data=Data, ModeMinCoancestrySpec=ModeMinCoancestrySpec, ModeMaxCriterionSpec=ModeMaxCriterionSpec)
          else if (present(CoancestryRate)) then
            call This%ModeOptSpec%SetTargets(CoancestryRate=CoancestryRate, &
                                             Data=Data, ModeMinCoancestrySpec=ModeMinCoancestrySpec, ModeMaxCriterionSpec=ModeMaxCriterionSpec)
          else if (present(MinCoancestryPct)) then
            call This%ModeOptSpec%SetTargets(MinCoancestryPct=MinCoancestryPct, &
                                             Data=Data, ModeMinCoancestrySpec=ModeMinCoancestrySpec, ModeMaxCriterionSpec=ModeMaxCriterionSpec)
          end if

          if (present(CoancestryWeightBelow)) then
            This%ModeOptSpec%CoancestryWeightBelow =      CoancestryWeightBelow
          else
            This%ModeOptSpec%CoancestryWeightBelow = This%CoancestryWeightBelow
          end if

          if      (This%TargetInbreedingGiven) then
            This%ModeOptSpec%ObjectiveInbreeding = .true.
            call This%ModeOptSpec%SetTargets(Inbreeding=This%TargetInbreeding, &
                                             Data=Data, ModeMinInbreedingSpec=ModeMinInbreedingSpec)
          else if (This%TargetInbreedingRateGiven) then
            This%ModeOptSpec%ObjectiveInbreeding = .true.
            call This%ModeOptSpec%SetTargets(InbreedingRate=This%TargetInbreedingRate, &
                                             Data=Data, ModeMinInbreedingSpec=ModeMinInbreedingSpec)
          else if (This%TargetMinInbreedingPctGiven) then
            This%ModeOptSpec%ObjectiveInbreeding = .true.
            call This%ModeOptSpec%SetTargets(MinInbreedingPct=This%TargetMinInbreedingPct, &
                                             Data=Data, ModeMinInbreedingSpec=ModeMinInbreedingSpec)
          end if

          if (present(InbreedingWeightBelow)) then
            This%ModeOptSpec%InbreedingWeightBelow =      InbreedingWeightBelow
          else
            This%ModeOptSpec%InbreedingWeightBelow = This%InbreedingWeightBelow
          end if

          call This%ModeSpec%Assign(In=This%ModeOptSpec)

        case ("Ran")
          !@todo do we need ModeRanSpec?
          call This%ModeRanSpec%Initialise(Name="Ran")
          !@todo???
          This%ModeRanSpec%CoancestryWeightBelow = This%CoancestryWeightBelow
          !@todo???
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
      real(real64) :: NANREAL64
      NANREAL64 = IEEE_Value(x=NANREAL64, class=IEEE_Quiet_NaN)
      This%Name = Name
      This%ObjectiveCriterion = .false.
      This%ObjectiveCoancestry = .false.
      This%ObjectiveInbreeding = .false.
      This%TargetDegree = NANREAL64
      This%TargetSelCriterion = NANREAL64
      This%TargetSelIntensity = NANREAL64
      This%TargetMaxCriterionPct = NANREAL64
      This%TargetCoancestry = NANREAL64
      This%TargetCoancestryRate = NANREAL64
      This%CoancestryWeightBelow = .false.
      This%TargetMinCoancestryPct = NANREAL64
      This%TargetInbreeding = NANREAL64
      This%TargetInbreedingRate = NANREAL64
      This%InbreedingWeightBelow = .false.
      This%TargetMinInbreedingPct = NANREAL64
      This%Degree = NANREAL64
      This%SelCriterion = NANREAL64
      This%SelIntensity = NANREAL64
      This%MaxCriterionPct = NANREAL64
      This%Coancestry = NANREAL64
      This%CoancestryRate = NANREAL64
      This%MinCoancestryPct = NANREAL64
      This%Inbreeding = NANREAL64
      This%InbreedingRate = NANREAL64
      This%MinInbreedingPct = NANREAL64
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
      Out%TargetSelIntensity = In%TargetSelIntensity
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
      Out%SelIntensity = In%SelIntensity
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
      Out%SelIntensity = In%SelIntensity
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
    pure subroutine SetTargetsAlphaMateModeSpec(This, Data, Degree, &
                                                SelCriterion, SelIntensity, MaxCriterionPct, &
                                                Coancestry, CoancestryRate, MinCoancestryPct, &
                                                Inbreeding, InbreedingRate, MinInbreedingPct, &
                                                ModeMinCoancestrySpec, ModeMinInbreedingSpec, ModeMaxCriterionSpec)
      implicit none
      class(AlphaMateModeSpec), intent(inout)       :: This                  !< @return AlphaMateModeSpec holder
      type(AlphaMateData), intent(in), optional     :: Data                  !< AlphaMateData holder
      real(real64), intent(in), optional            :: Degree                !< Targeted selection/Coancestry frontier degree
      real(real64), intent(in), optional            :: SelCriterion          !< Targeted selection criterion
      real(real64), intent(in), optional            :: SelIntensity          !< Targeted selection intensity
      real(real64), intent(in), optional            :: MaxCriterionPct       !< Targeted maximum criterion percentage (100 means the maximum possible selection intensity)
      real(real64), intent(in), optional            :: Coancestry            !< Targeted coancestry
      real(real64), intent(in), optional            :: CoancestryRate        !< Targeted coancestry rate
      real(real64), intent(in), optional            :: MinCoancestryPct      !< Targeted minimum coancestry percentage (100 means the minimum possible coancestry)
      real(real64), intent(in), optional            :: Inbreeding            !< Targeted inbreeding
      real(real64), intent(in), optional            :: InbreedingRate        !< Targeted inbreeding rate
      real(real64), intent(in), optional            :: MinInbreedingPct      !< Targeted minimum inbreeding percentage (100 means the minimum possible inbreeding)
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
          This%TargetSelIntensity = MaxCriterionPct2SelIntensity(MaxCriterionPct=This%TargetMaxCriterionPct, &
                                                                 MinSelIntensity=ModeMinCoancestrySpec%SelIntensity, &
                                                                 MaxSelIntensity=ModeMaxCriterionSpec%SelIntensity)
          if (present(Data)) then
            This%TargetCoancestry = CoancestryRate2Coancestry(CoancestryRate=This%TargetCoancestryRate, &
                                                              CurrentCoancestry=Data%CoancestryRanMate)
            This%TargetSelCriterion = SelIntensity2SelCriterion(SelIntensity=This%TargetSelIntensity, &
                                                                Mean=Data%SelCriterionStat%Mean, &
                                                                Sd=Data%SelCriterionStat%Sd)
          end if
        end if
      else if (present(SelCriterion)) then
        This%TargetSelCriterion = SelCriterion
        if (present(Data)) then
            This%TargetSelIntensity = SelCriterion2SelIntensity(SelCriterion=This%TargetSelCriterion, &
                                                                Mean=Data%SelCriterionStat%Mean, &
                                                                Sd=Data%SelCriterionStat%Sd)
        end if
        if (present(ModeMinCoancestrySpec) .and. present(ModeMaxCriterionSpec)) then
          This%TargetMaxCriterionPct = SelIntensity2MaxCriterionPct(SelIntensity=This%TargetSelIntensity, &
                                                                    MinSelIntensity=ModeMinCoancestrySpec%SelIntensity, &
                                                                    MaxSelIntensity=ModeMaxCriterionSpec%SelIntensity)
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
      else if (present(SelIntensity)) then
        This%TargetSelIntensity = SelIntensity
        if (present(Data)) then
            This%TargetSelCriterion = SelIntensity2SelCriterion(SelIntensity=This%TargetSelIntensity, &
                                                                Mean=Data%SelCriterionStat%Mean, &
                                                                Sd=Data%SelCriterionStat%Sd)
        end if
        if (present(ModeMinCoancestrySpec) .and. present(ModeMaxCriterionSpec)) then
          This%TargetMaxCriterionPct = SelIntensity2MaxCriterionPct(SelIntensity=This%TargetSelIntensity, &
                                                                    MinSelIntensity=ModeMinCoancestrySpec%SelIntensity, &
                                                                    MaxSelIntensity=ModeMaxCriterionSpec%SelIntensity)
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
          This%TargetSelIntensity = MaxCriterionPct2SelIntensity(MaxCriterionPct=This%TargetMaxCriterionPct, &
                                                                 MinSelIntensity=ModeMinCoancestrySpec%SelIntensity, &
                                                                 MaxSelIntensity=ModeMaxCriterionSpec%SelIntensity)
          if (present(Data)) then
            This%TargetCoancestry = CoancestryRate2Coancestry(CoancestryRate=This%TargetCoancestryRate, &
                                                              CurrentCoancestry=Data%CoancestryRanMate)
            This%TargetSelCriterion = SelIntensity2SelCriterion(SelIntensity=This%TargetSelIntensity, &
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
          This%TargetSelIntensity = MaxCriterionPct2SelIntensity(MaxCriterionPct=This%TargetMaxCriterionPct, &
                                                                 MinSelIntensity=ModeMinCoancestrySpec%SelIntensity, &
                                                                 MaxSelIntensity=ModeMaxCriterionSpec%SelIntensity)
          if (present(Data)) then
            This%TargetSelCriterion = SelIntensity2SelCriterion(SelIntensity=This%TargetSelIntensity, &
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
          This%TargetSelIntensity = MaxCriterionPct2SelIntensity(MaxCriterionPct=This%TargetMaxCriterionPct, &
                                                                 MinSelIntensity=ModeMinCoancestrySpec%SelIntensity, &
                                                                 MaxSelIntensity=ModeMaxCriterionSpec%SelIntensity)
          if (present(Data)) then
            This%TargetSelCriterion = SelIntensity2SelCriterion(SelIntensity=This%TargetSelIntensity, &
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
          This%TargetSelIntensity = MaxCriterionPct2SelIntensity(MaxCriterionPct=This%TargetMaxCriterionPct, &
                                                                 MinSelIntensity=ModeMinCoancestrySpec%SelIntensity, &
                                                                 MaxSelIntensity=ModeMaxCriterionSpec%SelIntensity)
          if (present(Data)) then
            This%TargetCoancestry = CoancestryRate2Coancestry(CoancestryRate=This%TargetCoancestryRate, &
                                                              CurrentCoancestry=Data%CoancestryRanMate)
            This%TargetSelCriterion = SelIntensity2SelCriterion(SelIntensity=This%TargetSelIntensity, &
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
                                                                        MaxCoancestryRate=1.0d0)
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
                                                                        MaxCoancestryRate=1.0d0)
        end if
      else if (present(MinInbreedingPct)) then
        This%TargetMinInbreedingPct = MinInbreedingPct
        if (present(ModeMinInbreedingSpec)) then
          This%TargetInbreedingRate = MinCoancestryPct2CoancestryRate(MinCoancestryPct=This%TargetMinInbreedingPct, &
                                                                      MinCoancestryRate=ModeMinInbreedingSpec%InbreedingRate, &
                                                                      MaxCoancestryRate=1.0d0)
          if (present(Data)) then
            This%TargetInbreeding = CoancestryRate2Coancestry(CoancestryRate=This%TargetInbreedingRate, &
                                                              CurrentCoancestry=Data%Inbreeding)
          end if
        end if
      end if
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Write an AlphaMate optimisation mode specifications
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   April 13, 2017
    !---------------------------------------------------------------------------
    subroutine WriteAlphaMateModeSpec(This, Unit)
      implicit none
      class(AlphaMateModeSpec), intent(in) :: This !< AlphaMateModeSpec holder
      integer(int32), intent(in)           :: Unit !< Unit to write to
      write(Unit, *) "Name: ", trim(This%Name)
      write(Unit, *) "ObjectiveCriterion: ", This%ObjectiveCriterion
      write(Unit, *) "ObjectiveCoancestry: ", This%ObjectiveCoancestry
      write(Unit, *) "ObjectiveInbreeding: ", This%ObjectiveInbreeding
      write(Unit, *) "TargetDegree: ", This%TargetDegree
      write(Unit, *) "TargetSelCriterion: ", This%TargetSelCriterion
      write(Unit, *) "TargetSelIntensity: ", This%TargetSelIntensity
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
      write(Unit, *) "SelIntensity: ", This%SelIntensity
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
    subroutine LogTargetsAlphaMateModeSpec(This, Spec, Unit)
      implicit none
      class(AlphaMateModeSpec), intent(in) :: This !< AlphaMateModeSpec holder
      type(AlphaMateSpec), intent(in)      :: Spec !< AlphaMateSpec holder
      integer(int32), intent(in)           :: Unit !< Unit to write to
      write(Unit, "(a)") "   Selection intensity / selection criterion"
      write(Unit, "(a)") "     @MinCoancestry: "//trim(Real2Char(Spec%ModeMinCoancestrySpec%SelIntensity, fmt=FMTREAL2CHAR))//&
                         " /"//trim(Real2Char(Spec%ModeMinCoancestrySpec%SelCriterion, fmt=FMTREAL2CHAR))
      write(Unit, "(a)") "     @MaxCriterion:  "//trim(Real2Char(Spec%ModeMaxCriterionSpec%SelIntensity, fmt=FMTREAL2CHAR))//&
                         " /"//trim(Real2Char(Spec%ModeMaxCriterionSpec%SelCriterion, fmt=FMTREAL2CHAR))
      write(Unit, "(a)") "     Pct of max:     "//trim(Real2Char(This%TargetMaxCriterionPct, fmt="(f7.1)"))
      write(Unit, "(a)") "     Target:         "//trim(Real2Char(This%TargetSelIntensity, fmt=FMTREAL2CHAR))//&
                         " /"//trim(Real2Char(This%TargetSelCriterion, fmt=FMTREAL2CHAR))
      write(Unit, "(a)") " "

      write(Unit, "(a)") "   Coancestry rate / coancestry coefficient"
      write(Unit, "(a)") "     @MinCoancestry: "//trim(Real2Char(Spec%ModeMinCoancestrySpec%CoancestryRate, fmt=FMTREAL2CHAR))//&
                         " /"//trim(Real2Char(Spec%ModeMinCoancestrySpec%Coancestry, fmt=FMTREAL2CHAR))
      write(Unit, "(a)") "     @MaxCriterion:  "//trim(Real2Char(Spec%ModeMaxCriterionSpec%CoancestryRate, fmt=FMTREAL2CHAR))//&
                         " /"//trim(Real2Char(Spec%ModeMaxCriterionSpec%Coancestry, fmt=FMTREAL2CHAR))
      write(Unit, "(a)") "     Pct of min:     "//trim(Real2Char(This%TargetMinCoancestryPct, fmt="(f7.1)"))
      write(Unit, "(a)") "     Target:         "//trim(Real2Char(This%TargetCoancestryRate, fmt=FMTREAL2CHAR))//&
                         " /"//trim(Real2Char(This%TargetCoancestry, fmt=FMTREAL2CHAR))
      write(Unit, "(a)") " "

      write(Unit, "(a)") "   Degree"
      write(Unit, "(a)") "     Target:     "//trim(Real2Char(This%TargetDegree, fmt="(f7.1)"))
      write(Unit, "(a)") " "

      write(Unit, "(a)") "   Inbreeding rate / inbreeding coefficient"
      write(Unit, "(a)") "     @MinInbreeding: "//trim(Real2Char(Spec%ModeMinInbreedingSpec%Inbreedingrate, fmt=FMTREAL2CHAR))//&
                         " /"//trim(Real2Char(Spec%ModeMinInbreedingSpec%Inbreeding, fmt=FMTREAL2CHAR))
      write(Unit, "(a)") "     MaxInbreeding:  "//trim(Real2Char(+1.0d0, fmt=FMTREAL2CHAR))//&
                         " / ???"
      write(Unit, "(a)") "     Pct of min:     "//trim(Real2Char(This%TargetMinInbreedingPct, fmt="(f7.1)"))
      write(Unit, "(a)") "     Target:         "//trim(Real2Char(This%TargetInbreedingRate, fmt=FMTREAL2CHAR))//&
                         " /"//trim(Real2Char(This%TargetInbreeding, fmt=FMTREAL2CHAR))
      write(Unit, "(a)") " "
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
      type(AlphaMateSpec), intent(inout) :: Spec      !< AlphaMateSpec holder (inout because we save some info also in Spec that is based on Data)
      logical, optional                  :: LogStdout !< Log process on stdout (default .false.)

      integer(int32) :: Ind, IndLoc, IndLoc2, nIndTmp, Mat, nMatTmp, GenderTmp, l, m, IndPair(2), Crit
      integer(int32) :: SelCriterionUnit, GenderUnit, GenericIndCritUnit, GenericMatCritUnit, SeedUnit
      integer(int32) :: CoancestrySummaryUnit, InbreedingSummaryUnit, CriterionSummaryUnit
      integer(int32) :: GenericIndCritSummaryUnit, GenericMatCritSummaryUnit

      real(real64) :: SelCriterionTmp, SelCriterionTmp2
      real(real64), allocatable :: GenericIndCritTmp(:), GenericMatCritTmp(:)

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

      if (LogStdoutInternal) then
        write(STDOUT, "(a)") " Number of individuals in the coancestry matrix file: "//trim(Int2Char(This%nInd))
      end if

      if (This%nInd .lt. Spec%nPar) then
        write(STDERR, "(a)") "ERROR: Number of individuals can not be smaller than number of parents!"
        write(STDERR, "(a)") "ERROR: Number of individuals: "//trim(Int2Char(This%nInd))
        write(STDERR, "(a)") "ERROR: Number of     parents: "//trim(Int2Char(Spec%nPar))
        write(STDERR, "(a)") " "
        stop 1
      end if

      ! --- Selection criterion ---

      allocate(This%SelCriterion(This%nInd))
      allocate(This%SelIntensity(This%nInd))
      if (Spec%PAGEPar) then
        allocate(This%SelCriterionPAGE(This%nInd))
        allocate(This%SelIntensityPAGE(This%nInd))
      end if

      if (.not. Spec%SelCriterionGiven) then
        This%SelCriterion = 0.0d0
        This%SelIntensity = 0.0d0
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
          IndLoc = FindLoc(IdCTmp, This%Coancestry%OriginalId(1:))
          if (IndLoc .eq. 0) then
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

      ! --- Gender ---

      allocate(This%Gender(This%nInd))
      This%Gender = 0
      if (Spec%GenderGiven) then
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
          IndLoc = FindLoc(IdCTmp, This%Coancestry%OriginalId(1:))
          if (IndLoc .eq. 0) then
            write(STDERR, "(a)") " ERROR: Individual "//trim(IdCTmp)//" from the gender file not present in the coancestry matrix file!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          This%Gender(IndLoc) = GenderTmp
        end do
        close(GenderUnit)

        write(STDOUT, "(a)") " Number of   males: "//trim(Int2Char(This%nMal))
        write(STDOUT, "(a)") " Number of females: "//trim(Int2Char(This%nFem))

        if (Spec%nPar1 > This%nMal) then
          write(STDERR, "(a)") " ERROR: Number of male parents can not be larger than number of males"
          write(STDERR, "(a)") " ERROR: Number of male parents: "//trim(Int2Char(Spec%nPar1))
          write(STDERR, "(a)") " ERROR: Number of        males: "//trim(Int2Char(This%nMal))
          write(STDERR, "(a)") " "
          stop 1
        end if
        if (Spec%nPar2 > This%nFem) then
          write(STDERR, "(a)") " ERROR: Number of female parents can not be larger than number of females"
          write(STDERR, "(a)") " ERROR: Number of female parents: "//trim(Int2Char(Spec%nPar2))
          write(STDERR, "(a)") " ERROR: Number of        females: "//trim(Int2Char(This%nFem))
          write(STDERR, "(a)") " "
          stop 1
        end if
      end if

      ! --- Define potential parents ---

      if (.not. Spec%GenderGiven) then
        This%nPotPar1 = This%nInd
        This%nPotPar2 = This%nInd
        allocate(This%IdPotPar1(This%nPotPar1))
        do Ind = 1, This%nInd
          This%IdPotPar1(Ind) = Ind
        end do
      else
        This%nPotPar1 = This%nMal
        This%nPotPar2 = This%nFem
        allocate(This%IdPotPar1(This%nPotPar1))
        allocate(This%IdPotPar2(This%nPotPar2))
        allocate(This%IdPotParSeq(This%nInd))
        block
          integer(int32) :: jMal, jFem
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
        end block
      end if

      ! --- Number of all potential matings ---

      if (Spec%GenderGiven) then
        This%nPotMat = This%nPotPar1 * This%nPotPar2
      else
        This%nPotMat = real(This%nPotPar1 * This%nPotPar1) / 2
        if (Spec%SelfingAllowed) then
          This%nPotMat = nint(This%nPotMat + real(This%nPotPar1) / 2)
        else
          This%nPotMat = nint(This%nPotMat - real(This%nPotPar1) / 2)
        end if
      end if

      if (Spec%nMat .gt. This%nPotMat) then
        ! @todo what about MOET, AI, JIVET, ... etc?
        write(STDERR, "(a)") " ERROR: Number of specified matings is larger than the number of all potential matings!"
        write(STDERR, "(a)") " ERROR: Number of all potential matings: "//trim(Int2Char(This%nPotMat))
        write(STDERR, "(a)") " ERROR: Number of     specified matings: "//trim(Int2Char(Spec%nMat))
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
        This%GenericIndCrit = 0.0d0
        open(newunit=GenericIndCritUnit, file=Spec%GenericIndCritFile, status="unknown")
        do Ind = 1, This%nInd
          read(GenericIndCritUnit, *) IdCTmp, GenericIndCritTmp
          IndLoc = FindLoc(IdCTmp, This%Coancestry%OriginalId(1:))
          if (IndLoc .eq. 0) then
            write(STDERR, "(a)") " ERROR: Individual "//trim(IdCTmp)//" from the generic individual criterion file not present in the coancestry matrix file!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          This%GenericIndCrit(IndLoc, :) = GenericIndCritTmp
        end do
        close(GenericIndCritUnit)
      end if

      ! --- Generic mating criterion ---

      if (Spec%GenericMatCritGiven) then
        nMatTmp = CountLines(Spec%GenericMatCritFile)
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " Number of matings in the generic mating criterion file: "//trim(Int2Char(nMatTmp))
        end if
        if (nMatTmp .ne. This%nPotMat) then
          write(STDERR, "(a)") " ERROR: Number of matings in the generic mating criterion file and the number of all potential matings is not the same!"
          write(STDERR, "(a)") " ERROR: Number of matings in the generic mating criterion file: "//trim(Int2Char(nMatTmp))
          write(STDERR, "(a)") " ERROR: Number of all potential matings:                        "//trim(Int2Char(This%nPotMat))
          write(STDERR, "(a)") " "
          stop 1
        end if
        allocate(This%GenericMatCrit(This%nPotPar1, This%nPotPar2, Spec%nGenericMatCrit))
        allocate(GenericMatCritTmp(Spec%nGenericMatCrit))
        This%GenericMatCrit = 0.0d0
        open(newunit=GenericMatCritUnit, file=Spec%GenericMatCritFile, status="unknown")
        do Mat = 1, This%nPotMat
          read(GenericMatCritUnit, *) IdCTmp, IdCTmp2, GenericMatCritTmp
          IndLoc = FindLoc(IdCTmp, This%Coancestry%OriginalId(1:))
          if (IndLoc .eq. 0) then
            write(STDERR, "(a)") " ERROR: Individual "//trim(IdCTmp)//" from the generic mating criterion file not present in the coancestry matrix file!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          IndLoc2 = FindLoc(IdCTmp2, This%Coancestry%OriginalId(1:))
          if (IndLoc2 .eq. 0) then
            write(STDERR, "(a)") " ERROR: Individual "//trim(IdCTmp2)//" from the generic mating criterion file not present in the coancestry matrix file!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          if (Spec%GenderGiven) then
            l = FindLoc(IndLoc, This%IdPotPar1)
            if (l .eq. 0) then
              write(STDERR, "(a)") " ERROR: Individual "//trim(IdCTmp)//" from the first column in the generic mating criterion file should be a male!"
              write(STDERR, "(a)") " ERROR: Generic mating criterion file:"
              write(STDERR, "(a)") " ERROR:   - line:                  "//trim(Int2Char(Mat))
              write(STDERR, "(a)") " ERROR:   - individual 1   (male): "//trim(IdCTmp)//" gender "//trim(Int2Char(This%Gender(IndLoc)))
              write(STDERR, "(a)") " ERROR:   - individual 2 (female): "//trim(IdCTmp2)//" gender "//trim(Int2Char(This%Gender(IndLoc2)))
              write(STDERR, "(a)") " "
              stop 1
            end if
            m = FindLoc(IndLoc2, This%IdPotPar2)
            if (m .eq. 0) then
              write(STDERR, "(a)") " ERROR: Individual "//trim(IdCTmp2)//" from the second column in the generic mating criterion file should be a female!"
              write(STDERR, "(a)") " ERROR: Generic mating criterion file:"
              write(STDERR, "(a)") " ERROR:   - line:                  "//trim(Int2Char(Mat))
              write(STDERR, "(a)") " ERROR:   - individual 1   (male): "//trim(IdCTmp)//" gender "//trim(Int2Char(This%Gender(IndLoc)))
              write(STDERR, "(a)") " ERROR:   - individual 2 (female): "//trim(IdCTmp2)//" gender "//trim(Int2Char(This%Gender(IndLoc2)))
              write(STDERR, "(a)") " "
              stop 1
            end if
            ! fill full-matrix
            ! - need to locate position in GenericMatCrit that pertains to a particular potential mating
            ! - l and m are respectively locations within IdPotPar1 and IdPotPar2
            ! - values in IdPotPar1 and IdPotPar2 are "joint" Id of males and females, i.e., they range from 1:n
            ! - values in IdPotParSeq are separate Id of males and females, i.e., one ranges from 1:nMal and the other from 1:nFem
            This%GenericMatCrit(This%IdPotParSeq(This%IdPotPar1(l)), This%IdPotParSeq(This%IdPotPar2(m)), :) = GenericMatCritTmp
          else
            ! fill lower-triangle (half-diallel)
            IndPair = [IndLoc, IndLoc2]
            This%GenericMatCrit(maxval(IndPair), minval(IndPair), :) = GenericMatCritTmp
          end if
        end do
        close(GenericMatCritUnit)
      end if

      ! --- Seed ---

      if (.not. (Spec%SeedGiven .or. Spec%SeedFileGiven)) then
        call SetSeed(Out=Spec%Seed)
      else
        if ((.not. Spec%SeedGiven) .and. Spec%SeedFileGiven) then
          open(newunit=SeedUnit, file=Spec%SeedFile, status="old")
          read(SeedUnit, *) Spec%Seed
          close(SeedUnit)
        end if
      end if
      open(newunit=SeedUnit, file="SeedUsed.txt", status="unknown")
      write(SeedUnit, *) Spec%Seed
      close(SeedUnit)
      if (LogStdoutInternal) then
        write(STDOUT, "(a)") " RNG seed: "//trim(Int2Char(Spec%Seed))
      end if

      ! --- Current coancestry summary ---

      if (LogStdoutInternal) then
        write(STDOUT, "(a)") " "
        write(STDOUT, "(a)") " Current coancestry summary (average identity of the four genome combinations of two individuals)"
      end if

      ! @todo Should we use DescStatMatrix or DescStatLowTriMatrix?
      This%CoancestryStat = DescStatMatrix(This%Coancestry%Value(1:, 1:))
      if (Spec%GenderGiven) then
        ! @todo Should we use DescStatMatrix or DescStatLowTriMatrix?
        This%CoancestryStatGender1    = DescStatMatrix(This%Coancestry%Value(This%IdPotPar1, This%IdPotPar1))
        ! @todo Should we use DescStatMatrix or DescStatLowTriMatrix?
        This%CoancestryStatGender2    = DescStatMatrix(This%Coancestry%Value(This%IdPotPar2, This%IdPotPar2))
        This%CoancestryStatGenderDiff = DescStatMatrix(This%Coancestry%Value(This%IdPotPar1, This%IdPotPar2))
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
        write(STDOUT, "(a)") "    Among   = coancestry among   individuals (including self-coancestry) = expected inbreeding in their progeny under random mating, including selfing"
        write(STDOUT, "(a)") "    Between = coancestry between individuals (excluding self-coancestry) = expected inbreeding in their progeny under random mating, excluding selfing"

        if (Spec%GenderGiven) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") "  - coancestry between males and females"
          write(STDOUT, "(a)") "    (=expected inbreeding in their progeny under random mating between genders)"
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

      ! Save means to a file
      open(newunit=CoancestrySummaryUnit, file="CoancestrySummary.txt", status="unknown")
      write(CoancestrySummaryUnit, "(a, f)") "Current (random mating),                 ",   This%CoancestryRanMate
      write(CoancestrySummaryUnit, "(a, f)") "Current (random mating, no selfing),     ",   This%CoancestryRanMateNoSelf
      if (Spec%GenderGiven) then
        write(CoancestrySummaryUnit, "(a, f)") "Current (random mating between genders), ", This%CoancestryGenderMate
      end if
      close(CoancestrySummaryUnit)

      ! --- Current inbreeding summary ---

      if (LogStdoutInternal) then
        write(STDOUT, "(a)") " "
        write(STDOUT, "(a)") " Current inbreeding summary (identity between the two genomes of an individual)"
      end if

      block
        type(InbVec) :: Inbreeding
        call This%Coancestry%Inbreeding(Out=Inbreeding, Nrm=.false.)
        This%InbreedingStat = DescStat(Inbreeding%Value(1:))
        This%Inbreeding = This%InbreedingStat%Mean
      end block

      if (LogStdoutInternal) then
        write(STDOUT, "(a)") "  - n:       "//trim( Int2Char(This%InbreedingStat%n,    fmt=FMTINT2CHAR))
        write(STDOUT, "(a)") "  - average: "//trim(Real2Char(This%InbreedingStat%Mean, fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "  - st.dev.: "//trim(Real2Char(This%InbreedingStat%Sd,   fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "  - minimum: "//trim(Real2Char(This%InbreedingStat%Min,  fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "  - maximum: "//trim(Real2Char(This%InbreedingStat%Max,  fmt=FMTREAL2CHAR))
      end if

      open(newunit=InbreedingSummaryUnit, file="InbreedingSummary.txt", status="unknown")
      write(InbreedingSummaryUnit, "(a, f)") "Current, ", This%Inbreeding
      close(InbreedingSummaryUnit)

      ! --- Current selection criterion summary ---

      if (Spec%SelCriterionGiven) then

        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") " Current selection criterion summary"
        end if

        This%SelCriterionStat = DescStat(This%SelCriterion)
        This%SelIntensity = SelCriterion2SelIntensity(SelCriterion=This%SelCriterion, &
                                                      Mean=This%SelCriterionStat%Mean, &
                                                      Sd=This%SelCriterionStat%Sd)
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") "  - n:       "//trim( Int2Char(This%SelCriterionStat%n,    fmt=FMTINT2CHAR))
          write(STDOUT, "(a)") "  - average: "//trim(Real2Char(This%SelCriterionStat%Mean, fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "  - st.dev.: "//trim(Real2Char(This%SelCriterionStat%Sd,   fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "  - minimum: "//trim(Real2Char(This%SelCriterionStat%Min,  fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "  - maximum: "//trim(Real2Char(This%SelCriterionStat%Max,  fmt=FMTREAL2CHAR))
        end if

        if (This%SelCriterionStat%Sd .eq. 0.0) then
          write(STDERR, "(a)") " ERROR: There is no variation in selection criterion!"
          write(STDERR, "(a)") " "
          stop 1
        end if

        open(newunit=CriterionSummaryUnit, file="SelCriterionSummary.txt", status="unknown")
        write(CriterionSummaryUnit, "(a, f)") "Mean, ", This%SelCriterionStat%Mean
        close(CriterionSummaryUnit)

        if (Spec%PAGEPar) then
          ! must have the same scale as selection criterion!!!!
          This%SelIntensityPAGE = SelCriterion2SelIntensity(SelCriterion=This%SelCriterionPAGE, &
                                                            Mean=This%SelCriterionStat%Mean, &
                                                            Sd=This%SelCriterionStat%Sd)
          ! only the PAGE bit of SelCriterion
          This%SelCriterionPAGE = This%SelCriterionPAGE - This%SelCriterion
          This%SelIntensityPAGE = This%SelIntensityPAGE - This%SelIntensity
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

          open(newunit=CriterionSummaryUnit, file="PAGESummary.txt", status="unknown")
          write(CriterionSummaryUnit, "(a, f)") "Mean, ", This%SelCriterionPAGEStat%Mean
          close(CriterionSummaryUnit)
        end if

      end if

      ! --- Current generic individual values ---

      if (Spec%GenericIndCritGiven) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") " Current generic individual selection criterion summary"
        end if

        open(newunit=GenericIndCritSummaryUnit, file="GenericIndCritSummary.txt", status="unknown")

        allocate(This%GenericIndCritStat(Spec%nGenericIndCrit))
        do Crit = 1, Spec%nGenericIndCrit
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") "  - criterion "//trim(Int2Char(Crit))
          This%GenericIndCritStat(Crit) = DescStat(This%GenericIndCrit(:, Crit))
          write(STDOUT, "(a)") "    - n:       "//trim( Int2Char(This%GenericIndCritStat(Crit)%n,    fmt=FMTINT2CHAR))
          write(STDOUT, "(a)") "    - average: "//trim(Real2Char(This%GenericIndCritStat(Crit)%Mean, fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(This%GenericIndCritStat(Crit)%Sd,   fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(This%GenericIndCritStat(Crit)%Min,  fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(This%GenericIndCritStat(Crit)%Max,  fmt=FMTREAL2CHAR))
          write(GenericIndCritSummaryUnit, "(a, f)") "Mean criterion "//trim(Int2Char(Crit)), This%GenericIndCritStat(Crit)%Mean
        end do

        close(GenericIndCritSummaryUnit)
      end if

      ! --- Generic mating values ---

      if (Spec%GenericMatCritGiven) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") " Generic mating selection criterion summary"
        end if

        open(newunit=GenericMatCritSummaryUnit, file="GenericMatCritSummary.txt", status="unknown")

        allocate(This%GenericMatCritStat(Spec%nGenericMatCrit))
        do Crit = 1, Spec%nGenericMatCrit
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") "  - criterion "//trim(Int2Char(Crit))
          if (Spec%GenderGiven) then
            This%GenericMatCritStat(Crit) = DescStatMatrix(This%GenericMatCrit(:, :, Crit))
            write(STDOUT, "(a)") "    - n:       "//trim( Int2Char(This%GenericMatCritStat(Crit)%All%n,    fmt=FMTINT2CHAR))
            write(STDOUT, "(a)") "    - average: "//trim(Real2Char(This%GenericMatCritStat(Crit)%All%Mean, fmt=FMTREAL2CHAR))
            write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(This%GenericMatCritStat(Crit)%All%Sd,   fmt=FMTREAL2CHAR))
            write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(This%GenericMatCritStat(Crit)%All%Min,  fmt=FMTREAL2CHAR))
            write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(This%GenericMatCritStat(Crit)%All%Max,  fmt=FMTREAL2CHAR))
            write(GenericMatCritSummaryUnit, "(a, f)") "Mean criterion "//trim(Int2Char(Crit)), This%GenericMatCritStat(Crit)%All%Mean
          else
            if (Spec%SelfingAllowed) then
              This%GenericMatCritStat(Crit) = DescStatLowTriMatrix(This%GenericMatCrit(:, :, Crit))
              write(STDOUT, "(a)") "    - n:       "//trim( Int2Char(This%GenericMatCritStat(Crit)%All%n,    fmt=FMTINT2CHAR))
              write(STDOUT, "(a)") "    - average: "//trim(Real2Char(This%GenericMatCritStat(Crit)%All%Mean, fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(This%GenericMatCritStat(Crit)%All%Sd,   fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(This%GenericMatCritStat(Crit)%All%Min,  fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(This%GenericMatCritStat(Crit)%All%Max,  fmt=FMTREAL2CHAR))
              write(GenericMatCritSummaryUnit, "(a, f)") "Mean criterion "//trim(Int2Char(Crit)), This%GenericMatCritStat(Crit)%All%Mean
            end if
              This%GenericMatCritStat(Crit) = DescStatLowTriMatrix(This%GenericMatCrit(:, :, Crit), Diag=.false.)
              write(STDOUT, "(a)") "    - n:       "//trim( Int2Char(This%GenericMatCritStat(Crit)%OffDiag%n,    fmt=FMTINT2CHAR))
              write(STDOUT, "(a)") "    - average: "//trim(Real2Char(This%GenericMatCritStat(Crit)%OffDiag%Mean, fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(This%GenericMatCritStat(Crit)%OffDiag%Sd,   fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(This%GenericMatCritStat(Crit)%OffDiag%Min,  fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(This%GenericMatCritStat(Crit)%OffDiag%Max,  fmt=FMTREAL2CHAR))
              write(GenericMatCritSummaryUnit, "(a, f)") "Mean criterion "//trim(Int2Char(Crit)), This%GenericMatCritStat(Crit)%OffDiag%Mean
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
          write(Unit, "(a, 2f)") This%Coancestry%OriginalId(Ind), This%SelCriterion(Ind), This%SelIntensity(Ind)
        end do
      else
        write(Unit, "(a)") " not allocated"
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " SelCriterionPAGE:"
      if (allocated(This%SelCriterionPAGE)) then
        do Ind = 1, This%nInd
          write(Unit, "(a, f)") This%Coancestry%OriginalId(Ind), This%SelCriterionPAGE(Ind), This%SelIntensityPAGE(Ind)
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
        write(Unit, "(a)") " Gender not allocated"
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " CoancestryStatGender2:"
      if (allocated(This%Gender)) then
        write(Unit, *) This%CoancestryStatGender2
      else
        write(Unit, "(a)") " Gender not allocated"
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " CoancestryStatGenderDiff:"
      if (allocated(This%Gender)) then
        write(Unit, *) This%CoancestryStatGenderDiff
      else
        write(Unit, "(a)") " Gender not allocated"
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
      write(Unit, "(a)") " SelCriterionStat:"
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
      write(Unit, "(a)") " nMal:"
      if (allocated(This%Gender)) then
        write(Unit, "(i)") This%nMal
      else
        write(Unit, "(a)") " Gender not allocated"
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " nFem:"
      if (allocated(This%Gender)) then
        write(Unit, "(i)") This%nFem
      else
        write(Unit, "(a)") " Gender not allocated"
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " IdPotPar1:"
      write(Unit, "(i)") This%IdPotPar1

      write(Unit, "(a)") " "
      write(Unit, "(a)") " IdPotPar2:"
      if (allocated(This%Gender)) then
        write(Unit, "(i)") This%IdPotPar2
      else
        write(Unit, "(a)") " Gender not allocated"
      end if

      write(Unit, "(a)") " "
      write(Unit, "(a)") " IdPotParSeq:"
      if (allocated(This%Gender)) then
        write(Unit, "(i)") This%IdPotParSeq
      else
        write(Unit, "(a)") " Gender not allocated"
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
        write(Unit, "(a)") " Gender not allocated"
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
    !> @brief  Initialise AlphaMate solution
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !---------------------------------------------------------------------------
    pure subroutine InitialiseAlphaMateSol(This, Chrom, Spec)
      implicit none

      ! Argument
      class(AlphaMateSol), intent(out)             :: This     !< @return AlphaMateSol holder
      real(real64), intent(in)                     :: Chrom(:) !< Provided initial solution
      class(AlphaEvolveSpec), intent(in), optional :: Spec     !< AlphaEvolveSpec --> AlphaMateSpec holder

      real(real64) :: NANREAL64
      NANREAL64 = IEEE_Value(x=NANREAL64, class=IEEE_Quiet_NaN)

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
          This%Penalty = 0.0d0
          This%PenaltyCoancestry = 0.0d0
          This%PenaltyInbreeding = 0.0d0
          This%PenaltySelfing = 0.0d0
          This%PenaltyLimitPar1 = 0.0d0
          This%PenaltyLimitPar2 = 0.0d0
          This%PenaltyGenericIndCrit = 0.0d0
          This%PenaltyGenericMatCrit = 0.0d0
          This%Degree = NANREAL64
          This%SelCriterion = 0.0d0
          This%SelIntensity = 0.0d0
          This%MaxCriterionPct = NANREAL64
          This%CoancestryRanMate = 0.0d0
          This%CoancestryRateRanMate = 0.0d0
          This%MinCoancestryPct = NANREAL64
          This%Inbreeding = 0.0d0
          This%InbreedingRate = 0.0d0
          This%MinInbreedingPct = NANREAL64
          if (Spec%GenericIndCritGiven) then
            allocate(This%GenericIndCrit(Spec%nGenericIndCrit))
            This%GenericIndCrit = 0.0d0
          end if
          if (Spec%GenericMatCritGiven) then
            allocate(This%GenericMatCrit(Spec%nGenericMatCrit))
            This%GenericMatCrit = 0.0d0
          end if
          This%Cost = 0.0d0
          allocate(This%nVec(Spec%nInd))
          This%nVec = 0
          allocate(This%xVec(Spec%nInd))
          This%xVec = 0.0d0
          allocate(This%MatingPlan(2, Spec%nMat))
          This%MatingPlan = 0
          if (Spec%PAGEPar) then
            allocate(This%GenomeEdit(Spec%nInd))
            This%GenomeEdit = 0.0d0
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
          Out%PenaltySelfing = In%PenaltySelfing
          Out%PenaltyLimitPar1 = In%PenaltyLimitPar1
          Out%PenaltyLimitPar2 = In%PenaltyLimitPar2
          Out%PenaltyGenericIndCrit = In%PenaltyGenericIndCrit
          Out%PenaltyGenericMatCrit = In%PenaltyGenericMatCrit
          Out%Degree = In%Degree
          Out%SelCriterion = In%SelCriterion
          Out%SelIntensity = In%SelIntensity
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
          Out%Cost = In%Cost
          allocate(Out%nVec(size(In%nVec)))
          Out%nVec = In%nVec
          allocate(Out%xVec(size(In%xVec)))
          Out%xVec = In%xVec
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
      real(real64) :: kR

      ! Updates
      kR = (dble(n) - 1.0d0) / n

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
          This%PenaltySelfing            = This%PenaltySelfing            * kR + Add%PenaltySelfing            / n
          This%PenaltyLimitPar1          = This%PenaltyLimitPar1          * kR + Add%PenaltyLimitPar1          / n
          This%PenaltyLimitPar2          = This%PenaltyLimitPar2          * kR + Add%PenaltyLimitPar2          / n
          This%PenaltyGenericIndCrit     = This%PenaltyGenericIndCrit     * kR + Add%PenaltyGenericIndCrit     / n
          This%PenaltyGenericMatCrit     = This%PenaltyGenericMatCrit     * kR + Add%PenaltyGenericMatCrit     / n
          This%Degree                    = This%Degree                    * kR + Add%Degree                    / n
          This%SelCriterion              = This%SelCriterion              * kR + Add%SelCriterion              / n
          This%SelIntensity              = This%SelIntensity              * kR + Add%SelIntensity              / n
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
          This%Cost                      = This%Cost                      * kR + Add%Cost                      / n
          This%nVec                      = This%nVec                      * kR + Add%nVec                      / n
          This%xVec                      = This%xVec                      * kR + Add%xVec                      / n
          This%MatingPlan                = This%MatingPlan                * kR + Add%MatingPlan                / n
          if (allocated(This%GenomeEdit)) then
            This%GenomeEdit              = This%GenomeEdit                * kR + Add%GenomeEdit                / n
          end if
          ! Constant
          ! This%TargetDegree
          ! This%TargetSelCriterion
          ! This%TargetSelIntensity
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
    !> @brief  AlphaMate evaluate function plus much MORE (this is the core!!!!)
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !---------------------------------------------------------------------------
    subroutine FixSolEtcMateAndEvaluateAlphaMateSol(This, Chrom, Spec, Data) ! not pure due to RNG
      implicit none
      ! Arguments
      class(AlphaMateSol), intent(inout)           :: This     !< @return AlphaMateSol holder(out as we fix solution)
      real(real64), intent(in)                     :: Chrom(:) !< A solution
      class(AlphaEvolveSpec), intent(in)           :: Spec     !< AlphaEvolveSpec --> AlphaMateSpec holder
      class(AlphaEvolveData), intent(in), optional :: Data     !< AlphaEvolveData --> AlphaMateData holder

      ! Other
      integer(int32) :: i, j, k, l, g, nCumMat, TmpMin, TmpMax, TmpI
      integer(int32), allocatable :: Rank(:), ChromInt(:), MatPar2(:), nVecPar1(:), nVecPar2(:)

      real(real64) :: TmpR, RanNum
      real(real64), allocatable :: TmpVec(:, :)

      select type (Spec)
        class default
          error stop " ERROR: FixSolEtcMateAndEvaluate works only with argument Spec being of type AlphaMateSpec!"
        class is (AlphaMateSpec)
          select type (Data)
            class default
              error stop " ERROR: FixSolEtcMateAndEvaluate works only with argument Data being of type AlphaMateData!"
            class is (AlphaMateData)

              call This%Initialise(Chrom=Chrom, Spec=Spec)
              This%Objective = 0.0d0

              allocate(Rank(Data%nInd))
              allocate(ChromInt(Data%nInd))
              allocate(MatPar2(Spec%nMat))
              allocate(nVecPar1(Data%nPotPar1))
              allocate(nVecPar2(Data%nPotPar2))
              allocate(TmpVec(Data%nInd, 1))

              ! The solution (based on the mate selection driver) has:
              ! - Data%nInd individual contributions
              !   - Data%nPotPar1 individual contributions for "parent1" (males   when GenderGiven, all ind when .not. GenderGiven)
              !   - Data%nPotPar2 individual contributions for "parent2" (females when GenderGiven, meaningful only when GenderGiven)
              ! - Spec%nMat rankings of parent1 1:Spec%nMat matings to pair with 1:Data%nPotPar2 "parent2" (see below)
              ! - Data%nInd edit indicators
              !   - Data%nPotPar1 edit indicators for "parent1" (males   when GenderGiven, all ind when .not. GenderGiven)
              !   - Data%nPotPar2 edit indicators for "parent2" (females when GenderGiven, present only when GenderGiven)
              !
              ! Say we have Chrom=(| 0, 2, 0, 1 | ... | 2.5, 1.5, 1.0 | 0, 1, 0, 0 | ...) then we:
              ! - mate male 2 with the first  available female (rank 2.5)
              ! - mate male 2 with the second available female (rank 1.5)
              ! - mate male 4 with the third  available female (rank 1.0)
              ! - edit male 2
              !
              ! @todo consider spliting the Chrom() vector internally into a type with
              !   separate vectors to simplify the code, e.g.,
              !   - Chrom2%ContPar1
              !   - Chrom2%ContPar2
              !   - Chrom2%MateRank
              !   - Chrom2%EditPar1
              !   - Chrom2%EditPar2
              !   and then at the end combine it back en exit - since we modify/fix some
              !   elements of a solution, we need to combine before exit!

              ! --- Parse the mate selection driver (=Is the solution valid?) ---

              ! The approach below assures that we have Spec%nMat contributions for each of
              ! the two parent sets. It does this by ranking internal solution values and
              ! traverses from top to the defined number of parents checking when the sum of
              ! interegrised values gives Spec%nMat. If values below 0.5 are found, they are
              ! changed to 1 contribution to achieve nMat. If this still does not give Spec%nMat,
              ! then we start adding on contributions to each parent (starting at those contributing
              ! the least and work up to those contributing most - to avoid local minima) until
              ! we reach Spec%nMat. How to treat values for the individuals that do not contribute
              ! is unlcear. None of the tested methods seemed to be very different. Intuitively,
              ! using properly ordered negative values should inform optim. alg. which individuals
              ! should less likely contribute, but this did not seem to be the case - better final
              ! solution was found when this strategy was not implemented - either zeroing
              ! values for those individuals (was the fastest) or giving random value (was
              ! marginally better, but slower). Potential advantage of not preserving the
              ! order is that this gives more randomness and more solutions being explored.

              ! "Parent1"
              if (Spec%GenderGiven) then
                g = 1
              else
                g = 2
              end if
              ! ... ranks to find the top contributors
              if (.not.(Spec%EqualizePar1 .and. (Spec%nPar1 .eq. Data%nPotPar1))) then
                Rank(1:Data%nPotPar1) = MrgRnk(This%Chrom(1:Data%nPotPar1))
                Rank(1:Data%nPotPar1) = Rank(Data%nPotPar1:1:-1) ! MrgRnk ranks small to large
              end if
              ! ... equal contributions
              if (Spec%EqualizePar1) then
                if (Spec%nPar1 .eq. Data%nPotPar1) then
                  ! ... set integers to all the values (no need for sorting here)
                  This%Chrom(1:Data%nPotPar1) = dble(Spec%nMat * g) / Spec%nPar1
                else
                  ! ... set integers to the top values
                  This%Chrom(Rank(1:Spec%nPar1)) = dble(Spec%nMat * g) / Spec%nPar1
                  ! @todo a better way to preserve the order of non contributing individuals? See below!
                  This%Chrom(Rank((Spec%nPar1 + 1):Data%nPotPar1)) = 0.0d0
                  ! This%Chrom(Rank((Spec%nPar1+1):Data%nPotPar1)) = -1.0d0
                end if
              else
                ! ... unequal contributions
                ! ... work for the defined number or parents
                nCumMat = 0
                do i = 1, Spec%nPar1
                  j = Rank(i)
                  ! ... these would have "zero" contributions and if we hit them then they need at least minimum contrib
                  if (This%Chrom(j) .lt. Spec%LimitPar1Min) then
                    This%Chrom(j) = Spec%LimitPar1Min ! set/fix to minimum usage @todo could consider penalising solution instead?
                  end if
                  ! ... but not above max allowed
                  if (This%Chrom(j) .gt. Spec%LimitPar1Max) then
                    This%Chrom(j) = Spec%LimitPar1Max ! set/fix to maximum usage @todo could consider penalising solution instead?
                  end if
                  ! ... accumulate and check if we reached Spec%nMat
                  nCumMat = nCumMat + nint(This%Chrom(j)) ! internally real, externally integer
                  if (nCumMat .ge. Spec%nMat * g) then
                    ! ... there should be exactly Spec%nMat contributions
                    if (nCumMat .gt. Spec%nMat * g) then
                      This%Chrom(j) = This%Chrom(j) - dble(nCumMat - Spec%nMat * g)
                      ! ... did we go below the minimum usage limit?
                      if (nint(This%Chrom(j)) .lt. Spec%LimitPar1Min) then
                        TmpR = Spec%LimitPar1MinWeight * (Spec%LimitPar1Min - nint(This%Chrom(j)))
                        This%Objective = This%Objective + TmpR
                        if (Spec%LimitPar1MinWeight .lt. 0.0d0) then
                          This%PenaltyLimitPar1 = This%PenaltyLimitPar1 + TmpR
                          This%Penalty          = This%Penalty          + TmpR
                        end if
                      end if
                      nCumMat = Spec%nMat * g
                    end if
                    exit
                  end if
                end do
                ! ... increment i if we have hit the exit, do loop would have ended with i=Spec%nPar1 + 1
                if (i .le. Spec%nPar1) then
                  i = i + 1
                end if
                ! ... the other individuals do not contribute
                if (i .le. Data%nPotPar1) then ! "equal" to capture incremented i + 1 on the do loop exit
                  ! ... zero (the same for all ind so no order)
                  This%Chrom(Rank(i:Data%nPotPar1)) = 0.0d0
                  ! ... negative (the same for all ind so no order)
                  ! This%Chrom(Rank(i:Data%nPotPar1)) = -1.0d0
                  ! ... negative (variable with partially preserving order)
                  !     Found faster convergence than with properly decreasing negative values?
                  !     I guess it adds some more randomness, i.e., it causes suffling of individuals that just did not make it onto the mating list.
                  ! This%Chrom(Rank(i:Data%nPotPar1)) = sign(This%Chrom(Rank(i:Data%nPotPar1)), -1.0d0)
                  ! ... negative and properly decreasing
                  ! TmpR = maxval(This%Chrom(Rank(i:Data%nPotPar1)))
                  ! if (TmpR .gt. 0.0d0) then
                  !     This%Chrom(Rank(i:Data%nPotPar1)) = This%Chrom(Rank(i:Data%nPotPar1)) - abs(TmpR)
                  ! end if
                  ! ... negative (random so no order)
                  ! do j = i, Data%nPotPar1
                  !   call random_number(RanNum)
                  !   This%Chrom(Rank(j)) = -1.0d0 * RanNum
                  ! end do
                end if
                ! ... Spec%nMat still not reached?
                do while (nCumMat .lt. Spec%nMat * g)
                  ! ... add more contributions
                  do i = Spec%nPar1, 1, -1 ! start with the lowest ranked individuals selected as parents (to avoid local optima)
                    j = Rank(i)
                    This%Chrom(j) = This%Chrom(j) + 1.0d0
                    ! ... accumulate and check if we reached Spec%nMat
                    nCumMat = nCumMat + 1
                    if (nCumMat .ge. Spec%nMat * g) then
                      ! To cater for real vs. integer issues
                      TmpI = sum(nint(This%Chrom(Rank(1:Spec%nPar1))))
                      if (TmpI .ne. Spec%nMat * g) then
                        if (TmpI .gt. Spec%nMat * g) then
                          This%Chrom(j) = dble(nint(This%Chrom(j)) - 1)
                        else
                          This%Chrom(j) = dble(nint(This%Chrom(j)) + 1)
                        end if
                      end if
                      exit
                    end if
                  end do
                end do
              end if

              ! "Parent2"
              if (Spec%GenderGiven) then
                ! ... find ranks to find the top values
                if (.not. (Spec%EqualizePar2 .and. (Spec%nPar2 .eq. Data%nPotPar2))) then
                  Rank(1:Data%nPotPar2) = MrgRnk(This%Chrom((Data%nPotPar1 + 1):(Data%nPotPar1 + Data%nPotPar2)))
                  Rank(1:Data%nPotPar2) = Rank(Data%nPotPar2:1:-1) ! MrgRnk ranks small to large
                end if
                ! ... handle cases with equalized contributions
                if (Spec%EqualizePar2) then
                  if (Spec%nPar2 .eq. Data%nPotPar2) then
                    ! ... set integers to all the values (no need for sorting here)
                    This%Chrom((Data%nPotPar1 + 1):(Data%nPotPar1 + Data%nPotPar2)) = dble(Spec%nMat) / Spec%nPar2
                  else
                    ! ... set integers to the top values
                    This%Chrom(Data%nPotPar1 + Rank(1:Spec%nPar2)) = dble(Spec%nMat) / Spec%nPar2
                    ! @todo anything better to preserve the order of non contributing individuals? See below!
                    This%Chrom(Data%nPotPar1 + Rank((Spec%nPar2 + 1):Data%nPotPar2)) = 0.0d0
                    ! This%Chrom(Data%nPotPar1+Rank((Spec%nPar2+1):Data%nPotPar2)) = -1.0d0
                  end if
                else
                  ! ... handle cases with unequal contributions
                  ! ... work for the defined number or parents
                  nCumMat = 0
                  do i = 1, Spec%nPar2
                    j = Data%nPotPar1 + Rank(i)
                    ! ... these would have "zero" contributions and if we hit them then they need at least minimum contrib
                    if (This%Chrom(j) .lt. Spec%LimitPar2Min) then
                      This%Chrom(j) = Spec%LimitPar2Min ! set/fix to minimum usage @todo could consider penalising solution instead?
                    end if
                    ! ... but not above max allowed
                    if (This%Chrom(j) .gt. Spec%LimitPar2Max) then
                      This%Chrom(j) = Spec%LimitPar2Max ! set/fix to maximum usage @todo could consider penalising solution instead?
                    end if
                    ! ... accumulate and check if we reached Spec%nMat
                    nCumMat = nCumMat + nint(This%Chrom(j)) ! internally real, externally integer
                    if (nCumMat .ge. Spec%nMat) then
                      ! ... there should be exactly Spec%nMat contributions
                      if (nCumMat .gt. Spec%nMat) then
                        This%Chrom(j) = This%Chrom(j) - dble(nCumMat - Spec%nMat)
                        ! ... did we go below the minimum usage limit?
                        if (nint(This%Chrom(j)) .lt. Spec%LimitPar2Min) then
                          TmpR = Spec%LimitPar2MinWeight * (Spec%LimitPar2Min - nint(This%Chrom(j)))
                          This%Objective = This%Objective + TmpR
                          if (Spec%LimitPar2MinWeight .lt. 0.0d0) then
                            This%PenaltyLimitPar2 = This%PenaltyLimitPar2 + TmpR
                            This%Penalty          = This%Penalty          + TmpR
                          end if
                        end if
                        nCumMat = Spec%nMat
                      end if
                      exit
                    end if
                  end do
                  ! ... increment i if we have hit the exit, do loop would have ended with i=Spec%nPar2+1
                  if (i .le. Spec%nPar2) then
                    i = i + 1
                  end if
                  ! ... the other individuals do not contribute
                  if (i .le. Data%nPotPar2) then ! equal to capture incremented i+1 on the do loop exit
                    ! ... zero (the same for all ind so no order)
                    This%Chrom(Data%nPotPar1 + Rank(i:Data%nPotPar2)) = 0.0d0
                    ! ... negative (the same for all ind so no order)
                    ! This%Chrom(Data%nPotPar1 + Rank(i:Data%nPotPar2)) = -1.0d0
                    ! ... negative (variable with partially preserving order, i.e., ~large positives become ~large negatives)
                    !     Found faster convergence than with properly decreasing negative values?
                    !     I guess it adds some more randomness, i.e., it causes suffling of individuals that just did not make it onto the mating list.
                    ! This%Chrom(Data%nPotPar1 + Rank(i:Data%nPotPar2)) = sign(This%Chrom(Data%nPotPar1 + Rank(i:Data%nPotPar2)), -1.0d0)
                    ! ... negative and properly decreasing
                    ! TmpR = maxval(This%Chrom(Data%nPotPar1 + Rank(i:Data%nPotPar2)))
                    ! if (TmpR .gt. 0.0d0) then
                    !     This%Chrom(Data%nPotPar1 + Rank(i:Data%nPotPar2)) = This%Chrom(Data%nPotPar1 + Rank(i:Data%nPotPar2)) - abs(TmpR)
                    ! end if
                    ! ... negative (random so no order)
                    ! do j = i, Data%nPotPar2 ! @todo really need this loop?
                    !   call random_number(RanNum)
                    !   This%Chrom(Data%nPotPar1 + Rank(j)) = -1.0d0 * RanNum
                    ! end do
                  end if
                  ! ... Spec%nMat still not reached?
                  do while (nCumMat .lt. Spec%nMat)
                    ! ... add more contributions
                    do i = Spec%nPar2, 1, -1 ! to bottom ranked selected individuals (to avoid local optima)
                      j = Data%nPotPar1 + Rank(i)
                      This%Chrom(j) = This%Chrom(j) + 1.0d0
                      ! ... accumulate and check if we reached Spec%nMat
                      nCumMat = nCumMat + 1
                      if (nCumMat .eq. Spec%nMat) then
                        ! To cater for real vs. integer issues
                        TmpI = sum(nint(This%Chrom(Data%nPotPar1 + Rank(1:Spec%nPar2))))
                        if (TmpI .ne. Spec%nMat) then
                          if (TmpI .gt. Spec%nMat) then
                            This%Chrom(j) = dble(nint(This%Chrom(j)) - 1)
                          else
                            This%Chrom(j) = dble(nint(This%Chrom(j)) + 1)
                          end if
                        end if
                        exit
                      end if
                    end do
                  end do
                end if
              end if

              ! --- Contributions (nVec & xVec) ---

              nVecPar1 = 0

              ! "Parent1"
              ! ... get integer values
              ChromInt(1:Data%nPotPar1) = nint(This%Chrom(1:Data%nPotPar1))
              ! ... remove negatives
              do i = 1, Data%nPotPar1
                if (ChromInt(i) .lt. 0) then
                  ChromInt(i) = 0
                end if
              end do
              ! ... map internal to external order
              nVecPar1 = ChromInt(1:Data%nPotPar1)
              if (.not. Spec%GenderGiven) then
                This%nVec = nVecPar1
              else
                This%nVec(Data%IdPotPar1) = nVecPar1
              end if

              ! "Parent2"
              if (Spec%GenderGiven) then
                nVecPar2 = 0
                ! ... get integer values
                ChromInt(1:Data%nPotPar2) = nint(This%Chrom((Data%nPotPar1 + 1):(Data%nPotPar1 + Data%nPotPar2)))
                ! ... remove negatives
                do i = 1, Data%nPotPar2
                  if (ChromInt(i) .lt. 0) then
                    ChromInt(i) = 0
                  end if
                end do
                ! ... map internal to external order
                nVecPar2 = ChromInt(1:Data%nPotPar2)
                This%nVec(Data%IdPotPar2) = nVecPar2
              end if

              This%xVec = dble(This%nVec) / (2 * Spec%nMat)

              ! --- PAGE ---

              if (Spec%PAGEPar) then
                if (.not. Spec%GenderGiven) then
                  Rank(1:Data%nInd) = MrgRnk(This%Chrom((Data%nPotPar1 + Spec%nMat + 1):(Data%nPotPar1 + Spec%nMat + Data%nInd)))
                  This%GenomeEdit(Rank(Data%nInd:(Data%nInd-Spec%PAGEPar1Max+1):-1)) = 1.0d0 ! MrgRnk ranks small to large
                else
                  if (Spec%PAGEPar1) then
                    Rank(1:Data%nPotPar1) = MrgRnk(This%Chrom((Data%nPotPar1 + Data%nPotPar2 + Spec%nMat + 1):(Data%nPotPar1 + Data%nPotPar2 + Spec%nMat + Data%nPotPar1)))
                    This%GenomeEdit(Data%IdPotPar1(Rank(Data%nPotPar1:(Data%nPotPar1 - Spec%PAGEPar1Max+1):-1))) = 1.0d0 ! MrgRnk ranks small to large
                  end if
                  if (Spec%PAGEPar2) then
                    Rank(1:Data%nPotPar2) = MrgRnk(This%Chrom((Data%nPotPar1 + Data%nPotPar2 + Spec%nMat + Data%nPotPar1 + 1):(Data%nPotPar1 + Data%nPotPar2 + Spec%nMat + Data%nPotPar1 + Data%nPotPar2)))
                    This%GenomeEdit(Data%IdPotPar2(Rank(Data%nPotPar2:(Data%nPotPar2 - Spec%PAGEPar2Max + 1):-1))) = 1.0d0 ! MrgRnk ranks small to large
                  end if
                end if
              end if

              ! --- Mate allocation ---

              MatPar2 = 0
              if (Spec%GenderGiven) then
                ! Distribute parent2 (=female) contributions into matings
                k = 0
                do i = 1, Data%nPotPar2 ! need to loop whole nVecPar2 as some entries are zero
                  do j = 1, nVecPar2(i)
                    k = k + 1
                    MatPar2(k) = Data%IdPotPar2(i)
                  end do
                end do
                ! Reorder parent2 contributions according to the rank of matings
                Rank(1:Spec%nMat) = MrgRnk(This%Chrom((Data%nPotPar1 + Data%nPotPar2 + 1):(Data%nPotPar1 + Data%nPotPar2 + Spec%nMat)))
                MatPar2 = MatPar2(Rank(1:Spec%nMat))
              else
                ! Distribute one half of contributions into matings
                k = 0
                do while (k .lt. Spec%nMat)
                  do i = 1, Data%nPotPar1 ! need to loop whole nVecPar1 as some entries are zero
                    l = nVecPar1(i) / 2
                    if (mod(nVecPar1(i), 2) .eq. 1) then
                      call random_number(RanNum)
                      if (RanNum .gt. 0.5) then
                        l = l + 1
                      end if
                    end if
                    do j = 1, l
                      if (k .eq. Spec%nMat) then
                        exit
                      end if
                      k = k + 1
                      MatPar2(k) = Data%IdPotPar1(i)
                      nVecPar1(i) = nVecPar1(i) - 1
                    end do
                  end do
                end do
                ! Reorder one half of contributions according to the rank of matings
                Rank(1:Spec%nMat) = MrgRnk(This%Chrom((Data%nPotPar1 + 1):(Data%nPotPar1 + Spec%nMat)))
                MatPar2 = MatPar2(Rank(1:Spec%nMat))
              end if

              ! Pair the contributions (=Mating plan)
              k = Spec%nMat ! MrgRnk ranks small to large
              if (Spec%GenderGiven .or. Spec%SelfingAllowed) then
                ! When gender matters selfing can not happen (we have two distinct sets of parents,
                ! unless the user adds individuals of one sex in both sets) and when SelfingAllowed
                ! we do not need to care about it - faster code
                do i = 1, Data%nPotPar1
                  do j = 1, nVecPar1(i)
                    !if (k<2) print*, k, i, j, nVecPar1(i), Spec%nMat, sum(nVecPar1)
                    This%MatingPlan(1, k) = Data%IdPotPar1(i)
                    This%MatingPlan(2, k) = MatPar2(k)
                    k = k - 1
                  end do
                end do
              else
                ! When gender does not matter, selfing can happen (we have one set of parents)
                ! and when selfing is not allowed we need to avoid it - slower code
                do i = 1, Data%nPotPar1
                  do j = 1, nVecPar1(i)
                    This%MatingPlan(1, k) = Data%IdPotPar1(i)
                    if (MatPar2(k) .eq. Data%IdPotPar1(i)) then
                      ! Try to avoid selfing by swapping the MatPar2 and Rank elements
                      do l = k, 1, -1
                        if (MatPar2(l) .ne. Data%IdPotPar1(i)) then
                          MatPar2([k, l]) = MatPar2([l, k])
                          This%Chrom(Data%nPotPar1 + Rank([k, l])) = This%Chrom(Data%nPotPar1 + Rank([l, k]))
                          exit
                        end if
                      end do
                      if (l .lt. 1) then ! Above loop ran out without finding a swap
                        This%Objective = This%Objective + Spec%SelfingWeight
                        if (Spec%SelfingWeight .lt. 0.0d0) then
                          This%PenaltySelfing = This%PenaltySelfing + Spec%SelfingWeight
                          This%Penalty        = This%Penalty        + Spec%SelfingWeight
                        end if
                      end if
                    end if
                    This%MatingPlan(2, k) = MatPar2(k)
                    k = k - 1
                  end do
                end do
              end if

              ! --- Selection criterion ---

              if (Spec%SelCriterionGiven) then
                This%SelIntensity = dot_product(This%xVec, Data%SelIntensity)
                if (Spec%PAGEPar) then
                  This%SelIntensity = This%SelIntensity + dot_product(This%xVec, Data%SelIntensityPAGE * This%GenomeEdit)
                end if
                ! Inlined SelIntensity2SelCriterion
                This%SelCriterion = This%SelIntensity * Data%SelCriterionStat%Sd + Data%SelCriterionStat%Mean
                ! Inlined SelIntensity2MaxCriterionPct
                This%MaxCriterionPct = (This%SelIntensity                      - Spec%ModeMinCoancestrySpec%SelIntensity) / &
                                       (Spec%ModeMaxCriterionSpec%SelIntensity - Spec%ModeMinCoancestrySpec%SelIntensity) * 100.0d0
                if (Spec%ModeSpec%ObjectiveCriterion) then
                  if      (trim(Spec%ModeSpec%Name) .eq. "MaxCriterion") then
                    This%Objective = This%Objective + This%SelIntensity
                  else if (trim(Spec%ModeSpec%Name) .eq. "Opt") then
                    ! This gives SelCriterion objective in the [0, 1] form of the TargetMaxCriterionPct (Kinghorn)
                    This%Objective = This%Objective + This%MaxCriterionPct / Spec%ModeSpec%TargetMaxCriterionPct
                  end if
                end if
                ! Handle beyond the nadir point case so that degree calculation will be meaningful
                if (This%MaxCriterionPct .lt. 0.0d0) then
                  This%MaxCriterionPct = 0.0d0
                end if
              end if

              ! --- Generic individual criterion ---

              if (Spec%GenericIndCritGiven) then
                do j = 1, Spec%nGenericIndCrit
                  TmpR = dot_product(This%xVec, Data%GenericIndCrit(:, j))
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

              ! @todo Enable running AlphaMate without coancestry matrix, i.e.,
              !       select most performant individuals, but consider some limits
              !       on the amount of use etc without coancestry among them.

              ! x'C
              do i = 1, Data%nInd
                TmpVec(i, 1) = dot_product(This%xVec, Data%Coancestry%Value(1:, i))
              end do
              ! @todo consider using matmul instead of repeated dot_product?
              ! @todo consider using BLAS/LAPACK - perhaps non-symmetric is more optimised?
              ! Matrix multiplication with a symmetric matrix using BLAS routine
              ! (it was ~5x slower than the above with 1.000 individuals, might be benefical with larger cases)
              ! http://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3.html#ga253c8edb8b21d1b5b1783725c2a6b692
              ! call dsymm(side="l", uplo="l", m=Data%nInd, n=1, alpha=1.0d0, A=CoaMtx, lda=Data%nInd, b=This%xVec, ldb=Data%nInd, beta=0, c=TmpVec, ldc=Data%nInd)
              ! call dsymm(     "l",      "l",   Data%nInd,   1,       1.0d0,   CoaMtx,     Data%nInd,   This%xVec,     Data%nInd,      0,   TmpVec,     Data%nInd)

              ! x'Cx
              This%CoancestryRanMate = dot_product(TmpVec(:, 1), This%xVec)

              ! Inlined Coancestry2CoancestryRate
              This%CoancestryRateRanMate = (This%CoancestryRanMate - Data%CoancestryRanMate) / &
                                                            (1.0d0 - Data%CoancestryRanMate)
              ! Inlined CoancestryRate2MinCoancestryPct
              This%MinCoancestryPct = (Spec%ModeMaxCriterionSpec%CoancestryRate - This%CoancestryRateRanMate) / &
                                      (Spec%ModeMaxCriterionSpec%CoancestryRate - Spec%ModeMinCoancestrySpec%CoancestryRate) * 100.d0
              ! Handle beyond the nadir point case so that degree calculation will be meaningful
              if (This%MinCoancestryPct .lt. 0.0d0) then
                This%MinCoancestryPct = 0.0d0
              end if

              ! Degree
              ! This calculation ASSUMES unit circular shape of selection/coancestry frontier
              ! - given results from ModeMinCoancestry and ModeMaxCriterion we can map the solution on the unit circle,
              !   i.e., the circle center is at (MaxCoancestryRate, MinSelIntensity) and solution is
              !   at (SolCoancestryRate, SolSelIntensity), which is mapped to (MinCoancestryPct, 1 - MaxCriterionPct)
              !   on unit circle with center at (1, 0) (MinCoancestryMode is at (0, 0) and MaxCriterionMode is at (1, 1))
              ! - then we can calculate degrees of the angle between the max-line (MaxCoancestryRate, MinSelIntensity)-(MaxCoancestryRate, MaxSelIntensity)
              !   and the sol-line (MaxCoancestryRate, MinSelIntensity)-(SolCoancestryRate, SolSelIntensity)
              ! - the min-line is (MaxCoancestryRate, MinSelIntensity)-(MinCoancestryRate, MinSelIntensity)
              ! - the solution y-coordinate on the max-line is adjacent to the angle and
              !   the solution x-coordinate on the min-line is opposite to the angle (if we
              !   put the min-line parallely up) so we use the arctangent function to compute
              !   angle degrees
              This%Degree = atan(This%MinCoancestryPct / This%MaxCriterionPct) * RAD2DEG

              !@todo ModeRan?

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
                  !     TmpR = 0.0d0
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
                  ! TmpR = 1.0d0 - This%MinCoancestryPct / Spec%ModeSpec%TargetMinCoancestryPct
                  ! if (This%MinCoancestryPct .lt. Spec%ModeSpec%TargetMinCoancestryPct) then
                  !   ! CANCEL SelCriterion from the objective (want to focus on coancestry for such solutions)
                  !   This%Objective = This%Objective - This%MaxCriterionPct / Spec%ModeSpec%TargetMaxCriterionPct
                  ! else ! if (This%MinCoancestryPct .ge. Spec%ModeSpec%TargetMinCoancestryPct) then
                  !   if (Spec%ModeSpec%CoancestryWeightBelow) then
                  !     TmpR = abs(TmpR)
                  !   else
                  !     TmpR = 0.0d0
                  !   end if
                  ! end if

                  ! Weight * (1 - Degree / TargetDegree) [negative = good for coancestry, positive = bad for criterion]
                  !   * larger Degree is better in terms of coancestry, but smaller degree is better in terms of criterion
                  !     1 - c(50, 60, 70)/70 =  0.28,  0.14, 0
                  !     1 - c(90, 80, 70)/70 = -0.28, -0.14, 0
                  !   * this gives coancestry objective in the [-Inf, Inf] form of the TargetDegree
                  !   * this penalty formulation works ~well with weights close to -1 and works
                  !     well on any part of the frontier (due to exploiting its geometry, albeit assuming circular shape)
                  TmpR = 1.0d0 - This%Degree / Spec%ModeSpec%TargetDegree
                  if (This%Degree .lt. Spec%ModeSpec%TargetDegree) then
                    ! CANCEL SelCriterion from the objective (want to focus on coancestry for such solutions)
                    This%Objective = This%Objective - This%MaxCriterionPct / Spec%ModeSpec%TargetMaxCriterionPct
                  else ! if (This%Degree .ge. Spec%ModeSpec%TargetDegree) then
                    if (Spec%ModeSpec%CoancestryWeightBelow) then
                      TmpR = abs(TmpR)
                    else
                      TmpR = 0.0d0
                    end if
                  end if

                  TmpR = Spec%CoancestryWeight * TmpR
                  This%Objective = This%Objective + TmpR
                  if (Spec%CoancestryWeight .lt. 0.0d0) then
                    This%PenaltyCoancestry = This%PenaltyCoancestry + TmpR
                    This%Penalty           = This%Penalty           + TmpR
                  end if
                end if
              end if

              ! --- Progeny inbreeding (=inbreeding of a mating) ---

              TmpR = 0.0d0
              do j = 1, Spec%nMat
                ! Lower triangle to speedup lookup
                TmpMax = maxval(This%MatingPlan(:, j))
                TmpMin = minval(This%MatingPlan(:, j))
                TmpR = TmpR + Data%Coancestry%Value(TmpMax, TmpMin)
              end do
              ! @Todo: different number of progeny per mating???
              This%Inbreeding = TmpR / Spec%nMat
              ! Inlined Coancestry2CoancestryRate
              This%InbreedingRate = (This%Inbreeding - Data%Inbreeding) / (1.0d0 - Data%Inbreeding)
              ! Inlined CoancestryRate2MinCoancestryPct
              This%MinInbreedingPct = (+1.0d0 - This%InbreedingRate) / &
                                      (+1.0d0 - Spec%ModeMinInbreedingSpec%InbreedingRate) * 100.d0

              ! @todo: Pareto front formulation for this objective too?
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
                  !     TmpR = 0.0d0
                  !   end if
                  ! end if

                  ! Weight * (1 - MinInbreedingPct / TargetMinInbreedingPct) [negative = good, positive = bad]
                  !   * bigger MinInbreedingPct is better
                  !     1 - c(50, 60, 70)/70 =  0.28,  0.14, 0
                  !     1 - c(90, 80, 70)/70 = -0.28, -0.14, 0
                  !     1 -            0 /70 =  1
                  !   * this gives inbreeding objective in the [-Inf, Inf] form of the TargetMinInbreedingPct
                  !   * this penalty formulation works ~well with weights close to -1
                  TmpR = 1.0d0 - This%MinInbreedingPct / Spec%ModeSpec%TargetMinInbreedingPct
                  if (This%MinInbreedingPct .lt. Spec%ModeSpec%TargetMinInbreedingPct) then
                    ! CANCEL SelCriterion from the objective (want to focus on inbreeding for such solutions)
                    This%Objective = This%Objective - This%MaxCriterionPct / Spec%ModeSpec%TargetMaxCriterionPct
                  else ! if (This%MinInbreedingPct .ge. Spec%ModeSpec%TargetMinInbreedingPct) then
                    if (Spec%ModeSpec%InbreedingWeightBelow) then
                      TmpR = abs(TmpR)
                    else
                      TmpR = 0.0d0
                    end if
                  end if

                  TmpR = Spec%InbreedingWeight * TmpR
                  This%Objective = This%Objective + TmpR
                  if (Spec%InbreedingWeight .lt. 0.0d0) then
                    This%PenaltyInbreeding = This%PenaltyInbreeding + TmpR
                    This%Penalty           = This%Penalty           + TmpR
                  end if
                end if
              end if

              ! --- Generic mating criterion ---

              if (Spec%GenericMatCritGiven) then
                do k = 1, Spec%nGenericMatCrit
                  TmpR = 0.0d0
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

              ! @todo how should we handle costs?
          end select
      end select
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Call various optimisations for AlphaMate
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !---------------------------------------------------------------------------
    ! @todo: the best random mating?
    ! @todo: return solutions?
    subroutine AlphaMateSearch(Spec, Data, LogStdout) ! not pure due to IO
      implicit none
      type(AlphaMateSpec), intent(inout) :: Spec      !< AlphaMateSpec holder (out because we set and reset some parameters for different search modes)
      type(AlphaMateData), intent(inout) :: Data      !< AlphaMateData holder (out because we set and reset some parameters for different search modes)
      logical, intent(in), optional      :: LogStdout !< Log process on stdout (default .false.)

      type(AlphaMateSol) :: SolMinCoancestry, SolMinInbreeding, SolMaxCriterion, Sol !< For frontier modes and random mating (no optimisation) mode
      type(AlphaMateSol), allocatable :: SolOpt(:) !< Optimal solutions

      integer(int32) :: nParam, Point, Target, FrontierUnit

      real(real64), allocatable :: InitEqual(:, :)

      logical :: LogStdoutInternal !, OptimOK

      character(len=FILELENGTH) :: LogFile, LogPopFile, ContribFile, MatingFile

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

      ! @todo if (Spec%MateAllocation) then
      nParam = nParam + Spec%nMat
      ! @todo end if

      if (Spec%PAGEPar) then
        nParam = nParam + Data%nInd
      end if

      ! --- Minimum future coancestry ---

      if (Spec%ModeMinCoancestry) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") " Optimise contributions for minimum future coancestry (ModeMinCoancestry) ..."
          write(STDOUT, "(a)") " "
        end if

        ! Setup
        call Spec%SetupMode(Mode="MinCoancestry")
        call SolMinCoancestry%SetupColNamesAndFormats(Spec=Spec)

        LogFile     = "OptimisationLogModeMinCoancestry.txt"
        LogPopFile  = "OptimisationLogPopModeMinCoancestry.txt"
        ContribFile = "ContributionsModeMinCoancestry.txt"
        MatingFile  = "MatingPlanModeMinCoancestry.txt"

        ! Search
        ! @todo Can we do this in a better way where we take "structure" of Chrom into account?
        allocate(InitEqual(nParam, nint(Spec%EvolAlgNSol * 0.1)))
        InitEqual = 1.0d0 ! A couple of solutions that would give equal contributions to everybody

        if (trim(Spec%EvolAlg) .eq. "DE") then
          call DifferentialEvolution(Spec=Spec, Data=Data, nParam=nParam, nSol=Spec%EvolAlgNSol, Init=InitEqual, &
            nIter=Spec%EvolAlgNIter, nIterBurnIn=Spec%DiffEvolNIterBurnIn, nIterStop=Spec%EvolAlgNIterStop, StopTolerance=Spec%EvolAlgStopTol, nIterPrint=Spec%EvolAlgNIterPrint, &
            LogStdout=LogStdoutInternal, LogFile=LogFile, LogPop=Spec%EvolAlgLogPop, LogPopFile=LogPopFile, &
            CRBurnIn=Spec%DiffEvolParamCrBurnIn, CRLate=Spec%DiffEvolParamCr, FBase=Spec%DiffEvolParamFBase, FHigh1=Spec%DiffEvolParamFHigh1, FHigh2=Spec%DiffEvolParamFHigh2, &
            BestSol=SolMinCoancestry)!, Status=OptimOK)
          ! if (.not. OptimOK) then
          !   write(STDERR, "(a)") " ERROR: Optimisation failed!"
          !   write(STDERR, "(a)") " "
          !   stop 1
          ! end if
        end if

        deallocate(InitEqual)

        ! Save
        !@todo Do these lines still make sense when we go above two objectives?
        SolMinCoancestry%Degree           =  90.0d0
        SolMinCoancestry%MinCoancestryPct = 100.0d0
        SolMinCoancestry%MinInbreedingPct =   0.0d0
        SolMinCoancestry%MaxCriterionPct  =   0.0d0
        call Spec%ModeMinCoancestrySpec%SaveSol2ModeSpec(In=SolMinCoancestry)
        call SolMinCoancestry%WriteContributions(Data, ContribFile)
        call SolMinCoancestry%WriteMatingPlan(Data, MatingFile)
      end if

      ! --- Minimum future inbreeding ---

      if (Spec%ModeMinInbreeding) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") " Optimise contributions for minimum future inbreeding (ModeMinInbreeding) ..."
          write(STDOUT, "(a)") " "
        end if

        ! Setup
        call Spec%SetupMode(Mode="MinInbreeding")
        call SolMinInbreeding%SetupColNamesAndFormats(Spec=Spec)

        LogFile     = "OptimisationLogModeMinInbreeding.txt"
        LogPopFile  = "OptimisationLogPopModeMinInbreeding.txt"
        ContribFile = "ContributionsModeMinInbreeding.txt"
        MatingFile  = "MatingPlanModeMinInbreeding.txt"

        ! Search
        ! @todo Can we do this in a better way where we take "structure" of Chrom into account?
        allocate(InitEqual(nParam, nint(Spec%EvolAlgNSol * 0.1)))
        InitEqual = 1.0d0 ! A couple of solutions that would give equal contributions to everybody

        if (trim(Spec%EvolAlg) .eq. "DE") then
          call DifferentialEvolution(Spec=Spec, Data=Data, nParam=nParam, nSol=Spec%EvolAlgNSol, Init=InitEqual, &
            nIter=Spec%EvolAlgNIter, nIterBurnIn=Spec%DiffEvolNIterBurnIn, nIterStop=Spec%EvolAlgNIterStop, StopTolerance=Spec%EvolAlgStopTol, nIterPrint=Spec%EvolAlgNIterPrint, &
            LogStdout=LogStdoutInternal, LogFile=LogFile, LogPop=Spec%EvolAlgLogPop, LogPopFile=LogPopFile, &
            CRBurnIn=Spec%DiffEvolParamCrBurnIn, CRLate=Spec%DiffEvolParamCr, FBase=Spec%DiffEvolParamFBase, FHigh1=Spec%DiffEvolParamFHigh1, FHigh2=Spec%DiffEvolParamFHigh2, &
            BestSol=SolMinInbreeding)!, Status=OptimOK)
          ! if (.not. OptimOK) then
          !   write(STDERR, "(a)") " ERROR: Optimisation failed!"
          !   write(STDERR, "(a)") " "
          !   stop 1
          ! end if
        end if

        deallocate(InitEqual)

        ! Save
        !@todo Do these lines still make sense when we go above two objectives?
        SolMinInbreeding%Degree           =  45.0d0
        SolMinInbreeding%MinCoancestryPct =   0.0d0
        SolMinInbreeding%MinInbreedingPct = 100.0d0
        SolMinInbreeding%MaxCriterionPct  =   0.0d0
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
        call SolMinInbreeding%WriteContributions(Data, ContribFile)
        call SolMinInbreeding%WriteMatingPlan(Data, MatingFile)
      end if

      ! --- Maximum future selection criterion ---

      if (Spec%ModeMaxCriterion) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") " Optimise contributions for maximum future selection criterion (ModeMaxCriterion) ..."
          write(STDOUT, "(a)") " "
        end if

        ! Setup
        call Spec%SetupMode(Mode="MaxCriterion")
        call SolMaxCriterion%SetupColNamesAndFormats(Spec=Spec)

        LogFile     = "OptimisationLogModeMaxCriterion.txt"
        LogPopFile  = "OptimisationLogPopModeMaxCriterion.txt"
        ContribFile = "ContributionsModeMaxCriterion.txt"
        MatingFile  = "MatingPlanModeMaxCriterion.txt"

        ! Search
        ! @todo add some clever initial values, say:
        !       - equal contributions for top 2/3 or 1/2 of BV distribution,
        !       - decreasing contributions with decreasing value
        !       - SDP solution, ...?
        if (trim(Spec%EvolAlg) .eq. "DE") then
          call DifferentialEvolution(Spec=Spec, Data=Data, nParam=nParam, nSol=Spec%EvolAlgNSol, nIter=Spec%EvolAlgNIter, nIterBurnIn=Spec%DiffEvolNIterBurnIn, &
            nIterStop=Spec%EvolAlgNIterStop, StopTolerance=Spec%EvolAlgStopTol, nIterPrint=Spec%EvolAlgNIterPrint, &
            LogStdout=LogStdoutInternal, LogFile=LogFile, LogPop=Spec%EvolAlgLogPop, LogPopFile=LogPopFile, &
            CRBurnIn=Spec%DiffEvolParamCrBurnIn, CRLate=Spec%DiffEvolParamCr, FBase=Spec%DiffEvolParamFBase, FHigh1=Spec%DiffEvolParamFHigh1, FHigh2=Spec%DiffEvolParamFHigh2, &
            BestSol=SolMaxCriterion)!, Status=OptimOK)
          ! if (.not. OptimOK) then
          !   write(STDERR, "(a)") " ERROR: Optimisation failed!"
          !   write(STDERR, "(a)") " "
          !   stop 1
          ! end if
        end if

        ! Save
        !@todo Do these lines still make sense when we go above two objectives?
        SolMaxCriterion%Degree           =   0.0d0
        SolMaxCriterion%MinCoancestryPct =   0.0d0
        SolMinInbreeding%MinInbreedingPct =  0.0d0
        SolMaxCriterion%MaxCriterionPct  = 100.0d0
        call Spec%ModeMaxCriterionSpec%SaveSol2ModeSpec(In=SolMaxCriterion)
        call SolMaxCriterion%WriteContributions(Data, ContribFile)
        call SolMaxCriterion%WriteMatingPlan(Data, MatingFile)
      end if

      ! --- Evaluate frontier ---
      ! ModeMinCoancestry, ModeMinInbreeding, and ModeMaxCriterion must be run prior to this!
      ! @todo This code will have to change in light of more than two objectives using the NBI or AWS method

      if (Spec%EvaluateFrontier) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") " Evaluate frontier ..."
        end if

        open(newunit=FrontierUnit, file="Frontier.txt", status="unknown")

        ! Setup
        call Sol%SetupColNamesAndFormats(Spec=Spec) ! need to do this here so that we can log previous results
        call Sol%LogHead(LogUnit=FrontierUnit, String="ModeOrPoint", StringNum=18)

        ! Add minimum coancestry solution to frontier (90 degress with two objectives)
        call SolMinCoancestry%Log(FrontierUnit, Iteration=-1, AcceptPct=-1.0d0, String="ModeMinCoancestry", StringNum=18)

        ! Add minimum inbreeding solution to frontier
        if (Spec%ModeMinInbreeding) then
          call SolMinInbreeding%Log(FrontierUnit, Iteration=-1, AcceptPct=-1.0d0, String="ModeMinInbreeding", StringNum=18)
        end if

        ! Frontier
        do Point = 1, size(TARGETDEGREEFRONTIER) ! 80, 70, ..., 10 degrees

          ! Setup
          call Spec%SetupMode(Mode="Opt", Data=Data, &
                              Degree=TARGETDEGREEFRONTIER(Point), &
                              ModeMinCoancestrySpec=Spec%ModeMinCoancestrySpec, CoancestryWeightBelow=.true., &
                              ModeMinInbreedingSpec=Spec%ModeMinInbreedingSpec, &
                              ModeMaxCriterionSpec=Spec%ModeMaxCriterionSpec)
                              ! CoancestryWeightBelow=.true. as we want to get ~exact degrees
          if (LogStdoutInternal) then
            write(STDOUT, "(a)") " "
            write(STDOUT, "(a)") "  Point "//Int2Char(Point)//" of "//Int2Char(size(TARGETDEGREEFRONTIER))
            write(STDOUT, "(a)") " "
            call Spec%ModeSpec%LogTargets(Unit=STDOUT, Spec=Spec)
          end if
          call Sol%SetupColNamesAndFormats(Spec=Spec) ! call again as InitialiseAlphaMateSol and AssignAlphaMateSol "nullify" the  above SetupColNamesAndFormats call (ugly, but works ...)

          LogFile     = "OptimisationLogModeFrontier"//trim(Int2Char(Point))//".txt"
          LogPopFile  = "OptimisationLogPopModeFrontier"//trim(Int2Char(Point))//".txt"
          ContribFile = "ContributionsModeFrontier"//trim(Int2Char(Point))//".txt"
          MatingFile  = "MatingPlanModeFrontier"//trim(Int2Char(Point))//".txt"

          ! Search
          if (trim(Spec%EvolAlg) .eq. "DE") then
            call DifferentialEvolution(Spec=Spec, Data=Data, nParam=nParam, nSol=Spec%EvolAlgNSol, nIter=Spec%EvolAlgNIter, nIterBurnIn=Spec%DiffEvolNIterBurnIn, &
              nIterStop=Spec%EvolAlgNIterStop, StopTolerance=Spec%EvolAlgStopTol, nIterPrint=Spec%EvolAlgNIterPrint, &
              LogStdout=LogStdoutInternal, LogFile=LogFile, LogPop=Spec%EvolAlgLogPop, LogPopFile=LogPopFile, &
              CRBurnIn=Spec%DiffEvolParamCrBurnIn, CRLate=Spec%DiffEvolParamCr, FBase=Spec%DiffEvolParamFBase, FHigh1=Spec%DiffEvolParamFHigh1, FHigh2=Spec%DiffEvolParamFHigh2, &
              BestSol=Sol)!, Status=OptimOK)
            ! if (.not. OptimOK) then
            !   write(STDERR, "(a)") " ERROR: Optimisation failed!"
            !   write(STDERR, "(a)") " "
            !   stop 1
            ! end if
          end if

          ! Save
          call Sol%Log(FrontierUnit, Iteration=-1, AcceptPct=-1.0d0, String=trim("ModeFrontier"//trim(Int2Char(Point))), StringNum=18)
          call Sol%WriteContributions(Data, ContribFile)
          call Sol%WriteMatingPlan(Data, MatingFile)

        end do

        ! Add maximum criterion solution to frontier (0 degress with two objectives)
        call SolMaxCriterion%Log(FrontierUnit, Iteration=-1, AcceptPct=-1.0d0, String="ModeMaxCriterion", StringNum=18)

        close(FrontierUnit)

      end if

      ! --- Random mating ---

      ! @todo: is this mode any good?
      if (Spec%ModeRan) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") " Evaluate random mating (ModeRan) ..."
          write(STDOUT, "(a)") " "
        end if

        ! Setup
        call Spec%SetupMode(Mode="Ran")
        call Sol%SetupColNamesAndFormats(Spec=Spec)

        LogFile = "OptimisationLogModeRan.txt"
        ! @todo: other reports from here - at least the best random solution?

        ! Search
        ! @todo Can we do this in a better way where we take "structure" of Chrom into account?
        allocate(InitEqual(nParam, nint(Spec%EvolAlgNSol * 0.1)))
        InitEqual = 1.0d0 ! A couple of solutions that would give equal contributions for everybody

        call RandomSearch(Mode="avg", Spec=Spec, Data=Data, nParam=nParam, Init=InitEqual, &
          nSamp=Spec%EvolAlgNSol*Spec%EvolAlgNIter*Spec%RanAlgStricter, nSampStop=Spec%EvolAlgNIterStop*Spec%RanAlgStricter, &
          StopTolerance=Spec%EvolAlgStopTol/Spec%RanAlgStricter, nSampPrint=Spec%EvolAlgNIterPrint, &
          LogStdout=LogStdoutInternal, LogFile=LogFile, BestSol=Sol)

        deallocate(InitEqual)
      end if

      ! --- Maximum future selection criterion with constraint on coancestry/inbreeding ---
      ! ModeMinCoancestry, ModeMinInbreeding, and ModeMaxCriterion should be run prior to this!

      ! if (Spec%ModeOpt) then
      !   if (LogStdoutInternal) then
      !     write(STDOUT, "(a)") " "
      !     write(STDOUT, "(a)") " Optimise contributions for maximum future selection criterion with constraint on coancestry/inbreeding (ModeOpt) ..."
      !     write(STDOUT, "(a)") " "
      !   end if

      !   ! Setup
      !   allocate(SolOpt(Spec%nTargets))
      !   do Target = 1, Spec%nTargets
      !     call Spec%SetupMode(Mode="Opt", Data=Data, ModeMinCoancestrySpec=Spec%ModeMinCoancestrySpec, ModeMaxCriterionSpec=Spec%ModeMaxCriterionSpec)
      !     if (LogStdoutInternal) then
      !       write(STDOUT, "(a)") "  Target "//Int2Char(Target)//" of "//Int2Char(Spec%nTargets)
      !       write(STDOUT, "(a)") " "
      !       call Spec%ModeSpec%LogTargets(Unit=STDOUT, Spec=Spec)
      !     end if
      !     call SolOpt(Target)%SetupColNamesAndFormats(Spec=Spec)

      !     LogFile     = "OptimisationLogModeOptTarget"//Int2Char(Target)//".txt"
      !     LogPopFile  = "OptimisationLogPopModeOptTarget"//Int2Char(Target)//".txt"
      !     ContribFile = "ContributionsModeOptTarget"//Int2Char(Target)//".txt"
      !     MatingFile  = "MatingPlanModeOptTarget"//Int2Char(Target)//".txt"

      !     ! Search
      !     ! @todo add some clever initial values, say:
      !     !       - equal contributions for top 2/3 or 1/2 of BV distribution,
      !     !       - decreasing contributions with decreasing value
      !     !       - SDP solution, ...?
      !     if (trim(Spec%EvolAlg) .eq. "DE") then
      !       call DifferentialEvolution(Spec=Spec, Data=Data, nParam=nParam, nSol=Spec%EvolAlgNSol, nIter=Spec%EvolAlgNIter, nIterBurnIn=Spec%DiffEvolNIterBurnIn, &
      !         nIterStop=Spec%EvolAlgNIterStop, StopTolerance=Spec%EvolAlgStopTol, nIterPrint=Spec%EvolAlgNIterPrint, &
      !         LogStdout=LogStdoutInternal, LogFile=LogFile, LogPop=Spec%EvolAlgLogPop, LogPopFile=LogPopFile, &
      !         CRBurnIn=Spec%DiffEvolParamCrBurnIn, CRLate=Spec%DiffEvolParamCr, FBase=Spec%DiffEvolParamFBase, FHigh1=Spec%DiffEvolParamFHigh1, FHigh2=Spec%DiffEvolParamFHigh2, &
      !         BestSol=SolOpt(Target))!, Status=OptimOK)
      !       ! if (.not. OptimOK) then
      !       !   write(STDERR, "(a)") " ERROR: Optimisation failed!"
      !       !   write(STDERR, "(a)") " "
      !       !   stop 1
      !       ! end if
      !     end if

      !     ! Save
      !     call SolOpt(Target)%WriteContributions(Data, ContribFile)
      !     call SolOpt(Target)%WriteMatingPlan(Data, MatingFile)
      !   end do
      ! end if

        ! if (Spec%ModeOpt) then
        !   call SolOpt%Log   (FrontierUnit, Iteration=-1, AcceptPct=-1.0d0, String="ModeOpt",    StringNum=15)
        ! end if
        ! if (Spec%ModeRan) then
        !   call SolRan%Log   (FrontierUnit, Iteration=-1, AcceptPct=-1.0d0, String="ModeRan",    StringNum=15)
        ! end if

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

      ! Solution components

      write(Unit, *) "Objective: ", This%Objective
      write(Unit, *) "nParam: ", This%nParam
      write(Unit, *) "Chrom: ", This%Chrom
      write(Unit, *) "Penalty: ", This%Penalty
      write(Unit, *) "PenaltyCoancestry: ", This%PenaltyCoancestry
      write(Unit, *) "PenaltyInbreeding: ", This%PenaltyInbreeding
      write(Unit, *) "PenaltySelfing: ", This%PenaltySelfing
      write(Unit, *) "PenaltyLimitPar1: ", This%PenaltyLimitPar1
      write(Unit, *) "PenaltyLimitPar2: ", This%PenaltyLimitPar2
      write(Unit, *) "PenaltyGenericIndCrit: ", This%PenaltyGenericIndCrit
      write(Unit, *) "PenaltyGenericMatCrit: ", This%PenaltyGenericMatCrit
      write(Unit, *) "SelCriterion: ", This%SelCriterion
      write(Unit, *) "SelIntensity: ", This%SelIntensity
      write(Unit, *) "CoancestryRanMate: ", This%CoancestryRanMate
      write(Unit, *) "CoancestryRateRanMate: ", This%CoancestryRateRanMate
      write(Unit, *) "Inbreeding: ", This%Inbreeding
      write(Unit, *) "InbreedingRate: ", This%InbreedingRate
      write(Unit, *) "MinInbreedingPct: ", This%MinInbreedingPct
      write(Unit, *) "GenericIndCrit: ", This%GenericIndCrit
      write(Unit, *) "GenericMatCrit: ", This%GenericMatCrit
      write(Unit, *) "Cost: ", This%Cost
      write(Unit, *) "nVec: ", This%nVec
      write(Unit, *) "xVec: ", This%xVec
      write(Unit, *) "Mating plan:"
      do Mat = 1, size(This%MatingPlan, dim=2)
        write(Unit, *) Mat, This%MatingPlan(:, Mat)
      end do
      write(Unit, *) "GenomeEdit: ", This%GenomeEdit

      ! Formats for logging

      write(Unit, *) "FmtLogStdoutHead: ",        trim(This%FmtLogStdoutHead)
      write(Unit, *) "FmtLogStdout: ",            trim(This%FmtLogStdout)
      if (allocated(This%ColnameLogStdout)) then
        write(Unit, *) "ColnameLogStdout: ",        This%ColnameLogStdout
      else
        write(Unit, *) "ColnameLogStdout: not allocated"
      end if
      write(Unit, *) "FmtLogUnitHead: ",          trim(This%FmtLogUnitHead)
      if (allocated(This%ColnameLogUnit)) then
        write(Unit, *) "ColnameLogUnit: ",          This%ColnameLogUnit
      else
        write(Unit, *) "ColnameLogUnit: not allocated"
      end if
      write(Unit, *) "FmtLogPopUnit: ",           trim(This%FmtLogPopUnit)
      if (allocated(This%ColnameLogPopUnit)) then
        write(Unit, *) "ColnameLogPopUnit: ",       This%ColnameLogPopUnit
      else
        write(Unit, *) "ColnameLogPopUnit: not allocated"
      end if
      write(Unit, *) "FmtContributionHead: ",     trim(This%FmtContributionHead)
      write(Unit, *) "FmtContributionHeadEdit: ", trim(This%FmtContributionHeadEdit)
      write(Unit, *) "FmtContribution: ",         trim(This%FmtContribution)
      write(Unit, *) "FmtContributionEdit: ",     trim(This%FmtContributionEdit)
      write(Unit, *) "FmtMatingHead: ",           trim(This%FmtMatingHead)
      write(Unit, *) "FmtMating: ",               trim(This%FmtMating)

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
    subroutine WriteContributions(This, Data, ContribFile) ! not pure due to IO
      implicit none

      ! Arguments
      class(AlphaMateSol), intent(in)        :: This        !< AlphaMateSol holder
      type(AlphaMateData), intent(in)        :: Data        !< AlphaMateData holder
      character(len=*), intent(in), optional :: ContribFile !< File to write individual contributions to (default STDOUT)

      ! Other
      integer(int32) :: IndRev, Ind, ContribUnit, Rank(Data%nInd)

      if (.not. present(ContribFile)) then
        ContribUnit = STDOUT
      end if

      ! @todo should we have constant output no matter which options are switched on?
      if (present(ContribFile)) then
        open(newunit=ContribUnit, file=ContribFile, status="unknown")
      end if

      Rank = MrgRnk(This%nVec)
      if (.not. allocated(This%GenomeEdit)) then
        !                                             1234567890123456789012
        write(ContribUnit, This%FmtContributionHead) "             Id", &
                                                     "         Gender", &
                                                     "   SelCriterion", &
                                                     "  AvgCoancestry", &
                                                     "   Contribution", &
                                                     "        nMating"
        do IndRev = Data%nInd, 1, -1 ! MrgRnk ranks small to large
          Ind = Rank(IndRev)
          write(ContribUnit, This%FmtContribution) Data%Coancestry%OriginalId(Ind), Data%Gender(Ind), Data%SelCriterion(Ind), &
                                                   mean(Data%Coancestry%Value(1:, Ind)), This%xVec(Ind), This%nVec(Ind)
        end do
      else
        !                                                 12345678901234567890123456789012
        write(ContribUnit, This%FmtContributionHeadEdit) "                             Id", &
                                                         "         Gender", &
                                                         "   SelCriterion", &
                                                         "  AvgCoancestry", &
                                                         "   Contribution", &
                                                         "        nMating", &
                                                         "     GenomeEdit", &
                                                         "  EditedSelCrit"
        do IndRev = Data%nInd, 1, -1 ! MrgRnk ranks small to large
          Ind = Rank(IndRev)
          write(ContribUnit, This%FmtContributionEdit) Data%Coancestry%OriginalId(Ind), Data%Gender(Ind), Data%SelCriterion(Ind), &
                                                       mean(Data%Coancestry%Value(1:, Ind)), This%xVec(Ind), This%nVec(Ind), &
                                                       nint(This%GenomeEdit(Ind)), Data%SelCriterion(Ind) + This%GenomeEdit(Ind) * Data%SelCriterionPAGE(Ind)
        end do
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
    subroutine WriteMatingPlan(This, Data, MatingFile) ! not pure due to IO
      implicit none

      ! Arguments
      class(AlphaMateSol), intent(in)        :: This       !< AlphaMateSol holder
      type(AlphaMateData), intent(in)        :: Data       !< AlphaMateData holder
      character(len=*), intent(in), optional :: MatingFile !< File to write mating plan to (default STDOUT)

      ! Other
      integer(int32) :: Mat, MatingUnit

      if (.not. present(MatingFile)) then
        MatingUnit = STDOUT
      end if

      if (present(MatingFile)) then
        open(newunit=MatingUnit, file=MatingFile, status="unknown")
      end if

      !                                      12345678901234567890123456789012
      write(MatingUnit, This%FmtMatingHead) "         Mating", &
                                            "                         Parent1", &
                                            "                         Parent2"
      do Mat = 1, size(This%MatingPlan, dim=2)
        write(MatingUnit, This%FmtMating) Mat, Data%Coancestry%OriginalId(This%MatingPlan(1, Mat)), Data%Coancestry%OriginalId(This%MatingPlan(2, Mat))
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
    pure subroutine SetupColNamesAndFormats(This, Spec)
      implicit none
      class(AlphaMateSol), intent(inout) :: This !< @return AlphaMateSol holder
      type(AlphaMateSpec), intent(in)    :: Spec !< AlphaMateSpec

      integer(int32) :: nCol, nColTmp, i
      character(len=CHARLENGTH) :: Tmp

      ! --- Optimisation log ---

      nCol = 14
      if (Spec%GenericIndCritGiven) then
        nCol = nCol + Spec%nGenericIndCrit
      end if
      if (Spec%GenericMatCritGiven) then
        nCol = nCol + Spec%nGenericMatCrit
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

      !                          1234567890123456789012
      This%ColnameLogUnit(1)  = "             Iteration"
      This%ColnameLogUnit(2)  = "             AcceptPct"
      This%ColnameLogUnit(3)  = "             Objective"
      This%ColnameLogUnit(4)  = "             Penalties"
      This%ColnameLogUnit(5)  = "        FrontierDegree"
      This%ColnameLogUnit(6)  = "          SelCriterion"
      This%ColnameLogUnit(7)  = "          SelIntensity"
      This%ColnameLogUnit(8)  = "       MaxSelCriterPct" ! @todo: Can we extend this to MaxSelCriterPct?
      This%ColnameLogUnit(9)  = "            Coancestry"
      This%ColnameLogUnit(10) = "        CoancestryRate"
      This%ColnameLogUnit(11) = "      MinCoancestryPct"
      This%ColnameLogUnit(12) = "            Inbreeding"
      This%ColnameLogUnit(13) = "        InbreedingRate"
      This%ColnameLogUnit(14) = "      MinInbreedingPct"

      This%ColnameLogStdout(1)  =         "  Iteration"
      This%ColnameLogStdout(2)  =         "  AcceptPct"
      This%ColnameLogStdout(3)  =         "  Objective"
      This%ColnameLogStdout(4)  =         "  Penalties"
      This%ColnameLogStdout(5)  =         "     Degree"
      This%ColnameLogStdout(6)  =         "  SelCriter"
      This%ColnameLogStdout(7)  =         "  SelIntens"
      This%ColnameLogStdout(8)  =         "     ...Pct"
      This%ColnameLogStdout(9)  =         " Coancestry"
      This%ColnameLogStdout(10) =         "    ...Rate"
      This%ColnameLogStdout(11) =         "     ...Pct"
      This%ColnameLogStdout(12) =         " Inbreeding"
      This%ColnameLogStdout(13) =         "    ...Rate"
      This%ColnameLogStdout(14) =         "     ...Pct"

      nColTmp = nCol
      if (Spec%GenericIndCritGiven) then
        do i = 1, Spec%nGenericIndCrit
          nColTmp = nColTmp + 1
          !                               12345678901
          This%ColnameLogUnit(nColTmp) = "GenIndCrit"//trim(Int2Char(i))
          This%ColnameLogUnit(nColTmp) = adjustr(This%ColnameLogUnit(nColTmp))
          Tmp = This%ColnameLogUnit(nColTmp)
          This%ColnameLogStdout(nColTmp) = Tmp(11:21)
          This%ColnameLogStdout(nColTmp) = adjustr(This%ColnameLogStdout(nColTmp))
        end do
      end if
      if (Spec%GenericMatCritGiven) then
        do i = 1, Spec%nGenericMatCrit
          nColTmp = nColTmp + 1
          !                               12345678901
          This%ColnameLogUnit(nColTmp) = "GenMatCrit"//trim(Int2Char(i))
          This%ColnameLogUnit(nColTmp) = adjustr(This%ColnameLogUnit(nColTmp))
          Tmp = This%ColnameLogUnit(nColTmp)
          This%ColnameLogStdout(nColTmp) = Tmp(11:21)
          This%ColnameLogStdout(nColTmp) = adjustr(This%ColnameLogStdout(nColTmp))
        end do
      end if
      This%ColnameLogPopUnit = This%ColnameLogUnit
      !                            1234567890123456789012
      This%ColnameLogPopUnit(2) = "              Solution"
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
    !> @brief  Write head of the AlphaMate log
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !> @return Output to file or standard output
    !---------------------------------------------------------------------------
    subroutine LogHeadAlphaMateSol(This, LogUnit, String, StringNum) ! not pure due to IO
      implicit none
      class(AlphaMateSol), intent(in)        :: This      !< AlphaMateSol holder
      integer(int32), intent(in), optional   :: LogUnit   !< Unit to write to (default STDOUT)
      character(len=*), intent(in), optional :: String    !< Additional string that will be written before the head
      integer(int32), optional               :: StringNum !< How much space is needed for the String
      character(len=10) :: StringFmt

      if (present(String)) then
        if (present(StringNum)) then
          StringFmt = "("//trim(Int2Char(StringNum))//"a)"
        else
          StringFmt = "(a)"
        end if
      end if
      if (present(LogUnit)) then
        if (present(String)) then
          write(LogUnit, StringFmt, Advance="No") trim(adjustl(String))
        end if
        write(LogUnit, This%FmtLogUnitHead)  This%ColnameLogUnit
      else
        if (present(String)) then
          write(STDOUT, StringFmt, Advance="No") trim(adjustl(String))
        end if
        write(STDOUT, This%FmtLogStdoutHead) This%ColnameLogStdout
      end if
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Write the AlphaMate log
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !> @return Output to file or standard output
    !---------------------------------------------------------------------------
    subroutine LogAlphaMateSol(This, LogUnit, Iteration, AcceptPct, String, StringNum) ! not pure due to IO
      implicit none
      class(AlphaMateSol), intent(in)              :: This      !< AlphaMateSol holder
      integer(int32), intent(in), optional         :: LogUnit   !< Unit to write to (default STDOUT)
      integer(int32), intent(in)                   :: Iteration !< Generation/Iteration
      real(real64), intent(in)                     :: AcceptPct !< Acceptance rate
      character(len=*), intent(in), optional       :: String    !< Additional string that will be written before the head
      integer(int32), optional                     :: StringNum !< How much space is needed for the String
      integer(int32) :: Unit
      character(len=CHARLENGTH) :: Fmt, StringFmt

      if (present(LogUnit)) then
        Unit = LogUnit
        Fmt = This%FmtLogUnit
      else
        Unit = STDOUT
        Fmt = This%FmtLogStdout
      end if
      if (present(String)) then
        if (present(StringNum)) then
          StringFmt = "("//trim(Int2Char(StringNum))//"a)"
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
                           This%SelCriterion, This%SelIntensity, This%MaxCriterionPct, &
                           This%CoancestryRanMate, This%CoancestryRateRanMate, This%MinCoancestryPct, &
                           This%Inbreeding, This%InbreedingRate, This%MinInbreedingPct, &
                           This%GenericIndCrit, This%GenericMatCrit
        else
          if (present(String)) then
            write(Unit, StringFmt, Advance="No") trim(adjustl(String))
          end if
          write(Unit, Fmt) Iteration, AcceptPct, This%Objective, This%Penalty, This%Degree, &
                           This%SelCriterion, This%SelIntensity, This%MaxCriterionPct, &
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
                           This%SelCriterion, This%SelIntensity, This%MaxCriterionPct, &
                           This%CoancestryRanMate, This%CoancestryRateRanMate, This%MinCoancestryPct, &
                           This%Inbreeding, This%InbreedingRate, This%MinInbreedingPct, &
                                                This%GenericMatCrit
        else
          if (present(String)) then
            write(Unit, StringFmt, Advance="No") trim(adjustl(String))
          end if
          write(Unit, Fmt) Iteration, AcceptPct, This%Objective, This%Penalty, This%Degree, &
                           This%SelCriterion, This%SelIntensity, This%MaxCriterionPct, &
                           This%CoancestryRanMate, This%CoancestryRateRanMate, This%MinCoancestryPct, &
                           This%Inbreeding, This%InbreedingRate, This%MinInbreedingPct
        end if
      end if
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Write head of the AlphaMate log - for the swarm/population of solutions
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !> @return Output to file or standard output
    !---------------------------------------------------------------------------
    subroutine LogPopHeadAlphaMateSol(This, LogPopUnit) ! not pure due to IO
      implicit none
      class(AlphaMateSol), intent(in)      :: This       !< AlphaMateSol holder
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
    !> @brief  Write the AlphaMate log - for the swarm/population of solutions
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !> @return Output to file or standard output
    !---------------------------------------------------------------------------
    subroutine LogPopAlphaMateSol(This, LogPopUnit, Iteration, i) ! not pure due to IO
      implicit none
      class(AlphaMateSol), intent(in)      :: This       !< AlphaMateSol holder
      integer(int32), intent(in), optional :: LogPopUnit !< population log file unit (default STDOUT)
      integer(int32), intent(in)           :: Iteration  !< generation/iteration
      integer(int32), intent(in)           :: i          !< solution id
      integer(int32) :: Unit

      if (present(LogPopUnit)) then
        Unit = LogPopUnit
      else
        Unit = STDOUT
      end if
      if (allocated(This%GenericIndCrit)) then
        if (allocated(This%GenericMatCrit)) then
          write(Unit, This%FmtLogPopUnit) Iteration, i, This%Objective, This%Penalty, This%Degree, &
                                          This%SelCriterion, This%SelIntensity, This%MaxCriterionPct, &
                                          This%CoancestryRanMate, This%CoancestryRateRanMate, This%MinCoancestryPct, &
                                          This%Inbreeding, This%InbreedingRate, This%MinInbreedingPct, &
                                          This%GenericIndCrit, This%GenericMatCrit
        else
          write(Unit, This%FmtLogPopUnit) Iteration, i, This%Objective, This%Penalty, This%Degree, &
                                          This%SelCriterion, This%SelIntensity, This%MaxCriterionPct, &
                                          This%CoancestryRanMate, This%CoancestryRateRanMate, This%MinCoancestryPct, &
                                          This%Inbreeding, This%InbreedingRate, This%MinInbreedingPct, &
                                          This%GenericIndCrit
        end if
      else
        if (allocated(This%GenericMatCrit)) then
          write(Unit, This%FmtLogPopUnit) Iteration, i, This%Objective, This%Penalty, This%Degree, &
                                          This%SelCriterion, This%SelIntensity, This%MaxCriterionPct, &
                                          This%CoancestryRanMate, This%CoancestryRateRanMate, This%MinCoancestryPct, &
                                          This%Inbreeding, This%InbreedingRate, This%MinInbreedingPct, &
                                                               This%GenericMatCrit
        else
          write(Unit, This%FmtLogPopUnit) Iteration, i, This%Objective, This%Penalty, This%Degree, &
                                          This%SelCriterion, This%SelIntensity, This%MaxCriterionPct, &
                                          This%CoancestryRanMate, This%CoancestryRateRanMate, This%MinCoancestryPct, &
                                          This%Inbreeding, This%InbreedingRate, This%MinInbreedingPct
        end if
      end if
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert selection criterion to selection intensity
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    April 18, 2017
    !---------------------------------------------------------------------------
    elemental function SelCriterion2SelIntensity(SelCriterion, Mean, Sd) result(SelIntensity)
      implicit none
      real(real64), intent(in) :: SelCriterion !< selection criterion
      real(real64), intent(in) :: Mean         !< mean of selection criterion
      real(real64), intent(in) :: Sd           !< standard deviation of selection criterion
      real(real64)             :: SelIntensity !< @return selection intensity
      SelIntensity = (SelCriterion - Mean) / Sd
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert selection intensity to selection criterion
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    April 18, 2017
    !---------------------------------------------------------------------------
    elemental function SelIntensity2SelCriterion(SelIntensity, Mean, Sd) result(SelCriterion)
      implicit none
      real(real64), intent(in) :: SelIntensity !< selection intensity
      real(real64), intent(in) :: Mean         !< mean of selection criterion
      real(real64), intent(in) :: Sd           !< standard deviation of selection criterion
      real(real64)             :: SelCriterion !< @return selection criterion
      SelCriterion = SelIntensity * Sd + Mean
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
      real(real64), intent(in) :: CurrentCoancestry !< Current coancestry
      real(real64), intent(in) :: FutureCoancestry  !< Future coancestry
      real(real64)             :: CoancestryRate  !< @return CoancestryRate
      CoancestryRate = (FutureCoancestry - CurrentCoancestry) / (1.0d0 - CurrentCoancestry)
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
      real(real64), intent(in) :: CoancestryRate    !< CoancestryRate
      real(real64), intent(in) :: CurrentCoancestry !< Current coancestry
      real(real64)             :: FutureCoancestry  !< @return Future coancestry
      FutureCoancestry = CoancestryRate + (1.0d0 - CoancestryRate) * CurrentCoancestry
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
      real(real64), intent(in) :: MinCoancestryPct !< percentage of minimum coancestry achieved (100 means we achieved the minimum possible coancestry)
      real(real64)             :: Degree           !< @return Frontier degree
      Degree = asin(MinCoancestryPct / 100.0d0) * RAD2DEG
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert Frontier degree to MinCoancestryPct
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    April 13, 2017
    !> @details Assuming unit circular selection/coancestry frontier and given the
    !!          the angle between the maximum selection line and the solution line
    !!          we can evaluate the percentage of maximum selection intensity achieved
    !!          (y-axis) by recognising that sin(angle) = opposite/hypothenuse, where
    !!          opposite = MinCoancestryPct/100 and hypothenuse = 1 (unit circle). Then
    !!          MinCoancestryPct = sin(angle) * 100.
    !---------------------------------------------------------------------------
    pure function Degree2MinCoancestryPct(Degree) result(MinCoancestryPct)
      implicit none
      real(real64), intent(in) :: Degree           !< Frontier degree
      real(real64)             :: MinCoancestryPct !< @return Percentage of minimum coancestry achieved (100 means we achieved the minimum possible coancestry)
      MinCoancestryPct = sin(Degree * DEG2RAD) * 100.0d0
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert MaxCriterionPct to Frontier degree
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    April 13, 2017
    !> @details Assuming unit circular selection/coancestry frontier and given the
    !!          percentage of maximum selection intensity achieved (y-axis) we can evaluate
    !!          the angle between the maximum selection line and the solution line
    !!          recognising that cos(angle) = adjacent/hypothenuse, where
    !!          adjacent = MaxCriterionPct/100 and hypothenuse = 1 (unit circle). Then
    !!          angle = acos(MaxCriterionPct/100).
    !---------------------------------------------------------------------------
    pure function MaxCriterionPct2Degree(MaxCriterionPct) result(Degree)
      implicit none
      real(real64), intent(in) :: MaxCriterionPct !< percentage of maximum criterion achieved (100 means we achieved the maximum possible selection intensity)
      real(real64)             :: Degree          !< @return Frontier degree
      Degree = acos(MaxCriterionPct / 100.0d0) * RAD2DEG
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert Frontier degree to MaxCriterionPct
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    April 13, 2017
    !> @details Assuming unit circular selection/coancestry frontier and given the
    !!          the angle between the maximum selection line and the solution line
    !!          we can evaluate the percentage of maximum selection intensity achieved
    !!          (y-axis) by recognising that cos(angle) = adjacent/hypothenuse, where
    !!          adjacent = MaxCriterionPct/100 and hypothenuse = 1 (unit circle). Then
    !!          MaxCriterionPct = cos(angle) * 100.
    !---------------------------------------------------------------------------
    pure function Degree2MaxCriterionPct(Degree) result(MaxCriterionPct)
      implicit none
      real(real64), intent(in) :: Degree          !< Frontier degree
      real(real64)             :: MaxCriterionPct !< @return Percentage of maximum criterion achieved (100 means we achieved the maximum possible selection intensity)
      MaxCriterionPct = cos(Degree * DEG2RAD) * 100.0d0
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
      real(real64), intent(in) :: MinCoancestryPct  !< MinCoancestryPct of a solution
      real(real64), intent(in) :: MinCoancestryRate !< Minimum possible coancestry rate
      real(real64), intent(in) :: MaxCoancestryRate !< Maximum possible coancestry rate
      real(real64)             :: CoancestryRate    !< @return Coancestry rate at a given MinCoancestryPct
      CoancestryRate = MinCoancestryRate + (100.0d0 - MinCoancestryPct) / 100.0d0 * (MaxCoancestryRate - MinCoancestryRate)
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
      real(real64), intent(in) :: CoancestryRate    !< Coancestry rate at a given MinCoancestryPct
      real(real64), intent(in) :: MinCoancestryRate !< Minimum possible coancestry rate
      real(real64), intent(in) :: MaxCoancestryRate !< Maximum possible coancestry rate
      real(real64)             :: MinCoancestryPct  !< @return MinCoancestryPct of a solution
      MinCoancestryPct = (MaxCoancestryRate - CoancestryRate) / (MaxCoancestryRate - MinCoancestryRate) * 100.d0
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert MaxCriterionPct to selection intensity
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    April 13, 2017
    !> @details Compute difference between Max and Min and take MaxCriterionPct of
    !!          the difference as we express MaxCriterionPct as 100% when solution
    !!          achieves Max criterion
    !---------------------------------------------------------------------------
    pure function MaxCriterionPct2SelIntensity(MaxCriterionPct, MinSelIntensity, MaxSelIntensity) result(SelIntensity)
      implicit none
      real(real64), intent(in) :: MaxCriterionPct !< MaxCriterionPct of a solution
      real(real64), intent(in) :: MinSelIntensity !< Minimum possible selection intensity
      real(real64), intent(in) :: MaxSelIntensity !< Maximum possible selection intensity
      real(real64)             :: SelIntensity    !< @return Selection intensity at a given MaxCriterionPct
      SelIntensity = MinSelIntensity + MaxCriterionPct / 100.0d0 * (MaxSelIntensity - MinSelIntensity)
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert selection intensity to MaxCriterionPct
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    April 18, 2017
    !> @details See how close SelIntensity is to MaxSelIntensity
    !---------------------------------------------------------------------------
    pure function SelIntensity2MaxCriterionPct(SelIntensity, MinSelIntensity, MaxSelIntensity) result(MaxCriterionPct)
      implicit none
      real(real64), intent(in) :: SelIntensity    !< Selection intensity of a solution
      real(real64), intent(in) :: MinSelIntensity !< Minimum possible selection intensity
      real(real64), intent(in) :: MaxSelIntensity !< Maximum possible selection intensity
      real(real64)             :: MaxCriterionPct !< @return MaxCriterionPct of a solution
      MaxCriterionPct = (SelIntensity - MinSelIntensity) / (MaxSelIntensity - MinSelIntensity) * 100.0d0
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert degree to selection intensity
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    April 18, 2017
    !---------------------------------------------------------------------------
    pure function Degree2SelIntensity(Degree, MinSelIntensity, MaxSelIntensity) result(SelIntensity)
      implicit none
      real(real64), intent(in) :: Degree          !< Degree
      real(real64), intent(in) :: MinSelIntensity !< Minimum possible selection intensity
      real(real64), intent(in) :: MaxSelIntensity !< Maximum possible selection intensity
      real(real64)             :: SelIntensity    !< @return Selection intensity
      SelIntensity = MaxCriterionPct2SelIntensity(MaxCriterionPct=Degree2MaxCriterionPct(Degree=Degree), &
                                                  MinSelIntensity=MinSelIntensity, &
                                                  MaxSelIntensity=MaxSelIntensity)
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert selection intensity to degree
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    April 18, 2017
    !---------------------------------------------------------------------------
    pure function SelIntensity2Degree(SelIntensity, MinSelIntensity, MaxSelIntensity) result(Degree)
      implicit none
      real(real64), intent(in) :: SelIntensity    !< Selection intensity
      real(real64), intent(in) :: MinSelIntensity !< Minimum possible selection intensity
      real(real64), intent(in) :: MaxSelIntensity !< Maximum possible selection intensity
      real(real64)             :: Degree          !< @return Degree
      Degree = MaxCriterionPct2Degree(MaxCriterionPct=SelIntensity2MaxCriterionPct(SelIntensity=SelIntensity, &
                                                                                   MinSelIntensity=MinSelIntensity, &
                                                                                   MaxSelIntensity=MaxSelIntensity))
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert degree to coancestry rate
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    April 18, 2017
    !---------------------------------------------------------------------------
    pure function Degree2CoancestryRate(Degree, MinCoancestryRate, MaxCoancestryRate) result(CoancestryRate)
      implicit none
      real(real64), intent(in) :: Degree            !< Degree
      real(real64), intent(in) :: MinCoancestryRate !< Minimum possible coancestry rate
      real(real64), intent(in) :: MaxCoancestryRate !< Maximum possible coancestry rate
      real(real64)             :: CoancestryRate    !< @return Coancestry rate
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
      real(real64), intent(in) :: CoancestryRate      !< Coancestry rate
      real(real64), intent(in) :: MinCoancestryRate !< Minimum possible coancestry rate
      real(real64), intent(in) :: MaxCoancestryRate !< Maximum possible coancestry rate
      real(real64)             :: Degree            !< @return Degree
      Degree = MinCoancestryPct2Degree(MinCoancestryPct=CoancestryRate2MinCoancestryPct(CoancestryRate=CoancestryRate, &
                                                                                        MinCoancestryRate=MinCoancestryRate, &
                                                                                        MaxCoancestryRate=MaxCoancestryRate))
    end function

    !###########################################################################

end module

!###############################################################################