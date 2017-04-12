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
  use ConstantModule, only : FILELENGTH, SPECOPTIONLENGTH, IDLENGTH
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

  ! Module parameters
  INTEGER,                   PARAMETER :: CHARLENGTH = 100
  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTREAL2CHAR = "(f11.5)"
  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTINT2CHAR  = "(i11)"

  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTLOGSTDOUTHEADA = "("
  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTLOGSTDOUTHEADB = "a15)"

  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTLOGSTDOUTA = "(i15, "
  CHARACTER(len=CHARLENGTH), PARAMETER :: FMTLOGSTDOUTB = "(4x, f11.5))"

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
    logical :: ModeMin, ModeMinPct, ModeMax, ModeMaxPct, ModeOpt, ModeRan, ModeFrontier
    integer(int32) :: nFrontierPoints

    ! Biological specifications
    logical :: NrmInsteadOfCoancestry
    real(real64) :: TargetDegree, TargetCoancestryRate, TargetInbreedingRate
    real(real64) :: MinPct, MaxPct
    real(real64), allocatable :: TargetDegreeFrontier(:), TargetCoancestryRateFrontier(:)
    real(real64) :: TargetDegreeWeight, TargetCoancestryRateWeight, TargetInbreedingRateWeight, SelfingWeight
    integer(int32) :: nInd, nMat, nPar, nPar1, nPar2 ! NOTE: nInd is here just for OO-flexibility (do not use it; the main one is in Data!!!)
    logical :: EqualizePar, EqualizePar1, EqualizePar2, LimitPar, LimitPar1, LimitPar2
    real(real64) :: LimitParMin, LimitPar1Min, LimitPar2Min, LimitParMax, LimitPar1Max, LimitPar2Max, LimitParMinWeight, LimitPar1MinWeight, LimitPar2MinWeight
    integer(int32) :: nGenericIndCrit, nGenericMatCrit
    real(real64), allocatable :: GenericIndCritWeight(:), GenericMatCritWeight(:)
    logical :: SelfingAllowed, TargetDegreeWeightBelow,  TargetCoancestryRateWeightBelow, TargetInbreedingRateWeightBelow
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

    ! Data&Spec derived quantities
    ! @todo (not sure if this should be in data though?)
    real(real64) :: MaxSelCriterion, MaxSelIntensity, MinSelCriterion, MinSelIntensity
    real(real64) :: MaxCoancestry, MaxCoancestryRate, MinCoancestry, MinCoancestryRate

    contains
      procedure :: Initialise => InitialiseAlphaMateSpec
      procedure :: Read       => ReadAlphaMateSpec
      procedure :: Write      => WriteAlphaMateSpec
  end type

  !> @brief AlphaMate data
  type, extends(AlphaEvolveData) :: AlphaMateData
    ! Raw data
    type(RelMat) :: Coancestry
    type(InbVec) :: Inbreeding
    real(real64), allocatable :: SelCriterion(:), SelCriterionStand(:), SelCriterionPAGE(:), SelCriterionPAGEStand(:)
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
    real(real64) :: CurrentCoancestryRanMate, CurrentCoancestryRanMateNoSelf, CurrentCoancestryGenderMate
    real(real64) :: CurrentInbreeding
! @todo Should this be in the spec type, given that these are quite different with different modes?
    real(real64) :: TargetCoancestryRanMate, TargetCoancestryRanMateNoSelf, TargetCoancestryGenderMate
    real(real64) :: TargetInbreeding
    contains
      procedure :: Read  => ReadAlphaMateData
      procedure :: Write => WriteAlphaMateData
  end type

  !> @brief AlphaMate solution
  type, extends(AlphaEvolveSol) :: AlphaMateSol
    real(real64)                :: Penalty
    real(real64)                :: PenaltyCoancestryRate
    real(real64)                :: PenaltyInbreedingRate
    real(real64)                :: PenaltySelfing
    real(real64)                :: PenaltyLimitPar1
    real(real64)                :: PenaltyLimitPar2
    real(real64)                :: PenaltyGenericIndCrit
    real(real64)                :: PenaltyGenericMatCrit
    real(real64)                :: SelCriterion
    real(real64)                :: SelIntensity
    real(real64)                :: MaxPct
    real(real64)                :: FutureCoancestryRanMate
    real(real64)                :: CoancestryRateRanMate
    real(real64)                :: MinPct
    real(real64)                :: Degree
    real(real64)                :: FutureInbreeding
    real(real64)                :: InbreedingRate
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
      write(STDOUT, "(a)") ""
      write(STDOUT, "(a)") "                            ***********************                           "
      write(STDOUT, "(a)") "                            *                     *                           "
      write(STDOUT, "(a)") "                            *      AlphaMate      *                           "
      write(STDOUT, "(a)") "                            *                     *                           "
      write(STDOUT, "(a)") "                            ***********************                           "
      write(STDOUT, "(a)") "                                                                              "
      write(STDOUT, "(a)") "         Software for optimizing contributions to the next generation         "
      write(STDOUT, "(a)") "                       http://AlphaGenes.Roslin.ed.ac.uk                      "
      write(STDOUT, "(a)") "                                 No liability                                 "
      write(STDOUT, "(a)") ""
      write(STDOUT, "(a)") "                       Commit:   "//TOSTRING(COMMIT)//"                       "
      write(STDOUT, "(a)") "                       Compiled: "//__DATE__//", "//__TIME__
      write(STDOUT, "(a)") ""
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

      This%ModeMin = .false.
      This%ModeMinPct = .false.
      This%ModeMax = .false.
      This%ModeMaxPct = .false.
      This%ModeOpt = .false.
      This%ModeRan = .false.
      This%ModeFrontier = .false.
      This%nFrontierPoints = 0

      ! Biological specifications

      This%nInd = 0
      This%nMat = 0
      This%nPar = 0
      This%nPar1 = 0
      This%nPar2 = 0

      This%TargetDegree = 45
      This%TargetCoancestryRate = 0.01d0
      This%TargetInbreedingRate = 0.01d0
      ! This%TargetDegreeFrontier ! allocatable so skip here
      ! This%TargetCoancestryRateFrontier ! allocatable so skip here
      This%TargetDegreeWeight = -100.0d0
      This%TargetCoancestryRateWeightBelow = .false.
      This%TargetCoancestryRateWeight = -100.0d0
      This%TargetCoancestryRateWeightBelow = .false.
      This%TargetInbreedingRateWeight =    0.0d0
      This%TargetInbreedingRateWeightBelow = .false.
      This%MinPct = 0
      This%MaxPct = 0
      This%SelfingAllowed = .false.
      This%SelfingWeight = 0.0d0
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
      This%LimitParMinWeight  = 0.0d0
      This%LimitPar1MinWeight = 0.0d0
      This%LimitPar2MinWeight = 0.0d0

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

      ! Data&Spec derived quantities

      This%MaxSelCriterion   = IEEE_Value(x=This%MaxSelCriterion,   class=IEEE_Quiet_NaN)
      This%MaxSelIntensity   = IEEE_Value(x=This%MaxSelIntensity,   class=IEEE_Quiet_NaN)
      This%MinSelCriterion   = IEEE_Value(x=This%MinSelCriterion,   class=IEEE_Quiet_NaN)
      This%MinSelIntensity   = IEEE_Value(x=This%MinSelIntensity,   class=IEEE_Quiet_NaN)
      This%MaxCoancestry     = IEEE_Value(x=This%MaxCoancestry,     class=IEEE_Quiet_NaN)
      This%MaxCoancestryRate = IEEE_Value(x=This%MaxCoancestryRate, class=IEEE_Quiet_NaN)
      This%MinCoancestry     = IEEE_Value(x=This%MinCoancestry,     class=IEEE_Quiet_NaN)
      This%MinCoancestryRate = IEEE_Value(x=This%MinCoancestryRate, class=IEEE_Quiet_NaN)
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

      write(Unit, *) "ModeMin: ",          This%ModeMin
      write(Unit, *) "ModeMinPct: ",       This%ModeMinPct
      write(Unit, *) "ModeMax: ",          This%ModeMax
      write(Unit, *) "ModeMaxPct: ",       This%ModeMaxPct
      write(Unit, *) "ModeOpt: ",          This%ModeOpt
      write(Unit, *) "ModeRan: ",          This%ModeRan
      write(Unit, *) "ModeFrontier: ",     This%ModeFrontier
      write(Unit, *) "nFrontierPoints: ",  This%nFrontierPoints

      ! Biological specifications

      write(Unit, *) "nInd: ",  This%nInd
      write(Unit, *) "nMat: ",  This%nMat
      write(Unit, *) "nPar: ",  This%nPar
      write(Unit, *) "nPar1: ", This%nPar1
      write(Unit, *) "nPar2: ", This%nPar2

      write(Unit, *) "TargetDegree: ",                    This%TargetDegree
      write(Unit, *) "TargetCoancestryRate: ",            This%TargetCoancestryRate
      write(Unit, *) "TargetInbreedingRate: ",            This%TargetInbreedingRate
      if (allocated(This%TargetDegreeFrontier)) then
        write(Unit, *) "TargetDegreeFrontier: ",          This%TargetDegreeFrontier
      else
        write(Unit, *) "TargetDegreeFrontier: not allocated"
      end if
      if (allocated(This%TargetCoancestryRateFrontier)) then
        write(Unit, *) "TargetCoancestryRateFrontier: ",  This%TargetCoancestryRateFrontier
      else
        write(Unit, *) "TargetCoancestryRateFrontier: not allocated"
      end if
      write(Unit, *) "TargetDegreeWeight: ",              This%TargetDegreeWeight
      write(Unit, *) "TargetDegreeWeightBelow: ",         This%TargetDegreeWeightBelow
      write(Unit, *) "TargetCoancestryRateWeight: ",      This%TargetCoancestryRateWeight
      write(Unit, *) "TargetCoancestryRateWeightBelow: ", This%TargetCoancestryRateWeightBelow
      write(Unit, *) "TargetInbreedingRateWeight: ",      This%TargetInbreedingRateWeight
      write(Unit, *) "TargetInbreedingRateWeightBelow: ", This%TargetInbreedingRateWeightBelow
      write(Unit, *) "MinPct: ",                          This%MinPct
      write(Unit, *) "MaxPct: ",                          This%MaxPct
      if (allocated(This%TargetDegreeFrontier)) then
        write(Unit, *) "TargetDegreeFrontier: ",          This%TargetDegreeFrontier
      else
        write(Unit, *) "TargetDegreeFrontier: not allocated"
      end if
      if (allocated(This%TargetCoancestryRateFrontier)) then
        write(Unit, *) "TargetCoancestryRateFrontier: ",  This%TargetCoancestryRateFrontier
      else
        write(Unit, *) "TargetCoancestryRateFrontier: not allocated"
      end if
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

      write(Unit, *) "MaxSelCriterion: ",   This%MaxSelCriterion
      write(Unit, *) "MaxSelIntensity: ",   This%MaxSelIntensity
      write(Unit, *) "MinSelCriterion: ",   This%MinSelCriterion
      write(Unit, *) "MinSelIntensity: ",   This%MinSelIntensity
      write(Unit, *) "MaxCoancestry: ",     This%MaxCoancestry
      write(Unit, *) "MaxCoancestryRate: ", This%MaxCoancestryRate
      write(Unit, *) "MinCoancestry: ",     This%MinCoancestry
      write(Unit, *) "MinCoancestryRate: ", This%MinCoancestryRate

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

      integer(int32) :: SpecUnit, Stat, nFrontierPoint, nGenericIndCrit, nGenericMatCrit

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
                write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                  write(STDERR, "(a)") ""
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a file for NrmMatrixFile, i.e., NrmMatrixFile, NrmMtx.txt"
                write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                  write(STDERR, "(a)") ""
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
                  write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                  write(STDERR, "(a)") ""
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
                  write(STDERR, "(a)") ""
                  stop 1
                end if
              end if

            ! Search mode specifications
            case ("modemin")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%ModeMin = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " ModeMin"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for ModeMin, i.e., ModeMin, Yes"
                write(STDERR, "(a)") ""
                stop 1
              end if

            case ("modeminpct")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%ModeMinPct = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " ModeMinPct"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for ModeMinPct, i.e., ModeMinPct, Yes"
                write(STDERR, "(a)") ""
                stop 1
              end if

            case ("modemax")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%ModeMax = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " ModeMax"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for ModeMax, i.e., ModeMax, Yes"
                write(STDERR, "(a)") ""
                stop 1
              end if

            case ("modemaxpct")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%ModeMaxPct = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " ModeMaxPct"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for ModeMaxPct, i.e., ModeMaxPct, Yes"
                write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
                stop 1
              end if

            case ("modefrontier")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%ModeFrontier = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " ModeFrontier"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for ModeFrontier, i.e., ModeFrontier, Yes"
                write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
                stop 1
              end if

            case ("targetdegree")
              if (allocated(Second)) then
                This%TargetDegree = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Targeted degree: "//trim(Real2Char(This%TargetDegree, fmt=FMTREAL2CHAR))
                end if
                if ((This%TargetDegree .lt. 0.0d0) .or. (This%TargetDegree .gt. 90.0d0)) then
                  write(STDERR, "(a)") "ERROR: TargetDegree must be between 0 and 90!"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetDegree, i.e., TargetDegree, 45"
                write(STDERR, "(a)") ""
                stop 1
              end if

            case ("targetdegreeweight")
              if (allocated(Second)) then
                This%TargetDegreeWeight = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Targeted degree - weight: "//trim(Real2Char(This%TargetDegreeWeight, fmt=FMTREAL2CHAR))
                end if
                if (This%TargetDegreeWeight .gt. 0.0d0) then
                  write(STDOUT, "(a)") " NOTE: Positive weight for targeted degree, i.e., encourage higher rate. Was this intended?"
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetDegreeWeight, i.e., TargetDegreeWeight, -100"
                write(STDERR, "(a)") ""
                stop 1
              end if

            case ("targetdegreeweightbelow")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%TargetDegreeWeightBelow = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Targeted degree - weight also values below the target"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for TargetDegreeWeightBelow, i.e., TargetDegreeWeightBelow, No"
                write(STDERR, "(a)") ""
                stop 1
              end if

            case ("targetcoancestryrate")
              if (allocated(Second)) then
                This%TargetCoancestryRate = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Targeted rate of coancestry: "//trim(Real2Char(This%TargetCoancestryRate, fmt=FMTREAL2CHAR))
                end if
                ! if (This%TargetInbreedingRate .eq. 0.0) then
                !   write(STDERR, "(a)") "ERROR: Can not work with the targeted rate of inbreeding exactly equal to zero - it is numerically unstable!"
                !   write(STDERR, "(a)") " "
                !   stop 1
                ! end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetedRateOfCoancestry, i.e., TargetedRateOfCoancestry, 0.01"
                write(STDERR, "(a)") ""
                stop 1
              end if

            case ("targetcoancestryrateweight")
              if (allocated(Second)) then
                This%TargetCoancestryRateWeight = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Targeted rate of coancestry - weight: "//trim(Real2Char(This%TargetCoancestryRateWeight, fmt=FMTREAL2CHAR))
                end if
                if (This%TargetCoancestryRateWeight .gt. 0.0d0) then
                  write(STDOUT, "(a)") " NOTE: Positive weight for targeted rate of coancestry, i.e., encourage higher rate. Was this intended?"
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetCoancestryRateWeight, i.e., TargetCoancestryRateWeight, -100"
                write(STDERR, "(a)") ""
                stop 1
              end if

            case ("targetcoancestryrateweightbelow")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%TargetCoancestryRateWeightBelow = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Targeted rate of coancestry - weight also values below the target"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for TargetCoancestryRateWeightBelow, i.e., TargetCoancestryRateWeightBelow, No"
                write(STDERR, "(a)") ""
                stop 1
              end if

            case ("targetinbreedingrate")
              if (allocated(Second)) then
                This%TargetInbreedingRate = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Targeted rate of inbreeding: "//trim(Real2Char(This%TargetInbreedingRate, fmt=FMTREAL2CHAR))
                end if
                ! if (This%TargetInbreedingRate .eq. 0.0) then
                !   write(STDERR, "(a)") "ERROR: Can not work with the targeted rate of inbreeding exactly equal to zero - it is numerically unstable!"
                !   write(STDERR, "(a)") " "
                !   stop 1
                ! end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetInbreedingRate, i.e., TargetInbreedingRate, 0.01"
                write(STDERR, "(a)") ""
                stop 1
              end if



            case ("targetinbreedingrateweight")
              if (allocated(Second)) then
                This%TargetInbreedingRateWeight = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Targeted rate of inbreeding - weight: "//trim(Real2Char(This%TargetInbreedingRateWeight, fmt=FMTREAL2CHAR))
                end if
                if (This%TargetInbreedingRateWeight .gt. 0.0d0) then
                  write(STDOUT, "(a)") " NOTE: Positive weight for targeted rate of inbreeding, i.e., encourage higher rate. Was this intended?"
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetInbreedingRateWeight, i.e., TargetInbreedingRateWeight, -100"
                write(STDERR, "(a)") ""
                stop 1
              end if

            case ("targetinbreedingrateweightbelow")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%TargetInbreedingRateWeightBelow = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Targeted rate of inbreeding - weight also values below the target"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for TargetInbreedingRateWeightBelow, i.e., TargetInbreedingRateWeightBelow, No"
                write(STDERR, "(a)") ""
                stop 1
              end if

            case ("minpct")
              if (allocated(Second)) then
                This%MinPct = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " ModeMin percentage: "//trim(Real2Char(This%MinPct, fmt=FMTREAL2CHAR))
                end if
                if (This%MinPct .lt. 0.0d0 .or. This%MinPct .gt. 100.0d0) then
                  write(STDERR, "(a)") "ERROR: MinPct must be between 0 and 100!"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for MinPct, i.e., MinPct, 10"
                write(STDERR, "(a)") ""
                stop 1
              end if

            case ("maxpct")
              if (allocated(Second)) then
                This%MaxPct = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " ModeMax percentage: "//trim(Real2Char(This%MaxPct, fmt=FMTREAL2CHAR))
                end if
                if (This%MaxPct .lt. 0.0d0 .or. This%MaxPct .gt. 100.0d0) then
                  write(STDERR, "(a)") "ERROR: MaxPct must be between 0 and 100!"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for MaxPct, i.e., MaxPct, 10"
                write(STDERR, "(a)") ""
                stop 1
              end if

            case ("frontiernumberofpoints")
              if (This%ModeFrontier) then
                if (allocated(Second)) then
                  This%nFrontierPoints = Char2Int(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Evaluate selection/coancestry frontier - number of points: "//trim(Int2Char(This%nFrontierPoints))
                  end if
                  allocate(This%TargetCoancestryRateFrontier(This%nFrontierPoints))
                  nFrontierPoint = 0
                else
                  write(STDERR, "(a)") " ERROR: Must specify a number for FrontierNumberOfPoints, i.e., FrontierNumberOfPoints, 3"
                  write(STDERR, "(a)") ""
                  stop 1
                end if
              end if

            case ("frontiertargetdegree")
              if (This%ModeFrontier) then
                if (allocated(Second)) then
                  nFrontierPoint = nFrontierPoint + 1
                  if (nFrontierPoint .le. This%nFrontierPoints) then
                    This%TargetDegreeFrontier(nFrontierPoint) = Char2Double(trim(adjustl(Second(1))))
                    if (LogStdoutInternal) then
                      write(STDOUT, "(a)") " Evaluate selection/coancestry frontier - degree ("//trim(Int2Char(nFrontierPoint))//"): "//trim(Real2Char(This%TargetDegreeFrontier(nFrontierPoint), fmt=FMTREAL2CHAR))
                    end if
                    if ((This%TargetDegreeFrontier(nFrontierPoint) .lt. 0.0d0) .or. (This%TargetDegreeFrontier(nFrontierPoint) .gt. 90.0d0)) then
                      write(STDERR, "(a)") "ERROR: TargetDegreeFrontier must be between 0 and 90!"
                      write(STDERR, "(a)") " "
                      stop 1
                    end if
                  else
                    write(STDOUT, "(a)") " NOTE: Specification '"//trim(Line)//"' was ignored - already read all frontier points!"
                    write(STDOUT, "(a)") " "
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify a value for FrontierTargetCoancestryRate, i.e., FrontierTargetCoancestryRate, 0.001"
                  write(STDERR, "(a)") ""
                  stop 1
                end if
              end if

            case ("frontiertargetcoancestryrate")
              if (This%ModeFrontier) then
                if (allocated(Second)) then
                  nFrontierPoint = nFrontierPoint + 1
                  if (nFrontierPoint .le. This%nFrontierPoints) then
                    This%TargetCoancestryRateFrontier(nFrontierPoint) = Char2Double(trim(adjustl(Second(1))))
                    if (LogStdoutInternal) then
                      write(STDOUT, "(a)") " Evaluate selection/coancestry frontier - coancestry rate ("//trim(Int2Char(nFrontierPoint))//"): "//trim(Real2Char(This%TargetCoancestryRateFrontier(nFrontierPoint), fmt=FMTREAL2CHAR))
                    end if
                    ! if (This%TargetCoancestryRateFrontier(nFrontierPoint) .eq. 0.0d0) then
                    !   write(STDERR, "(a)") "ERROR: Can not work with the targeted rate of coancestry exactly equal to zero - it is numerically unstable!"
                    !   write(STDERR, "(a)") " "
                    !   stop 1
                    ! end if
                  else
                    write(STDOUT, "(a)") " NOTE: Specification '"//trim(Line)//"' was ignored - already read all frontier points!"
                    write(STDOUT, "(a)") " "
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify a value for FrontierTargetCoancestryRate, i.e., FrontierTargetCoancestryRate, 0.001"
                  write(STDERR, "(a)") ""
                  stop 1
                end if
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
                write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                  write(STDERR, "(a)") ""
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
                  write(STDERR, "(a)") ""
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
                  write(STDERR, "(a)") " ERROR: Must specify a value for LimitContributionsMinWeight, i.e., LimitContributionsMinWeight, -0.01"
                  write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                  write(STDERR, "(a)") ""
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
                  write(STDERR, "(a)") ""
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
                  write(STDERR, "(a)") " ERROR: Must specify a value for LimitMaleContributionsMinWeight, i.e., LimitMaleContributionsMinWeight, -0.01"
                  write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                  write(STDERR, "(a)") ""
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
                  write(STDERR, "(a)") ""
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
                  write(STDERR, "(a)") " ERROR: Must specify a value for LimitFemaleContributionsMinWeight, i.e., LimitFemaleContributionsMinWeight, -0.01"
                  write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                  write(STDERR, "(a)") " ERROR: Must specify a value for SelfingWeight, i.e., SelfingWeight, -0.01"
                  write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                  write(STDERR, "(a)") ""
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
                  write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                  write(STDERR, "(a)") ""
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
                  write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                  write(STDERR, "(a)") ""
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
                  write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
                stop 1
              end if

            case ("evolalg")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .ne. "none") then
                  write(This%EvolAlg, *) trim(adjustl(Second(1)))
                  This%EvolAlg = adjustl(This%EvolAlg)
                  if (.not. (This%EvolAlg .eq. "DE")) then
                    write(STDERR, "(a)") " ERROR: Must specify a valid algorithm [DE/???] for EvolAlg, i.e., EvolAlg, DE"
                    write(STDERR, "(a)") ""
                    stop 1
                  end if
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Evolutionary algorithm: "//trim(This%EvolAlg)
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a algorithm for EvolAlg, i.e., EvolAlg, DE"
                write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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
                write(STDERR, "(a)") ""
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

      if (.not. (This%ModeMin .or. This%ModeMinPct .or. &
                 This%ModeMax .or. This%ModeMaxPct .or. &
                 This%ModeOpt .or. &
                 This%ModeRan .or. &
                 This%ModeFrontier)) then
        write(STDERR, "(a)") " ERROR: One of the modes must be activated!"
        write(STDERR, "(a)") " ERROR: ModeMin, ModeMinPct, ModeMax, ModeMaxPct, ModeOpt, ModeRan, or ModeFrontier"
        write(STDERR, "(a)") ""
        stop 1
      end if

      if (.not. This%SelCriterionGiven .and. &
          (This%ModeMinPct .or. This%ModeMax .or. This%ModeMaxPct .or. This%ModeOpt .or. This%ModeFrontier)) then
        write(STDERR, "(a)") " ERROR: Selection criterion is needed for modes: ModeMinPct, ModeMax, ModeMaxPct, ModeOpt, and ModeFrontier"
        write(STDERR, "(a)") ""
        stop 1
      end if

      if (This%ModeMinPct .or. This%ModeMaxPct .or. This%ModeFrontier) then
        This%ModeMin = .true.
        This%ModeMax = .true.
      end if

      if (.not. This%RelMtxGiven) then
        write(STDERR, "(a)") " ERROR: One of CoancestryMatrixFile or NrmMatrixFile must be specified!"
        write(STDERR, "(a)") ""
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
        This%LimitParMinWeight  = 0.0d0
        This%LimitPar1MinWeight = 0.0d0
        This%LimitPar2MinWeight = 0.0d0
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
      allocate(This%SelCriterionStand(This%nInd))
      if (Spec%PAGEPar) then
        allocate(This%SelCriterionPAGE(This%nInd))
        allocate(This%SelCriterionPAGEStand(This%nInd))
      end if

      if (.not. Spec%SelCriterionGiven) then
        This%SelCriterion = 0.0d0
        This%SelCriterionStand = 0.0d0
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

      ! Current
      This%CurrentCoancestryRanMate       = This%CoancestryStat%All%Mean
      This%CurrentCoancestryRanMateNoSelf = This%CoancestryStat%OffDiag%Mean
      if (Spec%GenderGiven) then
        This%CurrentCoancestryGenderMate  = This%CoancestryStatGenderDiff%All%Mean
      end if

      ! Obtain limit/target based on given rates
      ! F_t = DeltaF + (1 - DeltaF) * F_t-1
      This%TargetCoancestryRanMate       = Spec%TargetCoancestryRate + (1.0d0 - Spec%TargetCoancestryRate) * This%CurrentCoancestryRanMate
      This%TargetCoancestryRanMateNoSelf = Spec%TargetCoancestryRate + (1.0d0 - Spec%TargetCoancestryRate) * This%CurrentCoancestryRanMateNoSelf
      if (Spec%GenderGiven) then
        ! @todo Is this a good target? It does not take into account relatedness between male parents and between female parents
        This%TargetCoancestryGenderMate  = Spec%TargetCoancestryRate + (1.0d0 - Spec%TargetCoancestryRate) * This%CurrentCoancestryGenderMate
      end if

      if (LogStdoutInternal) then
        write(STDOUT, "(a)") " "
        write(STDOUT, "(a)") "  - coancestry among/between individuals"
        write(STDOUT, "(a)") "                     Among     Between"
        write(STDOUT, "(a)") "    - n:       "//trim( Int2Char(This%CoancestryStat%All%n,    fmt=FMTINT2CHAR)) //" "//trim( Int2Char(This%CoancestryStat%OffDiag%n,      fmt=FMTINT2CHAR))
        write(STDOUT, "(a)") "    - average: "//trim(Real2Char(This%CoancestryStat%All%Mean, fmt=FMTREAL2CHAR))//" "//trim(Real2Char(This%CoancestryStat%OffDiag%Mean,   fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(This%CoancestryStat%All%SD,   fmt=FMTREAL2CHAR))//" "//trim(Real2Char(This%CoancestryStat%OffDiag%SD,     fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(This%CoancestryStat%All%Min,  fmt=FMTREAL2CHAR))//" "//trim(Real2Char(This%CoancestryStat%OffDiag%Min,    fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(This%CoancestryStat%All%Max,  fmt=FMTREAL2CHAR))//" "//trim(Real2Char(This%CoancestryStat%OffDiag%Max,    fmt=FMTREAL2CHAR))
        ! @todo: how does this interplay with the different modes
        write(STDOUT, "(a)") "    - target:  "//trim(Real2Char(This%TargetCoancestryRanMate, fmt=FMTREAL2CHAR))//" "//trim(Real2Char(This%TargetCoancestryRanMateNoSelf, fmt=FMTREAL2CHAR))&
                                              //" (given rate of "//trim(Real2Char(Spec%TargetCoancestryRate, fmt=FMTREAL2CHAR))//")"
        write(STDOUT, "(a)") "    Among   = coancestry among   individuals (including self-coancestry) = expected inbreeding in their progeny under random mating, including selfing"
        write(STDOUT, "(a)") "    Between = coancestry between individuals (excluding self-coancestry) = expected inbreeding in their progeny under random mating, excluding selfing"

        if (Spec%GenderGiven) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") "  - coancestry between males and females"
          write(STDOUT, "(a)") "    (=expected inbreeding in their progeny under random mating between genders)"
          write(STDOUT, "(a)") "    - n:       "//trim( Int2Char(This%CoancestryStatGenderDiff%All%n,    fmt=FMTINT2CHAR))
          write(STDOUT, "(a)") "    - average: "//trim(Real2Char(This%CoancestryStatGenderDiff%All%Mean, fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(This%CoancestryStatGenderDiff%All%SD,   fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(This%CoancestryStatGenderDiff%All%Min,  fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(This%CoancestryStatGenderDiff%All%Max,  fmt=FMTREAL2CHAR))
          ! @todo: how does this interplay with the different modes
          ! @todo: once I figure our what the TargetCoancestryGenderMate should be
          ! write(STDOUT, "(a)") "    - target:  "//trim(Real2Char(This%TargetCoancestryGenderMate,        fmt=FMTREAL2CHAR))&
          !                                       //" (given rate of "//trim(Real2Char(Spec%TargetCoancestryRate, fmt=FMTREAL2CHAR))//")"

          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") "  - coancestry among/between males"
          write(STDOUT, "(a)") "                     Among     Between"
          write(STDOUT, "(a)") "    - n:       "//trim( Int2Char(This%CoancestryStatGender1%All%n,    fmt=FMTINT2CHAR)) //" "//trim( Int2Char(This%CoancestryStatGender1%OffDiag%n,    fmt=FMTINT2CHAR))
          write(STDOUT, "(a)") "    - average: "//trim(Real2Char(This%CoancestryStatGender1%All%Mean, fmt=FMTREAL2CHAR))//" "//trim(Real2Char(This%CoancestryStatGender1%OffDiag%Mean, fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(This%CoancestryStatGender1%All%SD,   fmt=FMTREAL2CHAR))//" "//trim(Real2Char(This%CoancestryStatGender1%OffDiag%SD,   fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(This%CoancestryStatGender1%All%Min,  fmt=FMTREAL2CHAR))//" "//trim(Real2Char(This%CoancestryStatGender1%OffDiag%Min,  fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(This%CoancestryStatGender1%All%Max,  fmt=FMTREAL2CHAR))//" "//trim(Real2Char(This%CoancestryStatGender1%OffDiag%Max,  fmt=FMTREAL2CHAR))

          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") "  - coancestry among/between females"
          write(STDOUT, "(a)") "                     Among     Between"
          write(STDOUT, "(a)") "    - n:       "//trim( Int2Char(This%CoancestryStatGender2%All%n,    fmt=FMTINT2CHAR)) //" "//trim( Int2Char(This%CoancestryStatGender2%OffDiag%n,    fmt=FMTINT2CHAR))
          write(STDOUT, "(a)") "    - average: "//trim(Real2Char(This%CoancestryStatGender2%All%Mean, fmt=FMTREAL2CHAR))//" "//trim(Real2Char(This%CoancestryStatGender2%OffDiag%Mean, fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(This%CoancestryStatGender2%All%SD,   fmt=FMTREAL2CHAR))//" "//trim(Real2Char(This%CoancestryStatGender2%OffDiag%SD,   fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(This%CoancestryStatGender2%All%Min,  fmt=FMTREAL2CHAR))//" "//trim(Real2Char(This%CoancestryStatGender2%OffDiag%Min,  fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(This%CoancestryStatGender2%All%Max,  fmt=FMTREAL2CHAR))//" "//trim(Real2Char(This%CoancestryStatGender2%OffDiag%Max,  fmt=FMTREAL2CHAR))
        end if
      end if

      ! Save means to a file
      open(newunit=CoancestrySummaryUnit, file="CoancestrySummary.txt", status="unknown")
      write(CoancestrySummaryUnit, "(a, f)") "Current (random mating),                 ",   This%CurrentCoancestryRanMate
      write(CoancestrySummaryUnit, "(a, f)") "Current (random mating, no selfing),     ",   This%CurrentCoancestryRanMateNoSelf
      ! @todo once I figure our what the TargetCoancestryGenderMate should be
      ! if (Spec%GenderGiven) then
      !   write(CoancestrySummaryUnit, "(a, f)") "Current (random mating between genders), ", This%CurrentCoancestryGenderMate
      ! end if
      write(CoancestrySummaryUnit, "(a, f)") "Target (random mating),                  ",   This%TargetCoancestryRanMate
      write(CoancestrySummaryUnit, "(a, f)") "Target (random mating, no selfing),      ",   This%TargetCoancestryRanMateNoSelf
      ! @todo once I figure our what the TargetCoancestryGenderMate should be
      ! if (Spec%GenderGiven) then
      !   write(CoancestrySummaryUnit, "(a, f)") "Target (random mating between genders),  ", This%TargetCoancestryGenderMate
      ! end if
      close(CoancestrySummaryUnit)

      ! --- Current inbreeding summary ---

      if (LogStdoutInternal) then
        write(STDOUT, "(a)") " "
        write(STDOUT, "(a)") " Current inbreeding summary (identity between the two genomes of an individual)"
      end if

      call This%Coancestry%Inbreeding(Out=This%Inbreeding, Nrm=.false.)

      ! Current
      This%InbreedingStat = DescStat(This%Inbreeding%Value(1:))
      This%CurrentInbreeding = This%InbreedingStat%Mean

      ! Obtain limit/target based on given rates
      ! F_t = DeltaF + (1 - DeltaF) * F_t-1
      This%TargetInbreeding = Spec%TargetInbreedingRate + (1.0d0 - Spec%TargetInbreedingRate) * This%CurrentInbreeding

      if (LogStdoutInternal) then
        write(STDOUT, "(a)") "  - n:       "//trim( Int2Char(This%InbreedingStat%n,    fmt=FMTINT2CHAR))
        write(STDOUT, "(a)") "  - average: "//trim(Real2Char(This%InbreedingStat%Mean, fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "  - st.dev.: "//trim(Real2Char(This%InbreedingStat%SD,   fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "  - minimum: "//trim(Real2Char(This%InbreedingStat%Min,  fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "  - maximum: "//trim(Real2Char(This%InbreedingStat%Max,  fmt=FMTREAL2CHAR))
        ! @todo: how does this interplay with the different modes
        write(STDOUT, "(a)") "  - target:  "//trim(Real2Char(This%TargetInbreeding,    fmt=FMTREAL2CHAR))&
                                            //" (given rate of "//trim(Real2Char(Spec%TargetInbreedingRate, fmt=FMTREAL2CHAR))//")"
      end if

      open(newunit=InbreedingSummaryUnit, file="InbreedingSummary.txt", status="unknown")
      write(InbreedingSummaryUnit, "(a, f)") "Current, ", This%CurrentInbreeding
      write(InbreedingSummaryUnit, "(a, f)") "Target,  ", This%TargetInbreeding
      close(InbreedingSummaryUnit)

      ! --- Current selection criterion summary ---

      if (Spec%SelCriterionGiven) then

        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") " Current selection criterion summary"
        end if

        This%SelCriterionStat = DescStat(This%SelCriterion)
        This%SelCriterionStand = (This%SelCriterion - This%SelCriterionStat%Mean) / This%SelCriterionStat%SD
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") "  - n:       "//trim( Int2Char(This%SelCriterionStat%n,    fmt=FMTINT2CHAR))
          write(STDOUT, "(a)") "  - average: "//trim(Real2Char(This%SelCriterionStat%Mean, fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "  - st.dev.: "//trim(Real2Char(This%SelCriterionStat%SD,   fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "  - minimum: "//trim(Real2Char(This%SelCriterionStat%Min,  fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "  - maximum: "//trim(Real2Char(This%SelCriterionStat%Max,  fmt=FMTREAL2CHAR))
        end if

        if (This%SelCriterionStat%SD .eq. 0.0) then
          write(STDERR, "(a)") " ERROR: There is no variation in selection criterion!"
          write(STDERR, "(a)") " "
          stop 1
        end if

        open(newunit=CriterionSummaryUnit, file="SelCriterionSummary.txt", status="unknown")
        write(CriterionSummaryUnit, "(a, f)") "Mean, ", This%SelCriterionStat%Mean
        close(CriterionSummaryUnit)

        if (Spec%PAGEPar) then
          ! must have the same scale as selection criterion!!!!
          This%SelCriterionPAGEStand = (This%SelCriterionPAGE - This%SelCriterionStat%Mean) / This%SelCriterionStat%SD
          ! only the PAGE bit of SelCriterion
          This%SelCriterionPAGE = This%SelCriterionPAGE - This%SelCriterion
          This%SelCriterionPAGEStand = This%SelCriterionPAGEStand - This%SelCriterionStand
          This%SelCriterionPAGEStat = DescStat(This%SelCriterionPAGE)
          if (LogStdoutInternal) then
            write(STDOUT, "(a)") " "
            write(STDOUT, "(a)") "Selection criterion increments with PAGE"
            write(STDOUT, "(a)") "  - n:       "//trim( Int2Char(This%SelCriterionPAGEStat%n,    fmt=FMTINT2CHAR))
            write(STDOUT, "(a)") "  - average: "//trim(Real2Char(This%SelCriterionPAGEStat%Mean, fmt=FMTREAL2CHAR))
            write(STDOUT, "(a)") "  - st.dev.: "//trim(Real2Char(This%SelCriterionPAGEStat%SD,   fmt=FMTREAL2CHAR))
            write(STDOUT, "(a)") "  - minimum: "//trim(Real2Char(This%SelCriterionPAGEStat%Min,  fmt=FMTREAL2CHAR))
            write(STDOUT, "(a)") "  - maximum: "//trim(Real2Char(This%SelCriterionPAGEStat%Max,  fmt=FMTREAL2CHAR))
          end if

          if (This%SelCriterionPAGEStat%SD .eq. 0.0) then
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
          write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(This%GenericIndCritStat(Crit)%SD,   fmt=FMTREAL2CHAR))
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
            write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(This%GenericMatCritStat(Crit)%All%SD,   fmt=FMTREAL2CHAR))
            write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(This%GenericMatCritStat(Crit)%All%Min,  fmt=FMTREAL2CHAR))
            write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(This%GenericMatCritStat(Crit)%All%Max,  fmt=FMTREAL2CHAR))
            write(GenericMatCritSummaryUnit, "(a, f)") "Mean criterion "//trim(Int2Char(Crit)), This%GenericMatCritStat(Crit)%All%Mean
          else
            if (Spec%SelfingAllowed) then
              This%GenericMatCritStat(Crit) = DescStatLowTriMatrix(This%GenericMatCrit(:, :, Crit))
              write(STDOUT, "(a)") "    - n:       "//trim( Int2Char(This%GenericMatCritStat(Crit)%All%n,    fmt=FMTINT2CHAR))
              write(STDOUT, "(a)") "    - average: "//trim(Real2Char(This%GenericMatCritStat(Crit)%All%Mean, fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(This%GenericMatCritStat(Crit)%All%SD,   fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(This%GenericMatCritStat(Crit)%All%Min,  fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(This%GenericMatCritStat(Crit)%All%Max,  fmt=FMTREAL2CHAR))
              write(GenericMatCritSummaryUnit, "(a, f)") "Mean criterion "//trim(Int2Char(Crit)), This%GenericMatCritStat(Crit)%All%Mean
            end if
              This%GenericMatCritStat(Crit) = DescStatLowTriMatrix(This%GenericMatCrit(:, :, Crit), Diag=.false.)
              write(STDOUT, "(a)") "    - n:       "//trim( Int2Char(This%GenericMatCritStat(Crit)%OffDiag%n,    fmt=FMTINT2CHAR))
              write(STDOUT, "(a)") "    - average: "//trim(Real2Char(This%GenericMatCritStat(Crit)%OffDiag%Mean, fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(This%GenericMatCritStat(Crit)%OffDiag%SD,   fmt=FMTREAL2CHAR))
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

      write(Unit, "(a)") ""
      write(Unit, "(a)") " --- Raw data ---"

      write(Unit, "(a)") ""
      write(Unit, "(a)") " Coancestry:"
      call This%Coancestry%Write(File=File)

      write(Unit, "(a)") ""
      write(Unit, "(a)") " Inbreeding:"
      call This%Inbreeding%Write(File=File)

      write(Unit, "(a)") ""
      write(Unit, "(a)") " Gender:"
      if (allocated(This%Gender)) then
        do Ind = 1, This%nInd
          write(Unit, "(a, i)") This%Coancestry%OriginalId(Ind), This%Gender(Ind)
        end do
      else
        write(Unit, "(a)") " not allocated"
      end if

      write(Unit, "(a)") ""
      write(Unit, "(a)") " SelCriterion:"
      if (allocated(This%SelCriterion)) then
        do Ind = 1, This%nInd
          write(Unit, "(a, 2f)") This%Coancestry%OriginalId(Ind), This%SelCriterion(Ind), This%SelCriterionStand(Ind)
        end do
      else
        write(Unit, "(a)") " not allocated"
      end if

      write(Unit, "(a)") ""
      write(Unit, "(a)") " SelCriterionPAGE:"
      if (allocated(This%SelCriterionPAGE)) then
        do Ind = 1, This%nInd
          write(Unit, "(a, f)") This%Coancestry%OriginalId(Ind), This%SelCriterionPAGE(Ind), This%SelCriterionPAGEStand(Ind)
        end do
      else
        write(Unit, "(a)") " not allocated"
      end if

      write(Unit, "(a)") ""
      write(Unit, "(a)") " GenericIndCrit:"
      if (allocated(This%GenericIndCrit)) then
        do Ind = 1, This%nInd
          write(Unit, "(a, f)") This%Coancestry%OriginalId(Ind), This%GenericIndCrit(:, Ind)
        end do
      else
        write(Unit, "(a)") " not allocated"
      end if

      write(Unit, "(a)") ""
      write(Unit, "(a)") " GenericMatCrit:"
      if (allocated(This%GenericMatCrit)) then
        write(Unit, "(a)") " @todo"
        ! do Mat = 1, This%nPotMat
        !   write(Unit, "(i, f)") Mat, This%GenericMatCrit(:, :)
        ! end do
      else
        write(Unit, "(a)") " not allocated"
      end if

      write(Unit, "(a)") ""
      write(Unit, "(a)") " --- Data summaries ---"

      write(Unit, "(a)") ""
      write(Unit, "(a)") " CoancestryStat:"
      write(Unit, *) This%CoancestryStat

      write(Unit, "(a)") ""
      write(Unit, "(a)") " CoancestryStatGender1:"
      if (allocated(This%Gender)) then
        write(Unit, *) This%CoancestryStatGender1
      else
        write(Unit, "(a)") " Gender not allocated"
      end if

      write(Unit, "(a)") ""
      write(Unit, "(a)") " CoancestryStatGender2:"
      if (allocated(This%Gender)) then
        write(Unit, *) This%CoancestryStatGender2
      else
        write(Unit, "(a)") " Gender not allocated"
      end if

      write(Unit, "(a)") ""
      write(Unit, "(a)") " CoancestryStatGenderDiff:"
      if (allocated(This%Gender)) then
        write(Unit, *) This%CoancestryStatGenderDiff
      else
        write(Unit, "(a)") " Gender not allocated"
      end if

      write(Unit, "(a)") ""
      write(Unit, "(a)") " InbreedingStat:"
      write(Unit, *) This%InbreedingStat

      write(Unit, "(a)") ""
      write(Unit, "(a)") " SelCriterionStat:"
      if (allocated(This%SelCriterion)) then
        write(Unit, *) This%SelCriterionStat
      else
        write(Unit, "(a)") " SelCriterion not allocated"
      end if

      write(Unit, "(a)") ""
      write(Unit, "(a)") " SelCriterionStat:"
      if (allocated(This%SelCriterionPAGE)) then
        write(Unit, *) This%SelCriterionPAGEStat
      else
        write(Unit, "(a)") " SelCriterionPAGE not allocated"
      end if

      write(Unit, "(a)") ""
      write(Unit, "(a)") " GenericIndCritStat:"
      if (allocated(This%GenericIndCrit)) then
        write(Unit, *) This%GenericIndCritStat
      else
        write(Unit, "(a)") " GenericIndCrit not allocated"
      end if

      write(Unit, "(a)") ""
      write(Unit, "(a)") " GenericMatCritStat:"
      if (allocated(This%GenericMatCrit)) then
        write(Unit, *) This%GenericMatCritStat
      else
        write(Unit, "(a)") " GenericMatCrit not allocated"
      end if

      write(Unit, "(a)") ""
      write(Unit, "(a)") " --- Derived data ---"

      write(Unit, "(a)") ""
      write(Unit, "(a)") " nInd:"
      write(Unit, "(i)") This%nInd

      write(Unit, "(a)") ""
      write(Unit, "(a)") " nPotMat:"
      write(Unit, "(i)") This%nPotMat

      write(Unit, "(a)") ""
      write(Unit, "(a)") " nPotPar1:"
      write(Unit, "(i)") This%nPotPar1

      write(Unit, "(a)") ""
      write(Unit, "(a)") " nPotPar2:"
      write(Unit, "(i)") This%nPotPar2

      write(Unit, "(a)") ""
      write(Unit, "(a)") " nMal:"
      if (allocated(This%Gender)) then
        write(Unit, "(i)") This%nMal
      else
        write(Unit, "(a)") " Gender not allocated"
      end if

      write(Unit, "(a)") ""
      write(Unit, "(a)") " nFem:"
      if (allocated(This%Gender)) then
        write(Unit, "(i)") This%nFem
      else
        write(Unit, "(a)") " Gender not allocated"
      end if

      write(Unit, "(a)") ""
      write(Unit, "(a)") " IdPotPar1:"
      write(Unit, "(i)") This%IdPotPar1

      write(Unit, "(a)") ""
      write(Unit, "(a)") " IdPotPar2:"
      if (allocated(This%Gender)) then
        write(Unit, "(i)") This%IdPotPar2
      else
        write(Unit, "(a)") " Gender not allocated"
      end if

      write(Unit, "(a)") ""
      write(Unit, "(a)") " IdPotParSeq:"
      if (allocated(This%Gender)) then
        write(Unit, "(i)") This%IdPotParSeq
      else
        write(Unit, "(a)") " Gender not allocated"
      end if

      write(Unit, "(a)") ""
      write(Unit, "(a)") " CurrentCoancestryRanMate:"
      write(Unit, "(f)") This%CurrentCoancestryRanMate

      write(Unit, "(a)") ""
      write(Unit, "(a)") " CurrentCoancestryRanMateNoSelf:"
      write(Unit, "(f)") This%CurrentCoancestryRanMateNoSelf

      write(Unit, "(a)") ""
      write(Unit, "(a)") " CurrentCoancestryGenderMate:"
      if (allocated(This%Gender)) then
        write(Unit, "(f)") This%CurrentCoancestryGenderMate
      else
        write(Unit, "(a)") " Gender not allocated"
      end if

      write(Unit, "(a)") ""
      write(Unit, "(a)") " CurrentInbreeding:"
      write(Unit, "(f)") This%CurrentInbreeding

      write(Unit, "(a)") ""
      write(Unit, "(a)") " TargetCoancestryRanMate:"
      write(Unit, "(f)") This%TargetCoancestryRanMate

      write(Unit, "(a)") ""
      write(Unit, "(a)") " TargetCoancestryRanMateNoSelf:"
      write(Unit, "(f)") This%TargetCoancestryRanMateNoSelf

      write(Unit, "(a)") ""
      write(Unit, "(a)") " TargetCoancestryGenderMate:"
      if (allocated(This%Gender)) then
        write(Unit, "(f)") This%TargetCoancestryGenderMate
      else
        write(Unit, "(a)") " Gender not allocated"
      end if

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

      ! Initialisation
      select type (Spec)
        class default
          ! @todo This will work with Fortran 2015 standard (or at least with new? ifort)
          ! error stop " ERROR: InitialiseAlphaMateSol works only with argument Spec being of type AlphaMateSpec!"
        class is (AlphaMateSpec)
          This%Objective = -huge(This%Objective)
          This%nParam = size(Chrom)
          allocate(This%Chrom(This%nParam))
          This%Chrom = Chrom
          This%Penalty = 0.0d0
          This%PenaltyCoancestryRate = 0.0d0
          This%PenaltyInbreedingRate = 0.0d0
          This%PenaltySelfing = 0.0d0
          This%PenaltyLimitPar1 = 0.0d0
          This%PenaltyLimitPar2 = 0.0d0
          This%PenaltyGenericIndCrit = 0.0d0
          This%PenaltyGenericMatCrit = 0.0d0
          This%SelCriterion = 0.0d0
          This%SelIntensity = 0.0d0
          This%MaxPct = IEEE_Value(x=This%MaxPct, class=IEEE_Quiet_NaN)
          This%FutureCoancestryRanMate = 0.0d0
          This%CoancestryRateRanMate = 0.0d0
          This%MinPct = IEEE_Value(x=This%MinPct, class=IEEE_Quiet_NaN)
          This%Degree = IEEE_Value(x=This%Degree, class=IEEE_Quiet_NaN)
          This%FutureInbreeding = 0.0d0
          This%InbreedingRate = 0.0d0
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
          ! @todo This will work with Fortran 2015 standard (or at least with new? ifort)
          ! error stop " ERROR: AssignAlphaMateSol works only with argument In being of type AlphaMateSol!"
        class is (AlphaMateSol)
          Out%Objective = In%Objective
          Out%nParam = In%nParam
          Out%Chrom = In%Chrom
          Out%Penalty = In%Penalty
          Out%PenaltyCoancestryRate = In%PenaltyCoancestryRate
          Out%PenaltyInbreedingRate = In%PenaltyInbreedingRate
          Out%PenaltySelfing = In%PenaltySelfing
          Out%PenaltyLimitPar1 = In%PenaltyLimitPar1
          Out%PenaltyLimitPar2 = In%PenaltyLimitPar2
          Out%PenaltyGenericIndCrit = In%PenaltyGenericIndCrit
          Out%PenaltyGenericMatCrit = In%PenaltyGenericMatCrit
          Out%SelCriterion = In%SelCriterion
          Out%SelIntensity = In%SelIntensity
          Out%MaxPct = In%MaxPct
          Out%FutureCoancestryRanMate = In%FutureCoancestryRanMate
          Out%CoancestryRateRanMate = In%CoancestryRateRanMate
          Out%MinPct = In%MinPct
          Out%Degree = In%Degree
          Out%FutureInbreeding = In%FutureInbreeding
          Out%InbreedingRate = In%InbreedingRate
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
          ! @todo This will work with Fortran 2015 standard (or at least with new? ifort)
          ! error stop " ERROR: UpdateMeanAlphaMateSol works only with argument Add being of type AlphaMateSol!"
        class is (AlphaMateSol)
          This%Objective                 = This%Objective                 * kR + Add%Objective                 / n
          ! This%nParam                    = This%nParam                    * kR + Add%nParam                    / n ! the same all the time
          ! This%Chrom                     = This%Chrom                     * kR + Add%Chrom                     / n ! hmm, do we really want to average over chromosomes?
          This%Penalty                   = This%Penalty                   * kR + Add%Penalty                   / n
          This%Penalty                   = This%PenaltyCoancestryRate     * kR + Add%PenaltyCoancestryRate     / n
          This%Penalty                   = This%PenaltyInbreedingRate     * kR + Add%PenaltyInbreedingRate     / n
          This%Penalty                   = This%PenaltySelfing            * kR + Add%PenaltySelfing            / n
          This%Penalty                   = This%PenaltyLimitPar1          * kR + Add%PenaltyLimitPar1          / n
          This%Penalty                   = This%PenaltyLimitPar2          * kR + Add%PenaltyLimitPar2          / n
          This%Penalty                   = This%PenaltyGenericIndCrit     * kR + Add%PenaltyGenericIndCrit     / n
          This%Penalty                   = This%PenaltyGenericMatCrit     * kR + Add%PenaltyGenericMatCrit     / n
          This%SelCriterion              = This%SelCriterion              * kR + Add%SelCriterion              / n
          This%SelIntensity              = This%SelIntensity              * kR + Add%SelIntensity              / n
          This%MaxPct                    = This%MaxPct                    * kR + Add%MaxPct                    / n
          This%FutureCoancestryRanMate   = This%FutureCoancestryRanMate   * kR + Add%FutureCoancestryRanMate   / n
          This%CoancestryRateRanMate     = This%CoancestryRateRanMate     * kR + Add%CoancestryRateRanMate     / n
          This%MinPct                    = This%MinPct                    * kR + Add%MinPct                    / n
          This%Degree                    = This%Degree                    * kR + Add%Degree                    / n
          This%FutureInbreeding          = This%FutureInbreeding          * kR + Add%FutureInbreeding          / n
          This%InbreedingRate            = This%InbreedingRate            * kR + Add%InbreedingRate            / n
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

      real(real64) :: TmpR, RanNum, Opposite, Adjacent
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
                  if (i .le. Data%nPotPar2) then ! "="" to capture incremented i+1 on the do loop exit
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
                This%SelIntensity = dot_product(This%xVec, Data%SelCriterionStand)
                if (Spec%PAGEPar) then
                  This%SelIntensity = This%SelIntensity + dot_product(This%xVec, Data%SelCriterionPAGEStand * This%GenomeEdit)
                end if
                This%SelCriterion = This%SelIntensity * Data%SelCriterionStat%SD + Data%SelCriterionStat%Mean
                This%MaxPct = (1.0d0 - This%SelIntensity / Spec%MaxSelIntensity) * 100.0d0
                if (Spec%ModeOpt .or. Spec%ModeMax) then
                  This%Objective = This%Objective + This%SelIntensity
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
              This%FutureCoancestryRanMate = dot_product(TmpVec(:, 1), This%xVec)

              ! dF = (F_t+1 - F_t) / (1 - F_t)
              This%CoancestryRateRanMate = (This%FutureCoancestryRanMate - Data%CurrentCoancestryRanMate) / (1.0d0 - Data%CurrentCoancestryRanMate)
              This%MinPct = (1.0d0 - This%CoancestryRateRanMate / Spec%MinCoancestryRate) * 100.0d0

              ! Degree
              ! This calculation ASSUMES unit circular shape of selection/coancestry frontier
              ! - given results from ModeMin and ModeMax we can map the solution on the unit circle,
              !   i.e., the circle center is at (MaxCoancestryRate, MinSelIntensity) and solution is
              !   at (SolCoancestryRate, SolSelIntensity), which is mapped to (Opposite, Adjacent)
              !   on unit circle with center at (1, 0) (MinMode is at (0,0 ) and MaxMode is at (1, 1))
              ! - then we can calculate degrees of the angle between the max-line (MaxCoancestryRate, MinSelIntensity)-(MaxCoancestryRate, MaxSelIntensity)
              !   and the sol-line (MaxCoancestryRate, MinSelIntensity)-(SolCoancestryRate, SolSelIntensity)
              ! - the min-line would be (MaxCoancestryRate, MinSelIntensity)-(MinCoancestryRate, MinSelIntensity)
              ! - the solution y-coordinate on the max-line is adjacent to the angle and
              !   the solution x-coordinate on the min-line is opposite to the angle (if we
              !   put the min-line parallely up) so we use the arctangent function
              Opposite = (Spec%MaxCoancestryRate - This%CoancestryRateRanMate) / (Spec%MaxCoancestryRate - Spec%MinCoancestryRate)
              Adjacent = (This%SelIntensity      - Spec%MinSelIntensity)       / (Spec%MaxSelIntensity   - Spec%MinSelIntensity)
              This%Degree = atan(Opposite / Adjacent)

!@todo ModeRan?
              TmpR = 0.0d0
              if      (Spec%ModeMin) then
                TmpR = Spec%TargetCoancestryRateWeight * This%CoancestryRateRanMate
              else if (Spec%ModeOpt) then
                TmpR = This%CoancestryRateRanMate - Spec%TargetCoancestryRate
                if (This%CoancestryRateRanMate .lt. Spec%TargetCoancestryRate) then
                  if (Spec%TargetCoancestryRateWeightBelow) then
                    TmpR = abs(TmpR)
                  else
                    TmpR = 0.0d0
                  end if
                end if
                TmpR = Spec%TargetCoancestryRateWeight * TmpR
              end if
              This%Objective = This%Objective + TmpR
              if (TmpR .lt. 0.0d0) then
                This%PenaltyCoancestryRate = This%PenaltyCoancestryRate + TmpR
                This%Penalty               = This%Penalty               + TmpR
              end if

              ! --- Progeny inbreeding (=inbreeding of a mating) ---

              TmpR = 0.0d0
              do j = 1, Spec%nMat
                ! Lower triangle to speedup lookup
                TmpMax = maxval(This%MatingPlan(:, j))
                TmpMin = minval(This%MatingPlan(:, j))
                TmpR = TmpR + Data%Coancestry%Value(TmpMax, TmpMin)
              end do
              ! TODO: different number of progeny per mating???
              This%FutureInbreeding = TmpR / Spec%nMat
              ! dF = (F_t+1 - F_t) / (1 - F_t)
              This%InbreedingRate = (This%FutureInbreeding - Data%CurrentInbreeding) / (1.0d0 - Data%CurrentInbreeding)
              TmpR = This%InbreedingRate - Spec%TargetInbreedingRate
              if (.not. Spec%TargetInbreedingRateWeightBelow .and. (This%InbreedingRate .lt. Spec%TargetInbreedingRate)) then
                TmpR = 0.0d0
              end if
              TmpR = Spec%TargetInbreedingRateWeight * TmpR
              ! ! We know the targeted rate of inbreeding so we can work with relative values,
              ! ! which makes the TargetInbreedingRateWeight generic for ~any scenario.
              ! TmpR = This%InbreedingRate / Spec%TargetInbreedingRate
              ! if (TmpR .gt. 1.0d0) then
              !   ! Rate of inbreeding for the solution is higher than the target
              !   TmpR = Spec%TargetInbreedingRateWeight * abs(1.0d0 - TmpR)
              ! else
              !   ! Rate of inbreeding for the solution is lower than the target
              !   if (Spec%TargetInbreedingRateWeightBelow) then
              !     TmpR = Spec%TargetInbreedingRateWeight * abs(1.0d0 - abs(TmpR)) ! the second abs is to handle negative inbreeding cases
              !   else
              !     TmpR = 0.0d0
              !   end if
              ! end if
              This%Objective = This%Objective + TmpR
              if (Spec%TargetInbreedingRateWeight .lt. 0.0d0) then
                This%PenaltyInbreedingRate = This%PenaltyInbreedingRate + TmpR
                This%Penalty               = This%Penalty               + TmpR
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
    ! TODO: Produce an output object @return object that holds SolMin, SolRan, SolOpt, SolFrontier etc?
    ! TODO: the best random mating?
    subroutine AlphaMateSearch(Spec, Data, LogStdout) ! not pure due to IO
      implicit none
      type(AlphaMateSpec), intent(inout) :: Spec      !< AlphaMateSpec holder (out because we set and reset some parameters for different search modes)
      type(AlphaMateData), intent(inout) :: Data      !< AlphaMateData holder (out because we set and reset some parameters for different search modes)
      logical, intent(in), optional      :: LogStdout !< Log process on stdout (default .false.)

      integer(int32) :: nParam, Point, FrontierUnit

      real(real64) :: HoldTargetDegree
      real(real64) :: HoldTargetCoancestryRanMate, HoldTargetCoancestryRate
      real(real64), allocatable :: InitEqual(:, :)

      logical :: LogStdoutInternal, HoldTargetDegreeWeightBelow, HoldTargetCoancestryRateWeightBelow
      logical :: HoldModeMin, HoldModeMax, HoldModeOpt, HoldModeRan, HoldModeFrontier

      character(len=FILELENGTH) :: LogFile, LogPopFile, ContribFile, MatingFile

      type(AlphaMateSol) :: SolMin, SolMax, SolMinPct, SolMaxPct, SolOpt, SolRan, SolFrontier

      HoldModeMin = .false.
      HoldModeMax = .false.
      HoldModeOpt = .false.
      HoldModeRan = .false.
      HoldModeFrontier = .false.

      if (present(LogStdout)) then
        LogStdoutInternal = LogStdout
      else
        LogStdoutInternal = .false.
      end if

      ! --- Number of parameters to optimise ---

      if (Spec%GenderGiven) then
        nParam = Data%nPotPar1 + Data%nPotPar2 + Spec%nMat
      else
        nParam = Data%nPotPar1 + Spec%nMat
      end if

      if (Spec%PAGEPar) then
        nParam = nParam + Data%nInd
      end if

      ! --- Minimum future coancestry/inbreeding ---

      if (Spec%ModeMin) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") " Optimise contributions for minimum future coancestry/inbreeding"
          write(STDOUT, "(a)") " "
        end if

        ! Setup
        ! Hold modes so that the Evaluate method is not confused what we are optimising for
        HoldModeMax = Spec%ModeMax
        Spec%ModeMax = .false.
        HoldModeOpt = Spec%ModeOpt
        Spec%ModeOpt = .false.
        HoldModeRan = Spec%ModeRan
        Spec%ModeRan = .false.
        HoldModeFrontier = Spec%ModeFrontier
        Spec%ModeFrontier = .false.

        ! We want to minimise rate of coancestry
        HoldTargetCoancestryRateWeightBelow = Spec%TargetCoancestryRateWeightBelow
        Spec%TargetCoancestryRateWeightBelow = .false.

        LogFile     = "OptimisationLogModeMin.txt"
        LogPopFile  = "OptimisationLogPopModeMin.txt"
        ContribFile = "ContributionsModeMin.txt"
        MatingFile  = "MatingPlanModeMin.txt"

        call SolMin%SetupColNamesAndFormats(Spec=Spec)

        ! Search
        ! @todo Can we do this in a better way where we take "structure" of Chrom into account?
        allocate(InitEqual(nParam, nint(Spec%EvolAlgNSol * 0.1)))
        InitEqual = 1.0d0 ! A couple of solutions that would give equal contributions to everybody

        if (trim(Spec%EvolAlg) .eq. "DE") then
          call DifferentialEvolution(Spec=Spec, Data=Data, nParam=nParam, nSol=Spec%EvolAlgNSol, Init=InitEqual, &
            nIter=Spec%EvolAlgNIter, nIterBurnIn=Spec%DiffEvolNIterBurnIn, nIterStop=Spec%EvolAlgNIterStop, StopTolerance=Spec%EvolAlgStopTol, nIterPrint=Spec%EvolAlgNIterPrint, &
            LogStdout=LogStdoutInternal, LogFile=LogFile, LogPop=Spec%EvolAlgLogPop, LogPopFile=LogPopFile, &
            CRBurnIn=Spec%DiffEvolParamCrBurnIn, CRLate=Spec%DiffEvolParamCr, FBase=Spec%DiffEvolParamFBase, FHigh1=Spec%DiffEvolParamFHigh1, FHigh2=Spec%DiffEvolParamFHigh2, &
            BestSol=SolMin)
        end if

        deallocate(InitEqual)

        call SolMin%WriteContributions(Data, ContribFile)
        call SolMin%WriteMatingPlan(Data, MatingFile)

        ! Save for later use
        Spec%MinCoancestry     = SolMin%FutureCoancestryRanMate
        Spec%MinCoancestryRate = SolMin%CoancestryRateRanMate
        Spec%MinSelCriterion   = SolMin%SelCriterion
        Spec%MinSelIntensity   = SolMin%SelIntensity

        ! Reset
        Spec%ModeMax = HoldModeMax
        Spec%ModeOpt = HoldModeOpt
        Spec%ModeRan = HoldModeRan
        Spec%ModeFrontier = HoldModeFrontier
        Spec%TargetCoancestryRateWeightBelow = HoldTargetCoancestryRateWeightBelow
      end if

      ! --- Maximum future selection criterion ---

      if (Spec%ModeMax) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") " Optimise contributions for maximum future selection criterion"
          write(STDOUT, "(a)") " "
        end if

        ! Setup
        ! Hold modes so that the Evaluate method is not confused what we are optimising for
        HoldModeMin = Spec%ModeMin
        Spec%ModeMin = .false.
        HoldModeOpt = Spec%ModeOpt
        Spec%ModeOpt = .false.
        HoldModeRan = Spec%ModeRan
        Spec%ModeRan = .false.
        HoldModeFrontier = Spec%ModeFrontier
        Spec%ModeFrontier = .false.

        LogFile     = "OptimisationLogModeMax.txt"
        LogPopFile  = "OptimisationLogPopModeMax.txt"
        ContribFile = "ContributionsModeMax.txt"
        MatingFile  = "MatingPlanModeMax.txt"

        call SolMax%SetupColNamesAndFormats(Spec=Spec)

        ! Search
        ! @todo add some clever initial values, say:
        !       - equal contributions for top 2/3 or 1/2 of BV distribution,
        !       - decreasing contributions with decreasing value
        !       - SDP solution, ...?
        if (trim(Spec%EvolAlg) .eq. "DE") then
          call DifferentialEvolution(Spec=Spec, Data=Data, nParam=nParam, nSol=Spec%EvolAlgNSol, nIter=Spec%EvolAlgNIter, nIterBurnIn=Spec%DiffEvolNIterBurnIn, &
            nIterStop=Spec%EvolAlgNIterStop, StopTolerance=Spec%EvolAlgStopTol, nIterPrint=Spec%EvolAlgNIterPrint,&
            LogStdout=LogStdoutInternal, LogFile=LogFile, LogPop=Spec%EvolAlgLogPop, LogPopFile=LogPopFile, &
            CRBurnIn=Spec%DiffEvolParamCrBurnIn, CRLate=Spec%DiffEvolParamCr, FBase=Spec%DiffEvolParamFBase, FHigh1=Spec%DiffEvolParamFHigh1, FHigh2=Spec%DiffEvolParamFHigh2, &
            BestSol=SolMax)
        end if

        ! Output
        call SolMax%WriteContributions(Data, ContribFile)
        call SolMax%WriteMatingPlan(Data, MatingFile)

        ! Save for later use
        Spec%MaxCoancestry     = SolMax%FutureCoancestryRanMate
        Spec%MaxCoancestryRate = SolMax%CoancestryRateRanMate
        Spec%MaxSelCriterion   = SolMax%SelCriterion
        Spec%MaxSelIntensity   = SolMax%SelIntensity

        ! Reset
        Spec%ModeMin = HoldModeMin
        Spec%ModeOpt = HoldModeOpt
        Spec%ModeRan = HoldModeRan
        Spec%ModeFrontier = HoldModeFrontier
      end if

      ! --- Percentage above minimum future coancestry/inbreeding ---

      if (Spec%ModeMinPct) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") " Optimise contributions for a 'percentage' above minimum future coancestry/inbreeding"
          write(STDOUT, "(a)") " "
        end if

        ! Setup
        ! Hold modes so that the Evaluate method is not confused what we are optimising for
        HoldModeMin = Spec%ModeMin
        Spec%ModeMin = .false.
        HoldModeMax = Spec%ModeMax
        Spec%ModeMax = .false.
        HoldModeOpt = Spec%ModeOpt
        Spec%ModeOpt = .true. ! we will balance future sel. criterion and coancestry
        HoldModeRan = Spec%ModeRan
        Spec%ModeRan = .false.
        HoldModeFrontier = Spec%ModeFrontier
        Spec%ModeFrontier = .false.

        ! Hold old results
        HoldTargetCoancestryRanMate = Data%TargetCoancestryRanMate
        HoldTargetCoancestryRate = Spec%TargetCoancestryRate

        ! This is the rate we target
        ! How do we obtain it?
        ! - given MinPct we simply compute the target coancestry rate as the
        !   minimum possible coancestry rate plus the percentage of difference
        !   between the maximum and the minimum possible coancestry rates
        Spec%TargetCoancestryRate = SolMin%CoancestryRateRanMate + Spec%MinPct / 100 * (SolMax%CoancestryRateRanMate - SolMin%CoancestryRateRanMate)
        ! F_t = DeltaF + (1 - DeltaF) * F_t-1
        Data%TargetCoancestryRanMate = Spec%TargetCoancestryRate + (1.0d0 - Spec%TargetCoancestryRate) * Data%CurrentCoancestryRanMate
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " Minimum  rate of coancestry: "//trim(Real2Char(SolMin%CoancestryRateRanMate, fmt=FMTREAL2CHAR))//&
                               " (=coancestry "//trim(Real2Char(SolMin%FutureCoancestryRanMate, fmt=FMTREAL2CHAR))//")"
          write(STDOUT, "(a)") " Maximum  rate of coancestry: "//trim(Real2Char(SolMax%CoancestryRateRanMate, fmt=FMTREAL2CHAR))//&
                               " (=coancestry "//trim(Real2Char(SolMax%FutureCoancestryRanMate, fmt=FMTREAL2CHAR))//")"
          write(STDOUT, "(a)") " Percentage above minimum:     "//trim(Real2Char(Spec%MinPct, fmt="(f7.2)"))
          write(STDOUT, "(a)") " Targeted rate of coancestry: "//trim(Real2Char(Spec%TargetCoancestryRate, fmt=FMTREAL2CHAR))//&
                               " (=coancestry "//trim(Real2Char(Data%TargetCoancestryRanMate, fmt=FMTREAL2CHAR))//")"
          write(STDOUT, "(a)") ""
        end if

        LogFile     = "OptimisationLogModeMinPct.txt"
        LogPopFile  = "OptimisationLogPopModeMinPct.txt"
        ContribFile = "ContributionsModeMinPct.txt"
        MatingFile  = "MatingPlanModeMinPct.txt"

        call SolMinPct%SetupColNamesAndFormats(Spec=Spec)

        ! Search
        ! @todo add some clever initial values, say:
        !       - equal contributions for top 2/3 or 1/2 of BV distribution,
        !       - decreasing contributions with decreasing value
        !       - SDP solution, ...?
        if (trim(Spec%EvolAlg) .eq. "DE") then
          call DifferentialEvolution(Spec=Spec, Data=Data, nParam=nParam, nSol=Spec%EvolAlgNSol, nIter=Spec%EvolAlgNIter, nIterBurnIn=Spec%DiffEvolNIterBurnIn, &
            nIterStop=Spec%EvolAlgNIterStop, StopTolerance=Spec%EvolAlgStopTol, nIterPrint=Spec%EvolAlgNIterPrint,&
            LogStdout=LogStdoutInternal, LogFile=LogFile, LogPop=Spec%EvolAlgLogPop, LogPopFile=LogPopFile, &
            CRBurnIn=Spec%DiffEvolParamCrBurnIn, CRLate=Spec%DiffEvolParamCr, FBase=Spec%DiffEvolParamFBase, FHigh1=Spec%DiffEvolParamFHigh1, FHigh2=Spec%DiffEvolParamFHigh2, &
            BestSol=SolMinPct)
        end if

        ! Output
        call SolMinPct%WriteContributions(Data, ContribFile)
        call SolMinPct%WriteMatingPlan(Data, MatingFile)

        ! Reset
        Spec%ModeMin = HoldModeMin
        Spec%ModeMax = HoldModeMax
        Spec%ModeOpt = HoldModeOpt
        Spec%ModeRan = HoldModeRan
        Spec%ModeFrontier = HoldModeFrontier
        Spec%TargetCoancestryRateWeightBelow = HoldTargetCoancestryRateWeightBelow
        Data%TargetCoancestryRanMate = HoldTargetCoancestryRanMate
      end if

      ! --- Percentage below maximum future selection criterion ---

      if (Spec%ModeMaxPct) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") " Optimise contributions for a 'percentage' below maximum future selection criterion"
          write(STDOUT, "(a)") " "
        end if

        ! Setup
        ! Hold modes so that the Evaluate method is not confused what we are optimising for
        HoldModeMin = Spec%ModeMin
        Spec%ModeMin = .false.
        HoldModeMax = Spec%ModeMax
        Spec%ModeMax = .false.
        HoldModeOpt = Spec%ModeOpt
        Spec%ModeOpt = .true. ! we will balance future sel. criterion and coancestry
        HoldModeRan = Spec%ModeRan
        Spec%ModeRan = .false.
        HoldModeFrontier = Spec%ModeFrontier
        Spec%ModeFrontier = .false.

        ! Hold old results
        HoldTargetCoancestryRanMate = Data%TargetCoancestryRanMate
        HoldTargetCoancestryRate = Spec%TargetCoancestryRate

        ! These are the criterion and rate we target
        ! How do we obtain them?
        ! - given MaxPct we compute targeted selection criterion (TargetSelCritRat)
        ! - we need to find target coancestry rate (TargetCoaRate) that will give us TargetSelCritRat
        ! - we have
        !   - distance from MinSelCrit    to MaxSelCrit    of 1
        !   - distance from MinSelCrit    to TargetSelCrit of p
        !   - distance from TargetSelCrit to MaxSelCrit    of 1 - p (MaxPct)
        !   - distance from MinCoaRate    to MaxCoaRate    of 1
        !   - distance from MaxCoaRate    to TargetCoaRate of x
        !   - distance from MinCoaRate    to TargetCoaRate of y = 1 - x (we need to find this)
        !   - a triangle with sides:
        !       x = (TargetCoaRate, MinSelCrit)    to (MaxCoaRate,    MinSelCrit)
        !       p = (TargetCoaRate, MinSelCrit)    to (TargetCoaRate, TargetSelCrit)
        !       r = (TargetCoaRate, TargetSelCrit) to (MaxCoaRate,    MinSelCrit)
        !  - ASSUMING unit circular shape of selection/coancestry frontier the radius (r) is 1
        !  - using Pythagorean theorem we have
        !      r^2 = x^2 + p^2
        !      x^2 = r^2 - p^2
        !          = 1   - p^2
        !      x   = sqrt(1 - p^2)
        !      y = 1 - x
        !        = 1 - sqrt(1 - p^2)
        block
          real(real64) :: TargetSelCritRat, TargetCoaRat
          TargetSelCritRat = 1.0d0 - Spec%MaxPct / 100
          TargetCoaRat = 1.0d0 - sqrt(1.0d0 - (TargetSelCritRat * TargetSelCritRat))
          Spec%TargetCoancestryRate = SolMin%CoancestryRateRanMate + TargetCoaRat * (SolMax%CoancestryRateRanMate - SolMin%CoancestryRateRanMate)
          ! F_t = DeltaF + (1 - DeltaF) * F_t-1
          Data%TargetCoancestryRanMate = Spec%TargetCoancestryRate + (1.0d0 - Spec%TargetCoancestryRate) * Data%CurrentCoancestryRanMate
          if (LogStdoutInternal) then
            write(STDOUT, "(a)") " Minimum  selection intensity: "//trim(Real2Char(SolMin%SelIntensity, fmt=FMTREAL2CHAR))//&
                                " (=selection criterion"//trim(Real2Char(SolMin%SelCriterion, fmt=FMTREAL2CHAR))//")"
            write(STDOUT, "(a)") " Maximum  selection intensity: "//trim(Real2Char(SolMax%SelIntensity, fmt=FMTREAL2CHAR))//&
                                " (=selection criterion"//trim(Real2Char(SolMax%SelCriterion, fmt=FMTREAL2CHAR))//")"
            write(STDOUT, "(a)") " Percentage below maximum:      "//trim(Real2Char(Spec%MaxPct, fmt="(f7.2)"))
            write(STDOUT, "(a)") " Targeted selection intensity: "//trim(Real2Char(SolMax%SelIntensity * TargetSelCritRat, fmt=FMTREAL2CHAR))//&
                                " (=selection criterion"//trim(Real2Char(SolMax%SelCriterion * TargetSelCritRat, fmt=FMTREAL2CHAR))//")"
            write(STDOUT, "(a)") ""

            write(STDOUT, "(a)") " Minimum  rate of coancestry: "//trim(Real2Char(SolMin%CoancestryRateRanMate, fmt=FMTREAL2CHAR))//&
                                " (=coancestry "//trim(Real2Char(SolMin%FutureCoancestryRanMate, fmt=FMTREAL2CHAR))//")"
            write(STDOUT, "(a)") " Maximum  rate of coancestry: "//trim(Real2Char(SolMax%CoancestryRateRanMate, fmt=FMTREAL2CHAR))//&
                                " (=coancestry "//trim(Real2Char(SolMax%FutureCoancestryRanMate, fmt=FMTREAL2CHAR))//")"
            write(STDOUT, "(a)") " Percentage above minimum:     "//trim(Real2Char(TargetCoaRat * 100.0d0, fmt="(f7.2)"))
            write(STDOUT, "(a)") " Targeted rate of coancestry: "//trim(Real2Char(Spec%TargetCoancestryRate, fmt=FMTREAL2CHAR))//&
                                " (=coancestry "//trim(Real2Char(Data%TargetCoancestryRanMate, fmt=FMTREAL2CHAR))//")"
            write(STDOUT, "(a)") ""
          end if
        end block

        LogFile     = "OptimisationLogModeMaxPct.txt"
        LogPopFile  = "OptimisationLogPopModeMaxPct.txt"
        ContribFile = "ContributionsModeMaxPct.txt"
        MatingFile  = "MatingPlanModeMaxPct.txt"

        call SolMaxPct%SetupColNamesAndFormats(Spec=Spec)

        ! Search
        ! @todo add some clever initial values, say:
        !       - equal contributions for top 2/3 or 1/2 of BV distribution,
        !       - decreasing contributions with decreasing value
        !       - SDP solution, ...?
        if (trim(Spec%EvolAlg) .eq. "DE") then
          call DifferentialEvolution(Spec=Spec, Data=Data, nParam=nParam, nSol=Spec%EvolAlgNSol, nIter=Spec%EvolAlgNIter, nIterBurnIn=Spec%DiffEvolNIterBurnIn, &
            nIterStop=Spec%EvolAlgNIterStop, StopTolerance=Spec%EvolAlgStopTol, nIterPrint=Spec%EvolAlgNIterPrint,&
            LogStdout=LogStdoutInternal, LogFile=LogFile, LogPop=Spec%EvolAlgLogPop, LogPopFile=LogPopFile, &
            CRBurnIn=Spec%DiffEvolParamCrBurnIn, CRLate=Spec%DiffEvolParamCr, FBase=Spec%DiffEvolParamFBase, FHigh1=Spec%DiffEvolParamFHigh1, FHigh2=Spec%DiffEvolParamFHigh2, &
            BestSol=SolMaxPct)
        end if

        ! Output
        call SolMaxPct%WriteContributions(Data, ContribFile)
        call SolMaxPct%WriteMatingPlan(Data, MatingFile)

        ! Reset
        Spec%ModeMin = HoldModeMin
        Spec%ModeMax = HoldModeMax
        Spec%ModeOpt = HoldModeOpt
        Spec%ModeRan = HoldModeRan
        Spec%ModeFrontier = HoldModeFrontier
        Spec%TargetCoancestryRateWeightBelow = HoldTargetCoancestryRateWeightBelow
        Data%TargetCoancestryRanMate = HoldTargetCoancestryRanMate
      end if

      ! --- Maximum future selection criterion with constraint on coancestry/inbreeding ---

      if (Spec%ModeOpt) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") " Optimise contributions for maximum future selection criterion with constraint on coancestry/inbreeding"
          write(STDOUT, "(a)") " "
        end if

        ! Setup
        ! Hold modes so that the Evaluate method is not confused what we are optimising for
        HoldModeMin = Spec%ModeMin
        Spec%ModeMin = .false.
        HoldModeMax = Spec%ModeMax
        Spec%ModeMax = .false.
        HoldModeRan = Spec%ModeRan
        Spec%ModeRan = .false.
        HoldModeFrontier = Spec%ModeFrontier
        Spec%ModeFrontier = .false.

        LogFile     = "OptimisationLogModeOpt.txt"
        LogPopFile  = "OptimisationLogPopModeOpt.txt"
        ContribFile = "ContributionsModeOpt.txt"
        MatingFile  = "MatingPlanModeOpt.txt"

        call SolOpt%SetupColNamesAndFormats(Spec=Spec)

        ! Search
        ! @todo add some clever initial values, say:
        !       - equal contributions for top 2/3 or 1/2 of BV distribution,
        !       - decreasing contributions with decreasing value
        !       - SDP solution, ...?
        if (trim(Spec%EvolAlg) .eq. "DE") then
          call DifferentialEvolution(Spec=Spec, Data=Data, nParam=nParam, nSol=Spec%EvolAlgNSol, nIter=Spec%EvolAlgNIter, nIterBurnIn=Spec%DiffEvolNIterBurnIn, &
            nIterStop=Spec%EvolAlgNIterStop, StopTolerance=Spec%EvolAlgStopTol, nIterPrint=Spec%EvolAlgNIterPrint,&
            LogStdout=LogStdoutInternal, LogFile=LogFile, LogPop=Spec%EvolAlgLogPop, LogPopFile=LogPopFile, &
            CRBurnIn=Spec%DiffEvolParamCrBurnIn, CRLate=Spec%DiffEvolParamCr, FBase=Spec%DiffEvolParamFBase, FHigh1=Spec%DiffEvolParamFHigh1, FHigh2=Spec%DiffEvolParamFHigh2, &
            BestSol=SolOpt)
        end if

        ! Output
        call SolOpt%WriteContributions(Data, ContribFile)
        call SolOpt%WriteMatingPlan(Data, MatingFile)

        ! Reset
        Spec%ModeMin = HoldModeMin
        Spec%ModeMax = HoldModeMax
        Spec%ModeRan = HoldModeRan
        Spec%ModeFrontier = HoldModeFrontier
      end if

      ! --- Random mating ---

      ! @todo: is this mode any good?
      if (Spec%ModeRan) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") " Evaluate random mating"
          write(STDOUT, "(a)") " "
        end if

        ! Setup
        ! Hold modes so that the Evaluate method is not confused what we are optimising for
        HoldModeMin = Spec%ModeMin
        Spec%ModeMin = .false.
        HoldModeMax = Spec%ModeMax
        Spec%ModeMax = .false.
        HoldModeOpt = Spec%ModeOpt
        Spec%ModeOpt = .false.
        HoldModeFrontier = Spec%ModeFrontier
        Spec%ModeFrontier = .false.

        LogFile = "OptimisationLogModeRan.txt"
        ! @todo: other reports from here - at least the best random solution?

        call SolRan%SetupColNamesAndFormats(Spec=Spec)

        ! Search
        allocate(InitEqual(nParam, nint(Spec%EvolAlgNSol * 0.1)))
        InitEqual = 1.0d0 ! A couple of solutions that would give equal contributions for everybody

        call RandomSearch(Mode="avg", Spec=Spec, Data=Data, nParam=nParam, Init=InitEqual, &
          nSamp=Spec%EvolAlgNSol*Spec%EvolAlgNIter*Spec%RanAlgStricter, nSampStop=Spec%EvolAlgNIterStop*Spec%RanAlgStricter, &
          StopTolerance=Spec%EvolAlgStopTol/Spec%RanAlgStricter, nSampPrint=Spec%EvolAlgNIterPrint, &
          LogStdout=LogStdoutInternal, LogFile=LogFile, BestSol=SolRan)

        deallocate(InitEqual)

        ! Reset
        Spec%ModeMin = HoldModeMin
        Spec%ModeMax = HoldModeMax
        Spec%ModeOpt = HoldModeOpt
        Spec%ModeFrontier = HoldModeFrontier
      end if

      ! --- Evaluate frontier ---
      ! this mode needs to be after ModeMin and ModeMax!!!

      if (Spec%ModeFrontier) then
        if (LogStdoutInternal) then
          write(STDOUT, "(a)") " "
          write(STDOUT, "(a)") " Evaluate frontier"
        end if

        open(newunit=FrontierUnit, file="Frontier.txt", status="unknown")

        ! Setup
        call SolFrontier%SetupColNamesAndFormats(Spec=Spec) ! so we can log any previous results
        call SolFrontier%LogHead(LogUnit=FrontierUnit, String="ModeOrPoint", StringNum=15)

        ! Add any previous results to frontier
        if (Spec%ModeMin) then
          call SolMin%Log(FrontierUnit, Iteration=-1, AcceptRate=-1.0d0, String="ModeMin", StringNum=15)
        end if
        if (Spec%ModeMax) then
          call SolMax%Log(FrontierUnit, Iteration=-1, AcceptRate=-1.0d0, String="ModeMax", StringNum=15)
        end if
        if (Spec%ModeOpt) then
          call SolOpt%Log(FrontierUnit, Iteration=-1, AcceptRate=-1.0d0, String="ModeOpt", StringNum=15)
        end if
        if (Spec%ModeRan) then
          call SolRan%Log(FrontierUnit, Iteration=-1, AcceptRate=-1.0d0, String="ModeRan", StringNum=15)
        end if

        ! Hold modes so that the Evaluate method is not confused what we are optimising for
        HoldModeMin = Spec%ModeMin
        Spec%ModeMin = .false.
        HoldModeMax = Spec%ModeMax
        Spec%ModeMax = .false.
        HoldModeOpt = Spec%ModeOpt
        Spec%ModeOpt = .true. ! we will balance future sel. criterion and coancestry
        HoldModeRan = Spec%ModeRan
        Spec%ModeRan = .false.

        ! We want to target exact rates of coancestry
        HoldTargetCoancestryRateWeightBelow = Spec%TargetCoancestryRateWeightBelow
        Spec%TargetCoancestryRateWeightBelow = .true.

        ! Hold old results
        HoldTargetCoancestryRanMate = Data%TargetCoancestryRanMate
        HoldTargetCoancestryRate = Spec%TargetCoancestryRate

        ! Evaluate
        do Point = 1, Spec%nFrontierPoints
          ! These are the rates
          Spec%TargetCoancestryRate = Spec%TargetCoancestryRateFrontier(Point)
          ! F_t = DeltaF + (1 - DeltaF) * F_t-1
          Data%TargetCoancestryRanMate = Spec%TargetCoancestryRate + (1.0d0 - Spec%TargetCoancestryRate) * Data%CurrentCoancestryRanMate
          if (LogStdoutInternal) then
            write(STDOUT, "(a)") ""
            write(STDOUT, "(a)") " Point "//trim(Int2Char(Point))//" out of "//trim(Int2Char(Spec%nFrontierPoints))
            write(STDOUT, "(a)") " Targeted rate of coancestry "//trim(Real2Char(Spec%TargetCoancestryRate, fmt=FMTREAL2CHAR))//&
                                " (=coancestry "//trim(Real2Char(Data%TargetCoancestryRanMate, fmt=FMTREAL2CHAR))//")"
            write(STDOUT, "(a)") ""
          end if

          LogFile     = "OptimisationLogModeFrontier"//trim(Int2Char(Point))//".txt"
          LogPopFile  = "OptimisationLogPopModeFrontier"//trim(Int2Char(Point))//".txt"
          ContribFile = "ContributionsModeFrontier"//trim(Int2Char(Point))//".txt"
          MatingFile  = "MatingPlanModeFrontier"//trim(Int2Char(Point))//".txt"

          call SolFrontier%SetupColNamesAndFormats(Spec=Spec) ! call again as InitialiseAlphaMateSol and AssignAlphaMateSol "nullify" the  above SetupColNamesAndFormats call (ugly, but works ...)

          ! Search
          if (trim(Spec%EvolAlg) .eq. "DE") then
            call DifferentialEvolution(Spec=Spec, Data=Data, nParam=nParam, nSol=Spec%EvolAlgNSol, nIter=Spec%EvolAlgNIter, nIterBurnIn=Spec%DiffEvolNIterBurnIn, &
              nIterStop=Spec%EvolAlgNIterStop, StopTolerance=Spec%EvolAlgStopTol, nIterPrint=Spec%EvolAlgNIterPrint,&
              LogStdout=LogStdoutInternal, LogFile=LogFile, LogPop=Spec%EvolAlgLogPop, LogPopFile=LogPopFile, &
              CRBurnIn=Spec%DiffEvolParamCrBurnIn, CRLate=Spec%DiffEvolParamCr, FBase=Spec%DiffEvolParamFBase, FHigh1=Spec%DiffEvolParamFHigh1, FHigh2=Spec%DiffEvolParamFHigh2, &
              BestSol=SolFrontier)
          end if

          ! Output
          call SolFrontier%Log(FrontierUnit, Iteration=-1, AcceptRate=-1.0d0, String=trim("ModeFrontier"//trim(Int2Char(Point))), StringNum=15)
          call SolFrontier%WriteContributions(Data, ContribFile)
          call SolFrontier%WriteMatingPlan(Data, MatingFile)

          if ((Spec%TargetCoancestryRate - SolFrontier%CoancestryRateRanMate) .gt. 0.01d0) then
            if (LogStdoutInternal) then
              write(STDOUT, "(a)") ""
              write(STDOUT, "(a)") "NOTE: Could not achieve the rate of coancestry "//trim(Real2Char(Spec%TargetCoancestryRate, fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "NOTE: Stopping frontier evaluation."
              write(STDOUT, "(a)") ""
            end if
            exit
          end if
        end do

        ! Reset
        Spec%ModeMin = HoldModeMin
        Spec%ModeMax = HoldModeMax
        Spec%ModeOpt = HoldModeOpt
        Spec%ModeRan = HoldModeRan
        Spec%TargetCoancestryRate = HoldTargetCoancestryRate
        Data%TargetCoancestryRanMate = HoldTargetCoancestryRanMate
        Spec%TargetCoancestryRateWeightBelow = HoldTargetCoancestryRateWeightBelow

        close(FrontierUnit)

      end if
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

      write(Unit, *) "Objective: ", This%Objective
      write(Unit, *) "nParam: ", This%nParam
      write(Unit, *) "Chrom: ", This%Chrom
      write(Unit, *) "Penalty: ", This%Penalty
      write(Unit, *) "PenaltyCoancestryRate: ", This%PenaltyCoancestryRate
      write(Unit, *) "PenaltyInbreedingRate: ", This%PenaltyInbreedingRate
      write(Unit, *) "PenaltySelfing: ", This%PenaltySelfing
      write(Unit, *) "PenaltyLimitPar1: ", This%PenaltyLimitPar1
      write(Unit, *) "PenaltyLimitPar2: ", This%PenaltyLimitPar2
      write(Unit, *) "PenaltyGenericIndCrit: ", This%PenaltyGenericIndCrit
      write(Unit, *) "PenaltyGenericMatCrit: ", This%PenaltyGenericMatCrit
      write(Unit, *) "SelCriterion: ", This%SelCriterion
      write(Unit, *) "SelIntensity: ", This%SelIntensity
      write(Unit, *) "FutureCoancestryRanMate: ", This%FutureCoancestryRanMate
      write(Unit, *) "CoancestryRateRanMate: ", This%CoancestryRateRanMate
      write(Unit, *) "FutureInbreeding: ", This%FutureInbreeding
      write(Unit, *) "InbreedingRate: ", This%InbreedingRate
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

      nCol = 10
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
      This%ColnameLogUnit(2)  = "            AcceptRate"
      This%ColnameLogUnit(3)  = "             Objective"
      This%ColnameLogUnit(4)  = "             Penalties"
      This%ColnameLogUnit(5)  = "          SelCriterion"
      This%ColnameLogUnit(6)  = "          SelIntensity"
      This%ColnameLogUnit(7)  = "            Coancestry"
      This%ColnameLogUnit(8)  = "        CoancestryRate"
      This%ColnameLogUnit(9)  = "            Inbreeding"
      This%ColnameLogUnit(10) = "        InbreedingRate"
      nColTmp = nCol
      if (Spec%GenericIndCritGiven) then
        do i = 1, Spec%nGenericIndCrit
          nColTmp = nColTmp + 1
          This%ColnameLogUnit(nColTmp) = "GenIndCrit"//trim(Int2Char(i))
          This%ColnameLogUnit(nColTmp) = adjustr(This%ColnameLogUnit(nColTmp))
        end do
      end if
      if (Spec%GenericMatCritGiven) then
        do i = 1, Spec%nGenericMatCrit
          nColTmp = nColTmp + 1
          This%ColnameLogUnit(nColTmp) = "GenMatCrit"//trim(Int2Char(i))
          This%ColnameLogUnit(nColTmp) = adjustr(This%ColnameLogUnit(nColTmp))
        end do
      end if
      do i = 1, nCol
        Tmp = This%ColnameLogUnit(i)
        This%ColnameLogStdout(i) = Tmp(9:22)
        This%ColnameLogStdout(i) = adjustr(This%ColnameLogStdout(i))
      end do
      This%ColnameLogPopUnit = This%ColnameLogUnit
      !                            1234567890123456789012
      This%ColnameLogPopUnit(2) = "              Solution"
      This%FmtLogStdoutHead  = trim(FMTLOGSTDOUTHEADA) //trim(Int2Char(nCol))  //trim(FMTLOGSTDOUTHEADB)
      This%FmtLogStdout      = trim(FMTLOGSTDOUTA)     //trim(Int2Char(nCol-1))//trim(FMTLOGSTDOUTB)
      This%FmtLogUnitHead    = trim(FMTLOGUNITHEADA)   //trim(Int2Char(nCol)  )//trim(FMTLOGUNITHEADB)
      This%FmtLogUnit        = trim(FMTLOGUNITA)       //trim(Int2Char(nCol-1))//trim(FMTLOGUNITB)
      This%FmtLogPopUnit     = trim(FMTLOGPOPUNITA)    //trim(Int2Char(nCol-2))//trim(FMTLOGUNITB)

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
    subroutine LogAlphaMateSol(This, LogUnit, Iteration, AcceptRate, String, StringNum) ! not pure due to IO
      implicit none
      class(AlphaMateSol), intent(in)              :: This       !< AlphaMateSol holder
      integer(int32), intent(in), optional         :: LogUnit    !< Unit to write to (default STDOUT)
      integer(int32), intent(in)                   :: Iteration  !< Generation/Iteration
      real(real64), intent(in)                     :: AcceptRate !< Acceptance rate
      character(len=*), intent(in), optional       :: String     !< Additional string that will be written before the head
      integer(int32), optional                     :: StringNum  !< How much space is needed for the String
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
          write(Unit, Fmt) Iteration, AcceptRate, This%Objective, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate, This%GenericIndCrit, This%GenericMatCrit
        else
          if (present(String)) then
            write(Unit, StringFmt, Advance="No") trim(adjustl(String))
          end if
          write(Unit, Fmt) Iteration, AcceptRate, This%Objective, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate, This%GenericIndCrit
        end if
      else
        if (allocated(This%GenericMatCrit)) then
          if (present(String)) then
            write(Unit, StringFmt, Advance="No") trim(adjustl(String))
          end if
          write(Unit, Fmt) Iteration, AcceptRate, This%Objective, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate,                      This%GenericMatCrit
        else
          if (present(String)) then
            write(Unit, StringFmt, Advance="No") trim(adjustl(String))
          end if
          write(Unit, Fmt) Iteration, AcceptRate, This%Objective, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate
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
          write(Unit, This%FmtLogPopUnit) Iteration, i, This%Objective, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate, This%GenericIndCrit, This%GenericMatCrit
        else
          write(Unit, This%FmtLogPopUnit) Iteration, i, This%Objective, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate, This%GenericIndCrit
        end if
      else
        if (allocated(This%GenericMatCrit)) then
          write(Unit, This%FmtLogPopUnit) Iteration, i, This%Objective, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate,                      This%GenericMatCrit
        else
          write(Unit, This%FmtLogPopUnit) Iteration, i, This%Objective, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate
        end if
      end if
    end subroutine

    !###########################################################################
end module

!###############################################################################