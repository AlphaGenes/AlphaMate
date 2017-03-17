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
  use IFPort, only : SystemQQ
  use OrderPackModule, only : MrgRnk
  use AlphaHouseMod, only : CountLines, Int2Char, Real2Char, RandomOrder, SetSeed, ToLower, FindLoc
  use AlphaStatMod, only : DescStat, DescStatReal64, DescStatMatrix, DescStatMatrixReal64, DescStatSymMatrix, DescStatLowTriMatrix
  use AlphaEvolveModule, only : AlphaEvolveSol, DifferentialEvolution, RandomSearch
  use AlphaRelateModule

  implicit none

  private
  ! Types
  public :: AlphaMateSpec, AlphaMateData, AlphaMateSol
  ! Functions
  public :: AlphaMateTitle, ReadSpecAndDataForAlphaMate, SetupColNamesAndFormats, AlphaMateSearch

  !> @brief AlphaMate specifications
  type AlphaMateSpec
    ! Biological specifications
    logical :: NrmInsteadOfCoancestry
    real(real64) :: TargetCoancestryRate
    real(real64) :: TargetInbreedingRate
    real(real64) :: CoancestryWeight, InbreedingWeight, SelfingWeight
    logical :: SelCriterionAvailable, GenderMatters, EqualizePar1, EqualizePar2
    real(real64) :: LimitPar1Min, LimitPar1Max, LimitPar2Min, LimitPar2Max, LimitPar1Weight, LimitPar2Weight
    logical :: GenericIndValAvailable, GenericMatValAvailable
    real(real64), allocatable :: GenericIndValWeight(:), GenericMatValWeight(:)
    logical :: SelfingAllowed, CoancestryWeightBellow, InbreedingWeightBellow
    logical :: PAGE, PAGEPar1, PAGEPar2
    integer(int32) :: PAGEPar1Max, PAGEPar2Max
    real(real64) :: PAGEPar1Cost, PAGEPar2Cost
    ! Search specifications
    logical :: ModeMin, ModeRan, ModeOpt
    logical :: EvaluateFrontier
    integer(int32) :: nFrontierSteps
    real(real64), allocatable :: TargetCoancestryRateFrontier(:)
    ! Algorithm specifications
    integer(int32) :: EvolAlgNSol, EvolAlgNIter, EvolAlgNIterBurnIn, EvolAlgNIterStop, EvolAlgNIterPrint, RanAlgStricter
    real(real64) :: EvolAlgStopTol, EvolAlgCRBurnIn, EvolAlgCRLate, EvolAlgFBase, EvolAlgFHigh1, EvolAlgFHigh2
    logical :: EvolAlgLogPop
    contains
      procedure :: Init => InitAlphaMateSpec
      procedure :: Read => ReadAlphaMateSpec
  end type

  !> @brief AlphaMate data
  type AlphaMateData
    ! Raw data
    type(RelMat) :: Coancestry
    type(InbVec) :: Inbreeding
    real(real64), allocatable :: SelCriterion(:), SelCriterionStand(:), SelCriterionPAGE(:), SelCriterionPAGEStand(:)
    integer(int32), allocatable :: Gender(:)
    real(real64), allocatable :: GenericIndVal(:, :), GenericMatVal(:, :, :)
    ! Data summaries
    type(DescStatReal64) :: InbreedingStat, SelCriterionStat, SelCriterionPAGEStat
    type(DescStatReal64), allocatable :: GenericIndValStat(:)
    type(DescStatMatrixReal64) :: CoancestryStat, CoancestryStatGenderDiff, CoancestryStatGender1, CoancestryStatGender2
    type(DescStatMatrixReal64), allocatable :: GenericMatValStat(:)
    ! Derived data
    integer(int32) :: nInd, nMat, nPotMat, nPar, nPotPar1, nPotPar2, nMal, nFem, nPar1, nPar2
    integer(int32), allocatable :: IdPotPar1(:), IdPotPar2(:), IdPotParSeq(:)
    integer(int32) :: nGenericIndVal, nGenericMatVal
    real(real64) :: CurrentCoancestryRanMate, CurrentCoancestryRanMateNoSelf, CurrentCoancestryGenderMate
    real(real64) :: CurrentInbreeding
    real(real64) :: TargetCoancestryRanMate, TargetCoancestryRanMateNoSelf, TargetCoancestryGenderMate
    real(real64) :: TargetInbreeding
    contains
      procedure :: Read => ReadAlphaMateData
  end type

  !> @brief AlphaMate solution
  type, extends(AlphaEvolveSol) :: AlphaMateSol
    real(real64)                :: Penalty
    real(real64)                :: SelCriterion
    real(real64)                :: SelIntensity
    real(real64)                :: FutureCoancestryRanMate
    real(real64)                :: CoancestryRateRanMate
    real(real64)                :: FutureInbreeding
    real(real64)                :: InbreedingRate
    real(real64), allocatable   :: GenericIndVal(:)
    real(real64), allocatable   :: GenericMatVal(:)
    real(real64)                :: Cost
    integer(int32), allocatable :: nVec(:)
    real(real64), allocatable   :: xVec(:)
    integer(int32), allocatable :: MatingPlan(:, :)
    real(real64), allocatable   :: GenomeEdit(:)
    contains
      procedure         :: Initialise => InitialiseAlphaMateSol
      procedure         :: Assign     => AssignAlphaMateSol
      procedure         :: UpdateMean => UpdateMeanAlphaMateSol
      procedure         :: Evaluate   => FixSolEtcMateAndEvaluate
      procedure, nopass :: LogHead    => LogHeadAlphaMateSol
      procedure         :: Log        => LogAlphaMateSol
      procedure, nopass :: LogPopHead => LogPopHeadAlphaMateSol
      procedure         :: LogPop     => LogPopAlphaMateSol
      procedure         :: Write      => WriteAlphaMateSol
  end type

  type(AlphaMateSpec) :: Spec
  type(AlphaMateData) :: Data

  CHARACTER(len=100), PARAMETER :: FMTREAL2CHAR = "(f11.5)"

  CHARACTER(len=100), PARAMETER  :: FMTLOGSTDOUTHEADA = "("
  CHARACTER(len=100), PARAMETER  :: FMTLOGSTDOUTHEADB = "a15)"
  CHARACTER(len=100)             :: FMTLOGSTDOUTHEAD
  CHARACTER(len=100), PARAMETER  :: FMTLOGSTDOUTA = "(i15, "
  CHARACTER(len=100), PARAMETER  :: FMTLOGSTDOUTB = "(4x, f11.5))"
  CHARACTER(len=100)             :: FMTLOGSTDOUT
  CHARACTER(len=15), ALLOCATABLE :: COLNAMELOGSTDOUT(:)

  CHARACTER(len=100), PARAMETER  :: FMTLOGUNITHEADA = "("
  CHARACTER(len=100), PARAMETER  :: FMTLOGUNITHEADB = "a22)"
  CHARACTER(len=100)             :: FMTLOGUNITHEAD
  CHARACTER(len=100), PARAMETER  :: FMTLOGUNITA = "(i22, "
  CHARACTER(len=100), PARAMETER  :: FMTLOGUNITB = "(1x, es21.13e3))"
  CHARACTER(len=100)             :: FMTLOGUNIT
  CHARACTER(len=22), ALLOCATABLE :: COLNAMELOGUNIT(:)

  CHARACTER(len=100), PARAMETER  :: FMTLOGPOPUNITA = "(2i22, "
  CHARACTER(len=100)             :: FMTLOGPOPUNIT
  CHARACTER(len=22), ALLOCATABLE :: COLNAMELOGPOPUNIT(:)

  CHARACTER(len=100), PARAMETER  :: FMTCONTRIBUTIONHEAD = "(6a12)"
  CHARACTER(len=100), PARAMETER  :: FMTCONTRIBUTIONHEADEDIT = "(8a12)"
  CHARACTER(len=100), PARAMETER  :: FMTCONTRIBUTION = "(a12, 1x, i11, 3(1x, f11.5), 1x, i11)"
  CHARACTER(len=100), PARAMETER  :: FMTCONTRIBUTIONEDIT = "(a12, 1x, i11, 3(1x, f11.5), 2(1x, i11), 1x, f11.5)"

  CHARACTER(len=100), PARAMETER  :: FMTMATINGHEAD = "(3a12)"
  CHARACTER(len=100), PARAMETER  :: FMTMATING = "(i12, 2(1x, a11))"

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

    !-------------------------------------------------------------------------
    !> @brief  Print AlphaMate title
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !-------------------------------------------------------------------------
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

    !-------------------------------------------------------------------------
    !> @brief  Initialize AlphaMate specifications
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !-------------------------------------------------------------------------
    pure subroutine InitAlphaMateSpec(This)
      implicit none
      class(AlphaMateSpec), intent(out) :: This !< AlphaMateSpec holder

      ! Inputs
      This%NrmInsteadOfCoancestry = .false.
      This%SelCriterionAvailable = .false.
      This%GenericIndValAvailable = .false.
      This%GenericMatValAvailable = .false.

      ! Biological specifications
      This%TargetCoancestryRate = 0.01d0
      This%TargetInbreedingRate = 0.01d0
      ! This%TargetCoancestryRateFrontier(:) ! allocatable so skip here
      This%CoancestryWeight = 0.5d0
      This%CoancestryWeightBellow = .false.
      This%InbreedingWeight =  0.5d0
      This%InbreedingWeightBellow = .false.
      This%SelfingWeight = 0.0d0
      ! This%GenericIndValWeight(:) ! allocatable so skip here
      ! This%GenericMatValWeight(:) ! allocatable so skip here

      This%GenderMatters = .false.
      This%EqualizePar1 = .false.
      This%EqualizePar2 = .false.
      This%SelfingAllowed = .false.

      This%LimitPar1Min = 1.0d0
      This%LimitPar1Max = huge(Spec%LimitPar1Max) - 1.0d0
      This%LimitPar2Min = 1.0d0
      This%LimitPar2Max = huge(Spec%LimitPar1Max) - 1.0d0
      This%LimitPar1Weight = 0.0d0
      This%LimitPar2Weight = 0.0d0

      This%PAGE = .false.
      This%PAGEPar1 = .false.
      This%PAGEPar2 = .false.
      This%PAGEPar1Max = 0
      This%PAGEPar2Max = 0
      This%PAGEPar1Cost = 0.0d0
      This%PAGEPar2Cost = 0.0d0

      ! Search mode specifications
      This%ModeMin = .false.
      This%ModeRan = .false.
      This%ModeOpt = .false.
      This%EvaluateFrontier = .false.
      This%nFrontierSteps = 0

      ! Search algorithm specifications
      This%EvolAlgNSol = 100
      This%EvolAlgNIter = 10000
      This%EvolAlgNIterBurnIn = 500
      This%EvolAlgNIterStop = 1000
      This%EvolAlgNIterPrint = 100
      This%RanAlgStricter = 10
      This%EvolAlgStopTol = 0.001d0
      This%EvolAlgCRBurnIn = 0.4d0
      This%EvolAlgCRLate = 0.2d0
      This%EvolAlgFBase = 0.1d0
      This%EvolAlgFHigh1 = 1.0d0
      This%EvolAlgFHigh2 = 4.0d0
      This%EvolAlgLogPop = .false.
    end subroutine

    !###########################################################################

    !-------------------------------------------------------------------------
    !> @brief  Read AlphaMate specifications from a file
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !-------------------------------------------------------------------------
    subroutine ReadAlphaMateSpec(This, SpecFile) ! not pure due to IO
      implicit none
      class(AlphaMateSpec), intent(out) :: This     !< AlphaMateSpec holder
      character(len=*), intent(in)      :: SpecFile !< Spec file; when missing, a stub with defaults is created

      call This%Init
      ! TODO
    end subroutine

    !###########################################################################

    !-------------------------------------------------------------------------
    !> @brief  Read AlphaMate data from a file
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !-------------------------------------------------------------------------
    subroutine ReadAlphaMateData(This) ! not pure due to IO
      implicit none
      class(AlphaMateData), intent(inout) :: This !< AlphaMateData holder
      ! TODO
    end subroutine

    !###########################################################################

    subroutine ReadSpecAndDataForAlphaMate ! not pure due to IO
      implicit none

      integer(int32) :: i, j, k, l, m, DumI, jMal, jFem, nIndTmp, GenderTmp, Seed
      integer(int32) :: SpecUnit, SelCriterionUnit, GenderUnit, AvgCoancestryUnit, AvgInbreedingUnit
      integer(int32) :: GenericIndValUnit, GenericMatValUnit
      ! integer(int32), allocatable :: Order(:)

      real(real64) :: SelCriterionTmp, SelCriterionTmp2
      real(real64), allocatable :: GenericIndValTmp(:), GenericMatValTmp(:)

      character(len=1000) :: CoaMtxFile, SelCriterionFile, GenderFile, SeedFile
      character(len=1000) :: GenericIndValFile, GenericMatValFile
      character(len=100) :: DumC, IdCTmp, IdCTmp2

      write(STDOUT, "(a)") "--- Specifications ---"
      write(STDOUT, "(a)") " "

      open(newunit=SpecUnit, file="AlphaMateSpec.txt", status="old")
      write(STDOUT, "(a)") "SpecFile: AlphaMateSpec.txt"

      ! --- Mode ---

      Spec%ModeMin = .false.
      Spec%ModeRan = .false.
      Spec%ModeOpt = .false.

      read(SpecUnit, *) DumC, DumC
      if (index(ToLower(trim(DumC)), "min") > 0) then
        Spec%ModeMin = .true.
        write(STDOUT, "(a)") "Mode: Min"
      end if
      if (index(ToLower(trim(DumC)), "ran") > 0) then
        Spec%ModeRan = .true.
        write(STDOUT, "(a)") "Mode: Ran"
      end if
      if (index(ToLower(trim(DumC)), "opt") > 0) then
        Spec%ModeOpt = .true.
        write(STDOUT, "(a)") "Mode: Opt"
      end if
      if (.not.Spec%ModeMin .and. .not.Spec%ModeRan .and. .not.Spec%ModeOpt) then
        write(STDERR, "(a)") "ERROR: Mode must be: Min, Ran, Opt, or a combination of the three, e.g., MinOpt, RanOpt, ...!"
        write(STDERR, "(a)") " "
        stop 1
      end if

      ! --- CoancestryMatrixFile ---

      ! read(SpecUnit, *) DumC, CoaMtxFile
      ! write(STDOUT, "(a)") "CoancestryMatrixFile: "//trim(CoaMtxFile)
      ! Spec%NrmInsteadOfCoancestry = .false.

      ! --- NrmMatrixFile ---

      read(SpecUnit, *) DumC, CoaMtxFile
      write(STDOUT, "(a)") "NrmMatrixFile: "//trim(CoaMtxFile)
      Spec%NrmInsteadOfCoancestry = .true.

      ! --- SelCriterionFile ---

      read(SpecUnit, *) DumC, SelCriterionFile
      if (ToLower(trim(SelCriterionFile)) /= "none") then
        Spec%SelCriterionAvailable = .true.
        write(STDOUT, "(a)") "SelCriterionFile: "//trim(SelCriterionFile)
      else
        Spec%SelCriterionAvailable = .false.
      end if

      ! --- GenderFile ---

      read(SpecUnit, *) DumC, GenderFile
      if (ToLower(trim(GenderFile)) /= "none") then
        Spec%GenderMatters = .true.
        write(STDOUT, "(a)") "GenderFile: "//trim(GenderFile)
      else
        Spec%GenderMatters = .false.
      end if

      ! --- NumberOfIndividuals ---

      read(SpecUnit, *) DumC, Data%nInd
      write(STDOUT, "(a)") "NumberOfIndividuals: "//trim(Int2Char(Data%nInd))

      ! --- NumberOfMatings ---

      read(SpecUnit, *) DumC, Data%nMat
      write(STDOUT, "(a)") "NumberOfMatings: "//trim(Int2Char(Data%nMat))

      ! --- NumberOfParents ---

      read(SpecUnit, *) DumC, Data%nPar
      write(STDOUT, "(a)") "NumberOfParents: "//trim(Int2Char(Data%nPar))
      if (Data%nPar > Data%nInd) then
        write(STDERR, "(a)") "ERROR: The number of parents can not be larger than the number of individuals!"
        write(STDERR, "(a)") "ERROR: Number of     parents: "//trim(Int2Char(Data%nPar))
        write(STDERR, "(a)") "ERROR: Number of individuals: "//trim(Int2Char(Data%nInd))
        write(STDERR, "(a)") " "
        stop 1
      end if
      if (Data%nMat > Data%nPar) then
        write(STDOUT, "(a)") "NOTE: The number of matings is larger than the number of parents! Was this really the intention?"
        write(STDOUT, "(a)") "NOTE: Number of matings: "//trim(Int2Char(Data%nMat))
        write(STDOUT, "(a)") "NOTE: Number of parents: "//trim(Int2Char(Data%nPar))
        write(STDOUT, "(a)") " "
      end if

      ! --- NumberOfMaleParents ---

      read(SpecUnit, *) DumC, Data%nPar1

      ! --- NumberOfFemaleParents ---

      read(SpecUnit, *) DumC, Data%nPar2

      if (.not.Spec%GenderMatters) then
        Data%nPar1 = Data%nPar
      else
        write(STDOUT, "(a)") "NumberOfMaleParents:   "//trim(Int2Char(Data%nPar1))
        write(STDOUT, "(a)") "NumberOfFemaleParents: "//trim(Int2Char(Data%nPar2))
      end if
      if (Spec%GenderMatters) then
        if ((Data%nPar1 + Data%nPar2) > Data%nInd) then
          write(STDERR, "(a)") "ERROR: The number of parents can not be larger than the number of individuals!"
          write(STDERR, "(a)") "ERROR: Number of        parents: "//trim(Int2Char(Data%nPar1 + Data%nPar2))
          write(STDERR, "(a)") "ERROR: Number of   male parents: "//trim(Int2Char(Data%nPar1))
          write(STDERR, "(a)") "ERROR: Number of female parents: "//trim(Int2Char(Data%nPar2))
          write(STDERR, "(a)") "ERROR: Number of    individuals: "//trim(Int2Char(Data%nInd))
          write(STDERR, "(a)") " "
          stop 1
        end if
        if ((Data%nPar1 + Data%nPar2) /= Data%nPar) then
          write(STDOUT, "(a)") "NOTE: The number of male and female parents does not match with the total number of parents - redefined!"
          write(STDOUT, "(a)") "NOTE: Number of   male parents: "//trim(Int2Char(Data%nPar1))
          write(STDOUT, "(a)") "NOTE: Number of female parents: "//trim(Int2Char(Data%nPar2))
          write(STDOUT, "(a)") "NOTE: Number of        parents: "//trim(Int2Char(Data%nPar))//" (defined)"
          Data%nPar = Data%nPar1 + Data%nPar2
          write(STDOUT, "(a)") "NOTE: Number of        parents: "//trim(Int2Char(Data%nPar))//" (redefined)"
          write(STDOUT, "(a)") " "
        end if
        if ((Data%nMat > Data%nPar1) .and. (Data%nMat > Data%nPar2)) then
          write(STDOUT, "(a)") "NOTE: The number of matings is larger than the number of male and female parents! Was this really the intention?"
          write(STDOUT, "(a)") "NOTE: Number of        matings: "//trim(Int2Char(Data%nMat))
          write(STDOUT, "(a)") "NOTE: Number of   male parents: "//trim(Int2Char(Data%nPar1))
          write(STDOUT, "(a)") "NOTE: Number of female parents: "//trim(Int2Char(Data%nPar2))
          write(STDOUT, "(a)") " "
        end if
      end if

      ! --- EqualizeParentContributions ---

      read(SpecUnit, *) DumC, DumC
      if (.not.Spec%GenderMatters) then
        if      (ToLower(trim(DumC)) == "yes") then
          Spec%EqualizePar1 = .true.
          write(STDOUT, "(a)") "EqualizeParentContributions: yes"
        else if (ToLower(trim(DumC)) == "no") then
          Spec%EqualizePar1 = .false.
          !write(STDOUT, "(a)") "EqualizeParentContributions: no"
        else
          write(STDERR, "(a)") "ERROR: EqualizeParentContributions must be: Yes or no!"
          write(STDERR, "(a)") " "
          stop 1
        end if
      end if

      ! --- EqualizeMaleParentContributions ---

      read(SpecUnit, *) DumC, DumC
      if (Spec%GenderMatters) then
        if      (ToLower(trim(DumC)) == "yes") then
          Spec%EqualizePar1 = .true.
          write(STDOUT, "(a)") "EqualizeMaleParentContributions: yes"
          if (mod(Data%nMat, Data%nPar1) /= 0) then
            ! @todo might consider handling this better at some point
            write(STDERR, "(a)") "ERROR: When contributions are equalized the number of matings needs to divide into the number of parents"
            write(STDERR, "(a)") "ERROR: Number of       matings: "//trim(Int2Char(Data%nMat))
            write(STDERR, "(a)") "ERROR: Number of  male parents: "//trim(Int2Char(Data%nPar1))
            write(STDERR, "(a)") "ERROR: Modulo (should be zero): "//trim(Int2Char(mod(Data%nMat, Data%nPar1)))
            write(STDERR, "(a)") " "
            stop 1
          end if
        else if (ToLower(trim(DumC)) == "no") then
          Spec%EqualizePar1 = .false.
          !write(STDOUT, "(a)") "EqualizeMaleParentContributions: no"
        else
          write(STDERR, "(a)") "ERROR: EqualizeMaleParentContributions must be: Yes or No!"
          write(STDERR, "(a)") " "
          stop 1
        end if
      end if

      ! --- EqualizeFemaleParentContributions ---

      read(SpecUnit, *) DumC, DumC
      if (Spec%GenderMatters) then
        if      (ToLower(trim(DumC)) == "yes") then
          Spec%EqualizePar2 = .true.
          write(STDOUT, "(a)") "EqualizeFemaleParentContributions: yes"
          if (mod(Data%nMat, Data%nPar2) /= 0) then
            ! @todo might consider handling this better at some point
            write(STDERR, "(a)") "ERROR: When contributions are equalized the number of matings needs to divide into the number of parents"
            write(STDERR, "(a)") "ERROR: Number of        matings: "//trim(Int2Char(Data%nMat))
            write(STDERR, "(a)") "ERROR: Number of female parents: "//trim(Int2Char(Data%nPar2))
            write(STDERR, "(a)") "ERROR: Modulo  (should be zero): "//trim(Int2Char(mod(Data%nMat, Data%nPar2)))
            write(STDERR, "(a)") " "
            stop 1
          end if
        else if (ToLower(trim(DumC)) == "no") then
          Spec%EqualizePar2 = .false.
          !write(STDOUT, "(a)") "EqualizeFemaleParentContributions: no"
        else
          write(STDERR, "(a)") "ERROR: EqualizeFemaleParentContributions must be: Yes or No!"
          write(STDERR, "(a)") " "
          stop 1
        end if
      end if

      ! --- LimitParentContributions ---

      read(SpecUnit, *) DumC, DumC
      Spec%LimitPar1Min = 1.0d0
      Spec%LimitPar1Max = huge(Spec%LimitPar1Max) - 1.0d0
      Spec%LimitPar1Weight = 0.0d0
      if (.not.Spec%GenderMatters) then
        if      (ToLower(trim(DumC)) == "yes") then
          if (Spec%EqualizePar1) then
            write(STDOUT, "(a)") "LimitParentContributions: no"
            write(STDOUT, "(a)") "NOTE: Limit parent contributions option ignored when contributions are to be equalized!"
            write(STDOUT, "(a)") " "
          else
            backspace(SpecUnit)
            read(SpecUnit, *) DumC, DumC, Spec%LimitPar1Min, Spec%LimitPar1Max, Spec%LimitPar1Weight
            write(STDOUT, "(a)") "LimitParentContributions: yes, min "//trim(Int2Char(nint(Spec%LimitPar1Min)))//", max "//&
              trim(Int2Char(nint(Spec%LimitPar1Max)))//", penalty weight "//trim(Real2Char(Spec%LimitPar1Weight, fmt=FMTREAL2CHAR))
            if (Spec%LimitPar1Weight > 0.0) then
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
      if (Spec%GenderMatters) then
        if      (ToLower(trim(DumC)) == "yes") then
          if (Spec%EqualizePar1) then
            write(STDOUT, "(a)") "LimitMaleParentContributions: no"
            write(STDOUT, "(a)") "NOTE: Limit male parent contributions option ignored when contributions are to be equalized!"
            write(STDOUT, "(a)") " "
          else
            backspace(SpecUnit)
            read(SpecUnit, *) DumC, DumC, Spec%LimitPar1Min, Spec%LimitPar1Max, Spec%LimitPar1Weight
            write(STDOUT, "(a)") "LimitMaleParentContributions: yes, min "//trim(Int2Char(nint(Spec%LimitPar1Min)))//", max "//&
              trim(Int2Char(nint(Spec%LimitPar1Max)))//", penalty weight "//trim(Real2Char(Spec%LimitPar1Weight, fmt=FMTREAL2CHAR))
            if (Spec%LimitPar1Weight > 0.0) then
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
      Spec%LimitPar2Min = 1.0d0
      Spec%LimitPar2Max = huge(Spec%LimitPar2Max) - 1.0d0
      Spec%LimitPar2Weight = 0.0d0
      if (Spec%GenderMatters) then
        if      (ToLower(trim(DumC)) == "yes") then
          if (Spec%EqualizePar2) then
            write(STDOUT, "(a)") "LimitFemaleParentContributions: no"
            write(STDOUT, "(a)") "NOTE: Limit female parent contributions option ignored when contributions are to be equalized!"
            write(STDOUT, "(a)") " "
          else
            backspace(SpecUnit)
            read(SpecUnit, *) DumC, DumC, Spec%LimitPar2Min, Spec%LimitPar2Max, Spec%LimitPar2Weight
            write(STDOUT, "(a)") "LimitFemaleParentContributions: yes, min "//trim(Int2Char(nint(Spec%LimitPar2Min)))//", max "//&
              trim(Int2Char(nint(Spec%LimitPar2Max)))//", penalty weight "//trim(Real2Char(Spec%LimitPar2Weight, fmt=FMTREAL2CHAR))
            if (Spec%LimitPar2Weight > 0.0) then
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
        Spec%SelfingAllowed = .true.
        write(STDOUT, "(a)") "AllowSelfing: Yes"
        if (Spec%GenderMatters) then
          write(STDOUT, "(a)") "ERROR: When gender matters, AlphaMate can not perform selfing! See the manual for a solution."
          write(STDOUT, "(a)") " "
          stop 1
        end if
      else if (ToLower(trim(DumC)) == "no") then
        Spec%SelfingAllowed = .false.
        if (.not.Spec%GenderMatters) then
          backspace(SpecUnit)
          read(SpecUnit, *) DumC, DumC, Spec%SelfingWeight
          write(STDOUT, "(a)") "AllowSelfing: no, penalty weight "//trim(Real2Char(Spec%SelfingWeight, fmt=FMTREAL2CHAR))
          if (Spec%SelfingWeight > 0.0) then
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
      if (.not.Spec%GenderMatters) then
        if      (ToLower(trim(DumC)) == "yes") then
          Spec%PAGEPar1 = .true.
          backspace(SpecUnit)
          read(SpecUnit, *) DumC, DumC, Spec%PAGEPar1Max!, Spec%PAGEPar1Cost
          write(STDOUT, "(a)") "PAGE: yes, no. of individuals "//trim(Int2Char(Spec%PAGEPar1Max))//", cost "//trim(Real2Char(Spec%PAGEPar1Cost, fmt=FMTREAL2CHAR))
          if (.not.Spec%SelCriterionAvailable) then
            write(STDERR, "(a)") "ERROR: Can not use PAGE when selection criterion file is not given!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          if (Spec%PAGEPar1Max > Data%nPar) then
            write(STDERR, "(a)") "ERROR: The max number of individuals to edit must not be greater than the total number of parents!"
            write(STDERR, "(a)") "ERROR: Number of             parents: "//trim(Int2Char(Data%nPar))
            write(STDERR, "(a)") "ERROR: Number of individuals to edit: "//trim(Int2Char(Spec%PAGEPar1Max))
            write(STDERR, "(a)") " "
          end if
        else if (ToLower(trim(DumC)) == "no") then
          Spec%PAGEPar1 = .false.
          !write(STDOUT, "(a)") "PAGE: no"
        else
          write(STDERR, "(a)") "ERROR: PAGE must be: Yes or No!"
          write(STDERR, "(a)") " "
          stop 1
        end if
      end if

      ! --- PAGEMales ---

      read(SpecUnit, *) DumC, DumC
      if (Spec%GenderMatters) then
        if      (ToLower(trim(DumC)) == "yes") then
          Spec%PAGEPar1 = .true.
          backspace(SpecUnit)
          read(SpecUnit, *) DumC, DumC, Spec%PAGEPar1Max!, Spec%PAGEPar1Cost
          write(STDOUT, "(a)") "PAGEMales: yes, no. of individuals "//trim(Int2Char(Spec%PAGEPar1Max))//", cost "//trim(Real2Char(Spec%PAGEPar1Cost, fmt=FMTREAL2CHAR))
          if (.not.Spec%SelCriterionAvailable) then
            write(STDERR, "(a)") "ERROR: Can not use PAGE when selection criterion file is not given!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          if (Spec%PAGEPar1Max > Data%nPar1) then
            write(STDERR, "(a)") "ERROR: The max number of male individuals to edit must not be greater than the total number of male parents!"
            write(STDERR, "(a)") "ERROR: Number of male             parents: "//trim(Int2Char(Data%nPar1))
            write(STDERR, "(a)") "ERROR: Number of male individuals to edit: "//trim(Int2Char(Spec%PAGEPar1Max))
            write(STDERR, "(a)") " "
          end if
        else if (ToLower(trim(DumC)) == "no") then
          Spec%PAGEPar1 = .false.
          !write(STDOUT, "(a)") "PAGEMales: no"
        else
          write(STDERR, "(a)") "ERROR: PAGEMales must be: Yes or No!"
          write(STDERR, "(a)") " "
          stop 1
        end if
      end if

      ! --- PAGEFemales ---

      read(SpecUnit, *) DumC, DumC
      if (Spec%GenderMatters) then
        if      (ToLower(trim(DumC)) == "yes") then
          Spec%PAGEPar2 = .false.
          backspace(SpecUnit)
          read(SpecUnit, *) DumC, DumC, Spec%PAGEPar2Max!, Spec%PAGEPar2Cost
          write(STDOUT, "(a)") "PAGEFemales: yes, no. of individuals "//trim(Int2Char(Spec%PAGEPar2Max))//", cost "//trim(Real2Char(Spec%PAGEPar2Cost, fmt=FMTREAL2CHAR))
          if (.not.Spec%SelCriterionAvailable) then
            write(STDERR, "(a)") "ERROR: Can not use PAGE when selection criterion file is not given!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          if (Spec%PAGEPar2Max > Data%nPar2) then
            write(STDERR, "(a)") "ERROR: The max number of female individuals to edit must not be greater than the total number of female parents!"
            write(STDERR, "(a)") "ERROR: Number of female             parents: "//trim(Int2Char(Data%nPar2))
            write(STDERR, "(a)") "ERROR: Number of female individuals to edit: "//trim(Int2Char(Spec%PAGEPar2Max))
            write(STDERR, "(a)") " "
          end if
        else if (ToLower(trim(DumC)) == "no") then
          Spec%PAGEPar2 = .false.
          !write(STDOUT, "(a)") "PAGEFemales: no"
        else
          write(STDERR, "(a)") "ERROR: PAGEFemales must be: Yes or No!"
          write(STDERR, "(a)") " "
          stop 1
        end if
      end if

      if (Spec%PAGEPar1 .or. Spec%PAGEPar2) then
        Spec%PAGE = .true.
      else
        Spec%PAGE = .false.
      end if

      ! --- TargetedRateOfCoancestry ---

      read(SpecUnit, *) DumC, Spec%TargetCoancestryRate, Spec%CoancestryWeight, DumC
      write(STDOUT, "(a)") "TargetedRateOfCoancestry: "//trim(Real2Char(Spec%TargetCoancestryRate, fmt=FMTREAL2CHAR))//&
        ", penalty weight "//trim(Real2Char(Spec%CoancestryWeight, fmt=FMTREAL2CHAR))//", mode "//trim(DumC)
      if      (ToLower(trim(DumC)) == "above") then
        Spec%CoancestryWeightBellow = .false.
      else if (ToLower(trim(DumC)) == "aboveandbellow") then
        Spec%CoancestryWeightBellow = .true.
      else
        write(STDERR, "(a)") "ERROR: CoancestryWeightMode must be: Above or AboveAndBellow!"
        write(STDERR, "(a)") " "
        stop 1
      end if
      if (Spec%TargetCoancestryRate == 0.0) then
        write(STDERR, "(a)") "ERROR: Can not work with the targeted rate of coancestry exactly equal to zero - it is numerically unstable!"
        write(STDERR, "(a)") " "
        stop 1
      end if
      if (Spec%CoancestryWeight > 0.0) then
        write(STDERR, "(a)") "ERROR: Penalty weight for the targeted rate of coancestry should be zero or negative!"
        write(STDERR, "(a)") " "
        stop 1
      end if

      ! --- TargetedRateOfInbreeding (=inbreeding of a mating) ---

      read(SpecUnit, *) DumC, Spec%TargetInbreedingRate, Spec%InbreedingWeight, DumC
      write(STDOUT, "(a)") "TargetedRateOfInbreeding: "//trim(Real2Char(Spec%TargetInbreedingRate, fmt=FMTREAL2CHAR))//&
        ", penalty weight "//trim(Real2Char(Spec%InbreedingWeight, fmt=FMTREAL2CHAR))//", mode "//trim(DumC)
      if      (ToLower(trim(DumC)) == "above") then
        Spec%InbreedingWeightBellow = .false.
      else if (ToLower(trim(DumC)) == "aboveandbellow") then
        Spec%InbreedingWeightBellow = .true.
      else
        write(STDERR, "(a)") "ERROR: InbreedingWeightMode must be: Above or AboveAndBellow!"
        write(STDERR, "(a)") " "
        stop 1
      end if
      if (Spec%TargetInbreedingRate == 0.0) then
        write(STDERR, "(a)") "ERROR: Can not work with the targeted rate of inbreeding exactly equal to zero - it is numerically unstable!"
        write(STDERR, "(a)") " "
        stop 1
      end if
      if (Spec%InbreedingWeight > 0.0) then
        write(STDERR, "(a)") "ERROR: Penalty weight for the targeted rate of inbreeding should be zero or negative!"
        write(STDERR, "(a)") " "
        stop 1
      end if

      ! --- EvaluateFrontier ---

      ! @todo Should this be part of TargetedCoancestryRate??
      read(SpecUnit, *) DumC, DumC
      if      (ToLower(trim(DumC)) == "no") then
        Spec%EvaluateFrontier = .false.
        !write(STDOUT, "(a)") "EvaluateFrontier: no"
      else if (ToLower(trim(DumC)) == "yes") then
        Spec%EvaluateFrontier = .true.
        backspace(SpecUnit)
        read(SpecUnit, *) DumC, DumC, Spec%nFrontierSteps
        allocate(Spec%TargetCoancestryRateFrontier(Spec%nFrontierSteps))
        backspace(SpecUnit)
        read(SpecUnit, *) DumC, DumC, Spec%nFrontierSteps, Spec%TargetCoancestryRateFrontier(:)
        write(STDOUT, "("//Int2Char(1+Spec%nFrontierSteps)//"a)") "EvaluateFrontier: yes, #steps: "//trim(Int2Char(Spec%nFrontierSteps))//&
          ", rates of coancestry: ", (trim(Real2Char(Spec%TargetCoancestryRateFrontier(i), fmt=FMTREAL2CHAR)), i = 1, Spec%nFrontierSteps)
        if (any(Spec%TargetCoancestryRateFrontier(:) == 0.0)) then
          write(STDERR, "(a)") "ERROR: Can not work with TargetCoancestryRateFrontier equal to zero - it is numerically unstable!"
          write(STDERR, "(a)") " "
          stop 1
        end if
      else
        write(STDERR, "(a)") "ERROR: EvaluateFrontier must be: Yes or No!"
        write(STDERR, "(a)") " "
        stop 1
      end if

      ! --- EvolutionaryAlgorithmIterations ---

      read(SpecUnit, *) DumC, Spec%EvolAlgNSol, Spec%EvolAlgNIter, Spec%EvolAlgNIterBurnIn, Spec%EvolAlgNIterStop, Spec%EvolAlgStopTol, Spec%EvolAlgNIterPrint, Spec%EvolAlgLogPop

      ! --- EvolutionaryAlgorithmParameters ---

      read(SpecUnit, *) DumC, Spec%EvolAlgCRBurnIn, Spec%EvolAlgCRLate, Spec%EvolAlgFBase, Spec%EvolAlgFHigh1, Spec%EvolAlgFHigh2

      ! --- RandomSearchIterationsStricter ---

      read(SpecUnit, *) DumC, Spec%RanAlgStricter

      ! --- Seed ---

      read(SpecUnit, *) DumC, DumC
      SeedFile = "Seed.txt"
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

      Data%nGenericIndVal = 0
      read(SpecUnit, *) DumC, GenericIndValFile
      if (ToLower(trim(GenericIndValFile)) == "none") then
        Spec%GenericIndValAvailable = .false.
        read(SpecUnit, *) DumC
      else
        Spec%GenericIndValAvailable = .true.
        backspace(SpecUnit)
        read(SpecUnit, *) DumC, GenericIndValFile, Data%nGenericIndVal
        allocate(Spec%GenericIndValWeight(Data%nGenericIndVal))
        read(SpecUnit, *) DumC, Spec%GenericIndValWeight(:)
        write(STDOUT, "("//Int2Char(1 + Data%nGenericIndVal)//"a)") "GenericIndividualValuesFile: "//trim(GenericIndValFile)//&
          ", weight(s): ", (trim(Real2Char(Spec%GenericIndValWeight(i), fmt=FMTREAL2CHAR)), i = 1, Data%nGenericIndVal)
      end if

      ! --- GenericMatingValuesFile ---
      ! --- GenericMatingValuesWeight ---

      Data%nGenericMatVal = 0
      read(SpecUnit, *) DumC, GenericMatValFile
      if (ToLower(trim(GenericMatValFile)) == "none") then
        Spec%GenericMatValAvailable = .false.
        read(SpecUnit, *) DumC
      else
        Spec%GenericMatValAvailable = .true.
        backspace(SpecUnit)
        read(SpecUnit, *) DumC, GenericMatValFile, Data%nGenericMatVal
        allocate(Spec%GenericMatValWeight(Data%nGenericMatVal))
        read(SpecUnit, *) DumC, Spec%GenericMatValWeight(:)
        write(STDOUT, "("//Int2Char(1 + Data%nGenericMatVal)//"a)") "GenericMatingValuesFile: "//trim(GenericMatValFile)//&
          ", weight(s): ", (trim(Real2Char(Spec%GenericMatValWeight(i), fmt=FMTREAL2CHAR)), i = 1, Data%nGenericMatVal)
      end if

      write(STDOUT, "(a)") " "
      close(SpecUnit)

      write(STDOUT, "(a)") "--- Data ---"
      write(STDOUT, "(a)") " "

      ! --- Coancestry ---

      write(STDOUT, "(a)") "Coancestry (average identity of the four genome combinations of two individuals)"

      call Data%Coancestry%Read(File=CoaMtxFile)
      if (Data%Coancestry%nInd < Data%nInd) then
        write(STDERR, "(a)") "ERROR: The coancestry matrix file has less rows than there are defined number of individuals!"
        write(STDERR, "(a)") "ERROR: Number of defined individuals:                       "//trim(Int2Char(Data%nInd))
        write(STDERR, "(a)") "ERROR: Number of individuals in the coancestry matrix file: "//trim(Int2Char(Data%Coancestry%nInd))
        write(STDERR, "(a)") " "
      end if
      if (Spec%NrmInsteadOfCoancestry) then
        call Data%Coancestry%Nrm2Coancestry
      end if

      Data%CoancestryStat = DescStatSymMatrix(Data%Coancestry%Value(1:, 1:))
      if (Spec%GenderMatters) then
        ! TODO
        Data%CoancestryStatGenderDiff = DescStatSymMatrix(Data%Coancestry%Value(1:, 1:))
        Data%CoancestryStatGender1    = DescStatSymMatrix(Data%Coancestry%Value(1:, 1:))
        Data%CoancestryStatGender2    = DescStatSymMatrix(Data%Coancestry%Value(1:, 1:))
      end if

      ! Current
      Data%CurrentCoancestryRanMate       = Data%CoancestryStat%All%Mean
      Data%CurrentCoancestryRanMateNoSelf = Data%CoancestryStat%OffDiag%Mean
      if (Spec%GenderMatters) then
        Data%CurrentCoancestryGenderMate  = Data%CoancestryStatGenderDiff%All%Mean
      end if

      ! Obtain limit/target based on given rates
      ! F_t = DeltaF + (1 - DeltaF) * F_t-1
      Data%TargetCoancestryRanMate       = Spec%TargetCoancestryRate + (1.0d0 - Spec%TargetCoancestryRate) * Data%CurrentCoancestryRanMate
      Data%TargetCoancestryRanMateNoSelf = Spec%TargetCoancestryRate + (1.0d0 - Spec%TargetCoancestryRate) * Data%CurrentCoancestryRanMateNoSelf
      if (Spec%GenderMatters) then
        Data%TargetCoancestryGenderMate  = Spec%TargetCoancestryRate + (1.0d0 - Spec%TargetCoancestryRate) * Data%CurrentCoancestryGenderMate
      end if

      write(STDOUT, "(a)") "  - coancestry among individuals (including self-coancestry)"
      write(STDOUT, "(a)") "    (expected inbreeding under random mating, including selfing)"
      write(STDOUT, "(a)") "    - average: "//trim(Real2Char(Data%CoancestryStat%All%Mean, fmt=FMTREAL2CHAR))//", limit/target: "//trim(Real2Char(Data%TargetCoancestryRanMate, fmt=FMTREAL2CHAR))
      write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(Data%CoancestryStat%All%SD,   fmt=FMTREAL2CHAR))
      write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(Data%CoancestryStat%All%Min,  fmt=FMTREAL2CHAR))
      write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(Data%CoancestryStat%All%Max,  fmt=FMTREAL2CHAR))

      write(STDOUT, "(a)") ""
      write(STDOUT, "(a)") "  - coancestry between individuals"
      write(STDOUT, "(a)") "    (expected inbreeding under random mating, excluding selfing)"
      write(STDOUT, "(a)") "    - average: "//trim(Real2Char(Data%CoancestryStat%OffDiag%Mean, fmt=FMTREAL2CHAR))//", limit/target: "//trim(Real2Char(Data%TargetCoancestryRanMateNoSelf, fmt=FMTREAL2CHAR))
      write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(Data%CoancestryStat%OffDiag%SD,   fmt=FMTREAL2CHAR))
      write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(Data%CoancestryStat%OffDiag%Min,  fmt=FMTREAL2CHAR))
      write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(Data%CoancestryStat%OffDiag%Max,  fmt=FMTREAL2CHAR))

      if (Spec%GenderMatters) then
        write(STDOUT, "(a)") ""
        write(STDOUT, "(a)") "  - coancestry between individuals of different gender"
        write(STDOUT, "(a)") "    (expected inbreeding under random mating, excluding selfing and equal-gender mating)"
        write(STDOUT, "(a)") "    - average: "//trim(Real2Char(Data%CoancestryStatGenderDiff%All%Mean, fmt=FMTREAL2CHAR))//", limit/target: "//trim(Real2Char(Data%TargetCoancestryGenderMate, fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(Data%CoancestryStatGenderDiff%All%SD,   fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(Data%CoancestryStatGenderDiff%All%Min,  fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(Data%CoancestryStatGenderDiff%All%Max,  fmt=FMTREAL2CHAR))

        write(STDOUT, "(a)") ""
        write(STDOUT, "(a)") "  - coancestry between individuals of gender 1"
        write(STDOUT, "(a)") "    - average: "//trim(Real2Char(Data%CoancestryStatGender1%All%Mean, fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(Data%CoancestryStatGender1%All%SD,   fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(Data%CoancestryStatGender1%All%Min,  fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(Data%CoancestryStatGender1%All%Max,  fmt=FMTREAL2CHAR))

        write(STDOUT, "(a)") ""
        write(STDOUT, "(a)") "  - coancestry between individuals of gender 2"
        write(STDOUT, "(a)") "    - average: "//trim(Real2Char(Data%CoancestryStatGender2%All%Mean, fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(Data%CoancestryStatGender2%All%SD,   fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(Data%CoancestryStatGender2%All%Min,  fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(Data%CoancestryStatGender2%All%Max,  fmt=FMTREAL2CHAR))
      end if

      ! Report
      open(newunit=AvgCoancestryUnit, file="AverageCoancestry.txt", status="unknown")
      write(AvgCoancestryUnit, "(a, f)") "Current (random mating),                  ", Data%CurrentCoancestryRanMate
      write(AvgCoancestryUnit, "(a, f)") "Current (random mating, no selfing),      ", Data%CurrentCoancestryRanMateNoSelf
      write(AvgCoancestryUnit, "(a, f)") "Limit/Target (random mating),             ", Data%TargetCoancestryRanMate
      write(AvgCoancestryUnit, "(a, f)") "Limit/Target (random mating, no selfing), ", Data%TargetCoancestryRanMateNoSelf
      close(AvgCoancestryUnit)

      ! --- Inbreeding ---

      write(STDOUT, "(a)") " "
      write(STDOUT, "(a)") "Inbreeding (identity between the two genomes of an individual)"

      call Data%Coancestry%Inbreeding(Out=Data%Inbreeding, Nrm=.false.)

      ! Current
      Data%InbreedingStat = DescStat(Data%Inbreeding%Value(1:))
      Data%CurrentInbreeding = Data%InbreedingStat%Mean

      ! Obtain limit/target based on given rates
      ! F_t = DeltaF + (1 - DeltaF) * F_t-1
      Data%TargetInbreeding = Spec%TargetInbreedingRate + (1.0d0 - Spec%TargetInbreedingRate) * Data%CurrentInbreeding

      write(STDOUT, "(a)") "  - average: "//trim(Real2Char(Data%InbreedingStat%Mean, fmt=FMTREAL2CHAR))//", limit/target: "//trim(Real2Char(Data%TargetInbreeding, fmt=FMTREAL2CHAR))
      write(STDOUT, "(a)") "  - st.dev.: "//trim(Real2Char(Data%InbreedingStat%SD,   fmt=FMTREAL2CHAR))
      write(STDOUT, "(a)") "  - minimum: "//trim(Real2Char(Data%InbreedingStat%Min,  fmt=FMTREAL2CHAR))
      write(STDOUT, "(a)") "  - maximum: "//trim(Real2Char(Data%InbreedingStat%Max,  fmt=FMTREAL2CHAR))
      write(STDOUT, "(a)") " "

      open(newunit=AvgInbreedingUnit, file="AverageInbreeding.txt", status="unknown")
      write(AvgInbreedingUnit, "(a, f)") "Current,      ", Data%CurrentInbreeding
      write(AvgInbreedingUnit, "(a, f)") "Limit/Target, ", Data%TargetInbreeding
      close(AvgInbreedingUnit)

      ! --- Selection criterion ---

      allocate(Data%SelCriterion(Data%nInd))
      allocate(Data%SelCriterionStand(Data%nInd))
      if (Spec%PAGE) then
        allocate(Data%SelCriterionPAGE(Data%nInd))
        allocate(Data%SelCriterionPAGEStand(Data%nInd))
      end if

      if (.not.Spec%SelCriterionAvailable) then
        Data%SelCriterion(:) = 0.0d0
        Data%SelCriterionStand(:) = 0.0d0
      else
        write(STDOUT, "(a)") "Selection criterion"
        nIndTmp = CountLines(SelCriterionFile)
        if (nIndTmp /= Data%nInd) then
          write(STDERR, "(a)") "ERROR: Number of individuals in the selection criterion file and the coancestry matrix file is not the same!"
          write(STDERR, "(a)") "ERROR: Number of individuals in the coancestry matrix file:   "//trim(Int2Char(Data%nInd))
          write(STDERR, "(a)") "ERROR: Number of individuals in the selection criterion file: "//trim(Int2Char(nIndTmp))
          write(STDERR, "(a)") " "
          stop 1
        end if
        open(newunit=SelCriterionUnit, file=trim(SelCriterionFile), status="old")
        do i = 1, Data%nInd
          if (Spec%PAGE) then
            read(SelCriterionUnit, *) IdCTmp, SelCriterionTmp, SelCriterionTmp2
          else
            read(SelCriterionUnit, *) IdCTmp, SelCriterionTmp
          end if
          j = FindLoc(IdCTmp, Data%Coancestry%OriginalId(1:))
          if (j == 0) then
            write(STDERR, "(a)") "ERROR: Individual "//trim(IdCTmp)//" from the selection criterion file not present in the coancestry matrix file!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          Data%SelCriterion(j) = SelCriterionTmp
          if (Spec%PAGE) then
            Data%SelCriterionPAGE(j) = SelCriterionTmp2
          end if
        end do
        close(SelCriterionUnit)

        Data%SelCriterionStat = DescStat(Data%SelCriterion)
        Data%SelCriterionStand(:) = (Data%SelCriterion(:) - Data%SelCriterionStat%Mean) / Data%SelCriterionStat%SD
        write(STDOUT, "(a)") "  - average: "//trim(Real2Char(Data%SelCriterionStat%Mean, fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "  - st.dev.: "//trim(Real2Char(Data%SelCriterionStat%SD,   fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "  - minimum: "//trim(Real2Char(Data%SelCriterionStat%Min,  fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") "  - maximum: "//trim(Real2Char(Data%SelCriterionStat%Max,  fmt=FMTREAL2CHAR))
        write(STDOUT, "(a)") " "

        if (Data%SelCriterionStat%SD == 0.0) then
          write(STDERR, "(a)") "ERROR: There is no variation in values!"
          write(STDERR, "(a)") "ERROR: Is this intentional?"
          write(STDERR, "(a)") " "
          stop 1
        end if

        if (Spec%PAGE) then
          ! must have the same scaling as selection criterion!!!!
          Data%SelCriterionPAGEStand(:) = (Data%SelCriterionPAGE(:) - Data%SelCriterionStat%Mean) / Data%SelCriterionStat%SD
          ! only the PAGE bit of SelCriterion
          Data%SelCriterionPAGE(:) = Data%SelCriterionPAGE(:) - Data%SelCriterion(:)
          Data%SelCriterionPAGEStand(:) = Data%SelCriterionPAGEStand(:) - Data%SelCriterionStand(:)
          Data%SelCriterionPAGEStat = DescStat(Data%SelCriterionPAGE)
          write(STDOUT, "(a)") "Genome editing increments"
          write(STDOUT, "(a)") "  - average: "//trim(Real2Char(Data%SelCriterionPAGEStat%Mean, fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "  - st.dev.: "//trim(Real2Char(Data%SelCriterionPAGEStat%SD,   fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "  - minimum: "//trim(Real2Char(Data%SelCriterionPAGEStat%Min,  fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "  - maximum: "//trim(Real2Char(Data%SelCriterionPAGEStat%Max,  fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") " "
        end if

        if (Data%SelCriterionPAGEStat%SD == 0.0) then
          write(STDERR, "(a)") "ERROR: There is no variation in values!"
          write(STDERR, "(a)") "ERROR: Is this intentional?"
          write(STDERR, "(a)") " "
          stop 1
        end if

      end if

      ! --- Gender ---

      allocate(Data%Gender(Data%nInd))

      Data%Gender(:) = 0
      if (Spec%GenderMatters) then
        write(STDOUT, "(a)") "Gender"
        nIndTmp = CountLines(GenderFile)
        if (nIndTmp /= Data%nInd) then
          write(STDERR, "(a)") "ERROR: Number of individuals in the gender file and the coancestry matrix file is not the same!"
          write(STDERR, "(a)") "ERROR: Number of individuals in the coancestry matrix file: "//trim(Int2Char(Data%nInd))
          write(STDERR, "(a)") "ERROR: Number of individuals in the gender file:            "//trim(Int2Char(nIndTmp))
          write(STDERR, "(a)") " "
          stop 1
        end if

        Data%nMal = 0
        Data%nFem = 0

        open(newunit=GenderUnit, file=trim(GenderFile), status="old")
        do i = 1, Data%nInd
          read(GenderUnit, *) IdCTmp, GenderTmp
          if      (GenderTmp == 1) then
            Data%nMal = Data%nMal + 1
          else if (GenderTmp == 2) then
            Data%nFem = Data%nFem + 1
          else
            write(STDERR, "(a)") "ERROR: Gender code must be either 1 for male individuals or 2 for female individuals!"
            write(STDERR, "(a)") "ERROR: "//trim(Int2Char(i))//" "//trim(IdCTmp)//" "//trim(Int2Char(GenderTmp))
            write(STDERR, "(a)") " "
            stop 1
          end if
          j = FindLoc(IdCTmp, Data%Coancestry%OriginalId(1:))
          if (j == 0) then
            write(STDERR, "(a)") "ERROR: Individual "//trim(IdCTmp)//" from the gender file not present in the coancestry matrix file!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          Data%Gender(j) = GenderTmp
        end do
        close(GenderUnit)

        write(STDOUT, "(a)") "  - number of   males in data: "//trim(Int2Char(Data%nMal))
        write(STDOUT, "(a)") "  - number of females in data: "//trim(Int2Char(Data%nFem))
        write(STDOUT, "(a)") " "

        if (Data%nPar1 > Data%nMal) then
          write(STDERR, "(a)") "ERROR: The number of male parents can not be larger than the number of males"
          write(STDERR, "(a)") "ERROR: Number of male parents: "//trim(Int2Char(Data%nPar1))
          write(STDERR, "(a)") "ERROR: Number of        males: "//trim(Int2Char(Data%nMal))
          write(STDERR, "(a)") " "
          stop 1
        end if
        if (Data%nPar2 > Data%nFem) then
          write(STDERR, "(a)") "ERROR: The number of female parents can not be larger than the number of females"
          write(STDERR, "(a)") "ERROR: Number of female parents: "//trim(Int2Char(Data%nPar2))
          write(STDERR, "(a)") "ERROR: Number of        females: "//trim(Int2Char(Data%nFem))
          write(STDERR, "(a)") " "
          stop 1
        end if
      end if

      ! --- Define potential parents ---

      if (.not.Spec%GenderMatters) then
        Data%nPotPar1 = Data%nInd
        Data%nPotPar2 = Data%nInd
        allocate(Data%IdPotPar1(Data%nPotPar1))
        do i = 1, Data%nInd
          Data%IdPotPar1(i) = i
        end do
      else
        Data%nPotPar1 = Data%nMal
        Data%nPotPar2 = Data%nFem
        allocate(Data%IdPotPar1(Data%nPotPar1))
        allocate(Data%IdPotPar2(Data%nPotPar2))
        allocate(Data%IdPotParSeq(Data%nInd))
        jMal = 0
        jFem = 0
        do i = 1, Data%nInd
          if (Data%Gender(i) == 1) then
            jMal = jMal + 1
            Data%IdPotPar1(jMal) = i
            Data%IdPotParSeq(i) = jMal
          else
            jFem = jFem + 1
            Data%IdPotPar2(jFem) = i
            Data%IdPotParSeq(i) = jFem
          end if
        end do
        if (Data%nPar1 > Data%nPotPar1) then
          Data%nPar1 = Data%nPotPar1
          write(STDOUT, "(a)") "NOTE: The number of male parents reduced to the number of male individuals!"
          write(STDOUT, "(a)") " "
        end if
        if (Data%nPar2 > Data%nPotPar2) then
          Data%nPar2 = Data%nPotPar2
          write(STDOUT, "(a)") "NOTE: The number of female parents reduced to the number of female individuals!"
          write(STDOUT, "(a)") " "
        end if
      end if

      ! --- Number of all potential matings ---

      if (Spec%GenderMatters) then
        Data%nPotMat = Data%nPotPar1 * Data%nPotPar2
      else
        Data%nPotMat = real(Data%nPotPar1 * Data%nPotPar1) / 2
        if (Spec%SelfingAllowed) then
          Data%nPotMat = nint(Data%nPotMat + real(Data%nPotPar1) / 2)
        else
          Data%nPotMat = nint(Data%nPotMat - real(Data%nPotPar1) / 2)
        end if
      end if
      if (Data%nMat > Data%nPotMat) then
        write(STDOUT, "(a)") "NOTE: The number of matings is larger than the number of all potential matings!"
        write(STDOUT, "(a)") "NOTE: Number of all potential matings: "//trim(Int2Char(Data%nPotMat))
        if (Spec%GenderMatters) then
          write(STDOUT, "(a)") "NOTE: = all males with all females"
          write(STDOUT, "(a)") "NOTE: = no. of males ("//trim(Int2Char(Data%nPotPar1))//") * no. of females ("//trim(Int2Char(Data%nPotPar2))//")"
        else
          if (Spec%SelfingAllowed) then
            write(STDOUT, "(a)") "NOTE: = half-diallel including selfing"
            write(STDOUT, "(a)") "NOTE: = no. of individuals * no. of individuals / 2 + individuals / 2"
          else
            write(STDOUT, "(a)") "NOTE: = half-diallel excluding selfing"
            write(STDOUT, "(a)") "NOTE: = no. of individuals * no. of individuals / 2 - individuals / 2"
          end if
          write(STDOUT, "(a)") "NOTE:   (no. of individuals = "//trim(Int2Char(Data%nPotPar1))//")"
        end if
        write(STDOUT, "(a)") "NOTE: Number of              matings: "//trim(Int2Char(Data%nMat))
        write(STDOUT, "(a)") " "
      end if

      ! --- Generic individual values ---

      if (Spec%GenericIndValAvailable) then
        write(STDOUT, "(a)") "Generic individual values"
        nIndTmp = CountLines(GenericIndValFile)
        if (nIndTmp /= Data%nInd) then
          write(STDERR, "(a)") "ERROR: Number of individuals in the generic individual values file and the coancestry matrix file is not the same!"
          write(STDERR, "(a)") "ERROR: Number of individuals in the coancestry matrix file:         "//trim(Int2Char(Data%nInd))
          write(STDERR, "(a)") "ERROR: Number of individuals in the generic individual values file: "//trim(Int2Char(nIndTmp))
          write(STDERR, "(a)") " "
          stop 1
        end if
        allocate(Data%GenericIndVal(Data%nInd, Data%nGenericIndVal))
        allocate(GenericIndValTmp(Data%nGenericIndVal))
        Data%GenericIndVal(:, :) = 0.0d0
        open(newunit=GenericIndValUnit, file=GenericIndValFile, status="unknown")
        do i = 1, Data%nInd
          read(GenericIndValUnit, *) IdCTmp, GenericIndValTmp(:)
          j = FindLoc(IdCTmp, Data%Coancestry%OriginalId(1:))
          if (j == 0) then
            write(STDERR, "(a)") "ERROR: Individual "//trim(IdCTmp)//" from the generic individual values file not present in the coancestry matrix file!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          Data%GenericIndVal(j, :) = GenericIndValTmp(:)
        end do
        close(GenericIndValUnit)

        allocate(Data%GenericIndValStat(Data%nGenericIndVal))
        do j = 1, Data%nGenericIndVal
          write(STDOUT, "(a)") "  - column "//trim(Int2Char(j))
          Data%GenericIndValStat(i) = DescStat(Data%GenericIndVal(:, j))
          write(STDOUT, "(a)") "    - average: "//trim(Real2Char(Data%GenericIndValStat(i)%Mean, fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(Data%GenericIndValStat(i)%SD,   fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(Data%GenericIndValStat(i)%Min,  fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(Data%GenericIndValStat(i)%Max,  fmt=FMTREAL2CHAR))
        end do
        write(STDOUT, "(a)") " "
      end if

      ! --- Generic mating values ---

      if (Spec%GenericMatValAvailable) then
        write(STDOUT, "(a)") "Generic mating values"
        DumI = CountLines(GenericMatValFile)
        if (DumI /= Data%nPotMat) then
          write(STDERR, "(a)") "ERROR: Number of matings in the generic mating values file and the number of all potential matings is not the same!"
          write(STDERR, "(a)") "ERROR: Number of all potential matings:                         "//trim(Int2Char(Data%nPotMat))
          write(STDERR, "(a)") "ERROR: Number of individuals in the generic mating values file: "//trim(Int2Char(DumI))
          write(STDERR, "(a)") " "
          stop 1
        end if
        allocate(Data%GenericMatVal(Data%nPotPar1, Data%nPotPar2, Data%nGenericMatVal))
        allocate(GenericMatValTmp(Data%nGenericMatVal))
        Data%GenericMatVal(:, :, :) = 0.0d0
        open(newunit=GenericMatValUnit, file=GenericMatValFile, status="unknown")
        do i = 1, Data%nPotMat
          read(GenericMatValUnit, *) IdCTmp, IdCTmp2, GenericMatValTmp(:)
          j = FindLoc(IdCTmp, Data%Coancestry%OriginalId(1:))
          if (j == 0) then
            write(STDERR, "(a)") "ERROR: Individual "//trim(IdCTmp)//" from the generic mating values file not present in the coancestry matrix file!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          k = FindLoc(IdCTmp2, Data%Coancestry%OriginalId(1:))
          if (k == 0) then
            write(STDERR, "(a)") "ERROR: Individual "//trim(IdCTmp2)//" from the generic mating values file not present in the coancestry matrix file!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          if (Spec%GenderMatters) then
            l = FindLoc(j, Data%IdPotPar1)
            if (l == 0) then
              write(STDERR, "(a)") "ERROR: Individual "//trim(IdCTmp)//" from the first column in the generic mating values file should be a male!"
              write(STDERR, "(a)") "ERROR: Generic mating values file (line "//trim(Int2Char(i))//"): "//trim(IdCTmp)//" "//trim(IdCTmp2)
              write(STDERR, "(a)") " "
              stop 1
            end if
            m = FindLoc(k, Data%IdPotPar2)
            if (l == 0) then
              write(STDERR, "(a)") "ERROR: Individual "//trim(IdCTmp2)//" from the second column in the generic mating values file should be a female!"
              write(STDERR, "(a)") "ERROR: Generic mating values file (line "//trim(Int2Char(i))//"): "//trim(IdCTmp)//" "//trim(IdCTmp2)
              write(STDERR, "(a)") " "
              stop 1
            end if
            ! fill full-matrix
            ! - l and m tell respectively the position within IdPotPar1 and IdPotPar2
            ! - values in IdPotPar1 and IdPotPar2 are "joint" Id of males and females
            ! - values in IdPotParSeq are separate Id of males and females (need these to find matching row and column)
            Data%GenericMatVal(Data%IdPotParSeq(Data%IdPotPar1(l)), Data%IdPotParSeq(Data%IdPotPar2(m)), :) = GenericMatValTmp(:)
          else
            ! fill lower-triangle (half-diallel)
            Data%GenericMatVal(maxval([j, k]), minval([j, k]), :) = GenericMatValTmp(:)
          end if
        end do
        close(GenericMatValUnit)

        allocate(Data%GenericMatValStat(Data%nGenericMatVal))
        do k = 1, Data%nGenericMatVal
          write(STDOUT, "(a)") "  - column "//trim(Int2Char(k))
          if (Spec%GenderMatters) then
            Data%GenericMatValStat(k) = DescStatMatrix(Data%GenericMatVal(:, :, k))
            write(STDOUT, "(a)") "    - average: "//trim(Real2Char(Data%GenericMatValStat(k)%All%Mean, fmt=FMTREAL2CHAR))
            write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(Data%GenericMatValStat(k)%All%SD,   fmt=FMTREAL2CHAR))
            write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(Data%GenericMatValStat(k)%All%Min,  fmt=FMTREAL2CHAR))
            write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(Data%GenericMatValStat(k)%All%Max,  fmt=FMTREAL2CHAR))
          else
            if (Spec%SelfingAllowed) then
              Data%GenericMatValStat(k) = DescStatLowTriMatrix(Data%GenericMatVal(:, :, k))
              write(STDOUT, "(a)") "    - average: "//trim(Real2Char(Data%GenericMatValStat(k)%All%Mean, fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(Data%GenericMatValStat(k)%All%SD,   fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(Data%GenericMatValStat(k)%All%Min,  fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(Data%GenericMatValStat(k)%All%Max,  fmt=FMTREAL2CHAR))
            end if
              Data%GenericMatValStat(k) = DescStatLowTriMatrix(Data%GenericMatVal(:, :, k), Diag=.false.)
              write(STDOUT, "(a)") "    - average: "//trim(Real2Char(Data%GenericMatValStat(k)%OffDiag%Mean, fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(Data%GenericMatValStat(k)%OffDiag%SD,   fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(Data%GenericMatValStat(k)%OffDiag%Min,  fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(Data%GenericMatValStat(k)%OffDiag%Max,  fmt=FMTREAL2CHAR))
          end if
        end do
        write(STDOUT, "(a)") " "
      end if
    end subroutine

    !###########################################################################

    !-------------------------------------------------------------------------
    !> @brief  Setup colnames and formats for output
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !-------------------------------------------------------------------------
    subroutine SetupColNamesAndFormats ! not pure due to setting module-wise variables
      implicit none
      integer(int32) :: nCol, nColTmp, i

      ! --- Optimisation log ---

      nCol = 10
      if (Spec%GenericIndValAvailable) then
        nCol = nCol + Data%nGenericIndVal
      end if
      if (Spec%GenericMatValAvailable) then
        nCol = nCol + Data%nGenericMatVal
      end if
      allocate(COLNAMELOGUNIT(nCol))
      allocate(COLNAMELOGSTDOUT(nCol))
      allocate(COLNAMELOGPOPUNIT(nCol))
      !                     1234567890123456789012
      COLNAMELOGUNIT(1)  = "             Iteration"
      COLNAMELOGUNIT(2)  = "            AcceptRate"
      COLNAMELOGUNIT(3)  = "          OptCriterion"
      COLNAMELOGUNIT(4)  = "             Penalties"
      COLNAMELOGUNIT(5)  = "          SelCriterion"
      COLNAMELOGUNIT(6)  = "          SelIntensity"
      COLNAMELOGUNIT(7)  = "            Coancestry"
      COLNAMELOGUNIT(8)  = "        CoancestryRate"
      COLNAMELOGUNIT(9)  = "            Inbreeding"
      COLNAMELOGUNIT(10) = "        InbreedingRate"
      nColTmp = nCol
      if (Spec%GenericIndValAvailable) then
        do i = 1, Data%nGenericIndVal
          nColTmp = nColTmp + 1
          COLNAMELOGUNIT(nColTmp) = "GenIndVal"//trim(Int2Char(i))
          COLNAMELOGUNIT(nColTmp) = adjustr(COLNAMELOGUNIT(nColTmp))
        end do
      end if
      if (Spec%GenericMatValAvailable) then
        do i = 1, Data%nGenericMatVal
          nColTmp = nColTmp + 1
          COLNAMELOGUNIT(nColTmp) = "GenMatVal"//trim(Int2Char(i))
          COLNAMELOGUNIT(nColTmp) = adjustr(COLNAMELOGUNIT(nColTmp))
        end do
      end if
      do i = 1, nCol
        COLNAMELOGSTDOUT(i) = COLNAMELOGUNIT(i)(9:22)
        COLNAMELOGSTDOUT(i) = adjustr(COLNAMELOGSTDOUT(i))
      end do
      COLNAMELOGPOPUNIT = COLNAMELOGUNIT
      !                       1234567890123456789012
      COLNAMELOGPOPUNIT(2) = "              Solution"
      FMTLOGSTDOUTHEAD  = trim(FMTLOGSTDOUTHEADA) //trim(Int2Char(nCol))  //trim(FMTLOGSTDOUTHEADB)
      FMTLOGSTDOUT      = trim(FMTLOGSTDOUTA)     //trim(Int2Char(nCol-1))//trim(FMTLOGSTDOUTB)
      FMTLOGUNITHEAD    = trim(FMTLOGUNITHEADA)   //trim(Int2Char(nCol)  )//trim(FMTLOGUNITHEADB)
      FMTLOGUNIT        = trim(FMTLOGUNITA)       //trim(Int2Char(nCol-1))//trim(FMTLOGUNITB)
      FMTLOGPOPUNIT     = trim(FMTLOGPOPUNITA)    //trim(Int2Char(nCol-2))//trim(FMTLOGUNITB)
    end subroutine

    !###########################################################################

    !-------------------------------------------------------------------------
    !> @brief  Call various optimisations for AlphaMate
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !-------------------------------------------------------------------------
    subroutine AlphaMateSearch ! not pure due to IO

      implicit none

      integer(int32) :: nParam, k, FrontierUnit

      real(real64) :: HoldTargetCoancestryRanMate, HoldTargetCoancestryRate
      real(real64), allocatable :: InitEqual(:, :)

      character(len=1000) :: LogFile, LogPopFile, ContribFile, MatingFile
      character(len=100) :: DumC

      type(AlphaMateSol) :: SolMin, SolRan, SolOpt, Sol

      ! --- Number of parameters to optimise ---

      if (Spec%GenderMatters) then
        nParam = Data%nPotPar1 + Data%nPotPar2 + Data%nMat
      else
        nParam = Data%nPotPar1 + Data%nMat
      end if

      if (Spec%PAGE) then
        nParam = nParam + Data%nInd
      end if

      ! --- Optimise contributions for minimum future coancestry/inbreeding ---

      if (Spec%ModeMin) then
        write(STDOUT, "(a)") "--- Optimise contributions for minimum future coancestry/inbreeding --- "
        write(STDOUT, "(a)") " "

        LogFile     = "OptimisationLogMinimumInbreeding.txt"
        LogPopFile  = "OptimisationLogPopMinimumInbreeding.txt"
        ContribFile = "IndividualResultsMinimumInbreeding.txt"
        MatingFile  = "MatingResultsMinimumInbreeding.txt"

        allocate(InitEqual(nParam, nint(Spec%EvolAlgNSol * 0.1)))
        InitEqual(:, :) = 1.0d0 ! A couple of solutions that would give equal contributions to everybody

        call DifferentialEvolution(nParam=nParam, nSol=Spec%EvolAlgNSol, Init=InitEqual, &
          nIter=Spec%EvolAlgNIter, nIterBurnIn=Spec%EvolAlgNIterBurnIn, nIterStop=Spec%EvolAlgNIterStop, StopTolerance=Spec%EvolAlgStopTol, nIterPrint=Spec%EvolAlgNIterPrint, &
          LogFile=LogFile, LogPop=Spec%EvolAlgLogPop, LogPopFile=LogPopFile, &
          CritType="min", CRBurnIn=Spec%EvolAlgCRBurnIn, CRLate=Spec%EvolAlgCRLate, FBase=Spec%EvolAlgFBase, FHigh1=Spec%EvolAlgFHigh1, FHigh2=Spec%EvolAlgFHigh2, &
          BestSol=SolMin)

        deallocate(InitEqual)

        call SolMin%Write(Data, Spec, ContribFile, MatingFile)
      end if

      ! --- Evaluate random mating ---

      if (Spec%ModeRan) then
        write(STDOUT, "(a)") "--- Evaluate random mating --- "
        write(STDOUT, "(a)") " "

        LogFile = "OptimisationLogRandomMating.txt"

        allocate(InitEqual(nParam, nint(Spec%EvolAlgNSol * 0.1)))
        InitEqual(:, :) = 1.0d0 ! A couple of solutions that would give equal contributions for everybody

        call RandomSearch(Mode="avg", nParam=nParam, Init=InitEqual, nSamp=Spec%EvolAlgNSol*Spec%EvolAlgNIter*Spec%RanAlgStricter, nSampStop=Spec%EvolAlgNIterStop*Spec%RanAlgStricter, &
          StopTolerance=Spec%EvolAlgStopTol/Spec%RanAlgStricter, nSampPrint=Spec%EvolAlgNIterPrint, LogFile=LogFile, CritType="ran", BestSol=SolRan)

        deallocate(InitEqual)
      end if

      ! --- Optimise contributions for maximum value with constraint on future coancestry/inbreeding ---

      if (Spec%ModeOpt) then
        write(STDOUT, "(a)") "--- Optimise contributions for maximum value with constraint on future coancestry/inbreeding ---"
        write(STDOUT, "(a)") " "

        LogFile     = "OptimisationLogMaximumValue.txt"
        LogPopFile  = "OptimisationLogPopMaximumValue.txt"
        ContribFile = "IndividualResultsMaximumValue.txt"
        MatingFile  = "MatingResultsMaximumValue.txt"

        ! @todo add some clever initial values, say:
        !       - equal contributions for top 2/3 or 1/2 of BV distribution,
        !       - decreasing contributions with decreasing value
        !       - SDP solution, ...?
        call DifferentialEvolution(nParam=nParam, nSol=Spec%EvolAlgNSol, nIter=Spec%EvolAlgNIter, nIterBurnIn=Spec%EvolAlgNIterBurnIn, &
          nIterStop=Spec%EvolAlgNIterStop, StopTolerance=Spec%EvolAlgStopTol, nIterPrint=Spec%EvolAlgNIterPrint,&
          LogFile=LogFile, LogPop=Spec%EvolAlgLogPop, LogPopFile=LogPopFile, &
          CritType="opt", CRBurnIn=Spec%EvolAlgCRBurnIn, CRLate=Spec%EvolAlgCRLate, FBase=Spec%EvolAlgFBase, FHigh1=Spec%EvolAlgFHigh1, FHigh2=Spec%EvolAlgFHigh2, &
          BestSol=SolOpt)

        call SolOpt%Write(Data, Spec, ContribFile, MatingFile)
      end if

      ! --- Evaluate the full frontier ---

      if (Spec%EvaluateFrontier) then
        write(STDOUT, "(a)") "--- Evaluate the frontier ---"
        write(STDOUT, "(a)") " "

        Spec%CoancestryWeightBellow = .true. ! we want to target certain rates of coancestry

        open(newunit=FrontierUnit, file="Frontier.txt", status="unknown")
        call Sol%LogHead(FrontierUnit, String="Bla", StringNum=10)
        if (Spec%ModeMin) then
          DumC = "Min"
          call SolMin%Log(FrontierUnit, Iteration=-1, AcceptRate=-1.0d0, String=DumC, StringNum=10)
        end if
        if (Spec%ModeRan) then
          DumC = "Ran"
          call SolRan%Log(FrontierUnit, Iteration=-1, AcceptRate=-1.0d0, String=DumC, StringNum=10)
        end if
        if (Spec%ModeOpt) then
          DumC = "Opt"
          call SolOpt%Log(FrontierUnit, Iteration=-1, AcceptRate=-1.0d0, String=DumC, StringNum=10)
        end if

        ! Hold old results
        HoldTargetCoancestryRanMate = Data%TargetCoancestryRanMate
        HoldTargetCoancestryRate = Spec%TargetCoancestryRate

        ! Evaluate
        do k = 1, Spec%nFrontierSteps
          Spec%TargetCoancestryRate = Spec%TargetCoancestryRateFrontier(k)
          ! F_t = DeltaF + (1 - DeltaF) * F_t-1
          Data%TargetCoancestryRanMate = Spec%TargetCoancestryRate + (1.0d0 - Spec%TargetCoancestryRate) * Data%CurrentCoancestryRanMate
          write(STDOUT, "(a)") "Step "//trim(Int2Char(k))//" out of "//trim(Int2Char(Spec%nFrontierSteps))//&
                               " for the rate of coancestry "//trim(Real2Char(Spec%TargetCoancestryRate, fmt=FMTREAL2CHAR))//&
                               " (=targeted coancestry "//trim(Real2Char(Data%TargetCoancestryRanMate, fmt=FMTREAL2CHAR))//")"
          write(STDOUT, "(a)") ""

          LogFile     = "OptimisationLogFrontier"//trim(Int2Char(k))//".txt"
          LogPopFile  = "OptimisationLogPopFrontier"//trim(Int2Char(k))//".txt"
          ContribFile = "IndividualResultsFrontier"//trim(Int2Char(k))//".txt"
          MatingFile  = "MatingResultsFrontier"//trim(Int2Char(k))//".txt"

          call DifferentialEvolution(nParam=nParam, nSol=Spec%EvolAlgNSol, nIter=Spec%EvolAlgNIter, nIterBurnIn=Spec%EvolAlgNIterBurnIn, &
            nIterStop=Spec%EvolAlgNIterStop, StopTolerance=Spec%EvolAlgStopTol, nIterPrint=Spec%EvolAlgNIterPrint,&
            LogFile=LogFile, LogPop=Spec%EvolAlgLogPop, LogPopFile=LogPopFile, &
            CritType="opt", CRBurnIn=Spec%EvolAlgCRBurnIn, CRLate=Spec%EvolAlgCRLate, FBase=Spec%EvolAlgFBase, FHigh1=Spec%EvolAlgFHigh1, FHigh2=Spec%EvolAlgFHigh2, &
            BestSol=Sol)

          DumC = "Frontier"//trim(Int2Char(k))
          call Sol%Log(FrontierUnit, Iteration=-1, AcceptRate=-1.0d0, String=DumC, StringNum=10)
          call Sol%Write(Data, Spec, ContribFile, MatingFile)

          if ((Spec%TargetCoancestryRate - Sol%CoancestryRateRanMate) > 0.01d0) then
            write(STDOUT, "(a)") "NOTE: Could not achieve the rate of coancestry "//trim(Real2Char(Spec%TargetCoancestryRate, fmt=FMTREAL2CHAR))
            write(STDOUT, "(a)") "NOTE: Stopping the frontier evaluation."
            write(STDOUT, "(a)") ""
            exit
          end if
        end do

        ! Put back old results
        Spec%TargetCoancestryRate = HoldTargetCoancestryRate
        Data%TargetCoancestryRanMate = HoldTargetCoancestryRanMate

        close(FrontierUnit)

      end if
    end subroutine

    !###########################################################################

    !-------------------------------------------------------------------------
    !> @brief  Write AlphaMate solution to files or standard output
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !-------------------------------------------------------------------------
    subroutine WriteAlphaMateSol(This, Data, Spec, ContribFile, MatingFile) ! not pure due to IO
      implicit none

      ! Arguments
      class(AlphaMateSol)          :: This        !< AlphaMateSol holder
      type(AlphaMateData)          :: Data        !< AlphaMateData holder
      type(AlphaMateSpec)          :: Spec        !< AlphaMateSpec holder
      character(len=*), optional   :: ContribFile !< File to write individual contributions to (default STDOUT)
      character(len=*), optional   :: MatingFile  !< File to write mating plan to (default STDOUT)

      ! Other
      integer(int32) :: i, j, ContribUnit, MatingUnit, Rank(Data%nInd)
      if (.not. present(ContribFile)) then
        ContribUnit = STDOUT
      end if
      if (.not. present(MatingFile)) then
        MatingUnit = STDOUT
      end if

      ! @todo should we have constant output no matter which options are switched on?
      if (present(ContribFile)) then
        open(newunit=ContribUnit, file=ContribFile, status="unknown")
      end if
      Rank = MrgRnk(This%nVec + Data%SelCriterionStand / 100.0d0) ! @todo is this really good sorting?
      if (.not.Spec%PAGE) then
        !                                        1234567890123456789012
        write(ContribUnit, FMTCONTRIBUTIONHEAD) "          Id", &
                                                "      Gender", &
                                                "SelCriterion", &
                                                "AvCoancestry", &
                                                "Contribution", &
                                                "     nMating"
        do i = Data%nInd, 1, -1 ! MrgRnk ranks small to large
          j = Rank(i)
          write(ContribUnit, FMTCONTRIBUTION) Data%Coancestry%OriginalId(j), Data%Gender(j), Data%SelCriterion(j), &
                                              sum(Data%Coancestry%Value(1:, j)) / Data%nInd, &
                                              This%xVec(j), This%nVec(j)
        end do
      else
        !                                            1234567890123456789012
        write(ContribUnit, FMTCONTRIBUTIONHEADEDIT) "          Id", &
                                                    "      Gender", &
                                                    "SelCriterion", &
                                                    "AvCoancestry", &
                                                    "Contribution", &
                                                    "     nMating", &
                                                    "  GenomeEdit", &
                                                    " EditedValue"
        do i = Data%nInd, 1, -1 ! MrgRnk ranks small to large
          j = Rank(i)
          write(ContribUnit, FMTCONTRIBUTIONEDIT) Data%Coancestry%OriginalId(j), Data%Gender(j), Data%SelCriterion(j), &
                                                  sum(Data%Coancestry%Value(1:, j)) / Data%nInd, &
                                                  This%xVec(j), This%nVec(j), &
                                                  nint(This%GenomeEdit(j)), Data%SelCriterion(j) + This%GenomeEdit(j) * Data%SelCriterionPAGE(j)
        end do
      end if
      if (present(ContribFile)) then
        close(ContribUnit)
      end if

      if (present(MatingFile)) then
        open(newunit=MatingUnit, file=MatingFile, status="unknown")
      end if
      !                                 1234567890123456789012
      write(MatingUnit, FMTMATINGHEAD) "      Mating", &
                                       "     Parent1", &
                                       "     Parent2"
      do i = 1, Data%nMat
        write(MatingUnit, FMTMATING) i, Data%Coancestry%OriginalId(This%MatingPlan(1, i)), Data%Coancestry%OriginalId(This%MatingPlan(2, i))
      end do
      if (present(MatingFile)) then
        close(MatingUnit)
      end if
    end subroutine

    !###########################################################################

    !-------------------------------------------------------------------------
    !> @brief  Initialize AlphaMate solution
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !-------------------------------------------------------------------------
    pure subroutine InitialiseAlphaMateSol(This)
      implicit none

      ! Argument
      class(AlphaMateSol), intent(out) :: This

      ! Initialisation
      This%Criterion = 0.0d0
      This%Penalty = 0.0d0
      This%SelCriterion = 0.0d0
      This%SelIntensity = 0.0d0
      This%FutureCoancestryRanMate = 0.0d0
      This%CoancestryRateRanMate = 0.0d0
      This%FutureInbreeding = 0.0d0
      This%InbreedingRate = 0.0d0
      if (Spec%GenericIndValAvailable) then
        allocate(This%GenericIndVal(Data%nGenericIndVal))
        This%GenericIndVal(:) = 0.0d0
      else
        allocate(This%GenericIndVal(0))
      end if
      if (Spec%GenericMatValAvailable) then
        allocate(This%GenericMatVal(Data%nGenericMatVal))
        This%GenericMatVal(:) = 0.0d0
      else
        allocate(This%GenericMatVal(0))
      end if
      This%Cost = 0.0d0
      allocate(This%nVec(Data%nInd))
      This%nVec(:) = 0
      allocate(This%xVec(Data%nInd))
      This%xVec(:) = 0.0d0
      allocate(This%MatingPlan(2, Data%nMat))
      This%MatingPlan(:, :) = 0
      if (Spec%PAGE) then
        allocate(This%GenomeEdit(Data%nInd))
        This%GenomeEdit(:) = 0.0d0
      else
        allocate(This%GenomeEdit(0))
      end if
    end subroutine

    !###########################################################################

    !-------------------------------------------------------------------------
    !> @brief  Assign one AlphaMate solution to another
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !-------------------------------------------------------------------------
    pure subroutine AssignAlphaMateSol(Out, In)
      implicit none

      ! Arguments
      class(AlphaMateSol), intent(out)   :: Out
      class(AlphaEvolveSol), intent(in)  :: In

      ! Assignments
      ! (Need to go via the select type stuff as all but the first arguments must
      !  be the same as in the base class/type)
      select type (In)
        class is (AlphaMateSol)
          Out%Criterion = In%Criterion
          Out%Penalty = In%Penalty
          Out%SelCriterion = In%SelCriterion
          Out%SelIntensity = In%SelIntensity
          Out%FutureCoancestryRanMate = In%FutureCoancestryRanMate
          Out%CoancestryRateRanMate = In%CoancestryRateRanMate
          Out%FutureInbreeding = In%FutureInbreeding
          Out%InbreedingRate = In%InbreedingRate
          if (allocated(In%GenericIndVal)) then
            allocate(Out%GenericIndVal(size(In%GenericIndVal)))
            Out%GenericIndVal = In%GenericIndVal
          end if
          if (allocated(In%GenericMatVal)) then
            allocate(Out%GenericMatVal(size(In%GenericMatVal)))
            Out%GenericMatVal = In%GenericMatVal
          end if
          Out%Cost = In%Cost
          if (allocated(In%nVec)) then
            allocate(Out%nVec(size(In%nVec)))
            Out%nVec = In%nVec
          end if
          if (allocated(In%xVec)) then
            allocate(Out%xVec(size(In%xVec)))
            Out%xVec = In%xVec
          end if
          if (allocated(In%MatingPlan)) then
            allocate(Out%MatingPlan(size(In%MatingPlan, dim=1), size(In%MatingPlan, dim=2)))
            Out%MatingPlan = In%MatingPlan
          end if
          if (allocated(In%GenomeEdit)) then
            allocate(Out%GenomeEdit(size(In%GenomeEdit)))
            Out%GenomeEdit = In%GenomeEdit
          end if
      end select
    end subroutine

    !###########################################################################

    !-------------------------------------------------------------------------
    !> @brief  Update mean of AlphaMate solution (when performing random search)
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !-------------------------------------------------------------------------
    pure subroutine UpdateMeanAlphaMateSol(This, Add, n)
      implicit none

      ! Arguments
      class(AlphaMateSol), intent(inout) :: This
      class(AlphaEvolveSol), intent(in)  :: Add
      integer(int32), intent(in)         :: n

      ! Other
      real(real64) :: kR

      ! Updates
      kR = (dble(n) - 1.0d0) / n

      ! (Need to go via the select type stuff as all but the first arguments must
      !  be the same as in the base class/type)
      select type (Add)
        class is (AlphaMateSol)
          This%Criterion                 = This%Criterion                 * kR + Add%Criterion                 / n
          This%Penalty                   = This%Penalty                   * kR + Add%Penalty                   / n
          This%SelCriterion              = This%SelCriterion              * kR + Add%SelCriterion              / n
          This%SelIntensity              = This%SelIntensity              * kR + Add%SelIntensity              / n
          This%FutureCoancestryRanMate   = This%FutureCoancestryRanMate   * kR + Add%FutureCoancestryRanMate   / n
          This%CoancestryRateRanMate     = This%CoancestryRateRanMate     * kR + Add%CoancestryRateRanMate     / n
          This%FutureInbreeding          = This%FutureInbreeding          * kR + Add%FutureInbreeding          / n
          This%InbreedingRate            = This%InbreedingRate            * kR + Add%InbreedingRate            / n
          if (allocated(This%GenericIndVal)) then
            This%GenericIndVal           = This%GenericIndVal             * kR + Add%GenericIndVal             / n
          end if
          if (allocated(This%GenericMatVal)) then
            This%GenericMatVal           = This%GenericMatVal             * kR + Add%GenericMatVal             / n
          end if
          This%Cost                      = This%Cost                      * kR + Add%Cost                      / n
          if (allocated(This%nVec)) then
            This%nVec                    = This%nVec                      * kR + Add%nVec                      / n
          end if
          if (allocated(This%xVec)) then
            This%xVec                    = This%xVec                      * kR + Add%xVec                      / n
          end if
          if (allocated(This%MatingPlan)) then
            This%MatingPlan              = This%MatingPlan                * kR + Add%MatingPlan                / n
          end if
          if (allocated(This%GenomeEdit)) then
            This%GenomeEdit              = This%GenomeEdit                * kR + Add%GenomeEdit                / n
          end if
      end select
    end subroutine

    !###########################################################################

    !-------------------------------------------------------------------------
    !> @brief  AlphaMate evaluate function plus much MORE (this is the core!!!!)
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !-------------------------------------------------------------------------
    subroutine FixSolEtcMateAndEvaluate(This, Chrom, CritType)
      implicit none
      ! Arguments
      class(AlphaMateSol), intent(inout)     :: This      ! Solution
      real(real64), intent(inout), optional  :: Chrom(:)  ! Internal representation of the solution
      character(len=*), intent(in), optional :: CritType  ! Type of criterion (Min, Ran, Opt)

      ! Other
      integer(int32) :: i, j, k, l, g, nCumMat, Rank(Data%nInd), ChromInt(Data%nInd), MatPar2(Data%nMat)
      integer(int32) :: nVecPar1(Data%nPotPar1), nVecPar2(Data%nPotPar2), TmpMin, TmpMax, TmpI

      real(real64) :: TmpVec(Data%nInd, 1), TmpR, RanNum

      ! Initialize the solution
      call This%Initialise()

      ! The solution (based on the mate selection driver) has:
      ! - Data%nInd individual contributions
      !   - Data%nPotPar1 individual contributions for "parent1" (males   when GenderMatters, all ind when .not. GenderMatters)
      !   - Data%nPotPar2 individual contributions for "parent2" (females when GenderMatters, meaningful only when GenderMatters)
      ! - Data%nMat     rankings of parent1 1:Data%nMat matings to pair with 1:Data%nPotPar2 "parent2" (see bellow)
      ! - Data%nInd edit indicators
      !   - Data%nPotPar1 edit indicators for "parent1" (males   when GenderMatters, all ind when .not. GenderMatters)
      !   - Data%nPotPar2 edit indicators for "parent2" (females when GenderMatters, present only when GenderMatters)

      ! Say we have Chrom=(| 0, 2, 0, 1 | ... | 2.5, 1.5, 1.0 | 0, 1, 0, 0 | ...) then we:
      ! - mate male 2 with the first  available female (rank 2.5)
      ! - mate male 2 with the second available female (rank 1.5)
      ! - mate male 4 with the third  available female (rank 1.0)
      ! - edit male 2

      ! @todo consider spliting the Chrom() vector internally into a type with
      !       separate vectors to simplify the code, e.g.,
      ! Chrom2%ContPar1
      ! Chrom2%ContPar2
      ! Chrom2%MateRank
      ! Chrom2%EditPar1
      ! Chrom2%EditPar2
      !       and then at the end combine it back. Since I modify some elements
      !       it would have to be put back.

      ! --- Parse the mate selection driver (=Is the solution valid?) ---

      ! The approach below assures that we have Data%nMat contributions for each of
      ! the two parent sets. It does this by ranking internal solution values and
      ! traverses from top to the defined number of parents checking when the sum of
      ! interegrised values gives Data%nMat. If values bellow 0.5 are found, they are
      ! changed to 1 contribution. If this still does not give Data%nMat, then we start

      ! @todo least contributing? Those that already contribute or thos that do not contribute at all?

      ! adding on contribution to each parent (starting from least contributing
      ! parents to avoid local minima) until we reach Data%nMat. How to treat values
      ! for the individuals that do not contribute is unlcear. None of the tested
      ! methods seemed to be very different. Intuitively, using properly ordered
      ! negative values should inform optim. alg. which individuals should less
      ! likely contribute, but this did not seem to be the case - better final
      ! solution was found when this strategy was not implemented - either zeroing
      ! values for those individuals (was the fastest) or giving random value (was
      ! marginally better, but slower). Potential advantage of not preserving the
      ! order is that this gives more randomness and more solutions being explored.

      ! "Parent1"
      if (Spec%GenderMatters) then
        g = 1
      else
        g = 2
      end if
      ! ... find ranks to find the top values
      if (.not.(Spec%EqualizePar1 .and. (Data%nPar1 == Data%nPotPar1))) then
        Rank(1:Data%nPotPar1) = MrgRnk(Chrom(1:Data%nPotPar1))
        Rank(1:Data%nPotPar1) = Rank(Data%nPotPar1:1:-1) ! MrgRnk ranks small to large
        !@todo decreasing option in MrgRnk?
      end if
      ! ... handle cases with equalized contributions
      if (Spec%EqualizePar1) then
        if (Data%nPar1 == Data%nPotPar1) then
          ! ... set integers to all the values (no need for sorting here)
          Chrom(1:Data%nPotPar1) = dble(Data%nMat * g) / Data%nPar1
        else
          ! ... set integers to the top values
          Chrom(Rank(1:Data%nPar1)) = dble(Data%nMat * g) / Data%nPar1
          !@todo anything better to preserve the order of non contributing individuals? See below!
          Chrom(Rank((Data%nPar1+1):Data%nPotPar1)) = 0.0d0
          ! Chrom(Rank((Data%nPar1+1):Data%nPotPar1)) = -1.0d0
        end if
      else
        ! ... handle cases with unequal contributions
        ! ... work for the defined number or parents
        nCumMat = 0
        do i = 1, Data%nPar1
          j = Rank(i)
          ! ... these would have "zero" contributions and if we hit them then they need at least minimum contrib
          if (Chrom(j) < Spec%LimitPar1Min) then
            Chrom(j) = Spec%LimitPar1Min
          end if
          ! ... but not above max allowed
          if (Chrom(j) > Spec%LimitPar1Max) then
            Chrom(j) = Spec%LimitPar1Max
          end if
          ! ... accumulate and check if we reached Data%nMat
          nCumMat = nCumMat + nint(Chrom(j)) ! internally real, externally integer
          if (nCumMat >= Data%nMat * g) then
            ! ... there should be exactly Data%nMat contributions
            if (nCumMat > Data%nMat * g) then
              Chrom(j) = Chrom(j) - dble(nCumMat - Data%nMat * g)
              if (nint(Chrom(j)) < Spec%LimitPar1Min) then
                TmpR = Spec%LimitPar1Weight * (Spec%LimitPar1Min - nint(Chrom(j)))
                This%Criterion = This%Criterion + TmpR
                if (Spec%LimitPar1Weight < 0.0d0) then
                  This%Penalty = This%Penalty + abs(TmpR)
                end if
              end if
              nCumMat = Data%nMat * g
            end if
            exit
          end if
        end do
        ! ... increment i if we have hit the exit, do loop would have ended with i=Data%nPar1 + 1
        if (i <= Data%nPar1) then
          i = i + 1
        end if
        ! ... the other individuals do not contribute
        if (i <= Data%nPotPar1) then ! "=" to capture incremented i+1 on the do loop exit
          ! ... zero (the same for all ind so no order)
          Chrom(Rank(i:Data%nPotPar1)) = 0.0d0
          ! ... negative (the same for all ind so no order)
          ! Chrom(Rank(i:Data%nPotPar1)) = -1.0d0
          ! ... negative (variable with partially preserving order)
          !     Found faster convergence than with properly decreasing negative values?
          !     I guess it adds some more randomness, i.e., it causes suffling of individuals that just did not make it onto the mating list.
          ! Chrom(Rank(i:Data%nPotPar1)) = sign(Chrom(Rank(i:Data%nPotPar1)), -1.0d0)
          ! ... negative and properly decreasing
          ! TmpR = maxval(Chrom(Rank(i:Data%nPotPar1)))
          ! if (TmpR > 0.0d0) then
          !     Chrom(Rank(i:Data%nPotPar1)) = Chrom(Rank(i:Data%nPotPar1)) - abs(TmpR)
          ! end if
          ! ... negative (random so no order)
          ! do j = i, Data%nPotPar1 ! @todo really need this loop?
          !   call random_number(RanNum)
          !   Chrom(Rank(j)) = -1.0d0 * RanNum
          ! end do
        end if
        ! ... Data%nMat still not reached?
        do while (nCumMat < Data%nMat * g)
          ! ... add more contributions
          do i = Data%nPar1, 1, -1 ! to bottom ranked selected individuals (to avoid local optima)
            j = Rank(i)
            Chrom(j) = Chrom(j) + 1.0d0
            ! ... accumulate and check if we reached Data%nMat
            nCumMat = nCumMat + 1
            if (nCumMat >= Data%nMat * g) then
              ! To cater for real vs. integer issues
              TmpI = sum(nint(Chrom(Rank(1:Data%nPar1))))
              if (TmpI /= Data%nMat * g) then
                if (TmpI > Data%nMat * g) then
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
      if (Spec%GenderMatters) then
        ! ... find ranks to find the top values
        if (.not.(Spec%EqualizePar2 .and. (Data%nPar2 == Data%nPotPar2))) then
          Rank(1:Data%nPotPar2) = MrgRnk(Chrom((Data%nPotPar1+1):(Data%nPotPar1+Data%nPotPar2)))
          Rank(1:Data%nPotPar2) = Rank(Data%nPotPar2:1:-1) ! MrgRnk ranks small to large
        end if
        ! ... handle cases with equalized contributions
        if (Spec%EqualizePar2) then
          if (Data%nPar2 == Data%nPotPar2) then
            ! ... set integers to all the values (no need for sorting here)
            Chrom((Data%nPotPar1+1):(Data%nPotPar1+Data%nPotPar2)) = dble(Data%nMat) / Data%nPar2
          else
            ! ... set integers to the top values
            Chrom(Data%nPotPar1+Rank(1:Data%nPar2)) = dble(Data%nMat) / Data%nPar2
            ! @todo anything better to preserve the order of non contributing individuals? See below!
            Chrom(Data%nPotPar1+Rank((Data%nPar2+1):Data%nPotPar2)) = 0.0d0
            ! Chrom(Data%nPotPar1+Rank((Data%nPar2+1):Data%nPotPar2)) = -1.0d0
          end if
        else
          ! ... handle cases with unequal contributions
          ! ... work for the defined number or parents
          nCumMat = 0
          do i = 1, Data%nPar2
            j = Data%nPotPar1 + Rank(i)
            ! ... these would have "zero" contributions and if we hit them then they need at least minimum contrib
            if (Chrom(j) < Spec%LimitPar2Min) then
              Chrom(j) = Spec%LimitPar2Min
            end if
            ! ... but not above max allowed
            if (Chrom(j) > Spec%LimitPar2Max) then
              Chrom(j) = Spec%LimitPar2Max
            end if
            ! ... accumulate and check if we reached Data%nMat
            nCumMat = nCumMat + nint(Chrom(j)) ! internally real, externally integer
            if (nCumMat >= Data%nMat) then
              ! ... there should be exactly Data%nMat contributions
              if (nCumMat > Data%nMat) then
                Chrom(j) = Chrom(j) - dble(nCumMat - Data%nMat)
                if (nint(Chrom(j)) < Spec%LimitPar2Min) then
                  TmpR = Spec%LimitPar2Weight * (Spec%LimitPar2Min - nint(Chrom(j)))
                  This%Criterion = This%Criterion + TmpR
                  if (Spec%LimitPar2Weight < 0.0d0) then
                    This%Penalty = This%Penalty + abs(TmpR)
                  end if
                end if
                nCumMat = Data%nMat
              end if
              exit
            end if
          end do
          ! ... increment i if we have hit the exit, do loop would have ended with i=Data%nPar2+1
          if (i <= Data%nPar2) then
            i = i + 1
          end if
          ! ... the other individuals do not contribute
          if (i <= Data%nPotPar2) then ! "="" to capture incremented i+1 on the do loop exit
            ! ... zero (the same for all ind so no order)
            Chrom(Data%nPotPar1+(Rank(i:Data%nPotPar2))) = 0.0d0
            ! ... negative (the same for all ind so no order)
            ! Chrom(Data%nPotPar1+(Rank(i:Data%nPotPar2))) = -1.0d0
            ! ... negative (variable with partially preserving order, i.e., ~large positives become ~large negatives)
            !     Found faster convergence than with properly decreasing negative values?
            !     I guess it adds some more randomness, i.e., it causes suffling of individuals that just did not make it onto the mating list.
            ! Chrom(Data%nPotPar1+(Rank(i:Data%nPotPar2))) = sign(Chrom(Data%nPotPar1+(Rank(i:Data%nPotPar2))), -1.0d0)
            ! ... negative and properly decreasing
            ! TmpR = maxval(Chrom(Data%nPotPar1+(Rank(i:Data%nPotPar2))))
            ! if (TmpR > 0.0d0) then
            !     Chrom(Data%nPotPar1+(Rank(i:Data%nPotPar2))) = Chrom(Data%nPotPar1+(Rank(i:Data%nPotPar2))) - abs(TmpR)
            ! end if
            ! ... negative (random so no order)
            ! do j = i, Data%nPotPar2 ! @todo really need this loop?
            !   call random_number(RanNum)
            !   Chrom(Data%nPotPar1+Rank(j)) = -1.0d0 * RanNum
            ! end do
          end if
          ! ... Data%nMat still not reached?
          do while (nCumMat < Data%nMat)
            ! ... add more contributions
            do i = Data%nPar2, 1, -1 ! to bottom ranked selected individuals (to avoid local optima)
              j = Data%nPotPar1 + Rank(i)
              Chrom(j) = Chrom(j) + 1.0d0
              ! ... accumulate and check if we reached Data%nMat
              nCumMat = nCumMat + 1
              if (nCumMat == Data%nMat) then
                ! To cater for real vs. integer issues
                TmpI = sum(nint(Chrom(Data%nPotPar1+Rank(1:Data%nPar2))))
                if (TmpI /= Data%nMat) then
                  if (TmpI > Data%nMat) then
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

      ! --- Contributions (nVec & xVec) ---

      nVecPar1(:) = 0

      ! "Parent1"
      ! ... get integer values
      ChromInt(1:Data%nPotPar1) = nint(Chrom(1:Data%nPotPar1))
      ! ... remove negatives
      do i = 1, Data%nPotPar1
        if (ChromInt(i) < 0) then
          ChromInt(i) = 0
        end if
      end do
      ! ... map internal to external order
      nVecPar1(:) = ChromInt(1:Data%nPotPar1)
      if (.not.Spec%GenderMatters) then
        This%nVec(:) = nVecPar1(:)
      else
        This%nVec(Data%IdPotPar1) = nVecPar1(:)
      end if

      ! "Parent2"
      if (Spec%GenderMatters) then
        nVecPar2(:) = 0
        ! ... get integer values
        ChromInt(1:Data%nPotPar2) = nint(Chrom((Data%nPotPar1+1):(Data%nPotPar1+Data%nPotPar2)))
        ! ... remove negatives
        do i = 1, Data%nPotPar2
          if (ChromInt(i) < 0) then
            ChromInt(i) = 0
          end if
        end do
        ! ... map internal to external order
        nVecPar2(:) = ChromInt(1:Data%nPotPar2)
        This%nVec(Data%IdPotPar2) = nVecPar2(:)
      end if

      This%xVec(:) = dble(This%nVec(:)) / (2 * Data%nMat)

      ! --- PAGE ---

      if (Spec%PAGE) then
        if (.not.Spec%GenderMatters) then
          Rank(1:Data%nInd) = MrgRnk(Chrom((Data%nPotPar1+Data%nMat+1):(Data%nPotPar1+Data%nMat+Data%nInd)))
          This%GenomeEdit(Rank(Data%nInd:(Data%nInd-Spec%PAGEPar1Max+1):-1)) = 1.0d0 ! MrgRnk ranks small to large
        else
          if (Spec%PAGEPar1) then
            Rank(1:Data%nPotPar1) = MrgRnk(Chrom((Data%nPotPar1+Data%nPotPar2+Data%nMat+1):(Data%nPotPar1+Data%nPotPar2+Data%nMat+Data%nPotPar1)))
            This%GenomeEdit(Data%IdPotPar1(Rank(Data%nPotPar1:(Data%nPotPar1-Spec%PAGEPar1Max+1):-1))) = 1.0d0 ! MrgRnk ranks small to large
          end if
          if (Spec%PAGEPar2) then
            Rank(1:Data%nPotPar2) = MrgRnk(Chrom((Data%nPotPar1+Data%nPotPar2+Data%nMat+Data%nPotPar1+1):(Data%nPotPar1+Data%nPotPar2+Data%nMat+Data%nPotPar1+Data%nPotPar2)))
            This%GenomeEdit(Data%IdPotPar2(Rank(Data%nPotPar2:(Data%nPotPar2-Spec%PAGEPar2Max+1):-1))) = 1.0d0 ! MrgRnk ranks small to large
          end if
        end if
      end if

      ! --- Mate allocation ---

      MatPar2(:) = 0
      if (Spec%GenderMatters) then
        ! Distribute parent2 (=female) contributions into matings
        k = 0
        do i = 1, Data%nPotPar2 ! need to loop whole nVecPar2 as some entries are zero
          do j = 1, nVecPar2(i)
            k = k + 1
            MatPar2(k) = Data%IdPotPar2(i)
          end do
        end do
        ! Reorder parent2 contributions according to the rank of matings
        Rank(1:Data%nMat) = MrgRnk(Chrom((Data%nPotPar1+Data%nPotPar2+1):(Data%nPotPar1+Data%nPotPar2+Data%nMat)))
        MatPar2(:) = MatPar2(Rank(1:Data%nMat))
      else
        ! Distribute one half of contributions into matings
        k = 0
        do while (k < Data%nMat)
          do i = 1, Data%nPotPar1 ! need to loop whole nVecPar1 as some entries are zero
            l = nVecPar1(i) / 2
            if (mod(nVecPar1(i), 2) == 1) then
              call random_number(RanNum)
              if (RanNum > 0.5) then
                l = l + 1
              end if
            end if
            do j = 1, l
              if (k == Data%nMat) then
                exit
              end if
              k = k + 1
              MatPar2(k) = Data%IdPotPar1(i)
              nVecPar1(i) = nVecPar1(i) - 1
            end do
          end do
        end do
        ! Reorder one half of contributions according to the rank of matings
        Rank(1:Data%nMat) = MrgRnk(Chrom((Data%nPotPar1+1):(Data%nPotPar1+Data%nMat)))
        MatPar2(:) = MatPar2(Rank(1:Data%nMat))
      end if

      ! Pair the contributions (=Mating plan)
      k = Data%nMat ! MrgRnk ranks small to large
      if (Spec%GenderMatters .or. Spec%SelfingAllowed) then
        ! When gender matters selfing can not happen (we have two distinct sets of parents,
        ! unless the user adds individuals of one sex in both sets) and when SelfingAllowed
        ! we do not need to care about it - faster code
        do i = 1, Data%nPotPar1
          do j = 1, nVecPar1(i)
            !if (k<2) print*, k, i, j, nVecPar1(i), Data%nMat, sum(nVecPar1(:))
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
            if (MatPar2(k) == Data%IdPotPar1(i)) then
              ! Try to avoid selfing by swapping the MatPar2 and Rank elements
              do l = k, 1, -1
                if (MatPar2(l) /= Data%IdPotPar1(i)) then
                  MatPar2([k, l]) = MatPar2([l, k])
                  Chrom(Data%nPotPar1+Rank([k, l])) = Chrom(Data%nPotPar1+Rank([l, k]))
                  exit
                end if
              end do
              if (l < 1) then ! Above loop ran out without finding a swap
                This%Criterion = This%Criterion + Spec%SelfingWeight
                if (Spec%SelfingWeight < 0.0d0) then
                  This%Penalty = This%Penalty + abs(Spec%SelfingWeight)
                end if
              end if
            end if
            This%MatingPlan(2, k) = MatPar2(k)
            k = k - 1
          end do
        end do
      end if

      ! --- Contribution value ---

      if (Spec%SelCriterionAvailable) then
        !@todo save SelCriterion mean and sd in the data object and then compute this dot_product only once and
        !      compute This%SelCriterion as This%SelCriterion = This%SelIntensity * SelCriterionSD + SelCriterionMean
        This%SelCriterion = dot_product(This%xVec, Data%SelCriterion)
        This%SelIntensity = dot_product(This%xVec, Data%SelCriterionStand)
        if (Spec%PAGE) then
          !@todo as above
          This%SelCriterion = This%SelCriterion + dot_product(This%xVec, Data%SelCriterionPAGE(:)      * This%GenomeEdit(:))
          This%SelIntensity = This%SelIntensity + dot_product(This%xVec, Data%SelCriterionPAGEStand(:) * This%GenomeEdit(:))
        end if
        if (CritType == "opt") then
          This%Criterion = This%Criterion + This%SelIntensity
        end if
      end if

      ! --- Generic individual values ---

      if (Spec%GenericIndValAvailable) then
        do j = 1, Data%nGenericIndVal
          TmpR = dot_product(This%xVec, Data%GenericIndVal(:, j))
          This%GenericIndVal(j) = TmpR
          TmpR = Spec%GenericIndValWeight(j) * This%GenericIndVal(j)
          This%Criterion = This%Criterion + TmpR
          if (Spec%GenericIndValWeight(j) .lt. 0.0) then
            This%Penalty = This%Penalty + abs(TmpR)
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

      if      (CritType == "min" .or. CritType == "ran") then
        This%Criterion = This%Criterion - This%FutureCoancestryRanMate
      else if (CritType == "opt") then
        ! We know the targeted rate of coancestry so we can work with relative values,
        ! which makes the CoancestryWeight generic for ~any scenario.
        TmpR = This%CoancestryRateRanMate / Spec%TargetCoancestryRate
        if (TmpR > 1.0d0) then
          ! Rate of coancestry for the solution is higher than the target
          TmpR = Spec%CoancestryWeight * abs(1.0d0 - TmpR)
        else
          ! Rate of coancestry for the solution is lower than the target
          if (Spec%CoancestryWeightBellow) then
            TmpR = Spec%CoancestryWeight * abs(1.0d0 - abs(TmpR)) ! the second abs is to handle negative coancestry cases
          else
            TmpR = 0.0d0
          end if
        end if
        This%Criterion = This%Criterion + TmpR
        if (Spec%CoancestryWeight < 0.0d0) then
          This%Penalty = This%Penalty + abs(TmpR)
        end if
      else
        write(STDERR, "(a)") "ERROR: Wrong CritType!!!"
        write(STDERR, "(a)") " "
        stop 1
      end if

      ! --- Progeny inbreeding (=inbreeding of a mating) ---

      !print*, "@todo need to check progeny inbreeding in light of genomic coancestry matrix!!!"
      TmpR = 0.0d0
      do j = 1, Data%nMat
        ! Lower triangle to speedup lookup
        TmpMax = maxval(This%MatingPlan(:, j))
        TmpMin = minval(This%MatingPlan(:, j))
        TmpR = TmpR + Data%Coancestry%Value(TmpMax, TmpMin)
      end do
      ! TODO: different number of progeny per mating???
      This%FutureInbreeding = TmpR / Data%nMat
      ! dF = (F_t+1 - F_t) / (1 - F_t)
      This%InbreedingRate = (This%FutureInbreeding - Data%CurrentInbreeding) / (1.0d0 - Data%CurrentInbreeding)
      ! We know the targeted rate of inbreeding so we can work with relative values,
      ! which makes the InbreedingWeight generic for ~any scenario.
      TmpR = This%InbreedingRate / Spec%TargetInbreedingRate
      if (TmpR > 1.0d0) then
        ! Rate of inbreeding for the solution is higher than the target
        TmpR = Spec%InbreedingWeight * abs(1.0d0 - TmpR)
      else
        ! Rate of inbreeding for the solution is lower than the target
        if (Spec%InbreedingWeightBellow) then
          TmpR = Spec%InbreedingWeight * abs(1.0d0 - abs(TmpR)) ! the second abs is to handle negative inbreeding cases
        else
          TmpR = 0.0d0
        end if
      end if
      This%Criterion = This%Criterion + TmpR
      if (Spec%InbreedingWeight < 0.0d0) then
        This%Penalty = This%Penalty + abs(TmpR)
      end if

      ! --- Generic mating values ---

      if (Spec%GenericMatValAvailable) then
        do k = 1, Data%nGenericMatVal
          TmpR = 0.0d0
          if (Spec%GenderMatters) then
            do j = 1, Data%nMat
              TmpR = TmpR + Data%GenericMatVal(Data%IdPotParSeq(This%MatingPlan(1, j)), &
                                               Data%IdPotParSeq(This%MatingPlan(2, j)), k)
            end do
          else
            do j = 1, Data%nMat
              ! Speedup lookup
              TmpMax = maxval(This%MatingPlan(:, j))
              TmpMin = minval(This%MatingPlan(:, j))
              TmpR = TmpR + Data%GenericMatVal(TmpMax, TmpMin, k)
            end do
          end if
          This%GenericMatVal(k) = TmpR / Data%nMat
          TmpR = Spec%GenericMatValWeight(k) * This%GenericMatVal(k)
          This%Criterion = This%Criterion + TmpR
          if (Spec%GenericMatValWeight(k) < 0.0) then
            This%Penalty = This%Penalty + abs(TmpR)
          end if
        end do
      end if

      ! @todo how should we handle costs?
    end subroutine

    !###########################################################################

    !-------------------------------------------------------------------------
    !> @brief  Write head of the AlphaMate log
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !-------------------------------------------------------------------------
    subroutine LogHeadAlphaMateSol(LogUnit, String, StringNum) ! not pure due to IO
      implicit none
      integer(int32), intent(in), optional   :: LogUnit !< Unit to write to (default STDOUT)
      character(len=*), intent(in), optional :: String     !< Additional string that will be written before the head
      integer(int32), optional               :: StringNum  !< How much space is needed for the String
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
        write(LogUnit, FMTLOGUNITHEAD)  COLNAMELOGUNIT(:)
      else
        if (present(String)) then
          write(STDOUT, StringFmt, Advance="No") trim(adjustl(String))
        end if
        write(STDOUT, FMTLOGSTDOUTHEAD) COLNAMELOGSTDOUT(:)
      end if
    end subroutine

    !###########################################################################

    !-------------------------------------------------------------------------
    !> @brief  Write the AlphaMate log
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !-------------------------------------------------------------------------
    subroutine LogAlphaMateSol(This, LogUnit, Iteration, AcceptRate, String, StringNum) ! not pure due to IO
      implicit none
      class(AlphaMateSol), intent(in)        :: This       !< AlphaMateSol holder
      integer(int32), intent(in), optional   :: LogUnit    !< Unit to write to (default STDOUT)
      integer(int32), intent(in)             :: Iteration  !< Generation/Iteration
      real(real64), intent(in)               :: AcceptRate !< Acceptance rate
      character(len=*), intent(in), optional :: String     !< Additional string that will be written before the head
      integer(int32), optional               :: StringNum  !< How much space is needed for the String
      integer(int32) :: Unit
      character(len=100) :: Fmt, StringFmt
      if (present(LogUnit)) then
        Unit = LogUnit
        Fmt = FMTLOGUNIT
      else
        Unit = STDOUT
        Fmt = FMTLOGSTDOUT
      end if
      if (present(String)) then
        if (present(StringNum)) then
          StringFmt = "("//trim(Int2Char(StringNum))//"a)"
        else
          StringFmt = "(a)"
        end if
      end if
      if (Spec%GenericIndValAvailable) then
        if (Spec%GenericMatValAvailable) then
          if (present(String)) then
            write(Unit, StringFmt, Advance="No") trim(adjustl(String))
          end if
          write(Unit, Fmt) Iteration, AcceptRate, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate, This%GenericIndVal, This%GenericMatVal
        else
          if (present(String)) then
            write(Unit, StringFmt, Advance="No") trim(adjustl(String))
          end if
          write(Unit, Fmt) Iteration, AcceptRate, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate, This%GenericIndVal
        end if
      else
        if (Spec%GenericMatValAvailable) then
          if (present(String)) then
            write(Unit, StringFmt, Advance="No") trim(adjustl(String))
          end if
          write(Unit, Fmt) Iteration, AcceptRate, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate,                     This%GenericMatVal
        else
          if (present(String)) then
            write(Unit, StringFmt, Advance="No") trim(adjustl(String))
          end if
          write(Unit, Fmt) Iteration, AcceptRate, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate
        end if
      end if
    end subroutine

    !###########################################################################

    !-------------------------------------------------------------------------
    !> @brief  Write head of the AlphaMate log - for the swarm/population of solutions
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !-------------------------------------------------------------------------
    subroutine LogPopHeadAlphaMateSol(LogPopUnit) ! not pure due to IO
      implicit none
      integer(int32), intent(in), optional :: LogPopUnit  !< log file unit (default STDOUT)
      integer(int32) :: Unit
      if (present(LogPopUnit)) then
        Unit = LogPopUnit
      else
        Unit = STDOUT
      end if
      write(Unit, FMTLOGUNITHEAD) COLNAMELOGPOPUNIT(:)
    end subroutine

    !###########################################################################

    !-------------------------------------------------------------------------
    !> @brief  Write the AlphaMate log - for the swarm/population of solutions
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !-------------------------------------------------------------------------
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
      if (Spec%GenericIndValAvailable) then
        if (Spec%GenericMatValAvailable) then
          write(Unit, FMTLOGPOPUNIT) Iteration, i, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate, This%GenericIndVal, This%GenericMatVal
        else
          write(Unit, FMTLOGPOPUNIT) Iteration, i, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate, This%GenericIndVal
        end if
      else
        if (Spec%GenericMatValAvailable) then
          write(Unit, FMTLOGPOPUNIT) Iteration, i, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate,                     This%GenericMatVal
        else
          write(Unit, FMTLOGPOPUNIT) Iteration, i, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate
        end if
      end if
    end subroutine

    !###########################################################################
end module

!###############################################################################