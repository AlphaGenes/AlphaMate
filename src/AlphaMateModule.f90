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
  use ConstantModule, only : FILELENGTH, SPECOPTIONLENGTH, IDLENGTH, IDINTLENGTH

! TODO: use IDLENGTH and IDINTLENGTH

  use OrderPackModule, only : MrgRnk
  use AlphaHouseMod, only : CountLines, Char2Int, Char2Double, Int2Char, Real2Char, &
                            RandomOrder, SetSeed, ToLower, FindLoc, &
                            ParseToFirstWhitespace, SplitLineIntoTwoParts
  use AlphaStatMod, only : DescStat, DescStatReal64, DescStatMatrix, DescStatMatrixReal64, &
                           DescStatSymMatrix, DescStatLowTriMatrix
  use AlphaEvolveModule, only : AlphaEvolveSol, DifferentialEvolution, RandomSearch
  use AlphaRelateModule

  implicit none

  private
  ! Types
  public :: AlphaMateSpec, AlphaMateData, AlphaMateSol
  ! Functions
  public :: AlphaMateTitle, ReadAlphaMateSpec, ReadAlphaMateData, SetupColNamesAndFormats, AlphaMateSearch

  !> @brief AlphaMate specifications
  type AlphaMateSpec
    ! Files
    character(len=FILELENGTH) :: SpecFile, RelMtxFile, SelCriterionFile, GenderFile, SeedFile
    character(len=FILELENGTH) :: GenericIndCritFile, GenericMatCritFile
    character(len=FILELENGTH) :: OutputBasename
    logical :: RelMtxFileGiven, SelCriterionFileGiven, GenderFileGiven, SeedFileGiven, GenericIndCritFileGiven, GenericMatCritFileGiven

    ! Biological specifications
    logical :: NrmInsteadOfCoancestry
    real(real64) :: TargetCoancestryRate
    real(real64) :: TargetInbreedingRate
    real(real64) :: TargetCoancestryRateWeight, TargetInbreedingRateWeight, SelfingWeight
    integer(int32) :: nMat, nPar, nPar1, nPar2
    logical :: EqualizePar, EqualizePar1, EqualizePar2, LimitPar, LimitPar1, LimitPar2
    real(real64) :: LimitParMin, LimitPar1Min, LimitPar2Min, LimitParMax, LimitPar1Max, LimitPar2Max, LimitParMinWeight, LimitPar1MinWeight, LimitPar2MinWeight
    integer(int32) :: nGenericIndCrit, nGenericMatCrit
    logical :: GenericIndCritGiven, GenericMatCritGiven
    real(real64), allocatable :: GenericIndCritWeight(:), GenericMatCritWeight(:)
    logical :: SelfingAllowed, TargetCoancestryRateWeightBelow, TargetInbreedingRateWeightBelow
    logical :: PAGEPar, PAGEPar1, PAGEPar2
    integer(int32) :: PAGEParMax, PAGEPar1Max, PAGEPar2Max
    real(real64) :: PAGEParCost, PAGEPar1Cost, PAGEPar2Cost

    ! Search specifications
    logical :: SeedGiven
    logical :: ModeMin, ModeRan, ModeOpt
    logical :: EvaluateFrontier
    integer(int32) :: Seed, nFrontierPoints
    real(real64), allocatable :: TargetCoancestryRateFrontier(:)

    ! Algorithm specifications
    integer(int32) :: EvolAlgNSol, EvolAlgNIter, EvolAlgNIterBurnIn, EvolAlgNIterStop, EvolAlgNIterPrint, RanAlgStricter
    real(real64) :: EvolAlgStopTol, EvolAlgParamCrBurnIn, EvolAlgParamCr, EvolAlgParamFBase, EvolAlgParamFHigh1, EvolAlgParamFHigh2
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
    real(real64), allocatable :: GenericIndCrit(:, :), GenericMatCrit(:, :, :)
    ! Data summaries
    type(DescStatReal64) :: InbreedingStat, SelCriterionStat, SelCriterionPAGEStat
    type(DescStatReal64), allocatable :: GenericIndCritStat(:)
    type(DescStatMatrixReal64) :: CoancestryStat, CoancestryStatGenderDiff, CoancestryStatGender1, CoancestryStatGender2
    type(DescStatMatrixReal64), allocatable :: GenericMatCritStat(:)
    ! Derived data
    integer(int32) :: nInd, nPotMat, nPotPar1, nPotPar2, nMal, nFem
    integer(int32), allocatable :: IdPotPar1(:), IdPotPar2(:), IdPotParSeq(:)
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
    real(real64), allocatable   :: GenericIndCrit(:)
    real(real64), allocatable   :: GenericMatCrit(:)
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

    !---------------------------------------------------------------------------
    !> @brief  Print AlphaMate title
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
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
    !> @brief  Initialize AlphaMate specifications
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !---------------------------------------------------------------------------
    pure subroutine InitAlphaMateSpec(This)
      implicit none
      class(AlphaMateSpec), intent(out) :: This !< AlphaMateSpec holder

      ! Inputs

      This%SpecFile = ""
      This%RelMtxFile = ""
      This%SelCriterionFile = ""
      This%GenderFile = ""
      This%GenericIndCritFile = ""
      This%nGenericIndCrit = 0
      This%GenericMatCritFile = ""
      This%nGenericMatCrit = 0
      This%SeedFile = ""

      This%RelMtxFileGiven = .false.
      This%SelCriterionFileGiven = .false.
      This%GenderFileGiven = .false.
      This%GenericIndCritGiven = .false.
      This%GenericMatCritGiven = .false.
      This%SeedFileGiven = .false.
      This%SeedGiven = .false.
      This%NrmInsteadOfCoancestry = .false.

      ! Biological specifications

      This%TargetCoancestryRate = 0.01d0
      This%TargetInbreedingRate = 0.01d0
      ! This%TargetCoancestryRateFrontier(:) ! allocatable so skip here
      This%TargetCoancestryRateWeight = 0.5d0
      This%TargetCoancestryRateWeightBelow = .false.
      This%TargetInbreedingRateWeight =  0.5d0
      This%TargetInbreedingRateWeightBelow = .false.
      This%SelfingWeight = 0.0d0
      ! This%GenericIndCritWeight(:) ! allocatable so skip here
      ! This%GenericMatCritWeight(:) ! allocatable so skip here

      This%SelfingAllowed = .false.

      This%EqualizePar  = .false.
      This%EqualizePar1 = .false.
      This%EqualizePar2 = .false.

      This%LimitPar  = .false.
      This%LimitPar1 = .false.
      This%LimitPar2 = .false.
      This%LimitParMin  = 1.0d0
      This%LimitPar1Min = 1.0d0
      This%LimitPar2Min = 1.0d0
      This%LimitParMax  = huge(Spec%LimitParMax)  - 1.0d0
      This%LimitPar1Max = huge(Spec%LimitPar1Max) - 1.0d0
      This%LimitPar2Max = huge(Spec%LimitPar2Max) - 1.0d0
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

      ! Search mode specifications

      This%ModeMin = .false.
      This%ModeRan = .false.
      This%ModeOpt = .false.
      This%EvaluateFrontier = .false.
      This%nFrontierPoints = 0

      ! Search algorithm specifications

      This%EvolAlgNSol = 100
      This%EvolAlgNIter = 10000
      This%EvolAlgNIterBurnIn = 500
      This%EvolAlgNIterStop = 1000
      This%EvolAlgNIterPrint = 100
      This%EvolAlgStopTol = 0.001d0
      This%EvolAlgParamCrBurnIn = 0.4d0
      This%EvolAlgParamCr = 0.2d0
      This%EvolAlgParamFBase = 0.1d0
      This%EvolAlgParamFHigh1 = 1.0d0
      This%EvolAlgParamFHigh2 = 4.0d0
      This%EvolAlgLogPop = .false.
      This%RanAlgStricter = 10
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Read AlphaMate specifications from a file
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !---------------------------------------------------------------------------
    subroutine ReadAlphaMateSpec(This, SpecFile, LogStdout) ! not pure due to IO
      implicit none
      class(AlphaMateSpec), intent(out) :: This      !< AlphaMateSpec holder
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
      call This%Init

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

            case ("coancestrymatrixfile")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .ne. "none") then
                  write(This%RelMtxFile, *) trim(adjustl(Second(1)))
                  This%RelMtxFile = adjustl(This%RelMtxFile)
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Coancestry matrix file: "//trim(This%RelMtxFile)
                  end if
                  This%RelMtxFileGiven = .true.
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
                  This%RelMtxFileGiven = .true.
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
                  This%SelCriterionFileGiven = .true.
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
                  This%GenderFileGiven = .true.
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a file for GenderFile, i.e., GenderFile, Gender.txt"
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
              if (This%GenericIndCritFileGiven) then
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
              if (This%GenericIndCritFileGiven) then
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
              if (This%GenericMatCritFileGiven) then
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
              if (This%GenericMatCritFileGiven) then
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
                  write(STDOUT, "(a)") " Number of targeted parents: "//trim(Int2Char(This%nPar))
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
                  write(STDOUT, "(a)") " Number of targeted male parents: "//trim(Int2Char(This%nPar1))
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
                  write(STDOUT, "(a)") " Number of targeted female parents: "//trim(Int2Char(This%nPar2))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfFemaleParents, i.e., NumberOfFemaleParents, 10"
                write(STDERR, "(a)") ""
                stop 1
              end if

            case ("targetcoancestryrate")
              if (allocated(Second)) then
                This%TargetCoancestryRate = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Targeted rate of coancestry: "//trim(Real2Char(This%TargetCoancestryRate, fmt=FMTREAL2CHAR))
                end if
                if (This%TargetCoancestryRate .eq. 0.0) then
                  write(STDERR, "(a)") "ERROR: Can not work with the targeted rate of coancestry exactly equal to zero - it is numerically unstable!"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
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
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetCoancestryRateWeight, i.e., TargetCoancestryRateWeight, 0.01"
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
                if (This%TargetInbreedingRate .eq. 0.0) then
                  write(STDERR, "(a)") "ERROR: Can not work with the targeted rate of inbreeding exactly equal to zero - it is numerically unstable!"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
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
                write(STDERR, "(a)") " ERROR: Must specify a value for TargetInbreedingRateWeight, i.e., TargetInbreedingRateWeight, 0.01"
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

            case ("evaluatefrontier")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  This%EvaluateFrontier = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Evaluate selection/coancestry frontier"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify Yes or No for EvaluateFrontier, i.e., EvaluateFrontier, Yes"
                write(STDERR, "(a)") ""
                stop 1
              end if

            case ("evaluatefrontiernumberofpoints")
              if (This%EvaluateFrontier) then
                if (allocated(Second)) then
                  This%nFrontierPoints = Char2Int(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Evaluate selection/coancestry frontier - number of points: "//trim(Int2Char(This%nFrontierPoints))
                  end if
                  allocate(This%TargetCoancestryRateFrontier(This%nFrontierPoints))
                  nFrontierPoint = 0
                else
                  write(STDERR, "(a)") " ERROR: Must specify a number for EvaluateFrontierNumberOfPoints, i.e., EvaluateFrontierNumberOfPoints, 3"
                  write(STDERR, "(a)") ""
                  stop 1
                end if
              end if

            case ("evaluatefrontiertargetcoancestryrate")
              if (This%EvaluateFrontier) then
                if (allocated(Second)) then
                  nFrontierPoint = nFrontierPoint + 1
                  This%TargetCoancestryRateFrontier(nFrontierPoint) = Char2Double(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Evaluate selection/coancestry frontier - coancestry rate ("//trim(Int2Char(nFrontierPoint))//"): "//trim(Real2Char(This%TargetCoancestryRateFrontier(nFrontierPoint), fmt=FMTREAL2CHAR))
                  end if
                  if (This%TargetCoancestryRateFrontier(nFrontierPoint) .eq. 0.0) then
                    write(STDERR, "(a)") "ERROR: Can not work with the targeted rate of coancestry exactly equal to zero - it is numerically unstable!"
                    write(STDERR, "(a)") " "
                    stop 1
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify a value for EvaluateFrontierTargetCoancestryRate, i.e., EvaluateFrontierTargetCoancestryRate, 0.001"
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

            case ("evolalgnumberofiterationsburnin")
              if (allocated(Second)) then
                This%EvolAlgNIterBurnIn = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Evolutionary algorithm - number of warming iterations (burn-in): "//trim(Int2Char(This%EvolAlgNIterBurnIn))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for EvolAlgNumberOfIterationsBurnin, i.e., EvolAlgNumberOfIterationsBurnin, 1000"
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

            case ("evolalgparametercrburnin")
              if (allocated(Second)) then
                This%EvolAlgParamCrBurnIn = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Evolutionary algorithm - cross-over parameter for warmup (burn-in): "//trim(Real2Char(This%EvolAlgParamCrBurnIn, fmt=FMTREAL2CHAR))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for EvolAlgParameterCrBurnin, i.e., EvolAlgParameterCrBurnin, 0.4"
                write(STDERR, "(a)") ""
                stop 1
              end if

            case ("evolalgparametercr")
              if (allocated(Second)) then
                This%EvolAlgParamCr = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Evolutionary algorithm - cross-over parameter: "//trim(Real2Char(This%EvolAlgParamCr, fmt=FMTREAL2CHAR))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for EvolAlgParameterCr, i.e., EvolAlgParameterCr, 0.2"
                write(STDERR, "(a)") ""
                stop 1
              end if

            case ("evolalgparameterfbase")
              if (allocated(Second)) then
                This%EvolAlgParamFBase = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Evolutionary algorithm - parameter F (base value): "//trim(Real2Char(This%EvolAlgParamFBase, fmt=FMTREAL2CHAR))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for EvolAlgParameterFBase, i.e., EvolAlgParameterFBase, 0.1"
                write(STDERR, "(a)") ""
                stop 1
              end if

            case ("evolalgparameterfhigh1")
              if (allocated(Second)) then
                This%EvolAlgParamFHigh1 = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Evolutionary algorithm - parameter F (high value 1): "//trim(Real2Char(This%EvolAlgParamFHigh1, fmt=FMTREAL2CHAR))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for EvolAlgParameterFHigh1, i.e., EvolAlgParameterFHigh1, 1.0"
                write(STDERR, "(a)") ""
                stop 1
              end if

            case ("evolalgparameterfhigh2")
              if (allocated(Second)) then
                This%EvolAlgParamFHigh2 = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Evolutionary algorithm - parameter F (high value 2): "//trim(Real2Char(This%EvolAlgParamFHigh2, fmt=FMTREAL2CHAR))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for EvolAlgParameterFHigh2, i.e., EvolAlgParameterFHigh2, 4.0"
                write(STDERR, "(a)") ""
                stop 1
              end if

            case ("randomsearchstricter")
              if (allocated(Second)) then
                This%RanAlgStricter = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Random search - perform k times more iterations that with the evolutionary algorithm: k="//trim(Int2Char(This%RanAlgStricter))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for RandomSearchStricter, i.e., RandomSearchStricter, 10"
                write(STDERR, "(a)") ""
                stop 1
              end if

            case ("stop")
              if (LogStdoutInternal) then
                write(STDOUT, "(3a)") " NOTE: Encountered Stop specification - the rest of specifications will be ignored"
                write(STDOUT, "(a)") " "
              end if
              exit

            case default
              if (LogStdoutInternal) then
                write(STDOUT, "(3a)") " NOTE: Specification '", trim(Line), "' was ignored"
                write(STDOUT, "(a)") " "
              end if
          end select
        end if
      end do ReadSpec
      close(SpecUnit)

      if (.not. This%RelMtxFileGiven) then
        write(STDERR, "(a)") " ERROR: Must specify CoancestryMatrixFile or NrmMatrixFile!"
        write(STDERR, "(a)") ""
        stop 1
      end if

      if (This%LimitPar .and. This%EqualizePar) then
        write(STDOUT, "(a)") " NOTE: The specification Equalize*Contributions has priority over Limit*Contributions."
        write(STDOUT, "(a)") " "
        ! ... therefore reset all limit specifications to default values
        This%LimitPar  = .false.
        This%LimitPar1 = .false.
        This%LimitPar2 = .false.
        This%LimitParMin  = 1.0d0
        This%LimitPar1Min = 1.0d0
        This%LimitPar2Min = 1.0d0
        This%LimitParMax  = huge(Spec%LimitParMax)  - 1.0d0
        This%LimitPar1Max = huge(Spec%LimitPar1Max) - 1.0d0
        This%LimitPar2Max = huge(Spec%LimitPar2Max) - 1.0d0
        This%LimitParMinWeight  = 0.0d0
        This%LimitPar1MinWeight = 0.0d0
        This%LimitPar2MinWeight = 0.0d0
      end if

      if (.not. This%GenderFileGiven) then
        This%nPar1 = This%nPar
        This%EqualizePar1 = This%EqualizePar
        This%LimitPar1 = This%LimitPar
        This%PAGEPar1 = This%PAGEPar
      end if

      if (This%GenderFileGiven) then
        This%nPar = This%nPar1 + This%nPar2
      end if

      if (This%GenderFileGiven .and. This%SelfingAllowed) then
        write(STDERR, "(a)") " ERROR: When gender matters, AlphaMate can not perform selfing! See the manual for a solution."
        ! @todo: what is the solution? Provide the same individual both as male and a female?
        write(STDERR, "(a)") " "
        stop 1
      end if

      if ((.not. This%SelCriterionFileGiven) .and. This%PAGEPar) then
        write(STDERR, "(a)") " ERROR: Can not use PAGE when selection criterion file is not given!"
        ! @todo: what about using the GenericIndCrit values?
        write(STDERR, "(a)") " "
        stop 1
      end if

      if (This%SeedFileGiven .and. This%SeedGiven) then
        write(STDOUT, "(a)") " NOTE: The specification Seed has priority over SeedFile."
        write(STDOUT, "(a)") " "
        This%SeedFile = ""
        This%SeedFileGiven = .false.
      end if
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Read AlphaMate data from a file and summarize it for further use
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !---------------------------------------------------------------------------
    subroutine ReadAlphaMateData(This, Spec, LogStdout) ! not pure due to IO
      implicit none
      class(AlphaMateData), intent(out)  :: This      !< AlphaMateData holder
      type(AlphaMateSpec), intent(inout) :: Spec      !< AlphaMateSpec holder
      logical, optional                  :: LogStdout !< Log process on stdout (default .false.)

      integer(int32) :: i, j, nIndTmp, GenderTmp
      integer(int32) :: SelCriterionUnit, GenderUnit, SeedUnit

      real(real64) :: SelCriterionTmp, SelCriterionTmp2

      character(len=IDLENGTH) :: IdCTmp

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

      if (LogStdoutInternal) then
        write(STDOUT, "(a)") " Number of individuals in the coancestry matrix file: "//trim(Int2Char(This%nInd))
      end if

      ! --- Selection criterion ---

      allocate(This%SelCriterion(This%nInd))
      allocate(This%SelCriterionStand(This%nInd))
      if (Spec%PAGEPar) then
        allocate(This%SelCriterionPAGE(This%nInd))
        allocate(This%SelCriterionPAGEStand(This%nInd))
      end if

      if (.not. Spec%SelCriterionFileGiven) then
        This%SelCriterion(:) = 0.0d0
        This%SelCriterionStand(:) = 0.0d0
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
        do i = 1, This%nInd
          if (Spec%PAGEPar) then
            read(SelCriterionUnit, *) IdCTmp, SelCriterionTmp, SelCriterionTmp2
          else
            read(SelCriterionUnit, *) IdCTmp, SelCriterionTmp
          end if
          j = FindLoc(IdCTmp, This%Coancestry%OriginalId(1:))
          if (j .eq. 0) then
            write(STDERR, "(a)") " ERROR: Individual "//trim(IdCTmp)//" from the selection criterion file not present in the coancestry matrix file!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          This%SelCriterion(j) = SelCriterionTmp
          if (Spec%PAGEPar) then
            This%SelCriterionPAGE(j) = SelCriterionTmp2
          end if
        end do
        close(SelCriterionUnit)
      end if

        ! --- Gender ---

        allocate(This%Gender(This%nInd))
        This%Gender(:) = 0
        if (Spec%GenderFileGiven) then
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
          do i = 1, This%nInd
            read(GenderUnit, *) IdCTmp, GenderTmp
            if      (GenderTmp .eq. 1) then
              This%nMal = This%nMal + 1
            else if (GenderTmp .eq. 2) then
              This%nFem = This%nFem + 1
            else
              write(STDERR, "(a)") " ERROR: Gender code must be either 1 for male individuals or 2 for female individuals!"
              write(STDERR, "(a)") " ERROR: "//trim(Int2Char(i))//" "//trim(IdCTmp)//" "//trim(Int2Char(GenderTmp))
              write(STDERR, "(a)") " "
              stop 1
            end if
            j = FindLoc(IdCTmp, This%Coancestry%OriginalId(1:))
            if (j == 0) then
              write(STDERR, "(a)") " ERROR: Individual "//trim(IdCTmp)//" from the gender file not present in the coancestry matrix file!"
              write(STDERR, "(a)") " "
              stop 1
            end if
            This%Gender(j) = GenderTmp
          end do
          close(GenderUnit)

          write(STDOUT, "(a)") " Number of   males: "//trim(Int2Char(This%nMal))
          write(STDOUT, "(a)") " Number of females: "//trim(Int2Char(This%nFem))
          write(STDOUT, "(a)") " "

          if (Spec%nPar1 > This%nMal) then
            write(STDERR, "(a)") " ERROR: The number of male parents can not be larger than the number of males"
            write(STDERR, "(a)") " ERROR: Number of male parents: "//trim(Int2Char(Spec%nPar1))
            write(STDERR, "(a)") " ERROR: Number of        males: "//trim(Int2Char(This%nMal))
            write(STDERR, "(a)") " "
            stop 1
          end if
          if (Spec%nPar2 > This%nFem) then
            write(STDERR, "(a)") " ERROR: The number of female parents can not be larger than the number of females"
            write(STDERR, "(a)") " ERROR: Number of female parents: "//trim(Int2Char(Spec%nPar2))
            write(STDERR, "(a)") " ERROR: Number of        females: "//trim(Int2Char(This%nFem))
            write(STDERR, "(a)") " "
            stop 1
          end if
        end if

! TODO

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

    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Setup colnames and formats for output
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !---------------------------------------------------------------------------
    subroutine SetupColNamesAndFormats ! not pure due to setting module-wise variables
      implicit none
      integer(int32) :: nCol, nColTmp, i

      ! --- Optimisation log ---

      nCol = 10
      if (Spec%GenericIndCritGiven) then
        nCol = nCol + Spec%nGenericIndCrit
      end if
      if (Spec%GenericMatCritGiven) then
        nCol = nCol + Spec%nGenericMatCrit
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
      if (Spec%GenericIndCritGiven) then
        do i = 1, Spec%nGenericIndCrit
          nColTmp = nColTmp + 1
          COLNAMELOGUNIT(nColTmp) = "GenIndCrit"//trim(Int2Char(i))
          COLNAMELOGUNIT(nColTmp) = adjustr(COLNAMELOGUNIT(nColTmp))
        end do
      end if
      if (Spec%GenericMatCritGiven) then
        do i = 1, Spec%nGenericMatCrit
          nColTmp = nColTmp + 1
          COLNAMELOGUNIT(nColTmp) = "GenMatCrit"//trim(Int2Char(i))
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

    !---------------------------------------------------------------------------
    !> @brief  Call various optimisations for AlphaMate
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !---------------------------------------------------------------------------
    subroutine AlphaMateSearch ! not pure due to IO

      implicit none

      integer(int32) :: nParam, k, FrontierUnit

      real(real64) :: HoldTargetCoancestryRanMate, HoldTargetCoancestryRate
      real(real64), allocatable :: InitEqual(:, :)

      character(len=1000) :: LogFile, LogPopFile, ContribFile, MatingFile
      character(len=100) :: DumC

      type(AlphaMateSol) :: SolMin, SolRan, SolOpt, Sol

      ! --- Number of parameters to optimise ---

      if (Spec%GenderFileGiven) then
        nParam = Data%nPotPar1 + Data%nPotPar2 + Spec%nMat
      else
        nParam = Data%nPotPar1 + Spec%nMat
      end if

      if (Spec%PAGEPar) then
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
          CritType="min", CRBurnIn=Spec%EvolAlgParamCrBurnIn, CRLate=Spec%EvolAlgParamCr, FBase=Spec%EvolAlgParamFBase, FHigh1=Spec%EvolAlgParamFHigh1, FHigh2=Spec%EvolAlgParamFHigh2, &
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
          CritType="opt", CRBurnIn=Spec%EvolAlgParamCrBurnIn, CRLate=Spec%EvolAlgParamCr, FBase=Spec%EvolAlgParamFBase, FHigh1=Spec%EvolAlgParamFHigh1, FHigh2=Spec%EvolAlgParamFHigh2, &
          BestSol=SolOpt)

        call SolOpt%Write(Data, Spec, ContribFile, MatingFile)
      end if

      ! --- Evaluate the full frontier ---

      if (Spec%EvaluateFrontier) then
        write(STDOUT, "(a)") "--- Evaluate the frontier ---"
        write(STDOUT, "(a)") " "

        Spec%TargetCoancestryRateWeightBelow = .true. ! we want to target certain rates of coancestry

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
        do k = 1, Spec%nFrontierPoints
          Spec%TargetCoancestryRate = Spec%TargetCoancestryRateFrontier(k)
          ! F_t = DeltaF + (1 - DeltaF) * F_t-1
          Data%TargetCoancestryRanMate = Spec%TargetCoancestryRate + (1.0d0 - Spec%TargetCoancestryRate) * Data%CurrentCoancestryRanMate
          write(STDOUT, "(a)") "Step "//trim(Int2Char(k))//" out of "//trim(Int2Char(Spec%nFrontierPoints))//&
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
            CritType="opt", CRBurnIn=Spec%EvolAlgParamCrBurnIn, CRLate=Spec%EvolAlgParamCr, FBase=Spec%EvolAlgParamFBase, FHigh1=Spec%EvolAlgParamFHigh1, FHigh2=Spec%EvolAlgParamFHigh2, &
            BestSol=Sol)

          DumC = "Frontier"//trim(Int2Char(k))
          call Sol%Log(FrontierUnit, Iteration=-1, AcceptRate=-1.0d0, String=DumC, StringNum=10)
          call Sol%Write(Data, Spec, ContribFile, MatingFile)

          if ((Spec%TargetCoancestryRate - Sol%CoancestryRateRanMate) .gt. 0.01d0) then
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

    !---------------------------------------------------------------------------
    !> @brief  Write AlphaMate solution to files or standard output
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !---------------------------------------------------------------------------
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
      if (.not.Spec%PAGEPar) then
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
      do i = 1, Spec%nMat
        write(MatingUnit, FMTMATING) i, Data%Coancestry%OriginalId(This%MatingPlan(1, i)), Data%Coancestry%OriginalId(This%MatingPlan(2, i))
      end do
      if (present(MatingFile)) then
        close(MatingUnit)
      end if
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Initialize AlphaMate solution
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !---------------------------------------------------------------------------
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
      if (Spec%GenericIndCritGiven) then
        allocate(This%GenericIndCrit(Spec%nGenericIndCrit))
        This%GenericIndCrit(:) = 0.0d0
      else
        allocate(This%GenericIndCrit(0))
      end if
      if (Spec%GenericMatCritGiven) then
        allocate(This%GenericMatCrit(Spec%nGenericMatCrit))
        This%GenericMatCrit(:) = 0.0d0
      else
        allocate(This%GenericMatCrit(0))
      end if
      This%Cost = 0.0d0
      allocate(This%nVec(Data%nInd))
      This%nVec(:) = 0
      allocate(This%xVec(Data%nInd))
      This%xVec(:) = 0.0d0
      allocate(This%MatingPlan(2, Spec%nMat))
      This%MatingPlan(:, :) = 0
      if (Spec%PAGEPar) then
        allocate(This%GenomeEdit(Data%nInd))
        This%GenomeEdit(:) = 0.0d0
      else
        allocate(This%GenomeEdit(0))
      end if
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
          if (allocated(In%GenericIndCrit)) then
            allocate(Out%GenericIndCrit(size(In%GenericIndCrit)))
            Out%GenericIndCrit = In%GenericIndCrit
          end if
          if (allocated(In%GenericMatCrit)) then
            allocate(Out%GenericMatCrit(size(In%GenericMatCrit)))
            Out%GenericMatCrit = In%GenericMatCrit
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

    !---------------------------------------------------------------------------
    !> @brief  Update mean of AlphaMate solution (when performing random search)
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !---------------------------------------------------------------------------
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
          if (allocated(This%GenericIndCrit)) then
            This%GenericIndCrit          = This%GenericIndCrit            * kR + Add%GenericIndCrit            / n
          end if
          if (allocated(This%GenericMatCrit)) then
            This%GenericMatCrit          = This%GenericMatCrit            * kR + Add%GenericMatCrit            / n
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

    !---------------------------------------------------------------------------
    !> @brief  AlphaMate evaluate function plus much MORE (this is the core!!!!)
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !---------------------------------------------------------------------------
! TODO: how will I pass Data and Spec here??????
    subroutine FixSolEtcMateAndEvaluate(This, Chrom, CritType)
      implicit none
      ! Arguments
      class(AlphaMateSol), intent(inout)     :: This      ! Solution
      real(real64), intent(inout), optional  :: Chrom(:)  ! Internal representation of the solution
      character(len=*), intent(in), optional :: CritType  ! Type of criterion (Min, Ran, Opt)

      ! Other
      integer(int32) :: i, j, k, l, g, nCumMat, Rank(Data%nInd), ChromInt(Data%nInd), MatPar2(Spec%nMat)
      integer(int32) :: nVecPar1(Data%nPotPar1), nVecPar2(Data%nPotPar2), TmpMin, TmpMax, TmpI

      real(real64) :: TmpVec(Data%nInd, 1), TmpR, RanNum

      ! Initialize the solution
      call This%Initialise()

      ! The solution (based on the mate selection driver) has:
      ! - Data%nInd individual contributions
      !   - Data%nPotPar1 individual contributions for "parent1" (males   when GenderFileGiven, all ind when .not. GenderFileGiven)
      !   - Data%nPotPar2 individual contributions for "parent2" (females when GenderFileGiven, meaningful only when GenderFileGiven)
      ! - Spec%nMat     rankings of parent1 1:Spec%nMat matings to pair with 1:Data%nPotPar2 "parent2" (see below)
      ! - Data%nInd edit indicators
      !   - Data%nPotPar1 edit indicators for "parent1" (males   when GenderFileGiven, all ind when .not. GenderFileGiven)
      !   - Data%nPotPar2 edit indicators for "parent2" (females when GenderFileGiven, present only when GenderFileGiven)

      ! Say we have Chrom=(| 0, 2, 0, 1 | ... | 2.5, 1.5, 1.0 | 0, 1, 0, 0 | ...) then we:
      ! - mate male 2 with the first  available female (rank 2.5)
      ! - mate male 2 with the second available female (rank 1.5)
      ! - mate male 4 with the third  available female (rank 1.0)
      ! - edit male 2

      ! @todo consider spliting the Chrom() vector internally into a type with
      !   separate vectors to simplify the code, e.g.,
      !   - Chrom2%ContPar1
      !   - Chrom2%ContPar2
      !   - Chrom2%MateRank
      !   - Chrom2%EditPar1
      !   - Chrom2%EditPar2
      !   and then at the end combine it back en extit - since we modify/fix some
      !   elements of a solution, we need to combine before exit!

      ! --- Parse the mate selection driver (=Is the solution valid?) ---

      ! The approach below assures that we have Spec%nMat contributions for each of
      ! the two parent sets. It does this by ranking internal solution values and
      ! traverses from top to the defined number of parents checking when the sum of
      ! interegrised values gives Spec%nMat. If values below 0.5 are found, they are
      ! changed to 1 contribution. If this still does not give Spec%nMat, then we start

      ! @todo least contributing? Those that already contribute or thos that do not contribute at all?

      ! adding on contribution to each parent (starting from least contributing
      ! parents to avoid local minima) until we reach Spec%nMat. How to treat values
      ! for the individuals that do not contribute is unlcear. None of the tested
      ! methods seemed to be very different. Intuitively, using properly ordered
      ! negative values should inform optim. alg. which individuals should less
      ! likely contribute, but this did not seem to be the case - better final
      ! solution was found when this strategy was not implemented - either zeroing
      ! values for those individuals (was the fastest) or giving random value (was
      ! marginally better, but slower). Potential advantage of not preserving the
      ! order is that this gives more randomness and more solutions being explored.

      ! "Parent1"
      if (Spec%GenderFileGiven) then
        g = 1
      else
        g = 2
      end if
      ! ... find ranks to find the top values
      if (.not.(Spec%EqualizePar1 .and. (Spec%nPar1 .eq. Data%nPotPar1))) then
        Rank(1:Data%nPotPar1) = MrgRnk(Chrom(1:Data%nPotPar1))
        Rank(1:Data%nPotPar1) = Rank(Data%nPotPar1:1:-1) ! MrgRnk ranks small to large
        !@todo decreasing option in MrgRnk?
      end if
      ! ... handle cases with equalized contributions
      if (Spec%EqualizePar1) then
        if (Spec%nPar1 .eq. Data%nPotPar1) then
          ! ... set integers to all the values (no need for sorting here)
          Chrom(1:Data%nPotPar1) = dble(Spec%nMat * g) / Spec%nPar1
        else
          ! ... set integers to the top values
          Chrom(Rank(1:Spec%nPar1)) = dble(Spec%nMat * g) / Spec%nPar1
          !@todo anything better to preserve the order of non contributing individuals? See below!
          Chrom(Rank((Spec%nPar1+1):Data%nPotPar1)) = 0.0d0
          ! Chrom(Rank((Spec%nPar1+1):Data%nPotPar1)) = -1.0d0
        end if
      else
        ! ... handle cases with unequal contributions
        ! ... work for the defined number or parents
        nCumMat = 0
        do i = 1, Spec%nPar1
          j = Rank(i)
          ! ... these would have "zero" contributions and if we hit them then they need at least minimum contrib
          if (Chrom(j) .lt. Spec%LimitPar1Min) then
            Chrom(j) = Spec%LimitPar1Min ! set/fix to minimum usage @todo could consider penalising solution instead?
          end if
          ! ... but not above max allowed
          if (Chrom(j) .gt. Spec%LimitPar1Max) then
            Chrom(j) = Spec%LimitPar1Max ! set/fix to maximum usage @todo could consider penalising solution instead?
          end if
          ! ... accumulate and check if we reached Spec%nMat
          nCumMat = nCumMat + nint(Chrom(j)) ! internally real, externally integer
          if (nCumMat .ge. Spec%nMat * g) then
            ! ... there should be exactly Spec%nMat contributions
            if (nCumMat .gt. Spec%nMat * g) then
              Chrom(j) = Chrom(j) - dble(nCumMat - Spec%nMat * g)
              ! ... did we go below the minimum usage limit?
              if (nint(Chrom(j)) .lt. Spec%LimitPar1Min) then
                TmpR = Spec%LimitPar1MinWeight * (Spec%LimitPar1Min - nint(Chrom(j)))
                This%Criterion = This%Criterion + TmpR
                if (Spec%LimitPar1MinWeight .lt. 0.0d0) then
                  This%Penalty = This%Penalty + abs(TmpR)
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
        if (i .le. Data%nPotPar1) then ! "=" to capture incremented i+1 on the do loop exit
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
          ! if (TmpR .gt. 0.0d0) then
          !     Chrom(Rank(i:Data%nPotPar1)) = Chrom(Rank(i:Data%nPotPar1)) - abs(TmpR)
          ! end if
          ! ... negative (random so no order)
          ! do j = i, Data%nPotPar1 ! @todo really need this loop?
          !   call random_number(RanNum)
          !   Chrom(Rank(j)) = -1.0d0 * RanNum
          ! end do
        end if
        ! ... Spec%nMat still not reached?
        do while (nCumMat .lt. Spec%nMat * g)
          ! ... add more contributions
          do i = Spec%nPar1, 1, -1 ! start with the lowest ranked individuals selected as parents (to avoid local optima)
            j = Rank(i)
            Chrom(j) = Chrom(j) + 1.0d0
            ! ... accumulate and check if we reached Spec%nMat
            nCumMat = nCumMat + 1
            if (nCumMat .ge. Spec%nMat * g) then
              ! To cater for real vs. integer issues
              TmpI = sum(nint(Chrom(Rank(1:Spec%nPar1))))
              if (TmpI .ne. Spec%nMat * g) then
                if (TmpI .gt. Spec%nMat * g) then
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
      if (Spec%GenderFileGiven) then
        ! ... find ranks to find the top values
        if (.not.(Spec%EqualizePar2 .and. (Spec%nPar2 .eq. Data%nPotPar2))) then
          Rank(1:Data%nPotPar2) = MrgRnk(Chrom((Data%nPotPar1+1):(Data%nPotPar1+Data%nPotPar2)))
          Rank(1:Data%nPotPar2) = Rank(Data%nPotPar2:1:-1) ! MrgRnk ranks small to large
        end if
        ! ... handle cases with equalized contributions
        if (Spec%EqualizePar2) then
          if (Spec%nPar2 .eq. Data%nPotPar2) then
            ! ... set integers to all the values (no need for sorting here)
            Chrom((Data%nPotPar1+1):(Data%nPotPar1+Data%nPotPar2)) = dble(Spec%nMat) / Spec%nPar2
          else
            ! ... set integers to the top values
            Chrom(Data%nPotPar1+Rank(1:Spec%nPar2)) = dble(Spec%nMat) / Spec%nPar2
            ! @todo anything better to preserve the order of non contributing individuals? See below!
            Chrom(Data%nPotPar1+Rank((Spec%nPar2+1):Data%nPotPar2)) = 0.0d0
            ! Chrom(Data%nPotPar1+Rank((Spec%nPar2+1):Data%nPotPar2)) = -1.0d0
          end if
        else
          ! ... handle cases with unequal contributions
          ! ... work for the defined number or parents
          nCumMat = 0
          do i = 1, Spec%nPar2
            j = Data%nPotPar1 + Rank(i)
            ! ... these would have "zero" contributions and if we hit them then they need at least minimum contrib
            if (Chrom(j) .lt. Spec%LimitPar2Min) then
              Chrom(j) = Spec%LimitPar2Min ! set/fix to minimum usage @todo could consider penalising solution instead?
            end if
            ! ... but not above max allowed
            if (Chrom(j) .gt. Spec%LimitPar2Max) then
              Chrom(j) = Spec%LimitPar2Max ! set/fix to maximum usage @todo could consider penalising solution instead?
            end if
            ! ... accumulate and check if we reached Spec%nMat
            nCumMat = nCumMat + nint(Chrom(j)) ! internally real, externally integer
            if (nCumMat .ge. Spec%nMat) then
              ! ... there should be exactly Spec%nMat contributions
              if (nCumMat .gt. Spec%nMat) then
                Chrom(j) = Chrom(j) - dble(nCumMat - Spec%nMat)
                ! ... did we go below the minimum usage limit?
                if (nint(Chrom(j)) .lt. Spec%LimitPar2Min) then
                  TmpR = Spec%LimitPar2MinWeight * (Spec%LimitPar2Min - nint(Chrom(j)))
                  This%Criterion = This%Criterion + TmpR
                  if (Spec%LimitPar2MinWeight .lt. 0.0d0) then
                    This%Penalty = This%Penalty + abs(TmpR)
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
            Chrom(Data%nPotPar1+(Rank(i:Data%nPotPar2))) = 0.0d0
            ! ... negative (the same for all ind so no order)
            ! Chrom(Data%nPotPar1+(Rank(i:Data%nPotPar2))) = -1.0d0
            ! ... negative (variable with partially preserving order, i.e., ~large positives become ~large negatives)
            !     Found faster convergence than with properly decreasing negative values?
            !     I guess it adds some more randomness, i.e., it causes suffling of individuals that just did not make it onto the mating list.
            ! Chrom(Data%nPotPar1+(Rank(i:Data%nPotPar2))) = sign(Chrom(Data%nPotPar1+(Rank(i:Data%nPotPar2))), -1.0d0)
            ! ... negative and properly decreasing
            ! TmpR = maxval(Chrom(Data%nPotPar1+(Rank(i:Data%nPotPar2))))
            ! if (TmpR .gt. 0.0d0) then
            !     Chrom(Data%nPotPar1+(Rank(i:Data%nPotPar2))) = Chrom(Data%nPotPar1+(Rank(i:Data%nPotPar2))) - abs(TmpR)
            ! end if
            ! ... negative (random so no order)
            ! do j = i, Data%nPotPar2 ! @todo really need this loop?
            !   call random_number(RanNum)
            !   Chrom(Data%nPotPar1+Rank(j)) = -1.0d0 * RanNum
            ! end do
          end if
          ! ... Spec%nMat still not reached?
          do while (nCumMat .lt. Spec%nMat)
            ! ... add more contributions
            do i = Spec%nPar2, 1, -1 ! to bottom ranked selected individuals (to avoid local optima)
              j = Data%nPotPar1 + Rank(i)
              Chrom(j) = Chrom(j) + 1.0d0
              ! ... accumulate and check if we reached Spec%nMat
              nCumMat = nCumMat + 1
              if (nCumMat .eq. Spec%nMat) then
                ! To cater for real vs. integer issues
                TmpI = sum(nint(Chrom(Data%nPotPar1+Rank(1:Spec%nPar2))))
                if (TmpI .ne. Spec%nMat) then
                  if (TmpI .gt. Spec%nMat) then
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
        if (ChromInt(i) .lt. 0) then
          ChromInt(i) = 0
        end if
      end do
      ! ... map internal to external order
      nVecPar1(:) = ChromInt(1:Data%nPotPar1)
      if (.not.Spec%GenderFileGiven) then
        This%nVec(:) = nVecPar1(:)
      else
        This%nVec(Data%IdPotPar1) = nVecPar1(:)
      end if

      ! "Parent2"
      if (Spec%GenderFileGiven) then
        nVecPar2(:) = 0
        ! ... get integer values
        ChromInt(1:Data%nPotPar2) = nint(Chrom((Data%nPotPar1+1):(Data%nPotPar1+Data%nPotPar2)))
        ! ... remove negatives
        do i = 1, Data%nPotPar2
          if (ChromInt(i) .lt. 0) then
            ChromInt(i) = 0
          end if
        end do
        ! ... map internal to external order
        nVecPar2(:) = ChromInt(1:Data%nPotPar2)
        This%nVec(Data%IdPotPar2) = nVecPar2(:)
      end if

      This%xVec(:) = dble(This%nVec(:)) / (2 * Spec%nMat)

      ! --- PAGE ---

      if (Spec%PAGEPar) then
        if (.not.Spec%GenderFileGiven) then
          Rank(1:Data%nInd) = MrgRnk(Chrom((Data%nPotPar1+Spec%nMat+1):(Data%nPotPar1+Spec%nMat+Data%nInd)))
          This%GenomeEdit(Rank(Data%nInd:(Data%nInd-Spec%PAGEPar1Max+1):-1)) = 1.0d0 ! MrgRnk ranks small to large
        else
          if (Spec%PAGEPar1) then
            Rank(1:Data%nPotPar1) = MrgRnk(Chrom((Data%nPotPar1+Data%nPotPar2+Spec%nMat+1):(Data%nPotPar1+Data%nPotPar2+Spec%nMat+Data%nPotPar1)))
            This%GenomeEdit(Data%IdPotPar1(Rank(Data%nPotPar1:(Data%nPotPar1-Spec%PAGEPar1Max+1):-1))) = 1.0d0 ! MrgRnk ranks small to large
          end if
          if (Spec%PAGEPar2) then
            Rank(1:Data%nPotPar2) = MrgRnk(Chrom((Data%nPotPar1+Data%nPotPar2+Spec%nMat+Data%nPotPar1+1):(Data%nPotPar1+Data%nPotPar2+Spec%nMat+Data%nPotPar1+Data%nPotPar2)))
            This%GenomeEdit(Data%IdPotPar2(Rank(Data%nPotPar2:(Data%nPotPar2-Spec%PAGEPar2Max+1):-1))) = 1.0d0 ! MrgRnk ranks small to large
          end if
        end if
      end if

      ! --- Mate allocation ---

      MatPar2(:) = 0
      if (Spec%GenderFileGiven) then
        ! Distribute parent2 (=female) contributions into matings
        k = 0
        do i = 1, Data%nPotPar2 ! need to loop whole nVecPar2 as some entries are zero
          do j = 1, nVecPar2(i)
            k = k + 1
            MatPar2(k) = Data%IdPotPar2(i)
          end do
        end do
        ! Reorder parent2 contributions according to the rank of matings
        Rank(1:Spec%nMat) = MrgRnk(Chrom((Data%nPotPar1+Data%nPotPar2+1):(Data%nPotPar1+Data%nPotPar2+Spec%nMat)))
        MatPar2(:) = MatPar2(Rank(1:Spec%nMat))
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
        Rank(1:Spec%nMat) = MrgRnk(Chrom((Data%nPotPar1+1):(Data%nPotPar1+Spec%nMat)))
        MatPar2(:) = MatPar2(Rank(1:Spec%nMat))
      end if

      ! Pair the contributions (=Mating plan)
      k = Spec%nMat ! MrgRnk ranks small to large
      if (Spec%GenderFileGiven .or. Spec%SelfingAllowed) then
        ! When gender matters selfing can not happen (we have two distinct sets of parents,
        ! unless the user adds individuals of one sex in both sets) and when SelfingAllowed
        ! we do not need to care about it - faster code
        do i = 1, Data%nPotPar1
          do j = 1, nVecPar1(i)
            !if (k<2) print*, k, i, j, nVecPar1(i), Spec%nMat, sum(nVecPar1(:))
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
                  Chrom(Data%nPotPar1+Rank([k, l])) = Chrom(Data%nPotPar1+Rank([l, k]))
                  exit
                end if
              end do
              if (l .lt. 1) then ! Above loop ran out without finding a swap
                This%Criterion = This%Criterion + Spec%SelfingWeight
                if (Spec%SelfingWeight .lt. 0.0d0) then
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

      if (Spec%SelCriterionFileGiven) then
        !@todo save SelCriterion mean and sd in the data object and then compute this dot_product only once and
        !      compute This%SelCriterion as This%SelCriterion = This%SelIntensity * SelCriterionSD + SelCriterionMean
        This%SelCriterion = dot_product(This%xVec, Data%SelCriterion)
        This%SelIntensity = dot_product(This%xVec, Data%SelCriterionStand)
        if (Spec%PAGEPar) then
          !@todo as above
          This%SelCriterion = This%SelCriterion + dot_product(This%xVec, Data%SelCriterionPAGE(:)      * This%GenomeEdit(:))
          This%SelIntensity = This%SelIntensity + dot_product(This%xVec, Data%SelCriterionPAGEStand(:) * This%GenomeEdit(:))
        end if
        if (CritType .eq. "opt") then
          This%Criterion = This%Criterion + This%SelIntensity
        end if
      end if

      ! --- Generic individual values ---

      if (Spec%GenericIndCritGiven) then
        do j = 1, Spec%nGenericIndCrit
          TmpR = dot_product(This%xVec, Data%GenericIndCrit(:, j))
          This%GenericIndCrit(j) = TmpR
          TmpR = Spec%GenericIndCritWeight(j) * This%GenericIndCrit(j)
          This%Criterion = This%Criterion + TmpR
          if (Spec%GenericIndCritWeight(j) .lt. 0.0) then
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

      if      (CritType .eq. "min" .or. CritType .eq. "ran") then
        This%Criterion = This%Criterion - This%FutureCoancestryRanMate
      else if (CritType .eq. "opt") then
        ! We know the targeted rate of coancestry so we can work with relative values,
        ! which makes the TargetCoancestryRateWeight generic for ~any scenario.
        TmpR = This%CoancestryRateRanMate / Spec%TargetCoancestryRate
        if (TmpR .gt. 1.0d0) then
          ! Rate of coancestry for the solution is higher than the target
          TmpR = Spec%TargetCoancestryRateWeight * abs(1.0d0 - TmpR)
        else
          ! Rate of coancestry for the solution is lower than the target
          if (Spec%TargetCoancestryRateWeightBelow) then
            TmpR = Spec%TargetCoancestryRateWeight * abs(1.0d0 - abs(TmpR)) ! the second abs is to handle negative coancestry cases
          else
            TmpR = 0.0d0
          end if
        end if
        This%Criterion = This%Criterion + TmpR
        if (Spec%TargetCoancestryRateWeight .lt. 0.0d0) then
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
      ! We know the targeted rate of inbreeding so we can work with relative values,
      ! which makes the TargetInbreedingRateWeight generic for ~any scenario.
      TmpR = This%InbreedingRate / Spec%TargetInbreedingRate
      if (TmpR .gt. 1.0d0) then
        ! Rate of inbreeding for the solution is higher than the target
        TmpR = Spec%TargetInbreedingRateWeight * abs(1.0d0 - TmpR)
      else
        ! Rate of inbreeding for the solution is lower than the target
        if (Spec%TargetInbreedingRateWeightBelow) then
          TmpR = Spec%TargetInbreedingRateWeight * abs(1.0d0 - abs(TmpR)) ! the second abs is to handle negative inbreeding cases
        else
          TmpR = 0.0d0
        end if
      end if
      This%Criterion = This%Criterion + TmpR
      if (Spec%TargetInbreedingRateWeight .lt. 0.0d0) then
        This%Penalty = This%Penalty + abs(TmpR)
      end if

      ! --- Generic mating values ---

      if (Spec%GenericMatCritGiven) then
        do k = 1, Spec%nGenericMatCrit
          TmpR = 0.0d0
          if (Spec%GenderFileGiven) then
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
          This%Criterion = This%Criterion + TmpR
          if (Spec%GenericMatCritWeight(k) .lt. 0.0) then
            This%Penalty = This%Penalty + abs(TmpR)
          end if
        end do
      end if

      ! @todo how should we handle costs?
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Write head of the AlphaMate log
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !---------------------------------------------------------------------------
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

    !---------------------------------------------------------------------------
    !> @brief  Write the AlphaMate log
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !---------------------------------------------------------------------------
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
      if (Spec%GenericIndCritGiven) then
        if (Spec%GenericMatCritGiven) then
          if (present(String)) then
            write(Unit, StringFmt, Advance="No") trim(adjustl(String))
          end if
          write(Unit, Fmt) Iteration, AcceptRate, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate, This%GenericIndCrit, This%GenericMatCrit
        else
          if (present(String)) then
            write(Unit, StringFmt, Advance="No") trim(adjustl(String))
          end if
          write(Unit, Fmt) Iteration, AcceptRate, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate, This%GenericIndCrit
        end if
      else
        if (Spec%GenericMatCritGiven) then
          if (present(String)) then
            write(Unit, StringFmt, Advance="No") trim(adjustl(String))
          end if
          write(Unit, Fmt) Iteration, AcceptRate, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate,                      This%GenericMatCrit
        else
          if (present(String)) then
            write(Unit, StringFmt, Advance="No") trim(adjustl(String))
          end if
          write(Unit, Fmt) Iteration, AcceptRate, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate
        end if
      end if
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Write head of the AlphaMate log - for the swarm/population of solutions
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !---------------------------------------------------------------------------
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

    !---------------------------------------------------------------------------
    !> @brief  Write the AlphaMate log - for the swarm/population of solutions
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
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
      if (Spec%GenericIndCritGiven) then
        if (Spec%GenericMatCritGiven) then
          write(Unit, FMTLOGPOPUNIT) Iteration, i, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate, This%GenericIndCrit, This%GenericMatCrit
        else
          write(Unit, FMTLOGPOPUNIT) Iteration, i, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate, This%GenericIndCrit
        end if
      else
        if (Spec%GenericMatCritGiven) then
          write(Unit, FMTLOGPOPUNIT) Iteration, i, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate,                      This%GenericMatCrit
        else
          write(Unit, FMTLOGPOPUNIT) Iteration, i, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate
        end if
      end if
    end subroutine

    !###########################################################################
end module

!###############################################################################