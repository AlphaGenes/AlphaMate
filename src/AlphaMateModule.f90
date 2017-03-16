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
      real(real64) :: TargetCoancestryRate
      real(real64) :: TargetInbreedingRate
  end type

  !> @brief AlphaMate data
  type AlphaMateData
    ! Raw data
    type(RelMat) :: Coancestry
    type(InbVec) :: Inbreeding
    ! Data summaries
    type(DescStatReal64) :: InbreedingStat, SelCriterionStat, SelCriterionPAGEStat
    type(DescStatReal64), allocatable :: GenericIndValStat(:)
    type(DescStatMatrixReal64) :: CoancestryStat, CoancestryStatGenderDiff, CoancestryStatGender1, CoancestryStatGender2
    type(DescStatMatrixReal64), allocatable :: GenericMatValStat(:)
    ! Derived data
    real(real64) :: CurrentCoancestryRanMate, CurrentCoancestryRanMateNoSelf, CurrentCoancestryGenderMate
    real(real64) :: CurrentInbreeding
    real(real64) :: TargetCoancestryRanMate, TargetCoancestryRanMateNoSelf, TargetCoancestryGenderMate
    real(real64) :: TargetInbreeding
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
  end type

  type(AlphaMateSpec) :: Spec
  type(AlphaMateData) :: Data

  integer(int32) :: nInd, nMat, nPotMat, nPar, nPotPar1, nPotPar2, nMal, nFem, nPar1, nPar2, nFrontierSteps
  integer(int32) :: EvolAlgNSol, EvolAlgNGen, EvolAlgNGenBurnIn, EvolAlgNGenStop, EvolAlgNGenPrint, RanAlgStricter
  integer(int32) :: PAGEPar1Max, PAGEPar2Max, nGenericIndVal, nGenericMatVal
  integer(int32), allocatable :: Gender(:), IdPotPar1(:), IdPotPar2(:), IdPotParSeq(:)

  real(real64) :: LimitPar1Min, LimitPar1Max, LimitPar2Min, LimitPar2Max
  real(real64) :: EvolAlgStopTol, EvolAlgCRBurnIn, EvolAlgCRLate, EvolAlgFBase, EvolAlgFHigh1, EvolAlgFHigh2
  real(real64) :: CoancestryWeight, InbreedingWeight, SelfingWeight, LimitPar1Weight, LimitPar2Weight
  real(real64), allocatable :: GenericIndValWeight(:), GenericMatValWeight(:)
  real(real64) :: PAGEPar1Cost, PAGEPar2Cost
  real(real64), allocatable :: SelCriterion(:), SelCriterionStand(:), SelCriterionPAGE(:), SelCriterionPAGEStand(:)
  real(real64), allocatable :: TargetCoancestryRateFrontier(:), GenericIndValTmp(:), GenericMatValTmp(:)
  real(real64), allocatable :: GenericIndVal(:, :), GenericMatVal(:, :, :)

  logical :: NrmInsteadOfCoancestry
  logical :: ModeMin, ModeRan, ModeOpt, SelCriterionAvailable, GenderMatters, EqualizePar1, EqualizePar2
  logical :: SelfingAllowed, CoancestryWeightBellow, InbreedingWeightBellow, EvolAlgLogPop, EvaluateFrontier
  logical :: PAGE, PAGEPar1, PAGEPar2, GenericIndValAvailable, GenericMatValAvailable

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

  CHARACTER(len=100), PARAMETER  :: FMTFRONTEIRHEAD = "(a12, 8a22)"
  CHARACTER(len=100), PARAMETER  :: FMTFRONTEIR = "(a12, 8(1x, es21.13e3))"

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
    pure subroutine InitAlphaMateSpec
! TODO
    end subroutine

    !###########################################################################

    !-------------------------------------------------------------------------
    !> @brief  Read AlphaMate specifications from a file
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !-------------------------------------------------------------------------
    subroutine ReadAlphaMateSpec ! not pure due to IO
! TODO
    end subroutine

    !###########################################################################

    !-------------------------------------------------------------------------
    !> @brief  Read AlphaMate data from a file
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !-------------------------------------------------------------------------
    subroutine ReadAlphaMateData ! not pure due to IO
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

      character(len=1000) :: CoaMtxFile, SelCriterionFile, GenderFile, SeedFile
      character(len=1000) :: GenericIndValFile, GenericMatValFile
      character(len=100) :: DumC, IdCTmp, IdCTmp2

      write(STDOUT, "(a)") "--- Specifications ---"
      write(STDOUT, "(a)") " "

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

      ! --- CoancestryMatrixFile ---

      ! read(SpecUnit, *) DumC, CoaMtxFile
      ! write(STDOUT, "(a)") "CoancestryMatrixFile: "//trim(CoaMtxFile)
      ! NrmInsteadOfCoancestry = .false.

      ! --- NrmMatrixFile ---

      read(SpecUnit, *) DumC, CoaMtxFile
      write(STDOUT, "(a)") "NrmMatrixFile: "//trim(CoaMtxFile)
      NrmInsteadOfCoancestry = .true.

      ! --- SelCriterionFile ---

      read(SpecUnit, *) DumC, SelCriterionFile
      if (ToLower(trim(SelCriterionFile)) /= "none") then
        SelCriterionAvailable = .true.
        write(STDOUT, "(a)") "SelCriterionFile: "//trim(SelCriterionFile)
      else
        SelCriterionAvailable = .false.
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
            ! @todo might consider handling this better at some point
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
            ! @todo might consider handling this better at some point
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
          if (.not.SelCriterionAvailable) then
            write(STDERR, "(a)") "ERROR: Can not use PAGE when selection criterion file is not given!"
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
          if (.not.SelCriterionAvailable) then
            write(STDERR, "(a)") "ERROR: Can not use PAGE when selection criterion file is not given!"
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
          if (.not.SelCriterionAvailable) then
            write(STDERR, "(a)") "ERROR: Can not use PAGE when selection criterion file is not given!"
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

      ! --- TargetedRateOfCoancestry ---

      read(SpecUnit, *) DumC, Spec%TargetCoancestryRate, CoancestryWeight, DumC
      write(STDOUT, "(a)") "TargetedRateOfCoancestry: "//trim(Real2Char(Spec%TargetCoancestryRate, fmt=FMTREAL2CHAR))//&
        ", penalty weight "//trim(Real2Char(CoancestryWeight, fmt=FMTREAL2CHAR))//", mode "//trim(DumC)
      if      (ToLower(trim(DumC)) == "above") then
        CoancestryWeightBellow = .false.
      else if (ToLower(trim(DumC)) == "aboveandbellow") then
        CoancestryWeightBellow = .true.
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
      if (CoancestryWeight > 0.0) then
        write(STDERR, "(a)") "ERROR: Penalty weight for the targeted rate of coancestry should be zero or negative!"
        write(STDERR, "(a)") " "
        stop 1
      end if

      ! --- TargetedRateOfInbreeding (=inbreeding of a mating) ---

      read(SpecUnit, *) DumC, Spec%TargetInbreedingRate, InbreedingWeight, DumC
      write(STDOUT, "(a)") "TargetedRateOfInbreeding: "//trim(Real2Char(Spec%TargetInbreedingRate, fmt=FMTREAL2CHAR))//&
        ", penalty weight "//trim(Real2Char(InbreedingWeight, fmt=FMTREAL2CHAR))//", mode "//trim(DumC)
      if      (ToLower(trim(DumC)) == "above") then
        InbreedingWeightBellow = .false.
      else if (ToLower(trim(DumC)) == "aboveandbellow") then
        InbreedingWeightBellow = .true.
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
      if (InbreedingWeight > 0.0) then
        write(STDERR, "(a)") "ERROR: Penalty weight for the targeted rate of inbreeding should be zero or negative!"
        write(STDERR, "(a)") " "
        stop 1
      end if

      ! --- EvaluateFrontier ---

      ! @todo Should this be part of TargetedCoancestryRate??
      read(SpecUnit, *) DumC, DumC
      if      (ToLower(trim(DumC)) == "no") then
        EvaluateFrontier = .false.
        !write(STDOUT, "(a)") "EvaluateFrontier: no"
      else if (ToLower(trim(DumC)) == "yes") then
        EvaluateFrontier = .true.
        backspace(SpecUnit)
        read(SpecUnit, *) DumC, DumC, nFrontierSteps
        allocate(TargetCoancestryRateFrontier(nFrontierSteps))
        backspace(SpecUnit)
        read(SpecUnit, *) DumC, DumC, nFrontierSteps, TargetCoancestryRateFrontier(:)
        write(STDOUT, "("//Int2Char(1+nFrontierSteps)//"a)") "EvaluateFrontier: yes, #steps: "//trim(Int2Char(nFrontierSteps))//&
          ", rates of coancestry: ", (trim(Real2Char(TargetCoancestryRateFrontier(i), fmt=FMTREAL2CHAR)), i = 1, nFrontierSteps)
        if (any(TargetCoancestryRateFrontier(:) == 0.0)) then
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

      read(SpecUnit, *) DumC, EvolAlgNSol, EvolAlgNGen, EvolAlgNGenBurnIn, EvolAlgNGenStop, EvolAlgStopTol, EvolAlgNGenPrint, EvolAlgLogPop

      ! --- EvolutionaryAlgorithmParameters ---

      read(SpecUnit, *) DumC, EvolAlgCRBurnIn, EvolAlgCRLate, EvolAlgFBase, EvolAlgFHigh1, EvolAlgFHigh2

      ! --- RandomSearchIterationsStricter ---

      read(SpecUnit, *) DumC, RanAlgStricter

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

      ! --- Coancestry ---

      write(STDOUT, "(a)") "Coancestry (average identity of the four genome combinations of two individuals)"

      call Data%Coancestry%Read(File=CoaMtxFile)
      if (Data%Coancestry%nInd < nInd) then
        write(STDERR, "(a)") "ERROR: The coancestry matrix file has less rows than there are defined number of individuals!"
        write(STDERR, "(a)") "ERROR: Number of defined individuals:                       "//trim(Int2Char(nInd))
        write(STDERR, "(a)") "ERROR: Number of individuals in the coancestry matrix file: "//trim(Int2Char(Data%Coancestry%nInd))
        write(STDERR, "(a)") " "
      end if
      if (NrmInsteadOfCoancestry) then
        call Data%Coancestry%Nrm2Coancestry
      end if

      Data%CoancestryStat = DescStatSymMatrix(Data%Coancestry%Value(1:, 1:))
      if (GenderMatters) then
! TODO
        Data%CoancestryStatGenderDiff = DescStatSymMatrix(Data%Coancestry%Value(1:, 1:))
        Data%CoancestryStatGender1    = DescStatSymMatrix(Data%Coancestry%Value(1:, 1:))
        Data%CoancestryStatGender2    = DescStatSymMatrix(Data%Coancestry%Value(1:, 1:))
      end if

      ! Current
      Data%CurrentCoancestryRanMate       = Data%CoancestryStat%All%Mean
      Data%CurrentCoancestryRanMateNoSelf = Data%CoancestryStat%OffDiag%Mean
      if (GenderMatters) then
        Data%CurrentCoancestryGenderMate  = Data%CoancestryStatGenderDiff%All%Mean
      end if

      ! Obtain limit/target based on given rates
      ! F_t = DeltaF + (1 - DeltaF) * F_t-1
      Data%TargetCoancestryRanMate       = Spec%TargetCoancestryRate + (1.0d0 - Spec%TargetCoancestryRate) * Data%CurrentCoancestryRanMate
      Data%TargetCoancestryRanMateNoSelf = Spec%TargetCoancestryRate + (1.0d0 - Spec%TargetCoancestryRate) * Data%CurrentCoancestryRanMateNoSelf
      if (GenderMatters) then
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

      if (GenderMatters) then
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

      allocate(SelCriterion(nInd))
      allocate(SelCriterionStand(nInd))
      if (PAGE) then
        allocate(SelCriterionPAGE(nInd))
        allocate(SelCriterionPAGEStand(nInd))
      end if

      if (.not.SelCriterionAvailable) then
        SelCriterion(:) = 0.0d0
        SelCriterionStand(:) = 0.0d0
      else
        write(STDOUT, "(a)") "Selection criterion"
        nIndTmp = CountLines(SelCriterionFile)
        if (nIndTmp /= nInd) then
          write(STDERR, "(a)") "ERROR: Number of individuals in the selection criterion file and the coancestry matrix file is not the same!"
          write(STDERR, "(a)") "ERROR: Number of individuals in the coancestry matrix file:   "//trim(Int2Char(nInd))
          write(STDERR, "(a)") "ERROR: Number of individuals in the selection criterion file: "//trim(Int2Char(nIndTmp))
          write(STDERR, "(a)") " "
          stop 1
        end if
        open(newunit=SelCriterionUnit, file=trim(SelCriterionFile), status="old")
        do i = 1, nInd
          if (PAGE) then
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
          SelCriterion(j) = SelCriterionTmp
          if (PAGE) then
            SelCriterionPAGE(j) = SelCriterionTmp2
          end if
        end do
        close(SelCriterionUnit)

        Data%SelCriterionStat = DescStat(SelCriterion)
        SelCriterionStand(:) = (SelCriterion(:) - Data%SelCriterionStat%Mean) / Data%SelCriterionStat%SD
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

        if (PAGE) then
          ! must have the same scaling as selection criterion!!!!
          SelCriterionPAGEStand(:) = (SelCriterionPAGE(:) - Data%SelCriterionStat%Mean) / Data%SelCriterionStat%SD
          ! only the PAGE bit of SelCriterion
          SelCriterionPAGE(:) = SelCriterionPAGE(:) - SelCriterion(:)
          SelCriterionPAGEStand(:) = SelCriterionPAGEStand(:) - SelCriterionStand(:)
          Data%SelCriterionPAGEStat = DescStat(SelCriterionPAGE)
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

      allocate(Gender(nInd))

      Gender(:) = 0
      if (GenderMatters) then
        write(STDOUT, "(a)") "Gender"
        nIndTmp = CountLines(GenderFile)
        if (nIndTmp /= nInd) then
          write(STDERR, "(a)") "ERROR: Number of individuals in the gender file and the coancestry matrix file is not the same!"
          write(STDERR, "(a)") "ERROR: Number of individuals in the coancestry matrix file: "//trim(Int2Char(nInd))
          write(STDERR, "(a)") "ERROR: Number of individuals in the gender file:            "//trim(Int2Char(nIndTmp))
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
          j = FindLoc(IdCTmp, Data%Coancestry%OriginalId(1:))
          if (j == 0) then
            write(STDERR, "(a)") "ERROR: Individual "//trim(IdCTmp)//" from the gender file not present in the coancestry matrix file!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          Gender(j) = GenderTmp
        end do
        close(GenderUnit)

        write(STDOUT, "(a)") "  - number of   males in data: "//trim(Int2Char(nMal))
        write(STDOUT, "(a)") "  - number of females in data: "//trim(Int2Char(nFem))
        write(STDOUT, "(a)") " "

        if (nPar1 > nMal) then
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
        nPotMat = real(nPotPar1 * nPotPar1) / 2
        if (SelfingAllowed) then
          nPotMat = nint(nPotMat + real(nPotPar1) / 2)
        else
          nPotMat = nint(nPotMat - real(nPotPar1) / 2)
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
          write(STDERR, "(a)") "ERROR: Number of individuals in the generic individual values file and the coancestry matrix file is not the same!"
          write(STDERR, "(a)") "ERROR: Number of individuals in the coancestry matrix file:         "//trim(Int2Char(nInd))
          write(STDERR, "(a)") "ERROR: Number of individuals in the generic individual values file: "//trim(Int2Char(nIndTmp))
          write(STDERR, "(a)") " "
          stop 1
        end if
        allocate(GenericIndVal(nInd, nGenericIndVal))
        allocate(GenericIndValTmp(nGenericIndVal))
        GenericIndVal(:, :) = 0.0d0
        open(newunit=GenericIndValUnit, file=GenericIndValFile, status="unknown")
        do i = 1, nInd
          read(GenericIndValUnit, *) IdCTmp, GenericIndValTmp(:)
          j = FindLoc(IdCTmp, Data%Coancestry%OriginalId(1:))
          if (j == 0) then
            write(STDERR, "(a)") "ERROR: Individual "//trim(IdCTmp)//" from the generic individual values file not present in the coancestry matrix file!"
            write(STDERR, "(a)") " "
            stop 1
          end if
          GenericIndVal(j, :) = GenericIndValTmp(:)
        end do
        close(GenericIndValUnit)

        allocate(Data%GenericIndValStat(nGenericIndVal))
        do j = 1, nGenericIndVal
          write(STDOUT, "(a)") "  - column "//trim(Int2Char(j))
          Data%GenericIndValStat(i) = DescStat(GenericIndVal(:, j))
          write(STDOUT, "(a)") "    - average: "//trim(Real2Char(Data%GenericIndValStat(i)%Mean, fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(Data%GenericIndValStat(i)%SD,   fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(Data%GenericIndValStat(i)%Min,  fmt=FMTREAL2CHAR))
          write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(Data%GenericIndValStat(i)%Max,  fmt=FMTREAL2CHAR))
        end do
        write(STDOUT, "(a)") " "
      end if

      ! --- Generic mating values ---

      if (GenericMatValAvailable) then
        write(STDOUT, "(a)") "Generic mating values"
        DumI = CountLines(GenericMatValFile)
        if (DumI /= nPotMat) then
          write(STDERR, "(a)") "ERROR: Number of matings in the generic mating values file and the number of all potential matings is not the same!"
          write(STDERR, "(a)") "ERROR: Number of all potential matings:                         "//trim(Int2Char(nPotMat))
          write(STDERR, "(a)") "ERROR: Number of individuals in the generic mating values file: "//trim(Int2Char(DumI))
          write(STDERR, "(a)") " "
          stop 1
        end if
        allocate(GenericMatVal(nPotPar1, nPotPar2, nGenericMatVal))
        allocate(GenericMatValTmp(nGenericMatVal))
        GenericMatVal(:, :, :) = 0.0d0
        open(newunit=GenericMatValUnit, file=GenericMatValFile, status="unknown")
        do i = 1, nPotMat
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
          if (GenderMatters) then
            l = FindLoc(j, IdPotPar1)
            if (l == 0) then
              write(STDERR, "(a)") "ERROR: Individual "//trim(IdCTmp)//" from the first column in the generic mating values file should be a male!"
              write(STDERR, "(a)") "ERROR: Generic mating values file (line "//trim(Int2Char(i))//"): "//trim(IdCTmp)//" "//trim(IdCTmp2)
              write(STDERR, "(a)") " "
              stop 1
            end if
            m = FindLoc(k, IdPotPar2)
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
            GenericMatVal(IdPotParSeq(IdPotPar1(l)), IdPotParSeq(IdPotPar2(m)), :) = GenericMatValTmp(:)
          else
            ! fill lower-triangle (half-diallel)
            GenericMatVal(maxval([j, k]), minval([j, k]), :) = GenericMatValTmp(:)
          end if
        end do
        close(GenericMatValUnit)

        allocate(Data%GenericMatValStat(nGenericMatVal))
        do k = 1, nGenericMatVal
          write(STDOUT, "(a)") "  - column "//trim(Int2Char(k))
          if (GenderMatters) then
            Data%GenericMatValStat(k) = DescStatMatrix(GenericMatVal(:, :, k))
            write(STDOUT, "(a)") "    - average: "//trim(Real2Char(Data%GenericMatValStat(k)%All%Mean, fmt=FMTREAL2CHAR))
            write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(Data%GenericMatValStat(k)%All%SD,   fmt=FMTREAL2CHAR))
            write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(Data%GenericMatValStat(k)%All%Min,  fmt=FMTREAL2CHAR))
            write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(Data%GenericMatValStat(k)%All%Max,  fmt=FMTREAL2CHAR))
          else
            if (SelfingAllowed) then
              Data%GenericMatValStat(k) = DescStatLowTriMatrix(GenericMatVal(:, :, k))
              write(STDOUT, "(a)") "    - average: "//trim(Real2Char(Data%GenericMatValStat(k)%All%Mean, fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - st.dev.: "//trim(Real2Char(Data%GenericMatValStat(k)%All%SD,   fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - minimum: "//trim(Real2Char(Data%GenericMatValStat(k)%All%Min,  fmt=FMTREAL2CHAR))
              write(STDOUT, "(a)") "    - maximum: "//trim(Real2Char(Data%GenericMatValStat(k)%All%Max,  fmt=FMTREAL2CHAR))
            end if
              Data%GenericMatValStat(k) = DescStatLowTriMatrix(GenericMatVal(:, :, k), Diag=.false.)
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
    subroutine SetupColNamesAndFormats
      implicit none
      integer(int32) :: nCol, nColTmp, i

      ! --- Optimisation log ---

      nCol = 10
      if (GenericIndValAvailable) then
        nCol = nCol + nGenericIndVal
      end if
      if (GenericMatValAvailable) then
        nCol = nCol + nGenericMatVal
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
      nColTmp = 10
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
    subroutine AlphaMateSearch

      implicit none

      integer(int32) :: nParam, k, FrontierUnit

      real(real64) :: HoldTargetCoancestryRanMate, HoldTargetCoancestryRate
      real(real64), allocatable :: InitEqual(:, :)

      character(len=1000) :: LogFile, LogPopFile, ContribFile, MatingFile
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

      ! --- Optimise contributions for minimal future coancestry/inbreeding ---

      if (ModeMin) then
        write(STDOUT, "(a)") "--- Optimise contributions for minimal future coancestry/inbreeding --- "
        write(STDOUT, "(a)") " "

        LogFile     = "OptimisationLogMinimalInbreeding.txt"
        LogPopFile  = "OptimisationLogPopMinimalInbreeding.txt"
        ContribFile = "IndividualResultsMinimalInbreeding.txt"
        MatingFile  = "MatingResultsMinimalInbreeding.txt"

        allocate(InitEqual(nParam, nint(EvolAlgNSol * 0.1)))
        InitEqual(:, :) = 1.0d0 ! A couple of solutions that would give equal contributions for everybody

        call DifferentialEvolution(nParam=nParam, nSol=EvolAlgNSol, Init=InitEqual, nGen=EvolAlgNGen, nGenBurnIn=EvolAlgNGenBurnIn, &
          nGenStop=EvolAlgNGenStop, StopTolerance=EvolAlgStopTol, nGenPrint=EvolAlgNGenPrint, LogFile=LogFile, LogPop=EvolAlgLogPop, LogPopFile=LogPopFile, &
          CritType="min", CRBurnIn=EvolAlgCRBurnIn, CRLate=EvolAlgCRLate, FBase=EvolAlgFBase, FHigh1=EvolAlgFHigh1, FHigh2=EvolAlgFHigh2, &
          BestSol=SolMin)

        deallocate(InitEqual)

        call SaveSolution(SolMin, ContribFile, MatingFile)
      end if

      ! --- Evaluate random mating ---

      if (ModeRan) then
        write(STDOUT, "(a)") "--- Evaluate random mating --- "
        write(STDOUT, "(a)") " "

        LogFile = "OptimisationLogRandomMating.txt"

        allocate(InitEqual(nParam, nint(EvolAlgNSol * 0.1)))
        InitEqual(:, :) = 1.0d0 ! A couple of solutions that would give equal contributions for everybody

        call RandomSearch(Mode="avg", nParam=nParam, Init=InitEqual, nSamp=EvolAlgNSol*EvolAlgNGen*RanAlgStricter, nSampStop=EvolAlgNGenStop*RanAlgStricter, &
          StopTolerance=EvolAlgStopTol/RanAlgStricter, nSampPrint=EvolAlgNGenPrint, LogFile=LogFile, CritType="ran", BestSol=SolRan)

        deallocate(InitEqual)
      end if

      ! --- Optimise contributions for maximal value with constraint on future coancestry/inbreeding ---

      if (ModeOpt) then
        write(STDOUT, "(a)") "--- Optimise contributions for maximal value with constraint on future coancestry/inbreeding ---"
        write(STDOUT, "(a)") " "

        LogFile     = "OptimisationLogMaximalValue.txt"
        LogPopFile  = "OptimisationLogPopMaximalValue.txt"
        ContribFile = "IndividualResultsMaximalValue.txt"
        MatingFile  = "MatingResultsMaximalValue.txt"

        ! @todo add some clever initial values, say:
        !       - equal contributions for top 2/3 or 1/2 of BV distribution,
        !       - decreasing contributions with decreasing value
        !       - SDP solution, ...?
        call DifferentialEvolution(nParam=nParam, nSol=EvolAlgNSol, nGen=EvolAlgNGen, nGenBurnIn=EvolAlgNGenBurnIn, &
          nGenStop=EvolAlgNGenStop, StopTolerance=EvolAlgStopTol, nGenPrint=EvolAlgNGenPrint, LogFile=LogFile, LogPop=EvolAlgLogPop, LogPopFile=LogPopFile, &
          CritType="opt", CRBurnIn=EvolAlgCRBurnIn, CRLate=EvolAlgCRLate, FBase=EvolAlgFBase, FHigh1=EvolAlgFHigh1, FHigh2=EvolAlgFHigh2, &
          BestSol=SolOpt)

        call SaveSolution(SolOpt, ContribFile, MatingFile)
      end if

      ! --- Evaluate the full frontier ---

      if (EvaluateFrontier) then
        write(STDOUT, "(a)") "--- Evaluate the frontier ---"
        write(STDOUT, "(a)") " "

        CoancestryWeightBellow = .true. ! we want to target certain rates of coancestry

        open(newunit=FrontierUnit, file="Frontier.txt", status="unknown")
        !                                     1234567890123456789012
        write(FrontierUnit, FMTFRONTEIRHEAD) "   Iteration", &
                                             "             Criterion", &
                                             "             Penalties", &
                                             "          SelCriterion", &
                                             "          SelIntensity", &
                                             "            Coancestry", &
                                             "        CoancestryRate", &
                                             "            Inbreeding", &
                                             "        InbreedingRate"
! @todo add the generic stuff from the log? Just call This%Log?
        if (ModeMin) then
          DumC = "Min"
          write(FrontierUnit, FMTFRONTEIR) adjustl(DumC), SolMin%Criterion, SolMin%Penalty, SolMin%SelCriterion, SolMin%SelIntensity, SolMin%FutureCoancestryRanMate, SolMin%CoancestryRateRanMate, SolMin%FutureInbreeding, SolMin%InbreedingRate
        end if
        if (ModeRan) then
          DumC = "Ran"
          write(FrontierUnit, FMTFRONTEIR) adjustl(DumC), SolRan%Criterion, SolRan%Penalty, SolRan%SelCriterion, SolRan%SelIntensity, SolRan%FutureCoancestryRanMate, SolRan%CoancestryRateRanMate, SolRan%FutureInbreeding, SolRan%InbreedingRate
        end if
        if (ModeOpt) then
          DumC = "Opt"
          write(FrontierUnit, FMTFRONTEIR) adjustl(DumC), SolOpt%Criterion, SolOpt%Penalty, SolOpt%SelCriterion, SolOpt%SelIntensity, SolOpt%FutureCoancestryRanMate, SolOpt%CoancestryRateRanMate, SolOpt%FutureInbreeding, SolOpt%InbreedingRate
        end if

        ! Hold old results
        HoldTargetCoancestryRanMate = Data%TargetCoancestryRanMate
        HoldTargetCoancestryRate = Spec%TargetCoancestryRate

        ! Evaluate
        do k = 1, nFrontierSteps
          Spec%TargetCoancestryRate = TargetCoancestryRateFrontier(k)
          ! F_t = DeltaF + (1 - DeltaF) * F_t-1
! TODO
          Data%TargetCoancestryRanMate = Spec%TargetCoancestryRate + (1.0d0 - Spec%TargetCoancestryRate) * Data%CurrentCoancestryRanMate
          write(STDOUT, "(a)") "Step "//trim(Int2Char(k))//" out of "//trim(Int2Char(nFrontierSteps))//&
                               " for the rate of coancestry "//trim(Real2Char(Spec%TargetCoancestryRate, fmt=FMTREAL2CHAR))//&
                               " (=targeted coancestry "//trim(Real2Char(Data%TargetCoancestryRanMate, fmt=FMTREAL2CHAR))//")"
          write(STDOUT, "(a)") ""

          LogFile     = "OptimisationLogFrontier"//trim(Int2Char(k))//".txt"
          LogPopFile  = "OptimisationLogPopFrontier"//trim(Int2Char(k))//".txt"
          ContribFile = "IndividualResultsFrontier"//trim(Int2Char(k))//".txt"
          MatingFile  = "MatingResultsFrontier"//trim(Int2Char(k))//".txt"

          call DifferentialEvolution(nParam=nParam, nSol=EvolAlgNSol, nGen=EvolAlgNGen, nGenBurnIn=EvolAlgNGenBurnIn, &
            nGenStop=EvolAlgNGenStop, StopTolerance=EvolAlgStopTol, nGenPrint=EvolAlgNGenPrint, LogFile=LogFile, LogPop=EvolAlgLogPop, LogPopFile=LogPopFile, &
            CritType="opt", CRBurnIn=EvolAlgCRBurnIn, CRLate=EvolAlgCRLate, FBase=EvolAlgFBase, FHigh1=EvolAlgFHigh1, FHigh2=EvolAlgFHigh2, &
            BestSol=Sol)

          ! @todo add the generic stuff from the log? Just call This%Log?
          DumC = "Frontier"//trim(Int2Char(k))
          write(FrontierUnit, FMTFRONTEIR) adjustl(DumC), Sol%Criterion, Sol%Penalty, Sol%SelCriterion, Sol%SelIntensity, Sol%FutureCoancestryRanMate, Sol%CoancestryRateRanMate, Sol%FutureInbreeding, Sol%InbreedingRate

          call SaveSolution(Sol, ContribFile, MatingFile)

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
    subroutine SaveSolution(Sol, ContribFile, MatingFile) ! not pure due to IO

! TODO: add stdout when ContribFile and MatingFile are not given
      implicit none

      ! Arguments
      type(AlphaMateSol) :: Sol
      character(len=*)   :: ContribFile
      character(len=*)   :: MatingFile

      ! Other
      integer(int32) :: i, j, ContribUnit, MatingUnit, Rank(nInd)

      ! @todo should we have constant output no matter which options are switched on?
      open(newunit=ContribUnit, file=ContribFile, status="unknown")
      Rank = MrgRnk(Sol%nVec + SelCriterionStand / 100.0d0) ! @todo is this really good sorting?
      if (.not.PAGE) then
        !                                        1234567890123456789012
        write(ContribUnit, FMTCONTRIBUTIONHEAD) "          Id", &
                                                "      Gender", &
                                                "SelCriterion", &
                                                "AvCoancestry", &
                                                "Contribution", &
                                                "     nMating"
        do i = nInd, 1, -1 ! MrgRnk ranks small to large
          j = Rank(i)
          write(ContribUnit, FMTCONTRIBUTION) Data%Coancestry%OriginalId(j), Gender(j), SelCriterion(j), &
                                              sum(Data%Coancestry%Value(1:, j)) / nInd, &
                                              Sol%xVec(j), Sol%nVec(j)
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
        do i = nInd, 1, -1 ! MrgRnk ranks small to large
          j = Rank(i)
          write(ContribUnit, FMTCONTRIBUTIONEDIT) Data%Coancestry%OriginalId(j), Gender(j), SelCriterion(j), &
                                                  sum(Data%Coancestry%Value(1:, j)) / nInd, &
                                                  Sol%xVec(j), Sol%nVec(j), &
                                                  nint(Sol%GenomeEdit(j)), SelCriterion(j) + Sol%GenomeEdit(j) * SelCriterionPAGE(j)
        end do
      end if
      close(ContribUnit)

      open(newunit=MatingUnit, file=MatingFile, status="unknown")
      !                                 1234567890123456789012
      write(MatingUnit, FMTMATINGHEAD) "      Mating", &
                                       "     Parent1", &
                                       "     Parent2"
      do i = 1, nMat
        write(MatingUnit, FMTMATING) i, Data%Coancestry%OriginalId(Sol%MatingPlan(1, i)), Data%Coancestry%OriginalId(Sol%MatingPlan(2, i))
      end do
      close(MatingUnit)
    end subroutine

    !###########################################################################

    !-------------------------------------------------------------------------
    !> @brief  Initialize AlphaMate solution
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !-------------------------------------------------------------------------
    subroutine InitialiseAlphaMateSol(This)
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
      This%MatingPlan(:, :) = 0
      if (PAGE) then
        allocate(This%GenomeEdit(nInd))
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
          Out%Criterion = In%Criterion
          Out%Penalty = In%Penalty
          Out%SelCriterion = In%SelCriterion
          Out%SelIntensity = In%SelIntensity
          Out%FutureCoancestryRanMate = In%FutureCoancestryRanMate
          Out%CoancestryRateRanMate = In%CoancestryRateRanMate
          Out%FutureInbreeding = In%FutureInbreeding
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
        class default
          write(STDERR, "(a)") "ERROR: Both Out and In must be of the AlphaMateSol class for AlphaMate!"
          write(STDERR, "(a)") " "
          stop 1
      end select
    end subroutine

    !###########################################################################

! TODO: keep the best random solution?
    !-------------------------------------------------------------------------
    !> @brief  Update mean of AlphaMate solution (when performing random search)
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !-------------------------------------------------------------------------
    subroutine UpdateMeanAlphaMateSol(This, Add, n)
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
        class default
          write(STDERR, "(a)") "ERROR: Both This and Add must be of the AlphaMateSol class for AlphaMate!"
          write(STDERR, "(a)") " "
          stop 1
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
      integer(int32) :: i, j, k, l, g, nCumMat, Rank(nInd), ChromInt(nInd), MatPar2(nMat)
      integer(int32) :: nVecPar1(nPotPar1), nVecPar2(nPotPar2), TmpMin, TmpMax, TmpI

      real(real64) :: TmpVec(nInd, 1), TmpR, RanNum

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

      ! The approach below assures that we have nMat contributions for each of
      ! the two parent sets. It does this by ranking internal solution values and
      ! traverses from top to the defined number of parents checking when the sum of
      ! interegrised values gives nMat. If values bellow 0.5 are found, they are
      ! changed to 1 contribution. If this still does not give nMat, then we start

      ! @todo least contributing? Those that already contribute or thos that do not contribute at all?

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
        !@todo decreasing option in MrgRnk?
      end if
      ! ... handle cases with equalized contributions
      if (EqualizePar1) then
        if (nPar1 == nPotPar1) then
          ! ... set integers to all the values (no need for sorting here)
          Chrom(1:nPotPar1) = dble(nMat * g) / nPar1
        else
          ! ... set integers to the top values
          Chrom(Rank(1:nPar1)) = dble(nMat * g) / nPar1
          !@todo anything better to preserve the order of non contributing individuals? See below!
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
          ! do j = i, nPotPar1 ! @todo really need this loop?
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
            Chrom((nPotPar1+1):(nPotPar1+nPotPar2)) = dble(nMat) / nPar2
          else
            ! ... set integers to the top values
            Chrom(nPotPar1+Rank(1:nPar2)) = dble(nMat) / nPar2
            ! @todo anything better to preserve the order of non contributing individuals? See below!
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
            ! do j = i, nPotPar2 ! @todo really need this loop?
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

      ! --- Contributions (nVec & xVec) ---

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

      This%xVec(:) = dble(This%nVec(:)) / (2 * nMat)

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
            This%MatingPlan(1, k) = IdPotPar1(i)
            This%MatingPlan(2, k) = MatPar2(k)
            k = k - 1
          end do
        end do
      else
        ! When gender does not matter, selfing can happen (we have one set of parents)
        ! and when selfing is not allowed we need to avoid it - slower code
        do i = 1, nPotPar1
          do j = 1, nVecPar1(i)
            This%MatingPlan(1, k) = IdPotPar1(i)
            if (MatPar2(k) == IdPotPar1(i)) then
              ! Try to avoid selfing by swapping the MatPar2 and Rank elements
              do l = k, 1, -1
                if (MatPar2(l) /= IdPotPar1(i)) then
                  MatPar2([k, l]) = MatPar2([l, k])
                  Chrom(nPotPar1+Rank([k, l])) = Chrom(nPotPar1+Rank([l, k]))
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
            This%MatingPlan(2, k) = MatPar2(k)
            k = k - 1
          end do
        end do
      end if

      ! --- Contribution value ---

      if (SelCriterionAvailable) then
        !@todo save SelCriterion mean and sd in the data object and then compute this dot_product only once and
        !      compute This%SelCriterion as This%SelCriterion = This%SelIntensity * SelCriterionSD + SelCriterionMean
        This%SelCriterion = dot_product(This%xVec, SelCriterion)
        This%SelIntensity = dot_product(This%xVec, SelCriterionStand)
        if (PAGE) then
          !@todo as above
          This%SelCriterion = This%SelCriterion + dot_product(This%xVec, SelCriterionPAGE(:)      * This%GenomeEdit(:))
          This%SelIntensity = This%SelIntensity + dot_product(This%xVec, SelCriterionPAGEStand(:) * This%GenomeEdit(:))
        end if
        if (CritType == "opt") then
          This%Criterion = This%Criterion + This%SelIntensity
        end if
      end if

      ! --- Generic individual values ---

      if (GenericIndValAvailable) then
        do j = 1, nGenericIndVal
          TmpR = dot_product(This%xVec, GenericIndVal(:, j))
          This%GenericIndVal(j) = TmpR
          TmpR = GenericIndValWeight(j) * This%GenericIndVal(j)
          This%Criterion = This%Criterion + TmpR
          if (GenericIndValWeight(j) .lt. 0.0) then
            This%Penalty = This%Penalty + abs(TmpR)
          end if
        end do
      end if

      ! --- Coancestry ---

      ! x'C
      do i = 1, nInd
        TmpVec(i, 1) = dot_product(This%xVec, Data%Coancestry%Value(1:, i))
      end do
      ! @todo consider using matmul instead of repeated dot_product?
      ! @todo consider using BLAS/LAPACK - perhaps non-symmetric is more optimised?
      ! Matrix multiplication with a symmetric matrix using BLAS routine
      ! (it was ~5x slower than the above with 1.000 individuals, might be benefical with larger cases)
      ! http://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3.html#ga253c8edb8b21d1b5b1783725c2a6b692
      ! call dsymm(side="l", uplo="l", m=nInd, n=1, alpha=1.0d0, A=CoaMtx, lda=nInd, b=This%xVec, ldb=nInd, beta=0, c=TmpVec, ldc=nInd)
      ! call dsymm(     "l",      "l",   nInd,   1,       1.0d0,   CoaMtx,     nInd,   This%xVec,     nInd,      0,   TmpVec,     nInd)

      ! x'Cx
!TODO
      This%FutureCoancestryRanMate = dot_product(TmpVec(:, 1), This%xVec)

      ! dF = (F_t+1 - F_t) / (1 - F_t)
!TODO
      This%CoancestryRateRanMate = (This%FutureCoancestryRanMate - Data%CurrentCoancestryRanMate) / (1.0d0 - Data%CurrentCoancestryRanMate)

      if      (CritType == "min" .or. CritType == "ran") then
        This%Criterion = This%Criterion - This%FutureCoancestryRanMate
      else if (CritType == "opt") then
        ! We know the targeted rate of coancestry so we can work with relative values,
        ! which makes the CoancestryWeight generic for ~any scenario.
        TmpR = This%CoancestryRateRanMate / Spec%TargetCoancestryRate
        if (TmpR > 1.0d0) then
          ! Rate of coancestry for the solution is higher than the target
          TmpR = CoancestryWeight * abs(1.0d0 - TmpR)
        else
          ! Rate of coancestry for the solution is lower than the target
          if (CoancestryWeightBellow) then
            TmpR = CoancestryWeight * abs(1.0d0 - abs(TmpR)) ! the second abs is to handle negative coancestry cases
          else
            TmpR = 0.0d0
          end if
        end if
        This%Criterion = This%Criterion + TmpR
        if (CoancestryWeight < 0.0d0) then
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
      do j = 1, nMat
        ! Lower triangle to speedup lookup
        TmpMax = maxval(This%MatingPlan(:, j))
        TmpMin = minval(This%MatingPlan(:, j))
        TmpR = TmpR + 0.5d0 * Data%Coancestry%Value(TmpMax, TmpMin)
      end do
      This%FutureInbreeding = TmpR / nMat
      This%InbreedingRate = (This%FutureInbreeding - Data%CurrentInbreeding) / (1.0d0 - Data%CurrentInbreeding)
      ! We know the targeted rate of inbreeding so we can work with relative values,
      ! which makes the InbreedingWeight generic for ~any scenario.
      TmpR = This%InbreedingRate / Spec%TargetInbreedingRate
      if (TmpR > 1.0d0) then
        ! Rate of inbreeding for the solution is higher than the target
        TmpR = InbreedingWeight * abs(1.0d0 - TmpR)
      else
        ! Rate of inbreeding for the solution is lower than the target
        if (InbreedingWeightBellow) then
          TmpR = InbreedingWeight * abs(1.0d0 - abs(TmpR)) ! the second abs is to handle negative inbreeding cases
        else
          TmpR = 0.0d0
        end if
      end if
      This%Criterion = This%Criterion + TmpR
      if (InbreedingWeight < 0.0d0) then
        This%Penalty = This%Penalty + abs(TmpR)
      end if

      ! --- Generic mating values ---

      if (GenericMatValAvailable) then
        do k = 1, nGenericMatVal
          TmpR = 0.0d0
          if (GenderMatters) then
            do j = 1, nMat
              TmpR = TmpR + GenericMatVal(IdPotParSeq(This%MatingPlan(1, j)), &
                                          IdPotParSeq(This%MatingPlan(2, j)), k)
            end do
          else
            do j = 1, nMat
              TmpMax = maxval(This%MatingPlan(:, j))
              TmpMin = minval(This%MatingPlan(:, j))
              TmpR = TmpR + GenericMatVal(TmpMax, TmpMin, k)
            end do
          end if
          This%GenericMatVal(k) = TmpR / nMat
          TmpR = GenericMatValWeight(k) * This%GenericMatVal(k)
          This%Criterion = This%Criterion + TmpR
          if (GenericMatValWeight(k) < 0.0) then
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
    subroutine LogHeadAlphaMateSol(LogUnit)
      implicit none
      integer(int32), intent(in), optional :: LogUnit
      write(STDOUT, FMTLOGSTDOUTHEAD)   COLNAMELOGSTDOUT(:)
      if (present(LogUnit)) then
        write(LogUnit, FMTLOGUNITHEAD) COLNAMELOGUNIT(:)
      end if
    end subroutine

    !###########################################################################

! TODO: should this be a method for AlphaMateSol?
    !-------------------------------------------------------------------------
    !> @brief  Write the AlphaMate log
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !-------------------------------------------------------------------------
    subroutine LogAlphaMateSol(This, LogUnit, Gen, AcceptRate)
      implicit none
      class(AlphaMateSol), intent(in)      :: This
      integer(int32), intent(in), optional :: LogUnit
      integer(int32), intent(in)           :: Gen
      real(real64), intent(in)             :: AcceptRate
      if (GenericIndValAvailable) then
        if (GenericMatValAvailable) then
          write(STDOUT, FMTLOGSTDOUT)  Gen, AcceptRate, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate, This%GenericIndVal, This%GenericMatVal
        else
          write(STDOUT, FMTLOGSTDOUT)  Gen, AcceptRate, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate, This%GenericIndVal
        end if
      else
        if (GenericMatValAvailable) then
          write(STDOUT, FMTLOGSTDOUT)  Gen, AcceptRate, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate,                     This%GenericMatVal
        else
          write(STDOUT, FMTLOGSTDOUT)  Gen, AcceptRate, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate
        end if
      end if
      if (present(LogUnit)) then
        if (GenericIndValAvailable) then
          if (GenericMatValAvailable) then
            write(LogUnit, FMTLOGUNIT) Gen, AcceptRate, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate, This%GenericIndVal, This%GenericMatVal
          else
            write(LogUnit, FMTLOGUNIT) Gen, AcceptRate, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate, This%GenericIndVal
          end if
        else
          if (GenericMatValAvailable) then
            write(LogUnit, FMTLOGUNIT) Gen, AcceptRate, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate,                     This%GenericMatVal
          else
            write(LogUnit, FMTLOGUNIT) Gen, AcceptRate, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate
          end if
        end if
      end if
    end subroutine

    !###########################################################################

    !-------------------------------------------------------------------------
    !> @brief  Write head of the AlphaMate log - for the swarm/population of solutions
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !-------------------------------------------------------------------------
    subroutine LogPopHeadAlphaMateSol(LogPopUnit)
      implicit none
      integer(int32), intent(in) :: LogPopUnit
      write(LogPopUnit, FMTLOGUNITHEAD) COLNAMELOGPOPUNIT(:)
    end subroutine

    !###########################################################################

! TODO: should this be a method for AlphaMateSol?
    !-------------------------------------------------------------------------
    !> @brief  Write the AlphaMate log - for the swarm/population of solutions
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 16, 2017
    !-------------------------------------------------------------------------
    subroutine LogPopAlphaMateSol(This, LogPopUnit, Gen, i)
      implicit none
      class(AlphaMateSol), intent(in) :: This
      integer(int32), intent(in)      :: LogPopUnit
      integer(int32), intent(in)      :: Gen
      integer(int32), intent(in)      :: i

      if (GenericIndValAvailable) then
        if (GenericMatValAvailable) then
          write(LogPopUnit, FMTLOGPOPUNIT) Gen, i, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate, This%GenericIndVal, This%GenericMatVal
        else
          write(LogPopUnit, FMTLOGPOPUNIT) Gen, i, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate, This%GenericIndVal
        end if
      else
        if (GenericMatValAvailable) then
          write(LogPopUnit, FMTLOGPOPUNIT) Gen, i, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate,                     This%GenericMatVal
        else
          write(LogPopUnit, FMTLOGPOPUNIT) Gen, i, This%Criterion, This%Penalty, This%SelCriterion, This%SelIntensity, This%FutureCoancestryRanMate, This%CoancestryRateRanMate, This%FutureInbreeding, This%InbreedingRate
        end if
      end if
    end subroutine

    !###########################################################################
end module

!###############################################################################