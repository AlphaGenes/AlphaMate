
#ifdef OS_UNIX
#define MKDIR "mkdir -p"
#define RMDIR "rm -r"
#else
#define MKDIR "md"
#define RMDIR "rmdir /S"
#endif

!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     AlphaMate.f90
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

program AlphaMate
  use ISO_Fortran_Env, STDIN => input_unit, STDOUT => output_unit, STDERR => error_unit
  use IFPort, only : SystemQQ
  use AlphaHouseMod, only : Int2Char, PrintElapsedTime
  use AlphaMateModule

  implicit none

  real(real32) :: StartTime, EndTime

  logical :: Success

  call cpu_time(StartTime)
  call AlphaMateTitle

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

  call cpu_time(EndTime)
  call AlphaMateTitle
  call PrintElapsedTime(StartTime, EndTime)
end program

!###############################################################################

!@todo make spec type
!@todo make data type
!@todo etc
!  use ConstantModule, only : FILELENGTH
!  integer(int32) :: nArg
!  character(len=FILELENGTH) :: SpecFile
!  type(AlphaRelateSpec) :: Spec
!  type(AlphaRelateData) :: Data

!  write(STDOUT, "(a)") ""
!  write(STDOUT, "(a)") " Processing specifications ..."
!  nArg = command_argument_count()
!  if (nArg > 0) then
!    call get_command_argument(1, SpecFile)
!  else
!    SpecFile = "AlphaRelateSpec.txt"
!  end if
!  write(STDOUT, "(2a)") " Using specification file: ", trim(SpecFile)
!  call Spec%Read(SpecFile=trim(SpecFile))

!  write(STDOUT, "(a)") ""
!  write(STDOUT, "(a)") " Processing data ..."
!  call Data%Read(Spec=Spec)
!  if (Spec%PedigreeGiven) then
!    call Data%RecPed%Write(File=trim(Spec%OutputBasename)//trim(Spec%PedigreeFile)//"_Recoded.txt")
!  end if

