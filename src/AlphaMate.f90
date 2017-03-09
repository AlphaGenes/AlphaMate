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

  call ReadSpecAndDataForAlphaMate
  call SetupColNamesAndFormats
  call AlphaMateSearch

  call cpu_time(EndTime)
  call AlphaMateTitle
  call PrintElapsedTime(StartTime, EndTime)
end program

!###############################################################################