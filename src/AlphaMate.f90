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
!> @details  Optimize contributions of individuals to the next generation and
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
  use ISO_Fortran_env, STDIN => input_unit, STDOUT => output_unit, STDERR => error_unit
  use ConstantModule, only : FILELENGTH
  use AlphaHouseMod, only : PrintElapsedTime
  use AlphaMateModule

  implicit none

  integer(int32) :: nArg
  real(real32) :: StartTime, EndTime
  character(len=FILELENGTH) :: SpecFile
  type(AlphaMateSpec) :: Spec
  type(AlphaMateData) :: Data

  call cpu_time(StartTime)
  call AlphaMateTitle

  write(STDOUT, "(a)") ""
  write(STDOUT, "(a)") " Specifications ..."
  nArg = command_argument_count()
  if (nArg > 0) then
    call get_command_argument(1, SpecFile)
  else
    SpecFile = "AlphaMateSpec.txt"
  end if
  write(STDOUT, "(2a)") " Using specification file: ", trim(SpecFile)
  call Spec%Read(SpecFile=SpecFile, LogStdout=.true.)

  write(STDOUT, "(a)") ""
  write(STDOUT, "(a)") " Data ..."
  call Data%Read(Spec=Spec, LogStdout=.true.)

  call SetupColNamesAndFormats
  ! call AlphaMateSearch(Spec=Spec, Data=Data, LogStdout=.true.)

  call cpu_time(EndTime)
  call AlphaMateTitle
  call PrintElapsedTime(StartTime, EndTime)
end program

!###############################################################################