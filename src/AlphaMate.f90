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

program AlphaMate
  use Iso_Fortran_Env, STDOUT => output_unit
  use ConstantModule, only : FILELENGTH
  use AlphaHouseMod, only : PrintCpuTime, PrintElapsedTime, PrintDateTime
  use AlphaMateModule

  implicit none

  integer(int32) :: nArg, ClockRate, ClockMax, ClockStartCount, ClockEndCount
  real(real32) :: CpuStartTime, CpuEndTime
  character(len=FILELENGTH) :: SpecFile
  type(AlphaMateSpec) :: Spec
  type(AlphaMateData) :: Data

  write(STDOUT, "(a)") ""
  call AlphaMateTitle
  write(STDOUT, "(a)") ""
  call PrintDateTime
  call Cpu_Time(CpuStartTime)
  call System_Clock(count_rate=ClockRate, count_max=ClockMax)
  call System_Clock(count=ClockStartCount)

  write(STDOUT, "(a)") ""
  write(STDOUT, "(a)") " --- Specifications ---"
  nArg = command_argument_count()
  if (nArg .gt. 0) then
    call Get_Command_Argument(1, SpecFile)
  else
    SpecFile = "AlphaMateSpec.txt"
  end if
  write(STDOUT, "(a)") ""
  write(STDOUT, "(2a)") " Using specification file: ", trim(SpecFile)
  call Spec%Read(SpecFile=SpecFile, LogStdout=.true.)

  write(STDOUT, "(a)") ""
  write(STDOUT, "(a)") " --- System ---"
  call AlphaMateSystem(Spec=Spec, LogStdout=.true.)

  write(STDOUT, "(a)") ""
  write(STDOUT, "(a)") " --- Data ---"
  call Data%Read(Spec=Spec, LogStdout=.true.)

  write(STDOUT, "(a)") ""
  write(STDOUT, "(a)") " --- Optimisation ---"
  call AlphaMateSearch(Spec=Spec, Data=Data, LogStdout=.true.)

  write(STDOUT, "(a)") ""
  write(STDOUT, "(a)") " --- End ---"

  write(STDOUT, "(a)") ""
  call PrintDateTime
  write(STDOUT, "(a)") ""
  call Cpu_Time(CpuEndTime)
  call PrintCpuTime(CpuStartTime, CpuEndTime)
  call System_Clock(count=ClockEndCount)
  call PrintElapsedTime(Start=ClockStartCount, End=ClockEndCount, Rate=ClockRate, Max=ClockMax)
  write(STDOUT, "(a)") ""
  call AlphaMateTitle
  write(STDOUT, "(a)") ""
end program

!###############################################################################
