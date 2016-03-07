
!###############################################################################

module AlphaSuiteModule

  use ISO_Fortran_Env

  implicit none

  private
  public :: CountLines,Int2Char,RandomOrder,SetSeed,ToLower

  ! List of characters for case conversion in ToLower
  CHARACTER(*),PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
  CHARACTER(*),PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  contains

    !###########################################################################

    subroutine CountLines(FileName,nLines)

      implicit none

      ! Arguments
      character(len=*),intent(in) :: FileName
      integer(int32),intent(out)  :: nLines

      ! Other
      integer(int32) :: f,Unit

      character(len=300) :: DumC

      nLines=0
      open(newunit=Unit,file=trim(FileName),status="old")
      do
        read(Unit,*,iostat=f) DumC
        nLines=nLines+1
        if (f /= 0) then
          nLines=nLines-1
          exit
        end if
      end do
      close(Unit)
    end subroutine

    !###########################################################################

    function Int2Char(i) result(Res)
      ! From http://stackoverflow.com/questions/1262695/converting-integers-to-strings-in-fortran

      implicit none

      integer(int32),intent(in) :: i

      character(:),allocatable :: Res
      character(range(i)+2) :: Tmp

      write(Tmp,'(i0)') i
      Res=trim(Tmp)
    end function

    !###########################################################################

    subroutine RandomOrder(order,n)
      ! Generate a random ordering of the integers 1 ... n.

      implicit none

      ! Arguments
      integer(int32),intent(in)  :: n
      integer(int32),intent(out) :: order(n)

      ! Other
      integer(int32) :: i,j,k

      real(real64) :: wk(n-1)

      do i=1,n
        order(i)=i
      end do

      ! Starting at the end, swap the current last indicator with one
      ! randomly chosen from those preceeding it.
      call random_number(wk)
      do i=n,2,-1
        j=1 + i * wk(i)
        if (j < i) then
          k=order(i)
          order(i)=order(j)
          order(j)=k
        end if
      end do

      return
    end subroutine

    !###########################################################################

    subroutine SetSeed(Seed,SeedFile,Out)

      implicit none

      ! Arguments
      integer(int32),intent(in),optional  :: Seed     ! A number to initialize RNG
      character(len=*),optional           :: SeedFile ! File to save the seed in
      integer(int32),intent(out),optional :: Out      ! Make the seed value available outside

      ! Other
      integer(int32) :: Size,Unit
      integer(int32),allocatable :: SeedList(:)

      ! Get the size of seed array
      call random_seed(size=Size)
      allocate(SeedList(Size))

      ! Set seed
      if (present(Seed)) then ! using the given value
        SeedList(1)=Seed
        SeedList(2:Size)=1
        call random_seed(put=SeedList)
      else                    ! using system/compiler value
        call random_seed
        call random_seed(get=SeedList)
      end if

      ! Save to a file
      if (present(SeedFile)) then
        open(newunit=Unit,file=trim(SeedFile),status="unknown")
        write(Unit,*) SeedList(1)
        close(Unit)
      end if

      ! Output
      if (present(Out)) then
          Out=SeedList(1)
      end if
      deallocate(SeedList)
    end subroutine

    !###########################################################################

    function ToLower(StringIn) result(StringOut)
      ! From https://groups.google.com/forum/#!topic/comp.lang.fortran/CKx1L2Ahkxg
      character(len=*),intent(in) :: StringIn
      character(len(StringIn))    :: StringOut
      integer :: i,n

      ! Copy input string
      StringOut=StringIn

      ! Convert case character by character
      do i=1,len(StringOut)
        n=index(UPPER_CASE,StringOut(i:i))
        if (n /= 0) StringOut(i:i)=LOWER_CASE(n:n)
      end do
    end function

    !###########################################################################
end module

!###############################################################################
