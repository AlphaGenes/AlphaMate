! TODO: PAGE
! TODO: variance within generated crosses - another file with some values per
!       each possible mating or perhaps just matings that one would consider
! TODO: Manual

#ifdef BINARY
#define BINFILE ,form="unformatted"
#else
#define BINFILE
#endif

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#ifdef OS_UNIX
#define DASH "/"
#define COPY "cp"
#define MKDIR "mkdir -p"
#define RMDIR "rm -r"
#define RM "rm"
#define RENAME "mv"
#else
#define DASH "\"
#define COPY "copy"
#define MKDIR "md"
#define RMDIR "rmdir /S"
#define RM "del"
#define RENAME "move /Y"
#endif

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

      real(real64) :: wk

      do i=1,n
        order(i)=i
      end do

      ! Starting at the end, swap the current last indicator with one
      ! randomly chosen from those preceeding it.
      do i=n,2,-1
        call random_number(wk)
        j=1 + i * wk
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
      character(len=*)                    :: SeedFile ! File to save the seed in
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
      open(newunit=Unit,file=trim(SeedFile),status="unknown")
      write(Unit,*) SeedList(1)
      close(Unit)

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

module OrderPack
  ! From http://www.fortran-2000.com/rank/MrgRnk.f90 (2016-02-15)
  ! Modified to work with ISO_Fortran_Env types

  use ISO_Fortran_Env

  implicit none

  private
  public :: MrgRnk

  interface MrgRnk
    module procedure D_MrgRnk,R_MrgRnk,I_MrgRnk
  end interface

  contains

    !###########################################################################

    Subroutine D_MrgRnk (XDONT, IRNGT)
      ! __________________________________________________________
      !   MRGRNK = Merge-sort ranking of an array
      !   For performance reasons, the first 2 passes are taken
      !   out of the standard loop, and use dedicated coding.
      ! __________________________________________________________
      ! __________________________________________________________
        Real(real64), Dimension (:), Intent (In) :: XDONT
        Integer(int32), Dimension (:), Intent (Out) :: IRNGT
      ! __________________________________________________________
        Real(real64) :: XVALA, XVALB
      !
        Integer(int32), Dimension (SIZE(IRNGT)) :: JWRKT
        Integer(int32) :: LMTNA, LMTNC, IRNG1, IRNG2
        Integer(int32) :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
      !
        NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
        Select Case (NVAL)
        Case (:0)
           Return
        Case (1)
           IRNGT (1) = 1
           Return
        Case Default
           Continue
        End Select
      !
      !  Fill-in the index array, creating ordered couples
      !
        Do IIND = 2, NVAL, 2
           If (XDONT(IIND-1) <= XDONT(IIND)) Then
              IRNGT (IIND-1) = IIND - 1
              IRNGT (IIND) = IIND
           Else
              IRNGT (IIND-1) = IIND
              IRNGT (IIND) = IIND - 1
           End If
        End Do
        If (Modulo(NVAL, 2) /= 0) Then
           IRNGT (NVAL) = NVAL
        End If
      !
      !  We will now have ordered subsets A - B - A - B - ...
      !  and merge A and B couples into     C   -   C   - ...
      !
        LMTNA = 2
        LMTNC = 4
      !
      !  First iteration. The length of the ordered subsets goes from 2 to 4
      !
        Do
           If (NVAL <= 2) Exit
      !
      !   Loop on merges of A and B into C
      !
           Do IWRKD = 0, NVAL - 1, 4
              If ((IWRKD+4) > NVAL) Then
                 If ((IWRKD+2) >= NVAL) Exit
      !
      !   1 2 3
      !
                 If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
      !
      !   1 3 2
      !
                 If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                    IRNG2 = IRNGT (IWRKD+2)
                    IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                    IRNGT (IWRKD+3) = IRNG2
      !
      !   3 1 2
      !
                 Else
                    IRNG1 = IRNGT (IWRKD+1)
                    IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                    IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                    IRNGT (IWRKD+2) = IRNG1
                 End If
                 Exit
              End If
      !
      !   1 2 3 4
      !
              If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
      !
      !   1 3 x x
      !
              If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                 IRNG2 = IRNGT (IWRKD+2)
                 IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                 If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
      !   1 3 2 4
                    IRNGT (IWRKD+3) = IRNG2
                 Else
      !   1 3 4 2
                    IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                    IRNGT (IWRKD+4) = IRNG2
                 End If
      !
      !   3 x x x
      !
              Else
                 IRNG1 = IRNGT (IWRKD+1)
                 IRNG2 = IRNGT (IWRKD+2)
                 IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                 If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                    IRNGT (IWRKD+2) = IRNG1
                    If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
      !   3 1 2 4
                       IRNGT (IWRKD+3) = IRNG2
                    Else
      !   3 1 4 2
                       IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                       IRNGT (IWRKD+4) = IRNG2
                    End If
                 Else
      !   3 4 1 2
                    IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                    IRNGT (IWRKD+3) = IRNG1
                    IRNGT (IWRKD+4) = IRNG2
                 End If
              End If
           End Do
      !
      !  The Cs become As and Bs
      !
           LMTNA = 4
           Exit
        End Do
      !
      !  Iteration loop. Each time, the length of the ordered subsets
      !  is doubled.
      !
        Do
           If (LMTNA >= NVAL) Exit
           IWRKF = 0
           LMTNC = 2 * LMTNC
      !
      !   Loop on merges of A and B into C
      !
           Do
              IWRK = IWRKF
              IWRKD = IWRKF + 1
              JINDA = IWRKF + LMTNA
              IWRKF = IWRKF + LMTNC
              If (IWRKF >= NVAL) Then
                 If (JINDA >= NVAL) Exit
                 IWRKF = NVAL
              End If
              IINDA = 1
              IINDB = JINDA + 1
      !
      !   Shortcut for the case when the max of A is smaller
      !   than the min of B. This line may be activated when the
      !   initial set is already close to sorted.
      !
      !          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
      !
      !  One steps in the C subset, that we build in the final rank array
      !
      !  Make a copy of the rank array for the merge iteration
      !
              JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
      !
              XVALA = XDONT (JWRKT(IINDA))
              XVALB = XDONT (IRNGT(IINDB))
      !
              Do
                 IWRK = IWRK + 1
      !
      !  We still have unprocessed values in both A and B
      !
                 If (XVALA > XVALB) Then
                    IRNGT (IWRK) = IRNGT (IINDB)
                    IINDB = IINDB + 1
                    If (IINDB > IWRKF) Then
      !  Only A still with unprocessed values
                       IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                       Exit
                    End If
                    XVALB = XDONT (IRNGT(IINDB))
                 Else
                    IRNGT (IWRK) = JWRKT (IINDA)
                    IINDA = IINDA + 1
                    If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                    XVALA = XDONT (JWRKT(IINDA))
                 End If
      !
              End Do
           End Do
      !
      !  The Cs become As and Bs
      !
           LMTNA = 2 * LMTNA
        End Do
      !
        Return
      !
    End Subroutine

    !###########################################################################

    Subroutine R_MrgRnk (XDONT, IRNGT)
      ! __________________________________________________________
      !   MRGRNK = Merge-sort ranking of an array
      !   For performance reasons, the first 2 passes are taken
      !   out of the standard loop, and use dedicated coding.
      ! __________________________________________________________
      ! _________________________________________________________
        Real(real32), Dimension (:), Intent (In) :: XDONT
        Integer(int32), Dimension (:), Intent (Out) :: IRNGT
      ! __________________________________________________________
        Real(real32) :: XVALA, XVALB
      !
        Integer(int32), Dimension (SIZE(IRNGT)) :: JWRKT
        Integer(int32) :: LMTNA, LMTNC, IRNG1, IRNG2
        Integer(int32) :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
      !
        NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
        Select Case (NVAL)
        Case (:0)
           Return
        Case (1)
           IRNGT (1) = 1
           Return
        Case Default
           Continue
        End Select
      !
      !  Fill-in the index array, creating ordered couples
      !
        Do IIND = 2, NVAL, 2
           If (XDONT(IIND-1) <= XDONT(IIND)) Then
              IRNGT (IIND-1) = IIND - 1
              IRNGT (IIND) = IIND
           Else
              IRNGT (IIND-1) = IIND
              IRNGT (IIND) = IIND - 1
           End If
        End Do
        If (Modulo(NVAL, 2) /= 0) Then
           IRNGT (NVAL) = NVAL
        End If
      !
      !  We will now have ordered subsets A - B - A - B - ...
      !  and merge A and B couples into     C   -   C   - ...
      !
        LMTNA = 2
        LMTNC = 4
      !
      !  First iteration. The length of the ordered subsets goes from 2 to 4
      !
        Do
           If (NVAL <= 2) Exit
      !
      !   Loop on merges of A and B into C
      !
           Do IWRKD = 0, NVAL - 1, 4
              If ((IWRKD+4) > NVAL) Then
                 If ((IWRKD+2) >= NVAL) Exit
      !
      !   1 2 3
      !
                 If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
      !
      !   1 3 2
      !
                 If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                    IRNG2 = IRNGT (IWRKD+2)
                    IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                    IRNGT (IWRKD+3) = IRNG2
      !
      !   3 1 2
      !
                 Else
                    IRNG1 = IRNGT (IWRKD+1)
                    IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                    IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                    IRNGT (IWRKD+2) = IRNG1
                 End If
                 Exit
              End If
      !
      !   1 2 3 4
      !
              If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
      !
      !   1 3 x x
      !
              If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                 IRNG2 = IRNGT (IWRKD+2)
                 IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                 If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
      !   1 3 2 4
                    IRNGT (IWRKD+3) = IRNG2
                 Else
      !   1 3 4 2
                    IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                    IRNGT (IWRKD+4) = IRNG2
                 End If
      !
      !   3 x x x
      !
              Else
                 IRNG1 = IRNGT (IWRKD+1)
                 IRNG2 = IRNGT (IWRKD+2)
                 IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                 If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                    IRNGT (IWRKD+2) = IRNG1
                    If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
      !   3 1 2 4
                       IRNGT (IWRKD+3) = IRNG2
                    Else
      !   3 1 4 2
                       IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                       IRNGT (IWRKD+4) = IRNG2
                    End If
                 Else
      !   3 4 1 2
                    IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                    IRNGT (IWRKD+3) = IRNG1
                    IRNGT (IWRKD+4) = IRNG2
                 End If
              End If
           End Do
      !
      !  The Cs become As and Bs
      !
           LMTNA = 4
           Exit
        End Do
      !
      !  Iteration loop. Each time, the length of the ordered subsets
      !  is doubled.
      !
        Do
           If (LMTNA >= NVAL) Exit
           IWRKF = 0
           LMTNC = 2 * LMTNC
      !
      !   Loop on merges of A and B into C
      !
           Do
              IWRK = IWRKF
              IWRKD = IWRKF + 1
              JINDA = IWRKF + LMTNA
              IWRKF = IWRKF + LMTNC
              If (IWRKF >= NVAL) Then
                 If (JINDA >= NVAL) Exit
                 IWRKF = NVAL
              End If
              IINDA = 1
              IINDB = JINDA + 1
      !
      !   Shortcut for the case when the max of A is smaller
      !   than the min of B. This line may be activated when the
      !   initial set is already close to sorted.
      !
      !          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
      !
      !  One steps in the C subset, that we build in the final rank array
      !
      !  Make a copy of the rank array for the merge iteration
      !
              JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
      !
              XVALA = XDONT (JWRKT(IINDA))
              XVALB = XDONT (IRNGT(IINDB))
      !
              Do
                 IWRK = IWRK + 1
      !
      !  We still have unprocessed values in both A and B
      !
                 If (XVALA > XVALB) Then
                    IRNGT (IWRK) = IRNGT (IINDB)
                    IINDB = IINDB + 1
                    If (IINDB > IWRKF) Then
      !  Only A still with unprocessed values
                       IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                       Exit
                    End If
                    XVALB = XDONT (IRNGT(IINDB))
                 Else
                    IRNGT (IWRK) = JWRKT (IINDA)
                    IINDA = IINDA + 1
                    If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                    XVALA = XDONT (JWRKT(IINDA))
                 End If
      !
              End Do
           End Do
      !
      !  The Cs become As and Bs
      !
           LMTNA = 2 * LMTNA
        End Do
      !
        Return
      !
    End Subroutine

    !###########################################################################

    Subroutine I_MrgRnk (XDONT, IRNGT)
      ! __________________________________________________________
      !   MRGRNK = Merge-sort ranking of an array
      !   For performance reasons, the first 2 passes are taken
      !   out of the standard loop, and use dedicated coding.
      ! __________________________________________________________
      ! __________________________________________________________
        Integer(int32), Dimension (:), Intent (In)  :: XDONT
        Integer(int32), Dimension (:), Intent (Out) :: IRNGT
      ! __________________________________________________________
        Integer(int32) :: XVALA, XVALB
      !
        Integer(int32), Dimension (SIZE(IRNGT)) :: JWRKT
        Integer(int32) :: LMTNA, LMTNC, IRNG1, IRNG2
        Integer(int32) :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
      !
        NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
        Select Case (NVAL)
        Case (:0)
           Return
        Case (1)
           IRNGT (1) = 1
           Return
        Case Default
           Continue
        End Select
      !
      !  Fill-in the index array, creating ordered couples
      !
        Do IIND = 2, NVAL, 2
           If (XDONT(IIND-1) <= XDONT(IIND)) Then
              IRNGT (IIND-1) = IIND - 1
              IRNGT (IIND) = IIND
           Else
              IRNGT (IIND-1) = IIND
              IRNGT (IIND) = IIND - 1
           End If
        End Do
        If (Modulo(NVAL, 2) /= 0) Then
           IRNGT (NVAL) = NVAL
        End If
      !
      !  We will now have ordered subsets A - B - A - B - ...
      !  and merge A and B couples into     C   -   C   - ...
      !
        LMTNA = 2
        LMTNC = 4
      !
      !  First iteration. The length of the ordered subsets goes from 2 to 4
      !
        Do
           If (NVAL <= 2) Exit
      !
      !   Loop on merges of A and B into C
      !
           Do IWRKD = 0, NVAL - 1, 4
              If ((IWRKD+4) > NVAL) Then
                 If ((IWRKD+2) >= NVAL) Exit
      !
      !   1 2 3
      !
                 If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
      !
      !   1 3 2
      !
                 If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                    IRNG2 = IRNGT (IWRKD+2)
                    IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                    IRNGT (IWRKD+3) = IRNG2
      !
      !   3 1 2
      !
                 Else
                    IRNG1 = IRNGT (IWRKD+1)
                    IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                    IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                    IRNGT (IWRKD+2) = IRNG1
                 End If
                 Exit
              End If
      !
      !   1 2 3 4
      !
              If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
      !
      !   1 3 x x
      !
              If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                 IRNG2 = IRNGT (IWRKD+2)
                 IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                 If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
      !   1 3 2 4
                    IRNGT (IWRKD+3) = IRNG2
                 Else
      !   1 3 4 2
                    IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                    IRNGT (IWRKD+4) = IRNG2
                 End If
      !
      !   3 x x x
      !
              Else
                 IRNG1 = IRNGT (IWRKD+1)
                 IRNG2 = IRNGT (IWRKD+2)
                 IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                 If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                    IRNGT (IWRKD+2) = IRNG1
                    If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
      !   3 1 2 4
                       IRNGT (IWRKD+3) = IRNG2
                    Else
      !   3 1 4 2
                       IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                       IRNGT (IWRKD+4) = IRNG2
                    End If
                 Else
      !   3 4 1 2
                    IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                    IRNGT (IWRKD+3) = IRNG1
                    IRNGT (IWRKD+4) = IRNG2
                 End If
              End If
           End Do
      !
      !  The Cs become As and Bs
      !
           LMTNA = 4
           Exit
        End Do
      !
      !  Iteration loop. Each time, the length of the ordered subsets
      !  is doubled.
      !
        Do
           If (LMTNA >= NVAL) Exit
           IWRKF = 0
           LMTNC = 2 * LMTNC
      !
      !   Loop on merges of A and B into C
      !
           Do
              IWRK = IWRKF
              IWRKD = IWRKF + 1
              JINDA = IWRKF + LMTNA
              IWRKF = IWRKF + LMTNC
              If (IWRKF >= NVAL) Then
                 If (JINDA >= NVAL) Exit
                 IWRKF = NVAL
              End If
              IINDA = 1
              IINDB = JINDA + 1
      !
      !   Shortcut for the case when the max of A is smaller
      !   than the min of B. This line may be activated when the
      !   initial set is already close to sorted.
      !
      !          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
      !
      !  One steps in the C subset, that we build in the final rank array
      !
      !  Make a copy of the rank array for the merge iteration
      !
              JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
      !
              XVALA = XDONT (JWRKT(IINDA))
              XVALB = XDONT (IRNGT(IINDB))
      !
              Do
                 IWRK = IWRK + 1
      !
      !  We still have unprocessed values in both A and B
      !
                 If (XVALA > XVALB) Then
                    IRNGT (IWRK) = IRNGT (IINDB)
                    IINDB = IINDB + 1
                    If (IINDB > IWRKF) Then
      !  Only A still with unprocessed values
                       IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                       Exit
                    End If
                    XVALB = XDONT (IRNGT(IINDB))
                 Else
                    IRNGT (IWRK) = JWRKT (IINDA)
                    IINDA = IINDA + 1
                    If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                    XVALA = XDONT (JWRKT(IINDA))
                 End If
      !
              End Do
           End Do
      !
      !  The Cs become As and Bs
      !
           LMTNA = 2 * LMTNA
        End Do
      !
        Return
      !
    End Subroutine

    !###########################################################################
end module

!###############################################################################

module AlphaEvolAlg

  use ISO_Fortran_Env, STDIN=>input_unit,STDOUT=>output_unit,STDERR=>error_unit

  implicit none

  type :: EvolAlgCrit
    ! TODO: How do I get this done generically? I have created an extended type
    !       in AlphaMateModule of base EvolAlgCrit here, but got bitten through
    !       using polymorphic (class() stuff) objects, where I am not allowed to
    !       do a=b etc.
    real(real64) :: Value
    real(real64) :: Penalty
    real(real64) :: Gain
    real(real64) :: GainStand
    real(real64) :: PopInb
    real(real64) :: RatePopInb
    real(real64) :: PopInb2
    real(real64) :: IndInb
  end type

  private
  public :: EvolAlgDE,EvolAlgCrit

  contains

    !###########################################################################

    subroutine EvolAlgDE(nParam,nSol,Init,nGen,nGenBurnIn,nGenStop,&
      StopTolerance,nGenPrint,File,CritType,CRBurnIn,CRLate,FBase,FHigh1,FHigh2,&
      CalcCriterion,LogHeader,Log,BestCriterion)

      implicit none

      ! Arguments
      integer(int32),intent(in)          :: nParam         ! No. of parameters in a solution
      integer(int32),intent(in)          :: nSol           ! No. of solutions to test each generation
      real(real64),intent(in),optional   :: Init(:,:)      ! Initial solutions to start with
      integer(int32),intent(in)          :: nGen           ! No. of generations to run
      integer(int32),intent(in)          :: nGenBurnIn     ! No. of generations with more
      integer(int32),intent(in)          :: nGenStop       ! Stop after no progress for nGenerationStop
      real(real64),intent(in)            :: StopTolerance  ! Stopping tolerance
      integer(int32),intent(in)          :: nGenPrint      ! Print changed solution every nGenerationPrint
      character(len=*),intent(in)        :: File           ! Which file to write to
      character(len=*),intent(in)        :: CritType       ! Passed to CalcCriterion
      real(real64),intent(in),optional   :: CRBurnIn       ! Crossover rate for nGenBurnIn
      real(real64),intent(in),optional   :: CRLate         ! Crossover rate
      real(real64),intent(in),optional   :: FBase          ! F is multiplier of difference used to mutate
      real(real64),intent(in),optional   :: FHigh1         ! F is multiplier of difference used to mutate
      real(real64),intent(in),optional   :: FHigh2         ! F is multiplier of difference used to mutate
      type(EvolAlgCrit),intent(out)      :: BestCriterion  ! Criterion for the best solution found

      interface
        subroutine CalcCriterion(Sol,CritType,Criterion)
          use ISO_Fortran_Env
          import :: EvolAlgCrit
          real(real64),intent(inout)     :: Sol(:)
          character(len=*),intent(in)    :: CritType
          type(EvolAlgCrit),intent(out) :: Criterion
        end subroutine

        subroutine LogHeader(LogUnit)
          use ISO_Fortran_Env
          integer(int32),intent(in) :: LogUnit
        end subroutine

        subroutine Log(LogUnit,Gen,AcceptRate,Criterion)
          use ISO_Fortran_Env
          import :: EvolAlgCrit
          integer(int32),intent(in)    :: LogUnit
          integer(int32),intent(in)    :: Gen
          real(real64),intent(in)      :: AcceptRate
          type(EvolAlgCrit),intent(in) :: Criterion
        end subroutine
      end interface

      ! Other
      integer(int32) :: nInit,Param,ParamLoc,Sol,Gen,LastGenPrint
      integer(int32) :: SolA,SolB,SolC,BestSol
      integer(int32) :: Unit

      real(real64) :: RanNum,FInt,FBaseInt,FHigh1Int,FHigh2Int,CRInt,CRBurnInInt,CRLateInt
      real(real64) :: AcceptRate,OldChrom(nParam,nSol),NewChrom(nParam,nSol),Chrom(nParam)

      logical :: DiffOnly,BestSolChanged

      type(EvolAlgCrit) :: Criterion(nSol),CriterionHold,BestCriterionStop

      ! --- Initialize ---

      LastGenPrint=0
      BestCriterion%Value=-huge(RanNum)
      BestCriterionStop%Value=-huge(RanNum)

      ! --- Printout ---

      open(newunit=Unit,file=trim(File),status="unknown")
      call LogHeader(LogUnit=Unit)

      ! --- Set parameters ---

      ! Crossover rate
      ! ... for later climbs
      if (present(CRLate)) then
        CRLateInt=CRLate
      else
        CRLateInt=0.1d0
      end if
      ! ... for first few generations (burn-in)
      if (present(CRBurnIn)) then
        CRBurnInInt=CRBurnIn
      else
        CRBurnInInt=2.0d0*CRLateInt
      end if

      ! F is multiplier of difference used to mutate
      ! Typically between 0.2 and 2.0
      ! (if alleles should be integer, keep F as integer)
      ! ... conservative moves
      if (present(FBase)) then
        FBaseInt=FBase
      else
        FBaseInt=0.1d0
      end if
      ! ... adventurous moves
      if (present(FHigh1)) then
        FHigh1Int=FHigh1
      else
        FHigh1Int=10.0d0*FBaseInt
      end if
      if (present(FHigh2)) then
        FHigh2Int=FHigh2
      else
        FHigh2Int=4.0d0*FHigh1Int
      end if

      ! --- Initialise foundation population of solutions ---

      if (present(Init)) then
        nInit=size(Init(:,:),2)
        do Sol=1,nInit
          OldChrom(:,Sol)=Init(:,Sol)
          call CalcCriterion(OldChrom(:,Sol),CritType,Criterion(Sol))
        end do
      else
        nInit=1
      end if
      do Sol=nInit,nSol
        call random_number(OldChrom(:,Sol))
        call CalcCriterion(OldChrom(:,Sol),CritType,Criterion(Sol))
      end do

      ! --- Evolve ---

      do Gen=1,nGen

        ! Vary differential and non-differential mutation to escape valleys
        if (mod(Gen,3) == 0) then
          DiffOnly=.true.
        else
          DiffOnly=.false.
        end if

        ! Burn-in
        if (Gen < nGenBurnIn) then
          CRInt=CRBurnInInt
        else
          CRInt=CRLateInt
        end if

        ! Vary mutation rate every few generations
        if (mod(Gen,4) == 0) then
          FInt=FHigh1Int
        else
          FInt=FBaseInt
        end if

        if (mod(Gen,7) == 0) then
          FInt=FHigh2Int
        else
          FInt=FBaseInt
        end if

        ! --- Generate competitors ---

        ! TODO: Paralelize this loop? Is it worth it?
        !       Can we do it given the saving to OldChrom&NewChrom?
        BestSolChanged=.false.
        AcceptRate=0.0d0
        do Sol=1,nSol

          ! --- Mutate and recombine ---

          ! Get three different solutions
          SolA=Sol
          do while (SolA == Sol)
            call random_number(RanNum)
            SolA=int(RanNum*nSol)+1
          end do
          SolB=Sol
          do while ((SolB == Sol) .or. (SolB == SolA))
            call random_number(RanNum)
            SolB=int(RanNum*nSol)+1
          end do
          SolC=Sol
          do while ((SolC == Sol) .or. (SolC == SolA) .or. (SolC == SolB))
            call random_number(RanNum)
            SolC=int(RanNum*nSol)+1
          end do

          ! Mate the solutions
          call random_number(RanNum)
          Param=int(RanNum*nParam)+1 ! Cycle through parameters starting at a random point
          do ParamLoc=1,nParam
            call random_number(RanNum)
            if ((RanNum < CRInt) .or. (ParamLoc == nParam)) then
              ! Recombine
              call random_number(RanNum)
              if ((RanNum < 0.8d0) .or. DiffOnly) then
                ! Differential mutation (with prob 0.8 or 1)
                Chrom(Param)=OldChrom(Param,SolC) + FInt*(OldChrom(Param,SolA)-OldChrom(Param,SolB))
              else
                ! Non-differential mutation (to avoid getting stuck)
                call random_number(RanNum)
                if (RanNum < 0.5d0) then
                  call random_number(RanNum)
                  Chrom(Param)=OldChrom(Param,SolC) * (0.9d0+0.2d0*RanNum)
                else
                  call random_number(RanNum)
                  Chrom(Param)=OldChrom(Param,SolC) + 0.01d0*FInt*(OldChrom(Param,SolA)+0.01d0)*(RanNum-0.5d0)
                end if
              end if
            else
              ! Do not recombine
              Chrom(Param)=OldChrom(Param,Sol)
            end if
            Param=Param+1
            if (Param > nParam) then
              Param=Param-nParam
            end if
          end do

          ! --- Evaluate and Select ---

          call CalcCriterion(Chrom,CritType,CriterionHold)      ! Merit of competitor
          if (CriterionHold%Value >= Criterion(Sol)%Value) then ! If competitor is better or equal, keep it
            NewChrom(:,Sol)=Chrom(:)                            !   ("equal" to force evolution)
            Criterion(Sol)=CriterionHold
            AcceptRate=AcceptRate+1.0d0
          else
            NewChrom(:,Sol)=OldChrom(:,Sol)                     ! Else keep the old solution
          end if
        end do ! Sol

        AcceptRate=AcceptRate/dble(nSol)

        ! --- New parents ---

        OldChrom(:,:)=NewChrom(:,:)

        ! --- Find the best solution in this generation ---

        BestSol=maxloc(Criterion(:)%Value,dim=1)
        if (Criterion(BestSol)%Value > BestCriterion%Value) then
          BestSolChanged=.true.
          BestCriterion=Criterion(BestSol)
        end if

        ! --- Monitor ---

        if (BestSolChanged) then
          if ((Gen == 1) .or. ((Gen - LastGenPrint) >= nGenPrint)) then
            LastGenPrint=Gen
            call Log(Unit,Gen,AcceptRate,BestCriterion)
          end if
        end if

        ! --- Test if solution is improving to stop early ---

        if (mod(Gen,nGenStop) == 0) then
          if ((BestCriterion%Value-BestCriterionStop%Value) > StopTolerance) then
            BestCriterionStop=BestCriterion
          else
            write(STDOUT,"(a,f9.5,a,i9,a)") "NOTE: Evolutionary algorithm did not improve objective for ",StopTolerance, " in the last ",nGenStop," generations. Stopping."
            write(STDOUT,"(a)") " "
            exit
          end if
        end if

      end do ! Gen

      ! --- The winner ---

      call Log(Unit,Gen,AcceptRate,BestCriterion)
      write(STDOUT,"(a)") " "
      close(Unit)
    end subroutine

    !###########################################################################
end module

!###############################################################################

module AlphaMateModule

  use ISO_Fortran_Env, STDIN=>input_unit,STDOUT=>output_unit,STDERR=>error_unit
  use IFPort,only : SystemQQ
  use OrderPack,only : MrgRnk
  use AlphaEvolAlg,only : EvolAlgDE,EvolAlgCrit
  use AlphaSuiteModule,only : CountLines,Int2Char,RandomOrder,SetSeed,ToLower

  implicit none

  integer(int32) :: nInd,nMat,nPar,nPotPar1,nPotPar2,nMal,nFem,nPar1,nPar2,nFrontierSteps
  integer(int32) :: EvolAlgNSol,EvolAlgNGen,EvolAlgNGenBurnIn,EvolAlgNGenStop,EvolAlgNGenPrint
  integer(int32),allocatable :: Gender(:),IdPotPar1(:),IdPotPar2(:),nVecPar1(:),nVecPar2(:),nVec(:),Mate(:,:)

  real(real64) :: LimitPar1Min,LimitPar1Max,LimitPar2Min,LimitPar2Max
  real(real64) :: EvolAlgStopTol,EvolAlgCRBurnIn,EvolAlgCRLate,EvolAlgFBase,EvolAlgFHigh1,EvolAlgFHigh2
  real(real64) :: PopInbOld,PopInbTarget,RatePopInbTarget,GainMinStand
  real(real64) :: PopInbPenalty,IndInbPenalty,SelfingPenalty,LimitPar1Penalty,LimitPar2Penalty
  real(real64),allocatable :: Bv(:),BvStand(:),RelMtx(:,:),RatePopInbFrontier(:),xVec(:)

  character(len=300),allocatable :: IdC(:)

  logical :: MinThenOpt,GenderMatters,EqualizePar1,EqualizePar2
  logical :: SelfingAllowed,InferPopInbOld,EvaluateFrontier

  private
  public :: AlphaMateTitles,ReadSpecAndDataForAlphaMate,SetInbreedingParameters
  public :: AlphaMateSearch,EvolAlgLogHeaderForAlphaMate,EvolAlgLogForAlphaMate
  public :: FixSolMateAndCalcCrit

  contains

    !###########################################################################

    subroutine AlphaMateTitles
      write(STDOUT,"(a)") ""
      write(STDOUT,"(a30,a,a30)") " ","**********************"," "
      write(STDOUT,"(a30,a,a30)") " ","*                    *"," "
      write(STDOUT,"(a30,a,a30)") " ","*     AlphaMate      *"," "
      write(STDOUT,"(a30,a,a30)") " ","*                    *"," "
      write(STDOUT,"(a30,a,a30)") " ","**********************"
      write(STDOUT,"(a30,a,a30)") " ","VERSION:"//TOSTRING(VERS)," "
      write(STDOUT,"(a15,a)")     " ","Software for optimizing contributions to the next generation"
      write(STDOUT,"(a)") ""
      write(STDOUT,"(a35,a)")     " ","No Liability"
      write(STDOUT,"(a25,a)")     " ","Bugs to John.Hickey@roslin.ed.ac.uk"
      write(STDOUT,"(a)") ""
    end subroutine

    !###########################################################################

    subroutine ReadSpecAndDataForAlphaMate

      implicit none

      integer(int32) :: i,j,DumI,jMal,jFem,nIndTmp,GenderTmp,Seed
      integer(int32) :: UnitSpec,UnitRelMtx,UnitBv,UnitGender
      integer(int32),allocatable :: Order(:)

      real(real64) :: BvTmp

      logical :: Success

      character(len=1000) :: DumC,IdCTmp
      character(len=1000) :: RelMtxFile,BvFile,GenderFile,SeedFile

      ! --- Spec file ---

      Success=SystemQQ(COPY//" AlphaMateSpec.txt AlphaMateResults"//DASH//"AlphaMateSpec.txt")
      if (.not.Success) then
        write(STDERR,"(a)") "ERROR: Failed to copy the AlphaMateSpec.txt file in the output folder!"
        write(STDERR,"(a)") " "
        stop 1
      end if

      open(newunit=UnitSpec,file="AlphaMateSpec.txt",status="old")
      write(STDOUT,"(a)") "SpecFile: AlphaMateSpec.txt"

      ! Mode
      read(UnitSpec,*) DumC,DumC
      if (ToLower(trim(DumC)) == "minthenopt") then
        MinThenOpt=.true.
        write(STDOUT,"(a,l)") "Mode: MinThenOpt"
      else
        MinThenOpt=.false.
        write(STDOUT,"(a,l)") "Mode: Opt"
      end if

      ! RelationshipMatrixFile
      read(UnitSpec,*) DumC,RelMtxFile
      write(STDOUT,"(a,a)") "RelationshipMatrixFile: ",trim(RelMtxFile)

      ! BreedingValueFile
      read(UnitSpec,*) DumC,BvFile
      write(STDOUT,"(a,a)") "BreedingValueFile: ",trim(BvFile)

      call CountLines(RelMtxFile,nInd)
      call CountLines(BvFile,nIndTmp)

      if (nIndTmp /= nInd) then
        write(STDERR,"(a)") "ERROR: Number of individuals in ebv file and Relationship Matrix file is not the same!"
        write(STDERR,"(a)") " "
        stop 1
      end if

      ! GenderFile
      read(UnitSpec,*) DumC,GenderFile
      if (ToLower(trim(GenderFile)) /= "none") then
        GenderMatters=.true.
        write(STDOUT,"(a,a)") "GenderFile: ",trim(GenderFile)
      else
        GenderMatters=.false.
      end if

      ! NumberOfIndividuals
      read(UnitSpec,*) DumC,nInd
      write(STDOUT,"(a,i9)") "NumberOfIndividuals: ",nInd

      ! NumberOfMatings
      read(UnitSpec,*) DumC,nMat
      write(STDOUT,"(a,i9)") "NumberOfMatings: ",nMat
      ! TODO: In animals one would not be able to generate more matings than there is individuals, i.e.,
      !       10 males and 10 females can give 10 matings only (females are a bottleneck). But if we do
      !       plants or collect ova from females, then we could technically do 10*10=100 matings (each male
      !       with each female).
      if (nMat > nInd) then
        write(STDOUT,"(a)") "NOTE: The number of matings is larger than the number of all individuals! Was this really the intention?"
        write(STDOUT,"(a,i9)") "NOTE: Number of     matings: ",nMat
        write(STDOUT,"(a,i9)") "NOTE: Number of individuals: ",nInd
        write(STDOUT,"(a)") " "
      end if

      ! NumberOfParents
      read(UnitSpec,*) DumC,nPar
      write(STDOUT,"(a,i9)") "NumberOfParents: ",nPar
      if (nPar > nInd) then
        write(STDERR,"(a)") "ERROR: The number of parents can not be larger than the number of individuals!"
        write(STDERR,"(a,i9)") "ERROR: Number of     parents: ",nPar
        write(STDERR,"(a,i9)") "ERROR: Number of individuals: ",nInd
        write(STDERR,"(a)") " "
        stop 1
      end if
      if (nMat > nPar) then
        write(STDOUT,"(a)") "NOTE: The number of matings is larger than the number of parents! Was this really the intention?"
        write(STDOUT,"(a,i9)") "NOTE: Number of matings: ",nMat
        write(STDOUT,"(a,i9)") "NOTE: Number of parents: ",nPar
        write(STDOUT,"(a)") " "
      end if

      ! NumberOfMaleParents
      read(UnitSpec,*) DumC,nPar1

      ! NumberOfFemaleParents
      read(UnitSpec,*) DumC,nPar2

      if (.not.GenderMatters) then
        nPar1=nPar
      else
        write(STDOUT,"(a,i9)") "NumberOfMaleParents: ",nPar1
        write(STDOUT,"(a,i9)") "NumberOfFemaleParents: ",nPar2
      end if
      if (GenderMatters) then
        if ((nPar1+nPar2) > nInd) then
          write(STDERR,"(a)") "ERROR: The number of parents can not be larger than the number of individuals!"
          write(STDERR,"(a,i9)") "ERROR: Number of        parents: ",nPar1+nPar2
          write(STDERR,"(a,i9)") "ERROR: Number of   male parents: ",nPar1
          write(STDERR,"(a,i9)") "ERROR: Number of female parents: ",nPar2
          write(STDERR,"(a,i9)") "ERROR: Number of    individuals: ",nInd
          write(STDERR,"(a)") " "
          stop 1
        end if
        if ((nPar1+nPar2) /= nPar) then
          write(STDOUT,"(a)") "NOTE: The number of male and female parents does not match with the total number of parents - redefined!"
          write(STDOUT,"(a,i9)")   "NOTE: Number of   male parents: ",nPar1
          write(STDOUT,"(a,i9)")   "NOTE: Number of female parents: ",nPar2
          write(STDOUT,"(a,i9,a)") "NOTE: Number of        parents: ",nPar," (defined)"
          nPar=nPar1+nPar2
          write(STDOUT,"(a,i9,a)") "NOTE: Number of        parents: ",nPar," (redefined)"
          write(STDOUT,"(a)") " "
        end if
        if ((nMat > nPar1) .and. (nMat > nPar2)) then
          write(STDOUT,"(a)") "NOTE: The number of matings is larger than the number of male and female parents! Was this really the intention?"
          write(STDOUT,"(a,i9)") "NOTE: Number of        matings: ",nMat
          write(STDOUT,"(a,i9)") "NOTE: Number of   male parents: ",nPar1
          write(STDOUT,"(a,i9)") "NOTE: Number of female parents: ",nPar2
          write(STDOUT,"(a)") " "
        end if
      end if

      ! EqualizeParentContributions
      read(UnitSpec,*) DumC,DumC
      if (ToLower(trim(DumC)) == "yes") then
        EqualizePar1=.true.
        write(STDOUT,"(a)") "EqualizeParentContributions: yes"
      else
        EqualizePar1=.false.
        write(STDOUT,"(a)") "EqualizeParentContributions: no"
      end if

      ! EqualizeMaleParentContributions
      read(UnitSpec,*) DumC,DumC
      if (ToLower(trim(DumC)) == "yes") then
        if (GenderMatters) then
          EqualizePar1=.true.
          if (mod(nMat,nPar1) /= 0) then
            ! TODO: might consider handling this better at some point
            write(STDERR,"(a)") "ERROR: When contributions are equalized the number of matings needs to divide into the number of parents"
            write(STDERR,"(a,i9)") "ERROR: Number of       matings: ",nMat
            write(STDERR,"(a,i9)") "ERROR: Number of  male parents: ",nPar1
            write(STDERR,"(a,i9)") "ERROR: Modulo (should be zero): ",mod(nMat,nPar1)
            write(STDERR,"(a)") " "
            stop 1
          end if
          write(STDOUT,"(a)") "EqualizeMaleParentContributions: yes"
        else
          EqualizePar1=.false.
          write(STDOUT,"(a)") "EqualizeMaleParentContributions: no"
          write(STDOUT,"(a)") "NOTE: Equalize male parent contributions option ignored when the gender file is not given!"
          write(STDOUT,"(a)") " "
        end if
      else
        if (GenderMatters) then
          EqualizePar1=.false.
          write(STDOUT,"(a)") "EqualizeMaleParentContributions: no"
        end if
      end if

      ! EqualizeFemaleParentContributions
      read(UnitSpec,*) DumC,DumC
      if (ToLower(trim(DumC)) == "yes") then
        if (GenderMatters) then
          EqualizePar2=.true.
          if (mod(nMat,nPar2) /= 0) then
            ! TODO: might consider handling this better at some point
            write(STDERR,"(a)") "ERROR: When contributions are equalized the number of matings needs to divide into the number of parents"
            write(STDERR,"(a,i9)") "ERROR: Number of        matings: ",nMat
            write(STDERR,"(a,i9)") "ERROR: Number of female parents: ",nPar2
            write(STDERR,"(a,i9)") "ERROR: Modulo  (should be zero): ",mod(nMat,nPar2)
            write(STDERR,"(a)") " "
            stop 1
          end if
          write(STDOUT,"(a)") "EqualizeFemaleParentContributions: yes"
        else
          EqualizePar2=.false.
          write(STDOUT,"(a)") "EqualizeFemaleParentContributions: no"
          write(STDOUT,"(a)") "NOTE: Equalize female parent contributions option ignored when the gender file is not given!"
          write(STDOUT,"(a)") " "
        end if
      else
        if (GenderMatters) then
          EqualizePar2=.false.
          write(STDOUT,"(a)") "EqualizeFemaleParentContributions: no"
        end if
      end if

      ! LimitParentContributions
      read(UnitSpec,*) DumC,DumC
      LimitPar1Min=1.0d0
      LimitPar1Max=huge(LimitPar1Max)-1.0d0
      LimitPar1Penalty=0.0d0
      if (ToLower(trim(DumC)) == "yes") then
        if (EqualizePar1) then
          write(STDOUT,"(a)") "LimitParentContributions: no"
          write(STDOUT,"(a)") "NOTE: Limit parent contributions option ignored when contributions are to be equalized!"
          write(STDOUT,"(a)") " "
        else
          backspace(UnitSpec)
          read(UnitSpec,*) DumC,DumC,LimitPar1Min,LimitPar1Max,LimitPar1Penalty
          write(STDOUT,"(a,i9,a,i9,a,f9.5)") "LimitParentContributions: yes, min ",&
            nint(LimitPar1Min),", max ",nint(LimitPar1Max),", penalty ",LimitPar1Penalty
        end if
      else
        write(STDOUT,"(a)") "LimitParentContributions: no"
      end if

      ! LimitMaleParentContributions
      read(UnitSpec,*) DumC,DumC
      if (ToLower(trim(DumC)) == "yes") then
        if (GenderMatters) then
          if (EqualizePar1) then
            write(STDOUT,"(a)") "LimitMaleParentContributions: no"
            write(STDOUT,"(a)") "NOTE: Limit male parent contributions option ignored when contributions are to be equalized!"
            write(STDOUT,"(a)") " "
          else
            backspace(UnitSpec)
            read(UnitSpec,*) DumC,DumC,LimitPar1Min,LimitPar1Max,LimitPar1Penalty
            write(STDOUT,"(a,i9,a,i9,a,f9.5)") "LimitMaleParentContributions: yes, min ",&
              nint(LimitPar1Min),", max ",nint(LimitPar1Max),", penalty ",LimitPar1Penalty
          end if
        else
          write(STDOUT,"(a)") "LimitMaleParentContributions: no"
          write(STDOUT,"(a)") "NOTE: Limit male parent contributions option ignored when the gender file is not given!"
          write(STDOUT,"(a)") " "
        end if
      else
        if (GenderMatters) then
          write(STDOUT,"(a)") "LimitMaleParentContributions: no"
        end if
      end if

      ! LimitFemaleParentContributions
      read(UnitSpec,*) DumC,DumC
      LimitPar2Min=1.0d0
      LimitPar2Max=huge(LimitPar2Max)-1.0d0
      LimitPar2Penalty=0.0d0
      if (ToLower(trim(DumC)) == "yes") then
        if (GenderMatters) then
          if (EqualizePar2) then
            write(STDOUT,"(a)") "LimitFemaleParentContributions: no"
            write(STDOUT,"(a)") "NOTE: Limit female parent contributions option ignored when contributions are to be equalized!"
            write(STDOUT,"(a)") " "
          else
            backspace(UnitSpec)
            read(UnitSpec,*) DumC,DumC,LimitPar2Min,LimitPar2Max,LimitPar2Penalty
            write(STDOUT,"(a,i9,a,i9,a,f9.5)") "LimitFemaleParentContributions: yes, min ",&
              nint(LimitPar2Min),", max ",nint(LimitPar2Max),", penalty ",LimitPar2Penalty
          end if
        else
          write(STDOUT,"(a)") "LimitFemaleParentContributions: no"
          write(STDOUT,"(a)") "NOTE: Limit female parent contributions option ignored when the gender file is not given!"
          write(STDOUT,"(a)") " "
        end if
      else
        if (GenderMatters) then
          write(STDOUT,"(a)") "LimitFemaleParentContributions: no"
        end if
      end if

      ! AllowSelfing
      read(UnitSpec,*) DumC,DumC
      if (ToLower(trim(DumC)) == "yes") then
        SelfingAllowed=.true.
        if (.not.GenderMatters) then
          write(STDOUT,"(a)") "AllowSelfing: Yes"
          ! TODO: make it clear in the manual that when GenderMatters, we avoid selfing
          !       by design. In that case selfing can be achieved if a user includes the
          !       same individual as male and female individual in the data!
        end if
      else
        SelfingAllowed=.false.
        backspace(UnitSpec)
        read(UnitSpec,*) DumC,DumC,SelfingPenalty
        write(STDOUT,"(a,f9.5)") "AllowSelfing: no, penalty ",SelfingPenalty
      end if

      ! OldCoancestry
      read(UnitSpec,*) DumC,DumC
      if (ToLower(trim(DumC)) == "unknown") then
        InferPopInbOld=.true.
        write(STDOUT,"(a,a)") "OldCoancestry: unknown"
      else
        InferPopInbOld=.false.
        backspace(UnitSpec)
        read(UnitSpec,*) DumC,PopInbOld
        write(STDOUT,"(a,f9.5)") "OldCoancestry: ",PopInbOld
      end if

      ! TargetedRateOfPopulationInbreeding
      read(UnitSpec,*) DumC,RatePopInbTarget,PopInbPenalty
      write(STDOUT,"(a,f9.5,a,f9.5)") "TargetedRateOfPopulationInbreeding: ",RatePopInbTarget,", penalty ",PopInbPenalty

      ! IndividualInbreedingPenalty
      read(UnitSpec,*) DumC,IndInbPenalty
      write(STDOUT,"(a,f9.5)") "IndividualInbreedingPenalty: ",IndInbPenalty

      ! EvaluateFrontier
      read(UnitSpec,*) DumC,DumC
      if (ToLower(trim(DumC)) == "no") then
        EvaluateFrontier=.false.
        write(STDOUT,"(a)") "EvaluateFrontier: no"
      else
        EvaluateFrontier=.true.
        backspace(UnitSpec)
        read(UnitSpec,*) DumC,DumC,nFrontierSteps
        allocate(RatePopInbFrontier(nFrontierSteps))
        backspace(UnitSpec)
        read(UnitSpec,*) DumC,DumC,nFrontierSteps,RatePopInbFrontier(:)
        write(STDOUT,"(a,i9)") "EvaluateFrontier: yes, #steps: ",nFrontierSteps
      end if

      ! EvolutionaryAlgorithmIterations
      read(UnitSpec,*) DumC,EvolAlgNSol,EvolAlgNGen,EvolAlgNGenBurnIn,EvolAlgNGenStop,EvolAlgStopTol,EvolAlgNGenPrint

      ! EvolutionaryAlgorithmParameters
      read(UnitSpec,*) DumC,EvolAlgCRBurnIn,EvolAlgCRLate,EvolAlgFBase,EvolAlgFHigh1,EvolAlgFHigh2

      ! Seed
      read(UnitSpec,*) DumC,DumC
      SeedFile="AlphaMateResults"//DASH//"Seed.txt"
      if ((ToLower(trim(DumC)) == "unknown") .or. (ToLower(trim(DumC)) == "none")) then
        call SetSeed(SeedFile=SeedFile,Out=Seed)
      else
        backspace(UnitSpec)
        read(UnitSpec,*) DumC,DumI
        call SetSeed(Seed=DumI,SeedFile=SeedFile,Out=Seed)
      end if
      write(STDOUT,"(a,i)") "Seed: ",Seed
      write(STDOUT,"(a)") " "

      close(UnitSpec)

      allocate(Bv(nInd))
      allocate(BvStand(nInd))
      allocate(RelMtx(nInd,nInd))
      allocate(IdC(nInd))
      allocate(Gender(nInd))
      allocate(xVec(nInd))
      allocate(nVec(nInd))
      allocate(Mate(nMat,2))

      ! --- Breeding values and relationship matrix ---

      open(newunit=UnitRelMtx,file=trim(RelMtxFile),status="old")
      open(newunit=UnitBv,file=trim(BvFile),status="old")

      do i=1,nInd
        read(UnitRelMtx,*) IdC(i),RelMtx(:,i)
      end do
      close(UnitRelMtx)

      do i=1,nInd
        read(UnitBv,*) IdCTmp,BvTmp
        do j=1,nInd
          if (trim(IdCTmp) == trim(IdC(j))) then
            Bv(j)=BvTmp
            exit
          end if
        end do
      end do
      close(UnitBv)

      ! --- Gender ---

      if (.not.GenderMatters) then
        Gender(:)=0
      else
        nMal=0
        nFem=0
        open(newunit=UnitGender,file=trim(GenderFile),status="old")

        do i=1,nInd
          read(UnitGender,*) IdCTmp,GenderTmp
          if (GenderTmp == 1) then
            nMal=nMal+1
          elseif (GenderTmp == 2) then
            nFem=nFem+1
          else
            write(STDERR,"(a)") "ERROR: Gender code must be either 1 for male individuals or 2 for female individuals!"
            write(STDERR,"(a,i6,a,i3)") "ERROR: ",i,IdCTmp,GenderTmp
            write(STDERR,"(a)") " "
            stop 1
          end if
          do j=1,nInd
            if (trim(IdCTmp) == trim(IdC(j))) then
              Gender(j)=GenderTmp
              exit
            end if
          end do
        end do
        close(UnitGender)
        write(STDOUT,"(a)") "NOTE: The following number of male and female individuals found in the data:"
        write(STDOUT,"(a,i)") "NOTE: Number of   males:",nMal
        write(STDOUT,"(a,i)") "NOTE: Number of females:",nFem
        write(STDOUT,"(a)") " "
        if (nPar1 > nMat) then
          write(STDERR,"(a)") "ERROR: The number of male parents can not be larger than the number of males:"
          write(STDERR,"(a,i)") "ERROR: Number of male parents:",nPar1
          write(STDERR,"(a,i)") "ERROR: Number of        males:",nMal
          write(STDERR,"(a)") " "
          stop 1
        end if
        if (nPar2 > nFem) then
          write(STDERR,"(a)") "ERROR: The number of female parents can not be larger than the number of females:"
          write(STDERR,"(a,i)") "ERROR: Number of female parents:",nPar2
          write(STDERR,"(a,i)") "ERROR: Number of        females:",nFem
          write(STDERR,"(a)") " "
          stop 1
        end if
      end if

      ! --- Shuffle the data ---

      ! To avoid having good animals together - better for cross-overs in Evolutionary algorithms

      allocate(Order(nInd))
      call RandomOrder(Order,nInd)
      IdC(:)=IdC(Order)
      Bv(:)=Bv(Order)
      Gender(:)=Gender(Order)
      RelMtx(:,:)=RelMtx(Order,Order)
      deallocate(Order)

      ! --- Define ---

      if (.not.GenderMatters) then
        nPotPar1=nInd
        allocate(IdPotPar1(nPotPar1))
        allocate(nVecPar1(nPotPar1))
        do i=1,nInd
          IdPotPar1(i)=i
        end do
      else
        nPotPar1=nMal
        nPotPar2=nFem
        allocate(IdPotPar1(nPotPar1))
        allocate(IdPotPar2(nPotPar2))
        allocate(nVecPar1(nPotPar1))
        allocate(nVecPar2(nPotPar2))
        jMal=0
        jFem=0
        do i=1,nInd
          if (Gender(i) == 1) then
            jMal=jMal+1
            IdPotPar1(jMal)=i
          else
            jFem=jFem+1
            IdPotPar2(jFem)=i
          end if
        end do
        if (nPar1 > nPotPar1) then
          nPar1=nPotPar1
          write(STDOUT,"(a)") "NOTE: The number of male parents reduced to the number of male individuals!"
          write(STDOUT,"(a)") " "
        end if
        if (nPar2 > nPotPar2) then
          nPar2=nPotPar2
          write(STDOUT,"(a)") "NOTE: The number of female parents reduced to the number of female individuals!"
          write(STDOUT,"(a)") " "
        end if
      end if
    end subroutine

    !###########################################################################

    subroutine SetInbreedingParameters

      implicit none

      integer(int32) :: i,UnitInbree

      real(real64) :: Tmp

      ! Old inbreeding
      if (InferPopInbOld) then
        PopInbOld=0.0d0
        do i=1,nInd
          Tmp=RelMtx(i,i)-1.0d0
          if (Tmp < 0.0d0) then
            write(STDERR,"(a)") "ERROR: Relationship matrix must have diagonals equal or more than 1.0!"
            write(STDERR,"(a)") " "
            stop 1
          end if
          PopInbOld=PopInbOld+Tmp
        end do
        PopInbOld=PopInbOld/dble(nInd)
      end if

      ! New inbreeding
      PopInbTarget=PopInbOld*(1.0d0-RatePopInbTarget)+RatePopInbTarget

      ! Report
      write(STDOUT,"(a,f9.5)") "Old coancestry: ",PopInbOld
      write(STDOUT,"(a,f9.5)") "Targeted rate of 'population' inbreeding: ",RatePopInbTarget
      write(STDOUT,"(a,f9.5)") "Targeted 'population' inbreeding: ",PopInbTarget
      write(STDOUT,"(a)") " "

      open(newunit=UnitInbree,file="AlphaMateResults"//DASH//"ConstraintPopulationInbreeding.txt",status="unknown")
      write(UnitInbree,"(a,f)") "Old_coancestry_defined, ",PopInbOld
      write(UnitInbree,"(a,f)") "Targeted_rate_of_population_inbreeding_defined, ",RatePopInbTarget
      write(UnitInbree,"(a,f)") "Targeted_population_inbreeding_defined, ",PopInbTarget
      close(UnitInbree)
    end subroutine

    !###########################################################################

    subroutine AlphaMateSearch

      implicit none

      integer(int32) :: i,j,UnitInbree,UnitMating,UnitContri,UnitFrontier,nTmp,Rank(nInd)

      real(real64) :: BvMean,BvStdDev,PopInbTargetHold,RatePopInbTargetHold
      real(real64),allocatable :: InitEqual(:,:)

      character(len=300) :: EvolAlgLogFile

      type(EvolAlgCrit) :: CritMin,CritOpt,Crit

      ! --- Scale breeding values ---

      BvMean=sum(Bv)/dble(nInd)
      BvStand(:)=Bv(:)-BvMean
      BvStdDev=sqrt(sum(BvStand(:)*BvStand(:))/dble(nInd))
      BvStand(:)=(Bv(:)-BvMean)/BvStdDev

      ! --- Optimise for minimum inbreeding ---

      if (MinThenOpt) then
        write(STDOUT,"(a)") "Optimise for minimum inbreeding:"
        write(STDOUT,"(a)") " "

        EvolAlgLogFile="AlphaMateResults"//DASH//"OptimisationLog1MinimumInbreeding.txt"
        if (GenderMatters) then
          nTmp=nPotPar1+nPotPar2+nMat ! TODO: add PAGE dimension
        else
          nTmp=nPotPar1+nMat          ! TODO: add PAGE dimension
        end if
        allocate(InitEqual(nTmp,1))
        InitEqual(:,1)=1.0d0 ! A solution that would give equal contribution for everybody
        call EvolAlgDE(nParam=nTmp,nSol=EvolAlgNSol,Init=InitEqual,nGen=EvolAlgNGen,nGenBurnIn=EvolAlgNGenBurnIn,&
                       nGenStop=EvolAlgNGenStop,StopTolerance=EvolAlgStopTol,&
                       nGenPrint=EvolAlgNGenPrint,File=EvolAlgLogFile,CritType="Min",&
                       CRBurnIn=EvolAlgCRBurnIn,CRLate=EvolAlgCRLate,&
                       FBase=EvolAlgFBase,FHigh1=EvolAlgFHigh1,FHigh2=EvolAlgFHigh2,&
                       CalcCriterion=FixSolMateAndCalcCrit,&
                       LogHeader=EvolAlgLogHeaderForAlphaMate,Log=EvolAlgLogForAlphaMate,&
                       BestCriterion=CritMin)
        deallocate(InitEqual)
        GainMinStand=CritMin%GainStand

        open(newunit=UnitContri,file="AlphaMateResults"//DASH//"ContribAndMatingsPerIndivMinimumInbreeding.txt",status="unknown")
        !                           12345678901   12345678901   12345678901   12345678901   12345678901   12345678901
        write(UnitContri,"(6a11)") "     OrigId","     Gender","      Merit"," AvgCoances"," Contribute","   nMatings"
        call MrgRnk(nVec,Rank)
        do i=nInd,1,-1 ! MrgRnk ranks small to large
          j=Rank(i)
          write(UnitContri,"(a11,i11,3f11.5,i11)") trim(IdC(j)),Gender(j),Bv(j),0.5d0*sum(RelMtx(:,j))/dble(nInd),xVec(j),nVec(j)
        end do
        close(UnitContri)

        open(newunit=UnitMating,file="AlphaMateResults"//DASH//"MatingListMinimumInbreeding.txt",status="unknown")
        !                           12345678901   12345678901   12345678901
        write(UnitMating,"(3a11)") "     Mating","    Parent1","    Parent2"
        do i=1,nMat
          write(UnitMating,"(3i11)") i,Mate(i,:)
        end do
        close(UnitMating)

        if (PopInbOld > CritMin%PopInb) then
          write(STDOUT,"(a)") "NOTE: Old coancestry is higher than the minimum group coancestry (x'Ax/2) under no selection."
          write(STDOUT,"(a)") "NOTE: Resetting the old coancestry to the minimum group coancestry under no selection and"
          write(STDOUT,"(a)") "NOTE:   recomputing the targeted 'population' inbreeding."
          write(STDOUT,"(a)") " "
          PopInbOld=CritMin%PopInb
          PopInbTarget=(1.0d0-RatePopInbTarget)*PopInbOld+RatePopInbTarget
          CritMin%RatePopInb=0.0d0
          write(STDOUT,"(a,f9.5)") "Old coancestry: ",PopInbOld
          write(STDOUT,"(a,f9.5)") "Targeted rate of 'population' inbreeding: ",RatePopInbTarget
          write(STDOUT,"(a,f9.5)") "Targeted 'population' inbreeding:",PopInbTarget
          write(STDOUT,"(a)") " "
        end if

        if (PopInbTarget < CritMin%PopInb) then
          write(STDERR,"(a)") "ERROR: Targeted 'population' inbreeding is lower than the group coancestry (x'Ax/2) under no selection."
          write(STDERR,"(a)") "ERROR: Can not optimise!"
          write(STDERR,"(a)") " "
          stop 1
        end if

        open(newunit=UnitInbree,file="AlphaMateResults"//DASH//"ConstraintPopulationInbreeding.txt",status="old")
        write(UnitInbree,"(a,f)") "Old_coancestry_redefined, ",PopInbOld
        write(UnitInbree,"(a,f)") "Targeted_rate_of_population_inbreeding_redefined, ",RatePopInbTarget
        write(UnitInbree,"(a,f)") "Targeted_population_inbreeding_redefined, ",PopInbTarget
        close(UnitInbree)
      else
        GainMinStand=0.0d0
      end if

      ! --- Optimise for maximum gain with constraint on inbreeding ---

      write(STDOUT,"(a)") "Optimise for maximum gain with constraint on inbreeding:"
      write(STDOUT,"(a)") " "

      EvolAlgLogFile="AlphaMateResults"//DASH//"OptimisationLog2OptimumGain.txt"
      if (GenderMatters) then
        nTmp=nPotPar1+nPotPar2+nMat ! TODO: add PAGE dimension
      else
        nTmp=nPotPar1+nMat          ! TODO: add PAGE dimension
      end if
      call EvolAlgDE(nParam=nTmp,nSol=EvolAlgNSol,nGen=EvolAlgNGen,nGenBurnIn=EvolAlgNGenBurnIn,&
                     nGenStop=EvolAlgNGenStop,StopTolerance=EvolAlgStopTol,&
                     nGenPrint=EvolAlgNGenPrint,File=EvolAlgLogFile,CritType="Opt",&
                     CRBurnIn=EvolAlgCRBurnIn,CRLate=EvolAlgCRLate,&
                     FBase=EvolAlgFBase,FHigh1=EvolAlgFHigh1,FHigh2=EvolAlgFHigh2,&
                     CalcCriterion=FixSolMateAndCalcCrit,&
                     LogHeader=EvolAlgLogHeaderForAlphaMate,Log=EvolAlgLogForAlphaMate,&
                     BestCriterion=CritOpt)

      open(newunit=UnitContri,file="AlphaMateResults"//DASH//"ContribAndMatingNbPerIndivOptimumGain.txt",status="unknown")
      !                           12345678901   12345678901   12345678901   12345678901   12345678901   12345678901
      write(UnitContri,"(6a11)") "     OrigId","     Gender","      Merit"," AvgCoances"," Contribute","   nMatings"
      call MrgRnk(nVec,Rank)
      do i=nInd,1,-1 ! MrgRnk ranks small to large
        j=Rank(i)
        write(UnitContri,"(a11,i11,3f11.5,i11,2f11.5)") trim(IdC(j)),Gender(j),Bv(j),0.5d0*sum(RelMtx(:,j))/dble(nInd),xVec(j),nVec(j)
      end do
      close(UnitContri)

      open(newunit=UnitMating,file="AlphaMateResults"//DASH//"MatingListOptimumGain.txt",status="unknown")
      !                           12345678901   12345678901   12345678901
      write(UnitMating,"(3a11)") "     Mating","    Parent1","    Parent2"
      do i=1,nMat
        write(UnitMating,"(3i11)") i,Mate(i,:)
      end do
      close(UnitMating)

      ! --- Evaluate the full frontier ---

      if (EvaluateFrontier) then
        write(STDOUT,"(a)") "Evaluate the full frontier (this might take quite some time!):"
        write(STDOUT,"(a)") " "

        open(newunit=UnitFrontier,file="AlphaMateResults"//DASH//"Frontier.txt",status="unknown")
        !                             12345678901   12345678901   12345678901   12345678901   12345678901   12345678901   12345678901   12345678901
        write(UnitFrontier,"(8a11)") "       Step","  Objective","       Gain","  GainStand"," PopInbreed"," RatePopInb"," PopInbree2"," IndInbreed"
        j=0
        if (MinThenOpt) then
          j=j+1
          write(UnitFrontier,"(i11,7f11.5)") j,CritMin%Value,CritMin%Gain,CritMin%GainStand,CritMin%PopInb,CritMin%RatePopInb,CritMin%PopInb2,CritMin%IndInb
        end if
        j=j+1
        write(UnitFrontier,"(i11,7f11.5)")   j,CritOpt%Value,CritOpt%Gain,CritOpt%GainStand,CritOpt%PopInb,CritOpt%RatePopInb,CritMin%PopInb2,CritOpt%IndInb

        ! Hold old results
        PopInbTargetHold=PopInbTarget
        RatePopInbTargetHold=RatePopInbTarget

        ! Evaluate
        do i=1,nFrontierSteps
          RatePopInbTarget=RatePopInbFrontier(i)
          PopInbTarget=(1.0d0-RatePopInbTarget)*PopInbOld+RatePopInbTarget
          write(STDOUT,"(a,i3,a,i3,a,f8.5)") "Step ",i," out of ",nFrontierSteps, " for the rate of inbreeding of",RatePopInbTarget
          write(STDOUT,"(a)") ""
          EvolAlgLogFile="AlphaMateResults"//DASH//"OptimisationLog"//Int2Char(i)//".txt"
          if (GenderMatters) then
            nTmp=nPotPar1+nPotPar2+nMat ! TODO: add PAGE dimension
          else
            nTmp=nPotPar1+nMat          ! TODO: add PAGE dimension
          end if
          call EvolAlgDE(nParam=nTmp,nSol=EvolAlgNSol,nGen=EvolAlgNGen,nGenBurnIn=EvolAlgNGenBurnIn,&
                         nGenStop=EvolAlgNGenStop,StopTolerance=EvolAlgStopTol,&
                         nGenPrint=EvolAlgNGenPrint,File=EvolAlgLogFile,CritType="Opt",&
                         CRBurnIn=EvolAlgCRBurnIn,CRLate=EvolAlgCRLate,&
                         FBase=EvolAlgFBase,FHigh1=EvolAlgFHigh1,FHigh2=EvolAlgFHigh2,&
                         CalcCriterion=FixSolMateAndCalcCrit,&
                         LogHeader=EvolAlgLogHeaderForAlphaMate,Log=EvolAlgLogForAlphaMate,&
                         BestCriterion=Crit)
          write(UnitFrontier,"(i11,7f11.5)") j+i,Crit%Value,Crit%Gain,Crit%GainStand,Crit%PopInb,Crit%RatePopInb,Crit%PopInb2,Crit%IndInb
          if ((RatePopInbTarget-Crit%RatePopInb) > 0.01d0) then
            write(STDOUT,"(a,f8.5)") "NOTE: Could not achieve the rate of 'population' inbreeding of ",RatePopInbTarget
            write(STDOUT,"(a)") "NOTE: Stopping the frontier evaluation."
            write(STDOUT,"(a)") ""
            exit
          end if
        end do

        ! Put back old results
        PopInbTarget=PopInbTargetHold
        RatePopInbTarget=RatePopInbTargetHold

        close(UnitFrontier)

      end if
    end subroutine

    !###########################################################################

    subroutine EvolAlgLogHeaderForAlphaMate(LogUnit)
      implicit none
      integer(int32),intent(in) :: LogUnit
      !                         12345678901   12345678901   12345678901   12345678901   12345678901   12345678901   12345678901   12345678901   12345678901   12345678901
      write(STDOUT, "(12a11)") "       Step"," AcceptRate","  Criterion","  Penalties","       Gain","  GainStand"," PopInbreed"," RatePopInb"," PopInbree2"," IndInbreed"
      write(LogUnit,"(12a11)") "       Step"," AcceptRate","  Criterion","  Penalties","       Gain","  GainStand"," PopInbreed"," RatePopInb"," PopInbree2"," IndInbreed"
    end subroutine

    !###########################################################################

    subroutine EvolAlgLogForAlphaMate(LogUnit,Gen,AcceptRate,Criterion)
      implicit none
      integer(int32),intent(in)    :: LogUnit
      integer(int32),intent(in)    :: Gen
      real(real64),intent(in)      :: AcceptRate
      type(EvolAlgCrit),intent(in) :: Criterion
      write(STDOUT, "(i11,11f11.5)") Gen,AcceptRate,Criterion%Value,Criterion%Penalty,Criterion%Gain,Criterion%GainStand,Criterion%PopInb,Criterion%RatePopInb,Criterion%PopInb2,Criterion%IndInb
      write(LogUnit,"(i11,11f11.5)") Gen,AcceptRate,Criterion%Value,Criterion%Penalty,Criterion%Gain,Criterion%GainStand,Criterion%PopInb,Criterion%RatePopInb,Criterion%PopInb2,Criterion%IndInb
    end subroutine

    !###########################################################################

    subroutine InitialiseAlphaMateCrit(This,Value)
      implicit none
      type(EvolAlgCrit),intent(out)    :: This
      real(real64),intent(in),optional :: Value
      if (present(Value)) then
        This%Value=Value
      else
        This%Value=0.0d0
      end if
      This%Penalty=0.0d0
      This%Gain=0.0d0
      This%GainStand=0.0d0
      This%PopInb=0.0d0
      This%RatePopInb=0.0d0
      This%PopInb2=0.0d0
      This%IndInb=0.0d0
    end subroutine

    !###########################################################################

    subroutine FixSolMateAndCalcCrit(Sol,CritType,Criterion)

      implicit none

      ! Arguments
      real(real64),intent(inout)    :: Sol(:)    ! Solution
      character(len=*),intent(in)   :: CritType  ! Type of criterion (Min,Opt)
      type(EvolAlgCrit),intent(out) :: Criterion ! Criterion of the solution

      ! Other
      integer(int32) :: i,j,k,l,g,nCumMat,RankSol(nInd),SolInt(nInd),MatPar2(nMat)
      integer(int32) :: TmpMin,TmpMax

      real(real64) :: TmpVec(nInd,1),TmpR!,RanNum

      ! Criterion base
      call InitialiseAlphaMateCrit(Criterion)

      ! The solution (based on the mate selection driver) has:
      ! - nPotPar1 individual contributions for "parent1" (males when GenderMatters)
      ! - nPotPar2 individual contributions for "parent2" (females when GenderMatters, present only when GenderMatters)
      ! - nMat     rankings of parent1 1:nMat matings to pair with 1:nPotPar2 "parent2" (see bellow)
      ! - nInd     rankings for genome editing

      ! Say we have Sol=(|0,2,0,1|...|2.5,1.5,1.0) then we mate
      ! - male 2 with the first  available female (rank 2.5)
      ! - male 2 with the second available female (rank 1.5)
      ! - male 4 with the third  available female (rank 1.0)

      ! --- Parse the mate selection driver (=Is the solution valid?) ---

      ! The approach bellow assures that we have nMat contributions for each of
      ! the two parent sets. It does this by ranking internal solution values and
      ! traverses from top to the defined number of parents checking when the sum of
      ! interegrised values gives nMat. If values bellow 0.5 are found, they are
      ! changed to 1 contribution. If this still does not give nMat, then we start
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
        g=1
      else
        g=2
      end if
      ! ... find ranks to find the top values
      if (.not.(EqualizePar1 .and. (nPar1 == nPotPar1))) then
        call MrgRnk(Sol(1:nPotPar1),RankSol(1:nPotPar1))
        RankSol(1:nPotPar1)=RankSol(nPotPar1:1:-1) ! MrgRnk ranks small to large
      end if
      ! ... handle cases with equalized contributions
      if (EqualizePar1) then
        if (nPar1 == nPotPar1) then
          ! ... set integers to all the values (no need for sorting here)
          Sol(1:nPotPar1)=dble(nMat*g)/dble(nPar1)
        else
          ! ... set integers to the top values
          Sol(RankSol(1:nPar1))=dble(nMat*g)/dble(nPar1)
          ! TODO: anything better to preserve the order of non contributing individuals? See below!
          Sol(RankSol((nPar1+1):nPotPar1))=0.0d0
          ! Sol(RankSol((nPar1+1):nPotPar1))=-1.0d0
        end if
      else
        ! ... handle cases with unequal contributions
        ! ... work for the defined number or parents
        nCumMat=0
        do i=1,nPar1
          j=RankSol(i)
          ! ... these would have "zero" contributions and if we hit them then they need at least minimum contrib
          if (Sol(j) < LimitPar1Min) then
            Sol(j)=LimitPar1Min
          end if
          ! ... but not above max allowed
          if (Sol(j) > LimitPar1Max) then
            Sol(j)=LimitPar1Max
          end if
          ! ... accumulate and check if we reached nMat
          nCumMat=nCumMat+nint(Sol(j)) ! internally real, externally integer
          if (nCumMat >= nMat*g) then
            ! ... there should be exactly nMat contributions
            if (nCumMat > nMat*g) then
              Sol(j)=Sol(j)-dble(nCumMat-nMat*g)
              if (nint(Sol(j)) < LimitPar1Min) then
                TmpR=LimitPar1Penalty*(LimitPar1Min-nint(Sol(j)))
                Criterion%Value=Criterion%Value-TmpR
                Criterion%Penalty=Criterion%Penalty+TmpR
              end if
              nCumMat=nMat*g
            end if
            exit
          end if
        end do
        ! ... increment i if we have hit the exit, do loop would have ended with i=nPar1+1
        if (i <= nPar1) then
          i=i+1
        end if
        ! ... the other individuals do not contribute
        if (i <= nPotPar1) then ! "=" to capture incremented i+1 on the do loop exit
          ! ... zero (the same for all ind so no order)
          Sol(RankSol(i:nPotPar1))=0.0d0
          ! ... negative (the same for all ind so no order)
          ! Sol(RankSol(i:nPotPar1))=-1.0d0
          ! ... negative (variable with partially preserving order)
          !     Found faster convergence than with properly decreasing negative values?
          !     I guess it adds some more randomness, i.e., it causes suffling of individuals that just did not make it onto the mating list.
          ! Sol(RankSol(i:nPotPar1))=sign(Sol(RankSol(i:nPotPar1)),-1.0d0)
          ! ... negative and properly decreasing
          ! TmpR=maxval(Sol(RankSol(i:nPotPar1)))
          ! if (TmpR > 0.0d0) then
          !     Sol(RankSol(i:nPotPar1))=Sol(RankSol(i:nPotPar1))-abs(TmpR)
          ! end if
          ! ... negative (random so no order)
          ! do j=i,nPotPar1 ! TODO: really need this loop?
          !   call random_number(RanNum)
          !   Sol(RankSol(j))=-1.0d0*RanNum
          ! end do
        end if
        ! ... nMat still not reached?
        do while (nCumMat < nMat*g)
          ! ... add more contributions
          do i=nPar1,1,-1 ! to low ranked individuals (to avoid local optima)
            j=RankSol(i)
            Sol(j)=Sol(j)+1.0d0
            ! ... accumulate and check if we reached nMat
            nCumMat=nCumMat+1
            if (nCumMat >= nMat*g) then
              exit
            end if
          end do
        end do
      end if

      ! "Parent2"
      if (GenderMatters) then
        ! ... find ranks to find the top values
        if (.not.(EqualizePar2 .and. (nPar2 == nPotPar2))) then
          call MrgRnk(Sol((nPotPar1+1):(nPotPar1+nPotPar2)),RankSol(1:nPotPar2))
          RankSol(1:nPotPar2)=RankSol(nPotPar2:1:-1) ! MrgRnk ranks small to large
        end if
        ! ... handle cases with equalized contributions
        if (EqualizePar2) then
          if (nPar2 == nPotPar2) then
            ! ... set integers to all the values (no need for sorting here)
            Sol((nPotPar1+1):(nPotPar1+nPotPar2))=dble(nMat)/dble(nPar2)
          else
            ! ... set integers to the top values
            Sol(nPotPar1+RankSol(1:nPar2))=dble(nMat)/dble(nPar2)
            ! TODO: anything better to preserve the order of non contributing individuals? See below!
            Sol(nPotPar1+RankSol((nPar2+1):nPotPar2))=0.0d0
            ! Sol(nPotPar1+RankSol((nPar2+1):nPotPar2))=-1.0d0
          end if
        else
          ! ... handle cases with unequal contributions
          ! ... work for the defined number or parents
          nCumMat=0
          do i=1,nPar2
            j=nPotPar1+RankSol(i)
            ! ... these would have "zero" contributions and if we hit them then they need at least minimum contrib
            if (Sol(j) < LimitPar2Min) then
              Sol(j)=LimitPar2Min
            end if
            ! ... but not above max allowed
            if (Sol(j) > LimitPar2Max) then
              Sol(j)=LimitPar2Max
            end if
            ! ... accumulate and check if we reached nMat
            nCumMat=nCumMat+nint(Sol(j)) ! internally real, externally integer
            if (nCumMat >= nMat) then
              ! ... there should be exactly nMat contributions
              if (nCumMat > nMat) then
                Sol(j)=Sol(j)-dble(nCumMat-nMat)
                if (nint(Sol(j)) < LimitPar2Min) then
                  TmpR=LimitPar2Penalty*(LimitPar2Min-nint(Sol(j)))
                  Criterion%Value=Criterion%Value-TmpR
                  Criterion%Penalty=Criterion%Penalty+TmpR
                end if
                nCumMat=nMat
              end if
              exit
            end if
          end do
          ! ... increment i if we have hit the exit, do loop would have ended with i=nPar2+1
          if (i <= nPar2) then
            i=i+1
          end if
          ! ... the other individuals do not contribute
          if (i <= nPotPar2) then ! "="" to capture incremented i+1 on the do loop exit
            ! ... zero (the same for all ind so no order)
            Sol(nPotPar1+(RankSol(i:nPotPar2)))=0.0d0
            ! ... negative (the same for all ind so no order)
            ! Sol(nPotPar1+(RankSol(i:nPotPar2)))=-1.0d0
            ! ... negative (variable with partially preserving order, i.e., ~large positives become ~large negatives)
            !     Found faster convergence than with properly decreasing negative values?
            !     I guess it adds some more randomness, i.e., it causes suffling of individuals that just did not make it onto the mating list.
            ! Sol(nPotPar1+(RankSol(i:nPotPar2)))=sign(Sol(nPotPar1+(RankSol(i:nPotPar2))),-1.0d0)
            ! ... negative and properly decreasing
            ! TmpR=maxval(Sol(nPotPar1+(RankSol(i:nPotPar2))))
            ! if (TmpR > 0.0d0) then
            !     Sol(nPotPar1+(RankSol(i:nPotPar2)))=Sol(nPotPar1+(RankSol(i:nPotPar2)))-abs(TmpR)
            ! end if
            ! ... negative (random so no order)
            ! do j=i,nPotPar2 ! TODO: really need this loop?
            !   call random_number(RanNum)
            !   Sol(nPotPar1+RankSol(j))=-1.0d0*RanNum
            ! end do
          end if
          ! ... nMat still not reached?
          do while (nCumMat < nMat)
            ! ... add more contributions
            do i=nPar2,1,-1 ! to low ranked individuals (to avoid local optima)
              j=nPotPar1+RankSol(i)
              Sol(j)=Sol(j)+1.0d0
              ! ... accumulate and check if we reached nMat
              nCumMat=nCumMat+1
              if (nCumMat == nMat) then
                exit
              end if
            end do
          end do
        end if
      end if

      ! --- Genetic contributions (nVec) ---

      nVecPar1(:)=0

      ! "Parent1"
      ! ... get integer values
      SolInt(1:nPotPar1)=nint(Sol(1:nPotPar1))
      ! ... remove negatives
      do i=1,nPotPar1
        if (SolInt(i) < 0) then
          SolInt(i)=0
        end if
      end do
      ! ... map internal to external order
      nVecPar1(:)=SolInt(1:nPotPar1)
      if (.not.GenderMatters) then
        nVec(:)=nVecPar1(:)
      else
        nVec(IdPotPar1)=nVecPar1(:)
      end if

      ! "Parent2"
      if (GenderMatters) then
        nVecPar2(:)=0
        ! ... get integer values
        SolInt(1:nPotPar2)=nint(Sol((nPotPar1+1):(nPotPar1+nPotPar2)))
        ! ... remove negatives
        do i=1,nPotPar2
          if (SolInt(i) < 0) then
            SolInt(i)=0
          end if
        end do
        ! ... map internal to external order
        nVecPar2(:)=SolInt(1:nPotPar2)
        nVec(IdPotPar2)=nVecPar2(:)
      end if

      xVec(:)=dble(nVec(:))/(dble(2*nMat))

      ! --- Genetic gain ---

      Criterion%Gain=dot_product(xVec,Bv)
      Criterion%GainStand=dot_product(xVec,BvStand)

      TmpR=Criterion%GainStand-GainMinStand
      if (ToLower(trim(CritType)) == "min") then
        Criterion%Value=Criterion%Value+1.0d-6*TmpR
      else
        Criterion%Value=Criterion%Value+TmpR
      end if

      ! --- Future population inbreeding (=selected group coancestry) ---

      ! xA
      do i=1,nInd
        TmpVec(i,1)=dot_product(xVec,RelMtx(:,i))
      end do
      ! xAx
      Criterion%PopInb=0.5d0*dot_product(TmpVec(:,1),xVec)

      ! TODO: Should group coancestry account for self-coancestries or not (the code
      !       bellow does not include them), but then should we account for coancestries
      !       between say males (and equally between females) - we will never mate males
      !       with males etc.
      ! Perhaps useful if we do not allow selfing. How does this work with theory of mate allocations?
      !
      ! ! xA
      do i=1,nInd
          TmpVec(i,1)=dot_product(xVec,RelMtx(:,i))
          TmpVec(i,1)=TmpVec(i,1)-xVec(i)*RelMtx(i,i) ! remove self-coancestry 1+Fi
          ! TmpVec(i,1)=TmpVec(i,1)-xVec(i)*1.0d0       ! remove self-coancestry just the 1
      end do
      ! xAx
      Criterion%PopInb2=0.5d0*dot_product(TmpVec(:,1),xVec)

      ! Matrix multiplication with symmetric matrix using BLAS routine
      ! (it was ~5x slower than the above with 1.000 individuals, might be
      !  benefical with larger cases so kept in commented.)
      ! http://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3.html#ga253c8edb8b21d1b5b1783725c2a6b692
      ! Ax
      ! call dsymm(side="l",uplo="l",m=nInd,n=1,alpha=1.0d0,A=RelMtx,lda=nInd,b=xVec,ldb=nInd,beta=0,c=TmpVec,ldc=nInd)
      ! call dsymm(     "l",     "l",  nInd,  1,      1.0d0,  RelMtx,    nInd,  xVec,    nInd,     0,  TmpVec,    nInd)
      ! xAx
      ! PopInb=0.5d0*dot_product(xVec,TmpVec(:,1))
      ! print*,xVec,TmpVec,PopInb
      ! stop 1

      Criterion%RatePopInb=(Criterion%PopInb-PopInbOld)/(1.0d0-PopInbOld)

      if (ToLower(trim(CritType)) == "min") then
        Criterion%Value=Criterion%Value-Criterion%PopInb
      else
        TmpR=PopInbPenalty*abs(RatePopInbTarget-Criterion%RatePopInb)
        Criterion%Value=Criterion%Value-TmpR
        Criterion%Penalty=Criterion%Penalty+TmpR
      endif

      ! --- Mate allocation ---

      if (GenderMatters) then
        ! Distribute parent2 (=female) contributions into matings
        k=0
        do i=1,nPotPar2 ! need to loop whole nVecPar2 as some entries are zero
          do j=1,nVecPar2(i)
            k=k+1
            MatPar2(k)=i
          end do
        end do
        ! Reorder parent2 contributions according to the rank of matings
        call MrgRnk(Sol((nPotPar1+nPotPar2+1):(nPotPar1+nPotPar2+nMat)),RankSol(1:nMat))
        MatPar2(:)=MatPar2(RankSol(1:nMat))
      else
        ! Distribute one half of contributions into matings
        k=0
        do while (k < nMat)
          do i=1,nPotPar1 ! need to loop whole nVecPar1 as some entries are zero
            l=nVecPar1(i)/2
            do j=1,l
              if (k == nMat) then
                exit
              end if
              k=k+1
              MatPar2(k)=i
              nVecPar1(i)=nVecPar1(i)-1
            end do
          end do
        end do
        ! Reorder one half of contributions according to the rank of matings
        call MrgRnk(Sol((nPotPar1+1):(nPotPar1+nMat)),RankSol(1:nMat))
        MatPar2(:)=MatPar2(RankSol(1:nMat))
      end if

      ! Pair the contributions
      k=nMat ! MrgRnk ranks small to large
      Mate(:,:)=0
      if (GenderMatters .or. SelfingAllowed) then
        ! When GenderMatters selfing can not happen (we have two distinct sets of parents,
        ! unless the user adds individuals of one sex in both sets) and when SelfingAllowed
        ! we do not need to care about it - faster code
        do i=1,nPotPar1
          do j=1,nVecPar1(i)
            ! print*,k,i,j
            Mate(k,1)=i
            Mate(k,2)=MatPar2(k)
            k=k-1
          end do
        end do
      else
        ! When .not. GenderMatters selfing can happen (we have one set of parents)
        ! and when .not. SelfingAllowed we do need to care about it - slower code
        do i=1,nPotPar1
          do j=1,nVecPar1(i)
            Mate(k,1)=i
            if (MatPar2(k) == i) then
              ! Try to avoid selfing by swapping the MatPar2 and Rank elements
              do l=k,1,-1
                if (MatPar2(l) /= i) then
                  MatPar2([k,l])=MatPar2([l,k])
                  Sol(nPotPar1+RankSol([k,l]))=Sol(nPotPar1+RankSol([l,k]))
                  exit
                end if
              end do
              if (l < 1) then ! Above loop ran out without finding a swap
                Criterion%Value=Criterion%Value-SelfingPenalty
                Criterion%Penalty=Criterion%Penalty+SelfingPenalty
              end if
            end if
            Mate(k,2)=MatPar2(k)
            k=k-1
          end do
        end do
      end if

      ! --- Individual (=progeny) inbreeding ---

      TmpR=0.0d0
      do i=1,nMat
        ! Try to speed-up retrieval by targeting lower triangle
        TmpMin=minval([Mate(i,1),Mate(i,2)])
        TmpMax=maxval([Mate(i,1),Mate(i,2)])
        TmpR=TmpR+0.5d0*RelMtx(TmpMax,TmpMin)
      end do
      Criterion%IndInb=TmpR/dble(nMat)
      TmpR=IndInbPenalty*Criterion%IndInb
      Criterion%Value=Criterion%Value-TmpR
      Criterion%Penalty=Criterion%Penalty+TmpR
    end subroutine

    !###########################################################################
end module

!###############################################################################

program AlphaMate

  use ISO_Fortran_Env, STDIN=>input_unit,STDOUT=>output_unit,STDERR=>error_unit
  use IFPort,only : SystemQQ
  use AlphaMateModule

  implicit none

  real(real32) :: Start,Finish

  logical :: Success

  call cpu_time(Start)
  call AlphaMateTitles

  Success=SystemQQ(RMDIR//" AlphaMateResults")
  if (.not.Success) then
    write(STDERR,"(a)") "ERROR: Failed to remove the old output folder (AlphaMateResults)!"
    write(STDERR,"(a)") " "
    stop 1
  end if
  Success=SystemQQ(MKDIR//" AlphaMateResults")
  if (.not.Success) then
    write(STDERR,"(a)") "ERROR: Failed to make the output folder (AlphaMateResults)!"
    write(STDERR,"(a)") " "
    stop 1
  end if

  call ReadSpecAndDataForAlphaMate
  call SetInbreedingParameters
  call AlphaMateSearch
  call cpu_time(Finish)

  write(STDOUT,"(a,f20.4,a)") "Time duration of AlphaMate: ",Finish-Start," seconds"
  write(STDOUT,"(a)") " "
end program

!###############################################################################