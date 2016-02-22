#ifdef BINARY
#define BINFILE ,form="unformatted"
#else
#DEFINE BINFILE
#endif

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#ifdef OS_UNIX

#DEFINE DASH "/"
#DEFINE COPY "cp"
#DEFINE MKDIR "mkdir -p"
#DEFINE RMDIR "rm -r"
#DEFINE RM "rm"
#DEFINE RENAME "mv"

#else
#DEFINE DASH "\"
#DEFINE COPY "copy"
#DEFINE MKDIR "md"
#DEFINE RMDIR "RMDIR /S"
#DEFINE RM "del"
#DEFINE RENAME "MOVE /Y"
#endif

!###############################################################################

module AlphaSuiteModule

    contains

        !#######################################################################

        subroutine SetSeed(Seed,SeedFile,Out)
            implicit none
            ! Arguments
            integer,intent(in),optional  :: Seed     ! A number to initialize RNG
            character(len=*)             :: SeedFile ! File to save the seed in
            integer,intent(out),optional :: Out      ! Make the seed value available outside

            ! Other
            integer,allocatable :: SeedList(:)
            integer :: Size,Unit

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
            endif
            open(newunit=Unit,file=trim(SeedFile),status="unknown")
            write(Unit,*) SeedList(1)
            close(Unit)

            ! Output
            if (present(Out)) then
                Out=SeedList(1)
            endif
            deallocate(SeedList)
        end subroutine SetSeed

        !#######################################################################

        subroutine RandomOrder(order,n)
            implicit none

            !     Generate a random ordering of the integers 1 ... n.

            integer,intent(in)  :: n
            integer,intent(out) :: order(n)

            !     Local variables

            integer :: i,j,k
            double precision :: wk

            do i=1,n
                order(i)=i
            enddo

            !     Starting at the end, swap the current last indicator with one
            !     randomly chosen from those preceeding it.

            do i=n,2,-1
                call random_number(wk)
                j=1 + i * wk
                if (j < i) then
                    k=order(i)
                    order(i)=order(j)
                    order(j)=k
                endif
            enddo

            return
        end subroutine RandomOrder

        !#######################################################################

        subroutine CountLines(FileName,nLines)
            implicit none

            character(len=*) :: FileName
            character(len=300) :: DumC
            integer :: nLines,f,Unit

            nLines=0
            open(newunit=Unit,file=trim(FileName),status="old")
            do
                read(Unit,*,iostat=f) DumC
                nLines=nLines+1
                if (f /= 0) then
                    nLines=nLines-1
                    exit
                endif
            enddo
            close(Unit)
        end subroutine CountLines

        !#######################################################################

        function Int2Char(i) result(Res)
            ! TODO: source of this?
            implicit none
            character(:),allocatable :: Res
            integer,intent(in) :: i
            character(range(i)+2) :: Tmp
            write(Tmp,'(i0)') i
            Res=trim(Tmp)
        end function

        !#######################################################################
end module AlphaSuiteModule

!###############################################################################

Module m_MrgRnk
    ! http://www.fortran-2000.com/rank/MrgRnk.f90 (2016-02-15)
    Integer, Parameter :: kdp = selected_real_kind(15)
    public :: MrgRnk
    private :: kdp
    private :: R_MrgRnk, I_MrgRnk, D_MrgRnk
    interface MrgRnk
      module procedure D_MrgRnk, R_MrgRnk, I_MrgRnk
    end interface MrgRnk

    contains

        !#######################################################################

        Subroutine D_MrgRnk (XDONT, IRNGT)
            ! __________________________________________________________
            !   MRGRNK = Merge-sort ranking of an array
            !   For performance reasons, the first 2 passes are taken
            !   out of the standard loop, and use dedicated coding.
            ! __________________________________________________________
            ! __________________________________________________________
              Real (kind=kdp), Dimension (:), Intent (In) :: XDONT
              Integer, Dimension (:), Intent (Out) :: IRNGT
            ! __________________________________________________________
              Real (kind=kdp) :: XVALA, XVALB
            !
              Integer, Dimension (SIZE(IRNGT)) :: JWRKT
              Integer :: LMTNA, LMTNC, IRNG1, IRNG2
              Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
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
        End Subroutine D_MrgRnk

        !#######################################################################

        Subroutine R_MrgRnk (XDONT, IRNGT)
            ! __________________________________________________________
            !   MRGRNK = Merge-sort ranking of an array
            !   For performance reasons, the first 2 passes are taken
            !   out of the standard loop, and use dedicated coding.
            ! __________________________________________________________
            ! _________________________________________________________
              Real, Dimension (:), Intent (In) :: XDONT
              Integer, Dimension (:), Intent (Out) :: IRNGT
            ! __________________________________________________________
              Real :: XVALA, XVALB
            !
              Integer, Dimension (SIZE(IRNGT)) :: JWRKT
              Integer :: LMTNA, LMTNC, IRNG1, IRNG2
              Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
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
        End Subroutine R_MrgRnk

        !#######################################################################

        Subroutine I_MrgRnk (XDONT, IRNGT)
            ! __________________________________________________________
            !   MRGRNK = Merge-sort ranking of an array
            !   For performance reasons, the first 2 passes are taken
            !   out of the standard loop, and use dedicated coding.
            ! __________________________________________________________
            ! __________________________________________________________
              Integer, Dimension (:), Intent (In)  :: XDONT
              Integer, Dimension (:), Intent (Out) :: IRNGT
            ! __________________________________________________________
              Integer :: XVALA, XVALB
            !
              Integer, Dimension (SIZE(IRNGT)) :: JWRKT
              Integer :: LMTNA, LMTNC, IRNG1, IRNG2
              Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
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
        End Subroutine I_MrgRnk

        !#######################################################################
end module m_MrgRnk

!###############################################################################

module AlphaMateModule
#ifdef f2003
    use,intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                             stdout=>output_unit, &
                                             stderr=>error_unit
#else
#define stdin  5
#define stdout 6
#define stderr 0
#endif
    use IFPort,only : SystemQQ,SortQQ
    use m_MrgRnk,only : MrgRnk
    use AlphaSuiteModule,only : SetSeed,CountLines,Int2Char

    implicit none

    integer :: nInd,nMat,nPar,nPotPar1,nPotPar2,nMal,nFem,nMalPar,nFemPar,nFrontierSteps
    integer :: EvolAlgNSol,EvolAlgNGen,EvolAlgNGenBurnIn,EvolAlgNGenStop,EvolAlgNGenPrint
    integer,allocatable :: Gender(:),IdPotPar1(:),IdPotPar2(:),nVec(:),nVecPar1(:),nVecPar2(:),Mate(:,:)

    ! TODO: can we not work with single precision to speed up?
    double precision :: EvolAlgStopTol
    double precision :: Gain,GainScaled,GainMinInb,GainMinInbScaled,GainOpt,GainOptScaled
    double precision :: RateInbTarget,RateInbSol,RateInbMinInb,RateInbOpt
    double precision :: InbOld,InbTarget,InbTargetRebased,InbSol,InbSolRebased,InbMinInb
    double precision :: InbMinInbRebased,InbOpt,InbOptRebased,IndInbSol,IndInbSolRebased
    double precision :: IndInbTargetRebased,RateIndInbSol,RateIndInbTarget,IndInbMinInb
    double precision :: IndInbOpt,RateIndInbMinInb,RateIndInbOpt,PopInbPenalty,IndInbPenalty
    double precision :: ValueHold,ValueHoldMinInb,ValueHoldOpt

    double precision,allocatable :: Bv(:),BvScaled(:),xVec(:),RelMtx(:,:),RateInbFrontier(:)

    character(len=300),allocatable :: IdC(:)

    logical :: GenderMatters,EqualizeMales,EqualizeFemales,InferInbOld,EvaluateFrontier

    contains

        !#######################################################################

        subroutine AlphaMateTitles
            write(stdout,"(a)") ""
            write(stdout,"(a30,a,a30)") " ","**********************"," "
            write(stdout,"(a30,a,a30)") " ","*                    *"," "
            write(stdout,"(a30,a,a30)") " ","*     AlphaMate      *"," "
            write(stdout,"(a30,a,a30)") " ","*                    *"," "
            write(stdout,"(a30,a,a30)") " ","**********************"
            write(stdout,"(a30,a,a30)") " ","VERSION:"//TOSTRING(VERS)," "
            write(stdout,"(a15,a)")     " ","Software for optimizing contributions to the next generation"
            write(stdout,"(a)") ""
            write(stdout,"(a35,a)")     " ","No Liability"
            write(stdout,"(a25,a)")     " ","Bugs to John.Hickey@roslin.ed.ac.uk"
            write(stdout,"(a)") ""
        end subroutine AlphaMateTitles

        !#######################################################################

        subroutine ReadSpecAndDataForAlphaMate
            implicit none

            integer :: i,j,DumI,jMal,jFem,nIndTmp,GenderTmp,Seed
            integer :: UnitSpec,UnitRelMtx,UnitBv,UnitGender
            integer,allocatable :: BvRank(:)

            double precision :: BvTmp

            character(len=1000) :: DumC,IdCTmp
            character(len=1000) :: RelMtxFile,BvFile,GenderFile,SeedFile

            ! --- Spec file ---

            open(newunit=UnitSpec,file="AlphaMateSpec.txt",status="old")

            read(UnitSpec,*) DumC,RelMtxFile
            read(UnitSpec,*) DumC,BvFile

            call CountLines(RelMtxFile,nInd)
            call CountLines(BvFile,nIndTmp)

            if (nIndTmp /= nInd) then
                write(stderr,"(a)") "ERROR: Number of individuals in ebv file and Relationship Matrix file is not the same!"
                stop 1
            endif

            read(UnitSpec,*) DumC,GenderFile
            if (trim(GenderFile) /= "None") then
                GenderMatters=.true.
            else
                GenderMatters=.false.
            endif

            read(UnitSpec,*) DumC,nInd
            read(UnitSpec,*) DumC,nMat
            if (nMat > nInd) then
                write(stderr,"(a)") "ERROR: Number of matings can not be larger than the number or all individuals!"
                write(stderr,"(a,i)") "ERROR: Number of individuals: ",nInd
                write(stderr,"(a,i)") "ERROR: Number of     matings: ",nMat
                stop 1
            endif
            read(UnitSpec,*) DumC,nPar
            if (nPar > nInd) then
                write(stderr,"(a)") "ERROR: Number of parents can not be larger than the number or all individuals!"
                write(stderr,"(a,i)") "ERROR: Number of individuals: ",nInd
                write(stderr,"(a,i)") "ERROR: Number of     parents: ",nPar
                stop 1
            endif
            read(UnitSpec,*) DumC,nMalPar
            read(UnitSpec,*) DumC,nFemPar
            if (GenderMatters .and. ((nMalPar+nFemPar) /= nPar)) then
                write(stderr,"(a)") "ERROR: Number of male and female parents does not match with the number of parents!"
                write(stderr,"(a,i)") "ERROR: Number of        parents: ",nPar
                write(stderr,"(a,i)") "ERROR: Number of   male parents: ",nMalPar
                write(stderr,"(a,i)") "ERROR: Number of female parents: ",nFemPar
                stop 1
            endif

            ! TODO: should there be a parameter about how many progeny to produce and what is litter size?

            read(UnitSpec,*) DumC,DumC
            if (GenderMatters .and. (trim(DumC) == "Yes")) then
                EqualizeMales=.true.
            else
                EqualizeMales=.false.
            endif
            read(UnitSpec,*) DumC,DumC
            if (GenderMatters .and. (trim(DumC) == "Yes")) then
                EqualizeFemales=.true.
            else
                EqualizeFemales=.false.
            endif

            read(UnitSpec,*) DumC,DumC
            if (trim(DumC) == "Unknown") then
                InferInbOld=.true.
            else
                InferInbOld=.false.
                backspace(UnitSpec)
                read(UnitSpec,*) DumC,InbOld
            endif
            read(UnitSpec,*) DumC,RateInbTarget,PopInbPenalty
            read(UnitSpec,*) DumC,RateIndInbTarget,IndInbPenalty

            read(UnitSpec,*) DumC,DumC
            if (trim(DumC) == "No") then
                EvaluateFrontier=.false.
            else
                EvaluateFrontier=.true.
                backspace(UnitSpec)
                read(UnitSpec,*) DumC,DumC,nFrontierSteps
                allocate(RateInbFrontier(nFrontierSteps))
                backspace(UnitSpec)
                read(UnitSpec,*) DumC,DumC,DumC,RateInbFrontier(:)
            endif

            read(UnitSpec,*) DumC,EvolAlgNSol,EvolAlgNGen,EvolAlgNGenBurnIn,EvolAlgNGenStop,EvolAlgStopTol,EvolAlgNGenPrint

            read(UnitSpec,*) DumC,DumC
            SeedFile="AlphaMateResults"//DASH//"Seed.txt"
            if ((trim(DumC) == "Unknown") .or. (trim(DumC) == "None")) then
                call SetSeed(SeedFile=SeedFile,Out=Seed)
            else
                backspace(UnitSpec)
                read(UnitSpec,*) DumC,DumI
                call SetSeed(Seed=DumI,SeedFile=SeedFile,Out=Seed)
            endif
            write(stdout,"(a,i)") "Used seed: ",Seed
            write(stdout,"(a)") " "

            close(UnitSpec)

            allocate(Bv(nInd))
            allocate(RelMtx(nInd,nInd))
            allocate(IdC(nInd))
            allocate(nVec(nInd))
            allocate(xVec(nInd))
            allocate(Mate(nMat,2))

            ! --- Breeding values and relationship matrix ---

            open(newunit=UnitRelMtx,file=trim(RelMtxFile),status="old")
            open(newunit=UnitBv,file=trim(BvFile),status="old")

            do i=1,nInd
                read(UnitRelMtx,*) IdC(i),RelMtx(:,i)
            enddo
            close(UnitRelMtx)

            do i=1,nInd
                read(UnitBv,*) IdCTmp,BvTmp
                do j=1,nInd
                    if (trim(IdCTmp) == trim(IdC(j))) then
                        Bv(j)=BvTmp
                        exit
                    endif
                enddo
            enddo
            close(UnitBv)

            ! --- Gender ---

            allocate(Gender(nInd))
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
                        write(stderr,"(a)") "ERROR: Gender code must be either 1 for males or 2 for females!"
                        write(stderr,"(a,i6,a,i3)") "ERROR: ",i,IdCTmp,GenderTmp
                        stop 1
                    endif
                    do j=1,nInd
                        if (trim(IdCTmp) == trim(IdC(j))) then
                            Gender(j)=GenderTmp
                            exit
                        endif
                    enddo
                enddo
                close(UnitGender)
            endif

            ! --- Order the data by breeding values ---

            allocate(BvRank(nInd))
            call MrgRnk(Bv,BvRank)
            BvRank(:)=BvRank(nInd:1:-1) ! MrgRnk ranks small to large
            Bv(:)=Bv(BvRank)
            Gender(:)=Gender(BvRank)
            RelMtx(:,:)=RelMtx(BvRank,BvRank)
            deallocate(BvRank)

            if (.not.GenderMatters) then
                allocate(IdPotPar1(nInd))
                allocate(IdPotPar2(nInd))
                allocate(nVecPar1(nInd))
                allocate(nVecPar2(nInd))
                nPotPar1=nInd
                nPotPar2=nInd
                IdPotPar1(:)=[1:nPotPar1]
                IdPotPar2(:)=[1:nPotPar2]
            else
                nPotPar1=nMal
                nPotPar2=nFem
                allocate(IdPotPar1(nMal))
                allocate(IdPotPar2(nFem))
                allocate(nVecPar1(nMal))
                allocate(nVecPar2(nFem))
                jMal=0
                jFem=0
                do i=1,nInd
                    if (Gender(i) == 1) then
                        jMal=jMal+1
                        IdPotPar1(jMal)=i
                    else
                        jFem=jFem+1
                        IdPotPar2(jFem)=i
                    endif
                enddo
            endif
        end subroutine ReadSpecAndDataForAlphaMate

        !#######################################################################

        subroutine SetInbreedingParameters
            implicit none
            integer :: i,UnitInbree
            double precision :: Tmp

            ! Old inbreeding
            if (InferInbOld) then
                InbOld=0.0d0
                do i=1,nInd
                    Tmp=RelMtx(i,i)-1.0d0
                    if (Tmp < 0.0d0) then
                        write(stderr,"(a)") "ERROR: Relationship matrix must have diagonals equal or more than 1.0!"
                        stop 1
                    endif
                    InbOld=InbOld+Tmp
                enddo
                InbOld=InbOld/dble(nInd)
            endif

            ! New inbreeding
            InbTarget=InbOld*(1.0d0-RateInbTarget)+RateInbTarget
            InbTargetRebased=(InbTarget-InbOld)/(1.0d0-InbOld)

            ! Report
            write(stdout,"(a,f)") "Previous coancestry (=old inbreeding): ",InbOld
            write(stdout,"(a,f)") "Targeted rate of inbreeding: ",RateInbTarget
            write(stdout,"(a,f)") "Targeted inbreeding: ",InbTarget
            write(stdout,"(a)") " "

            open(newunit=UnitInbree,file="AlphaMateResults"//DASH//"ConstraintPopulationInbreeding.txt",status="unknown")
            write(UnitInbree,"(a,f)") "Old_inbreeding_defined, ",InbOld
            write(UnitInbree,"(a,f)") "Targeted_rate_of_inbreeding_defined, ",RateInbTarget
            write(UnitInbree,"(a,f)") "Targeted_inbreeding_defined, ",InbTarget
            close(UnitInbree)
        end subroutine SetInbreedingParameters

        !#######################################################################

        subroutine AlphaMateSearch
            implicit none

            integer :: i,UnitInbree,UnitMating,UnitContri,UnitFrontier,nTmp

            double precision :: StdDev,Mean,RateInbFrontierStep
            double precision :: InbTargetRebasedHold,InbTargetHold,RateInbTargetHold

            character(len=300) :: EvolAlgLogFile

            ! --- Scale breeding values ---

            allocate(BvScaled(nInd))
            Mean=sum(Bv)/dble(nInd)
            BvScaled(:)=Bv(:)-Mean
            StdDev=sqrt(sum(BvScaled(:)*BvScaled(:))/dble(nInd))
            BvScaled(:)=(Bv(:)-Mean)/StdDev

            ! --- Optimise for minimum inbreeding ---

            write(stdout,"(a)") "Optimise for minimum inbreeding:"
            write(stdout,"(a)") " "

            EvolAlgLogFile="AlphaMateResults"//DASH//"OptimisationLog1MinimumInbreeding.txt"
            nTmp=nPotPar1+nPotPar2+nMat ! TODO: add PAGE dimension
            call EvolAlgForAlphaMate(nParam=nTmp,nSol=EvolAlgNSol,nGen=EvolAlgNGen,nGenBurnIn=EvolAlgNGenBurnIn,&
                                     nGenStop=EvolAlgNGenStop,StopTolerance=EvolAlgStopTol,&
                                     nGenPrint=EvolAlgNGenPrint,File=EvolAlgLogFile,CritType="MinInb")
            GainMinInb=Gain
            GainMinInbScaled=GainScaled
            InbMinInb=InbSol
            InbMinInbRebased=InbSolRebased
            RateInbMinInb=RateInbSol
            IndInbMinInb=IndInbSol
            RateIndInbMinInb=RateIndInbSol
            ValueHoldMinInb=ValueHold

            open(newunit=UnitContri,file="AlphaMateResults"//DASH//"ContribAndMatingsPerIndivMinimumInbreeding.txt",status="unknown")
            !                           12345678901   12345678901   12345678901   12345678901   12345678901
            write(UnitContri,"(5a11)") "         Id","     OrigId","     Gender"," Contribute","   nMatings"
            do i=1,nInd
                write(UnitContri,"(i11,a11,i11,f11.4,i11)") i,trim(IdC(i)),Gender(i),xVec(i),nVec(i)
            enddo
            close(UnitContri)

            open(newunit=UnitMating,file="AlphaMateResults"//DASH//"MatingListMinimumInbreeding.txt",status="unknown")
            !                           12345678901   12345678901
            write(UnitMating,"(2a11)") "    Parent1","    Parent2"
            do i=1,nMat
                write(UnitMating,"(2i11)") Mate(i,:)
            enddo
            close(UnitMating)

            ! --- Optimise for maximum gain with constraint on inbreeding ---

            write(stdout,"(a)") "Optimise for maximum gain with constraint on inbreeding:"
            write(stdout,"(a)") " "

            if (InbOld > InbSol) then
                write(stdout,"(a)") "NOTE: Old inbreeding is higher than the minimum group coancestry (x'Ax/2) under no selection."
                write(stdout,"(a)") "NOTE: Resetting the old inbreeding to the minimum group coancestry under no selection and"
                write(stdout,"(a)") "NOTE:   recomputing the targeted inbreeding."
                InbOld=InbSol
                InbTarget=InbOld*(1.0d0-RateInbTarget)+RateInbTarget
                InbTargetRebased=(InbTarget-InbOld)/(1.0d0-InbOld)
                InbMinInbRebased=0.0d0
                RateInbMinInb=0.0d0
                write(stdout,"(a,f)") "Previous coancestry (=old inbreeding): ",InbOld
                write(stdout,"(a,f)") "Targeted rate of inbreeding: ",RateInbTarget
                write(stdout,"(a,f)") "Targeted inbreeding:",InbTarget
                write(stdout,"(a)") " "
            endif

            if (InbSol > InbTarget) then
                write(stderr,"(a)") "ERROR: Targeted inbreeding is lower than the group coancestry (x'Ax/2) under no selection."
                write(stderr,"(a)") "ERROR: Can not optimise!"
                stop 1
            endif

            open(newunit=UnitInbree,file="AlphaMateResults"//DASH//"ConstraintPopulationInbreeding.txt",status="old")
            write(UnitInbree,"(a,f)") "Old_inbreeding_redefined, ",InbOld
            write(UnitInbree,"(a,f)") "Targeted_rate_of_inbreeding_redefined, ",RateInbTarget
            write(UnitInbree,"(a,f)") "Targeted_inbreeding_redefined, ",InbTarget
            close(UnitInbree)

            EvolAlgLogFile="AlphaMateResults"//DASH//"OptimisationLog2OptimumGain.txt"
            nTmp=nPotPar1+nPotPar2+nMat ! TODO: add PAGE dimension
            call EvolAlgForAlphaMate(nParam=nTmp,nSol=EvolAlgNSol,nGen=EvolAlgNGen,nGenBurnIn=EvolAlgNGenBurnIn,&
                                     nGenStop=EvolAlgNGenStop,StopTolerance=EvolAlgStopTol,&
                                     nGenPrint=EvolAlgNGenPrint,File=EvolAlgLogFile,CritType="OptGain")
            GainOpt=Gain
            GainOptScaled=GainScaled
            InbOpt=InbSol
            InbOptRebased=InbSolRebased
            RateInbOpt=RateInbSol
            IndInbOpt=IndInbSol
            RateIndInbOpt=RateIndInbSol
            ValueHoldOpt=ValueHold

            open(newunit=UnitContri,file="AlphaMateResults"//DASH//"ContribAndMatingNbPerIndivOptimumGain.txt",status="unknown")
            !                           12345678901   12345678901   12345678901   12345678901   12345678901
            write(UnitContri,"(5a11)") "         Id","     OrigId","     Gender"," Contribute","   nMatings"
            do i=1,nInd
                write(UnitContri,"(i11,a11,i11,f11.4,i11)") i,trim(IdC(i)),Gender(i),xVec(i),nVec(i)
            enddo
            close(UnitContri)

            open(newunit=UnitMating,file="AlphaMateResults"//DASH//"MatingListOptimumGain.txt",status="unknown")
            !                           12345678901   12345678901
            write(UnitMating,"(2a11)") "    Parent1","    Parent2"
            do i=1,nMat
                write(UnitMating,"(2i11)") i,Mate(i,:)
            enddo
            close(UnitMating)

            ! --- Evaluate the full frontier ---

            if (EvaluateFrontier) then

                write(stdout,"(a)") "Evaluate the full frontier (this might take quite some time!):"
                write(stdout,"(a)") " "

                open(newunit=UnitFrontier,file="AlphaMateResults"//DASH//"Frontier.txt",status="unknown")
                !                             12345678901   12345678901   12345678901   12345678901   12345678901   12345678901   12345678901   12345678901
                write(UnitFrontier,"(8a11)") "       Step","       Gain"," GainScaled"," PopInbreed"," RatePopInb"," IndInbreed"," RateIndInb","  Objective"
                write(UnitFrontier,"(i11,7f11.4)") 1,GainMinInb,GainMinInbScaled,InbMinInb,RateInbMinInb,IndInbMinInb,RateIndInbMinInb,ValueHoldMinInb
                write(UnitFrontier,"(i11,7f11.4)") 2,GainOpt,   GainOptScaled,   InbOpt,   RateInbOpt,   IndInbOpt,   RateIndInbOpt,   ValueHoldOpt

                ! Hold old results
                InbTargetRebasedHold=InbTargetRebased
                InbTargetHold=InbTarget
                RateInbTargetHold=RateInbTarget
                RateInbTarget=RateInbMinInb

                ! Evaluate
                do i=1,nFrontierSteps
                    RateInbTarget=RateInbFrontier(i)
                    InbTargetRebased=RateInbTarget ! due to rebasing F=RateInb
                    InbTarget=InbOld*(1.0d0-RateInbTarget)+RateInbTarget
                    write(stdout,"(a,i3,a,i3,a,f7.4)") "Step ",i," out of ",nFrontierSteps, " for the rate of inbreeding of",RateInbTarget
                    write(stdout,"(a)") ""
                    EvolAlgLogFile="AlphaMateResults"//DASH//"OptimisationLog"//Int2Char(i)//".txt"
                    nTmp=nPotPar1+nPotPar2+nMat ! TODO: add PAGE dimension
                    call EvolAlgForAlphaMate(nParam=nTmp,nSol=EvolAlgNSol,nGen=EvolAlgNGen,nGenBurnIn=EvolAlgNGenBurnIn,&
                                             nGenStop=EvolAlgNGenStop,StopTolerance=EvolAlgStopTol,&
                                             nGenPrint=EvolAlgNGenPrint,File=EvolAlgLogFile,CritType="OptGain")
                    write(UnitFrontier,"(i11,7f11.4)") i+2,Gain,GainScaled,InbSol,RateInbSol,IndInbSol,RateIndInbSol,ValueHold
                    if ((RateInbTarget-RateInbSol) > 0.01d0) then
                        write(stdout,"(a,f)") "NOTE: Could not achieve the rate of inbreeding of ",RateInbTarget
                        write(stdout,"(a,f)") "NOTE: Stopping the frontier evaluation."
                        write(stdout,"(a)") ""
                        exit
                    endif
                enddo

                ! Put back old results
                InbTargetRebased=InbTargetRebasedHold
                InbTarget=InbTargetHold
                RateInbTarget=RateInbTargetHold

                close(UnitFrontier)

            endif

            deallocate(BvScaled)
        end subroutine AlphaMateSearch

        !#######################################################################

        subroutine EvolAlgForAlphaMate(nParam,nSol,nGen,nGenBurnIn,nGenStop,&
            StopTolerance,nGenPrint,File,CritType)
            implicit none

            ! Arguments
            integer,intent(in)          :: nParam        ! No. of parameters in a solution
            integer,intent(in)          :: nSol          ! No. of solutions to test each generation
            integer,intent(in)          :: nGen          ! No. of generations to run
            integer,intent(in)          :: nGenBurnIn    ! No. of generations with more
            integer,intent(in)          :: nGenStop      ! Stop after no progress for nGenerationStop
            double precision,intent(in) :: StopTolerance ! Stopping tolerance
            integer,intent(in)          :: nGenPrint     ! Print changed solution every nGenerationPrint
            character(len=*),intent(in) :: File          ! Which file to write to
            character(len=*),intent(in) :: CritType      ! Passed to FixSolMateAndCalcCrit

            ! Other
            integer :: Param,ParamLoc,Sol,Gen,LastGenPrint
            integer :: SolA,SolB,SolC,BestSol,BestSolOld
            integer :: Unit

            double precision :: RanNum,F,FHold,FHigh1,FHigh2,CR,CRLow,CRHigh
            double precision :: BestValue,BestValueOld,BestValueStop,AcceptRate
            double precision,allocatable :: OldChrom(:,:),NewChrom(:,:),Chrom(:),Value(:)!,MiVal(:),MaVal(:)

            logical :: DiffOnly,BestSolChanged

            LastGenPrint=0
            BestSolOld=0
            BestValueOld=-999999.0d0
            BestValueStop=BestValueOld

            allocate(OldChrom(nParam,nSol))
            allocate(NewChrom(nParam,nSol))
            allocate(Chrom(nParam))
            allocate(Value(nSol))
            ! allocate(MiVal(nParam))
            ! allocate(MaVal(nParam))

            ! --- Printout ---

            ! TODO: make a subroutine for this to make evol alg code generic?
            open(newunit=Unit,file=trim(File),status="unknown")
            !                        12345678901   12345678901   12345678901   12345678901   12345678901   12345678901   12345678901   12345678901   12345678901
            write(stdout,"(9a11)") " SearchMode","       Step","       Gain"," PopInbreed"," RatePopInb"," IndInbreed"," RateIndInb","  Criterion"," AcceptRate"
            write(Unit,  "(9a11)") " SearchMode","       Step","       Gain"," PopInbreed"," RatePopInb"," IndInbreed"," RateIndInb","  Criterion"," AcceptRate"

            ! --- Set parameters ---

            ! TODO: make arguments for this to make evol alg code generic?
            ! Crossover rate
            CRHigh=0.4d0 ! For first few generations (burn-in)
            CRLow=0.2d0  ! For later climbs

            ! F is multiplier of difference used to mutate
            ! Typically between 0.2 and 2.0
            ! (if alleles should be integer, keep F as integer)
            FHold=0.4d0  ! Conservative moves
            FHigh1=1.0d0 ! Adventurous moves
            FHigh2=4.0d0 ! Adventurous moves

            ! Constrain parameters
            ! MiVal=0.0d0
            ! MaVal=1.0d0

            ! --- Initialise foundation population of solutions ---

            ! A solution with equal contributions ! TODO: make arguments for this to make evol alg code generic?
            !OldChrom(:,1)=1.0d0 ! TODO: how should this be setup for the mate selection driver?
            !Value(1)=FixSolMateAndCalcCrit(nParam,OldChrom(:,1),CritType)

            ! Solutions with varied contributions
            !do Sol=2,nSol
            do Sol=1,nSol
                do Param=1,nParam                   ! Set a wide range for each parameter
                    call random_number(RanNum)      ! May need integer values for some problems
                    OldChrom(Param,Sol)=RanNum*(dble(nMat)/2.0d0) ! TODO: what would a good range be?
                enddo
                Value(Sol)=FixSolMateAndCalcCrit(nParam,OldChrom(:,Sol),CritType)
            enddo

            ! --- Evolve ---

            do Gen=1,nGen

                ! Vary differential and non-differential mutation to escape valleys
                if (mod(Gen,3) == 0) then
                    DiffOnly=.true.
                else
                    DiffOnly=.false.
                endif

                ! Burn-in
                if (Gen < nGenBurnIn) then
                    CR=CRHigh
                else
                    CR=CRLow
                endif

                ! Vary mutation rate every few generations
                if (mod(Gen,4) == 0) then
                    F=FHigh1
                else
                    F=FHold
                endif

                if (mod(Gen,7) == 0) then
                    F=FHigh2
                else
                    F=FHold
                endif

                ! --- Generate competitors ---

                ! TODO: paralelize this loop? Is it worth it?
                BestSolChanged=.false.
                AcceptRate=0.0d0
                do Sol=1,nSol

                    ! --- Mutate and recombine ---

                    ! Get three different solutions
                    SolA=Sol
                    do while (SolA == Sol)
                        call random_number(RanNum)
                        SolA=int(RanNum*nSol)+1
                    enddo
                    SolB=Sol
                    do while ((SolB == Sol) .or. (SolB == SolA))
                        call random_number(RanNum)
                        SolB=int(RanNum*nSol)+1
                    enddo
                    SolC=Sol
                    do while ((SolC == Sol) .or. (SolC == SolA) .or. (SolC == SolB))
                        call random_number(RanNum)
                        SolC=int(RanNum*nSol)+1
                    enddo

                    ! Mate the solutions
                    call random_number(RanNum)
                    Param=int(RanNum*nParam)+1 ! Cycle through parameters starting at a random point
                    do ParamLoc=1,nParam
                        call random_number(RanNum)
                        if ((RanNum < CR) .or. (ParamLoc == nParam)) then
                            ! Recombine
                            call random_number(RanNum)
                            if ((RanNum < 0.8d0) .or. DiffOnly) then
                                ! Differential mutation (with prob 0.8 or 1)
                                Chrom(Param)=OldChrom(Param,SolC) + F*(OldChrom(Param,SolA)-OldChrom(Param,SolB))
                            else
                                ! Non-differential mutation (to avoid getting stuck)
                                call random_number(RanNum)
                                if (RanNum < 0.5d0) then
                                    call random_number(RanNum)
                                    Chrom(Param)=OldChrom(Param,SolC) * (0.9d0 + 0.2d0 * RanNum)
                                else
                                    call random_number(RanNum)
                                    Chrom(Param)=OldChrom(Param,SolC) + 0.01d0 * F * (OldChrom(Param,SolA) + 0.01d0) * (RanNum - 0.5d0)
                                endif
                            endif
                        else
                            ! Do not recombine
                            Chrom(Param)=OldChrom(Param,Sol)
                        endif
                        Param=Param+1
                        if (Param > nParam) Param=Param-nParam
                    enddo

                    ! --- Constrain parameters ---

                    ! TODO: make a subroutine/condition for this to make evol alg code generic?
                    ! - this is preferablly not done as it slows convergence
                    ! - it is preferable for criterion function to handle constraints
                    ! - for AlphaMate GG found that it is advised to keep parameters
                    !   unconstrained and set negative values to zero in criterion
                    ! do Param=1,nParam
                    !    if (Chrom(Param) < MiVal(Param)) Chrom(Param)=MiVal(Param)
                    !    if (Chrom(Param) > MaVal(Param)) Chrom(Param)=MaVal(Param)
                    ! enddo

                    ! --- Evaluate and Select ---

                    ValueHold=FixSolMateAndCalcCrit(nParam,Chrom,CritType) ! Merit of competitor
                    if (ValueHold >= Value(Sol)) then                      ! If competitor is better or equal, keep it
                        NewChrom(:,Sol)=Chrom(:)                           !   ("equal" to force evolution)
                        Value(Sol)=ValueHold
                        AcceptRate=AcceptRate+1.0d0
                    else
                        NewChrom(:,Sol)=OldChrom(:,Sol)                    ! Else keep the old solution
                    endif
                enddo ! nSol

                AcceptRate=AcceptRate/dble(nSol)

                ! --- New parents ---

                ! TODO: could we reduce memory requirement by having only OldChrom and Chrom
                !       and replace solutions within OldChrom directly when a new better solution
                !       is found - would this be too greedy? The logic would be that in-silico
                !       we can already use new solution for "mating" as we do not need to follow
                !       discrete generations logic - make generations continous. But this might
                !       become very GREEDY scheme that gets stuck in local optima. An option to
                !       alternate between the two?
                do Sol=1,nSol
                    OldChrom(:,Sol)=NewChrom(:,Sol)
                enddo

                ! --- Find the best solution in this generation ---

                BestSol=maxloc(Value,dim=1)
                BestValue=Value(BestSol)
                if (BestValue > BestValueOld) BestSolChanged=.true.
                BestSolOld=BestSol
                BestValueOld=BestValue

                ! --- Test if solution is improving to stop early ---

                if (mod(Gen,nGenStop) == 0) then
                    if ((BestValue-BestValueStop) > StopTolerance) then
                        BestValueStop=BestValue
                    else
                        write(stdout,"(a,f,a,i,a)") "NOTE: Evolutionary algorithm did not improve objective for ",StopTolerance, " in the last ",nGenStop," generations. Stopping."
                        exit
                    endif
                endif

                ! --- Monitor ---

                if (BestSolChanged) then
                    if (((Gen - LastGenPrint) >= nGenPrint)) then
                        LastGenPrint=Gen
                        ! TODO: make a subroutine for this to make evol alg code generic?
                        ValueHold=FixSolMateAndCalcCrit(nParam,NewChrom(:,BestSol),CritType)
                        write(stdout,"(a11,i11,7f11.4)") CritType,Gen,Gain,InbSol,RateInbSol,IndInbSol,RateIndInbSol,ValueHold,AcceptRate
                        write(Unit,  "(a11,i11,7f11.4)") CritType,Gen,Gain,InbSol,RateInbSol,IndInbSol,RateIndInbSol,ValueHold,AcceptRate
                    endif
                endif

            enddo ! Gen

            ! --- Evaluate the winner ---

            ! TODO: make a subroutine for this to make evol alg code generic?
            ! TODO: add ind. inbreeding to print-out
            ValueHold=FixSolMateAndCalcCrit(nParam,NewChrom(:,BestSol),CritType)
            write(stdout,"(a11,i11,7f11.4)") CritType,Gen,Gain,InbSol,RateInbSol,IndInbSol,RateIndInbSol,ValueHold,AcceptRate
            write(Unit,  "(a11,i11,7f11.4)") CritType,Gen,Gain,InbSol,RateInbSol,IndInbSol,RateIndInbSol,ValueHold,AcceptRate
            write(stdout,"(a)") " "

            close(Unit)

            deallocate(OldChrom)
            deallocate(NewChrom)
            deallocate(Chrom)
            deallocate(Value)
            ! deallocate(MiVal)
            ! deallocate(MaVal)
        end subroutine EvolAlgForAlphaMate

        !#######################################################################

        function FixSolMateAndCalcCrit(nParam,Sol,CritType) result(Criterion)
            implicit none
            ! Arguments
            integer,intent(in)             :: nParam      ! No. of parameters
            double precision,intent(inout) :: Sol(nParam) ! Solution
            character(len=*),intent(in)    :: CritType    ! Type of criterion (MinInb,OptGain)
            ! Other
            integer :: i,j,k,m,nCumMat,RankSol(nInd),SolInt(nInd),MatPar2(nMat),TmpMin,TmpMax
            double precision :: TmpVec(nInd,1),Criterion

            ! Solution has:
            ! - nPotPar1 individual contributions for "parent1" (males when GenderMatters)
            ! - nPotPar2 individual contributions for "parent2" (females when GenderMatters)
            ! - nMat     rankings of parent1 1:nMat matings to pair with 1:nPotPar2 "parent2" (see bellow)
            ! - nInd     rankings for genome editing

            ! Say we have Sol=(|0,2,0,1|...|2.5,1.5,1.0) then we mate
            ! - male 2 with the first  available female (rank 2.5)
            ! - male 2 with the second available female (rank 1.5)
            ! - male 4 with the third  available female (rank 1.0)

            ! --- Is the solution valid? ---

            ! TODO: ALLOW SELFFERTILISATION - AN OPTION? In mate allocations?
            !       will we have to modify solution if this is not achieved?
            !       progeny inbreeding would help in lowering objective, but that
            !       might not be desired in say plants

            ! TODO: This implementation adds one contribution to "trailing" individuals (have values < 1)
            !       until nMat is reached. Should we add these contributions to the first trailing
            !       ind and then jump back to top ind and start adding contributions to top inds sequentially
            !       until we reach next trailing ind and then we repeat this process?

            ! Assure we have nMat contributions for "parent1"
            ! ... find ranks to find the top values
            !write(stdout,"(a)") " "
            call MrgRnk(Sol(1:nPotPar1),RankSol(1:nPotPar1))
            RankSol(1:nPotPar1)=RankSol(nPotPar1:1:-1) ! MrgRnk ranks small to large
            ! ... accumulate the top values as integers (towards nMat)
            !write(stdout,"("//Int2Char(nPotPar1)//"f6.1)") Sol(RankSol(1:nPotPar1))
            nCumMat=0
            do i=1,nPotPar1
                if (Sol(RankSol(i)) < 1.0d0) then
                    ! ... fix by setting to 1 mating
                    Sol(RankSol(i))=1.0d0
                endif
                nCumMat=nCumMat+nint(Sol(RankSol(i))) ! internally real, externally integer
                if (nCumMat >= nMat) then
                    if (nCumMat > nMat) then
                        ! ... make sure there are exactly nMat contributions
                        Sol(RankSol(i))=Sol(RankSol(i))-dble(nCumMat-nMat)
                    endif
                    ! ... make sure all other values are negative
                    if (i < nPotPar1) then
                        Sol(RankSol((i+1):nPotPar1))=sign(Sol(RankSol((i+1):nPotPar1)),-1.0d0)
                    endif
                    exit
                endif
            enddo
            !write(stdout,"(a2,2i4,"//Int2Char(nPotPar1)//"f6.1)") "1",nCumMat,nMat,Sol(RankSol(1:nPotPar1))
            if (nCumMat < nMat) then
                ! TODO
                write(stderr,"(a)") "TODO: nCumMat < nMat!"
                stop 1
            endif

            !write(stdout,"(a)") " "
            ! Assure we have nMat contributions for "parent2"
            ! ... find ranks to find the top values
            call MrgRnk(Sol((nPotPar1+1):(nPotPar1+nPotPar2)),RankSol(1:nPotPar2))
            RankSol(1:nPotPar2)=RankSol(nPotPar2:1:-1) ! MrgRnk ranks small to large
            !write(stdout,"("//Int2Char(nPotPar2)//"f6.1)") Sol(nPotPar1+RankSol(1:nPotPar2))
            ! ... accumulate the top values as integers (towards nMat)
            nCumMat=0
            do i=1,nPotPar2
                if (Sol(nPotPar1+RankSol(i)) < 1.0d0) then
                    ! ... fix by setting to 1 mating
                    Sol(nPotPar1+RankSol(i))=1.0d0
                endif
                nCumMat=nCumMat+nint(Sol(nPotPar1+RankSol(i))) ! internally real, externally integer
                if (nCumMat >= nMat) then
                    if (nCumMat > nMat) then
                        ! ... make sure there are only nMat contributions
                        Sol(nPotPar1+RankSol(i))=Sol(nPotPar1+RankSol(i))-dble(nCumMat-nMat)
                    endif
                    ! ... make sure all other values are negative
                    if (i < nPotPar2) then
                        Sol(nPotPar1+(RankSol((i+1):nPotPar2)))=sign(Sol(nPotPar1+(RankSol((i+1):nPotPar2))),-1.0d0)
                    endif
                    exit
                endif
            enddo
            !write(stdout,"(a2,2i4,"//Int2Char(nPotPar2)//"f6.1)") "2",nCumMat,nMat,Sol(nPotPar1+RankSol(1:nPotPar2))
            if (nCumMat < nMat) then
                ! TODO
                write(stderr,"(a)") "TODO: nCumMat < nMat!"
                stop 1
            endif

            ! TODO: how do we equalize contributions? - distribute equal contributions to nPar positions
            ! TODO: how do we limit the number of parents? - take the top nPar (or nMalPar and nFemPar) positions
            !       (will have to add missing contributions only to these parents!!!)
            ! TODO: how to limit min and max usage

            ! do i=1,nInd
            !    write(stdout,"(i,3f)") i,Sol(i),Sol(nPotPar1+i),Sol(nPotPar1+nPotPar2+i)
            ! enddo

            ! --- Genetic contributions (nVec) ---

            ! "Parent1"
            ! ... get integer values
            SolInt(1:nPotPar1)=nint(Sol(1:nPotPar1))
            ! ... remove negatives
            do i=1,nPotPar1
                if (SolInt(i) < 0) then
                    SolInt(i)=0
                endif
            enddo
            ! ... map internal to external order
            nVecPar1(:)=SolInt(1:nPotPar1)
            if (.not.GenderMatters) then
                nVec(:)=nVecPar1(:)
            else
                nVec(IdPotPar1)=nVecPar1(:)
            endif

            ! "Parent2"
            ! ... get integer values
            SolInt(1:nPotPar2)=nint(Sol((nPotPar1+1):(nPotPar1+nPotPar2)))
            ! ... remove negatives
            do i=1,nPotPar2
                if (SolInt(i) < 0) then
                    SolInt(i)=0
                endif
            enddo
            ! ... map internal to external order
            nVecPar2(:)=SolInt(1:nPotPar2)
            if (.not.GenderMatters) then
                nVec(:)=nVec(:)+nVecPar2(:)
            else
                nVec(IdPotPar2)=nVecPar2(:)
            endif

            xVec(:)=dble(nVec(:))/(2.0d0*dble(nMat))

            ! if (.not.GenderMatters) then
            !     if (nPar == nInd) then
            !         ! Take all individuals
            !         nVec(:)=nint(SolIn(1:nInd))!/sum(SolIn(1:nInd))
            !     else ! nPar < nInd
            !         ! Take the top nPar individuals
            !         SolInSorted(:)=SolIn(1:nInd)
            !         nTmp=nInd
            !         call SortQQ(Loc(SolInSorted),nTmp,SRT$REAL8) ! https://software.intel.com/en-us/node/526803 (2016-02-15)
            !         if (nTmp < nInd) then
            !             write(stderr,"(a)") "Sorting failed"
            !             stop 1
            !         endif
            !         Threshold=SolInSorted(nInd-nPar+1) ! SolInSorted is sorted small to large
            !         !ThresholdSum=sum(SolInSorted((nInd-nPar+1):nInd))
            !         do i=1,nInd
            !             if (SolIn(i) < Threshold) then
            !                 nVec(i)=0!.0d0
            !             else
            !                 nVec(i)=nint(SolIn(i))!/ThresholdSum
            !             endif
            !         enddo
            !     endif
            ! else
            !     if (nMalPar == nMal) then
            !         ! Take all males
            !         if (EqualizeMales) then
            !             nVec(IdPotPar1)=0.5d0/dble(nMal)
            !         else
            !             nVec(IdPotPar1)=0.5d0*SolIn(IdPotPar1)/sum(SolIn(IdPotPar1))
            !         endif
            !     else ! nMalPar < nMal
            !         ! Take the top nMalPar males
            !         SolInSorted(1:nMal)=SolIn(IdPotPar1)
            !         nTmp=nMal
            !         call SortQQ(Loc(SolInSorted(1:nMal)),nTmp,SRT$REAL8) ! https://software.intel.com/en-us/node/526803 (2016-02-15)
            !         if (nTmp < nMal) then
            !             write(stderr,"(a)") "Sorting failed"
            !             stop 1
            !         endif
            !         Threshold=SolInSorted(nMal-nMalPar+1) ! SolInSorted is sorted small to large
            !         ThresholdSum=sum(SolInSorted((nMal-nMalPar+1):nMal))
            !         do i=1,nMal
            !             if (SolIn(IdPotPar1(i)) < Threshold) then
            !                 nVec(IdPotPar1(i))=0.0d0
            !             else
            !                 if (EqualizeMales) then
            !                     nVec(IdPotPar1(i))=0.5d0/dble(nMalPar)
            !                 else
            !                     nVec(IdPotPar1(i))=0.5d0*SolIn(IdPotPar1(i))/ThresholdSum
            !                 endif
            !             endif
            !         enddo
            !     endif
            !     if (nFemPar == nFem) then
            !         ! Take all females
            !         if (EqualizeFemales) then
            !             nVec(IdPotPar2)=0.5d0/dble(nFem)
            !         else
            !             nVec(IdPotPar2)=0.5d0*SolIn(IdPotPar2)/sum(SolIn(IdPotPar2))
            !         endif
            !     else ! nFemPar < nMal
            !         ! Take the top nFemPar females
            !         SolInSorted(1:nFem)=SolIn(IdPotPar2)
            !         nTmp=nFem
            !         call SortQQ(Loc(SolInSorted(1:nFem)),nTmp,SRT$REAL8) ! https://software.intel.com/en-us/node/526803 (2016-02-15)
            !         if (nTmp < nFem) then
            !             write(stderr,"(a)") "Sorting failed"
            !             stop 1
            !         endif
            !         Threshold=SolInSorted(nFem-nFemPar+1) ! SolInSorted is sorted small to large
            !         ThresholdSum=sum(SolInSorted((nFem-nFemPar+1):nFem))
            !         do i=1,nFem
            !             if (SolIn(IdPotPar2(i)) < Threshold) then
            !                 nVec(IdPotPar2(i))=0.0d0
            !             else
            !                 if (EqualizeFemales) then
            !                     nVec(IdPotPar2(i))=0.5d0/dble(nFemPar)
            !                 else
            !                     nVec(IdPotPar2(i))=0.5d0*SolIn(IdPotPar2(i))/ThresholdSum
            !                 endif
            !             endif
            !         enddo
            !     endif
            ! endif
            ! nVec(:)=nVec(:)/(2.0d0*dble(nMat))

            ! --- Mate allocation ---

            ! Distribute parent2 (=females) into matings
            k=0
            do i=1,nPotPar2
                do j=1,nVecPar2(i)
                    k=k+1
                    MatPar2(k)=i
                enddo
            enddo

            ! Reorder them according to the rank of matings
            call MrgRnk(Sol((nPotPar1+nPotPar2+1):(nPotPar1+nPotPar2+nMat)),RankSol(1:nMat))
            RankSol(1:nMat)=RankSol(nMat:1:-1) ! MrgRnk ranks small to large
            MatPar2(:)=MatPar2(RankSol(1:nMat))

            ! Pair parent1 (=males) and parent2 (=females)
            m=0
            Mate(:,:)=0
            do i=1,nPotPar1
                do j=1,nVecPar1(i)
                    m=m+1
                    Mate(m,1)=i
                    Mate(m,2)=MatPar2(m)
                enddo
            enddo

            ! print*,""
            ! do i=1,nMat
            !     write(stdout,"(3i20)") i,Mate(i,:)
            ! enddo
            ! do i=1,nPotPar1
            !     write(stdout,"(2i6)") i,nVecPar1(i)
            ! enddo
            ! do i=1,nPotPar2
            !     write(stdout,"(2i6)") i,nVecPar2(i)
            ! enddo

            ! --- Group coancestry (=future population inbreeding) ---

            ! xA
            do i=1,nInd
                TmpVec(i,1)=dot_product(xVec,RelMtx(:,i))
            enddo
            ! xAx
            InbSol=0.5d0*dot_product(TmpVec(:,1),xVec)
            ! print*,xVec,TmpVec,InbSol

            ! Matrix multiplication with symmetric matrix using BLAS routine
            ! (it was ~5x slower than the above with 1.000 individuals, might be
            !  benefical with larger cases so kept in as commented.)
            ! http://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3.html#ga253c8edb8b21d1b5b1783725c2a6b692
            ! Ax
            ! call dsymm(side="l",uplo="l",m=nInd,n=1,alpha=1.0d0,A=RelMtx,lda=nInd,b=xVec,ldb=nInd,beta=0,c=TmpVec,ldc=nInd)
            ! call dsymm(     "l",     "l",  nInd,  1,      1.0d0,  RelMtx,    nInd,  xVec,    nInd,     0,  TmpVec,    nInd)
            ! xAx
            ! InbSol=0.5d0*dot_product(xVec,TmpVec(:,1))
            ! print*,xVec,TmpVec,InbSol
            ! stop 1

            InbSolRebased=(InbSol-InbOld)/(1.0d0-InbOld)
            RateInbSol=InbSolRebased

            ! --- Progeny inbreeding ---

            IndInbSol=0.0d0
            do i=1,nMat
                ! Try to speed-up retrieval by targeting lower triangle
                TmpMin=minval([Mate(i,1),Mate(i,2)])
                TmpMax=maxval([Mate(i,1),Mate(i,2)])
                IndInbSol=IndInbSol+0.5d0*RelMtx(TmpMax,TmpMin)
            enddo
            IndInbSol=IndInbSol/dble(nMat)
            IndInbSolRebased=(IndInbSol-InbOld)/(1.0d0-InbOld)
            !if (IndInbSolRebased < 0.0) then
            !    IndInbSolRebased=0.0d0
            !endif
            RateIndInbSol=IndInbSolRebased
            IndInbTargetRebased=RateIndInbTarget

            ! --- Genetic gain ---

            Gain=dot_product(xVec,Bv)
            GainScaled=dot_product(xVec,BvScaled)

            ! --- Criterion ---

            Criterion=0.0d0

            ! Gain
            if (trim(CritType) == "OptGain") then
                Criterion=Criterion+(GainScaled-GainMinInbScaled)
            endif

            ! Population inbreeding
            if (RateInbSol <= RateInbTarget) then
                Criterion=Criterion-(InbSolRebased-InbTargetRebased)
            else
                ! Say IndInbPenalty=100, this would be a penalty of 1 SD for 0.01
                ! increase in inbreeding above the target. Note that gain is scaled
                ! and inbreeding is rebased so this should be a "stable" soft constraint.
                Criterion=Criterion-PopInbPenalty*(InbSolRebased-InbTargetRebased)
            endif

            ! Individual inbreeding
            if (RateIndInbSol <= RateIndInbTarget) then
                Criterion=Criterion-(IndInbSolRebased-IndInbTargetRebased)
            else
                ! Say IndInbPenalty=100, this would be a penalty of 1 SD for 0.01
                ! increase in inbreeding above the target. Note that gain is scaled
                ! and inbreeding is rebased so this should be a "stable" soft constraint.
                Criterion=Criterion-IndInbPenalty*(IndInbSolRebased-IndInbTargetRebased)
            endif

            return
        end function FixSolMateAndCalcCrit

        !#######################################################################
end module AlphaMateModule

!###############################################################################

program AlphaMate
    use AlphaMateModule

    implicit none
    real :: Start,Finish
    logical :: Success

    call cpu_time(Start)
    call AlphaMateTitles

    Success=SystemQQ(RMDIR//" AlphaMateResults")
    if (.not.Success) then
        write(stderr,"(a)") "ERROR: Failure to remove the old output folder!"
        stop 1
    endif
    Success=SystemQQ(MKDIR//" AlphaMateResults")
    if (.not.Success) then
        write(stderr,"(a)") "ERROR: Failure to make output folder!"
        stop 1
    endif

    call ReadSpecAndDataForAlphaMate
    call SetInbreedingParameters
    call AlphaMateSearch
    call cpu_time(Finish)

    write(stdout,"(a,f20.4,a)") "Time duration of AlphaMate: ",Finish-Start," seconds"
    write(stdout,"(a)") " "
end program AlphaMate

!###############################################################################