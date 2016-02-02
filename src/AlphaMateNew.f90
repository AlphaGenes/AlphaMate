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
#DEFINE MD "mkdir"
#DEFINE RMDIR "rm -r"
#DEFINE RM "rm"
#DEFINE RENAME "mv"

#else
#DEFINE DASH "\"
#DEFINE COPY "copy"
#DEFINE MD "md"
#DEFINE RMDIR "RMDIR /S"
#DEFINE RM "del"
#DEFINE RENAME "MOVE /Y"
#endif

!################################################################################################################################################################

program AlphaMate
    implicit none

    real :: start, finish
    character(len=1000) :: ebvFile, GFile, genderFile
    double precision,allocatable,dimension(:) :: ebv
    double precision,allocatable,dimension(:,:) :: A
    logical :: genderMatters
    integer,allocatable,dimension(:) :: idMale, idFemale, newSire, newDam
    integer :: maxIter=40, maxGen, popSize, nInd, nMatings, nMale, nFemale, idum
    double precision :: startK=40, multK=2, divK=10, dCoA, tol=0.001, genAlgTol=0.001


    call cpu_time(start)
    call Titles

    call ReadSpecFile(GFile, ebvFile, genderMatters, genderFile, nMale, nFemale, &
                      nInd, nMatings, dCoA, popSize, maxGen)
    allocate(ebv(nInd))
    call ReadEbv(ebvFile, ebv)
    allocate(A(nInd,nInd))
    call ReadA(AFile, A)

    if(genderMatters) then
        allocate(idMale(nMale))
        allocate(idFemale(nFemale))
        call ReadGender(genderFile, idMale, idFemale)
    endif

    call InitiateSeed(idum)

    call RunAlphaMate(ebv, G, dCoA, genderMatters, idMale, idFemale, nMatings, startK, multK, &
                      divK, tol, maxIter, popSize, maxGen, genAlgTol, newSire, newDam, idum)

    call cpu_time(finish)
    print *," "
    print '("  Time duration of AlphaMate = ",f20.4," seconds.")',finish-start
    print *," "

end program AlphaMate

!################################################################################################################################################################

subroutine ReadSpecFile(AFile, ebvFile, genderMatters, genderFile, nMale, nFemale, nInd, &
                        nMatings, dCoA, popSize, maxGen)
    implicit none

    character(len=1000), INTENT(OUT) :: ebvFile, AFile, genderFile
    logical, INTENT(OUT) :: genderMatters
    integer, INTENT(OUT) :: nMale=0, nFemale=0, nInd=0, nMatings=0, popSize=0, maxGen=0
    double precision, INTENT(OUT) :: dCoA
    integer :: nSpecLines, nIndTmp
    character(len=3) :: genderYesNo
    character(len=300) :: dumC

    open (unit=1,file="AlphaMateSpec.txt",status="old")
    call CountLines("AlphaMateSpec.txt",nSpecLines)

    read (1,*) dumC,AFile
    read (1,*) dumC,ebvFile

    call CountLines(AFile,nInd)
    call CountLines(ebvFile,nIndTmp)

    if (nIndTmp/=nInd) then
        print*, "Number of individuals in EBV file and A Matrix file are not the same"
        stop
    endif

    read (1,*) dumC, genderYesNo
    if (trim(genderYesNo)=='Yes') then
        genderMatters=.true.
        read (1,*) dumC,genderFile
        read (1,*) dumC,nMale
        read (1,*) dumC,nFemale
        nInd=nMale+nFemale
    else if(trim(genderYesNo)=='No') then
        genderMatters=.false.
        read (1,*) dumC,nInd
    else
        print*, "PerformMatingConsideringGender must be 'Yes' or 'No'"
        stop
    endif

    if (nIndTmp/=nInd) then
        print*, "Number of individuals in SpecFile doesn't match lines in EBV file"
        stop
    endif

    read (1,*) dumC,nMatings

    read (1,*) dumC,dCoA
    read (1,*) dumC,popSize
    read (1,*) dumC,maxGen

    close(1)

    return

end subroutine ReadSpecFile

!################################################################################################################################################################

subroutine ReadEbv(ebvFile, ebv)
    implicit none

    character(len=1000),intent(IN) :: ebvFile
    integer :: nInd, i
    character(len=1000) :: idCTmp
    double precision,allocatable,intent(INOUT) :: ebv(:)

    if(allocated(ebv)) then
        nInd = ubound(ebv, dim=1)
    else
        print*, "ebv not allocated"
        stop
    endif

    open (unit=1,file=trim(ebvFile),status="old")

    do i=1,nInd
        read (1,*) idCTmp,ebv(i)
    enddo

    close(1)

    return

end subroutine ReadEbv

!################################################################################################################################################################

subroutine ReadA(AFile, A)
    implicit none

    character(len=1000),intent(IN) :: AFile
    integer :: nInd, i
    character(len=1000) :: idCTmp
    double precision,allocatable,intent(INOUT) :: A(:,:)

    if(allocated(A)) then
        nInd = ubound(A, dim=1)
    else
        print*, "A not allocated"
        stop
    endif

    open (unit=1,file=trim(ebvFile),status="old")

    do i=1,nInd
        read (1,*) idCTmp,A(i,:)
    enddo

    close(1)

    return

end subroutine ReadA

!################################################################################################################################################################

subroutine ReadGender(genderFile, idMale, idFemale)
    implicit none

    character(len=1000),intent(IN) :: genderFile
    integer :: nInd, nMale, nFemale, i, tmp, jMale=0, jFemale=0
    character(len=1000) :: idCTmp
    integer,allocatable,intent(INOUT) :: idMale(:),idFemale(:)

    if(allocated(idMale)) then
        nMale = ubound(idMale, dim=1)
    else
        print*, "idMale not allocated"
        stop
    endif

    if(allocated(idFemale)) then
        nFemale = ubound(idFemale, dim=1)
    else
        print*, "idFemale not allocated"
        stop
    endif

    nInd = nMale + nFemale

    open (unit=1,file=trim(genderFile),status="old")


    do i=1,nInd
        read (1,*) idCTmp,genderTmp
        if(genderTmp==1) then
            jMale = jMale+1
            idMale(jMale) = i
        else if(genderTmp==2) then
            jFemale = jFemale+1
            idFemale(jFemale) = i
        else
            print*, "Gender must be 1 for male and 2 for female"
            stop
        endif
    enddo

    return

end subroutine ReadGender

!################################################################################################################################################################

subroutine RunAlphaMate(ebv, G, dCoA, genderMatters, idMale, idFemale, nMatings, startK, multK, &
                        divK, tol, maxIter, popSize, maxGen, genAlgTol, newSire, newDam, idum)
    implicit none

    double precision,allocatable,dimension(:),intent(IN) :: ebv
    double precision,allocatable,dimension(:,:),intent(IN) :: G
    logical,intent(IN) :: genderMatters
    integer,allocatable,dimension(:),intent(IN) :: idMale, idFemale
    integer,intent(IN) :: maxIter, maxGen, popSize, nMatings
    double precision,intent(IN) :: startK, multK, divK, dCoA, tol, genAlgTol
    integer,intent(INOUT) :: idum
    integer,allocatable,dimension(:),intent(OUT) :: newSire, newDam
    integer :: nInd,i,k,ExitCount
    logical :: success=.FALSE.
    character(len=7) :: mode
    character(len=300) :: specifierForCriterion
    double precision,allocatable,dimension(:) :: x

    !integer :: nInd,idum,CurrentLambdaId,nMatings,nMal,nFem
    !integer, allocatable, dimension(:) :: Gender,IdMal,IdFem,nMatingPerInd
    !integer,allocatable,dimension(:,:) :: Matings
    !double precision :: xa,xAx,DeltaFSpec,DeltaF,CurrentLambdaVal,GlobalMinxAx,CurrentxAx,OldLamdaVal
    !double precision :: InitialLambda,LambdaJump,LambdaExitCount,DeltaFPrecision
    !double precision,allocatable,dimension(:) :: Ebv,XVec
    !double precision,allocatable,dimension(:,:) :: Gmat
    !character(len=300) :: dumC,SpecifierForCriterion,SexualType
    !character(len=300),allocatable,dimension(:) :: IdC

    if(allocated(ebv)) then
        nInd = ubound(ebv, dim=1)
    else
        print*, "ebv not allocated"
        stop
    endif
    allocate(x(nInd))
    allocate(newSire(nMatings))
    allocate(newDam(nMatings))

#ifdef OS_UNIX
    success = systemQQ("mkdir -p AlphaMateResults")
#endif
#ifdef OS_WIN
    success = systemQQ("md AlphaMateResults")
#endif
    if(success == .false.) then
        print *, "Failure to create output directory"
        stop
    endif

#ifdef OS_UNIX
    open (unit=101,file="./AlphaMateResults/Frontier.txt",status="unknown")
    open (unit=102,file="./AlphaMateResults/MatingListZero.txt",status="unknown")
    open (unit=103,file="./AlphaMateResults/ContribAndMatingNbPerIndivZero.txt",status="unknown")
#endif
#ifdef OS_WIN
    open (unit=101,file=".\AlphaMateResults\Frontier.txt",status="unknown")
    open (unit=102,file=".\AlphaMateResults\MatingListZero.txt",status="unknown")
    open (unit=103,file=".\AlphaMateResults\ContribAndMatingNbPerIndivZero.txt",status="unknown")
#endif

    !                          1234567890   1234567890   12345678901234567890  12345678901234567890   12345678901234567890   12345678901234567890   12345678901234567890   12345678901234567890
    write (101,'(2a10,6a20)') "SearchMode","      Step","                x'a","                x'Ax","              DeltaF"," |DeltaF-DeltaFSpec|","              Lambda","     x'a+Lambda*x'Ax"

    mode="Coarse"

    ! --- Find lower bound ---
    ! This gives GlobalMinxAx (=F0, minimum possible inbreeding without optimizing genetic gain)

    specifierForCriterion='Inbreeding'
    call GeneticAlgorithm
    call PerformMating


end subroutine RunAlphaMate

!###########################################################################################################################################################

subroutine GeneticAlgorithm(ebv, A, genderMatters, idMale, idFemale, x)

    implicit none
    double precision,allocatable,dimension(:),intent(IN) :: ebv
    double precision,allocatable,dimension(:,:),intent(IN) :: G
    logical,intent(IN) :: genderMatters
    integer,allocatable,dimension(:),intent(IN) :: idMale, idFemale
    integer,intent(IN) :: maxIter, maxGen, popSize, nMatings
    double precision,intent(IN) :: startK, multK, divK, dCoA, tol, genAlgTol
    integer,intent(INOUT) :: idum
    integer,allocatable,dimension(:),intent(OUT) :: newSire, newDam
    integer :: nInd,i,k,ExitCount
    logical :: success=.FALSE.
    character(len=7) :: mode
    character(len=300) :: specifierForCriterion
    double precision,allocatable,dimension(:) :: x
    integer ::  i,j,k,kold=0,a,b,c,generation,kseed,DiffOnly
    double precision :: F,Fhold,CR,value_hold,rand_num,criterion
    double precision,allocatable :: parent_allele(:,:),progeny_allele(:,:),value(:),allele(:),MiVAL(:),MaVAL(:)
    integer,allocatable :: seed(:)
    logical :: solutionchanged

    ! Set some EA parameters ...
    CR=(loci/2.-1.)/(loci-1.) ! 'Crossover Rate' set about (P/2-1)/(P-1) to get .5 effectively
    Fhold=0.5                 ! F is multiplier of difference used to mutate (Fhold holds original F). Set about 0.4

    !  housekeeping ...
    allocate (parent_allele(popsize,loci))
    allocate (progeny_allele(popsize,loci))
    allocate (allele(loci))
    allocate (value(popsize))
    allocate (MiVAL(loci))
    allocate (MaVAL(loci))

    call RANDOM_SEED(size=kseed)

    allocate (seed(kseed))

    do i=1,kseed
        seed(i)=3
    enddo

    call RANDOM_SEED(put=seed)

    MiVal=0
    MaVal=1

    !initialise foundation population of solutions ...
    do i=1, popsize
        do j=1,loci
            call RANDOM_NUMBER(rand_num) ! uniform random number 0 to 1 range
            parent_allele(i,j)=rand_num  ! Set a wide range here appropriate to each parameter
        enddo                            ! May need to generate integer values depending
        allele(:)=parent_allele(i,:)     ! on the nature of the parameters involved
        value(i)=criterion(loci,allele)
    enddo

    !  run over generations

    do generation=1, max_gens

        !  this bit varies 'mutation' every few generations
        if (MOD(generation,4).eq.0) then
            F=1
        else
            F=Fhold
        endif

        if (MOD(generation,7).eq.0) then
            F=5
        else
            F=Fhold
        endif

        if (MOD(generation,3).eq.0) then
            DiffOnly=1
        else
            DiffOnly=0
        endif

        ! This loop generates a competitor for each animal in the population, and replaces it if better ...
        solutionchanged=.false.
        do i=1,popsize

            ! Mutate and recombine

1           call RANDOM_NUMBER(rand_num)
            a=INT(rand_num*popsize) +1    ! choose a different solution to be a template for competitor
            if(a.eq.i) go to 1
2           call RANDOM_NUMBER(rand_num)
            b=INT(rand_num*popsize) +1    ! choose 2 more solutions.  The difference between these drives mutation
            if(b.eq.i .or. b.eq.a) go to 2
3           call RANDOM_NUMBER(rand_num)
            c=INT(rand_num*popsize) +1    ! b and c are these 2 solutions.
            if(c.eq.i .or. c.eq.a .or. c.eq.b) go to 3   !    a, b and c must differ

            call RANDOM_NUMBER(rand_num)
            j=INT(rand_num*loci) + 1   ! cycle through loci starting at a random point
            do k=1,loci
                call RANDOM_NUMBER(rand_num)
                if (rand_num.lt.CR .or. k.eq.loci) then
                    call RANDOM_NUMBER(rand_num)
                    if (rand_num < 0.8 .Or. DiffOnly == 1) Then
                        ! differential mutation ... If alleles should be integer, keep F as integer.
                        allele(j)=parent_allele(c,j) + F*(parent_allele(a,j)-parent_allele(b,j))
                    else
                        ! non - differential mutation
                        call RANDOM_NUMBER(rand_num)
                        allele(j)=parent_allele(c,j) * (0.9 + 0.2 * rand_num)

                    !         call RANDOM_NUMBER(rand_num)
                    !         allele(j)=allele(j) +  0.01*(rand_num - 0.5)  ! because of danger of getting stuck at zero.  Sensible coefficient would help
                    endif
                else
                    allele(j)=parent_allele(i,j)                        ! or no change depending on CR
                end if
                j=(j+1)
                if(j.gt.loci) j=j-loci
            end do
            do j = 1,loci
                If (allele(j) < MiVal(j)) allele(j) = MiVal(j)
                If (allele(j) > MaVal(j)) allele(j) = MaVal(j)
            enddo

            ! Evaluate and Select */

            ! Ensure alleles are integers here if needed. Eg.  "allele(j)=int(allele(j)+.1)"

            value_hold=criterion(loci,allele)     ! merit of competitor
            if(value_hold>=value(i)) then         ! if competitor is better, replace previous
                do j=1,loci                       ! generation solution 'TITLE HOLDER' for this
                    progeny_allele(i,j)=allele(j) ! position in the population
                end do
                if(i==kold)solutionchanged=.true.
                value(i)=value_hold
            else
                do j=1,loci
                    progeny_allele(i,j)=parent_allele(i,j) ! else the current 'title holder' keeps this position
                end do
            endif
        end do  ! popsize

        ! End of population loop, so breed by swapping arrays ...
        do i=1,popsize
            do j=1,loci
                parent_allele(i,j)=progeny_allele(i,j)
            enddo
        enddo

        ! Find best solution in this generation
        value_hold=-9999999999.0
        do i=1,popsize
            if (value(i).gt.value_hold) then
                k=i                 ! k gives current idea of best soln.
                value_hold=value(i) ! " best value_hold
            endif
        enddo
        if(kold/=k)solutionchanged=.true.
        kold=k

    ! print results every 10 generations...
    ! if (MOD(generation,1).eq.0 .and. solutionchanged) then
    !   write(*,'(i5,2f8.3,f20.13)') generation,(progeny_allele(k,i),i=1,loci),value_hold
    ! endif

    end do ! generation

    do i=1,loci
        allele(i)=progeny_allele(k,i)
    end do

    value_hold=criterion(loci,allele)

    if (trim(SpecifierForCriterion)=='Inbreeding') then
        GlobalMinxAx=xAx
        DeltaF=0
    else
        DeltaF=xAx-GlobalMinxAx
    endif
    write (101,'(a10,i10,6f20.10)') Mode,CurrentLambdaId,xa,xAx,DeltaF,abs(DeltaF-DeltaFSpec),CurrentLambdaVal,value_hold

end subroutine GeneticAlgorithm

!###########################################################################################################################################################

subroutine PerformMating
    use Global
    implicit none

    integer :: i,j,k,l,MatingIdVec(nMatings*2),ShuffleMatingIdVec(nMatings*2),nTmp
    double precision :: SortedXvec(nInd),ran1

    SortedXvec(:)=Xvec(:)
    nMatingPerInd=0
    MatingIdVec=0
    l=0
    do i=1,nInd
        k=maxloc(SortedXvec(:),dim=1)
        SortedXvec(k)=(minval(SortedXvec(:)))-1

        nTmp=int(Xvec(k)*(nMatings*2))
        nMatingPerInd(k)=nTmp
        do j=1,nTmp
            l=l+1 ! l evolves from 1 to nMatings independently of i.
            MatingIdVec(l)=k
            if (l==(nMatings*2)) exit
        enddo
        if (l==(nMatings*2)) exit
    enddo

    if (sum(nMatingPerInd(:))<(nMatings*2)) then ! if all of the remaining matings have not received individuals
        do
            k=int(ran1(idum)*nInd)+1
            if (ran1(idum)<Xvec(k)) then ! ids are selected for these matings according to their contribution.
                l=l+1
                MatingIdVec(l)=k
                nMatingPerInd(k)=nMatingPerInd(k)+1
            endif
            if (l==(nMatings*2)) exit
        enddo
    endif

    call RandomOrder(ShuffleMatingIdVec,(nMatings*2),idum)

    l=0
    do i=1,nMatings
        l=l+1
        Matings(i,1)=MatingIdVec(ShuffleMatingIdVec(l))
        l=l+1
        Matings(i,2)=MatingIdVec(ShuffleMatingIdVec(l))
    enddo

end subroutine PerformMating

!###########################################################################################################################################################

subroutine PerformMatingGenderIncl
    use Global
    implicit none

    integer :: i,j,k,l,nTmp,MatingMalIdVec(nMatings),MatingFemIdVec(nMatings),nMatingPerMal(nMal),nMatingPerFem(nFem)
    integer :: ShuffleMatingMalVec(nMatings),ShuffleMatingFemVec(nMatings)
    double precision :: ran1,MalXvec(nMal),FemXvec(nFem),MalSortedXvec(nMal),FemSortedXvec(nFem)

    ! Select male parents

    MatingMalIdVec=0
    nMatingPerMal=0

    l=0
    do i=1,nMal
        MalXvec(i)=Xvec(IdMal(i))
    enddo

    MalSortedXvec(:)=MalXvec(:)

    do i=1,nMal
        k=maxloc(MalSortedXvec,dim=1) ! k is the position (id) in vector x having the highest contribution.
        MalSortedXvec(k)=(minval(MalSortedXvec))-1

        nTmp=int(MalXvec(k)*(nMatings)) ! changed! Xvec(SortedIdVec(i))=Xvec(k)
        nMatingPerMal(k)=nTmp
        do j=1,nTmp
            l=l+1 ! l evolves from 1 to nMatings independently of i.
            MatingMalIdVec(l)=IdMal(k) ! MatingIdVec takes the value of k, or SortedIdVec(i), at positions 1 to nTmp.
            if (l==(nMatings)) exit
        enddo
        if (l==(nMatings)) exit
    enddo

    if (sum(nMatingPerMal)<(nMatings)) then ! if all of the remaining matings have not received individuals.
        do
            k=int(ran1(idum)*nMal)+1
            if (ran1(idum)<MalXvec(k)) then ! ids are selected for these matings according to their contribution.
                l=l+1
                MatingMalIdVec(l)=IdMal(k)
                nMatingPerMal(k)=nMatingPerMal(k)+1
            endif
            if (l==(nMatings)) exit
        enddo
    endif

    ! Select female parents

    MatingFemIdVec=0
    nMatingPerFem=0

    l=0
    do i=1,nFem
        FemXvec(i)=Xvec(IdFem(i))
    enddo

    FemSortedXvec(:)=FemXvec(:)

    do i=1,nFem
        k=maxloc(FemSortedXvec,dim=1)
        FemSortedXvec(k)=(minval(FemSortedXvec))-1

        nTmp=int(FemXvec(k)*(nMatings))
        nMatingPerFem(k)=nTmp
        do j=1,nTmp
            l=l+1
            MatingFemIdVec(l)=IdFem(k)
            if (l==(nMatings)) exit
        enddo
        if (l==(nMatings)) exit
    enddo

    if (sum(nMatingPerFem)<(nMatings)) then
        do
            k=int(ran1(idum)*nFem)+1
            if (ran1(idum)<FemXvec(k)) then
                l=l+1
                MatingFemIdVec(l)=IdFem(k)
                nMatingPerFem(k)=nMatingPerFem(k)+1
            endif
            if (l==(nMatings)) exit
        enddo
    endif

    ! Perform matings

    call RandomOrder(ShuffleMatingMalVec,nMatings,idum)
    call RandomOrder(ShuffleMatingFemVec,nMatings,idum)

    l=0
    do i=1,nMatings
        l=l+1
        Matings(i,1)=MatingMalIdVec(ShuffleMatingMalVec(l))
        Matings(i,2)=MatingFemIdVec(ShuffleMatingFemVec(l))
    enddo

    ! Fill nMatingPerInd

    nMatingPerInd=0

    do i=1,nMal
        nMatingPerInd(IdMal(i))=nMatingPerMal(i)
    enddo
    do i=1,nFem
        nMatingPerInd(IdFem(i))=nMatingPerFem(i)
    enddo

end subroutine PerformMatingGenderIncl

!###########################################################################################################################################################

subroutine InitiateSeed(idum)
    implicit none
    integer :: edum,idum
    DOUBLE PRECISION :: W(1),GASDEV

    open (unit=3,file="Seed.txt",status="old")
    !READ AND WRITE SEED BACK TO FILE
    READ (3,*) idum
    W(1)=GASDEV(idum)
    !Code to write new seed to file
    IF (idum>=0) THEN
        edum=(-1*idum)
    ELSE
        edum=idum
    END IF
    REWIND (3)
    WRITE (3,*) edum
    close(3)
    idum=edum

end subroutine InitiateSeed

!#############################################################################################################################################################################################################################

subroutine RandomOrder(order,n,idum)
    implicit none

    !     Generate a random ordering of the integers 1 ... n.

    integer, INTENT(IN)  :: n
    integer, INTENT(OUT) :: order(n)
    integer :: idum
    double precision ran1

    !     Local variables

    integer :: i, j, k
    double precision    :: wk

    do i = 1, n
        order(i) = i
    end do

    !     Starting at the end, swap the current last indicator with one
    !     randomly chosen from those preceeding it.

    do i = n, 2, -1
        wk=ran1(idum)
        j = 1 + i * wk
        if (j < i) then
            k = order(i)
            order(i) = order(j)
            order(j) = k
        end if
    end do

    RETURN
end subroutine RandomOrder

!#############################################################################################################################################################################################################################

FUNCTION gasdev(idum)
    IMPLICIT NONE
    !C USES ran1
    !Returns a normally distributed deviate with zero mean and unit variance, using ran1(idum)
    !as the source of uniform deviates.

    INTEGER idum
    DOUBLE PRECISION :: gasdev
    INTEGER iset
    DOUBLE PRECISION fac,gset,rsq,v1,v2,ran1
    SAVE iset,gset
    DATA iset/0/
    if (idum.lt.0) iset=0
    if (iset.eq.0) then
1       v1=2.*ran1(idum)-1.
        v2=2.*ran1(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
    else
        gasdev=gset
        iset=0
    endif
    return
END

!#############################################################################################################################################################################################################################

! This Function returns a uniform random deviate between 0.0 and 1.0.
! Set IDUM to any negative value to initialize or reinitialize the sequence.
!MODIFIED FOR REAL

FUNCTION ran1(idum)
    IMPLICIT NONE
    INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
    DOUBLE PRECISION ran1,AM,EPS,RNMX
    PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
    INTEGER j,k,iv(NTAB),iy
    SAVE iv,iy
    DATA iv /NTAB*0/, iy /0/
    IF (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        DO 11 j=NTAB+8,1,-1
            k=idum/IQ
            idum=IA*(idum-k*IQ)-IR*k
            IF (idum.lt.0) idum=idum+IM
            IF (j.le.NTAB) iv(j)=idum

11      CONTINUE
        iy=iv(1)
    END IF
    k=idum/IQ
    idum=IA*(idum-k*IQ)-IR*k
    IF (idum.lt.0) idum=idum+IM
    j=1+iy/NDIV
    iy=iv(j)
    iv(j)=idum
    ran1=min(AM*iy,RNMX)
    RETURN
END

!  (C) Copr. 1986-92 Numerical Recipes Software 6

!#############################################################################################################################################################################################################################

subroutine CountLines(FileName,nLines)
    implicit none

    character(len=*) :: FileName
    character(len=300) :: dumC
    integer :: nLines,f

    nLines=0
    open (unit=101,file=trim(FileName),status="old")
    do
        read (101,*,iostat=f) dumC
        nLines=nLines+1
        if (f/=0) then
            nLines=nLines-1
            exit
        endif
    enddo
    close(101)

end subroutine CountLines

!#############################################################################################################################################################################################################################

   subroutine Titles
       print*, ""
       write(*,'(a30,a,a30)') " ","**********************"," "
       write(*,'(a30,a,a30)') " ","*                    *"," "
       write(*,'(a30,a,a30)') " ","*     AlphaMate      *"," "
       write(*,'(a30,a,a30)') " ","*                    *"," "
       write(*,'(a30,a,a30)') " ","**********************"
       write(*,'(a30,a,a30)') " ","VERSION:"//TOSTRING(VERS)," "
       print*, "              Software for optimizing contributions to the next generation       "
       print*, ""
       print*, "                                    No Liability              "
       print*, "                          Bugs to John.Hickey@roslin.ed.ac.uk"
       print*, ""
   end subroutine Titles

!#############################################################################################################################################################################################################################
