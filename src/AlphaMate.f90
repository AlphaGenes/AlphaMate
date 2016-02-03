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

module AlphaMateModule
    use ifport  ! required for systemQQ
    implicit none

    integer :: nInd,idum,currentLambdaId,nMatings,nMal,nFem,popSize,maxGen,maxIter=20
    integer, allocatable, dimension(:) :: gender,idMal,idFem,nMatingPerInd
    integer,allocatable,dimension(:,:) :: matings
    double precision :: xa,xAx,dFSpec,dF,currentLambda,xAxMin,xAxCurrent
    double precision :: lambdaDown=10,lambdaUp=10,dFTol=0.001,initialLambda=10 !May add to SpecFile
    double precision :: genAlgTol=0.0001 !Minimum require change in objective
    double precision,allocatable,dimension(:) :: ebv,xVec,ebvS
    double precision,allocatable,dimension(:,:) :: Amat
    character(len=300) :: dumC,criterionType
    character(len=300),allocatable,dimension(:) :: idC
    character(len=7) :: mode
    logical :: genderMatters

CONTAINS

    subroutine DeallocateAlphaMate
        !Intended for AlphaSim
        implicit none
        if(allocated(nMatingPerInd)) deallocate(nMatingPerInd)
        if(allocated(gender)) deallocate(gender)
        if(allocated(idMal)) deallocate(idMal)
        if(allocated(idFem)) deallocate(idFem)
        if(allocated(matings)) deallocate(matings)
        if(allocated(ebv)) deallocate(ebv)
        if(allocated(ebvS)) deallocate(ebvS)
        if(allocated(xVec)) deallocate(xVec)
        if(allocated(Amat)) deallocate(Amat)
        if(allocated(idC)) deallocate(idC)
    end subroutine DeallocateAlphaMate

    !##################################################################################

    function CalcCriterion(nInd,allele) result(criterion)
        implicit none

        integer :: i,jMal,jFem,nInd
        double precision :: criterion,tmpVec(nInd),totMal,totFem,xaS
        double precision, intent (IN) :: allele(nInd)

        if (genderMatters) then
            totMal=sum(allele(1:nMal))
            totFem=sum(allele((nMal+1):nInd))

            jMal=0
            jFem=0
            do i=1,nInd
                if (gender(i)==1) then
                    jMal=jMal+1
                    xVec(i)=allele(jMal)
                else
                    jFem=jFem+1
                    xVec(i)=allele(nMal+jFem)
                endif
            enddo

            do i=1,nInd
                if (gender(i)==1) then
                    xVec(i)=xVec(i)*0.5/totMal
                else
                    xVec(i)=xVec(i)*0.5/totFem
                endif
            enddo
        else
            do i=1,nInd
                xVec(i)=allele(i)/sum(allele(:))
            enddo
        endif

        do i=1,nInd
            tmpVec(i)=dot_product(xVec,Amat(i,:)) ! xA
        enddo
        xAx=dot_product(tmpVec,xVec)              ! xAx
        xAx = xAx/2

        xa=dot_product(xVec,ebv)
        xaS=dot_product(xVec,ebvS) !Scaled genetic gain

        xAxCurrent=xAx
        if (trim(criterionType)=='Both') then
            criterion=(xaS - (currentLambda*xAx))
        endif

        if (trim(criterionType)=='Inbreeding') then
            criterion=-1.0*(xAx)
        endif

        ! can set criterion=-9999 if vector 'allele' is
        ! found to be illegal or breaking a constraint

        return

    end function CalcCriterion

    !##################################################################################

    subroutine SearchLambda
        implicit none

        integer :: i,k
        double precision :: highLambda,lowLambda,stdDev=0,mean,xAxLow,xAxHigh
        double precision :: slope, intercept
        logical :: success

        !Create output folder
#ifdef OS_UNIX
        success = systemQQ("mkdir -p AlphaMateResults")
#endif
#ifdef OS_WIN
        success = systemQQ("md AlphaMateResults")
#endif
        if(success == .false.) print *, "Failure to make output folder"

        !Scale ebv
        allocate(ebvS(nInd))
        mean = sum(ebv)/size(ebv)
        do i=1,nInd
            stdDev = stdDev + (ebv(i) - mean)**2
        enddo
        stdDev = SQRT(stdDev)
        do i=1,nInd
            ebvS(i) = (ebv(i) - mean) / stdDev
        enddo

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
        write (101,'(2a10,6a20)') "SearchMode","      Step","                x'a","              x'Ax/2","              dF","   dFSpec-dF","              lambda","   x'a-lambda*x'Ax/2"

        mode="Coarse"

        ! --- Find lower bound ---
        ! This gives xAxMin (=F0, minimum possible inbreeding without optimizing genetic gain)

        criterionType='Inbreeding'
        currentLambda=0.0
        call GeneticAlgorithm
        xAxMin=xAx

        if (genderMatters) then
            call PerformMatingGenderIncl
        else
            call PerformMating
        endif

        do i=1,nMatings
            write (102,*) matings(i,:)
        enddo
        if (genderMatters) then
            write (103,'(a82)') "  Id          AlphaDropId     gender              Contribution       nMatings"
            do i=1,nInd
                write (103,'(i20,a20,i10,f20.10,i20)') i,trim(idC(i)),gender(i),xVec(i),nMatingPerInd(i) ! "MatingPerIndividual"i0,".txt"
            enddo
        else
            write (103,'(a70)') "  Id         AlphaDropId        Contribution            nMatings"
            do i=1,nInd
                write (103,'(i20,a20,f20.10,i20)') i,trim(idC(i)),xVec(i),nMatingPerInd(i) ! "MatingPerIndividual"i0,".txt"
            enddo
        endif

        close(102)
        close(103)

        ! --- Explore frontier in rough steps ---
        ! Start from the upper bound down to the region of targeted inbreeding (xAxMin+dFSpec)
        k=1
        currentLambdaId=k
        currentLambda=initialLambda
        criterionType='Both'
        call GeneticAlgorithm

        if((dF>dFSpec)) then
            !Step lambda value up
            do i=1,100
                if(i==100) then
                    print*, "Lambda search failed"
                    stop
                endif
                if((dF<=dFSpec)) then
                    exit
                endif
                k=k+1
                currentLambdaId=k
                lowLambda = currentLambda
                xAxLow = xAxCurrent
                highLambda = currentLambda*lambdaUp
                currentLambda = highLambda
                call GeneticAlgorithm
                xAxHigh = xAxCurrent
            enddo
        else
            !Step lambda value down
            do i=1,100
                if(i==100) then
                    print*, "Lambda search failed"
                    stop
                endif
                if((dF>=(dFSpec-dFTol))) then
                    exit
                endif
                k=k+1
                currentLambdaId=k
                lowLambda = currentLambda/lambdaDown
                highLambda = currentLambda
                xAxHigh = xAxCurrent
                currentLambda = lowLambda
                call GeneticAlgorithm
                xAxLow = xAxCurrent
            enddo
        endif

        mode="Fine"
        do i=1,maxIter
            if((dF<dFSpec) .AND. (dF>=(dFSpec-dFTol))) then
                exit
            endif
            k=k+1
            currentLambdaId=k
            if(mod(i,2)==1) then
                slope = (xAxHigh-xAxLow)/(highLambda-lowLambda)
                intercept = xAxLow - lowLambda * slope
                currentLambda = (dFSpec - xAxMin - intercept) / slope
            else
                currentLambda = (highLambda + lowLambda)/2
            endif
            call GeneticAlgorithm
            if(i==maxIter) then
                print*, "Fine search reached maximum iterations"
            endif
            if(dF>dFSpec) then
                lowLambda = currentLambda
                xAxLow = xAxCurrent
            else
                highLambda = currentLambda
                xAxHigh = xAxCurrent
            endif
        enddo

        if (genderMatters) then
           call PerformMatingGenderIncl
        else
           call PerformMating
        endif
#ifdef OS_UNIX
        open (unit=102,file="./AlphaMateResults/MatingListFinal.txt",status="unknown")
        open (unit=103,file="./AlphaMateResults/ContribAndMatingNbPerIndivFinal.txt",status="unknown")
        open (unit=104,file="./AlphaMateResults/LambdaValueFinal.txt",status="unknown")
#endif
#ifdef OS_WIN
        open (unit=102,file=".\AlphaMateResults\MatingListFinal.txt",status="unknown")
        open (unit=103,file=".\AlphaMateResults\ContribAndMatingNbPerIndivFinal.txt",status="unknown")
        open (unit=104,file=".\AlphaMateResults\LambdaValueFinal.txt",status="unknown")
#endif

        do i=1,nMatings
            write (102,*) matings(i,:)
        enddo
        if (genderMatters) then
            write (103,'(a82)') "  Id          AlphaDropId     gender              Contribution       nMatings"
            do i=1,nInd
                write (103,'(i20,a20,i10,f20.10,i20)') i,trim(idC(i)),gender(i),xVec(i),nMatingPerInd(i) ! "MatingPerIndividual"i0,".txt"
            enddo
        else
            write (103,'(a70)') "  Id         AlphaDropId        Contribution            nMatings"
            do i=1,nInd
                write (103,'(i20,a20,f20.10,i20)') i,trim(idC(i)),xVec(i),nMatingPerInd(i) ! "MatingPerIndividual"i0,".txt"
            enddo
        endif
        write (104,*) currentLambda

        close(101)
        close(102)
        close(103)
        close(104)

    end subroutine SearchLambda

    !##################################################################################

    subroutine GeneticAlgorithm
        !Called within SearchLambda
        implicit none

        integer ::  i,j,k,kold=0,a,b,c,generation,kseed,DiffOnly
        double precision :: F,Fhold,CR,value_hold,rand_num,criterion,oldBest
        double precision,allocatable :: parent_allele(:,:),progeny_allele(:,:),value(:),allele(:),MiVAL(:),MaVAL(:)
        integer,allocatable :: seed(:)
        logical :: solutionchanged

        ! Set some EA parameters ...
        CR=(nInd/2.-1.)/(nInd-1.) ! 'Crossover Rate' set about (P/2-1)/(P-1) to get .5 effectively
        Fhold=0.5                 ! F is multiplier of difference used to mutate (Fhold holds original F). Set about 0.4

        !  housekeeping ...
        allocate (parent_allele(popSize,nInd))
        allocate (progeny_allele(popSize,nInd))
        allocate (allele(nInd))
        allocate (value(popSize))
        allocate (MiVAL(nInd))
        allocate (MaVAL(nInd))

        call RANDOM_SEED(size=kseed)

        allocate (seed(kseed))

        do i=1,kseed
            seed(i)=3
        enddo

        call RANDOM_SEED(put=seed)

        MiVal=0
        MaVal=1

        !initialise foundation population of solutions ...
        do j=1,nInd
            parent_allele(1,j)=1/nInd
        enddo
        allele(:)=parent_allele(1,:)
        value(1)=CalcCriterion(nInd,allele)
        do i=2, popSize
            do j=1,nInd
                call RANDOM_NUMBER(rand_num) ! uniform random number 0 to 1 range
                parent_allele(i,j)=rand_num  ! Set a wide range here appropriate to each parameter
            enddo                            ! May need to generate integer values depending
            allele(:)=parent_allele(i,:)     ! on the nature of the parameters involved
            value(i)=CalcCriterion(nInd,allele)
        enddo

        !  run over generations

        do generation=1, maxGen

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
            do i=1,popSize

                ! Mutate and recombine

1               call RANDOM_NUMBER(rand_num)
                a=INT(rand_num*popSize) +1    ! choose a different solution to be a template for competitor
                if(a.eq.i) go to 1
2               call RANDOM_NUMBER(rand_num)
                b=INT(rand_num*popSize) +1    ! choose 2 more solutions.  The difference between these drives mutation
                if(b.eq.i .or. b.eq.a) go to 2
3               call RANDOM_NUMBER(rand_num)
                c=INT(rand_num*popSize) +1    ! b and c are these 2 solutions.
                if(c.eq.i .or. c.eq.a .or. c.eq.b) go to 3   !    a, b and c must differ

                call RANDOM_NUMBER(rand_num)
                j=INT(rand_num*nInd) + 1   ! cycle through nInd starting at a random point
                do k=1,nInd
                    call RANDOM_NUMBER(rand_num)
                    if (rand_num.lt.CR .or. k.eq.nInd) then
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
                    if(j.gt.nInd) j=j-nInd
                end do
                do j = 1,nInd
                    If (allele(j) < MiVal(j)) allele(j) = MiVal(j)
                    If (allele(j) > MaVal(j)) allele(j) = MaVal(j)
                enddo

                ! Evaluate and Select */

                ! Ensure alleles are integers here if needed. Eg.  "allele(j)=int(allele(j)+.1)"

                value_hold=CalcCriterion(nInd,allele)     ! merit of competitor
                if(value_hold>=value(i)) then         ! if competitor is better, replace previous
                    do j=1,nInd                       ! generation solution 'TITLE HOLDER' for this
                        progeny_allele(i,j)=allele(j) ! position in the population
                    end do
                    if(i==kold)solutionchanged=.true.
                    value(i)=value_hold
                else
                    do j=1,nInd
                        progeny_allele(i,j)=parent_allele(i,j) ! else the current 'title holder' keeps this position
                    end do
                endif
            end do  ! popSize

            ! End of population loop, so breed by swapping arrays ...
            do i=1,popSize
                do j=1,nInd
                    parent_allele(i,j)=progeny_allele(i,j)
                enddo
            enddo

            ! Find best solution in this generation
            value_hold=-9999999999.0
            do i=1,popSize
                if (value(i).gt.value_hold) then
                    k=i                 ! k gives current idea of best soln.
                    value_hold=value(i) ! " best value_hold
                endif
            enddo
            if(kold/=k)solutionchanged=.true.
            kold=k

            !Test if solution is improving every 100 generations
            !If it isn't stop optimization
            if(generation==1) then
                oldBest = value_hold
            else if(mod(generation,100)==0) then
                if((value_hold-oldBest)>genAlgTol) then
                    oldBest = value_hold
                else
                    exit
                endif
            endif

        ! print results every 10 generations...
        ! if (MOD(generation,1).eq.0 .and. solutionchanged) then
        !   write(*,'(i5,2f8.3,f20.13)') generation,(progeny_allele(k,i),i=1,nInd),value_hold
        ! endif

        end do ! generation

        do i=1,nInd
            allele(i)=progeny_allele(k,i)
        end do

        value_hold=CalcCriterion(nInd,allele)

        if (trim(criterionType)=='Inbreeding') then
            dF=0
        else
            dF=xAx-xAxMin
        endif
        write (101,'(a10,i10,6f20.10)') mode,currentLambdaId,xa,xAx,dF,(dFSpec-dF),currentLambda,(xa-currentLambda*xAx)

    end subroutine GeneticAlgorithm

    !##################################################################################

    subroutine PerformMating
        !Called within SearchLambda
        implicit none

        integer :: i,j,k,l,MatingIdVec(nMatings*2),ShuffleMatingIdVec(nMatings*2),nTmp
        double precision :: SortedXvec(nInd),ran1

        SortedXvec(:)=xVec(:)
        nMatingPerInd=0
        MatingIdVec=0
        l=0
        do i=1,nInd
            k=maxloc(SortedXvec(:),dim=1)
            SortedXvec(k)=(minval(SortedXvec(:)))-1

            nTmp=int(xVec(k)*(nMatings*2))
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
                if (ran1(idum)<xVec(k)) then ! ids are selected for these matings according to their contribution.
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
            matings(i,1)=MatingIdVec(ShuffleMatingIdVec(l))
            l=l+1
            matings(i,2)=MatingIdVec(ShuffleMatingIdVec(l))
        enddo

    end subroutine PerformMating

    !##################################################################################

    subroutine PerformMatingGenderIncl
        !Called within SearchLambda
        implicit none

        integer :: i,j,k,l,nTmp,MatingMalIdVec(nMatings),MatingFemIdVec(nMatings),nMatingPerMal(nMal),nMatingPerFem(nFem)
        integer :: ShuffleMatingMalVec(nMatings),ShuffleMatingFemVec(nMatings)
        double precision :: ran1,MalXvec(nMal),FemXvec(nFem),MalSortedXvec(nMal),FemSortedXvec(nFem)

        ! Select male parents

        MatingMalIdVec=0
        nMatingPerMal=0

        l=0
        do i=1,nMal
            MalXvec(i)=xVec(idMal(i))
        enddo

        MalSortedXvec(:)=MalXvec(:)

        do i=1,nMal
            k=maxloc(MalSortedXvec,dim=1) ! k is the position (id) in vector x having the highest contribution.
            MalSortedXvec(k)=(minval(MalSortedXvec))-1

            nTmp=int(MalXvec(k)*(nMatings)) ! changed! xVec(SortedIdVec(i))=xVec(k)
            nMatingPerMal(k)=nTmp
            do j=1,nTmp
                l=l+1 ! l evolves from 1 to nMatings independently of i.
                MatingMalIdVec(l)=idMal(k) ! MatingIdVec takes the value of k, or SortedIdVec(i), at positions 1 to nTmp.
                if (l==(nMatings)) exit
            enddo
            if (l==(nMatings)) exit
        enddo

        if (sum(nMatingPerMal)<(nMatings)) then ! if all of the remaining matings have not received individuals.
            do
                k=int(ran1(idum)*nMal)+1
                if (ran1(idum)<MalXvec(k)) then ! ids are selected for these matings according to their contribution.
                    l=l+1
                    MatingMalIdVec(l)=idMal(k)
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
            FemXvec(i)=xVec(idFem(i))
        enddo

        FemSortedXvec(:)=FemXvec(:)

        do i=1,nFem
            k=maxloc(FemSortedXvec,dim=1)
            FemSortedXvec(k)=(minval(FemSortedXvec))-1

            nTmp=int(FemXvec(k)*(nMatings))
            nMatingPerFem(k)=nTmp
            do j=1,nTmp
                l=l+1
                MatingFemIdVec(l)=idFem(k)
                if (l==(nMatings)) exit
            enddo
            if (l==(nMatings)) exit
        enddo

        if (sum(nMatingPerFem)<(nMatings)) then
            do
                k=int(ran1(idum)*nFem)+1
                if (ran1(idum)<FemXvec(k)) then
                    l=l+1
                    MatingFemIdVec(l)=idFem(k)
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
            matings(i,1)=MatingMalIdVec(ShuffleMatingMalVec(l))
            matings(i,2)=MatingFemIdVec(ShuffleMatingFemVec(l))
        enddo

        ! Fill nMatingPerInd

        nMatingPerInd=0

        do i=1,nMal
            nMatingPerInd(idMal(i))=nMatingPerMal(i)
        enddo
        do i=1,nFem
            nMatingPerInd(idFem(i))=nMatingPerFem(i)
        enddo

    end subroutine PerformMatingGenderIncl

end module AlphaMateModule

!################################################################################################################################################################

program AlphaMate
    use AlphaMateModule

    implicit none
    real :: start, finish

    call cpu_time(start)
    call Titles
    call ReadParametersData !Read in SpecFile
    call InitiateSeed(idum)
    call SearchLambda !Runs AlphaMate
    call cpu_time(finish)

    print *," "
    print '("  Time duration of AlphaMate = ",f20.4," seconds.")',finish-start
    print *," "

end program AlphaMate

!################################################################################################################################################################

subroutine ReadParametersData

    use AlphaMateModule
    implicit none

    integer :: i,j,jMal,jFem,nIndTmp,GenderTmp,nSpecLines
    double precision :: EbvTmp
    character(len=1000) :: RelMatFile,EbvFile,GenderFile,IdCTmp,GenderMattersChar


    open (unit=1,file="AlphaMateSpec.txt",status="old")
    call CountLines("AlphaMateSpec.txt",nSpecLines)


    read (1,*) dumC,RelMatFile
    read (1,*) dumC,EbvFile

    call CountLines(RelMatFile,nInd)
    call CountLines(EbvFile,nIndTmp)

    if (nIndTmp/=nInd) then
        print*, "Number of individuals in ebv file and Relationship Matrix file is not the same"
        stop
    endif

    read (1,*) dumC,GenderMattersChar
    if (trim(GenderMattersChar)=='Yes') then
        genderMatters=.true.
        read (1,*) dumC,GenderFile
        read (1,*) dumC,nMal
        read (1,*) dumC,nFem
        nInd=nMal+nFem
    else
        genderMatters=.false.
        read (1,*) dumC,nInd
    endif

    read (1,*) dumC,nMatings

    read (1,*) dumC,dFSpec
    read (1,*) dumC,popSize
    read (1,*) dumC,maxGen

    close (1)

    allocate(matings(nMatings,2))
    allocate(nMatingPerInd(nInd))

    allocate(Amat(nInd,nInd))
    allocate(ebv(nInd))
    allocate(idC(nInd))
    allocate(xVec(nInd))

    open (unit=11,file=trim(RelMatFile),status="old")
    open (unit=12,file=trim(EbvFile),status="old")

    do i=1,nInd
        read (11,*) idC(i),Amat(i,:)
    enddo

    do i=1,nInd
        read (12,*) IdCTmp,EbvTmp
        do j=1,nInd
            if (trim(IdCTmp)==trim(idC(j))) then
                ebv(j)=EbvTmp
                exit
            endif
        enddo
    enddo

    if (genderMatters) then

        allocate(gender(nInd))
        allocate(idMal(nMal))
        allocate(idFem(nFem))

        open (unit=13,file=trim(GenderFile),status="old")

        do i=1,nInd
            read (13,*) IdCTmp,GenderTmp
            do j=1,nInd
                if (trim(IdCTmp)==trim(idC(j))) then
                    gender(j)=GenderTmp
                    exit
                endif
            enddo
        enddo

        jMal=0
        jFem=0
        do i=1,nInd
            if (gender(i)==1) then
                jMal=jMal+1
                idMal(jMal)=i
            else
                jFem=jFem+1
                idFem(jFem)=i
            endif
        enddo

    endif

    close (11)
    close (12)
    close (13)

end subroutine ReadParametersData

!###########################################################################################################################################################

subroutine InitiateSeed(idum)
    !subroutine by John Hickey October 2009
    implicit none
    integer :: idum,edum
    double precision :: W(1),ran1

    open(unit=3,file="Seed.txt",status="old")
    open(unit=4,file="SeedOld.txt",status="unknown")
    !READ AND write SEED BACK TO file

    read(3,*) idum
    write(4,*) idum

    W(1)=ran1(idum)
    !Code to write new seed to file
    if (idum>=0) then
        edum=(-1*idum)
    ELSE
        edum=idum
    end if
    REWIND (3)
    write(3,*) edum
    idum=edum

end subroutine InitiateSeed

!#############################################################################################################################################################################################################################

subroutine RandomOrder(order,n,idum)
    !Called within GeneticAlgorithm
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
! Set idum to any negative value to initialize or reinitialize the sequence.
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
    !Called within ReadParametersData
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
