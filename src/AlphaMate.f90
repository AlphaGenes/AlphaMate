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
end module AlphaSuiteModule

!###############################################################################

module AlphaMateModule
    use ifport  ! required for systemQQ

#ifdef f2003
    use,intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                             stdout=>output_unit, &
                                             stderr=>error_unit
#else
#define stdin  5
#define stdout 6
#define stderr 0
#endif

    implicit none

    integer :: nInd,idum,nMatings,nMal,nFem
    integer :: EvolAlgNSol,EvolAlgNGen,EvolAlgNGenBurnIn,EvolAlgNGenStop,EvolAlgNGenPrint
    integer,allocatable,dimension(:) :: Gender,IdMal,IdFem,nMatingPerInd
    integer,allocatable,dimension(:,:) :: Matings

    double precision :: EvolAlgStopTol,Gain,DeltaFTarget,DeltaFCurrent
    double precision :: FOld,FTarget,FTargetRebase,FCurrent,FCurrentRebase
    double precision,allocatable,dimension(:) :: Bv,BvScaled,XVec
    double precision,allocatable,dimension(:,:) :: RelMat

    character(len=300),allocatable,dimension(:) :: IdC

    logical :: GenderMatters,InferFOld,Verbose=.false.

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

        subroutine ReadParametersAndDataForAlphaMate
            use AlphaSuiteModule,only : SetSeed,CountLines
            implicit none

            integer :: DumI,i,j,jMal,jFem,nIndTmp,GenderTmp,UnitSpec,UnitRelMat,UnitBv,UnitGender,Seed

            double precision :: BvTmp

            character(len=1000) :: DumC,RelMatFile,BvFile,GenderFile,SeedFile,IdCTmp

            open(newunit=UnitSpec,file="AlphaMateSpec.txt",status="old")

            read(UnitSpec,*) DumC,RelMatFile
            read(UnitSpec,*) DumC,BvFile

            call CountLines(RelMatFile,nInd)
            call CountLines(BvFile,nIndTmp)

            if (nIndTmp /= nInd) then
                write(stderr,"(a)") "ERROR: Number of individuals in ebv file and Relationship Matrix file is not the same!"
                stop 1
            endif

            read(UnitSpec,*) DumC,DumC
            if (trim(DumC) == "Yes") then
                GenderMatters=.true.
                read(UnitSpec,*) DumC,GenderFile
                read(UnitSpec,*) DumC,nMal
                read(UnitSpec,*) DumC,nFem
                nInd=nMal+nFem
            else
                GenderMatters=.false.
                read(UnitSpec,*) DumC,nInd
            endif

            read(UnitSpec,*) DumC,nMatings

            read(UnitSpec,*) DumC,DumC
            if (trim(DumC) == "Unknown") then
                InferFOld=.true.
            else
                InferFOld=.false.
                backspace(UnitSpec)
                read(UnitSpec,*) DumC,FOld
            endif
            read(UnitSpec,*) DumC,DeltaFTarget

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

            close(UnitSpec)

            allocate(Matings(nMatings,2))
            allocate(nMatingPerInd(nInd))

            allocate(RelMat(nInd,nInd))
            allocate(Bv(nInd))
            allocate(IdC(nInd))
            allocate(XVec(nInd))

            open(newunit=UnitRelMat,file=trim(RelMatFile),status="old")
            open(newunit=UnitBv,file=trim(BvFile),status="old")

            do i=1,nInd
                read(UnitRelMat,*) idC(i),RelMat(:,i)
            enddo
            close(UnitRelMat)

            do i=1,nInd
                read(UnitBv,*) IdCTmp,BvTmp
                do j=1,nInd
                    if (trim(IdCTmp) == trim(idC(j))) then
                        Bv(j)=BvTmp
                        exit
                    endif
                enddo
            enddo
            close(UnitBv)

            allocate(Gender(nInd))
            if (.not.GenderMatters) then
                Gender(:)=0
            else
                allocate(IdMal(nMal))
                allocate(IdFem(nFem))

                open(unit=UnitGender,file=trim(GenderFile),status="old")

                do i=1,nInd
                    read(UnitGender,*) IdCTmp,GenderTmp
                    do j=1,nInd
                        if (trim(IdCTmp) == trim(idC(j))) then
                            Gender(j)=GenderTmp
                            exit
                        endif
                    enddo
                enddo
                close(UnitGender)

                jMal=0
                jFem=0
                do i=1,nInd
                    if (Gender(i) == 1) then
                        jMal=jMal+1
                        IdMal(jMal)=i
                    else
                        jFem=jFem+1
                        IdFem(jFem)=i
                    endif
                enddo

            endif
        end subroutine ReadParametersAndDataForAlphaMate

        !#######################################################################

        subroutine SetInbreedingParameters
            implicit none
            integer :: i,UnitInbree
            double precision :: Tmp

            ! Old inbreeding
            if (InferFOld) then
                FOld=0.0
                do i=1,nInd
                    Tmp=RelMat(i,i)-1.0d0
                    if (Tmp < 0.0) then
                        write(stderr,"(a)") "ERROR: Relationship matrix must have diagonals equal or more than 1.0!"
                        stop 2
                    endif
                    FOld=FOld+Tmp
                enddo
                FOld=FOld/dble(nInd)
            endif

            ! New inbreeding
            FTarget=FOld*(1-DeltaFTarget)+DeltaFTarget
            FTargetRebase=(FTarget-FOld)/(1.0d0-FOld)

            ! Report
            write(stdout,"(a,f)") "Old inbreeding: ",FOld
            write(stdout,"(a,f)") "Targeted rate of inbreeding: ",DeltaFTarget
            write(stdout,"(a,f)") "Targeted inbreeding: ",FTarget

            open(newunit=UnitInbree,file="AlphaMateResults"//DASH//"ConstraintInbreeding.txt",status="unknown")
            write(UnitInbree,"(a,f)") "Old_inbreeding, ",FOld
            write(UnitInbree,"(a,f)") "Targeted_rate_of_inbreeding, ",DeltaFTarget
            write(UnitInbree,"(a,f)") "Targeted_inbreeding, ",FTarget
            close(UnitInbree)
        end subroutine SetInbreedingParameters

        !#######################################################################

        subroutine AlphaMateSearch
            implicit none

            integer :: i,UnitInbree,UnitMating,UnitContri

            double precision :: StdDev,Mean

            character(len=300) :: EvolAlgSearchFile

            ! --- Scale breeding values ---

            allocate(BvScaled(nInd))
            Mean=sum(Bv)/dble(nInd)
            BvScaled(:)=Bv(:)-Mean
            StdDev=sqrt(sum(BvScaled(:)*BvScaled(:))/dble(nInd))
            BvScaled(:)=(Bv(:)-Mean)/StdDev

            ! --- Minimize inbreeding ---

            EvolAlgSearchFile="AlphaMateResults"//DASH//"SearchMinimumInbreeding.txt"
            open(newunit=UnitMating,file="AlphaMateResults"//DASH//"MatingListMinimumInbreeding.txt",status="unknown")
            open(newunit=UnitContri,file="AlphaMateResults"//DASH//"ContribAndMatingNbPerIndivMinimumInbreeding.txt",status="unknown")

            call EvolAlgForAlphaMate(nParam=nInd,nSolution=EvolAlgNSol,nGeneration=EvolAlgNGen,nGenerationBurnIn=EvolAlgNGenBurnIn,&
                                     nGenerationStop=EvolAlgNGenStop,StopTolerance=EvolAlgStopTol,&
                                     nGenerationPrint=EvolAlgNGenPrint,File=EvolAlgSearchFile,CriterionType="MinInb")

            if (GenderMatters) then
                call PerformMatingGenderIncl
            else
                call PerformMating
            endif

            ! TODO header; make pedigree?
            do i=1,nMatings
                write(UnitMating,*) Matings(i,:) ! TODO inefficient
            enddo

            !                           1234567890   1234567890   1234567890   1234567890   1234567890
            write(UnitContri,"(5a10)") "        Id","    OrigId","    Gender","Contribute","  nMatings"
            do i=1,nInd
                write(UnitContri,"(i10,a10,i10,f10.7,i10)") i,trim(IdC(i)),Gender(i),XVec(i),nMatingPerInd(i)
            enddo

            close(UnitMating)
            close(UnitContri)

            ! --- Maximize genetic gain with constraint on inbreeding ---

            if (FCurrent < FOld) then
                write(stdout,"(a)") " "
                write(stdout,"(a)") "NOTE: Old inbreeding is higher than the minimum group coancestry (x'Ax/2) under no selection."
                write(stdout,"(a)") "NOTE: Resetting the old inbreeding to the minimum group coancestry under no selection and"
                write(stdout,"(a)") "NOTE:   recomputing the targeted inbreeding."
                FOld=FCurrent
                FTarget=FOld*(1-DeltaFTarget)+DeltaFTarget
                FTargetRebase=(FTarget-FOld)/(1.0d0-FOld)
                write(stdout,"(a,f)") "Old inbreeding: ",FOld
                write(stdout,"(a,f)") "Targeted rate of inbreeding: ",DeltaFTarget
                write(stdout,"(a,f)") "Targeted inbreeding:",FTarget
                write(stdout,"(a)") " "
            endif

            if (FCurrent > FTarget) then
                write(stderr,"(a)") "ERROR: Targeted inbreeding is lower than the group coancestry (x'Ax/2) under no selection."
                write(stderr,"(a)") "ERROR: Can not optimise!"
                stop 3
            endif

            open(newunit=UnitInbree,file="AlphaMateResults"//DASH//"ConstraintInbreeding.txt",status="old")
            write(UnitInbree,"(a,f)") "Old_inbreeding, ",FOld
            write(UnitInbree,"(a,f)") "Targeted_rate_of_inbreeding, ",DeltaFTarget
            write(UnitInbree,"(a,f)") "Targeted_inbreeding, ",FTarget
            close(UnitInbree)

            EvolAlgSearchFile="AlphaMateResults"//DASH//"SearchMaximumGain.txt"
            open(newunit=UnitMating,file="AlphaMateResults"//DASH//"MatingListMaximumGain.txt",status="unknown")
            open(newunit=UnitContri,file="AlphaMateResults"//DASH//"ContribAndMatingNbPerIndivMaximumGain.txt",status="unknown")

            call EvolAlgForAlphaMate(nParam=nInd,nSolution=EvolAlgNSol,nGeneration=EvolAlgNGen,nGenerationBurnIn=EvolAlgNGenBurnIn,&
                                     nGenerationStop=EvolAlgNGenStop,StopTolerance=EvolAlgStopTol,&
                                     nGenerationPrint=EvolAlgNGenPrint,File=EvolAlgSearchFile,CriterionType="MaxGain")

            if (GenderMatters) then
               call PerformMatingGenderIncl
            else
               call PerformMating
            endif

            do i=1,nMatings
                write(UnitMating,*) Matings(i,:) ! TODO inefficient
            enddo
            !                           1234567890   1234567890   1234567890   1234567890   1234567890
            write(UnitContri,"(5a10)") "        Id","    OrigId","    Gender","Contribute","  nMatings"
            do i=1,nInd
                write(UnitContri,"(i10,a10,i10,f10.7,i10)") i,trim(IdC(i)),Gender(i),XVec(i),nMatingPerInd(i)
            enddo

            close(UnitMating)
            close(UnitContri)

            deallocate(BvScaled)
        end subroutine AlphaMateSearch

        !#######################################################################

        subroutine EvolAlgForAlphaMate(nParam,nSolution,nGeneration,nGenerationBurnIn,&
            nGenerationStop,StopTolerance,nGenerationPrint,File,CriterionType)
            implicit none

            ! Arguments
            integer,intent(in)          :: nParam            ! Number of parameters in a solution
            integer,intent(in)          :: nSolution         ! Number of solutions to test each generation
            integer,intent(in)          :: nGeneration       ! Number of generations to run
            integer,intent(in)          :: nGenerationBurnIn ! Number of generations with more
            integer,intent(in)          :: nGenerationStop   ! Stop after no progress for nGenerationStop
            double precision,intent(in) :: StopTolerance     ! Stopping tolerance
            integer,intent(in)          :: nGenerationPrint  ! Print changed solution every nGenerationPrint
            character(len=*),intent(in) :: File              ! Which file to write to
            character(len=*),intent(in) :: CriterionType     ! Passed to CalcCriterion

            ! Other
            integer :: Param,ParamLoc,Solution,Generation,LastGenerationPrint
            integer :: SolutionA,SolutionB,SolutionC,BestSolution,BestSolutionOld
            integer :: DiffOnly,Unit

            double precision :: RanNum,F,FHold,FHigh1,FHigh2,CR,CRLow,CRHigh
            double precision :: ValueHold,BestValue,BestValueOld,BestValueStop,AcceptRate
            double precision,allocatable :: ParentChrom(:,:),ProgenyChrom(:,:),Chrom(:),MiVal(:),MaVal(:),Value(:)

            logical :: BestSolutionChanged

            LastGenerationPrint=0
            BestSolutionOld=0
            BestValueOld=-999999.0
            BestValueStop=BestValueOld

            allocate(ParentChrom(nParam,nSolution))
            allocate(ProgenyChrom(nParam,nSolution))
            allocate(Chrom(nParam))
            allocate(Value(nSolution))
            allocate(MiVal(nParam))
            allocate(MaVal(nParam))

            ! --- Printout ---

            open(newunit=Unit,file=trim(File),status="unknown")
            write(stdout,"(a)") " "
            !                       12345678901   12345678901   12345678901   12345678901   12345678901   12345678901   12345678901   12345678901   12345678901
            write(stdout,"(9a11)") " SearchMode","       Step","       Gain"," Inbreeding"," ...-Target","  RateOfInb"," ...-Target","  Criterion"," AcceptRate"
            write(Unit,  "(9a11)") " SearchMode","       Step","       Gain"," Inbreeding"," ...-Target","  RateOfInb"," ...-Target","  Criterion"," AcceptRate"

            ! --- Set parameters ---

            ! Crossover rate
            CRHigh=0.4 ! For first few generations (burn-in)
            CRLow=0.2  ! For later climbs

            ! F is multiplier of difference used to mutate
            ! Typically between 0.2 and 2.0
            ! (if alleles should be integer, keep F as integer)
            FHold=0.2  ! Conservative moves
            FHigh1=0.4 ! Adventurous moves
            FHigh2=1.0 ! Adventurous moves

            ! Constrain parameters
            MiVal=0.0
            MaVal=1.0

            ! --- Initialise foundation population of solutions ---

            ParentChrom(:,1)=1.0
            Value(1)=CalcCriterion(nParam,ParentChrom(:,1),CriterionType)
            do Solution=2,nSolution
                do Param=1,nParam                      ! Set a wide range for each parameter
                    call random_number(RanNum)         ! May need integer values for some problems
                    ParentChrom(Param,Solution)=RanNum
                enddo
                Value(Solution)=CalcCriterion(nParam,ParentChrom(:,Solution),CriterionType)
            enddo

            ! --- Evolve ---

            do Generation=1,nGeneration

                ! Vary differential and non-differential mutation to escape valleys
                if (mod(Generation,3) == 0) then
                    DiffOnly=1
                else
                    DiffOnly=0
                endif

                ! Burn-in
                if (Generation < nGenerationBurnIn) then
                    CR=CRHigh
                else
                    CR=CRLow
                endif

                ! Vary mutation rate every few generations
                if (mod(Generation,4) == 0) then
                    F=FHigh1
                else
                    F=FHold
                endif

                if (mod(Generation,7) == 0) then
                    F=FHigh2
                else
                    F=FHold
                endif

                ! --- Generate competitors ---

                ! TODO: paralelize this loop? Is it worth it?
                BestSolutionChanged=.false.
                AcceptRate=0.0
                do Solution=1,nSolution

                    ! --- Mutate and recombine ---

                    ! Get three different solutions
                    SolutionA=Solution
                    do while (SolutionA == Solution)
                        call random_number(RanNum)
                        SolutionA=int(RanNum*nSolution)+1
                    enddo
                    SolutionB=Solution
                    do while ((SolutionB == Solution) .or. (SolutionB == SolutionA))
                        call random_number(RanNum)
                        SolutionB=int(RanNum*nSolution)+1
                    enddo
                    SolutionC=Solution
                    do while ((SolutionC == Solution) .or. (SolutionC == SolutionA) .or. (SolutionC == SolutionB))
                        call random_number(RanNum)
                        SolutionC=int(RanNum*nSolution)+1
                    enddo

                    ! Mate the solutions
                    call random_number(RanNum)
                    Param=int(RanNum*nParam)+1 ! Cycle through parameters starting at a random point
                    do ParamLoc=1,nParam
                        call random_number(RanNum)
                        if ((RanNum < CR) .or. (ParamLoc == nParam)) then
                            ! Recombine
                            call random_number(RanNum)
                            if ((RanNum < 0.8) .or. (DiffOnly == 1)) then
                                ! Differential mutation (with prob 0.8 or 1)
                                Chrom(Param)=ParentChrom(Param,SolutionC) + F*(ParentChrom(Param,SolutionA)-ParentChrom(Param,SolutionB))
                            else
                                ! Non-differential mutation (to avoid getting stuck)
                                call random_number(RanNum)
                                if (RanNum < 0.5) then
                                    call random_number(RanNum)
                                    Chrom(Param)=ParentChrom(Param,SolutionC) * (0.9d0 + 0.2d0 * RanNum)
                                else
                                    call random_number(RanNum)
                                    Chrom(Param)=ParentChrom(Param,SolutionC) + 0.01d0 * F * (ParentChrom(Param,SolutionA) + 0.01d0) * (RanNum - 0.5d0)
                                endif
                            endif
                        else
                            ! Do not recombine
                            Chrom(Param)=ParentChrom(Param,Solution)
                        endif
                        Param=Param+1
                        if (Param > nParam) Param=Param-nParam
                    enddo

                    ! Constrain parameters
                    ! - this is preferablly not done as it slows convergence
                    ! - it is preferable for criterion function to handle constraints
                    ! - for AlphaMate GG found that it is advised to keep parameters
                    !   unconstrained and set negative values to zero in criterion
                    ! do Param=1,nParam
                    !    if (Chrom(Param) < MiVal(Param)) Chrom(Param)=MiVal(Param)
                    !    if (Chrom(Param) > MaVal(Param)) Chrom(Param)=MaVal(Param)
                    ! enddo

                    ! --- Evaluate and Select ---

                    ValueHold=CalcCriterion(nParam,Chrom,CriterionType)  ! Merit of competitor
                    if (ValueHold >= Value(Solution)) then               ! If competitor is better or equal, keep it
                        ProgenyChrom(:,Solution)=Chrom(:)                ! ("equal" to force evolution)
                        Value(Solution)=ValueHold
                        AcceptRate=AcceptRate+1.0
                    else
                        ProgenyChrom(:,Solution)=ParentChrom(:,Solution) ! else keep the old solution
                    endif
                enddo ! nSolution

                AcceptRate=AcceptRate/dble(nSolution)

                ! --- New parents ---

                do Solution=1,nSolution
                    ParentChrom(:,Solution)=ProgenyChrom(:,Solution)
                enddo

                ! --- Find the best solution in this generation ---

                BestSolution=maxloc(Value,dim=1)
                BestValue=Value(BestSolution)
                if (BestValue > BestValueOld) BestSolutionChanged=.true.
                BestSolutionOld=BestSolution
                BestValueOld=BestValue

                ! --- Test if solution is improving to stop early ---

                if (mod(Generation,nGenerationStop) == 0) then
                    if ((BestValue-BestValueStop) > StopTolerance) then
                        BestValueStop=BestValue
                    else
                        write(stdout,"(a,f,a,i,a)") "NOTE: Evolutionary algorithm did not improve objective for ",StopTolerance, " in the last ",nGenerationStop," generations. Stopping."
                        exit
                    endif
                endif

                ! --- Monitor ---

                if (BestSolutionChanged) then
                    if (((Generation - LastGenerationPrint) >= nGenerationPrint)) then
                        LastGenerationPrint=Generation
                        ValueHold=CalcCriterion(nParam,ProgenyChrom(:,BestSolution),CriterionType)
                        write(stdout,"(a11,i11,7f11.4)") CriterionType,Generation,Gain,FCurrent,(FCurrent-FTarget),DeltaFCurrent,(DeltaFCurrent-DeltaFTarget),ValueHold,AcceptRate
                        write(Unit,  "(a11,i11,7f11.4)") CriterionType,Generation,Gain,FCurrent,(FCurrent-FTarget),DeltaFCurrent,(DeltaFCurrent-DeltaFTarget),ValueHold,AcceptRate
                    endif
                endif

            enddo ! Generation

            ! --- Evaluate the winner once more ---

            ValueHold=CalcCriterion(nParam,ProgenyChrom(:,BestSolution),CriterionType)
            write(stdout,"(a11,i11,7f11.4)") CriterionType,Generation,Gain,FCurrent,(FCurrent-FTarget),DeltaFCurrent,(DeltaFCurrent-DeltaFTarget),ValueHold,AcceptRate
            write(Unit,  "(a11,i11,7f11.4)") CriterionType,Generation,Gain,FCurrent,(FCurrent-FTarget),DeltaFCurrent,(DeltaFCurrent-DeltaFTarget),ValueHold,AcceptRate

            close(Unit)

            deallocate(ParentChrom)
            deallocate(ProgenyChrom)
            deallocate(Chrom)
            deallocate(Value)
            deallocate(MiVal)
            deallocate(MaVal)
        end subroutine EvolAlgForAlphaMate

        !#######################################################################

        function CalcCriterion(nInd,Solution,CriterionType) result(Criterion)
            implicit none
            ! Arguments
            integer,intent(in)          :: nInd           ! No. of individuals
            double precision,intent(in) :: Solution(nInd) ! Solution
            character(len=*),intent(in) :: CriterionType  ! Type of criterion (MinInb,MaxGain)
            ! Other
            integer :: i,jMal,jFem
            double precision :: MiVal,Criterion,TmpVec(nInd,1),TotMal,TotFem,GainScaled,SolutionScaled(nInd)

            ! --- Negative-to-zero mapping, i.e., [-Inf,0,Inf] --> [0,0,Inf] ---

            ! (this one allows unconstrained parameters within evolutionary algorithm,
            !  and speeds up convergence of AlphaMate)
            SolutionScaled(:)=Solution(:)
            do i=1,nInd ! TODO: handle nInd here when PAGE comes in
                if (SolutionScaled(i) < 0.0) then
                    SolutionScaled(i)=0.0
                endif
            enddo

            ! --- Probability mapping, i.e., [-Inf,0,Inf] --> [0,1/2,1] ---

            ! For some reason, not limiting within the evolutionary algorithm and
            ! using the logit mapping gave much worse outcomes. Perhaps due to
            ! non-linearity of the link function.
            ! SolutionScaled=Solution(:)
            ! SolutionScaled(:)=exp(Solution(:))/(1.0d0 + exp(Solution(:)))

            ! --- Genetic contributions (XVec) ---

            if (.not.GenderMatters) then
                XVec(1:nInd)=SolutionScaled(1:nInd)/sum(SolutionScaled(1:nInd))
            else
                ! WARNING: the order of males and females differs in Solution and XVec
                TotMal=sum(SolutionScaled(1:nMal))
                TotFem=sum(SolutionScaled((nMal+1):nInd))
                jMal=0
                jFem=0
                do i=1,nInd
                    if (Gender(i) == 1) then
                        jMal=jMal+1
                        XVec(i)=0.5d0*SolutionScaled(jMal)/TotMal
                    else
                        jFem=jFem+1
                        XVec(i)=0.5d0*SolutionScaled(nMal+jFem)/TotFem
                    endif
                enddo
            endif

            ! --- Group coancestry (future inbreeding) ---

            ! xA
            do i=1,nInd
                TmpVec(i,1)=dot_product(XVec,RelMat(:,i))
            enddo
            ! xAx
            FCurrent=0.5d0*dot_product(TmpVec(:,1),XVec)
            ! print*,XVec,TmpVec,FCurrent

            ! Matrix multiplication with symmetric matrix using BLAS routine
            ! (it was slower than the above with 1.000 individuals, might be
            !  benefical with larger cases so kept in for now.)
            ! http://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3.html#ga253c8edb8b21d1b5b1783725c2a6b692
            ! Ax
            ! call dsymm(side="l",uplo="l",m=nInd,n=1,alpha=1.0d0,A=RelMat,lda=nInd,b=XVec,ldb=nInd,beta=0,c=TmpVec,ldc=nInd)
            ! call dsymm(     "l",     "l",  nInd,  1,      1.0d0,  RelMat,    nInd,  XVec,    nInd,     0,  TmpVec,    nInd)
            ! xAx
            ! FCurrent=0.5d0*dot_product(XVec,TmpVec(:,1))
            ! print*,XVec,TmpVec,FCurrent
            ! stop

            FCurrentRebase=(FCurrent-FOld)/(1.0d0-FOld)

            ! --- Genetic gain ---

            Gain=dot_product(XVec,Bv)
            GainScaled=dot_product(XVec,BvScaled)

            ! --- Criterion ---

            if (trim(CriterionType) == "MinInb") then
                DeltaFCurrent=(FCurrent-FOld)/(1.0d0-FOld)
                !Criterion=-1.0d0*(FCurrentRebase-FTargetRebase)
                Criterion=FTargetRebase-FCurrentRebase
            endif

            if (trim(CriterionType) == "MaxGain") then
                DeltaFCurrent=(FCurrent-FOld)/(1.0d0-FOld)
                if (DeltaFCurrent <= DeltaFTarget) then
                    Criterion=GainScaled-(FCurrentRebase-FTargetRebase)
                else
                    ! This gives a penalty of 1 SD for 0.01 increase in inbreeding above the target.
                    ! Note that gain is scaled and inbreeding is rebased so this should be a "stable"
                    ! soft constraint.
                    Criterion=GainScaled-100.0d0*(FCurrentRebase-FTargetRebase)
                    !write(stderr,"(a,8f10.5)") "Illegal",Gain,FCurrent,FTarget,DeltaFCurrent,DeltaFTarget,FCurrentRebase,FTargetRebase,Criterion
                endif
            endif

            return
        end function CalcCriterion

        !#######################################################################

        subroutine PerformMating
            use AlphaSuiteModule,only : RandomOrder
            implicit none

            integer :: i,j,k,l,MatingIdVec(nMatings*2),ShuffleMatingIdVec(nMatings*2),nTmp
            double precision :: SortedXvec(nInd),RanNum

            SortedXvec(:)=XVec(:)
            nMatingPerInd=0
            MatingIdVec=0
            l=0
            do i=1,nInd
                k=maxloc(SortedXvec(:),dim=1)
                SortedXvec(k)=(minval(SortedXvec(:)))-1

                nTmp=int(XVec(k)*(nMatings*2))
                nMatingPerInd(k)=nTmp
                do j=1,nTmp
                    l=l+1 ! l evolves from 1 to nMatings independently of i.
                    MatingIdVec(l)=k
                    if (l == (nMatings*2)) exit
                enddo
                if (l == (nMatings*2)) exit
            enddo

            if (sum(nMatingPerInd(:)) < (nMatings*2)) then ! if all of the remaining matings have not received individuals
                do
                    call random_number(RanNum)
                    k=int(RanNum*nInd)+1
                    call random_number(RanNum)
                    if (RanNum < XVec(k)) then ! ids are selected for these matings according to their contribution.
                        l=l+1
                        MatingIdVec(l)=k
                        nMatingPerInd(k)=nMatingPerInd(k)+1
                    endif
                    if (l == (nMatings*2)) exit
                enddo
            endif

            call RandomOrder(ShuffleMatingIdVec,(nMatings*2))

            l=0
            do i=1,nMatings
                l=l+1
                Matings(i,1)=MatingIdVec(ShuffleMatingIdVec(l))
                l=l+1
                Matings(i,2)=MatingIdVec(ShuffleMatingIdVec(l))
            enddo
        end subroutine PerformMating

        !#######################################################################

        subroutine PerformMatingGenderIncl
            use AlphaSuiteModule,only : RandomOrder
            implicit none

            integer :: i,j,k,l,nTmp,MatingMalIdVec(nMatings),MatingFemIdVec(nMatings),nMatingPerMal(nMal),nMatingPerFem(nFem)
            integer :: ShuffleMatingMalVec(nMatings),ShuffleMatingFemVec(nMatings)
            double precision :: MalXvec(nMal),FemXvec(nFem),MalSortedXvec(nMal),FemSortedXvec(nFem),RanNum

            ! Select male parents

            MatingMalIdVec=0
            nMatingPerMal=0

            l=0
            do i=1,nMal
                MalXvec(i)=XVec(IdMal(i))
            enddo

            MalSortedXvec(:)=MalXvec(:)

            do i=1,nMal
                k=maxloc(MalSortedXvec,dim=1) ! k is the position (id) in vector x having the highest contribution.
                MalSortedXvec(k)=(minval(MalSortedXvec))-1

                nTmp=int(MalXvec(k)*(nMatings)) ! changed! XVec(SortedIdVec(i))=XVec(k)
                nMatingPerMal(k)=nTmp
                do j=1,nTmp
                    l=l+1 ! l evolves from 1 to nMatings independently of i.
                    MatingMalIdVec(l)=IdMal(k) ! MatingIdVec takes the value of k, or SortedIdVec(i), at positions 1 to nTmp.
                    if (l == (nMatings)) exit
                enddo
                if (l == (nMatings)) exit
            enddo

            if (sum(nMatingPerMal) < (nMatings)) then ! if all of the remaining matings have not received individuals.
                do
                    call random_number(RanNum)
                    k=int(RanNum*nMal)+1
                    call random_number(RanNum)
                    if (RanNum < MalXvec(k)) then ! ids are selected for these matings according to their contribution.
                        l=l+1
                        MatingMalIdVec(l)=IdMal(k)
                        nMatingPerMal(k)=nMatingPerMal(k)+1
                    endif
                    if (l == (nMatings)) exit
                enddo
            endif

            ! Select female parents

            MatingFemIdVec=0
            nMatingPerFem=0

            l=0
            do i=1,nFem
                FemXvec(i)=XVec(IdFem(i))
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
                    if (l == (nMatings)) exit
                enddo
                if (l == (nMatings)) exit
            enddo

            if (sum(nMatingPerFem) < (nMatings)) then
                do
                    call random_number(RanNum)
                    k=int(RanNum*nFem)+1
                    call random_number(RanNum)
                    if (RanNum < FemXvec(k)) then
                        l=l+1
                        MatingFemIdVec(l)=IdFem(k)
                        nMatingPerFem(k)=nMatingPerFem(k)+1
                    endif
                    if (l == (nMatings)) exit
                enddo
            endif

            ! Perform matings

            call RandomOrder(ShuffleMatingMalVec,nMatings)
            call RandomOrder(ShuffleMatingFemVec,nMatings)

            l=0
            do i=1,nMatings
                l=l+1
                matings(i,1)=MatingMalIdVec(ShuffleMatingMalVec(l))
                matings(i,2)=MatingFemIdVec(ShuffleMatingFemVec(l))
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

    ! Create output folder
    Success=systemQQ(RMDIR//" AlphaMateResults")
    if (.not.Success) then
        write(stderr,"(a)") "ERROR: Failure to remove old output folder!"
    endif
    Success=systemQQ(MKDIR//" AlphaMateResults")
    if (.not.Success) then
        write(stderr,"(a)") "ERROR: Failure to make output folder!"
    endif

    call ReadParametersAndDataForAlphaMate
    call SetInbreedingParameters
    call AlphaMateSearch
    call cpu_time(Finish)

    write(stdout,"(a)") " "
    write(stdout,"(a,f20.4,a)") "Time duration of AlphaMate: ",Finish-Start," seconds"
    write(stdout,"(a)") " "
end program AlphaMate

!###############################################################################
