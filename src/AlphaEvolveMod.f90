
!###############################################################################

module AlphaEvolveModule

  use ISO_Fortran_Env, STDIN=>input_unit,STDOUT=>output_unit,STDERR=>error_unit
  use AlphaSuiteModule, only : Int2Char,Real2Char

  implicit none

  type :: EvolveCrit
    ! TODO: How do I get this done generically? This is type specific for AlphaMate
    !       I have created an extended type in AlphaMateModule of base EvolveCrit here,
    !       but got bitten through using polymorphic (class() stuff) objects, where I
    !       am not allowed to do a=b etc.
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
  public :: EvolAlgDE,EvolveCrit

  contains

    !###########################################################################

    subroutine EvolAlgDE(nParam,nSol,Init,nGen,nGenBurnIn,nGenStop,&
      StopTolerance,nGenPrint,File,CritType,CRBurnIn,CRLate,FBase,FHigh1,FHigh2,&
      CalcCriterion,LogHeader,Log,BestCriterion)

      implicit none

      ! Arguments
      integer(int32),intent(in)         :: nParam         ! No. of parameters in a solution
      integer(int32),intent(in)         :: nSol           ! No. of solutions to test each generation
      real(real64),intent(in),optional  :: Init(:,:)      ! Initial solutions to start with
      integer(int32),intent(in)         :: nGen           ! No. of generations to run
      integer(int32),intent(in)         :: nGenBurnIn     ! No. of generations with more
      integer(int32),intent(in)         :: nGenStop       ! Stop after no progress for nGenerationStop
      real(real64),intent(in)           :: StopTolerance  ! Stopping tolerance
      integer(int32),intent(in)         :: nGenPrint      ! Print changed solution every nGenerationPrint
      character(len=*),intent(in)       :: File           ! Which file to write to
      character(len=*),intent(in)       :: CritType       ! Passed to CalcCriterion
      real(real64),intent(in),optional  :: CRBurnIn       ! Crossover rate for nGenBurnIn
      real(real64),intent(in),optional  :: CRLate         ! Crossover rate
      real(real64),intent(in),optional  :: FBase          ! F is multiplier of difference used to mutate
      real(real64),intent(in),optional  :: FHigh1         ! F is multiplier of difference used to mutate
      real(real64),intent(in),optional  :: FHigh2         ! F is multiplier of difference used to mutate
      type(EvolveCrit),intent(out)      :: BestCriterion  ! Criterion for the best solution found

      interface
        subroutine CalcCriterion(Sol,CritType,Criterion)
          use ISO_Fortran_Env
          import :: EvolveCrit
          real(real64),intent(inout)   :: Sol(:)
          character(len=*),intent(in)  :: CritType
          type(EvolveCrit),intent(out) :: Criterion
        end subroutine

        subroutine LogHeader(LogUnit)
          use ISO_Fortran_Env
          integer(int32),intent(in) :: LogUnit
        end subroutine

        subroutine Log(LogUnit,Gen,AcceptRate,Criterion)
          use ISO_Fortran_Env
          import :: EvolveCrit
          integer(int32),intent(in)   :: LogUnit
          integer(int32),intent(in)   :: Gen
          real(real64),intent(in)     :: AcceptRate
          type(EvolveCrit),intent(in) :: Criterion
        end subroutine
      end interface

      ! Other
      integer(int32) :: nInit,Param,ParamLoc,Sol,Gen,LastGenPrint
      integer(int32) :: SolA,SolB,SolC,BestSol
      integer(int32) :: Unit

      real(real64) :: RanNum,FInt,FBaseInt,FHigh1Int,FHigh2Int,CRInt,CRBurnInInt,CRLateInt
      real(real64) :: AcceptRate,OldChrom(nParam,nSol),NewChrom(nParam,nSol),Chrom(nParam)

      logical :: DiffOnly,BestSolChanged

      character(len=100) :: DumC,DumC2

      type(EvolveCrit) :: Criterion(nSol),CriterionHold,BestCriterionStop

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

        ! TODO: Paralelize this loop? Is it worth it? Would be possible if
        !       CalcCriterion would be pure (no I/O and having side effects
        !       outside of the function/subroutine, like module variables)
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
            DumC=Real2Char(StopTolerance)!,fmt="(f9.5)")
            DumC2=Int2Char(nGenStop)
            write(STDOUT,"(5a)") "NOTE: Objective did not improve for ",trim(adjustl(DumC))," in the last ",trim(adjustl(DumC2))," generations. Stopping the optimisation."
            write(STDOUT,"(a)") " "
            exit
          end if
        end if

      end do ! Gen

      ! --- The winner ---

      ! Re-evaluate the winner to fill any potential global objects!!!
      ! (if CalcCriterion does not include any global objects, then this would not be needed)
      call CalcCriterion(NewChrom(:,BestSol),CritType,BestCriterion)
      call Log(Unit,Gen,AcceptRate,BestCriterion)
      write(STDOUT,"(a)") " "
      close(Unit)
    end subroutine

    !###########################################################################
end module
