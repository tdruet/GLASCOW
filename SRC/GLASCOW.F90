program GLMM
  !$ use OMP_LIB
  !Copyright 2012, Frederic Farnir, Francois Guillaume, Tom Druet

  ! This program is free software: you can redistribute it and/or modify
  ! it under the terms of the GNU General Public License as published by
  ! the Free Software Foundation, either version 3 of the License, or
  ! (at your option) any later version.
  !
  ! This program is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ! GNU General Public License for more details.
  !
  ! You should have received a copy of the GNU General Public License
  ! along with this program.  If not, see <http://www.gnu.org/licenses/>.

  implicit none
  !$Id: GLASCOW.F90 24 2012-09-14 09:39:00Z fguillaume $

  ! This program implements various GLMM 
  ! Fixed effects can be any combination of covariates and class effects
  ! Random effects are 'individual' effect and 'haplotype' effects:
  !
  !     eta = Xb + Za + Hh + e
  !
  ! with V(a)=A*va, V(h)=I*vh and Ve=I*ve
  ! 
  ! The analysis is performed in two steps
  ! a) eta = Xb + Za + e => e = eta - Xb - Zu
  ! b) e = mu + Hh + e
  !
  ! The (b) step is repeated for each position (where H changes) and
  ! a score value is computed in each position
  ! Finally, a score p-value is computed using permutations
  !
  !========================================================
  ! Parameters
  !========================================================
  integer, parameter::dp = selected_real_kind(12,60)
  !========================================================
  ! Variables declarations
  !========================================================
  ! Type declarations
  type parameters
     integer::nb_rec                                   ! Number of records
     integer::nb_fix                                   ! Number of fixed effects
     integer::nb_rnd                                   ! Number of random effects
     integer::nb_pos                                   ! Number of genomic positions
     character(30)::aname                              ! Name of the analysis
     character(100)::dfile                             ! Data file
     character(100)::hfile                             ! Haplotypes file
     integer::datid                                    ! Position of id in data file
     integer::hapid                                    ! Position of id in haplos file
     integer,allocatable::dim_fix(:)                   ! Dimensions of fixed effects
     integer,allocatable::dim_rnd(:)                   ! Dimensions of random effects
     integer,allocatable::pos_fix(:)                   ! Positions of fixed effects
     integer,allocatable::pos_rnd(:)                   ! Positions of random effects
     character(30),allocatable::name_fix(:)            ! Names of fixed effects
     character(30),allocatable::name_rnd(:)            ! Names of random effects
     character(100)::rfiles(2)
     character(100),dimension(2)::rifiles              ! Inverse correlation files
     character(6)::linkf                               ! Link function
     character(6)::distr                               ! Trait distribution function
     character(30)::name_trt                           ! Name of trait
     integer::pos_trt                                  ! Position of trait
     real(kind=dp)::start_tau2                         ! Initial value for error variance
     real(kind=dp),dimension(2)::start_sigma2          ! Initial values for random effects variances
     character(100)::start_b                           ! File with initial values for fixed effects solutions
     character(100)::start_u                           ! File with initial values for random effects solutions
     character(2)::reml_method                         ! Method used to obtain REML estimates of VC
     integer::reml_maxiter                             ! Max number of allowed REML iterations 
     real(kind=dp)::reml_maxtol                        ! Max allowed diff between estimates in REML
     integer,dimension(2)::scores_pos                  ! Positions for which scores should be computed
     character(5)::score_method                        ! Method used to compute score
     integer::nb_scores                                ! Number of scores to be computed
     integer::slice                                    ! Number of positions to be read in at a time
     integer::showp                                    ! Control display of positions on screen
     integer::nb_perms                                 ! Number of permutations to be computed
  end type parameters
  type(parameters)::p_p
  ! Correlation and inverse correlation matrices
  type correl
     integer::n                                        ! Size of the correlation matrices
     logical::is_direct                                ! Check presence
     logical::is_inverse                               ! Check presence
     logical::is_diagonal                              ! Check whether diagonal
     logical::is_identity                              ! Check whether identity
     real(kind=dp),allocatable::dir_cor(:,:)           ! Direct correlation matrix
     real(kind=dp),allocatable::inv_cor(:,:)           ! Inverse correlation matrix
  end type correl
  type(correl)::corr
  ! Model variables
  real(kind=dp),allocatable::X(:,:),Z(:,:),y(:),b(:),u(:)
  integer,allocatable::pheno_id(:),start_indices(:)
  ! Declarations for variances
  real(kind=dp)::tau2                                ! Residual variance 
  real(kind=dp),dimension(2)::sigma2                 ! Random effects variances 
  ! Declarations for scores
  real(kind=dp),dimension(:,:),allocatable::scores     ! Scores
  real(kind=dp),dimension(:,:),allocatable::p_values   ! Scores p-values
  real(kind=dp),dimension(:),allocatable::resid,eta,mu ! Utilities for scores
  ! Work variables
  integer::rc,ok,nbr,i,j,k,is,nb_fix_tot,nb_rnd_tot,base,info,first_fact,io 
  logical::blankp
  character(3)::version
  character(100)::prmfile,aname
  integer,dimension(:),allocatable::use_fix,use_rnd
  real(kind=dp),allocatable::ci(:,:)
  !---------------------------------------------------------------------!
  integer::narg,cptArg   !Check number of argument
  Logical::ExFile !Check if file exist
  Logical::DbgMode=.FALSE.,Nocompute=.FALSE.,ExportSol=.FALSE.
  Logical::SkipREML=.FALSE.,exportcorr=.FALSE.,ExportPerm=.FALSE.
  Logical::exportresidual=.FALSE.,DoReadMap=.FALSE.,WaitNumThread=.FALSE.
  !$ integer::NumOfThreads=1
  integer,dimension(:)::racine(12)
  character(len=100)::filename
  !---------------------------------------------------------------------!
  real::MCsg,VCsg,CsgThre
  logical::SkipTwin=.FALSE.,PossTwin=.FALSE.
  !----------------------------------------------------------------------!
  real :: temps, t_cpu_0, t_cpu_1, t_cpu
  integer::ir,t1,t0
  !----------------------------------------------------------------------!
  real(kind=dp),allocatable::R(:)
  character(len=30),allocatable::GenoFileName(:)
  integer,allocatable::pardim(:,:),IPIV(:)


  !========================================================
  ! Code
  !========================================================
  rc=0
  !---------------------------------------------------------------------!
  ! Obtain the name of the analysis                                     !
  !---------------------------------------------------------------------!
  narg=command_argument_count()
  if(narg<1)THEN
     write(*,*)"At least one argument should follow glascow invocation : Process will stop"
     stop
  ELSE
     call get_command_argument(1,aname)
  endif
  !---------------------------------------------------------------------!
  ! Obtain the name of the parameters file                              !
  !---------------------------------------------------------------------!
  prmfile=trim(aname)
  if(narg==2)call get_command_argument(2,aname)
  inquire(file=prmfile,EXIST=ExFile) !check if file exist
  if(.not.ExFile.and.(adjustl(aname).ne."-init"))THEN
     write(*,*)"File ",prmfile," doesn't exist : Process will stop"
     stop
  endif
  !---------------------------------------------------------------------!
  ! If more arg exist read them
  !---------------------------------------------------------------------!
  if(narg>1)Then
     !loop on additionnal arguments 
     LoopArg: do cptArg=2,narg
        call get_command_argument(cptArg,aname)    
        !If OMP--------------------------
        if((len_trim(adjustl(aname))==2).and.(trim(adjustl(aname))=="-n"))then
           WaitNumThread=.TRUE.
           cycle LoopArg
        endif
        if(WaitNumThread)then
           WaitNumThread=.FALSE.
           !$ read(aname,*)NumOfThreads
           !$ cycle LoopArg
           write(*,*)"This version was not compiled with OPENMP enabled"
           write(*,*)"The -n option will have no effect"
           cycle LoopArg
        endif
        !Other arguments
        select case(adjustl(aname))
        case("-version") 
           version='1.0'
           call banner(version)
           STOP 
        case("-dbg")
           write(*,'(3a)')"Option ",adjustl(aname)," : ON"
           write(*,*)"Additional output will be written"
           DbgMode=.TRUE.
        case("-init")
           call InitParam(prmfile)
        case("-map")
           write(*,'(a)')"Name of the map file ?"
           read(*,*)filename
           inquire(file=filename,EXIST=ExFile) !check if file exist
           if(.not.ExFile)THEN
              write(*,*)"File ",filename," doesn't exist : Process will stop"
              stop
           endif
           DoReadMap=.TRUE.
        case("-seed")
           racine(2)=56
           write(*,*)"seed ?"
           read(*,*)racine(1)
           call random_seed(put=racine)
        case("-nocompute")
           Nocompute=.TRUE.
        case("-exportperm")
           ExportPerm=.TRUE.
        case("-exportsol")
           Exportsol=.TRUE.
        case("-exportcorr")
           exportcorr=.TRUE.
        case("-exportresidual")
           exportResidual=.TRUE.
        case("-skipreml")
           write(*,*)"No REML estimation"
           SkipREML=.TRUE.
           Nocompute=.TRUE.
        case("-allowtwins")
           SkipTwin=.TRUE.
        case default 
           write(*,'(3a)')"Option ",trim(adjustl(aname))," unknown"
        end select
     end do LoopArg
  endif
  !Set OMP NUM THREADS
  !$ call omp_set_num_threads(NumOfThreads)   
  !---------------------------------------------------------------------!
  ! Get mixed model dimensions (nb_fix,nb_rnd)                          !
  !---------------------------------------------------------------------!
  rc=Get_Dimensions(p_p,prmfile)
  if (rc.ne.0) goto 996
  !---------------------------------------------------------------------!
  ! Allocate space for fixed and random effects data                    !
  !---------------------------------------------------------------------!
  allocate(p_p%name_fix(p_p%nb_fix),p_p%name_rnd(p_p%nb_rnd),& 
       &p_p%pos_fix(p_p%nb_fix),p_p%pos_rnd(p_p%nb_rnd),&     
       &p_p%dim_fix(p_p%nb_fix),p_p%dim_rnd(p_p%nb_rnd),stat=ok)
  if (ok.ne.0)then
     write(*,*)'Error: cannot allocate parts of parameters structure'
  endif

  !---------------------------------------------------------------------!
  ! Obtain parameters from parameters file                              !
  !---------------------------------------------------------------------!
  rc=Get_Parameters(p_p,prmfile)
  if (rc.ne.0) goto 999
  nbr=p_p%nb_rec
  !---------------------------------------------------------------------!
  ! Obtain data from data file (y)                             
  !---------------------------------------------------------------------!
  allocate(start_indices(p_p%nb_fix+p_p%nb_rnd))
  allocate(y(nbr),pheno_id(nbr),stat=ok)
  if (ok.ne.0) goto 990
  rc=Get_Data(p_p%dfile,p_p%pos_trt,y,p_p%datid,pheno_id,nbr)
  if (rc.gt.0) goto 995
  !---------------------------------------------------------------------!
  ! Create design matrices (X, Z, b, u)
  !---------------------------------------------------------------------!
  ! Fixed effects incidence matrix
  nb_fix_tot=0;first_fact=0
  do i=1,p_p%nb_fix
     if (p_p%dim_fix(i).gt.0) then
        nb_fix_tot=nb_fix_tot+p_p%dim_fix(i)
        if(first_fact==0)then
           first_fact=i
        else
           nb_fix_tot=nb_fix_tot-1
        endif
     else
        nb_fix_tot=nb_fix_tot+1   ! Covariate
     endif
  enddo
  allocate(X(nbr,p_p%nb_fix),b(nb_fix_tot),stat=ok)
  if (ok.ne.0) goto 990
  X=0
  ! Random effects incidence
  nb_rnd_tot=0
  do i=1,1 ! only for POLYG
     nb_rnd_tot=nb_rnd_tot+p_p%dim_rnd(i)
  enddo
  allocate(Z(nbr,1),u(nb_rnd_tot),stat=ok)
  if (ok.ne.0) goto 990
  Z=0  

  !---------------------------------------------------------------------!
  ! Get design matrices (X) from data file
  !---------------------------------------------------------------------!
  base=0
  blankp=.true.
  do is=1,p_p%nb_fix
     rc=Get_Design(p_p%dfile,is,p_p%pos_fix(is), &
          p_p%dim_fix(is),X,nbr,p_p%nb_fix,base,blankp)
     if (rc.gt.0) goto 988
     if(p_p%dim_fix(is).gt.0 .and. is.gt.first_fact)p_p%dim_fix(is)=p_p%dim_fix(is)-1
     base=base+p_p%dim_fix(is)
  enddo

  !---------------------------------------------------------------------!
  ! Get rank of X matrix
  !---------------------------------------------------------------------!
  if (rc.gt.0) goto 992

  !---------------------------------------------------------------------!
  ! Get design matrices (Z) from data file for non-haplo effects
  ! In this version, this design matrix corresponds to polygenic effect
  !---------------------------------------------------------------------!
  blankp=.true.
  rc=Get_Design(p_p%dfile,1,p_p%pos_rnd(1), &
       p_p%dim_rnd(1),Z,nbr,p_p%nb_rnd,0,blankp)
  if (rc.gt.0) goto 988

  !---------------------------------------------------------------------!
  ! Initialize VC matrix sigma2
  !---------------------------------------------------------------------!
  sigma2=0  
  !---------------------------------------------------------------------!
  ! Allocate (inverse) correlation data for POLYG
  !---------------------------------------------------------------------!
  is=p_p%dim_rnd(1)
  ok=0
  if (p_p%rfiles(1).ne.'I') then
     allocate(corr%dir_cor(is,is),stat=ok)
     if (ok.ne.0) goto 990
     corr%dir_cor=0
  endif
  if (p_p%rifiles(1).ne.'I') then
     allocate(corr%inv_cor(is,is),stat=ok)
     if (ok.ne.0) goto 990
     corr%inv_cor=0
  endif
  ! Only allocate diagonals, else. 
  if ((p_p%rfiles(1).eq.'I').and.(p_p%rifiles(1).eq.'I')) then
     allocate(corr%dir_cor(is,1),stat=ok)
     if (ok.ne.0) goto 990
     corr%dir_cor=0
     allocate(corr%inv_cor(is,1),stat=ok)
     if (ok.ne.0) goto 990
     corr%inv_cor=0
  endif
  !---------------------------------------------------------------------!
  ! Obtain (inverse) POLYG correlation data
  !---------------------------------------------------------------------!

  rc=Get_Corr(p_p,corr,1) 
  if (rc.gt.0) goto 994

  !---------------------------------------------------------------------!
  ! Step (a): solve without haplotype effect
  !---------------------------------------------------------------------!
  !---------------------------------------------------------------------!
  ! Set initial values for parameters (VC + b + u)
  !---------------------------------------------------------------------!
  tau2=p_p%start_tau2                                 ! Initial value
  sigma2(1:p_p%nb_rnd)=p_p%start_sigma2(1:p_p%nb_rnd) ! Initial values
  rc=Initialize_Solutions(p_p,b,nb_fix_tot,u,nb_rnd_tot, &
       p_p%start_b,p_p%start_u)
  if (rc.gt.0) goto 987                        
  !---------------------------------------------------------------------!
  ! Check availability of inverse correlation matrices
  !---------------------------------------------------------------------!
  if ((p_p%rfiles(1).ne.'I').and.(p_p%rifiles(1).eq.'I')) then
     is=p_p%dim_rnd(1)
     allocate(corr%inv_cor(is,is),stat=ok)
     if (ok.ne.0) goto 990
     corr%inv_cor=corr%dir_cor
     write(*,*)"=> Computing inverse correlation matrix for effect POLYG"
     !======================================================================
     Twins:if(.not.SkipTwin)then
        !Check if no individuals are twins
        MCsg=0 ;VCsg=0 ;CsgThre=0
        if(allocated(R))THEN
           do i=1,is 
              MCsg=R(IAd(i,i))+MCsg
              VCsg=R(IAd(i,i))**2+VCsg
           end do
           Mcsg=Mcsg/is
           VCsg=sqrt((VCsg/is)-(Mcsg)**2)
           write(*,'(a,f6.3,a,f6.3)')"Average relationship gii :",Mcsg," std ",VCsg
           CsgThre=MCsg-(1.96*VCsg)
           Mcsg=(Sum(R)-(Mcsg*is))/(Iad(is,is)-is)     
           write(*,*)"Average relationship between individuals gij:",Mcsg," std ",VCsg 
           do i=1,is-1
              do j=i+1,is
                 if(R(IAd(i,j))>CsgThre)THEN
                    write(*,'(a,i9,a,i9,a,f7.4,a,f7.4,a,f7.4,a)')"Relationship between ",i," and ",j," is ",R(IAd(i,j))&
                         &," whereas IC95 for average inbreeding : [",CsgThre," : ",CsgThre+(3.92*VCsg),&
                         &"]"
                    PossTwin=.TRUE.
                 endif
              end do
           end do
        ELSE
           do i=1,is 
              MCsg=corr%dir_cor(i,i)+MCsg
              VCsg=(corr%dir_cor(i,i)**2)+VCsg
           end do
           Mcsg=Mcsg/is
           VCsg=sqrt((VCsg/is)-(Mcsg)**2)
           write(*,'(a,f7.4,a,f7.4)')"Average relationship gii :",Mcsg," std ",VCsg 
           CsgThre=MCsg-(1.96*VCsg)
           Mcsg=(Sum(corr%dir_cor)-(Mcsg*is))/(is*(is-1))
           write(*,*)"Average relationship between individuals gij:",Mcsg," std ",VCsg 
           do i=1,is-1
              do j=i+1,is
                 if(corr%dir_cor(i,j)>CsgThre)THEN
                    write(*,'(a,i9,a,i9,a,f7.4,a,f7.4,a,f7.4,a)')"Relationship between ",i," and ",j," is ",corr%dir_cor(i,j)&
                         &," whereas IC95 for average inbreeding : [",CsgThre," : ",CsgThre+(3.92*VCsg),&
                         &"]"
                    PossTwin=.TRUE.                 
                 endif
              end do
           end do
        endif
        !If we suspect the existence of twin stop (REML would either not
        !converge or lead to dubious results)
        if(PossTwin)THEN
           write(*,*)"Some individuals appear to be higly similar this may cause numerical problem"
           write(*,*)"You should treat the latter as only one individual (e.g. remove one individual of the couple)"
           write(*,*)"GLASCOW will now stop"
           if(DbgMode)write(*,*)"To force computation despite this warning use option -allowtwins "
           stop
        endif
     endif Twins

     !======================================================================

#if LAPACK
     !Si R bien stocke ou pas
     if(allocated(R))THEN
        call DPPTRF('U',is,R,info)      
        if(info/=0)write(*,*)"Pb LAPACK ",info
        call DPPTRI('U',is,R,info)
        if(info/=0)write(*,*)"Pb LAPACK ",info        
     else           
        call DPOTRF('U',is,corr%dir_cor,is,info)
        call DPOTRI('U',is,corr%dir_cor,is,info)
     endif
     do j=1,is
        do k=1,j
           if(allocated(R))THEN!passage stockage triangulaire en carre
              corr%dir_cor(k,j)=R(IAd(k,j))
           endif
           corr%dir_cor(j,k)=corr%dir_cor(k,j)
        enddo
     enddo
#else
     if(allocated(R))THEN
        do j=1,is
           do k=1,j
              if(allocated(R))THEN!passage stockage triangulaire en carre
                 corr%dir_cor(k,j)=R(IAd(k,j))
              endif
              corr%dir_cor(j,k)=corr%dir_cor(k,j)
           enddo
        enddo
     endif
     call invert(corr%dir_cor,is)
#endif
     corr%is_inverse=.true.
     ! Store upper triangle of inverse correlation matrix
     if(exportcorr)open(20,file=trim(p_p%aname)//'.invcorr.polyg',err=982)
     do k=1,is
        do j=k,is
           corr%inv_cor(k,j)=corr%dir_cor(k,j)
           corr%inv_cor(j,k)=corr%dir_cor(k,j)
           if(exportcorr.and.(abs(corr%inv_cor(k,j)).gt.0.0)) write(20,*) k,j,corr%inv_cor(k,j)
        enddo
     enddo
     if(exportcorr)close(20)
  endif

  !---------------------------------------------------------------------!
  ! Compute REML estimates of tau2 and sigma2(1)
  !---------------------------------------------------------------------!
  write(*,*)' => Computing REML estimates' !
  allocate(use_fix(p_p%nb_fix),use_rnd(p_p%nb_rnd),stat=ok)
  if (ok.ne.0) goto 990
  use_fix=1     ! Include all fixed effects
  use_rnd(1)=1  ! Include polygenic effect
  use_rnd(2)=0  ! Exclude haplotypic effect
  !If we want to run REML
  CdtReml:if(.not.SkipREML)then
     write(*,*)" => Computing REML estimates" 
     rc=Compute_REML(y,X,Z,corr,use_rnd,nbr,nb_fix_tot,nb_rnd_tot, &
          p_p%nb_fix,p_p%nb_rnd,p_p%dim_fix,p_p%dim_rnd, &
          p_p%linkf,p_p%distr,b,u,tau2,sigma2, &
          p_p%reml_maxiter,p_p%reml_maxtol,p_p%reml_method)
     if (rc.gt.0) goto 989                  
     if (rc.lt.0) print *,' => Warning: REML convergence was not achieved'

     !---------------------------------------------------------------------!
     ! Report REML results
     !---------------------------------------------------------------------!
     if(ExportSol)then 
        ! Reinitialization file
        open(10,file=trim(p_p%aname)//'.fix.res')
        base=0
        do i=1,p_p%nb_fix
           if(p_p%dim_fix(i)==0 .and. sum(X(:,base+1))==0.00)write(10,*) 'F ',p_p%name_fix(i),1,0.00 
           if(p_p%dim_fix(i)==0 .and. sum(X(:,base+1))/=0.00)write(10,*) 'F ',p_p%name_fix(i),1,b(base+1)
           if(i>first_fact .and. p_p%dim_fix(i)>0)write(10,*) 'F ',p_p%name_fix(i),1,0.00 
           do j=1,p_p%dim_fix(i)
              if(i==first_fact)write(10,*) 'F ',p_p%name_fix(i),j,b(base+j)
              if(i>first_fact)write(10,*) 'F ',p_p%name_fix(i),j+1,b(base+j)
           enddo
           base=base+p_p%dim_fix(i)
           if(p_p%dim_fix(i)==0 .and. sum(X(:,i))/=0.00)base=base+1
        enddo
        close(10)
        open(10,file=trim(p_p%aname)//'.rnd.res')
        base=0
        do i=1,p_p%nb_rnd
           if (use_rnd(i).gt.0) then
              do j=1,p_p%dim_rnd(i)
                 write(10,*) p_p%name_rnd(i),j,u(base+j) 
              enddo
           endif
           base=base+p_p%dim_rnd(i)
        enddo
        close(10)
     end if
  endif CdtReml
  deallocate(use_rnd,use_fix)

  !-------------------------------------!
  ! Start looping on scores positions   !
  !-------------------------------------!
  if (p_p%nb_scores.gt.0) then
     write(*,*)' => Now working on scores...'
     !-------------------------------------!
     !In case, we skipped the construction of G^{-1}
     !-------------------------------------!
     if(allocated(GenoFileName))THEN
        if(pardim(1,3).ne.1)then
           !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(pardim,GenoFileName)
           do i=1,size(GenoFileName,1)
              if(DbgMode)write(*,*)"check dimension of ",GenoFileName(i)
              call LecDimTypage(GenoFileName(i),pardim(i,1:4))
           end do
           !$OMP END PARALLEL DO                         
           do i=2,size(GenoFileName,1)
              pardim(i,3)=pardim(i,3)+pardim(i-1,4)
              pardim(i,4)=pardim(i,4)+pardim(i-1,4)
           end do
        endif
     endif
     !---------------------!
     ! Computing residuals !
     !---------------------!
     allocate(resid(p_p%nb_rec),eta(p_p%nb_rec),mu(p_p%nb_rec),stat=ok)
     if (ok.ne.0) goto 990
     eta=0.00
     CdtNoREML:if(.not.SkipREML)then     
        do i=1,p_p%nb_rec
           do j=1,p_p%nb_fix
              if(p_p%dim_fix(j)>0 .and. X(i,j)==0)cycle 
              if(p_p%dim_fix(j)==0)then
                 k=start_indices(j)+1
                 eta(i)=eta(i)+X(i,j)*b(k)
              else
                 k=start_indices(j)+X(i,j)
                 eta(i)=eta(i)+b(k)
              endif
           enddo
           do j=1,1 ! only for polyg
              k=Z(i,j)
              eta(i)=eta(i)+1.0*u(k)
           enddo
        enddo
        rc=get_mu(p_p%linkf,eta,mu,p_p%nb_rec)  
        if (rc.ne.0) goto 983

        resid=y-mu
        !---------------------!
        ! Reporting residuals !
        !---------------------!
        if(exportResidual)then
           open(10,file=trim(p_p%aname)//'.residuals',err=981,action="WRITE")
           do i=1,size(resid,1)
              write(10,*) i,resid(i)
           enddo
           close(10)
        endif
     else 
        !---------------------!
        ! Reading residuals   !
        !---------------------!
        inquire(file=trim(p_p%aname)//'.residuals',exist=Exfile)
        if(Exfile)THEN
           open(10,file=trim(p_p%aname)//'.residuals',err=981,action="read")
        ELSE 
           write(*,'(2a)')"No residual file named ",trim(p_p%aname)//'.residuals'
        endif
        do
           read(10,*,iostat=io) i,resid(i)
           if(io/=0)exit
           if(i>size(resid,1))then
              write(*,*)"Error number of residuals is bigger than number of phenotypes"
              stop
           endif
        enddo
        close(10)
        tau2=p_p%start_tau2
     endif CdtNoREML
     !---------------------!
     ! Computing scores    !
     !---------------------!
     allocate(scores(p_p%nb_scores,3),p_values(p_p%nb_scores,2),stat=ok)
     if (ok.ne.0) goto 990
     rc=Get_Scores(p_p,resid,pheno_id,tau2,scores,p_values)
     if (rc.ne.0) goto 985
     deallocate(p_values,scores)
     deallocate(mu,eta,resid)
  endif

  goto 1000
981 continue
  print *,'Error: problem while trying to create residuals file'
  goto 1000
982 continue
  print *,'Error: problems while creating a correlation file'
  goto 1000
983 continue
  print *,'Error: problems while computing mu'
  goto 1000
985 continue
  print *,'Error: problem while computing scores'
  goto 1000
987 continue
  print *,'Error: problems while initializing solutions'
  goto 1000
988 continue
  print *,'Error: problems while obtaining design matrices'
  goto 1000
989 continue
  print *,'Error: problems while computing REML'
  goto 1000
990 continue
  print *,'Error: problems while allocating'
  goto 1000
992 continue
  print *,'Error: problems while computing the rank of X'
  goto 1000
994 continue
  print *,'Error: problem occurred while obtaining correlations'
  goto 1000
995 continue
  print *,'Error: problems while reading data file'
  goto 1000
996 continue
  print *,'Error: problem occurred while obtaining dimensions'
  goto 1000
999 continue
  print *,'Error: problems while reading parameters from PRM file'
  goto 1000
  !---------------------------------------------------------------------!
  ! That's all, folks
  !---------------------------------------------------------------------!
1000 continue
  ! Deallocating design matrices
  if (allocated(corr%dir_cor)) deallocate(corr%dir_cor,stat=ok)
  if (ok.ne.0) write(*,*) 'Error: problem while deallocating dir_cor'
  if (allocated(corr%inv_cor)) deallocate(corr%inv_cor,stat=ok)
  if (ok.ne.0) write(*,*) 'Error: problem while deallocating inv_cor'
  if (allocated(u)) deallocate(u,stat=ok)
  if (ok.ne.0) write(*,*) 'Error: problem while deallocating u'
  if (allocated(Z)) deallocate(Z,stat=ok)
  if (ok.ne.0) write(*,*) 'Error: problem while deallocating Z'
  if (allocated(b)) deallocate(b,stat=ok)
  if (ok.ne.0) write(*,*) 'Error: problem while deallocating b'
  if (allocated(X)) deallocate(X,stat=ok)
  if (ok.ne.0) write(*,*) 'Error: problem while deallocating X'
  if (allocated(pheno_id)) deallocate(pheno_id,stat=ok)
  if (ok.ne.0) write(*,*) 'Error: problem while deallocating pheno_id'
  if (allocated(y)) deallocate(y,stat=ok)
  if (ok.ne.0) write(*,*) 'Error: problem while deallocating y'
  ! Deallocate parameters structure
  if (allocated(p_p%dim_rnd)) deallocate(p_p%dim_rnd,stat=ok)
  if (ok.ne.0) write(*,*) 'Error: problem while deallocating dim_rnd'
  if (allocated(p_p%dim_fix)) deallocate(p_p%dim_fix,stat=ok)
  if (ok.ne.0) write(*,*) 'Error: problem while deallocating dim_fix'
  if (allocated(p_p%pos_rnd)) deallocate(p_p%pos_rnd,stat=ok)
  if (ok.ne.0) write(*,*) 'Error: problem while deallocating pos_rnd'
  if (allocated(p_p%pos_fix)) deallocate(p_p%pos_fix,stat=ok)
  if (ok.ne.0) write(*,*) 'Error: problem while deallocating pos_fix'
  if (allocated(p_p%name_rnd)) deallocate(p_p%name_rnd,stat=ok)
  if (ok.ne.0) write(*,*) 'Error: problem while deallocating name_rnd'
  if (allocated(p_p%name_fix)) deallocate(p_p%name_fix,stat=ok)
  if (ok.ne.0) write(*,*) 'Error: problem while deallocating name_fix'

contains

  !======================================================================
  ! SUBROUTINES AND FUNCTIONS
  !======================================================================

  function Permute(y,n,yp,m)
    !---------------------------------------------------------------------!
    ! Obtain m permutations of the observations                           !
    !---------------------------------------------------------------------!
    ! y : vector of size n containing phenotypes                          !
    ! n : # of observations                                               !
    ! yp: Matrice of m permutations of y                                  !
    ! m : # of permutations                                               !
    !---------------------------------------------------------------------!
    integer,intent(in)::m,n
    real(kind=dp),intent(in)::y(n)
    real(kind=dp),intent(out)::yp(n,m)
    integer::Permute
    !
    real(kind=dp)::h(n),x
    integer::ip,i,j
    Permute=0
    do ip=1,m
       yp(:,ip)=y
       call random_number(h)
       do i=n,1,-1
          j=int(h(i)*i)+1
          x=yp(i,ip)
          yp(i,ip)=yp(j,ip)
          yp(j,ip)=x
       enddo
    enddo
  end function Permute

  function get_W(t2,v,d,y,n)
    !----------------------------------------------------------------------
    ! W=diag(Wi)=inv[var(i)*g'^{2}]
    ! tag provides the link function
    ! y=dg_dmu(x) 
    !--------------------------------------- 
    !
    integer,intent(in)::n
    real(kind=dp),intent(in)::v(n),d(n),t2
    real(kind=dp),intent(out)::y(n)
    integer::get_W
    !
    integer::i
    real(kind=dp)::tol
    !
    tol=1.0d-10
    get_W=0
    y=t2*v*d*d
    do i=1,n
       if (abs(y(i)).lt.tol) get_W=1
    enddo
    y=1/y   
  end function get_W

  function get_delta(linkf,x,y,n)
    !--------------------------------------- 
    ! Derivative functions
    !  linkf provides the link function
    !  y=dg_dmu(x) 
    !--------------------------------------- 
    !
    character(6),intent(in)::linkf
    integer,intent(in)::n
    real(kind=dp),intent(in)::x(n)
    real(kind=dp),intent(out)::y(n)
    integer::get_delta
    !
    get_delta=0
    select case(linkf)
    case('LOGIT ')
       y=1/(1-x)/x
    case('LOG   ')
       y=1/x
    case('IDENT ')
       y=1
    case default
       write(*,'(a,a6)')"Error: unexpected LINK function ",linkf  
       get_delta=1
    end select
  end function get_delta
  !---------------------------------------------------
  ! Quick sort algorithm (2 subroutines)
  !---------------------------------------------------
  recursive subroutine qsort(A,n)
    !
    real(kind=dp), intent(inout), dimension(:) :: A
    integer, intent(in) :: n
    integer :: ip
    !
    if (size(A).gt.1) then
       ip=int(size(A)/2) 
       call partition(A,n,ip)
       call qsort(A(1:(ip-1)),(ip-1))
       call qsort(A((ip+1):n),(n-ip))
    endif
  end subroutine qsort

  subroutine partition(A,n,ip)
    real(kind=dp),dimension(:)::A
    integer::n,ip

    real(kind=dp),dimension(n)::B
    real(kind=dp)::pvt
    integer::il,ih
    pvt=A(ip)
    il=1
    ih=n
    do i=1,n
       if (i.eq.ip) cycle
       if (A(i).lt.pvt) then
          B(il)=A(i)
          il=il+1
       else
          B(ih)=A(i)      
          ih=ih-1
       endif
    enddo
    ip=(il+ih)/2
    B(ip)=pvt
    A=B
  end subroutine partition

  function Get_Stat(y,n,mean,var) 
    !--------------------------------------- 
    ! Getting mean and variance
    !--------------------------------------- 
    real(kind=dp),intent(in)::y(:)
    integer,intent(in)::n
    real(kind=dp),intent(out)::mean,var
    integer::Get_Stat
    !
    integer::i
    real(kind=dp)::sx,sx2
    !
    Get_Stat=0
    if (n.eq.0)then
       write(*,*)"Error: empty table => no mean nor variance"
       Get_Stat=1
    else
       sx=0 ; sx2=0 !MAZ
       sx=sum(y(:))
       do i=1,n
          sx2=sx2+y(i)*y(i)
       enddo
       mean=sx/n
       if (n.eq.1) then
          var=0.0
       else
          var=(sx2-sx*sx/n)/(n-1)
       endif
    end if
  end function Get_Stat

  function Get_Scores(p,r,p2g,tau2,scores,pvalues)
    !======================================================================
    !
    ! Compute score and associated p-values
    !
    !======================================================================
    !P       : Object containing parameters
    !r       : Residuals
    !p2g     : Correspondance with id
    !tau2    : Tau2 estimation
    !scores  : Score 
    !pvalues : pvalues
    !======================================================================
    type(parameters)::p

    real(kind=dp),intent(in),dimension(:)::r
    integer,intent(in),dimension(:)::p2g
    real(kind=dp),intent(in)::tau2
    real(kind=dp),intent(out),dimension(:,:)::scores
    real(kind=dp),intent(out),dimension(:,:)::pvalues
    integer::Get_Scores,deb,fin,lg
    ! Permuted scores table
    real(kind=dp),allocatable::pscores(:),gpscores(:)
    integer,allocatable::bscores(:),gbscores(:)
    ! Permuted phenotypes
    real(kind=dp),allocatable::rp(:,:)
    ! Work variables
    real(kind=dp),allocatable::eta(:),delta(:),var(:),W(:),dW(:),rdWZ(:)
    real(kind=dp),allocatable::p_scores(:),a_gamma(:),b_gamma(:)
    integer(kind=1),allocatable::D(:,:,:)
    real(kind=dp)::SX,SX2,EPS,NU
    integer::ipos,i,j,ip,nb_scores,ifrom,ito,nb_pos,iscore,itoNext,ifromNext
    character(len=1000)::Misc !in order to write extra info 
    integer(kind=8)::BP,NSnp,NumFicGeno,PosWithinRun
    integer::CptPos,NbPos,NextLoad,LgLecture,PosRnd
    character(len=3)::CHR
    character(len=100)::Outputname,hFile
    character(len=60)::NameSNP
    Integer::test,OMP_GET_THREAD_NUM,NbPerm,NbRec,PosFin,PosDebut,slice,Nhap,Nphase,Ani
    real(kind=dp)::UnSurTau4
    !
    Get_Scores=0
    allocate(a_gamma(p_p%nb_scores),b_gamma(p_p%nb_scores),bscores(p%nb_perms),pscores(p%nb_perms),&
         &p_scores(p%nb_perms),gbscores(p%nb_perms),gpscores(p%nb_perms),stat=ok)
    scores=0.0   ;    pscores=0.0 ; bscores=0 ;    
    gpscores=1.0 ;     gbscores=0 ; 
    a_gamma=0.0  ;  b_gamma=0.0 ;UnSurTau4=1.0_dp/(tau2*tau2)
    !
    nb_scores=p%scores_pos(2)-p%scores_pos(1)+1
    ! First compute parts that are independent of the position
    allocate(eta(p%nb_rec),delta(p%nb_rec),&
         &var(p%nb_rec),W(p%nb_rec),&
         &dW(p%nb_rec),stat=ok)
    !***********************************************
    ! Assuming beta=0, the GLMM reduces to a GLM g(y) = m + e
    ! ===> Include correction for m here ? <===
    ! rc=Get_Stat(r,p%nb_rec,mean,v2) 
    ! if (rc.ne.0) goto 998
    ! eta=mean
    ! t2=v2
    ! print *,' ===> Residual mean     in residuals model: ',mean
    ! print *,' ===> Residual variance in residuals model: ',t2
    ! ===> Then, new model must be REMLed  <===
    !***********************************************
    !We rename variable (unaccepted as shared by openMP)
    PosDebut=p%scores_pos(1); PosFin=p%scores_pos(2) ;slice=p%slice
    NbPerm=p%nb_perms ;  NbRec=p%nb_rec ; hfile=p%hfile ; PosRnd=p%pos_rnd(2)
    ! Permute phenotypes-----------------------------------------------
    if (p%nb_perms.gt.0) then
       allocate(rp(p%nb_rec,p%nb_perms),stat=ok)       
       rc=Permute(r,p%nb_rec,rp,p%nb_perms)
    endif
    ! Now loop on positions---------------------------------------------
    NbPos=0 ;NumFicGeno=1;NextLoad=p%scores_pos(1);CptPos=PosDebut-1 !MAZ
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(pardim,PosFin,PosDebut,NbPos,NumFicGeno,GenoFileName,CptPos,NextLoad,p2g,r,rp,p_scores,&
    !$OMP SX,SX2,scores,bscores,pscores,a_gamma,b_gamma,gpscores,nbrec,LgLecture,slice,UnSurTau4,DbgMode,nbperm,PosRnd,hfile)

    Bpos:do 
       !$OMP SINGLE      
       CptPos=CptPos+1 ; NbPos=NbPos+1 !Count of really tested position    
       !$OMP END SINGLE
       !$OMP BARRIER
       if(CptPos>PosFin)exit Bpos
       CdtLoad:if(CptPos==NextLoad)then !we must (re)read genotypes     
          !$OMP BARRIER
          !Load data=============================================================
          !If multiple geno file slices should be change accordingly 
          if(allocated(GenoFileName))THEN
             !Find the correct file to read
             !$OMP MASTER             
             FindFile:do
                if(CptPos>=pardim(NumFicGeno,3).and.(CptPos<=pardim(NumFicGeno,4)))exit FindFile
                !$if(DbgMode)write(*,*)"Avant ",pardim(NumFicGeno,3:4),CptPos,NumFicGeno,OMP_GET_THREAD_NUM()
                NumFicGeno=NumFicGeno+1                
                if(NumFicGeno>size(GenoFileName,1))then 
                   write(*,*)"Position ",CptPos," seems incorrect ",NumFicGeno-1
                   stop
                endif
             end do FindFile
             !$OMP FLUSH(NumFicGeno)
             LgLecture=min(slice,(1+pardim(NumFicGeno,4)-CptPos)) !We will read a block of length slice or till the end of file.
             !$OMP CRITICAL
             NextLoad=CptPos+LgLecture !Nextposition to load
             !$OMP END CRITICAL
             !$OMP FLUSH (NextLoad,LgLecture,NumFicGeno)
             !$OMP END MASTER
             !$OMP BARRIER
             if(allocated(D))deallocate(D)   !The trick is to stop all thread to deallocate/allocate them in parallel
             allocate(D(nbrec,LgLecture,2),stat=ok)
             D=0
             !$if(ok.ne.0)write(*,*)" Pb During allocation of D ",ok," thread ",OMP_GET_THREAD_NUM()
             !$OMP BARRIER
             deb=CptPos-pardim(NumFicGeno,3)+1  !Compute beginning end and length of the slice
             fin=NextLoad-1
             lg=1+pardim(NumFicGeno,4)-pardim(NumFicGeno,3)
             !$OMP SINGLE
             rc=Get_Design_Haplos2(GenoFileName(NumFicGeno),deb,fin,lg,D)   !One thread will load data and copy it's value to the other
             !$OMP END SINGLE COPYPRIVATE(D)
             !$OMP BARRIER
             PosWithinRun=0 !Intra slice position set to 0 
          else 
             if ((mod(NbPos,slice).eq.1))then
                !$ if(OMP_GET_THREAD_NUM()==0)write(*,*)' ----> Obtaining new haplotypes from file...'
                ifrom=PosDebut+nbpos-1
                if ((ifrom+slice-1).gt.PosFin) then
                   ito=PosFin
                else
                   ito=ifrom+slice-1
                endif
                LgLecture=1+ito-ifrom
                !$OMP BARRIER
                if(allocated(D))deallocate(D)
                allocate(D(nbrec,LgLecture,2),stat=ok)
                !$OMP BARRIER 
                !$OMP SINGLE 
                rc=Get_Design_Haplos(hfile,PosRnd,ifrom,ito,D)
                !$OMP END SINGLE COPYPRIVATE(D)
                !$OMP BARRIER                
                PosWithinRun=0
                NextLoad=CptPos+LgLecture
             end if
          end if
       endif CdtLoad
       PosWithinRun=PosWithinRun+1
       ! $ OMP BARRIER
       if(Allocated(rdWZ))Deallocate(rdWZ)
       allocate(rdWZ(maxval(D(:,PosWithinRun,:))),stat=ok)
       ! $ OMP BARRIER
       if(DbgMode)write(*,*)"Test of position ",PosWithinRun,size(r,1),size(rdWZ,1)
       ! Compute score
       rdWZ=0
       !$OMP MASTER   
       do Ani=1,size(r,1)
          do Nphase=1,2
             Nhap=D(p2g(Ani),PosWithinRun,Nphase)
             if(Nhap==0 .or.Nhap >size(rdWZ,1))then
                write(*,*)"Pb with ",Ani,p2g(Ani),PosWithinRun,Nphase,Nhap
                stop
             endif
             rdWZ(Nhap)=rdWZ(NHap)+r(Ani)!*dW(i)
          end do
       enddo
       scores(NbPos,1)=0.5*UnSurTau4*dot_product(rdWZ,rdWZ)! /(tau2*tau2)
       !$OMP END MASTER
       !$OMP BARRIER

       !----------------------------------------------------------------------
       ! Compute permutations
       !----------------------------------------------------------------------
       if(DbgMode.and.(mod(NbPos,10)==1))write(*,*)NbPos,"Score = ",scores(NbPos,1)," debut perm"
       !$ if(DbgMode)write(*,*)"thread ",OMP_GET_THREAD_NUM(),size(p_scores,1),PosWithinRun,p2g(1),D(1,1,1)
       !$OMP DO 
       LPerm:do ip=1,size(rp,2)
          rdWZ=0
          do i=1,size(rp,1)
             do Nphase=1,2
                Nhap=D(p2g(i),PosWithinRun,Nphase)
                if(Nhap==0)cycle
                rdWZ(Nhap)=rdWZ(Nhap)+rp(i,ip)!*dW(i)
             end do
          enddo
          !Store solution for permutation
          p_scores(ip)=0.5*UnSurTau4*dot_product(rdWZ,rdWZ)
          !fin
       enddo LPerm
       !$OMP END DO
       ! $ OMP BARRIER
       !$OMP single
       SX=0    ;  SX2=0 !MAZ

       !$OMP FLUSH(SX,SX2)
       !$OMP END SINGLE 
       !$OMP DO REDUCTION (+:SX) REDUCTION (+:SX2)  
       do  ip=1,size(p_scores,1)
          SX=SX+p_scores(ip)
          SX2=SX2+p_scores(ip)*p_scores(ip)
       end do
       !$OMP END DO

       !$OMP WORKSHARE
       scores(NbPos,3)=count(p_scores(:).gt.scores(NbPos,1))
       !$OMP END WORKSHARE

       !$OMP WORKSHARE
       where(pscores<p_scores)bscores=CptPos
       where(pscores<p_scores)pscores=p_scores
       !$OMP END WORKSHARE
       if(DBGmode.and.mod(NbPos,1000)==0)write(*,*)"Debut Gamma ",NbPos,CptPos,SX,nbperm,unsurtau4
       !----------------------------------------------------------------------
       ! Compute Gamma parameters
       !----------------------------------------------------------------------
       if (nbperm.gt.1) then
          !$OMP MASTER
          EPS=SX/nbperm
          NU=(SX2-SX*SX/nbperm)/(nbperm-1)
          a_gamma(NbPos)=EPS*EPS/NU          
          b_gamma(NbPos)=NU/EPS          
          scores(NbPos,2)=1_dp-gammad(scores(NbPos,1)/b_gamma(NbPos),a_gamma(NbPos))
          if(scores(NbPos,2)<10.0**(-16))scores(NbPos,2)=10.0**(-16)
          scores(NbPos,3)=(1.00+scores(NbPos,3))/(1.00*nbperm)
          !$OMP END MASTER 
          !$OMP BARRIER
          !$OMP FLUSH (a_gamma,b_gamma,scores,p_scores)          
          !$OMP DO
          do ip=1,nbperm
             p_scores(ip)=1_dp-gammad(p_scores(ip)/b_gamma(NbPos),a_gamma(NbPos))
             if(p_scores(ip)<10.0**(-16))p_scores(ip)=10.0**(-16)
          enddo
          !$OMP END DO

          !$OMP WORKSHARE
          gpscores(nbperm)=minval(p_scores)
          !$OMP END WORKSHARE
       endif
       if(DBGmode)write(*,*)"END of Gamma ",NbPos,CptPos
    enddo Bpos
    if(DBGmode)write(*,*)"Sortie Boucle position"
    if(allocated(D))deallocate(D)      
    !$if(DBGmode)write(*,*)"Deallocation D : Thread ",OMP_GET_THREAD_NUM()
    !$OMP END PARALLEL

    !--------------------------------------------------   
    ! Compute p_values
    !--------------------------------------------------   
    if(DBGmode)write(*,*)" P values computation"
    if(ExportPerm)then
       open(25,file=trim(p_p%aname)//'.perms')
       write(25,fmt='("Before sorting scores:")')
       write(25,fmt='("Perms Score            Bestp")')
       do i=1,p%nb_perms
          write(25,fmt='(i5,1x,f16.5,1x,i7)') i,pscores(i),bscores(i)
       enddo
    endif
    call qsort(pscores,p%nb_perms)
    if(ExportPerm)then
       write(25,fmt='("After sorting:")')
       do i=p%nb_perms,1,-1
          write(25,fmt='(i5,1x,f16.5)') i,pscores(i)
       enddo
       write(25,fmt='("Before sorting gamma scores:")')
       write(25,fmt='("Perms Score            Bestp")')
       do i=1,p%nb_perms
          write(25,fmt='(i5,1x,f17.14,1x,i7)') i,gpscores(i),gbscores(i)
       enddo
    endif
    call qsort(gpscores,p%nb_perms)
    if(ExportPerm)then
       write(25,fmt='("After sorting:")')
       do i=p%nb_perms,1,-1
          write(25,fmt='(i5,1x,f17.14)') i,gpscores(i)
       enddo
       close(25)
    endif
    open(21,file=trim(p_p%aname)//'.out')
    do ipos=p%scores_pos(1),p%scores_pos(2)
       iscore=ipos-p%scores_pos(1)+1
       pvalues(iscore,1)=0.0
       !       do i=p%nb_perms,1,-1
       !          if (scores(ipos,1).gt.pscores(i)) exit
       !          pvalues(ipos,1)=pvalues(ipos,1)+1.0/p%nb_perms ! rank of score corrected for multiple testing - incorrect score is local
       !       enddo
       pvalues(iscore,2)=1.0/p%nb_perms
       do i=1,p%nb_perms
          if (scores(iscore,2).lt.gpscores(i))exit
          pvalues(iscore,2)=pvalues(iscore,2)+1.0/p%nb_perms ! rank of gamma p-value corrected for multiple-testing
       enddo
       write(21,'(i12,3(1x,f18.6),3(3x,e12.4))')ipos,a_gamma(iscore),&
            &b_gamma(iscore),scores(iscore,1),scores(iscore,2),scores(iscore,3),pvalues(iscore,2)                  
    enddo
    close(21)
    !Creation of a .assoc output 
    if(DoReadMap)THEN
       write(Outputname,'(a,".assoc")')adjustl(trim(p_p%aname))
       !======================================================================
       open(12,file=Outputname,action="WRITE")
       open(11,file=filename,action="READ")
       read(11,'(a1000)')Misc
       if(adjustl(Misc(1:1))=="#")THEN
          write(12,'(a)')adjustl(trim(Misc))
       else 
          write(12,'(a)')"# CHR SNP BP SCORE P P2 P3"
          rewind(11)
       endif
       !======================================================================
       NSnp=0
       do 
          read(11,'(a,a,i10)',iostat=io)CHR,NameSNP,BP
          if(io/=0)exit
          NSnp=NSnp+1
          write(12,'(a,a,i0,f18.6,3(3x,e12.4),a1000)')CHR,adjustl(trim(NameSNP)),BP,scores(NSNP,1:3),pvalues(NSNP,2)
       end do
       close(12); close(11)
    endif
    goto 1000
1000 continue
    if(allocated(rp))deallocate(rp)
    if(allocated(dW))deallocate(dW)
    if(allocated(W))deallocate(W)
    if(allocated(var))deallocate(var)
    if(allocated(delta))deallocate(delta)
    if(allocated(eta))deallocate(eta)
    if(allocated(D))deallocate(D)
    if(allocated(gpscores))deallocate(gpscores)
    if(allocated(gbscores))deallocate(gbscores)
    if(allocated(p_scores))deallocate(p_scores)
    if(allocated(pscores))deallocate(pscores)
    if(allocated(bscores))deallocate(bscores)  

  end function Get_Scores

  function get_D(linkf,mu,D,n)
    !--------------------------------------- 
    ! Getting D (diagonal) matrix
    !  tag provides the link function
    !--------------------------------------- 
    character(len=6),intent(in)::linkf
    integer,intent(in)::n
    real(kind=dp),intent(in)::mu(n)
    real(kind=dp),intent(out)::D(n)
    integer::get_D
    !
    get_D=0
    select case(linkf)
    case('LOGIT ')
       D=mu*(1-mu)
    case('LOG   ')
       D=mu
    case('IDENT ')
       D=1.0
    case Default
       write(*,'(a,a6)')"Error: unexpected LINK function ",linkf  
       get_D=1 
    end select
  end function get_D

  function get_var(distr,mu,t2,var,n)
    !--------------------------------------- 
    ! Variance functions
    !  distr provides the distribution function
    !--------------------------------------- 
    character(6),intent(in)::distr
    integer,intent(in)::n
    real(kind=dp),intent(in)::mu(n)
    real(kind=dp),intent(inout)::t2
    real(kind=dp),intent(out)::var(n)
    integer::get_var
    !
    get_var=0
    select case(distr)
    case('BINOM')
       t2=1.0
       var=mu*(1-mu)
    case('POISS ')
       t2=1.0
       var=mu
    case('NORMAL')
       var=1.0
    case Default
       write(*,'(a,a6)')"Error: unexpected DISTRIBUTION function ",distr  
       get_var=1
    end select
  end function get_var

  function get_mu(linkf,x,y,n)
    !--------------------------------------- 
    ! Probability functions
    ! (i.e. inverse of link function)
    !  tag provides the link function
    !  y=inv_link(x) 
    !--------------------------------------- 
    !
    character(6),intent(in)::linkf
    integer,intent(in)::n
    real(kind=dp),intent(in)::x(n)
    real(kind=dp),intent(out)::y(n)
    integer::get_mu
    !
    get_mu=0
    select case(linkf)
    case('LOGIT ')
       y=exp(x)/(1+exp(x))
    case('LOG   ')
       y=exp(x)
    case('IDENT ')
       y=x
    case Default
       write(*,'(a,a6)')"Error: unexpected LINK function ",linkf  
       get_mu=1
    end select
  end function get_mu


  function Compute_REML(y,X,Z,corr,use_rnd,n,nft,nrt,nf,nr,dim_fix,dim_rnd, &
       linkf,distr,b,u,tau2,sigma2,max_iter,max_tol,method)
    !--------------------------------------- 
    ! REML estimation of VC and of b and u 
    ! This version is limited to 1 random effect
    !--------------------------------------- 
    character(6),intent(in)::linkf,distr
    integer,intent(in)::use_rnd(:)
    integer::dim_fix(:),dim_rnd(:)
    integer,intent(in)::n,nft,nrt,nf,nr,max_iter
    real(kind=dp),intent(in)::y(n),X(n,nft),Z(n,nrt),max_tol
    real(kind=dp),intent(inout)::b(nft),u(nrt),tau2,sigma2(nr)
    type(correl),intent(in)::corr
    character(2),intent(in)::method
    integer::Compute_REML
    ! Work variables
    real(kind=dp)::eta(n),mu(n),var(n),D(n),R(n),Ri(n),e(n),f(n),g(n)
    real(kind=dp)::ys(n)
    real(kind=dp),allocatable::MMEl(:,:),MMEr(:),MMEf(:),MMEg(:)
    real(kind=dp),allocatable,dimension(:,:)::MMElp
    real(kind=dp),allocatable::oldsol(:),sol2(:),solf(:),solg(:)
    real(kind=dp)::sol(nft+nrt)
    real(kind=dp)::new_tau2,glmm_tau2  ! Error variance
    real(kind=dp)::new_sigma2(nr),glmm_sigma2  ! Random effects variances 
    real(kind=dp)::det,trace,AI(2,2),EM(2,2),MIX(2,2),dlr,dlg
    logical::reml_done,reml_ok,new_glmm
    real(kind=dp)::reml_diff,elt1,elt2,glmm_diff,sol_diff
    integer::reml_iter,n1,n2,info,glmm_iter,tot_iter
    integer::i,ntr,nftr,nrtr,base,irnd,mixing
    !
    Compute_REML=0
    !---------------------------------------------------------------------!
    ! Allocate space for MME
    !---------------------------------------------------------------------!
    ntr=0
    do i=1,nf
       if(dim_fix(i)>0)ntr=ntr+dim_fix(i)
       if(dim_fix(i)==0 .and. sum(X(:,i))/=0.00)ntr=ntr+1 ! regression
    enddo
    nftr=ntr
    ntr=0
    do i=1,nr
       if (use_rnd(i).eq.1) then
          ntr=ntr+dim_rnd(i)
          irnd=i
       endif
    enddo
    nrtr=ntr
    ntr=nftr+nrtr 
    allocate(MMElp(ntr,ntr),MMEl(ntr,ntr),MMEr(ntr),MMEf(ntr),MMEg(ntr),&
         &oldsol(ntr),sol2(ntr),solf(ntr),solg(ntr),stat=ok)
    if (ok.ne.0 ) goto 106
    MMElp = MMEl
    !---------------------------------------------------------------------!
    ! Set initial values b(0), u(0)
    !---------------------------------------------------------------------!
    !sol=0
    sol(1:nft)=b                     ! Initial values for fixed effects
    sol((nft+1):(nft+nrt))=u         ! Initial values for random effects
    !---------------------------------------------------------------------!
    ! Start REML iterations
    !---------------------------------------------------------------------!
    reml_done=.FALSE.
    reml_iter=0;glmm_iter=0;tot_iter=0;new_glmm=.true.
    glmm_sigma2=sigma2(1);glmm_tau2=tau2
    REML1:do while (.not. reml_done)
       if(new_glmm)then
          glmm_iter=glmm_iter+1
          reml_iter=1
          tot_iter=tot_iter+1
       else
          reml_iter=reml_iter+1
          tot_iter=tot_iter+1
       endif
       if (reml_iter.gt.max_iter) then
          reml_done=.TRUE.
          reml_ok=.FALSE.
       else
          !---------------------------------------------------------------------!
          ! 0) Save solutions (i.e. b and u)
          !---------------------------------------------------------------------!
          NewGlmm:if(new_glmm)then ! only when new_glmm
             oldsol=sol
             !---------------------------------------------------------------------!
             ! 1) Compute eta
             !---------------------------------------------------------------------!
             base=0
             do i=1,nf
                start_indices(i)=base
                if(dim_fix(i)==0 .and. sum(X(:,i))/=0.00)base=base+1 ! regression only for fixed effects
                if(dim_fix(i)>0)base=base+dim_fix(i)
             enddo
             do i=1,nr
                if(use_rnd(i)==0)cycle
                start_indices(nf+i)=base
                base=base+dim_rnd(i)
             enddo

             eta=0.00
             do i=1,n
                do j=1,nf
                   if(dim_fix(j)>0 .and. X(i,j)==0)cycle 
                   k=start_indices(j)+X(i,j)
                   if(dim_fix(j)==0)k=start_indices(j)+1
                   elt1=1.00
                   if(dim_fix(j)==0)elt1=X(i,j)
                   !write(*,*)"elt1 solk",elt1,sol(k)
                   eta(i)=eta(i)+elt1*sol(k)
                enddo
                do j=1,nr
                   if(use_rnd(j)==0)cycle
                   k=start_indices(nf+j)+Z(i,j)
                   elt1=1.00 ! no regression for random effects
                   !write(*,*)" solk",sol(k)
                   eta(i)=eta(i)+elt1*sol(k)
                enddo
             enddo
             !---------------------------------------------------------------------!
             ! 2) Compute mu=inv_g(eta)
             !---------------------------------------------------------------------!
             rc=get_mu(linkf,eta,mu,n)
             if (rc.gt.0) goto 102
             !---------------------------------------------------------------------!
             ! 3) Compute var function (without tau squared term)
             !---------------------------------------------------------------------!
             rc=get_var(distr,mu,tau2,var,n)
             if (rc.gt.0) goto 103
             !---------------------------------------------------------------------!
             ! 4) Compute work matrices
             !---------------------------------------------------------------------!
             !     => D(k) = d[mu(k)]/d[eta(k)]
             rc=get_D(linkf,mu,D,n) 
             if (rc.gt.0) goto 104
             !   !     => R(k) = tau2(k)*var(k)
             !   R=tau2*var
             !   !     => inv(R(k))
             !   Ri=1/R
             !   !     => Ys(k) = Y - mu(k) +D(k)*eta(k)
             Ys=(Y-mu+D*eta)/sqrt(var)
             new_glmm=.false.

          endif NewGlmm ! new_glmm
          !     => R(k) = tau2(k)*var(k)
          R=tau2*var
          !     => inv(R(k))
          Ri=1/R
          !---------------------------------------------------------------------!
          ! 5) Compute other variances using MME
          !---------------------------------------------------------------------!
          MMEl=0.0;MMEr=0.0
          do k=1,n
             ! X'*Ri*X
             do i=1,nf
                if(dim_fix(i)>0 .and. X(k,i)==0)cycle 
                do j=1,nf
                   if(dim_fix(j)>0 .and. X(k,j)==0)cycle 
                   n1=start_indices(i)+X(k,i)
                   n2=start_indices(j)+X(k,j)
                   elt1=1.00;elt2=1.00
                   if(dim_fix(i)==0)then
                      n1=start_indices(i)+1
                      elt1=X(k,i)
                   endif
                   if(dim_fix(j)==0)then
                      n2=start_indices(j)+1
                      elt2=X(k,j)
                   endif
                   !       MMEl(n1,n2)=MMEl(n1,n2)+elt1*elt2*D(k)*D(k)*Ri(k)
                   MMEl(n1,n2)=MMEl(n1,n2)+elt1*elt2*D(k)/sqrt(var(k))*D(k)/sqrt(var(k))*1/tau2
                enddo
                !     MMEr(n1)=MMEr(n1)+elt1*D(k)*Ri(k)*Ys(k)
                MMEr(n1)=MMEr(n1)+elt1*D(k)/sqrt(var(k))*1/tau2*Ys(k)
             enddo
             ! X'*Ri*Z
             do i=1,nf
                if(dim_fix(i)>0 .and. X(k,i)==0)cycle 
                do j=1,nr
                   if(use_rnd(j)==0)cycle
                   n1=start_indices(i)+X(k,i)
                   n2=start_indices(nf+j)+Z(k,j)
                   elt1=1.00;elt2=1.00
                   if(dim_fix(i)==0)then
                      n1=start_indices(i)+1
                      elt1=X(k,i)
                   endif
                   MMEl(n1,n2)=MMEl(n1,n2)+elt1*elt2*D(k)/sqrt(var(k))*D(k)/sqrt(var(k))*1/tau2
                   MMEl(n2,n1)=MMEl(n2,n1)+elt1*elt2*D(k)/sqrt(var(k))*D(k)/sqrt(var(k))*1/tau2
                enddo
             enddo
             ! Z'*Ri*Z
             do i=1,nr
                if(use_rnd(i)==0)cycle
                do j=1,nr
                   if(use_rnd(j)==0)cycle
                   n1=start_indices(nf+i)+Z(k,i)
                   n2=start_indices(nf+j)+Z(k,j)
                   elt1=1.00;elt2=1.00
                   MMEl(n1,n2)=MMEl(n1,n2)+elt1*elt2*D(k)/sqrt(var(k))*D(k)/sqrt(var(k))*1/tau2 
                enddo
                elt1=1.00
                MMEr(n1)=MMEr(n1)+elt1*D(k)/sqrt(var(k))*1/tau2*Ys(k)
             enddo
          enddo

          ! Adding inverse variance matrices
          base=nftr
          do i=1,nr
             if (use_rnd(i).eq.0) cycle
             if (corr%is_diagonal) then
                do j=1,dim_rnd(i)
                   MMEl(base+j,base+j)=MMEl(base+j,base+j)+corr%inv_cor(j,1)/sigma2(i)
                enddo
             else
                MMEl((base+1):(base+dim_rnd(i)),(base+1):(base+dim_rnd(i))) = &
                     MMEl((base+1):(base+dim_rnd(i)),(base+1):(base+dim_rnd(i))) + &
                     corr%inv_cor/sigma2(i)
             endif
             base=base+dim_rnd(i)
          enddo

          ! Find a (generalized) inverse of MMEl
#if LAPACK
          call DPOTRF('U',ntr,MMEl,ntr,info)
          call DPOTRI('U',ntr,MMEl,ntr,info)
          do j=2,ntr
             do k=1,j
                MMEl(j,k)=MMEl(k,j)
             enddo
          enddo
#else   
          call invert(MMEl,ntr)
#endif 
          ! Compute solutions vector
          sol=matmul(MMEl,MMEr)

          !-------------------------------------------
          ! Obtain new estimates of the VC using method
          !-------------------------------------------
          CdtAI:if (method.eq.'AI') then
             !-------------------------------------------
             !(AI-REML)
             ! newv = v - inv(AI)*score
             ! With BINOM & POISS, ve=1
             !      and newva = va - inv(AI11)*score_a
             ! With NORMAL, both va and ve are recomputed
             !-------------------------------------------

             ! compute residuals

             e=Ys
             do i=1,n
                do j=1,nf
                   if(dim_fix(j)>0 .and. X(i,j)==0)cycle 
                   k=start_indices(j)+X(i,j)
                   elt1=1.00
                   if(dim_fix(j)==0)then
                      k=start_indices(j)+1
                      elt1=X(i,j)
                   endif
                   !       e(i)=e(i)-elt1*D(i)*sol(k)
                   e(i)=e(i)-elt1*D(i)/sqrt(var(i))*sol(k)
                enddo
                elt1=1.00 ! no regression for random effects 
                do j=1,nr
                   if(use_rnd(j)==0)cycle
                   k=start_indices(nf+j)+Z(i,j)
                   !       e(i)=e(i)-elt1*D(i)*sol(k)
                   e(i)=e(i)-elt1*D(i)/sqrt(var(i))*sol(k)
                enddo
             enddo

             ! compute dlr

             dlr=0.00
             do k=1,n
                ! X*Cii*X'
                do i=1,nf
                   if(dim_fix(i)>0 .and. X(k,i)==0)cycle 
                   do j=1,nf
                      if(dim_fix(j)>0 .and. X(k,j)==0)cycle 
                      n1=start_indices(i)+X(k,i)
                      n2=start_indices(j)+X(k,j)
                      elt1=1.00;elt2=1.00
                      if(dim_fix(i)==0)then
                         n1=start_indices(i)+1
                         elt1=X(k,i)
                      endif
                      if(dim_fix(j)==0)then
                         n2=start_indices(j)+1
                         elt2=X(k,j)
                      endif
                      dlr=dlr+elt1*elt2*D(k)/sqrt(var(k))*D(k)/sqrt(var(k))*MMEl(n1,n2)
                   enddo
                enddo
                ! X*Cii*Z'
                do i=1,nf
                   if(dim_fix(i)>0 .and. X(k,i)==0)cycle 
                   do j=1,nr
                      if(use_rnd(j)==0)cycle
                      n1=start_indices(i)+X(k,i)
                      n2=start_indices(nf+j)+Z(k,j)
                      elt1=1.00;elt2=1.00
                      if(dim_fix(i)==0)then
                         n1=start_indices(i)+1
                         elt1=X(k,i)
                      endif
                      dlr=dlr+elt1*elt2*D(k)/sqrt(var(k))*D(k)/sqrt(var(k))*MMEl(n1,n2)
                      dlr=dlr+elt1*elt2*D(k)/sqrt(var(k))*D(k)/sqrt(var(k))*MMEl(n1,n2)
                   enddo
                enddo
                ! Z*Cii*Z'
                do i=1,nr
                   if(use_rnd(i)==0)cycle
                   do j=1,nr
                      if(use_rnd(j)==0)cycle
                      n1=start_indices(nf+i)+Z(k,i)
                      n2=start_indices(nf+j)+Z(k,j)
                      elt1=1.00;elt2=1.00
                      dlr=dlr+elt1*elt2*D(k)/sqrt(var(k))*D(k)/sqrt(var(k))*MMEl(n1,n2)
                   enddo
                enddo
             enddo

             dlr=-(dlr+dot_product(e,e))/(tau2**2)+n/tau2

             ! compute dlg

             trace=0
             if (corr%is_diagonal) then
                do i=1,nrtr
                   trace=trace+corr%inv_cor(i,1)*MMEl(nftr+i,nftr+i)
                enddo
             else
                do i=1,nrtr
                   do j=1,nrtr
                      trace=trace+corr%inv_cor(i,j)*MMEl(nftr+j,nftr+i)
                   enddo
                enddo
             endif

             dlg=nrtr/sigma2(irnd)-(trace+dot_product(sol(nftr+1:nftr+nrtr), &
                  matmul(corr%inv_cor,sol(nftr+1:nftr+nrtr)))) &
                  /(sigma2(irnd)**2)

             ! compute f vector for random effect

             f=0.00
             do i=1,n
                do j=1,nr
                   if(use_rnd(j)==0)cycle
                   k=start_indices(nf+j)+Z(i,j)
                   !                   write(*,*)" var et sol",var(i),sol(k)
                   f(i)=f(i)+D(i)/sqrt(var(i))*sol(k)
                   !      f(i)=f(i)+D(k)*sol(k)
                enddo
             enddo

             f=f/sigma2(irnd)

             ! use f vector to compute RHS

             MMEf=0.00
             do i=1,n
                do j=1,nf
                   if(dim_fix(j)>0 .and. X(i,j)==0)cycle 
                   k=start_indices(j)+X(i,j)
                   elt1=1.00
                   if(dim_fix(j)==0)then
                      k=start_indices(j)+1
                      elt1=X(i,j)
                   endif
                   !     MMEr(k)=MMEr(k)+elt1*D(i)*Ri(i)*f(i)
                   MMEf(k)=MMEf(k)+elt1*D(i)/sqrt(var(i))*f(i)
                enddo
                do j=1,nr
                   if(use_rnd(j)==0)cycle
                   k=start_indices(nf+j)+Z(i,j)
                   elt1=1.00 ! no regression for random effects
                   !     MMEr(k)=MMEr(k)+elt1*D(i)*Ri(i)*f(i)
                   MMEf(k)=MMEf(k)+elt1*D(i)/sqrt(var(i))*f(i)
                enddo
             enddo
             MMEf=MMEf/tau2

             ! obtain solutions

             solf=matmul(MMEl,MMEf)

             ! obtain f vector for residals => g

             g=e/tau2

             ! construct new RHS

             MMEg=0.00
             do i=1,n
                do j=1,nf
                   if(dim_fix(j)>0 .and. X(i,j)==0)cycle
                   k=start_indices(j)+X(i,j)
                   elt1=1.00
                   if(dim_fix(j)==0)then
                      k=start_indices(j)+1
                      elt1=X(i,j)
                   endif
                   !      MMEr(k)=MMEr(k)+elt1*D(i)*Ri(i)*g(i)
                   MMEg(k)=MMEg(k)+elt1*D(i)/sqrt(var(i))*g(i)
                enddo
                do j=1,nr
                   if(use_rnd(j)==0)cycle
                   k=start_indices(nf+j)+Z(i,j)
                   elt1=1.00 ! no regression for random effects
                   !      MMEr(k)=MMEr(k)+elt1*D(i)*Ri(i)*g(i)
                   MMEg(k)=MMEg(k)+elt1*D(i)/sqrt(var(i))*g(i)
                enddo
             enddo
             MMEg=MMEg/tau2

             ! obtain solutions

             solg=matmul(MMEl,MMEg)

             ! compute AI

             AI=0.d00

             do k=1,ntr
                AI(1,1)=AI(1,1)-solf(k)*MMEf(k)
                AI(1,2)=AI(1,2)-solg(k)*MMEf(k)
                AI(2,1)=AI(2,1)-solf(k)*MMEg(k)
                AI(2,2)=AI(2,2)-solg(k)*MMEg(k)
             enddo

             AI(1,1)=AI(1,1)+dot_product(f,f)/tau2
             AI(1,2)=AI(1,2)+dot_product(f,g)/tau2
             AI(2,1)=AI(2,1)+dot_product(g,f)/tau2
             AI(2,2)=AI(2,2)+dot_product(g,g)/tau2

             EM=0.d00

             EM(1,1)=(1.0*nrtr)/(sigma2(1)**2)
             EM(2,2)=(1.0*n)/(tau2**2)

             mixing=0;new_tau2=-1.00;new_sigma2(irnd)=-1.00
             do while(new_tau2<0.00 .or. new_sigma2(irnd)<0.00)
                MIX=(1.00-0.001*mixing)*AI+0.001*mixing*EM
                mixing=mixing+1                
                select case(distr)
                case('BINOM ','POISS')
                   new_tau2=tau2     ! tau2 remains fixed to 1
                   new_sigma2(irnd)=sigma2(irnd)-dlg/MIX(1,1)
                   !write(*,*)"Sigma2",sigma2(irnd),dlg,MIX(1,1)
                case('NORMAL')
                   det=MIX(1,1)*MIX(2,2)-MIX(1,2)*MIX(1,2)
                   new_sigma2(irnd)=sigma2(irnd)-(MIX(2,2)*dlg-MIX(1,2)*dlr)/det
                   new_tau2=tau2-(MIX(1,1)*dlr-MIX(1,2)*dlg)/det
                case Default
                   Compute_REML=4
                   write(*,'(2a)')"Error: unknown distribution function ",trim(adjustl(distr))
                   goto 1000
                end select
             enddo

             if(new_sigma2(irnd)<0.000001)new_sigma2(irnd)=0.000001
             if(new_tau2<0.000001)new_tau2=0.000001

          endif CdtAI ! à retirer sinon   

          ! Check convergence

          sol_diff=sqrt(dot_product(sol(1:nftr)-oldsol(1:nftr),sol(1:nftr)-oldsol(1:nftr)))
          reml_diff=((new_tau2-tau2)**2+(new_sigma2(1)-sigma2(1))**2)/(new_tau2**2+new_sigma2(1)**2)
          glmm_diff=((new_tau2-glmm_tau2)**2+(new_sigma2(1)-glmm_sigma2)**2)/(new_tau2**2+new_sigma2(1)**2)
          !write(*,*)new_tau2,tau2,sigma2(1),new_sigma2(1),sigma2(1),glmm_sigma2
          ! print*,reml_diff,glmm_diff          
          if (reml_diff.lt.max_tol) then
             if(glmm_diff.lt.max_tol .and. sol_diff<0.00001) then
                !     if(glmm_iter.eq.20) then
                reml_done=.TRUE.
                reml_ok=.TRUE.
             else ! start new glmm iteration
                new_glmm=.true.
                glmm_sigma2=new_sigma2(1);glmm_tau2=new_tau2
             endif
          endif
          ! New VC estimates
          tau2=new_tau2
          sigma2=new_sigma2    
          write(*,'(a,3i6," : ",2g16.9,i6)')"=> Convergence at the end of iter",glmm_iter,&
               &reml_iter,tot_iter,glmm_diff,reml_diff,(mixing-1)
          write(*,'("   Residual variance = ",2f20.5)')tau2
          do i=1,nr
             if (use_rnd(i).eq.1) print '("   Random variance ",i1," = ",f20.5)',i,sigma2(i)
          enddo
       endif ! end AI/EM
    enddo REML1 ! end REML

    if (.not. reml_ok) Compute_REML=-1 
    ! Report solutions
    b=sol(1:nftr)
    u=sol(nftr+1:nftr+nrtr)

    print '("=> REML estimates after ",2i3," iterations ")',glmm_iter,reml_iter
    do i=1,nr
       if (use_rnd(i).eq.1) print '("   variance random effect ",i1," = ",f10.5)',i,new_sigma2(i)
    enddo
    print '("   residual variance = ",f20.5)',new_tau2
    goto 1000
102 continue
    Compute_REML=1
    print *,'Error: problems while computing probability functions'
    goto 1000
103 continue
    Compute_REML=2
    print *,'Error: problems while computing variance functions'
    goto 1000
104 continue
    Compute_REML=3
    print *,'Error: problems while computing D matrices'
    goto 1000
106 continue
    Compute_REML=5
    print *,'Error: problem while allocating space'
    goto 1000
1000 continue
    if (allocated(sol2)) deallocate(sol2)
    if (allocated(oldsol)) deallocate(oldsol)
    if (allocated(MMEr)) deallocate(MMEr)
    if (allocated(MMEl)) deallocate(MMEl)
  end function Compute_REML
  !---------------------
  ! Own inversion matrix
  !---------------------
  subroutine invert(A,n)
    !
    integer,intent(in)::n 
    real(kind=dp),dimension(:,:)::A
    !
    real(kind=dp)::B(n,n),tol,pivot,sol
    integer::index(n),i,j,ipvt,itemp
    ! Initialize B
    B=0
    do i=1,n
       B(i,i)=1
    enddo
    ! Tolerance
    tol=1.0d-10
    ! Initialize index (no permuted rows)
    do i=1,n
       index(i)=i
    enddo
    ! Gaussian elimination
    do j=1,n-1
       ! Look for largest pivot
       ipvt=j
       do i=j+1,n
          if (abs(A(index(i),j)).gt.abs(A(index(ipvt),j))) ipvt=i
       enddo
       ! Check for null pivot
       if (abs(A(index(ipvt),j)).lt.tol) print *,'Matrix is singular !'
       ! Pivot lines index(j) and index(ipvt)
       itemp=index(j)
       index(j)=index(ipvt)
       index(ipvt)=itemp
       ! Compute new elements
       do i=j+1,n
          pivot=A(index(i),j)/A(index(j),j)
          A(index(i),:)=A(index(i),:)-pivot*A(index(j),:)
          B(index(i),:)=B(index(i),:)-pivot*B(index(j),:)      
       enddo
    enddo
    ! A is now triangularized. Back solve    
    do i=n,1,-1
       do j=1,n
          sol=B(index(i),j)
          do k=i+1,n
             sol=sol-A(index(i),k)*B(index(k),j)
          enddo
          B(index(i),j)=sol/A(index(i),i)
       enddo
    enddo
    A=B
    return
  end subroutine invert
  !--------------------------------------- 
  ! Obtain initial values for b and u
  !--------------------------------------- 
  function Initialize_Solutions(p,b,nft,u,nrt,fb,fu)
    !
    type(parameters)::p
    integer,intent(in)::nrt,nft
    real(kind=dp), intent(out)::b(nft),u(nrt)
    character(100),intent(in)::fb,fu
    !
    character(30)::fname
    integer::Initialize_Solutions
    integer::n,base,ifound
    real(kind=dp)::v
    !
    Initialize_Solutions=0
    b=0
    u=0
    if (fb.ne.'none') then
       open(11,file=fb,err=100)
10     read(11,*,end=20) fname,n,v
       if (n.lt.0) goto 102
       base=0
       ifound=0
       do i=1,p%nb_fix
          if (p%name_fix(i).eq.fname) then
             ifound=1
          else
             base=base+p%dim_fix(i)
          endif
          if (ifound.eq.1) exit
       enddo
       if (ifound.eq.0) goto 101
       if (n.gt.p%dim_fix(i)) goto 102
       b(base+n)=v
       goto 10
20     continue
       close(11)
       goto 200
100    continue
       print '("Error: problem while opening initial values file ",a30)',fb
       Initialize_Solutions=1
       goto 200
101    continue
       print '("Error: cannot find fixed effect ",a30)',fname
       Initialize_Solutions=2
       goto 200
102    continue
       print '("Error: invalid index ",i5,"for fixed effect ",a30)',n,fname
       Initialize_Solutions=3
       goto 200
200    continue  
    endif
    if (fu.ne.'none') then
       open(11,file=fu,err=300)
30     read(11,*,end=40) fname,n,v
       if (n.lt.0) goto 302
       base=0
       ifound=0
       do i=1,p%nb_rnd
          if (p%name_rnd(i).eq.fname) then
             ifound=1
          else
             base=base+p%dim_rnd(i)
          endif
          if (ifound.eq.1) exit
       enddo
       if ((ifound.eq.0).and.(fname.ne.'haplo')) goto 301
       if ((ifound.eq.0).and.(fname.eq.'haplo')) goto 30
       if (n.gt.p%dim_rnd(i)) goto 302
       u(base+n)=v
       goto 30
40     continue
       close(11)
       goto 400
300    continue
       print '("Error: problem while opening initial values file ",a30)',fu
       Initialize_Solutions=4
       goto 400
301    continue
       print '("Error: cannot find random effect ",a30)',fname
       Initialize_Solutions=5
       goto 400
302    continue
       print '("Error: invalid index ",i5,"for random effect ",a30)',n,fname
       Initialize_Solutions=6
       goto 400
400    continue  
    endif
  end function Initialize_Solutions
  !---------------------------------------------------------------------!
  ! Compute correlation matrix 
  ! This routine is invoked to compute (inv) correlation matrices
  ! for the polygenic random effect, either using haplotypes (to infer R)
  ! or using pedigree (to infer R or Rinv)
  !  dir => 1 = direct, -1 = inverse
  !  tag => 'HAPLOTYPES': use haplotypic information to infer correlation
  !         'PEDIGREE  ': use pedigree information to infer correlation
  !---------------------------------------------------------------------!
  function Compute_Inv_Corr(p,R,dir,tag)
    !
    type(parameters)::p
    integer,intent(in)::dir
    real(kind=dp),dimension(:,:)::R
    character(10),intent(in)::tag
    integer::Compute_Inv_Corr
    !
    integer(kind=1),allocatable::haplos(:,:,:)
    integer,allocatable::dummy(:)
    integer::ok,i,j,k,nb,istart,io 
    real(kind=dp)::Rmin
    !
    Compute_Inv_Corr=0
    if (tag.eq.'HAPLOTYPES') then
       if (dir.eq.-1) goto 998
       allocate(haplos(p%nb_rec,2,p%nb_pos),stat=ok)
       if (ok.ne.0) goto 997
       ! Look for starting position in file
       istart=p%pos_rnd(2)
       ! Read haplotypes file
       allocate(dummy(istart-1),stat=ok)
       if (ok.ne.0) goto 997
       open(20,file=p%hfile,err=995)
       nb=0
       do 
          read(20,*,iostat=io)(dummy(i),i=1,istart-1), &
               (haplos(nb+1,1,i),haplos(nb+1,2,i),i=1,p%nb_pos)
          if(io/=0)exit
          nb=nb+1
       end do
       close(20)
       deallocate(dummy)
       ! Use haplotypes to build correlation matrix
       Rmin=1.0
       do i=1,nb
          do j=i,nb
             R(i,j)=0.0
             do k=1,p%nb_pos
                if (haplos(i,1,k).eq.haplos(j,1,k)) R(i,j)=R(i,j)+0.25
                if (haplos(i,1,k).eq.haplos(j,2,k)) R(i,j)=R(i,j)+0.25
                if (haplos(i,2,k).eq.haplos(j,1,k)) R(i,j)=R(i,j)+0.25
                if (haplos(i,2,k).eq.haplos(j,2,k)) R(i,j)=R(i,j)+0.25
             enddo
             R(i,j)=R(i,j)/p%nb_pos
             R(j,i)=R(i,j)
             if (R(i,j).lt.Rmin) Rmin=R(i,j)
          enddo
       enddo
       ! Rescale
       do i=1,nb
          do j=i,nb
             R(i,j)=(R(i,j)-Rmin)/(1.0-Rmin)
             R(j,i)=R(i,j)
          enddo
       enddo
       ! Deallocate
       deallocate(haplos)
       ! Store upper triangle of correlation matrix
       if(exportcorr)THEN 
          open(20,file=trim(p_p%aname)//'.corr.polyg',err=992)
          do i=1,nb
             do j=i,nb
                if (R(i,j).gt.0.0) write(20,*) i,j,R(i,j)
             enddo
          enddo
          close(20)
       endif
    else
       print *,'Error: unknown tag to compute (inverse) correlation matrix: ',tag
       Compute_Inv_Corr=1
    endif
    goto 1000
992 continue
    print *,'Error: cannot create correlation output file...'
    Compute_Inv_Corr=8
    goto 1000
995 continue
    print *,'Error: cannot open haplotypes file ',p%hfile
    Compute_Inv_Corr=5
    goto 1000
997 continue
    print *,'Error: problem while allocating in Compute_Inv_Corr'
    Compute_Inv_Corr=3
    goto 1000
998 continue
    print *,'Error: haplotypes can only be used to build direct correlation matrices'
    Compute_Inv_Corr=2
    goto 1000
    !999 continue

    !    goto 1000
1000 continue
  end function Compute_Inv_Corr
  !--------------------------------------- 
  ! Reading correlation matrices 
  !--------------------------------------- 
  subroutine Read_Corr(pc,fname)
    real(kind=dp),dimension(:,:)::pc     
    character(100),intent(in)::fname
    integer::ii,jj,io,line
    real(kind=dp)::xx
    logical::ExFile
    !
    inquire(file=fname,exist=ExFile)
    if(.not.ExFile)THEN
       write(*,'(3a)')"Error: file ",trim(adjustl(fname))," does not exist "
       stop
    else 
       write(*,'(2a)')"Reading ",trim(adjustl(fname))
    endif
    line=0
    open(20,file=fname,iostat=io,action="READ")
    do 
       read(20,*,iostat=io) ii,jj,xx
       if(io/=0)exit
       line=line+1
       if(max(ii,jj)>size(pc,1))then
          write(*,'(3a,i0,a,2(i0," "),a,2(i0,a))')"Error : ",trim(adjustl(fname)),&
               &" : line ",line," have element ",ii,jj," that are out of bounds [",&
               &size(pc,1),":",size(pc,2),"]"
          stop 
       else 
          pc(ii,jj)=xx
          pc(jj,ii)=xx
       end if
    end do
  end subroutine Read_Corr

  !---------------------------------------------------------------------!
  ! Obtain random effects (inverse) correlation data from data files
  ! Since 3 possibilities exist for RFILE and RIFILE, 9 situations exist
  ! RF = I       and RIF = I         => no correlation
  !    = Compute     RIF = I         => Compute R (later, RI=Ginv(R))
  !    = <fname>     RIF = I         => Obtain R (later, RI=Ginv(R))
  ! RF = I       and RIF = Compute   => Compute RI directly (no R needed)
  !    = I       and RIF = <fname>   => Obtain RI directly (no R needed)
  ! Other situations do not make sense.
  !---------------------------------------------------------------------!
  function Get_Corr(p,corr,i)
    !
    integer,intent(in)::i
    type(parameters)::p
    type(correl),intent(inout)::corr
    integer::Get_Corr
    !
    integer::is,rc,CptElt,Ani1,Ani2
    real(kind=dp)::InvRmin,Rmin
    !
    Get_Corr=0
    rc=0
    corr%is_direct=.false.
    corr%is_inverse=.false.
    corr%is_diagonal=.false.
    corr%is_identity=.false.
    is=p%dim_rnd(i)
    CdtRfile:if (p%rfiles(i).eq.'I') then
       if (p%rifiles(i).eq.'I') then
          corr%inv_cor(:,1)=1.0
          corr%dir_cor(:,1)=1.0
          corr%is_direct=.true.
          corr%is_inverse=.true.
          corr%is_diagonal=.true.
          corr%is_identity=.true.
          !elseif (p%rifiles(i).eq.'PEDIGREE') then
          !   !Check wether corresponding file already exists
          !   if(NoCompute)then 
          !      inquire(file=trim(p%aname)//".invcorr.polyg.txt",exist=is_file)
          !      if (is_file) then 
          !         the_file=trim(p%aname)//".invcorr.polyg.txt" 
          !         call Read_Corr(corr%inv_cor,the_file,'txt')
          !         corr%is_inverse=.true.
          !      else
          !         inquire(file=trim(p%aname)//".invcorr.polyg.bin",exist=is_file)
          !         if (is_file) then
          !            the_file=trim(p%aname)//".invcorr.polyg.bin" 
          !            call Read_Corr(corr%inv_cor,the_file,'bin')
          !            corr%is_inverse=.true.
          !         else 
          !            NoCompute=.FALSE.
          !         endif
          !      endif
          !   endif
          !   ! If external files does not exist, or we have to compute
          !   if (.not. NoCompute) then
          !      rc=Compute_Inv_Corr(p,corr%inv_cor,-1,'PEDIGREE  ')
          !      if (rc.gt.0) goto 991
          !      corr%is_inverse=.true.
          !   endif
       elseif (p%rifiles(i).eq.'COMPUTE') then
          ! Check wether corresponding file already exists
          !inquire(file=trim(p%aname)//".invcorr.polyg.txt",exist=is_file)
          !if (is_file) then
          !   the_file=trim(p%aname)//".invcorr.polyg.txt" 
          !   call Read_Corr(corr%inv_cor,the_file,'txt')
          !   corr%is_inverse=.true.
          !else
          !   inquire(file=trim(p%aname)//".invcorr.polyg.bin",exist=is_file)
          !   if (is_file) then
          !      the_file=trim(p%aname)//".invcorr.polyg.bin" 
          !      call Read_Corr(corr%inv_cor,the_file,'bin')
          !      corr%is_inverse=.true.
          !   endif
          !endif
          !! If external files do not exist, cannot compute
          !!if (.not. is_file) goto 993
       else ! Read it
          call Read_Corr(corr%inv_cor,p%rifiles(i))
          corr%is_inverse=.true.
       endif
    elseif  (p%rfiles(i).eq.'PEDIGREE') then
       ! Check wether corresponding file already exists
       !inquire(file=trim(p%aname)//".corr.polyg.txt",exist=is_file)
       !if (is_file) then
       !   the_file=trim(p%aname)//".corr.polyg.txt" 
       !   call Read_Corr(corr%dir_cor,the_file,'txt')
       !else
       !   inquire(file=trim(p%aname)//".corr.polyg.bin",exist=is_file)
       !   if (is_file) then
       !      the_file=trim(p%aname)//".corr.polyg.bin" 
       !      call Read_Corr(corr%dir_cor,the_file,'bin')
       !   endif
       !endif
       !! If external files do not exist, compute
       !if (.not. is_file) then
       print *,'=> Computing correlation matrix for effect POLYG'
       rc=Compute_Inv_Corr(p,corr%dir_cor,+1,'PEDIGREE  ')
       if (rc.gt.0) goto 991
       !endif
       corr%is_direct=.true.
       if (p%rifiles(i).ne.'I') goto 996
    elseif  (p%rfiles(i).eq.'COMPUTE') then
       ! Check wether corresponding file already exists
       !inquire(file=trim(p%aname)//".corr.polyg.txt",exist=is_file)
       !if (is_file) then
       !   the_file=trim(p%aname)//".corr.polyg.txt" 
       !   call Read_Corr(corr%dir_cor,the_file,'txt')
       !else
       !   inquire(file=trim(p%aname)//".corr.polyg.bin",exist=is_file)
       !   if (is_file) then
       !      the_file=trim(p%aname)//".corr.polyg.bin" 
       !      call Read_Corr(corr%dir_cor,the_file,'bin')
       !   endif
       !endif
       !! If external files do not exist, compute
       !if (.not. is_file) then
       write(*,*)'=> Computing correlation matrix for effect polyg'
       if(allocated(GenoFileName))then
          CptElt=size(corr%dir_cor,1) 
          allocate(R(IAd(CptElt,CptElt)));R=0.0
          !---------------------------------------------------------------------
          ! START of the OPENMP section
          !---------------------------------------------------------------------
          !$OMP PARALLEL DEFAULT(PRIVATE)shared(pardim,GenoFileName,R,InvRmin,Rmin,pheno_id)
          !$OMP DO SCHEDULE(RUNTIME)
          do CptElt=1,size(GenoFileName,1)                
             call   LecFicTypage(GenoFileName(CptElt),pardim(CptElt,1:4),R,pheno_id)
          end do
          !$OMP END DO 
          !$OMP SINGLE
          do CptElt=2,size(GenoFileName,1)
             pardim(CptElt,3)=pardim(CptElt,3)+pardim(CptElt-1,4)
             pardim(CptElt,4)=pardim(CptElt,4)+pardim(CptElt-1,4)
          end do
          if(DbgMode)then !Control boundaries---------------------
             do CptElt=1,size(GenoFileName,1)
                write(*,*)"File ",GenoFileName(CptElt)," Pos ",pardim(CptElt,3)," to "&
                     &,pardim(CptElt,4)
             end do
          endif
          !$OMP END SINGLE
          !$OMP FLUSH
          !$OMP WORKSHARE
          InvRmin=1.0_dp/maxval(pardim(:,4))
          R=InvRmin*R
          Rmin=minval(R(:))
          InvRmin=1.0_dp/(1.0_dp-Rmin)
          R=InvRmin*(R-Rmin)
          !$OMP END WORKSHARE
          !$OMP END PARALLEL
          !---------------------------------------------------------------------
          ! END of the OPENMP section
          !---------------------------------------------------------------------
          if(exportcorr)THEN 
             open(20,file=trim(p_p%aname)//'.corr.polyg')
             do Ani1=1,size(corr%dir_cor,1)
                do Ani2=Ani1,size(corr%dir_cor,1)
                   if (R(IAd(Ani1,Ani2)).gt.0.0) write(20,*)Ani1,Ani2,R(IAd(Ani1,Ani2))
                enddo
             enddo
             close(20)
          endif
          !Version longue !!!
       else 
          rc=Compute_Inv_Corr(p,corr%dir_cor,+1,'HAPLOTYPES')
       endif
       if (rc.gt.0) goto 991
       !endif
       corr%is_direct=.true.
       if (p%rifiles(i).ne.'I') goto 996
    else 
       if (p%rifiles(i).eq.'I') then
          call Read_Corr(corr%dir_cor,p%rfiles(i))
          corr%is_direct=.true.
       else
          goto 996
       endif
    endif CdtRfile

    goto 1000
991 continue
    print *,'Error: problem while computing correlation files'
    Get_Corr=2
    goto 1000
996 continue
    print *,'Error: when corr is a file name or needs to be computed, inv corr should be I'
    Get_Corr=6
    goto 1000
1000 continue
  end function Get_Corr
  !---------------------------------------------------------------------- 
  ! Getting rank of matrix A 
  !---------------------------------------------------------------------- 
  function Get_Rank(A,n,m,r)
    !
    integer,intent(in)::n,m
    real(kind=dp),intent(in)::A(n,m)
    integer,intent(out)::r
    integer::Get_Rank
    !
    real(kind=dp)::B(n,m),max_tol,pivot
    integer::indx(n),itemp,i,j,k,ipvt
    !
    Get_Rank=0
    max_tol=1.0d-10
    B=A
    r=0
    do i=1,n
       indx(i)=i
    enddo
    do j=1,m
       ! Look for pivot
       ipvt=j
       do i=j+1,n
          if (abs(B(indx(i),j)).gt.abs(B(indx(ipvt),j))) ipvt=i
       enddo
       if (abs(B(indx(ipvt),j)).lt.max_tol) cycle
       r=r+1 
       ! Pivot lines  
       itemp=indx(j)
       indx(j)=indx(ipvt)
       indx(ipvt)=itemp
       ! Substract k*pivot line
       do i=j+1,n
          pivot=B(indx(i),j)/B(indx(j),j)
          do k=j+1,m
             B(indx(i),k)=B(indx(i),k)-B(indx(j),k)*pivot
          enddo
       enddo
    enddo
  end function Get_Rank
  !--------------------------------------- 
  ! Obtain design matrices for haplotypes
  !--------------------------------------- 
  function Get_Design_Haplos(hfile,posr,ifrom,ito,D)
    character(100),intent(in)::hfile
    integer,intent(in)::posr,ifrom,ito
    integer(kind=1),intent(out)::D(:,:,:)
    integer::Get_Design_Haplos
    !
    integer::np,i,ir,io
    real(kind=dp)::dummy
    !
    Get_Design_Haplos=0
    D=0
    np=ito-ifrom+1
    open(20,file=hfile,iostat=io,action="READ")
    if(io/=0)then
       write(*,'("Error: problem while reading haplotypes file ",a)')trim(adjustl(hfile))
       Get_Design_Haplos=1
    endif
    ir=1
    do 
       read(20,*,iostat=io)(dummy,i=1,posr+2*ifrom-3),(D(ir,i,1),D(ir,i,2),i=1,np)
       if(io/=0)exit
       ir=ir+1
    end do
    close(20)
  end function Get_Design_Haplos
  !--------------------------------------- 
  ! Obtain design matrices for haplotypes
  !--------------------------------------- 
  function Get_Design_Haplos2(hfile,deb,fin,Taille,D)
    character(*),intent(in)::hfile
    integer,intent(in)::deb,fin,Taille
    integer(kind=1),intent(out)::D(:,:,:)
    integer(kind=1),allocatable::temp(:)
    integer::Get_Design_Haplos2
    !
    integer::np,i,j,ir,io
    real(kind=dp)::dummy
    !
    Get_Design_Haplos2=0
    D=0
    open(20,file=hfile,iostat=io,action="READ")
    if(io/=0)then
       write(*,'("Error: problem while reading haplotypes file ",a)')trim(adjustl(hfile))
       Get_Design_Haplos2=1
    endif
    ir=1 ; np=size(D,2) ;
    allocate(temp(Taille));temp=0
    if(DbgMode)write(*,*)deb,fin
    do 
       read(20,*,iostat=io)i,j,temp
       if(io/=0)exit
       D(ir,:,1)=temp(deb:fin)
       read(20,*,iostat=io)i,j,temp
       if(io/=0)exit
       D(ir,:,2)=temp(deb:fin)
       ir=ir+1
    end do
    close(20)
    Deallocate(temp)
  end function Get_Design_Haplos2

  !--------------------------------------- 
  ! Obtain design matrices from data file
  !---------------------------------------
  function Get_Design(dfile,nf,posf,dimf,D,n,p,base,blankp)
    !
    character(100),intent(in)::dfile
    integer,intent(in)::nf,posf,dimf,n,p,base
    logical,intent(in)::blankp
    real(kind=dp),intent(inout)::D(n,p)
    integer::Get_Design
    !
    integer::i,j,io
    real(kind=dp)::dummy,dum1
    !
    Get_Design=0
    if (blankp) D(:,nf)=0 
    !
    open(11,file=dfile,iostat=io,action="read")
    if(io/=0)then
       write(*,'("Error: problem while reading haplotypes file ",a)')trim(adjustl(dfile))
       Get_Design=1
    endif
    !
    do i=1,n
       read(11,*,iostat=io) (dummy,j=1,posf-1),dum1
       if (dimf.gt.0) then 
          D(i,nf)=int(dum1)
          if(nf>first_fact .and. base>0)D(i,nf)=D(i,nf)-1 ! remove arbitrarily one level, base =0 for random effects
       else
          D(i,nf)=dum1
       endif
    enddo
    if(i<n)then
       write(*,*)"Error: premature end of file in design matrix creation "
       stop
    endif
    close(11)  
  end function Get_Design
  !--------------------------------------- 
  ! Obtain data from data file
  !--------------------------------------- 
  function Get_Data(dfile,pos_y,y,pos_id,id,n)
    !
    character(100),intent(in)::dfile
    integer,intent(in)::pos_y,pos_id,n
    real(kind=dp),intent(out)::y(n)
    integer,intent(out)::id(n)
    integer::Get_Data
    !
    integer::nr,i,io
    real(kind=dp)::dummy(max(pos_id,pos_y))
    !
    Get_Data=0
    open(unit=11,file=dfile,action="READ")
    nr=0
    do 
       read(11,*,iostat=io) (dummy(i),i=1,max(pos_id,pos_y)) 
       if(io/=0)exit
       nr=nr+1
       y(nr)=dummy(pos_y)
       id(nr)=int(dummy(pos_id))
    end do
    close(unit=11)
    if (nr.ne.n)then 
       print '("Error: number of records and number of reads dont match (",i4,":",i4)',n,nr
       Get_Data=1
    endif

  end function Get_Data
  !--------------------------------------- 
  ! Obtain parameters from parameters file
  !--------------------------------------- 
  function Get_Parameters(p,pfile)
    integer Get_Parameters
    type(parameters)::p
    ! Work variables
    logical::is_file
    character(len=5)::key
    integer::i,j,nft,nrt,nb_levels,oldid,io
    real(kind=dp),allocatable::dum1(:),dum2(:)
    character(len=30000)::line,lines2
    character(len=100)::pfile
    character(len=30)::rname
    character(len=30)::method
    real(kind=dp)::start
    !Multigenofiles
    integer::NbFicTyp
    ! Start code
    Get_Parameters=0
    ! Init parameters
    p%aname=''   
    p%dfile=''   
    p%hfile=''
    p%datid=1    
    p%hapid=1
    p%nb_fix=0 ! already computed, but reused here
    p%nb_rnd=0 ! already computed, but reused here
    nft=0        
    nrt=0
    p%linkf='IDENT '
    p%distr=''
    p%pos_trt=0
    p%rfiles(1)='COMPUTE'
    p%rfiles(2)='I'
    p%rifiles(1:2)='I'
    p%start_tau2=1.0
    p%start_sigma2=1.0 
    p%start_b='none'
    p%start_u='none'
    p%reml_method='AI'
    p%reml_maxiter=1000
    p%reml_maxtol=0.0000000001
    p%scores_pos(1:2)=0
    p%score_method='TZENG'
    p%nb_scores=0
    p%nb_perms=1000
    p%slice=10000
    p%showp=1000
    !======================================================================
    ! Part 1 Reading parameter file
    !======================================================================
    ! Read all parameters
    open(unit=10,file=pfile,err=999)
    !----------------------------------------------------------------------------
    ParReadLoop:do 
       read(10,'(a30000)',iostat=io)line
       if(io/=0)exit
       lines2=trim(adjustl(line))
       line=trim(lines2)
       !skip empty or commentary lines
       if((lines2(1:1)=="#").or.(len_trim(lines2)==0))cycle
       read(lines2,'(a,a30000)')key,lines2
       !Check for key length correctness
       if(len_trim(adjustl(key)).ne.5)THEN
          write(*,*)"mispelled  key ",key," Keys should be 5 characters long"
          write(*,*)"glascow will now stop"
          stop
       End if
       select case(adjustl(key))
       case('PFILE')
          ! Format is: 'PFILE datafile'
          read(lines2,*,iostat=io) p%dfile
          if(io/=0)then
             write(*,'(5a)')"Error in file : '",trim(adjustl(pfile))," ' line : ' ",&
                  &trim(adjustl(line))," ' is incorrectly written "
             stop
          endif
          !Test if the file exist
          inquire(file=p%dfile,exist=is_file)
          if(.not.is_file)then 
             write(*,*)"Error: file ",trim(adjustl(p%dfile))," does not exist "
             stop            
          endif
          if (p%hfile.eq.'') p%hfile=p%dfile  ! Default
       case('GFILE') 
          ! Format is: 'GFILE haplos_file'
          ! haplos_file can be the same as datafile         
          call Nb_Element_Chaine(lines2,NbFicTyp)
          NbFicTyp=NbFicTyp
          if(NbFicTyp>1)then
             allocate(GenoFileName(NbFicTyp),pardim(NbFicTyp,4))
             read(lines2,*,iostat=io)GenoFileName
             pardim(NbFicTyp,4)=0
          else  
             read(lines2,*,iostat=io) p%hfile
             inquire(file=p%hfile,exist=is_file)
             if(.not.is_file)then 
                write(*,*)"Error: file ",trim(adjustl(p%hfile))," does not exist "
                stop            
             endif
          endif
          if(io/=0)then
             write(*,'(5a)')"Error in file : '",trim(adjustl(pfile))," ' line : ' ",&
                  &trim(adjustl(line))," ' is incorrectly written "
             stop
          endif
       case('ANAME')
          ! Format is: 'ANAME analysis_name'
          read(lines2,*,iostat=io) p%aname
          if(io/=0)then
             write(*,'(5a)')"Error in file : '",trim(adjustl(pfile))," ' line : ' ",&
                  &trim(adjustl(line))," ' is incorrectly written "
             stop
          endif
          !case('GENID')
          !   ! Format is: 'GENID id'
          !   ! Position of id field in haplotypes file
          !   ! In the PROCESS phase, the value should be 1
          !   read(lines2,*,iostat=io) p%hapid
          !   if(io/=0)then
          !      write(*,'(5a)')"Error in file : '",trim(adjustl(pfile))," ' line : ' ",&
          !           &trim(adjustl(line))," ' is incorrectly written "
          !      stop
          !   endif
          !case('PHEID')
          !   ! Format is: 'PHEID id'
          !   ! Position of id field in data file
          !   ! In the PROCESS phase, the value should be 1
          !   read(lines2,*,iostat=io) p%datid
          !   if(io/=0)then
          !      write(*,'(5a)')"Error in file : '",trim(adjustl(pfile))," ' line : ' ",&
          !           &trim(adjustl(line))," ' is incorrectly written "
          !      stop
          !   endif
       case('FIXED')
          p%nb_fix=p%nb_fix+1
          ! Format is: 'FIXED name nb_levels position'
          ! nb_levels is 0 for covariates  
          read(lines2,*,iostat=io) p%name_fix(p%nb_fix),p%dim_fix(p%nb_fix),p%pos_fix(p%nb_fix)
          if(io/=0)then
             write(*,'(5a)')"Error in file : '",trim(adjustl(pfile))," ' line : ' ",&
                  &trim(adjustl(line))," ' is incorrectly written "
             stop
          endif
          nft=nft+p%dim_fix(p%nb_fix)
       case('POLYG')
          p%nb_rnd=p%nb_rnd+1
          ! Format is: 'POLYG nb_indiv pos'
          ! nb_indivs is the number of individuals
          ! pos is the position of the animal id
          read(lines2,*,iostat=io) p%dim_rnd(1),p%pos_rnd(1)
          if(io/=0)then
             write(*,'(5a)')"Error in file : '",trim(adjustl(pfile))," ' line : ' ",&
                  &trim(adjustl(line))," ' is incorrectly written "
             stop
          endif
          p%name_rnd(1)='polyg'
          nrt=nrt+p%dim_rnd(1)
       case('HAPLO')
          p%nb_rnd=p%nb_rnd+1
          ! Format is: 'HAPLO nb_haplos nb_pos pos1'
          ! nb_haplos is the (max) number of haplotypes in each position
          ! nb_pos is the number of positions along the tested genomic region  
          ! pos1 is the position of the first haplotype for the first position
          !  (successive haplotypes are supposed to be stored on successive positions)
          read(lines2,*,iostat=io) p%dim_rnd(2),p%nb_pos,p%pos_rnd(2)
          if(io/=0)then
             write(*,'(5a)')"Error in file : '",trim(adjustl(pfile))," ' line : ' ",&
                  &trim(adjustl(line))," ' is incorrectly written "
             stop
          endif
          p%name_rnd(2)='haplo'
          nrt=nrt+p%dim_rnd(2)  ! Note that this is a maximum, some positions might show less haplotypes
          !set pscore_pos to default values if not already set
          if(p%scores_pos(1)==0)then
             p%scores_pos(1)=1;p%scores_pos(2)=p%nb_pos
             p%nb_scores=p%scores_pos(2)-p%scores_pos(1)+1
          endif
          !case('DISTR')
          !   ! Format is: 'DISTR NORMAL|BINOM |POISS '  
          !   read(lines2,*,iostat=io)tmpdist
          !   if(io/=0)then
          !      write(*,'(5a)')"Error in file : '",trim(adjustl(pfile))," ' line : ' ",&
          !           &trim(adjustl(line))," ' is incorrectly written "
          !      stop
          !   endif
          !   if((p%distr.ne.'').and.(p%distr.ne.tmpdist))THEN
          !      write(*,'(4a)')"Error :",trim(adjustl(tmpdist)),"defined in line : ",trim(adjustl(line))
          !      write(*,'(3a)')"While ",p%distr," was expected based on phenotypes info"
          !      stop
          !   endif
          !case('LINKF')
          !   ! Format is: 'LINKF link_function'  
          !   read(lines2,*,iostat=io)p%linkf
          !   if(io/=0)then
          !      write(*,'(5a)')"Error in file : '",trim(adjustl(pfile))," ' line : ' ",&
          !           &trim(adjustl(line))," ' is incorrectly written "
          !      stop
          !   endif
       case('PHENO')
          ! Format is: 'PHENO name nb_levels position'  
          read(lines2,*,iostat=io) p%name_trt,nb_levels,p%pos_trt
          if(io/=0)then
             write(*,'(5a)')"Error in file : '",trim(adjustl(pfile))," ' line : ' ",&
                  &trim(adjustl(line))," ' is incorrectly written "
             stop
          endif
          select case(nb_levels)
          case(1)
             if (p%distr.eq.''.or.p%distr.eq.'NORMAL') then
                p%distr='NORMAL'
             else 
                write(*,*)"Normal phenotypes as indicate in line :",trim(lines2)
                write(*,*)"Are in disagrement with the previously define distribution ",p%distr
                stop
             endif
          case(2)
             if (p%distr.eq.''.or.p%distr.eq.'BINOM') then
                p%distr='BINOM '
             else 
                write(*,'(2a)')"Binary phenotypes as indicate in line :",trim(lines2)
                write(*,'(2a)')"Are in disagrement with the previously define distribution ",p%distr
                STOP 
             endif
          case default
             if(nb_levels>2)then
                if (p%distr.eq.''.or.p%distr.eq.'POISS ') then
                   p%distr='POISS '
                else 
                   write(*,'(2a)')"Count phenotypes as indicate in line :",trim(lines2)
                   write(*,'(2a)')"Are in disagrement with the previously define distribution",p%distr
                   STOP 
                endif
             else 
                write(*,'(2a)')"Phenotypes type is incorrect in line :",trim(lines2)
                stop
             endif
          end select
       case('RFILE')
          ! Format is: 'RFILE random_name -1/+1 file_name|method '
          ! The 1/-1 indicates whether we have a direct (1) or inverse(-1) corr matrix
          ! <file_name> should contain the non-null coefficients of the (inv) corr matrix
          !             For HAPLO, the name will be expanded with the position
          !                        f.e. corr.hap => corr.hap.<pos>
          ! <method> can either be 'PEDIGREE' (for <random_name> = 'POLYG') 
          !                     or 'COMPUTE' (for <random_name> = 'POLYG') 
          read(lines2,*,iostat=io)i,method
          if(io/=0)then
             write(*,'(5a)')"Error in file : '",trim(adjustl(pfile))," ' line : ' ",&
                  &trim(adjustl(line))," ' is incorrectly written "
             stop
          endif
          !Check if specified correlation file exist
          if(adjustl(trim(method)).ne.'COMPUTE')THEN
             inquire(file=adjustl(trim(method)),exist=ExFile)
             if(.not.ExFile)THEN
                write(*,*)"File ",adjustl(trim(method))," doesn't exist : Process will stop"
                stop
             endif
          endif
          j=1
          !set r(i)files according to i
          if (i.eq.1) then
             p%rfiles(j)=method
             p%rifiles(j)='I'
          else
             p%rfiles(j)='I'
             p%rifiles(j)=method
          endif
       case('REMLT')
          ! Format is: 'REMLT max_reml_tolerance'  
          read(lines2,*,iostat=io) p%reml_maxtol
          if(io/=0)then
             write(*,'(5a)')"Error in file : '",trim(adjustl(pfile))," ' line : ' ",&
                  &trim(adjustl(line))," ' is incorrectly written "
             stop
          endif
       case('REMLI')
          ! Format is: 'REMLI max_reml_iterations'  
          read(lines2,*,iostat=io) p%reml_maxiter
          if(io/=0)then
             write(*,'(5a)')"Error in file : '",trim(adjustl(pfile))," ' line : ' ",&
                  &trim(adjustl(line))," ' is incorrectly written "
             stop
          endif
       case('REMLS')
          ! Format is: 'REMLS random_name start_sigma_n'  
          read(lines2,*,iostat=io) rname,start
          if(io/=0)then
             write(*,'(5a)')"Error in file : '",trim(adjustl(pfile))," ' line : ' ",&
                  &trim(adjustl(line))," ' is incorrectly written "
             stop
          endif
          select case(trim(adjustl(rname)))
          case('RESIDUAL')
             p%start_tau2=start
          case('POLYG')
             p%start_sigma2(1)=start 
          case('HAPLO')
             p%start_sigma2(2)=start 
          case default 
             write(*,'(3a)')"Error : unknown Variance TYPE",trim(adjustl(rname))," in line ",&
                  &trim(adjustl(line))
             stop 
          end select
       case('SCORE')
          ! Format is: 'SCORE pos_1 pos_n method'
          read(lines2,*,iostat=io) p%scores_pos(1),p%scores_pos(2)!,p%score_method
          if(io/=0)then
             write(*,'(5a)')"Error in file : '",trim(adjustl(pfile))," ' line : ' ",&
                  &trim(adjustl(line))," ' is incorrectly written "
             stop
          endif
          p%nb_scores=p%scores_pos(2)-p%scores_pos(1)+1
       case('SLICE')
          ! Format is: 'SLICE slice_size'
          read(lines2,*,iostat=io) p%slice
          if(io/=0)then
             write(*,'(5a)')"Error in file : '",trim(adjustl(pfile))," ' line : ' ",&
                  &trim(adjustl(line))," ' is incorrectly written "
             stop
          endif
       case('SHOWP')
          ! Format is: 'SHOWP nb_pos'
          read(lines2,*,iostat=io) p%showp
          if(io/=0)then
             write(*,'(5a)')"Error in file : '",trim(adjustl(pfile))," ' line : ' ",&
                  &trim(adjustl(line))," ' is incorrectly written "
             stop
          endif
       case('PERMS')
          ! Format is: 'PERMS nb_perms'
          read(lines2,*,iostat=io) p%nb_perms
          if(io/=0)then
             write(*,'(5a)')"Error in file : '",trim(adjustl(pfile))," ' line : ' ",&
                  &trim(adjustl(line))," ' is incorrectly written "
             stop
          endif
       case default  
          write(*,*)"unknow key ",key," in parameter file"
          write(*,*)"glascow will now stop"
          stop
       end select
    end do ParReadLoop
    close(unit=10)
    !======================================================================
    ! Part 2 checking parameter correctness
    !======================================================================
    !Check Distribution and link function----------------------------------------------
    select case(trim(adjustl(p%distr)))
    case("NORMAL") 
       p%linkf='IDENT'
    case('BINOM ')
       p%linkf='LOGIT'
    case('POISS ')
       p%linkf='LOG'
    end select
    !---------------------------------------------------------------------------------
    if (p%aname.eq.'') goto 973
    if (p%pos_trt.eq.0) goto 977
    if (p%name_rnd(2).eq.'') goto 979
    ! Check whether pheno and geno records match
    allocate(dum1(p%datid),dum2(p%hapid),stat=ok)
    if (ok.ne.0) goto 983
    if (p%hfile.ne.p%dfile) then
       !the file have previsously inquired
       open(11,file=p%dfile,iostat=io)
       if(io/=0)then 
          write(*,*)'Error: cannot open data file ',p%dfile
          Get_Parameters=3
          return 
       end if
       open(12,file=p%hfile,err=994)
       oldid=0
55     continue
       read(11,*,end=56) (dum1(i),i=1,p%datid)
       print*,'Datid ::',p%datid,(dum1(p%datid))
       read(12,*,end=993) (dum2(i),i=1,p%hapid)
       print*,'Hapid ::',p%hapid,(dum2(p%hapid))
57     continue
       if (int(dum1(p%datid)).lt.int(dum2(p%hapid))) then
          if (int(dum1(p%datid)).eq.oldid) then
             write(*,*) 'Warning: several phenotyped individuals (twins ?) share same id ',int(dum1(p%datid))
             read(11,*,end=56) (dum1(i),i=1,p%datid)
             goto 57
          else
             goto 992
          endif
       elseif (int(dum1(p%datid)).gt.int(dum2(p%hapid))) then
          goto 974
       end if
       oldid=int(dum1(p%datid))
       goto 55
56     continue
       read(12,*,end=60)
       goto 985
60     continue
       close(12)
       close(11)
    else
       if (p%datid.ne.p%hapid) goto 984
    endif
    deallocate(dum2,dum1)
    ! Obtain number of records from data file
    ! Read from data file
    ! Count number of records
    p%nb_rec=0
    open(unit=11,file=p%dfile,action="READ",iostat=io)
    if(io/=0)then 
       write(*,*)'Error: cannot open data file ',p%dfile
       Get_Parameters=3
       return 
    end if
    do 
       read(11,*,iostat=io)line
       if(io/=0)exit
       p%nb_rec=p%nb_rec+1
    end do

    close(unit=11)
    ! Report parameters
    print '("Nb fixed effects:  ",i2)',p%nb_fix
    print '("Nb random effects: ",i2)',p%nb_rnd
    print '("Nb records:        ",i5)',p%nb_rec
    goto 1000
973 continue
    print *,'Error: analysis name (ANAME) should be specified'
    Get_Parameters=26
    goto 1000
974 continue
    print *,'Error: genotyped id should have a phenotypic record',int(dum2(p%hapid))
    Get_Parameters=25
    goto 1000
977 continue
    print *,'Error: trait should be specified...'
    Get_Parameters=23
    goto 1000
979 continue
    print *,'Error: haplotype effect should be included in PRM file'
    Get_Parameters=21
    goto 1000
983 continue
    print *,'Error: problem while allocating dummy structures'
    Get_Parameters=17
    goto 1000
984 continue
    print *,'Error: HAPID and DATID in the same file should be the same'
    Get_Parameters=16
    goto 1000
985 continue
    print *,'Error: data file ',p%dfile,' end prematurely'
    Get_Parameters=15
    goto 1000
992 continue
    print *,'Error: phenotyped individuals should have an haplotype: ',int(dum1(p%datid))
    Get_Parameters=7
    goto 1000
993 continue
    print *,'Error: haplotypes file ',p%hfile,' end prematurely'
    Get_Parameters=7
    goto 1000
994 continue
    print *,'Error: cannot open haplotypes file ',p%hfile
    Get_Parameters=6
    goto 1000
999 continue
    print *,'Error: cannot open parameters file'
    Get_Parameters=1
    goto 1000
    ! End of the subroutine
1000 continue
  end function Get_Parameters

  !--------------------------------------- 
  ! Obtain dimensions from parameters file
  ! Input:   p   parameters
  ! Output:  nb_fix     set
  !          nb_rnd(=2) set and checked
  !--------------------------------------- 
  function Get_Dimensions(p,pfile)
    integer::Get_Dimensions
    type(parameters)::p
    ! Work variables
    logical::is_polyg,is_haplo
    character(len=300)::line,lines2
    character(len=5)::key
    integer::i,j,k,io
    character(len=30)::rname
    character(len=100)::pfile
    ! Start code
    Get_Dimensions=0
    ! Read from parameters file
    open(unit=10,file=pfile,iostat=io)
    if(io/=0)Then
       write(*,*)'Error: cannot open parameters file',pfile
       Get_Dimensions=1
    else
       ! Count fixed and random effects in parameters file
       p%nb_fix=0       ; p%nb_rnd=0
       is_polyg=.false. ; is_haplo=.false.          
       LoopDim:do 
          read(10,'(a300)',iostat=io)line
          if(io/=0)exit
          lines2=adjustl(line)
          line=trim(lines2)
          !skip empty or commentary lines
          if((lines2(1:1)=="#").or.(len_trim(lines2)==0))cycle LoopDim
          !Split line in key + remaining line
          read(line,'(a,a300)')key,lines2
          if(len_trim(adjustl(key)).ne.5)THEN
             write(*,*)"unknow key ",key
             write(*,*)"glascow will now stop"
             stop
          end if
          select case(adjustl(key))
          case('FIXED')!--------------------------------------
             p%nb_fix=p%nb_fix+1
             read(lines2,*,iostat=io) rname,i,j
             if(io/=0)then
                write(*,'(5a)')"Error in file : '",trim(adjustl(pfile))," ' line : ' ",&
                     &trim(adjustl(line))," ' is incorrectly written "
                stop
             endif
          case ('POLYG')!--------------------------------------
             is_polyg=.true.
             p%nb_rnd=p%nb_rnd+1
             read(lines2,*,iostat=io) i,j
             if(io/=0)then
                write(*,'(5a)')"Error in file : '",trim(adjustl(pfile))," ' line : ' ",&
                     &trim(adjustl(line))," ' is incorrectly written "
                stop
             endif
          case('HAPLO')!--------------------------------------
             is_haplo=.true.
             p%nb_rnd=p%nb_rnd+1
             read(lines2,*,iostat=io) i,j,k
             if(io/=0)then
                write(*,'(5a)')"Error in file : '",trim(adjustl(pfile))," ' line : ' ",&
                     &trim(adjustl(line))," ' is incorrectly written "
                stop
             endif
          case default!--------------------------------------
             read(lines2,*)
          end select
       end do LoopDim
       !---------------------------------------------------------------------------
       close(10)
       if (.not. (is_polyg.and.is_haplo))then
          print *,'Error: HAPLO and POLYG effects should be specified'
          Get_Dimensions=4
       endif
       if (p%nb_rnd.ne.2)then
          print *,'Error: only 2 random effects should be specified'
          Get_Dimensions=5
       end if
    endif
    ! End of the subroutine
  end function Get_Dimensions

  !======================================================================
  ! Gamma density function
  !======================================================================
  FUNCTION gammad(x, p) RESULT(fn_val)

    !      ALGORITHM AS239  APPL. STATIST. (1988) VOL. 37, NO. 3

    !      Computation of the Incomplete Gamma Integral

    !      Auxiliary functions required: ALOGAM = logarithm of the gamma
    !      function, and ALNORM = algorithm AS66

    ! ELF90-compatible version by Alan Miller
    ! Latest revision - 27 October 2000

    ! N.B. Argument IFAULT has been removed

    IMPLICIT NONE
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN) :: x, p
    REAL (dp)             :: fn_val

    ! Local variables
    REAL (dp)             :: pn1, pn2, pn3, pn4, pn5, pn6, arg, c, rn, a, b, an
    REAL (dp), PARAMETER  :: zero = 0.d0, one = 1.d0, two = 2.d0, &
         oflo = 1.d+37, three = 3.d0, nine = 9.d0, &
         tol = 1.d-20, xbig = 1.d+8, plimit = 1000.d0, &
         elimit = -88.d0
    INTEGER               :: ifault
    ! EXTERNAL alogam, alnorm

    fn_val = zero

    !      Check that we have valid values for X and P

    IF (p <= zero .OR. x < zero) THEN
       WRITE(*, *) 'AS239: Either p <= 0 or x < 0'
       RETURN
    END IF
    IF (x == zero) RETURN

    !      Use a normal approximation if P > PLIMIT

    IF (p > plimit) THEN
       pn1 = three * SQRT(p) * ((x / p) ** (one / three) + one /(nine * p) - one)
       fn_val = alnorm(pn1, .false.)
       RETURN
    END IF

    !      If X is extremely large compared to P then set fn_val = 1

    IF (x > xbig) THEN
       fn_val = one
       RETURN
    END IF

    IF (x <= one .OR. x < p) THEN

       !      Use Pearson's series expansion.
       !      (Note that P is not large enough to force overflow in ALOGAM).
       !      No need to test IFAULT on exit since P > 0.

       arg = p * LOG(x) - x - alogam(p + one, ifault)
       c = one
       fn_val = one
       a = p
40     a = a + one
       c = c * x / a
       fn_val = fn_val + c
       IF (c > tol) GO TO 40
       arg = arg + LOG(fn_val)
       fn_val = zero
       IF (arg >= elimit) fn_val = EXP(arg)

    ELSE

       !      Use a continued fraction expansion

       arg = p * LOG(x) - x - alogam(p, ifault)
       a = one - p
       b = a + x + one
       c = zero
       pn1 = one
       pn2 = x
       pn3 = x + one
       pn4 = x * b
       fn_val = pn3 / pn4
60     a = a + one
       b = b + two
       c = c + one
       an = a * c
       pn5 = b * pn3 - an * pn1
       pn6 = b * pn4 - an * pn2
       IF (ABS(pn6) > zero) THEN
          rn = pn5 / pn6
          IF (ABS(fn_val - rn) <= MIN(tol, tol * rn)) GO TO 80
          fn_val = rn
       END IF

       pn1 = pn3
       pn2 = pn4
       pn3 = pn5
       pn4 = pn6
       IF (ABS(pn5) >= oflo) THEN

          !      Re-scale terms in continued fraction if terms are largeg

          pn1 = pn1 / oflo
          pn2 = pn2 / oflo
          pn3 = pn3 / oflo
          pn4 = pn4 / oflo
       END IF
       GO TO 60
80     arg = arg + LOG(fn_val)
       fn_val = one
       IF (arg >= elimit) fn_val = one - EXP(arg)
    END IF
    RETURN
  END FUNCTION gammad

  function alogam ( x, ifault )

    !*****************************************************************************80
    !
    !! ALOGAM computes the logarithm of the Gamma function.
    !
    !  Modified:
    !
    !    28 March 1999
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Malcolm Pike, David Hill.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Malcolm Pike, David Hill,
    !    Algorithm 291:
    !    Logarithm of Gamma Function,
    !    Communications of the ACM,
    !    Volume 9, Number 9, September 1966, page 684.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the argument of the Gamma function.
    !    X should be greater than 0.
    !
    !    Output, integer ( kind = 4 ) IFAULT, error flag.
    !    0, no error.
    !    1, X <= 0.
    !
    !    Output, real ( kind = 8 ) ALOGAM, the logarithm of the Gamma
    !    function of X.
    !
    implicit none
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
    real    ( kind = dp ) alogam
    real    ( kind = dp ) f
    integer :: ifault
    real    ( kind = dp ) x
    real    ( kind = dp ) y
    real    ( kind = dp ) z

    if ( x <= 0.0D+00 ) then
       ifault = 1
       alogam = 0.0D+00
       return
    end if

    ifault = 0
    y = x

    if ( x < 7.0D+00 ) then

       f = 1.0D+00
       z = y

       do while ( z < 7.0D+00 )
          f = f * z
          z = z + 1.0D+00
       end do

       y = z
       f = - log ( f )

    else

       f = 0.0D+00

    end if

    z = 1.0D+00 / y / y

    alogam = f + ( y - 0.5D+00 ) * log ( y ) - y &
         + 0.918938533204673D+00 + &
         ((( &
         - 0.000595238095238D+00   * z &
         + 0.000793650793651D+00 ) * z &
         - 0.002777777777778D+00 ) * z &
         + 0.083333333333333D+00 ) / y

    return
  end function alogam

  !  Algorithm AS66 Applied Statistics (1973) vol.22, no.3

  !  Evaluates the tail area of the standardised normal curve
  !  from x to infinity if upper is .true. or
  !  from minus infinity to x if upper is .false.

  ! ELF90-compatible version by Alan Miller
  ! Latest revision - 29 November 2001

  function alnorm ( x, upper )

    !*****************************************************************************80
    !
    !! ALNORM computes the cumulative density of the standard normal distribution.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    17 January 2008
    !
    !  Author:
    !
    !    Original FORTRAN77 version by David Hill
    !    MATLAB version by John Burkardt
    !
    !  Reference:
    !
    !    David Hill,
    !    Algorithm AS 66:
    !    The Normal Integral,
    !    Applied Statistics,
    !    Volume 22, Number 3, 1973, pages 424-427.
    !
    !  Parameters:
    !
    !    Input, real X, is one endpoint of the semi-infinite interval
    !    over which the integration takes place.
    !
    !    Input, logical UPPER, determines whether the upper or lower
    !    interval is to be integrated:
    !    1 => integrate from X to + Infinity;
    !    0 => integrate from - Infinity to X.
    !
    !    Output, real VALUE, the integral of the standard normal
    !    distribution over the desired interval.
    !
    implicit none
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
    real(kind=dp), parameter ::  a1 = 5.75885480458
    real(kind=dp), parameter ::  a2 = 2.62433121679
    real(kind=dp), parameter ::  a3 = 5.92885724438
    real(kind=dp), parameter ::  b1 = -29.8213557807
    real(kind=dp), parameter ::  b2 = 48.6959930692
    real(kind=dp), parameter ::  c1 = -0.000000038052
    real(kind=dp), parameter ::  c2 = 0.000398064794
    real(kind=dp), parameter ::  c3 = -0.151679116635
    real(kind=dp), parameter ::  c4 = 4.8385912808
    real(kind=dp), parameter ::  c5 = 0.742380924027
    real(kind=dp), parameter ::  c6 = 3.99019417011
    real(kind=dp), parameter ::  con = 1.28
    real(kind=dp), parameter ::  d1 = 1.00000615302
    real(kind=dp), parameter ::  d2 = 1.98615381364
    real(kind=dp), parameter ::  d3 = 5.29330324926
    real(kind=dp), parameter ::  d4 = -15.1508972451
    real(kind=dp), parameter ::  d5 = 30.789933034
    real(kind=dp), parameter ::  ltone = 7.0
    real(kind=dp), parameter ::  p = 0.39894228044
    real(kind=dp), parameter ::  q = 0.39990348504
    real(kind=dp), parameter ::  r = 0.398942280385
    real(kind=dp), parameter ::  utzero = 18.66
    real(kind=dp) ::x,y,z,alnorm
    logical ::up,upper

    up = upper
    z = x

    if (z .lt. 0.0D+00 )then
       up = .not. up
       z = - z
    end if

    if ( z .gt. ltone .and. ( ( .not. up ) .or. utzero .lt. z ) ) then
       if ( up ) then
          alnorm = 0.0D+00
       else
          alnorm = 1.0D+00
       end if

       return

    end if

    y = 0.5D+00 * z * z

    if ( z .le. con ) then

       alnorm = 0.5D+00 - z * (p-q*y/(y+a1+b1/(y+a2+b2/( y + a3 ))))

    else

       alnorm = r * exp ( - y )/(z+c1+d1/(z+c2+d2/(z+c3+d3/(z+c4+d4/(z+c5+d5/(z+c6))))))

    end if

    if ( .not. up ) then
       alnorm = 1.0D+00 - alnorm
    end if

    return
  end function alnorm

  !---------------------------------------------------------------------!
  ! Display a banner of the program                                     !
  !---------------------------------------------------------------------!
  subroutine banner(version)
    character(3),intent(in)::version
#define  date __DATE__
#define time __TIME__
    write(*,fmt="('*--------------------------*')")
    write(*,fmt="('*    GLASCOW version ',a3,'   *')") version
    write(*,fmt="('*     February 2012        *')")
    write(*,fmt="('*--------------------------*')")
    write(*,fmt="('*-Compilation-date-:--',a,' ',a,'   *')")date,time
  end subroutine banner

  !---------------------------------------------------------------------!
  ! Help user to built a parameter file
  !---------------------------------------------------------------------!

  subroutine InitParam(file)
    !======================================================================
    !This subroutine will questions the user in order to create 
    ! A correct minimal paramater file
    !======================================================================
    ! file : Name of the parameter file to be output
    !======================================================================
    implicit NONE

    character(len=*),intent(in)::file
    character(len=5000)::StringAnswers
    integer::NumericAnswers,Nfix

    open(10,file=file)
    write(*,*)"Project name ?"
    read(*,*)StringAnswers
    write(10,'(2a)')"ANAME ",adjustl(trim(StringAnswers))
    write(*,*)"Name of the phenotype file ?"
    read(*,*)StringAnswers
    write(10,'(2a)')"PFILE ",adjustl(trim(StringAnswers))
    write(*,*)"Name of the phenotype ?"
    read(*,*)StringAnswers
    WRITE(10,'(2(a,1x))',ADVANCE="NO")"PHENO",adjustl(trim(StringAnswers))
    write(*,*)"Type ? 1 : binary  2 :Multinomial 3 : Normal" 
    read(*,*)NumericAnswers
    if(NumericAnswers==2)THEN
       write(*,*)"Number of levels ?"
       read(*,*)NumericAnswers
    elseif(NumericAnswers==1)then 
       NumericAnswers=2
    else 
       NumericAnswers=1
    endif
    WRITE(10,'(i0,1x)',ADVANCE="NO")NumericAnswers 
    write(*,*)"Position (field number) of ",adjustl(trim(StringAnswers))," in the file ?"
    read(*,*)NumericAnswers
    WRITE(10,'(i0,1x)',ADVANCE="YES")NumericAnswers 
    write(*,*)" Number of fixed effect"
    read(*,*)nFix
    do i=1,nFix
       write(*,*)"Name of the effect ?"
       read(*,*)StringAnswers
       WRITE(10,'(2(a,1x))',ADVANCE="NO")"FIXED",adjustl(trim(StringAnswers))
       write(*,*)"Number of level for the effect (1 for a covariate) ?"
       read(*,*)NumericAnswers
       WRITE(10,'(2(i0,1x))',ADVANCE="NO")NumericAnswers
       write(*,*)"Position (field number) of effect ",adjustl(trim(StringAnswers))," in the file ?"
       read(*,*)NumericAnswers
       WRITE(10,'(2(i0,1x))',ADVANCE="YES")NumericAnswers
    end do
    write(*,*)"Name(s) of the genotype file(s) ?"
    read(*,*)StringAnswers
    WRITE(10,'(2(a,1x))',ADVANCE="YES")"GFILE",adjustl(trim(StringAnswers))
    WRITE(10,'(a,1x)',ADVANCE="NO")"POLYG"
    write(*,*)"Number of individuals ?"
    read(*,*)NumericAnswers
    WRITE(10,'(2(i0,1x))',ADVANCE="YES")NumericAnswers,1
    WRITE(10,'(a,1x)',ADVANCE="NO")"HAPLO"
    write(*,*)"Highest level of the allele/hidden states ?" 
    read(*,*)NumericAnswers
    WRITE(10,'(2(i0,1x))',ADVANCE="NO")NumericAnswers
    write(*,*)"Number of markers/hidden states ?"
    read(*,*)NumericAnswers
    WRITE(10,'(2(i0,1x))',ADVANCE="NO")NumericAnswers,2
    write(*,*)"The parameter's file has been successfully created"
    write(*,*)"Normal analyse will now start"
    write(*,*)"You can tweak it by adding some extra keywords"
    write(*,*)"Have a look into the manual ! "
    close(10)
  end subroutine InitParam

  !
  subroutine LecFicTypage(FicGenot,pardim,R,Id)
    !=============================================================================!
    ! Read a genotype file
    !=============================================================================!
    !sous soutine permettant de lire un fichier de typages, de déterminer son type
    !de vérifier sa cohérence
    !des données de typage
    !-----------------------------------------------------------------------------!
    !FicGenot   : Un fichier de typage (phase ou typage)
    !pardim     : Un vecteur de dimensions des typages
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
    !$USE OMP_LIB
    implicit none
    character(len=*),intent(in)       ::FicGenot   
    integer,dimension(4),intent(inout)::pardim
    real(kind=dp),intent(inout)       ::R(:)
    integer,dimension(:),intent(in),optional::Id(:) 

    !Variables internes
    Logical                 ::Ex_Genot
    character(len=1000000)  ::Ligne
    integer                 ::nligne,LigRef,MaxIrcGen,io,Ilec
    integer                 ::IrcCourant,NbElt,NbEltPrec,NbAnim,NvxCtrl,OMP_GET_THREAD_NUM 
    integer(kind=1),allocatable::pheno(:,:,:),Allcar(:)
    integer(kind=1)::Numal
    integer::ircTyp,NumPhase,NligneTyp,PhaseForce,CptElt,NumMarq
    integer::SumEqui,Ani1,Ani2,NbAnimGardes
    !vérifications d'usage
    inquire(file=FicGenot,EXIST=Ex_Genot)
    if(.not.(Ex_Genot))then
       write(*,*)"No genotype file called ",trim(FicGenot)
       write(*,*)"Program will stop"
       stop
    else
       if(DbgMode)write(*,*)"File named ",trim(FicGenot)," exist"
    endif
    Ilec=0 ;
    !ajouter une exécution cdtionnelle en openMP
    !$ Ilec=OMP_GET_THREAD_NUM()
    !ouverture pour première lecture
    Ilec=30+Ilec ;nligne=0 ; MaxIrcGen=0;pardim=1;NbAnimGardes=0!MAZ
    open(Ilec,file=FicGenot,action="read")
    !$ if(DbgMode)write(*,*)"Execution by thread :",Ilec-30,"opening file ",FicGenot," on unit : ",Ilec
    BLect1:do 
       read(Ilec,'(a1000000)',iostat=io)Ligne
       if(io/=0)exit
       nligne=nligne+1
       call Nb_Element_Chaine(Ligne,NbElt)
       read(Ligne,*)IrcCourant
       !Test du type de fichier lu-phases-ou-typages-(fait-une-seule-fois)------------------------
       if(nligne==2)then
          if(IrcCourant==MaxIrcGen)then
             pardim(2)=2
             NbElt=NbElt-2
             if(DbgMode)write(*,*)"Genotype File ",trim(FicGenot)," contains phased data"
          else
             pardim(2)=1
             if(DbgMode)write(*,*)"Genotype File ",trim(FicGenot)," contains unphased data"
             NbElt=(NbElt-1)/2
          end if
          !Nb of SNP in the file
          pardim(4)=NbElt
       endif!------------------------------------------------------------------------------------
       if(IrcCourant>MaxIrcGen)MaxIrcGen=IrcCourant
       !Si on veut renumeroter les identifiants---------------------------
       if(present(Id))then
          if(IrcCourant>size(Id))cycle BLect1 !cas ou l'animal n'est pas dans le fichier de perf
          if(Id(IrcCourant)>0)then
             NbAnimGardes=NbAnimGardes+1
          else
             cycle BLect1 !on passe a un autre enregistrement          
          endif
       end if
    end do BLect1
    !======================================================================
    !End if first read
    !======================================================================
    !Calcul du nombre d individus genotype
    pardim(1)=MaxIrcGen
    !Si les animaux sont recodes
    if(present(Id))pardim(1)=NbAnimGardes/pardim(2)
    !Allocate a matrice to store phase (1=>SNP;2=>Phase;3=>Individual)
    allocate(pheno(pardim(4),2,MaxIrcGen),allcar(pardim(4)))
    pheno=-1 ; NligneTyp=0

    !======================================================================
    !Now start the second read 
    !======================================================================
    rewind(Ilec)
    BLect:do 
       if(pardim(2)==2)then !lecture fichier phase------------------------------------------------
          read(Ilec,*,iostat=io)ircTyp,NumPhase,AllCar(:)
          if(io/=0)exit BLect
          !if(ircTyp>size(pheno,3))exit BLect
          NligneTyp=NligneTyp+1
          !d'une ligne à l'autre numphase devient 1 puis 2
          if(NumPhase==0)then
             NumPhase=mod(PhaseForce,2)+1
             PhaseForce=NumPhase
          endif
       else   !lecture d un fichier de genotypes
          read(Ilec,*,iostat=io)ircTyp,AllCar(:)
          if(io/=0)exit BLect
          !if(ircTyp>size(pheno,1))exit BLect
          NligneTyp=NligneTyp+1
       end if
       !Si on veut renumeroter les identifiants---------------------------
       if(present(Id))then
          if(ircTyp>size(Id))cycle BLect !cas ou l'animal n'est pas connu (/a renumeroter)
          if(Id(ircTyp)>0)then
             ircTyp=Id(ircTyp)
          else
             cycle BLect !on passe à un autre enregistrement          
          endif
       end if
       !----------------------------------------------------------------------
       !on lit la ligne de typage elt par elt on determine le numero marqueur et la phase---
       do CptElt=1,size(AllCar)
          if(pardim(2)==2)then !si phase
             NumMarq=pardim(3)-1+CptElt
          else !si typage
             NumMarq=pardim(3)-1+((CptElt+1)/2)
             if(mod(CptElt,2)==1)then
                NumPhase=1
             else
                NumPhase=2
             endif
          endif
          Numal=AllCar(CptElt)
          !Et on stocke---------------------------------------------------
          if(ircTyp>size(pheno,3).or.(NumMarq)>size(pheno,1))then
             write(*,*)"Problem Individual",ircTyp," file ",FicGenot,"line",NligneTyp,"thread ",ilec-30
             cycle BLect
          endif
          if(size(pheno,2)==2)then
             if(pheno(NumMarq,NumPhase,ircTyp)==-1)then
                pheno(NumMarq,NumPhase,ircTyp)=Numal
             else
                !$ write(*,*)"Pb Thread",OMP_GET_THREAD_NUM()
                write(*,*)"-Problem : marker ",NumMarq, "phase ",NumPhase
                write(*,*)"-of individual ",ircTyp," already present"
                write(*,*)"-Program stop at line :",NligneTyp," of ",trim(FicGenot)
                stop
             endif
          end if
       end do
    end do BLect
    deallocate(AllCar)
    close(Ilec)    
    !-Calcul de R
    do Ani1=1,size(pheno,3)
       do Ani2=Ani1,size(pheno,3)
          SumEqui=0 ! Mise a zero nb hits
          SumEqui=count(pheno(:,1,Ani1)==pheno(:,1,Ani2))+&
               &count(pheno(:,1,Ani1)==pheno(:,2,Ani2))+&
               &count(pheno(:,2,Ani1)==pheno(:,1,Ani2))+&
               &count(pheno(:,2,Ani1)==pheno(:,2,Ani2))
          !$OMP ATOMIC
          R(IAd(Ani1,Ani2))=R(IAd(Ani1,Ani2))+(0.25_dp*SumEqui)
       enddo
    enddo
  end subroutine LecFicTypage

  subroutine LecDimTypage(FicGenot,pardim,Id)
    !=============================================================================!
    ! Calcul des dimensions des fichiers de typages
    !=============================================================================!
    !sous soutine permettant de lire un fichier de typages, de déterminer son type
    !de vérifier sa cohérence
    !-----------------------------------------------------------------------------!
    !FicGenot   : Un fichier de typage (phase ou typage)
    !pardim     : Un vecteur de dimensions des typages
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
    !$USE OMP_LIB
    implicit none
    character(len=*),intent(in)     ::FicGenot   
    integer,dimension(4),intent(inout)::pardim
    integer,dimension(:),intent(in),optional::Id 

    !Variables internes
    Logical                 ::Ex_Genot
    character(len=1000000)  ::Ligne
    integer                 ::nligne,MaxIrcGen,io,Ilec,OMP_GET_THREAD_NUM
    integer                 ::IrcCourant,NbElt,NbAnimGardes
    integer(kind=1),allocatable::pheno(:,:,:),Allcar(:)
    integer::NligneTyp
    !vérifications d'usage
    inquire(file=FicGenot,EXIST=Ex_Genot)
    if(.not.(Ex_Genot))then
       write(*,*)"No genotype file named ",trim(FicGenot)
       write(*,*)"GLASCOW will now stop"
       stop
    else
       if(DbgMode)write(*,*)"Genotype file named ",trim(FicGenot)," exist"
    endif
    Ilec=0 ;
    !ajouter une exécution cdtionnelle en openMP
    !$ Ilec=OMP_GET_THREAD_NUM()
    !ouverture pour première lecture
    Ilec=30+Ilec ;nligne=0 ; MaxIrcGen=0;pardim=1;NbAnimGardes=0!MAZ
    open(Ilec,file=FicGenot,action="read")
    !$ if(DbgMode)write(*,*)"Execution by thread :",Ilec-30,"openning file ",FicGenot," unit : ",Ilec
    BLect1:do 
       read(Ilec,'(a1000000)',iostat=io)Ligne
       if(io/=0)exit
       nligne=nligne+1
       call Nb_Element_Chaine(Ligne,NbElt)
       read(Ligne,*)IrcCourant
       !Test du type de fichier lu-phases-ou-typages-(fait-une-seule-fois)------------------------
       if(nligne==2)then
          if(IrcCourant==MaxIrcGen)then
             pardim(2)=2
             NbElt=NbElt-2
             if(DbgMode)write(*,*)"Genotype File ",trim(FicGenot)," contains phased data"
          else
             pardim(2)=1
             if(DbgMode)write(*,*)"Genotype File ",trim(FicGenot)," contains unphased data"
             NbElt=(NbElt-1)/2
          end if
          !Nb of SNP in the file
          pardim(4)=NbElt
       endif!------------------------------------------------------------------------------------
       if(IrcCourant>MaxIrcGen)MaxIrcGen=IrcCourant
       !Si on veut renumeroter les identifiants---------------------------
       if(present(Id))then
          if(IrcCourant>size(Id))cycle BLect1 !cas ou l'animal n'est pas dans le fichier de perf
          if(Id(IrcCourant)>0)then
             NbAnimGardes=NbAnimGardes+1
          else
             cycle BLect1 !on passe à un autre enregistrement          
          endif
       end if
    end do BLect1
    !Calcul du nombre d individus genotype
    pardim(1)=MaxIrcGen
    !Si les animaux sont recodes
    if(present(Id))pardim(1)=NbAnimGardes/pardim(2)
    !Allocate a matrice to store phase (1=>SNP;2=>Phase;3=>Individual)
    allocate(pheno(pardim(4),2,MaxIrcGen),allcar(pardim(4)))
    pheno=-1 ; NligneTyp=0

  end subroutine LecDimTypage


  !=============================================================================!
  subroutine Nb_Element_Chaine(Chaine,NbElt)
    !=============================================================================!
    ! Nombre d'éléments d'une chaine de caracteres
    !=============================================================================!
    ! cette sous routine sert a calculer le nombre d'elements presents dans
    !une chaine de caracteres
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
    !Chaine  :Une chaine de caractere 
    !NbElt   :Le nombre d'elements contenus dans la chaine
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
    character (len=*)::Chaine
    character (len=len_trim(Chaine))::Chaine_Red
    integer::InChaine,NbElt,EspCstf

    !on va debarasser Chaine de ses espaces inutiles a droite et a gauche
    Chaine_Red=trim(adjustl(Chaine))
    !verification
    if(len_trim(Chaine_Red)==0)then
       write(*,*)"Empty line"
    else 
       !init 
       NbElt=1 ;EspCstf=0
       do InChaine=1,len_trim(Chaine_Red)
          if((Chaine_Red(InChaine:InChaine))==" ")then!si on est sur 1 blanc
             if(EspCstf==0)then ! et si la case precedente n'etait pas blanche
                NbElt=NbElt+1 ! on avait une variable en plus
                EspCstf=1
             else
                EspCstf=EspCstf+1
             endif
          else
             ! si on est pas sur blanc
             EspCstf=0
          endif
       end do
    endif
  end subroutine Nb_Element_Chaine
  !======================================================================
  integer function IAd(i,j)
    !----------------------------------------------------------------------
    !function to compute address of element(i,j)
    !of a symetric matrix stored in a upper triangular
    !packed form
    !----------------------------------------------------------------------
    integer::i,j
    IAd=i+((j-1)*j)/2
  end function IAD
  !======================================================================

  subroutine ConstructPedAInv(FicPed,R)
    !----------------------------------------------------------------------
    ! Compute A inverse based on a pedigree file
    !----------------------------------------------------------------------
    ! o Ficped : Pedigree file
    ! o R      : The inverse of relationship matrice
    !----------------------------------------------------------------------
    implicit none
    !----------------------------------------------------------------------
    real,allocatable,intent(inout)::R(:,:)
    character(len=*),intent(in)::FicPed
    !----------------------------------------------------------------------
    integer::p(3)
    integer::MaxId,io,Nani,Ani,k,l
    real::w(3),res(4),mendel
    !Define coefficient to apply
    w=(/1., -.5, -.5/)
    res=(/2., 4/3., 1., 0./)
    !Define name of the pedigree file


    Nani=0 ;MaxId=0!MAZ
    !Count Number of individuals
    open(10,file=FicPed)
    do 
       read(10,*,iostat=io)p(1:3)
       if(io/=0)exit
       Nani=Nani+1
       MaxId=max(MaxId,maxval(p))
    end do
    p=0
    allocate(R(MaxId,MaxId)); R=0.0
    rewind(10)
    do 
       read(10,*,iostat=io)p(1:3)
       if(io/=0)exit
       mendel=res(3-count(p(2:3)>0))
       do k=1,3
          do l=1,3	
             if (p(k) /=0 .and. p(l) /=0) then
                R(p(k),p(l))=R(p(k),p(l))+w(k)*w(l)*mendel
             endif
          enddo
       enddo
    end do

    !do Nani=1,size(R,1)
    !   do k=1,size(R,1)
    !      write(*,'(f7.4,1x)',ADVANCE="NO")R(Nani,k)
    !   end do
    !   write(*,'(1x)',ADVANCE="YES")
    !end do

  end subroutine ConstructPedAInv

end program GLMM
