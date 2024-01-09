!!!************************************************************
!!!*						                	                *
!!!*           K Epsilon 1D, Hwang and Lin Model    	         	    *
!!!*							        	                    *
!!!*              Author: Dr. Andrea Boghi  		        *
!!!*							        	                    *
!!!*								                            *
!!!************************************************************
    
Program main_K_epsilon
    implicit none
    real*8, allocatable :: y(:),U(:),kt(:),eps(:),detady(:),d2etady2(:)
    real*8, allocatable :: U0(:),kt0(:),eps0(:)
    real*8, allocatable :: nut(:),dnutdy(:),dUdy(:)
    real*8, allocatable :: tau_mu(:),tau_R(:),Pk(:),Tk(:),Dk(:),Pik(:),Pieps(:),sigmak(:),sigmae(:)
    real*8, allocatable :: Peps(:),Teps(:),Deps(:),epseps(:)
    real*8, allocatable :: T(:),Th2(:),T0(:),Th20(:),lambda(:),dlambdadT(:),d2lambdadT2(:)
    real*8, allocatable :: dTdy(:),d2Tdy2(:),dTh2dy(:),d2Th2dy2(:)
    real*8, allocatable :: q_lam(:),q_R(:),q_new(:),P_Th2(:),eps_Th2(:),T_Th2(:),D_Th2(:),H_Th2(:)
    real*8, allocatable :: dsigmakdy(:),dsigmaedy(:),eps_hat(:),K_K(:),K_eps(:),dK_Kdy(:),dK_epsdy(:)
    integer j,ny,nhy,iter,niter
    real*8 Re_tau,Pr,Bp,Cp,Dp,Ep,dy_min,Cmu,Ce1,Ce2,f1,f2,alphaU,alphaKt,alphaeps
    real*8 resU,resK,resE,resT,resTh2,deta,aU_w,aU_e,sU,aK_w,aK_e,sK,aE_w,aE_e,sE, conv_fac
    real*8 resU_old,resK_old,resE_old,resT_old,resTh2_old
    real*8 aT_e,aT_w,sT,aTh2_e,aTh2_w,sTh2
    real*8 sigmaT,sigmaTh2,alphaT,alphaTh2,U_max
    CHARACTER(len=80)::fname_ke
    CHARACTER(len=80)::fname_th
    CHARACTER(len=80)::fname_res
    logical flag

    open(1,file='imp_ke_th_hl.dat')
    read(1,*) flag
    read(1,*) nhy
    read(1,*) niter
    read(1,*) Re_tau
    read(1,*) dy_min
    read(1,*) Pr
    read(1,*) Bp
    read(1,*) Cp
    read(1,*) Dp
    read(1,*) Ep
    read(1,*) alphaU
    read(1,*) alphaKt
    read(1,*) alphaeps 
    read(1,*) alphaT
    read(1,*) alphaTh2
    close(1)

    ny = nhy*2

    print*, ' niter =', niter, ' ny =', ny 
    print*, ' Re_tau =', Re_tau, ' flag =', flag 

    allocate(y(1:ny),U(1:ny),kt(1:ny),eps(1:ny),nut(1:ny),dnutdy(1:ny),dUdy(1:ny))
    allocate(U0(1:ny),kt0(1:ny),eps0(1:ny),detady(1:ny),d2etady2(1:ny))
    allocate(T(1:ny),Th2(1:ny),T0(1:ny),Th20(1:ny),lambda(1:ny),dlambdadT(1:ny),d2lambdadT2(1:ny))
    allocate(dTdy(1:ny),d2Tdy2(1:ny),dTh2dy(1:ny),d2Th2dy2(1:ny),sigmak(1:ny),sigmae(1:ny))
    allocate(q_lam(1:ny),q_R(1:ny),q_new(1:ny),P_Th2(1:ny),eps_Th2(1:ny),T_Th2(1:ny),D_Th2(1:ny),H_Th2(1:ny))
    allocate(dsigmakdy(1:ny),dsigmaedy(1:ny),Pik(1:ny),Pieps(1:ny),eps_hat(1:ny))
    allocate(K_K(1:ny),K_eps(1:ny),dK_Kdy(1:ny),dK_epsdy(1:ny))

    call initialization(flag,ny,Re_tau,Pr,Bp,Cp,dy_min,y,detady,d2etady2,U,Kt,eps,T,Th2,deta)

    call hwang_lin_k_epsilon_constants(Ce1,Ce2,Cmu,f1,f2)

    sigmaT = 1.d0
    sigmaTh2 = 1.d0

    conv_fac = 1.d0

    write(fname_res,110)'residualsNT_var_Re_tau=',int(Re_tau),'_log10(Pr)=',int(log10(Pr)),'.csv'
    110 format(a23,i3,a11,i2,a4)

    open(11,file=fname_res)
    write(11,*) '"iter","resU","resK","resE","resT","resTh2"'

    do iter=1,niter

        resU_old = resU
        resK_old = resK
        resE_old = resE
        resT_old = resT
        resTh2_old = resTh2

        U0 = U
        Kt0 = kt
        eps0 = eps
        T0 = T
        Th20 = Th2

        do j=1,ny
            call thernal_diffusivity(lambda(j),dlambdadT(j),d2lambdadT2(j),T(j),Bp,Cp,Dp,Ep)
            if(lambda(j) < 0) then
                print*, 'Warning: at j= ',j,' lambda<0'
                lambda(j) = -lambda(j)
            endif
        enddo

        call hwang_lin_k_epsilon_functions(nut,sigmak,sigmae,eps_hat,K_K,K_eps,dK_Kdy,dK_epsdy,ny,y,kt,eps,detady,d2etady2,deta,Cmu)
        
        call relevant_derivatives(dnutdy,dsigmakdy,dsigmaedy,dUdy,dTdy,dTh2dy,d2Tdy2,d2Th2dy2,nut,sigmak,sigmae,U,T,Th2, &
            detady,d2etady2,deta,ny)

        U_max = maxval(U, dim=1, mask=(U>0))

        call solve_u(U,nut,dnutdy,detady,d2etady2,deta,Re_tau,ny)    
        call solve_Kt(Kt,eps,dUdy,nut,dnutdy,detady,d2etady2,deta,sigmak,dsigmakdy,K_eps,dK_epsdy,eps_hat,ny)  
        !do j=1,10
        !    print*, ' j= ',j,'; U(j) = ',U(j),'; Kt(j) = ',Kt(j),'; eps(j) = ',eps(j)
            !print*, ' j= ',j,'; nut(j) = ',nut(j),'; dnutdy(j) = ',dnutdy(j),'; sigmak(j) = ',sigmak(j),'; dsigmakdy(j) = ', &
            !dsigmakdy(j),'; sigmae(j) = ',sigmae(j),'; dsigmaedy(j) = ',dsigmaedy(j)
        !enddo  
        call solve_eps(Kt,eps,dUdy,nut,dnutdy,detady,d2etady2,deta,sigmae,dsigmaedy,K_K,dK_Kdy,eps_hat,ce1,ce2,f1,f2,ny)
        call solve_T(T,nut,dnutdy,lambda,dlambdadT,d2lambdadT2,dTh2dy,d2Th2dy2,dTdy,Pr,sigmaT,deta,d2etady2,detady,ny)
        call solve_Th2(Th2,nut,dnutdy,lambda,dlambdadT,d2lambdadT2,dTdy,d2Tdy2,Pr,sigmaTh2,deta,d2etady2,detady,eps,Kt,ny)

        call residuals(ny,U,U0,resU)
        call residuals(ny,Kt,Kt0,resK)
        call residuals(ny,eps,eps0,resE)
        call residuals(ny,T,T0,resT)
        call residuals(ny,Th2,Th20,resTh2)
        write(11,102) conv_fac*iter,',',resU,',',resK,',',resE,',',resT,',',resTh2

        U = alphaU*U +(1.d0-alphaU)*U0
        Kt = alphaKt*Kt +(1.d0-alphaKt)*Kt0
        eps = alphaeps*eps +(1.d0-alphaeps)*eps0
        !!! T can be negative
        T = alphaT*T +(1.d0-alphaT)*T0
        Th2 = alphaTh2*Th2 +(1.d0-alphaTh2)*Th20
        print*, ' completed =', 100*real(iter)/real(niter), ' resU = ', resU, ' resK = ', resK, ' resE = ', resE, &
        ' resT = ', resT, ' resTh2 = ', resTh2

        !if(iter > 1) then
        !    call correct_residuals(alphaU,resU,resU_old)
        !    call correct_residuals(alphaKt,resK,resK_old)
        !    call correct_residuals(alphaeps,resE,resE_old)
        !    call correct_residuals(alphaT,resT,resT_old)
        !    call correct_residuals(alphaTh2,resTh2,resTh2_old)
        !endif
        
    enddo
    close(11)

    open(512,file='point_ke_var.dat',form='unformatted')
    write(512) y,detady,d2etady2,U,Kt,eps,nut,f2
    close(512)

    open(512,file='point_T_var.dat',form='unformatted')
    write(512) y,T,Th2
    close(512)

    allocate(tau_mu(1:ny),tau_R(1:ny),Pk(1:ny),Tk(1:ny),Dk(1:ny))
    allocate(Peps(1:ny),Teps(1:ny),Deps(1:ny),epseps(1:ny))

    call output_fields_k_eps(ny,U,Kt,eps,nut,f1,f2,deta,sigmaK,sigmaE,Ce1,Ce2,tau_mu,tau_R,Pk,Tk,Dk,Peps,Teps,Deps,epseps, &
        dsigmakdy,dsigmaedy, detady, d2etady2, Pik, Pieps)

    write(fname_ke,111)'momentumNT_var_Re_tau=',int(Re_tau),'_log10(Pr)=',int(log10(Pr)),'.csv'
    111 format(a22,i3,a11,i2,a4)

    open(14,file=fname_ke,form='formatted')
    write(14,*) '"y","U","k","epsilon","nut","tau_mu","tau_R","Pk","Tk","Dk","Peps","Teps","Deps","epsEps","Pik","Pieps","epsk"'
    do j=1,ny
       write(14,101) y(j),',',U(j),',',Kt(j),',',eps(j),',',nut(j),',',tau_mu(j),',',tau_R(j),',',Pk(j),',',Tk(j),',',Dk(j),',', &
       Peps(j),',',Teps(j),',',Deps(j),',',epsEps(j),',',Pik(j),',',Pieps(j),',',-(eps(j)+eps_hat(j))
    enddo
    close(14)

    call output_fields_thermal(ny,T,Th2,lambda,dlambdadT,d2lambdadT2,Kt,eps,nut,detady,d2etady2,deta,sigmaT,sigmaTh2,Pr, &
    q_lam,q_R,q_new,P_Th2,eps_Th2,T_Th2,D_Th2,H_Th2)

    write(fname_th,112)'thermalNT_var_Re_tau=',int(Re_tau),'_log10(Pr)=',int(log10(Pr)),'.csv'
    112 format(a21,i3,a11,i2,a4)

    open(15,file=fname_th,form='formatted')
    write(15,*) '"y","T","Th2","q_lam","q_R","q_new","P_Th2","eps_Th2","T_Th2","D_Th2","H_Th2","lambda"'
    do j=1,ny
       write(15,103) y(j),',',T(j),',',Th2(j),',',q_lam(j),',',q_R(j),',',q_new(j),',',P_Th2(j),',',eps_Th2(j),',',T_Th2(j), & 
       ',',D_Th2(j),',',H_Th2(j),',',lambda(j)
    enddo
    close(15)

    101 format(e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10, &
    A,e18.10,A,e18.10,A,e18.10,A,e18.10)

    102 format(e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10)

    103 format(e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10)

    end

!!!*************************************************
!!!*						         	             *
!!!*            initialization                        *
!!!*								                 *
!!!*************************************************

subroutine initialization(flag,ny,Re_tau,Pr,Bp,Cp,dy_min,y,detady,d2etady2,U,Kt,eps,T,Th2,deta)
    implicit none
    logical, intent(in) :: flag
    integer, intent(in) :: ny
    real*8, intent(in) :: Re_tau,dy_min,Pr,Bp,Cp
    real*8, intent(out) :: y(1:ny),detady(1:ny),d2etady2(1:ny),U(1:ny),kt(1:ny),T(1:ny),Th2(1:ny),eps(1:ny),deta
    integer j
    real*8 Kappa, Cmu,Ce1,Ce2,nut(1:ny),f2(1:ny),y_mid,dUdeta(1:ny),dUdy(1:ny),uk(1:ny),dukdeta(1:ny)

    Kappa = 4.d-1
    Cmu = 9.d-2
    Ce1=1.45d0
    Ce2=1.9d0
    y_mid = 11.635

    !initial conditions
    if (flag) then
        
        print *,'continuation'
  
        open(511,file='point_ke_var.dat',form='unformatted')
        read(511) y(1:ny),detady(1:ny),d2etady2(1:ny),U(1:ny),Kt(1:ny),eps(1:ny),nut(1:ny),f2(1:ny)
        close(511)        
        
        open(511,file='point_T_var.dat',form='unformatted')
        read(511) y(1:ny),T(1:ny),Th2(1:ny)
        close(511)

        deta = 2.d0*Re_tau/(ny-1)

        !T = T - T(1)
  
    else
  
        print *,'new run'

        call grid(ny,dy_min,Re_tau,y,detady,d2etady2,deta)

        do j=1,ny/2
            if(y(j)<=y_mid) then
                U(j) = y(j)
                eps(j) = (y(j)*Kappa -y(j)*y(j)*Kappa/Re_tau)/(Kappa*y_mid)**2.d0 ! 0.1335d0 / ( 1.d0 + ( ( y(j) - 15.515d0 )**2.d0 ) / 166.7634d0 ) !
                Kt(j) = dsqrt( (y(j)*Kappa -y(j)*y(j)*Kappa/Re_tau)*eps(j)/Cmu ) ! 0.019678d0 * y(j) * y(j) / ( 1.d0 + ( ( y(j) - 7.28d0 )**2.d0 ) / 88.263d0 ) !
            else
                U(j) = (1.d0/Kappa)*dlog(y(j)) +5.5d0 
                eps(j) = (y(j)*Kappa -y(j)*y(j)*Kappa/Re_tau)/(Kappa*y(j))**2.d0 !
                Kt(j) = dsqrt( (y(j)*Kappa -y(j)*y(j)*Kappa/Re_tau)*eps(j)/Cmu ) ! 0.019678d0 * y(j) * y(j) / ( 1.d0 + ( ( y(j) - 7.28d0 )**2.d0 ) / 88.263d0 ) !
            endif
        enddo

        !Kt = 0.01d0*y**2.d0
        !eps = Kt

        call ddeta(ny,U,dUdeta,deta)
        dUdy = dUdeta*detady

        do j=ny/2+1,ny
            U(j) = U(ny+1-j)
            eps(j) = eps(ny+1-j)
            Kt(j) = Kt(ny+1-j)
        enddo

        !!! Hwang and Lin correction
        uk = dsqrt(dabs(Kt))
        call ddeta(ny,uk,dukdeta,deta)

        !eps = eps - 2.d0 * ( dukdeta * detady )**2.d0

        !do j=1,ny/4
        !    print*, ' j= ',j,'; U(j) = ',U(j),'; Kt(j) = ',Kt(j),'; eps(j) = ',eps(j)
        !enddo

        do j=1,ny
            Th2(j) = Kt(j)
        enddo

        !!!********************************************************
        !!!
        !!! dy^+ = -K^+(T)dT^+
        !!!
        !!!********************************************************
  
        do j=1,ny/2
            T(j) = - Pr*( U(j) - U(ny/2) )
        enddo

        do j=ny/2+1,ny
            T(j) = -T(ny+1-j)
        enddo

        !T = T - T(1)

    endif

    end

!!!*************************************************
!!!*						         	             *
!!!*                    ddeta                        *
!!!*								                 *
!!!*************************************************
    
subroutine  ddeta(ny,A,DA,deta)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: A(1:ny),deta
    real*8, intent(out) :: DA(1:ny)
    integer j

    DA(1) =  -(3.d0*A(1) -4.d0*A(2) +A(3))/(2.d0*deta)

    do j=2,ny-1
        DA(j) = (A(j+1) -A(j-1))/(2.d0*deta)
    enddo

    DA(ny) = (3.d0*A(ny) -4.d0*A(ny-1) +A(ny-2))/(2.d0*deta)

    end

!!!*************************************************
!!!*						         	             *
!!!*                    d2deta2                        *
!!!*								                 *
!!!*************************************************

subroutine  d2deta2(ny,A,D2A,deta)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: A(1:ny),deta
    real*8, intent(out) :: D2A(1:ny)
    real*8 deta2
    integer j

    deta2 = deta*deta
    
    D2A(1) =  (12.d0*a(1) -30.d0*a(2) +24.d0*a(3) -6.d0*a(4))/(6.d0*deta2)
    
    do j=2,ny-1
        D2A(j) = (a(j+1) - 2.d0*a(j) + a(j-1))/deta2
    enddo
    
    D2A(ny) = (12.d0*a(ny) -30.d0*a(ny-1) +24.d0*a(ny-2) -6.d0*a(ny-3))/(6.d0*deta2)
    
    end

!!!***************************************************
!!!*						         	               *
!!!*       Hwand and Lin K - Epsilon Constants 	       	   *
!!!*								                   *
!!!***************************************************
    
subroutine  hwang_lin_k_epsilon_constants(Ce1,Ce2,Cmu,f1,f2)
    implicit none
    real*8, intent(out) :: Ce1,Ce2,Cmu,f1,f2

    Ce1=1.44d0
    Ce2=1.92d0
    Cmu=0.09d0
    f1=1.d0
    f2=1.d0
    
    end

!!!***************************************************
!!!*						         	               *
!!!*       Hwang and Lin K - Epsilon Functions 	       	   *
!!!*								                   *
!!!***************************************************
    
subroutine  hwang_lin_k_epsilon_functions(nut,sigmak,sigmae,eps_hat,K_K,K_eps,dK_Kdy,dK_epsdy,ny,y,kt,eps,detady,d2etady2,deta,Cmu)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: y(1:ny),Kt(1:ny),eps(1:ny),detady(1:ny),d2etady2(1:ny),Cmu,deta
    real*8, intent(out) :: nut(1:ny),sigmak(1:ny),sigmae(1:ny),eps_hat(1:ny),K_K(1:ny),K_eps(1:ny), dK_Kdy(1:ny),dK_epsdy(1:ny)
    real*8 y_lambda(1:ny), fmu(1:ny), uk(1:ny), dukdeta(1:ny),dukdy(1:ny),depsdeta(1:ny)
    real*8 dKtdeta(1:ny),d2Ktdeta2(1:ny),deps_hatdeta(1:ny),d2eps_hatdeta2(1:ny),Kt_min,eps_min
    integer j

    Kt_min = 1.d-60

    do j=1,ny/2
        if(dabs(Kt(j))<=Kt_min) then
            y_lambda(j)= y(j) * dsqrt( dabs(eps(j) / Kt_min) )
        else
            y_lambda(j)= y(j) * dsqrt( dabs(eps(j) / Kt(j)) )
        endif
    enddo
    !print*, 'y_lambda(1) = ', y_lambda(1), '; y_lambda(2) = ', y_lambda(2), '; y_lambda(3) = ', y_lambda(3), &
    !'; y_lambda(4) = ', y_lambda(4)

    do j=ny/2+1,ny
        if(dabs(Kt(j))<=Kt_min) then
            y_lambda(j)= y(ny+1-j) * dsqrt( dabs(eps(j) / Kt_min) )
        else
            y_lambda(j)= y(ny+1-j) * dsqrt( dabs(eps(j) / Kt(j)) )
        endif
    enddo

    do j=1,ny
        fmu(j)= 1.d0 -dexp(-0.01*y_lambda(j) -0.008d0*y_lambda(j)**3.d0)
    enddo

    eps_min = Kt_min
    do j=1,ny
        if(dabs(eps(j))<=eps_min) then
            nuT(j)= Cmu*fmu(j)*(Kt(j)**2.d0)/eps_min
        else 
            nuT(j)= Cmu*fmu(j)*(Kt(j)**2.d0)/eps(j)
        endif
    enddo

    do j=1,ny
        sigmak(j)= 1.4d0 -1.1d0*dexp(-y_lambda(j)/10.d0)
    enddo

    do j=1,ny
        sigmae(j)= 1.3d0 -1.d0*dexp(-y_lambda(j)/10.d0)
    enddo

    do j=1,ny
        uk(j)= dsqrt(dabs(Kt(j)))
    enddo

    call ddeta(ny,uk,dukdeta,deta)
    dukdy = dukdeta*detady

    eps_hat = 2.d0 * (dukdy)**2.d0

    call ddeta(ny,eps_hat,deps_hatdeta,deta)
    call d2deta2(ny,eps_hat,d2eps_hatdeta2,deta)
    call ddeta(ny,Kt,dKtdeta,deta)
    call d2deta2(ny,Kt,d2Ktdeta2,deta)

    call ddeta(ny,eps,depsdeta,deta)

    do j=1,ny
        if(dabs(eps_hat(j) + eps(j))<=eps_min) then
            K_eps(j) = deps_hatdeta(j) * detady(j) / eps_min
            dK_epsdy(j) = ( d2eps_hatdeta2(j) * detady(j)**2.d0 + deps_hatdeta(j) * d2etady2(j) ) / eps_min &
            - depsdeta(j) * deps_hatdeta(j) * ( detady(j) / eps_min )**2.d0
        else 
            K_eps(j) = deps_hatdeta(j) * detady(j) / (eps_hat(j) + eps(j))
            dK_epsdy(j) = ( d2eps_hatdeta2(j) * detady(j)**2.d0 + deps_hatdeta(j) * d2etady2(j) ) /eps_min &
            - depsdeta(j) * deps_hatdeta(j) * ( detady(j) / eps(j) )**2.d0
        endif

        if(dabs(Kt(j))<=Kt_min) then
            K_K(j) = dKtdeta(j) * detady(j) / Kt_min
            dK_Kdy(j) = ( d2Ktdeta2(j) * detady(j)**2.d0 + dKtdeta(j) * d2etady2(j) ) / Kt_min &
            - ( dKtdeta(j) * detady(j) / Kt_min )**2.d0
        else
            K_K(j) = dKtdeta(j) * detady(j) / Kt(j)
            dK_Kdy(j) = ( d2Ktdeta2(j) * detady(j)**2.d0 + dKtdeta(j) * d2etady2(j) ) / Kt(j) &
            - ( dKtdeta(j) * detady(j) / Kt(j) )**2.d0
        endif
    enddo

    end

!!!***************************************************
!!!*						         	               *
!!!*                Relevant Derivatives	       	   *
!!!*								                   *
!!!***************************************************

subroutine relevant_derivatives(dnutdy,dsigmakdy,dsigmaedy,dUdy,dTdy,dTh2dy,d2Tdy2,d2Th2dy2,nut,sigmak,sigmae,U,T,Th2, &
    detady,d2etady2,deta,ny)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: nut(1:ny),sigmak(1:ny),sigmae(1:ny),U(1:ny),T(1:ny),Th2(1:ny),detady(1:ny),d2etady2(1:ny),deta
    real*8, intent(out) :: dnutdy(1:ny),dsigmakdy(1:ny),dsigmaedy(1:ny),dUdy(1:ny)
    real*8, intent(out) :: dTdy(1:ny),dTh2dy(1:ny),d2Tdy2(1:ny),d2Th2dy2(1:ny)
    real*8 dnutdeta(1:ny),dsigmakdeta(1:ny),dsigmaedeta(1:ny),dUdeta(1:ny),dTdeta(1:ny),dTh2deta(1:ny)
    real*8 d2Tdeta2(1:ny),d2Th2deta2(1:ny)

    call ddeta(ny,nut,dnutdeta,deta)
    dnutdy = dnutdeta*detady
    call ddeta(ny,U,dUdeta,deta)
    dUdy = dUdeta*detady
    call ddeta(ny,sigmak,dsigmakdeta,deta)
    dsigmakdy = dsigmakdeta*detady
    call ddeta(ny,sigmae,dsigmaedeta,deta)
    dsigmaedy = dsigmaedeta*detady
    call ddeta(ny,T,dTdeta,deta)
    dTdy = dTdeta*detady
    call ddeta(ny,Th2,dTh2deta,deta)
    dTh2dy = dTh2deta*detady
    call d2deta2(ny,T,d2Tdeta2,deta)
    d2Tdy2 = d2Tdeta2*detady**2.d0 + dTdeta*d2etady2
    call d2deta2(ny,Th2,d2Th2deta2,deta)
    d2Th2dy2 = d2Th2deta2*detady**2.d0 + dTh2deta*d2etady2

    end

!!!***************************************************
!!!*						         	               *
!!!*                U coefficients	       	   *
!!!*								                   *
!!!***************************************************
subroutine  u_coefficients(aU_w,aU_e,sU,nut,dnutdy,deta,Re_tau,d2etady2,detady)
    implicit none
    real*8, intent(in) :: nut,dnutdy,deta,Re_tau,d2etady2,detady
    real*8, intent(out) :: aU_w,aU_e,sU
    real*8 dev

    dev = deta*( (1.d0+nut)*d2etady2 + dnutdy*detady )/(4.d0*(1.d0+nut)*(detady)**2.d0)

    aU_w = 5.d-1 - dev
    aU_e = 5.d-1 + dev
    sU = (deta*deta)/(2.d0*Re_tau*(1.d0+nut)*(detady)**2.d0)

    end

!!!***************************************************
!!!*						         	               *
!!!*                K coefficients	       	   *
!!!*								                   *
!!!***************************************************
subroutine  K_coefficients(aK_w,aK_e,sK,eps,nut,dnutdy,dUdy,deta,sigmak,dsigmakdy,K_eps,dK_epsdy,eps_hat,d2etady2,detady)
    implicit none
    real*8, intent(in) :: eps,nut,dnutdy,dUdy,deta,sigmak,d2etady2,detady,dsigmakdy,K_eps,dK_epsdy,eps_hat
    real*8, intent(out) :: aK_w,aK_e,sK
    real*8 dev

    dev = deta*( (sigmak+nut)*d2etady2 + ( dnutdy -(nut/sigmak)*dsigmakdy - sigmak * K_eps / 2.d0 )*detady ) &
        /(4.d0*(sigmak+nut)*detady**2.d0)

    aK_w = ( 5.d-1 - dev ) / ( 1.d0 + (deta**2.d0) * dK_epsdy / ( 4.d0 * detady**2.d0) )
    aK_e = ( 5.d-1 + dev ) / ( 1.d0 + (deta**2.d0) * dK_epsdy / ( 4.d0 * detady**2.d0) )
    sK = ( (nut*dUdy*dUdy - eps_hat - eps)*(deta*deta)/(2.d0*(1.d0+nut/sigmak)*detady**2.d0) ) &
        / ( 1.d0 + (deta**2.d0) * dK_epsdy / ( 4.d0 * detady**2.d0) )

    end

!!!***************************************************
!!!*						         	               *
!!!*                E coefficients	       	   *
!!!*								                   *
!!!***************************************************
subroutine  E_coefficients(aE_w,aE_e,sE,eps,Kt,nut,dnutdy,dUdy,deta,sigmae,dsigmaedy,K_K,dK_Kdy,Ce1,f1,Ce2,f2,d2etady2,detady)
    implicit none
    real*8, intent(in) :: eps,Kt,nut,dnutdy,dUdy,deta,sigmae,Ce1,f1,Ce2,f2,d2etady2,detady,dsigmaedy,K_K,dK_Kdy
    real*8, intent(out) :: aE_w,aE_e,sE
    real*8 K_min, Kb, dev
    logical method1

    dev = deta*( (sigmae+nut)*d2etady2 + ( dnutdy -(nut/sigmae)*dsigmaedy - sigmae * K_K )*detady )/(4.d0*(sigmae+nut)*detady**2.d0)

    K_min = 1.d-60

    Kb = (Ce1*f1*nut*dUdy*dUdy -Ce2*f2*eps)*(deta*deta/(2.d0*(1.d0+nut/sigmae)*detady**2.d0))

    method1 = .true.
    if (method1) then
        aE_w = ( 5.d-1 - dev ) / ( 1.d0 + (deta**2.d0) * dK_Kdy / ( 2.d0 * detady**2.d0) )
        aE_e = ( 5.d-1 + dev ) / ( 1.d0 + (deta**2.d0) * dK_Kdy / ( 2.d0 * detady**2.d0) )
        if (dabs(Kt)<=K_min) then
            sE = Kb*eps/K_min
        else
            sE = Kb*eps/Kt
        endif
        sE = sE / ( 1.d0 + (deta**2.d0) * dK_Kdy / ( 2.d0 * detady**2.d0) )
    else
        aE_w = ( (5.d-1 - dev)/(1.d0 - Kb/Kt) ) / ( 1.d0 + (deta**2.d0) * dK_Kdy / ( 2.d0 * detady**2.d0) )
        aE_e = ( (5.d-1 + dev)/(1.d0 - Kb/Kt) ) / ( 1.d0 + (deta**2.d0) * dK_Kdy / ( 2.d0 * detady**2.d0) )
        sE = 0.d0 
    endif
    
    end

!!!***************************************************
!!!*						         	             *
!!!*                T coefficients	       	         *
!!!*								                 *
!!!***************************************************

subroutine  T_coefficients(aT_w,aT_e,sT,nut,dnutdy,lambda,dlambdadT,d2lambdadT2,dTh2dy,d2Th2dy2,dTdy, &
    Pr,sigmaT,deta,d2etady2,detady)
    implicit none
    real*8, intent(in) :: nut,dnutdy,lambda,dlambdadT,d2lambdadT2,dTh2dy,d2Th2dy2,dTdy,Pr,sigmaT,deta,d2etady2,detady
    real*8, intent(out) :: aT_w,aT_e,sT
    real*8 A1,A2,A3, dev

    A1 = lambda / Pr + nut / sigmaT
    A2 = ( dlambdadT * dTdy + d2lambdadT2 * dTh2dy ) / Pr + dnutdy / sigmaT
    A3 =  dlambdadT * d2Th2dy2 / Pr

    !print*, ' A1 = ',A1,' A2 = ',A2,' A3 = ',A3

    dev = deta * ( d2etady2 / detady + A2 / A1 ) / ( 4.d0 * detady )

    aT_w = ( 5.d-1 - dev )
    aT_e = ( 5.d-1 + dev )
    sT = ( ( deta * deta ) / ( 2.d0 * detady**2.d0 ) ) * A3 / A1

    end

!!!***************************************************
!!!*						         	             *
!!!*                Th2 coefficients	       	     *
!!!*								                 *
!!!***************************************************

subroutine  Th2_coefficients(aTh2_w,aTh2_e,sTh2,nut,dnutdy,lambda,dlambdadT,d2lambdadT2,dTdy,d2Tdy2, &
    Pr,sigmaTh2,deta,d2etady2,detady,eps,Kt)
    implicit none
    real*8, intent(in) :: nut,dnutdy,lambda,dlambdadT,d2lambdadT2,dTdy,d2Tdy2,Pr,sigmaTh2,deta,d2etady2,detady
    real*8, intent(in) :: eps,Kt
    real*8, intent(out) :: aTh2_w,aTh2_e,sTh2
    real*8 A1,A2,A3,A4,dev, den

    A1 = lambda / Pr + nut / sigmaTh2
    A2 = 2.d0 * dlambdadT * dTdy / Pr + dnutdy / sigmaTh2
    A3 =  ( nut / sigmaTh2 ) * dTdy**2.d0
    A4 = ( ( eps / Kt ) * lambda - 2.d0 * dlambdadT * d2Tdy2 - 2.d0 *d2lambdadT2 * dTdy**2.d0 ) / Pr 

    dev = deta * ( d2etady2 / detady + A2 / A1 ) / ( 4.d0 * detady )

    den = ( 2.d0 * detady**2.d0 + deta * deta * ( A4 / A1 ) )

    aTh2_w = ( 5.d-1  - dev ) * ( 2.d0 * detady**2.d0 ) / den
    aTh2_e = ( 5.d-1  + dev ) * ( 2.d0 * detady**2.d0 ) / den
    sTh2 = deta * deta * ( A3 / A1 ) / den

    end

!!!*************************************************
!!!*						         	           *
!!!*          output fields k-eps                  *
!!!*								               *
!!!*************************************************

subroutine  output_fields_k_eps(ny,U,Kt,eps,nut,f1,f2,deta,sigmaK,sigmaE,Ce1,Ce2,tau_mu,tau_R,Pk,Tk,Dk,Peps,Teps,Deps,epseps, &
    dsigmakdy, dsigmaedy, detady, d2etady2, Pik, Pieps)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: deta,sigmaK(1:ny),sigmaE(1:ny),Ce1,Ce2,f1,f2,Kt(1:ny),eps(1:ny),nut(1:ny),U(1:ny)
    real*8, intent(in) :: detady(1:ny),d2etady2(1:ny),dsigmakdy(1:ny),dsigmaedy(1:ny)
    real*8, INTENT(OUT) :: tau_mu(1:ny),tau_R(1:ny),Pk(1:ny),Tk(1:ny),Dk(1:ny)
    real*8, INTENT(OUT) :: Peps(1:ny),Teps(1:ny),Deps(1:ny),epseps(1:ny),Pik(1:ny),Pieps(1:ny)
    real*8 dUdy(1:ny),d2Ktdeta2(1:ny),d2epsdeta2(1:ny),dKtdeta(1:ny),depsdeta(1:ny),dnutdy(1:ny),dUdeta(1:ny),dnutdeta(1:ny)
    real*8 uk(1:ny),tau(1:ny),omega(1:ny),eps_hat(1:ny),dtaudeta(1:ny),domegadeta(1:ny),dukdeta(1:ny),deps_hatdeta(1:ny)
    real*8 dukdy(1:ny),d2eps_hatdeta2(1:ny),Kt_min,eps_min
    integer j

    call ddeta(ny,U,dUdeta,deta)
    dUdy = dUdeta*detady
    call ddeta(ny,nut,dnutdeta,deta)
    dnutdy = dnutdeta*detady

    tau_mu = dUdy
    tau_R = nut*dUdy
    Pk = tau_R*dUdy

    call d2deta2(ny,Kt,d2Ktdeta2,deta)
    call ddeta(ny,Kt,dKtdeta,deta)

    Dk = d2Ktdeta2*detady**2.d0 + dKtdeta*d2etady2
    Tk = (nut/sigmaK)*Dk + (dKtdeta*detady/sigmaK)*dnutdy - (nut/sigmaK**2.d0) * dKtdeta * detady * dsigmakdy

    Peps(2:ny) = f1*Ce1*(eps(2:ny)/Kt(2:ny))*Pk(2:ny)
    Peps(1) = Peps(2)

    call d2deta2(ny,eps,D2epsdeta2,deta)
    call ddeta(ny,eps,depsdeta,deta)

    Deps = d2epsdeta2*detady**2.d0 + depsdeta*d2etady2
    Teps = (nut/sigmaE)*Deps + (depsdeta*detady/sigmaE)*dnutdy - (nut/sigmaE**2.d0) * depsdeta * detady * dsigmaedy

    Kt_min = 1.d-60
    do j=1,ny
        if(dabs(Kt(j))<=Kt_min) then
            epsEps(j) = -f2*Ce2*(eps(j)/Kt_min)*eps(j)
        else
            epsEps(j) = -f2*Ce2*(eps(j)/Kt(j))*eps(j)
        endif
    enddo

        do j=1,ny
        uk(j)= dsqrt(dabs(Kt(j)))
    enddo

    call ddeta(ny,uk,dukdeta,deta)
    dukdy = dukdeta*detady

    eps_hat = 2.d0 * (dukdy)**2.d0

    do j=1,ny
        if(dabs(eps_hat(j) + eps(j))<=eps_min) then
            tau(j) = Kt(j) / eps_min
        else 
            tau(j) = Kt(j) / (eps_hat(j) + eps(j))
        endif
    enddo
    do j=1,ny
        if(dabs(Kt(j))<=Kt_min) then
            omega(j) = eps(j) / Kt_min
        else
            omega(j) = eps(j) / Kt(j)
        endif
    enddo
    call ddeta(ny,tau,dtaudeta,deta)
    call ddeta(ny,omega,domegadeta,deta)
    call d2deta2(ny,eps_hat,d2eps_hatdeta2,deta)
    call ddeta(ny,eps_hat,deps_hatdeta,deta)
    !!! This the Problematic term.
    Pik = -0.5d0 * ( tau * (d2eps_hatdeta2*detady**2.d0 + deps_hatdeta*d2etady2) + dtaudeta * deps_hatdeta * (detady)**2.d0 )

    call ddeta(ny,Kt,dKtdeta,deta)
    call d2deta2(ny,Kt,d2Ktdeta2,deta)
    
    Pieps= - ( omega * (d2Ktdeta2*detady**2.d0 + dKtdeta*d2etady2) + domegadeta * dKtdeta * (detady)**2.d0 )
    
    end

!!!*************************************************
!!!*						         	           *
!!!*          output fields thermal                  *
!!!*								               *
!!!*************************************************

subroutine  output_fields_thermal(ny,T,Th2,lambda,dlambdadT,d2lambdadT2,Kt,eps,nut,detady,d2etady2,deta,sigmaT,sigmaTh2,Pr, &
    q_lam,q_R,q_new,P_Th2,eps_Th2,T_Th2,D_Th2,H_Th2)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: deta,Kt(1:ny),eps(1:ny),nut(1:ny),detady(1:ny),d2etady2(1:ny)
    real*8, intent(in) :: T(1:ny),Th2(1:ny),lambda(1:ny),dlambdadT(1:ny),d2lambdadT2(1:ny),sigmaT,sigmaTh2,Pr
    real*8, INTENT(OUT) :: q_lam(1:ny),q_R(1:ny),q_new(1:ny),P_Th2(1:ny),eps_Th2(1:ny),T_Th2(1:ny),D_Th2(1:ny),H_Th2(1:ny)
    real*8 dTdy(1:ny),dTdeta(1:ny),d2Tdy2(1:ny),d2Tdeta2(1:ny),dTh2dy(1:ny),dTh2deta(1:ny),d2Th2dy2(1:ny),d2Th2deta2(1:ny)
    real*8 dnutdy(1:ny),dnutdeta(1:ny)

    call ddeta(ny,T,dTdeta,deta)
    dTdy = dTdeta*detady
    call ddeta(ny,Th2,dTh2deta,deta)
    dTh2dy = dTh2deta*detady
    call ddeta(ny,nut,dnutdeta,deta)
    dnutdy = dnutdeta*detady
    call d2deta2(ny,T,d2Tdeta2,deta)
    d2Tdy2 = d2Tdeta2*detady**2.d0 + dTdeta*d2etady2
    call d2deta2(ny,Th2,d2Th2deta2,deta)
    d2Th2dy2 = d2Th2deta2*detady**2.d0 + dTh2deta*d2etady2

    q_lam = ( lambda / Pr ) * dTdy
    q_R = ( nut / sigmaT ) * dTdy
    q_new = ( dlambdadT / Pr ) * dTh2dy

    P_Th2 = ( nut / sigmaT ) * dTdy**2.d0
    eps_Th2(2:ny-1) = ( eps(2:ny-1) / Kt(2:ny-1) ) * ( lambda(2:ny-1) / Pr ) * Th2(2:ny-1)
    eps_Th2(1) = eps_Th2(2)
    eps_Th2(ny) = eps_Th2(ny-1)

    T_Th2 = ( nut / sigmaTh2 )*d2Th2dy2 + ( dTh2dy / sigmaTh2 ) * dnutdy
    D_Th2 = ( lambda / Pr )*d2Th2dy2 + ( dTh2dy / Pr ) * dlambdadT * dTdy

    H_Th2 = ( 2.d0 / Pr ) * d2lambdadT2 * Th2 * dTdy**2.d0 + ( 2.d0 / Pr ) * dlambdadT * Th2 * d2Tdy2 &
    + ( 1.d0 / Pr ) * dlambdadT * dTh2dy * dTdy
    
    end


!!!*************************************************
!!!*						         	           *
!!!*               residuals                   *
!!!*								               *
!!!*************************************************

subroutine  residuals(ny,A,A0,resA)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: A(1:ny),A0(1:ny)
    real*8, INTENT(OUT) :: resA
    real*8 sumN, sumD
    integer j

    sumN = 0
    do j=1,ny
        sumN = sumN +dabs(A(j)- A0(j))
    enddo

    sumD = 0
    do j=1,ny
        sumD = sumD + dabs(A0(j))
    enddo

    resA = sumN/sumD

    end

!!!*************************************************
!!!*						         	             *
!!!*            grid                       *
!!!*								                 *
!!!*************************************************

subroutine grid(ny,dy_min,Re_tau,y,detady,d2etady2,deta)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: Re_tau,dy_min
    real*8, intent(out) :: y(1:ny),detady(1:ny),d2etady2(1:ny),deta
    integer j
    real*8 a, b, c, d, e, eta, y_max

    y_max = 2*Re_tau

    deta = y_max/(ny-1)

    a = (deta - dy_min)/(Re_tau*deta - deta*deta)
    b = (dy_min*Re_tau - deta*deta)/(Re_tau*deta - deta*deta)

    do j=1,ny/2
        eta = deta*(j-1)
        y(j) = a*eta**2.d0 + b*eta
        detady(j) = 1.d0/(2.d0*a*eta + b)
        d2etady2(j) = -2.d0*a/(2.d0*a*eta + b)**3.d0
    enddo

    c = -a
    d = 4.d0*Re_tau*a + b
    e = 2.d0*Re_tau*( 1.d0 - b -2.d0*a*Re_tau )

    do j=ny/2+1,ny
        eta = deta*(j-1)
        y(j) = c*eta**2.d0 + d*eta + e
        detady(j) = 1.d0/(2.d0*c*eta + d)
        d2etady2(j) = -2.d0*c/(2.d0*c*eta + d)**3.d0
    enddo

    print*, ' dy_max =', y(ny/2)-y(ny/2-1), ' dy_min =', y(2)-y(1), ' ratio =', (y(ny/2)-y(ny/2-1))/(y(2)-y(1))

    end

!!!*************************************************
!!!*						         	           *
!!!*            Thermal Diffusivity                *
!!!*								               *
!!!*************************************************
    
subroutine  thernal_diffusivity(lambda,dlambdadT,d2lambdadT2,T,Bp,Cp,Dp,Ep)
    implicit none
    real*8, intent(in) :: T,Bp,Cp,Dp,Ep
    real*8, intent(out) :: lambda,dlambdadT,d2lambdadT2
    integer j

    lambda = 1.d0 + Bp * T + Cp * T**2.d0 + Dp * T**3.d0 + Ep * T**4.d0
    dlambdadT = Bp + 2.d0 * Cp * T + 3.d0 * Dp * T**2.d0 + 4.d0 * Ep * T**3.d0
    d2lambdadT2 = 2.d0 * Cp + 6.d0 * Dp * T + 12.d0 * Ep * T**2.d0

    end

!!!*************************************************
!!!*						         	           *
!!!*            Thomas Algorithm                *
!!!*								               *
!!!*************************************************
    
subroutine  thomas_algorithm(var,ny,a_e,a_w,a_p,S)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(inout) :: var(1:ny),a_e,a_w,a_p,S
    real*8 A(1:ny),C_apex(1:ny)
    integer j

    A(1) = 0.d0
    A(ny) = 0.d0
    C_apex(1) = var(1)
    C_apex(ny) = var(ny)

    do j=2,ny-1
        A(j) = a_e / ( a_p - a_w * A(j-1) )
        C_apex(j) = ( a_w * C_apex(j-1) + S ) / ( a_p - a_w * A(j-1) )
        var(j) = A(j) * var(j+1) + C_apex(j)
    enddo

    end

!!!*************************************************
!!!*						         	           *
!!!*            Correct Residuals               *
!!!*								               *
!!!*************************************************
    
subroutine  correct_residuals(alpha,res,res_old)
    implicit none
    real*8, intent(in) :: res,res_old
    real*8, intent(inout) :: alpha
    real*8 increment

    increment = 1.d-06

    if(res < res_old) then
        alpha = alpha * ( 1.d0 + increment )
    elseif(res > res_old) then
        alpha = alpha * ( 1.d0 + increment )
    else
        alpha = alpha
    endif

    alpha = min(1.d0,max(alpha,increment))

    end

!!!*************************************************
!!!*						         	           *
!!!*                 Solve U                        *
!!!*								               *
!!!*************************************************
    
subroutine  solve_u(U,nut,dnutdy,detady,d2etady2,deta,Re_tau,ny)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: nut(1:ny),dnutdy(1:ny),detady(1:ny),d2etady2(1:ny),deta,Re_tau
    real*8, intent(inout) :: U(1:ny)
    real*8 aU_w,aU_e,sU
    real*8 A(1:ny),C_apex(1:ny),denominator
    integer j

    U(1) = 0.d0
    call u_coefficients(aU_w,aU_e,sU,nut(2),dnutdy(2),deta,Re_tau,d2etady2(2),detady(2))
    A(2) = aU_e
    C_apex(2) = sU + aU_w * U(1)
    do j =3,ny-1
        call u_coefficients(aU_w,aU_e,sU,nut(j),dnutdy(j),deta,Re_tau,d2etady2(j),detady(j))
        denominator = ( 1.d0 - aU_w * A(j-1) )
        A(j) = aU_e / denominator
        C_apex(j) = ( aU_w * C_apex(j-1) + sU ) / denominator
    enddo
    !call u_coefficients(aU_w,aU_e,sU,nut(ny),dnutdy(ny),deta,Re_tau,d2etady2(ny),detady(ny))
    !denominator = ( 1.d0 - aU_w * A(ny-1) )
    !A(ny) = aU_e / denominator
    !C_apex(ny) = ( aU_w * C_apex(ny-1) + sU ) / denominator
    A(ny-1) = 0.d0
    U(ny) = 0.d0

    !U(1) = 0.d0
    !do j =2,ny-1
    !    call u_coefficients(aU_w,aU_e,sU,nut(j),dnutdy(j),deta,Re_tau,d2etady2(j),detady(j))
    !    U(j) =  sU + aU_e*U(j+1) + aU_w*U(j-1)
    !    !!!U(j) = A(j) * U(j+1) + C_apex(j)
    !enddo
    !U(ny) = 0.d0

    do j =ny-1,2,-1
        U(j) = A(j) * U(j+1) + C_apex(j)
    enddo


    end

!!!*************************************************
!!!*						         	           *
!!!*                 Solve Kt                        *
!!!*								               *
!!!*************************************************
    
subroutine  solve_Kt(Kt,eps,dUdy,nut,dnutdy,detady,d2etady2,deta,sigmak,dsigmakdy,K_eps,dK_epsdy,eps_hat,ny)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: eps(1:ny),dUdy(1:ny),nut(1:ny),dnutdy(1:ny),detady(1:ny),d2etady2(1:ny),deta,sigmak(1:ny)
    real*8, intent(in) :: dsigmakdy(1:ny),K_eps(1:ny),dK_epsdy(1:ny),eps_hat(1:ny)
    real*8, intent(inout) :: Kt(1:ny)
    real*8 aK_w,aK_e,sK
    real*8 A(1:ny),C_apex(1:ny),denominator
    integer j

    Kt(1) = 0.d0
    call K_coefficients(aK_w,aK_e,sK,eps(2),nut(2),dnutdy(2),dUdy(2),deta,sigmak(2),dsigmakdy(2),K_eps(2),dK_epsdy(2), &
        eps_hat(2), d2etady2(2),detady(2))
    A(2) = aK_e
    C_apex(2) = sK + aK_w * Kt(1)
    do j =3,ny-1
        call K_coefficients(aK_w,aK_e,sK,eps(j),nut(j),dnutdy(j),dUdy(j),deta,sigmak(j),dsigmakdy(j),K_eps(j),dK_epsdy(j), &
            eps_hat(j),d2etady2(j),detady(j))
        denominator = ( 1.d0 - aK_w * A(j-1) )
        A(j) = aK_e / denominator
        C_apex(j) = ( aK_w * C_apex(j-1) + sK ) / denominator
        !print*, ' A(j) = ', A(j),' C_apex(j) = ', C_apex(j)
    enddo
    !call K_coefficients(aK_w,aK_e,sK,eps(ny),nut(ny),dnutdy(ny),dUdy(ny),deta,sigmak(ny),dsigmakdy(ny),K_eps(ny),dK_epsdy(ny),eps_hat(ny),d2etady2(ny),detady(ny))
    !denominator = ( 1.d0 - aK_w * A(ny-1) )
    !A(ny) = aK_e / denominator
    !C_apex(ny) = ( aK_w * C_apex(ny-1) + sK ) / denominator
    A(ny-1) = 0.d0
    Kt(ny) = 0.d0
    
    do j =2,ny-1
        call K_coefficients(aK_w,aK_e,sK,eps(j),nut(j),dnutdy(j),dUdy(j),deta,sigmak(j),dsigmakdy(j),K_eps(j),dK_epsdy(j), &
            eps_hat(j), d2etady2(j),detady(j))
        Kt(j) = sK + aK_e*Kt(j+1) + aK_w*Kt(j-1)
    enddo
    
    !do j =ny-1,2,-1
    !    Kt(j) = A(j) * Kt(j+1) + C_apex(j)
    !enddo

    end

!!!*************************************************
!!!*						         	           *
!!!*                 Solve Eps                        *
!!!*								               *
!!!*************************************************
    
subroutine  solve_eps(Kt,eps,dUdy,nut,dnutdy,detady,d2etady2,deta,sigmae,dsigmaedy,K_K,dK_Kdy,eps_hat,ce1,ce2,f1,f2,ny)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: Kt(1:ny),dUdy(1:ny),nut(1:ny),dnutdy(1:ny),detady(1:ny),d2etady2(1:ny),deta,sigmae(1:ny)
    real*8, intent(in) :: ce1,ce2,f1,f2,dsigmaedy(1:ny),K_K(1:ny),dK_Kdy(1:ny),eps_hat(1:ny)
    real*8, intent(inout) :: eps(1:ny)
    real*8 aE_w,aE_e,sE
    real*8 A(1:ny),C_apex(1:ny),denominator
    integer j

    eps(1) = 0.d0
    call E_coefficients(aE_w,aE_e,sE,eps(2),Kt(2),nut(2),dnutdy(2),dUdy(2),deta,sigmae(2),dsigmaedy(2),K_K(2),dK_Kdy(2), &
        Ce1,f1,Ce2,f2, d2etady2(2), detady(2))
    A(2) = aE_e
    C_apex(2) = sE + aE_w * eps(1)
    do j =3,ny-1
        call E_coefficients(aE_w,aE_e,sE,eps(j),Kt(j),nut(j),dnutdy(j),dUdy(j),deta,sigmae(j),dsigmaedy(j),K_K(j),dK_Kdy(j), &
            Ce1,f1,Ce2,f2,d2etady2(j), detady(j))
        denominator = ( 1.d0 - aE_w * A(j-1) )
        A(j) = aE_e / denominator
        C_apex(j) = ( aE_w * C_apex(j-1) + sE ) / denominator
    enddo
    !call E_coefficients(aE_w,aE_e,sE,eps(ny),Kt(ny),nut(ny),dnutdy(ny),dUdy(ny),deta,sigmae(ny),dsigmaedy(ny),K_K(ny),dK_Kdy(ny), &
    !   Ce1,f1,Ce2,f2,d2etady2(ny), detady(ny))
    !denominator = ( 1.d0 - aE_w * A(ny-1) )
    A(ny-1) = 0.d0
    eps(ny) = 0.d0
    
    do j =2,ny-1
        call E_coefficients(aE_w,aE_e,sE,eps(j),Kt(j),nut(j),dnutdy(j),dUdy(j),deta,sigmae(j),dsigmaedy(j),K_K(j),dK_Kdy(j), &
            Ce1,f1,Ce2,f2,d2etady2(j), detady(j))
        !if(j==2 .OR. j==ny-1) then
        !    print*, eps(j),Kt(j),nut(j),dnutdy(j),dUdy(j),deta,sigmae(j),dsigmaedy(j),Pieps(j),d2etady2(j), detady(j)
        !endif
        !PRINT*, ' j = ',j,'; sE = ',sE,'; aE_w = ',aE_w,',aE_e = ',aE_e
        eps(j) = sE + aE_e*eps(j+1) + aE_w*eps(j-1)
    enddo

    !do j =ny-1,2,-1
    !    eps(j) = A(j) * eps(j+1) + C_apex(j)
    !enddo    

    end

!!!*************************************************
!!!*						         	           *
!!!*                 Solve T                       *
!!!*								               *
!!!*************************************************
    
subroutine  solve_T(T,nut,dnutdy,lambda,dlambdadT,d2lambdadT2,dTh2dy,d2Th2dy2,dTdy,Pr,sigmaT,deta,d2etady2,detady,ny)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: dTdy(1:ny),nut(1:ny),dnutdy(1:ny),detady(1:ny),d2etady2(1:ny),deta,sigmaT
    real*8, intent(in) :: Pr,lambda(1:ny),dlambdadT(1:ny),d2lambdadT2(1:ny),dTh2dy(1:ny),d2Th2dy2(1:ny)
    real*8, intent(inout) :: T(1:ny)
    real*8 aT_w,aT_e,sT
    real*8 A(1:ny),C_apex(1:ny),denominator
    integer j

    call T_coefficients(aT_w,aT_e,sT,nut(1),dnutdy(1),lambda(1),dlambdadT(1),d2lambdadT2(1),dTh2dy(1),d2Th2dy2(1),dTdy(1), &
        Pr,sigmaT,deta,d2etady2(1),detady(1))
    T(1) = sT + aT_e*T(2) + aT_w*( T(2) - 2.d0 * deta * Pr / ( lambda(1) * detady(1) ) )
    A(1) = aT_e + aT_w
    C_apex(1) = sT - aT_w * 2.d0 * deta * Pr / ( lambda(1) * detady(1) ) 
    do j =2,ny
        call T_coefficients(aT_w,aT_e,sT,nut(j),dnutdy(j),lambda(j),dlambdadT(j),d2lambdadT2(j),dTh2dy(j),d2Th2dy2(j),dTdy(j), &
        Pr,sigmaT,deta,d2etady2(j),detady(j))
        denominator = ( 1.d0 - aT_w * A(j-1) )
        A(j) = aT_e / denominator
        C_apex(j) = ( aT_w * C_apex(j-1) + sT ) / denominator
    enddo
    !call T_coefficients(aT_w,aT_e,sT,nut(ny),dnutdy(ny),lambda(ny),dlambdadT(ny),d2lambdadT2(ny),dTh2dy(ny),d2Th2dy2(ny), & 
    !dTdy(ny), Pr,sigmaT,deta,d2etady2(ny),detady(ny))
    !denominator = ( 1.d0 - aT_w * A(ny-1) )
    A(ny) = 0.d0
    !C_apex(ny) = ( aT_w * C_apex(ny-1) + sT ) / denominator

    !do j =2,ny-1
    !    call T_coefficients(aT_w,aT_e,sT,nut(j),dnutdy(j),lambda(j),dlambdadT(j),d2lambdadT2(j),dTh2dy(j),d2Th2dy2(j),dTdy(j), &
    !    Pr,sigmaT,deta,d2etady2(j),detady(j))
    !    T(j) = sT + aT_e*T(j+1) + aT_w*T(j-1)
    !enddo
    call T_coefficients(aT_w,aT_e,sT,nut(ny),dnutdy(ny),lambda(ny),dlambdadT(ny),d2lambdadT2(ny),dTh2dy(ny),d2Th2dy2(ny), & 
    dTdy(ny), Pr,sigmaT,deta,d2etady2(ny),detady(ny))
    T(ny) = sT + aT_w*T(ny-1) + aT_e*( T(ny-1) + 2.d0 * deta * Pr / ( lambda(ny) * detady(ny) ) ) 
    
    do j =ny-1,1,-1
        T(j) = A(j) * T(j+1) + C_apex(j)
    enddo

    if (mod(ny,2).eq.0) then
        T = T - ( T(ny/2) + T(ny/2+1) ) / 2.d0
    else
        T = T - T((ny+1)/2)
    endif
    !T = T - T(1)

    end

!!!*************************************************
!!!*						         	           *
!!!*                 Solve Th2                       *
!!!*								               *
!!!*************************************************
    
subroutine  solve_Th2(Th2,nut,dnutdy,lambda,dlambdadT,d2lambdadT2,dTdy,d2Tdy2,Pr,sigmaTh2,deta,d2etady2,detady,eps,Kt,ny)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: dTdy(1:ny),nut(1:ny),dnutdy(1:ny),detady(1:ny),d2etady2(1:ny),deta,sigmaTh2
    real*8, intent(in) :: Pr,lambda(1:ny),dlambdadT(1:ny),d2lambdadT2(1:ny),d2Tdy2(1:ny),Kt(1:ny),eps(1:ny)
    real*8, intent(inout) :: Th2(1:ny)
    real*8 aTh2_w,aTh2_e,sTh2
    real*8 A(1:ny),C_apex(1:ny),denominator
    integer j

    Th2(1) = 0.d0
    call Th2_coefficients(aTh2_w,aTh2_e,sTh2,nut(2),dnutdy(2),lambda(2),dlambdadT(2),d2lambdadT2(2),dTdy(2),d2Tdy2(2), &
        Pr,sigmaTh2,deta,d2etady2(2),detady(2),eps(2),Kt(2))
    A(2) = aTh2_e
    C_apex(2) = sTh2 + aTh2_w * Th2(1)
    do j =3,ny-1
        call Th2_coefficients(aTh2_w,aTh2_e,sTh2,nut(j),dnutdy(j),lambda(j),dlambdadT(j),d2lambdadT2(j),dTdy(j),d2Tdy2(j), &
        Pr,sigmaTh2,deta,d2etady2(j),detady(j),eps(j),Kt(j))
        denominator = ( 1.d0 - aTh2_w * A(j-1) )
        A(j) = aTh2_e / denominator
        C_apex(j) = ( aTh2_w * C_apex(j-1) + sTh2 ) / denominator
    enddo
    !call Th2_coefficients(aTh2_w,aTh2_e,sTh2,nut(ny),dnutdy(ny),lambda(ny),dlambdadT(ny),d2lambdadT2(ny),dTdy(ny),d2Tdy2(ny), &
    !    Pr,sigmaTh2,deta,d2etady2(ny),detady(ny),eps(ny),Kt(ny))
    !denominator = ( 1.d0 - aTh2_w * A(ny-1) )
    A(ny-1) = 0.d0
    !C_apex(ny) = ( aTh2_w * C_apex(ny-1) + sTh2 ) / denominator
    Th2(ny) = 0.d0
    
    !do j =2,ny-1
    !    call Th2_coefficients(aTh2_w,aTh2_e,sTh2,nut(j),dnutdy(j),lambda(j),dlambdadT(j),d2lambdadT2(j),dTdy(j),d2Tdy2(j), &
    !    Pr,sigmaTh2,deta,d2etady2(j),detady(j),eps(j),Kt(j))
    !    Th2(j) = sTh2 + aTh2_e*Th2(j+1) + aTh2_w*Th2(j-1)
    !enddo

    do j =ny-1,2,-1
        Th2(j) = A(j) * Th2(j+1) + C_apex(j)
    enddo
    

    end
