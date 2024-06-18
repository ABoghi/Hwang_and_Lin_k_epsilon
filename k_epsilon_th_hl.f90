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
    real*8, allocatable :: y(:),U(:),kt(:),eps(:),detady(:),d2etady2(:),d3etady3(:)
    real*8, allocatable :: U0(:),kt0(:),eps0(:)
    real*8, allocatable :: nut(:),dnutdy(:),dUdy(:)
    real*8, allocatable :: tau_mu(:),tau_R(:),Pk(:),Tk(:),Dk(:),Pik(:),Pieps(:),sigmak(:),sigmae(:)
    real*8, allocatable :: Peps(:),Teps(:),Deps(:),epseps(:),deps_hatdy(:)
    real*8, allocatable :: T(:),Th2(:),T0(:),Th20(:),lambda(:),dlambdadT(:),d2lambdadT2(:)
    real*8, allocatable :: dTdy(:),d2Tdy2(:),dTh2dy(:),d2Th2dy2(:)
    real*8, allocatable :: q_lam(:),q_R(:),q_new(:),P_Th2(:),eps_Th2(:),T_Th2(:),D_Th2(:),H_Th2(:)
    real*8, allocatable :: dsigmakdy(:),dsigmaedy(:),eps_hat(:),Pi_K(:),Pi_eps(:),upkt(:),upeps(:)
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

    ny = nhy*2 !+ 1

    print*, ' niter =', niter, ' ny =', ny 
    print*, ' Re_tau =', Re_tau, ' flag =', flag 
    print*, ' alphaU =', alphaU, ' alphaKt =', alphaKt, ' alphaeps =', alphaeps

    allocate(y(1:ny),U(1:ny),kt(1:ny),eps(1:ny),nut(1:ny),dnutdy(1:ny),dUdy(1:ny))
    allocate(U0(1:ny),kt0(1:ny),eps0(1:ny),detady(1:ny),d2etady2(1:ny),deps_hatdy(1:ny))
    allocate(T(1:ny),Th2(1:ny),T0(1:ny),Th20(1:ny),lambda(1:ny),dlambdadT(1:ny),d2lambdadT2(1:ny))
    allocate(dTdy(1:ny),d2Tdy2(1:ny),dTh2dy(1:ny),d2Th2dy2(1:ny),sigmak(1:ny),sigmae(1:ny))
    allocate(q_lam(1:ny),q_R(1:ny),q_new(1:ny),P_Th2(1:ny),eps_Th2(1:ny),T_Th2(1:ny),D_Th2(1:ny),H_Th2(1:ny))
    allocate(dsigmakdy(1:ny),dsigmaedy(1:ny),Pik(1:ny),Pieps(1:ny),eps_hat(1:ny))
    allocate(Pi_K(1:ny),Pi_eps(1:ny),upkt(1:ny),upeps(1:ny),d3etady3(1:ny))

    call initialization(flag,ny,Re_tau,Pr,Bp,Cp,dy_min,y,detady,d2etady2,d3etady3,U,Kt,eps,T,Th2,deta)

    call hwang_lin_k_epsilon_constants(Ce1,Ce2,Cmu,f1,f2)

    sigmaT = 1.d0
    sigmaTh2 = 1.d0

    conv_fac = 1.d0

    write(fname_res,110)'residualsNT_var_Re_tau=',int(Re_tau),'_log10(Pr)=',int(log10(Pr)),'.csv'
    110 format(a23,i3,a11,i2,a4)

    open(11,file=fname_res)
    write(11,*) '"iter","resU","resK","resE","resT","resTh2"'

    call define_eps_hat(eps_hat,kt,detady,d2etady2,deta,ny)

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

        call hwang_lin_k_epsilon_functions(nut,sigmak,sigmae,ny,y,kt,eps,detady,d2etady2,deta,Cmu)
        !if (mod(iter,10000).eq.0) then !1500
        !if (iter > 10000) then
        !    print*, "Update eps_hat"
            call define_eps_hat(eps_hat,kt,detady,d2etady2,deta,ny)
        !endif
        
        call relevant_derivatives(dnutdy,dsigmakdy,dsigmaedy,dUdy,dTdy,dTh2dy,d2Tdy2,d2Th2dy2,nut,sigmak,sigmae,U,T,Th2, &
            detady,d2etady2,deta,ny)

        U_max = maxval(U, dim=1, mask=(U>0))

        call solve_u(U,nut,dnutdy,detady,d2etady2,deta,Re_tau,ny)    
        if (mod(iter,100).eq.0) then !1500
        call pressure_speed_Kt(upkt,Kt,eps,eps_hat,detady,d2etady2,d3etady3,deta,y,ny)
        call ddy(ny,eps_hat,deps_hatdy,detady,deta)
        upkt = deps_hatdy / ( 2.d0 * ( eps + eps_hat ) )
        do j=1,16
            call quadratic_variable_smoother(upkt,upkt,y,ny)
        enddo
        endif
        call solve_Kt_b(Kt,eps,dUdy,nut,dnutdy,detady,d2etady2,d3etady3,deta,sigmak,dsigmakdy,eps_hat,ny,upkt,y)  
        !do j=1,10
        !    print*, ' j= ',j,'; U(j) = ',U(j),'; Kt(j) = ',Kt(j),'; eps(j) = ',eps(j)
            !print*, ' j= ',j,'; nut(j) = ',nut(j),'; dnutdy(j) = ',dnutdy(j),'; sigmak(j) = ',sigmak(j),'; dsigmakdy(j) = ', &
            !dsigmakdy(j),'; sigmae(j) = ',sigmae(j),'; dsigmaedy(j) = ',dsigmaedy(j)
        !enddo  
        call solve_eps(Kt,eps,dUdy,nut,dnutdy,detady,d2etady2,deta,sigmae,dsigmaedy,eps_hat,ce1,ce2,f1,f2,ny,upeps)
        call solve_T(T,nut,dnutdy,lambda,dlambdadT,d2lambdadT2,dTh2dy,d2Th2dy2,dTdy,Pr,sigmaT,deta,d2etady2,detady,ny)
        call solve_Th2(Th2,nut,dnutdy,lambda,dlambdadT,d2lambdadT2,dTdy,d2Tdy2,Pr,sigmaTh2,deta,d2etady2,detady,eps,Kt,ny)

        U = alphaU*U +(1.d0-alphaU)*U0
        Kt = alphaKt*Kt +(1.d0-alphaKt)*Kt0
        eps = alphaeps*eps +(1.d0-alphaeps)*eps0
        !!! T can be negative
        T = alphaT*T +(1.d0-alphaT)*T0
        Th2 = alphaTh2*Th2 +(1.d0-alphaTh2)*Th20

        call residuals(ny,U,U0,resU)
        call residuals(ny,Kt,Kt0,resK)
        call residuals(ny,eps,eps0,resE)
        call residuals(ny,T,T0,resT)
        call residuals(ny,Th2,Th20,resTh2)

        write(11,102) conv_fac*iter,',',resU,',',resK,',',resE,',',resT,',',resTh2

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

    call output_fields_k_eps(ny,U,Kt,eps,nut,f1,f2,deta,sigmaK,sigmaE,Ce1,Ce2,tau_mu,tau_R,Pk,Pi_K,Tk,Dk,Peps,Pi_eps,Teps, &
        Deps,epseps, eps_hat, y, dsigmakdy,dsigmaedy, detady, d2etady2, d3etady3)

    write(fname_ke,111)'momentumNT_var_Re_tau=',int(Re_tau),'_log10(Pr)=',int(log10(Pr)),'.csv'
    111 format(a22,i3,a11,i2,a4)

    open(14,file=fname_ke,form='formatted')
    write(14,*) '"y","U","k","epsilon","nut","tau_mu","tau_R","Pk","Tk","Dk","Peps","Teps","Deps","epsEps","Pik","Pieps","epsk", &
        "Uk","Ueps","eps_hat"'
    do j=1,ny
       write(14,101) y(j),',',U(j),',',Kt(j),',',eps(j),',',nut(j),',',tau_mu(j),',',tau_R(j),',',Pk(j),',',Tk(j),',',Dk(j),',', &
       Peps(j),',',Teps(j),',',Deps(j),',',epsEps(j),',',Pi_k(j),',',Pi_eps(j),',',-(eps(j)+eps_hat(j)),',',upkt(j),',',upeps(j), &
       ',',eps_hat(j)
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
    A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10)

    102 format(e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10)

    103 format(e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10)

    end

!!!*************************************************
!!!*						         	             *
!!!*            initialization                        *
!!!*								                 *
!!!*************************************************

subroutine initialization(flag,ny,Re_tau,Pr,Bp,Cp,dy_min,y,detady,d2etady2,d3etady3,U,Kt,eps,T,Th2,deta)
    implicit none
    logical, intent(in) :: flag
    integer, intent(in) :: ny
    real*8, intent(in) :: Re_tau,dy_min,Pr,Bp,Cp
    real*8, intent(out) :: y(1:ny),detady(1:ny),d2etady2(1:ny),d3etady3(1:ny),U(1:ny)
    real*8, intent(out) :: kt(1:ny),T(1:ny),Th2(1:ny),eps(1:ny),deta
    integer j,m
    real*8 Kappa, Cmu,Ce1,Ce2,nut(1:ny),f2(1:ny),y_mid,dUdeta(1:ny),dUdy(1:ny),uk(1:ny),dukdeta(1:ny)
    real*8 U_ti(1:ny),Kt_ti(1:ny),eps_ti(1:ny)
    real*8 y_min, y_max, a1, b1, Bco, c1, c2

    Kappa = 4.d-1
    Bco = 5.5d0
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

        call grid(ny,dy_min,Re_tau,y,detady,d2etady2,d3etady3,deta)

        do j=1,ny/2
            if(y(j)<=y_mid) then
                eps(j) = (y(j)*Kappa -y(j)*y(j)*Kappa/Re_tau)/(Kappa*y_mid)**2.d0 ! 0.1335d0 / ( 1.d0 + ( ( y(j) - 15.515d0 )**2.d0 ) / 166.7634d0 ) !
                Kt(j) = dsqrt( (y(j)*Kappa -y(j)*y(j)*Kappa/Re_tau)*eps(j)/Cmu ) ! 0.019678d0 * y(j) * y(j) / ( 1.d0 + ( ( y(j) - 7.28d0 )**2.d0 ) / 88.263d0 ) !
            else
                eps(j) = (y(j)*Kappa -y(j)*y(j)*Kappa/Re_tau)/(Kappa*y(j))**2.d0 !
                Kt(j) = dsqrt( (y(j)*Kappa -y(j)*y(j)*Kappa/Re_tau)*eps(j)/Cmu ) ! 0.019678d0 * y(j) * y(j) / ( 1.d0 + ( ( y(j) - 7.28d0 )**2.d0 ) / 88.263d0 ) !
            endif
            eps(j) = (3.3234d-03 * y(j)**2.d0)/(1.d0 + 3.371626d-02 * y(j) + 4.99172d-03 * y(j)**2.d0)**2.d0 !!!0.1d0*((y(j)/12.5d0)**2.d0)*dexp(2.d0*(1.d0 - (y(j)/12.5d0)))
            Kt(j) = (8.13345d-02 * y(j)**2.d0)/(1.d0 + 8.25131877d-03 * y(j) + 5.46599692d-03 * y(j)**2.d0 &
            -5.35853011d-05 * y(j)**3.d0 + 2.27911529d-07 * y(j)**4.d0)**2.d0 !!!3.6d0*((y(j)/14d0)**2.d0)*dexp(2.d0*(1.d0 - (y(j)/14d0)))
            eps(j) = (0.00163893d0 * y(j)**2.d0)/(1.d0 + -0.03708228d0 * y(j) + 0.00672567d0 * y(j)**2.d0 &
            + ( (2.d0/(y(ny)/2.d0)**2.d0 - 0.00672567d0)/(y(ny)/2.d0) ) * y(j)**3.d0 &
             + ( ( 0.00672567d0 - 3.d0/(y(ny)/2.d0)**2.d0 )/(3.d0*(y(ny)/2.d0)**2.d0) ) * y(j)**4.d0 )**2.d0
            Kt(j) = (0.07866248d0 * y(j)**2.d0)/(1.d0 + 0.01064684d0 * y(j) + 0.00477908d0 * y(j)**2.d0 &
            + ( (2.d0/(y(ny)/2.d0)**2.d0 - 0.00477908d0)/(y(ny)/2.d0) ) * y(j)**3.d0 &
             + ( ( 0.00477908d0 - 3.d0/(y(ny)/2.d0)**2.d0 )/(3.d0*(y(ny)/2.d0)**2.d0) ) * y(j)**4.d0 )**2.d0
        enddo

        !m = 5
        !call parabolic_system(ny,Kt_ti,Kt,detady,d2etady2,deta,ny/2+1-m)
        !call parabolic_system(ny,eps_ti,eps,detady,d2etady2,deta,ny/2+1-m)

        !do j=ny/2+1-m,ny/2
        !    Kt(j) = Kt_ti(j) 
        !    eps(j) = eps_ti(j)
        !enddo

        y_min = 7.d0
        y_max = 20.d0
        c1 = -(1.d0 - 1.d0 / (Kappa*y_max) ) / (y_max - y_min)
        c2 = ( (1.d0 / Kappa) * dlog(y_max) + Bco - (y_max - y_min) - y_min ) / (y_max - y_min)**2.d0
        a1 = ( c1 - 2* c2 )/ (y_max - y_min)
        b1 = - c1 + 3* c2
        do j=1,ny/2
            if(y(j)<=y_min) then
                U(j) = y(j)
            elseif(y(j)>=y_min .and. y(j) < y_max) then
                U(j) = a1 * (y(j) - y_min)**3.d0 + b1 * (y(j) - y_min)**2.d0 + (y(j) - y_min) +y_min
            else
                U(j) = (1.d0/Kappa)*dlog(y(j)) + Bco
            endif
        enddo
        !do j=1,ny/2
        !    if(y(j)<=y_mid) then
        !        U(j) = y(j)
        !    else
        !        U(j) = (1.d0/Kappa)*dlog(y(j)) + ( y_mid - (1.d0/Kappa)*dlog(y_mid) )!+ Bco
        !    endif
        !enddo
        !Kt = 0.01d0*y**2.d0
        !eps = Kt

        do j=ny/2+1,ny
            U(j) = U(ny+1-j)
            eps(j) = eps(ny+1-j)
            Kt(j) = Kt(ny+1-j)
        enddo

        if (mod(ny,2).eq.0) then
            print*, "is even"
        else
            eps(ny/2+1) = ( eps(ny/2) + eps(ny/2+2) )/2.d0
            Kt(ny/2+1) = ( Kt(ny/2) + Kt(ny/2+2) )/2.d0
            U(ny/2+1) = ( U(ny/2) + U(ny/2+2) )/2.d0
        endif

        U_ti = U
        Kt_ti = Kt
        eps_ti = eps

        m = 3
        !do j=ny/2+1-m,ny/2+1+m
        !    U_ti(j) = ( U(j-1) + 1.d0*U(j) + U(j+1) )/3.d0
        !    Kt_ti(j) = ( Kt(j-1) + 1.d0*Kt(j) + Kt(j+1) )/3.d0
        !    eps_ti(j) = ( eps(j-1) + 1.d0*eps(j) + eps(j+1) )/3.d0
        !enddo

        !U = U_ti
        !Kt = Kt_ti
        !eps = eps_ti

        !!! Hwang and Lin correction
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

        print*, " ny/2 + 1 = ",ny/2 + 1

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

    !DA(1) = (-25.d0*a(1) +48.d0*a(2) -36.d0*a(3) +16.d0*a(4) -3.d0*a(5))/(12.d0*deta)
    !Da(2) = (-3.d0*a(1) -10.d0*a(2) +18.d0*a(3) -6.d0*a(4) +a(5))/(12.d0*deta)
    !do j=3,ny-2
    !    DA(j)= (-a(j+2) + 8.d0*a(j+1) - 8.d0*a(j-1) +a(j-2))/(12.d0*deta)
    !enddo
    !Da(ny-1) = (-a(ny-4) +6.d0*a(ny-3) -18.d0*a(ny-2) +10.d0*a(ny-1) +3.d0*a(ny))/(12.d0*deta)
    !DA(ny)= (25.d0*a(ny) -48.d0*a(ny-1) +36.d0*a(ny-2) -16.d0*a(ny-3) +3.d0*a(ny-4))/(12.d0*deta)

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
    
    D2A(1) =  (2.d0*a(1) -5.d0*a(2) +4.d0*a(3) -a(4))/deta2
    do j=2,ny-1
        D2A(j) = (a(j+1) - 2.d0*a(j) + a(j-1))/deta2
    enddo
    D2A(ny) = (2.d0*a(ny) -5.d0*a(ny-1) +4.d0*a(ny-2) -a(ny-3))/deta2

    !D2A(1) = (45.d0*a(1) -154.d0*a(2) +214.d0*a(3) -156.d0*a(4) + 61.d0*a(5) - 10.d0*a(6))/(12.d0*deta2)
    !D2A(2) = (10.d0*a(1) -15.d0*a(2) -4.d0*a(3) +14.d0*a(4) - 6.d0*a(5) + 1.d0*a(6))/(12.d0*deta2)
    !do j=3,ny-2
    !    D2A(j)= (-a(j+2) +16.d0*a(j+1) - 30.d0*a(j) + 16.d0*a(j-1) -a(j-2))/(12.d0*deta2)
    !enddo
    !D2A(ny-1) = (10.d0*a(ny) -15.d0*a(ny-1) -4.d0*a(ny-2) +14.d0*a(ny-3) -6.d0*a(ny-4) + 1.d0*a(ny-5))/(12.d0*deta2)
    !D2A(ny) = (45.d0*a(ny) -154.d0*a(ny-1) +214.d0*a(ny-2) -156.d0*a(ny-3) + 61.d0*a(ny-4) - 10.d0*a(ny-5))/(12.d0*deta2)
    
    end

!!!*************************************************
!!!*						         	             *
!!!*                    d3deta3                        *
!!!*								                 *
!!!*************************************************

subroutine  d3deta3(ny,A,D3A,deta)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: A(1:ny),deta
    real*8, intent(out) :: D3A(1:ny)
    real*8 deta3
    integer j

    deta3 = deta**3.d0
    
    D3A(1) = (-3.d0*a(5) + 14.d0*a(4) - 24.d0*a(3) + 18.d0*a(2) - 5.d0*a(1))/(2.d0*deta3)
    D3A(2) = (-a(5) + 6.d0*a(4) - 12.d0*a(3) + 10.d0*a(2) - 3.d0*a(1))/(2.d0*deta3)
    do j=3,ny-2
        D3A(j) = (a(j+2) - 2.d0*a(j+1) + 2.d0*a(j-1) - a(j-2))/(2.d0*deta3)
    enddo
    D3A(ny-1) = (a(ny-4) - 6.d0*a(ny-3) + 12.d0*a(ny-2) - 10.d0*a(ny-1) + 3.d0*a(ny))/(2.d0*deta3)
    D3A(ny) = (3.d0*a(ny-4) - 14.d0*a(ny-3) + 24.d0*a(ny-2) - 18.d0*a(ny-1) + 5.d0*a(ny))/(2.d0*deta3)

    !D3A(1) = (-49.d0*a(1) +232.d0*a(2) -461.d0*a(3) +496.d0*a(4) - 307.d0*a(5) + 104.d0*a(6) - 15.d0*a(7))/(8.d0*deta3)
    !D3A(2) = (-15.d0*a(1) +56.d0*a(2) -83.d0*a(3) +64.d0*a(4) - 29.d0*a(5) + 8.d0*a(6) - a(7))/(8.d0*deta3)
    !D3A(3) = (-a(1) - 8.d0*a(2) + 35.d0*a(3) -48.d0*a(4) + 29.d0*a(5) - 8.d0*a(6) + a(7))/(8.d0*deta3)
    !do j=4,ny-3
    !    D3A(j) = (-a(j+3) + 8.d0*a(j+2) - 13.d0*a(j+1) + 13.d0*a(j-1) - 8.d0*a(j-2) + a(j-3))/(8.d0*deta3)
    !enddo
    !D3A(ny-2) = (a(ny) + 8.d0*a(ny-1) - 35.d0*a(ny-2) + 48.d0*a(ny-3) - 29.d0*a(ny-4) + 8.d0*a(ny-5) - a(ny-6))/(8.d0*deta3)
    !D3A(ny-1) = (15.d0*a(ny) -56.d0*a(ny-1) +83.d0*a(ny-2) -64.d0*a(ny-3) + 29.d0*a(ny-4) - 8.d0*a(ny-5) + a(ny-6))/(8.d0*deta3)
    !D3A(ny) = (49.d0*a(ny) -232.d0*a(ny-1) +461.d0*a(ny-2) -496.d0*a(ny-3) + 307.d0*a(ny-4) - 104.d0*a(ny-5) + 15.d0*a(ny-6)) &
    !    /(8.d0*deta3)
    
    end

!!!*************************************************
!!!*						         	             *
!!!*                    ddy                        *
!!!*								                 *
!!!*************************************************

subroutine  ddy(ny,A,DADY,detady,deta)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: A(1:ny),detady(1:ny),deta
    real*8, intent(out) :: DADY(1:ny)
    real*8 DADeta(1:ny)
    
    call ddeta(ny,A,DADeta,deta)

    DADY = DADeta * detady
    
    end

!!!*************************************************
!!!*						         	             *
!!!*                    d2dy2                        *
!!!*								                 *
!!!*************************************************

subroutine  d2dy2(ny,A,D2ADY2,detady,d2etady2,deta)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: A(1:ny),detady(1:ny),d2etady2(1:ny),deta
    real*8, intent(out) :: D2ADY2(1:ny)
    real*8 DADeta(1:ny),D2ADeta2(1:ny)
    
    call ddeta(ny,A,DADeta,deta)
    call d2deta2(ny,A,D2ADeta2,deta)

    D2ADY2 = D2ADeta2 * (detady)**2.d0 + DADeta * d2etady2
    
    end

!!!*************************************************
!!!*						         	             *
!!!*                    d3dy3                        *
!!!*								                 *
!!!*************************************************

subroutine  d3dy3(ny,A,D3ADY3,detady,d2etady2,d3etady3,deta)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: A(1:ny),detady(1:ny),d2etady2(1:ny),d3etady3(1:ny),deta
    real*8, intent(out) :: D3ADY3(1:ny)
    real*8 DADeta(1:ny),D2ADeta2(1:ny),D3ADeta3(1:ny)
    
    call ddeta(ny,A,DADeta,deta)
    call d2deta2(ny,A,D2ADeta2,deta)
    call d3deta3(ny,A,D3ADeta3,deta)

    D3ADY3 = D3ADeta3 * (detady)**3.d0 + 3.d0 * D2ADeta2 * detady * d2etady2 + DADeta * d3etady3
    
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
    
subroutine  hwang_lin_k_epsilon_functions(nut,sigmak,sigmae,ny,y,kt,eps,detady,d2etady2,deta,Cmu)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: y(1:ny),Kt(1:ny),eps(1:ny),detady(1:ny),d2etady2(1:ny),Cmu,deta
    real*8, intent(out) :: nut(1:ny),sigmak(1:ny),sigmae(1:ny)
    real*8 y_lambda(1:ny), fmu(1:ny)
    real*8 Kt_min,eps_min
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

    end

!!!***************************************************
!!!*						         	               *
!!!*                   Epsilon Hat	       	   *
!!!*								                   *
!!!***************************************************
    
subroutine  define_eps_hat(eps_hat,kt,detady,d2etady2,deta,ny)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: Kt(1:ny),detady(1:ny),d2etady2(1:ny),deta
    real*8, intent(out) :: eps_hat(1:ny)
    real*8 uk(1:ny), dukdeta(1:ny),dukdy(1:ny), dKtdy(1:ny), d2Ktdy2(1:ny)
    integer j

    do j=1,ny
        uk(j)= dsqrt(dabs(Kt(j)))
    enddo

    call ddy(ny,uk,dukdy,detady,deta)

    call ddy(ny,Kt,dKtdy,detady,deta)
    call d2dy2(ny,Kt,d2Ktdy2,detady,d2etady2,deta)

    eps_hat(1) = d2Ktdy2(1) ! 2.d0 * (dukdy(1))**2.d0 !
    do j=2,ny-1
        if(Kt(j) /= 0.d0) then
            eps_hat(j) = (1.d0/(2.d0*Kt(j)))*(dKtdy(j))**2.d0
        else
            eps_hat(j) = d2Ktdy2(j)
        endif
    enddo
    eps_hat(ny) = d2Ktdy2(ny) ! 2.d0 * (dukdy(ny))**2.d0 !

    eps_hat = 2.d0 * (dukdy)**2.d0

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

    call ddy(ny,nut,dnutdy,detady,deta)
    call ddy(ny,U,dUdy,detady,deta)
    call ddy(ny,sigmak,dsigmakdy,detady,deta)
    call ddy(ny,sigmae,dsigmaedy,detady,deta)
    call ddy(ny,T,dTdy,detady,deta)
    call ddy(ny,Th2,dTh2dy,detady,deta)
    call d2dy2(ny,T,d2Tdy2,detady,d2etady2,deta)
    call d2dy2(ny,Th2,d2Th2dy2,detady,d2etady2,deta)

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
subroutine  K_coefficients(aK_w,aK_e,sK,eps,nut,dnutdy,nup,dnupdy,dUdy,deta,sigmak,dsigmakdy,Uw,Up,Ue,eps_hat,d2etady2,detady, &
Pi_k,j)
    implicit none
    real*8, intent(in) :: eps,nut,dnutdy,nup,dnupdy,dUdy,deta,sigmak,d2etady2,detady,dsigmakdy,Uw,Up,Ue,eps_hat,Pi_k
    real*8, intent(out) :: aK_w,aK_e,sK
    integer, intent(in) :: j
    real*8 diff, conv, b_w, b_e, b_p, D_w, D_e, F_w, F_e, F_p, DF

    diff = (1.d0+nut/sigmak+nup)*(1.d0/deta**2.d0)*detady**2.d0
    conv = (1.d0/(2.d0*sigmak))*( (sigmak+nut+sigmak*nup)*d2etady2 + detady*(dnutdy -(nut/sigmak)*dsigmakdy +dnupdy) )/deta
    D_w = diff - conv
    D_e = diff + conv
    !F_w = Uw*detady/deta
    !F_e = Ue*detady/deta
    !F_p = Up*detady/deta

    !b_w = D_w + F_w/2.d0
    !b_e = D_e - F_e/2.d0 
    !b_p = D_w + D_e

    F_w = (Uw+Up)*detady/(2.d0*deta)
    F_e = (Ue+Up)*detady/(2.d0*deta)

    b_w = D_w + F_w/2.d0
    b_e = D_e - F_e/2.d0 
    b_p = b_w + b_e + F_e - F_w
    !if(j >= 235 .and. j <=255) then
    !    print*, "b_w = ",b_w,"; b_e = ",b_e,"; b_p = ",b_p,"; j=",j
    !endif

    !if ( (b_e < 0.d0) .and. (b_w >= 0.d0) ) then
    !    !print*, '(b_e < 0.d0) .and. (b_w >= 0.d0)'
    !    b_e = D_e
    !    b_w = D_w + F_w
    !    b_p = D_e + D_w + F_p
    !else if( (b_w < 0.d0) .and. (b_e >= 0.d0) ) then
    !    !print*, '(b_w < 0.d0) .and. (b_e >= 0.d0)'
    !    b_w = D_w
    !    b_e = D_e - F_e
    !    b_p = D_w + D_e - F_p
    !else if( (b_w < 0.d0) .and. (b_e < 0.d0) ) then
    !    b_w = D_w + F_p/2.d0
    !    b_e = D_e - F_p/2.d0  
    !    b_p = D_w + D_e +(F_e - F_w)/2.d0
    !endif

    !b_w = max(max(D_w + F_w/2.d0,F_w),0.d0)
    !b_e = max(max(D_e - F_e/2.d0,-F_e),0.d0)
    !b_w = D_w * max(0.d0,(1.d0 - 0.1d0 * dabs(F_w/D_w))**5.d0) + max(F_w,0.d0)
    !b_e = D_e * max(0.d0,(1.d0 - 0.1d0 * dabs(F_e/D_e))**5.d0) + max(-F_e,0.d0)
    b_w = D_w + max(F_w,0.d0)
    b_e = D_e + max(-F_e,0.d0) 
    !b_p = b_w + b_e + F_e - F_w
    !b_w = D_w + max(F_w,0.d0)
    !b_e = D_e + max(-F_e,0.d0) 
    b_p = b_w + b_e + F_e - F_w
    !if (max(F_w,0.d0) == F_w) then
    !    DF = F_p - F_w
    !elseif ( max(-F_e,0.d0) == -F_e) then
    !    DF = F_e - F_p
    !else
    !    DF = 0.d0
    !endif
    !DF = 0.d0
    !print*, "DF = ",DF
        
    !b_p = b_w + b_e + DF !+ max(F_p - F_w,F_e - F_p)
    !if (b_p < 0.d0) then
    !    print*, "b_p = ",b_p
    !endif
    !if(F_w > 0.d0 .and. F_e > 0.d0) then
    !    b_w = D_w + F_w
    !    b_e = D_e
    !    b_p = b_w + b_e + F_p - F_w
    !elseif(F_w < 0.d0 .and. F_e < 0.d0) then
    !    b_w = D_w
    !    b_e = D_e - F_e
    !    b_p = b_w + b_e + F_e - F_p
    !else
    !    b_w = D_w + F_w/2.d0
    !    b_e = D_e - F_e/2.d0 
    !    b_p = D_w + D_e
    !endif

    !if(b_p+F_e-F_w < 0.d0 .or. b_p+F_p-F_w < 0.d0 .or. b_p+F_e-F_p < 0.d0) then
    !    print*, "b_w + b_e = ",b_w+b_e,"; F_p-F_w = ",F_p-F_w,"; F_e-F_p = ",F_e-F_p,"; F_e-F_w = ",F_e-F_w
    !endif

    !if(b_e < 0.d0 .or. b_w < 0.d0) then
    !    print*, "b_w = ",b_w,"; b_e = ",b_e
    !endif

    aK_w = b_w / b_p
    aK_e = b_e / b_p
    sK = (nut*dUdy*dUdy + Pi_k - eps_hat - eps) / b_p

    end

!!!***************************************************
!!!*						         	               *
!!!*                K coefficients	       	   *
!!!*								                   *
!!!***************************************************
subroutine  K_coefficients_a(aK_w,aK_e,sK,eps,nut,nup,dnupdy,dUdy,deta,eps_hat,d2etady2,detady,sa)
    implicit none
    real*8, intent(in) :: eps,nut,nup,dnupdy,dUdy,deta,d2etady2,detady,eps_hat,sa
    real*8, intent(out) :: aK_w,aK_e,sK
    real*8 diff, conv, b_w, b_e, b_p

    diff = (nup/deta**2.d0)*detady**2.d0
    conv = ( nup*d2etady2 + detady*dnupdy )/(2.d0*deta)
    b_w = diff - conv
    b_e = diff + conv
    b_p = b_w + b_e

    aK_w = b_w / b_p
    aK_e = b_e / b_p
    sK = (nut*dUdy*dUdy + sa - eps_hat - eps) / b_p

    end

!!!***************************************************
!!!*						         	               *
!!!*                E coefficients	       	   *
!!!*								                   *
!!!***************************************************
subroutine  E_coefficients(aE_w,aE_e,sE,eps,Kt,nut,dnutdy,dUdy,deta,sigmae,dsigmaedy,Uw,Up,Ue,Ce1,f1,Ce2,f2,d2etady2,detady)
    implicit none
    real*8, intent(in) :: eps,Kt,nut,dnutdy,dUdy,deta,sigmae,Ce1,f1,Ce2,f2,d2etady2,detady,dsigmaedy,Uw,Up,Ue
    real*8, intent(out) :: aE_w,aE_e,sE
    real*8 K_min, eps_b, diff, conv, b_w, b_e, b_p, D_w, D_e, F_w, F_e, F_p
    logical method1

    diff = (1.d0+nut/sigmae)*(1.d0/deta**2.d0)*detady**2.d0
    conv = (1.d0/(2.d0*sigmae))*( (sigmae+nut)*d2etady2 + detady*(dnutdy -(nut/sigmae)*dsigmaedy) )/deta
    D_w = diff - conv
    D_e = diff + conv
    F_w = Uw*detady/deta
    F_e = Ue*detady/deta
    F_p = Up*detady/deta

    b_w = D_w + F_w/2.d0
    b_e = D_e - F_e/2.d0  
    b_p = D_w + D_e

    !if ( (b_e < 0.d0) .and. (b_w >= 0.d0) ) then
    !    !print*, '(b_e < 0.d0) .and. (b_w >= 0.d0)'
    !    b_e = D_e
    !    b_w = D_w + F_w
    !    b_p = D_e + D_w + F_p
    !else if( (b_w < 0.d0) .and. (b_e >= 0.d0) ) then
    !    !print*, '(b_w < 0.d0) .and. (b_e >= 0.d0)'
    !    b_w = D_w
    !    b_e = D_e - F_e
    !    b_p = D_w + D_e - F_p
    !else if( (b_w < 0.d0) .and. (b_e < 0.d0) ) then
    !    b_w = D_w + F_p/2.d0
    !    b_e = D_e - F_p/2.d0  
    !    b_p = D_w + D_e +(F_e - F_w)/2.d0
    !endif

    !b_w = max(max(D_w + F_w/2.d0,F_w),0.d0)
    !b_e = max(max(D_e - F_e/2.d0,-F_e),0.d0) 
    !b_p = b_w + b_e + (F_e - F_w)/2.d0

    eps_b = (Ce1*f1*nut*dUdy*dUdy -Ce2*f2*eps)

    K_min = 1.d-60
    method1 = .true.
    if (method1) then
        aE_w = b_w / b_p
        aE_e = b_e / b_p
        if (dabs(Kt)<=K_min) then
            sE = eps_b*eps / ( b_p * K_min )
        else
            sE = eps_b*eps / ( b_p * Kt )
        endif
    else
        aE_w = ( b_w / ( b_p - eps_b / Kt) ) 
        aE_e = ( b_w / ( b_p - eps_b / Kt) ) 
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

subroutine  output_fields_k_eps(ny,U,Kt,eps,nut,f1,f2,deta,sigmaK,sigmaE,Ce1,Ce2,tau_mu,tau_R,Pk,Pi_k,Tk,Dk,Peps,Pi_eps,Teps, &
    Deps,epseps, eps_hat, y, dsigmakdy, dsigmaedy, detady, d2etady2, d3etady3)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: deta,sigmaK(1:ny),sigmaE(1:ny),Ce1,Ce2,f1,f2,Kt(1:ny),eps(1:ny),nut(1:ny),U(1:ny)
    real*8, intent(in) :: detady(1:ny),d2etady2(1:ny),d3etady3(1:ny),dsigmakdy(1:ny),dsigmaedy(1:ny), eps_hat(1:ny), y(1:ny)
    real*8, INTENT(OUT) :: tau_mu(1:ny),tau_R(1:ny),Pk(1:ny),Tk(1:ny),Dk(1:ny)
    real*8, INTENT(OUT) :: Peps(1:ny),Teps(1:ny),Deps(1:ny),epseps(1:ny),Pi_K(1:ny),Pi_eps(1:ny)
    real*8 dUdy(1:ny),dnutdy(1:ny),dKtdy(1:ny),depsdy(1:ny),Kt_min,eps_min
    integer j

    call ddy(ny,U,dUdy,detady,deta)
    call ddy(ny,nut,dnutdy,detady,deta)

    tau_mu = dUdy
    tau_R = nut*dUdy
    Pk = tau_R*dUdy

    call d2dy2(ny,Kt,Dk,detady,d2etady2,deta)
    call ddy(ny,Kt,dKtdy,detady,deta)
    Tk = (nut/sigmaK)*Dk + (dKtdy/sigmaK)*dnutdy - (nut/sigmaK**2.d0) * dKtdy * dsigmakdy

    Peps = f1*Ce1*(eps/Kt)*Pk
    Peps(1) = Peps(2)
    Peps(ny) = Peps(ny-1)

    call d2dy2(ny,eps,Deps,detady,d2etady2,deta)
    call ddy(ny,eps,depsdy,detady,deta)
    Teps = (nut/sigmaE)*Deps + (depsdy/sigmaE)*dnutdy - (nut/sigmaE**2.d0) * depsdy * dsigmaedy

    Kt_min = 1.d-60
    do j=1,ny
        if(dabs(Kt(j))<=Kt_min) then
            epsEps(j) = -f2*Ce2*(eps(j)/Kt_min)*eps(j)
        else
            epsEps(j) = -f2*Ce2*(eps(j)/Kt(j))*eps(j)
        endif
    enddo

    call calculate_Pi_K(Pi_k,Kt,eps,eps_hat,deta,detady,d2etady2,d3etady3,ny)

    Pi_eps = - (eps / Kt) * Dk - ( depsdy - eps * dKtdy / Kt ) * dKtdy / Kt

    !!! Wall approximation
    Pi_eps(1) = (Pi_eps(2) * y(3) - Pi_eps(3) * y(2)) / (y(3) - y(2))
    Pi_eps(ny) = Pi_eps(ny-1) + (Pi_eps(ny-1) - Pi_eps(ny-2))*(y(ny)-y(ny-1))/(y(ny-1)-y(ny-2))
    !print*, 'Pi_eps(1) = ',Pi_eps(1), ' Pi_eps(2) = ',Pi_eps(2),' Pi_eps(ny-1) = ',Pi_eps(ny-1),' Pi_eps(ny) = ',Pi_eps(ny)
    
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
    real*8 dTdy(1:ny),d2Tdy2(1:ny),dTh2dy(1:ny),d2Th2dy2(1:ny),dnutdy(1:ny)

    call ddy(ny,T,dTdy,detady,deta)
    call ddy(ny,Th2,dTh2dy,detady,deta)
    call ddy(ny,nut,dnutdy,detady,deta)
    call d2dy2(ny,T,d2Tdy2,detady,d2etady2,deta)
    call d2dy2(ny,Th2,d2Th2dy2,detady,d2etady2,deta)

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

subroutine grid(ny,dy_min,Re_tau,y,detady,d2etady2,d3etady3,deta)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: Re_tau,dy_min
    real*8, intent(out) :: y(1:ny),detady(1:ny),d2etady2(1:ny),d3etady3(1:ny),deta
    integer j
    real*8 a, b, c, d, e, eta, y_max

    y_max = 2.d0*Re_tau

    deta = y_max/(ny-1)

    a = (deta - dy_min)/(Re_tau*deta - deta*deta)
    b = (dy_min*Re_tau - deta*deta)/(Re_tau*deta - deta*deta)

    do j=1,ny/2
        eta = deta*(j-1)
        y(j) = a*eta**2.d0 + b*eta
        detady(j) = 1.d0/(2.d0*a*eta + b)
        d2etady2(j) = -2.d0*a/(2.d0*a*eta + b)**3.d0
        d3etady3(j) = ( 12.d0*a**2.d0 ) / (2.d0*a*eta +b)**5.d0
    enddo

    c = -a
    d = 4.d0*Re_tau*a + b
    e = 2.d0*Re_tau*( 1.d0 - b -2.d0*a*Re_tau )

    do j=ny/2+1,ny
        eta = deta*(j-1)
        y(j) = c*eta**2.d0 + d*eta + e
        detady(j) = 1.d0/(2.d0*c*eta + d)
        d2etady2(j) = -2.d0*c/(2.d0*c*eta + d)**3.d0
        d3etady3(j) = ( 12.d0*c**2.d0 ) / (2.d0*c*eta +d)**5.d0
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

    increment = 1.d-01

    if(res < res_old) then
        alpha = alpha * ( 1.d0 + increment )
    elseif(res > res_old) then
        alpha = alpha * ( 1.d0 - increment )
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
    
subroutine  solve_Kt(Kt,eps,dUdy,nut,dnutdy,detady,d2etady2,d3etady3,deta,sigmak,dsigmakdy,eps_hat,ny,upk,y)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: eps(1:ny),dUdy(1:ny),nut(1:ny),dnutdy(1:ny),detady(1:ny),d2etady2(1:ny),deta,sigmak(1:ny)
    real*8, intent(in) :: dsigmakdy(1:ny),eps_hat(1:ny),d3etady3(1:ny),y(1:ny),upk(1:ny)
    real*8, intent(inout) :: Kt(1:ny)
    !real*8, intent(out) :: upk(1:ny)
    real*8 aK_w,aK_e,sK
    real*8 dupkdy(1:ny),nup(1:ny),dnupdy(1:ny),depsdy(1:ny),D3KtDY3(1:ny),D2KtDY2(1:ny),dKtdy(1:ny),Pi_k(1:ny)
    real*8 uk(1:ny),dukdy(1:ny),d2ukdy2(1:ny),dlnEdy(1:ny),Kt_star(1:ny),deps_hatdy(1:ny),nup_p(1:ny),sa(1:ny)
    integer j
    real*8 A(1:ny),C_apex(1:ny),denominator

    do j=1,ny
        uk(j)= dsqrt(dabs(Kt(j)))
    enddo

    call ddy(ny,uk,dukdy,detady,deta)
    call d2dy2(ny,uk,d2ukdy2,detady,d2etady2,deta)
    call d2dy2(ny,Kt,d2Ktdy2,detady,d2etady2,deta)
    call d3dy3(ny,Kt,D3KtDY3,detady,d2etady2,d3etady3,deta)
    call ddy(ny,Kt,dKtdy,detady,deta)
    call ddy(ny,eps,depsdy,detady,deta)
    call ddy(ny,eps_hat,deps_hatdy,detady,deta)

    !nup = - uk * dukdy / ( eps + eps_hat )
    nup = 0*( eps_hat - D2KtDY2 ) / ( 2.d0 * ( eps + eps_hat ) )
    !!! for a
    nup_p = 1.d0 + nut / sigmak + eps_hat / ( 2.d0 * ( eps + eps_hat ) ) 
    nup = nup_p + ( deps_hatdy + depsdy ) * dKtdy / ( 2.d0 * ( eps + eps_hat )**2.d0 )

    call ddy(ny,dlog(dabs(eps_hat+eps)),dlnEdy,detady,deta)
    call ddy(ny,nup_p,dnupdy,detady,deta)
    !sa =  - D3KtDY3 * dKtdy / ( 2.d0 * ( eps + eps_hat ) )
    call d2dy2(ny,eps_hat*Kt,sa,detady,d2etady2,deta)
    sa = - sa / ( 2.d0 * ( eps + eps_hat) )
    sa = sa/10.d0
    !dnupdy = 0.d0 !- nup * dlnEdy -eps_hat / ( 2.d0 * ( eps + eps_hat ) ) - uk * d2ukdy2 / ( eps + eps_hat )
    !dnupdy = dnupdy !/ 75.d0
    do j=1,16
        call quadratic_variable_smoother(sa,sa,y,ny)
    enddo

    !nup = ( eps_hat - d2Ktdy2 ) / ( 2.d0 * ( eps + eps_hat ) )
    !nup = 0.d0
    
    !dnupdy = - nup * ( depsdy / ( eps + eps_hat ) + 0*(1.d0 - 2.d0 *nup ) * dKtdy / Kt ) &
    !- 0*d3Ktdy3 / ( 2.d0 * ( eps + eps_hat ) )
    !print*, 'nup(1) = ',nup(1),'; nup(ny) = ',nup(ny)
    !print*, 'nup(2) = ',nup(2),'; nup(ny-1) = ',nup(ny-1),'; dnupdy(2) = ',dnupdy(2),'; dnupdy(ny-1) = ',dnupdy(ny-1)

    !dnupdy(1) = 0.d0
    !dnupdy(ny) = 0.d0
    !call linear_variable_smoother(dnupdy,dnupdy,0.d0*dnupdy,ny)

    !!! This is the problem
    !call pressure_speed_Kt(upk,Kt,eps,eps_hat,detady,d2etady2,d3etady3,deta,y,ny)
    !call ddy(ny,eps_hat,deps_hatdy,detady,deta)
    !upk = deps_hatdy / ( 2.d0 * ( eps + eps_hat ) )
    !do j=1,16
    !    call quadratic_variable_smoother(upk,upk,y,ny)
    !enddo
    !call calculate_Pi_K(Pi_k,Kt,eps,eps_hat,deta,detady,d2etady2,d3etady3,ny)
    !do j=1,16
    !    call quadratic_variable_smoother(Pi_k,Pi_k,y,ny)
    !enddo
    !Pi_k = 0.d0
    Kt(1) = 0.d0
    do j =2,ny-1
        !call K_coefficients(aK_w,aK_e,sK,eps(j),nut(j),dnutdy(j),nup(j),dnupdy(j),dUdy(j),deta,sigmak(j),dsigmakdy(j), &
        !    upk(j-1), upk(j), upk(j+1), eps_hat(j), d2etady2(j),detady(j), Pi_k(j), j)
        call K_coefficients_a(aK_w,aK_e,sK,eps(j),nut(j),nup(j),dnupdy(j),dUdy(j),deta,eps_hat(j),d2etady2(j),detady(j),sa(j))
        if (aK_e < 0.d0) then
            print*, "aK_e=",aK_e, " j=",j
        endif
        if (aK_w < 0.d0) then
            print*, "aK_w=",aK_w, " j=",j
        endif
        !if(j >= 235 .and. j <=255) then
        !    print*, "aK_w = ",aK_w,"; aK_e = ",aK_e,"; sK = ",sK,"; j=",j
        !endif
        !Kt_star(j) = sK
        Kt(j) = sK + aK_e*Kt(j+1) + aK_w*Kt(j-1)
        !if(j > 235 .and. j < 255) then
        !    print*, ' aK_w = ',aK_w,'; aK_e = ',aK_e,'; sK = ',sK, '; Kt(j) = ',Kt(j)
        !    !print*, ' Kt(j-1) = ',Kt(j-1),'; Kt(j) = ',Kt(j),'; Kt(j+1) = ',Kt(j+1),'; j = ',j
        !endif
    enddo
    Kt(ny) = 0.d0
    !Kt(1) = 0.d0
    !call K_coefficients(aK_w,aK_e,sK,eps(2),nut(2),dnutdy(2),nup(2),dnupdy(2),dUdy(2),deta,sigmak(2),dsigmakdy(2), &
    !        upk(1),upk(2), upk(3), eps_hat(2), d2etady2(2),detady(2), Pi_k(2))
    !A(2) = aK_e
    !C_apex(2) = sK + aK_w * Kt(1)
    !do j =3,ny-1
    !    call K_coefficients(aK_w,aK_e,sK,eps(j),nut(j),dnutdy(j),nup(j),dnupdy(j),dUdy(j),deta,sigmak(j),dsigmakdy(j), &
    !        upk(j-1),upk(j), upk(j+1), eps_hat(j), d2etady2(j),detady(j), Pi_k(j))
    !    denominator = ( 1.d0 - aK_w * A(j-1) )
    !    A(j) = aK_e / denominator
    !    C_apex(j) = ( aK_w * C_apex(j-1) + sK ) / denominator
    !enddo

    !A(ny-1) = 0.d0
    !Kt(ny) = 0.d0

    !do j =ny-1,2,-1
    !    Kt(j) = A(j) * Kt(j+1) + C_apex(j)
    !enddo

    !do j=1,10
    !call quadratic_variable_smoother(Kt,Kt,y,ny)
    !enddo
    !upk = detady !- D2KtDY2 / ( 2.d0 * ( eps + eps_hat ) ) !Kt_star !nup !dKtdy * ( d2Ktdy2 - eps_hat ) / ( 2.d0 * Kt * (eps + eps_hat))

    end

!!!*************************************************
!!!*						         	           *
!!!*                 Solve Eps                        *
!!!*								               *
!!!*************************************************
    
subroutine  solve_eps(Kt,eps,dUdy,nut,dnutdy,detady,d2etady2,deta,sigmae,dsigmaedy,eps_hat,ce1,ce2,f1,f2,ny,upeps)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: Kt(1:ny),dUdy(1:ny),nut(1:ny),dnutdy(1:ny),detady(1:ny),d2etady2(1:ny),deta,sigmae(1:ny)
    real*8, intent(in) :: ce1,ce2,f1,f2,dsigmaedy(1:ny),eps_hat(1:ny)
    real*8, intent(inout) :: eps(1:ny)
    real*8, intent(out) :: upeps(1:ny)
    real*8 aE_w,aE_e,sE
    real*8 dKtdeta(1:ny),duepsdy(1:ny),d2Ktdeta2(1:ny)
    real*8 A(1:ny),C_apex(1:ny),denominator
    integer j

    call ddeta(ny,Kt,dKtdeta,deta)
    call d2deta2(ny,Kt,d2Ktdeta2,deta)

    upeps = dKtdeta * detady / Kt
    !!! wall correction
    upeps(1) = (5.d0*upeps(2) -4.d0*upeps(3) +upeps(4))/2.d0
    upeps(ny) = (5.d0*upeps(ny-1) -4.d0*upeps(ny-2) +upeps(ny-3))/2.d0
    upeps(1) = 2.d0*upeps(2) 
    upeps(ny) = 2.d0*upeps(ny-1)

    !duepsdy = (1.d0 / Kt) * ( d2Ktdeta2 * ( detady )**2.d0 + dKtdeta * d2etady2 ) &
    !- ( dKtdeta * detady / Kt ) **2.d0
    !upeps = upeps / 1.d+01
    !duepsdy = 0.d0*duepsdy

    eps(1) = 0.d0
    do j =2,ny-1
        call E_coefficients(aE_w,aE_e,sE,eps(j),Kt(j),nut(j),dnutdy(j),dUdy(j),deta,sigmae(j),dsigmaedy(j),upeps(j-1), &
            upeps(j), upeps(j+1), Ce1,f1,Ce2,f2,d2etady2(j), detady(j))
        !if(j==2 .OR. j==ny-1) then
        !    print*, eps(j),Kt(j),nut(j),dnutdy(j),dUdy(j),deta,sigmae(j),dsigmaedy(j),Pieps(j),d2etady2(j), detady(j)
        !endif
        !PRINT*, ' j = ',j,'; sE = ',sE,'; aE_w = ',aE_w,',aE_e = ',aE_e
        !if (aE_e < 0.d0) then
        !    print*, "aE_e=",aE_e, " j=",j
        !endif
        !if (aE_w < 0.d0) then
        !    print*, "aE_w=",aE_w, " j=",j
        !endif
        eps(j) = sE + aE_e*eps(j+1) + aE_w*eps(j-1)
    enddo
    eps(ny) = 0.d0 

    !eps(1) = 0.d0
    !call E_coefficients(aE_w,aE_e,sE,eps(2),Kt(2),nut(2),dnutdy(2),dUdy(2),deta,sigmae(2),dsigmaedy(2),upeps(1), &
    !        upeps(2), upeps(3), Ce1,f1,Ce2,f2,d2etady2(2), detady(2))
    !A(2) = aE_e
    !C_apex(2) = sE + aE_w * eps(1)
    !do j =3,ny-1
    !    call E_coefficients(aE_w,aE_e,sE,eps(j),Kt(j),nut(j),dnutdy(j),dUdy(j),deta,sigmae(j),dsigmaedy(j),upeps(j-1), &
    !        upeps(j), upeps(j+1), Ce1,f1,Ce2,f2,d2etady2(j), detady(j))
    !    denominator = ( 1.d0 - aE_w * A(j-1) )
    !    A(j) = aE_e / denominator
    !    C_apex(j) = ( aE_w * C_apex(j-1) + sE ) / denominator
    !enddo

    !A(ny-1) = 0.d0
    !eps(ny) = 0.d0

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

!!!*************************************************
!!!*						         	           *
!!!*           Quadratic Variable Smoother                     *
!!!*								               *
!!!*************************************************
    
subroutine  quadratic_variable_smoother(phi_hat,phi,y,ny)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: phi(1:ny),y(1:ny)
    real*8, intent(out) :: phi_hat(1:ny)
    real*8 alpha,beta, gamma, delta
    integer j

    phi_hat = phi

    do j=2,ny-1
        delta = (y(j+1) - y(j)) / (y(j) - y(j-1))
        alpha = ( 2.d0 - delta )/6.d0
        gamma = ( 2.d0 - 1.d0 / delta )/6.d0
        beta = 1.d0 - alpha - gamma
        phi_hat(j) = alpha * phi(j-1) + beta * phi(j) + gamma * phi(j+1)
    enddo
    

    end

!!!*************************************************
!!!*						         	           *
!!!*           Linear Variable Smoother                     *
!!!*								               *
!!!*************************************************
    
subroutine  linear_variable_smoother(phi_hat,phi,y,ny)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: phi(1:ny),y(1:ny)
    real*8, intent(out) :: phi_hat(1:ny)
    integer j

    phi_hat = phi

    phi_hat(1) = 2.d0*phi(2) - phi(3)
    do j=2,ny-1
        !phi_hat(j) = ( phi(j-1) + 2.d0 * phi(j) + phi(j+1) ) / 4.d0
        phi_hat(j) = ( phi(j-1) + phi(j) + phi(j+1) ) / 3.d0
    enddo
    phi_hat(ny) = 2.d0*phi(ny-1) - phi(ny-2)
    

    end

!!!*************************************************
!!!*						         	           *
!!!*           Positive part                    *
!!!*								               *
!!!*************************************************
    
subroutine  positive_part(up,u)
    implicit none
    real*8, intent(in) :: u
    real*8, intent(out) :: up
    integer j

    up = ( u + dabs(u) ) / 2.d0
    
    end

!!!*************************************************
!!!*						         	           *
!!!*           Negative part                    *
!!!*								               *
!!!*************************************************
    
subroutine  negative_part(un,u)
    implicit none
    real*8, intent(in) :: u
    real*8, intent(out) :: un
    integer j

    un = ( u - dabs(u) ) / 2.d0
    
    end

!!!*************************************************
!!!*						         	           *
!!!*                 Pi k                       *
!!!*								               *
!!!*************************************************
    
subroutine  calculate_Pi_K(Pi_k,Kt,eps,eps_hat,deta,detady,d2etady2,d3etady3,ny)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: deta,detady(1:ny),d2etady2(1:ny),d3etady3(1:ny)
    real*8, intent(in) :: Kt(1:ny),eps(1:ny),eps_hat(1:ny)
    real*8, intent(out) :: Pi_k(1:ny)
    real*8 depsdy(1:ny),d3Ktdy3(1:ny),d2Ktdy2(1:ny),dKtdy(1:ny)
    integer j

    call ddy(ny,eps,depsdy,detady,deta)
    call ddy(ny,Kt,dKtdy,detady,deta)
    call d2dy2(ny,Kt,d2Ktdy2,detady,d2etady2,deta)
    call d3dy3(ny,Kt,d3Ktdy3,detady,d2etady2,d3etady3,deta)

    Pi_k = - ( dKtdy * d3Ktdy3 + ( d2Ktdy2 - eps_hat ) * ( ( d2Ktdy2 - eps_hat ) * &
        ( ( eps - eps_hat ) / ( eps + eps_hat ) ) - eps_hat - depsdy * dKtdy / ( eps + eps_hat ) ) &
        ) / ( 2.d0 * ( eps + eps_hat) )
    end

!!!***************************************************
!!!*						         	               *
!!!*                K coefficients	       	   *
!!!*								                   *
!!!***************************************************
subroutine  K_coefficients_2(aK_ww,aK_w,aK_e,aK_ee,sK,eps,nut,dnutdy,dUdy,deta,sigmak,dsigmakdy,Uww,Uw,Up,Ue,Uee,eps_hat, &
    d2etady2,detady)
    implicit none
    real*8, intent(in) :: eps,nut,dnutdy,dUdy,deta,sigmak,d2etady2,detady,dsigmakdy,Uww,Uw,Up,Ue,Uee,eps_hat
    real*8, intent(out) :: aK_ww,aK_w,aK_e,aK_ee,sK
    real*8 b_ww, b_w, b_e, b_ee, b_p, D_w, D_e, F_ww, F_w, F_e, F_ee, F_p
    real*8 A_k,B_k,C_k,E_k,sigma_K, beta,gamma, tollerance
    real*8 alpha_e, alpha_w

    tollerance = 1.d-02
    !DA_e(j) = (-2.d0*Uw*a(j-1) -3.d0*Up*a(j) + 6.d0*Ue*a(j+1) -Uee*a(j+2))/(6.d0*deta)
    !DA_w(j) = (Uww*a(j-2) - 6.d0*Uw*a(j-1) + 3.d0*Up*a(j) +2.d0*Ue*a(j+1))/(6.d0*deta)

    A_k = (1.d0+nut/sigmak)*detady**2.d0
    B_k = ( (1.d0+nut/sigmak)*d2etady2 + detady*(dnutdy -(nut/sigmak)*dsigmakdy)/sigmak )
    
    E_k = deta * B_k / ( 4.d0 * A_k )
    sigma_K = (nut*dUdy*dUdy - eps_hat - eps) * (deta**2.d0) / ( 2.d0 * A_k)

    D_w = 5.d-01 - E_k
    D_e = 5.d-01 + E_k

    goto 7373
    C_k = deta * detady / ( 4.d0 * A_k )
    F_ww = C_k*Uww
    F_w = C_k*Uw
    F_e = C_k*Ue
    F_ee = C_k*Uee
    F_p = C_k*Up

    !Targets.
    beta = 1.d0
    gamma = -1.d0
    !Initialize
    b_w = D_w - gamma * F_w
    b_e = D_e - beta * F_e
    b_p = 1.d0 - F_p * ( beta + gamma ) * 3.d0 / 4.d0
    b_ee = - F_ee * ( 4.d0 -3.d0 * beta +  gamma ) / 8.d0
    b_ww = F_ww * ( 4.d0 +3.d0 * gamma - beta ) / 8.d0

    !if(b_e<0 .or. b_w<0 .or. b_p < 0) then
    !    print*, 'D_w / F_w = ',D_w / F_w, '; D_e / F_e = ',D_e / F_e
    !    print*, 'F_ww = ',F_ww,'; F_p = ',F_p,'; F_ee = ',F_ee
    !endif
    !gamma < = D_w / F_w
    !beta < = D_e / F_e
    !beta + gamma < = 4.d0 / ( 3.d0 * F_p )
    !- F_ee * ( 4.d0 -3.d0 * beta +  gamma ) > =0
    !F_ww * ( 4.d0 +3.d0 * gamma - beta ) > =0

    if ( (b_e < 0.d0) .and. (b_w >= 0.d0) ) then
        print*, ' b_e <0'
        !print*, 'D_e / F_e = ',D_e / F_e,'; 4.d0* ( 1.d0 -  tollerance ) / ( 3.d0 * F_p ) - gamma = ',  &
        !4.d0* ( 1.d0 -  tollerance ) / ( 3.d0 * F_p ) - gamma
        beta = D_e / F_e
        !print*, ' beta = ',min(D_e / F_e,  4.d0* ( 1.d0 -  tollerance ) / ( 3.d0 * F_p ) - gamma)
        gamma = -1.d0
        b_p = 1.d0 - F_p * ( beta + gamma ) * 3.d0 / 4.d0
        if(b_p <= 0.d0) then
            print*, ' b_e = ',b_e,'; b_p = ',b_p
            gamma = 4.d0 * ( 1.d0 -  tollerance ) / ( 3.d0 * F_p) - beta
        endif
        b_ee = - F_ee * ( 4.d0 -3.d0 * beta +  gamma ) / 8.d0
        if(b_ee <= 0.d0) then
            print*, ' b_e = ',b_e,'; b_ee = ',b_ee
            !gamma = 3.d0 * beta -4.d0
        endif
        b_ww = F_ww * ( 4.d0 +3.d0 * gamma - beta ) / 8.d0
        if(b_ww <= 0.d0) then
            print*, ' b_e = ',b_e,'; b_ww = ',b_ww
            !gamma = ( beta -4.d0 ) / 3.d0
        endif
    else if( (b_w < 0.d0) .and. (b_e >= 0.d0) ) then
    print*, ' b_w<0'
        !print*, 'D_w / F_w = ',D_w / F_w,'; 4.d0* ( 1.d0 -  tollerance ) / ( 3.d0 * F_p ) - beta = ',  &
        !4.d0* ( 1.d0 -  tollerance ) / ( 3.d0 * F_p ) - beta
        gamma = D_w /F_w
        !print*, ' gamma = ',min( D_w /F_w, 4.d0 * ( 1.d0 -  tollerance ) / ( 3.d0 * F_p ) - beta)
        beta = 1.d0
        b_p = 1.d0 - F_p * ( beta + gamma ) * 3.d0 / 4.d0
        if(b_p <= 0.d0) then
            print*, ' b_w = ',b_w,'; b_p = ',b_p
            beta = 4.d0 * ( 1.d0 -  tollerance ) / ( 3.d0 * F_p) - gamma
        endif
        b_ee = - F_ee * ( 4.d0 -3.d0 * beta +  gamma ) / 8.d0
        if(b_ee <= 0.d0) then
            print*, ' b_w = ',b_w,'; b_ee = ',b_ee
            !beta = ( 4.d0 + gamma )/ 3.d0 
        endif
        b_ww = F_ww * ( 4.d0 +3.d0 * gamma - beta ) / 8.d0
        if(b_ww <= 0.d0) then
            print*, ' b_w = ',b_w,'; b_ww = ',b_ww
            !beta = 4.d0 +3.d0 * gamma
        endif
    else if( (b_w < 0.d0) .and. (b_e < 0.d0) .and. (b_p >= 0.d0) ) then
        print*,'both <0'
        beta = D_e / F_e
        gamma = D_w /F_w
        b_p = 1.d0 - F_p * ( beta + gamma ) * 3.d0 / 4.d0
        if(b_p <= 0.d0) then
            print*, 'both <0; b_p = ',b_p
        endif
        b_ee = - F_ee * ( 4.d0 -3.d0 * beta +  gamma ) / 8.d0
        if(b_ee <= 0.d0) then
            print*, 'both <0; b_ee = ',b_ee
        endif
        b_ww = F_ww * ( 4.d0 +3.d0 * gamma - beta ) / 8.d0
        if(b_ww <= 0.d0) then
            print*, 'both <0; b_ww = ',b_ww
        endif
    else
        beta = 1.d0
        gamma = -1.d0
    endif

    b_w = D_w - gamma * F_w
    b_e = D_e - beta * F_e
    b_p = 1.d0 - F_p * ( beta + gamma ) * 3.d0 / 4.d0
    b_ee = - F_ee * ( 4.d0 -3.d0 * beta +  gamma ) / 8.d0
    b_ww = F_ww * ( 4.d0 +3.d0 * gamma - beta ) / 8.d0

    aK_ww = b_ww / b_p
    aK_w = b_w / b_p
    aK_e = b_e / b_p
    aK_ee = b_ee / b_p
    sK = sigma_K / b_p

    !goto 8383
    7373 continue
    C_k = deta * detady / ( 2.d0 * A_k )
    F_w = C_k*( Uw + Up ) / 2.d0
    F_e = C_k*( Ue + Up ) / 2.d0

    if(F_w <= 0.d0) then
        alpha_w = 0.d0
    else
        alpha_w = 1.d0
    endif

    if(F_e <= 0.d0) then
        alpha_e = 0.d0
    else
        alpha_e = 1.d0
    endif

    b_ww = - (alpha_w / 8.d0) * F_w
    b_w = D_w + ( alpha_e / 8.d0 ) * F_e + ( 3.d0 * ( 1.d0 + alpha_w ) / 8.d0 ) * F_w
    b_e = D_e - ( 3.d0 * ( 2.d0 - alpha_e) / 8.d0 ) * F_e - ( ( 1.d0 - alpha_w ) / 8.d0 ) * F_w 
    b_ee = ( (1.d0 -alpha_e)/8.d0 ) * F_e
    b_p = b_ww + b_w + b_e + b_ee + F_e - F_w

    aK_ww = b_ww / b_p
    aK_w = b_w / b_p
    aK_e = b_e / b_p
    aK_ee = b_ee / b_p
    sK = sigma_K / b_p
    8383 continue

    end

!!!*************************************************
!!!*						         	           *
!!!*                 Solve Kt                        *
!!!*								               *
!!!*************************************************
    
subroutine  solve_Kt_2(Kt,eps,dUdy,nut,dnutdy,detady,d2etady2,d3etady3,deta,sigmak,dsigmakdy,eps_hat,ny,upk,y)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: eps(1:ny),dUdy(1:ny),nut(1:ny),dnutdy(1:ny),detady(1:ny),d2etady2(1:ny),deta,sigmak(1:ny)
    real*8, intent(in) :: dsigmakdy(1:ny),eps_hat(1:ny),d3etady3(1:ny),y(1:ny), upk(1:ny)
    real*8, intent(inout) :: Kt(1:ny)
    !real*8, intent(out) :: upk(1:ny)
    real*8 aK_ww,aK_w,aK_e,aK_ee,sK
    real*8 dupkdy(1:ny),nup(1:ny),dnupdy(1:ny)
    real*8 depsdy(1:ny),D3KtDY3(1:ny),D2KtDY2(1:ny),dKtdy(1:ny)
    real*8 uk(1:ny),dukdy(1:ny),d2ukdy2(1:ny),dlnEdy(1:ny),Kt_star(1:ny),deps_hatdy(1:ny)
    integer j

    !call pressure_speed_Kt(upk,Kt,eps,eps_hat,detady,d2etady2,d3etady3,deta,y,ny)
    !upk = Kt * deps_hatdy / ( 2.d0 * ( eps + eps_hat ) )
    !do j=1,16
    !    call quadratic_variable_smoother(upk,upk,y,ny)
    !enddo

    Kt(1) = 0.d0
    call K_coefficients_2(aK_ww,aK_w,aK_e,aK_ee,sK,eps(2),nut(2),dnutdy(2),dUdy(2),deta,sigmak(2), &
            dsigmakdy(2), 0.d0*upk(2), upk(1),upk(2), upk(3), upk(4), eps_hat(2), d2etady2(2),detady(2))
    Kt(2) = sK + aK_ee*Kt(4) + aK_e*Kt(3) + aK_w*Kt(1)
    do j =3,ny-2
        call K_coefficients_2(aK_ww,aK_w,aK_e,aK_ee,sK,eps(j),nut(j),dnutdy(j),dUdy(j),deta,sigmak(j), &
            dsigmakdy(j), upk(j-2), upk(j-1),upk(j), upk(j+1), upk(j+2), eps_hat(j), d2etady2(j),detady(j))
        if (aK_e < 0.d0) then
            print*, "aK_e=",aK_e, " j=",j
        endif
        if (aK_w < 0.d0) then
            print*, "aK_w=",aK_w, " j=",j
        endif
        !if (aK_ww < 0.d0) then
        !    print*, "aK_ww=",aK_ww, " j=",j
        !endif
        !if (aK_ee < 0.d0) then
        !    print*, "aK_ee=",aK_ee, " j=",j
        !endif
        Kt(j) = sK + aK_ee*Kt(j+2) + aK_e*Kt(j+1) + aK_w*Kt(j-1) + aK_ww*Kt(j-2)
    enddo
    call K_coefficients_2(aK_ww,aK_w,aK_e,aK_ee,sK,eps(ny-1),nut(ny-1),dnutdy(ny-1),dUdy(ny-1),deta,sigmak(ny-1), &
            dsigmakdy(ny-1), upk(ny-3), upk(ny-2),upk(ny-1), upk(ny), 0.d0*upk(ny-1), eps_hat(ny-1), d2etady2(ny-1),detady(ny-1))
    Kt(ny-1) = sK + aK_e*Kt(ny) + aK_w*Kt(ny-2) + aK_ww*Kt(ny-3)
    Kt(ny) = 0.d0

    !Kt = Kt_star
    !upk = ( d2Ktdy2 - eps_hat ) * dKtdy / ( 2.d0 * Kt * (eps + eps_hat))

    end

!!!***************************************************
!!!*						         	               *
!!!*                K coefficients	       	   *
!!!*								                   *
!!!***************************************************
subroutine  K_coefficients_tvd(aK_w,aK_e,sK,eps,nut,dnutdy,dUdy,deta,sigmak,dsigmakdy,Uw,Up,Ue,eps_hat, &
    d2etady2,detady,psi_w_m,psi_w_p,psi_e_m,psi_e_p,Kt_w,Kt_p,Kt_e)
    implicit none
    real*8, intent(in) :: eps,nut,dnutdy,dUdy,deta,sigmak,d2etady2,detady,dsigmakdy,Uw,Up,Ue,eps_hat
    real*8, intent(in) :: psi_w_m,psi_w_p,psi_e_m,psi_e_p,Kt_w,Kt_p,Kt_e
    real*8, intent(out) :: aK_w,aK_e,sK
    real*8 b_w, b_e, b_p, D_w, D_e, F_w, F_e, F_p
    real*8 A_k,B_k,C_k,E_k,sigma_K, beta,gamma
    real*8 alpha_e, alpha_w

    A_k = (1.d0+nut/sigmak)*detady**2.d0
    B_k = ( (1.d0+nut/sigmak)*d2etady2 + detady*(dnutdy -(nut/sigmak)*dsigmakdy)/sigmak )
    
    E_k = deta * B_k / ( 4.d0 * A_k )
    sigma_K = (nut*dUdy*dUdy - eps_hat - eps) * (deta**2.d0) / ( 2.d0 * A_k)

    D_w = 5.d-01 - E_k
    D_e = 5.d-01 + E_k

    C_k = deta * detady / ( 2.d0 * A_k )
    F_w = 2.d0 * C_k*( Uw + Up ) / 2.d0
    F_e = 2.d0 * C_k*( Ue + Uw ) / 2.d0

    if(F_w <= 0.d0) then
        alpha_w = 0.d0
    else
        alpha_w = 1.d0
    endif

    if(F_e <= 0.d0) then
        alpha_e = 0.d0
    else
        alpha_e = 1.d0
    endif

    !b_w = D_w + F_w * ( alpha_w * ( 1.d0 - psi_w_p / 2.d0 ) + ( 1.d0 - alpha_w ) * psi_w_m / 2.d0 )
    !b_e = D_e - F_e * ( alpha_e * psi_e_p / 2.d0 + ( 1.d0 - alpha_e ) * ( 1.d0 - psi_e_m / 2.d0 ) )
    !b_p = b_w + b_e + F_e - F_w
    b_w = D_w + max(F_w,0.d0) 
    b_e = D_e + max(- F_e,0.d0) 
    sigma_K = sigma_K + ( F_e / 2.d0 ) * ( ( 1.d0 - alpha_e ) * psi_e_m - alpha_e * psi_e_p ) * (Kt_e - Kt_p) &
    + ( F_w / 2.d0 ) * ( ( 1.d0 - alpha_w ) * psi_w_m - alpha_w * psi_w_p ) * (Kt_w - Kt_p)
    b_p = b_w + b_e + F_e - F_w
    !if ( b_w < 0.d0 .and. b_e >= 0.d0 ) then
    !    print*, 'b_w = ',b_w
    !    b_w = 0.d0
    !    b_p = b_w + b_e + F_e - F_w
    !elseif ( b_e < 0.d0 .and. b_w >= 0.d0  ) then
    !    print*, 'b_e = ',b_e
    !    b_e = 0.d0
    !    b_p = b_w + b_e + F_e - F_w
    !elseif ( b_w < 0.d0 .and. b_e < 0.d0  ) then
    !    print*, 'b_w = ',b_w,'; b_e = ',b_e
    !    b_w = 0.d0
    !    b_e = 0.d0
    !    b_p = b_w + b_e + F_e - F_w
    !endif
    
    aK_w = b_w / b_p
    aK_e = b_e / b_p
    sK = sigma_K / b_p

    end

!!!***************************************************
!!!*						         	               *
!!!*                Limiter Function	       	   *
!!!*								                   *
!!!***************************************************
subroutine  Limiter_Function(psi_w_m, psi_w_p, psi_e_m, psi_e_p, phi_ww, phi_w, phi_p, phi_e, phi_ee, method)
    implicit none
    real*8, intent(in) :: phi_ww, phi_w, phi_p, phi_e, phi_ee
    character(len = 10), intent(in) :: method
    real*8, intent(out) :: psi_w_m,psi_w_p,psi_e_m,psi_e_p
    real*8 r_w_m, r_w_p, r_e_m, r_e_p

    r_e_p = (phi_p - phi_w)/(phi_e - phi_p)
    r_w_p = (phi_w - phi_ww)/(phi_p - phi_w)
    r_e_m = (phi_ee - phi_e)/(phi_e - phi_p)
    r_w_m = (phi_e - phi_p)/(phi_p - phi_w)

    select case (method)
    case ("VanLeer")
        call Van_Leer(psi_e_p,r_e_p)
        call Van_Leer(psi_e_m,r_e_m)
        call Van_Leer(psi_w_p,r_w_p)
        call Van_Leer(psi_w_m,r_w_m)
    case ("VanAlbada")
        call Van_Albada(psi_e_p,r_e_p)
        call Van_Albada(psi_e_m,r_e_m)
        call Van_Albada(psi_w_p,r_w_p)
        call Van_Albada(psi_w_m,r_w_m)
    case ("MinMod")
        call Min_Mod(psi_e_p,r_e_p)
        call Min_Mod(psi_e_m,r_e_m)
        call Min_Mod(psi_w_p,r_w_p)
        call Min_Mod(psi_w_m,r_w_m)
    case ("SuperBee")
        call SuperBee(psi_e_p,r_e_p)
        call SuperBee(psi_e_m,r_e_m)
        call SuperBee(psi_w_p,r_w_p)
        call SuperBee(psi_w_m,r_w_m)
    case ("Quick")
        call QUICK(psi_e_p,r_e_p)
        call QUICK(psi_e_m,r_e_m)
        call QUICK(psi_w_p,r_w_p)
        call QUICK(psi_w_m,r_w_m)
    case ("Umist")
        call UMIST(psi_e_p,r_e_p)
        call UMIST(psi_e_m,r_e_m)
        call UMIST(psi_w_p,r_w_p)
        call UMIST(psi_w_m,r_w_m)
    case default
        call Van_Leer(psi_e_p,r_e_p)
        call Van_Leer(psi_e_m,r_e_m)
        call Van_Leer(psi_w_p,r_w_p)
        call Van_Leer(psi_w_m,r_w_m)
    end select

    end

!!!***************************************************
!!!*						         	               *
!!!*                Van Leer	       	   *
!!!*								                   *
!!!***************************************************
subroutine  Van_Leer(psi,r)
    implicit none
    real*8, intent(in) :: r
    real*8, intent(out) :: psi
    real*8 r_w_m, r_w_p, r_e_m, r_e_p

    psi = ( r + dabs(r) ) / ( 1.d0 + r )

    end

!!!***************************************************
!!!*						         	               *
!!!*                Van Albada	       	   *
!!!*								                   *
!!!***************************************************
subroutine  Van_Albada(psi,r)
    implicit none
    real*8, intent(in) :: r
    real*8, intent(out) :: psi
    real*8 r_w_m, r_w_p, r_e_m, r_e_p

    psi = ( r + r**2.d0 ) / ( 1.d0 + r**2.d0 )

    end

!!!***************************************************
!!!*						         	               *
!!!*                Min-Mod	       	   *
!!!*								                   *
!!!***************************************************
subroutine  Min_Mod(psi,r)
    implicit none
    real*8, intent(in) :: r
    real*8, intent(out) :: psi
    real*8 r_w_m, r_w_p, r_e_m, r_e_p

    if(r>0.d0) then
        psi = min(r,1.d0)
    else
        psi = 0.d0
    endif

    end

!!!***************************************************
!!!*						         	               *
!!!*                SuperBee       	   *
!!!*								                   *
!!!***************************************************
subroutine  SuperBee(psi,r)
    implicit none
    real*8, intent(in) :: r
    real*8, intent(out) :: psi
    real*8 r_w_m, r_w_p, r_e_m, r_e_p

    psi = max( 0.d0, min(2.d0*r,1.d0), min(r,2.d0) )

    end

!!!***************************************************
!!!*						         	               *
!!!*                QUICK       	   *
!!!*								                   *
!!!***************************************************
subroutine  QUICK(psi,r)
    implicit none
    real*8, intent(in) :: r
    real*8, intent(out) :: psi
    real*8 r_w_m, r_w_p, r_e_m, r_e_p

    psi = max( 0.d0, min( 2.d0*r, (3.d0+r)/4.d0, 2.d0 ) )

    end

!!!***************************************************
!!!*						         	               *
!!!*                UMIST      	   *
!!!*								                   *
!!!***************************************************
subroutine  UMIST(psi,r)
    implicit none
    real*8, intent(in) :: r
    real*8, intent(out) :: psi
    real*8 r_w_m, r_w_p, r_e_m, r_e_p

    psi = max( 0.d0, min( 2.d0*r, (1.d0+3.d0*r)/4.d0, (3.d0+r)/4.d0, 2.d0 ) )

    end

!!!*************************************************
!!!*						         	           *
!!!*                 Solve Kt                        *
!!!*								               *
!!!*************************************************
    
subroutine  solve_Kt_tvd(Kt,eps,dUdy,nut,dnutdy,detady,d2etady2,d3etady3,deta,sigmak,dsigmakdy,eps_hat,ny,upk,y)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: eps(1:ny),dUdy(1:ny),nut(1:ny),dnutdy(1:ny),detady(1:ny),d2etady2(1:ny),deta,sigmak(1:ny)
    real*8, intent(in) :: dsigmakdy(1:ny),eps_hat(1:ny),d3etady3(1:ny),y(1:ny)
    real*8, intent(inout) :: Kt(1:ny)
    real*8, intent(out) :: upk(1:ny)
    real*8 aK_ww,aK_w,aK_e,aK_ee,sK
    real*8 deps_hatdeta(1:ny), dupkdy(1:ny), depsdeta(1:ny),d2eps_hatdeta2(1:ny)
    real*8 nup(1:ny),d2Ktdeta2(1:ny),dKtdeta(1:ny),dnupdy(1:ny)
    real*8 depsdy(1:ny),D3KtDY3(1:ny),D2KtDY2(1:ny),dKtdy(1:ny)
    real*8 uk(1:ny),dukdy(1:ny),d2ukdy2(1:ny),dlnEdy(1:ny),Kt_star(1:ny)
    real*8 psi_w_m,psi_w_p,psi_e_m,psi_e_p
    character(len = 10) method
    integer j

    method = "Umist"

    call pressure_speed_Kt(upk,Kt,eps,eps_hat,detady,d2etady2,d3etady3,deta,y,ny)

    Kt(1) = 0.d0
    call Limiter_Function(psi_w_m, psi_w_p, psi_e_m, psi_e_p, Kt(1), Kt(1), Kt(2), Kt(3), Kt(4), method)
        call K_coefficients_tvd(aK_w,aK_e,sK,eps(2),nut(2),dnutdy(2),dUdy(2),deta,sigmak(2),dsigmakdy(2), &
            upk(1), upk(2), upk(3), eps_hat(2), d2etady2(2),detady(2),psi_w_m,psi_w_p,psi_e_m,psi_e_p, &
            Kt(1),Kt(2),Kt(3))
    Kt(2) = sK + aK_e*Kt(3) + aK_w*Kt(1)
    do j =3,ny-2
        call Limiter_Function(psi_w_m, psi_w_p, psi_e_m, psi_e_p, Kt(j-2), Kt(j-1), Kt(j), Kt(j+1), Kt(j+2), method)
        call K_coefficients_tvd(aK_w,aK_e,sK,eps(j),nut(j),dnutdy(j),dUdy(j),deta,sigmak(j),dsigmakdy(j), &
            upk(j-1), upk(j), upk(j+1), eps_hat(j), d2etady2(j),detady(j),psi_w_m,psi_w_p,psi_e_m,psi_e_p, &
            Kt(j-1),Kt(j),Kt(j+1))
        !if (aK_e < 0.d0) then
        !    print*, "aK_e=",aK_e, " j=",j
        !endif
        !if (aK_w < 0.d0) then
        !    print*, "aK_w=",aK_w, " j=",j
        !endif
        Kt(j) = sK + aK_e*Kt(j+1) + aK_w*Kt(j-1)
    enddo
    call Limiter_Function(psi_w_m, psi_w_p, psi_e_m, psi_e_p, Kt(ny-3), Kt(ny-2), Kt(ny-1), Kt(ny), Kt(ny), method)
        call K_coefficients_tvd(aK_w,aK_e,sK,eps(ny-1),nut(ny-1),dnutdy(ny-1),dUdy(ny-1),deta,sigmak(ny-1),dsigmakdy(ny-1), &
            upk(ny-2), upk(ny-1), upk(ny), eps_hat(ny-1), d2etady2(ny-1),detady(ny-1),psi_w_m,psi_w_p,psi_e_m,psi_e_p, &
            Kt(ny-2),Kt(ny-1),Kt(ny))
    Kt(ny-1) = sK + aK_e*Kt(ny) + aK_w*Kt(ny-2)
    Kt(ny) = 0.d0

    !Kt = Kt_star

    !upk = ( d2Ktdy2 - eps_hat ) * dKtdy / ( 2.d0 * Kt * (eps + eps_hat))
    !upk(1) = 0.d0*upk(2)
    !upk(ny) = 0.d0*upk(ny-1)

    end

!!!*************************************************
!!!*						         	           *
!!!*                 Pressue speed for K                        *
!!!*								               *
!!!*************************************************
    
subroutine  pressure_speed_Kt(upk,Kt,eps,eps_hat,detady,d2etady2,d3etady3,deta,y,ny)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: eps(1:ny),detady(1:ny),d2etady2(1:ny),y(1:ny),deta
    real*8, intent(in) :: Kt(1:ny),eps_hat(1:ny),d3etady3(1:ny)
    real*8, intent(out) :: upk(1:ny)
    real*8 deps_hatdy(1:ny), dupkdy(1:ny), depsdeta(1:ny),d2eps_hatdeta2(1:ny)
    real*8 d2Ktdeta2(1:ny),dKtdeta(1:ny)
    real*8 depsdy(1:ny),D3KtDY3(1:ny),D2KtDY2(1:ny),dKtdy(1:ny)
    real*8 uk(1:ny),dukdy(1:ny),d2ukdy2(1:ny),dlnEdy(1:ny),Kt_star(1:ny)
    real*8 psi_w_m,psi_w_p,psi_e_m,psi_e_p
    character(len = 10) method
        integer j

    do j=1,ny
        uk(j)= dsqrt(dabs(Kt(j)))
    enddo

    call ddy(ny,uk,dukdy,detady,deta)
    call ddy(ny,eps_hat,deps_hatdy,detady,deta)
    call d2dy2(ny,uk,d2ukdy2,detady,d2etady2,deta)
    call ddy(ny,dlog(dsqrt(dabs(eps_hat+eps))),dlnEdy,detady,deta)

    call d2dy2(ny,Kt,d2Ktdy2,detady,d2etady2,deta)
    call d3dy3(ny,Kt,D3KtDY3,detady,d2etady2,d3etady3,deta)
    call ddy(ny,Kt,dKtdy,detady,deta)
    call ddy(ny,eps,depsdy,detady,deta)

    !!! This is the problem
    !upk = 2.d0 * dukdy * d2ukdy2 / ( eps + eps_hat )
    upk = dlnEdy - depsdy / ( 2.d0 * ( eps + eps_hat ) )
    !upk = deps_hatdy / ( 2.d0 * ( eps + eps_hat ) )
    !upk = dKtdy * ( d2Ktdy2 - eps_hat ) / ( 2.d0 * Kt * (eps + eps_hat))
    !upk = 0.d0
    !upk(1) = 0.d0 ! This is verified
    !upk(ny) = 0.d0 ! This is verified
    !upk(1) = upk(2)
    !upk(ny) = upk(ny-1)

    !do j=1,100
    !call quadratic_variable_smoother(upk,upk,y,ny)
    !enddo

    end

!!!*************************************************
!!!*						         	             *
!!!*               q = D2Ady2                        *
!!!*								                 *
!!!*************************************************

subroutine  parabolic_system(ny,A_hat,A,detady,d2etady2,deta,m)
    implicit none
    integer, intent(in) :: ny,m
    real*8, intent(in) :: A(1:ny),detady(1:ny),d2etady2(1:ny),deta
    real*8, intent(out) :: A_hat(1:ny)
    real*8 deta2, S, aw, ae, aww, aee, awww, aeee, d2ady2(1:ny), q
    integer j

    call d2dy2(ny,A,D2ADY2,detady,d2etady2,deta)

    q = d2ady2(m)
    deta2 = deta**2.d0

    S = ( deta2 * q / (detady(1))**2.d0 ) / ( 2.d0 - 3.d0 * ( deta * d2etady2(1) / (2.d0*(detady(1))**2.d0) ) )
    ae = -( (-5.d0 + 4.d0 * ( deta * d2etady2(1) / (2.d0*(detady(1))**2.d0) ) ) &
    / ( 2.d0 - 3.d0 * ( deta * d2etady2(1) / (2.d0*(detady(1))**2.d0) ) ) )
    aee = -( ( 4.d0 - ( deta * d2etady2(1) / (2.d0*(detady(1))**2.d0) ) ) &
    / ( 2.d0 - 3.d0 * ( deta * d2etady2(1) / (2.d0*(detady(1))**2.d0) ) ) )
    aeee = -( - 1.d0 / ( 2.d0 - 3.d0 * ( deta * d2etady2(1) / (2.d0*(detady(1))**2.d0) ) ) )
    A_hat(1) = S + ae*A(2) + aee*A(3) + aeee*a(4) 
    
    do j=2,ny-1
        S = - deta2 * q / (2.d0*(detady(j))**2.d0)
        ae = (1.d0/2.d0 + deta * d2etady2(j) / (4.d0*(detady(j))**2.d0))
        aw = (1.d0/2.d0 - deta * d2etady2(j) / (4.d0*(detady(j))**2.d0))
        a_hat(j) = ae*a(j+1) + aw*a(j-1) + S
    enddo

    S = ( deta2 * q / (detady(ny))**2.d0) /  ( 2.d0 + 3.d0 * deta * d2etady2(ny) / (2.d0*(detady(ny))**2.d0) ) 
    aw = - ( (-5.d0 -4.d0 * deta * d2etady2(ny) / (2.d0*(detady(ny))**2.d0)) &
    / ( 2.d0 + 3.d0 * deta * d2etady2(ny) / (2.d0*(detady(ny))**2.d0) ) )
    aww = - ( (4.d0 + deta * d2etady2(ny) / (2.d0*(detady(ny))**2.d0)) &
    / ( 2.d0 + 3.d0 * deta * d2etady2(ny) / (2.d0*(detady(ny))**2.d0) ) )
    awww = - ( - 1.D0 / ( 2.d0 + 3.d0 * deta * d2etady2(ny) / (2.d0*(detady(ny))**2.d0) ) )
    a_hat(ny) = aw*a(ny-1) + aww*a(ny-2) + awww *a(ny-3) + S 
    
    end

!!!*************************************************
!!!*						         	           *
!!!*                 Solve Kt                        *
!!!*								               *
!!!*************************************************
    
subroutine  solve_Kt_b(Kt,eps,dUdy,nut,dnutdy,detady,d2etady2,d3etady3,deta,sigmak,dsigmakdy,eps_hat,ny,upk,y)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: eps(1:ny),dUdy(1:ny),nut(1:ny),dnutdy(1:ny),detady(1:ny),d2etady2(1:ny),deta,sigmak(1:ny)
    real*8, intent(in) :: dsigmakdy(1:ny),eps_hat(1:ny),d3etady3(1:ny),y(1:ny),upk(1:ny)
    real*8, intent(inout) :: Kt(1:ny)
    !real*8, intent(out) :: upk(1:ny)
    real*8 deltak,thetak,sK,aK_www,aK_ww,aK_w,aK_e,aK_ee,aK_eee
    real*8 dupkdy(1:ny),nup(1:ny),dnupdy(1:ny),depsdy(1:ny),D3KtDY3(1:ny),D2KtDY2(1:ny),dKtdy(1:ny),Pi_k(1:ny)
    real*8 uk(1:ny),dukdy(1:ny),d2ukdy2(1:ny),dlnEdy(1:ny),Kt_star(1:ny),deps_hatdy(1:ny),nup_p(1:ny),sa(1:ny)
    integer j

    call ddy(ny,Kt,dKtdy,detady,deta)
    call d2dy2(ny,Kt,d2Ktdy2,detady,d2etady2,deta)
    call d3dy3(ny,Kt,D3KtDY3,detady,d2etady2,d3etady3,deta)
    call ddy(ny,Kt,dKtdy,detady,deta)
    call ddy(ny,eps,depsdy,detady,deta)
    call ddy(ny,eps_hat,deps_hatdy,detady,deta)

    Kt(1) = 0.d0
    call K_coefficients_b(deltak,thetak,sK,dKtdy(2),d2Ktdy2(2),d3Ktdy3(2),depsdy(2),deps_hatdy(2),eps(2),eps_hat(2),nut(2), &
        sigmak(2), dnutdy(2), dsigmakdy(2),dUdy(2),detady(2),d2etady2(2),d3etady3(2),deta)
    !aK_ww,aK_w,aK_e,aK_ee
    aK_www = 0.d0
    aK_ww = 0.d0
    aK_w = ( 1.d0 / 2.d0 - deltak - 3.d0 * thetak / 2.d0 ) / ( 1.d0 - 5.d0 * thetak )
    aK_e = ( 1.d0 / 2.d0 + deltak - 6.d0 * thetak ) / ( 1.d0 - 5.d0 * thetak )
    aK_ee = 3.d0 * thetak / ( 1.d0 - 5.d0 * thetak )
    aK_eee = - ( thetak / 2.d0 ) / ( 1.d0 - 5.d0 * thetak )
    Kt(2) = ( sK  - ( thetak / 2.d0 ) * Kt(5) + 3.d0 * thetak * Kt(4) + ( 1.d0 / 2.d0 + deltak - 6.d0 * thetak ) * Kt(3) &
        + ( 1.d0 / 2.d0 - deltak - 3.d0 * thetak / 2.d0 ) * Kt(1) ) / ( 1.d0 - 5.d0 * thetak )
    do j =3,ny-2
        call K_coefficients_b(deltak,thetak,sK,dKtdy(j),d2Ktdy2(j),d3Ktdy3(j),depsdy(j),deps_hatdy(j),eps(j),eps_hat(j),nut(j), &
            sigmak(j), dnutdy(j), dsigmakdy(j),dUdy(j),detady(j),d2etady2(j),d3etady3(j),deta)
        aK_www = 0.d0
        aK_ww = - ( thetak / 2.d0 )
        aK_w = ( 1.d0 / 2.d0 - deltak + thetak )
        aK_e = ( 1.d0 / 2.d0 + deltak - thetak )
        aK_ee = ( thetak / 2.d0 )
        aK_eee = 0.d0
        if(aK_w < 0.d0) then 
            print*, "aK_w"
            print*, j,aK_w
        endif
        if(aK_e < 0.d0) then 
            print*, "aK_e"
            print*, j,aK_e
        endif
        Kt(j) = sK + aK_ee * Kt(j+2) + aK_e * Kt(j+1) + aK_w * Kt(j-1) + aK_ww * Kt(j-2)
    enddo
    call K_coefficients_b(deltak,thetak,sK,dKtdy(ny-1),d2Ktdy2(ny-1),d3Ktdy3(ny-1),depsdy(ny-1),deps_hatdy(ny-1),eps(ny-1), &
        eps_hat(ny-1), nut(ny-1),sigmak(ny-1), dnutdy(ny-1), dsigmakdy(ny-1),dUdy(ny-1),detady(ny-1),d2etady2(ny-1), &
        d3etady3(ny-1), deta)
    !aK_ww,aK_w,aK_e,aK_ee
    aK_www = ( thetak / 2.d0 ) / ( 1.d0 + 5.d0 * thetak )
    aK_ww = - 3.d0 * thetak / ( 1.d0 + 5.d0 * thetak )
    aK_w = ( 1.d0 / 2.d0 - deltak + 6.d0 * thetak ) / ( 1.d0 + 5.d0 * thetak )
    aK_e = ( 1.d0 / 2.d0 + deltak + 3.d0 * thetak / 2.d0 ) / ( 1.d0 + 5.d0 * thetak )
    aK_ee = 0.d0
    aK_eee = 0.d0
    Kt(ny-1) = ( sK  + ( 1.d0 / 2.d0 + deltak + 3.d0 * thetak / 2.d0 ) * Kt(ny) &
        + ( 1.d0 / 2.d0 - deltak + 6.d0 * thetak ) * Kt(ny-2) &
        - 3.d0 * thetak * Kt(ny-3) + ( thetak / 2.d0 ) * Kt(ny-4) ) / ( 1.d0 + 5.d0 * thetak )
    Kt(ny) = 0.d0

    end

!!!*************************************************
!!!*						         	           *
!!!*            Coeffcients Kt                        *
!!!*								               *
!!!*************************************************
    
subroutine  K_coefficients_b(deltak,thetak,sK,dKtdy,d2Ktdy2,d3Ktdy3,depsdy,deps_hatdy,eps,eps_hat,nut,sigmak,dnutdy, &
    dsigmakdy,dUdy,detady,d2etady2,d3etady3,deta)
    implicit none
    real*8, intent(in) :: dKtdy,d2Ktdy2,d3Ktdy3,depsdy,deps_hatdy,dUdy,detady,d2etady2,d3etady3,eps,eps_hat,nut,sigmak, &
        dnutdy,dsigmakdy,deta
    real*8, intent(out) :: deltak,thetak,sK
    real*8 Ak,Bk,Ck,Dk,alphak,betak,gammak

    Ak = - dKtdy / ( 2.d0 * ( eps_hat + eps ) )
    Bk = 1.d0 + nut / sigmak + ( eps_hat - d2Ktdy2 + dKtdy * ( depsdy + deps_hatdy ) / ( eps_hat + eps ) ) &
        / ( 2.d0 * ( eps_hat + eps ) ) 
    Ck = dnutdy / sigmak - nut * dsigmakdy / sigmak**2.d0 + deps_hatdy / ( 2.d0 * ( eps_hat + eps ) ) &
        - eps_hat * ( depsdy + deps_hatdy ) / ( 2.d0 * ( eps_hat + eps )**2.d0 )
    Dk = nut * ( dUdy )**2.d0 - eps_hat - eps
    alphak = Ak * detady**3.d0
    betak = 3.d0 * Ak * detady * d2etady2 + Bk * detady**2.d0
    gammak = Ak * d3etady3 + Bk * d2etady2 + Ck * detady

    deltak = gammak * deta / ( 4.d0 * betak )
    thetak = alphak / ( 2.d0 * betak * deta )
    Sk = deta * deta * Dk / ( 2.d0 * betak )

    end