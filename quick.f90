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
    real*8 A_k,B_k,C_k,E_k,sigma_K, beta,gamma

    !DA_e(j) = (-2.d0*Uw*a(j-1) -3.d0*Up*a(j) + 6.d0*Ue*a(j+1) -Uee*a(j+2))/(6.d0*deta)
    !DA_w(j) = (Uww*a(j-2) - 6.d0*Uw*a(j-1) + 3.d0*Up*a(j) +2.d0*Ue*a(j+1))/(6.d0*deta)

    A_k = (1.d0+nut/sigmak)*detady**2.d0
    B_k = ( (1.d0+nut/sigmak)*d2etady2 + detady*(dnutdy -(nut/sigmak)*dsigmakdy)/sigmak )
    C_k = deta * detady / ( 4.d0 * A_k )
    E_k = deta * B_k / ( 4.d0 * A_k )
    sigma_K = (nut*dUdy*dUdy - eps_hat - eps) * (deta**2.d0) / ( 2.d0 * A_k)

    D_w = 5.d-01 - E_k
    D_e = 5.d-01 + E_k
    F_ww = C_k*Uww
    F_w = C_k*Uw
    F_e = C_k*Ue
    F_ee = C_k*Uee
    F_p = C_k*Up

    beta = 1.d0
    gamma = -1.d0

    b_w = D_w - gamma * F_w
    b_e = D_e - beta * F_e
    b_p = 1.d0 - F_p * ( beta + gamma ) * 3.d0 / 4.d0
    b_ee = - F_ee * ( 4.d0 -3.d0 * beta +  gamma ) / 8.d0
    b_ww = F_ww * ( 4.d0 +3.d0 * gamma - beta ) / 8.d0

    if ( (b_e < 0.d0) .and. (b_w >= 0.d0) ) then
        beta = D_e / F_e
        gamma = -1.d0
    else if( (b_w < 0.d0) .and. (b_e >= 0.d0) ) then
        gamma = D_w /F_w
        beta = 1.d0
    else if( (b_w < 0.d0) .and. (b_e < 0.d0) ) then
        beta = D_e / F_e
        gamma = D_w /F_w
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
    sK = (nut*dUdy*dUdy - eps_hat - eps) / b_p

    end

!!!*************************************************
!!!*						         	           *
!!!*                 Solve Kt                        *
!!!*								               *
!!!*************************************************
    
subroutine  solve_Kt_2(Kt,eps,dUdy,nut,dnutdy,detady,d2etady2,d3etady3,deta,sigmak,dsigmakdy,eps_hat,ny,upk)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: eps(1:ny),dUdy(1:ny),nut(1:ny),dnutdy(1:ny),detady(1:ny),d2etady2(1:ny),deta,sigmak(1:ny)
    real*8, intent(in) :: dsigmakdy(1:ny),eps_hat(1:ny),d3etady3(1:ny)
    real*8, intent(inout) :: Kt(1:ny)
    real*8, intent(out) :: upk(1:ny)
    real*8 aK_ww,aK_w,aK_e,aK_ee,sK
    real*8 deps_hatdeta(1:ny), dupkdy(1:ny), depsdeta(1:ny),d2eps_hatdeta2(1:ny)
    real*8 nup(1:ny),d2Ktdeta2(1:ny),dKtdeta(1:ny),dnupdy(1:ny)
    real*8 depsdy(1:ny),D3KtDY3(1:ny),D2KtDY2(1:ny),dKtdy(1:ny)
    integer j

    call d2deta2(ny,Kt,d2Ktdeta2,deta)
    call ddeta(ny,Kt,dKtdeta,deta)
    call d2dy2(ny,Kt,d2Ktdy2,detady,d2etady2,deta)
    call d3dy3(ny,Kt,D3KtDY3,detady,d2etady2,d3etady3,deta)
    call ddeta(ny,eps,depsdeta,deta)
    call ddeta(ny,eps_hat,deps_hatdeta,deta)
    call d2deta2(ny,eps_hat,d2eps_hatdeta2,deta)
    call ddy(ny,Kt,dKtdy,detady,deta)
    call ddy(ny,eps,depsdy,detady,deta)

    !!! This is the problem
    upk = ( d2Ktdy2 - eps_hat ) * dKtdy / ( 2.d0 * Kt * (eps + eps_hat))
    upk(1) = 2.d0*upk(2)
    upk(ny) = 2.d0*upk(ny-1)
    !upk = 0.d0

    Kt(1) = 0.d0
    Kt(2) = Kt(3)/4.d0
    do j =3,ny-2
        call K_coefficients_2(aK_ww,aK_w,aK_e,aK_ee,sK,eps(j),nut(j),dnutdy(j),dUdy(j),deta,sigmak(j), &
            dsigmakdy(j), upk(j-2), upk(j-1),upk(j), upk(j+1), upk(j+2), eps_hat(j), d2etady2(j),detady(j))
        if (aK_e < 0.d0) then
            print*, "aK_e=",aK_e, " j=",j
        endif
        if (aK_w < 0.d0) then
            print*, "aK_w=",aK_w, " j=",j
        endif
        Kt(j) = sK + aK_ee*Kt(j+2) + aK_e*Kt(j+1) + aK_w*Kt(j-1) + aK_ww*Kt(j-2)
    enddo
    Kt(ny-1) = Kt(ny-2)/4.d0
    Kt(ny) = 0.d0

    !upk = d3Ktdy3

    end
