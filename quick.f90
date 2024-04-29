!!!***************************************************
!!!*						         	               *
!!!*                K coefficients	       	   *
!!!*								                   *
!!!***************************************************
subroutine  K_coefficients(aK_ww,aK_w,aK_e,aK_ee,sK,eps,nut,dnutdy,nup,dnupdy,dUdy,deta,sigmak,dsigmakdy,Uww,Uw,Up,Ue,Uee,eps_hat, &
    d2etady2,detady)
    implicit none
    real*8, intent(in) :: eps,nut,dnutdy,nup,dnupdy,dUdy,deta,sigmak,d2etady2,detady,dsigmakdy,Uww,Uw,Up,Ue,Uee,eps_hat
    real*8, intent(out) :: aK_ww,aK_w,aK_e,aK_ee,sK
    real*8 diff, conv, b_ww, b_w, b_e, b_ee, b_p, D_w, D_e, F_ww, F_w, F_e, F_ee, F_p

    !DA_e(j) = (-2.d0*Uw*a(j-1) -3.d0*Up*a(j) + 6.d0*Ue*a(j+1) -Uee*a(j+2))/(6.d0*deta)
    !DA_w(j) = (Uww*a(j-2) - 6.d0*Uw*a(j-1) + 3.d0*Up*a(j) +2.d0*Ue*a(j+1))/(6.d0*deta)

    diff = (1.d0+nut/sigmak+nup)*(1.d0/deta**2.d0)*detady**2.d0
    conv = (1.d0/(2.d0*sigmak))*( (sigmak+nut+sigmak*nup)*d2etady2 + detady*(dnutdy -(nut/sigmak)*dsigmakdy +dnupdy) )/deta
    D_w = diff - conv
    D_e = diff + conv
    F_ww = Uww*detady/deta
    F_w = Uw*detady/deta
    F_e = Ue*detady/deta
    F_p = Up*detady/deta
    F_ee = Uee*detady/deta

    b_w = D_w + F_w/2.d0
    b_e = D_e - F_e/2.d0 
    b_p = D_w + D_e
    b_ee = 0.d0
    b_ww = 0.d0

    if ( (b_e < 0.d0) .and. (b_w >= 0.d0) ) then
        !print*, '(b_e < 0.d0) .and. (b_w >= 0.d0)'
        b_ee =  0.d0
        b_e = D_e - F_e/3.d0
        b_w = D_w + F_w
        b_ww = - F_ww/6.d0
        b_p = D_e + D_w + F_p/2.d0
    else if( (b_w < 0.d0) .and. (b_e >= 0.d0) ) then
        !print*, '(b_w < 0.d0) .and. (b_e >= 0.d0)'
        b_ww = 0.d0
        b_w = D_w + F_w/3.d0
        b_e = D_e - F_e
        b_ee = F_ee/6.d0
        b_p = D_w + D_e - F_p/2.d0
    else if( (b_w < 0.d0) .and. (b_e < 0.d0) ) then
        b_w = D_w + F_p/2.d0
        b_e = D_e - F_p/2.d0  
        b_p = D_w + D_e +(F_e - F_w)/2.d0
    endif

    !b_w = max(max(D_w + F_w/2.d0,F_w),0.d0)
    !b_e = max(max(D_e - F_e/2.d0,-F_e),0.d0) 
    !b_p = b_w + b_e + (F_e - F_w)/2.d0

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
    
subroutine  solve_Kt(Kt,eps,dUdy,nut,dnutdy,detady,d2etady2,d3etady3,deta,sigmak,dsigmakdy,eps_hat,ny,upk)
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

    nup = ( eps_hat - d2Ktdy2 ) / ( 2.d0 * ( eps + eps_hat ) )
    
    dnupdy = - 0*nup * ( depsdy / ( eps + eps_hat ) + 0*(1.d0 - 2.d0 *nup ) * dKtdy / Kt ) &
    - 0*d3Ktdy3 / ( 2.d0 * ( eps + eps_hat ) )
    nup = 0.d0
    !print*, 'nup(1) = ',nup(1),'; nup(ny) = ',nup(ny)
    !print*, 'nup(2) = ',nup(2),'; nup(ny-1) = ',nup(ny-1),'; dnupdy(2) = ',dnupdy(2),'; dnupdy(ny-1) = ',dnupdy(ny-1)

    dnupdy(1) = 0.d0
    dnupdy(ny) = 0.d0
    !call linear_variable_smoother(dnupdy,dnupdy,0.d0*dnupdy,ny)

    !!! This is the problem
    upk = ( d2Ktdy2 - eps_hat ) * dKtdy / ( 2.d0 * Kt * (eps + eps_hat))
    upk(1) = 2.d0*upk(2)
    upk(ny) = 2.d0*upk(ny-1)
    !upk = 0.d0

    Kt(1) = 0.d0
    Kt(2) = Kt(3)/4.d0
    do j =3,ny-2
        call K_coefficients(aK_ww,aK_w,aK_e,aK_ee,sK,eps(j),nut(j),dnutdy(j),nup(j),dnupdy(j),dUdy(j),deta,sigmak(j), &
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
