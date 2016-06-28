

!! Here we have CO cooling in the parametrization of Neufeld et al., 
!! (based only on temperture and NCOtilde)
!! but using data constructed from DESPOTIC.
!!subroutine CO_cooling_neufeld(y, temperature, lambda, rpar)
subroutine co_cooling(yco, yh2, ye, temperature, nh, divv, lambda)

  implicit none

  real, intent(IN) :: yco, yh2, ye, temperature, nh, divv
  real, intent(OUT) :: lambda
  
  integer, parameter :: ntemp = 14
  integer, parameter :: nnco = 11
  
   
  real, dimension(ntemp) :: logtempv = (/ 0.69897,1.0,1.30103,1.47712125,1.69897, &
       1.90308999, 2.0, 2.30103, 2.47712125, 2.69897, 2.90308999, 3.0, 3.30103, 3.47712125 /)


  real, dimension(nnco) :: logncov = (/ 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, &
       12.0, 12.5, 13.0, 13.5, 14.0 /)
   
  real, dimension(ntemp) :: logL0v = (/ 24.721, 24.261, 24.0216, & 
       23.931, 23.836, 23.746, 23.699, 23.518, & 
       23.387, 23.206, 23.034, 22.950, 22.691, 22.562 /)
  



  ! here, columns are constant temperature
  ! rows are constant Ntilde


  real, dimension(ntemp,nnco) :: logLLTEv = reshape((/ 21.943, 21.085, 20.350, 19.944, 19.449, 19.006, 18.799, 18.170, 17.813, 17.376, 17.010, 16.864, 16.552, 16.451,  &
21.959 ,  21.093 ,  20.354 ,  19.947 ,  19.450 ,  19.007 ,  18.800 ,  18.171 ,  17.813 ,  17.376 ,  17.010 ,  16.864 ,  16.552 ,  16.451 ,  &
22.007 ,  21.118 ,  20.367 ,  19.955 ,  19.455 ,  19.010 ,  18.802 ,  18.172 ,  17.814 ,  17.376 ,  17.011 ,  16.864 ,  16.552 ,  16.451 ,  &
22.139 ,  21.192 ,  20.404 ,  19.980 ,  19.470 ,  19.019 ,  18.810 ,  18.176 ,  17.816 ,  17.378 ,  17.012 ,  16.865 ,  16.552 ,  16.451 ,  &
22.409 ,  21.370 ,  20.509 ,  20.053 ,  19.516 ,  19.048 ,  18.833 ,  18.187 ,  17.824 ,  17.382 ,  17.015 ,  16.867 ,  16.554 ,  16.453 ,  &
22.770 ,  21.673 ,  20.731 ,  20.227 ,  19.636 ,  19.129 ,  18.900 ,  18.222 ,  17.847 ,  17.397 ,  17.024 ,  16.875 ,  16.560 ,  16.457 ,  &
23.198 ,  22.045 ,  21.051 ,  20.514 ,  19.872 ,  19.314 ,  19.060 ,  18.318 ,  17.915 ,  17.439 ,  17.052 ,  16.900 ,  16.577 ,  16.471 ,  &
23.619 ,  22.443 ,  21.420 ,  20.862 ,  20.191 ,  19.602 ,  19.330 ,  18.522 ,  18.076 ,  17.551 ,  17.133 ,  16.971 ,  16.630 ,  16.514 ,  &
24.062 ,  22.874 ,  21.817 ,  21.243 ,  20.552 ,  19.944 ,  19.662 ,  18.817 ,  18.343 ,  17.774 ,  17.322 ,  17.150 ,  16.780 ,  16.640 ,  &
24.543 ,  23.302 ,  22.230 ,  21.645 ,  20.940 ,  20.317 ,  20.029 ,  19.161 ,  18.671 ,  18.082 ,  17.637 ,  17.480 ,  17.121 ,  16.948 ,  &
25.011 ,  23.757 ,  22.661 ,  22.062 ,  21.345 ,  20.712 ,  20.418 ,  19.533 ,  19.032 ,  18.439 ,  18.050 ,  17.929 ,  17.605 ,  17.423 /),(/ntemp,nnco/))





  real, dimension(ntemp,nnco) :: lognhalfv = reshape((/ 2.643 ,  2.918 ,  3.333 ,  3.628 ,  4.024 ,  4.380 ,  4.542 ,  4.982 ,  5.198 ,  5.441 ,  5.656 ,  5.753 ,  5.953 ,  5.977 , &  
2.621 ,  2.902 ,  3.324 ,  3.621 ,  4.020 ,  4.378 ,  4.540 ,  4.981 ,  5.198 ,  5.441 ,  5.656 ,  5.752 ,  5.953 ,  5.976 ,  &
2.553 ,  2.853 ,  3.294 ,  3.600 ,  4.007 ,  4.370 ,  4.534 ,  4.977 ,  5.195 ,  5.439 ,  5.655 ,  5.752 ,  5.953 ,  5.976 ,  &
2.365 ,  2.715 ,  3.207 ,  3.538 ,  3.967 ,  4.344 ,  4.513 ,  4.967 ,  5.188 ,  5.435 ,  5.652 ,  5.749 ,  5.951 ,  5.975 ,  &
1.974 ,  2.394 ,  2.980 ,  3.365 ,  3.853 ,  4.269 ,  4.452 ,  4.934 ,  5.165 ,  5.421 ,  5.643 ,  5.742 ,  5.947 ,  5.972 ,  &
1.480 ,  1.923 ,  2.563 ,  3.001 ,  3.574 ,  4.066 ,  4.281 ,  4.837 ,  5.098 ,  5.377 ,  5.615 ,  5.719 ,  5.935 ,  5.961 ,  &
0.981 ,  1.424 ,  2.071 ,  2.520 ,  3.123 ,  3.669 ,  3.918 ,  4.590 ,  4.912 ,  5.252 ,  5.531 ,  5.651 ,  5.896 ,  5.928 ,  &
0.481 ,  0.924 ,  1.571 ,  2.022 ,  2.628 ,  3.182 ,  3.437 ,  4.158 ,  4.533 ,  4.953 ,  5.311 ,  5.465 ,  5.782 ,  5.830 ,  &
-0.019 ,  0.424 ,  1.071 ,  1.522 ,  2.129 ,  2.684 ,  2.939 ,  3.666 ,  4.048 ,  4.493 ,  4.898 ,  5.084 ,  5.504 ,  5.583 ,  &
-0.519 ,  -0.076 ,  0.571 ,  1.022 ,  1.629 ,  2.184 ,  2.439 ,  3.167 ,  3.550 ,  3.998 ,  4.408 ,  4.599 ,  5.053 ,  5.150 ,  &
-1.019 ,  -0.576 ,  0.071 ,  0.522 ,  1.129 ,  1.684 ,  1.939 ,  2.667 ,  3.051 ,  3.499 ,  3.910 ,  4.102 ,  4.559 ,  4.656 /),(/ntemp,nnco/))

  





  real, dimension(ntemp,nnco) :: alphav = reshape((/ 0.651 ,  0.621 ,  0.586 ,  0.559 ,  0.526 ,  0.505 ,  0.507 ,  0.496 ,  0.494 ,  0.487 ,  0.472 ,  0.460 ,  0.402 ,  0.378 , &  
0.645 ,  0.618 ,  0.582 ,  0.555 ,  0.524 ,  0.503 ,  0.505 ,  0.495 ,  0.493 ,  0.487 ,  0.472 ,  0.460 ,  0.402 ,  0.378 ,  &
0.630 ,  0.608 ,  0.571 ,  0.545 ,  0.516 ,  0.497 ,  0.500 ,  0.491 ,  0.490 ,  0.485 ,  0.471 ,  0.458 ,  0.401 ,  0.377 ,   &
0.612 ,  0.589 ,  0.547 ,  0.522 ,  0.496 ,  0.481 ,  0.485 ,  0.481 ,  0.482 ,  0.480 ,  0.466 ,  0.455 ,  0.399 ,  0.375 ,   &
0.638 ,  0.577 ,  0.519 ,  0.490 ,  0.463 ,  0.451 ,  0.456 ,  0.459 ,  0.464 ,  0.465 ,  0.456 ,  0.445 ,  0.392 ,  0.370 ,  &
0.677 ,  0.593 ,  0.521 ,  0.484 ,  0.447 ,  0.427 ,  0.427 ,  0.426 ,  0.431 ,  0.437 ,  0.432 ,  0.425 ,  0.378 ,  0.359 , & 
0.693 ,  0.610 ,  0.536 ,  0.499 ,  0.460 ,  0.434 ,  0.429 ,  0.410 ,  0.406 ,  0.405 ,  0.400 ,  0.394 ,  0.356 ,  0.344 ,  &
0.710 ,  0.621 ,  0.548 ,  0.512 ,  0.475 ,  0.451 ,  0.447 ,  0.422 ,  0.412 ,  0.398 ,  0.384 ,  0.375 ,  0.341 ,  0.338 ,  &
0.718 ,  0.627 ,  0.555 ,  0.521 ,  0.487 ,  0.465 ,  0.461 ,  0.438 ,  0.427 ,  0.412 ,  0.395 ,  0.385 ,  0.352 ,  0.359 , & 
0.718 ,  0.631 ,  0.560 ,  0.527 ,  0.495 ,  0.475 ,  0.472 ,  0.450 ,  0.440 ,  0.425 ,  0.407 ,  0.397 ,  0.372 ,  0.388 ,  &
0.716 ,  0.629 ,  0.561 ,  0.530 ,  0.501 ,  0.483 ,  0.480 ,  0.460 ,  0.451 ,  0.436 ,  0.413 ,  0.401 ,  0.373 ,  0.393 /),(/ntemp,nnco/))





  real :: LLTE, L0, nH2, NCO, sigh2, ve
  real :: nhalf, alpha, tff, log_nco, LCO
  real :: log_temp, n_co, dvdz
  real :: t, u, y1, y2, y3, y4
  integer :: temp_index, nco_index
  real :: dx, slope
  
  
!! effective nH2 - see Glover et al. 2010
  sigh2 = 3.3e-16 / (temperature / 1000.0)**0.25
  ve = 1.03e4 * sqrt(temperature)
  nH2 = (yh2 + 1.414 * (2.3e-15 / sigh2) + (1.3e-8 / sigh2 / ve) * ye) * nh

  
  log_temp = log10(temperature)

  ! the NCO here (read: NCOtilde) must be the optical depth parameter
  ! introduced in Neufeld et al.
  ! NCO = n(CO) / |dv/dr|
  ! measured in cm^-2 per (cm/s) OR s/cm^3 (all CGS)

  NCO = nh * yco / divv
  log_nco = log10(NCO)
 
  
   
  call locate(logtempv, ntemp, log_temp, temp_index) 
  if(temp_index .le. 0) then
     temp_index = 1
     !log_temp = logtempv(temp_index)
  endif
  if(temp_index .ge. ntemp) then
     temp_index = ntemp-1
     !log_temp = logtempv(ntemp)
  endif
  call locate(logncov, nnco, log_nco, nco_index)
  if(nco_index .le. 0) then
     nco_index = 1
     log_nco = logncov(nco_index)
  endif
  if(nco_index .ge. nnco) then
     nco_index = nnco-1
     log_nco = logncov(nnco)
  endif
  t = (log_temp - logtempv(temp_index)) / ( logtempv(temp_index+1) - logtempv(temp_index))
  u = (log_nco - logncov(nco_index)) / ( logncov(nco_index+1) - logncov(nco_index))


  y1 = logLLTEv(temp_index, nco_index)
  y2 = logLLTEv(temp_index+1,nco_index)
  y3 = logLLTEv(temp_index+1,nco_index+1)
  y4 = logLLTEv(temp_index,nco_index+1)
  LLTE = 10.0**(-((1.0-t)*(1.0-u)*y1 + t*(1.0-u)*y2  + t*u*y3 + (1.0-t)*u*y4))


  y1 = alphav(temp_index, nco_index)
  y2 = alphav(temp_index+1,nco_index)
  y3 = alphav(temp_index+1,nco_index+1)
  y4 = alphav(temp_index,nco_index+1)
  alpha = (1.0-t)*(1.0-u)*y1 + t*(1.0-u)*y2 + t*u*y3 + (1.0-t)*u*y4

  y1 = lognhalfv(temp_index, nco_index)
  y2 = lognhalfv(temp_index+1,nco_index)
  y3 = lognhalfv(temp_index+1,nco_index+1)
  y4 = lognhalfv(temp_index,nco_index+1)
  nhalf = 10.0**((1.0-t)*(1.0-u)*y1 + t*(1.0-u)*y2  + t*u*y3 + (1.0-t)*u*y4)

  slope = (logL0v(temp_index+1) - logL0v(temp_index)) / ( logtempv(temp_index+1) - logtempv(temp_index))
  dx = log_temp - logtempv(temp_index)
  L0 = 10.0**(-(logL0v(temp_index) + slope * dx))

  
  LCO = 1.0 / (1.0 / L0 + nH2 / LLTE + (nH2 / nhalf)**alpha * (1.0 - nhalf*L0/LLTE) / L0)
  lambda =  nH2 * yco * nh * LCO


  return
end subroutine co_cooling




! see Lee et al. 1996, table 11
subroutine co_shielding(NH2, NCO, fh2, fco)

  implicit none
    
  real, intent(IN) :: NH2, NCO
  real, intent(OUT) :: fh2, fco
  
  integer, parameter :: n_h2 = 42
  integer, parameter :: n_co = 51
  
  integer :: i
  real :: dx, dydx
  
  real, dimension(n_h2) :: NH2v = (/ 2.666e13, 3.801e14, 6.634e15, 8.839e16, &
       9.268e17, 1.007e18, 2.021e18, 3.036e18, 4.051e18, 5.066e18, 6.082e18, &
       7.097e18, 8.112e18, 9.341e18, 1.014e19, 2.030e19, 3.045e19, 4.061e19, &
       5.076e19, 6.092e19, 7.107e19, 8.123e19, 9.353e19, 1.015e20, 2.031e20, & 
       3.047e20, 4.062e20, 5.078e20, 6.094e20, 7.109e20, 8.125e20, 9.355e20, &
       1.016e21, 2.031e21, 3.047e21, 4.063e21, 5.078e21, 6.094e21, 7.110e21, & 
       8.125e21, 9.355e21, 1.016e22 /)
 real, dimension(n_h2) :: fh2v = (/ 0.999, 0.9893, 0.9678, 0.9465, 0.9137, 0.9121, 0.8966, & 
       0.8862, 0.8781, 0.8716, 0.8660, 0.8612, 0.8569, 0.8524, 0.8497, 0.8262, 0.8118, & 
       0.8011, 0.7921, 0.7841, 0.7769, 0.7702, 0.7626, 0.7579, 0.7094, 0.6712, 0.6378, & 
       0.6074, 0.5791, 0.5524, 0.5271, 0.4977, 0.4793, 0.2837, 0.1526, 7.774e-2, 3.952e-2, & 
       2.093e-2, 1.199e-2, 7.666e-3, 5.333e-3, 4.666e-3  /)
  real, dimension(n_co) :: NCOv = (/ 1.0e12, 1.65e12, 2.995e12, 5.979e12, 1.313e13, &
       3.172e13, 8.429e13, 2.464e14, 7.923e14, 1.670e15, 2.595e15, 4.435e15, 6.008e15, &
       8.952e15, 1.334e16, 1.661e16, 2.274e16, 3.115e16, 4.266e16, 5.843e16, 8.002e16, & 
       1.096e17, 1.501e17, 2.055e17, 2.815e17, 4.241e17, 6.389e17, 9.625e17, 1.450e18, & 
       2.184e18, 3.291e18, 4.124e18, 5.685e18, 7.838e18, 1.080e19, 1.285e19, 1.681e19, & 
       2.199e19, 2.538e19, 3.222e19, 4.091e19, 5.193e19, 5.893e19, 7.356e19, 8.269e19, & 
       9.246e19, 1.031e20, 1.148e20, 1.277e20, 1.419e20, 1.578e20 /)
  real, dimension(n_co) :: fcov = (/ 0.999, 0.9981, 0.9961, 0.9912, 0.9815, 0.9601, 0.9113, & 
       0.8094,   0.6284,   0.4808,   0.3889,   0.2827,   0.2293,   0.1695,   0.1224, & 
       0.1017,   7.764e-2, 5.931e-2, 4.546e-2, 3.506e-2, 2.728e-2, 2.143e-2, 1.700e-2, & 
       1.360e-2, 1.094e-2, 8.273e-3, 6.283e-3, 4.773e-3, 3.611e-3, 2.704e-3, 1.986e-3, & 
       1.657e-3, 1.258e-3, 9.332e-4, 6.745e-4, 5.596e-4, 4.123e-4, 2.982e-4, 2.490e-4, & 
       1.827e-4, 1.324e-4, 9.473e-5, 7.891e-5, 5.668e-5, 4.732e-5, 3.967e-5, 3.327e-5, & 
       2.788e-5, 2.331e-5, 1.944e-5, 1.619e-5  /)
  
  !if(NH2 .ge. 1.016e22) then
  !   fh2 = 4.66e-3
  !else 

  if(NH2 .le. NH2v(1)) then
     fh2 = 1.0
  else
     call locate(NH2v, n_h2, NH2, i)

     !if(i .eq. 0 .or. i .ge. n_h2) then
     !   print*, "shouldn't be here. i for CO / H2 shielidng out of range", i, NH2
     !   call abort(1)
     !endif
     ! extrapolate if H2 column is larger than maximum
     if(i .ge. n_h2) then
        i = n_h2-1
     endif
     
     dx = NH2 - NH2v(i)
     dydx = (fh2v(i+1) - fh2v(i)) / (NH2v(i+1) - NH2v(i))
     fh2 = fh2v(i) + dx * dydx
  endif

  
  !if(NCO .ge. 1.578e20) then
  !   fco = 1.62e-5
  !else 

  if(NCO .le. 1.0e12) then
     fco = 1.0
  else
     call locate(NCOv, n_co, NCO, i)
     !if(i .eq. 0 .or. i .ge. n_co) then
     !   print*, "shouldn't be here. i for CO / CO shielidng out of range", i, NCO
     !   call abort(1)
     !endif

     ! extrapolate if CO column is larger than maximum
     if(i .ge. n_co) then
        i = n_co-1
     endif

     dx = NCO - NCOv(i)
     dydx = (fcov(i+1) - fcov(i)) / (NCOv(i+1) - NCOv(i))
     fco = fcov(i) + dx * dydx
  endif
  

  return
end subroutine co_shielding






subroutine fine_structure_cooling(yc, ycp, yo, yh, ye, yhp, yh2, temp, nh, NHtot, lambda_fs)

  implicit none
  

  
  real, intent(IN) :: yc, ycp, yo, yh, ye, yhp, yh2, temp, nh, NHtot
  real, intent(OUT) :: lambda_fs
  
  
  !! CII parameters
  real, parameter :: A10_CII = 2.3e-6
  real, parameter :: E10_CII = 1.27e-14
  real, parameter :: g1_CII = 4.0
  real, parameter :: g0_CII = 2.0
  
  !! OI parameters
  real, parameter :: A10_OI = 8.9e-5
  real, parameter :: A20_OI = 1.3e-10
  real, parameter :: A21_OI = 1.8e-5
  real, parameter :: E10_OI = 3.17e-14
  real, parameter :: E20_OI = 4.55e-14
  real, parameter :: E21_OI = 1.35e-14
  real, parameter :: g0_OI = 5.0
  real, parameter :: g1_OI = 3.0
  real, parameter :: g2_OI = 1.0
  
  !! CI parameters
  real, parameter :: A10_CI = 7.9e-8
  real, parameter :: A20_CI = 2.1e-14
  real, parameter :: A21_CI = 2.7e-7
  real, parameter :: E10_CI = 3.31e-15
  real, parameter :: E20_CI = 8.69e-15
  real, parameter :: E21_CI = 5.38e-15
  real, parameter :: g0_CI = 1.0
  real, parameter :: g1_CI = 3.0
  real, parameter :: g2_CI = 5.0

  real :: f0, f1, f2, w
  real :: beta, ngamma
  real :: C10, C21, C20
  real :: lambda_Cplus, Lambda_C, lambda_O
  
  real :: ng20,ng21,ng10,t20,t21,t10
  real :: b20,b21,b10
  logical :: error
  
  real :: f0_old, f1_old, f2_old, Q
  integer, parameter :: maxiter = 50
  real, parameter :: tolerance = 1.0e-4
  logical :: converged
  integer :: i
  
  real :: ortho, para, opratio, fortho, fpara, T2, ln_temp
  real :: t22, t23, t24, t25, t26, t27, t28, t29, t30
  real :: t31, t32, t33, t34, t35, t36, t37, t38, t39, t40
  real :: t41, t42, t43, t44, t45, t46, t47, t48, t49, t50
  
  logical :: fine_structure_line_opacity = .false.
  
  real :: TCMB = 2.725
  real, parameter :: kb = 1.38e-16

  
  opratio = 2.4
  fortho = opratio / (1.0 + opratio)
  fpara = 1.0 / (1.0 + opratio)
  T2 = 1.0e-2 * temp
  ln_temp = log(temp)
  ! CII:
    
    ! collisions with e-
    if(temp .le. 2.0e3) then
       t22 = 3.86e-7 * T2**(-0.5)
    else
       t22 = 2.43e-7 * T2**(-0.345)
    end if

    ! collisions with H
    if(temp .le. 2.0e3) then
       t23 = 8.0e-10 * T2**(0.07) 
    else
       t23 = 3.1e-10 * T2**(0.385)
    end if
    
    ! collisions with H2
    if(temp .le. 250.0) then
       ortho = (4.7e-10 + 4.6e-13 * temp) * fortho
       para = 2.5e-10 * temp**0.12 * fpara
    else
       !ortho = 5.85e-10 * temp**0.07 * fortho
       !para = 4.85e-10 * temp**0.07 * fpara
       ortho = 3.975e-10 * temp**0.07 * fortho
       para = 3.295e-10 * temp**0.07 * fpara
    end if
    t24 = ortho + para



  
  ! OI:
    
    ! collisions with H2
    
    
    ! 1->0
    ortho = 2.7e-11 * temp**0.362 * fortho
    para = 3.46e-11 * temp**0.316 * fpara
    t25 = ortho + para
    
    ! 2->0
    ortho = 5.49e-11 * temp**0.317 * fortho 
    para = 7.07e-11 * temp**0.268 * fpara
    t26 = ortho + para

    ! 2->1
    ortho = 2.74e-14 * temp**1.060 * fortho
    para = 3.33e-15 * temp**1.360 * fpara
    t27 = ortho + para
    

    ! collisions with H

    ! 1->0
    t28 = 9.2e-11 * T2**0.67
    
    ! 2->0
    t29 = 4.3e-11 * T2**0.80

    ! 2->1
    t30 = 1.1e-10 * T2**0.44
    
    
    ! collisions with e-

    ! 1->0
    t31 = 5.12e-10 * temp**(-0.075)
    
    ! 2->0
    t32 = 4.86e-10 * temp**(-0.026)

    ! 2->1
    t33 = 1.08e-14 * temp**(0.926)


    ! collisions with H+

    ! 1->0
    if(temp .le. 194.0) then
       t34 = 6.38e-11 * temp**0.4
    else if(temp .le. 3686.0) then
       t34 = 7.75e-12 * temp**0.8 
    else
       t34 = 2.65e-10 * temp**0.37
    end if
    
    ! 2->0
    if(temp .le. 511.0) then
       t35 = 6.10e-13 * temp**1.1
    else if(temp .le. 7150.0) then
       t35 = 2.12e-12 * temp**0.9
    else
       t35 = 49e-10 * temp**0.3
    end if
 

    ! 2->1
    if(temp .le. 2090.0) then
       t36 = 2.03e-11 * temp**0.56
    else
       t36 = 3.43e-10 * temp**0.19
    end if




  ! CI:
    
    ! collisions with H2
    
    ! 1->0
    ortho = 8.7e-11 - 6.6e-11 * exp(-temp / 218.3) + 6.6e-11 * exp(-2.0 * temp / 218.3) * fortho
    para = 7.9e-11 - 8.7e-11 * exp(-temp / 126.4) + 1.3e-10 * exp(-2.0 * temp / 126.4) * fpara
    t39 = ortho + para
    
    ! 2->0
    ortho = 1.2e-10 - 6.1e-11 * exp(-temp / 387.3) * fortho
    para = 1.1e-10 - 8.6e-11 * exp(-temp / 223.0) + 8.7e-11 * exp(-2.0 * temp / 223.0) * fpara
    t40 = ortho + para

    ! 2->1
    ortho = 2.9e-10 - 1.9e-10 * exp(-temp / 348.9) * fortho
    para = 2.7e-10 - 2.6e-10 * exp(-temp / 250.7) + 1.8e-10 * exp(-2.0 * temp / 250.7) * fpara
    t41 = ortho + para
    

    ! collisions with H

    ! 1->0
    t42 = 1.6e-10 * T2**0.14
    
    ! 2->0
    t43 = 9.2e-11 * T2**0.26

    ! 2->1
    t44 = 2.9e-10 * T2**0.26
    
    
    ! collisions with e-

    ! 1->0
    if(temp .le. 1000.0) then
       t45 = 2.88e-6 * temp**(-0.5) * exp(-9.25141 - 7.73782e-1 * ln_temp + & 
            3.61184e-1 * ln_temp**2 - 1.50892e-2 * ln_temp**3 - 6.56325e-4 * ln_temp**4)
    else if(temp .le. 1.0e4) then
       !! slightly wrong expression given in Glover & Jappsen, corrected here --- see Johnson et al. 1987
       t45 = 2.88e-6 * temp**(-0.5) * exp(4.446e2 - 2.27913e2 * ln_temp + & 
            4.2595e1 * ln_temp**2 - 3.47620 * ln_temp**3 + 1.0508e-1 * ln_temp**4)
    else
       t45 = 0.0
    end if
    
    
    ! 2->0
    if(temp .le. 1000.0) then
       t46 = 1.73e-6 * temp**(-0.5) * exp(-7.69735 - 1.30743 * ln_temp + & 
            0.697638 * ln_temp**2 - 0.111338 * ln_temp**3 + 0.705277e-2 * ln_temp**4)
    else if(temp .le. 1.0e4) then
       t46 = 1.73e-6 * temp**(-0.5) * exp(3.50609e2 - 1.87474e2 * ln_temp + & 
            3.61803e1 * ln_temp**2 - 3.03283 * ln_temp**3 + 9.38138e-2 * ln_temp**4)
    else
       t46 = 0.0
    end if
    
    ! 2->1

    if(temp .le. 1000.0) then
       t47 = 1.73e-6 * temp**(-0.5) * exp(-7.4387 - 0.57443 * ln_temp + & 
            0.358264 * ln_temp**2 - 4.18166e-2 * ln_temp**3 + 2.35272e-3 * ln_temp**4)
    else if(temp .le. 1.0e4) then
       t47 = 1.73e-6 * temp**(-0.5) * exp(3.86186e2 - 2.02193e2 * ln_temp + & 
            3.85049e1 * ln_temp**2 - 3.19268 * ln_temp**3 + 9.78573e-2 * ln_temp**4)
    else
       t47 = 0.0
    end if


    ! collisions with H+

    ! 1->0
    if(temp .le. 5000.0) then
       t48 = temp**0.45 * (9.6e-11 - 1.8e-14 * temp + 1.9e-18 * temp**2)
    else
       t48 = 8.9e-10 * temp**0.117
    end if
    
    ! 2->0
    if(temp .le. 5000.0) then
       t49 = temp * (3.1e-12 - 6.0e-16 * temp + 3.9e-20 * temp**2)
    else
       t49 = 2.3e-9 * temp**0.0965
    end if
 

    ! 2->1
    if(temp .le. 5000.0) then
       t50 = temp**0.7 * (1.0e-10 - 2.2e-14 * temp + 1.7e-18 * temp**2)
    else
       t50 = 9.2e-9 * temp**0.0535
    end if



  !! C+, two level system:
  
  lambda_Cplus = 0.0
  if(ycp .gt. 0.0) then
     beta = 1.0 
     ng10 = 1.0 / (exp(E10_CII / kb / TCMB) - 1.0)
     
     C10 = (  ye  * t22  & 
          + yh  * t23 & 
          + yh2 * t24) * nh
     
     call two_level_system(temp,A10_CII,C10,beta,ng10,g1_CII/g0_CII,E10_CII,f1,f0) 
     
     if(fine_structure_line_opacity) then
     do i=1,maxiter
        
        f1_old = f1
        f0_old = f0
        
        call calculate_beta(ycp,12.0,temp,f1,f0,E10_CII,A10_CII,g1_CII/g0_CII,NHtot,beta)
        call two_level_system(temp,A10_CII,C10,beta,ng10,g1_CII/g0_CII,E10_CII,f1,f0) 
        
        ! agressively force convergence:
        w = 1.0 - real(i) / real(maxiter)
        f0 = (1.0 - w) * f0_old + w * f0
        f1 = (1.0 - w) * f1_old + w * f1
        
        ! check for convergence
        converged = .false.
        converged = (abs(f1_old - f1) / (f1 * tolerance) .lt. 1.0) & 
             .and.  (abs(f0_old - f0) / (f0 * tolerance) .lt. 1.0)
        
        if(converged) exit
        
        
     enddo
     
     if(.not. converged) then
        !print*, "C+ did not converge! wtf", temp
        !print*, y
        stop
     endif
     endif
     
     lambda_Cplus = beta * ((1.0 + ng10) * f1 - (g1_CII / g0_CII) * ng10 * f0) * A10_CII * E10_CII * ycp
  endif
  

  ! neutral C - three level system:

  lambda_C = 0.0
  if(yc .gt. 0.0) then
  C10 = (ye * t45 + & 
       yh * t42 + & 
       yhp * t48 + & 
       yh2 * t39) * nh
  C20 = (ye * t46 + & 
       yh * t43 + & 
       yhp * t49 + & 
       yh2 * t40) * nh
  C21 = (ye * t47 + & 
       yh * t44 + & 
       yhp * t50 + & 
       yh2 * t41) * nh
  
  b20 = 1.0
  b21 = 1.0
  b10 = 1.0

  ng20 = 1.0 / (exp(E20_CI / kb / TCMB) - 1.0)
  ng21 = 1.0 / (exp(E21_CI / kb / TCMB) - 1.0)
  ng10 = 1.0 / (exp(E10_CI / kb / TCMB) - 1.0)

  
  t20 = exp(-E20_CI / (kb * temp))
  t21 = exp(-E21_CI / (kb * temp))
  t10 = exp(-E10_CI / (kb * temp))
  
  
  error = .false.
  call three_level_system(A20_CI,A21_CI,A10_CI,C20,C21,C10,b20,b21,b10,&
       ng20,ng21,ng10,g2_CI/g0_CI,g2_CI/g1_CI,g1_CI/g0_CI,t20,t21,t10,f0,f1,f2,error)
  if(error) then
     !print*, "error out of three level system CI. y=", y
     print*, "temp = ", temp
     stop
  endif
  
  


  if(fine_structure_line_opacity) then
  do i=1,maxiter
     
     f2_old = f2
     f1_old = f1
     f0_old = f0
     
     call calculate_beta(yc,12.0,temp,f2,f0,E20_CI,A20_CI,g2_CI/g0_CI,NHtot,b20)
     call calculate_beta(yc,12.0,temp,f2,f1,E21_CI,A21_CI,g2_CI/g1_CI,NHtot,b21)
     call calculate_beta(yc,12.0,temp,f1,f0,E10_CI,A10_CI,g1_CI/g0_CI,NHtot,b10)
     error = .false.
     call three_level_system(A20_CI,A21_CI,A10_CI,C20,C21,C10,b20,b21,b10,&
          ng20,ng21,ng10,g2_CI/g0_CI,g2_CI/g1_CI,g1_CI/g0_CI,t20,t21,t10,f0,f1,f2,error)
     if(error) then
        !print*, "error out of three level system. y=", y
        !print*, "temp = ", temp
        stop
     endif

          
     ! agressively force convergence:
     w = 1.0 - real(i) / real(maxiter)
     f0 = (1.0 - w) * f0_old + w * f0
     f1 = (1.0 - w) * f1_old + w * f1
     f2 = (1.0 - w) * f2_old + w * f2

     ! check for convergence
     converged = .false.
     converged = (abs(f1_old - f1) / (f1 * tolerance) .lt. 1.0) & 
          .and.  (abs(f0_old - f0) / (f0 * tolerance) .lt. 1.0) & 
          .and.  (abs(f2_old - f2) / (f2 * tolerance) .lt. 1.0) 
     
     if(converged) exit
     
     
  enddo

  if(.not. converged) then
     print*, "CI did not converge! wtf"
     !print*, "y = ", y
     print*, "temp = ", temp
     stop
  endif
  endif

  lambda_C = (b20 * ((1.0 + ng20) * f2 - (g2_CI / g0_CI) * ng20 * f0) * A20_CI * E20_CI & 
            + b21 * ((1.0 + ng21) * f2 - (g2_CI / g1_CI) * ng21 * f1) * A21_CI * E21_CI  &
            + b10 * ((1.0 + ng10) * f1 - (g1_CI / g0_CI) * ng10 * f0) * A10_CI * E10_CI) * yc
    
  endif


  ! neutral O - three level system:
  
  
  C10 = (ye * t31 + & 
       yh * t28 + & 
       yhp * t34 + & 
       yh2 * t25) * nh
  C20 = (ye * t32 + & 
       yh * t29 + & 
       yhp * t35 + & 
       yh2 * t26) * nh
  C21 = (ye * t33 + & 
       yh * t30 + & 
       yhp * t36 + & 
       yh2 * t27) * nh
  
  b20 = 1.0
  b21 = 1.0
  b10 = 1.0

  ng20 = 1.0 / (exp(E20_OI / kb / TCMB) - 1.0)
  ng21 = 1.0 / (exp(E21_OI / kb / TCMB) - 1.0)
  ng10 = 1.0 / (exp(E10_OI / kb / TCMB) - 1.0)
  
  t20 = exp(-E20_OI / (kb * temp))
  t21 = exp(-E21_OI / (kb * temp))
  t10 = exp(-E10_OI / (kb * temp))
  
  error = .false.
  call three_level_system(A20_OI,A21_OI,A10_OI,C20,C21,C10,b20,b21,b10,&
       ng20,ng21,ng10,g2_OI/g0_OI,g2_OI/g1_OI,g1_OI/g0_OI,t20,t21,t10,f0,f1,f2,error)
  if(error) then
     !print*, "error out of three level system OI. y=", y
     !print*, "temp = ", temp
     stop
  endif


  if(fine_structure_line_opacity) then
  do i=1,maxiter
     
     f2_old = f2
     f1_old = f1
     f0_old = f0
     
     call calculate_beta(yo,16.0,temp,f2,f0,E20_OI,A20_OI,g2_OI/g0_OI,NHtot,b20)
     call calculate_beta(yo,16.0,temp,f2,f1,E21_OI,A21_OI,g2_OI/g1_OI,NHtot,b21)
     call calculate_beta(yo,16.0,temp,f1,f0,E10_OI,A10_OI,g1_OI/g0_OI,NHtot,b10)
     error = .false.
     call three_level_system(A20_OI,A21_OI,A10_OI,C20,C21,C10,b20,b21,b10,&
          ng20,ng21,ng10,g2_OI/g0_OI,g2_OI/g1_OI,g1_OI/g0_OI,t20,t21,t10,f0,f1,f2,error)
     if(error) then
        !print*, "error out of three level system. y=", y
        !print*, "temp = ", temp
        stop
     endif


          
     ! agressively force convergence:
     w = 1.0 - real(i) / real(maxiter)
     f0 = (1.0 - w) * f0_old + w * f0
     f1 = (1.0 - w) * f1_old + w * f1
     f2 = (1.0 - w) * f2_old + w * f2

     ! check for convergence
     converged = .false.
     converged = (abs(f1_old - f1) / (f1 * tolerance) .lt. 1.0) & 
          .and.  (abs(f0_old - f0) / (f0 * tolerance) .lt. 1.0) & 
          .and.  (abs(f2_old - f2) / (f2 * tolerance) .lt. 1.0) 
     
     if(converged) exit
     
     !print*, i
     
  enddo

  if(.not. converged) then
     !print*, "OI did not converge! wtf"
     stop
  endif
  endif

  lambda_O = (b20 * ((1.0 + ng20) * f2 - (g2_OI / g0_OI) * ng20 * f0) * A20_OI * E20_OI  & 
           + b21 * ((1.0 + ng21) * f2 - (g2_OI / g1_OI) * ng21 * f1) * A21_OI * E21_OI  &
           + b10 * ((1.0 + ng10) * f1 - (g1_OI / g0_OI) * ng10 * f0) * A10_OI * E10_OI) * yo
 

  
  lambda_fs = (lambda_Cplus + lambda_C + lambda_O) * nh

 
 return
 
end subroutine fine_structure_cooling



subroutine calculate_beta(y,mu,temp,f1,f0,e,a,gfac,NHtot,beta)

  implicit none

  real, intent(IN) :: y, mu, temp, f1, f0, e, a, gfac, NHtot
  real, intent(OUT) :: beta

  real, parameter :: hc = 1.99e-16
  real, parameter :: mh = 1.67e-24
  real, parameter :: kb = 1.38e-16
  
  real :: tauline, cs, sigma_tot, lambda
  
  if(NHtot .le. 0.0) then
     beta = 0.0
     return
  endif
  
  lambda = hc / e
  
  ! assume no non-thermal component
  cs = sqrt(kb * temp / (mh * mu))
  sigma_tot = cs
  
  tauline = gfac * a * lambda**3 / sigma_tot  * y * NHtot * f0 * (1.0 - f1 / f0 / gfac) / 63.0
  
  if(tauline .le. 1.0e-5) then
     beta = 1.0
  else
     beta = (1.0 - exp(-3.0 * tauline)) / (3.0 * tauline)
  endif
  
  return
end subroutine calculate_beta





subroutine two_level_system(t,a,c,b,ng,gfac,e,f1,f0)


  implicit none

  real, intent(IN) :: t,a,c,b,gfac,e,ng
  real, intent(OUT) :: f1,f0
  
  real :: alpha, gamma
  real :: kb = 1.38e-16
  
  alpha = a * b * (1.0 + ng) + c
  gamma = gfac * (c * exp(-e / (kb * t)) + a * ng * b)
  
  f1 = gamma / (alpha + gamma)
  f0 = 1.0 - f1
  
  if(f0 .lt. 0.0) then
     print*, "problem in two level system. negative level pop", f1, f0
     stop
  endif
     
  return
end subroutine two_level_system




! use cramers rule to solve three level systems
subroutine three_level_system(a20,a21,a10,c20,c21,c10,b20,b21,b10, & 
                    ng20,ng21,ng10,g20,g21,g10,t20,t21,t10,f0,f1,f2,error)



  real, intent(IN) :: a20,a21,a10,c20,c21,c10,b20,b21,b10
  real, intent(IN) :: ng20,ng21,ng10,g20,g21,g10,t20,t21,t10
  real, intent(OUT) :: f0,f1,f2
  logical, intent(INOUT) :: error
  
  real :: q1, q2, q3, q4, q5, q6
  real :: d, dx, dy, dz
  
  
  q1 = c21 + a20 * b20 * (1.0 + ng20) + a21 * b21 * (1.0 + ng21)
  q2 = c21 * g21 * t21 + g21 * a21 * ng21 * b21
  q3 = c20 * g20 * t20 + g20 * a20 * ng20 * b20
  q4 = c10 + c21 * g21 * t21 + b10 * a10 * (1.0 + ng10) + g21 * a21 * ng21 * b21
  q5 = b21 * a21 * c21 + a21 * ng21 * b21
  q6 = g10 * a10 * ng10 * b10 + g10 * c10 * t10
  
  d =   q3 * (-q4 - q5) & 
      - q2 * (q6 - q5) & 
      - q1 * (q6 + q4)
  
  dx = - q2 * (-q5) & 
       - q1 * (q4)
  
  dy = q3 * (-q5) & 
       -q1 * (q6)
  
  dz = -q3 * q4 & 
       -q2 * q6
  
  f0 = dx/d
  f1 = dy/d
  f2 = dz/d

  
  if(f0 .lt. 0.0 .or. f1 .lt. 0.0 .or. f2 .lt. 0.0) then
     print*, "problem in three level system. negative level pop", f0, f1, f2
     error = .true.
  endif
  
  
  return
end subroutine three_level_system






SUBROUTINE locate(xx,n,x,j)
  INTEGER j,n
  REAL x,xx(n)
  INTEGER jl,jm,ju
  jl=0
  ju=n+1
10 if(ju-jl.gt.1)then
     jm=(ju+jl)/2
     if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
        jl=jm
     else
        ju=jm
     endif
     goto 10
  endif
  j=jl
  return
END SUBROUTINE locate
