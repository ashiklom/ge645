	PROGRAM DISORD1D
*
* this simple program illustrates numerical solution of the 1D RTE
* with the discrete ordinate finite kernel method as described in
* chapter 4 of GG 645 (http://cybele.bu.edu/courses/gg645/main.html)
*
*  - assumes PLANOPHILE leaf normal inclination distribution
*  - assumes UNIFORM leaf normal azimuthal distribution
* ___________________________________________________________________
*
*   standard constants
*
	PARAMETER (PI = 3.141592654)
	PARAMETER (degtorad = 3.141592654/180.0)
	PARAMETER (radtodeg = 180.0/3.141592654)
	PARAMETER (ng = 8)
	PARAMETER (Nlayers = 100)
	PARAMETER (Epsilon = 0.0001)
*
        REAL xg(ng), wg(ng), gL(ng), hL(ng), theta_o
	REAL phi_o,  mu_o,   Ftot,   fdir,   I_o	
	REAL I_d,    rho_Ld, tau_Ld, LAI,    DeltaL
        REAL R_s,    Gdir,   Gdif(ng,ng),    dummy(ng,ng)
*
	REAL Gamma_d_dir(ng,ng),      Gamma_d_dif(ng,ng,ng,ng)
	REAL I_o_uc_d(Nlayers),       I_o_uc_u(Nlayers,ng,ng)
	REAL I_d_uc_d(Nlayers,ng,ng), I_d_uc_u(Nlayers,ng,ng)
        REAL F_o_uc_d_soil,           F_d_uc_d_soil
	REAL Q(Nlayers,ng,ng),        S(Nlayers,ng,ng)
	REAL Ic(Nlayers+1,ng,ng),     Ic_old(ng,ng)
*
	REAL HR_uc,   HR_c,  HR
	REAL HT_uc,   HT_c,  HT
	REAL AB_uc,   AB_c,  AB
	REAL theta_v, phi_v, RF
*
	LOGICAL convergence
* ___________________________________________________________________
*
*
* ng            : gauss quadrature order (4,6,8,10 or 12)
* Nlayers       : number of spatial nodes (100)
* Epsilon       : convergence criterion for iteration of the 
*                 scattering integral
* xg            : gauss ordinates (between -1 to 1)
* wg            : gauss weights (between 0 and 1)
* gL            : pdf of leaf normal inclination (assumed planophile)
* hL            : pdf of lear normal azimuths (assumed uniform)
* theta_o       : polar angle of the sun (between 90 and 180 degrees)
* phi_o         : azimuthal angle of the sun (between 0 and 360 
*                  degrees)
* mu_o          : cosine of theta_o
* Ftot          : total incident flux density (1 W/m2)
* fdir          : fraction of Ftot in the direct solar beam 
* I_o           : intensity of the direct beam
* I_d           : intensity of diffuse sky light (assumed isotropic)
* rho_LD        : leaf hemispherical reflectance (diffuse internal 
*                 scattering)
* tau_LD        : leaf hemispherical transmittance (diffuse internal 
*                 scattering)
* LAI           : one-sided leaf area per unit ground area
* DeltaL        : thickness of sptatial cells (LAI/Nlayers)
* R_s           : soil hemipsherical reflectance (assumed Lambertian)
* Gdir          : Geometry factor for direct solar radiation
* Gdif          : Geometry factor for scattered radiation field
* dummy         : well, a dummy variable
* Gamma_d_dir   : area scattering phase function for direct solar 
*                 radiation
* Gamma_d_dif   : area scattering phase function for scattered 
*                 radiation field
* I_o_uc_d      : downward uncollided direct solar radiation
* I_o_uc_u      : upward uncollided direct solar radiation
* I_d_uc_d      : downward uncollided diffuse sky radiation
* I_d_uc_u      : upward uncollided diffuse sky radiation
* F_o_uc_d_soil : donward uncollided flux density of direct solar 
*                 radiation incident on the ground below the canopy
* F_d_uc_d_soil : donward uncollided flux density of diffuse sky
*                 radiation incident on the ground below the canopy
* Q             : first collision source 
* S             : multiple collision source
* Ic            : collided intensity field
* Ic_old        : canopy leaving collided intensity field from  
*                 previous iteration
* HR_uc         : uncollided hemispherical reflectance
* HR_c          : collided hemispherical reflectance
* HT_uc         : uncollided hemispherical transmittance
* HT_c          : uncollided hemispherical reflectance
* AB_uc         : uncollided canopy absortance
* AB_c          : collided canopy absorptance
* HR            : hemispherical reflectance (DHR if fdir=1, else BHR)
* HT            : hemispherical transmittance
* AB            : canopy absorptance
* convergence   : logical flag to test for convergence of the  
*                 iteration on the multiple collision source
* theta_v       : view polar angle
* phi_v         : view azimuthal angle
* RF            : reflectance factor (BRF if fdir=1, else HDRF)
* ___________________________________________________________________
*
* BEGIN INPUTS
*
	theta_o = 120.0
	phi_o   = 0.0
	theta_o = theta_o*degtorad
	phi_o   = phi_o*degtorad
	mu_o    = cos(theta_o)
     	Ftot    = 1.0
	fdir    = 0.7
	I_o     = Ftot*(fdir/(abs(mu_o)))
        I_d     = Ftot*(1-fdir)/PI 
	LAI     = 3.0
	DeltaL  = LAI/Nlayers
	rho_Ld  = 0.45
	tau_Ld  = 0.45
	R_s     = 0.3
*
* END INPUTS
* ___________________________________________________________________
*
* get cross sections
*
	CALL XSECTIONS (ng,mu_o,phi_o,rho_LD,tau_LD,xg,wg,gL,hL,dummy,
     *                  Gdir,Gdif,Gamma_d_dir,Gamma_d_dif)
* ___________________________________________________________________
*
* Evaluate Uncollided direct solar radiation (I_o_uc_d, I_o_uc_u)
*
	CALL I_o_uncol_down (DeltaL,Nlayers,Gdir,I_o,mu_o,I_o_uc_d)
	CALL I_o_uncol_up   (DeltaL,Nlayers,Gdir, Gdif,I_o,mu_o,ng,
     *                       xg,R_s,F_o_uc_d_soil,I_o_uc_u)
* ___________________________________________________________________
*
* Evaluate Uncollided diffuse sky radiation (I_d_uc_d, I_d_uc_u)
*
	CALL I_d_uncol_down (DeltaL,Nlayers,Gdif,I_d,ng,xg,I_d_uc_d)
	CALL I_d_uncol_up   (DeltaL,Nlayers,Gdif,I_d,ng,xg,wg,R_s,
     *                       F_d_uc_d_soil,I_d_uc_u)
* ___________________________________________________________________
*
* Evaluate First-Collision Source Q
*
	CALL FCS (Nlayers,ng,xg,wg,Gamma_d_dir,Gamma_d_dif,
     *            I_o_uc_d, I_o_uc_u,I_d_uc_d, I_d_uc_u,Q)
* ___________________________________________________________________
*
* Iterate on the Multiple-Collision Source S
*
	do ims = 1, 100
*
* add multiple-collision source to first-collision source
*
	   do i = 1, ng
	   do j = 1, ng
	   do k = 1, Nlayers
	      S(k,j,i) = Q(k,j,i)+S(k,j,i)
	   end do
	   end do
	   end do
*
* sweep downwards plus handle the bottom boundary condition
*
	   CALL SWEEP_DOWN (Nlayers,ng,xg,wg,Gdif,DeltaL,S,R_s,Ic)
*
* sweep upwards and check for convergence
*
	   CALL SWEEP_UP (Nlayers,ng,xg,wg,Gdif,DeltaL,S,Ic,Ic_old,
     *                    epsilon,convergence)
*
* evaluate Multiple-Collision Source
*
	   if (convergence) goto 100
	   CALL  MULTI_COLL_S (Nlayers,ng,xg,wg,Gamma_d_dif,Ic,S)

	end do
*
100     continue
* ___________________________________________________________________
*
* do energy balance
*
	CALL ENERGY_BAL (Nlayers,ng,xg,wg,mu_o,Q,S,DeltaL,R_s,rho_LD,
     *                   tau_LD,Gdir,Gdif,F_o_uc_d_soil, F_d_uc_d_soil,
     *                   I_o_uc_d,I_d_uc_d, I_o_uc_u,I_d_uc_u,Ic,
     *                   HR_uc,HR_c,HT_uc,HT_c,AB_uc,AB_c)
* ___________________________________________________________________
*
* write output
*
	HR = HR_uc + HR_c
	HT = HT_uc + HT_c
	AB = AB_uc + AB_c
*
	write(*,101) HR_uc
	write(*,102) HR_c
	write(*,103) HR
	write(*,104) HT_uc
	write(*,105) HT_c
	write(*,106) HT
	write(*,107) AB_uc
	write(*,108) AB_c
	write(*,109) AB
	write(*,110) (1.-R_S)*HT_uc
	write(*,111) (1.-R_S)*HT_c
	write(*,112) (1.-R_S)*HT
	write(*,113) HR + AB + (1.-R_S)*HT
*
	write(*,*)
	do i = (ng/2)+1, ng
	   theta_v = xg(i)*radtodeg
	   do j = 1, ng
	      phi_v = xg(j)*PI + PI
	      phi_v = phi_v*radtodeg
	      RF    = Ic(1,j,i)*PI
	      write(*,114) theta_v, phi_v, RF	   
	   end do
	end do
*
* ___________________________________________________________________
*
* format statements
*
101	format(' Uncollided Hemispherical Reflectance     : ', F8.6)
102	format(' Collided Hemispherical Reflectance       : ', F8.6)
103	format(' Hemispherical Reflectance                : ', F8.6)
104	format(' Uncollided Hemispherical Transmittance   : ', F8.6)
105	format(' Collided Hemispherical Transmittance     : ', F8.6)
106	format(' Hemispherical Transmittance              : ', F8.6)
107	format(' Uncollided Canopy Absorptance            : ', F8.6)
108	format(' Collided Canopy Absorptance              : ', F8.6)
109	format(' Canopy Absorptance                       : ', F8.6)
110	format(' Uncollided Soil Absorptance              : ', F8.6)
111	format(' Collided Soil Absorptance                : ', F8.6)
112	format(' Soil Absorptance                         : ', F8.6)
113	format(' Energy Balance (=1.0)                    : ', F8.6)
114     format(' Theta_view, Phi_view, Reflectance Factor : ',
     *           F9.6,3X,F10.6,3X,F8.6)
* ___________________________________________________________________
*
* Done
* 
	STOP
	END
*------------------------------------------------------------------------------*


*------------------------------------------------------------------------------*
*
	SUBROUTINE XSECTIONS (ng,mu_o,phi_o,rho_LD,tau_LD,xg,wg,
     *                        gL,hL,dummy,Gdir,Gdif,Gamma_d_dir,
     *                        Gamma_d_dif)
*
* this routine evaluates cross sections for the discrete ordinates code
*
*  - assumes PLANOPHILE leaf normal inclination distribution
*  - assumes UNIFORM leaf normal azimuthal distribution
*
*   declarations
*
	PARAMETER (PI = 3.141592654)
*
        REAL xg(ng), wg(ng), gL(ng), hL(ng)
	REAL phi_o, mu_o, rho_Ld, tau_Ld
        REAL Gdir, Gdif(ng,ng), dummy(ng,ng)
	REAL Gamma_d_dir(ng,ng), Gamma_d_dif(ng,ng,ng,ng)
*
	INTEGER ng
*
* get quadrature
*
	CALL gauss_quad (ng,xg,wg)
*
* check the quadrature
*
	CALL check_quad (ng,xg,wg)
*
* example integral
*
	CALL example_integral (ng,xg,wg)
*
* get pdf of leaf normal orientation gL(thetaL) and hL(phiL)
*
	CALL leaf_normal_pdf (ng,xg,wg,gL,hL)
*
* get G_FUNCTION for direct solar radiation
*
	CALL G_dir_function (ng,xg,wg,gL,hL,mu_o,phi_o,Gdir)
*
* get G_FUNCTION for all quadrature directions (phi,mu)
*
       CALL G_dif_function (ng,xg,wg,gL,hL,Gdif)
*
* get GAMMA_d_dir (phi_o,mu_o -> phi,mu)
*
	CALL GAMMA_d_dir_function (ng,xg,wg,gL,hL,rho_Ld,tau_Ld, 
     *                             mu_o,phi_o,Gamma_d_dir) 
	CALL CHECK_Gamma_d_dir    (ng,xg,wg,rho_Ld,tau_Ld,Gdir,
     *                             Gamma_d_dir)
*
* get GAMMA_d_dif (phiprime,muprime -> phi,mu) 
*
	CALL GAMMA_d_dif_function (ng,xg,wg,gL,hL,rho_Ld,tau_Ld, 
     *                             Gdif,Gamma_d_dif,dummy)
* Done
* 
	RETURN
	END
*------------------------------------------------------------------------------*



*------------------------------------------------------------------------------*
*
        SUBROUTINE gauss_quad  (ng,xg,wg)
*
* program for obtaining gauss quadrature of order ng
*
* inputs:
*   ng: quadrature order (integer)
*   ng must be 4, 6, 8, 10 OR 12
*
* outputs:
*   xg: ordinates (real)
*   wg: weights   (real)
*
*
* begin declarations
*
        REAL xg(12), wg(12), xx(20), ww(20)
        INTEGER ishift(6), ng, ng2, i
*
        DATA xx/-0.861136312,-0.339981044,-0.9324695,-0.6612094,
     *          -0.2386192,-0.960289856,-0.796666477,-0.525532410,
     *          -0.183434642,
     *          -0.973906529,-0.865063367,-0.679409568,-0.433395394,
     *          -0.148874339,
     *          -0.981560634,-0.904117256,-0.769902674,-0.587317954,
     *          -0.367831499,-0.125233409/
*
        DATA ww/ 0.347854845, 0.652145155,0.1713245, 0.3607616,
     *           0.4679139, 0.101228536, 0.222381034, 0.313706646,
     *           0.362683783,
     *           0.066671344,0.149451349,0.219086363,0.269266719,
     *           0.295524225,
     *           0.047175336,0.106939326,0.160078329,0.203167427,
     *           0.233492537,0.249147046/
*
        DATA ishift/0,0,2,5,9,14/
*
* end declarations
*
*
* set values for xg and wg
*
 	ng2 = ng/2
*
        DO i = 1, ng2
           xg(i)  = xx(i+ishift(ng2))
           wg(i)  = ww(i+ishift(ng2))
        END DO
*
        DO i = ng2+1, ng
           xg(i)  = -xg(ng+1-i)
           wg(i)  =  wg(ng+1-i)
        END DO
* 
* done
*
        RETURN
        END
*------------------------------------------------------------------------------*


*------------------------------------------------------------------------------*
*
	SUBROUTINE check_quad (ng, xg, wg)
*
* this routine -  
*
*  (1) checks if the quadrature weights sum to 2.0
*  (2) checks if [ int_0^1 dx x ] is equal to 0.5
* 
* begin declarations
*
*
        REAL    xg(ng), wg(ng), sum
        INTEGER ng, i
*
* end declarations
*
*
* check if the weights sum to 2.0
*
	sum = 0.0
	DO i = 1, ng
	   sum = sum + wg(i)
	END DO
c	WRITE (*,*) " Qwts check (=2.0?): ", sum
*
* check if [ int_0^1 dx x ] is equal to 0.5
*
	sum = 0.0
	DO i = (ng/2)+1, ng
	   sum = sum + (xg(i)*wg(i))
	END DO
c	WRITE (*,*) " Qord check (=0.5?): ", sum
* 
* done
*
        RETURN
        END
*------------------------------------------------------------------------------*


*------------------------------------------------------------------------------*
*
	SUBROUTINE example_integral (ng,xg,wg)
*
* this routine illustrates -
*
* how to do the integral whose lower and upper bounds are A and B
* for example [ int_A^B dx |x| ] , A = -1 and B = +1
*
* begin declarations
*
*
        REAL      xg(ng), wg(ng), sum, conv1, conv2, neword
	REAL      upperlimit, lowerlimit
*
	INTEGER   i
*
* end declarations
*
*  step 1: define conversion factors conv1 and conv2
*
	upperlimit =  1.0
	lowerlimit = -1.0
  	conv1 = (upperlimit-lowerlimit)/2.0
	conv2 = (upperlimit+lowerlimit)/2.0
*
*  step 2: do the integral by making sure the ordinates run between 
*          the upperlimit and the lowerlimit
*
	sum = 0.0
	do i = 1, ng
	   neword = conv1*xg(i) + conv2
	   sum    = sum + abs(neword*wg(i))
	END DO
*
*  step 3: make sure not to forget to apply the conversion factor 1 again
*
	sum = sum*conv1
c	WRITE (*,*) " Intg check (=1.0?): ", sum
	sum = 0.0
* 
* done
*
        RETURN
        END
*------------------------------------------------------------------------------*



*------------------------------------------------------------------------------*
*
	SUBROUTINE leaf_normal_pdf(ng,xg,wg,gL,hL)
* 
* this routine evaluates PLANOPHILE leaf normal inclination pdf (gL)
* and UNIFORM leaf normal azimuthal pdf (hL)
* 
* begin declarations
*
	PARAMETER (PI = 3.141592654)
*
        REAL     xg(ng), wg(ng), gL(ng), hL(ng)
	REAL     upperlimit, lowerlimit, conv1, conv2, sum
        REAL     neword
*
	INTEGER  i
*
* end declarations
*
* set hL = 1.0
*
	DO i = 1, ng
	   hL(i) = 1.0
	END DO
*
* obtain the planophile gL
*
*   gL(thetaL) = (2/pi) (1+cos(2thetaL))
*
* and at the same time check if it satisfies the required condition
* of normalization, that is, 
*
*   [ int_0^(pi/2) dthetaL gL(thetaL) = 1.0 ]
*
*  step 1: define limits of integration and the convertion factors
*
	upperlimit = PI/2.0
	lowerlimit = 0.0
  	conv1 = (upperlimit-lowerlimit)/2.0
	conv2 = (upperlimit+lowerlimit)/2.0
*
*  step 2: do the integral by making sure the ordinates run between 
*          the upperlimit and the lowerlimit
*
	sum = 0.0
	do i = 1, ng
	   neword = conv1*xg(i) + conv2
	   gL(i)  = (2.0/PI)*(1.0+COS(2.0*neword))
	   sum    = sum + gL(i)*wg(i)
	END DO
*
*  step 3: make sure not to forget to apply the conversion factor 1 again
*
	sum = sum*conv1
c	WRITE (*,*) " LNO  check (=1.0?): ", sum
* 
* done
*
        RETURN
        END
*------------------------------------------------------------------------------*


*------------------------------------------------------------------------------*
*
	SUBROUTINE G_dir_function (ng,xg,wg,gL,hL,
     *                              muprime,phiprime,Gdir)
* 
* this routine calculates the G_FUNCTION given a direction of
* photon travel OMEGA^PRIME and the leaf normal pdf
* 
* begin declarations
*
	PARAMETER (PI = 3.141592654)
*
        REAL    xg(ng), wg(ng), gL(ng), hL(ng)
        REAL    muprime, phiprime
	REAL    mu_tp, sin_tp
	REAL    upperlimit_tL, lowerlimit_tL
        REAL    conv1_tL, conv2_tL, sum_tL
	REAL    upperlimit_pL, lowerlimit_pL
  	REAL    conv1_pL, conv2_pL, sum_pL
 	REAL    neword_tL, mu_tL, sin_tL
	REAL    neword_pL, dotproduct, Gdir
*
	INTEGER i, j
*
* end declarations
*
	mu_tp  = muprime
	sin_tp = sqrt(1.0 - muprime*muprime)
*
* define limits of integration and the convertion factors for integration
* over thetaL (note the tL suffix!)
*
	upperlimit_tL = PI/2.0
	lowerlimit_tL = 0.0
  	conv1_tL = (upperlimit_tL-lowerlimit_tL)/2.0
	conv2_tL = (upperlimit_tL+lowerlimit_tL)/2.0
*
* define limits of integration and the convertion factors for integration
* over phiL (note the pL suffix!)
*
	upperlimit_pL = 2.0*PI
	lowerlimit_pl = 0.0
  	conv1_pL = (upperlimit_pL-lowerlimit_pL)/2.0
	conv2_pL = (upperlimit_pL+lowerlimit_pL)/2.0
*
* integral over theta_L
*
	sum_tL = 0.0
	do i = 1, ng
	   neword_tL = conv1_tL*xg(i) + conv2_tL
	   mu_tL     = cos(neword_tL)
	   sin_tL    = sin(neword_tL)
*
* integral over phi_L
*
	   sum_pL = 0.0
	   do j = 1, ng
	      neword_pL  = conv1_pL*xg(j) + conv2_pL
              dotproduct = abs (mu_tL*mu_tp + 
     *                           sin_tL*sin_tp*cos(neword_pL-phiprime))
              sum_pL     = sum_pL + wg(j)*hL(j)/(2.0*PI)*dotproduct
           end do
* 
* finish the phi_L integral
*
	   sum_pL = sum_pL*conv1_pL
           sum_tL = sum_tL + wg(i)*gL(i)*sum_pL
	end do
* 
* finish the theta_L integral
*
	sum_tL = sum_tL*conv1_tL
	Gdir  = sum_tL
*
* done
* 
	RETURN
	END
*------------------------------------------------------------------------------*



*------------------------------------------------------------------------------*
*
	SUBROUTINE G_dif_function (ng,xg,wg,gL,hL,Gdif)
*
* this routine evaluates the G function in all quadrature directions 
* and checks for normalization
* 
* begin declarations
*
	PARAMETER (PI = 3.141592654)
*
        REAL xg(ng), wg(ng)
	REAL gL(ng), hL(ng)
        REAL phiprime, muprime
	REAL upperlimit_pp, lowerlimit_pp, conv1_pp, conv2_pp
	REAL Gdif(ng,ng), sum_tp, sum_pp
*
	INTEGER i, j
*
* end declarations
*
*  conversion factors to have the ordinates simulate phiprime
*
	upperlimit_pp = 2.0*PI
	lowerlimit_pp = 0.0
  	conv1_pp = (upperlimit_pp-lowerlimit_pp)/2.0
	conv2_pp = (upperlimit_pp+lowerlimit_pp)/2.0
*
*  now get the Gdif matrix direction by direction
*
	DO i = 1, ng
	   muprime = xg(i)
	   DO j = 1, ng
	      phiprime = conv1_pp*xg(j) + conv2_pp
	      CALL G_dir_function (ng,xg,wg,gL,hL,muprime,phiprime,
     *                             Gdif(j,i))
	   END DO
	END DO
*
*  check for normalization
*
*  (1/4PI) int_0^2PI dphi^prime int_1^1 dmu^prime G(OMEGA^prime) = 0.5
*
	sum_tp = 0.0
	DO i = 1, ng
           sum_pp = 0.0
	   DO j = 1, ng
	      sum_pp = sum_pp + wg(j)*Gdif(j,i)
	   END DO
	   sum_pp = sum_pp*conv1_pp
	   sum_tp = sum_tp + wg(i)*sum_pp
	END DO
	sum_tp = sum_tp/(4.0*PI)
c	WRITE (*,*) " Gfun check (=0.5?): ", sum_tp
*
* done
* 
	RETURN
	END
*------------------------------------------------------------------------------*




*------------------------------------------------------------------------------*
*
	SUBROUTINE GAMMA_d_dir_function (ng,xg,wg,gL,hL,rho_Ld,tau_Ld, 
     *                                   muprime,phiprime,Gamma_d_dir)
*
* this routine evaluates the GAMMA_d function (muprime,phiprime ->
* mu,phi) where (mu,phi) are quadrature directions  and checks for 
* normalization
* 
* begin declarations
*
	PARAMETER (PI = 3.141592654)
*
        REAL xg(ng), wg(ng)
	REAL gL(ng), hL(ng)
        REAL muprime, phiprime, rho_Ld, tau_Ld
	REAL upperlimit_pp, lowerlimit_pp, conv1_pp, conv2_pp
	REAL mu, phi, Gamma_d_dir(ng,ng)
*
	INTEGER i, j
*
* end declarations
*
*  conversion factors to have the ordinates simulate phiprime
*
	upperlimit_pp = 2.0*PI
	lowerlimit_pp = 0.0
  	conv1_pp = (upperlimit_pp-lowerlimit_pp)/2.0
	conv2_pp = (upperlimit_pp+lowerlimit_pp)/2.0
*
*  now get the Gamma_d_dir matrix direction by direction
*
	DO i = 1, ng
	   mu = xg(i)
	   DO j = 1, ng
	      phi = conv1_pp*xg(j) + conv2_pp
	      CALL GAMMA_d_function (ng,xg,wg,gL,hL,rho_Ld,tau_Ld, 
     *                               muprime,phiprime,mu,phi,
     *                               Gamma_d_dir(j,i))
	   END DO
	END DO
*
* done
* 
	RETURN
	END
*------------------------------------------------------------------------------*


*------------------------------------------------------------------------------*
*
	SUBROUTINE GAMMA_d_dif_function (ng,xg,wg,gL,hL,rho_Ld,tau_Ld, 
     *                                   Gdif,Gamma_d_dif,dummy)
*
* this routine evaluates the Gamma_d_dif function for scattering
* from all quadrature directions to all exit quadrature directions 
* 
* begin declarations
*
	PARAMETER (PI = 3.141592654)
*
        REAL xg(ng), wg(ng), gL(ng), hL(ng), rho_Ld, tau_Ld
        REAL phiprime, muprime, Gdif(ng,ng)
	REAL upperlimit_pp, lowerlimit_pp, conv1_pp, conv2_pp
	REAL Gamma_d_dif(ng,ng,ng,ng), dummy(ng,ng)
*
	INTEGER ng, i, j, n, m
*
* end declarations
*
*  conversion factors to have the ordinates simulate phiprime
*
	upperlimit_pp = 2.0*PI
	lowerlimit_pp = 0.0
  	conv1_pp = (upperlimit_pp-lowerlimit_pp)/2.0
	conv2_pp = (upperlimit_pp+lowerlimit_pp)/2.0
*
*  now get the Gamma_d_dif matrix direction by direction
*
	DO i = 1, ng
	   muprime = xg(i)
	   DO j = 1, ng
	      phiprime = conv1_pp*xg(j) + conv2_pp
              CALL GAMMA_d_dir_function (ng,xg,wg,gL,hL,
     *                                   rho_Ld,tau_Ld, 
     *                                   muprime,phiprime,
     *                                   dummy)
*
	      CALL CHECK_Gamma_d_dir (ng,xg,wg,rho_Ld,tau_Ld,
     *                                Gdif(j,i),dummy)

	      DO m = 1, ng
	      DO n = 1, ng
	         Gamma_d_dif(j,i,n,m) = dummy(n,m)
	      END DO
	      END DO
	   END DO
	END DO
*
* done
* 
	RETURN
	END
*------------------------------------------------------------------------------*


*------------------------------------------------------------------------------*
*
	SUBROUTINE GAMMA_d_function (ng,xg,wg,gL,hL,rho_Ld,tau_Ld, 
     *                               muprime,phiprime,mu,phi,
     *                               Gamma_d)
* 
* this routine calculates the G_FUNCTION given a direction of
* photon travel OMEGA^PRIME and the leaf normal pdf
* 
* begin declarations
*
	PARAMETER (PI = 3.141592654)
*
        REAL    xg(ng), wg(ng), gL(ng), hL(ng), rho_Ld, tau_Ld
        REAL    muprime, phiprime, mu, phi
	REAL    mu_t, mu_tp, sin_t, sin_tp
	REAL    upperlimit_tL, lowerlimit_tL
        REAL    conv1_tL, conv2_tL, sum_tL
	REAL    upperlimit_pL, lowerlimit_pL
  	REAL    conv1_pL, conv2_pL, sum_pL
 	REAL    neword_tL, mu_tL, sin_tL
	REAL    neword_pL, dotproduct1, dotproduct2, Gamma_d
*
	INTEGER i, j
*
* end declarations
*
	mu_t   = mu
	mu_tp  = muprime
	sin_t  = sqrt(1.0 - mu*mu)
	sin_tp = sqrt(1.0 - muprime*muprime)
*
* define limits of integration and the convertion factors for integration
* over thetaL (note the tL suffix!)
*
	upperlimit_tL = PI/2.0
	lowerlimit_tL = 0.0
  	conv1_tL = (upperlimit_tL-lowerlimit_tL)/2.0
	conv2_tL = (upperlimit_tL+lowerlimit_tL)/2.0
*
* define limits of integration and the convertion factors for integration
* over phiL (note the pL suffix!)
*
	upperlimit_pL = 2.0*PI
	lowerlimit_pl = 0.0
  	conv1_pL = (upperlimit_pL-lowerlimit_pL)/2.0
	conv2_pL = (upperlimit_pL+lowerlimit_pL)/2.0
*
* integral over theta_L
*
	sum_tL = 0.0
	do i = 1, ng
	   neword_tL = conv1_tL*xg(i) + conv2_tL
	   mu_tL     = cos(neword_tL)
	   sin_tL    = sin(neword_tL)
*
* integral over phi_L
*
	   sum_pL = 0.0
	   do j = 1, ng
	      neword_pL  = conv1_pL*xg(j) + conv2_pL
              dotproduct1 = ( mu_tL*mu_tp + 
     *                        sin_tL*sin_tp*cos(neword_pL-phiprime) )
              dotproduct2 = ( mu_tL*mu_t  + 
     *                        sin_tL*sin_t *cos(neword_pL-phi     ) )
	      If (dotproduct1*dotproduct2 .LE. 0.0) then
                  sum_pL = sum_pL + rho_Ld*wg(j)*hL(j)/(2.0*PI)*
     *                              abs(dotproduct1*dotproduct2)
	      else
                  sum_pL = sum_pL + tau_Ld*wg(j)*hL(j)/(2.0*PI)*
     *                              abs(dotproduct1*dotproduct2)
	      endif
           end do
* 
* finish the phi_L integral
*
	   sum_pL = sum_pL*conv1_pL
           sum_tL = sum_tL + wg(i)*gL(i)*sum_pL
	end do
* 
* finish the theta_L integral
*
	sum_tL   = sum_tL*conv1_tL
	Gamma_d  = sum_tL
*
* done
* 
	RETURN
	END
*------------------------------------------------------------------------------*


*------------------------------------------------------------------------------*
*
	SUBROUTINE CHECK_Gamma_d_dir (ng,xg,wg,rho_Ld,tau_Ld,
     *                                Gdir,Gamma_d_dir)

*
* this routine checks the Gamma_d_dir for normalization 
* 
* begin declarations
*
	PARAMETER (PI = 3.141592654)
*
        REAL xg(ng), wg(ng), rho_Ld, tau_Ld, Gdir, Gamma_d_dir(ng,ng)
	REAL upperlimit_pp, lowerlimit_pp, conv1_pp
	REAL sum_tp, sum_pp
*
	INTEGER i, j
*
* end declarations
*
*  conversion factors to have the ordinates simulate phiprime
*
	upperlimit_pp = 2.0*PI
	lowerlimit_pp = 0.0
  	conv1_pp = (upperlimit_pp-lowerlimit_pp)/2.0
*
*  check for normalization
*
*  (1/PI) int_0^2PI dphi int_(-1)^1 dmu 
*          Gamma_d_dir(muprime,phiprime -> mu,phi) = 
*           Gdir*(rho_Ld+tau_Ld)    
*
	sum_tp = 0.0
	DO i = 1, ng
           sum_pp = 0.0
	   DO j = 1, ng
	      sum_pp = sum_pp + wg(j)*Gamma_d_dir(j,i)
	   END DO
	   sum_pp = sum_pp*conv1_pp
	   sum_tp = sum_tp + wg(i)*sum_pp
	END DO
	sum_tp = sum_tp/(PI)
	sum_tp = sum_tp/(Gdir*(rho_Ld+tau_Ld))
c	WRITE (*,*) " Gamma_d chk(=1.0?): ", sum_tp
*
* done
* 
	RETURN
	END
*------------------------------------------------------------------------------*


*------------------------------------------------------------------------------*
*
	SUBROUTINE I_o_uncol_down (DeltaL,Nlayers,Gdir,I_o,mu_o, 
     *                             I_o_uc_d)
*
* this routine evaluates the dpwnward uncollided direct solar radiation  
* 
* begin declarations
*
	PARAMETER (PI = 3.141592654)
*
	REAL DeltaL, Gdir, I_o, mu_o
	REAL L1, L2, Prob, I_o_uc_d(Nlayers)
*
	INTEGER Nlayers, k
*
* end declarations
*
*
* evaluate downward uncollided direct solar intensity layer by layer
*
	L1 = 0.0
	do k = 1, Nlayers
	   L2 = float(k)*DeltaL - (0.5*DeltaL)
	   Prob = exp(- (1/abs(mu_o)) * Gdir * (L2-L1) )
           I_o_uc_d(k) = I_o * Prob
	end do
* 
* done
*
        RETURN
        END
*------------------------------------------------------------------------------*


*------------------------------------------------------------------------------*
*
	SUBROUTINE I_o_uncol_up (DeltaL,Nlayers,Gdir,Gdif,I_o,mu_o,
     *                           ng,xg,R_s,F_o_uc_d_soil,I_o_uc_u)
*
* this routine evaluates the upward uncollided direct solar radiation  
* 
* begin declarations
*
	PARAMETER (PI = 3.141592654)
*
	REAL DeltaL, Gdir, Gdif(ng,ng), I_o
	REAL mu_o, R_s, xg(ng), L1
	REAL L2, Prob, I_o_uc_u_soil, F_o_uc_d_soil
        REAL I_o_uc_u(Nlayers,ng,ng)
*
	INTEGER Nlayers, ng, i, j, k
*
* end declarations
*
*
* the uncollided direct solar radiation incident on the ground
*
	L1 = 0.0
	L2 = Nlayers*DeltaL
	F_o_uc_d_soil = abs(mu_o)*I_o*exp(-(1/abs(mu_o))*Gdir*(L2-L1))
*
* upward uncollided intensity from reflection by soil of I_o_uc_d_soil
*
	I_o_uc_u_soil = (R_s/PI)*F_o_uc_d_soil
*
* evaluate upward uncollided direct solar intensity layer by layer in all
* upward directions
*
	do i = (ng/2)+1, ng
	do j = 1, ng
	   L2 = Nlayers*DeltaL
	   do k = Nlayers, 1, -1
	      L1 = float(k)*DeltaL - (0.5*DeltaL)
	      Prob = exp(- (1/abs(xg(i))) * Gdif(j,i) * (L2-L1) )
              I_o_uc_u(k,j,i) = I_o_uc_u_soil * Prob
	   end do
	end do
	end do
* 
* done
*
        RETURN
        END
*------------------------------------------------------------------------------*


*------------------------------------------------------------------------------*
*
	SUBROUTINE I_d_uncol_down (DeltaL,Nlayers,Gdif,I_d,ng,xg,
     *                             I_d_uc_d)
*
* this routine evaluates the downward uncollided diffuse sky radiation  
* 
* begin declarations
*
	PARAMETER (PI = 3.141592654)
*
	REAL DeltaL, Gdif(ng,ng), I_d, xg(ng)
	REAL L1, L2, Prob
        REAL I_d_uc_d(Nlayers,ng,ng)
*
	INTEGER Nlayers, ng, i, j, k
*
* end declarations
*
* evaluate downward uncollided diffuse sky intensity layer by layer in all
* downward directions
*
	do i = 1, ng/2
	do j = 1, ng
	   L1 = 0.0
	   do k = 1, Nlayers
	      L2 = float(k)*DeltaL - (0.5*DeltaL)
	      Prob = exp(- (1/abs(xg(i))) * Gdif(j,i) * (L2-L1) )
              I_d_uc_d(k,j,i) = I_d * Prob
	   end do
	end do
	end do
* 
* done
*
        RETURN
        END
*------------------------------------------------------------------------------*


*------------------------------------------------------------------------------*
*
	SUBROUTINE I_d_uncol_up (DeltaL,Nlayers,Gdif,I_d,ng,xg,wg,R_s,
     *                           F_d_uc_d_soil,I_d_uc_u)
*
* this routine evaluates the upward uncollided diffuse sky radiation  
* 
* begin declarations
*
	PARAMETER (PI = 3.141592654)
*
	REAL DeltaL, Gdif(ng,ng), I_d, xg(ng), wg(ng), R_s
	REAL upperlimit, lowerlimit, conv, sum1, sum2
	REAL L1, L2, Prob
        REAL I_d_uc_d_soil, I_d_uc_u_soil, F_d_uc_d_soil
 	REAL I_d_uc_u(Nlayers,ng,ng)
*
	INTEGER Nlayers, ng, i, j, k
*
* end declarations
*
*
* evaluate diffuse sky flux density incident on the soil below canopy
*
	L1 = 0.0
	L2 = Nlayers*DeltaL
	sum1 = 0.0
	do i = 1, ng/2
	   upperlimit = 2.0*PI
	   lowerlimit = 0.0
  	   conv = (upperlimit-lowerlimit)/2.0
	   sum2 = 0.0
	   do j = 1, ng
	      Prob = exp(- (1/abs(xg(i))) * Gdif(j,i) * (L2-L1) )
	      I_d_uc_d_soil = I_d * Prob
	      sum2 = sum2 + wg(j)*I_d_uc_d_soil
	   end do
           sum1 = sum1 + wg(i)*abs(xg(i))*sum2*conv
	end do
	F_d_uc_d_soil = sum1
*
* evaluate upward uncollided intensity due to reflection from soil
*
	I_d_uc_u_soil = F_d_uc_d_soil*(R_s/PI)
*
* evaluate downward uncollided diffuse sky intensity layer by layer in all
* downward directions
*
	do i = (ng/2)+1, ng
	do j = 1, ng
	   L2 = Nlayers*DeltaL
	   do k = Nlayers, 1, -1
	      L1 = float(k)*DeltaL - (0.5*DeltaL)
	      Prob = exp(- (1/abs(xg(i))) * Gdif(j,i) * (L2-L1) )
              I_d_uc_u(k,j,i) = I_d_uc_u_soil * Prob
	   end do
	end do
	end do
* 
* done
*
        RETURN
        END
*------------------------------------------------------------------------------*



*------------------------------------------------------------------------------*
*
	SUBROUTINE FCS (Nlayers,ng,xg,wg,Gamma_d_dir,Gamma_d_dif,
     *                  I_o_uc_d, I_o_uc_u,I_d_uc_d, I_d_uc_u,
     *                  Q)
*
* this routine evaluates the first-collision source Q
* 
* begin declarations
*
	PARAMETER (PI = 3.141592654)
*
	REAL xg(ng), wg(ng)
	REAL Gamma_d_dir(ng,ng),Gamma_d_dif(ng,ng,ng,ng)
	REAL I_o_uc_d(Nlayers), I_o_uc_u(Nlayers,ng,ng)
	REAL I_d_uc_d(Nlayers,ng,ng), I_d_uc_u(Nlayers,ng,ng)
	REAL upperlimit1, lowerlimit1, conv11, sum1
	REAL upperlimit2, lowerlimit2, conv21, sum2
 	REAL Q(Nlayers,ng,ng)
*
	INTEGER Nlayers, ng, i, j, k
*
* end declarations
*
*
* evaluate first-collision source due to I_o_uc_d
*
*
	do i = 1, ng
	do j = 1, ng
	do k = 1, Nlayers
	   Q(k,j,i) = (1.0/PI)*Gamma_d_dir(j,i)*I_o_uc_d(k)
	end do
	end do
	end do
*
* evaluate first-collision source due to I_o_uc_u, I_d_uc_d, I_d_uc_u
*
*
	do i = 1, ng
	do j = 1, ng
	do k = 1, Nlayers
*
	   upperlimit1 =  1.0
	   lowerlimit1 = -1.0
  	   conv11 = (upperlimit1-lowerlimit1)/2.0
	   sum1 = 0.0
	   do n = 1, ng
	      upperlimit2 = 2.0*PI
	      lowerlimit2 = 0.0
  	      conv21 = (upperlimit2-lowerlimit2)/2.0
	      sum2 = 0.0
	      do m = 1, ng
	         sum2 = sum2 + wg(m)*(1.0/PI)*Gamma_d_dif(m,n,j,i)*
     *           (I_o_uc_u(k,m,n)+I_d_uc_d(k,m,n)+I_d_uc_u(k,m,n))
	      end do
              sum1 = sum1 + wg(n)*sum2*conv21
	   end do
	   Q(k,j,i) = Q(k,j,i) + sum1*conv11
*
	end do
	end do
	end do
* 
* done
*
        RETURN
        END
*------------------------------------------------------------------------------*


*------------------------------------------------------------------------------*
*
	SUBROUTINE SWEEP_DOWN (Nlayers,ng,xg,wg,Gdif,
     *                         DeltaL,JJ,R_s,
     *                         Ic)
*
* this routine sweeps downwards in the phase-space mesh and handle the
* bottom boundary condition and evaluate the upward Ic at the ground
* 
* begin declarations
*
	PARAMETER (PI = 3.141592654)
*
	REAL xg(ng), wg(ng)
	REAL Gdif(ng,ng), DeltaL, JJ(Nlayers,ng,ng), R_s
	REAL fij, aij, bij
	REAL Ic(Nlayers+1,ng,ng)
	REAL upperlimit, lowerlimit, conv, sum1, sum2
	REAL Fc_soil
*
	INTEGER Nlayers, ng, i, j, k
*
*
* end declarations
*
*
* sweep downwards
*
	do i = 1, ng/2
	do j = 1, ng
	   fij = (xg(i)/DeltaL) - (0.5*Gdif(j,i))
           aij = ((0.5*Gdif(j,i)) + (xg(i)/DeltaL)) / fij
	   bij = 1.0/fij
	   do k = 1, Nlayers
	      Ic(k+1,j,i) = aij*Ic(k,j,i) - bij*JJ(k,j,i)
	   end do
	end do
	end do
*
* evaluate flux density incident on the ground
*
	sum1 = 0.0
	do i = 1, ng/2
	   upperlimit = 2.0*PI
	   lowerlimit = 0.0
  	   conv = (upperlimit-lowerlimit)/2.0
	   sum2 = 0.0
	   do j = 1, ng
	      sum2 = sum2 + wg(j)*Ic(Nlayers+1,j,i)
	   end do
           sum1 = sum1 + wg(i)*abs(xg(i))*sum2*conv
	end do
	Fc_soil = sum1
*
* evluate Ic upward at the ground
*
	do i = (ng/2)+1, ng
	do j = 1, ng
	   Ic(Nlayers+1,j,i) = (Fc_soil*R_s)/PI
	end do
	end do
* 
* done
*
        RETURN
        END
*------------------------------------------------------------------------------*



*------------------------------------------------------------------------------*
*
	SUBROUTINE SWEEP_UP (Nlayers,ng,xg,wg,Gdif,DeltaL,JJ,
     *                       Ic,Ic_old,epsilon,convergence)
*
* this routine sweeps upwards in the phase-space mesh and checks for
* convergence
* 
* begin declarations
*
	PARAMETER (PI = 3.141592654)
*
	REAL xg(ng), wg(ng), Gdif(ng,ng), DeltaL, JJ(Nlayers,ng,ng)
	REAL fij, cij, dij
	REAL Ic(Nlayers+1,ng,ng), epsilon, Ic_old(ng,ng)
*
	INTEGER Nlayers, ng, i, j, k
*
	LOGICAL convergence
*
*
* end declarations
*
*
* save earlier iterate of Ic(1,j,i) for convergence checking
*
	do i = (ng/2)+1, ng
	do j = 1, ng
	   Ic_old(j,i) = Ic(1,j,i)
	end do
	end do
*
* sweep upwards
*
	do i = (ng/2)+1, ng
	do j = 1, ng
	   fij = ((xg(i)/DeltaL)+(0.5*Gdif(j,i)))
           cij = ((xg(i)/DeltaL)-(0.5*Gdif(j,i))) / fij
	   dij = 1.0/fij
	   do k = Nlayers, 1, -1
	      Ic(k,j,i) = cij*Ic(k+1,j,i) + dij*JJ(k,j,i)
	   end do
	end do
	end do
*
* check for convergence
*
	convergence = .true.
	do i = (ng/2)+1, ng
	do j = 1, ng
           if (abs(Ic(1,j,i)-Ic_old(j,i)) .gt. epsilon) then
	        convergence = .false.
                go to 100
	   endif
	end do
	end do
* 
* done
*
100     RETURN
        END
*------------------------------------------------------------------------------*

*------------------------------------------------------------------------------*
*
	SUBROUTINE MULTI_COLL_S (Nlayers,ng,xg,wg,Gamma_d_dif,Ic,
     *                           S)
*
* this routine evaluates the multiple-collision source S
* 
* begin declarations
*
	PARAMETER (PI = 3.141592654)
*
	REAL xg(ng), wg(ng)
	REAL Gamma_d_dif(ng,ng,ng,ng)
	REAL upperlimit1, lowerlimit1, conv11, sum1
	REAL upperlimit2, lowerlimit2, conv21, sum2
 	REAL Ic(Nlayers+1,ng,ng), Ic_cell_center
 	REAL S(Nlayers,ng,ng)
*
	INTEGER Nlayers, ng, i, j, k, m, n
*
* end declarations
*
*
* evaluate multiple-collision source due to Ic
*
	do i = 1, ng
	do j = 1, ng
	do k = 1, Nlayers
*
	   upperlimit1 =  1.0
	   lowerlimit1 = -1.0
  	   conv11 = (upperlimit1-lowerlimit1)/2.0
	   sum1 = 0.0
	   do n = 1, ng
	      upperlimit2 = 2.0*PI
	      lowerlimit2 = 0.0
  	      conv21 = (upperlimit2-lowerlimit2)/2.0
	      sum2 = 0.0
	      do m = 1, ng
	         Ic_cell_center = 0.5*(Ic(k,m,n)+Ic(k+1,m,n))
	         sum2 = sum2 + wg(m)*(1.0/PI)*Gamma_d_dif(m,n,j,i)*
     *                         Ic_cell_center
	      end do
              sum1 = sum1 + wg(n)*sum2*conv21
	   end do
	   S(k,j,i) =  sum1*conv11
*
	end do
	end do
	end do
* 
* done
*
        RETURN
        END
*------------------------------------------------------------------------------*

*------------------------------------------------------------------------------*
*
	SUBROUTINE ENERGY_BAL (Nlayers,ng,xg,wg,mu_o,Q,S,DeltaL,R_s,
     *                         rho_LD,tau_LD,Gdir,Gdif,
     *                         F_o_uc_d_soil,F_d_uc_d_soil,
     *                         I_o_uc_d,I_d_uc_d,
     *                         I_o_uc_u,I_d_uc_u,
     *                         Ic,
     *                         HR_uc,HR_c,HT_uc,HT_c,AB_uc,AB_c)

*
* this routine does energy balance
* 
* begin declarations
*
	PARAMETER (PI = 3.141592654)
*
	REAL xg(ng), wg(ng), mu_o
	REAL rho_LD, tau_LD, Gdif(ng,ng), DeltaL, R_s
	REAL Q(Nlayers,ng,ng), S(Nlayers,ng,ng)
	REAL L1, L2, Prob, I_uc_u
	REAL I_o_uc_d(Nlayers), I_d_uc_d(Nlayers,ng,ng)
	REAL I_o_uc_u(Nlayers,ng,ng), I_d_uc_u(Nlayers,ng,ng)
	REAL F_o_uc_d_soil, F_d_uc_d_soil
 	REAL Ic(Nlayers+1,ng,ng)
	REAL upperlimit, lowerlimit, conv, sum1, sum2
	REAL HR_uc, HR_c, HT_uc, HT_c
	REAL AB_c,AB_uc,AB_o_uc_d,AB_o_uc_u,AB_d_uc_d,AB_d_uc_u

*
	INTEGER Nlayers, ng, i, j
*
*
* end declarations
*
	upperlimit = 2.0*PI
	lowerlimit = 0.0
  	conv = (upperlimit-lowerlimit)/2.0
*
* evaluate uncollided hemispherical transmittance
*
	HT_uc = F_o_uc_d_soil + F_d_uc_d_soil
*
* evaluate collided hemispherical transmittance
*
	sum1 = 0.0
	do i = 1, ng/2
	   sum2 = 0.0
	   do j = 1, ng
	      sum2 = sum2 + wg(j)*Ic(Nlayers+1,j,i)
	   end do
           sum1 = sum1 + wg(i)*abs(xg(i))*sum2*conv
	end do
	HT_c = sum1
*
* evaluate uncollided hemispherical reflectance
*
	L1 = 0.0
	L2 = Nlayers*DeltaL
	sum1 = 0.0
	do i = (ng/2)+1, ng
	   sum2 = 0.0
	   do j = 1, ng
	      Prob = exp(-(1/abs(xg(i)))*Gdif(j,i)*(L2-L1))
              I_uc_u = ((F_o_uc_d_soil+F_d_uc_d_soil)*R_s/PI)*Prob
	      sum2 = sum2 + wg(j)*I_uc_u
	   end do
           sum1 = sum1 + wg(i)*abs(xg(i))*sum2*conv
	end do
	HR_uc = sum1
*
* evaluate collided hemispherical reflectance
*
	sum1 = 0.0
	do i = (ng/2)+1, ng
	   sum2 = 0.0
	   do j = 1, ng
	      sum2 = sum2 + wg(j)*Ic(1,j,i)
	   end do
           sum1 = sum1 + wg(i)*abs(xg(i))*sum2*conv
	end do
	HR_c = sum1
*
* evaluate canopy absorption from I_o_uc_d
*
	sum1 = 0.0
	do k = 1, Nlayers
	   sum1 = sum1 + I_o_uc_d(k)*Gdir
	end do
	AB_o_uc_d = sum1*(1.0-(rho_LD+tau_LD))*DeltaL
*
* evaluate canopy absorption from I_o_uc_u
*
	sum1 = 0.0
	do k = 1, Nlayers
	do i = (ng/2)+1, ng
	   sum2 = 0.0
	   do j = 1, ng
	      sum2 = sum2 + wg(j)*I_o_uc_u(k,j,i)*Gdif(j,i)
	   end do
	   sum 1 = sum1 + wg(i)*sum2*conv
	end do
	end do
	AB_o_uc_u = sum1*(1.0-(rho_LD+tau_LD))*DeltaL
*
* evaluate canopy absorption from I_d_uc_d
*
	sum1 = 0.0
	do k = 1, Nlayers
	do i = 1, ng/2
	   sum2 = 0.0
	   do j = 1, ng
	      sum2 = sum2 + wg(j)*I_d_uc_d(k,j,i)*Gdif(j,i)
	   end do
	   sum 1 = sum1 + wg(i)*sum2*conv
	end do
	end do
	AB_d_uc_d = sum1*(1.0-(rho_LD+tau_LD))*DeltaL
*
* evaluate canopy absorption from I_d_uc_u
*
	sum1 = 0.0
	do k = 1, Nlayers
	do i = (ng/2)+1, ng
	   sum2 = 0.0
	   do j = 1, ng
	      sum2 = sum2 + wg(j)*I_d_uc_u(k,j,i)*Gdif(j,i)
	   end do
	   sum 1 = sum1 + wg(i)*sum2*conv
	end do
	end do
	AB_d_uc_u = sum1*(1.0-(rho_LD+tau_LD))*DeltaL
*
* evaluate canopy absorption from I_c
*
	sum1 = 0.0
	do k = 1, Nlayers
	do i = 1, ng
	   sum2 = 0.0
	   do j = 1, ng
	      sum2 = sum2 + wg(j)*Ic(k,j,i)*Gdif(j,i)
	   end do
	   sum 1 = sum1 + wg(i)*sum2*conv
	end do
	end do
	AB_c = sum1*(1.0-(rho_LD+tau_LD))*DeltaL
*
	AB_uc = AB_o_uc_d + AB_o_uc_u + AB_d_uc_d + AB_d_uc_u
* 
* done
*
        RETURN
        END
*------------------------------------------------------------------------------*

