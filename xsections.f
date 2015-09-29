	PROGRAM XSECTIONS
*
* this simple program illustrates 
*
* (1) obtain quadrature [ng,xg,wg]
* (2) obtain leaf normal orientation pdfs [gL,hL]
* (3) evaluate the G function [G(mu,phi)]
* (4) evaluate the Gamma_d function [Gamma_d(phip,mup->phi,mu]
*
* 
* begin declarations
*
	PARAMETER (PI = 3.141592654)
        PARAMETER (ng = 4)
	PARAMETER (degtorad = 3.141592654/180.0)
	PARAMETER (radtodeg = 180.0/3.141592654)
*
        REAL xg(ng), wg(ng), gL(ng), hL(ng), rho_Ld, tau_Ld
        REAL thetaprime, phiprime, muprime, Gdir, Gdif(ng,ng)
	REAL Gamma_d_dir(ng,ng), Gamma_d_dif(ng,ng,ng,ng)
*
*
* end declarations
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
* get G_FUNCTION for a direction OMEGA^prime
*
	thetaprime = 120.0
	phiprime   = 0.0
	thetaprime = thetaprime*degtorad
	phiprime   = phiprime*degtorad
	muprime    = cos(thetaprime)
*
	CALL G_dir_function (ng,xg,wg,gL,hL,muprime,phiprime,Gdir)
*
* get G_FUNCTION for all quadrature directions
*
*  Gdif(phi,mu)
* 
       CALL G_dif_function (ng,xg,wg,gL,hL,Gdif)
*
* get GAMMA_d (phiprime, muprime -> phi, mu) where (mu,phi) are all
* quadrature directions
*
*  Gamma_d_dir(phi,mu)
*
	rho_Ld = 0.05
	tau_Ld = 0.25
*
	CALL GAMMA_d_dir_function (ng,xg,wg,gL,hL,rho_Ld,tau_Ld, 
     *                             muprime,phiprime,Gamma_d_dir) 
*
	CALL CHECK_Gamma_d_dir (ng,xg,wg,rho_Ld,tau_Ld,
     *                          Gdir,Gamma_d_dir)

*
* get GAMMA_d (muprime,phiprime -> mu,phi) where both the incident
* and exit directions are all quadrature directions
*
*  Gamma_d_dif(phiprime,muprime,phi,mu)
*
	CALL GAMMA_d_dif_function (ng,xg,wg,gL,hL,rho_Ld,tau_Ld, 
     *                             Gdif,Gamma_d_dif)
*
* 
	STOP
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
	WRITE (*,*) " Qwts check (=2.0?): ", sum
*
* check if [ int_0^1 dx x ] is equal to 0.5
*
	sum = 0.0
	DO i = (ng/2)+1, ng
	   sum = sum + (xg(i)*wg(i))
	END DO
	WRITE (*,*) " Qord check (=0.5?): ", sum
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
	WRITE (*,*) " Intg check (=1.0?): ", sum
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
	WRITE (*,*) " LNO  check (=1.0?): ", sum
* 
* done
*
        RETURN
        END
*------------------------------------------------------------------------------*


*------------------------------------------------------------------------------*
*
	SUBROUTINE G_dir_ function (ng,xg,wg,gL,hL,
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
	      CALL G_dir_ function (ng,xg,wg,gL,hL,muprime,phiprime,
     *                              Gdif(j,i))
	   END DO
	END DO
*
*  check for normalization
*
*  (1/2PI) int_0^2PI dphi^prime int_0^1 dmu^prime G(OMEGA^prime) = 0.5
*
	sum_tp = 0.0
	DO i = (ng*0.5)+1, ng
           sum_pp = 0.0
	   DO j = 1, ng
	      sum_pp = sum_pp + wg(j)*Gdif(j,i)
	   END DO
	   sum_pp = sum_pp*conv1_pp
	   sum_tp = sum_tp + wg(i)*sum_pp
	END DO
	sum_tp = sum_tp/(2.0*PI)
	WRITE (*,*) " Gfun check (=0.5?): ", sum_tp
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
	REAL upperlimit_p, lowerlimit_p, conv1_p, conv2_p
	REAL mu, phi, Gamma_d_dir(ng,ng)
*
	INTEGER i, j
*
* end declarations
*
*  conversion factors to have the ordinates simulate phi
*
	upperlimit_p = 2.0*PI
	lowerlimit_p = 0.0
  	conv1_p = (upperlimit_p-lowerlimit_p)/2.0
	conv2_p = (upperlimit_p+lowerlimit_p)/2.0
*
*  now get the Gamma_d_dir matrix direction by direction
*
	DO i = 1, ng
	   mu = xg(i)
	   DO j = 1, ng
	      phi = conv1_p*xg(j) + conv2_p
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
     *                                   Gdif,Gamma_d_dif)
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
	REAL Gamma_d_dif(ng,ng,ng,ng), dummy(12,12)
*
	INTEGER i, j, n, m
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
	REAL upperlimit_pp, lowerlimit_pp, conv1_pp, conv2_pp
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
	conv2_pp = (upperlimit_pp+lowerlimit_pp)/2.0
*
*  check for normalization
*
*  (1/PI) int_0^2PI dphi int_1^1 dmu 
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
	WRITE (*,*) " Gamma_d chk(=1.0?): ", sum_tp
*
* done
* 
	RETURN
	END
*------------------------------------------------------------------------------*

