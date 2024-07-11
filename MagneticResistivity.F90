!!!!****if* source/physics/materialProperties/MagneticResistivity/MagneticResistivityMain/AmbD/MagneticResistivity_fullState
!!!!
!!!! NAME
!!!!  MagneticResistivity_fullState
!!!!
!!!! SYNOPSIS
!!!!  call MagneticResistivity_fullState(real(in)    :: solnVec(NUNK_VARS),
!!!!                                     real(out)   :: eta_paral,
!!!!                                     real(out)   :: eta_perp,
!!!!                                     real(out)   :: eta_hall)
!!!! DESCRIPTION
!!!!
!!!! Computes the conductivities for all species (ions & electrons if
!!!! the inelastic_flag is false and all species - i, e and grains if
!!!! if inelastic flag is true). From conductivities --> resistivities 
!!!! as Kunz & Mouschovias 2009 did it
!!!! 
!!!!
!!!!     eta_v = eta_v_ad + eta_v_ohm =>                   
!!!!                                                       
!!!!                    sigma_v                            
!!!!     eta_v = ------------------                        
!!!!             sigma_v^2+sigma_h^2                       
!!!!                                                       
!!!!                      1                                
!!!!     eta_//= -------------------                       
!!!!                  sigma_//                             
!!!!
!!!!                  sigma_h                              
!!!!     eta_h = -------------------                       
!!!!             sigma_v^2+sigma_h^2                       
!!!!                                                       
!!!!  The rest of the resistivities may be computed from   
!!!!      these relations:                                 
!!!!                                                       
!!!!    i.) eta_//_ad=0.        ii.) eta_//_ohm = eta_//   
!!!!  iii.) eta_v_ohm=eta_//    iv.) eta_v_ad=eta_v-eta_// 
!!!!                                                       
!!!!
!!!!  Returns Magnetic Resistivity, parallel, perpendicular and hall components.
!!!!
!!!! ARGUMENTS
!!!!
!!!!   solnVec  :   solution state, a vector from UNK with all variables
!!!!   eta_paral:   parallel component of Magnetic Resistivity
!!!!   eta_perp :   perpendicular component of Magnetic Resistivity
!!!!   eta_hall :   hall component of Magnetic Resistivity
!!!!
!!!!***

#include "Flash.h"
#include "constants.h" 


subroutine MagneticResistivity(temp,dens, magx, magy, magz, curx, cury, curz, vise, crze, xn, &
                               eta_paral, eta_hall, eta_perp, sigma_paral, sigma_hall, sigma_perp)

  use MagneticResistivity_data
  use Driver_data , ONLY : dr_simTime
 
  implicit none
  
#include "Flash_mpi.h"

  !!-------------------------------------------------------!!
  !!      what I want to exit this piece of code. Either   !!
  !!            this is a matrix or a single number        !!
  real, intent(IN) :: temp, dens, magx, magy, magz, curx, cury, curz, vise, crze
  real, intent(IN), dimension(NSPECIES) :: xn              !!
  real, intent(OUT):: eta_paral, eta_hall, eta_perp, sigma_paral, sigma_hall, sigma_perp
                                                           !!
  !!-------------------------------------------------------!!

  !!-------------------------------------------------------!!
  !!   More parameters that I will definitely use and are  !!
  !!                   defined here                        !!
  integer :: nsp, nsp2, gnb, cnt, exitCounter              !!
                                                           !!
  real :: freeelectrons_magres, gdratio                    !!
  real :: mcr_grains, ionsH, elecH                         !!
  real :: etaPa, etaHa, etaPe                              !!      ** used for while loop **
  real :: mcri, mcre, mcrg, G0, ass, rdm                   !!
  real :: etaPaL, etaHaL, etaPeL, sigmaPaL, sigmaHaL, sigmaPeL  !  ** These are computed based on the Langevin approximation **
                                                           !!          for the actual eta_paral, eta_hall, eta_perp and   
  !!       Parameters coming elsewhere from the code       !!      corresponding sigmas drift velocities are taken into account
  real :: rho_magres, temperature_magres                   !!
  real :: Bx_magres, By_magres, Bz_magres, B_magres        !!
  real :: Jx, Jy, Jz                                       !!
  !!     Some more dummy parameters to define a grain      !!
  !!            distribution for the grains                !!
  real :: rhoG_tot, ng_tot, ng_loc, rad_loc, gmass_loc     !!
  !!    More parameters in case a simple resistivity       !!
  !!            calculation is requested                   !!
  real :: va, n_i_tot, numDens, rho_H2, res_constant!, f, g_H, g_He, g_HCO, g_H3O, k, l, m
  integer :: av_iproc, av_err, av_meshComm, av_numProcs
  real :: gamma, consty, Av, zeta
  !!-------------------------------------------------------!!

  !!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!!
  !!     All the conductivities/resistivities will be      !!
  !!                 Computed here                         !!
  !!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!!

  !!-------------------------------------------------------!!
  !! Take density magnetic field in all directions and     !!
  !!   abundances of all molecules for the specific gp     !!
  abund_magres(:) = 0. ; nmol_magres(:) = 0. ; cyclo_freq(:) = 0.
  sigma_indivS(:) = 0. ; coll_times(:) = 0. ; species_masses(:) = 0.
                                                           !!
  rho_magres = dens                                        !!
                                                           !!
  temperature_magres = temp                                !!
                                                           !!
  Bx_magres = magx * SQRT(4*PI)                            !!
  By_magres = magy * SQRT(4*PI)                            !!
  Bz_magres = magz * SQRT(4*PI)                            !!
                                                           !!
  B_magres=sqrt(Bx_magres**2+By_magres**2+Bz_magres**2)    !!
       
     rho_H2=2*rho_magres/res_mu
     va = B_magres/sqrt(rho_H2*4*PI)                       !!
     numDens = rho_magres/res_amu/res_mu
     
     
     call Driver_getComm(GLOBAL_COMM, av_meshComm)
                                                                 !!
     !n_i_tot = max(3.2e-3*numDens**(-0.5), 2.24e-03/((numDens/3.e+4)**(-0.45) * (1.+(numDens/3.e+4)**(-0.4)))) !LowAv
     !n_i_tot = max(2.8e-3*numDens**(-0.5), 2.23e-03/((numDens/3.e+4)**(-0.5) * (1.+(numDens/3.e+4)**(-0.4)))) ! Fidu/MF/Delayed
     !n_i_tot = max(1.93e-3*numDens**(-0.5), 1.53e-03/((numDens/3.e+4)**(-0.5) * (1.+(numDens/3.e+4)**(-0.4)))) !LowZ
     !n_i_tot = max(6.02e-3*numDens**(-0.5), 2.87e-03/((numDens/3.e+4)**(-0.2) * (1.+(numDens/3.e+4)**(-0.48)))) ! HighZ
                                                           !!
     
     if (saved_time .ne. dr_simTime) then
      saved_time=dr_simTime
      Call mpi_allreduce(0, allcounter, 1, FLASH_INTEGER, MPI_MIN, av_meshComm, av_err)
      Call mpi_allreduce(counter, allcounter, 1, FLASH_INTEGER, MPI_SUM, av_meshComm, av_err)
      
      temporary_logmean=temporary_logmean * counter / allcounter
      Call mpi_allreduce(0.0, mass_dens_logmean, 1, FLASH_REAL, MPI_MIN, av_meshComm, av_err)
      Call mpi_allreduce(temporary_logmean, mass_dens_logmean, 1, FLASH_REAL, MPI_SUM, av_meshComm, av_err)

      mass_dens_logmean=10**(mass_dens_logmean)
      
      Call mpi_allreduce(0.0, mass_dens_int, 1, FLASH_REAL, MPI_MIN, av_meshComm, av_err)
      Call mpi_allreduce(temporary_max, mass_dens_int, 1, FLASH_REAL, MPI_MAX, av_meshComm, av_err)
      
      mass_dens_red=(mass_dens_int * mass_dens_logmean)**(0.5)
      
      temporary_max=0
      temporary_logmean=0
      counter=0
     end if
     
     
     Call InterpolateCONSTANTS(mass_dens_red, crze, vise, temperature_magres, numdens_0, alpha, beta, res_constant)
     n_i_tot=rho_H2*alpha*(rho_H2/mass_dens_red)**beta
     
     
     if (rho_H2>temporary_max) then
      temporary_max=rho_H2
     end if
     counter=counter+1
     temporary_logmean=temporary_logmean * (counter - 1) / counter + log10(rho_H2) / counter 
                                             
     
     eta_perp = res_constant * va**2 / n_i_tot           !!
                                                           
     eta_paral = 2.74e6 * rho_H2 / n_i_tot
     
     eta_hall = 0.0694 * B_magres / n_i_tot                                        !!
     return                                                !!


subroutine InterpolateCONSTANTS(mass_dens_int, z, Av, T, numdens_0, alpha, beta, res_constant) 


#include "Flash.h"
use MagneticResistivity_data , ONLY : IntMassD,A_Fid,B_Fid,A_Lowz,B_Lowz,A_Highz,B_Highz,A_LowAv,B_LowAv,A_MediumAv,B_MediumAv,&
                                               A_LowT,B_LowT,A_HighT,B_HighT,res_amu,A_HighDens,B_HighDens

implicit none
integer :: N, i, size=201
real :: mass_dens_int, z, Av, T, D, zx, Ax, Tx, Dx, rho1, rho2, MIa, MIb, rho_adj, numdens_0, DT
real :: alpha1, alpha2,  beta1, beta2, alpha, beta
real :: res_constant, FidConst, LowzConst, HighzConst, LowAvConst, MediumAvConst, HighTConst, LowTConst
real :: HighDensConst
real :: zConst, AvConst, TConst, ka1, kb1, ka2, kb2 
real :: la1, lb1, la2, lb2, y, expAv, expT
real, dimension(2) :: af, aT, aA, az, aD, bf, bz, bA, bT, bD

FidConst=7.5e-12
LowzConst=7.4e-12
HighzConst=7.0e-12
LowAvConst=7.3e-12
MediumAvConst=7.4e-12
LowTConst=7.8e-12
HighTConst=7.6e-12
HighDensConst=7.2e-12


rho_adj=mass_dens_int*(300/numdens_0 &
+(log10(mass_dens_int/(numdens_0*2*res_amu))/(numdens_0/300*3))**2)/ &
(1+(log10(mass_dens_int/(numdens_0*2*res_amu))/(numdens_0/300*3))**2)


IF (rho_adj<IntMassD(size)) THEN
do i=2,size
 if (IntMassD(i)>rho_adj) then
  N=i
  exit
 end if
end do

rho1=IntMassD(N-1)
rho2=IntMassD(N)

af(1)=A_Fid(N-1)
bf(1)=B_Fid(N-1)

af(2)=A_Fid(N)
bf(2)=B_Fid(N)

If (z>=1.) then
 zx=2.
 az(1)=A_Highz(N-1)
 bz(1)=B_Highz(N-1)
 
 az(2)=A_Highz(N)
 bz(2)=B_Highz(N)
 
 zConst=z2Const
Else
 zx=0.5
 az(1)=A_Lowz(N-1)
 bz(1)=B_Lowz(N-1)
 
 az(2)=A_Lowz(N)
 bz(2)=B_Lowz(N)
 
 zConst=Z05Const
End if


If (Av>=5) then
 Ax=5.
 aA(1)=A_MediumAv(N-1)
 bA(1)=B_MediumAv(N-1)
  
 aA(2)=A_MediumAv(N)
 bA(2)=B_MediumAv(N)
 
 AvConst=MediumAvConst
 
Else
 Ax=3.
 aA(1)=A_LowAv(N-1)
 bA(1)=B_LowAv(N-1)
  
 aA(2)=A_LowAv(N)
 bA(2)=B_LowAv(N)
 
 AvConst=LowAvConst
End if


If (T>=10) then
 Tx=15.
 DT=450.
 aT(1)=A_HighT(N-1)
 bT(1)=B_HighT(N-1)

 aT(2)=A_HighT(N)
 bT(2)=B_HighT(N)

 TConst=HighTConst
Else
 Tx=6.
 DT=180.
 aT(1)=A_LowT(N-1)
 bT(1)=B_LowT(N-1)

 aT(2)=A_LowT(N)
 bT(2)=B_LowT(N)

 TConst=LowTConst
End if


Dx=750.
aD(1)=A_HighDens(N-1)
bD(1)=B_HighDens(N-1)

aD(2)=A_HighDens(N)
bD(2)=B_HighDens(N)


expAv=(exp(-Av/1.)-exp(-10./1.))/(exp(-Ax/1.)-exp(-10./1.))


ka1=log10(aD(1)/af(1))/log10(Dx/300.)
la1=(log10(aT(1)/af(1))-ka1*log10(DT/300.))/log10(Tx/10)
alpha1=(af(1) + 300/numdens_0*(az(1)-af(1))*(z-1.)/(zx-1.))&
*(T/10)**la1 *(numdens_0/300.)**ka1 &
* (aA(1)/af(1))**expAv

kb1=log10(bD(1)/bf(1))/log10(Dx/300)
lb1=(log10(bT(1)/bf(1))-kb1*log10(DT/300.))/log10(Tx/10)
beta1=(bf(1) + 300/numdens_0*(bz(1)-bf(1))*(z-1.)/(zx-1.))&
*(T/10)**lb1 *(numdens_0/300.)**kb1 &
* (bA(1)/bf(1))**expAv

ka2=log10(aD(2)/af(2))/log10(Dx/300)
la2=(log10(aT(2)/af(2))-ka2*log10(DT/300.))/log10(Tx/10)
alpha2=(af(2) + 300/numdens_0*(az(2)-af(2))*(z-1.)/(zx-1.))&
*(T/10)**la2 *(numdens_0/300.)**ka2 &
* (aA(2)/af(2))**expAv

kb2=log10(bD(2)/bf(2))/log10(Dx/300)
lb2=(log10(bT(2)/bf(2))-kb2*log10(DT/300.))/log10(Tx/10)
beta2=(bf(2) + 300/numdens_0*(bz(2)-bf(2))*(z-1.)/(zx-1.))&
*(T/10)**lb2 *(numdens_0/300.)**kb2 &
* (bA(2)/bf(2))**expAv


alpha=alpha1*(rho_adj-rho1)/(rho2-rho1)+alpha2*(rho2-rho_adj)/(rho2-rho1)
beta=beta1*(rho_adj-rho1)/(rho2-rho1)+beta2*(rho2-rho_adj)/(rho2-rho1)

res_constant=FidConst+(TConst-FidConst)*(T-10.)/(Tx-10.)+(AvConst-FidConst)*(Av-10.)/(Ax-10.)+(zConst-FidConst)*(z-1.)/(zx-1.) & 
+(D750Const-FidConst)*(numdens_0-300.-(DT-300.)*(T-10.)/(Tx-10.))/(Dx-300.)

return

ELSE

rho1=IntMassD(size-5)
rho2=IntMassD(size)

af(1)=A_Fid(size-5)
bf(1)=B_Fid(size-5)

af(2)=A_Fid(size)
bf(2)=B_Fid(size)

If (z>=1.) then
 zx=2.
 az(1)=A_Highz(size-5)
 bz(1)=B_Highz(size-5)
 
 az(2)=A_Highz(size)
 bz(2)=B_Highz(size)
 
 zConst=HighzConst
Else
 zx=0.5
 az(1)=A_Lowz(size-5)
 bz(1)=B_Lowz(size-5)
 
 az(2)=A_Lowz(size)
 bz(2)=B_Lowz(size)
 
 zConst=LowzConst
End if


If (Av>=5) then
 Ax=5.
 aA(1)=A_MediumAv(size-5)
 bA(1)=B_MediumAb(size-5)
  
 aA(2)=A_MediumAv(size)
 bA(2)=B_MediumAv(size)
 
 AvConst=MediumAvConst
 
Else
 Ax=3.
 aA(1)=A_LowAv(size-5)
 bA(1)=B_LowAv(size-5)
  
 aA(2)=A_LowAv(size)
 bA(2)=B_LowAv(size)
 
 AvConst=LowAvConst
End if


If (T>=10) then
 Tx=15.
 DT=450.
 aT(1)=A_HighT(size-5)
 bT(1)=B_HighT(size-5)

 aT(2)=A_HighT(size)
 bT(2)=B_HighT(size)

 TConst=HighTConst
Else
 Tx=6.
 DT=180.
 aT(1)=A_LowT(size-5)
 bT(1)=B_LowT(size-5)

 aT(2)=A_LowT(size)
 bT(2)=B_LowT(size)

 TConst=LowTConst
End if


Dx=750.
aD(1)=A_HighDens(size-5)
bD(1)=B_HighDens(size-5)

aD(2)=A_HighDens(size)
bD(2)=B_HighDens(size)


expAv=(exp(-Av/1.)-exp(-10./1.))/(exp(-Ax/1.)-exp(-10./1.))

ka1=log10(aD(1)/af(1))/log10(Dx/300.)
la1=(log10(aT(1)/af(1))-ka1*log10(DT/300.))/log10(Tx/10)
alpha1=(af(1) + 300/numdens_0*(az(1)-af(1))*(z-1.)/(zx-1.))&
*(T/10)**la1 *(numdens_0/300.)**ka1 &
* (aA(1)/af(1))**expAv

kb1=log10(bD(1)/bf(1))/log10(Dx/300)
lb1=(log10(bT(1)/bf(1))-kb1*log10(DT/300.))/log10(Tx/10)
beta1=(bf(1) + 300/numdens_0*(bz(1)-bf(1))*(z-1.)/(zx-1.))&
*(T/10)**lb1 *(numdens_0/300.)**kb1 &
* (bA(1)/bf(1))**expAv


ka2=log10(aD(2)/af(2))/log10(Dx/300)
la2=(log10(aT(2)/af(2))-ka2*log10(DT/300.))/log10(Tx/10)
alpha2=(af(2) + 300/numdens_0*(az(2)-af(2))*(z-1.)/(zx-1.))&
*(T/10)**la2 *(numdens_0/300.)**ka2 &
* (aA(2)/af(2))**expAv

kb2=log10(bD(2)/bf(2))/log10(Dx/300)
lb2=(log10(bT(2)/bf(2))-kb2*log10(DT/300.))/log10(Tx/10)
beta2=(bf(2) + 300/numdens_0*(bz(2)-bf(2))*(z-1.)/(zx-1.))&
*(T/10)**lb2 *(numdens_0/300.)**kb2 &
* (bA(2)/bf(2))**expAv

MIa=log10(alpha2/alpha1)/log10(rho2/rho1)
MIb=log10(beta2/beta1)/log10(rho2/rho1)

alpha=alpha1*(rho_adj/rho1)**MIa
beta=beta1*(rho_adj/rho1)**MIb

res_constant=FidConst+(TConst-FidConst)*(T-10.)/(Tx-10.)+(AvConst-FidConst)*(Av-10.)/(Ax-10.)+(zConst-FidConst)*(z-1.)/(zx-1.) & 
+(HighDensConst-FidConst)*(numdens_0-300.-(DT-300.)*(T-10.)/(Tx-10.))/(Dx-300.)


return

END IF
End subroutine InterpolateCONSTANTS
