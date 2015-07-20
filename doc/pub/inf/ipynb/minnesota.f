real*8 function vmom_minnesota(p,q,r,s)
    USE single_particle_orbits
    USE constants
    use chiral_constants
  
    implicit none 
    INTEGER :: p,q,r,s, m1,m2,m3,m4, spin, iph, t1,t2,t3,t4, Tiso
    INTEGER :: nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3, nx4, ny4, nz4 
    REAL*8 :: k1(3), k2(3), k3(3), k4(3), kmean(3) 
    REAL*8 :: qtrans(3), prel(3), pprel(3)
    REAL*8 :: q2, p2, qabs, pp2, vdir, vexc, cg1, cg2 
    REAL*8 :: delta, nucleon_mass, relativity_factor
    REAL*8 :: v0r, v0t, v0s, kr,kt,ks, vr,vt,vs,vs_ex,vcentral, vsigma, spin_exc1, spin_exc2
  
    nx1 = all_orbit%nx(p)
    ny1 = all_orbit%ny(p)
    nz1 = all_orbit%nz(p)
    nx2 = all_orbit%nx(q)
    ny2 = all_orbit%ny(q)
    nz2 = all_orbit%nz(q)
    nx3 = all_orbit%nx(r)
    ny3 = all_orbit%ny(r)
    nz3 = all_orbit%nz(r)
    nx4 = all_orbit%nx(s)
    ny4 = all_orbit%ny(s)
    nz4 = all_orbit%nz(s)
    ! 
    ! Conservation of linear momentum
    !
    if ( nx1 + nx2 /= nx3 + nx4 ) return 
    if ( ny1 + ny2 /= ny3 + ny4 ) return 
    if ( nz1 + nz2 /= nz3 + nz4 ) return 
  
    k1(1) = all_orbit%kx(p)
    k1(2) = all_orbit%ky(p)
    k1(3) = all_orbit%kz(p)
    k2(1) = all_orbit%kx(q)
    k2(2) = all_orbit%ky(q)
    k2(3) = all_orbit%kz(q)
    k3(1) = all_orbit%kx(r)
    k3(2) = all_orbit%ky(r)
    k3(3) = all_orbit%kz(r)
    k4(1) = all_orbit%kx(s)
    k4(2) = all_orbit%ky(s)
    k4(3) = all_orbit%kz(s)
  
  
    ! 
    ! conservation of spin and isospin 
    !
    if ( all_orbit%szp(p) + all_orbit%szp(q) /= all_orbit%szp(r) + all_orbit%szp(s) ) return 
    if ( all_orbit%itzp(p) + all_orbit%itzp(q) /= all_orbit%itzp(r) + all_orbit%itzp(s) ) return 
  
    m1 = all_orbit%szp(p) 
    m2 = all_orbit%szp(q) 
    m3 = all_orbit%szp(r) 
    m4 = all_orbit%szp(s) 
  
    t1 = all_orbit%itzp(p) 
    t2 = all_orbit%itzp(q) 
    t3 = all_orbit%itzp(r) 
    t4 = all_orbit%itzp(s) 
  
    ! 
    ! RELATIVE MOMENTA <prel |v| pprel > 
    ! 
    prel  = 0.5d0*(k1-k2)
    pprel = 0.5d0*(k3-k4)
    !
    ! momentum transfer 
    !
    qtrans = prel - pprel
    q2 = sum(qtrans*qtrans) 
    
    
    v0r=200.0  ! MeV
    v0t=178.0  ! MeV
    v0s=91.85  ! MeV
    kr=1.487  ! fm**-2
    kt=0.639  ! fm**-2
    ks=0.465  ! fm**-2

  
    ! r-space 
    !vr=v0r*exp(-kr*rr**2)
    !vt=-v0t*exp(-kt*rr**2)
    !vs=-v0s*exp(-ks*rr**2)
    
    vr =  v0r * pi**1.5d0 * exp(-q2/(4.d0*kr) )/ (kr**1.5d0) 
    vt = -v0t * pi**1.5d0 * exp(-q2/(4.d0*kt) )/ (kt**1.5d0)
    vs = -v0s * pi**1.5d0 * exp(-q2/(4.d0*ks) )/ (ks**1.5d0)
   
    
    vr = vr * (3.d0*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4) - & 
         chp_tau_dot_tau(t1,t2,t3,t4)*delta(m1,m3)*delta(m2,m4) - & 
         chp_sigma_dot_sigma(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4) - & 
         chp_tau_dot_tau(t1,t2,t3,t4)*chp_sigma_dot_sigma(m1,m2,m3,m4) )/8.d0
    
    vs = vs * (3.d0*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4) + & 
         chp_tau_dot_tau(t1,t2,t3,t4)*delta(m1,m3)*delta(m2,m4) - & 
         3.d0*chp_sigma_dot_sigma(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4) - & 
         chp_tau_dot_tau(t1,t2,t3,t4)*chp_sigma_dot_sigma(m1,m2,m3,m4) )/16.d0
    
    vt = vt * (3.d0*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4) - & 
         3.d0*chp_tau_dot_tau(t1,t2,t3,t4)*delta(m1,m3)*delta(m2,m4) + & 
         chp_sigma_dot_sigma(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4) - & 
         chp_tau_dot_tau(t1,t2,t3,t4)*chp_sigma_dot_sigma(m1,m2,m3,m4) )/16.d0
    
   
    vdir = vs+vr+vt
    
    vmom_minnesota = vdir/volume
    

  end function vmom_minnesota


FUNCTION chp_tau_dot_tau_mtx(ms1,ms2,ms3,ms4) RESULT(res)
    
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8 :: res, delta
    complex*16 :: res1
    COMPLEX*16 :: chi1(2), chi2(2), chi3(2), chi4(2)
    INTEGER :: i1 
    
    res = 0.0D0
    chi1 = 0.d0; chi2 = 0.d0; chi3 = 0.d0; chi4 = 0.d0
    i1 = nint(1.5-0.5*ms1) 
    chi1(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms2) 
    chi2(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms3) 
    chi3(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms4) 
    chi4(i1) = 1.d0 
    
    
    res1 = dot_product(chi1, matmul( sigma_x, chi3)) * dot_product(chi2, matmul( sigma_x, chi4)) &
         + dot_product(chi1, matmul( sigma_y, chi3)) * dot_product(chi2, matmul( sigma_y, chi4)) &
         + dot_product(chi1, matmul( sigma_z, chi3)) * dot_product(chi2, matmul( sigma_z, chi4)) 
    res = res1
    
    
  END FUNCTION chp_tau_dot_tau_mtx
  
