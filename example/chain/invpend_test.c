/* this code is to check if the inverse dynamics, angular momentum and
   kinematic energy computation are correct, comparing with the
   analytical solutions on a cylindrical inverted pendulum.
 */
#include <roki/rk_chain.h>

/* the following parameters must be consistent with invpend.zkc */
#define H  0.3
#define L1 0.49
#define L2 0.195

#define I1 0.075
#define I2 0.0011403
#define M2 0.085

#define DT 0.001

zVec6D *invpend_vel(rkChain *ip, double t1, double t2, double dt1, double dt2, zVec6D *v)
{
  double s1, c1, s2, c2;

  zSinCos( t1, &s1, &c1 );
  zSinCos( t2, &s2, &c2 );
  return zVec6DCreate( v,
   -L1*dt1*s1+L2*(dt2*s1*c2+dt1*c1*s2),
    L1*dt1*c1-L2*(dt2*c1*c2-dt1*s1*s2),
             -L2* dt2   *s2,
    dt2*c1, dt2*s1, dt1 );
}

zVec6D *invpend_acc(rkChain *ip, double t1, double t2, double dt1, double dt2, double ddt1, double ddt2, zVec6D *a)
{
  double s1, c1, s2, c2;

  zSinCos( t1, &s1, &c1 );
  zSinCos( t2, &s2, &c2 );
  return zVec6DCreate( a,
    (-L1*s1+L2*c1*s2)*ddt1+L2*s1*c2*ddt2+(-L1*c1-L2*s1*s2)*dt1*dt1+2*L2*c1*c2*dt1*dt2-L2*s1*s2*dt2*dt2,
    ( L1*c1+L2*s1*s2)*ddt1-L2*c1*c2*ddt2+(-L1*s1+L2*c1*s2)*dt1*dt1+2*L2*s1*c2*dt1*dt2+L2*c1*s2*dt2*dt2,
     RK_G - L2*s2*ddt2 - L2*c2*dt2*dt2,
    c1*ddt2-s1*dt1*dt2, s1*ddt2+c1*dt1*dt2, ddt1 );
}

zVec3D *invpend_am(rkChain *ip, double t1, double t2, double dt1, double dt2, zVec3D *am)
{
  double s1, c1, s2, c2;

  zSinCos( t1, &s1, &c1 );
  zSinCos( t2, &s2, &c2 );
  return zVec3DCreate( am,
   -(M2*(H*(L1*c1+L2*s1*s2)+L1*L2*c1*c2+L2*L2*s1*s2*c2)+I2*s1*s2*c2)*dt1 + (M2*(H*L2*c1*c2-L1*L2*s1*s2+L2*L2*c1)+I2*c1)*dt2,
    (M2*(H*(-L1*s1+L2*c1*s2)-L1*L2*s1*c2+L2*L2*c1*s2*c2)+I2*c1*s2*c2)*dt1 + (M2*(H*L2*s1*c2+L1*L2*c1*s2+L2*L2*s1)+I2*s1)*dt2,
    (I1+M2*(L1*L1+L2*L2*s2*s2)+I2*s2*s2)*dt1 - M2*L1*L2*c2*dt2 );
}

double invpend_ke(rkChain *ip, double t1, double t2, double dt1, double dt2)
{
  double s1, c1, s2, c2;

  zSinCos( t1, &s1, &c1 );
  zSinCos( t2, &s2, &c2 );
  return 0.5*(I1+M2*L1*L1+(I2+M2*L2*L2)*s2*s2)*dt1*dt1
               - M2*L1*L2*c2*dt1*dt2
       + 0.5*(I2+M2*L2*L2)*dt2*dt2;
}

void invpend_trq(rkChain *ip, double t1, double t2, double dt1, double dt2, double ddt1, double ddt2, double *u1, double *u2)
{
  double s1, c1, s2, c2;

  zSinCos( t1, &s1, &c1 );
  zSinCos( t2, &s2, &c2 );
  *u1 = (I1+M2*L1*L1+(M2*L2*L2+I2)*s2*s2)*ddt1 - M2*L1*L2*c2*ddt2 + 2*(M2*L2*L2+I2)*s2*c2*dt1*dt2 + M2*L1*L2*s2*dt2*dt2;
  *u2 = -M2*L1*L2*c2*ddt1 + (M2*L2*L2+I2)*ddt2 - (M2*L2*L2+I2)*s2*c2*dt1*dt1 - M2*RK_G*L2*s2;
}


void set_joint_state(zVec t, zVec dt, zVec ddt)
{
  zVecSetElem( t, 0, zRandF(-zPI,zPI) );
  zVecSetElem( t, 1, zRandF(-zPI,zPI) );
  zVecSetElem( dt, 0, zRandF(-zPI,zPI) * 0.1 / DT );
  zVecSetElem( dt, 1, zRandF(-zPI,zPI) * 0.1 / DT );
  zVecSetElem( ddt, 0, zRandF(-zPI,zPI) * zSqr(0.1/DT) );
  zVecSetElem( ddt, 1, zRandF(-zPI,zPI) * zSqr(0.1/DT) );
}

int main(void)
{
  rkChain invpend;
  zVec t, dt, ddt;
  zVec6D v1, v2, a1, a2, e6;
  zVec3D am1, am2, e3;
  double k1, k2, u11, u12, u21, u22;

  zRandInit();
  rkChainReadFile( &invpend, "../model/invpend.zkc" );
  t   = zVecAlloc( rkChainJointSize(&invpend) );
  dt  = zVecAlloc( rkChainJointSize(&invpend) );
  ddt = zVecAlloc( rkChainJointSize(&invpend) );

  set_joint_state( t, dt, ddt );
  rkChainFK( &invpend, t );
  rkChainID( &invpend, dt, ddt );

  printf( "<velocity test>\n" );
  printf( " - analytical computation\n" );
  invpend_vel( &invpend, zVecElem(t,0), zVecElem(t,1), zVecElem(dt,0), zVecElem(dt,1), &v1 );
  zVec6DWrite( &v1 );
  printf( " - recursive computation\n" );
  zMulMatVec3D( rkChainLinkWldAtt(&invpend,2), rkChainLinkCOMVel(&invpend,2), zVec6DLin(&v2) );
  zMulMatVec3D( rkChainLinkWldAtt(&invpend,2), rkChainLinkAngVel(&invpend,2), zVec6DAng(&v2) );
  zVec6DWrite( &v2 );
  printf( "*** (error)\n" ); zVec6DWrite( zVec6DSub( &v1, &v2, &e6 ) );

  printf( "<acceleration test>\n" );
  printf( " - analytical computation\n" );
  invpend_acc( &invpend, zVecElem(t,0), zVecElem(t,1), zVecElem(dt,0), zVecElem(dt,1), zVecElem(ddt,0), zVecElem(ddt,1), &a1 );
  zVec6DWrite( &a1 );
  printf( " - recursive computation\n" );
  zMulMatVec3D( rkChainLinkWldAtt(&invpend,2), rkChainLinkCOMAcc(&invpend,2), zVec6DLin(&a2) );
  zMulMatVec3D( rkChainLinkWldAtt(&invpend,2), rkChainLinkAngAcc(&invpend,2), zVec6DAng(&a2) );
  zVec6DWrite( &a2 );
  printf( "*** (error)\n" ); zVec6DWrite( zVec6DSub( &a1, &a2, &e6 ) );

  printf( "<angular momentum test>\n" );
  printf( " - analytical computation\n" );
  invpend_am( &invpend, zVecElem(t,0), zVecElem(t,1), zVecElem(dt,0), zVecElem(dt,1), &am1 );
  zVec3DWrite( &am1 );
  printf( " - recursive computation\n" );
  zVec3DWrite( rkChainAM( &invpend, Z_ZEROVEC3D, &am2 ) );
  printf( "*** (error)\n" ); zVec3DWrite( zVec3DSub( &am1, &am2, &e3 ) );

  printf( "<kinematic energy test>\n" );
  printf( " - analytical computation\n" );
  printf( "%.16g\n", ( k1 = invpend_ke( &invpend, zVecElem(t,0), zVecElem(t,1), zVecElem(dt,0), zVecElem(dt,1) ) ) );
  printf( " - recursive computation\n" );
  printf( "%.16g\n", ( k2 = rkChainKE( &invpend ) ) );
  printf( "*** (error) = %.16g\n", k1 - k2 );

  printf( "<torque test>\n" );
  printf( " - analytical computation\n" );
  invpend_trq( &invpend, zVecElem(t,0), zVecElem(t,1), zVecElem(dt,0), zVecElem(dt,1), zVecElem(ddt,0), zVecElem(ddt,1), &u11, &u12 );
  printf( "%.16f %.16f\n", u11, u12 );
  printf( " - recursive computation\n" );
  rkChainLinkGetJointTrq( &invpend, 1, &u21 );
  rkChainLinkGetJointTrq( &invpend, 2, &u22 );
  printf( "%.16f %.16f\n", u21, u22 );
  printf( "*** (error) = %.16g, %.16g\n", u11 - u21, u12 - u22 );

  rkChainDestroy( &invpend );
  zVecFreeAO( 3, t, dt, ddt );
  return 0;
}
