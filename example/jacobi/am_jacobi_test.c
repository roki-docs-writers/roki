#include <roki/rk_jacobi.h>

#define N   6
#define TIP 3

void link_mp_rand(rkLink *l)
{
  double i11, i12, i13, i22, i23, i33;
  rkLinkSetMass( l, zRandF(0.01,0.1) );
  zVec3DCreate( rkLinkCOM(l), zRandF(-1,1), zRandF(-1,1), zRandF(-1,1) );
  i11 = zRandF(0.01,0.1);
  i12 =-zRandF(0.0001,0.001);
  i13 =-zRandF(0.0001,0.001);
  i22 = zRandF(0.01,0.1);
  i23 =-zRandF(0.00001,0.0001);
  i33 = zRandF(0.01,0.1);
  zMat3DCreate( rkLinkInertia(l), i11, i12, i13, i12, i22, i23, i13, i23, i33 );
}

void chain_init(rkChain *chain)
{
  register int i;
  char name[BUFSIZ];
  zVec3D aa;

  rkChainInit( chain );
  zArrayAlloc( &chain->link, rkLink, N );
  for( i=0; i<N; i++ ){
    sprintf( name, "link#%02d", i );
    rkLinkInit( rkChainLink(chain,i) );
    zNameSet( rkChainLink(chain,i), name );
    link_mp_rand( rkChainLink(chain,i) );
    rkChainMass(chain) += rkChainLinkMass(chain,i);
    zVec3DCreate( &aa, zRandF(-1,1), zRandF(-1,1), zRandF(-1,1) );
    zVec3DCreate( rkChainLinkOrgPos(chain,i), zRandF(-1,1), zRandF(-1,1), zRandF(-1,1) );
    zMat3DAA( rkChainLinkOrgAtt(chain,i), &aa );
  }
  rkLinkAddChild( rkChainLink(chain,0), rkChainLink(chain,1) );
  rkLinkAddChild( rkChainLink(chain,1), rkChainLink(chain,2) );
  rkLinkAddChild( rkChainLink(chain,2), rkChainLink(chain,3) );
  rkLinkAddChild( rkChainLink(chain,0), rkChainLink(chain,4) );
  rkLinkAddChild( rkChainLink(chain,4), rkChainLink(chain,5) );
  rkJointCreate( rkChainLinkJoint(chain,0), RK_JOINT_FLOAT );
  rkJointCreate( rkChainLinkJoint(chain,1), RK_JOINT_REVOL );
  rkJointCreate( rkChainLinkJoint(chain,2), RK_JOINT_CYLIN );
  rkJointCreate( rkChainLinkJoint(chain,3), RK_JOINT_FIXED );
  rkJointCreate( rkChainLinkJoint(chain,4), RK_JOINT_SPHER );
  rkJointCreate( rkChainLinkJoint(chain,5), RK_JOINT_FIXED );
/*
  zVec3DCreate( rkChainLinkOrgPos(chain,1), 1, 0, 0 );
  zVec3DCreate( rkChainLinkOrgPos(chain,2), 0, 1, 0 );
  zVec3DCreate( rkChainLinkOrgPos(chain,3), 1, 0, 0 );
  zVec3DCreate( rkChainLinkOrgPos(chain,4), 0, 1, 0 );
  zVec3DCreate( rkChainLinkOrgPos(chain,5), 0, 2, 0 );
*/

  rkChainSetOffset( chain );
  rkChainUpdateFK( chain );
  rkChainUpdateID( chain );
}

void link_am_test(rkChain *chain, zMat jacobi, zVec3D *v)
{
  zVec3D tp;

  eprintf( "link angular momentum test\n" );
  rkChainLinkAMJacobi( chain, TIP, ZVEC3DZERO, jacobi );
  zXfer3DInv( rkChainLinkWldFrame(chain,TIP), ZVEC3DZERO, &tp );
  rkLinkAM( rkChainLink(chain,TIP), &tp, v );
  zMulMatVec3DDRC( rkChainLinkWldAtt(chain,TIP), v );
}

void am_test(rkChain *chain, zMat jacobi, zVec3D *v)
{
  eprintf( "chain angular momentum test\n" );
  rkChainAMJacobi( chain, ZVEC3DZERO, jacobi );
  rkChainAM( chain, ZVEC3DZERO, v );
}

int main(int argc, char *argv[])
{
  rkChain chain;
  zVec dis, vel, acc, ev;
  zMat jacobi;
  zVec3D v, err;

  /* initialization */
  zRandInit();
  chain_init( &chain );
  dis = zVecAlloc( rkChainJointSize(&chain) );
  vel = zVecAlloc( rkChainJointSize(&chain) );
  acc = zVecAlloc( rkChainJointSize(&chain) ); /* dummy */
  ev = zVecAlloc( 3 );
  jacobi = zMatAlloc( 3, rkChainJointSize(&chain) );
  rkChainGetJointDisAll( &chain, dis );

  zVecRandUniform( dis, -10, 10 );
  zVecRandUniform( vel, -10, 10 );
  rkChainFK( &chain, dis );
  rkChainID( &chain, vel, acc );
  argc > 1 && atoi(argv[1]) == 0 ?
    link_am_test( &chain, jacobi, &v ) : am_test( &chain, jacobi, &v );

  zMulMatVec( jacobi, vel, ev );
  zVec3DSub( (zVec3D*)zVecBuf(ev), &v, &err );
  zVec3DDataWrite( (zVec3D*)zVecBuf(ev) );
  zVec3DDataWrite( &v );
  zVec3DDataNLWrite( &err );
  printf( " ... %s.\n", zVec3DIsTiny( &err ) ? "OK" : "BUG probably in Jacobian matrix computation" );

  /* termination */
  zVecFree( dis );
  zVecFree( vel );
  zVecFree( acc );
  zVecFree( ev );
  zMatFree( jacobi );
  rkChainDestroy( &chain );
  return 0;
}
