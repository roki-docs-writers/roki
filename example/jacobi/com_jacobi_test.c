#include <roki/rk_jacobi.h>

#define N 8

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
    zVec3DCreate( rkChainLinkCOM(chain,i), zRandF(-1,1), zRandF(-1,1), zRandF(-1,1) );
    rkLinkSetMass( rkChainLink(chain,i), zRandF( 0.1, 1.0 ) );
    rkChainMass( chain ) += rkChainLinkMass( chain, i );
    zVec3DCreate( &aa, zRandF(-1,1), zRandF(-1,1), zRandF(-1,1) );
    zVec3DCreate( rkChainLinkOrgPos(chain,i), zRandF(-1,1), zRandF(-1,1), zRandF(-1,1) );
    zMat3DAA( rkChainLinkOrgAtt(chain,i), &aa );
    zNameSet( rkChainLink(chain,i), name );
  }
  rkLinkAddChild( rkChainLink(chain,0), rkChainLink(chain,1) );
  rkLinkAddChild( rkChainLink(chain,1), rkChainLink(chain,2) );
  rkLinkAddChild( rkChainLink(chain,2), rkChainLink(chain,3) );
  rkLinkAddChild( rkChainLink(chain,3), rkChainLink(chain,4) );
  rkLinkAddChild( rkChainLink(chain,0), rkChainLink(chain,5) );
  rkLinkAddChild( rkChainLink(chain,5), rkChainLink(chain,6) );
  rkLinkAddChild( rkChainLink(chain,6), rkChainLink(chain,7) );
  rkJointCreate( rkChainLinkJoint(chain,0), RK_JOINT_FLOAT );
  rkJointCreate( rkChainLinkJoint(chain,1), RK_JOINT_SPHER );
  rkJointCreate( rkChainLinkJoint(chain,2), RK_JOINT_REVOL );
  rkJointCreate( rkChainLinkJoint(chain,3), RK_JOINT_CYLIN );
  rkJointCreate( rkChainLinkJoint(chain,4), RK_JOINT_REVOL );
  rkJointCreate( rkChainLinkJoint(chain,5), RK_JOINT_PRISM );
  rkJointCreate( rkChainLinkJoint(chain,6), RK_JOINT_HOOKE );
  rkJointCreate( rkChainLinkJoint(chain,7), RK_JOINT_FIXED );

  rkChainSetOffset( chain );
  rkChainUpdateFK( chain );
  rkChainUpdateID( chain );
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

  zVecRandUniform( dis, -10.0, 10.0 );
  zVecRandUniform( vel, -10.0, 10.0 );
  rkChainFK( &chain, dis );
  rkChainID( &chain, vel, acc );

  rkChainCOMJacobi( &chain, jacobi );
  zVec3DCopy( rkChainCOMVel(&chain), &v );
  zMulMatVec( jacobi, vel, ev );

  zVec3DSub( (zVec3D*)zVecArray(ev), &v, &err );
  printf( "%g %g %g ",  zVecElem(ev,0), zVecElem(ev,1), zVecElem(ev,2) );
  printf( "%g %g %g ",  zVec3DElem(&v,0), zVec3DElem(&v,1), zVec3DElem(&v,2) );
  printf( "%g %g %g\n", zVec3DElem(&err,0), zVec3DElem(&err,1), zVec3DElem(&err,2) );
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
