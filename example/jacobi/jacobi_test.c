#include <roki/rk_jacobi.h>

#define N   8
#define TIP 4
#define ANO 7

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
    zNameSet( rkChainLink(chain,i), name );
    zVec3DCreate( &aa, zRandF(-1,1), zRandF(-1,1), zRandF(-1,1) );
    zVec3DCreate( rkChainLinkOrgPos(chain,i), zRandF(-1,1), zRandF(-1,1), zRandF(-1,1) );
    zMat3DAA( rkChainLinkOrgAtt(chain,i), &aa );
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
  rkChainSetMass( chain, 1.0 ); /* dummy */
  rkChainUpdateFK( chain );
  rkChainUpdateID( chain );
}

void world_ang_test(rkChain *chain, zMat jacobi, zVec3D *v)
{
  eprintf( "World (angular) test\n" );
  rkChainLinkWldAngJacobi( chain, TIP, jacobi );
  zMulMatVec3D( rkChainLinkWldAtt(chain,TIP), rkChainLinkAngVel(chain,TIP), v );
}

void world_lin_test(rkChain *chain, zMat jacobi, zVec3D *v)
{
  eprintf( "World (linear) test\n" );
  rkChainLinkWldLinJacobi( chain, TIP, Z_ZEROVEC3D, jacobi );
  zMulMatVec3D( rkChainLinkWldAtt(chain,TIP), rkChainLinkLinVel(chain,TIP), v );
}

void world_com_test(rkChain *chain, zMat jacobi, zVec3D *v)
{
  eprintf( "World (COM) test\n" );
  rkChainLinkWldLinJacobi( chain, TIP, rkChainLinkCOM(chain,TIP), jacobi );
  zMulMatVec3D( rkChainLinkWldAtt(chain,TIP), rkChainLinkCOMVel(chain,TIP), v );
}

void l2l_ang_test(rkChain *chain, zMat jacobi, zVec3D *v)
{
  zVec3D av;

  eprintf( "Link->Link (angular) test\n" );
  rkChainLinkToLinkAngJacobi( chain, ANO, TIP, jacobi );
  zMulMatVec3D( rkChainLinkWldAtt(chain,TIP), rkChainLinkAngVel(chain,TIP), v );
  zMulMatVec3D( rkChainLinkWldAtt(chain,ANO), rkChainLinkAngVel(chain,ANO), &av );
  zVec3DSubDRC( v, &av );
}

void l2l_lin_test(rkChain *chain, zMat jacobi, zVec3D *v)
{
  zVec3D av, vr, tmp;

  eprintf( "Link->Link (linear) test\n" );
  rkChainLinkToLinkLinJacobi( chain, ANO, TIP, Z_ZEROVEC3D, jacobi );
  zMulMatVec3D( rkChainLinkWldAtt(chain,TIP), rkChainLinkLinVel(chain,TIP), v );
  zMulMatVec3D( rkChainLinkWldAtt(chain,ANO), rkChainLinkLinVel(chain,ANO), &av );
  zVec3DSubDRC( v, &av );

  zMulMatVec3D( rkChainLinkWldAtt(chain,ANO), rkChainLinkAngVel(chain,ANO), &vr );
  zVec3DSub( rkChainLinkWldPos(chain,TIP), rkChainLinkWldPos(chain,ANO), &av );
  zVec3DOuterProd( &vr, &av, &tmp );
  zVec3DSubDRC( v, &tmp );
}

void (*testfunc[])(rkChain*,zMat,zVec3D*) = {
  world_ang_test,
  world_lin_test,
  world_com_test,
  l2l_ang_test,
  l2l_lin_test,
  NULL,
};

int main(int argc, char *argv[])
{
  rkChain chain;
  zVec dis, vel, acc, ev;
  zMat jacobi;
  zVec3D v, err;
  void (*test)(rkChain*,zMat,zVec3D*);

  if( !( test = testfunc[argc>1 ? atoi(argv[1]) : 0] ) ){
    ZRUNERROR( "invalid test indicator: %d", atoi(argv[1]) );
    return EXIT_FAILURE;
  }
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
zFrame3DWrite(rkChainLinkWldFrame(&chain,TIP));
  test( &chain, jacobi, &v );
  zMulMatVec( jacobi, vel, ev );
  zVec3DSub( (zVec3D*)zVecBuf(ev), &v, &err );
  zVec3DWrite( (zVec3D*)zVecBuf(ev) );
  zVec3DWrite( &v );
  zVec3DWrite( &err );
/*
  printf( "%g %g %g ",  zVecElem(ev,0), zVecElem(ev,1), zVecElem(ev,2) );
  printf( "%g %g %g ",  zVec3DElem(&v,0), zVec3DElem(&v,1), zVec3DElem(&v,2) );
  printf( "%g %g %g\n", zVec3DElem(&err,0), zVec3DElem(&err,1), zVec3DElem(&err,2) );
*/
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
