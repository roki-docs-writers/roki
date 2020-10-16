#include <roki/rk_chain.h>

#define CHAIN_FILE "../model/humanoid.zkc"

#define DT 0.001
#define DTHETA zDeg2Rad(20.0)
#define STEP 100

rkChain chain;
zVec dis;

zVec3D com, vel, acc;

void update_dis(int st)
{
  register int i;
  double val;
  zMat3D r;

  val = DTHETA * ( 1 - cos(2*zPI*st/STEP) );
  zMat3DZYX( &r, val, val, val );
  zMat3DToAA( &r, (zVec3D*)&zVecElem(dis,3) );
  for( i=6; i<_zVecSize(dis); i++ )
    zVecSetElem( dis, i, val );
}

void world_com_test(void)
{
  zVec3D v, a;

  printf( "%.10f %.10f %.10f ", zVec3DElem(rkChainWldCOM(&chain),zX), zVec3DElem(rkChainWldCOM(&chain),zY), zVec3DElem(rkChainWldCOM(&chain),zZ) );
  printf( "%.10f %.10f %.10f ", zVec3DElem(rkChainCOMVel(&chain),zX), zVec3DElem(rkChainCOMVel(&chain),zY), zVec3DElem(rkChainCOMVel(&chain),zZ) );
  printf( "%.10f %.10f %.10f ", zVec3DElem(rkChainCOMAcc(&chain),zX), zVec3DElem(rkChainCOMAcc(&chain),zY), zVec3DElem(rkChainCOMAcc(&chain),zZ) );

  zVec3DSub( rkChainWldCOM(&chain), &com, &v );
  zVec3DDivDRC( &v, DT );
  printf( "%.10f %.10f %.10f ", zVec3DElem(&v,zX), zVec3DElem(&v,zY), zVec3DElem(&v,zZ) );
  zVec3DSub( &v, &vel, &a );
  zVec3DDivDRC( &a, DT );
  printf( "%.10f %.10f %.10f\n", zVec3DElem(&a,zX), zVec3DElem(&a,zY), zVec3DElem(&a,zZ) );

  zVec3DCopy( rkChainWldCOM(&chain), &com );
  zVec3DCopy( &v, &vel );
}

int main(void)
{
  register int i;

  rkChainReadFile( &chain, CHAIN_FILE );
  rkChainUpdateID( &chain );
  zVec3DCopy( rkChainWldCOM(&chain), &com );
  zVec3DCopy( rkChainCOMVel(&chain), &vel );
  dis = zVecAlloc( rkChainJointSize( &chain ) );
  rkChainGetJointDisAll( &chain, dis );
  for( i=0; i<=STEP; i++ ){
    update_dis( i );
    rkChainFKCNT( &chain, dis, DT );
    world_com_test();
  }
  zVecFree( dis );
  rkChainDestroy( &chain );
  return 0;
}
