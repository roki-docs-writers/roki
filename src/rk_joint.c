/* RoKi - Robot Kinetics library
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * rk_joint - joint structure
 */

#include <roki/rk_joint.h>

/* ********************************************************** */
/* joint type
 * ********************************************************** */
static char *__rkjointtypename[] = {
  "fix", "revolute", "prism", "cylinder", "hooke", "sphere", "float",
  NULL
};

/* rkJointTypeExpr
 * - convert joint type to string.
 */
char *rkJointTypeExpr(byte type)
{
  return __rkjointtypename[zLimit(type,RK_JOINT_FIXED,RK_JOINT_FLOAT)];
}

/* rkJointTypeFromStr
 * - convert string to joint type.
 */
byte rkJointTypeFromStr(char *str)
{
  char **jp;
  byte type;

  for( type=RK_JOINT_FIXED, jp=__rkjointtypename; *jp; jp++, type++ )
    if( !strcmp( str, *jp ) ) return type;
  return RK_JOINT_FIXED;
}

/* ********************************************************** */
/* CLASS: rkJoint
 * joint class
 * ********************************************************** */

static rkJoint *(* rk_joint_create[])(rkJoint*) = {
  rkJointCreateFixed,
  rkJointCreateRevol,
  rkJointCreatePrism,
  rkJointCreateCylin,
  rkJointCreateHooke,
  rkJointCreateSpher,
  rkJointCreateFloat,
};

/* rkJointCreate
 * - create joint object.
 */
rkJoint *rkJointCreate(rkJoint *j, byte type)
{
  if( type < RK_JOINT_FIXED || type > RK_JOINT_FLOAT ){
    ZRUNERROR( "invalid joint type specified - %d", type );
    return NULL;
  }
  rkJointInit( j );
  if( !rk_joint_create[( (j)->type = type )]( j ) ){
    ZRUNERROR( "cannot create joint instance" );
    rkJointDestroy( j );
    return NULL;
  }
  rkJointNeutral( j );
  return j;
}

/* rkJointDestroy
 * - destroy joint object.
 */
void rkJointDestroy(rkJoint *j)
{
  zFree( j->prp );
  rkJointInit( j );
}

/* rkJointNeutral
 * - neutralize joint displacement.
 */
void rkJointNeutral(rkJoint *j)
{
  double dis[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  rkJointSetDis( j, dis );
}

/* rkJointIsNeutral
 * - check if joint displacement is neutral.
 */
bool rkJointIsNeutral(rkJoint *j)
{
  double dis[6];
  register int i;

  rkJointGetDis( j, dis );
  for( i=0; i<rkJointSize(j); i++ )
    if( !zIsTiny( dis[i] ) ) return false;
  return true;
}

/* rkJointCopyState
 * - copy joint state.
 */
rkJoint *rkJointCopyState(rkJoint *src, rkJoint *dst)
{
  double val[6];

  rkJointGetDis( src, val ); rkJointSetDis( dst, val );
  rkJointGetVel( src, val ); rkJointSetVel( dst, val );
  rkJointGetAcc( src, val ); rkJointSetAcc( dst, val );
  rkJointGetTrq( src, val ); rkJointSetTrq( dst, val );
  return dst;
}

/* rkJointIncRate
 * - increment motion rate due to the joint rate.
 */
void rkJointIncRate(rkJoint *j, zVec3D *w, zVec6D *vel, zVec6D *acc)
{
  rkJointIncVel( j, w, vel, acc );
  rkJointIncAcc( j, acc );
}

/* NOTE: The following macros and functions are for sharing
 * some operation codes. Do not use them in users programs. */
zVec3D *_rkJointAxisNull(void *prp, zFrame3D *f, zVec3D *a){
  return NULL;
}

zVec3D *_rkJointAxisZ(void *prp, zFrame3D *f, zVec3D *a){
  zVec3DCopy( zMat3DVec(zFrame3DAtt(f),zZ), a );
  return a;
}

/* joint torsion */

double rkJointTorsionDisRevol(zFrame3D *dev, zVec6D *t)
{
  zMat3D rm;
  zVec3D aa;
  double l, angle;

  zVec3DCreate( &aa, -zMat3DElem(zFrame3DAtt(dev),1,2), zMat3DElem(zFrame3DAtt(dev),0,2), 0 );
  l = sqrt( zSqr(zVec3DElem(&aa,zX))+zSqr(zVec3DElem(&aa,zY)) );
  angle = atan2( l, zMat3DElem(zFrame3DAtt(dev),2,2) );
  zIsTiny( angle ) ?
    zVec3DClear( &aa ) : zVec3DMulDRC( &aa, angle/l );
  zMulMatTVec3D( zFrame3DAtt(dev), &aa, zVec6DAng(t) );
  /* intermediate attitude */
  zMat3DAA( &rm, &aa );
  /* joint displacement */
  return 0.5 *
    ( zVec3DAngle( zMat3DVec(&rm,zX), zMat3DVec(zFrame3DAtt(dev),zX), zMat3DVec(&rm,zZ) )
    + zVec3DAngle( zMat3DVec(&rm,zY), zMat3DVec(zFrame3DAtt(dev),zY), zMat3DVec(&rm,zZ) ) );
}

double rkJointTorsionDisPrism(zFrame3D *dev, zVec6D *t)
{
  double q;

  zMulMatTVec3D( zFrame3DAtt(dev), zFrame3DPos(dev), zVec6DLin(t) );
  /* joint displacement */
  q = zVec6DElem( t, zZ );
  zVec6DElem(t,zZ) = 0;
  return q;
}
