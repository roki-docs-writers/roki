/* RoKi - Robot Kinetics library
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * rk_joint_hooke - joint structure: universal joint
 */

#include <roki/rk_joint.h>

static double _rkJointLimDis1Hooke(void *prp, int i, double testval);
static void _rkJointLimDisHooke(void *prp, double *testval, double *limval);
static void _rkJointSetDis1Hooke(void *prp, int i, double val);
static void _rkJointSetDisHooke(void *prp, double *val);
static void _rkJointSetVelHooke(void *prp, double *val);
static void _rkJointSetAccHooke(void *prp, double *val);
static void _rkJointSetTrqHooke(void *prp, double *val);
static void _rkJointGetDisHooke(void *prp, double *val);
static void _rkJointGetVelHooke(void *prp, double *val);
static void _rkJointGetAccHooke(void *prp, double *val);
static void _rkJointGetTrqHooke(void *prp, double *val);
static void _rkJointGetMotorHooke(void *prp, rkMotor **m);
static void _rkJointCatDisHooke(void *prp, double *dis, double k, double *val);
static void _rkJointSubDisHooke(void *prp, double *dis, double *sdis);
static void _rkJointSetDisCNTHooke(void *prp, double *val, double dt);

static zFrame3D *_rkJointXferHooke(void *prp, zFrame3D *fo, zFrame3D *f);
static void _rkJointIncVelHooke(void *prp, zVec3D *w, zVec6D *vel, zVec6D *acc);
static void _rkJointIncAccHooke(void *prp, zVec6D *acc);

static void _rkJointCalcTrqHooke(void *prp, zVec6D *f);
static void _rkJointTorsionHooke(zFrame3D *dev, zVec6D *t, double dis[]);

static void _rkJointSetFricHooke(void *prp, double *val);
static void _rkJointGetFricHooke(void *prp, double *val);
static void _rkJointGetSFricHooke(void *prp, double *val);
static void _rkJointGetKFricHooke(void *prp, double *val);

static void _rkJointSetRefHooke(void *prp, rkJointRef *ref);
static void _rkJointGetRefHooke(void *prp, rkJointRef *ref);

static zVec3D *_rkJointAngAxisHooke1(void *prp, zFrame3D *f, zVec3D *a);
static zVec3D *_rkJointAngAxisHooke2(void *prp, zFrame3D *f, zVec3D *a);

static bool _rkJointQueryFReadHooke(FILE *fp, char *buf, void *prp, rkMotor *marray, int nm);
static void _rkJointFWriteHooke(FILE *fp, void *prp, char *name);

#define _rkc(p) ((rkJointPrpHooke *)p)

/* limit joint displacement */
double _rkJointLimDis1Hooke(void *prp, int i, double testval){
  testval = zPhaseNormalize( testval );
  return zLimit( testval, _rkc(prp)->min[i], _rkc(prp)->max[i] );
}

void _rkJointLimDisHooke(void *prp, double *testval, double *limval){
  limval[0] = _rkJointLimDis1Hooke( prp, 0, testval[0] );
  limval[1] = _rkJointLimDis1Hooke( prp, 1, testval[1] );
}

/* joint displacement set function */
void _rkJointSetDis1Hooke(void *prp, int i, double val){
  _rkc(prp)->dis[i] = _rkJointLimDis1Hooke( prp, i, val );
  zSinCos( _rkc(prp)->dis[i], &_rkc(prp)->_s[i], &_rkc(prp)->_c[i] );
}

void _rkJointSetDisHooke(void *prp, double *val){
  _rkJointSetDis1Hooke( prp, 0, val[0] );
  _rkJointSetDis1Hooke( prp, 1, val[1] );
}

void _rkJointSetVelHooke(void *prp, double *val){
  memcpy( _rkc(prp)->vel, val, sizeof(double)*2 );
}

void _rkJointSetAccHooke(void *prp, double *val){
  memcpy( _rkc(prp)->acc, val, sizeof(double)*2 );
}

void _rkJointSetTrqHooke(void *prp, double *val){
  memcpy( _rkc(prp)->trq, val, sizeof(double)*2 );
}

/* get joint displacement, velocity, acceleration and torque */
void _rkJointGetDisHooke(void *prp, double *val){
  memcpy( val, _rkc(prp)->dis, sizeof(double)*2 );
}

void _rkJointGetVelHooke(void *prp, double *val){
  memcpy( val, _rkc(prp)->vel, sizeof(double)*2 );
}

void _rkJointGetAccHooke(void *prp, double *val){
  memcpy( val, _rkc(prp)->acc, sizeof(double)*2 );
}

void _rkJointGetTrqHooke(void *prp, double *val){
  memcpy( val, _rkc(prp)->trq, sizeof(double)*2 );
}

/* motor */
void _rkJointGetMotorHooke(void *prp, rkMotor **m){
  *m = &_rkc(prp)->m;
}

void _rkJointCatDisHooke(void *prp, double *dis, double k, double *val)
{
  dis[0] += val[0] * k;
  dis[1] += val[1] * k;
}

void _rkJointSubDisHooke(void *prp, double *dis, double *sdis)
{
  dis[0] -= sdis[0];
  dis[1] -= sdis[1];
}

/* continuously update joint displacement */
void _rkJointSetDisCNTHooke(void *prp, double *val, double dt)
{
  double olddis[2], oldvel[2];

  _rkJointGetDisHooke( prp, olddis );
  _rkJointGetVelHooke( prp, oldvel );
  _rkJointSetDisHooke( prp, val );
  _rkc(prp)->vel[0] = ( val[0] - olddis[0] ) / dt;
  _rkc(prp)->vel[1] = ( val[1] - olddis[1] ) / dt;
  _rkc(prp)->acc[0] = ( _rkc(prp)->vel[0] - oldvel[0] ) / dt;
  _rkc(prp)->acc[1] = ( _rkc(prp)->vel[1] - oldvel[1] ) / dt;
}

/* joint frame transfer function */
zFrame3D *_rkJointXferHooke(void *prp, zFrame3D *fo, zFrame3D *f)
{
  zMat3D m;

  /* position */
  zVec3DCopy( zFrame3DPos(fo), zFrame3DPos(f) );
  /* joint displacements correspond to the rotation angle about
   * z-axis and y-axis, respectively */
  zMat3DZYXSC( &m, _rkc(prp)->_s[0], _rkc(prp)->_c[0], _rkc(prp)->_s[1], _rkc(prp)->_c[1], 0, 1 );
  zMulMatMat3D( zFrame3DAtt(fo), &m, zFrame3DAtt(f) );
  return f;
}

/* joint motion rate transfer function */
void _rkJointIncVelHooke(void *prp, zVec3D *w, zVec6D *vel, zVec6D *acc)
{
  rkJointPrpHooke *p;
  zVec3D v1, v2;
  double dq2;

  p = prp;
  zVec3DCreate( &v1, -p->_s[1]*p->vel[0], p->vel[1], p->_c[1]*p->vel[0] );
  zVec3DOuterProd( w, &v1, &v2 );
  zVec3DAddDRC( zVec6DAng(vel), &v1 );
  zVec3DAddDRC( zVec6DAng(acc), &v2 );
  dq2 = p->vel[0] * p->vel[1];
  zVec3DCreate( &v1, -p->_c[1]*dq2, 0, -p->_s[1]*dq2 );
  zVec3DAddDRC( zVec6DAng(acc), &v1 );
}

void _rkJointIncAccHooke(void *prp, zVec6D *acc)
{
  zVec3D v1;
  rkJointPrpHooke *p;

  p = prp;
  zVec3DCreate( &v1, -p->_s[1]*p->acc[0], p->acc[1], p->_c[1]*p->acc[0] );
  zVec3DAddDRC( zVec6DAng(acc), &v1 );
}

/* joint torque transfer function */
void _rkJointCalcTrqHooke(void *prp, zVec6D *f)
{
  rkJointPrpHooke *p;

  p = prp;
  p->trq[0] =-p->_s[1]*zVec6DElem(f,zXA)+p->_c[1]*zVec6DElem(f,zZA);
  p->trq[1] = zVec6DElem(f,zYA);
}

void _rkJointSetFricHooke(void *prp, double *val)
{
  _rkc(prp)->tf[0] = val[0];
  _rkc(prp)->tf[1] = val[1];
}

void _rkJointGetFricHooke(void *prp, double *val){
  val[0] = _rkc(prp)->tf[0];
  val[1] = _rkc(prp)->tf[1];
}

void _rkJointGetSFricHooke(void *prp, double *val)
{
  val[0] = _rkc(prp)->sf[0];
  val[1] = _rkc(prp)->sf[1];
}

void _rkJointGetKFricHooke(void *prp, double *val)
{
  val[0] = _rkJointRestTrq( _rkc(prp)->stiff[0], _rkc(prp)->viscos[0], _rkc(prp)->coulomb[0], _rkc(prp)->dis[0], _rkc(prp)->vel[0] );
  val[1] = _rkJointRestTrq( _rkc(prp)->stiff[1], _rkc(prp)->viscos[1], _rkc(prp)->coulomb[1], _rkc(prp)->dis[1], _rkc(prp)->vel[1] );
}

void _rkJointSetRefHooke(void *prp, rkJointRef *ref){
  ref[0] = _rkc(prp)->_ref[0];
  ref[1] = _rkc(prp)->_ref[1];
}

void _rkJointGetRefHooke(void *prp, rkJointRef *ref){
  _rkc(prp)->_ref[0] = ref[0];
  _rkc(prp)->_ref[1] = ref[1];
}

/* inverse computation of joint torsion and displacement */
void _rkJointTorsionHooke(zFrame3D *dev, zVec6D *t, double dis[])
{
  zMat3D *r;
  double qx, s, c;

  r = zFrame3DAtt(dev);
  zMulMatTVec3D( r, zFrame3DPos(dev), zVec6DLin(t) );
  qx = atan2( zMat3DElem(r,2,1), zMat3DElem(r,1,1) );
  zMat3DRow( r, 0, zVec6DAng(t) );
  zVec3DMulDRC( zVec6DAng(t), qx );
  zSinCos( qx, &s, &c );
  dis[0] = atan2( -zMat3DElem(r,0,1), c*zMat3DElem(r,1,1)+s*zMat3DElem(r,2,1) );
  c = cos(dis[0]);
  dis[1] = atan2( c*zMat3DElem(r,0,2), c*zMat3DElem(r,0,0) );
}

/* joint axis function */
zVec3D *_rkJointAngAxisHooke1(void *prp, zFrame3D *f, zVec3D *a)
{
  zVec3DMul( zMat3DVec(zFrame3DAtt(f),zX),-_rkc(prp)->_s[1], a );
  return zVec3DCatDRC( a, _rkc(prp)->_c[1], zMat3DVec(zFrame3DAtt(f),zZ) );
}

zVec3D *_rkJointAngAxisHooke2(void *prp, zFrame3D *f, zVec3D *a){
  zVec3DCopy( zMat3DVec(zFrame3DAtt(f),zY), a );
  return a;
}

/* query joint properties */
bool _rkJointQueryFReadHooke(FILE *fp, char *buf, void *prp, rkMotor *marray, int nm)
{
  rkMotor *mp;

  if( strcmp( buf, "dis" ) == 0 ){
    _rkJointSetDis1Hooke( prp, 0, zDeg2Rad(zFDouble(fp)) );
    _rkJointSetDis1Hooke( prp, 1, zDeg2Rad(zFDouble(fp)) );
  } else
  if( strcmp( buf, "min" ) == 0 ){
    _rkc(prp)->min[0] = zDeg2Rad(zFDouble(fp));
    _rkc(prp)->min[1] = zDeg2Rad(zFDouble(fp));
  } else
  if( strcmp( buf, "max" ) == 0 ){
    _rkc(prp)->max[0] = zDeg2Rad(zFDouble(fp));
    _rkc(prp)->max[1] = zDeg2Rad(zFDouble(fp));
  } else
  if( strcmp( buf, "stiff" ) == 0 ){
    _rkc(prp)->stiff[0] = zFDouble(fp);
    _rkc(prp)->stiff[1] = zFDouble(fp);
  } else
  if( strcmp( buf, "viscos" ) == 0 ){
    _rkc(prp)->viscos[0] = zFDouble(fp);
    _rkc(prp)->viscos[1] = zFDouble(fp);
  } else
  if( strcmp( buf, "coulomb" ) == 0 ){
    _rkc(prp)->coulomb[0] = zFDouble(fp);
    _rkc(prp)->coulomb[1] = zFDouble(fp);
  } else
  if( strcmp( buf, "staticfriction" ) == 0 ){
    _rkc(prp)->sf[0] = zFDouble(fp);
    _rkc(prp)->sf[1] = zFDouble(fp);
  }
  if( strcmp( buf, "motor" ) == 0 ){
    zFToken( fp, buf, BUFSIZ );
    zNameFind( marray, nm, buf, mp );
    if( !mp ){
      ZRUNERROR( "invalid motor name %s detected", buf );
      return true;
    }
    if( rkMotorSize(mp) != 2 ){
      ZRUNERROR( "unmatched motor size" );
      return true;
    }
    rkMotorClone( mp, &_rkc(prp)->m );
  } else
  if( !rkMotorQueryFRead( fp, buf, &_rkc(prp)->m ) )
    return false;
  return true;
}

void _rkJointFWriteHooke(FILE *fp, void *prp, char *name)
{
  rkJointPrpHooke *v;

  v = prp;
  if( !zIsTiny( v->dis[0] ) || !zIsTiny( v->dis[1] ) )
    fprintf( fp, "%s: %.10f %.10f\n", name, zRad2Deg(v->dis[0]), zRad2Deg(v->dis[1]) );
  fprintf( fp, "min: %.10f %.10f\n", zRad2Deg(v->min[0]), zRad2Deg(v->min[1]) );
  fprintf( fp, "max: %.10f %.10f\n", zRad2Deg(v->max[0]), zRad2Deg(v->max[1]) );
  fprintf( fp, "stiff: %.10f %.10f\n", zDeg2Rad(v->stiff[0]), zDeg2Rad(v->stiff[1]) );
  fprintf( fp, "viscos: %.10f %.10f\n", zDeg2Rad(v->viscos[0]), zDeg2Rad(v->viscos[1]) );
  fprintf( fp, "coulomb: %.10f %.10f\n", zDeg2Rad(v->coulomb[0]), zDeg2Rad(v->coulomb[1]) );
}

static zVec3D* (*_rk_joint_axis_hooke_ang[])(void*,zFrame3D*,zVec3D*) = {
  _rkJointAngAxisHooke1,
  _rkJointAngAxisHooke2,
};
static zVec3D* (*_rk_joint_axis_hooke_lin[])(void*,zFrame3D*,zVec3D*) = {
  _rkJointAxisNull,
  _rkJointAxisNull,
};
static rkJointCom rk_joint_hooke = {
  2,
  _rkJointLimDisHooke,
  _rkJointSetDisHooke,
  _rkJointSetVelHooke,
  _rkJointSetAccHooke,
  _rkJointSetTrqHooke,
  _rkJointGetDisHooke,
  _rkJointGetVelHooke,
  _rkJointGetAccHooke,
  _rkJointGetTrqHooke,
  _rkJointGetMotorHooke,
  _rkJointCatDisHooke,
  _rkJointSubDisHooke,
  _rkJointSetDisCNTHooke,
  _rkJointXferHooke,
  _rkJointIncVelHooke,
  _rkJointIncAccHooke,
  _rkJointCalcTrqHooke,
  _rkJointTorsionHooke,
  _rkJointSetFricHooke,
  _rkJointGetFricHooke,
  _rkJointGetSFricHooke,
  _rkJointGetKFricHooke,
  _rkJointSetRefHooke,
  _rkJointGetRefHooke,
  _rk_joint_axis_hooke_ang,
  _rk_joint_axis_hooke_lin,
  _rkJointQueryFReadHooke,
  _rkJointFWriteHooke,
};

/* motor */
static byte _rkJointMotorHooke(void *prp);
static void _rkJointMotorSetInputHooke(void *prp, double *val);
static void _rkJointMotorInertiaHooke(void *prp, double *val);
static void _rkJointMotorInputTrqHooke(void *prp, double *val);
static void _rkJointMotorResistanceHooke(void *prp, double *val);
static void _rkJointMotorDrivingTrqHooke(void *prp, double *val);

byte _rkJointMotorHooke(void *prp){
  return rkMotorType( &_rkc(prp)->m );
}
void _rkJointMotorSetInputHooke(void *prp, double *val){
  rkMotorSetInput( &_rkc(prp)->m, val );
}
void _rkJointMotorInertiaHooke(void *prp, double *val){
  rkMotorInertia( &_rkc(prp)->m, val );
}
void _rkJointMotorInputTrqHooke(void *prp, double *val){
  rkMotorInputTrq( &_rkc(prp)->m, val );
}
void _rkJointMotorResistanceHooke(void *prp, double *val){
  rkMotorRegistance( &_rkc(prp)->m, _rkc(prp)->dis, _rkc(prp)->vel, val );
}
void _rkJointMotorDrivingTrqHooke(void *prp, double *val){
  rkMotorDrivingTrq( &_rkc(prp)->m, _rkc(prp)->dis, _rkc(prp)->vel, _rkc(prp)->acc, val );
}

static rkJointMotorCom rk_joint_motor_hooke = {
  _rkJointMotorHooke,
  _rkJointMotorSetInputHooke,
  _rkJointMotorInertiaHooke,
  _rkJointMotorInputTrqHooke,
  _rkJointMotorResistanceHooke,
  _rkJointMotorDrivingTrqHooke,
};

/* ABI */
static void _rkJointABIAxisInertiaHooke(void *prp, zMat6D *m, zMat h);
static void _rkJointABIAddAbiBiosHooke(void *prp, zMat6D *I, zMat6D *J, zVec6D *b, zMat h, zMat6D *pi, zVec6D *pb);
static void _rkJointABIQAccHooke(void *prp, zMat3D *R, zMat6D *I, zVec6D *b, zVec6D *jac, zMat h, zVec6D *acc);

void _rkJointABIAxisInertiaHooke(void *prp, zMat6D *m, zMat h)
{
  rkJointPrpHooke *p;
  zMat3D *m22;

  p = prp;
  m22 = zMat6DMat3D(m, 1, 1);
  zMatElem(h, 0, 0) = zMat3DElem(m22, 0, 0)*p->_s[1]*p->_s[1] - (zMat3DElem(m22, 0, 2) + zMat3DElem(m22, 2, 0))*p->_s[1]*p->_c[1] + zMat3DElem(m22, 2, 2)*p->_c[1]*p->_c[1];
  zMatElem(h, 1, 0) = zMat3DElem(m22, 1, 2)*p->_c[1] - zMat3DElem(m22, 1, 0)*p->_s[1];
  zMatElem(h, 0, 1) = zMat3DElem(m22, 2, 1)*p->_c[1] - zMat3DElem(m22, 0, 1)*p->_s[1];
  zMatElem(h, 1, 1) = zMat3DElem(m22, 1, 1);
}

void _rkJointABIAddAbiBiosHooke(void *prp, zMat6D *I, zMat6D *J, zVec6D *b, zMat h, zMat6D *pi, zVec6D *pb)
{
  eprintf("under construction error: abi update for hooke joint\n");
}

void _rkJointABIQAccHooke(void *prp, zMat3D *R, zMat6D *I, zVec6D *b, zVec6D *jac, zMat h, zVec6D *acc){}

static rkJointABICom rk_joint_abi_hooke = {
  _rkJointABIAxisInertiaHooke,
  _rkJointABIAddAbiBiosHooke,
  _rkJointABIQAccHooke,
};

/* rkJointCreateHooke
 * - create universal joint instance.
 */
rkJoint *rkJointCreateHooke(rkJoint *j)
{
  if( !( j->prp = zAlloc( rkJointPrpHooke, 1 ) ) )
    return NULL;
  _rkc(j->prp)->max[0] = zPI;
  _rkc(j->prp)->min[0] =-zPI;
  _rkc(j->prp)->max[1] = zPI;
  _rkc(j->prp)->min[1] =-zPI;
  rkMotorCreate( &_rkc(j->prp)->m, RK_MOTOR_NONE );
  j->com = &rk_joint_hooke;
  j->mcom = &rk_joint_motor_hooke;
  j->acom = &rk_joint_abi_hooke;
  return j;
}

#undef _rkc
