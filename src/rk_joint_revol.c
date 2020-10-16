/* RoKi - Robot Kinetics library
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * rk_joint_revol - joint structure: revolutional joint
 */

#include <roki/rk_joint.h>

static void _rkJointLimDisRevol(void *prp, double *testval, double *limval);
static void _rkJointSetDisRevol(void *prp, double *val);
static void _rkJointSetVelRevol(void *prp, double *val);
static void _rkJointSetAccRevol(void *prp, double *val);
static void _rkJointSetTrqRevol(void *prp, double *val);
static void _rkJointGetDisRevol(void *prp, double *val);
static void _rkJointGetVelRevol(void *prp, double *val);
static void _rkJointGetAccRevol(void *prp, double *val);
static void _rkJointGetTrqRevol(void *prp, double *val);
static void _rkJointGetMotorRevol(void *prp, rkMotor **m);
static void _rkJointCatDisRevol(void *prp, double *dis, double k, double *val);
static void _rkJointSubDisRevol(void *prp, double *dis, double *sdis);
static void _rkJointSetDisCNTRevol(void *prp, double *val, double dt);
static zFrame3D *_rkJointXferRevol(void *prp, zFrame3D *fo, zFrame3D *f);
static void _rkJointIncVelRevol(void *prp, zVec3D *w, zVec6D *vel, zVec6D *acc);
static void _rkJointIncAccRevol(void *prp, zVec6D *acc);
static void _rkJointCalcTrqRevol(void *prp, zVec6D *f);
static void _rkJointTorsionRevol(zFrame3D *dev, zVec6D *t, double dis[]);
static void _rkJointSetFricRevol(void *prp, double *val);
static void _rkJointGetFricRevol(void *prp, double *val);
static void _rkJointGetSFricRevol(void *prp, double *val);
static void _rkJointGetKFricRevol(void *prp, double *val);
static void _rkJointSetRefRevol(void *prp, rkJointRef *ref);
static void _rkJointGetRefRevol(void *prp, rkJointRef *ref);
static bool _rkJointQueryFReadRevol(FILE *fp, char *buf, void *prp, rkMotor *marray, int nm);
static void _rkJointFWriteRevol(FILE *fp, void *prp, char *name);

#define _rkc(p) ((rkJointPrpRevol *)p)

/* limit joint displacement */
void _rkJointLimDisRevol(void *prp, double *testval, double *limval){
  double angle;

  angle = zPhaseNormalize( *testval );
  *limval = zLimit( angle, _rkc(prp)->min, _rkc(prp)->max );
}

/* set joint displacement, velocity, acceleration and torque */
void _rkJointSetDisRevol(void *prp, double *val){
  _rkJointLimDisRevol( prp, val, &_rkc(prp)->dis );
  zSinCos( _rkc(prp)->dis, &_rkc(prp)->_s, &_rkc(prp)->_c );
}

void _rkJointSetVelRevol(void *prp, double *val){
  _rkc(prp)->vel = *val;
}

void _rkJointSetAccRevol(void *prp, double *val){
  _rkc(prp)->acc = *val;
}

void _rkJointSetTrqRevol(void *prp, double *val){
  _rkc(prp)->trq = *val;
}

/* get joint displacement, velocity, acceleration and torque */
void _rkJointGetDisRevol(void *prp, double *val){
  *val = _rkc(prp)->dis;
}

void _rkJointGetVelRevol(void *prp, double *val){
  *val = _rkc(prp)->vel;
}

void _rkJointGetAccRevol(void *prp, double *val){
  *val = _rkc(prp)->acc;
}

void _rkJointGetTrqRevol(void *prp, double *val){
  *val = _rkc(prp)->trq;
}

/* motor */
void _rkJointGetMotorRevol(void *prp, rkMotor **m){
  *m = &_rkc(prp)->m;
}

void _rkJointCatDisRevol(void *prp, double *dis, double k, double *val)
{
  *dis += k * *val;
}

void _rkJointSubDisRevol(void *prp, double *dis, double *sdis)
{
  *dis -= *sdis;
}

/* continuously update joint displacement */
void _rkJointSetDisCNTRevol(void *prp, double *val, double dt)
{
  double olddis, oldvel;

  _rkJointGetDisRevol( prp, &olddis );
  _rkJointGetVelRevol( prp, &oldvel );
  _rkJointSetDisRevol( prp, val );
  _rkc(prp)->vel = ( *val - olddis ) / dt;
  _rkc(prp)->acc = ( _rkc(prp)->vel - oldvel ) / dt;
}

/* joint frame transfer function */
zFrame3D *_rkJointXferRevol(void *prp, zFrame3D *fo, zFrame3D *f)
{
  /* position */
  zVec3DCopy( zFrame3DPos(fo), zFrame3DPos(f) );
  /* attitude */
  zVec3DMul( zMat3DVec(zFrame3DAtt(fo),0), _rkc(prp)->_c, zMat3DVec(zFrame3DAtt(f),0) );
  zVec3DCatDRC( zMat3DVec(zFrame3DAtt(f),0), _rkc(prp)->_s, zMat3DVec(zFrame3DAtt(fo),1) );
  zVec3DMul( zMat3DVec(zFrame3DAtt(fo),0),-_rkc(prp)->_s, zMat3DVec(zFrame3DAtt(f),1) );
  zVec3DCatDRC( zMat3DVec(zFrame3DAtt(f),1), _rkc(prp)->_c, zMat3DVec(zFrame3DAtt(fo),1) );
  zVec3DCopy( zMat3DVec(zFrame3DAtt(fo),2), zMat3DVec(zFrame3DAtt(f),2) );
  return f;
}

/* joint velocity transfer function */
void _rkJointIncVelRevol(void *prp, zVec3D *w, zVec6D *vel, zVec6D *acc)
{
  zVec6DElem(vel,zZA) += _rkc(prp)->vel;
  zVec6DElem(acc,zXA) += _rkc(prp)->vel * zVec3DElem(w,zY);
  zVec6DElem(acc,zYA) -= _rkc(prp)->vel * zVec3DElem(w,zX);
}

/* joint acceleration transfer function */
void _rkJointIncAccRevol(void *prp, zVec6D *acc)
{
  zVec6DElem(acc,zZA) += _rkc(prp)->acc;
}

/* joint torque transfer function */
void _rkJointCalcTrqRevol(void *prp, zVec6D *f)
{
  _rkc(prp)->trq = zVec6DElem(f,zZA);
}

/* inverse computation of joint torsion and displacement */
void _rkJointTorsionRevol(zFrame3D *dev, zVec6D *t, double dis[])
{
  zMulMatTVec3D( zFrame3DAtt(dev), zFrame3DPos(dev), zVec6DLin(t) );
  dis[0] = rkJointTorsionDisRevol( dev, t );
}

void _rkJointSetFricRevol(void *prp, double *val){
  _rkc(prp)->tf = *val;
}
void _rkJointGetFricRevol(void *prp, double *val){
  *val = _rkc(prp)->tf;
}

void _rkJointGetSFricRevol(void *prp, double *val){
  *val = _rkc(prp)->sf;
}
void _rkJointGetKFricRevol(void *prp, double *val){
  *val = _rkJointRestTrq( _rkc(prp)->stiff, _rkc(prp)->viscos, _rkc(prp)->coulomb, _rkc(prp)->dis, _rkc(prp)->vel );;
}

void _rkJointSetRefRevol(void *prp, rkJointRef *ref){
  _rkc(prp)->_ref = *ref;
}
void _rkJointGetRefRevol(void *prp, rkJointRef *ref){
  *ref = _rkc(prp)->_ref;
}

/* query joint properties */
bool _rkJointQueryFReadRevol(FILE *fp, char *buf, void *prp, rkMotor *marray, int nm)
{
  double val;
  rkMotor *mp;

  if( strcmp( buf, "dis" ) == 0 ){
    val = zDeg2Rad( zFDouble(fp) );
    _rkJointSetDisRevol( prp, &val );
  } else
  if( strcmp( buf, "min" ) == 0 )
    _rkc(prp)->min = zDeg2Rad(zFDouble(fp));
  else
  if( strcmp( buf, "max" ) == 0 )
    _rkc(prp)->max = zDeg2Rad(zFDouble(fp));
  else
  if( strcmp( buf, "stiff" ) == 0 )
    _rkc(prp)->stiff = zFDouble(fp);
  else
  if( strcmp( buf, "viscos" ) == 0 )
    _rkc(prp)->viscos = zFDouble(fp);
  else
  if( strcmp( buf, "coulomb" ) == 0 )
    _rkc(prp)->coulomb = zFDouble(fp);
  else
  if( strcmp( buf, "staticfriction" ) == 0 )
    _rkc(prp)->sf = zFDouble(fp);
  else
  if( strcmp( buf, "motor" ) == 0 ){
    zFToken( fp, buf, BUFSIZ );
    zNameFind( marray, nm, buf, mp );
    if( !mp ){
      ZRUNERROR( "invalid motor name %s detected", buf );
      return true;
    }
    if( rkMotorSize(mp) != 1 ){
      ZRUNERROR( "unmatched motor size" );
      return true;
    }
    rkMotorClone( mp, &_rkc(prp)->m );
  } else
  if( !rkMotorQueryFRead( fp, buf, &_rkc(prp)->m ) )
    return false;
  return true;
}

void _rkJointFWriteRevol(FILE *fp, void *prp, char *name)
{
  rkJointPrpRevol *v;

  v = prp;
  if( !zIsTiny( v->dis ) )
    fprintf( fp, "%s: %.10f\n", name, zRad2Deg(v->dis) );
  fprintf( fp, "min: %.10f\n", zRad2Deg(v->min) );
  fprintf( fp, "max: %.10f\n", zRad2Deg(v->max) );
  fprintf( fp, "stiff: %.10f\n", zDeg2Rad(v->stiff) );
  fprintf( fp, "viscos: %.10f\n", zDeg2Rad(v->viscos) );
}

static zVec3D* (*_rk_joint_axis_revol_ang[])(void*,zFrame3D*,zVec3D*) = {
  _rkJointAxisZ,
};
static zVec3D* (*_rk_joint_axis_revol_lin[])(void*,zFrame3D*,zVec3D*) = {
  _rkJointAxisNull,
};
static rkJointCom rk_joint_revol = {
  1,
  _rkJointLimDisRevol,
  _rkJointSetDisRevol,
  _rkJointSetVelRevol,
  _rkJointSetAccRevol,
  _rkJointSetTrqRevol,
  _rkJointGetDisRevol,
  _rkJointGetVelRevol,
  _rkJointGetAccRevol,
  _rkJointGetTrqRevol,
  _rkJointGetMotorRevol,
  _rkJointCatDisRevol,
  _rkJointSubDisRevol,
  _rkJointSetDisCNTRevol,
  _rkJointXferRevol,
  _rkJointIncVelRevol,
  _rkJointIncAccRevol,
  _rkJointCalcTrqRevol,
  _rkJointTorsionRevol,
  _rkJointSetFricRevol,
  _rkJointGetFricRevol,
  _rkJointGetSFricRevol,
  _rkJointGetKFricRevol,
  _rkJointSetRefRevol,
  _rkJointGetRefRevol,
  _rk_joint_axis_revol_ang,
  _rk_joint_axis_revol_lin,
  _rkJointQueryFReadRevol,
  _rkJointFWriteRevol,
};

/* motor */
static byte _rkJointMotorTypeRevol(void *prp);
static void _rkJointMotorSetInputRevol(void *prp, double *val);
static void _rkJointMotorInertiaRevol(void *prp, double *val);
static void _rkJointMotorInputTrqRevol(void *prp, double *val);
static void _rkJointMotorResistanceRevol(void *prp, double *val);
static void _rkJointMotorDrivingTrqRevol(void *prp, double *val);

byte _rkJointMotorTypeRevol(void *prp){
  return rkMotorType( &_rkc(prp)->m );
}
void _rkJointMotorSetInputRevol(void *prp, double *val){
  rkMotorSetInput( &_rkc(prp)->m, val );
}
void _rkJointMotorInertiaRevol(void *prp, double *val){
  rkMotorInertia( &_rkc(prp)->m, val );
}
void _rkJointMotorInputTrqRevol(void *prp, double *val){
  rkMotorInputTrq( &_rkc(prp)->m, val );
}
void _rkJointMotorResistanceRevol(void *prp, double *val){
  rkMotorRegistance( &_rkc(prp)->m, &_rkc(prp)->dis, &_rkc(prp)->vel, val );
}
void _rkJointMotorDrivingTrqRevol(void *prp, double *val){
  rkMotorDrivingTrq( &_rkc(prp)->m, &_rkc(prp)->dis, &_rkc(prp)->vel, &_rkc(prp)->acc, val );
}

static rkJointMotorCom rk_joint_motor_revol = {
  _rkJointMotorTypeRevol,
  _rkJointMotorSetInputRevol,
  _rkJointMotorInertiaRevol,
  _rkJointMotorInputTrqRevol,
  _rkJointMotorResistanceRevol,
  _rkJointMotorDrivingTrqRevol,
};

/* ABI */
static void _rkJointABIAxisInertiaRevol(void *prp, zMat6D *m, zMat h);
static void _rkJointABIAddAbiBiosRevol(void *prp, zMat6D *I, zMat6D *J, zVec6D *b, zMat h, zMat6D *pi, zVec6D *pb);
static void _rkJointABIQAccRevol(void *prp, zMat3D *R, zMat6D *I, zVec6D *b, zVec6D *jac, zMat h, zVec6D *acc);

void _rkJointABIAxisInertiaRevol(void *prp, zMat6D *m, zMat h)
{
  zMatElem(h, 0, 0) = zMat3DElem(zMat6DMat3D(m, 1, 1), zZ, zZ);
}

void _rkJointABIAddAbiBiosRevol(void *prp, zMat6D *I, zMat6D *J, zVec6D *b, zMat h, zMat6D *pI, zVec6D *pb)
{
  zVec6D tempv, tempv2;
  zMat6D tempm, tempm2;
  double u = 0.0, val;

  /* I */
  zMat3DCol(zMat6DMat3D(I, 0, 1), zZ, zVec6DLin(&tempv));
  zMat3DCol(zMat6DMat3D(I, 1, 1), zZ, zVec6DAng(&tempv));
  zMat3DRow(zMat6DMat3D(I, 1, 0), zZ, zVec6DLin(&tempv2));
  zMat3DRow(zMat6DMat3D(I, 1, 1), zZ, zVec6DAng(&tempv2));
  zVec6DMulDRC(&tempv, -1.0*zMatElem(h, 0, 0));
  zMat6DDyad(&tempv, &tempv2, &tempm);
  zMat6DAddDRC(&tempm, I);

  zMulMatMat6D(J, &tempm, &tempm2);
  zMulMatMatT6D(&tempm2, J, &tempm);

  zMat6DAddDRC(pI, &tempm);

  /* b */
  _rkJointMotorInputTrqRevol( prp, &val );
  u += val;
  _rkJointMotorResistanceRevol( prp, &val );
  u -= val;
  u += _rkc(prp)->tf;
  zVec6DMulDRC(&tempv, u - zVec6DElem(b, zZA));
  zVec6DSub(b, &tempv, &tempv2);

  zMulMat6DVec6D(J, &tempv2, &tempv);
  zVec6DAddDRC(pb, &tempv);
}

void _rkJointABIQAccRevol(void *prp, zMat3D *R, zMat6D *I, zVec6D *b, zVec6D *jac, zMat h, zVec6D *acc)
{
  zVec6D tempv;
  double u, val;

  zMat3DRow(zMat6DMat3D(I, 1, 0), zZ, zVec6DLin(&tempv));
  zMat3DRow(zMat6DMat3D(I, 1, 1), zZ, zVec6DAng(&tempv));
  /* u */
  _rkJointMotorInputTrqRevol( prp, &val );
  u = val;
  _rkJointMotorResistanceRevol( prp, &val );
  u -= val;
  u += _rkc(prp)->tf;
  _rkc(prp)->acc = zMatElem(h, 0, 0)*(u - zVec6DInnerProd(&tempv, jac) - zVec6DElem(b, zZA));

  /* acc */
  zVec6DCopy( jac, acc );
  zVec6DElem( acc, zZA ) += _rkc(prp)->acc;
}

static rkJointABICom rk_joint_abi_revol = {
  _rkJointABIAxisInertiaRevol,
  _rkJointABIAddAbiBiosRevol,
  _rkJointABIQAccRevol,
};

/* rkJointCreateRevol
 * - create revolutional joint instance.
 */
rkJoint *rkJointCreateRevol(rkJoint *j)
{
  if( !( j->prp = zAlloc( rkJointPrpRevol, 1 ) ) )
    return NULL;
  _rkc(j->prp)->max = zPI;
  _rkc(j->prp)->min =-zPI;
  rkMotorCreate( &_rkc(j->prp)->m, RK_MOTOR_NONE );
  j->com = &rk_joint_revol;
  j->mcom = &rk_joint_motor_revol;
  j->acom = &rk_joint_abi_revol;
  return j;
}

#undef _rkc
