/* RoKi - Robot Kinetics library
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * rk_joint_cylin - joint structure: cylindrical joint
 */

#include <roki/rk_joint.h>

static void _rkJointLimDisCylin(void *prp, double *testval, double *limval);
static void _rkJointSetDisCylin(void *prp, double *val);
static void _rkJointSetVelCylin(void *prp, double *val);
static void _rkJointSetAccCylin(void *prp, double *val);
static void _rkJointSetTrqCylin(void *prp, double *val);
static void _rkJointGetDisCylin(void *prp, double *val);
static void _rkJointGetVelCylin(void *prp, double *val);
static void _rkJointGetAccCylin(void *prp, double *val);
static void _rkJointGetTrqCylin(void *prp, double *val);
static void _rkJointGetMotorCylin(void *prp, rkMotor **m);
static void _rkJointCatDisCylin(void *prp, double *dis, double k, double *val);
static void _rkJointSubDisCylin(void *prp, double *dis, double *sdis);
static void _rkJointSetDisCNTCylin(void *prp, double *val, double dt);

static zFrame3D *_rkJointXferCylin(void *prp, zFrame3D *fo, zFrame3D *f);
static void _rkJointIncVelCylin(void *prp, zVec3D *w, zVec6D *vel, zVec6D *acc);
static void _rkJointIncAccCylin(void *prp, zVec6D *acc);

static void _rkJointCalcTrqCylin(void *prp, zVec6D *f);
static void _rkJointTorsionCylin(zFrame3D *dev, zVec6D *t, double dis[]);

static void _rkJointSetFricCylin(void *prp, double *val);
static void _rkJointGetFricCylin(void *prp, double *val);
static void _rkJointGetSFricCylin(void *prp, double *val);
static void _rkJointGetKFricCylin(void *prp, double *val);

static void _rkJointSetRefCylin(void *prp, rkJointRef *ref);
static void _rkJointGetRefCylin(void *prp, rkJointRef *ref);

static bool _rkJointQueryFReadCylin(FILE *fp, char *buf, void *prp, rkMotor *marray, int nm);
static void _rkJointFWriteCylin(FILE *fp, void *prp, char *name);

#define _rkc(p) ((rkJointPrpCylin *)p)

/* limit joint displacement */
void _rkJointLimDisCylin(void *prp, double *testval, double *limval){
  double angle;

  /* 0: prismatic */
  limval[0] = zLimit( testval[0], _rkc(prp)->min[0], _rkc(prp)->max[0] );
  /* 1: revolutional */
  angle = zPhaseNormalize( testval[1] );
  limval[1] = zLimit( angle, _rkc(prp)->min[1], _rkc(prp)->max[1] );
}

/* joint displacement set function */
void _rkJointSetDisCylin(void *prp, double *val){
  _rkJointLimDisCylin( prp, val, _rkc(prp)->dis );
  zSinCos( _rkc(prp)->dis[1], &_rkc(prp)->_s, &_rkc(prp)->_c );
}

void _rkJointSetVelCylin(void *prp, double *val){
  memcpy( _rkc(prp)->vel, val, sizeof(double)*2 );
}

void _rkJointSetAccCylin(void *prp, double *val){
  memcpy( _rkc(prp)->acc, val, sizeof(double)*2 );
}

void _rkJointSetTrqCylin(void *prp, double *val){
  memcpy( _rkc(prp)->trq, val, sizeof(double)*2 );
}

/* get joint displacement, velocity, acceleration and torque */
void _rkJointGetDisCylin(void *prp, double *val){
  memcpy( val, _rkc(prp)->dis, sizeof(double)*2 );
}

void _rkJointGetVelCylin(void *prp, double *val){
  memcpy( val, _rkc(prp)->vel, sizeof(double)*2 );
}

void _rkJointGetAccCylin(void *prp, double *val){
  memcpy( val, _rkc(prp)->acc, sizeof(double)*2 );
}

void _rkJointGetTrqCylin(void *prp, double *val){
  memcpy( val, _rkc(prp)->trq, sizeof(double)*2 );
}

/* motor */
void _rkJointGetMotorCylin(void *prp, rkMotor **m){
  *m = &_rkc(prp)->m;
}

void _rkJointCatDisCylin(void *prp, double *dis, double k, double *val){
  dis[0] += val[0] * k;
  dis[1] += val[1] * k;
}

void _rkJointSubDisCylin(void *prp, double *dis, double *sdis){
  dis[0] -= sdis[0];
  dis[1] -= sdis[1];
}

/* continuously update joint displacement */
void _rkJointSetDisCNTCylin(void *prp, double *val, double dt)
{
  double olddis[2], oldvel[2];

  _rkJointGetDisCylin( prp, olddis );
  _rkJointGetVelCylin( prp, oldvel );
  _rkJointSetDisCylin( prp, val );
  _rkc(prp)->vel[0] = ( val[0] - olddis[0] ) / dt;
  _rkc(prp)->vel[1] = ( val[1] - olddis[1] ) / dt;
  _rkc(prp)->acc[0] = ( _rkc(prp)->vel[0] - oldvel[0] ) / dt;
  _rkc(prp)->acc[1] = ( _rkc(prp)->vel[1] - oldvel[1] ) / dt;
}

/* joint frame transfer function */
zFrame3D *_rkJointXferCylin(void *prp, zFrame3D *fo, zFrame3D *f)
{
  /* rotation */
  zVec3DMul( zMat3DVec(zFrame3DAtt(fo),0), _rkc(prp)->_c, zMat3DVec(zFrame3DAtt(f),0) );
  zVec3DCatDRC( zMat3DVec(zFrame3DAtt(f),0), _rkc(prp)->_s, zMat3DVec(zFrame3DAtt(fo),1) );
  zVec3DMul( zMat3DVec(zFrame3DAtt(fo),0),-_rkc(prp)->_s, zMat3DVec(zFrame3DAtt(f),1) );
  zVec3DCatDRC( zMat3DVec(zFrame3DAtt(f),1), _rkc(prp)->_c, zMat3DVec(zFrame3DAtt(fo),1) );
  zVec3DCopy( zMat3DVec(zFrame3DAtt(fo),2), zMat3DVec(zFrame3DAtt(f),2) );
  /* slide */
  zVec3DCat( zFrame3DPos(fo),
    _rkc(prp)->dis[0], zMat3DVec(zFrame3DAtt(fo),zZ), zFrame3DPos(f) );
  return f;
}

/* joint velocity transfer function */
void _rkJointIncVelCylin(void *prp, zVec3D *w, zVec6D *vel, zVec6D *acc)
{
  double xa, ya;

  xa = zVec3DElem(w,zX);
  ya = zVec3DElem(w,zY);
  zVec6DElem(vel,zZ) += _rkc(prp)->vel[0];
  zVec6DElem(vel,zZA) += _rkc(prp)->vel[1];
  zVec6DElem(acc,zX) += 2 * _rkc(prp)->vel[0] * ya;
  zVec6DElem(acc,zY) -= 2 * _rkc(prp)->vel[0] * xa;
  zVec6DElem(acc,zXA) += _rkc(prp)->vel[1] * ya;
  zVec6DElem(acc,zYA) -= _rkc(prp)->vel[1] * xa;
}

/* joint acceleration transfer function */
void _rkJointIncAccCylin(void *prp, zVec6D *acc)
{
  zVec6DElem(acc,zZ) += _rkc(prp)->acc[0];
  zVec6DElem(acc,zZA) += _rkc(prp)->acc[1];
}

/* joint torque transfer function */
void _rkJointCalcTrqCylin(void *prp, zVec6D *f)
{
  _rkc(prp)->trq[0] = zVec6DElem(f,zZ);
  _rkc(prp)->trq[1] = zVec6DElem(f,zZA);
}

/* inverse computation of joint torsion and displacement */
void _rkJointTorsionCylin(zFrame3D *dev, zVec6D *t, double dis[])
{
  dis[0] = rkJointTorsionDisPrism( dev, t );
  dis[1] = rkJointTorsionDisRevol( dev, t );
}

void _rkJointSetFricCylin(void *prp, double *val)
{
  _rkc(prp)->tf[0] = val[0];
  _rkc(prp)->tf[1] = val[1];
}
void _rkJointGetFricCylin(void *prp, double *val){
  val[0] = _rkc(prp)->tf[0];
  val[1] = _rkc(prp)->tf[1];
}

void _rkJointGetSFricCylin(void *prp, double *val)
{
  val[0] = _rkc(prp)->sf[0];
  val[1] = _rkc(prp)->sf[1];
}

void _rkJointGetKFricCylin(void *prp, double *val)
{
  val[0] = _rkJointRestTrq( _rkc(prp)->stiff[0], _rkc(prp)->viscos[0], _rkc(prp)->coulomb[0], _rkc(prp)->dis[0], _rkc(prp)->vel[0] );
  val[1] = _rkJointRestTrq( _rkc(prp)->stiff[1], _rkc(prp)->viscos[1], _rkc(prp)->coulomb[1], _rkc(prp)->dis[1], _rkc(prp)->vel[1] );
}

void _rkJointSetRefCylin(void *prp, rkJointRef *ref){
  ref[0] = _rkc(prp)->_ref[0];
  ref[1] = _rkc(prp)->_ref[1];
}
void _rkJointGetRefCylin(void *prp, rkJointRef *ref){
  _rkc(prp)->_ref[0] = ref[0];
  _rkc(prp)->_ref[1] = ref[1];
}

/* query joint properties */
bool _rkJointQueryFReadCylin(FILE *fp, char *buf, void *prp, rkMotor *marray, int nm)
{
  rkMotor *mp;
  double dis[2];

  if( strcmp( buf, "dis" ) == 0 ){
    dis[0] = zFDouble(fp);
    dis[1] = zDeg2Rad(zFDouble(fp));
    _rkJointSetDisCylin( prp, dis );
  } else
  if( strcmp( buf, "min" ) == 0 ){
    _rkc(prp)->min[0] = zFDouble(fp);
    _rkc(prp)->min[1] = zDeg2Rad(zFDouble(fp));
  } else
  if( strcmp( buf, "max" ) == 0 ){
    _rkc(prp)->max[0] = zFDouble(fp);
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

void _rkJointFWriteCylin(FILE *fp, void *prp, char *name)
{
  rkJointPrpCylin *v;

  v = prp;
  if( !zIsTiny( v->dis[0] ) || !zIsTiny( v->dis[1] ) )
    fprintf( fp, "%s: %.10f %.10f\n", name, v->dis[0], zRad2Deg(v->dis[1]) );
  fprintf( fp, "min: %.10f %.10f\n", v->min[0], zRad2Deg(v->min[1]) );
  fprintf( fp, "max: %.10f %.10f\n", v->max[0], zRad2Deg(v->max[1]) );
  fprintf( fp, "stiff: %.10f %.10f\n", v->stiff[0], zDeg2Rad(v->stiff[1]) );
  fprintf( fp, "viscos: %.10f %.10f\n", v->viscos[0], zDeg2Rad(v->viscos[1]) );
  fprintf( fp, "coulomb: %.10f %.10f\n", v->coulomb[0], zDeg2Rad(v->coulomb[1]) );
}

static zVec3D* (*_rk_joint_axis_cylin_ang[])(void*,zFrame3D*,zVec3D*) = {
  _rkJointAxisNull,
  _rkJointAxisZ,
};
static zVec3D* (*_rk_joint_axis_cylin_lin[])(void*,zFrame3D*,zVec3D*) = {
  _rkJointAxisZ,
  _rkJointAxisNull,
};
static rkJointCom rk_joint_cylin = {
  2,
  _rkJointLimDisCylin,
  _rkJointSetDisCylin,
  _rkJointSetVelCylin,
  _rkJointSetAccCylin,
  _rkJointSetTrqCylin,
  _rkJointGetDisCylin,
  _rkJointGetVelCylin,
  _rkJointGetAccCylin,
  _rkJointGetTrqCylin,
  _rkJointGetMotorCylin,
  _rkJointCatDisCylin,
  _rkJointSubDisCylin,
  _rkJointSetDisCNTCylin,
  _rkJointXferCylin,
  _rkJointIncVelCylin,
  _rkJointIncAccCylin,
  _rkJointCalcTrqCylin,
  _rkJointTorsionCylin,
  _rkJointSetFricCylin,
  _rkJointGetFricCylin,
  _rkJointGetSFricCylin,
  _rkJointGetKFricCylin,
  _rkJointSetRefCylin,
  _rkJointGetRefCylin,
  _rk_joint_axis_cylin_ang,
  _rk_joint_axis_cylin_lin,
  _rkJointQueryFReadCylin,
  _rkJointFWriteCylin,
};

/* motor */
static byte _rkJointMotorCylin(void *prp);
static void _rkJointMotorSetInputCylin(void *prp, double *val);
static void _rkJointMotorInertiaCylin(void *prp, double *val);
static void _rkJointMotorInputTrqCylin(void *prp, double *val);
static void _rkJointMotorResistanceCylin(void *prp, double *val);
static void _rkJointMotorDrivingTrqCylin(void *prp, double *val);

byte _rkJointMotorCylin(void *prp){
  return rkMotorType( &_rkc(prp)->m );
}
void _rkJointMotorSetInputCylin(void *prp, double *val){
  rkMotorSetInput( &_rkc(prp)->m, val );
}
void _rkJointMotorInertiaCylin(void *prp, double *val){
  rkMotorInertia( &_rkc(prp)->m, val );
}
void _rkJointMotorInputTrqCylin(void *prp, double *val){
  rkMotorInputTrq( &_rkc(prp)->m, val );
}
void _rkJointMotorResistanceCylin(void *prp, double *val){
  rkMotorRegistance( &_rkc(prp)->m, _rkc(prp)->dis, _rkc(prp)->vel, val );
}
void _rkJointMotorDrivingTrqCylin(void *prp, double *val){
  rkMotorDrivingTrq( &_rkc(prp)->m, _rkc(prp)->dis, _rkc(prp)->vel, _rkc(prp)->acc, val );
}

static rkJointMotorCom rk_joint_motor_cylin = {
  _rkJointMotorCylin,
  _rkJointMotorSetInputCylin,
  _rkJointMotorInertiaCylin,
  _rkJointMotorInputTrqCylin,
  _rkJointMotorResistanceCylin,
  _rkJointMotorDrivingTrqCylin,
};

/* ABI */
static void _rkJointABIAxisInertiaCylin(void *prp, zMat6D *m, zMat h);
static void _rkJointABIAddAbiBiosCylin(void *prp, zMat6D *I, zMat6D *J, zVec6D *b, zMat h, zMat6D *pi, zVec6D *pb);
static void _rkJointABIQAccCylin(void *prp, zMat3D *R, zMat6D *I, zVec6D *b, zVec6D *jac, zMat h, zVec6D *acc);

void _rkJointABIAxisInertiaCylin(void *prp, zMat6D *m, zMat h)
{
  zMatElem(h, 0, 0) = zMat3DElem(zMat6DMat3D(m, 0, 0), zZ, zZ);
  zMatElem(h, 1, 0) = zMat3DElem(zMat6DMat3D(m, 1, 0), zZ, zZ);
  zMatElem(h, 0, 1) = zMat3DElem(zMat6DMat3D(m, 0, 1), zZ, zZ);
  zMatElem(h, 1, 1) = zMat3DElem(zMat6DMat3D(m, 1, 1), zZ, zZ);
}

void _rkJointABIAddAbiBiosCylin(void *prp, zMat6D *I, zMat6D *J, zVec6D *b, zMat h, zMat6D *pI, zVec6D *pb)
{
  zVec6D v13, v31, v16, v61, tempv;
  zMat6D tempm, tempm2;
  double u[2], val[2];

  /* I */
  zMat3DCol(zMat6DMat3D(I, 0, 0), zZ, zVec6DLin(&v13));
  zMat3DCol(zMat6DMat3D(I, 1, 0), zZ, zVec6DAng(&v13));
  zMat3DCol(zMat6DMat3D(I, 1, 0), zZ, zVec6DLin(&v16));
  zMat3DCol(zMat6DMat3D(I, 1, 1), zZ, zVec6DAng(&v16));
  zMat3DRow(zMat6DMat3D(I, 0, 0), zZ, zVec6DLin(&v31));
  zMat3DRow(zMat6DMat3D(I, 0, 1), zZ, zVec6DAng(&v31));
  zMat3DRow(zMat6DMat3D(I, 1, 0), zZ, zVec6DLin(&v61));
  zMat3DRow(zMat6DMat3D(I, 1, 1), zZ, zVec6DAng(&v61));

  zVec6DMul(&v13, zMatElem(h, 0, 0), &tempv);
  zMat6DDyad(&tempv, &v31, &tempm);
  zVec6DMul(&v13, zMatElem(h, 0, 1), &tempv);
  zMat6DDyad(&tempv, &v61, &tempm2);
  zMat6DAddDRC(&tempm, &tempm2);
  zVec6DMul(&v16, zMatElem(h, 1, 0), &tempv);
  zMat6DDyad(&tempv, &v31, &tempm2);
  zMat6DAddDRC(&tempm, &tempm2);
  zVec6DMul(&v16, zMatElem(h, 1, 1), &tempv);
  zMat6DDyad(&tempv, &v61, &tempm2);
  zMat6DAddDRC(&tempm, &tempm2);

  zMulMatMat6D(J, &tempm, &tempm2);
  zMulMatMatT6D(&tempm2, J, &tempm);

  zMat6DAddDRC(pI, &tempm);

  /* b */
  _rkJointMotorInputTrqCylin( prp, val );
  u[0] = val[0];
  u[1] = val[1];
  _rkJointMotorResistanceCylin( prp, val );
  u[0] += val[0];
  u[1] += val[1];
  u[0] += _rkc(prp)->tf[0];
  u[1] += _rkc(prp)->tf[1];
  zVec6DCat(b,         (u[0] - zVec6DElem(b, 2))*zMatElem(h, 0, 0) + (u[1] - zVec6DElem(b, 5))*zMatElem(h, 0, 1), &v13, &tempv);
  zVec6DCatDRC(&tempv, (u[0] - zVec6DElem(b, 2))*zMatElem(h, 1, 0) + (u[1] - zVec6DElem(b, 5))*zMatElem(h, 1, 1), &v16);

  zMulMat6DVec6D(J, &tempv, &v13); /* v13:temp */
  zVec6DAddDRC(pb, &v13);
}

void _rkJointABIQAccCylin(void *prp, zMat3D *R, zMat6D *I, zVec6D *b, zVec6D *jac, zMat h, zVec6D *acc)
{
  double u[2], val[2];
  zVec6D v31, v61;

  zMat3DRow(zMat6DMat3D(I, 0, 0), zZ, zVec6DLin(&v31));
  zMat3DRow(zMat6DMat3D(I, 0, 1), zZ, zVec6DAng(&v31));
  zMat3DRow(zMat6DMat3D(I, 1, 0), zZ, zVec6DLin(&v61));
  zMat3DRow(zMat6DMat3D(I, 1, 1), zZ, zVec6DAng(&v61));
  /* u */
  _rkJointMotorInputTrqCylin( prp, val );
  u[0] = val[0];
  u[1] = val[1];
  _rkJointMotorResistanceCylin( prp, val );
  u[0] += val[0];
  u[1] += val[1];
  u[0] += _rkc(prp)->tf[0];
  u[1] += _rkc(prp)->tf[1];

  u[0] -= zVec6DInnerProd( &v31, jac ) + zVec6DElem( b, zZ );
  u[1] -= zVec6DInnerProd( &v61, jac ) + zVec6DElem( b, zZA );

  _rkc(prp)->acc[0] = u[0]*zMatElem(h, 0, 0) + u[1]*zMatElem(h, 0, 1);
  _rkc(prp)->acc[1] = u[0]*zMatElem(h, 1, 0) + u[1]*zMatElem(h, 1, 1);

  /* acc */
  zVec6DCopy( jac, acc );
  zVec6DElem( acc, zZ )  += _rkc(prp)->acc[0];
  zVec6DElem( acc, zZA ) += _rkc(prp)->acc[1];
}

static rkJointABICom rk_joint_abi_cylin = {
  _rkJointABIAxisInertiaCylin,
  _rkJointABIAddAbiBiosCylin,
  _rkJointABIQAccCylin,
};

/* rkJointCreateCylin
 * - create cylindrical joint instance.
 */
rkJoint *rkJointCreateCylin(rkJoint *j)
{
  if( !( j->prp = zAlloc( rkJointPrpCylin, 1 ) ) )
    return NULL;
  _rkc(j->prp)->max[0] = HUGE_VAL;
  _rkc(j->prp)->min[0] =-HUGE_VAL;
  _rkc(j->prp)->max[1] = zPI;
  _rkc(j->prp)->min[1] =-zPI;
  rkMotorCreate( &_rkc(j->prp)->m, RK_MOTOR_NONE );
  j->com = &rk_joint_cylin;
  j->mcom = &rk_joint_motor_cylin;
  j->acom = &rk_joint_abi_cylin;
  return j;
}

#undef _rkc
