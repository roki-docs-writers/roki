/* RoKi - Robot Kinetics library
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * rk_fd - forward dynamics computation
 * contributer: 2014-2015 Naoki Wakisaka
 */

/* ********************************************************** */
/* CLASS: rkFD
   forward dynamics class
 * ********************************************************** */

#include <roki/rk_fd.h>

#define RK_FD_DT_DEFAULT           0.001
#define RK_FD_CDT_DEFAULT          0.002
#define RK_FD_JOINT_COMP_K_DEFAULT 100
#define RK_FD_JOINT_COMP_L_DEFAULT 0.01
#define RK_FD_KINETIC_FRIC_WEIGHT_DEFAULT 100
#define RK_FD_FRIC_PYRAMID_ORDER_DEFAULT 8

#define rkChainUpdateRateGrav(c) \
  rkLinkUpdateRate( rkChainRoot(c), ZVEC6DZERO, ZVEC6DZERO )

void _rkFDCellDatFree(rkFDCellDat *ld)
{
  rkChainDestroy( &ld->chain );
}

rkFD *rkFDCreate(rkFD *fd)
{
  zListInit( &fd->list );
  zArrayInit( &fd->ci );
  fd->size = 0;
  fd->dis = NULL;
  fd->vel = NULL;
  fd->acc = NULL;
  fd->_comp_k = RK_FD_JOINT_COMP_K_DEFAULT;
  fd->_comp_l = RK_FD_JOINT_COMP_L_DEFAULT;
  fd->_kf_weight = RK_FD_KINETIC_FRIC_WEIGHT_DEFAULT;
  fd->_pyramid = RK_FD_FRIC_PYRAMID_ORDER_DEFAULT;
  /* collision detector */
  rkCDCreate( &fd->cd );
  /* ODE solver */
  rkFDODE2Assign( fd, Regular );
  rkFDODE2AssignRegular( fd, RKG );
  fd->_ode_step = 0;
  fd->dt = RK_FD_DT_DEFAULT;
  fd->t = 0.0;
  /* default collision solver */
  rkFDSetSolver( fd, Vert );
  return fd;
}

void rkFDDestroy(rkFD *fd)
{
  rkFDCell *lc;

  zVecFreeAO( 3, fd->dis, fd->vel, fd->acc );
  rkCDDestroy( &fd->cd );
  rkContactInfoPoolDestroy( &fd->ci );
  while( !zListIsEmpty( &fd->list ) ){
    zListDeleteHead( &fd->list, &lc );
    _rkFDCellDatFree( &lc->data );
    zFree( lc );
  }
}

void _rkFDCellDatSetOffset(rkFD *fd, rkFDCellDat *ld, int offset)
{
  ld->_offset = offset;
  zVecBuf(&ld->_dis) = &zVecElem(fd->dis,offset);
  zVecBuf(&ld->_vel) = &zVecElem(fd->vel,offset);
}

bool _rkFDAllocJointState(rkFD *fd, rkFDCell *rlc)
{
  rkFDCell *lc;
  zVec pdis=NULL, pvel=NULL;
  int offset = 0;
  int js, rjs;

  if( fd->dis && fd->vel )
    if( !(pdis = zVecClone( fd->dis )) || !(pvel = zVecClone( fd->vel )) )
      goto ERROR;

  zVecFreeAO( 3, fd->dis, fd->vel, fd->acc );

  rjs = rkChainJointSize( &rlc->data.chain );

  if( !( fd->dis = zVecAlloc( fd->size + rjs ) ) ||
      !( fd->vel = zVecAlloc( fd->size + rjs ) ) ||
      !( fd->acc = zVecAlloc( fd->size + rjs ) ) )
    goto ERROR;

  zListForEach( &fd->list, lc ){
    if( lc ==  rlc ) break;
    js = rkChainJointSize( &lc->data.chain );
    if( pdis ) zRawVecCopy( &zVecElem(pdis,offset), &zVecElem(fd->dis,offset), js );
    if( pvel ) zRawVecCopy( &zVecElem(pvel,offset), &zVecElem(fd->vel,offset), js );
    _rkFDCellDatSetOffset( fd, &lc->data, offset );
    offset += js;
  }
  zRawVecClear( &zVecElem(fd->dis,offset), rjs );
  zRawVecClear( &zVecElem(fd->vel,offset), rjs );
  zVecSetSize( &lc->data._dis, rjs );
  zVecSetSize( &lc->data._vel, rjs );
  _rkFDCellDatSetOffset( fd, &lc->data, offset );
  offset += rjs;
  for( lc=zListCellNext(lc); lc!=zListRoot(&fd->list); lc=zListCellNext(lc) ){
    js = rkChainJointSize( &lc->data.chain );
    if( pdis ) zRawVecCopy( &zVecElem(pdis,offset), &zVecElem(fd->dis,offset+rjs), js );
    if( pvel ) zRawVecCopy( &zVecElem(pvel,offset), &zVecElem(fd->vel,offset+rjs), js );
    _rkFDCellDatSetOffset( fd, &lc->data, offset + rjs );
    offset += js;
  }

  zVecFreeAO( 2, pdis, pvel );
  return true;
 ERROR:
  zVecFreeAO( 4, pdis, pvel, fd->dis, fd->vel );
  return false;
}

rkFDCellDat *_rkFDCellDatJointRefInit(rkFDCellDat *ld)
{
  register int i, j;
  rkJoint *joint;
  double val[6];
  rkJointRef jref[6];

  for( i=0; i<rkChainNum(&ld->chain); i++ ){
    joint = rkChainLinkJoint(&ld->chain,i);
    rkJointGetDis( joint, val );
    for( j=0; j<rkJointSize(joint); j++ ){
      jref[j].ref_dis = val[j];
      jref[j].ref_trq = 0.0;
      jref[j].type = RK_CONTACT_SF;
    }
    rkJointSetRef( joint, jref );
  }
  return ld;
}

/* for braking joint */
void _rkFDChainSetBreakJointFuncFix(rkChain *chain)
{
  register int i;

  for( i=0; i<rkChainNum(chain); i++ ){
    if( rkLinkBreakPrp(rkChainLink(chain,i)) ){
      rkJointSetFuncFix( rkChainLinkJoint(chain,i) );
    }
  }
}

rkFDCellDat *_rkFDCellDatInit(rkFDCellDat *ld)
{
  rkChainABIInit( &ld->chain );
  _rkFDChainSetBreakJointFuncFix( &ld->chain );
  ld->_dis.size = ld->_vel.size = ld->_acc.size = rkChainJointSize( &ld->chain );
  ld->_dis.elem = ld->_vel.elem = ld->_acc.elem = NULL;
  _rkFDCellDatJointRefInit( ld );
  return ld;
}

rkFDCell *_rkFDCellPush(rkFD *fd, rkFDCell *lc)
{
  rkCDPair *cdp;
  _rkFDCellDatInit( &lc->data );
  zListInsertHead( &fd->list, lc );
  fd->size += rkChainJointSize( &lc->data.chain );
  if( !_rkFDAllocJointState( fd, lc ) ){
    rkFDDestroy( fd );
    return NULL;
  }
  /* cd reg */
  rkCDChainReg( &fd->cd, &lc->data.chain, true );
  /* set contact info */
  zListForEach( &fd->cd.plist, cdp ){
    if( cdp->data.ci == NULL ){
      cdp->data.ci = rkContactInfoPoolAssoc( &fd->ci, rkLinkStuff(cdp->data.cell[0]->data.link),
                                                      rkLinkStuff(cdp->data.cell[1]->data.link) );
      if( cdp->data.ci == NULL )
        cdp->data.ci = &fd->_ci_def;
    }
  }
  return lc;
}

rkFDCell *rkFDChainReg(rkFD *fd, rkChain *chain)
{
  rkFDCell *lc;

  if( !chain ) return NULL;
  lc = zAlloc( rkFDCell, 1 );
  if( !rkChainClone( chain, &lc->data.chain ) ){
    zFree( lc );
    return NULL;
  }
  return _rkFDCellPush( fd, lc ) ? lc : NULL;
}

rkFDCell *rkFDChainRegFile(rkFD *fd, char filename[])
{
  rkFDCell *lc;

  lc = zAlloc( rkFDCell, 1 );
  if( !rkChainReadFile( &lc->data.chain, filename ) ){
    _rkFDCellDatFree( &lc->data );
    zFree( lc );
    return NULL;
  }
  return _rkFDCellPush( fd, lc ) ? lc : NULL;
}

bool rkFDChainUnreg(rkFD *fd, rkFDCell *cell)
{
  rkFDCell *lc;

  zListForEach( &fd->list, lc ){
    if( lc != cell ) continue;
    if( !_rkFDAllocJointState( fd, lc ) ){
      rkFDDestroy( fd );
      return false;
    }
    zListPurge( &fd->list, lc );
    _rkFDCellDatFree( &lc->data );
    zFree( lc );
    return true;
  }
  return false;
}

/******************************************************************************/
/* read contact information */
bool rkFDContactInfoReadFile(rkFD *fd, char filename[]){
  rkCDPair *cdp;

  if( zArrayNum( &fd->ci ) != 0)
    rkContactInfoPoolDestroy( &fd->ci );
  if( !rkContactInfoPoolReadFile( &fd->ci, filename ) ) return false;
  /* set contact info */
  zListForEach( &fd->cd.plist, cdp ){
    cdp->data.ci = rkContactInfoPoolAssoc( &fd->ci, rkLinkStuff(cdp->data.cell[0]->data.link),
                                                    rkLinkStuff(cdp->data.cell[1]->data.link) );
    if( cdp->data.ci == NULL )
      cdp->data.ci = &fd->_ci_def;
  }
  return true;
}

/******************************************************************************/
/* set state function */
void rkFDChainSetDis(rkFDCell *lc, zVec dis)
{
  zVecCopy( dis, &lc->data._dis );
  rkChainSetJointDisAll( &lc->data.chain, dis );
}

void rkFDChainSetVel(rkFDCell *lc, zVec vel)
{
  zVecCopy( vel, &lc->data._vel );
  rkChainSetJointVelAll( &lc->data.chain, vel );
}

/* connect each cell state to the total state */
void _rkFDConnectJointState(rkFD *fd, zVec dis, zVec vel, zVec acc)
{
  rkFDCell *lc;

  zListForEach( &fd->list, lc ){
    lc->data._dis.elem = &zVecElem(dis,lc->data._offset);
    lc->data._vel.elem = &zVecElem(vel,lc->data._offset);
    lc->data._acc.elem = &zVecElem(acc,lc->data._offset);
    rkChainFK( &lc->data.chain, &lc->data._dis );
    rkChainSetJointVelAll( &lc->data.chain, &lc->data._vel );
    rkChainUpdateVel( &lc->data.chain );
  }
}

/******************************************************************************/
/* ODE solver default function */
zVec rkFDODECatDefault(zVec x, double k, zVec v, zVec xnew, void *util)
{
  rkFDCell *lc;
  zVecStruct lv, lxn;

  zVecCopyNC( x, xnew );
  zListForEach( &((rkFD *)util)->list, lc ){
    lv.size  = rkChainJointSize( &lc->data.chain );
    lxn.size = rkChainJointSize( &lc->data.chain );
    lv.elem  = &zVecElem(v,   lc->data._offset);
    lxn.elem = &zVecElem(xnew,lc->data._offset);
    rkChainCatJointDisAll( &lc->data.chain, &lxn, k, &lv );
  }
  return xnew;
}

zVec rkFDODESubDefault(zVec x1, zVec x2, zVec dx, void *util)
{
  rkFDCell *lc;
  zVecStruct lx2, ldx;

  zVecCopyNC( x1, dx );
  zListForEach( &((rkFD *)util)->list, lc ){
    lx2.size = rkChainJointSize( &lc->data.chain );
    ldx.size = rkChainJointSize( &lc->data.chain );
    lx2.elem = &zVecElem(x2,lc->data._offset);
    ldx.elem = &zVecElem(dx,lc->data._offset);
    rkChainSubJointDisAll( &lc->data.chain, &ldx, &lx2 );
  }
  return dx;
}

/******************************************************************************/
/* forward kinematics */
/* NOTE:
 * this function doesn't update the total state in rkFD
 * it may be better to change the specification
 */
void rkFDFK(rkFD *fd, zVec dis)
{
  rkFDCell *lc;
  zVecStruct ldis;

  zListForEach( &fd->list, lc ){
    ldis.size = rkChainJointSize( &lc->data.chain );
    ldis.elem = &zVecElem(dis,lc->data._offset);
    rkChainFK( &lc->data.chain, &ldis );
  }
}

void rkFDUpdateRate(rkFD *fd, zVec vel, zVec acc)
{
  rkFDCell *lc;
  zVecStruct lvel, lacc;

  zListForEach( &fd->list, lc ){
    lvel.size = lacc.size = rkChainJointSize( &lc->data.chain );
    lvel.elem = &zVecElem(vel,lc->data._offset);
    lacc.elem = &zVecElem(acc,lc->data._offset);
    rkChainSetJointRateAll( &lc->data.chain, &lvel, &lacc );
    rkChainUpdateRateGrav( &lc->data.chain );
  }
}

void _rkFDChainUpdateFKRate(rkFDCell *lc)
{
  rkChainFK( &lc->data.chain, &lc->data._dis );
  rkChainSetJointRateAll( &lc->data.chain, &lc->data._vel, &lc->data._acc );
  rkChainUpdateRateGrav( &lc->data.chain );
}

void rkFDUpdateFKRate(rkFD *fd)
{
  rkFDCell *lc;

  zListForEach( &fd->list, lc ){
    _rkFDChainUpdateFKRate( lc );
  }
}

/******************************************************************************/
/* for a fake-crawler */
void rkFDCDCellSetSlideMode(rkCDCell *cell, bool mode)
{
  cell->data.slide_mode = mode;
}

void rkFDCDCellSetSlideVel(rkCDCell *cell, double vel)
{
  cell->data.slide_vel = vel;
}

void rkFDCDCellSetSlideAxis(rkCDCell *cell, zVec3D *axis)
{
  zVec3DCopy( axis, &cell->data.slide_axis );
}

rkCDCell *rkFDShape3DGetCDCell(rkFD *fd, zShape3D *shape)
{
  rkCDCell *cell;

  zListForEach( &fd->cd.clist, cell )
    if( cell->data.shape == shape ) return cell;
  return NULL;
}

rkCDCell *rkFDShape3DSetSlideMode(rkFD *fd, zShape3D *shape, bool mode)
{
  rkCDCell *cell;

  if( ( cell = rkFDShape3DGetCDCell( fd, shape ) ) == NULL )
    return NULL;
  rkFDCDCellSetSlideMode( cell, mode );
  return cell;
}

rkCDCell *rkFDShape3DSetSlideVel(rkFD *fd, zShape3D *shape, double vel)
{
  rkCDCell *cell;

  if( ( cell = rkFDShape3DGetCDCell( fd, shape ) ) == NULL )
    return NULL;
  rkFDCDCellSetSlideVel( cell, vel );
  return cell;
}

rkCDCell *rkFDShape3DSetSlideAxis(rkFD *fd, zShape3D *shape, zVec3D *axis)
{
  rkCDCell *cell;

  if( ( cell = rkFDShape3DGetCDCell( fd, shape ) ) == NULL )
    return NULL;
  rkFDCDCellSetSlideAxis( cell, axis );
  return cell;
}

/******************************************************************************/
/* relative velocity */
zVec3D *_rkFDLinkPointWldVel(rkLink *link, zVec3D *p, zVec3D *v)
{
  zVec6D v6;
  zVec3D tempv;

  zMulMatVec6D( rkLinkWldAtt(link), rkLinkVel(link), &v6 );
  zVec3DSub( p, rkLinkWldPos(link) ,&tempv );
  zVec3DOuterProd( zVec6DAng(&v6), &tempv, &tempv );
  zVec3DAdd( zVec6DLin(&v6), &tempv, v );
  return v;
}

void _rkFDLinkAddSlideVel(rkCDCell *cell, zVec3D *p, zVec3D *n, zVec3D *v)
{
  zVec3D sv, tmpv;

  zVec3DSub( p, rkLinkWldPos(cell->data.link), &tmpv );
  zMulMatVec3D( rkLinkWldAtt(cell->data.link), &cell->data.slide_axis, &sv );
  zVec3DOuterProd( &sv, &tmpv, &sv );
  zVec3DCatDRC( &sv, -zVec3DInnerProd( &sv, n ), n );
  if( zIsTiny( zVec3DNorm( &sv ) ) ){
    zVec3DClear( &sv );
  }else{
    zVec3DMulDRC( &sv, cell->data.slide_vel/zVec3DNorm( &sv ) );
  }
  zVec3DAddDRC( v, &sv );
}

zVec3D *_rkFDChainPointRelativeVel(rkCDPair *cdp, zVec3D *p, zVec3D *n, rkCDCell *cell, zVec3D *v)
{
  zVec3D vv[2];
  register int i;

  for( i=0; i<2; i++ ){
    if( cdp->data.cell[i]->data.type == RK_CD_CELL_STAT )
      zVec3DClear( &vv[i] );
    else
      _rkFDLinkPointWldVel( cdp->data.cell[i]->data.link, p, &vv[i] );
    if( cdp->data.cell[i]->data.slide_mode )
      _rkFDLinkAddSlideVel( cdp->data.cell[i], p, n, &vv[i] );
  }
  if( cdp->data.cell[0] == cell )
    zVec3DSub( &vv[0], &vv[1], v );
  else
    zVec3DSub( &vv[1], &vv[0], v );
  return v;
}

zVec6D *_rkFDLinkPointWldVel6D(rkLink *link, zVec3D *p, zVec6D *v)
{
  zVec3D tempv;

  zMulMatVec6D( rkLinkWldAtt(link), rkLinkVel(link), v );
  zVec3DSub( p, rkLinkWldPos(link) ,&tempv );
  zVec3DOuterProd( zVec6DAng(v), &tempv, &tempv );
  zVec3DAddDRC( zVec6DLin(v), &tempv );
  return v;
}

zVec6D *_rkFDChainPointRelativeVel6D(rkCDPair *cdp, zVec3D *p, zVec3D *n, zVec6D *v)
{
  zVec6D vv[2];
  register int i;

  for( i=0; i<2; i++ ){
    if( cdp->data.cell[i]->data.type == RK_CD_CELL_STAT )
      zVec6DClear( &vv[i] );
    else
      _rkFDLinkPointWldVel6D( cdp->data.cell[i]->data.link, p, &vv[i] );
    /* NOTE: this version of the volume-based method does not support slide-mode */
    /* if( cdp->data.cell[i]->data.slide_mode ) */
    /*  _rkFDLinkAddSlideVel( cdp->data.cell[i], p, n, &vv[i] ); */
  }
  zVec6DSub( &vv[0], &vv[1], v );
  return v;
}

/* relative acceleration */
zVec3D *_rkFDLinkPointWldAcc(rkLink *link, zVec3D *p, zVec3D *a)
{
  zVec3D vp;

  zVec3DSub( p, rkLinkWldPos(link), &vp );
  zMulMatTVec3DDRC( rkLinkWldAtt(link), &vp );
  rkLinkPointAcc( link, &vp, a );
  zMulMatVec3DDRC( rkLinkWldAtt(link), a );
  return a;
}

zVec3D *_rkFDChainPointRelativeAcc(rkCDPair *cdp, zVec3D *p, rkCDCell *cell, zVec3D *a)
{
  zVec3D av[2];
  register int i;

  for( i=0; i<2; i++ )
    if( cdp->data.cell[i]->data.type == RK_CD_CELL_STAT )
      zVec3DClear( &av[i] );
    else
      _rkFDLinkPointWldAcc( cdp->data.cell[i]->data.link, p, &av[i] );
  if( cdp->data.cell[0] == cell )
    zVec3DSub( &av[0], &av[1], a );
  else
    zVec3DSub( &av[1], &av[0], a );
  return a;
}

zVec6D *_rkFDLinkPointWldAcc6D(rkLink *link, zVec3D *p, zVec6D *a)
{
  zVec3D vp;

  zVec3DSub( p, rkLinkWldPos(link), &vp );
  zMulMatTVec3DDRC( rkLinkWldAtt(link), &vp );
  rkLinkPointAcc( link, &vp, zVec6DLin(a) );
  zMulMatVec3DDRC( rkLinkWldAtt(link), zVec6DLin(a) );
  zMulMatVec3D( rkLinkWldAtt(link), rkLinkAngAcc(link), zVec6DAng(a) );
  return a;
}

zVec6D *_rkFDChainPointRelativeAcc6D(rkCDPair *cdp, zVec3D *p, zVec6D *a)
{
  zVec6D av[2];
  register int i;

  for( i=0; i<2; i++ )
    if( cdp->data.cell[i]->data.type == RK_CD_CELL_STAT )
      zVec6DClear( &av[i] );
    else{
      _rkFDLinkPointWldAcc6D( cdp->data.cell[i]->data.link, p, &av[i] );
    }
  zVec6DSub( &av[0], &av[1], a );
  return a;
}

/******************************************************************************/
double _rkFDKineticFrictionWeight(rkFD *fd, double fs){
  return 1.0 - exp( -1.0 * fd->_kf_weight * fs );
}

/* joint friction and joint reference */
/* NOTE:
 * before this function is called,
 * the reference value at the one step before must be set
 */
void _rkFDJointFriction(rkFD *fd, bool doUpRef)
{
  rkFDCell *lc;
  rkFDCellDat *ld;
  rkJoint *joint;
  register int i, j;
  double A, b, kf[6], sf, v[6], val, tf;
  rkJointRef jref;

  zListForEach( &fd->list, lc ){
    ld = &lc->data;
    if( rkChainJointSize( &ld->chain ) == 0 ) continue;
    for( i=0; i<rkChainNum(&ld->chain); i++ ){
      joint = rkChainLinkJoint(&ld->chain,i);
      /* 1DoF joint and DC motor only */
      if( rkJointSize(joint) != 1 || rkJointMotorType(joint) != RK_MOTOR_DC ){
        rkJointGetKFric( joint, kf );
        rkJointGetVel( joint, v );
        for( j=0; j<rkJointSize(joint); j++ )
          kf[j] *= _rkFDKineticFrictionWeight( fd, fabs(v[j]) );
        rkJointSetFric( joint, kf );
        continue;
      }
      /* for DC motor model */
      rkJointMotorInertia( joint, &val );
      A = 1.0 / val;
      rkJointMotorInputTrq( joint, &b );
      rkJointMotorRegistance( joint, &val );
      rkJointGetRef( joint, &jref );
      b -= val + jref.ref_trq;
      b *= A * fd->dt;
      rkJointGetVel( joint, v );
      b += v[0];
      rkJointGetDis( joint, &val );
      b += fd->_comp_k * ( val - jref.ref_dis );
      tf = -( A * b / ( zSqr( A ) + fd->_comp_l ) ) / fd->dt;

      rkJointGetSFric( joint, &sf );
      if( fabs( tf ) > sf ){
        rkJointGetKFric( joint, &tf );
        tf *= _rkFDKineticFrictionWeight( fd, fabs(v[0]) );
        if( doUpRef ){
          jref.ref_dis = val;
          jref.type = RK_CONTACT_KF;
          rkJointSetRef( joint, &jref );
        }
      } else{
        if( doUpRef ){
          jref.type = RK_CONTACT_SF;
          rkJointSetRef( joint, &jref );
        }
      }
      rkJointSetFric( joint, &tf );
    }
  }
}

void _rkFDUpdateJointRefDrivingTrq(rkFD *fd)
{
  rkFDCell *lc;
  rkFDCellDat *ld;
  rkJoint *joint;
  register int i, j;
  double val[6], tf[6];
  rkJointRef jref[6];

  zListForEach( &fd->list, lc ){
    ld = &lc->data;
    for( i=0; i<rkChainNum(&ld->chain); i++ ){
      joint = rkChainLinkJoint(&ld->chain,i);
      rkJointGetRef( joint, jref );
      rkJointGetFric( joint, tf );
      rkJointMotorDrivingTrq( joint, val );
      for( j=0; j<rkJointSize(joint); j++ ){
        jref[j].ref_trq = val[j] + tf[j];
      }
      rkJointSetRef( joint, jref );
    }
  }
}

/******************************************************************************/
/* destroy temporary wrench list */
void _rkFDChainExtWrenchDestroy(rkChain *chain)
{
  register int i;

  for( i=0; i<rkChainNum(chain); i++ )
    rkWrenchListDestroy( &rkLinkABIPrp(rkChainLink(chain,i))->wlist );
}

void _rkFDExtWrenchDestroy(rkFD *fd)
{
  rkFDCell *lc;

  zListForEach( &fd->list, lc )
    _rkFDChainExtWrenchDestroy( &lc->data.chain );
}

/******************************************************************************/
/* acceleration update function */
/* this has not be used due to high computation cost */
void _rkFDUpdateAcc(rkFD *fd)
{
  rkFDCell *lc;

  zListForEach( &fd->list, lc ){
    rkChainABI( &lc->data.chain, &lc->data._dis, &lc->data._vel, &lc->data._acc );
  }
}

/* update and set ABIPrp */
void _rkFDUpdateAccBias(rkFD *fd)
{
  rkFDCell *lc;

  zListForEach( &fd->list, lc ){
    if( lc->data.chain._col_flag ){
      rkChainABIUpdate( &lc->data.chain );
      rkChainABIPushPrpAccBias( &lc->data.chain );
      _rkFDChainExtWrenchDestroy( &lc->data.chain );
    }
  }
}

void _rkFDChainUpdateAccAddExForceTwo(rkCDPair *cdp, rkWrench *w[])
{
  if( cdp->data.cell[0]->data.chain == cdp->data.cell[1]->data.chain ) /* self collision */
    rkChainABIUpdateAddExForceTwo( cdp->data.cell[0]->data.chain, cdp->data.cell[0]->data.link, w[0], cdp->data.cell[1]->data.link, w[1] );
  else {
    rkChainABIUpdateAddExForceTwo( cdp->data.cell[0]->data.chain, cdp->data.cell[0]->data.link, w[0], NULL, NULL );
    rkChainABIUpdateAddExForceTwo( cdp->data.cell[1]->data.chain, cdp->data.cell[1]->data.link, w[1], NULL, NULL );
  }
}

void _rkFDChainABIPopPrpExForceTwo(rkCDPair *cdp, rkWrench *w[])
{
  if( cdp->data.cell[0]->data.chain == cdp->data.cell[1]->data.chain )
    rkChainABIPopPrpAccBiasAddExForceTwo( cdp->data.cell[0]->data.chain, cdp->data.cell[0]->data.link, cdp->data.cell[1]->data.link );
  else {
    rkChainABIPopPrpAccBiasAddExForceTwo( cdp->data.cell[0]->data.chain, cdp->data.cell[0]->data.link, NULL );
    rkChainABIPopPrpAccBiasAddExForceTwo( cdp->data.cell[1]->data.chain, cdp->data.cell[1]->data.link, NULL );
  }
}

void _rkFDUpdateAccAddExForce(rkFD *fd)
{
  rkFDCell *lc;

  zListForEach( &fd->list, lc ){
    if( lc->data.chain._col_flag ){
      rkChainABIUpdateAddExForce( &lc->data.chain );
      lc->data.chain._col_flag = false;
    } else
      rkChainABIUpdate( &lc->data.chain );
    rkChainGetJointAccAll( &lc->data.chain, &lc->data._acc );
  }
}

/* update and compute link wrench */
void _rkFDUpdateAccRef(rkFD *fd)
{
  rkFDCell *lc;

  zListForEach( &fd->list, lc ){
    if( lc->data.chain._col_flag ){
      rkChainABIUpdateAddExForceGetWrench( &lc->data.chain );
      lc->data.chain._col_flag = false;
    } else
      rkChainABIUpdateGetWrench( &lc->data.chain );
    rkChainGetJointAccAll( &lc->data.chain, &lc->data._acc );
  }
}

/******************************************************************************/
/* for breaking joint */
bool _rkFDLinkIsBroken(rkLink *link)
{
  switch( rkLinkJointType(link) ){
  case RK_JOINT_REVOL:
    return ( zVec3DInnerProd( zMat3DVec(rkLinkWldAtt(link),zZ), zVec6DAng(rkLinkWrench(link))) > rkLinkBreakPrp(link)->ep_t );
  case RK_JOINT_FLOAT:
    return ( zVec3DNorm( zVec6DLin(rkLinkWrench(link)) ) > rkLinkBreakPrp(link)->ep_f ) ||
           ( zVec3DNorm( zVec6DAng(rkLinkWrench(link)) ) > rkLinkBreakPrp(link)->ep_t );
  default:
    return false;
  }
}

void _rkFDLinkJoinSetFunc(rkLink *link)
{
  switch( rkLinkJointType(link) ){
  case RK_JOINT_REVOL:
    rkJointSetFuncRevol( rkLinkJoint(link) );
    break;
  case RK_JOINT_FLOAT:
    rkJointSetFuncFloat( rkLinkJoint(link) );
    break;
  default:
    break;
  }
}

void _rkFDUpdateLinkBreak(rkFD *fd)
{
  rkFDCell *lc;
  rkBreakPrp *bprp;
  register int i;

  zListForEach( &fd->list, lc )
    for( i=0; i<rkChainNum(&lc->data.chain); i++ ){
      if( (bprp = rkLinkBreakPrp(rkChainLink(&lc->data.chain,i))) &&
          !bprp->is_broken && _rkFDLinkIsBroken( rkChainLink(&lc->data.chain,i) ) ){
        _rkFDLinkJoinSetFunc( rkChainLink(&lc->data.chain,i) );
        bprp->is_broken = true;
      }
    }
}

/******************************************************************************/
/* solver function */
void _rkFDSolveJointContact(rkFD *fd, bool doUpRef)
{
  _rkFDExtWrenchDestroy( fd );
  _rkFDJointFriction( fd, doUpRef );
  rkFDSolveContact( fd, doUpRef );
}

zVec _rkFDUpdate(double t, zVec dis, zVec vel, void *fd, zVec acc)
{
  zVecClear( acc );
  _rkFDConnectJointState( fd, dis, vel, acc );
  _rkFDSolveJointContact( fd, false );
  _rkFDUpdateAccAddExForce( fd );
  return acc;
}

void _rkFDUpdateRef(rkFD *fd)
{
  zVecClear( fd->acc );
  _rkFDConnectJointState( fd, fd->dis, fd->vel, fd->acc );
  _rkFDSolveJointContact( fd, true );
  _rkFDUpdateAccRef( fd );
  _rkFDUpdateJointRefDrivingTrq( fd );
  _rkFDUpdateLinkBreak( fd );
}

rkFD *rkFDSolve(rkFD *fd)
{
  /* reference */
  _rkFDUpdateRef( fd );
  /* integration */
  zODE2Init( &fd->_ode, zVecSize( fd->dis ), fd->_ode_step, _rkFDUpdate );
  zODE2Update( &fd->_ode, fd->t, fd->dis, fd->vel, fd->dt, fd );
  zODE2Destroy( &fd->_ode );
  fd->t += fd->dt;
  /* calc acc */
  _rkFDUpdate( fd->t, fd->dis, fd->vel, fd, fd->acc );
  return fd;
}

void rkFDUpdateInit(rkFD *fd)
{
  zODE2Init( &fd->_ode, zVecSize( fd->dis ), fd->_ode_step, _rkFDUpdate );
  _rkFDUpdateRef( fd );
}

rkFD *rkFDUpdate(rkFD *fd)
{
  zODE2Update( &fd->_ode, fd->t, fd->dis, fd->vel, fd->dt, fd );
  fd->t += fd->dt;
  _rkFDUpdateRef( fd );
  return fd;
}

void rkFDUpdateDestroy(rkFD *fd)
{
  zODE2Destroy( &fd->_ode );
}

/******************************************************************************/
/* other function */
/* update coontact flag */
/* this is called after collision detection */
void _rkFDUpdateRigidColFlag(rkFD *fd)
{
  rkFDCell *lc;
  rkCDPair *cdp;

  zListForEach( &fd->list, lc )
    lc->data.chain._col_flag = false;
  zListForEach( &fd->cd.plist, cdp )
    if( rkFDCDPairIsContactInfoType( cdp, RK_CONTACT_RIGID ) ){
      cdp->data.cell[0]->data.chain->_col_flag = true;
      cdp->data.cell[1]->data.chain->_col_flag = true;
    }
}

/* ************************************************************************** */
/* the forward dynamics based on vertex contact
 * ************************************************************************** */
void rkFDSetContactInfoDefault_Vert(rkFD *fd)
{
  rkContactInfoInit( &fd->_ci_def );
  rkContactInfoSetType( &fd->_ci_def, RK_CONTACT_RIGID );
  rkContactInfoSetK( &fd->_ci_def, 1000.0 );
  rkContactInfoSetL( &fd->_ci_def, 1.0 );
  rkContactInfoSetSF( &fd->_ci_def, 0.5 );
  rkContactInfoSetKF( &fd->_ci_def, 0.3 );
}

void _rkFDContactUpdateRefSlide_Vert(rkCDPair *cdp, rkCDVert *cdv, double dt)
{
  register int i;
  zVec3D tmpv, sv;

  /* for slide mode */
  for( i=0; i<2; i++ ){
    if( cdp->data.cell[i]->data.slide_mode ){
      zVec3DSub( cdv->data.vert, rkLinkWldPos(cdp->data.cell[i]->data.link), &tmpv );
      zMulMatVec3D( rkLinkWldAtt(cdp->data.cell[i]->data.link), &cdp->data.cell[i]->data.slide_axis, &sv );
      zVec3DOuterProd( &sv, &tmpv, &sv );
      zVec3DCatDRC( &sv, -zVec3DInnerProd( &sv, &cdv->data.norm ), &cdv->data.norm );
      if( !zIsTiny( zVec3DNorm( &sv ) ) ){
        zVec3DMulDRC( &sv, (cdp->data.cell[i]==cdv->data.cell?-1.0:1.0) * dt * cdp->data.cell[i]->data.slide_vel/zVec3DNorm( &sv ) );
        zMulMatTVec3DDRC( rkLinkWldAtt(cdp->data.cell[(cdp->data.cell[i]==cdv->data.cell?1:0)]->data.link), &sv );
        zVec3DAddDRC( &cdv->data._ref, &sv );
      }
    }
  }
}

void _rkFDContactModForceVert_Vert(rkFD *fd, rkCDPair *cdp, rkCDVert *cdv, zVec3D v, double dt, bool doUpRef)
{
  double fn, fs, vs;

  fn = zVec3DInnerProd( &cdv->data.f, &cdv->data.axis[0] );
  fs = sqrt( zSqr( zVec3DInnerProd( &cdv->data.f, &cdv->data.axis[1] ) ) + zSqr( zVec3DInnerProd( &cdv->data.f, &cdv->data.axis[2] ) ) );
  if( !zIsTiny( fs ) &&
      fs > ( cdv->data.type == RK_CONTACT_SF ? rkContactInfoSF(cdp->data.ci) : rkContactInfoKF(cdp->data.ci) )*fn ){
    /* kinetic friction */
    zVec3DCatDRC( &v, -zVec3DInnerProd( &v, &cdv->data.axis[0] ), &cdv->data.axis[0] );
    vs = zVec3DNorm( &v );
    zVec3DMul( &cdv->data.axis[0], fn, &cdv->data.f );
    if( !zIsTiny( vs ) ){
      zVec3DNormalizeNCDRC( &v );
      zVec3DCatDRC( &cdv->data.f, - _rkFDKineticFrictionWeight( fd, vs ) * rkContactInfoKF(cdp->data.ci) * fn, &v );
    }
    if( doUpRef ){
      cdv->data.type = RK_CONTACT_KF;
      zVec3DCopy( &cdv->data._pro, &cdv->data._ref );
    }
  } else {
    /* static friction */
    if( doUpRef ){
      cdv->data.type = RK_CONTACT_SF;
      _rkFDContactUpdateRefSlide_Vert( cdp, cdv, dt );
    }
  }
}

void _rkFDContactSetWrench_Vert(rkCDPair *cdp, rkCDVert *cdv)
{
  rkWrench *w;
  register int i;

  for( i=0; i<2; i++){
    w = zAlloc( rkWrench, 1 );
    rkWrenchInit( w );
    zXfer3DInv( rkLinkWldFrame(cdp->data.cell[i]->data.link), cdv->data.vert, rkWrenchPos(w) );
    zMulMatTVec3D( rkLinkWldAtt(cdp->data.cell[i]->data.link), &cdv->data.f, rkWrenchForce(w) );
    if( cdp->data.cell[i] != cdv->data.cell )
      zVec3DRevDRC( rkWrenchForce(w) );
    rkWrenchListPush( &rkLinkABIPrp(cdp->data.cell[i]->data.link)->wlist, w );
  }
}

/**************************************/
/* penalty forces */
void _rkFDContactPenalty_Vert(rkFD *fd, bool doUpRef)
{
  rkCDPair *cdp;
  rkCDVert *cdv;
  zVec3D d, vr;

  fd->colnum_e = 0;
  zListForEach( &fd->cd.plist, cdp ){
    if( !rkFDCDPairIsContactInfoType( cdp, RK_CONTACT_ELASTIC ) ) continue;
    zListForEach( &cdp->data.vlist, cdv ){
      fd->colnum_e++;
      zVec3DSub( cdv->data.vert, &cdv->data.ref, &d );
      _rkFDChainPointRelativeVel( cdp, cdv->data.vert, &cdv->data.norm, cdv->data.cell, &vr );

      /* penalty force */
      zVec3DMul( &d, - rkContactInfoE(cdp->data.ci), &cdv->data.f );
      zVec3DCatDRC( &cdv->data.f, -1.0*( rkContactInfoV(cdp->data.ci) + rkContactInfoE(cdp->data.ci) * rkFDDT(fd) ), &vr );
      if( zVec3DInnerProd( &cdv->data.f, &cdv->data.axis[0] ) < 0.0 ) continue;

      /* modify */
      _rkFDContactModForceVert_Vert( fd, cdp, cdv, vr, rkFDDT(fd), doUpRef );
      /* set wrench */
      _rkFDContactSetWrench_Vert( cdp, cdv );
    }
  }
}

/**************************************/
/* rigid contact */
void _rkFDContactBiasAcc_Vert(rkFD *fd, zVec acc)
{
  rkCDPair *cdp;
  rkCDVert *cdv;
  register int i;
  int offset = 0;
  zVec3D av;

  zListForEach( &fd->cd.plist, cdp ){
    if( !rkFDCDPairIsContactInfoType( cdp, RK_CONTACT_RIGID ) ) continue;
    zListForEach( &cdp->data.vlist, cdv ){
      _rkFDChainPointRelativeAcc( cdp, cdv->data.vert, cdv->data.cell, &av );
      for( i=0; i<3; i++ )
        zVecElem(acc,offset+i) = zVec3DInnerProd( &cdv->data.axis[i], &av );
      offset += 3;
    }
  }
}

void _rkFDContactRelativeAcc_Vert(rkFD *fd, rkCDPair *cp, zVec b, zVec a)
{
  rkCDPair *cdp;
  rkCDVert *cdv;
  int offset = 0;
  zVec3D av;
  register int i;

  zListForEach( &fd->cd.plist, cdp ){
    if( !rkFDCDPairIsContactInfoType( cdp, RK_CONTACT_RIGID ) ) continue;
    if( cdp->data.cell[0]->data.chain != cp->data.cell[0]->data.chain &&
        cdp->data.cell[0]->data.chain != cp->data.cell[1]->data.chain &&
        cdp->data.cell[1]->data.chain != cp->data.cell[0]->data.chain &&
        cdp->data.cell[1]->data.chain != cp->data.cell[1]->data.chain ){
      for( i=0; i<zListNum(&cdp->data.vlist); i++ ){
        zVec3DClear( (zVec3D *)&zVecElem(a,offset) );
        offset += 3;
      }
    } else {
      zListForEach( &cdp->data.vlist, cdv ){
        _rkFDChainPointRelativeAcc( cdp, cdv->data.vert, cdv->data.cell, &av );
        for( i=0; i<3; i++ )
          zVecElem(a,offset+i) = zVec3DInnerProd( &cdv->data.axis[i], &av ) - zVecElem(b,offset+i);
        offset += 3;
      }
    }
  }
}

void _rkFDContactRelationAccForce_Vert(rkFD *fd, zMat a, zVec b)
{
  rkCDPair *cdp;
  rkCDVert *cdv;
  rkWrench *w[2];
  register int i, j;
  int offset = 0;
  zVec t;

  t = zVecAlloc( zVecSize( b ) );
  w[0] = zAlloc( rkWrench, 1 );
  w[1] = zAlloc( rkWrench, 1 );
  /* b */
  _rkFDUpdateAccBias( fd );
  _rkFDContactBiasAcc_Vert( fd, b );
  /* a */
  zListForEach( &fd->cd.plist, cdp ){
    if( !rkFDCDPairIsContactInfoType( cdp, RK_CONTACT_RIGID ) ) continue;
    zListForEach( &cdp->data.vlist, cdv ){
      for( i=0; i<3; i++ ){
        for( j=0; j<2; j++ ){
          rkWrenchInit( w[j] );
          if( cdp->data.cell[j]->data.type == RK_CD_CELL_STAT ) continue;
          zXfer3DInv( rkLinkWldFrame(cdp->data.cell[j]->data.link), cdv->data.vert, rkWrenchPos(w[j]) );
          zMulMatTVec3D( rkLinkWldAtt(cdp->data.cell[j]->data.link), &cdv->data.axis[i], rkWrenchForce(w[j]) );
          if( cdp->data.cell[j] != cdv->data.cell )
            zVec3DRevDRC( rkWrenchForce(w[j]) );
        }
        _rkFDChainUpdateAccAddExForceTwo( cdp, w );
        _rkFDContactRelativeAcc_Vert( fd, cdp, b, t );
        zMatSetCol( a, offset+i, t );

        /* restore ABIPrp */
        _rkFDChainABIPopPrpExForceTwo( cdp, w );
      }
      offset += 3;
    }
  }

  zFree( w[0] );
  zFree( w[1] );
  zVecFree( t );
}

void _rkFDContactBiasVel_Vert(rkFD *fd, zVec b)
{
  rkCDPair *cdp;
  rkCDVert *cdv;
  int offset = 0;
  zVec3D vr;

  zVecMulDRC( b, fd->dt );
  zListForEach( &fd->cd.plist, cdp ){
    if( !rkFDCDPairIsContactInfoType( cdp, RK_CONTACT_RIGID ) ) continue;
    zListForEach( &cdp->data.vlist, cdv ){
      _rkFDChainPointRelativeVel( cdp, cdv->data.vert, &cdv->data.norm, cdv->data.cell, &vr );
      zVecElem(b,offset  ) += zVec3DInnerProd( &vr, &cdv->data.axis[0] );
      zVecElem(b,offset+1) += zVec3DInnerProd( &vr, &cdv->data.axis[1] );
      zVecElem(b,offset+2) += zVec3DInnerProd( &vr, &cdv->data.axis[2] );
      offset += 3;
    }
  }
}

void _rkFDContactCompensateDepth_Vert(rkFD *fd, zVec b, zVec c)
{
  rkCDPair *cdp;
  rkCDVert *cdv;
  int offset = 0;
  zVec3D d;

  zListForEach( &fd->cd.plist, cdp ){
    if( !rkFDCDPairIsContactInfoType( cdp, RK_CONTACT_RIGID ) ) continue;
    zListForEach( &cdp->data.vlist, cdv ){
      zVec3DSub( cdv->data.vert, &cdv->data.ref, &d );
      zVecElem(c,offset  ) = zVecElem(b,offset  ) + rkContactInfoK(cdp->data.ci)                                 * zVec3DInnerProd( &d, &cdv->data.axis[0] );
      zVecElem(c,offset+1) = zVecElem(b,offset+1) + rkContactInfoK(cdp->data.ci) * rkContactInfoKF(cdp->data.ci) * zVec3DInnerProd( &d, &cdv->data.axis[1] );
      zVecElem(c,offset+2) = zVecElem(b,offset+2) + rkContactInfoK(cdp->data.ci) * rkContactInfoKF(cdp->data.ci) * zVec3DInnerProd( &d, &cdv->data.axis[2] );
      offset += 3;
    }
  }
}

void _rkFDContactSetForce_Vert(rkFD *fd, zVec f)
{
  rkCDPair *cdp;
  rkCDVert *cdv;
  int offset = 0;
  register int i;

  zListForEach( &fd->cd.plist, cdp ){
    if( !cdp->data.is_col || rkContactInfoType(cdp->data.ci) != RK_CONTACT_RIGID ||
       ( cdp->data.cell[0]->data.type == RK_CD_CELL_STAT &&
         cdp->data.cell[1]->data.type == RK_CD_CELL_STAT ) ) continue;
    zListForEach( &cdp->data.vlist, cdv ){
      zVec3DClear( &cdv->data.f );
      for( i=0; i<3; i++ ) zVec3DCatDRC( &cdv->data.f, zVecElem(f,offset+i), &cdv->data.axis[i] );
      offset += 3;
    }
  }
}

void _rkFDContactSolveQP_Vert(rkFD *fd, zMat a, zVec b)
{
  zMat q, nf;
  zVec c, d, f;
  rkCDPair *cdp;
  rkCDVert *cdv;
  register int i;
  int offset = 0;

  nf = zMatAlloc( fd->colnum_r, 3*fd->colnum_r );
  d = zVecAlloc( fd->colnum_r );
  q = zMatAllocSqr( 3*fd->colnum_r );
  c = zVecAlloc( 3*fd->colnum_r );
  f = zVecAlloc( 3*fd->colnum_r );
  /* unilateral constraint */
  for( i=0; i<fd->colnum_r; i++ ){
    zMatElem(nf,i,3*i) = 1.0;
  }
  /* QP */
  zMulMatTMat( a, a, q );
  zListForEach( &fd->cd.plist, cdp ){
    if( !rkFDCDPairIsContactInfoType( cdp, RK_CONTACT_RIGID ) ) continue;
    zListForEach( &cdp->data.vlist, cdv ){
      for( i=0; i<3; i++ )
        zMatElem(q,offset+i,offset+i) += rkContactInfoL(cdp->data.ci);
      offset += 3;
    }
  }
  _rkFDContactCompensateDepth_Vert( fd, b, c );
  zMulMatTVecDRC( a, c );

  /* solve */
  zQPSolveASM( q, c, nf, d, f, NULL, NULL, NULL );
  zVecDivDRC( f, fd->dt );
  _rkFDContactSetForce_Vert( fd, f );

  zMatFreeAO( 2, q, nf );
  zVecFreeAO( 3, c, d, f );
}

void _rkFDContactRgidForce_Vert(rkFD *fd)
{
  zMat a;
  zVec b;

  a = zMatAllocSqr( 3*(fd->colnum_r) );
  b = zVecAlloc( 3*(fd->colnum_r) );
  _rkFDContactRelationAccForce_Vert( fd, a, b );
  _rkFDContactBiasVel_Vert( fd, b );
  _rkFDContactSolveQP_Vert( fd, a, b );
  zMatFree( a );
  zVecFree( b );
}

void _rkFDContactModRigidForce_Vert(rkFD *fd, bool doUpRef)
{
  rkCDPair *cdp;
  rkCDVert *cdv;
  zVec3D vr;

  zListForEach( &fd->cd.plist, cdp ){
    if( !rkFDCDPairIsContactInfoType( cdp, RK_CONTACT_RIGID ) ) continue;
    zListForEach( &cdp->data.vlist, cdv ){
      _rkFDChainPointRelativeVel( cdp, cdv->data.vert, &cdv->data.norm, cdv->data.cell, &vr );
      /* modify */
      _rkFDContactModForceVert_Vert( fd, cdp, cdv, vr, rkFDDT(fd), doUpRef );
      /* set wrench */
      _rkFDContactSetWrench_Vert( cdp, cdv );
    }
  }
}

void _rkFDContactRgid_Vert(rkFD *fd, bool doUpRef)
{
  _rkFDUpdateRigidColFlag( fd );
  _rkFDContactRgidForce_Vert( fd );
  _rkFDContactModRigidForce_Vert( fd, doUpRef );
}

void rkFDSolveContact_Vert(rkFD *fd, bool doUpRef)
{
  /* update CD */
  zEchoOff();
  rkCDColChkVert( &fd->cd );
  zEchoOn();
  if( fd->cd.colnum == 0 ) return;
  _rkFDContactPenalty_Vert( fd, doUpRef );
  if( (fd->colnum_r = fd->cd.colnum - fd->colnum_e) == 0 ) return;
  _rkFDContactRgid_Vert( fd, doUpRef );
  return;
}

/* ************************************************************************** */
/* the forward dynamics based on contact volume
 * ************************************************************************** */
void rkFDSetContactInfoDefault_Volume(rkFD *fd)
{
  rkContactInfoInit( &fd->_ci_def );
  rkContactInfoSetType( &fd->_ci_def, RK_CONTACT_RIGID );
  rkContactInfoSetK( &fd->_ci_def, 1000.0 );
  rkContactInfoSetL( &fd->_ci_def, 0.001 );
  rkContactInfoSetSF( &fd->_ci_def, 0.5 );
  rkContactInfoSetKF( &fd->_ci_def, 0.3 );
}

/**************************************/
/* penalty forces */
void _rkFDContactPenalty_Volume(rkFD *fd, bool doUpRef)
{
  fd->colnum_e = 0;
}

/**************************************/
/* rigid contact */
/* ok */
void _rkFDContactBiasAcc_Volume(rkFD *fd, zVec acc)
{
  rkCDPair *cdp;
  int offset = 0;
  zVec6D av;

  zListForEach( &fd->cd.plist, cdp ){
    if( !rkFDCDPairIsContactInfoType( cdp, RK_CONTACT_RIGID ) ) continue;
    _rkFDChainPointRelativeAcc6D( cdp, &cdp->data.center, &av );
    zVec6DCopy( &av, (zVec6D *)&zVecElem(acc,offset) );
    offset += 6;
  }
}

void _rkFDContactRelativeAcc_Volume(rkFD *fd, rkCDPair *cp, zVec b, zVec a)
{
  rkCDPair *cdp;
  int offset = 0;
  zVec6D av;

  zListForEach( &fd->cd.plist, cdp ){
    if( !rkFDCDPairIsContactInfoType( cdp, RK_CONTACT_RIGID ) ) continue;
    if( cdp->data.cell[0]->data.chain != cp->data.cell[0]->data.chain &&
        cdp->data.cell[0]->data.chain != cp->data.cell[1]->data.chain &&
        cdp->data.cell[1]->data.chain != cp->data.cell[0]->data.chain &&
        cdp->data.cell[1]->data.chain != cp->data.cell[1]->data.chain ){
      zVec6DClear( (zVec6D *)&zVecElem(a,offset) );
    } else {
      _rkFDChainPointRelativeAcc6D( cdp, &cdp->data.center, &av );
      zVec6DSub( &av, (zVec6D *)&zVecElem(b,offset), (zVec6D *)&zVecElem(a,offset) );
    }
    offset += 6;
  }
}

void _rkFDContactRelationAccForce_Volume(rkFD *fd, zMat a, zVec b)
{
  rkCDPair *cdp;
  rkWrench *w[2];
  register int i, j;
  int offset = 0;
  zVec t;
  zVec3D pos[2];

  t = zVecAlloc( zVecSize(b) );
  w[0] = zAlloc( rkWrench, 1 );
  w[1] = zAlloc( rkWrench, 1 );
  /* b */
  _rkFDUpdateAccBias( fd );
  _rkFDContactBiasAcc_Volume( fd, b );
  /* A */
  zListForEach( &fd->cd.plist, cdp ){
    if( !rkFDCDPairIsContactInfoType( cdp, RK_CONTACT_RIGID ) ) continue;
    /* let the origin be the center of contact volume */
    for( j=0; j<2; j++ )
      if( cdp->data.cell[j]->data.type != RK_CD_CELL_STAT )
        zXfer3DInv( rkLinkWldFrame(cdp->data.cell[j]->data.link), &cdp->data.center, &pos[j] );
    for( i=0; i<6; i++ ){
      for( j=0; j<2; j++ ){
        rkWrenchInit( w[j] );
        if( cdp->data.cell[j]->data.type == RK_CD_CELL_STAT ) continue;
        rkWrenchSetPos( w[j], &pos[j] );
        if( i<3 ) /* force */
          zMat3DRow( rkLinkWldAtt(cdp->data.cell[j]->data.link), i, rkWrenchForce(w[j]) );
        else /* torque */
          zMat3DRow( rkLinkWldAtt(cdp->data.cell[j]->data.link), i-3, rkWrenchTorque(w[j]) );
        if( j != 0 )
          zVec6DRevDRC( rkWrenchW(w[j]) );
      }
      _rkFDChainUpdateAccAddExForceTwo( cdp, w );
      _rkFDContactRelativeAcc_Volume( fd, cdp, b, t );
      zMatSetCol( a, offset+i, t );
      /* restore ABIPrp */
      _rkFDChainABIPopPrpExForceTwo( cdp, w );
    }
    offset += 6;
  }
  zFree( w[0] );
  zFree( w[1] );
  zVecFree( t );
}

void _rkFDContactBiasVel_Volume(rkFD *fd, zVec b)
{
  rkCDPair *cdp;
  int offset = 0;
  zVec6D vr;

  zVecMulDRC( b, fd->dt );
  zListForEach( &fd->cd.plist, cdp ){
    if( !rkFDCDPairIsContactInfoType( cdp, RK_CONTACT_RIGID ) ) continue;
    _rkFDChainPointRelativeVel6D( cdp, &cdp->data.center, &cdp->data.norm, &vr );
    zVec6DAddDRC( (zVec6D *)&zVecElem(b,offset), &vr );
    offset += 6;
  }
}

/******************/
void _rkFDContactConstraintMidDepth_Volume(double h[], rkContactInfo *ci, double s, double hm[], double *hc)
{
  double k;

  k = rkContactInfoK(ci) * s / 6.0;
  hm[0] = k * ( h[0] + h[1]        );
  hm[1] = k * (        h[1] + h[2] );
  hm[2] = k * ( h[0]        + h[2] );
  *hc   = k * ( h[0] + h[1] + h[2] ) * 2;
}

void _rkFDContactConstraintMidPoint_Volume(zVec3D p[], zVec3D pm[])
{
  zVec3DMid( &p[0], &p[1], &pm[0] );
  zVec3DMid( &p[1], &p[2], &pm[1] );
  zVec3DMid( &p[2], &p[0], &pm[2] );
}

void _rkFDContactConstraintAvePoint_Volume(zVec3D p[], double s, zVec3D *pc)
{
  double k;

  k = s / 3.0;
  zVec3DElem(pc,0) = k * ( zVec3DElem(&p[0],0) + zVec3DElem(&p[1],0) + zVec3DElem(&p[2],0) );
  zVec3DElem(pc,1) = k * ( zVec3DElem(&p[0],1) + zVec3DElem(&p[1],1) + zVec3DElem(&p[2],1) );
  zVec3DElem(pc,2) = k * ( zVec3DElem(&p[0],2) + zVec3DElem(&p[1],2) + zVec3DElem(&p[2],2) );
}

void _rkFDContactConstraintAveMat_Volume(zMat3D m[], double s, zMat3D *mc)
{
  double k;

  k = s / 3.0;
  zMat3DElem(mc,0,0) = k * ( zMat3DElem(&m[0],0,0) + zMat3DElem(&m[1],0,0) + zMat3DElem(&m[2],0,0) );
  zMat3DElem(mc,0,1) = k * ( zMat3DElem(&m[0],0,1) + zMat3DElem(&m[1],0,1) + zMat3DElem(&m[2],0,1) );
  zMat3DElem(mc,0,2) = k * ( zMat3DElem(&m[0],0,2) + zMat3DElem(&m[1],0,2) + zMat3DElem(&m[2],0,2) );
  zMat3DElem(mc,1,0) = k * ( zMat3DElem(&m[0],1,0) + zMat3DElem(&m[1],1,0) + zMat3DElem(&m[2],1,0) );
  zMat3DElem(mc,1,1) = k * ( zMat3DElem(&m[0],1,1) + zMat3DElem(&m[1],1,1) + zMat3DElem(&m[2],1,1) );
  zMat3DElem(mc,1,2) = k * ( zMat3DElem(&m[0],1,2) + zMat3DElem(&m[1],1,2) + zMat3DElem(&m[2],1,2) );
  zMat3DElem(mc,2,0) = k * ( zMat3DElem(&m[0],2,0) + zMat3DElem(&m[1],2,0) + zMat3DElem(&m[2],2,0) );
  zMat3DElem(mc,2,1) = k * ( zMat3DElem(&m[0],2,1) + zMat3DElem(&m[1],2,1) + zMat3DElem(&m[2],2,1) );
  zMat3DElem(mc,2,2) = k * ( zMat3DElem(&m[0],2,2) + zMat3DElem(&m[1],2,2) + zMat3DElem(&m[2],2,2) );
}

double _rkFDContactConstrainArea_Volume(zVec3D p[])
{
  zVec3D e1, e2;

  zVec3DSub( &p[1], &p[0], &e1 );
  zVec3DSub( &p[2], &p[0], &e2 );
  return 0.5 * zVec3DOuterProdNorm( &e1, &e2 );
}

void _rkFDContactConstraintAddQ_Volume(zVec3D p[], zVec3D pm[], double s, zMat6D *q)
{
  register int i;
  zVec3D pc;
  zMat3D mm[3], tmpm;

  _rkFDContactConstraintAvePoint_Volume( p, s, &pc );
  zMat3DMul( ZMAT3DIDENT, s, &tmpm );
  zMat3DAddDRC( zMat6DMat3D(q,0,0), &tmpm );
  zVec3DOuterProdMat3D( &pc, &tmpm );
  zMat3DAddDRC( zMat6DMat3D(q,1,0), &tmpm );
  zMat3DSubDRC( zMat6DMat3D(q,0,1), &tmpm );
  for( i=0; i<3; i++ ) zVec3DOuterProd2Mat3D( &pm[i], &pm[i], &mm[i] );
  _rkFDContactConstraintAveMat_Volume( mm, s, &tmpm );
  zMat3DSubDRC( zMat6DMat3D(q,1,1), &tmpm );
}

void _rkFDContactConstraintDepth_Volume(zVec3D pm[], double h[], rkContactInfo *ci, double s, zVec3D norm, zVec6D *c)
{
  register int i;
  double hm[3], hc;
  zVec3D tmpv;

  _rkFDContactConstraintMidDepth_Volume( h, ci, s, hm, &hc );
  zVec3DMul( &norm, -hc, zVec6DLin(c) );
  zVec3DClear( zVec6DAng(c) );
  for( i=0; i<3; i++ ){
    zVec3DMulDRC( &pm[i], hm[i] );
    zVec3DOuterProd( &norm, &pm[i], &tmpv );
    zVec3DAddDRC( zVec6DAng(c), &tmpv );
  }
}

void _rkFDContactConstraintSignDepth_Volume(double h[], int *st, int stp[])
{
  register int i;

  *st = 0;
  stp[0] = stp[1] = stp[2] = 0;
  for( i=0; i<3; i++ ){
    if( h[i] > zTOL ){
      *st += 1<<(i*2);
      stp[1] = i;
    } else if( h[i] < -zTOL ){
      *st += 1<<(i*2+1);
      stp[2] = i;
    } else{
      stp[0] = i;
    }
  }
}

void _rkFDContactConstraintInnerPoint_Volume(zVec3D *p1, zVec3D *p2, double h1, double h2, zVec3D *pp)
{
  /* for safety */
  if( zIsTiny( h1 ) ){
    zVec3DCopy( p1, pp );
    return;
  }
  if( zIsTiny( h2 ) ){
    zVec3DCopy( p2, pp );
    return;
  }
  if( zIsTiny( h2 - h1 ) ){
    zVec3DMid( p1, p2, pp );
    return;
  }
  zVec3DMul( p1, h2/(h2 - h1), pp );
  zVec3DCatDRC( pp, h1/(h1 - h2), p2 );
}

void _rkFDContactSetContactPlane_Volume(rkCDPair *cdp, zVec3D *p, zVec3D *norm)
{
  rkCDPlane *cdpl, *cdpl2;
  zVec3D tmpv;

  zVec3DCat( norm, -zVec3DInnerProd( norm, &cdp->data.norm ), &cdp->data.norm, &tmpv );
  if( zVec3DIsTiny( &tmpv ) ) return;
  cdpl = zAlloc( rkCDPlane, 1 );
  zVec3DCopy( p, &cdpl->data.v );
  zVec3DDiv( &tmpv, -zVec3DNorm( &tmpv ) , &cdpl->data.norm );

  /* to eliminate identical conditions */
  zListForEach( &cdp->data.cplane, cdpl2 ){
    if( zVec3DIsTol( zVec3DSub( &cdpl->data.norm, &cdpl2->data.norm, &tmpv ), 1e-8 ) &&
        zIsTol( zVec3DInnerProd( &cdpl->data.norm, zVec3DSub( &cdpl2->data.v, &cdpl->data.v, &tmpv ) ), 1e-8 ) ){
      zVec3DOuterProd( &cdpl->data.norm, &tmpv, &tmpv );
      if( zVec3DInnerProd( &cdp->data.norm, &tmpv ) > 0.0 ){
        zVec3DCopy( &cdpl->data.v, &cdpl2->data.v );
      }
      zFree( cdpl );
      return;
    }
  }
  zListInsertHead( &cdp->data.cplane, cdpl );
}

int __rk_fd_plane_cmp(void *p1, void *p2, void *priv)
{
  zVec3D *a, *n, tmp;
  double th1, th2;

  n = (zVec3D *)priv;
  a = (zVec3D *)priv + 1;
  zVec3DOuterProd( a, &((rkCDPlane *)p1)->data.norm, &tmp );
  if( zVec3DInnerProd( &tmp, n ) > 0 )
    th1 = atan2( -zVec3DNorm( &tmp ), zVec3DInnerProd( a, &((rkCDPlane *)p1)->data.norm ) );
  else
    th1 = atan2(  zVec3DNorm( &tmp ), zVec3DInnerProd( a, &((rkCDPlane *)p1)->data.norm ) );
  zVec3DOuterProd( a, &((rkCDPlane *)p2)->data.norm, &tmp );
  if( zVec3DInnerProd( &tmp, n ) > 0 )
    th2 = atan2( -zVec3DNorm( &tmp ), zVec3DInnerProd( a, &((rkCDPlane *)p2)->data.norm ) );
  else
    th2 = atan2(  zVec3DNorm( &tmp ), zVec3DInnerProd( a, &((rkCDPlane *)p2)->data.norm ) );
  if( zIsTiny( th1 - th2 ) ) return 0;
  return th1 > th2 ? 1: -1;
}

void _rkFDContactConstraint_Volume(rkFD *fd, rkCDPair *cdp, rkContactInfo *ci, zMat6D *q, zVec6D *c)
{
  register int i, j;
  zTri3D *face;
  double h[3], s;
  zVec3D pf[3], p[3], pm[3], pp[2];
  zVec6D cc;
  int stp[3], st;

  zMat6DClear( q );
  zVec6DClear( c );
  for( i=0; i<zPH3DFaceNum(&cdp->data.colvol); i++ ){
    face = zPH3DFace(&cdp->data.colvol,i);
    for( j=0; j<3; j++ ){
      zVec3DSub( zTri3DVert(face,j), &cdp->data.center, &pf[j] );
      h[j] = zVec3DInnerProd( &cdp->data.norm, &pf[j] );
      zVec3DCat( &pf[j], -h[j], &cdp->data.norm, &p[j] );
    }
    _rkFDContactConstraintMidPoint_Volume( p, pm );
    s = _rkFDContactConstrainArea_Volume( p );
    _rkFDContactConstraintAddQ_Volume( p, pm, s, q ); /* q */
    _rkFDContactConstraintDepth_Volume( pm, h, ci, s, cdp->data.norm, &cc ); /* c */

    /* a process dependent on signum of h */
    _rkFDContactConstraintSignDepth_Volume( h, &st, stp );
    switch( st ){
    case 0x01: case 0x04: case 0x10: /* two points are on the plane and one point is above it */
    case 0x05: case 0x11: case 0x14: /* one point is on the plane and two points are above it */
      _rkFDContactSetContactPlane_Volume( cdp, &pf[stp[0]], zTri3DNorm(face) );
    case 0x15: /* all points are above the plane */
      zVec6DAddDRC( c, &cc );
      continue;
    case 0x02: case 0x08: case 0x20: /* two points are on the line and one point is below it */
    case 0x0a: case 0x22: case 0x28: /* one point is on the line and two points are below it */
      _rkFDContactSetContactPlane_Volume( cdp, &pf[stp[0]], zTri3DNorm(face) );
    case 0x2a: /* all points are below the plane */
      zVec6DSubDRC( c, &cc );
      continue;
    case 0x24: /* 0:on    1:above 2:below */
    case 0x12: /* 0:below 1:on    2:above */
    case 0x09: /* 0:above 1:below 2:on    */
      zVec6DAddDRC( c, &cc );
      _rkFDContactConstraintInnerPoint_Volume( &pf[stp[1]], &pf[stp[2]], h[stp[1]], h[stp[2]], &p[stp[1]] );
      h[stp[1]] = 0.0;
      zVec3DCopy( &p[stp[0]], &pp[0] );
      zVec3DCopy( &p[stp[1]], &pp[1] );
      break;
    case 0x06: /* 0:below 1:above 2:on    */
    case 0x21: /* 0:above 1:on    2:below */
    case 0x18: /* 0:on    1:below 2:above */
      zVec6DAddDRC( c, &cc );
      _rkFDContactConstraintInnerPoint_Volume( &pf[stp[1]], &pf[stp[2]], h[stp[1]], h[stp[2]], &p[stp[2]] );
      h[stp[2]] = 0.0;
      zVec3DCopy( &p[stp[2]], &pp[0] );
      zVec3DCopy( &p[stp[0]], &pp[1] );
      break;
    case 0x16:case 0x19:case 0x25: /* two points are above the plane and one point is below it */
      stp[0] = (stp[2]+1) % 3;
      stp[1] = (stp[0]+1) % 3;
      zVec6DAddDRC( c, &cc );
      _rkFDContactConstraintInnerPoint_Volume( &pf[stp[2]], &pf[stp[0]], h[stp[2]], h[stp[0]], &p[stp[0]] );
      _rkFDContactConstraintInnerPoint_Volume( &pf[stp[2]], &pf[stp[1]], h[stp[2]], h[stp[1]], &p[stp[1]] );
      h[stp[0]] = h[stp[1]] = 0.0;
      zVec3DCopy( &p[stp[0]], &pp[0] );
      zVec3DCopy( &p[stp[1]], &pp[1] );
      break;
    case 0x1a:case 0x26:case 0x29: /* two points are below the plane and one point is above it */
      stp[0] = (stp[1]+1) % 3;
      stp[2] = (stp[0]+1) % 3;
      zVec6DSubDRC( c, &cc );
      _rkFDContactConstraintInnerPoint_Volume( &pf[stp[1]], &pf[stp[0]], h[stp[1]], h[stp[0]], &p[stp[0]] );
      _rkFDContactConstraintInnerPoint_Volume( &pf[stp[1]], &pf[stp[2]], h[stp[1]], h[stp[2]], &p[stp[2]] );
      h[stp[0]] = h[stp[2]] = 0.0;
      zVec3DCopy( &p[stp[2]], &pp[0] );
      zVec3DCopy( &p[stp[0]], &pp[1] );
      break;
    default: /* all points are on the plane */
      continue;
    }
    s = _rkFDContactConstrainArea_Volume( p );
    _rkFDContactConstraintMidPoint_Volume( p, pm );
    _rkFDContactConstraintDepth_Volume( pm, h, ci, s, cdp->data.norm, &cc );
    switch( st ){
    case 0x1a: case 0x26: case 0x29:
      zVec6DCatDRC( c, 2.0, &cc );
      break;
    default:
      zVec6DCatDRC( c, -2.0, &cc );
    }
    /* register pp */
    _rkFDContactSetContactPlane_Volume( cdp, &pp[0], zTri3DNorm(face) );
  }
  /* sort plane list */
  rkCDPlaneListQuickSort( &cdp->data.cplane, __rk_fd_plane_cmp, cdp->data.axis );
}

/******************/
void _rkFDContactSetForce_Volume(rkFD *fd, zVec f)
{
  rkCDPair *cdp;
  int offset =  0;

  zListForEach( &fd->cd.plist, cdp ){
    if( !rkFDCDPairIsContactInfoType( cdp, RK_CONTACT_RIGID ) ) continue;
    if( zListNum(&cdp->data.cplane) == 0 ){
      zVec6DClear( zVec6DLin(&cdp->data.f) );
      continue;
    }
    zVec6DCopy( (zVec6D *)&zVecElem(f,offset), &cdp->data.f );
    if( zVec3DIsTiny( zVec6DLin(&cdp->data.f) ) || zVec3DInnerProd( zVec6DLin(&cdp->data.f), &cdp->data.norm ) < zTOL ){
      zVec6DClear( zVec6DLin(&cdp->data.f) );
    }
    offset += 6;
  }
}

zVec _rkFDContactSolveConstraintInit_Volume(zMat a, zVec b, zVec ans, void *util){
  rkCDPair *cdp;
  register int offset = 0;

  zVecClear( ans );
  zListForEach( &((rkFD *)util)->cd.plist, cdp ){
    if( !rkFDCDPairIsContactInfoType( cdp, RK_CONTACT_RIGID ) ) continue;
    /* normal forces initial offset */
    zVec3DMul( &cdp->data.norm, 1.0, (zVec3D *)&zVecElem(ans,offset) );
    offset += 6;
  }
  return ans;
}

void _rkFDContactSolveConstraint_Volume(rkFD *fd, zMat a, zVec b)
{
  zMat q, nf;
  zVec c, df, f;
  zMat6D qv;
  zVec6D cv, tmpv;
  rkCDPair *cdp;
  rkCDPlane *cdpl;
  register int i, j;
  int offset = 0, pvertn = 0;

  /* evaluation function */
  q  = zMatAllocSqr( 6*fd->colnum_r );
  c  = zVecAlloc( 6*fd->colnum_r );
  zListForEach( &fd->cd.plist, cdp ){
    if( !rkFDCDPairIsContactInfoType( cdp, RK_CONTACT_RIGID ) ) continue;
    _rkFDContactConstraint_Volume( fd, cdp, cdp->data.ci, &qv, &cv );
    /* q */
    for( i=0; i<6; i++ )
      for( j=0; j<6; j++ )
        zRawMatCatDyad( zMatBuf(q), zMat6DElem(&qv,i,j), zMatRowBuf(a,offset+i), zMatColSizeNC(a), zMatRowBuf(a,offset+j), zMatColSizeNC(a) );
    /* c */
    zMulMat6DVec6D( &qv, (zVec6D *)&zVecElem(b,offset), &tmpv );
    zVec6DAddDRC( &cv, &tmpv );
    for( i=0; i<6; i++ )
      zRawVecCatDRC( zVecBuf(c), zVec6DElem(&cv,i), zMatRowBuf(a,offset+i), zMatColSizeNC(a) );
    offset += 6;
    pvertn += zListNum(&cdp->data.cplane);
  }

  /* constraint */
  nf = zMatAlloc( fd->colnum_r+pvertn, 6*fd->colnum_r );
  df = zVecAlloc( fd->colnum_r+pvertn );
  f  = zVecAlloc( 6*fd->colnum_r );
  i = j = 0;
  zListForEach( &fd->cd.plist, cdp ){
    if( !rkFDCDPairIsContactInfoType( cdp, RK_CONTACT_RIGID ) ) continue;
    zVec3DCopy( &cdp->data.norm, (zVec3D *)&zMatElem(nf,i++,j) );
    zListForEach( &cdp->data.cplane, cdpl ){
      zVec3DMul( &cdp->data.norm, -zVec3DInnerProd( &cdpl->data.norm, &cdpl->data.v ), (zVec3D *)&zMatElem(nf,i,j) );
      zVec3DMul( &cdp->data.axis[1],  zVec3DInnerProd( &cdpl->data.norm, &cdp->data.axis[2] ), (zVec3D *)&zMatElem(nf,i,j+3) );
      zVec3DCatDRC( (zVec3D *)&zMatElem(nf,i,j+3), -zVec3DInnerProd( &cdpl->data.norm, &cdp->data.axis[1] ), &cdp->data.axis[2] );
      i++;
    }
    j += 6;
  }
  /* relaxation */
  offset = 0;
  zListForEach( &fd->cd.plist, cdp ){
    if( !rkFDCDPairIsContactInfoType( cdp, RK_CONTACT_RIGID ) ) continue;
    for( i=0; i<6; i++ )
      zMatElem(q,offset+i,offset+i) += rkContactInfoL(cdp->data.ci);
    offset += 6;
  }

  /* solve qp */
  zQPSolveASM( q, c, nf, df, f, NULL, _rkFDContactSolveConstraintInit_Volume, fd );
  zVecDivDRC( f, fd->dt );
  _rkFDContactSetForce_Volume( fd, f );

  zMatFreeAO( 2, q, nf );
  zVecFreeAO( 3, c, df, f );
}

void _rkFDContactRigidForce_Volume(rkFD *fd)
{
  zMat a;
  zVec b;

  a = zMatAllocSqr( 6*(fd->colnum_r) );
  b = zVecAlloc( 6*(fd->colnum_r) );
  _rkFDContactRelationAccForce_Volume( fd, a, b );
  _rkFDContactBiasVel_Volume( fd, b );
  _rkFDContactSolveConstraint_Volume( fd, a, b );
  zMatFree( a );
  zVecFree( b );
}

/******************/
/* for safety */
/* NOTE:
 * the ASM solver does not work well occasionally.
 * in that case, the solution does not meet the constraints with respect to normal force.
 * before the friction modification, the solution is projected inside the constraints with a small margin.
 */
void _rkFDContactModNormForceCenterTrq_Volume(rkCDPair *cdp, zVec3D *r, double fn)
{
  zVec3DMul( &cdp->data.norm, zVec3DInnerProd( &cdp->data.norm, zVec6DAng(&cdp->data.f) ), zVec6DAng(&cdp->data.f) );
  zVec3DCatDRC( zVec6DAng(&cdp->data.f),  fn * zVec3DInnerProd( &cdp->data.axis[2], r ), &cdp->data.axis[1] );
  zVec3DCatDRC( zVec6DAng(&cdp->data.f), -fn * zVec3DInnerProd( &cdp->data.axis[1], r ), &cdp->data.axis[2] );
}

void _rkFDContactModNormForceCenter_Volume(rkFD *fd)
{
  rkCDPair *cdp;
  rkCDPlane *cdpl[4];
  zVec3D r0, r, dir, tmp;
  double fn, d, s;
  bool flag = false;

  zListForEach( &fd->cd.plist, cdp ){
    if( !rkFDCDPairIsContactInfoType( cdp, RK_CONTACT_RIGID ) ) continue;
    /* unilateral constraint */
    fn = zVec3DInnerProd( &cdp->data.norm, zVec6DLin(&cdp->data.f) );
    if( fn < zTOL ) continue;

    zVec3DMul( &cdp->data.axis[1], -zVec3DInnerProd( &cdp->data.axis[2], zVec6DAng(&cdp->data.f) ) / fn, &r0 );
    zVec3DCatDRC( &r0, zVec3DInnerProd( &cdp->data.axis[1], zVec6DAng(&cdp->data.f) ) / fn, &cdp->data.axis[2] );
		flag = false;
    cdpl[2] = zListHead( &cdp->data.cplane );
    cdpl[1] = zListCellPrev( cdpl[2] );
    cdpl[0] = zListCellPrev( cdpl[1] );
    for( cdpl[3]=zListTail(&cdp->data.cplane); cdpl[3]!=zListRoot(&cdp->data.cplane);
         cdpl[0]=cdpl[1],cdpl[1]=cdpl[2],cdpl[2]=cdpl[3],cdpl[3]=zListCellNext(cdpl[3]) ){
      zVec3DSub( &cdpl[2]->data.v, &cdpl[1]->data.v, &dir );
      d = zVec3DSqrNorm( &dir );
      if( zIsTiny( d ) ) continue;
      zVec3DSub( &r0, &cdpl[1]->data.v, &tmp );
      if( zVec3DInnerProd( &tmp, &cdpl[1]->data.norm ) > zTOL ) continue;
      s = zVec3DInnerProd( &dir, &tmp ) / d;
      if( s < zTOL ){
        if( flag ) break;
        /* r: cdpl[1]->data.v */
        zVec3DSub( &cdpl[0]->data.v, &cdpl[1]->data.v, &tmp );
        zVec3DAddDRC( &tmp, &dir );
        zVec3DCat( &cdpl[1]->data.v, zTOL / zVec3DNorm( &tmp ), &tmp, &r );
        _rkFDContactModNormForceCenterTrq_Volume( cdp, &r, fn );
        break;
      } else if ( s < 1.0-zTOL ){
        /* r: projected position */
        zVec3DCat( &cdpl[1]->data.v, s, &dir, &r );
        zVec3DCatDRC( &r, zTOL, &cdpl[1]->data.norm );
        _rkFDContactModNormForceCenterTrq_Volume( cdp, &r, fn );
        break;
      } else {
        /* r: cdpl[2]->data.v */
        zVec3DSub( &cdpl[3]->data.v, &cdpl[2]->data.v, &tmp );
        zVec3DSubDRC( &tmp, &dir );
        zVec3DCat( &cdpl[2]->data.v, zTOL / zVec3DNorm( &tmp ), &tmp, &r );
        _rkFDContactModNormForceCenterTrq_Volume( cdp, &r, fn );
        flag = true;
      }
    }
  }
}

/******************/
void _rkFDContactPlaneVertPos_Volume(rkCDPair *cdp, double *tl)
{
  rkCDPlane *cdpl;
  double rl;

  *tl = 0.0;
  zListForEach( &cdp->data.cplane, cdpl ){
    cdpl->data.r[0] = zVec3DInnerProd( &cdpl->data.v, &cdp->data.axis[1] );
    cdpl->data.r[1] = zVec3DInnerProd( &cdpl->data.v, &cdp->data.axis[2] );
    if( *tl < (rl = zVec2DNorm( cdpl->data.r )) )
      *tl = rl;
  }
}

void _rkFDContactModWrenchStaticConstraint_Volume(rkCDPair *cdp, int n, zVec6D *w, zMat a, zVec b)
{
  rkCDPlane *cdpl;
  double th;
  register int i;
  int offset = 0;

  zListForEach( &cdp->data.cplane, cdpl ){
    for( i=0,th=0.0; i<n; i++ ){
      zMatElem(a,0,offset+i) = 1.0;
      zMatElem(a,1,offset+i) = cdpl->data.r[1];
      zMatElem(a,2,offset+i) = -cdpl->data.r[0];
      zMatElem(a,3,offset+i) = rkContactInfoSF(cdp->data.ci)*cos(th);
      zMatElem(a,4,offset+i) = rkContactInfoSF(cdp->data.ci)*sin(th);
      zMatElem(a,5,offset+i) = -(zMatElem(a,2,offset+i)*zMatElem(a,4,offset+i) + zMatElem(a,1,offset+i)*zMatElem(a,3,offset+i));
      th += 2.0 * zPI / (double)n;
    }
    offset += n;
  }
  zVecElem(b,0) = zVec6DElem(w,0);
  zVecElem(b,1) = zVec6DElem(w,4);
  zVecElem(b,2) = zVec6DElem(w,5);
  zVecElem(b,3) = zVec6DElem(w,1);
  zVecElem(b,4) = zVec6DElem(w,2);
  zVecElem(b,5) = zVec6DElem(w,3);
}

bool _rkFDContactModWrenchStatic_Volume(rkFD *fd, rkCDPair *cdp, zVec6D *w, zMat a, zVec b, zVec f)
{
  bool ret;

  /* constraint */
  _rkFDContactModWrenchStaticConstraint_Volume( cdp, fd->_pyramid, w, a, b );
  zEchoOff();
  ret = zLPFeasibleBase( a, b, f );
  zEchoOn();
  return ret;
}

void _rkFDContactModWrenchKineticConstraint_Volume(rkCDPair *cdp, zVec6D *w, zMat a, zVec b)
{
  rkCDPlane *cdpl;
  int offset = 0;

  zListForEach( &cdp->data.cplane, cdpl ){
    zMatElem(a,0,offset) = 1.0;
    zMatElem(a,1,offset) = cdpl->data.r[1];
    zMatElem(a,2,offset) = -cdpl->data.r[0];
    offset++;
  }
  zVecElem(b,0) = zVec6DElem(w,0);
  zVecElem(b,1) = zVec6DElem(w,4);
  zVecElem(b,2) = zVec6DElem(w,5);
}

void _rkFDContactPlaneVertSlideDir_Volume(rkFD *fd, rkCDPair *cdp, rkCDPlane *cdpl)
{
  zVec3D p, v;
  double nv, w;

  zVec3DAdd( &cdp->data.center, &cdpl->data.v, &p );
  _rkFDChainPointRelativeVel( cdp, &p, &cdp->data.norm, cdp->data.cell[0], &v );
  zVec3DCatDRC( &v, -zVec3DInnerProd( &cdp->data.norm, &v ), &cdp->data.norm );
  nv = zVec3DNorm( &v );
  if( zIsTiny( nv ) ){
    cdpl->data.s[0] = 0.0;
    cdpl->data.s[1] = 0.0;
  } else {
    w = _rkFDKineticFrictionWeight( fd, nv ) * rkContactInfoKF(cdp->data.ci) / nv;
    cdpl->data.s[0] = -w * zVec3DInnerProd( &v, &cdp->data.axis[1] );
    cdpl->data.s[1] = -w * zVec3DInnerProd( &v, &cdp->data.axis[2] );
  }
}

void _rkFDContactModWrenchKineticTotal_Volume(rkCDPair *cdp, zVec f, zVec6D *w)
{
  rkCDPlane *cdpl;
  int offset = 0;
  zVec2D fs;

  zVec3DClear( (zVec3D *)&zVec6DElem(w,1) );
  zListForEach( &cdp->data.cplane, cdpl ){
    zVec2DMul( cdpl->data.s, zVecElem(f,offset), fs );
    zVec6DElem(w,1) += fs[0];
    zVec6DElem(w,2) += fs[1];
    zVec6DElem(w,3) += cdpl->data.r[0] * fs[1] - cdpl->data.r[1] * fs[0];
    offset++;
  }
}

void _rkFDContactModWrenchKinetic_Volume(rkFD *fd, rkCDPair *cdp, zVec6D *w, zMat a, zVec b, zVec f)
{
  rkCDPlane *cdpl;
  zVec c;
  int offset = 0;
  double wn[3];
  register int i;

  /* constraint */
  _rkFDContactModWrenchKineticConstraint_Volume( cdp, w, a, b );
  /* evaluation function */
  c = zVecAlloc( zListNum(&cdp->data.cplane) );
  for( i=0; i<3; i++ )
    if( !zIsTiny( zVec6DElem(w,i+1) ) )
      wn[i] = 1.0 / zVec6DElem(w,i+1);
    else
      wn[i] = 0.0;
  zListForEach( &cdp->data.cplane, cdpl ){
    _rkFDContactPlaneVertSlideDir_Volume( fd, cdp, cdpl );
    zVecElem(c,offset) = -wn[0] * cdpl->data.s[0] - wn[1] * cdpl->data.s[1] -
      wn[2] * ( cdpl->data.r[0] * cdpl->data.s[1] - cdpl->data.r[1] * cdpl->data.s[0] );
    offset++;
  }

  zEchoOff();
  if( !zLPSolveSimplex( a, b, c, f, NULL ) ){
    /* for safety */
    zMatSetSize( a, 1, zListNum(&cdp->data.cplane) );
    zVecSetSize( b, 1 );
    for( i=0; i<2; i++ )
      if( !zIsTiny( zVec6DElem(w,i+3) ) )
        wn[i] = 1.0 / zVec6DElem(w,i+3);
      else
        wn[i] = 0.0;
    offset = 0;
    zListForEach( &cdp->data.cplane, cdpl ){
      _rkFDContactPlaneVertSlideDir_Volume( fd, cdp, cdpl );
      zVecElem(c,offset) += wn[0] * cdpl->data.r[0] - wn[1] * cdpl->data.r[1];
      offset++;
    }
    zLPSolveSimplex( a, b, c, f, NULL );
  }
  /* zLPSolvePDIP_PC( a, b, c, f, NULL ); */
  zEchoOn();

  _rkFDContactModWrenchKineticTotal_Volume( cdp, f, w );
  zVecFree( c );
}

void _rkFDContactModWrenchKineticCenter_Volume(rkFD *fd, rkCDPair *cdp, zVec6D *w)
{
  zVec3D v;
  double nv, tmpd;

  _rkFDChainPointRelativeVel( cdp, &cdp->data.center, &cdp->data.norm, cdp->data.cell[0], &v );
  zVec3DCatDRC( &v, -zVec3DInnerProd( &cdp->data.norm, &v ), &cdp->data.norm );
  nv = zVec3DNorm( &v );
  if( zIsTiny( nv ) ){
    zVec6DElem(w,1) = 0.0;
    zVec6DElem(w,2) = 0.0;
  } else {
    tmpd = _rkFDKineticFrictionWeight( fd, nv ) * rkContactInfoKF(cdp->data.ci) * zVec6DElem(w,0) / nv;
    zVec6DElem(w,1) = -tmpd * zVec3DInnerProd( &v, &cdp->data.axis[1] );
    zVec6DElem(w,2) = -tmpd * zVec3DInnerProd( &v, &cdp->data.axis[2] );
  }
}

/* if this function's overhead is large, this function should be expanded */
void _rkFDContactModWrenchSetKinetic_Volume(rkCDPair *cdp, bool doUpRef)
{
  if( doUpRef ){
    cdp->data.type = RK_CONTACT_KF;
    /* this version doesn't use the reference frame */
    /* zFrame3DCopy( rkLinkWldFrame(cdp->data.cell[0]->data.link), &cdp->data.ref[0] ); */
    /* zFrame3DCopy( rkLinkWldFrame(cdp->data.cell[1]->data.link), &cdp->data.ref[1] ); */
  }
}

void _rkFDContactModWrenchSetStatic_Volume(rkCDPair *cdp, bool doUpRef)
{
  if( doUpRef ){
    cdp->data.type = RK_CONTACT_SF;
  }
}

void _rkFDContactModWrenchSetForce_Volume(rkCDPair *cdp, zVec6D *w)
{
  register int i;

  zVec6DClear( &cdp->data.f );
  for( i=0; i<3; i++ ){
    zVec3DCatDRC( zVec6DLin(&cdp->data.f), zVec6DElem(w,i  ), &cdp->data.axis[i] );
    zVec3DCatDRC( zVec6DAng(&cdp->data.f), zVec6DElem(w,i+3), &cdp->data.axis[i] );
  }
}

void _rkFDContactModWrench_Volume(rkFD *fd, bool doUpRef)
{
  rkCDPair *cdp;
  zVec6D w;
  double fn, fs, tl;
  register int i;
  zMat a;
  zVec b, f;

  zListForEach( &fd->cd.plist, cdp ){
    if( !rkFDCDPairIsContactInfoType( cdp, RK_CONTACT_RIGID ) ) continue;
    if( zListIsEmpty(&cdp->data.cplane) ) continue;
    if( zIsTiny( zVec3DInnerProd( zVec6DLin(&cdp->data.f), &cdp->data.axis[0] ) ) ) continue;
    for( i=0; i<3; i++ ){
      zVec6DElem(&w,i  ) = zVec3DInnerProd( zVec6DLin(&cdp->data.f), &cdp->data.axis[i] );
      zVec6DElem(&w,i+3) = zVec3DInnerProd( zVec6DAng(&cdp->data.f), &cdp->data.axis[i] );
    }
    /* friction */
    fn = zVec6DElem(&w,0);
    fs = sqrt( zSqr( zVec6DElem(&w,1) ) + zSqr( zVec6DElem(&w,2) ) );
    _rkFDContactPlaneVertPos_Volume( cdp, &tl );
    if( zIsTiny( tl ) ){
      /* the plane is too small */
      zVec3DClear( zVec6DAng(&w) );
      if( !zIsTiny( fs ) && fs > rkContactInfoSF(cdp->data.ci) * fn ){
        _rkFDContactModWrenchKineticCenter_Volume( fd, cdp, &w );
        _rkFDContactModWrenchSetKinetic_Volume( cdp, doUpRef );
      } else
        _rkFDContactModWrenchSetStatic_Volume( cdp, doUpRef );
      _rkFDContactModWrenchSetForce_Volume( cdp, &w );
      continue;
    }
    /* rectangular constraint: for reducing allocated memory */
    if( (!zIsTiny( fs ) && fs > rkContactInfoSF(cdp->data.ci) * fn) ||
        fabs(zVec6DElem(&w,3)) > tl * zVec6DElem(&w,0) ){
      /* workspace */
      a = zMatAlloc( 3, zListNum(&cdp->data.cplane) );
      b = zVecAlloc( 3 );
      f = zVecAlloc( zListNum(&cdp->data.cplane) );
      _rkFDContactModWrenchKinetic_Volume( fd, cdp, &w, a, b, f );
      _rkFDContactModWrenchSetForce_Volume( cdp, &w );
      _rkFDContactModWrenchSetKinetic_Volume( cdp, doUpRef );
    } else {
      /* workspace */
      a = zMatAlloc( 6, fd->_pyramid * zListNum(&cdp->data.cplane) );
      b = zVecAlloc( 6 );
      f = zVecAlloc( fd->_pyramid * zListNum(&cdp->data.cplane) );
      if( !_rkFDContactModWrenchStatic_Volume( fd, cdp, &w, a, b, f ) ){
        /* reuse momory of mat a, vec b, f */
        zMatSetSize( a, 3, zListNum(&cdp->data.cplane) );
        zVecSetSize( b, 3 );
        zVecSetSize( f, zListNum(&cdp->data.cplane) );
        _rkFDContactModWrenchKinetic_Volume( fd, cdp, &w, a, b, f );
        _rkFDContactModWrenchSetKinetic_Volume( cdp, doUpRef );
        _rkFDContactModWrenchSetForce_Volume( cdp, &w );
      } else
        _rkFDContactModWrenchSetStatic_Volume( cdp, doUpRef );
    }

    zMatFree( a );
    zVecFreeAO( 2, b, f );
  }
}

/******************/
void _rkFDContactSetWrench_Volume(rkFD *fd)
{
  rkCDPair *cdp;
  rkWrench *w;
  register int i;

  zListForEach( &fd->cd.plist, cdp ){
    if( !cdp->data.is_col || rkContactInfoType(cdp->data.ci) != RK_CONTACT_RIGID ) continue;
    for( i=0; i<2; i++ ){
      w = zAlloc( rkWrench, 1 );
      rkWrenchInit( w );
      zXfer3DInv( rkLinkWldFrame(cdp->data.cell[i]->data.link), &cdp->data.center, rkWrenchPos(w) );
      zMulMatTVec6D( rkLinkWldAtt(cdp->data.cell[i]->data.link), &cdp->data.f, rkWrenchW(w) );
      if( i == 1 )
        zVec6DRevDRC( rkWrenchW(w) );
      rkWrenchListPush( &rkLinkABIPrp(cdp->data.cell[i]->data.link)->wlist, w );
    }
  }
}

void _rkFDContactRgid_Volume(rkFD *fd, bool doUpRef)
{
  _rkFDUpdateRigidColFlag( fd );
  _rkFDContactRigidForce_Volume( fd );
  _rkFDContactModNormForceCenter_Volume( fd );
  _rkFDContactModWrench_Volume( fd, doUpRef );
  _rkFDContactSetWrench_Volume( fd );
}

void rkFDSolveContact_Volume(rkFD *fd, bool doUpRef)
{
  /* update CD */
  zEchoOff();
  rkCDColVolBREP( &fd->cd );
  zEchoOn();
  if( fd->cd.colnum == 0 ) return;
  _rkFDContactPenalty_Volume( fd, doUpRef );
  if( (fd->colnum_r = fd->cd.colnum - fd->colnum_e) == 0 ) return;
  _rkFDContactRgid_Volume( fd, doUpRef );
}

/******************************************************************************/
/* for debug */
void rkFDWrite(rkFD *fd)
{
  rkFDCell *lc;

  zListForEach( &fd->list, lc )
    rkChainWrite( &lc->data.chain );
}
