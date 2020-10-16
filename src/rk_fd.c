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
    if(!(pdis = zVecClone( fd->dis )) || !(pvel = zVecClone( fd->vel )) )
      goto ERROR;

  zVecFreeAO( 3, fd->dis, fd->vel, fd->acc );

  rjs = rkChainJointSize( &rlc->data.chain );

  if(!( fd->dis = zVecAlloc( fd->size + rjs ) ) ||
     !( fd->vel = zVecAlloc( fd->size + rjs ) ) ||
     !( fd->acc = zVecAlloc( fd->size + rjs ) ) )
    goto ERROR;

  zListForEach( &fd->list, lc ){
    if( lc ==  rlc ) break;
    js = rkChainJointSize( &lc->data.chain );
    if( pdis ) zRawVecCopy( &zVecElem( pdis, offset ), &zVecElem( fd->dis, offset ), js );
    if( pvel ) zRawVecCopy( &zVecElem( pvel, offset ), &zVecElem( fd->vel, offset ), js );
    _rkFDCellDatSetOffset( fd, &lc->data, offset );
    offset += js;
  }
  zRawVecClear( &zVecElem( fd->dis, offset ), rjs );
  zRawVecClear( &zVecElem( fd->vel, offset ), rjs );
  zVecSetSize( &lc->data._dis, rjs );
  zVecSetSize( &lc->data._vel, rjs );
  _rkFDCellDatSetOffset( fd, &lc->data, offset );
  offset += rjs;
  for( lc=zListCellNext(lc); lc!=zListRoot(&fd->list); lc=zListCellNext(lc) ){
    js = rkChainJointSize( &lc->data.chain );
    if( pdis ) zRawVecCopy( &zVecElem( pdis, offset ), &zVecElem( fd->dis, offset + rjs ), js );
    if( pvel ) zRawVecCopy( &zVecElem( pvel, offset ), &zVecElem( fd->vel, offset + rjs ), js );
    _rkFDCellDatSetOffset( fd, &lc->data, offset + rjs );
    offset += js;
  }

  zVecFreeAO( 2, pdis, pvel );
  return true;
 ERROR:
  zVecFreeAO( 4, pdis, pvel, fd->dis, fd->vel );
  return false;
}

rkFDCellDat *_rkFDCellDatInit(rkFDCellDat *ld)
{
  rkChainABIInit( &ld->chain );
  ld->_dis.size = ld->_vel.size = ld->_acc.size = rkChainJointSize( &ld->chain );
  ld->_dis.elem = ld->_vel.elem = ld->_acc.elem = NULL;
  return ld;
}

rkFDCell *_rkFDCellPush(rkFD *fd, rkFDCell *lc)
{
  _rkFDCellDatInit( &lc->data );
  zListInsertHead( &fd->list, lc );
  fd->size += rkChainJointSize( &lc->data.chain );
  if( !_rkFDAllocJointState( fd, lc ) ){
    rkFDDestroy( fd );
    return NULL;
  }
  /* cd reg */
  rkCDChainReg( &fd->cd, &lc->data.chain, true );
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
    if( lc == cell ){
      if( !_rkFDAllocJointState( fd, lc ) ){
        rkFDDestroy( fd );
        return false;
      }
      zListPurge( &fd->list, lc );
      _rkFDCellDatFree( &lc->data );
      zFree( lc );
      return true;
    }
  }
  return false;
}

/******************************************************************************/
void rkFDChainSetDis(rkFDCell *lc, zVec dis)
{
  zVecCopy( dis, &lc->data._dis );
  rkChainSetJointDisAll( &lc->data.chain, dis );
}

/******************************************************************************/
void rkFDChainFK(rkFD *fd, zVec dis)
{
  rkFDCell *lc;
  zVecStruct ldis;

  zListForEach( &fd->list, lc ){
    ldis.size = rkChainJointSize( &lc->data.chain );
    ldis.elem = &zVecElem( dis, lc->data._offset );
    rkChainFK( &lc->data.chain, &ldis );
  }
}

void rkFDChainUpdateRate(rkFD *fd, zVec vel, zVec acc)
{
  rkFDCell *lc;
  zVecStruct lvel, lacc;

  zListForEach( &fd->list, lc ){
    lvel.size = lacc.size = rkChainJointSize( &lc->data.chain );
    lvel.elem = &zVecElem( vel, lc->data._offset );
    lacc.elem = &zVecElem( acc, lc->data._offset );
    rkChainSetJointRateAll( &lc->data.chain, &lvel, &lacc );
    rkChainUpdateRateGrav( &lc->data.chain );
  }
}

void _rkFDChainUpdateFKRateCellDat(rkFDCellDat *ld)
{
  rkChainFK( &ld->chain, &ld->_dis );
  rkChainSetJointRateAll( &ld->chain, &ld->_vel, &ld->_acc );
  rkChainUpdateRateGrav( &ld->chain );
}

void rkFDChainUpdateFKRate(rkFD *fd)
{
  rkFDCell *lc;

  zListForEach( &fd->list, lc ){
    _rkFDChainUpdateFKRateCellDat( &lc->data );
  }
}

/******************************************************************************/
bool rkFDContactInfoReadFile(rkFD *fd, char filename[]){
  if( zArrayNum( &fd->ci ) != 0)
    rkContactInfoPoolDestroy( &fd->ci );
  return rkContactInfoPoolReadFile( &fd->ci, filename );
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
zVec rkFDODECatDefault(zVec x, double k, zVec v, zVec xnew, void *util)
{
  rkFDCell *lc;
  zVecStruct lv, lxn;

  zVecCopyNC( x, xnew );
  zListForEach( &((rkFD *)util)->list, lc ){
    lv.size  = rkChainJointSize( &lc->data.chain );
    lxn.size = rkChainJointSize( &lc->data.chain );
    lv.elem  = &zVecElem( v,    lc->data._offset );
    lxn.elem = &zVecElem( xnew, lc->data._offset );
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
    lx2.elem = &zVecElem( x2, lc->data._offset );
    ldx.elem = &zVecElem( dx, lc->data._offset );
    rkChainSubJointDisAll( &lc->data.chain, &ldx, &lx2 );
  }
  return dx;
}

/******************************************************************************/
void _rkFDChainConnectJointState(rkFD *fd, zVec dis, zVec vel, zVec acc)
{
  rkFDCell *lc;

  zListForEach( &fd->list, lc ){
    lc->data._dis.elem = &zVecElem( dis, lc->data._offset );
    lc->data._vel.elem = &zVecElem( vel, lc->data._offset );
    lc->data._acc.elem = &zVecElem( acc, lc->data._offset );
    _rkFDChainUpdateFKRateCellDat( &lc->data );
  }
}

/******************************************************************************/
void _rkFDJointCalcFriction(rkFD *fd, bool doUpRef)
{
  rkFDCell *lc;
  rkFDCellDat *ld;
  rkJoint *joint;
  register int i;
  double A, b, kf[6], sf, val, tf;
  rkJointRef jref;

  zListForEach( &fd->list, lc ){
    ld = &lc->data;
    if( rkChainJointSize( &ld->chain ) == 0 ) continue;
    for( i=0; i<rkChainNum( &ld->chain ); i++ ){
      joint = rkChainLinkJoint( &ld->chain, i );
      /* 1DoF joint and DC motor only */
      if( rkJointSize(joint) != 1 || rkJointMotorType(joint) != RK_MOTOR_DC ){
        rkJointGetKFric( joint, kf );
        rkJointSetFric( joint, kf );
        continue;
      }
      rkJointMotorInertia( joint, &val );
      A = 1.0 / val;
      rkJointMotorInputTrq( joint, &b );
      rkJointMotorRegistance( joint, &val );
      rkJointGetRef( joint, &jref );
      b -= val + jref.ref_trq;
      b *= A * fd->dt;
      rkJointGetVel( joint, &val );
      b += val;
      rkJointGetDis( joint, &val );
      b += fd->_comp_k * ( val - jref.ref_dis );
      tf = -( A * b / ( zSqr( A ) + fd->_comp_l ) ) / fd->dt;

      rkJointGetSFric( joint, &sf );

      if( fabs( tf ) > sf ){
        rkJointGetKFric( joint, kf );
        tf = kf[0];
        if( doUpRef ){
          jref.ref_dis = val;
          jref.type = RK_CONTACT_SLIP;
        }
      } else{
        if( doUpRef ){
          jref.type = RK_CONTACT_STICK;
        }
      }
      rkJointSetFric( joint, &tf );
      rkJointSetRef( joint, &jref );
    }
  }
}

/******************************************************************************/
void _rkFDUpdateAcc(rkFD *fd)
{
  rkFDCell *lc;

  zListForEach( &fd->list, lc )
    rkChainABI( &lc->data.chain, &lc->data._dis, &lc->data._vel, &lc->data._acc );
}

void _rkFDContactRelativeAcc(rkFD *fd, zVec acc)
{
  rkCDPair *cdp;
  rkCDVert *cdv;
  register int i;
  int offset = 0;
  zVec3D av[2], vp, lvw, tempv;
  zVec6D la;
  rkCDCellDat *celld[2];

  _rkFDUpdateAcc( fd );
  zListForEach( &fd->cd.plist, cdp ){
    if( !cdp->data.is_col ) continue;
    zListForEach( &cdp->data.vlist, cdv ){
      celld[0] = &cdv->data.cell->data;
      celld[1] = cdv->data.cell != cdp->data.cell[0] ? &cdp->data.cell[0]->data : &cdp->data.cell[1]->data;
      /* compute acc */
      for( i=0; i<2; i++){
        if( celld[i]->type == RK_CD_CELL_STAT ){
          zVec3DClear( &av[i] );
          continue;
        }
        zMulMatVec3D( rkLinkWldAtt( celld[i]->link ), zVec6DAng( rkLinkVel( celld[i]->link ) ), &lvw );
        zMulMatVec6D( rkLinkWldAtt( celld[i]->link ), rkLinkAcc( celld[i]->link ), &la );
        zVec3DSub( cdv->data.vert, rkLinkWldPos( celld[i]->link ), &vp );
        zVec3DTripleProd( &lvw, &lvw, &vp, &av[i] );
        zVec3DOuterProd( zVec6DAng( &la ), &vp, &tempv );
        zVec3DAddDRC( &av[i], &tempv );
        zVec3DAddDRC( &av[i], zVec6DLin( &la ) );
      }
      zVec3DSub( &av[0], &av[1], &tempv );

      for( i=0; i<3; i++ ){
        zVecElem( acc, offset+i ) = zVec3DInnerProd( &cdv->data.axis[i], &tempv );
      }
      offset += 3;
    }
  }
}

void _rkFDChainExtWrenchDestroy(rkFD *fd)
{
  rkFDCell *lc;

  zListForEach( &fd->list, lc )
    rkChainExtWrenchDestroy( &lc->data.chain );
}

void _rkFDContactRelationAccForce(rkFD *fd, zMat A, zVec b)
{
  rkCDPair *cdp;
  rkCDVert *cdv;
  rkWrench *w;
  register int i, j;
  int offset = 0;
  rkCDCellDat *celld[2];
  zVec t;

  t = zVecAlloc( zVecSize( b ) );
  /* b */
  _rkFDChainExtWrenchDestroy( fd );
  _rkFDContactRelativeAcc( fd, b );
  /* A */
  zListForEach( &fd->cd.plist, cdp ){
    if( !cdp->data.is_col ) continue;
    zListForEach( &cdp->data.vlist, cdv ){
      celld[0] = &cdv->data.cell->data;
      celld[1] = cdv->data.cell != cdp->data.cell[0] ? &cdp->data.cell[0]->data : &cdp->data.cell[1]->data;
      for( i=0; i<3; i++ ){
        for( j=0; j<2; j++ ){
          if( celld[j]->type == RK_CD_CELL_STAT ) continue;
          w = zAlloc( rkWrench, 1 );
          rkWrenchInit( w );
          zXfer3DInv( rkLinkWldFrame( celld[j]->link ), cdv->data.vert, rkWrenchPos( w ) );
          zMulMatTVec3D( rkLinkWldAtt( celld[j]->link ), &cdv->data.axis[i], rkWrenchForce( w ) );
          zVec3DMulDRC( rkWrenchForce( w ), j==0 ? 1.0 : -1.0 );
          rkLinkExtWrenchPush( celld[j]->link, w );
        }
        /* calc acc */
        _rkFDContactRelativeAcc( fd, t );
        zVecSubDRC( t, b );
        zMatSetCol( A, offset+i, t );
        /* destroy */
        _rkFDChainExtWrenchDestroy( fd );
      }
      offset += 3;
    }
  }
  zVecFree( t );
}

static rkContactInfo _rkFDContactInfoDefault={
  { NULL, NULL },
  RK_CONTACT_RIGID,
  { { 1000.0, 1.0 } },
  0.5, 0.3,
};

void _rkFDContactCompensateDepth(rkFD *fd, zMat A, zVec b)
{
  rkCDPair *cdp;
  rkCDVert *cdv;
  rkContactInfo *ci;
  register int i;
  int offset = 0;
  rkCDCellDat *celld[2];
  zVec3D d, vr, vv[2], tempv, sv;
  zVec6D v6;

  zVecMulDRC( b, fd->dt );
  _rkFDUpdateAcc( fd );
  zListForEach( &fd->cd.plist, cdp ){
    if( !cdp->data.is_col ) continue;
    ci = rkContactInfoPoolAssoc( &fd->ci, rkLinkStuff(cdp->data.cell[0]->data.link), rkLinkStuff(cdp->data.cell[1]->data.link) );
    if( !ci ) ci = &_rkFDContactInfoDefault;
    zListForEach( &cdp->data.vlist, cdv ){
      celld[0] = &cdv->data.cell->data;
      celld[1] = cdv->data.cell != cdp->data.cell[0] ? &cdp->data.cell[0]->data : &cdp->data.cell[1]->data;
      /* d */
      zVec3DSub( cdv->data.vert, &cdv->data.ref, &d );
      /* v */
      for( i=0; i<2; i++){
        if( celld[i]->type == RK_CD_CELL_STAT ) zVec3DClear( &vv[i] );
        else{
          zMulMatVec6D( rkLinkWldAtt( celld[i]->link ), rkLinkVel( celld[i]->link ), &v6 );
          zVec3DSub( cdv->data.vert, rkLinkWldPos( celld[i]->link ) ,&tempv );
          zVec3DOuterProd( zVec6DAng( &v6 ), &tempv, &tempv );
          zVec3DAdd( zVec6DLin( &v6 ), &tempv, &vv[i] );
        }
        /* for a fake-crawler */
        if( celld[i]->slide_mode ){
          zVec3DSub( cdv->data.vert, rkLinkWldPos( celld[i]->link ), &tempv );
          zMulMatVec3D( rkLinkWldAtt( celld[i]->link ), &celld[i]->slide_axis, &sv );
          zVec3DOuterProd( &sv, &tempv, &sv );
          zVec3DCatDRC( &sv, -zVec3DInnerProd( &sv, &cdv->data.norm ), &cdv->data.norm );
          if( !zVec3DIsTiny( &sv ) )
            zVec3DCatDRC( &vv[i], celld[i]->slide_vel/zVec3DNorm(&sv), &sv );
        }
      }
      zVec3DSub( &vv[0], &vv[1], &vr );

      /* compensation */
      zVecElem( b, offset   ) += zVec3DInnerProd( &vr, &cdv->data.axis[0] ) + rkContactInfoK( ci )                         * zVec3DInnerProd( &d, &cdv->data.axis[0] );
      zVecElem( b, offset+1 ) += zVec3DInnerProd( &vr, &cdv->data.axis[1] ) + rkContactInfoK( ci ) * rkContactInfoKF( ci ) * zVec3DInnerProd( &d, &cdv->data.axis[1] );
      zVecElem( b, offset+2 ) += zVec3DInnerProd( &vr, &cdv->data.axis[2] ) + rkContactInfoK( ci ) * rkContactInfoKF( ci ) * zVec3DInnerProd( &d, &cdv->data.axis[2] );
      offset += 3;
    }
  }
}

void _rkFDContactSetForce(rkFD *fd, zVec f)
{
  rkCDPair *cdp;
  rkCDVert *cdv;
  int offset = 0;

  zListForEach( &fd->cd.plist, cdp ){
    if( !cdp->data.is_col ||
       ( cdp->data.cell[0]->data.type == RK_CD_CELL_STAT &&
         cdp->data.cell[1]->data.type == RK_CD_CELL_STAT ) ) continue;
    zListForEach( &cdp->data.vlist, cdv ){
      zVec3DCopy( (zVec3D *)&zVecElem( f, offset ), &cdv->data.f );
      offset += 3;
    }
  }
}

void _rkFDContactSolveQP(rkFD *fd, zMat A, zVec b)
{
  zMat Q, C;
  zVec c, d, f;
  register int i;
  rkCDPair *cdp;
  rkCDVert *cdv;
  rkContactInfo *ci;
  int offset = 0;

  /* unilateral constraint */
  C = zMatAlloc( fd->cd.colnum, 3*fd->cd.colnum );
  d = zVecAlloc( fd->cd.colnum );
  zMatClear( C );
  zVecClear( d );
  for( i=0; i<fd->cd.colnum; i++ ){
    zMatElem( C, i, 3*i ) = 1.0;
  }

  /* QP */
  Q = zMatAllocSqr( 3*fd->cd.colnum );
  c = zVecAlloc( 3*fd->cd.colnum );
  zMulMatTMat( A, A, Q );

  zListForEach( &fd->cd.plist, cdp ){
    if( !cdp->data.is_col ) continue;
    ci = rkContactInfoPoolAssoc( &fd->ci, rkLinkStuff(cdp->data.cell[0]->data.link), rkLinkStuff(cdp->data.cell[1]->data.link) );
    if( !ci ) ci = &_rkFDContactInfoDefault;
    zListForEach( &cdp->data.vlist, cdv ){
      for( i=0; i<3; i++ )
        zMatElem( Q, offset+i, offset+i ) += rkContactInfoL( ci );
      offset += 3;
    }
  }
  zMulMatTVec( A, b, c );

  /* solve */
  f = zVecAlloc( 3*fd->cd.colnum );
  zQPSolveASM( Q, c, C, d, f, NULL, NULL, NULL );

  zVecDivDRC( f, fd->dt );
  _rkFDContactSetForce( fd, f );

  zMatFreeAO( 2, Q, C );
  zVecFreeAO( 3, c, d, f );
}

void _rkFDContactCalcForce(rkFD *fd)
{
  zMat A;
  zVec b;

  A = zMatAllocSqr( 3*(fd->cd.colnum) );
  b = zVecAlloc( 3*(fd->cd.colnum) );
  _rkFDContactRelationAccForce( fd, A, b );
  _rkFDContactCompensateDepth( fd, A, b );
  _rkFDContactSolveQP( fd, A, b );
  zMatFree( A );
  zVecFree( b );
}

void _rkFDContactModForce(rkFD *fd, bool doUpRef)
{
  rkCDPair *cdp;
  rkCDVert *cdv;
  rkContactInfo *ci;
  double fn, fs;
  rkCDCellDat *celld[2];
  zVec3D f, vr, vv[2], tempv, sv;
  zVec6D v6;
  rkWrench *w;
  register int i;

  _rkFDChainExtWrenchDestroy( fd );
  zListForEach( &fd->cd.plist, cdp ){
    if( !cdp->data.is_col ) continue;
    ci = rkContactInfoPoolAssoc( &fd->ci, rkLinkStuff(cdp->data.cell[0]->data.link), rkLinkStuff(cdp->data.cell[1]->data.link) );
    if( !ci ) ci = &_rkFDContactInfoDefault;
    zListForEach( &cdp->data.vlist, cdv ){
      celld[0] = &cdv->data.cell->data;
      celld[1] = cdv->data.cell != cdp->data.cell[0] ? &cdp->data.cell[0]->data : &cdp->data.cell[1]->data;
      /* v */
      for( i=0; i<2; i++){
        if( celld[i]->type == RK_CD_CELL_STAT ){
          zVec3DClear( &vv[i] );
        } else{
          zMulMatVec6D( rkLinkWldAtt( celld[i]->link ), rkLinkVel( celld[i]->link ), &v6 );
          zVec3DSub( cdv->data.vert, rkLinkWldPos( celld[i]->link ) ,&tempv );
          zVec3DOuterProd( zVec6DAng( &v6 ), &tempv, &tempv );
          zVec3DAdd( zVec6DLin( &v6 ), &tempv, &vv[i] );
        }
      }
      zVec3DSub( &vv[0], &vv[1], &vr );
      zVec3DCatDRC( &vr, -zVec3DInnerProd( &vr, &cdv->data.axis[0] ), &cdv->data.axis[0] );
      if( zVec3DIsTiny( &vr ) ) zVec3DClear( &vr );
      else zVec3DNormalizeDRC( &vr );

      fn = cdv->data.f.e[zX];
      fs = sqrt( zSqr(cdv->data.f.e[zY]) + zSqr(cdv->data.f.e[zZ]) );

      /* friction corn */
      if( !zIsTiny( fs ) &&
          fs > ( cdv->data.type == RK_CONTACT_STICK ? rkContactInfoSF(ci) : rkContactInfoKF(ci) )*fn ){
        zVec3DMul( &vr, -rkContactInfoKF(ci)*fn, &f );
        zVec3DCatDRC( &f, fn, &cdv->data.axis[0] );
        cdv->data.f.e[zY] = zVec3DInnerProd( &f, &cdv->data.axis[1] );
        cdv->data.f.e[zZ] = zVec3DInnerProd( &f, &cdv->data.axis[2] );
        if( doUpRef ){
          cdv->data.type = RK_CONTACT_SLIP;
          zVec3DCopy( &cdv->data._pro, &cdv->data._ref );
        }
      } else{
        zVec3DClear( &f );
        for( i=zX; i<=zZ; i++ )
          zVec3DCatDRC( &f, cdv->data.f.e[i], &cdv->data.axis[i] );
        if( doUpRef ){
          cdv->data.type = RK_CONTACT_STICK;
          /* for a fake-crawler */
          for( i=0; i<2; i++ ){
            if( celld[i]->slide_mode ){
              zVec3DSub( cdv->data.vert, rkLinkWldPos( celld[i]->link ), &tempv );
              zMulMatVec3D( rkLinkWldAtt( celld[i]->link ), &celld[i]->slide_axis, &sv );
              zVec3DOuterProd( &sv, &tempv, &sv );
              zVec3DCatDRC( &sv, -zVec3DInnerProd( &sv, &cdv->data.norm ), &cdv->data.norm );
              if( !zVec3DIsTiny( &sv ) ){
                zVec3DMulDRC( &sv, (i==0?-1.0:1.0)*rkFDDT(fd)*celld[i]->slide_vel/zVec3DNorm(&sv) );
                zMulMatTVec3DDRC( rkLinkWldAtt( celld[1]->link ), &sv );
                zVec3DAddDRC( &cdv->data._ref, &sv );
              }
            }
          }
        }
      }

      /* set wrench */
      for( i=0; i<2; i++){
        w = zAlloc( rkWrench, 1 );
        rkWrenchInit( w );
        zXfer3DInv( rkLinkWldFrame(celld[i]->link), cdv->data.vert, rkWrenchPos(w) );
        zMulMatTVec3D( rkLinkWldAtt(celld[i]->link), &f, rkWrenchForce(w) );
        rkLinkExtWrenchPush( celld[i]->link, w );
        zVec3DRevDRC( &f );
      }
    }
  }
}

/******************************************************************************/
void _rkFDSolveJointContact(rkFD *fd, bool doUpRef)
{
  /* update CD */
  zEchoOff();
  rkCDColChkVert( &fd->cd );
  zEchoOn();
  /* compute joint friction */
  _rkFDJointCalcFriction( fd, doUpRef );
  if( fd->cd.colnum == 0 ) return;
  /* compute contact force */
  _rkFDContactCalcForce( fd );
  /* modify contact force */
  _rkFDContactModForce( fd, doUpRef );
  return;
}

/******************************************************************************/
zVec _rkFDUpdate(double t, zVec dis, zVec vel, void *fd, zVec acc)
{
  zVecClear( acc );
  _rkFDChainConnectJointState( fd, dis, vel, acc );
  _rkFDChainExtWrenchDestroy( fd );
  _rkFDSolveJointContact( fd, false );
  zVecClear( acc );
  _rkFDUpdateAcc( fd );
  return acc;
}

void _rkFDUpdateRefDrivingTorque(rkFD *fd)
{
  rkFDCell *lc;
  rkFDCellDat *ld;
  rkJoint *joint;
  register int i, j;
  double val[6], tf[6];
  rkJointRef jref[6];

  zListForEach( &fd->list, lc ){
    ld = &lc->data;
    for( i=0; i<rkChainNum( &ld->chain ); i++ ){
      joint = rkChainLinkJoint( &ld->chain, i );
      rkJointGetRef( joint, jref );
      rkJointGetFric( joint, tf );
      rkJointMotorDrivingTrq( joint, val );
      for( j=0; j<rkJointSize( joint ); j++ ){
        jref[j].ref_trq = val[j] + tf[j];
      }
    }
  }
}

void _rkFDUpdateRef(rkFD *fd)
{
  zVecClear( fd->acc );
  _rkFDChainConnectJointState( fd, fd->dis, fd->vel, fd->acc );
  _rkFDChainExtWrenchDestroy( fd );

  _rkFDSolveJointContact( fd, true );
  _rkFDUpdateAcc( fd );
  zVecClear( fd->acc );
  _rkFDUpdateRefDrivingTorque( fd );
}

/******************************************************************************/
void _rkFDJointRefInit(rkFD *fd)
{
  register int i, j;
  rkJoint *joint;
  double val[6];
  rkJointRef jref[6];
  rkFDCell *lc;

  zListForEach( &fd->list, lc ){
    for( i=0; i<rkChainNum( &lc->data.chain ); i++ ){
      joint = rkChainLinkJoint( &lc->data.chain, i );
      rkJointGetDis( joint, val );
      for( j=0; j<rkJointSize( joint ); j++ ){
        jref[j].ref_dis = val[j];
        jref[j].ref_trq = 0.0;
        jref[j].type = RK_CONTACT_STICK;
      }
      rkJointSetRef( joint, jref );
    }
  }
}

/* rkFDSolve_Vert
 * - solve the forward dynamics based on point contacts.
 */
rkFD *rkFDSolve_Vert(rkFD *fd)
{
  if( zIsTiny( fd->t ) ) _rkFDJointRefInit( fd );
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

void rkFDUpdateInit_Vert(rkFD *fd)
{
  zODE2Init( &fd->_ode, zVecSize( fd->dis ), fd->_ode_step, _rkFDUpdate );
  _rkFDJointRefInit( fd );
  _rkFDUpdateRef( fd );
}

rkFD *rkFDUpdate_Vert(rkFD *fd)
{
  zODE2Update( &fd->_ode, fd->t, fd->dis, fd->vel, fd->dt, fd );
  fd->t += fd->dt;
  _rkFDUpdateRef( fd );
  return fd;
}

void rkFDUpdateDestroy_Vert(rkFD *fd)
{
  zODE2Destroy( &fd->_ode );
}

/******************************************************************************/
/* the forward dynamics based on volumetric contacts
 */
#define RK_FD_SOLVE_ITER 8
void _rkFDContactRelativeAcc_Volume(rkFD *fd, zVec acc)
{
  rkCDPair *cdp;
  register int i;
  int offset = 0;
  zVec3D vp, lvw, tempv;
  zVec6D av[2], la;

  _rkFDUpdateAcc( fd );
  zListForEach( &fd->cd.plist, cdp ){
    if( !cdp->data.is_col ) continue;
    for( i=0; i<2; i++ ){
      if( cdp->data.cell[i]->data.type == RK_CD_CELL_STAT ){
        zVec6DClear( &av[i] );
        continue;
      }
      zMulMatVec3D( rkLinkWldAtt( cdp->data.cell[i]->data.link ), zVec6DAng( rkLinkVel( cdp->data.cell[i]->data.link ) ), &lvw );
      zMulMatVec6D( rkLinkWldAtt( cdp->data.cell[i]->data.link ), rkLinkAcc( cdp->data.cell[i]->data.link ), &la );
      zVec6DCopy( &la, &av[i] );
      zVec3DSub( ZVEC3DZERO, rkLinkWldPos( cdp->data.cell[i]->data.link ), &vp );
      zVec3DTripleProd( &lvw, &lvw, &vp, &tempv );

      zVec3DAddDRC( zVec6DLin(&av[i]), &tempv );
      zVec3DOuterProd( zVec6DAng(&la), &vp, &tempv );
      zVec3DAddDRC( zVec6DLin(&av[i]), &tempv );
    }
    zVec6DSubDRC( &av[0], &av[1] );

    zVec6DCopy( &av[0], (zVec6D *)&zVecElem(acc,offset) );
    offset += 6;
  }
}

void _rkFDContactRelationAccForce_Volume(rkFD *fd, zMat A, zVec b)
{
  rkCDPair *cdp;
  rkWrench *w;
  register int i, j;
  int offset = 0;
  zVec t;
  zVec3D pos[2];

  t = zVecAlloc( zVecSize(b) );
  /* b */
  _rkFDChainExtWrenchDestroy( fd );
  _rkFDContactRelativeAcc_Volume( fd, b );
  /* A */
  zListForEach( &fd->cd.plist, cdp ){
    if( !cdp->data.is_col ) continue;
    /* preprocess */
    for( j=0; j<2; j++ ){
      /* the original point of the inertial frame (with respect to the link
         frame) is regarded as the point of application of contact force. */
      if( cdp->data.cell[j]->data.type == RK_CD_CELL_STAT ) continue;
      zXfer3DInv( rkLinkWldFrame(cdp->data.cell[j]->data.link), ZVEC3DZERO, &pos[j] );
    }
    for( i=0; i<6; i++ ){
      for( j=0; j<2; j++ ){
        if( cdp->data.cell[j]->data.type == RK_CD_CELL_STAT ) continue;
        w = zAlloc( rkWrench, 1 );
        zListCellInit( w );
        rkWrenchSetPos( w, &pos[j] );

        rkWrenchSetW( w, ZVEC6DZERO );
        zVec6DElem( rkWrenchW(w), i ) = 1.0;
        if( i<3 ){ /* force */
          zMulMatTVec3DDRC( rkLinkWldAtt(cdp->data.cell[j]->data.link), rkWrenchForce(w) );
        } else{ /* torque */
          zMulMatTVec3DDRC( rkLinkWldAtt(cdp->data.cell[j]->data.link), rkWrenchTorque(w) );
        }
        zVec6DMulDRC( rkWrenchW(w), j==0 ? 1.0 : -1.0 );
        rkLinkExtWrenchPush( cdp->data.cell[j]->data.link, w );
      }
      /* calc acc */
      _rkFDContactRelativeAcc_Volume( fd, t );
      zVecSubDRC( t, b );
      zMatSetCol( A, offset+i, t );
      /* destroy */
      _rkFDChainExtWrenchDestroy( fd );
    }
    offset += 6;
  }
  zVecFree( t );
}

void _rkFDContactVelAfterInt_Volume(rkFD *fd, zMat A, zVec b)
{
  rkCDPair *cdp;
  register int i;
  int offset = 0;
  zVec3D tempv;
  zVec6D vv[2];

  zVecMulDRC( b, fd->dt );
  _rkFDUpdateAcc( fd );
  zListForEach( &fd->cd.plist, cdp ){
    if( !cdp->data.is_col ) continue;
    for( i=0; i<2; i++){
      if( cdp->data.cell[i]->data.type == RK_CD_CELL_STAT )
        zVec6DClear( &vv[i] );
      else{
        zMulMatVec6D( rkLinkWldAtt(cdp->data.cell[i]->data.link), rkLinkVel(cdp->data.cell[i]->data.link), &vv[i] );
        zVec3DSub( ZVEC3DZERO, rkLinkWldPos(cdp->data.cell[i]->data.link) ,&tempv );
        zVec3DOuterProd( zVec6DAng(&vv[i]), &tempv, &tempv );
        zVec3DAddDRC( zVec6DLin(&vv[i]), &tempv );
      }
    }
    zVec6DSubDRC( &vv[0], &vv[1] );
    zVec6DAddDRC( (zVec6D *)&zVecElem(b,offset), &vv[0] );
    offset += 6;
  }
}

void _rkFDContactConstraintGetMidDepth_Volume(double h[], rkContactInfo *ci, double s, double hm[], double *hc)
{
  double k;

  k = rkContactInfoK(ci) * s / 6.0;
  hm[0] = k * ( h[0] + h[1]        );
  hm[1] = k * (        h[1] + h[2] );
  hm[2] = k * ( h[0]        + h[2] );
  *hc   = k * ( h[0] + h[1] + h[2] ) * 2;
}

void _rkFDContactConstraintGetMidPoint_Volume(zVec3D p[], zVec3D pm[])
{
  zVec3DMid( &p[0], &p[1], &pm[0] );
  zVec3DMid( &p[1], &p[2], &pm[1] );
  zVec3DMid( &p[2], &p[0], &pm[2] );
}

void _rkFDContactConstraintDepth_Volume(zVec3D pm[], double hm[], double hc, zVec3D norm, zVec6D *c)
{
  register int i;
  zVec3D tempv;

  zVec3DMul( &norm, -hc, zVec6DLin(c) );
  zVec3DClear( zVec6DAng(c) );
  for( i=0; i<3; i++ ){
    zVec3DMulDRC( &pm[i], hm[i] );
    zVec3DOuterProd( &norm, &pm[i], &tempv );
    zVec3DAddDRC( zVec6DAng(c), &tempv );
  }
}

double _rkFDContactConstrainGetArea(zVec3D p[])
{
  zVec3D e1, e2;

  zVec3DSub( &p[1], &p[0], &e1 );
  zVec3DSub( &p[2], &p[0], &e2 );
  return 0.5 * zVec3DOuterProdNorm( &e1, &e2 );
}

bool _rkFDContactConstraintGetSlice_Volume(zVec3D p[], double h[], zVec3D pp[])
{
  double c[2];

  c[0] = h[1] - h[0];
  c[1] = h[2] - h[0];
  if( zIsTiny( c[0] ) || zIsTiny( c[1] ) ) return false;
  zVec3DMul( &p[0], h[1] / c[0], &pp[0] );
  zVec3DCatDRC( &pp[0], -h[0] / c[0], &p[1] );
  zVec3DMul( &p[0], h[2] / c[1], &pp[1] );
  zVec3DCatDRC( &pp[1], -h[0] / c[1], &p[2] );
  h[1] = h[2] = 0;
  return true;
}

void _rkFDContactConstraintGetInnerPoint_Volume(zVec3D *p1, zVec3D *p2, double h1, double h2, zVec3D *pp)
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
  zVec3D tempv;

  zVec3DCat( norm, -zVec3DInnerProd( norm, &cdp->data.norm ), &cdp->data.norm, &tempv );
  if( zVec3DIsTiny( &tempv ) ) return;
  cdpl = zAlloc( rkCDPlane, 1 );
  zVec3DCopy( p, &cdpl->data.v );
  zVec3DDiv( &tempv, -zVec3DNorm( &tempv ) , &cdpl->data.norm );

  /* to eliminate identical conditions */
  /* NOTE: this is O(n^2) and might be able to be O(n) */
  zListForEach( &cdp->data.cplane, cdpl2 ){
    if( zVec3DIsTol( zVec3DSub( &cdpl->data.norm, &cdpl2->data.norm, &tempv ), 1e-8 ) &&
        zIsTol( zVec3DInnerProd( &cdpl->data.norm, zVec3DSub( &cdpl2->data.v, &cdpl->data.v, &tempv ) ), 1e-8 ) ){
      zFree( cdpl );
      return;
    }
  }
  zListInsertHead( &cdp->data.cplane, cdpl );
}

void _rkFDContactConstraint_Volume(rkFD *fd, rkCDPair *cdp, rkContactInfo *ci, zMat6D *Q, zVec6D *C)
{
  register int i, j;
  zTri3D *face;
  double h[3], hm[3], hc, S;
  zVec3D p[3], pm[3], pp[2], pc, tempv;
  zMat3D mmo[3], tempm;
  zVec6D Cc;
  int stp[3], st;

  zMat6DClear( Q );
  zVec6DClear( C );

  for( i=0; i<zPH3DFaceNum(&cdp->data.colvol); i++ ){
    face = zPH3DFace( &cdp->data.colvol, i );
    for( j=0; j<3; j++ ){
      zVec3DSub( zTri3DVert(face,j), &cdp->data.center, &tempv );
      h[j] = zVec3DInnerProd( &cdp->data.norm, &tempv );
    }
    for( j=0; j<3; j++ ){
      zVec3DCat( zTri3DVert(face,j), -h[j], &cdp->data.norm, &p[j] );
    }
    S = _rkFDContactConstrainGetArea( p );
    _rkFDContactConstraintGetMidDepth_Volume( h, ci, S, hm, &hc );
    _rkFDContactConstraintGetMidPoint_Volume( p, pm );
    zVec3DAdd( &p[0], &p[1], &pc );
    zVec3DAddDRC( &pc, &p[2] );
    zVec3DDivDRC( &pc, 3.0 );

    /* Q : needs pm, pc and S */
    zMat3DMul( ZMAT3DIDENT, S, &tempm );
    zMat3DAddDRC( zMat6DMat3D( Q, 0, 0 ), &tempm );
    zVec3DMulDRC( &pc, S );
    zVec3DOuterProdMat3D( &pc, &tempm );
    zMat3DAddDRC( zMat6DMat3D(Q,1,0), &tempm );
    zMat3DSubDRC( zMat6DMat3D(Q,0,1), &tempm );
    for( j=0; j<3; j++ )
      zVec3DOuterProd2Mat3D( &pm[j], &pm[j], &mmo[j] );
    zMat3DAdd( &mmo[0], &mmo[1], &tempm );
    zMat3DAddDRC( &tempm, &mmo[2] );
    zMat3DMulDRC( &tempm, S / 3.0 );
    zMat3DSubDRC( zMat6DMat3D(Q,1,1), &tempm );
    /* C */
    _rkFDContactConstraintDepth_Volume( pm, hm, hc, cdp->data.norm, &Cc );
    /* a process dependent on signum of h */
    st = 0;
    stp[0] = stp[1] = stp[2] = 0;
    for( j=0; j<3; j++ ){
      if( h[j] > zTOL ){
        st += 1<<(j*2);
        stp[1] = j;
      } else if( h[j] < -zTOL ){
        st += 1<<(j*2+1);
        stp[2] = j;
      } else{
        stp[0] = j;
      }
    }
    switch( st ){
    case 0x01: case 0x04: case 0x10:  /* 00上 */
      _rkFDContactSetContactPlane_Volume( cdp, zTri3DVert(face,stp[0]), zTri3DNorm(face) );
    case 0x05: case 0x11: case 0x14: case 0x15: /* 0上上, 上上上 */
      zVec6DAddDRC( C, &Cc );
      continue;
    case 0x02: case 0x08: case 0x20: /* 00下 */
      _rkFDContactSetContactPlane_Volume( cdp, zTri3DVert(face,stp[0]), zTri3DNorm(face) );
    case 0x0a: case 0x22: case 0x28: case 0x2a: /* すべて負 */
      zVec6DSubDRC( C, &Cc );
      continue;
    case 0x24: /* 0上下 */
    case 0x12: /* 下0上 */
    case 0x09: /* 上下0 */
      zVec6DAddDRC( C, &Cc );
      _rkFDContactConstraintGetInnerPoint_Volume( zTri3DVert(face,stp[1]), zTri3DVert(face,stp[2]), h[stp[1]], h[stp[2]], &p[stp[1]] );
      h[stp[1]] = 0.0;
      zVec3DCopy( &p[stp[0]], &pp[0] );
      zVec3DCopy( &p[stp[1]], &pp[1] );
      break;
    case 0x06: /* 下上0 */
    case 0x21: /* 上0下 */
    case 0x18: /* 0下上 */
      zVec6DAddDRC( C, &Cc );
      _rkFDContactConstraintGetInnerPoint_Volume( zTri3DVert(face,stp[1]), zTri3DVert(face,stp[2]), h[stp[1]], h[stp[2]], &p[stp[2]] );
      h[stp[2]] = 0.0;
      zVec3DCopy( &p[stp[2]], &pp[0] );
      zVec3DCopy( &p[stp[0]], &pp[1] );
      break;
    case 0x16: /* 下上上 */
    case 0x19: /* 上下上 */
    case 0x25: /* 上上下 */
      stp[0] = (stp[2]+1) % 3;
      stp[1] = (stp[0]+1) % 3;
      zVec6DAddDRC( C, &Cc );
      _rkFDContactConstraintGetInnerPoint_Volume( zTri3DVert(face,stp[2]), zTri3DVert(face,stp[0]), h[stp[2]], h[stp[0]], &p[stp[0]] );
      _rkFDContactConstraintGetInnerPoint_Volume( zTri3DVert(face,stp[2]), zTri3DVert(face,stp[1]), h[stp[2]], h[stp[1]], &p[stp[1]] );
      h[stp[0]] = h[stp[1]] = 0.0;
      zVec3DCopy( &p[stp[0]], &pp[0] );
      zVec3DCopy( &p[stp[1]], &pp[1] );
      break;
    case 0x1a: /* 下下上 */
    case 0x26: /* 下上下 */
    case 0x29: /* 上下下 */
      stp[0] = (stp[1]+1) % 3;
      stp[2] = (stp[0]+1) % 3;
      zVec6DSubDRC( C, &Cc );
      _rkFDContactConstraintGetInnerPoint_Volume( zTri3DVert(face,stp[1]), zTri3DVert(face,stp[0]), h[stp[1]], h[stp[0]], &p[stp[0]] );
      _rkFDContactConstraintGetInnerPoint_Volume( zTri3DVert(face,stp[1]), zTri3DVert(face,stp[2]), h[stp[1]], h[stp[2]], &p[stp[2]] );
      h[stp[0]] = h[stp[2]] = 0.0;
      zVec3DCopy( &p[stp[2]], &pp[0] );
      zVec3DCopy( &p[stp[0]], &pp[1] );
      break;
    default: /* 面上 */
      continue;
    }
    S = _rkFDContactConstrainGetArea( p );
    _rkFDContactConstraintGetMidDepth_Volume( h, ci, S, hm, &hc );
    _rkFDContactConstraintGetMidPoint_Volume( p, pm );
    _rkFDContactConstraintDepth_Volume( pm, hm, hc, cdp->data.norm, &Cc );
    switch( st ){
    case 0x1a: case 0x26: case 0x29:
      zVec6DCatDRC( C, 2.0, &Cc );
      break;
    default:
      zVec6DCatDRC( C, -2.0, &Cc );
    }
    /* register pp */
    _rkFDContactSetContactPlane_Volume( cdp, &pp[0], zTri3DNorm(face) );
  }
}

void _rkFDContactModForce_Volume(rkFD *fd, bool doUpRef)
{
  rkCDPair *cdp;
  rkContactInfo *ci;
  zVec3D vr, vv[2], tempv, s;
  zVec6D v6;
  register int i;
  double fn, fs;

  zListForEach( &fd->cd.plist, cdp ){
    if( !cdp->data.is_col ) continue;
    ci = rkContactInfoPoolAssoc( &fd->ci, rkLinkStuff(cdp->data.cell[0]->data.link), rkLinkStuff(cdp->data.cell[1]->data.link) );
    if( !ci ) ci = &_rkFDContactInfoDefault;

    /* v */
    for( i=0; i<2; i++ ){
      if( cdp->data.cell[i]->data.type == RK_CD_CELL_STAT )
        zVec3DClear( &vv[i] );
      else{
        zMulMatVec6D( rkLinkWldAtt(cdp->data.cell[i]->data.link), rkLinkVel(cdp->data.cell[i]->data.link), &v6 );
        zVec3DSub( &cdp->data.r, rkLinkWldPos(cdp->data.cell[i]->data.link) ,&tempv );
        zVec3DOuterProd( zVec6DAng(&v6), &tempv, &tempv );
        zVec3DAdd( zVec6DLin(&v6), &tempv, &vv[i] );
      }
    }
    zVec3DSub( &vv[0], &vv[1], &vr );
    zVec3DCatDRC( &vr, -zVec3DInnerProd( &vr, &cdp->data.norm ), &cdp->data.norm );
    if( zVec3DIsTiny( &vr ) ) zVec3DClear( &vr );
    else zVec3DNormalizeNCDRC( &vr );

    /* friction */
    fn = zVec3DInnerProd( &cdp->data.f, &cdp->data.norm );
    zVec3DCat( &cdp->data.f, -fn, &cdp->data.norm, &s );
    fs = zVec3DNorm( &s );
    if( !zIsTiny( fs ) && fs > rkContactInfoSF(ci) * fn ){
      zVec3DMul( &vr, -rkContactInfoKF(ci)*fn, &cdp->data.f );
      zVec3DCatDRC( &cdp->data.f, fn, &cdp->data.norm );
    }
  }
}

void _rkFDContactSetForceF_Volume(rkFD *fd, zVec f)
{
  rkCDPair *cdp;
  int offset =  0;

  zListForEach( &fd->cd.plist, cdp ){
    if( !cdp->data.is_col ) continue;
    zVec3DCopy( (zVec3D *)&zVecElem(f,offset), &cdp->data.f );
    if( zVec3DIsTiny( &cdp->data.f ) || zIsTiny( zVec3DInnerProd( &cdp->data.f, &cdp->data.norm ) ) ){
      zVec3DClear( &cdp->data.f );
    }
    offset += 3;
  }
}

void _rkFDContactSetForceR_Volume(rkFD *fd, zVec r)
{
  rkCDPair *cdp;
  int offset =  0;

  zListForEach( &fd->cd.plist, cdp ){
    if( !cdp->data.is_col ) continue;
    zVec3DCopy( &cdp->data.center, &cdp->data.r );
    zVec3DCatDRC( &cdp->data.r, zVecElem(r,offset  ), &cdp->data.axis[1] );
    zVec3DCatDRC( &cdp->data.r, zVecElem(r,offset+1), &cdp->data.axis[2] );
    offset += 2;
  }
}

zVec _rkFDContactSolveConstraintInit_Volume(zMat a, zVec b, zVec ans, void *util){
  zVecClear( ans );
  return ans;
}

void _rkFDContactSolveConstraintQc4F_Volume(rkFD *fd, zMat Q, zVec c, zMat Qr, zVec cr)
{
  rkCDPair *cdp, *cdp2;
  zMat6D Qcell;
  zMat3D Qp, r1, r2, tempm;
  zVec3D tempv;
  register int i, j;
  int ii, jj;

  ii = 0;
  zListForEach( &fd->cd.plist, cdp ){
    if( !cdp->data.is_col ) continue;
    jj = 0;
    /* Qr */
    zListForEach( &fd->cd.plist, cdp2 ){
      if( !cdp2->data.is_col ) continue;
      for( i=0; i<6; i++ ){
        for( j=0; j<6; j++ ){
          zMat6DElem(&Qcell,i,j) = zMatElem(Q,6*ii+i,6*jj+j);
        }
      }
      zMat3DCopy( zMat6DMat3D(&Qcell,0,0), &Qp );
      zVec3DOuterProdMat3D( &cdp->data.r, &r1 );
      zMulMatMat3D( &r1, zMat6DMat3D(&Qcell,1,0), &tempm );
      zMat3DSubDRC( &Qp, &tempm );
      zVec3DOuterProdMat3D( &cdp2->data.r, &r2 );
      zMulMatMat3D( zMat6DMat3D(&Qcell,0,1), &r2, &tempm );
      zMat3DAddDRC( &Qp, &tempm );
      zMulMatMat3D( zMat6DMat3D(&Qcell,1,1), &r2, &tempm );
      zMulMatMat3D( &r1, &tempm, &tempm );
      zMat3DSubDRC( &Qp, &tempm );
      for( i=0; i<3; i++ )
        for( j=0; j<3; j++ )
          zMatElem(Qr,3*ii+i,3*jj+j) = zMat3DElem(&Qp,i,j);
      jj++;
    }
    /* cr */
    zVec3DOuterProd( (zVec3D *)&zVecElem(c,6*ii+3), &cdp->data.r, &tempv );
    zVec3DAdd( (zVec3D *)&zVecElem(c,6*ii), &tempv, (zVec3D *)&zVecElem(cr,3*ii) );
    ii++;
  }
}

void _rkFDContactSolveConstraintQc4R_Volume(rkFD *fd, zMat Q, zVec c, zMat Qr, zVec cr)
{
  rkCDPair *cdp, *cdp2;
  zMat3D Qcell, Qcell2;
  zVec3D d1[2], d2[2], tempv, tempv2;
  register int i, j;
  int ii, jj, iq, jq;

  ii = iq = 0;
  zVecClear( cr );
  zListForEach( &fd->cd.plist, cdp ){
  if( !cdp->data.is_col ) continue;
    if( zListNum(&cdp->data.cplane) == 0 ){
      ii++;
      continue;
    }
    jj = jq = 0;
    /* Qr */
    zVec3DOuterProd( &cdp->data.axis[1], &cdp->data.f, &d1[0] );
    zVec3DOuterProd( &cdp->data.axis[2], &cdp->data.f, &d1[1] );
    zVec3DMulDRC( &d1[0], fd->dt );
    zVec3DMulDRC( &d1[1], fd->dt );
    zListForEach( &fd->cd.plist, cdp2 ){
      if( !cdp2->data.is_col ) continue;
      if( zListNum(&cdp2->data.cplane) == 0 ){
        jj++;
        continue;
      }
      for( i=0; i<3; i++ )
        for( j=0; j<3; j++ )
          zMat3DElem(&Qcell,i,j) = zMatElem(Q,6*ii+3+i,6*jj+3+j);
      zVec3DOuterProd( &cdp2->data.axis[1], &cdp2->data.f, &d2[0] );
      zVec3DOuterProd( &cdp2->data.axis[2], &cdp2->data.f, &d2[1] );
      zVec3DMulDRC( &d2[0], fd->dt );
      zVec3DMulDRC( &d2[1], fd->dt );
      for( i=0; i<2; i++ )
        for( j=0; j<2; j++ ){
          zMulMatVec3D( &Qcell, &d2[j], &tempv );
          zMatElem(Qr,2*iq+i,2*jq+j) = zVec3DInnerProd( &d1[i], &tempv );
        }
      /* cr */
      for( i=0; i<3; i++ )
        for( j=0; j<3; j++ )
          zMat3DElem(&Qcell2,i,j) = zMatElem(Q,6*ii+3+i,6*jj+j);
      zVec3DOuterProd( &cdp2->data.center, &cdp2->data.f, &tempv );
      zVec3DMulDRC( &tempv, fd->dt );
      zMulMatVec3DDRC( &Qcell, &tempv );
      zMulMatVec3D( &Qcell2, &cdp2->data.f, &tempv2 );
      zVec3DMulDRC( &tempv2, fd->dt );
      zVec3DAddDRC( &tempv, &tempv2 );

      zVecElem(cr,2*iq  ) += zVec3DInnerProd( &tempv, &d1[0] );
      zVecElem(cr,2*iq+1) += zVec3DInnerProd( &tempv, &d1[1] );
      jj++;
      jq++;
    }
    /* cr */
    zVecElem(cr,2*iq  ) += zVec3DInnerProd( (zVec3D *)&zVecElem(c,6*ii+3), &d1[0] );
    zVecElem(cr,2*iq+1) += zVec3DInnerProd( (zVec3D *)&zVecElem(c,6*ii+3), &d1[1] );
    ii++;
    iq++;
  }
}

void _rkFDContactSolveConstraint_Volume(rkFD *fd, zMat A, zVec b)
{
  zMat Q, Qf, Qr, Nf, Nr;
  zVec c, cf, cr, df, dr, f, fp, r, rp;
  zMat6D Qv;
  zVec6D Cv, tempv;
  rkCDPair *cdp;
  rkCDPlane *cdpl;
  rkContactInfo *ci;
  register int i, j;
  int offset = 0, cnt = 0;
  int m, n;

  Q = zMatAllocSqr( 6*fd->cd.colnum );
  c = zVecAlloc( 6*fd->cd.colnum );
  zMatClear( Q );
  zVecClear( c );
  zListForEach( &fd->cd.plist, cdp ){
    if( !cdp->data.is_col ) continue;
    ci = rkContactInfoPoolAssoc( &fd->ci, rkLinkStuff(cdp->data.cell[0]->data.link), rkLinkStuff(cdp->data.cell[1]->data.link) );
    if( !ci ) ci = &_rkFDContactInfoDefault;
    _rkFDContactConstraint_Volume( fd, cdp, ci, &Qv, &Cv );
    /* Q */
    for( i=0; i<6; i++ )
      for( j=0; j<6; j++ )
        zRawMatCatDyad( zMatBuf(Q), zMat6DElem(&Qv,i,j), zMatRowBuf(A,offset+i), zMatColSizeNC(A), zMatRowBuf(A,offset+j), zMatColSizeNC(A) );
    /* c */
    zMulMat6DVec6D( &Qv, (zVec6D *)&zVecElem(b,offset), &tempv );
    zVec6DAddDRC( &Cv, &tempv );
    for( i=0; i<6; i++ )
      zRawVecCatDRC( zVecBuf(c), zVec6DElem(&Cv,i), zMatRowBuf(A,offset+i), zMatColSizeNC(A) );
    offset += 6;
  }
  /* iterative computation of force and moment arm r */
  /* preprocess */
  /* f */
  Qf = zMatAllocSqr( 3*fd->cd.colnum );
  cf = zVecAlloc( 3*fd->cd.colnum );
  Nf = zMatAlloc( fd->cd.colnum, 3*fd->cd.colnum );
  df = zVecAlloc( fd->cd.colnum );
  f  = zVecAlloc( 3*fd->cd.colnum );
  fp = zVecAlloc( 3*fd->cd.colnum );

  i = m = n = 0;
  zListForEach( &fd->cd.plist, cdp ){
    if( !cdp->data.is_col ) continue;
    zVec3DCopy( &cdp->data.center, &cdp->data.r );
    zVec3DCopy( &cdp->data.norm, (zVec3D *)&zMatElem(Nf,i,3*i) );
    i++;
    m += zListNum( &cdp->data.cplane );
    if( zListNum(&cdp->data.cplane) != 0 ) n++;
  }
  /* r */
  if( m != 0 ){
    /* 接触しているが断面が得られない場合の自由度は消す */
    Qr = zMatAllocSqr( 2*n );
    cr = zVecAlloc( 2*n );
    Nr = zMatAlloc( m, 2*n );
    dr = zVecAlloc( m );
    r  = zVecAlloc( 2*n );
    rp = zVecAlloc( 2*n );

    i = j = 0;
    zListForEach( &fd->cd.plist, cdp ){
      if( !cdp->data.is_col || zListNum(&cdp->data.cplane) == 0 ) continue;
      zListForEach( &cdp->data.cplane, cdpl ){
        zMatElem(Nr,i,2*j  ) = zVec3DInnerProd( &cdpl->data.norm, &cdp->data.axis[1] );
        zMatElem(Nr,i,2*j+1) = zVec3DInnerProd( &cdpl->data.norm, &cdp->data.axis[2] );
        zVecElem(dr,i) = zVec3DInnerProd( &cdpl->data.norm, &cdpl->data.v ) - zVec3DInnerProd( &cdpl->data.norm, &cdp->data.center );
        /* for safety */
        if( zVecElem(dr,i) > 0.0 ){
          zMatElem(Nr,i,2*j  ) *= -1.0;
          zMatElem(Nr,i,2*j+1) *= -1.0;
          zVecElem(dr,i) *= -1.0;
        }
        i++;
      }
      j++;
    }
  }
  /* iteration */
  while( 1 ){
    /* f */
    _rkFDContactSolveConstraintQc4F_Volume( fd, Q, c, Qf, cf );
    offset = 0;
    zListForEach( &fd->cd.plist, cdp ){
      if( !cdp->data.is_col ) continue;
      ci = rkContactInfoPoolAssoc( &fd->ci, rkLinkStuff(cdp->data.cell[0]->data.link), rkLinkStuff(cdp->data.cell[1]->data.link) );
      if( !ci ) ci = &_rkFDContactInfoDefault;
      for( i=0; i<3; i++ )
        zMatElem(Qf,offset+i,offset+i) += rkContactInfoL(ci);
      offset += 3;
    }
    zQPSolveASM( Qf, cf, Nf, df, f, NULL, NULL, NULL );
    zVecDivDRC( f, fd->dt );

    if( m == 0 || cnt++ > RK_FD_SOLVE_ITER ) break;
    if( zVecIsEqual( f, fp ) && zVecIsEqual( r, rp ) ) break;
    zVecCopy( f, fp );
    zVecCopy( r, rp );

    /* correction of friction force */
    _rkFDContactModForce_Volume( fd, false );
    _rkFDContactSetForceF_Volume( fd, f );

    /* r */
    _rkFDContactSolveConstraintQc4R_Volume( fd, Q, c, Qr, cr );
    zQPSolveASM( Qr, cr, Nr, dr, r, NULL, _rkFDContactSolveConstraintInit_Volume, NULL );
    _rkFDContactSetForceR_Volume( fd, r );
  }

  zMatFreeAO( 3, Q, Qf, Nf );
  zVecFreeAO( 5, c, cf, df, f, fp );
  if( m != 0 ){
    zMatFreeAO( 2, Qr, Nr );
    zVecFreeAO( 4, cr, dr, r, rp );
  }
}

void _rkFDContactCalcForce_Volume(rkFD *fd)
{
  zMat A;
  zVec b;

  A = zMatAllocSqr( 6*(fd->cd.colnum) );
  b = zVecAlloc( 6*(fd->cd.colnum) );
  _rkFDContactRelationAccForce_Volume( fd, A, b );
  _rkFDContactVelAfterInt_Volume( fd, A, b );
  _rkFDContactSolveConstraint_Volume( fd, A, b );
  zMatFree( A );
  zVecFree( b );
}

void _rkFDContactSetWrench(rkFD *fd)
{
  rkCDPair *cdp;
  rkWrench *w;
  register int i;

  zListForEach( &fd->cd.plist, cdp ){
    if( !cdp->data.is_col ) continue;
    for( i=0; i<2; i++ ){
      w = zAlloc( rkWrench, 1 );
      rkWrenchInit( w );
      zXfer3DInv( rkLinkWldFrame(cdp->data.cell[i]->data.link), &cdp->data.r, rkWrenchPos(w) );
      zMulMatTVec3D( rkLinkWldAtt(cdp->data.cell[i]->data.link), &cdp->data.f, rkWrenchForce(w) );
      rkLinkExtWrenchPush( cdp->data.cell[i]->data.link, w );
      zVec3DRevDRC( &cdp->data.f );
    }
  }
}

void _rkFDSolveJointContact_Volume(rkFD *fd, bool doUpRef)
{
  /* update CD */
  zEchoOff();
  rkCDColVolBREP( &fd->cd );
  zEchoOn();

  /* compute joint friction */
  _rkFDJointCalcFriction( fd, doUpRef );
  if( fd->cd.colnum == 0 ) return;
  /* compute contact force */
  _rkFDContactCalcForce_Volume( fd );
  /* modify contact force */
  _rkFDContactModForce_Volume( fd, doUpRef );
  _rkFDContactSetWrench( fd );
}

zVec _rkFDUpdate_Volume(double t, zVec dis, zVec vel, void *fd, zVec acc)
{
  zVecClear( acc );
  _rkFDChainConnectJointState( fd, dis, vel, acc );
  _rkFDChainExtWrenchDestroy( fd );

  _rkFDSolveJointContact_Volume( fd, false );
  zVecClear( acc );
  _rkFDUpdateAcc( fd );
  return acc;
}

void _rkFDUpdateRef_Volume(rkFD *fd)
{
  zVecClear( fd->acc );
  _rkFDChainConnectJointState( fd, fd->dis, fd->vel, fd->acc );
  _rkFDChainExtWrenchDestroy( fd );

  _rkFDSolveJointContact_Volume( fd, true );
  _rkFDUpdateAcc( fd );
  zVecClear( fd->acc );
  _rkFDUpdateRefDrivingTorque( fd );
}

/* rkFDSolve
 * - solve the forward dynamics.
 */
rkFD *rkFDSolve_Volume(rkFD *fd)
{
  if( zIsTiny( fd->t ) ) _rkFDJointRefInit( fd );
  /* reference */
  _rkFDUpdateRef_Volume( fd );

  /* integration */
  zODE2Init( &fd->_ode, zVecSize(fd->dis), fd->_ode_step, _rkFDUpdate_Volume );
  zODE2Update( &fd->_ode, fd->t, fd->dis, fd->vel, fd->dt, fd );
  zODE2Destroy( &fd->_ode );
  fd->t += fd->dt;
  /* calc acc */
  _rkFDUpdate_Volume( fd->t, fd->dis, fd->vel, fd, fd->acc );
  return fd;
}

void rkFDUpdateInit_Volume(rkFD *fd)
{
  zODE2Init( &fd->_ode, zVecSize(fd->dis), fd->_ode_step, _rkFDUpdate_Volume );
  _rkFDJointRefInit( fd );
  _rkFDUpdateRef_Volume( fd );
}

rkFD *rkFDUpdate_Volume(rkFD *fd)
{
  zODE2Update( &fd->_ode, fd->t, fd->dis, fd->vel, fd->dt, fd );
  fd->t += fd->dt;
  _rkFDUpdateRef_Volume( fd );
  return fd;
}

void rkFDUpdateDestroy_Volume(rkFD *fd)
{
  zODE2Destroy( &fd->_ode );
}

/* for debug */
void rkFDWrite(rkFD *fd)
{
  rkFDCell *lc;

  zListForEach( &fd->list, lc )
    rkChainWrite( &lc->data.chain );
}
