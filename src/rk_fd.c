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
  zVec3D d, vr, vv[2], tempv;
  zVec6D v6;

  zVecMulDRC( b, fd->dt );
  _rkFDUpdateAcc( fd );
  zListForEach( &fd->cd.plist, cdp ){
    if( !cdp->data.is_col ) continue;
    ci = rkContactInfoPoolAssoc( &fd->ci, rkLinkStuff( cdp->data.cell[0]->data.link ), rkLinkStuff( cdp->data.cell[1]->data.link ) );
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
    ci = rkContactInfoPoolAssoc( &fd->ci, rkLinkStuff( cdp->data.cell[0]->data.link ), rkLinkStuff( cdp->data.cell[1]->data.link ) );
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
  zVec3D f, vr, vv[2], tempv;
  zVec6D v6;
  rkWrench *w;
  register int i;

  _rkFDChainExtWrenchDestroy( fd );
  zListForEach( &fd->cd.plist, cdp ){
    if( !cdp->data.is_col ) continue;
    ci = rkContactInfoPoolAssoc( &fd->ci, rkLinkStuff( cdp->data.cell[0]->data.link ), rkLinkStuff( cdp->data.cell[1]->data.link ) );
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
      zVec3DCatDRC( &vr, -1.0*zVec3DInnerProd( &vr, &cdv->data.axis[0] ), &cdv->data.axis[0] );
      if( zIsTiny( zVec3DNorm( &vr ) ) ) zVec3DClear( &vr );
      else zVec3DNormalizeDRC( &vr );

      fn = zVec3DElem( &cdv->data.f, 0 );
      fs = sqrt( zSqr( zVec3DElem(&cdv->data.f,1) ) + zSqr( zVec3DElem( &cdv->data.f, 2 ) ) );

      /* friction corn */
      if( !zIsTiny( fs ) &&
          fs > ( cdv->data.type == RK_CONTACT_STICK ? rkContactInfoSF(ci) : rkContactInfoKF(ci) )*fn ){
        zVec3DMul( &vr, -1.0*rkContactInfoKF(ci)*fn, &f );
        zVec3DCatDRC( &f, fn, &cdv->data.axis[0] );
        if( doUpRef ){
          cdv->data.type = RK_CONTACT_SLIP;
          zVec3DCopy( &cdv->data._pro, &cdv->data._ref );
        }
      } else{
        zVec3DClear( &f );
        for( i=0; i<3; i++ ){
          zVec3DCatDRC( &f, zVec3DElem(&cdv->data.f,i), &cdv->data.axis[i] );
        }
        if( doUpRef ){
          cdv->data.type = RK_CONTACT_STICK;
        }
      }

      /* set wrench */
      for( i=0; i<2; i++){
        w = zAlloc( rkWrench, 1 );
        rkWrenchInit( w );
        zXfer3DInv( rkLinkWldFrame( celld[i]->link ), cdv->data.vert, rkWrenchPos( w ) );
        zMulMatTVec3D( rkLinkWldAtt( celld[i]->link ), &f, rkWrenchForce( w ) );
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

/* rkFDSolve
 * - solve the forward dynamics.
 */
rkFD *rkFDSolve(rkFD *fd)
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

void rkFDUpdateInit(rkFD *fd)
{
  zODE2Init( &fd->_ode, zVecSize( fd->dis ), fd->_ode_step, _rkFDUpdate );
  _rkFDJointRefInit( fd );
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

/* for debug */
void rkFDWrite(rkFD *fd)
{
  rkFDCell *lc;

  zListForEach( &fd->list, lc )
    rkChainWrite( &lc->data.chain );
}
