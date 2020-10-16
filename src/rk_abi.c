/* RoKi - Robot Kinetics library
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * rk_abi - Articulated Body Inertia Method.
 * contributer: 2014-2015 Naoki Wakisaka
 */

#include <roki/rk_abi.h>

/* (static)
 * _rkLinkABIInitInertia
 * - initialize invariant mass properties.
 */
void _rkLinkABIInitInertia(rkLink *link)
{
  rkABIPrp *ap;
  zMat3D pcros;

  ap = rkLinkABIPrp(link);
  zMat3DMul( ZMAT3DIDENT, rkLinkMass(link), zMat6DMat3D(&ap->m,0,0) );
  zVec3DOuterProdMat3D( rkLinkCOM(link), &pcros );
  zMat3DMul( &pcros, rkLinkMass(link), zMat6DMat3D(&ap->m,1,0) );
  zMat3DRev( zMat6DMat3D(&ap->m,1,0), zMat6DMat3D(&ap->m,0,1) );
  zMulMatMat3D( &pcros, zMat6DMat3D(&ap->m,0,1), zMat6DMat3D(&ap->m,1,1) );
  zMat3DAddDRC( zMat6DMat3D(&ap->m,1,1), rkLinkInertia(link) );
}

void rkLinkABIInit(rkLink *link)
{
  rkABIPrp *ap;

  ap = rkLinkABIPrp(link);
  memset( ap, 0, sizeof(rkABIPrp) );
  if( rkLinkJointSize(link) == 0 ){
    ap->axi = ap->iaxi = NULL;
  } else{
    ap->axi  = zMatAllocSqr( rkLinkJointSize(link) );
    ap->iaxi = zMatAllocSqr( rkLinkJointSize(link) );
  }
  _rkLinkABIInitInertia( link );
}

void rkChainABIInit(rkChain *chain)
{
  register int i;

  for( i=0; i<rkChainNum(chain); i++ )
    rkLinkABIInit( rkChainLink(chain,i) );
}

/* ABI method */
void rkLinkABIUpdateInit(rkLink *link, zVec6D *pvel)
{
  rkABIPrp *ap;
  zMat3D pcrs;
  zVec3D tempv;

  ap = rkLinkABIPrp(link);

  /*J*/
  zMat6DPutMat3D( &ap->j, 0, 0, rkLinkAdjAtt(link) );
  zMat6DPutMat3D( &ap->j, 0, 1, ZMAT3DZERO );
  zVec3DOuterProdMat3D( rkLinkAdjPos(link), &pcrs );
  zMulMatMat3D( &pcrs, rkLinkAdjAtt(link), zMat6DMat3D(&ap->j,1,0) );
  zMat6DPutMat3D( &ap->j, 1, 1, rkLinkAdjAtt(link) );

  /*vel update*/
  zXfer6DLin( rkLinkAdjFrame(link), pvel, rkLinkVel(link) );
  zMulMatTVec3D( rkLinkAdjAtt( link ), zVec6DAng( pvel ), &tempv );
  zVec6DClear( &ap->c );
  rkJointIncVel( rkLinkJoint( link ), &tempv, rkLinkVel( link ), &ap->c );

  /*I*/
  zMat6DCopy( &ap->m, &ap->i );

  /*b*/
  zVec3DTripleProd( rkLinkAngVel( link ), rkLinkAngVel( link ), rkLinkCOM( link ), &tempv);
  zVec3DMul( &tempv, rkLinkMass( link ), zVec6DLin( &ap->b ) );
  zMulMatVec3D( zMat6DMat3D( &ap->i, 1, 1 ), rkLinkAngVel( link ), &tempv );
  zVec3DOuterProd( rkLinkAngVel( link ), &tempv, zVec6DAng( &ap->b ) );

  /*c*/
  zVec3DTripleProd( zVec6DAng( pvel ), zVec6DAng( pvel ), rkLinkAdjPos( link ), &tempv );
  zMulMatTVec3DDRC( rkLinkAdjAtt( link ), &tempv );
  zVec3DAddDRC( zVec6DLin( &ap->c ), &tempv );

  /* recursive initialization of ABI */
  if( rkLinkSibl( link ) )
    rkLinkABIUpdateInit( rkLinkSibl( link ), pvel );
  if( rkLinkChild( link ) )
    rkLinkABIUpdateInit( rkLinkChild( link ), rkLinkVel( link ) );
}

void rkLinkABIUpdateBackward(rkLink *link)
{
  rkABIPrp *ap, *pap;
  zVec6D icb;
  rkWrench *w;
  zVec3D tempv;

  /* recursive update of ABI */
  if( rkLinkSibl( link ) )
    rkLinkABIUpdateBackward( rkLinkSibl(link) );
  if( rkLinkChild( link ) )
    rkLinkABIUpdateBackward( rkLinkChild(link) );

  ap = rkLinkABIPrp(link);
  /* external forces */
  zListForEach( rkLinkExtWrench( link ), w ){
    zVec6DSubDRC( &ap->b, rkWrenchW( w ) );
    zVec3DOuterProd( rkWrenchPos( w ), rkWrenchForce( w ), &tempv );
    zVec3DSubDRC( zVec6DAng( &ap->b ), &tempv );
  }

  /* gravity force */
  zVec3DClear( &tempv );
  zVec3DElem(&tempv,zZ) = -RK_G * rkLinkMass(link);
  zMulMatTVec3DDRC( rkLinkWldAtt(link), &tempv );
  zVec3DSubDRC( zVec6DLin(&ap->b), &tempv );
  zVec3DOuterProd( rkLinkCOM(link), &tempv, &tempv );
  zVec3DSubDRC( zVec6DAng(&ap->b), &tempv );

  /* IsIs */
  if( ap->axi ){
    rkJointABIAxisInertia( rkLinkJoint( link ), &ap->i, ap->axi );
    zMatClear( ap->iaxi );
    rkJointMotorInertia( rkLinkJoint( link ), zMatBuf( ap->iaxi ) );
    zMatAddDRC( ap->axi, ap->iaxi );
    zMatInv( ap->axi, ap->iaxi );
  }

  if( !rkLinkParent(link) ) return;

  /* add ABI and bias acceleration to parent prp */
  /* Is (sIs)-1 sT I */
  /* Is (sIs)-1 (Q - sT (Ic + b)) */
  /*Ic+b*/
  zMulMat6DVec6D( &ap->i, &ap->c, &icb );
  zVec6DAddDRC( &icb, &ap->b );

  pap = rkLinkABIPrp( rkLinkParent( link ) );
  /* todo friction */
  rkJointABIAddAbiBios( rkLinkJoint( link ), &ap->i, &ap->j, &icb, ap->iaxi, &pap->i, &pap->b );
}

void rkLinkABIUpdateForward(rkLink *link, zVec6D *pa)
{
  rkABIPrp *ap;
  zVec6D jac;
  zMat3D att;

  ap = rkLinkABIPrp( link );
  /*J^Ta+c*/
  zMulMat6DTVec6D( &ap->j, pa, &jac );
  zVec6DAddDRC( &jac, &ap->c );

  /* q, acc update */
  /* todo friction */
  zMulMatTMat3D(rkLinkOrgAtt(link), rkLinkAdjAtt(link), &att);
  rkJointABIQAcc( rkLinkJoint( link ), &att, &ap->i, &ap->b, &jac, ap->iaxi, rkLinkAcc( link ) );

  /* recursive forward computation of ABI */
  if( rkLinkSibl( link ) )
    rkLinkABIUpdateForward( rkLinkSibl( link ), pa );
  if( rkLinkChild( link ) )
    rkLinkABIUpdateForward( rkLinkChild( link ), rkLinkAcc( link ) );
}

void rkChainABIUpdate(rkChain *chain)
{
  if( rkChainJointSize(chain) == 0 ) return;
  rkChainUpdateFK( chain );
  rkChainABIUpdateInit( chain );
  rkChainABIUpdateBackward( chain );
  rkChainABIUpdateForward( chain );
}

zVec rkChainABI(rkChain *chain, zVec dis, zVec vel, zVec acc)
{
  if( rkChainJointSize(chain) == 0 ) return NULL;
  rkChainSetJointDisAll( chain, dis );
  rkChainSetJointVelAll( chain, vel );
  rkChainABIUpdate( chain );
  rkChainGetJointAccAll( chain, acc );
  return acc;
}
