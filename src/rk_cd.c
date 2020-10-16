/* RoKi - Robot Kinetics library
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * rk_cd - collision detection
 * contributer: 2014-2015 Ken'ya Tanaka
 */

#include <roki/rk_cd.h>

static void _rkCDPairDestroy(rkCDPair *pair);

/* (static)
 * _rkCDCellInit
 * - initialize a collision detection cell.
 */
void _rkCDCellInit(rkCDCell *cell)
{
  cell->data.shape = NULL;
  cell->data.link = NULL;
  cell->data.chain = NULL;
  cell->data.type = RK_CD_CELL_STAT;
  zAABox3DInit( &cell->data.aabb );
  zBox3DInit( &cell->data.obb );
  zPH3DInit( &cell->data.ph );
}

/* (static)
 * _rkCDCellCreate
 * - create a collision detection cell.
 */
rkCDCell *_rkCDCellCreate(rkCDCell *cell, rkChain *chain, rkLink *link, zShape3D *shape, rkCDCellType type)
{
  _rkCDCellInit( cell );
  cell->data.shape = shape;
  cell->data.link = link;
  cell->data.chain = chain;
  cell->data.type = type;
  cell->data._ph_update_flag = 0;
  cell->data._bb_update_flag = 0;
  /* convert the original shape to a polyhedron */
  if( zShape3DType(cell->data.shape) != ZSHAPE_PH )
    if( !zShape3DToPH( cell->data.shape ) ) return NULL;
  /* create the bounding box of the shape */
  if( zBox3DDepth(zShape3DBB(cell->data.shape)) == 0 &&
      zBox3DWidth(zShape3DBB(cell->data.shape)) == 0 &&
      zBox3DHeight(zShape3DBB(cell->data.shape)) == 0 )
    zOBB( zShape3DBB(cell->data.shape), zShape3DVertBuf(cell->data.shape), zShape3DVertNum(cell->data.shape) );

  if( !zPH3DClone( zShape3DPH(cell->data.shape), &cell->data.ph ) )
    return NULL;
  zPH3DXfer( zShape3DPH(cell->data.shape), rkLinkWldFrame(cell->data.link), &cell->data.ph );

  zBox3DXfer( zShape3DBB(cell->data.shape), rkLinkWldFrame(cell->data.link), &cell->data.obb );
  zBox3DToAABox3D( &cell->data.obb, &cell->data.aabb );
  return cell;
}

/* (static)
 * _rkCDCellDestroy
 * - destroy a collision detection cell.
 */
void _rkCDCellDestroy(rkCDCell *cell)
{
  zPH3DDestroy( &cell->data.ph );
  _rkCDCellInit( cell );
}

/* rkCDCellReg
 * - register a collision detection cell.
 */
rkCDCell *rkCDCellReg(rkCDCellList *clist, rkChain *chain, rkLink *link, zShape3D *shape, rkCDCellType type)
{
  rkCDCell *cell;

  if( !( cell = zAlloc( rkCDCell, 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  if( !_rkCDCellCreate( cell, chain, link, shape, type ) ) return NULL;
  zListInsertHead( clist, cell );
  return cell;
}

/* rkCDCellUpdateBB
 * - update the bounding box of a collision detection cell.
 */
void rkCDCellUpdateBB(rkCDCell *cell)
{
  if( cell->data.type == RK_CD_CELL_STAT ||
      cell->data._bb_update_flag == 1) return;
  zBox3DXfer( zShape3DBB(cell->data.shape), rkLinkWldFrame(cell->data.link), &cell->data.obb );
  zBox3DToAABox3D( &cell->data.obb, &cell->data.aabb );
  cell->data._bb_update_flag = 1;
}

/* rkCDCellUpdatePH
 * - update the polyhedron of a collision detection cell.
 */
void rkCDCellUpdatePH(rkCDCell *cell)
{
  if( cell->data.type == RK_CD_CELL_STAT ||
      cell->data._ph_update_flag == 1) return;
  zPH3DXfer( zShape3DPH(cell->data.shape), rkLinkWldFrame(cell->data.link), &cell->data.ph );
  cell->data._ph_update_flag = 1;
}

/* rkCDCellUpdate
 * - update a cellcollision detection cell.
 */
void rkCDCellUpdate(rkCDCell *cell)
{
  if( cell->data.type == RK_CD_CELL_STAT ) return;
  if( cell->data._bb_update_flag == 0 ){
    zBox3DXfer( zShape3DBB(cell->data.shape), rkLinkWldFrame(cell->data.link), &cell->data.obb );
    zBox3DToAABox3D( &cell->data.obb, &cell->data.aabb );
    cell->data._bb_update_flag = 1;
  }
  if( cell->data._ph_update_flag == 0 ){
    zPH3DXfer( zShape3DPH(cell->data.shape), rkLinkWldFrame(cell->data.link), &cell->data.ph );
    cell->data._ph_update_flag = 1;
  }
}

/* rkCDCreate
 * - create a collision detector.
 */
rkCD *rkCDCreate(rkCD *cd)
{
  zListInit( &cd->clist );
  zListInit( &cd->plist );
  cd->colnum = 0;
  return cd;
}

/* rkCDDestroy
 * - destroy a collision detector.
 */
void rkCDDestroy(rkCD *cd)
{
  rkCDCell *cp;
  rkCDPair *pair;

  while( !zListIsEmpty( &cd->clist ) ){
    zListDeleteHead( &cd->clist, &cp );
    _rkCDCellDestroy( cp );
    zFree( cp );
  }
  while( !zListIsEmpty( &cd->plist ) ){
    zListDeleteHead( &cd->plist, &pair );
    _rkCDPairDestroy( pair );
    zFree( pair );
  }
}

/* (static)
 * _rkCDPairDestroy
 * - destroy a pair of collision detection cells.
 */
void _rkCDPairDestroy(rkCDPair *pair)
{
  zListDestroy( rkCDVert, &pair->data.vlist );
  zPH3DDestroy( &pair->data.colvol );
}

/* (static)
 * _rkCDPairReg
 * - register a pair of collision detection cells.
 */
rkCDPair *_rkCDPairReg(rkCD *cd, rkCDCell *c1, rkCDCell *c2)
{
  rkCDPair *pair;

  if( !( pair = zAlloc( rkCDPair, 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  pair->data.cell[0] = c1;
  pair->data.cell[1] = c2;
  pair->data.is_col = false;
  zListInit( &pair->data.vlist );
  zListInsertHead( &cd->plist, pair );
  zPH3DInit( &pair->data.colvol );
  return pair;
}

/* (static)
 * _rkCDPairCellReg
 * - register all pairs of collision detection cells.
 */
rkCD *_rkCDPairCellReg(rkCD *cd, rkCDCell *cell)
{
  rkCDCell *cp;
  zVec3D v1, v2;

  zListForEach( &cd->clist, cp ){
    if( cp->data.type == RK_CD_CELL_STAT && cell->data.type == RK_CD_CELL_STAT )
      continue;
    if( cp->data.link == cell->data.link )
      continue;
    /* ignore a pair in permanent collision */
    if( cp->data.chain ==  cell->data.chain )
      if( zColChkAABox3D(&cp->data.aabb, &cell->data.aabb) )
        if( zColChkBox3D(&cp->data.obb, &cell->data.obb) )
          if( zColChkPH3D( &cp->data.ph, &cell->data.ph, &v1, &v2 ) )
            continue;
    if( !_rkCDPairReg( cd, cp, cell ) ) break;
  }
  return cd;
}

/* rkCDReset
 * - reset a collision detector.
 */
void rkCDReset(rkCD *cd)
{
  rkCDPair *pair;
  rkCDCell *cell;

  zListForEach( &cd->plist, pair ){
    pair->data.is_col = false;
    zPH3DDestroy( &pair->data.colvol );
  }
   zListForEach(&cd->clist, cell){
    cell->data._ph_update_flag = 0;
    cell->data._bb_update_flag = 0;
  }
}

/* rkCDPairReg
 * - register a pair of links in collision detector.
 */
rkCD *rkCDPairReg(rkCD *cd, rkLink *link1, rkLink *link2)
{
  rkCDCell *cp0, *cp1;
  rkCDPair *pair;

  zListForEach( &cd->plist, pair )
    if( ( pair->data.cell[0]->data.link == link1 &&
          pair->data.cell[1]->data.link == link2 ) ||
        ( pair->data.cell[0]->data.link == link2 &&
          pair->data.cell[1]->data.link == link1 ) )
      return cd;

  zListForEach( &cd->clist, cp0 )
    for( cp1=zListCellNext(cp0); cp1!=zListRoot(&cd->clist); cp1=zListCellNext(cp1) )
      if( ( cp0->data.link == link1 && cp1->data.link == link2 ) ||
          ( cp0->data.link == link2 && cp1->data.link == link1 ) )
        _rkCDPairReg( cd, cp0, cp1 );
  return cd;
}

/* rkCDPairUnreg
 * - unregister a pair of links from a collision detector.
 */
void rkCDPairUnreg(rkCD *cd, rkLink *link1, rkLink *link2)
{
  rkCDPair *cp, *temp;

  if( zListIsEmpty( &cd->plist ) ) return;
  zListForEach( &cd->plist, cp ){
    if( ( cp->data.cell[0]->data.link == link1 &&
          cp->data.cell[1]->data.link == link2 ) ||
        ( cp->data.cell[0]->data.link == link2 &&
          cp->data.cell[1]->data.link == link1 ) ){
      temp = zListCellPrev(cp);
      _rkCDPairDestroy(cp);
      zListPurge( &cd->plist, cp );
      zFree( cp );
      cp = temp;
    }
  }
}

void rkCDPairChainUnreg(rkCD *cd, rkChain *chain)
{
  rkCDPair *cp, *temp;

  if( zListIsEmpty( &cd->plist ) ) return;
  zListForEach( &cd->plist, cp ){
    if( ( cp->data.cell[0]->data.chain == chain &&
          cp->data.cell[1]->data.chain == chain ) ){
      temp = zListCellPrev(cp);
       _rkCDPairDestroy(cp);
      zListPurge( &cd->plist, cp );
      zFree( cp );
      cp = temp;
    }
  }
}

void rkCDPairWrite(rkCD *cd)
{
  rkCDPair *cp;
  register int i = 0;

  printf( "number : flag : chain1 shape1 : chain2 shape2\n" );
  zListForEach( &cd->plist, cp ){
    printf( "%d : %s : %s %s : %s %s\n", i,
      zBoolExpr(cp->data.is_col),
      zName(cp->data.cell[0]->data.chain), zName(cp->data.cell[0]->data.shape),
      zName(cp->data.cell[1]->data.chain), zName(cp->data.cell[1]->data.shape) );
    i++;
  }
}

#define __vec_label_write(label,v) \
  printf( "%s %g %g %g", label, zVec3DElem(v,zX), zVec3DElem(v,zY), zVec3DElem(v,zZ) )

void rkCDPairVertWrite(rkCD *cd)
{
  rkCDPair *cp;
  rkCDVert *v;
  register int i = 0;

  printf( "number : flag : chain1 shape1 : chain2 shape2\n" );
  zListForEach( &cd->plist, cp ){
    printf( "%d : %s : %s %s : %s %s\n", i,
      zBoolExpr(cp->data.is_col),
      zName(cp->data.cell[0]->data.chain), zName(cp->data.cell[0]->data.shape),
      zName(cp->data.cell[1]->data.chain), zName(cp->data.cell[1]->data.shape) );
    zListForEach( &cp->data.vlist, v ){
      printf( "%s %s",
        zName(v->data.cell->data.chain), zName(v->data.cell->data.link) );
      __vec_label_write( " vert", v->data.vert );
      __vec_label_write( " norm", &v->data.norm );
      __vec_label_write( " pro",  &v->data.pro );
      __vec_label_write( " ref",  &v->data.ref );
      printf( "\n" );
    }
    i++;
  }
}

/* rkCDChainReg
 * - register a chain to a collision detector.
 */
rkCD *rkCDChainReg(rkCD *cd, rkChain *chain, rkCDCellType type)
{
  zShapeListCell *sc;
  register int i;

  for( i=0; i<rkChainNum(chain); i++ )
    zListForEach( &rkChainLink(chain,i)->body.shapelist, sc ){
      if( !rkCDCellReg( &cd->clist, chain, rkChainLink(chain,i), sc->data, type ) )
        return NULL;
      if( zListNum(&cd->clist) > 1 )
        _rkCDPairCellReg( cd, zListHead(&cd->clist) );
    }
  return cd;
}

void rkCDChainUnreg(rkCD *cd, rkChain *chain)
{
  rkCDPair *pair, *ptemp;
  rkCDCell *cell, *ctemp;

  zListForEach( &cd->plist, pair ){
    if( ( pair->data.cell[0]->data.chain == chain ) ||
        ( pair->data.cell[1]->data.chain == chain ) ){
      ptemp = zListCellPrev(pair);
      _rkCDPairDestroy(pair);
      zListPurge( &cd->plist, pair );
      zFree( pair );
      pair = ptemp;
    }
  }
  zListForEach( &cd->clist, cell ){
    if( cell->data.chain == chain ){
      ctemp = zListCellPrev(cell);
      _rkCDCellDestroy(cell);
      zListPurge( &cd->clist, cell );
      zFree( cell );
      cell = ctemp;
    }
  }
}

void _rkCDColChkAABB(rkCD *cd)
{
  rkCDPair *cp;

  zListForEach( &cd->plist, cp ){
    /* TO BE MODIFIED */
    rkCDCellUpdateBB( cp->data.cell[0] );
    rkCDCellUpdateBB( cp->data.cell[1] );
    if( zColChkAABox3D( &cp->data.cell[0]->data.aabb, &cp->data.cell[1]->data.aabb ) )
      cp->data.is_col = true;
  }
}

void _rkCDColChkOBB(rkCD *cd)
{
  rkCDPair *cp;

  zListForEach( &cd->plist, cp )
    if( cp->data.is_col == true &&
        !zColChkBox3D( &cp->data.cell[0]->data.obb, &cp->data.cell[1]->data.obb ) )
      cp->data.is_col = false;
}

void _rkCDColChkGJK(rkCD *cd)
{
  rkCDPair *cp;
  zVec3D v1, v2;

  zListForEach( &cd->plist, cp ){
    if( cp->data.is_col == true ){
      /* TO BE MODIFIED */
      rkCDCellUpdatePH( cp->data.cell[0] );
      rkCDCellUpdatePH( cp->data.cell[1] );
      if( !zColChkPH3D( &cp->data.cell[0]->data.ph, &cp->data.cell[1]->data.ph,
          &v1, &v2 ) )
        cp->data.is_col = false;
    }
  }
}

void rkCDColChkAABB(rkCD *cd)
{
  if( zListIsEmpty( &cd->plist ) ) return;
  rkCDReset( cd );
  _rkCDColChkAABB( cd );
}

void rkCDColChkOBB(rkCD *cd)
{
  rkCDColChkAABB( cd );
  _rkCDColChkOBB( cd );
}

void rkCDColChkGJK(rkCD *cd)
{
  rkCDColChkOBB( cd );
  _rkCDColChkGJK( cd );
}

void rkCDColChkGJKOnly(rkCD *cd)
{
  rkCDPair *cp;
  zVec3D v1, v2;

  if( zListIsEmpty( &cd->plist ) ) return;
  rkCDReset( cd );
  zListForEach( &cd->plist, cp ){
    /* TO BE MODIFIED */
    rkCDCellUpdatePH( cp->data.cell[0] );
    rkCDCellUpdatePH( cp->data.cell[1] );
    if( zColChkPH3D( &cp->data.cell[0]->data.ph, &cp->data.cell[1]->data.ph,
                     &v1, &v2 ) )
      cp->data.is_col = true;
  }
}

void _rkCDVertNorm(zPH3D *ph, zVec3D *vert, zVec3D *norm, zVec3D *pro)
{
  zPH3DClosest( ph, vert, pro );
  zVec3DSub( pro, vert, norm );
  zVec3DNormalizeDRC( norm );
}

rkCDVert *_rkCDVertReg(rkCDPair *pair, rkCDVertList *vlist, rkCDCell *cell0, int v_id)
{
  rkCDVert *v, *cp;
  zVec3D pro, sub, vert;
  bool flag = false;
  rkCDCell *cell1;

  if( !( v = zAlloc( rkCDVert, 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  cell1 = pair->data.cell[ pair->data.cell[0] == cell0 ? 1 : 0 ];
  v->data.cell = cell0;
  v->data.vert = zPH3DVert(&cell0->data.ph,v_id);

  zListForEach( &pair->data.vlist, cp )
    if( cp->data.vert == zPH3DVert(&cell0->data.ph,v_id) ){
      zVec3DCopy( &cp->data._ref, &v->data._ref );
      zVec3DCopy( &cp->data._norm, &v->data._norm );
      zVec3DCopy( &cp->data._axis[0], &v->data._axis[0] );
      zVec3DCopy( &cp->data._axis[1], &v->data._axis[1] );
      zVec3DCopy( &cp->data._axis[2], &v->data._axis[2] );
      zXfer3DInv( rkLinkWldFrame(cell1->data.link), zPH3DVert(&cell0->data.ph, v_id), &vert );
      zVec3DSub( &cp->data._pro, &vert, &sub );
      zVec3DProj( &sub, &v->data._norm, &pro );
      zVec3DAdd( &vert, &pro, &v->data._pro );
      zMulMatVec3D( rkLinkWldAtt(cell1->data.link), &v->data._norm, &v->data.norm );
      zXfer3D( rkLinkWldFrame(cell1->data.link), &v->data._pro, &v->data.pro );
      zXfer3D( rkLinkWldFrame(cell1->data.link), &v->data._ref, &v->data.ref );
      zXfer3D( rkLinkWldFrame(cell1->data.link), &v->data._axis[0], &v->data.axis[0] );
      zXfer3D( rkLinkWldFrame(cell1->data.link), &v->data._axis[1], &v->data.axis[1] );
      zXfer3D( rkLinkWldFrame(cell1->data.link), &v->data._axis[2], &v->data.axis[2] );
      v->data.type = cp->data.type;
      flag = true;
    }
  if( zListIsEmpty( &pair->data.vlist ) || flag == false ){
    _rkCDVertNorm( &cell1->data.ph, zPH3DVert(&cell0->data.ph, v_id), &v->data.norm, &v->data.pro );
    zVec3DCopy( &v->data.pro, &v->data.ref );
    zVec3DCopy( &v->data.norm, &v->data.axis[0] );
    zVec3DOrthoSpace( &v->data.axis[0], &v->data.axis[1], &v->data.axis[2] );
    zMulMatTVec3D( rkLinkWldAtt(cell1->data.link), &v->data.norm, &v->data._norm );
    zXfer3DInv( rkLinkWldFrame(cell1->data.link), &v->data.pro, &v->data._pro );
    zXfer3DInv( rkLinkWldFrame(cell1->data.link), &v->data.ref, &v->data._ref );
    zXfer3DInv( rkLinkWldFrame(cell1->data.link), &v->data.axis[0], &v->data._axis[0] );
    zXfer3DInv( rkLinkWldFrame(cell1->data.link), &v->data.axis[1], &v->data._axis[1] );
    zXfer3DInv( rkLinkWldFrame(cell1->data.link), &v->data.axis[2], &v->data._axis[2] );
    v->data.type = RK_CONTACT_STICK;
  }
  zListInsertHead( vlist, v );
  return v;
}

void _rkCDColChkVert(rkCD *cd)
{
  rkCDPair *cp;
  rkCDVertList temp;
  register int i;

  cd->colnum = 0;
  zListForEach( &cd->plist, cp ){
    if( cp->data.is_col == true ){
      /* TO BE MODIFIED */
      rkCDCellUpdatePH( cp->data.cell[0] );
      rkCDCellUpdatePH( cp->data.cell[1] );
      zListInit( &temp );
      for( i=0; i<zPH3DVertNum(&cp->data.cell[0]->data.ph); i++ )
        if( zPH3DPointIsInside( &cp->data.cell[1]->data.ph, zPH3DVert(&cp->data.cell[0]->data.ph,i), false ) ){
          _rkCDVertReg( cp, &temp, cp->data.cell[0], i );
          cd->colnum++;
        }
      for( i=0; i<zPH3DVertNum(&cp->data.cell[1]->data.ph); i++ )
        if( zPH3DPointIsInside( &cp->data.cell[0]->data.ph, zPH3DVert(&cp->data.cell[1]->data.ph,i), false ) ){
          _rkCDVertReg( cp, &temp, cp->data.cell[1], i );
          cd->colnum++;
        }
      zListDestroy( rkCDVert, &cp->data.vlist );
      zListMove( &temp, &cp->data.vlist );
      if( zListIsEmpty( &cp->data.vlist ) )
        cp->data.is_col = false;
    } else{
      if( !zListIsEmpty( &cp->data.vlist ) )
        zListDestroy( rkCDVert, &cp->data.vlist );
    }
  }
}

void rkCDColChkVert(rkCD *cd)
{
  rkCDColChkOBB( cd );
  _rkCDColChkVert( cd );
}

void _rkCDColChkOBBVert(rkCD *cd)
{
  rkCDPair *cp;
  rkCDVertList temp;
  register int i;

  cd->colnum = 0;
  zListForEach( &cd->plist, cp ){
    if( cp->data.is_col == true ){
      /* TO BE MODIFIED */
      zBox3DToPH( &cp->data.cell[0]->data.obb, &cp->data.cell[0]->data.ph );
      zBox3DToPH( &cp->data.cell[1]->data.obb, &cp->data.cell[1]->data.ph );
      zListInit( &temp );
      for( i=0; i<zPH3DVertNum(&cp->data.cell[0]->data.ph); i++ )
        if( zPH3DPointIsInside( &cp->data.cell[1]->data.ph, zPH3DVert(&cp->data.cell[0]->data.ph,i), false ) ){
          _rkCDVertReg( cp, &temp, cp->data.cell[0], i );
          cd->colnum++;
        }
      for( i=0; i<zPH3DVertNum(&cp->data.cell[1]->data.ph); i++ )
        if( zPH3DPointIsInside( &cp->data.cell[0]->data.ph, zPH3DVert(&cp->data.cell[1]->data.ph,i), false ) ){
          _rkCDVertReg( cp, &temp, cp->data.cell[1], i );
          cd->colnum++;
        }
      zListDestroy( rkCDVert, &cp->data.vlist );
      zListMove( &temp, &cp->data.vlist );
      if( zListIsEmpty( &cp->data.vlist ) )
        cp->data.is_col = false;
    } else{
      if( !zListIsEmpty( &cp->data.vlist ) )
        zListDestroy( rkCDVert, &cp->data.vlist );
    }
  }
}

void rkCDColChkOBBVert(rkCD *cd)
{
  rkCDColChkOBB( cd );
  _rkCDColChkOBBVert( cd );
}

void rkCDColVol(rkCD *cd)
{
  rkCDPair *cp;

  zListForEach( &cd->plist, cp )
    if( cp->data.is_col == true )
    zIntersectPH3DFast( &cp->data.cell[0]->data.ph, &cp->data.cell[1]->data.ph, &cp->data.colvol );
}
