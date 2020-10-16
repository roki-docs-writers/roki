/* RoKi - Robot Kinetics library
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * rk_cd - collision detection
 * contributer: 2014-2015 Ken'ya Tanaka
 */

#ifndef __RK_CD_H__
#define __RK_CD_H__

#include <roki/rk_chain.h>

__BEGIN_DECLS

/*! \brief type to classify stationary or movable shape */
typedef enum{ RK_CD_CELL_STAT, RK_CD_CELL_MOVE } rkCDCellType;

/* ********************************************************** */
/*! \brief collision detection cell class.
 *//* ******************************************************* */
typedef struct{
  zShape3D *shape;   /*!< the shape to be checked */
  rkLink *link;      /*!< the link which the shape belongs to */
  rkChain *chain;    /*!< the chain which the shape belongs to */
  rkCDCellType type; /*!< stationary or movable shape */
  zAABox3D aabb;     /*!< axis aligned bounding box in the world frame */
  zBox3D obb;        /*!< oriented bounding box in the world frame */
  zPH3D ph;          /*!< polyhedron in the world frame */
  /*! \cond */
  bool _bb_update_flag; /* check if boundary boxes are updated */
  bool _ph_update_flag; /* check if polyhedra are updated */
  /*! \endcond */
} rkCDCellDat;
zListClass( rkCDCellList, rkCDCell, rkCDCellDat );

/*! \brief register collision detection cell
 */
__EXPORT rkCDCell *rkCDCellReg(rkCDCellList *clist, rkChain *chain, rkLink *link, zShape3D *shape, rkCDCellType type);

/*! \brief update the bounding box of a collision detection cell
 */
__EXPORT void rkCDCellUpdateBB(rkCDCell *cell);

/*! \brief update the polyhedron of a collision detection cell
 */
__EXPORT void rkCDCellUpdatePH(rkCDCell *cell);

/*! \brief update a collision detection cell
 */
__EXPORT void rkCDCellUpdate(rkCDCell *cell);

/* ********************************************************** */
/*! \brief collision detection vertex class
 *//* ******************************************************* */
typedef struct{
  rkCDCell *cell;  /*!< collision detection cell */
  zVec3D *vert;    /*!< the vertex to be checked */
  zVec3D norm;     /*!< normal vector at the vertex */
  zVec3D axis[3];  /*!< z-x-y */
  zVec3D pro;      /*!< projection point of the vertex */
  zVec3D ref;      /*!< referential point to measure displacement of the vertex */
  /*! \cond */
  zVec3D _norm;    /* normal vector at the vertex in the link frame */
  zVec3D _axis[3]; /* axes of contact frame */
  zVec3D _pro;     /* projection point of the vertex in the link frame */
  zVec3D _ref;     /* referential point to measure displacement of the vertex in the link frame */
  zVec3D f;        /* contact force */
  rkContactType type; /* type to classify stick/slip mode */
  /*! \endcond */
} rkCDVertDat;
zListClass( rkCDVertList, rkCDVert, rkCDVertDat );

/* ********************************************************** */
/*! \brief collision detection pair class.
 *//* ******************************************************* */
typedef struct{
  rkCDCell *cell[2];  /*!< a pair of the collision detection cell */
  bool is_col;        /*!< flag to check collision */
  rkCDVertList vlist; /*!< a list of vertices */
  zVec3D norm;        /*!< normal vector of the pair */
  zPH3D colvol;       /*!< collision volume */
} rkCDPairDat;
zListClass( rkCDPairList, rkCDPair, rkCDPairDat );

/* ********************************************************** */
/*! \brief collision detection class.
 *//* ******************************************************* */
typedef struct{
  rkCDCellList clist; /*!< a list of collision detection cells */
  rkCDPairList plist; /*!< a list of collision detection pairs */
  int colnum;         /*!< the number of collision pairs */
} rkCD;

__EXPORT rkCD *rkCDCreate(rkCD *cd);
__EXPORT void rkCDDestroy(rkCD *cd);
__EXPORT void rkCDReset(rkCD *cd);

__EXPORT rkCD *rkCDChainReg(rkCD *cd, rkChain *chain, rkCDCellType type);
__EXPORT void rkCDChainUnreg(rkCD *cd, rkChain *chain);

__EXPORT rkCD *rkCDPairReg(rkCD *cd, rkLink *link1, rkLink *link2);
__EXPORT void rkCDPairUnreg(rkCD *cd, rkLink *link1, rkLink *link2);
__EXPORT void rkCDPairChainUnreg(rkCD *cd, rkChain *chain);

__EXPORT void rkCDPairWrite(rkCD *cd);
__EXPORT void rkCDPairVertWrite(rkCD *cd);

__EXPORT void rkCDColChkAABB(rkCD *cd);    /* AABB */
__EXPORT void rkCDColChkOBB(rkCD *cd);     /* AABB->OBB */
__EXPORT void rkCDColChkGJK(rkCD *cd);     /* AABB->OBB->GJK */
__EXPORT void rkCDColChkVert(rkCD *cd);    /* AABB->OBB->Vert(PH) */
__EXPORT void rkCDColChkOBBVert(rkCD *cd); /* AABB->OBB->Vert(OBB) */

__EXPORT void rkCDColChkGJKOnly(rkCD *cd); /* GJK */

__EXPORT void rkCDColVol(rkCD *cd);

__END_DECLS

#endif /* __RK_CD_H__ */
