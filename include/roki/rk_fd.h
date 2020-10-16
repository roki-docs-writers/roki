/* RoKi - Robot Kinetics library
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * rk_fd - forward dynamics computation
 * contributer: 2014-2015 Naoki Wakisaka
 */

#ifndef __RK_FD_H__
#define __RK_FD_H__

#include <roki/rk_chain.h>
#include <roki/rk_abi.h>
#include <roki/rk_cd.h>

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: rkFD
   forward dynamics class
 * ********************************************************** */
typedef struct{
  rkChain chain;   /* chain */
  int _offset;     /* offset */
  zVecStruct _dis; /* displacement */
  zVecStruct _vel; /* velocity */
  zVecStruct _acc; /* acceleration */
} rkFDCellDat;

zListClass( rkFDCellList, rkFDCell, rkFDCellDat );

typedef struct _rkFD{
  rkFDCellList list;    /* chain list */
  rkContactInfoPool ci; /* contact information */
  rkCD cd;              /* contact maneger */
  zVec dis, vel;        /* total joint state */
  zVec acc;
  int size;             /* total joint size */
  double t, dt;

  /* working space */
  zODE2 _ode;    /* ode solver */
  int _ode_step;
  double _comp_k;
  double _comp_l;
  double _kf_weight;
  int _pyramid;
  int colnum_e; /* elastioc collision */
  int colnum_r; /* rigid collision */

  /* solver */
  void (*_solve_contact)(struct _rkFD*,bool);
  rkContactInfo _ci_def;
} rkFD;

#define rkFDTime(f) (f)->t
#define rkFDDT(f)   (f)->dt

__EXPORT rkFD *rkFDCreate(rkFD *fd);
__EXPORT void rkFDDestroy(rkFD *fd);

__EXPORT rkFDCell *rkFDChainReg(rkFD *fd, rkChain *chain);
__EXPORT rkFDCell *rkFDChainRegFile(rkFD *fd, char filename[]);
__EXPORT bool rkFDChainUnreg(rkFD *fd, rkFDCell *cell);

__EXPORT void rkFDFK(rkFD *fd, zVec dis);
__EXPORT void rkFDUpdateRate(rkFD *fd, zVec vel, zVec acc);
__EXPORT void rkFDUpdateFKRate(rkFD *fd);

__EXPORT void rkFDChainSetDis(rkFDCell *lc, zVec dis);
__EXPORT void rkFDChainSetVel(rkFDCell *lc, zVec vel);
__EXPORT bool rkFDContactInfoReadFile(rkFD *fd, char filename[]);

/* for a fake-crawler */
__EXPORT void rkFDCDCellSetSlideMode(rkCDCell *cell, bool mode);
__EXPORT void rkFDCDCellSetSlideVel(rkCDCell *cell, double vel);
__EXPORT void rkFDCDCellSetSlideAxis(rkCDCell *cell, zVec3D *axis);
__EXPORT rkCDCell *rkFDShape3DGetCDCell(rkFD *fd, zShape3D *shape);
__EXPORT rkCDCell *rkFDShape3DSetSlideMode(rkFD *fd, zShape3D *shape, bool mode);
__EXPORT rkCDCell *rkFDShape3DSetSlideVel(rkFD *fd, zShape3D *shape, double vel);
__EXPORT rkCDCell *rkFDShape3DSetSlideAxis(rkFD *fd, zShape3D *shape, zVec3D *axis);

/* utility */
#define rkFDCDPairIsContactInfoType(cp,type) ((cp)->data.is_col && rkContactInfoType((cp)->data.ci) == (type))

/* ODE updater */
__EXPORT zVec rkFDODECatDefault(zVec x, double k, zVec v, zVec xnew, void *util);
__EXPORT zVec rkFDODESubDefault(zVec x1, zVec x2, zVec dx, void *util);
#define rkFDODE2Assign(f,t)        zODE2Assign( &(f)->_ode, t, rkFDODECatDefault, NULL, rkFDODESubDefault, NULL )
#define rkFDODE2AssignRegular(f,t) zODE2AssignRegular( &(f)->_ode, t )

/* setting */
#define rkFDSetDT(f,t)  ( (f)->dt = (t) )
#define rkFDSetJointFricPrp(f,k,l) do{\
  (f)->_comp_k = (k);\
  (f)->_comp_l = (l);\
} while(0)
#define rkFDSetKFWeight(f,w) ( (f)->_kf_weight = (w) )
#define rkFDSetPyramidOrder(f,p) ( (f)->_pyramid = (p) )

#define rkFDSetSolver(f,type) do{\
  (f)->_solve_contact = rkFDSolveContact_##type;\
  rkFDSetContactInfoDefault_##type(f);\
} while(0)

__EXPORT rkFD *rkFDSolve(rkFD *fd);
__EXPORT void rkFDUpdateInit(rkFD *fd);
__EXPORT rkFD *rkFDUpdate(rkFD *fd);
__EXPORT void rkFDUpdateDestroy(rkFD *fd);
#define rkFDSolveContact(f,b) (f)->_solve_contact(f,b)

/* solver */
__EXPORT void rkFDSetContactInfoDefault_Vert(rkFD *fd);
__EXPORT void rkFDSolveContact_Vert(rkFD *fd, bool doUpRef);

__EXPORT void rkFDSetContactInfoDefault_Volume(rkFD *fd);
__EXPORT void rkFDSolveContact_Volume(rkFD *fd, bool doUpRef);

/* for debug */
__EXPORT void rkFDWrite(rkFD *fd);

__END_DECLS

#endif /* __RK_FD_H__ */
