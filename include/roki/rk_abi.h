/* RoKi - Robot Kinetics library
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * rk_abi - Articulated Body Inertia Method.
 * contributer: 2014-2015 Naoki Wakisaka
 */

#ifndef __RK_ABI_H__
#define __RK_ABI_H__

#include <roki/rk_chain.h>

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: rkABI
 * ********************************************************** */

__EXPORT void rkLinkABIInit(rkLink *link);
__EXPORT void rkChainABIInit(rkChain *chain);

__EXPORT void rkLinkABIDestroy(rkLink *link);
__EXPORT void rkChainABIDestroy(rkChain *chain);

__EXPORT void rkLinkABIUpdateInit(rkLink *link, zVec6D *pvel);
__EXPORT void rkLinkABIUpdateBackward(rkLink *link);
__EXPORT void rkLinkABIUpdateForward(rkLink *link, zVec6D *pa);

#define rkChainABIUpdateInit(c)     rkLinkABIUpdateInit( rkChainRoot(c), ZVEC6DZERO )
#define rkChainABIUpdateBackward(c) rkLinkABIUpdateBackward( rkChainRoot(c) )
#define rkChainABIUpdateForward(c)  rkLinkABIUpdateForward( rkChainRoot(c), ZVEC6DZERO )

__EXPORT void rkChainABIUpdate(rkChain *chain);
__EXPORT zVec rkChainABI(rkChain *chain, zVec dis, zVec vel, zVec acc);

__END_DECLS

#endif /* __RK_ABI_H__ */
