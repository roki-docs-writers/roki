/* RoKi - Robot Kinetics library
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * rk_motor - motor model
 * contributer: 2014-2015 Naoki Wakisaka
 */

#ifndef __RK_MOTOR_H__
#define __RK_MOTOR_H__

#include <zm/zm.h>
#include <roki/rk_errmsg.h>

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: rkMotor
 * geared DC motor model class
 * ********************************************************** */

typedef struct{
  char *typestr; /*!< \brief a string for type identification */
  byte size; /*!< \brief size of motor input signals */
  void (* _init)(void*);
  void *(* _alloc)(void);
  void (* _copy)(void*,void*);

  void (* _setinput)(void*,double*); /*!< \brief set feasible motor input */
  void (* _inertia)(void*,double*); /*!< \brief rotor inertia of motor + gear */
  void (* _inputtrq)(void*,double*); /*!< \brief actuation torque from motor input signal */
  void (* _regist)(void*,double*,double*,double*); /*!< \brief registance torque (e.g. counter electromotive) */
  void (* _dtrq)(void*,double*,double*,double*,double*); /*!< \brief driving torque */

  bool (* _query)(FILE*,char*,void*);
  void (* _fprint)(FILE*,void*);
} rkMotorCom;

typedef struct{
  Z_NAMED_CLASS
  void *prp;
  rkMotorCom *com;
} rkMotor;

#define rkMotorTypeStr(m) (m)->com->typestr
#define rkMotorSize(m)    (m)->com->size

#define rkMotorInit(m) do{\
  zNameSetPtr( m, NULL );\
  (m)->prp = NULL;\
  (m)->com = NULL;\
} while(0)
__EXPORT rkMotor *rkMotorAssign(rkMotor *m, rkMotorCom *com);
__EXPORT rkMotor *rkMotorQueryAssign(rkMotor *m, char *str);
__EXPORT void rkMotorDestroy(rkMotor *m);

__EXPORT rkMotor *rkMotorClone(rkMotor *org, rkMotor *cln);

/* ********************************************************** */
/* method
 * ********************************************************** */

#define rkMotorSetInput(m,i)           (m)->com->_setinput( (m)->prp, (i) )
#define rkMotorInertia(m,i)            (m)->com->_inertia( (m)->prp, (i) )
#define rkMotorInputTrq(m,i)           (m)->com->_inputtrq( (m)->prp, (i) )
#define rkMotorRegistance(m,d,v,r)     (m)->com->_regist( (m)->prp, (d), (v), (r) )
#define rkMotorDrivingTrq(m,d,v,a,t)   (m)->com->_dtrq( (m)->prp, (d), (v), (a), (t) )

#define rkMotorQueryFScan(f,k,m)       (m)->com->_query( f, k, (m)->prp )
__EXPORT void rkMotorFPrint(FILE *fp, rkMotor *m);
#define rkMotorPrint(m)                rkMotorFPrint( stdout, m )

/* ********************************************************** */
/* CLASS: rkMotorArray
 * ********************************************************** */

#define ZTK_TAG_RKMOTOR "motor"
zArrayClass( rkMotorArray, rkMotor );

__EXPORT rkMotorArray *rkMotorArrayClone(rkMotorArray *org);
__EXPORT rkMotorArray *rkMotorArrayFScan(FILE *fp, rkMotorArray *m);
__EXPORT void rkMotorArrayFPrint(FILE *fp, rkMotorArray *m);

__END_DECLS

#include <roki/rk_motor_none.h> /* for unactuated joints */
#include <roki/rk_motor_trq.h>  /* direct torque motor */
#include <roki/rk_motor_dc.h>   /* DC motor */

#endif /* __RK_MOTOR_H__ */
