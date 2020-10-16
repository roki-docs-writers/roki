#include <roki/rk_link.h>

void link_init(rkLink *l, int id, char *name, byte jtype)
{
  rkLinkInit( l );
  zNameSet( l, name );
  rkJointCreate( rkLinkJoint(l), jtype );
}

void link_write(rkLink *l)
{
  rkLinkPostureWrite( l );
  if( rkLinkChild(l) ) link_write( rkLinkChild(l) );
  if( rkLinkSibl(l) ) link_write( rkLinkSibl(l) );
}

int main(void)
{
  rkLink l[3];
  double a1, a2;

  link_init( &l[0], 0, "root", RK_JOINT_FIXED );
  link_init( &l[1], 1, "l1", RK_JOINT_REVOL );
  link_init( &l[2], 2, "l2", RK_JOINT_REVOL );
  zVec3DCreate( rkLinkWldPos(&l[0]),-1, 0, 1 );
  zVec3DCreate( rkLinkOrgPos(&l[1]), 0, 1, 0 );
  zVec3DCreate( rkLinkOrgPos(&l[2]), 0, 1, 0 );
  rkLinkAddChild( &l[0], &l[1] );
  rkLinkAddChild( &l[1], &l[2] );
  zFrame3DCopy( rkLinkOrgFrame(&l[1]), rkLinkAdjFrame(&l[1]) );
  zFrame3DCopy( rkLinkOrgFrame(&l[2]), rkLinkAdjFrame(&l[2]) );
  rkLinkConnectionWrite( &l[0], 0 );

  do{
    printf( "theta1=?" );
    a1 = zDeg2Rad( zFDouble( stdin ) );
    printf( "theta2=?" );
    a2 = zDeg2Rad( zFDouble( stdin ) );

    rkLinkSetJointDis( &l[1], &a1 );
    rkLinkSetJointDis( &l[2], &a2 );
    rkLinkUpdateFrame( &l[0], Z_IDENTFRAME3D );
    link_write( &l[0] );
  } while( a1!=0 || a2!=0 );

  rkLinkDestroy( &l[0] );
  rkLinkDestroy( &l[1] );
  rkLinkDestroy( &l[2] );
  return 0;
}
