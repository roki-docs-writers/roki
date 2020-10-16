#include <roki/rk_link.h>

void link_init(rkLink *l, int id, char *name)
{
  rkLinkInit( l );
  zNameSet( l, name );
}

int main(void)
{
  rkLink l[5];

  link_init( &l[0], 0, "root" );
  link_init( &l[1], 1, "b1" );
  link_init( &l[2], 2, "b2" );
  link_init( &l[3], 3, "b1-1" );
  link_init( &l[4], 4, "b2-1" );

  rkLinkAddChild( &l[0], &l[1] );
  rkLinkAddChild( &l[0], &l[2] );
  rkLinkAddChild( &l[1], &l[3] );
  rkLinkAddChild( &l[2], &l[4] );

  rkLinkConnectionWrite( &l[0], 0 );
  return 0;
}
