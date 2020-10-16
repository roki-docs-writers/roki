#include <roki/rk_force.h>

int main(void)
{
  rkWrenchList wl;
  rkWrench w[4], *wp;
  zVec6D w6, rw;

  zListInit( &wl );
  zVec6DCreate( rkWrenchW(&w[0]),-1, 0, 0, 0, 0, 0 );
  zVec3DCreate( rkWrenchPos(&w[0]), 1, 1, 0 );
  zVec6DCreate( rkWrenchW(&w[1]), 0,-1, 0, 0, 0, 0 );
  zVec3DCreate( rkWrenchPos(&w[1]),-1, 1, 0 );
  zVec6DCreate( rkWrenchW(&w[2]), 1, 0, 0, 0, 0, 0 );
  zVec3DCreate( rkWrenchPos(&w[2]),-1,-1, 0 );
  zVec6DCreate( rkWrenchW(&w[3]), 0, 1, 0, 0, 0, 0 );
  zVec3DCreate( rkWrenchPos(&w[3]), 1,-1, 0 );

  zListInsertTail( &wl, &w[0] );
  zListInsertTail( &wl, &w[1] );
  zListInsertTail( &wl, &w[2] );
  zListInsertTail( &wl, &w[3] );

  rkWrenchListNet( &wl, &rw );
  printf( "<net wrench>\n" );
  zVec6DWrite( &rw );
  zEndl();

  zVec6DClear( &rw );
  while( !zListIsEmpty( &wl ) ){
    zListDeleteTail( &wl, &wp );
    printf( "<force1 :%p>\n", wp );
    rkWrenchWrite( wp );
    rkWrenchXfer( wp, &w6 );
    printf( "(equivalent wrench)\n" );
    zVec6DWrite( &w6 );
    zVec6DAddDRC( &rw, &w6 );
  }
  zEndl();
  zVec6DWrite( &rw );
  return 0;
}
