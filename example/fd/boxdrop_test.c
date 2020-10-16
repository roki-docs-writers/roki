#include <roki/rk_fd.h>

#define T   2.0
#define DT  0.001
#define DTC 0.002

#define NMAX  9
#define N     5

int main(int argc, char *argv[])
{
  rkFD fd;
  int i, n;
  rkFDCell *cell[NMAX];
  zVec dis[NMAX];
  FILE *fp[NMAX];
  char name[BUFSIZ];

  zRandInit();
  rkFDCreate( &fd );
  rkFDContactInfoReadFile( &fd, "../model/cinfo.zci" );

  n = argc > 1 ? atoi( argv[1] ) : N;
  if( n > NMAX ) n = NMAX;
  for( i=0; i<n; i++ ){
    sprintf( name, "%d.zvs", i+1 );
    fp[i] = fopen( name, "w" );
    cell[i] = rkFDChainRegFile( &fd, "../model/box.zkc" );
    dis[i] = zVecAlloc( rkChainJointSize(&cell[i]->data.chain) );
    zVecElem(dis[i],0) = zRandF( -0.1, 0.1 );
    zVecElem(dis[i],1) = zRandF( -0.1, 0.1 );
    zVecElem(dis[i],2) = 0.1 + i*0.15;
    zVecElem(dis[i],3) = zDeg2Rad(zRandF(-90.0, 90.0));
    zVecElem(dis[i],4) = zDeg2Rad(zRandF(-90.0, 90.0));
    zVecElem(dis[i],5) = zDeg2Rad(zRandF(-90.0, 90.0));
    rkFDChainSetDis( cell[i], dis[i] );
    rkCDPairChainUnreg( &fd.cd, &cell[i]->data.chain );
  }
  rkFDChainRegFile( &fd, "../model/floor.zkc" );

  /* ode */
  rkFDODE2Assign( &fd, Regular );
  rkFDODE2AssignRegular( &fd, RKG );
  rkFDSetDT( &fd, DT );
  rkFDUpdateInit( &fd );

  while( rkFDTime(&fd) < T ){
    eprintf( "t = %f\n", rkFDTime(&fd) );
    rkFDUpdate( &fd );
    for( i=0; i<n; i++ ){
      rkChainGetJointDisAll( &cell[i]->data.chain, dis[i] );
      fprintf( fp[i], "%f ", rkFDDT(&fd) );
      zVecFWrite( fp[i], dis[i] );
    }
  }
  rkFDUpdateDestroy( &fd );

  for( i=0; i<n; i++ ){
    zVecFree( dis[i] );
    fclose( fp[i] );
  }
  rkFDDestroy( &fd );
  return 0;
}
