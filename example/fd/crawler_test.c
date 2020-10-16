#include <roki/rk_fd.h>
#include <roki/rk_cd.h>

#define T 3.0
#define DT 0.001

#define CRAW_VEL 0.3

int main(int argc, char *argv[])
{
  rkFD fd;
  rkFDCell *fd_cell;
  rkCDCell *cd_cell[2];
  zVec dis;
  FILE *fp;
  register int i;

  fp = fopen( "crawler.zvs", "w" );

  rkFDCreate( &fd );
  rkFDContactInfoReadFile( &fd, "../model/cinfo.zci" );
  fd_cell = rkFDChainRegFile( &fd, "../model/crawler.zkc" );
  rkFDChainRegFile( &fd, "../model/floor2.zkc" );

  /* init */
  dis = zVecAlloc(rkChainJointSize(&fd_cell->data.chain));
  zVecElem(dis,1) = 0.5;
  zVecElem(dis,2) = 0.01;
  zVecElem(dis,5) = zDeg2Rad(-90);
  rkFDChainSetDis( fd_cell, dis );

  /* slide mode */
  for( i=0; i<2; i++ ){
    cd_cell[i] = rkFDShape3DGetCDCell( &fd, zMShape3DShape( rkChainShape( &fd_cell->data.chain ), i+1 ) );
    rkFDCDCellSetSlideMode( cd_cell[i], true );
    rkFDCDCellSetSlideVel( cd_cell[i], CRAW_VEL );
    rkFDCDCellSetSlideAxis( cd_cell[i], ZVEC3DY );
  }

  /* no self collision */
  rkCDPairChainUnreg( &fd.cd, &fd_cell->data.chain );

  /* ode */
  rkFDODE2Assign( &fd, Regular );
  rkFDODE2AssignRegular( &fd, RKG );
  rkFDSetDT( &fd, DT );
  rkFDUpdateInit( &fd );

  while( rkFDTime(&fd) < T ){
    eprintf( "t = %f\n", rkFDTime(&fd) );
    rkFDUpdate( &fd );
    rkChainGetJointDisAll( &fd_cell->data.chain, dis );
    fprintf( fp, "%f ", rkFDDT(&fd) );
    zVecFWrite( fp, dis );
  }
  rkFDUpdateDestroy( &fd );

  zVecFree( dis );
  rkFDDestroy( &fd );
  fclose( fp );
  return 0;
}
