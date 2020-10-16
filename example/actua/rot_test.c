#include <roki/rk_g.h>
#include <roki/rk_actua.h>

#define DT 0.001
#define T 10.0
#define I  1.0
#define V 12.0

int main(void)
{
  rkMotor motor;
  double angle = 0, angleveloc = 0;
  double t;

  rkMotorCreate( &motor, rkN2Kgf(0.0525), 0.5, 100, 42.0,-42.0, "RE35" );
  for( t=0; t<T; t+=DT ){
    angle += angleveloc * DT;
    angleveloc += ( rkMotorTorque( &motor, V, angleveloc ) / I ) * DT;
    printf( "%.10f %.10f\n", angle, angleveloc );
  }
  return 0;
}
