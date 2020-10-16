#include <roki/rk_g.h>
#include <roki/rk_actua.h>

int main(void)
{
  rkMotor motor;

  rkMotorCreate( &motor, rkN2Kgf(0.0525), 0.5, 100, 42.0,-42.0, "RE35" );
  printf( "viscosity = %g\n", rkMotorViscos( &motor ) );
  return 0;
}
