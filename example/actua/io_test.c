#include <roki/rk_actua.h>

int main(void)
{
  rkMotor motor;

  rkMotorRead( &motor );
  rkMotorWrite( &motor );
  return 0;
}
