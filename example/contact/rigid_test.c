#include <roki/rk_contact.h>

int main(void)
{
  rkContactInfo ci;

  rkContactInfoRigidCreate( &ci, 0, 0.8, 0.5, "b1", "b2" );
  rkContactInfoWrite( &ci );
  return 0;
}
