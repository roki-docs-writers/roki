#include <roki/rk_contact.h>

int main(void)
{
  rkContactInfo ci;

  rkContactInfoElasticCreate( &ci, 5000, 100, 0.5, "b1", "b2" );
  rkContactInfoWrite( &ci );
  return 0;
}
