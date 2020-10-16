#include <roki/rk_contact.h>

int main(void)
{
  rkContactInfo ci;

  rkContactInfoRead( &ci );
  rkContactInfoWrite( &ci );
  return 0;
}
