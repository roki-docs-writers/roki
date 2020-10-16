#include <roki/rk_contact.h>

int main(void)
{
  rkContactInfoPool ci;

  rkContactInfoPoolReadFile( &ci, "test" );
  rkContactInfoPoolWrite( &ci );
  rkContactInfoPoolDestroy( &ci );
  return 0;
}
