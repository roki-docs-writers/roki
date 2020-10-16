#include <roki/rk_contact.h>

void assoc(rkContactInfoPool *ci, char *key1, char *key2)
{
  rkContactInfo *cp;

  cp = rkContactInfoPoolAssoc( ci, key1, key2 );
  printf( "associated from <%s,%s>\n", key1, key2 );
  cp ? rkContactInfoWrite( cp ) : printf( "failed.\n" );
}

int main(void)
{
  rkContactInfoPool ci;

  rkContactInfoPoolReadFile( &ci, "test" );
  assoc( &ci, "stf1", "stf2" );
  assoc( &ci, "stf4", "stf3" );
  assoc( &ci, "stf1", "hoge" );
  rkContactInfoPoolDestroy( &ci );
  return 0;
}
