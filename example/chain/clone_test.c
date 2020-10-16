#include <roki/rk_chain.h>

int main(int argc, char *argv[])
{
  rkChain chain, chain_c, *ret;

  if( !rkChainReadFile( &chain, "../model/puma.zkc" ) )
    exit( 1 );
  ret = rkChainClone( &chain, &chain_c );
  rkChainDestroy( &chain );

  if( ret ){
    eprintf( "cloning succeeded.\n" );
    rkChainWriteFile( &chain_c, "clone" );
  }
  rkChainDestroy( &chain_c );
  return 0;
}
