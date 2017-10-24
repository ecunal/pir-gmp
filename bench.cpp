#include "benchmark/benchmark.h"

#include <unordered_set>

static void BM_UnorderedSetSequentialInsert( benchmark::State &state )
{
  while ( state.KeepRunning() )
  {
    std::unordered_set<int> insertion_test;
    for ( int i = 0, i_end = state.range( 0 ); i < i_end; i++ )
    {
      insertion_test.insert( i );
    }
  }
}
BENCHMARK( BM_UnorderedSetSequentialInsert )->Unit(benchmark::kMillisecond)->Range( 1 << 16, 1 << 24 );

BENCHMARK_MAIN();
