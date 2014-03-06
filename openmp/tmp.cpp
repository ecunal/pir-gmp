#include <stdio.h>
 int main()
 {
 	#pragma omp parallel for
 	for(int n=0; n<10; ++n)
 	{
   		printf(" %d", n);
 	}
 	printf(".\n");
 }
