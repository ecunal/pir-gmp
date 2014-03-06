#define __STDC_FORMAT_MACROS
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <inttypes.h>
#include <queue>
#include <pthread.h>
#include <unistd.h>

using namespace std;

#define PAR 1

uint64_t rdtsc() {
	unsigned int lo, hi;
	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
	return ((uint64_t) hi << 32) | lo;
}

const int Q_SIZE = 512;

int shared_q[Q_SIZE];

void* produce(void* arg) {

	int d = (rand() % 256) + 256, sum = 0;

	for (int i = 0; i < (d - 1); i++) {
		int elmt = (rand() % 4) + 1;
		shared_q[i] = elmt;
		sum += elmt;
	}

	shared_q[d - 1] = 0;

}

void* consume(void* arg) {

	int i = 0, sum = 0;

	while (true) {

		int elmt = shared_q[i];

		if (elmt == 0) {
			break;
		} else if (elmt != -1) {
			sum += elmt;
			i++;
		}
	}

}

int main(int argc, char *argv[]) {

	if (argc != 2) {
		cout << "usage: " << argv[0]
				<< " <0 or 1> : 0 for serial, 1 for parallel" << endl;
		exit(0);
	}

	srand(time(0));
	int test_count = 10000, k_size = 10, warmup = 10;
	priority_queue<uint64_t> k_best;

	clock_t begin, end;
	double time_spent;

	begin = clock();

	for (int i = 0; i < test_count; i++) {

		for (int j = 0; j < Q_SIZE; j++) {
			shared_q[j] = -1;
		}

		/* Start timing */
		//uint64_t _start_ = rdtsc();


		if (atoi(argv[1])) {

			pthread_t t1, t2;

			pthread_create(&t1, NULL, produce, NULL);
			pthread_create(&t2, NULL, consume, NULL);

			pthread_join(t1, NULL);
			pthread_join(t2, NULL);
		} else {

			produce(NULL);
			consume(NULL);
		}



		//uint64_t _time_ = rdtsc() - _start_;

		/* End timing */

//		if (i >= warmup)
//			k_best.push(_time_);
//
//		while (k_best.size() > k_size) {
//			k_best.pop();
//		}

	}

	end = clock();

	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

	cout << time_spent << endl;


//	uint64_t average = 0;
//
//	while (!k_best.empty()) {
//
//		average += k_best.top();
//		printf("Best ones: %" PRIu64 "\n", k_best.top());
//		k_best.pop();
//	}
//
//	average = average / k_size;
//
//	printf("Average cycles: %" PRIu64 "\n", average);

	return 0;
}
