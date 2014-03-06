#define __STDC_FORMAT_MACROS
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <inttypes.h>
#include <queue>

using namespace std;

queue<int> shared_q;

void clear(queue<int> &q) {
	queue<int> empty;
	swap(q, empty);
}

uint64_t rdtsc() {
	unsigned int lo, hi;
	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
	return ((uint64_t) hi << 32) | lo;
}

void producer() {

	int d = (rand() % 256) + 256;

	for (int i = 0; i < d; i++) {
		int elmt = (rand() % 4) + 1;
		shared_q.push(elmt);
	}

	shared_q.push(0);
}

void consumer() {

	while (true) {

			if (!shared_q.empty()) {

				int front = shared_q.front();
				shared_q.pop();

				if (front == 0) {
					break;
				}
			}
		}
}

int main() {

	int test_count = 100000, k_size = 1000, warmup = 10000;

	priority_queue<uint64_t> k_best;

	srand(time(0));

	for (int i = 0; i < test_count; i++) {

		clear(shared_q);

		uint64_t _start_ = rdtsc();

		producer();
		consumer();

		uint64_t _time_ = rdtsc() - _start_;

//		printf("Total cycles: %" PRIu64 "\n", _time_);

		if(i >= warmup)
			k_best.push(_time_);

		while(k_best.size() > k_size) {
			k_best.pop();
		}

	}

	uint64_t average = 0;

	while( !k_best.empty() ) {

		average += k_best.top();
//		printf("Best ones: %" PRIu64 "\n", k_best.top());
		k_best.pop();
	}

	average = average / k_size;

	printf("Average cycles: %" PRIu64 "\n", average);

	return 0;
}
