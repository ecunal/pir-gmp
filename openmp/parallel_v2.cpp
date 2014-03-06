#define __STDC_FORMAT_MACROS
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
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

vector<int> shared_q;
int Q_SIZE;

void produce() {

	int sum = 0;

	for (int i = 0; i < (Q_SIZE - 1); i++) {
		int elmt = (rand() % 4) + 1;
		shared_q[i] = elmt;
		sum += elmt;
	}

	shared_q[Q_SIZE - 1] = 0;

}

void consume() {

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

	srand(time(0));

	ofstream myfile;
	myfile.open("results.txt", ios::out | ios::app);

	for (int x = 512; x < 20500; x += 512) {

		Q_SIZE = x;
		shared_q = vector<int>(Q_SIZE, -1);

		int test_count = 10000, k_size = 100, warmup = 1000;
		priority_queue<uint64_t, vector<uint64_t>, greater<uint64_t> > k_best;

		clock_t begin, end;
		double time_spent;

		begin = clock();

		for (int i = 0; i < test_count; i++) {

			for (int j = 0; j < Q_SIZE; j++) {
				shared_q[j] = -1;
			}

			/* Start timing */
			uint64_t _start_ = rdtsc();

			#pragma omp parallel for
			for (int k = 0; k < 2; k++) {

				if (k == 0) {
					produce();
				} else if (k == 1) {
					consume();
				}
			}

			uint64_t _time_ = rdtsc() - _start_;

			/* End timing */

			if (i >= warmup)
				k_best.push(_time_);

		}

		end = clock();

		time_spent = (double) (end - begin) / CLOCKS_PER_SEC;

		//cout << time_spent << endl;

		uint64_t average = 0;

		int c = 0;

		while (c < k_size) {

			average += k_best.top();
			//printf("Best ones: %" PRIu64 "\n", k_best.top());
			k_best.pop();
			c++;
		}

		average = average / k_size;

		//printf("Average cycles: %" PRIu64 "\n", average);

		myfile << Q_SIZE << "\t" << average << "\n";

	}

	myfile.close();
	return 0;
}
