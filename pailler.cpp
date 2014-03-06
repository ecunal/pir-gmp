#include <iostream>
#include <gmp.h>
#include <time.h>
#include <fcntl.h> 
#include <unistd.h>

#define PREC 1E9
#define NUM_TESTS 1000

using namespace std;

int main() {

	double gm_time = 0, rn_time = 0;
	timespec begin, end;

	for (int t = 0; t < NUM_TESTS; t++) {

		// ********** KEYGEN PART ************* //

		// Define & initialize variables

		mpz_t p, q, n, g, n_sq, n_sqr, r, m, result1, result2;
		mpz_inits(p, q, n, g, n_sq, n_sqr, r, m, result1, result2, NULL);

		// Random device for seeding random state

		int random_dev = open("/dev/urandom", O_RDONLY);
		long random_seed;
		size_t random_len = 0;

		while (random_len < sizeof random_seed) {

			ssize_t result = read(random_dev,
					((char*) &random_seed) + random_len,
					(sizeof random_seed) - random_len);
			if (result < 0) {
				cout << "unable to read from urandom" << endl;
			}
			random_len += result;
		}

		close(random_dev);

		gmp_randstate_t state;
		gmp_randinit_mt(state);
		gmp_randseed_ui(state, random_seed);

		// 1. generate p & q of 512 bits

		mpz_urandomb(p, state, 512);
		mpz_nextprime(p, p);

		mpz_urandomb(q, state, 512);
		mpz_nextprime(q, q);

		// 2. n = p*q

		mpz_mul(n, p, q);

		// 3. g = n + 1

		mpz_add_ui(g, n, 1);

		// 4. n_sq = n^3

		mpz_pow_ui(n_sq, n, 2);
		mpz_pow_ui(n_sqr, n, 1);

		// ******* ENCRYPTION PART *********** //

		// 0. generate random r and m of 1024 bits

		mpz_urandomb(r, state, 1024);
		mpz_urandomb(m, state, 1024);

		// 1. g modpow m, n_sq  <--- measured

		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &begin);
		mpz_powm(result1, g, m, n_sq);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);

		double time_spent = (double) ((end.tv_sec - begin.tv_sec)
				+ ((end.tv_nsec - begin.tv_nsec) / PREC));

		gm_time += time_spent;

		// 2. r modpow n, n_sq  <--- measured

		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &begin);
		mpz_powm(result2, r, n_sqr, n_sq);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);

		time_spent = (double) ((end.tv_sec - begin.tv_sec)
				+ (end.tv_nsec - begin.tv_nsec) / PREC);

		rn_time += time_spent;

	}

	cout << "gm: " << (gm_time/NUM_TESTS)*1000 << " ms" << endl;
	cout << "rn: " << (rn_time/NUM_TESTS)*1000 << " ms" << endl;

	return 0;
}
