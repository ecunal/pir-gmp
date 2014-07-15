#include "DamgardJurik.h"
#define BIT_LENGTH 1024
#define TEST_CASES 100

using namespace std;

int main() {

	mpz_t cipher;
	mpz_init(cipher);

	for(int s=1; s<14; s++) {

		double gm_time = 0, rns_time = 0;

		DamgardJurik dj(BIT_LENGTH, s);

		if(s == 1) {
			dj.get_random(cipher, BIT_LENGTH);
		}

		if(s > 10) {
			for(int i=0; i<TEST_CASES; i++) {

				mpz_t r, random;
				mpz_inits(r, random, NULL);
				dj.get_random(random, BIT_LENGTH);
			
				double start_time = omp_get_wtime();
				dj.encrypt_exp_1(r, cipher);
				double end_time = omp_get_wtime();

				gm_time += end_time - start_time;
			
				start_time = omp_get_wtime();
				dj.encrypt_exp_2(r, random);
				end_time = omp_get_wtime();

				rns_time += end_time - start_time;

				mpz_clears(r, random, NULL);
			}
		}

		mpz_t result;
		mpz_init(result);

		dj.encrypt(result, cipher);
		mpz_set(cipher, result);

		if(s > 10) {

			gm_time /= TEST_CASES;
			rns_time /= TEST_CASES;

			cout << "s = " << s << endl;
			cout << "g^m:\t" << gm_time << "\tr^n^s:\t" << rns_time << "\n" << endl; 
		}
	
	}


return 0;

}
