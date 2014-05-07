#include "Server.h"

using namespace std;

Server::Server(int b_length, int s, int file_size, TreeType t, mpz_t n,
		mpz_t g) {

	bit_length = b_length;
	f_size = file_size;
	max_s = s;
	tree = t;

	dj = new DamgardJurik(bit_length, 1, n, g);

	generate_files(false);
}

void Server::generate_files(bool debug) {

	files = new mpz_t[f_size];

	//cout << "file size " << f_size << endl;

	if (debug) { // generate meaningful files for testing

		for (int i = 0; i < f_size; i++) {
			mpz_init_set_ui(files[i], i);
		}

	} else { // generate 1024 bit files

		init_random();

		for (int i = 0; i < f_size; i++) {

			mpz_init(files[i]);
			mpz_urandomb(files[i], state, 1024); // 1024 bit random files
		}

	}
}
double Server::get_file(mpz_t result, mpz_t s_bits[], int s_length, int parallel, bool extra_prl) {

	double time = 0;

	if (tree == BINARY) {

		int depth = (int) log2(f_size);

		mpz_t *R = new mpz_t[f_size];
		mpz_t *temp;

		for (int i = 0; i < f_size; i++) {
			mpz_init_set(R[i], files[i]);
		}

		int r_size = f_size, temp_size;

		for (int i = 0; i < depth; i++) {

			temp_size = f_size / pow(2, i + 1);
			temp = new mpz_t[temp_size];

			for (int j = 0; j < temp_size; j++) {
				mpz_init(temp[j]);
			}

			dj->set_s(i+1);

			double start_time = omp_get_wtime();

			omp_set_nested(1);

			#pragma omp parallel for if(parallel)
			for (int j = 0; j < temp_size; j++) {

				mpz_t f0, f1, subf, R0, R1, R2, exp1, exp2;
				mpz_inits(f0, f1, subf, R0, R1, R2, exp1, exp2, NULL);

				mpz_set(f0, R[2*j]);
				mpz_set(f1, R[2*j + 1]);

				if(!extra_prl) {
					dj->encrypt(R0, f0);
				} else {

					#pragma omp parallel for
					for(int p=0; p<2; p++) {
						if(p == 0) {
							dj->encrypt_exp_1(exp1, f0);
						} else {
							dj->encrypt_exp_2(exp2);
						}
					}

					dj->encrypt_mult(R0, exp1, exp2);
				}


				mpz_sub(subf, f1, f0);
				mpz_powm(R1, s_bits[i], subf, dj->n_sp);

				mpz_mul(R2, R0, R1);
				mpz_mod(R2, R2, dj->n_sp);

				mpz_set(temp[j], R2);

				mpz_clears(f0, f1, subf, R0, R1, R2, exp1, exp2, NULL);
			}

			double end_time = omp_get_wtime();

			time += end_time - start_time;

			for (int j = 0; j < r_size; j++) {
				mpz_clear(R[j]);
			}

			delete[] R;

			r_size = temp_size;
			R = new mpz_t[r_size];

			for (int j = 0; j < r_size; j++) {
				mpz_init_set(R[j], temp[j]);
			}

			for(int j=0; j<temp_size; j++)  {
				mpz_clear(temp[j]);
			}

			delete [] temp;

		} // end of all depths

//		cout << "end of all depths, R size: " << r_size << endl;

		mpz_set(result, R[0]);

		delete [] R;

	} // end of binary case

	return time;
}

void Server::init_random() {

	int random_dev = open("/dev/urandom", O_RDONLY);
	long random_seed;
	size_t random_len = 0;

	while (random_len < sizeof random_seed) {

		ssize_t result = read(random_dev, ((char*) &random_seed) + random_len,
				(sizeof random_seed) - random_len);
		if (result < 0) {
			cout << "unable to read from urandom" << endl;
		}
		random_len += result;
	}

	close(random_dev);

	gmp_randinit_mt(state);
	gmp_randseed_ui(state, random_seed);

}
