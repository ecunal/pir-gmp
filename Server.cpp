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

		cout << "DEBUG!!!!1" << endl;

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


double Server::get_file_new_p(mpz_t result, mpz_t s_bits[], int s_length) {

	double time = 0;

	if (tree != BINARY) {
		cout << "other tree types are not supported yet" << endl;
		return -1;
	}

	int depth = (int) log2(f_size);

	if (depth < 3) {
		cout << "no need for this for small file sizes, just use other one."
				<< endl;
		return -1;
	}

	omp_set_nested(1);

	mpz_t *new_files = new mpz_t[CORE_SIZE];
	DamgardJurik *djs[CORE_SIZE];

	for (int p = 0; p < CORE_SIZE; p++) {
		djs[p] = new DamgardJurik(dj->bit_length, dj->s, dj->n, dj->g);
		mpz_init(new_files[p]);
	}

	double start = omp_get_wtime();

	#pragma omp parallel for
	for (int p = 0; p < CORE_SIZE; p++) {

		double local_time = 0;

		int r_size = f_size / CORE_SIZE, temp_size;
		int new_f_size = f_size / CORE_SIZE;

		mpz_t *R = new mpz_t[r_size];
		mpz_t *temp;

		for (int i = 0; i < r_size; i++) {
			mpz_init_set(R[i], files[p * r_size + i]);
		}

		for (int i = 0; i < depth - 2; i++) { // -2 also relative to core size

			temp_size = new_f_size / pow(2, i + 1);
			temp = new mpz_t[temp_size];

			for (int j = 0; j < temp_size; j++) {
				mpz_init(temp[j]);
			}

			djs[p]->set_s(i + 1);

			double start_time = omp_get_wtime();

			for (int j = 0; j < temp_size; j++) {


				mpz_t f0, f1, subf, R0, R1, R2;
				mpz_inits(f0, f1, subf, R0, R1, R2, NULL);

				mpz_set(f0, R[2 * j]);
				mpz_set(f1, R[2 * j + 1]);

				djs[p]->encrypt(R0, f0);
				mpz_sub(subf, f1, f0);
				mpz_powm(R1, s_bits[i], subf, djs[p]->n_sp);

				mpz_mul(R2, R0, R1);
				mpz_mod(R2, R2, djs[p]->n_sp);

				mpz_set(temp[j], R2);

				mpz_clears(f0, f1, subf, R0, R1, R2, NULL);

			}

			double end_time = omp_get_wtime();

			local_time += end_time - start_time;

			for (int j = 0; j < r_size; j++) {
				mpz_clear(R[j]);
			}

			delete[] R;

			r_size = temp_size;
			R = new mpz_t[r_size];

			for (int j = 0; j < r_size; j++) {
				mpz_init_set(R[j], temp[j]);
			}

			for (int j = 0; j < temp_size; j++) {
				mpz_clear(temp[j]);
			}

			delete[] temp;

		}

		mpz_set(new_files[p], R[0]);

		delete[] R;

//		cout << "for p = " << p << ", local time is: " << local_time << endl;
	}

	double end = omp_get_wtime();
	time += end - start;

//	cout << "global time: " << time << endl;

	// following file sizes are also equal to CORE_SIZE

	mpz_t *R = new mpz_t[CORE_SIZE];
	mpz_t *temp;

	for (int i = 0; i < CORE_SIZE; i++) {
		mpz_init_set(R[i], new_files[i]);
	}

	int r_size = CORE_SIZE, temp_size;

	for (int i = depth - 2; i < depth; i++) {

		temp_size = f_size / pow(2, i + 1);
		temp = new mpz_t[temp_size];

		for (int j = 0; j < temp_size; j++) {
			mpz_init(temp[j]);
		}

		dj->set_s(i + 1);

		double start_time = omp_get_wtime();

		#pragma omp parallel for
		for (int j = 0; j < temp_size; j++) {

			mpz_t f0, f1, subf, R0, R1, R2;
			mpz_inits(f0, f1, subf, R0, R1, R2, NULL);

			mpz_set(f0, R[2 * j]);
			mpz_set(f1, R[2 * j + 1]);

			#pragma omp parallel for
			for (int p = 0; p < 2; p++) {

				if (p == 0) {
					dj->encrypt(R0, f0);
				} else {
					mpz_sub(subf, f1, f0);
					mpz_powm(R1, s_bits[i], subf, dj->n_sp);
				}

			}

			mpz_mul(R2, R0, R1);
			mpz_mod(R2, R2, dj->n_sp);

			mpz_set(temp[j], R2);

			mpz_clears(f0, f1, subf, R0, R1, R2, NULL);
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

		for (int j = 0; j < temp_size; j++) {
			mpz_clear(temp[j]);
		}

		delete[] temp;

	} // end of all depths

	mpz_set(result, R[0]);

	delete[] R;

	return time;
}


double Server::get_file(mpz_t result, mpz_t s_bits[], int s_length,
		int parallel, int extra_prl) {

	omp_set_nested(1);

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

			dj->set_s(i + 1);

			double start_time = omp_get_wtime();

			#pragma omp parallel for if(parallel)
			for (int j = 0; j < temp_size; j++) {

				mpz_t f0, f1, subf, R0, R1, R2, exp1, exp2;
				mpz_inits(f0, f1, subf, R0, R1, R2, exp1, exp2, NULL);

				mpz_set(f0, R[2 * j]);
				mpz_set(f1, R[2 * j + 1]);

				#pragma omp parallel for if(extra_prl)
				for(int p=0; p<2; p++) {

					if(p == 0) {
						dj->encrypt(R0, f0);
					} else {
						mpz_sub(subf, f1, f0);
						mpz_powm(R1, s_bits[i], subf, dj->n_sp);
					}

				}

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

			for (int j = 0; j < temp_size; j++) {
				mpz_clear(temp[j]);
			}

			delete[] temp;

		} // end of all depths

		mpz_set(result, R[0]);

		delete[] R;


		/****** END OF BINARY *****/


	} else if (tree == QUAD) {


		int depth = (int) (log2(f_size) / log2(4));

		mpz_t *R = new mpz_t[f_size];
		mpz_t *temp;

		for (int i = 0; i < f_size; i++) {
			mpz_init_set(R[i], files[i]);
		}

		int r_size = f_size, temp_size;

		for (int i = 0; i < depth; i++) {

			temp_size = f_size / pow(4, i + 1);
			temp = new mpz_t[temp_size];

			for (int j = 0; j < temp_size; j++) {
				mpz_init(temp[j]);
			}

			dj->set_s(i+1);

			double start_time = omp_get_wtime();

			#pragma omp parallel for if(parallel)
			for(int j=0; j<temp_size; j++) {

				mpz_t f0, f1, f2, f3, R0, R1, R2, R3, R00;
				mpz_inits(f0, f1, f2, f3, R0, R1, R2, R3, R00, NULL);

				mpz_t tmps[3];

				for(int p=0; p<3; p++) {
					mpz_init(tmps[p]);
				}

				mpz_set(f0, R[4*j]);
				mpz_set(f1, R[4*j + 1]);
				mpz_set(f2, R[4*j + 2]);
				mpz_set(f3, R[4*j + 3]);

				// selection bits: c0: s_bits[1], c1: s_bits[0], c01: s_bits[2]


				#pragma omp parallel for if(extra_prl)
				for(int p=0; p<4; p++) {

					if(p == 0) {

						dj->encrypt(R0, f0);

					} else if(p == 1) {

						mpz_sub(tmps[0], f1, f0);
						mpz_powm(R1, s_bits[3*i + 1], tmps[0], dj->n_sp);

					} else if(p == 2) {

						mpz_sub(tmps[1], f2, f0);
						mpz_powm(R2, s_bits[3*i], tmps[1], dj->n_sp);

					} else if(p == 3) {

						mpz_add(tmps[2], f3, f0);
						mpz_sub(tmps[2], tmps[2], f2);
						mpz_sub(tmps[2], tmps[2], f1);
						mpz_powm(R3, s_bits[3*i + 2], tmps[2], dj->n_sp);

					}

				}

				mpz_mul(R00, R0, R1);
				mpz_mod(R00, R00, dj->n_sp);
				mpz_mul(R00, R00, R2);
				mpz_mod(R00, R00, dj->n_sp);
				mpz_mul(R00, R00, R3);
				mpz_mod(R00, R00, dj->n_sp);

				mpz_set(temp[j], R00);

				for(int p=0; p<3; p++) {
					mpz_clear(tmps[p]);
				}

				mpz_clears(f0, f1, f2, f3, R0, R1, R2, R3, R00, NULL);
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

			for (int j = 0; j < temp_size; j++) {
				mpz_clear(temp[j]);
			}

			delete[] temp;

		}

		mpz_set(result, R[0]);

		/******* END OF QUAD *******/




	} else if (tree == OCTO) {


		int depth = (int) (log2(f_size) / log2(8));

		mpz_t *R = new mpz_t[f_size];
		mpz_t *temp;

		for (int i = 0; i < f_size; i++) {
			mpz_init_set(R[i], files[i]);
		}

		int r_size = f_size, temp_size;

		for (int i = 0; i < depth; i++) {

			temp_size = f_size / pow(8, i + 1);
			temp = new mpz_t[temp_size];

			for (int j = 0; j < temp_size; j++) {
				mpz_init(temp[j]);
			}

			dj->set_s(i + 1);

			double start_time = omp_get_wtime();

			#pragma omp parallel for if(parallel)
			for (int j = 0; j < temp_size; j++) {

				mpz_t f[8], r[8], tmps[8];

				for (int p = 0; p < 8; p++) {
					mpz_init_set(f[p], R[8 * j + p]);
					mpz_inits(r[p], tmps[p], NULL);
				}

				#pragma omp parallel for if(extra_prl)
				for (int p = 0; p < 8; p++) {

					if (p == 0) {
						dj->encrypt(r[p], f[0]);
					} else if (p == 1) {

						mpz_sub(tmps[p - 1], f[1], f[0]);
						mpz_powm(r[p], s_bits[7 * i], tmps[p - 1], dj->n_sp);

					} else if (p == 2) {

						mpz_sub(tmps[p - 1], f[2], f[0]);
						mpz_powm(r[p], s_bits[7 * i + 1], tmps[p - 1],
								dj->n_sp);

					} else if (p == 3) {

						mpz_sub(tmps[p - 1], f[4], f[0]);
						mpz_powm(r[p], s_bits[7 * i + 2], tmps[p - 1],
								dj->n_sp);

					} else if (p == 4) {

						mpz_add(tmps[p - 1], f[3], f[0]);
						mpz_sub(tmps[p - 1], tmps[p - 1], f[2]);
						mpz_sub(tmps[p - 1], tmps[p - 1], f[1]);
						mpz_powm(r[p], s_bits[7 * i + 3], tmps[p - 1],
								dj->n_sp);

					} else if (p == 5) {

						mpz_add(tmps[p - 1], f[5], f[0]);
						mpz_sub(tmps[p - 1], tmps[p - 1], f[4]);
						mpz_sub(tmps[p - 1], tmps[p - 1], f[1]);
						mpz_powm(r[p], s_bits[7 * i + 4], tmps[p - 1],
								dj->n_sp);

					} else if (p == 6) {

						mpz_add(tmps[p - 1], f[6], f[0]);
						mpz_sub(tmps[p - 1], tmps[p - 1], f[2]);
						mpz_sub(tmps[p - 1], tmps[p - 1], f[4]);
						mpz_powm(r[p], s_bits[7 * i + 5], tmps[p - 1],
								dj->n_sp);

					} else if (p == 7) {

						mpz_add(tmps[p - 1], f[7], f[4]);
						mpz_add(tmps[p - 1], tmps[p - 1], f[2]);
						mpz_add(tmps[p - 1], tmps[p - 1], f[1]);

						mpz_sub(tmps[p - 1], tmps[p - 1], f[6]);
						mpz_sub(tmps[p - 1], tmps[p - 1], f[5]);
						mpz_sub(tmps[p - 1], tmps[p - 1], f[3]);
						mpz_sub(tmps[p - 1], tmps[p - 1], f[0]);

						mpz_powm(r[p], s_bits[7 * i + 6], tmps[p - 1],
								dj->n_sp);

					}

				}

				// NOW MULTIPLY THE RESULTS

				mpz_t RR;
				mpz_init_set_ui(RR, 1);

				for(int p = 0; p < 8; p++) {
					mpz_mul(RR, RR, r[p]);
					mpz_mod(RR, RR, dj->n_sp);
				}

				mpz_set(temp[j], RR);


				for (int p = 0; p < 8; p++) {
					mpz_clears(f[p], r[p], tmps[p], NULL);
				}

			} // end of for j

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

			for (int j = 0; j < temp_size; j++) {
				mpz_clear(temp[j]);
			}

			delete[] temp;

		}

		mpz_set(result, R[0]);


		/******* END OF OCTO ******/

	}


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
