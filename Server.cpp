#include "Server.h"
#include <assert.h>

using namespace std;

extern unsigned int CORE_SIZE;

int isPowerOfTwo(unsigned int x)
{
	return ((x != 0) && !(x & (x - 1)));
}

Server::Server(int b_length, int file_size, TreeType t, mpz_t n, mpz_t g) {

	bit_length = b_length;
	f_size = file_size;
	tree = t;
	subtree_size = f_size; // burdan kontrol edebiliriz scalable mı değil mi
	scale = 0;

	dj = new DamgardJurik(bit_length, 1, n, g);

	generate_files(false);
}

Server::Server(int b_length, int file_size, TreeType t, mpz_t n, mpz_t g,
		int subtree) {

	bit_length = b_length;
	f_size = file_size;
	tree = t;
	subtree_size = subtree;
	scale = 1;

	dj = new DamgardJurik(bit_length, 1, n, g);

	generate_files(false);
}

void Server::generate_files(bool debug) {

	files = new mpz_t[f_size];

	if (debug) { // generate small meaningful files for testing

		cout << "DEBUG!!" << endl;

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

double Server::get_file_scalable(mpz_t result, mpz_t s_bits[],
		mpz_t subt_bits[], int subt_length) {

	assert(isPowerOfTwo(CORE_SIZE));

	double time = 0;

	// First, process the files with subtree selection bits
	// new files will be like E(subtree bits)^f

	mpz_t sub_files[subtree_size]; // subtree_size = subtree file size

	dj->set_s(1);

	int files_per_core = subtree_size / CORE_SIZE;

	double start = omp_get_wtime();

#pragma omp parallel for
	for (int p = 0; p < CORE_SIZE; p++) {
		for (int i = 0; i < files_per_core; i++) {

			mpz_t new_f;
			mpz_init_set_ui(new_f, 1);

			for (int j = 0; j < subt_length; j++) {

				mpz_t temp;
				mpz_init_set(temp, subt_bits[j]);
				mpz_powm(temp, temp,
						files[i + p * files_per_core + j * subtree_size],
						dj->n_sp);

				mpz_mul(new_f, new_f, temp);
				mpz_mod(new_f, new_f, dj->n_sp);
				mpz_clear(temp);
			}

			mpz_init_set(sub_files[i + p * files_per_core], new_f);
			mpz_clear(new_f);
		}
	}

	double end = omp_get_wtime();

	time += end - start;

	// simdi sub_files'ta collapsed subtree olması lazım.
	// onu normal file haline getirelim ki eski versiyonu kullanabilelim
	// simdilik

	for (int i = 0; i < f_size; i++) {
		mpz_clear(files[i]);
	}

	delete[] files;

	files = new mpz_t[subtree_size];

	for (int i = 0; i < subtree_size; i++) {
		mpz_init_set(files[i], sub_files[i]);
	}

	f_size = subtree_size;

	time += get_file_new_p(result, s_bits);

	return time;
}

/************ NORMAL METHOD WITH PARALLELIZATION *************/

double Server::get_file(mpz_t result, mpz_t s_bits[], int parallel,
		int extra_prl) {

	omp_set_nested(1);
	omp_set_num_threads(CORE_SIZE);

	double time = 0;

	if (tree == BINARY) {

		int depth = (int) log2(f_size);

		/* R array stores the ciphertext from lower levels.
		 * It gets updated at the end of every level.
		 *
		 * Therefore the result that client gets is the only
		 * element of this array at the end of calculation.
		 *
		 * This is our tree in a way, but we don't store the
		 * ciphertexts that we already consumed, therefore
		 * this is more memory friendly.
		 * */

		mpz_t *R = new mpz_t[f_size];

		/* During the processing of a level, computed values
		 * are first stored in this temp array,
		 * then at the end of the level, they are transferred
		 * to R array
		 * */
		mpz_t *temp;

		/* These R and temp arrays are used in all of the methods
		 * utilizing binary trees, quadtrees, octrees and also
		 * parallelized ones.
		 * */

		/* Since at the bottom of the tree we find our original files,
		 * we instantiate the R array with those files.
		 * */
		for (int i = 0; i < f_size; i++) {
			mpz_init_set(R[i], files[i]);
		}

		int r_size = f_size, temp_size;

		// For each level, do the calculations in a bottom-up manner
		for (int i = 0; i < depth; i++) {

			temp_size = f_size / pow(2, i + 1);
			temp = new mpz_t[temp_size];

			for (int j = 0; j < temp_size; j++) {
				mpz_init(temp[j]);
			}

			// Change s according to the level
			dj->set_s(i + 1);

			double start_time = omp_get_wtime();

#pragma omp parallel for if(parallel)
			for (int j = 0; j < temp_size; j++) {

				mpz_t f0, f1, subf, R0, R1, R2, exp1, exp2;
				mpz_inits(f0, f1, subf, R0, R1, R2, exp1, exp2, NULL);

				mpz_set(f0, R[2 * j]);
				mpz_set(f1, R[2 * j + 1]);

#pragma omp parallel for if(extra_prl)
				for (int p = 0; p < 2; p++) {

					if (p == 0) {
						dj->encrypt(R0, f0); // encryption
					} else {
						mpz_sub(subf, f1, f0);
						mpz_powm(R1, s_bits[i], subf, dj->n_sp); // exponentiation
					}

				}

				// Multiply the encryption & exponentiation results
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

			dj->set_s(i + 1);

			double start_time = omp_get_wtime();

#pragma omp parallel for if(parallel)
			for (int j = 0; j < temp_size; j++) {

				mpz_t f0, f1, f2, f3, R0, R1, R2, R3, R00;
				mpz_inits(f0, f1, f2, f3, R0, R1, R2, R3, R00, NULL);

				mpz_t tmps[3];

				for (int p = 0; p < 3; p++) {
					mpz_init(tmps[p]);
				}

				mpz_set(f0, R[4 * j]);
				mpz_set(f1, R[4 * j + 1]);
				mpz_set(f2, R[4 * j + 2]);
				mpz_set(f3, R[4 * j + 3]);

				// selection bits: c0: s_bits[1], c1: s_bits[0], c01: s_bits[2]

#pragma omp parallel for if(extra_prl)
				for (int p = 0; p < 4; p++) {

					if (p == 0) {

						dj->encrypt(R0, f0);

					} else if (p == 1) {

						mpz_sub(tmps[0], f1, f0);
						mpz_powm(R1, s_bits[3 * i + 1], tmps[0], dj->n_sp);

					} else if (p == 2) {

						mpz_sub(tmps[1], f2, f0);
						mpz_powm(R2, s_bits[3 * i], tmps[1], dj->n_sp);

					} else if (p == 3) {

						mpz_add(tmps[2], f3, f0);
						mpz_sub(tmps[2], tmps[2], f2);
						mpz_sub(tmps[2], tmps[2], f1);
						mpz_powm(R3, s_bits[3 * i + 2], tmps[2], dj->n_sp);

					}

				}

				mpz_mul(R00, R0, R1);
				mpz_mod(R00, R00, dj->n_sp);
				mpz_mul(R00, R00, R2);
				mpz_mod(R00, R00, dj->n_sp);
				mpz_mul(R00, R00, R3);
				mpz_mod(R00, R00, dj->n_sp);

				mpz_set(temp[j], R00);

				for (int p = 0; p < 3; p++) {
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

			dj->set_s(i + 1 + scale);

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

				for (int p = 0; p < 8; p++) {
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

/****** ALL CORES GET A PART OF THE TREE METHOD ******/

double Server::get_file_new_p(mpz_t result, mpz_t s_bits[]) {

	double time = 0;

	omp_set_nested(1);
	omp_set_num_threads(CORE_SIZE);

	/******************** BIN **********************/

	if (tree == BINARY) {

		int depth = (int) log2(f_size);

		if (depth < 3) {
			cout << "no need for this for small file sizes, just use other one."
					<< endl;
			return -1;
		}

		assert(isPowerOfTwo(CORE_SIZE));

		int remaining_depth = (int) log2(CORE_SIZE);

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

			for (int i = 0; i < depth - remaining_depth; i++) { // -2 also relative to core size

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

		}

		double end = omp_get_wtime();
		time += end - start;

		// 4 thread's works are finished
		// now the last two levels will be processed
		// this part is like get_file method
		// with parallel and extra_prl both 1

		mpz_t *R = new mpz_t[CORE_SIZE];
		mpz_t *temp;

		for (int i = 0; i < CORE_SIZE; i++) {
			mpz_init_set(R[i], new_files[i]);
		}

		int r_size = CORE_SIZE, temp_size;

		for (int i = depth - remaining_depth; i < depth; i++) {

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
	}

	/******************** QUAD **********************/

	else if (tree == QUAD) {

		int depth = (int) (log2(f_size) / log2(4));
		int remaining_depth = (int) ceil(log2(CORE_SIZE) / log2(4));

		if (depth <= remaining_depth) {
			cout << "no need for this method, just use other one." << endl;
			return -1;
		}

		assert(isPowerOfTwo(CORE_SIZE));

		int last_lvl_size = (int) pow(4, remaining_depth);
		int file_per_core = (int) ( last_lvl_size / CORE_SIZE);

		mpz_t outputs[last_lvl_size];
		DamgardJurik *djs[CORE_SIZE];

		for (int p = 0; p < CORE_SIZE; p++) {

			djs[p] = new DamgardJurik(dj->bit_length, dj->s, dj->n, dj->g);
		}

		for(int INDEX=0; INDEX < last_lvl_size; INDEX++ ) {
			mpz_init(outputs[INDEX]);
		}

		double start = omp_get_wtime();

#pragma omp parallel for
		for (int p = 0; p < CORE_SIZE; p++) {

			for(int y=0; y<file_per_core; y++) {

				double local_time = 0;

				int r_size = f_size / last_lvl_size, temp_size;
				int new_f_size = f_size / last_lvl_size;

				mpz_t *R = new mpz_t[r_size];
				mpz_t *temp;

				for (int i = 0; i < r_size; i++) {
					mpz_init_set(R[i], files[p * file_per_core * r_size + y*r_size + i]);
				}

				for (int i = 0; i < depth - remaining_depth; i++) { // depth -1: last level is integration

					temp_size = new_f_size / pow(4, i + 1);
					temp = new mpz_t[temp_size];

					for (int j = 0; j < temp_size; j++) {
						mpz_init(temp[j]);
					}

					djs[p]->set_s(i + 1);

					double start_time = omp_get_wtime();

					for (int j = 0; j < temp_size; j++) {

						mpz_t f0, f1, f2, f3, R0, R1, R2, R3, R00;
						mpz_inits(f0, f1, f2, f3, R0, R1, R2, R3, R00, NULL);

						mpz_t tmps[3];

						for (int t = 0; t < 3; t++) {
							mpz_init(tmps[t]);
						}

						mpz_set(f0, R[4 * j]);
						mpz_set(f1, R[4 * j + 1]);
						mpz_set(f2, R[4 * j + 2]);
						mpz_set(f3, R[4 * j + 3]);

						djs[p]->encrypt(R0, f0);

						mpz_sub(tmps[0], f1, f0);
						mpz_powm(R1, s_bits[3 * i + 1], tmps[0], djs[p]->n_sp);

						mpz_sub(tmps[1], f2, f0);
						mpz_powm(R2, s_bits[3 * i], tmps[1], djs[p]->n_sp);

						mpz_add(tmps[2], f3, f0);
						mpz_sub(tmps[2], tmps[2], f2);
						mpz_sub(tmps[2], tmps[2], f1);
						mpz_powm(R3, s_bits[3 * i + 2], tmps[2], djs[p]->n_sp);

						mpz_mul(R00, R0, R1);
						mpz_mod(R00, R00, djs[p]->n_sp);
						mpz_mul(R00, R00, R2);
						mpz_mod(R00, R00, djs[p]->n_sp);
						mpz_mul(R00, R00, R3);
						mpz_mod(R00, R00, djs[p]->n_sp);

						mpz_set(temp[j], R00);

						for (int p = 0; p < 3; p++) {
							mpz_clear(tmps[p]);
						}

						mpz_clears(f0, f1, f2, f3, R0, R1, R2, R3, R00, NULL);
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

				mpz_set(outputs[p*file_per_core + y], R[0]);

				delete[] R;

				//cout << "for p = " << p << ", local time is: " << local_time
				//		<< endl;

			}

		}

		double end = omp_get_wtime();
		time += end - start;

		mpz_t *R = new mpz_t[last_lvl_size];
		mpz_t *temp;

		for (int i = 0; i < last_lvl_size; i++) {
			mpz_init_set(R[i], outputs[i]);
		}

		int r_size = last_lvl_size, temp_size;

		for (int i = depth - remaining_depth; i < depth; i++) {

			temp_size = f_size / pow(4, i + 1);
			temp = new mpz_t[temp_size];

			for (int j = 0; j < temp_size; j++) {
				mpz_init(temp[j]);
			}

			dj->set_s(i + 1);

			double start_time = omp_get_wtime();

#pragma omp parallel for
			for (int j = 0; j < temp_size; j++) {

				mpz_t f0, f1, f2, f3, R0, R1, R2, R3, R00;
				mpz_inits(f0, f1, f2, f3, R0, R1, R2, R3, R00, NULL);

				mpz_t tmps[3];

				for (int p = 0; p < 3; p++) {
					mpz_init(tmps[p]);
				}

				mpz_set(f0, R[4 * j]);
				mpz_set(f1, R[4 * j + 1]);
				mpz_set(f2, R[4 * j + 2]);
				mpz_set(f3, R[4 * j + 3]);

				// selection bits: c0: s_bits[1], c1: s_bits[0], c01: s_bits[2]

#pragma omp parallel for
				for (int p = 0; p < 4; p++) {

					if (p == 0) {

						dj->encrypt(R0, f0);

					} else if (p == 1) {

						mpz_sub(tmps[0], f1, f0);
						mpz_powm(R1, s_bits[3 * i + 1], tmps[0], dj->n_sp);

					} else if (p == 2) {

						mpz_sub(tmps[1], f2, f0);
						mpz_powm(R2, s_bits[3 * i], tmps[1], dj->n_sp);

					} else if (p == 3) {

						mpz_add(tmps[2], f3, f0);
						mpz_sub(tmps[2], tmps[2], f2);
						mpz_sub(tmps[2], tmps[2], f1);
						mpz_powm(R3, s_bits[3 * i + 2], tmps[2], dj->n_sp);

					}

				}

				mpz_mul(R00, R0, R1);
				mpz_mod(R00, R00, dj->n_sp);
				mpz_mul(R00, R00, R2);
				mpz_mod(R00, R00, dj->n_sp);
				mpz_mul(R00, R00, R3);
				mpz_mod(R00, R00, dj->n_sp);

				mpz_set(temp[j], R00);

				for (int p = 0; p < 3; p++) {
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

	}

	/******************** OCTO **********************/

	else if (tree == OCTO) {

		int depth = (int) (log2(f_size) / log2(8));
		int remaining_depth = (int) ceil(log2(CORE_SIZE) / log2(8));

//		cout << "depth: " << depth << "rem.depth: "<< remaining_depth << endl;

		if (depth <= remaining_depth) {
			cout << "no need for this method, just use other one." << endl;
			return -1;
		}

		assert(isPowerOfTwo(CORE_SIZE));

		int last_lvl_size = (int) pow(8, remaining_depth);
		int file_per_core = (int) ( last_lvl_size / CORE_SIZE);

//		cout << "last lvl size: " << last_lvl_size << "file_per_core: "<< file_per_core << endl;

		mpz_t outputs[last_lvl_size];
		DamgardJurik *djs[CORE_SIZE];

		for (int p = 0; p < CORE_SIZE; p++) {

			djs[p] = new DamgardJurik(dj->bit_length, dj->s, dj->n, dj->g);
		}

		for(int INDEX=0; INDEX < last_lvl_size; INDEX++ ) {
			mpz_init(outputs[INDEX]);
		}

		double start = omp_get_wtime();

#pragma omp parallel for
		for (int p = 0; p < CORE_SIZE; p++) {

			for(int y=0; y<file_per_core; y++) {

				double local_time = 0;

				int r_size = f_size / last_lvl_size, temp_size;
				int new_f_size = f_size / last_lvl_size;

				mpz_t *R = new mpz_t[r_size];
				mpz_t *temp;

				for (int i = 0; i < r_size; i++) {
					mpz_init_set(R[i], files[p * file_per_core * r_size + y*r_size + i]);
				}

				for (int i = 0; i < depth - remaining_depth; i++) {

					temp_size = new_f_size / pow(8, i + 1);
					temp = new mpz_t[temp_size];

					for (int j = 0; j < temp_size; j++) {
						mpz_init(temp[j]);
					}

					djs[p]->set_s(i + 1);

					double start_time = omp_get_wtime();

					for (int j = 0; j < temp_size; j++) {

						mpz_t f[8], r[8], tmps[8];

						for (int t = 0; t < 8; t++) {

							mpz_init_set(f[t], R[8 * j + t]);
							mpz_inits(r[t], tmps[t], NULL);
						}

						int x = 0;

						djs[p]->encrypt(r[x], f[0]);

						x++;

						mpz_sub(tmps[x], f[1], f[0]);
						mpz_powm(r[x], s_bits[7 * i], tmps[x], djs[p]->n_sp);

						x++;

						mpz_sub(tmps[x], f[2], f[0]);
						mpz_powm(r[x], s_bits[7 * i + 1], tmps[x], djs[p]->n_sp);

						x++;

						mpz_sub(tmps[x], f[4], f[0]);
						mpz_powm(r[x], s_bits[7 * i + 2], tmps[x], djs[p]->n_sp);

						x++;

						mpz_add(tmps[x], f[3], f[0]);
						mpz_sub(tmps[x], tmps[x], f[2]);
						mpz_sub(tmps[x], tmps[x], f[1]);
						mpz_powm(r[x], s_bits[7 * i + 3], tmps[x], djs[p]->n_sp);

						x++;

						mpz_add(tmps[x], f[5], f[0]);
						mpz_sub(tmps[x], tmps[x], f[4]);
						mpz_sub(tmps[x], tmps[x], f[1]);
						mpz_powm(r[x], s_bits[7 * i + 4], tmps[x], djs[p]->n_sp);

						x++;

						mpz_add(tmps[x], f[6], f[0]);
						mpz_sub(tmps[x], tmps[x], f[2]);
						mpz_sub(tmps[x], tmps[x], f[4]);
						mpz_powm(r[x], s_bits[7 * i + 5], tmps[x], djs[p]->n_sp);

						x++;

						mpz_add(tmps[x], f[7], f[4]);
						mpz_add(tmps[x], tmps[x], f[2]);
						mpz_add(tmps[x], tmps[x], f[1]);

						mpz_sub(tmps[x], tmps[x], f[6]);
						mpz_sub(tmps[x], tmps[x], f[5]);
						mpz_sub(tmps[x], tmps[x], f[3]);
						mpz_sub(tmps[x], tmps[x], f[0]);

						mpz_powm(r[x], s_bits[7 * i + 6], tmps[x], djs[p]->n_sp);

						// NOW MULTIPLY THE RESULTS

						mpz_t RR;
						mpz_init_set_ui(RR, 1);

						for (int t = 0; t < 8; t++) {
							mpz_mul(RR, RR, r[t]);
							mpz_mod(RR, RR, djs[p]->n_sp);
						}

						mpz_set(temp[j], RR);

						for (int t = 0; t < 8; t++) {
							mpz_clears(f[t], r[t], tmps[t], NULL);
						}

					} // end of for j

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

				mpz_set(outputs[p*file_per_core + y], R[0]);

				delete[] R;
				// cout << "for p = " << p << " local time is: " << local_time << endl;

			}

		}

		// COLLAB AGAIN

		double end = omp_get_wtime();
		time += end - start;

		mpz_t *R = new mpz_t[last_lvl_size];
		mpz_t *temp;

		for (int i = 0; i < last_lvl_size; i++) {
			mpz_init_set(R[i], outputs[i]);
		}

		int r_size = last_lvl_size, temp_size;

		for (int i = depth - remaining_depth; i < depth; i++) {

			temp_size = f_size / pow(8, i + 1);
			temp = new mpz_t[temp_size];

			for (int j = 0; j < temp_size; j++) {
				mpz_init(temp[j]);
			}

			dj->set_s(i + 1 + scale);

			double start_time = omp_get_wtime();

#pragma omp parallel for
			for (int j = 0; j < temp_size; j++) {

				mpz_t f[8], r[8], tmps[8];

				for (int p = 0; p < 8; p++) {
					mpz_init_set(f[p], R[8 * j + p]);
					mpz_inits(r[p], tmps[p], NULL);
				}

#pragma omp parallel for
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

				for (int p = 0; p < 8; p++) {
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

	} // end of octo

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
