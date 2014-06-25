#include "Client.h"

using namespace std;

Client::Client(int m_s, int b_length, TreeType t) {

	max_s = m_s;
	tree = t;
	bit_length = b_length;

	dj = new DamgardJurik(b_length, 1);
}

void Client::get_pub_keys(mpz_t n, mpz_t g) {

	mpz_set(n, dj->n);
	mpz_set(g, dj->g);

}

double Client::decr_file(mpz_t dec, mpz_t file) {

	mpz_t *decr_array = new mpz_t[max_s];
	int s = max_s;

	double start_time = omp_get_wtime();
	for (int i = 0; i < max_s; i++) {

		mpz_init(decr_array[i]);

		dj->set_s(s);

		if (i > 0) {
			dj->decrypt(decr_array[i], decr_array[i - 1]);
		} else {
			dj->decrypt(decr_array[i], file);
		}

		s--;
	}
	double end_time = omp_get_wtime();

	mpz_set(dec, decr_array[max_s - 1]);

	return end_time - start_time;
}

double Client::encrypt_s_bits(mpz_t result[], int result_length,
		unsigned char s_bits[], int s_bit_length) {

	int temp_s = max_s;
	double time_elapsed = 0;

	omp_set_nested(1);

	if (tree == BINARY) {

		DamgardJurik *djs[s_bit_length];

		for (int i = 0; i < s_bit_length; i++) {
			djs[i] = new DamgardJurik(dj->bit_length, temp_s, dj->n, dj->g);
			temp_s--;
		}

		double start_time = omp_get_wtime();

		#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < s_bit_length; i++) {

			mpz_t enc, m;
			mpz_inits(enc, m, NULL);

			mpz_set_ui(m, s_bits[i]);

			djs[i]->encrypt(enc, m);

			mpz_set(result[s_bit_length - i - 1], enc);

			mpz_clears(enc, m, NULL);

		}

		double end_time = omp_get_wtime();

		time_elapsed = end_time - start_time;

	} else if (tree == QUAD) {

		DamgardJurik *djs[s_bit_length/2];

		for (int i = 0; i < s_bit_length/2; i++) {
			djs[i] = new DamgardJurik(dj->bit_length, temp_s, dj->n, dj->g);
			temp_s--;
		}

		double start_time = omp_get_wtime();

		#pragma omp parallel for
		for (int i = 0; i < s_bit_length; i += 2) {

			mpz_t x0, x1, x01, enc0, enc1, enc01;
			mpz_inits(x0, x1, x01, enc0, enc1, enc01, NULL);

			mpz_set_ui(x1, s_bits[i]);
			mpz_set_ui(x0, s_bits[i + 1]);
			mpz_mul(x01, x0, x1);

			#pragma omp parallel for
			for (int p = 0; p < 3; p++) {

				if (p == 0) {
					djs[i/2]->encrypt(enc0, x0);
				} else if (p == 1) {
					djs[i/2]->encrypt(enc1, x1);
				} else if (p == 2) {
					djs[i/2]->encrypt(enc01, x01);
				}

			}

			mpz_set(result[result_length - (i/2) * 3 - 1], enc01);
			mpz_set(result[result_length - (i/2) * 3 - 2], enc0);
			mpz_set(result[result_length - (i/2) * 3 - 3], enc1);

			mpz_clears(x0, x1, x01, enc0, enc1, enc01, NULL);

		}

		double end_time = omp_get_wtime();

		time_elapsed = end_time - start_time;

	} else if (tree == OCTO) {

		DamgardJurik *djs[s_bit_length/3];

		for (int i = 0; i < s_bit_length/3; i++) {
			djs[i] = new DamgardJurik(dj->bit_length, temp_s, dj->n, dj->g);
			temp_s--;
		}

		double start_time = omp_get_wtime();

		#pragma omp parallel for
		for (int i = 0; i < s_bit_length; i += 3) {

			mpz_t *x = new mpz_t[7]; // 0, 1, 2, 3->01, 4->02, 5->12, 6->012
			mpz_t *enc = new mpz_t[7];

			for (int j = 0; j < 7; j++) {
				mpz_init(x[j]);
				mpz_init(enc[j]);
			}

			mpz_set_ui(x[2], s_bits[i]);
			mpz_set_ui(x[1], s_bits[i + 1]);
			mpz_set_ui(x[0], s_bits[i + 2]);

			mpz_mul(x[3], x[0], x[1]);
			mpz_mul(x[5], x[1], x[2]);
			mpz_mul(x[4], x[0], x[2]);

			mpz_mul(x[6], x[3], x[2]);

			#pragma omp parallel for
			for (int j = 6; j >= 0; j--) {
				djs[i/3]->encrypt(enc[j], x[j]);
				mpz_set(result[result_length - (i/3)*7 - (7-j)], enc[j]);
			}

			for (int j = 0; j < 7; j++) {
				mpz_clear(x[j]);
				mpz_clear(enc[j]);
			}
		}

		double end_time = omp_get_wtime();

		time_elapsed = end_time - start_time;

	}

	return time_elapsed;
}

