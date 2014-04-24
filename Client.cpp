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

void Client::decr_file(mpz_t dec, mpz_t file) {

	mpz_t *decr_array = new mpz_t[max_s];
	int s = max_s;

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

	mpz_set(dec, decr_array[max_s - 1]);
}

void Client::encrypt_s_bits(mpz_t result[], int result_length,
		unsigned char s_bits[], int s_bit_length) {

	int temp_s = max_s;

	if (tree == BINARY) {

		cout << "binary tree, encrypting selection bits" << endl;

		DamgardJurik *djs[s_bit_length];

		for (int i = 0; i < s_bit_length; i++) {
			djs[i] = new DamgardJurik(dj->bit_length, temp_s, dj->n, dj->g);
			temp_s--;
		}

//		clock_t begin, end;
//		double time_spent;
//
//		begin = clock();

		double start_time = omp_get_wtime();

		#pragma omp parallel for
		for (int i = 0; i < s_bit_length; i++) {

			mpz_t enc, m;
			mpz_inits(enc, m, NULL);

			mpz_set_ui(m, s_bits[i]);

			djs[i]->encrypt(enc, m);

			mpz_set(result[s_bit_length - i - 1], enc);

			mpz_clears(enc, m, NULL);

		}

		double end_time = omp_get_wtime();

		double time_elapsed = end_time - start_time;

		cout << "enc of selection bits: " << time_elapsed << endl;

	} else if (tree == QUAD) {

		cout << "quadtree, encrypting selection bits" << endl;

		int r_i = result_length - 1;

		for (int i = 0; i < s_bit_length; i += 2) {

			dj->set_s(temp_s);

			mpz_t x0, x1, x01, enc0, enc1, enc01;
			mpz_inits(x0, x1, x01, enc0, enc1, enc01, NULL);

			mpz_set_ui(x1, s_bits[i]);
			mpz_set_ui(x0, s_bits[i + 1]);
			mpz_mul(x01, x0, x1);

			dj->encrypt(enc0, x0);
			dj->encrypt(enc1, x1);
			dj->encrypt(enc01, x01);

			mpz_set(result[r_i], enc01);
			r_i--;
			mpz_set(result[r_i], enc0);
			r_i--;
			mpz_set(result[r_i], enc1);
			r_i--;

			mpz_clears(x0, x1, x01, enc0, enc1, enc01, NULL);

			temp_s--;
		}

	} else if (tree == OCTO) {

		cout << "octree, encrypting selection bits" << endl;

		int r_i = result_length - 1;

		for (int i = 0; i < s_bit_length; i += 3) {

			dj->set_s(temp_s);

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

			for (int j = 6; j >= 0; j--) {
				dj->encrypt(enc[j], x[j]);
				mpz_set(result[r_i], enc[j]);
				r_i--;
			}

			for (int j = 0; j < 7; j++) {
				mpz_clear(x[j]);
				mpz_clear(enc[j]);
			}

			temp_s--;
		}

	}
}

