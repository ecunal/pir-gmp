#include "Client.h"

using namespace std;

Client::Client(int m_s, int b_length, TreeType t) {

	max_s = m_s;
	tree = t;
	bit_length = b_length;

	// Initialize DamgardJurik

	dj = DamgardJurik(b_length, 1);
}

void Client::get_pub_keys(mpz_t n, mpz_t g) {

	mpz_set(n, dj.n);
	mpz_set(g, dj.g);

}

void Client::encrypt_s_bits(mpz_t result[], unsigned char s_bits[], int s_bit_length) {

	int temp_s = max_s;

	if (tree == BINARY) {

		cout << "binary tree, encrypting selection bits" << endl;

		for (int i = 0; i < s_bit_length; i++) {

			dj.set_s(temp_s);

			mpz_t enc, m;
			mpz_inits(enc, m, NULL);

			mpz_set_ui(m, s_bits[i]);

			dj.encrypt(enc, m);

			mpz_set(result[s_bit_length-i-1], enc);

			temp_s--;
		}

	}
}
