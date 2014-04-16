#include "Client.h"
#include <cmath>

using namespace std;

#define FILE_SIZE 512
#define BIT_LENGTH 1024

int main() {

//	int s = (int) log2(FILE_SIZE);
	int s = 3;

	Client c(s, BIT_LENGTH, OCTO);

	unsigned char s_bits[] = { 1, 1, 1, 1, 0, 1, 0, 0, 0 };
	mpz_t *results = new mpz_t[7 * s];

	for (int i = 0; i < 7 * s; i++) {
		mpz_init(results[i]);
	}

	c.encrypt_s_bits(results, 7 * s, s_bits, 3 * s);

	int r_i = 0;

	for (int i = 0; i < s; i++) {

		mpz_t *decr = new mpz_t[7];

		for (int j = 0; j < 7; j++) {
			mpz_init(decr[j]);
		}

		c.dj->set_s(i + 1);

		for (int j = 0; j < 7; j++) {

			c.dj->decrypt(decr[j], results[r_i]);
			cout << decr[j] << endl;
			r_i++;
		}

		cout << "decrypted" << endl;

		for (int j = 0; j < 7; j++) {
			mpz_clear(decr[j]);
		}
	}

	return 0;
}
