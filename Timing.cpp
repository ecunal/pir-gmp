/*
 * Timing.cpp
 *
 *  Created on: 2 May 2014
 *      Author: ecem
 */

#include "Client.h"
#include "Server.h"

#define BIT_LENGTH 1024

using namespace std;

int main() {

	// BINARY DUZ PARALLEL

	TreeType tree = BINARY;

	for (int s = 1; s <= 9; s++) {

		int file_size = (int) pow(2, s);

		int TEST_CASES = (s < 7) ? (int) pow(2, (11 - s)) : 20;

		cout << TEST_CASES << endl;

		double c_enc_time = 0, s_enc_time = 0, c_decr_time = 0;

		for (int j = 0; j < TEST_CASES; j++) {

			mpz_t n, g, cipher, decr;
			mpz_inits(n, g, cipher, decr, NULL);

			Client client(s, BIT_LENGTH, tree);
			client.get_pub_keys(n, g);

			Server server(BIT_LENGTH, s, file_size, tree, n, g);

			unsigned char s_bits[s];
			mpz_t *results = new mpz_t[s];

			for (int i = 0; i < s; i++) {
				s_bits[i] = 1;
				mpz_init(results[i]);
			}

			c_enc_time += client.encrypt_s_bits(results, s, s_bits, s);
			s_enc_time += server.get_file(cipher, results, s, 1, 1); // 1: parallel
			c_decr_time += client.decr_file(decr, cipher);

			mpz_clears(n, g, cipher, decr, NULL);

			for (int i = 0; i < s; i++) {
				mpz_clear(results[i]);
			}

		}

		cout << "file size " << file_size << endl;
		cout << "client encryption (prl): " << (c_enc_time / TEST_CASES)
				<< endl;
		cout << "server encryption (prl): " << (s_enc_time / TEST_CASES)
				<< endl;
		cout << "client decryption: " << (c_decr_time / TEST_CASES) << "\n"
				<< endl;
	}

	return 0;
}
