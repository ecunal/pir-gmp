/*
 * Timing.cpp
 *
 *  Created on: 2 May 2014
 *      Author: ecem
 */

#include "Client.h"
#include "Server.h"
#include <vector>
#include <algorithm>

#define BIT_LENGTH 1024

using namespace std;

int main() {

	// BINARY

	TreeType tree = BINARY;

	for (int s = 6; s <= 9; s++) {

		int file_size = (int) pow(2, s);

		int TEST_CASES = /*(s < 7) ? ((int) pow(2, (11-s))) : 25 */ 100;

		cout << "\nTest Cases: " << TEST_CASES << endl;

		double c_enc_time = 0, s_enc_time = 0, c_decr_time = 0;

		vector<double> server_times;

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
			}

			for (int i = 0; i < s; i++) {
				mpz_init(results[i]);
			}

			c_enc_time += client.encrypt_s_bits(results, s, s_bits, s);

			double s_time = server.get_file_new_p(cipher, results, s); // NEW METHOD!
			s_enc_time += s_time;
			server_times.push_back(s_time);

			c_decr_time += client.decr_file(decr, cipher);

			mpz_clears(n, g, cipher, decr, NULL);

		}

		sort(server_times.begin(), server_times.end());
		double smallest = server_times.front();
		double average = 0;
		int count = 0;

		cout << "Smallest time: " << smallest << endl;

		for(vector<double>::iterator it = server_times.begin() ; it != server_times.end(); ++it ) {

			if(*it < smallest * 1.1) {
				average += *it;
				count++;
			} else {
				cout << "Stop on: " << *it << "\nCount: " << count << "\nSize was: " << server_times.size() << endl;
				break;
			}
		}

		

		cout << "~~~~~~~~ results ~~~~~~~~~" << endl;
		cout << "file size " << file_size << endl;
		cout << "client encryption (prl): " << (c_enc_time / TEST_CASES) << endl;
		cout << "server encryption (new method): " << (s_enc_time / TEST_CASES) << endl;
		cout << "server encryption (new method - k best): " << (average / count) << endl;
		cout << "client decryption: " << (c_decr_time / TEST_CASES) << "\n" << endl;
		cout << "~~~~~~~~ end ~~~~~~~~~" << endl;
	} 

	// QUAD

/*	TreeType tree = QUAD;

	for (int s = 1; s <= 4; s++) {

		int file_size = (int) pow(4, s);

		int TEST_CASES = (int) pow(2, (11-s));

		cout << TEST_CASES << endl;

		double c_enc_time = 0, s_enc_time = 0, c_decr_time = 0;

		for (int j = 0; j < TEST_CASES; j++) {

			mpz_t n, g, cipher, decr;
			mpz_inits(n, g, cipher, decr, NULL);

			Client client(s, BIT_LENGTH, tree);
			client.get_pub_keys(n, g);

			Server server(BIT_LENGTH, s, file_size, tree, n, g);

			unsigned char s_bits[2*s];
			mpz_t *results = new mpz_t[3*s];

			for (int i = 0; i < 2*s; i++) {
				s_bits[i] = 1;				
			}

			for (int i = 0; i < 3*s; i++) {
				mpz_init(results[i]);
			}

			c_enc_time += client.encrypt_s_bits(results, 3*s, s_bits, 2*s);
			s_enc_time += server.get_file(cipher, results, s, 1, 1); // 1: parallel
			c_decr_time += client.decr_file(decr, cipher);

			mpz_clears(n, g, cipher, decr, NULL);

		}

		cout << "file size " << file_size << endl;
		cout << "client encryption (seqo): " << (c_enc_time / TEST_CASES) << endl;
		cout << "server encryption (seqo): " << (s_enc_time / TEST_CASES) << endl;
		cout << "client decryption: " << (c_decr_time / TEST_CASES) << "\n" << endl;
	}  

	// OCTO 

	TreeType tree = OCTO;

	for (int s = 1; s <= 3; s++) {

		int file_size = (int) pow(8, s);

		int TEST_CASES = (int) pow(2, (11-s));

		cout << TEST_CASES << endl;

		double c_enc_time = 0, s_enc_time = 0, c_decr_time = 0;

		for (int j = 0; j < TEST_CASES; j++) {

			mpz_t n, g, cipher, decr;
			mpz_inits(n, g, cipher, decr, NULL);

			Client client(s, BIT_LENGTH, tree);
			client.get_pub_keys(n, g);

			Server server(BIT_LENGTH, s, file_size, tree, n, g);

			unsigned char s_bits[3*s];
			mpz_t *results = new mpz_t[7*s];

			for (int i = 0; i < 3*s; i++) {
				s_bits[i] = 1;				
			}

			for (int i = 0; i < 7*s; i++) {
				mpz_init(results[i]);
			}

			c_enc_time += client.encrypt_s_bits(results, 7*s, s_bits, 3*s);
			s_enc_time += server.get_file(cipher, results, s, 1, 1); // 1: parallel
			c_decr_time += client.decr_file(decr, cipher);

			mpz_clears(n, g, cipher, decr, NULL);

		}

		cout << "file size " << file_size << endl;
		cout << "client encryption (seqo): " << (c_enc_time / TEST_CASES) << endl;
		cout << "server encryption (prl2): " << (s_enc_time / TEST_CASES) << endl;
		cout << "client decryption: " << (c_decr_time / TEST_CASES) << "\n" << endl;
	} */

	return 0;
}
