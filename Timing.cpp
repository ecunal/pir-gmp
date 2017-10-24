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

using namespace std;

#define BIT_LENGTH 1024

unsigned int CORE_SIZE;
unsigned int TEST_CASES;

/******* BINARIES *******/

void binary_new_p(int start_s, int end_s) {

	TreeType tree = BINARY;

	cout << "\n**BINARY NEW PARALLEL**\n" << endl;

	for (int s = start_s; s <= end_s; s++) {

		int file_size = (int) pow(2, s);

		int TEST_CASES = (s < 7) ? ((int) pow(2, (11-s))) : 25; //  100;
		cout << "\nTest Cases: " << TEST_CASES << endl;

		double c_enc_time = 0, s_enc_time = 0, c_decr_time = 0;

		vector<double> server_times;

		for (int j = 0; j < TEST_CASES; j++) {

			mpz_t n, g, cipher, decr;
			mpz_inits(n, g, cipher, decr, NULL);

			Client client(s, BIT_LENGTH, tree);
			client.get_pub_keys(n, g);

			Server server(BIT_LENGTH, file_size, tree, n, g);

			unsigned char s_bits[s];
			mpz_t *results = new mpz_t[s];

			for (int i = 0; i < s; i++) {
				s_bits[i] = 1;
			}

			for (int i = 0; i < s; i++) {
				mpz_init(results[i]);
			}

			c_enc_time += client.encrypt_s_bits(results, s, s_bits, s);

			double s_time = server.get_file_new_p(cipher, results);
			s_enc_time += s_time;
			server_times.push_back(s_time);
			c_decr_time += client.decr_file(decr, cipher);

			mpz_clears(n, g, cipher, decr, NULL);

			for (int i = 0; i < s; i++) {
				mpz_clear(results[i]);
			}

			delete[] results;

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

}

void binary_full_p(int start_s, int end_s) {


	cout << "**BINARY FULL PARALLEL**\n" << endl;

	TreeType tree = BINARY;

	for (int s = start_s; s <= end_s; s++) {

		int file_size = (int) pow(2, s);

		int TEST_CASES = (s < 7) ? ((int) pow(2, (11-s))) : 25; //  100;
		cout << "\nTest Cases: " << TEST_CASES << endl;

		double c_enc_time = 0, s_enc_time = 0, c_decr_time = 0;

		vector<double> server_times;

		for (int j = 0; j < TEST_CASES; j++) {

			mpz_t n, g, cipher, decr;
			mpz_inits(n, g, cipher, decr, NULL);

			Client client(s, BIT_LENGTH, tree);
			client.get_pub_keys(n, g);

			Server server(BIT_LENGTH, file_size, tree, n, g);

			unsigned char s_bits[s];
			mpz_t *results = new mpz_t[s];

			for (int i = 0; i < s; i++) {
				s_bits[i] = 1;
			}

			for (int i = 0; i < s; i++) {
				mpz_init(results[i]);
			}

			c_enc_time += client.encrypt_s_bits(results, s, s_bits, s);

			double s_time = server.get_file(cipher, results, 1, 1);
			s_enc_time += s_time;
			server_times.push_back(s_time);
			c_decr_time += client.decr_file(decr, cipher);

			mpz_clears(n, g, cipher, decr, NULL);

			for (int i = 0; i < s; i++) {
				mpz_clear(results[i]);
			}

			delete[] results;

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

}

void binary_serial(int start_s, int end_s) {


	cout << "**BINARY SEQUENTIAL**\n" << endl;

	TreeType tree = BINARY;

	for (int s = start_s; s <= end_s; s++) {

		int file_size = (int) pow(2, s);

		int TEST_CASES = (s < 7) ? ((int) pow(2, (11-s))) : 25; //  100;
		cout << "\nTest Cases: " << TEST_CASES << endl;

		double c_enc_time = 0, s_enc_time = 0, c_decr_time = 0;

		vector<double> server_times;

		for (int j = 0; j < TEST_CASES; j++) {

			mpz_t n, g, cipher, decr;
			mpz_inits(n, g, cipher, decr, NULL);

			Client client(s, BIT_LENGTH, tree);
			client.get_pub_keys(n, g);

			Server server(BIT_LENGTH, file_size, tree, n, g);

			unsigned char s_bits[s];
			mpz_t *results = new mpz_t[s];

			for (int i = 0; i < s; i++) {
				s_bits[i] = 1;
			}

			for (int i = 0; i < s; i++) {
				mpz_init(results[i]);
			}

			c_enc_time += client.encrypt_s_bits(results, s, s_bits, s);

			double s_time = server.get_file(cipher, results, 0, 0);
			s_enc_time += s_time;
			server_times.push_back(s_time);
			c_decr_time += client.decr_file(decr, cipher);

			mpz_clears(n, g, cipher, decr, NULL);

			for (int i = 0; i < s; i++) {
				mpz_clear(results[i]);
			}

			delete[] results;

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

}

/******* QUADS *******/

void quad_new_p(int start_s, int end_s) {

	TreeType tree = QUAD;

	cout << "**QUAD NEW PARALLEL**\n" << endl;

	for (int s = start_s; s <= end_s; s++) {

		int file_size = (int) pow(4, s);

		int TEST_CASES = (int) pow(2, (11-s));
		cout << TEST_CASES << endl;

		double c_enc_time = 0, s_enc_time = 0, c_decr_time = 0;

		vector<double> server_times;

		for (int j = 0; j < TEST_CASES; j++) {

			mpz_t n, g, cipher, decr;
			mpz_inits(n, g, cipher, decr, NULL);

			Client client(s, BIT_LENGTH, tree);
			client.get_pub_keys(n, g);

			Server server(BIT_LENGTH, file_size, tree, n, g);

			unsigned char s_bits[2*s];
			mpz_t *results = new mpz_t[3*s];

			for (int i = 0; i < 2*s; i++) {
				s_bits[i] = 1;
			}

			for (int i = 0; i < 3*s; i++) {
				mpz_init(results[i]);
			}

			c_enc_time += client.encrypt_s_bits(results, 3*s, s_bits, 2*s);

			double s_time = server.get_file_new_p(cipher, results);
			s_enc_time += s_time;
			server_times.push_back(s_time);

			c_decr_time += client.decr_file(decr, cipher);

			mpz_clears(n, g, cipher, decr, NULL);

			for (int i = 0; i < 3*s; i++) {
				mpz_clear(results[i]);
			}

			delete[] results;

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

}

void quad_full_p(int start_s, int end_s) {

	TreeType tree = QUAD;

	cout << "**QUAD FULL PARALLEL**\n" << endl;

	for (int s = start_s; s <= end_s; s++) {

		int file_size = (int) pow(4, s);

		int TEST_CASES = (int) pow(2, (11-s));
		cout << TEST_CASES << endl;

		double c_enc_time = 0, s_enc_time = 0, c_decr_time = 0;

		vector<double> server_times;

		for (int j = 0; j < TEST_CASES; j++) {

			mpz_t n, g, cipher, decr;
			mpz_inits(n, g, cipher, decr, NULL);

			Client client(s, BIT_LENGTH, tree);
			client.get_pub_keys(n, g);

			Server server(BIT_LENGTH, file_size, tree, n, g);

			unsigned char s_bits[2*s];
			mpz_t *results = new mpz_t[3*s];

			for (int i = 0; i < 2*s; i++) {
				s_bits[i] = 1;
			}

			for (int i = 0; i < 3*s; i++) {
				mpz_init(results[i]);
			}

			c_enc_time += client.encrypt_s_bits(results, 3*s, s_bits, 2*s);

			double s_time = server.get_file(cipher, results, 1, 1);
			s_enc_time += s_time;
			server_times.push_back(s_time);

			c_decr_time += client.decr_file(decr, cipher);

			mpz_clears(n, g, cipher, decr, NULL);

			for (int i = 0; i < 3*s; i++) {
				mpz_clear(results[i]);
			}

			delete[] results;

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
}

void quad_serial(int start_s, int end_s) {

	TreeType tree = QUAD;

	cout << "**QUAD SEQUENTIAL**\n" << endl;

	for (int s = start_s; s <= end_s; s++) {

		int file_size = (int) pow(4, s);

		int TEST_CASES = (int) pow(2, (11-s));
		cout << TEST_CASES << endl;

		double c_enc_time = 0, s_enc_time = 0, c_decr_time = 0;

		vector<double> server_times;

		for (int j = 0; j < TEST_CASES; j++) {

			mpz_t n, g, cipher, decr;
			mpz_inits(n, g, cipher, decr, NULL);

			Client client(s, BIT_LENGTH, tree);
			client.get_pub_keys(n, g);

			Server server(BIT_LENGTH, file_size, tree, n, g);

			unsigned char s_bits[2*s];
			mpz_t *results = new mpz_t[3*s];

			for (int i = 0; i < 2*s; i++) {
				s_bits[i] = 1;
			}

			for (int i = 0; i < 3*s; i++) {
				mpz_init(results[i]);
			}

			c_enc_time += client.encrypt_s_bits(results, 3*s, s_bits, 2*s);

			double s_time = server.get_file(cipher, results, 0, 0);
			s_enc_time += s_time;
			server_times.push_back(s_time);

			c_decr_time += client.decr_file(decr, cipher);

			mpz_clears(n, g, cipher, decr, NULL);

			for (int i = 0; i < 3*s; i++) {
				mpz_clear(results[i]);
			}

			delete[] results;

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
}

/******* OCTOS *******/

void octo_new_p(int start_s, int end_s) {

	TreeType tree = OCTO;

	cout << "**OCTO NEW PARALLEL**\n" << endl;

	for (int s = start_s; s <= end_s; s++) {

		int file_size = (int) pow(8, s);

		int TEST_CASES = (int) pow(2, (11-s));

		cout << TEST_CASES << endl;

		vector<double> server_times;
		double c_enc_time = 0, s_enc_time = 0, c_decr_time = 0;

		for (int j = 0; j < TEST_CASES; j++) {

			mpz_t n, g, cipher, decr;
			mpz_inits(n, g, cipher, decr, NULL);

			Client client(s, BIT_LENGTH, tree);
			client.get_pub_keys(n, g);

			Server server(BIT_LENGTH, file_size, tree, n, g);

			unsigned char s_bits[3*s];
			mpz_t *results = new mpz_t[7*s];

			for (int i = 0; i < 3*s; i++) {
				s_bits[i] = 1;
			}

			for (int i = 0; i < 7*s; i++) {
				mpz_init(results[i]);
			}

			c_enc_time += client.encrypt_s_bits(results, 7*s, s_bits, 3*s);

			double s_time = server.get_file_new_p(cipher, results);
			s_enc_time += s_time;
			server_times.push_back(s_time);

			c_decr_time += client.decr_file(decr, cipher);

			mpz_clears(n, g, cipher, decr, NULL);

			for (int i = 0; i < 7*s; i++) {
				mpz_clear(results[i]);
			}

			delete[] results;

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

}

void octo_full_p(int start_s, int end_s) {

	TreeType tree = OCTO;

	cout << "**OCTO FULL PARALLEL**\n" << endl;

	for (int s = start_s; s <= end_s; s++) {

		int file_size = (int) pow(8, s);

		int TEST_CASES = (int) pow(2, (11-s));

		cout << TEST_CASES << endl;

		vector<double> server_times;
		double c_enc_time = 0, s_enc_time = 0, c_decr_time = 0;

		for (int j = 0; j < TEST_CASES; j++) {

			mpz_t n, g, cipher, decr;
			mpz_inits(n, g, cipher, decr, NULL);

			Client client(s, BIT_LENGTH, tree);
			client.get_pub_keys(n, g);

			Server server(BIT_LENGTH, file_size, tree, n, g);

			unsigned char s_bits[3*s];
			mpz_t *results = new mpz_t[7*s];

			for (int i = 0; i < 3*s; i++) {
				s_bits[i] = 1;
			}

			for (int i = 0; i < 7*s; i++) {
				mpz_init(results[i]);
			}

			c_enc_time += client.encrypt_s_bits(results, 7*s, s_bits, 3*s);

			double s_time = server.get_file(cipher, results, 1, 1);
			s_enc_time += s_time;
			server_times.push_back(s_time);

			c_decr_time += client.decr_file(decr, cipher);

			mpz_clears(n, g, cipher, decr, NULL);

			for (int i = 0; i < 7*s; i++) {
				mpz_clear(results[i]);
			}

			delete[] results;

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
}

void octo_serial(int start_s, int end_s) {

	TreeType tree = OCTO;

	cout << "**OCTO SEQUENTIAL**\n" << endl;

	for (int s = start_s; s <= end_s; s++) {

		int file_size = (int) pow(8, s);

		int TEST_CASES = (int) pow(2, (11-s));

		cout << TEST_CASES << endl;

		vector<double> server_times;
		double c_enc_time = 0, s_enc_time = 0, c_decr_time = 0;

		for (int j = 0; j < TEST_CASES; j++) {

			mpz_t n, g, cipher, decr;
			mpz_inits(n, g, cipher, decr, NULL);

			Client client(s, BIT_LENGTH, tree);
			client.get_pub_keys(n, g);

			Server server(BIT_LENGTH, file_size, tree, n, g);

			unsigned char s_bits[3*s];
			mpz_t *results = new mpz_t[7*s];

			for (int i = 0; i < 3*s; i++) {
				s_bits[i] = 1;
			}

			for (int i = 0; i < 7*s; i++) {
				mpz_init(results[i]);
			}

			c_enc_time += client.encrypt_s_bits(results, 7*s, s_bits, 3*s);

			double s_time = server.get_file(cipher, results, 0, 0);
			s_enc_time += s_time;
			server_times.push_back(s_time);

			c_decr_time += client.decr_file(decr, cipher);

			mpz_clears(n, g, cipher, decr, NULL);

			for (int i = 0; i < 7*s; i++) {
				mpz_clear(results[i]);
			}

			delete[] results;

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

}

/******* SCALABLE *******/

void octo_scalable(int subtree_start, int subtree_end, int s) {

	TreeType tree = OCTO;

	unsigned int file_size = (int) pow(8, s);

	cout << "SCALABLE OCTAL - " << file_size << " files\n" << endl;

	for (int i = subtree_start; i <= subtree_end; i++) {

		int subt_size = (int) pow(8, i);
		int subt_length = file_size / subt_size;

		vector<double> server_times;
		double c_enc_time = 0, s_enc_time = 0, c_decr_time = 0;

		for (int j = 0; j < TEST_CASES; j++) {

			int s_length = 3*i;
			int result_length = 7*i;

			int selected_subtree = 3; // bu demek ki 011 ile başlıyo selection

			unsigned char s_bits[s_length];

			for (int x = 0; x < s_length; x++) {
				s_bits[x] = 0; // 011000...
			}

			mpz_t results[result_length];

			for (int x=0; x < result_length; x++) {
				mpz_init(results[x]);
			}

			mpz_t encr_subtree[subt_length];

			mpz_t n, g, cipher, decr;
			mpz_inits(n, g, cipher, decr, NULL);

			Client client(i, BIT_LENGTH, tree, true);
			client.get_pub_keys(n, g);

			c_enc_time += client.get_scalable_s_bits(results, result_length, s_bits, s_length, encr_subtree, subt_length, selected_subtree);

			Server server(BIT_LENGTH, file_size, tree, n, g, subt_size);

			double s_time = server.get_file_scalable(cipher, results, encr_subtree, subt_length);
			s_enc_time += s_time;
			server_times.push_back(s_time);

			c_decr_time += client.decr_file(decr, cipher);

			//cout << "file: " << decr << endl;

			for (int x=0; x < result_length; x++) {
				mpz_clear(results[x]);
			}

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


}


int main(int argc, char *argv[]) {

	if(argc > 2) {
		CORE_SIZE = atoi(argv[1]);
		TEST_CASES = atoi(argv[2]);
	} else {
		cout << "Usage: ./main CORE_SIZE TEST_CASES" << endl;
		CORE_SIZE = 1;
		TEST_CASES = 1;
	}

	binary_new_p(1, 10);
	quad_new_p(1, 6);
	octo_new_p(1, 4);

	return 0;
}
