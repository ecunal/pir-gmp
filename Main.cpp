#include "Client.h"
#include "Server.h"
#include <stdlib.h>
#include <vector>
#include <algorithm>

using namespace std;

#define BIT_LENGTH 1024

int CORE_SIZE;
int FILE_SIZE;

int main(int argc, char *argv[]) {

	if(argc > 2) {
		CORE_SIZE = atoi(argv[1]);
		FILE_SIZE = atoi(argv[2]);
	} else if (argc > 1) {
		CORE_SIZE = atoi(argv[1]);
	} else {
		cout << "Usage: ./main CORE_SIZE FILE_SIZE" << endl;
		CORE_SIZE = 1;
		FILE_SIZE = 64;
	}

	TreeType tree = OCTO;

	// scalable method denemee

	int file_size = 64;
	int subtree_size = 8;
	int subt_length = 8;

	int s_length = 3;
	int result_length = 7;

	int selected_subtree = 3; // bu demek ki 011 ile başlıyo selection

	unsigned char s_bits[s_length];

	s_bits[0] = 1;
	s_bits[1] = 0;
	s_bits[2] = 0; // yani 011100 = 28i istiyorum

	mpz_t results[result_length];

	for(int i=0; i<result_length; i++) {
		mpz_init(results[i]);
	}

	mpz_t encr_subtree[subt_length];

	mpz_t n, g, cipher, decr;
	mpz_inits(n, g, cipher, decr, NULL);

	Client client(1, BIT_LENGTH, tree, true);
	client.get_pub_keys(n, g);

	double client_enc_time = client.get_scalable_s_bits(results, result_length, s_bits, s_length, encr_subtree, subt_length, selected_subtree);

	Server server(BIT_LENGTH, file_size, tree, n, g, subtree_size);

	double server_time = server.get_file_scalable(cipher, results, encr_subtree, subt_length);

	double client_decr_time = client.decr_file(decr, cipher);

	cout << "Decrypted: " << decr << endl;

	cout << "client encryption (prl): " << client_enc_time << endl;
	cout << "server encryption (prl): " << server_time << endl;
	cout << "client decryption: " << client_decr_time << endl;

	mpz_clears(n, g, cipher, decr, NULL);


//	TreeType tree = QUAD;
//
//	for (int s = 2; s <= 5; s++) {
//
//		int file_size = (int) pow(4, s);
//
//		int TEST_CASES = 1; //(int) pow(2, (11-s));
//		cout << TEST_CASES << endl;
//
//		double c_enc_time = 0, s_enc_time = 0, c_decr_time = 0;
//
//		vector<double> server_times;
//
//		for (int j = 0; j < TEST_CASES; j++) {
//
//			mpz_t n, g, cipher, decr;
//			mpz_inits(n, g, cipher, decr, NULL);
//
//			Client client(s, BIT_LENGTH, tree);
//			client.get_pub_keys(n, g);
//
//			Server server(BIT_LENGTH, file_size, tree, n, g);
//
//			unsigned char s_bits[2*s];
//			mpz_t *results = new mpz_t[3*s];
//
//			for (int i = 0; i < 2*s; i++) {
//				s_bits[i] = 1;
//			}
//
//			for (int i = 0; i < 3*s; i++) {
//				mpz_init(results[i]);
//			}
//
//			c_enc_time += client.encrypt_s_bits(results, 3*s, s_bits, 2*s);
//
//			double s_time = server.get_file_new_p(cipher, results); // NEW METHOD!
//			s_enc_time += s_time;
//			server_times.push_back(s_time);
//
//			c_decr_time += client.decr_file(decr, cipher);
//
//			cout << "decr: " << decr << endl;
//
//			mpz_clears(n, g, cipher, decr, NULL);
//
//			for (int i = 0; i < 3*s; i++) {
//				mpz_clear(results[i]);
//			}
//
//			delete[] results;
//
//		}
//
//		sort(server_times.begin(), server_times.end());
//		double smallest = server_times.front();
//		double average = 0;
//		int count = 0;
//
//		cout << "Smallest time: " << smallest << endl;
//
//		for(vector<double>::iterator it = server_times.begin() ; it != server_times.end(); ++it ) {
//
//			if(*it < smallest * 1.1) {
//				average += *it;
//				count++;
//			} else {
//				cout << "Stop on: " << *it << "\nCount: " << count << "\nSize was: " << server_times.size() << endl;
//				break;
//			}
//		}
//
//		cout << "~~~~~~~~ results ~~~~~~~~~" << endl;
//		cout << "file size " << file_size << endl;
//		cout << "client encryption (prl): " << (c_enc_time / TEST_CASES) << endl;
//		cout << "server encryption (new method): " << (s_enc_time / TEST_CASES) << endl;
//		cout << "server encryption (new method - k best): " << (average / count) << endl;
//		cout << "client decryption: " << (c_decr_time / TEST_CASES) << "\n" << endl;
//		cout << "~~~~~~~~ end ~~~~~~~~~" << endl;
//	}
//
//	tree = OCTO;
//
//	for (int s = 2; s <= 2; s++) {
//
//		int file_size = (int) pow(8, s);
//
//		int TEST_CASES = 1; // (int) pow(2, (11-s));
//
//		cout << TEST_CASES << endl;
//
//		vector<double> server_times;
//		double c_enc_time = 0, s_enc_time = 0, c_decr_time = 0;
//
//		for (int j = 0; j < TEST_CASES; j++) {
//
//			mpz_t n, g, cipher, decr;
//			mpz_inits(n, g, cipher, decr, NULL);
//
//			Client client(s, BIT_LENGTH, tree);
//			client.get_pub_keys(n, g);
//
//			Server server(BIT_LENGTH, file_size, tree, n, g);
//
//			unsigned char s_bits[3*s];
//			mpz_t *results = new mpz_t[7*s];
//
//			for (int i = 0; i < 3*s; i++) {
//				s_bits[i] = 1;
//			}
//
//			for (int i = 0; i < 7*s; i++) {
//				mpz_init(results[i]);
//			}
//
//			c_enc_time += client.encrypt_s_bits(results, 7*s, s_bits, 3*s);
//
//			double s_time = server.get_file_new_p(cipher, results); // NEW METHOD!
//			s_enc_time += s_time;
//			server_times.push_back(s_time);
//
//			c_decr_time += client.decr_file(decr, cipher);
//
//			cout << "decr: " << decr << endl;
//
//			mpz_clears(n, g, cipher, decr, NULL);
//
//			for (int i = 0; i < 7*s; i++) {
//				mpz_clear(results[i]);
//			}
//
//			delete[] results;
//
//		}
//
//		sort(server_times.begin(), server_times.end());
//		double smallest = server_times.front();
//		double average = 0;
//		int count = 0;
//
//		cout << "Smallest time: " << smallest << endl;
//
//		for(vector<double>::iterator it = server_times.begin() ; it != server_times.end(); ++it ) {
//
//			if(*it < smallest * 1.1) {
//				average += *it;
//				count++;
//			} else {
//				cout << "Stop on: " << *it << "\nCount: " << count << "\nSize was: " << server_times.size() << endl;
//				break;
//			}
//		}
//
//		cout << "~~~~~~~~ results ~~~~~~~~~" << endl;
//		cout << "file size " << file_size << endl;
//		cout << "client encryption (prl): " << (c_enc_time / TEST_CASES) << endl;
//		cout << "server encryption (new method): " << (s_enc_time / TEST_CASES) << endl;
//		cout << "server encryption (new method - k best): " << (average / count) << endl;
//		cout << "client decryption: " << (c_decr_time / TEST_CASES) << "\n" << endl;
//		cout << "~~~~~~~~ end ~~~~~~~~~" << endl;
//	}
//
//
//
//	/**********/

	// KULLANCAZ BUNU SONRA






	/*******/









/*	for(int s=1; s<=1; s++) {

		int file_size = (int) pow(8, s);
		int TEST_CASES = 1;

		cout << "TEST_CASES: " << TEST_CASES << endl;

		double c_enc_time = 0, s_enc_time = 0, c_decr_time = 0;

		for (int j = 0; j < TEST_CASES; j++) {

			mpz_t n, g, cipher, decr;
			mpz_inits(n, g, cipher, decr, NULL);

			Client client(s, BIT_LENGTH, tree);
			client.get_pub_keys(n, g);

			Server server(BIT_LENGTH, s, file_size, tree, n, g);

			unsigned char s_bits[3*s];
			mpz_t *results = new mpz_t[7*s];

//			for (int i = 0; i < 3*s; i++) {
//				s_bits[i] = i;
//			}

			s_bits[0] = 0;
			s_bits[1] = 1;
			s_bits[2] = 1;

			for(int i=0; i<7*s; i++) {
				mpz_init(results[i]);
			}

			c_enc_time += client.encrypt_s_bits(results, 7*s, s_bits, 3*s);
			s_enc_time += server.get_file(cipher, results, s, 1, 1);
			c_decr_time += client.decr_file(decr, cipher);

			cout << "Decrypted: " << decr << endl;

			mpz_clears(n, g, cipher, decr, NULL);

			for (int i = 0; i < s; i++) {
				mpz_clear(results[i]);
			}

			delete[] results;

		}

		cout << "file size " << file_size << endl;
		cout << "client encryption (prl): " << (c_enc_time / TEST_CASES)
				<< endl;
		cout << "server encryption (prl): " << (s_enc_time / TEST_CASES)
				<< endl;
		cout << "client decryption: " << (c_decr_time / TEST_CASES) << "\n"
				<< endl;
	}

*/

/*	TreeType tree = BINARY;

		for (int s = 4; s <= 8; s++) {

			int file_size = (int) pow(2, s);

			int TEST_CASES = 5;

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
				s_enc_time += server.get_file_new_p(cipher, results, s); // 1: parallel
				c_decr_time += client.decr_file(decr, cipher);

				cout << decr << endl;

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


	double all_time = 0;

	 for(int i=0; i<100; i++) {

	 DamgardJurik dj(BIT_LENGTH, 6);

	 mpz_t exp1, exp2, msg, enc, dec, plain;
	 mpz_inits(exp1, exp2, enc, msg, dec, plain, NULL);

	 dj.init_random();
	 dj.get_random(msg, BIT_LENGTH);

	 cout << msg << endl;

	 double start_time = omp_get_wtime();

	 //#pragma omp parallel for
	 for(int i=0; i<2; i++) {
	 if(i == 0)
	 dj.encrypt_exp_1(exp1, msg);
	 else
	 dj.encrypt_exp_2(exp2);
	 }

	 dj.encrypt_mult(enc, exp1, exp2);

	 double end_time = omp_get_wtime();

	 all_time += end_time - start_time;

	 dj.decrypt(dec, enc);

	 cout << dec << endl;

	 }

	 cout << all_time / 100 << endl;

	 #pragma omp parallel for
	 for(int i=0; i<5; i++){

	 cout << i << " sth1" << endl;

	 #pragma omp task
	 {
	 cout << i << " sth2" << endl;
	 }
	 #pragma omp task
	 {
	 cout << i << " sth3" << endl;
	 }
	 #pragma omp taskwait

	 cout << i << " sth4" << endl;

	 }

	int s = (int) (log2(FILE_SIZE) / log2(8));
	TreeType tree = OCTO;
	mpz_t n, g, cipher1, cipher2, decr1, decr2;
	mpz_inits(n, g, cipher1, cipher2, decr1, decr2, NULL);

	Client client(s, BIT_LENGTH, tree);
	client.get_pub_keys(n, g);

	cout << "client created" << endl;

	Server server(BIT_LENGTH, s, FILE_SIZE, tree, n, g);

	cout << "server created" << endl;

	int s_l = 6, r_l = (s_l * 7) / 3;

	unsigned char s_bits[] = { 1, 1, 1, 0, 0, 0};
	mpz_t *results = new mpz_t[r_l];

	for (int i = 0; i < r_l; i++) {
		mpz_init(results[i]);
	}

	cout << "selection bits initalized" << endl;

	double client_time = 0;

	for(int i=0; i<20; i++) {

		client_time += client.encrypt_s_bits(results, r_l, s_bits, s_l);

	}

	cout << "selection bits generated: " << client_time/20 << endl;

	double time_1 = server.get_file(cipher1, results, s, 1, 1); // iterative
	cout << "octo seque: " << time_1 << endl;

//	double time_2 = server.get_file(cipher2, results, s, 1, true); // lvl 2
//	cout << "enc 2 lvl parallel: " << time_2 << endl;

	cout << "server files encrypted" << endl;

	client.decr_file(decr1, cipher1);
//	client.decr_file(decr2, cipher2);

	cout << "decrypted 1: " << decr1 << endl;
//	cout << "decrypted 2: " << decr2 << endl;

	mpz_clears(n, g, cipher1, cipher2, decr1, decr2, NULL); */

	return 0;
}
