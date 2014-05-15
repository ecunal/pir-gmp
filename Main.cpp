#include "Client.h"
#include "Server.h"

using namespace std;

#define FILE_SIZE 64
#define BIT_LENGTH 1024

int main() {

	/*	double all_time = 0;

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

	 } */

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

	mpz_clears(n, g, cipher1, cipher2, decr1, decr2, NULL);

	return 0;
}
