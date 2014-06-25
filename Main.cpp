#include "Client.h"
#include "Server.h"

using namespace std;

#define FILE_SIZE 64
#define BIT_LENGTH 1024

int main() {

/*	#pragma omp parallel for
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

	int s = (int) log2(FILE_SIZE);
//	int s = 7;
	TreeType tree = BINARY;
	mpz_t n, g, cipher1, cipher2, decr1, decr2;
	mpz_inits(n, g, cipher1, cipher2, decr1, decr2, NULL);

	Client client(s, BIT_LENGTH, tree);
	client.get_pub_keys(n, g);

	cout << "client created" << endl;

	Server server(BIT_LENGTH, s, FILE_SIZE, tree, n, g);

	cout << "server created" << endl;

	unsigned char s_bits[] = { 1, 0, 1, 1, 0, 1 };
	mpz_t *results = new mpz_t[s];

	for (int i = 0; i < s; i++) {
		mpz_init(results[i]);
	}

	cout << "selection bits initalized" << endl;


	double client_time = client.encrypt_s_bits(results, s, s_bits, s);
	cout << "selection bits generated: " << client_time << endl;


	server.get_file(cipher1, results, s, 0); // iterative
	server.get_file(cipher2, results, s, 1); // parallel

	cout << "server file encrypted" << endl;

	client.decr_file(decr1, cipher1);
	client.decr_file(decr2, cipher2);

	cout << "decrypted 1: " << decr1 << endl;
	cout << "decrypted 2: " << decr2 << endl;

	return 0;
}
