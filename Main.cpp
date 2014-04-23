#include "Client.h"
#include "Server.h"

using namespace std;

#define FILE_SIZE 64
#define BIT_LENGTH 1024

int main() {

	int s = (int) log2(FILE_SIZE);
//	int s = 7;
	TreeType tree = BINARY;
	mpz_t n, g, cipher, decr;
	mpz_inits(n, g, cipher, decr, NULL);

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

	client.encrypt_s_bits(results, s, s_bits, s);

	cout << "selection bits generated" << endl;

	server.get_file_it(cipher, results, s);

	cout << "server file encryption" << endl;

	client.decr_file(decr, cipher);

	cout << "decrypted: " << decr << endl;

	return 0;
}
