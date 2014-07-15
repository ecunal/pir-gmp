#include "Client.h"
#include "Server.h"

using namespace std;

#define BIT_LENGTH 1024


void binary_example() {

	cout << "\nIterative 128 file binary free case" << endl;

	TreeType tree = BINARY;
	int file_size = 128; // should be a power of 2
	int s = (int) log2(file_size); // maximum s value

	mpz_t n, g, cipher, decr;
	mpz_inits(n, g, cipher, decr, NULL);

	// Initialize the client, get the public keys
	Client client(s, BIT_LENGTH, tree);
	client.get_pub_keys(n, g);

	// Initialize the server using that public keys
	Server server(BIT_LENGTH, s, file_size, tree, n, g);

	unsigned char s_bits[s];
	mpz_t *results = new mpz_t[s];

	for (int i = 0; i < s; i++) {
		s_bits[i] = 1; // selection bits are all 1, but this can obviously change
		mpz_init(results[i]);
	}

	// Client encrypts selection bits
	double c_enc_time = client.encrypt_s_bits(results, s, s_bits, s);

	// Server returns "cipher" using those results.
	double s_enc_time = server.get_file(cipher, results, s, 0, 0); // 0, 0: no parallelization

	// Client decrypts cipher.
	double c_decr_time = client.decr_file(decr, cipher);

	cout << "Got file : " << decr << endl;
	cout << "Client encryption time: " << c_enc_time << " sec" << endl;
	cout << "Server computation time: " << s_enc_time << " sec" << endl;
	cout << "Client decryption time: " << c_decr_time << " sec" << endl;

	mpz_clears(n, g, cipher, decr, NULL);

	for (int i = 0; i < s; i++) {
		mpz_clear(results[i]);
	}

}

void quad_example() {

	cout << "\nAll parallel 64 file quadratic tree case" << endl;

	TreeType tree = QUAD;
	int file_size = 64; // should be power of 4 for quad
	int s = (int) (log2(file_size) / 2); // maximum s value

	mpz_t n, g, cipher, decr;
	mpz_inits(n, g, cipher, decr, NULL);

	// Initialize the client, get the public keys
	Client client(s, BIT_LENGTH, tree);
	client.get_pub_keys(n, g);

	// Initialize the server using that public keys
	Server server(BIT_LENGTH, s, file_size, tree, n, g);

	unsigned char s_bits[2*s];
	mpz_t *results = new mpz_t[3*s];

	for (int i = 0; i < 2*s; i++) {
		s_bits[i] = 1;	// selection bits are again all 1
	}

	for (int i = 0; i < 3*s; i++) {
		mpz_init(results[i]);
	}

	// Client encrypts selection bits
	double c_enc_time = client.encrypt_s_bits(results, 3*s, s_bits, 2*s);

	// Server returns "cipher" using those results.
	double s_enc_time = server.get_file(cipher, results, s, 1, 1); // 1, 1: all parallelized

	// Client decrypts cipher.
	double c_decr_time = client.decr_file(decr, cipher);

	cout << "Got file : " << decr << endl;
	cout << "Client encryption time: " << c_enc_time << " sec" << endl;
	cout << "Server computation time: " << s_enc_time << " sec" << endl;
	cout << "Client decryption time: " << c_decr_time << " sec" << endl;

	mpz_clears(n, g, cipher, decr, NULL);

	for (int i = 0; i < 3*s; i++) {
		mpz_clear(results[i]);
	}
}

void octal_example() {

	cout << "\n512 file octal tree case with new method (split tree into cores)" << endl;

	TreeType tree = OCTO;
	int file_size = 512; // should be power of 8 for octal
	int s = (int) (log2(file_size) / 3); // maximum s value

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

	// Client encrypts selection bits
	double c_enc_time = client.encrypt_s_bits(results, 7*s, s_bits, 3*s);

	// Server returns "cipher" using those results
	double s_enc_time = server.get_file_new_p(cipher, results, s); // new method

	// Client decrypts cipher
	double c_decr_time = client.decr_file(decr, cipher);

	cout << "Got file : " << decr << endl;
	cout << "Client encryption time: " << c_enc_time << " sec" << endl;
	cout << "Server computation time: " << s_enc_time << " sec" << endl;
	cout << "Client decryption time: " << c_decr_time << " sec" << endl;

	mpz_clears(n, g, cipher, decr, NULL);

	for (int i = 0; i < 7*s; i++) {
		mpz_clear(results[i]);
	}

}

int main() {

	/*** Usage examples ***/

	/** Iterative 128 file binary tree **/

	binary_example();

	/** All parallel 64 file quadratic tree **/

	quad_example();

	/** 512 file octal tree with new method (split tree into cores)**/

	octal_example();

	return 0;
}
