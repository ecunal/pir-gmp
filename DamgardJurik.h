#ifndef DAMGARDJURIK_H
#define DAMGARDJURIK_H

#include <gmp.h>
#include <iostream>
#include <fcntl.h> 
#include <unistd.h>
#include <inttypes.h>
#include <omp.h>

enum TreeType {
	BINARY, QUAD, OCTO
};

class DamgardJurik {
	
private:
	
	// Private keys
	mpz_t p, q, d, mu;
	gmp_randstate_t state;
	bool enc_only; 		// if true, the object does not have private keys
	
public:
	
	// Public keys
	mpz_t n, n_s, n_sp, g; // n_s = n^s, n_sp = n^(s+1)
	int bit_length, s;


	/** Constructors & Related methods **/

	// Normal constructor with decrypting abilities
	DamgardJurik(int bit_length, int s);

	// Server case - it only knows public keys
	DamgardJurik(int bit_length, int s, mpz_t n, mpz_t g);
	
	// Key generation - constructors call
	void keygen();
	
	// Also changes n_s, n_sp, and mu (if !enc_only)
	void set_s(int s);


	/** Encryption functions **/

	// Encrypts m and stores in result
	void encrypt(mpz_t result, const mpz_t m);

	// g^m part of the encryption
	void encrypt_exp_1(mpz_t result, const mpz_t m);

	// r^n_s part
	void encrypt_exp_2(mpz_t result);

	// multiply given m*g mod n_sp
	void encrypt_mult(mpz_t result, const mpz_t m, const mpz_t g);
	

	/** Decryption Methods**/

	void decrypt(mpz_t result, const mpz_t c);
	void find_i (mpz_t result, const mpz_t c);
	
	/** Utility functions **/

	void init_random();
	void get_random(mpz_t result, int bit_length);
	void get_random_prime(mpz_t result);
	
};

#endif
