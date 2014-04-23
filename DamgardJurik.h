#ifndef DAMGARDJURIK_H
#define DAMGARDJURIK_H

#include <gmp.h>
#include <iostream>
#include <fcntl.h> 
#include <unistd.h>

enum TreeType {
	BINARY, QUAD, OCTO
};

class DamgardJurik {
	
private:
	
	mpz_t p, q, d, mu;
	gmp_randstate_t state;
	bool enc_only;
	
public:
	
	mpz_t n, n_s, n_sp, g;
	int bit_length, s;
	DamgardJurik(int bit_length, int s);
	DamgardJurik(int bit_length, int s, mpz_t n, mpz_t g);
	
	void keygen();
	
	void set_s(int s);

	void encrypt(mpz_t result, const mpz_t m);
	
	void find_i (mpz_t result, const mpz_t c);
	void decrypt(mpz_t result, const mpz_t c);
	
	void init_random();
	void get_random_prime(mpz_t result);
	
};

#endif
