import java.math.BigInteger;


public class Deneme {

	/**
	 * @param args
	 */
	public static void main(String[] args) {

		
		/*int total = 0, runs = 1000 ;
		
		for (int i = 0; i < runs; i++) {
			
			DamgardJurik dj = new DamgardJurik(512, 1);
			total += dj.keyGenerationTime;
		}
		
		System.out.println("In average: " + total/runs);*/

		DamgardJurik dj = new DamgardJurik(1024, 1);
		BigInteger b = new BigInteger("32323"), d;
		
		b = dj.encryption(b, BigInteger.valueOf(2).pow(1024));
		System.out.println("\n" + b + "\n");

		dj.setS(2);
		b = dj.encryption(b, BigInteger.valueOf(2).pow(1024));
		System.out.println("\n" + b + "\n");

		b = dj.decryption(b);
		System.out.println("first decr: " + b);
		dj.setS(1);
		System.out.println("\nsecond decr: " + dj.decryption(b));

		//
		//
		// dj.setS(2);
		// b = dj.encryption(b);
		//
		//
		// d = dj.decryption(b);
		// dj.setS(1);
		// d = dj.decryption(d);
		//
		// System.out.println(d);
		//
		//
		// dj.setS(3);
		// b = dj.encryption(b);
		//
		//
		// d = dj.decryption(b);
		// dj.setS(2);
		// d = dj.decryption(d);
		// dj.setS(1);
		// d = dj.decryption(d);
		//
		// System.out.println(d);
		//
		// dj.setS(4);
		// b = dj.encryption(b);
		//
		// d = dj.decryption(b);
		// dj.setS(3);
		// d = dj.decryption(d);
		// dj.setS(2);
		// d = dj.decryption(d);
		// dj.setS(1);
		// d = dj.decryption(d);
		//
		// System.out.println(d);
	}

}
