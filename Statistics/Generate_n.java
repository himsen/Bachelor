import java.math.*;
import java.util.Random;
import java.io.*;

public class Generate_n {

	//Input the number of bits in smallest factor.
	public static void main(String[] args) {
		FileWriter writer;
		BufferedWriter bufWriter;
		try {
			Random random = new Random();
			writer = new FileWriter(args[0] + "bit_sampleset.txt");
			bufWriter = new BufferedWriter(writer);
			BigInteger p1,p2;
			int bitSize = Integer.parseInt(args[0]);
			int sizeN = 200;
			BigInteger n;
			for(int i = 0; i < 1000; i++){
				do{
					p1 = new BigInteger(bitSize, 2056, random);
				}
				while(!p1.isProbablePrime(2056));
				do{
					p2 = new BigInteger(sizeN-bitSize, 2056, random);
				}
				while(!p2.isProbablePrime(2056));
				n = p1.multiply(p2);
				bufWriter.write(n.toString());
				bufWriter.newLine();
			}
			bufWriter.close();
			writer.close();
		} 
		catch (IOException e) {
			e.printStackTrace();
		}

	}
}
