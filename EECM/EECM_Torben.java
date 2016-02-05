/*
 * ECM using Edwards curves.
 * Implements the ECM algorithm using Edwards curves. Part of the authors bachelor thesis. 
 * Version 1.0.
 * Torben Hansen
 * Aarhus University, Institute of Mathematics.
 */

import java.math.*;
import java.util.*;

public class EECM_Torben implements Runnable{

	/*
	 * TODO: (and ideas)
	 * - B1 and B2 are limited due to the implementation of Erasthosthenes sieve. 
	 * - The way we calculate the bounds we restrict ourself to pick optimal bounds up a 70-bit factor size. 
	 * - Hard coded curves.
	 * - Option to set optimal number of curves that should be used. 
	 * - Atkin-Morain and Montgomery construction
	 * - Implement a factoring strategy
	 */
	
	//Global fields
	static BigInteger n; //Number we wish to find a factor in. 
	static int numberOfThreads; //Number of threads running at the same time. Default 2 threads can run.
	static int numberOfThreadsRunning; //Number of threads running a.t.m.
	static int numberOfCurves; //Number of curves willing to try. Default the program exists when a total of Integer.MAX_VALUE curves have been used.
	static boolean isFactorFound; //Flag when a factor is found. 
	static BigInteger factorFound; //Hopefully in the end this will contain a factor of n.
	static int numberOfCurvesUsed; //Keep track of number of curves tried a.t.m.
	static long B1; //Bound for 1. stage. Default B1 = 10^5
	static long B2; //Bound for 2. stage. Default B2 = 10^6
	static boolean boundSet; //Flag if bound B1 or B2 is set manually. Later used to flag which stage found a factor.
	static ArrayList<Integer> primes; //Store primes needed for ECM stage 1 and 2.
	//static BigInteger trash1; //Trash BigInteger-register. 
	static final long max = new Long("1400000000"); //Max bound input.
	static int FS; //Factor size in digits (we search for a prime in n with SF digits). Manually set. A SF != 0 indicate it has been set. 
	static long intermediateTime1, intermidateTime2; //Holding intermediate times. 
	static long M, S, D, A, I; //Modular operations counts. 
	static Random random = new Random(); //random generator is used to pick a random Edwards curve with no predescribed properties.
	static int radix; //Base for representation of n
	static final byte version1 = 1, version2 = 0; //Version of program
	static int maxPrimeGap; //maximal gap between two consecutive primes below B2.
	static boolean Montgomery; //Flag if the Montgomery construction should be used.
	static boolean Atkin_Morain; //Flag if the Atkin_Morain construction should be used.
	//static double logSave = -1; // Save log n here
	static BigInteger luckyX, luckyY; //Initial point point for the curve that found a factor. 
	static BigInteger luckyD; //Curve that found a factor;
	static boolean trialDivide; //Flag if we found a factor while doing trial division.
	static double time; //Used to output time for statistics
	static final int[] smallPrimes = new int[]{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 
										 101, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 
										 211, 223, 227,229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 
										 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 
										 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 
										 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 
										 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 
										 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 
										 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 
										 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997}; //Used to test for very small prime factors. 
	
	
	public static void main(String[] args) {
		//Start time
		final long startTime = System.currentTimeMillis();
		
		//Set global fields;
		numberOfThreads = 2;
		numberOfThreadsRunning = 0;
		numberOfCurves = Integer.MAX_VALUE-1;
		isFactorFound = false;
		factorFound = BigInteger.ZERO;
		numberOfCurvesUsed = 0;
		B1 = 10000;
		B2 = 100000;
		boundSet = false;
		primes = new ArrayList<Integer>();
		FS = 0;
		intermediateTime1 = 0;
		intermidateTime2 = 0;
		M = 0;
		S = 0;
		D = 0;
		A = 0;
		I = 0;
		radix = 10;
		maxPrimeGap = 0;
		Montgomery = false;
		Atkin_Morain = false;
		luckyX = BigInteger.ZERO;
		luckyY = BigInteger.ZERO;
		luckyD = BigInteger.ZERO;
		trialDivide = false;
		time = 0;
		
		
		//--------------------------------Parse Input--------------------------------//
		//Ugly implementation.
		
		HashMap<String, Integer> inputOrganiser = new HashMap<String, Integer>();
		//Specify the 1. stage bound.
		inputOrganiser.put("-B1",1); 
		
		//Specify the 2. stage bound.
		inputOrganiser.put("-B2",2); 
		
		//Specifying how many threads run. This may depend on the hardware on which this algorithm is run.
		inputOrganiser.put("-NT",3); 
		
		//One may specify the maximal amount of curves that should be tried. 
		inputOrganiser.put("-NC", 4); 
		
		//Search factor. This will calculate a (heuristic) optimal value for B1 and B2. This value is restricted to and must be in radix 10. 
		//input: number of digits of the prime number we search for.
		//Only configure new bounds if the user has not manually pick he's own bounds. 
		inputOrganiser.put("-FS", 5); 
		
		//Specifying that the Montgomery construction should be used.
		inputOrganiser.put("-MO",6);
		
		//Specifying that the Atkin-Morain construction should be used.
		inputOrganiser.put("-AM",7);
		
		//Print an overview of the different options the user has to choose from.
		inputOrganiser.put("-HE", 8);
		
		String inputString;
		for(int i = 0; i < args.length; i++){
			inputString = args[i];
			if(inputString.substring(0,1).equals("-")){
				//inputString = -??
				switch(inputOrganiser.get(inputString.substring(0,3)).intValue()){
					case 1: B1 = Long.parseLong(inputString.substring(3,inputString.length()));
							if(B1 > max){
								System.out.println("Bound B1 should be smaller than " + max + ", and in radix 10. Exiting..");
								System.exit(1);
							}
							boundSet = true;
							break;
					case 2: B2 = Long.parseLong(inputString.substring(3,inputString.length()));
							if(B2 > max){
								System.out.println("Bound B2 should be smaller than " + max + " and in radix 10. Exiting..");
								System.exit(1);
							}
							boundSet = true;
							break;
					case 3: numberOfThreads = Integer.parseInt(inputString.substring(3,inputString.length()));
							break;
					case 4: numberOfCurves = Integer.parseInt(inputString.substring(3,inputString.length()));
							break;
					case 5: FS = Integer.parseInt(inputString.substring(3,inputString.length()));
							if(FS > 74){
								System.out.println("The factor size must be (stricly) below 75 digits. Exiting..");
								System.exit(1);
							}
							break;
					case 6: Montgomery = true;
							if(Atkin_Morain){
								System.out.println("You need to decided on either Atkin-Morain or Montgomery construction. Exiting..");
								System.exit(1);
							}
							break;
					case 7: Atkin_Morain = true;
							if(Montgomery){
								System.out.println("You need to decided on either Atkin-Morain or Montgomery construction. Exiting..");
								System.exit(1);
							}
							break;
					case 8: System.out.println("\nFormat: -option ; information.");
							System.out.println("\n-B1 ; Specify the B1 bound. Default is 10000.");
							System.out.println("\n-B2 ; Specify the B2 bound. Default is 100000.");
							System.out.println("\n-NT ; Specify the number of threads that should be used. Defualt is 2.");
							System.out.println("\n-NC ; Specify the number of curves that should be tried. Default is");
							System.out.println("      Integer.MAX_VALUE-1.");
							System.out.println("\n-FS ; Sets 'optimal' bounds for a search after the particular factor size in");
							System.out.println("      digits.");
							System.out.println("\n-HE ; Help.");
							System.exit(0);		
				}
			}
			else if(inputString.substring(0,1).equals("(")){
				//inputString = (??)xxxxxxxxxxxxxxxxxxxxxxxx  - (??) contains base information for the input number. 
				//n = inputString.substring(index+1, inputString.length())
				//radix = Integer.parseInt(inputString.substring(1, index))			
				int index = inputString.indexOf(")");
				radix = Integer.parseInt(inputString.substring(1, index));
				n = new BigInteger( inputString.substring(index+1, inputString.length()), Integer.parseInt(inputString.substring(1, index)) );
			}
			else{
				//inputString = xxxxxxxxxxxxxxxxxxxxxxx
				//n = inputString				
				n = new BigInteger(inputString); //Default radix = 10.
			}
		}
		
		System.out.println("\nECM software writing by Torben Hansen, version " + version1 + "." + version2 + ".");
		System.out.println("For an overview of options input option -HE");
		if(radix != 10){
			System.out.println("Input (base "+ radix + "): " + n.toString(radix));
		}
		System.out.println("Input (base 10): " + n);
		System.out.println("Is a " + n.bitLength() + " bit number.");
		
		
		//--------------------------------Initial work--------------------------------//
		intermediateTime1 = System.currentTimeMillis();
		System.out.println("\n----------Doing initial work----------");
		System.out.println("Trial divide....");
		if(!trialDivide()){
			//Check bounds or calculate if FS is set manually. 
			System.out.println("Configure bounds....");
			Bounds();
			
			//Find primes for 1. and 2. stage.
			System.out.println("Sieving for primes...."); 
			ErasthosthenesSieve();
			System.out.println("Found " + primes.size() + " primes.");
			
			//Compute maximal prime gap for primes below B2.
			System.out.println("Computing maximal primegap with respect to B2....");
			computeMaxPrimeGap(); 	
		}
	
		//Running time of this step
		intermidateTime2 = System.currentTimeMillis();
		System.out.println("Doing initial work took: " + (double)(intermidateTime2 - intermediateTime1)/1000 + "sec.");
		
		
		//--------------------------------ECM algorithm--------------------------------//
		if(!isFactorFound){
			intermediateTime1 = System.currentTimeMillis();
			System.out.println("\n----------ECM----------");
			
			//Output starting parameters
			System.out.println("ECM is initialized with: B1 = " + B1 + " and B2 = " + B2 + ".");
			System.out.println("Trying at most " + numberOfCurves + " curves, using " + numberOfThreads + " thread/s.");	
			
		
		//--------------------------------System controller--------------------------------//
			intermediateTime1 = System.currentTimeMillis();	
			System.out.println("\nRunning ECM algorithm....."); 
			
			boundSet = false; //Now used as a flag if a factor is found in stage 1.
			
			//Manage threads
			do{
				if(numberOfThreadsRunning < numberOfThreads){
					new Thread(new EECM_Torben()).start();	
					incrementNumberOfThreadsRunning();
					//numberOfThreadsRunning++;
					incrementNumberOfCurvesUsed();
					//numberOfCurvesUsed++; 
				}
			}
			while(!isFactorFound && (numberOfCurvesUsed < numberOfCurves+1));
			
			/*
			//Stall until threads are finished.
			String dot = ".";
			System.out.print("Waiting for remaining running threads to finish...."); 
			while(numberOfThreadsRunning > 0){
				//System.out.print("\b \b"); //Delete latest output
				//System.out.print(dot);
				try {
					Thread.sleep(3000);
				} 
				catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			System.out.println(); //New line.
			*/
			
			//Running time of this step
			intermidateTime2 = System.currentTimeMillis();
			System.out.println("Running ECM algorithm took: " + (double)(intermidateTime2 - intermediateTime1)/1000 + "sec.");
		}
	
		
		//--------------------------------Output result--------------------------------//
		System.out.println("\n----------Result----------");
		
		//Factor found?
		if(isFactorFound){
			if(trialDivide){
				System.out.println("Factor found doing trial division");
				System.out.println("Factor is: "+ factorFound + ". Is prime.");
				System.out.println("Is a " + factorFound.bitLength() + " bit number.");
				if(n.divide(factorFound).isProbablePrime(1028)){
					System.out.println("Cofactor is: " + n.divide(factorFound));
					System.out.println("Is probably a prime - with probability 1-1/(2^1028).");
					System.out.println("Is a " + n.divide(factorFound).bitLength() + " bit number.");
				}
				else{
					System.out.println("Cofactor is: " + n.divide(factorFound));
					System.out.println("Is not a prime.");
					System.out.println("Is a " + n.divide(factorFound).bitLength() + " bit number.");			
				}
			}
			else if(boundSet){ //BoundSet true if and only if factor was found in stage 1.
				System.out.println("Factor found in 1. stage");
				if(factorFound.isProbablePrime(1028)){
					System.out.println("Factor is: " + factorFound + ".");
					System.out.println("Is probably a prime - with probability 1-1/(2^1028).");
					System.out.println("Is a " + factorFound.bitLength() + " bit number.");
				}
				else{
					System.out.println("Factor is: " + factorFound + ".");
					System.out.println("Is not a prime.");
					System.out.println("Is a " + factorFound.bitLength() + " bit number.");
				}
				
				if(n.divide(factorFound).isProbablePrime(1028)){
					System.out.println("Cofactor is: " + n.divide(factorFound));
					System.out.println("Is probably a prime - with probability 1-1/(2^1028).");
					System.out.println("Is a " + n.divide(factorFound).bitLength() + " bit number.");
				}
				else{
					System.out.println("Cofactor is: " + n.divide(factorFound));
					System.out.println("Is not a prime.");
					System.out.println("Is a " + n.divide(factorFound).bitLength() + " bit number.");			
				}
			}
			else{
				System.out.println("Factor found in 2. stage");
				if(factorFound.isProbablePrime(1028)){
					System.out.println("Factor is: " + factorFound + ".");
					System.out.println("Is probably a prime - with probability 1-1/(2^1028).");
					System.out.println("Is a " + factorFound.bitLength() + " bit number.");
				}
				else{
					System.out.println("Factor is: " + factorFound + ".");
					System.out.println("Is not a prime.");
					System.out.println("Is a " + factorFound.bitLength() + " bit number.");
				}
				
				if(n.divide(factorFound).isProbablePrime(1028)){
					System.out.println("Cofactor is: " + n.divide(factorFound));
					System.out.println("Is probably a prime - with probability 1-1/(2^1028).");
					System.out.println("Is a " + n.divide(factorFound).bitLength() + " bit number.");
				}
				else{
					System.out.println("Cofactor is: " + n.divide(factorFound));
					System.out.println("Is not a prime.");
					System.out.println("Is a " + n.divide(factorFound).bitLength() + " bit number.");				
				}
			}
		}
		else{
			System.out.println("Failed to find a factor.");
			System.out.println("Try changing the bounds B1 and B2.");
			System.out.println("Don't give up!");
		}
		
		//Stop time
		final long endTime = System.currentTimeMillis();
		
		time = (double)(endTime - startTime)/1000;
		
		//Output statistic 
		System.out.println("\n----------Statistics----------");
		System.out.println("Running program took: " + time + "sec.");
		if(!trialDivide && isFactorFound){
			System.out.println("Number of curves tried: " + numberOfCurvesUsed);
			System.out.println("Number of multiplications: " + M);
			System.out.println("Number of squarings: " + S);
			System.out.println("Number of additions/subtractions: " + A);
			System.out.println("Number of inversions: " + I);
			System.out.println("Number of multiplications by d: " + D);
			System.out.println("A total of " + (M+S+A+D+I) + " modular operations.");
			if(luckyD.compareTo(BigInteger.ZERO) != 0 ){
				System.out.println("Curve: x^2+y^2 = 1+" + luckyD + "x^2y^2.");
				System.out.println("Initial point (x,y) = ( " + luckyX + " , " + luckyY + " ).");
			}
			else{
				System.out.println("Factor found during some setup.");
			}
		}
		System.out.print("\nWaiting for remaining running threads to finish....");
	}

	
	///////////////////////////////////////
	//                                   //
	//           ECM algorithm           //
	//                                   //
	///////////////////////////////////////
	
	public void run() {
		//--------------------------------Initialize curve--------------------------------//
		
		//Pick curve and point
		BigInteger temp1, temp2, temp3, temp4, x, y, d = BigInteger.ZERO;
		//if(Montgomery){
			
		//}
		//else if(Atkin_Morain){
			
		//}
		//else{

		do{
			x = (new BigInteger(n.bitLength(), random)).mod(n);
			y = (new BigInteger(n.bitLength(), random)).mod(n);
			temp1 = x.modPow(new BigInteger("2"), n);
			temp2 = y.modPow(new BigInteger("2"), n);
			temp3 = (temp1.add(temp2)).subtract(BigInteger.ONE).mod(n);
			temp4 = temp1.multiply(temp2).mod(n);
			BigInteger g = temp4.gcd(n);
			M += 1; S += 2; A += 2;
			if(g.compareTo(BigInteger.ONE) == 0 && !isFactorFound){
				temp4 = temp4.modInverse(n);
				d = temp3.multiply(temp4).mod(n);
				M += 1; I += 1;
			}	
			else if(!isFactorFound){
				//When n is big is it unlikely that g > 1.
				if(g.compareTo(n) != 0){ 
					factorIsFoundFirstStage(x, y, d, g);
					//DEBUG System.out.println("Found1");
					break;
				}
				else if(temp3.gcd(n).compareTo(n) != 0){
					factorIsFoundFirstStage(x, y, d, temp3.gcd(n));
					//DEBUG System.out.println("Found2");
					break;
				}
				else if(temp2.gcd(n).compareTo(n) != 0){
					factorIsFoundFirstStage(x, y, d, temp2.gcd(n));
					//DEBUG System.out.println("Found3");
					break;
				}
				else if(temp1.gcd(n).compareTo(n) != 0){
					factorIsFoundFirstStage(x, y, d, temp1.gcd(n));
					//DEBUG System.out.println("Found4");
					break;
				}
				else if(x.gcd(n).compareTo(n) != 0){
					factorIsFoundFirstStage(x, y, d, x.gcd(n));
					//DEBUG System.out.println("Found5");
					break;
				}
				else if(y.gcd(n).compareTo(n) != 0){
					factorIsFoundFirstStage(x, y, d, y.gcd(n));
					//DEBUG System.out.println("Found6");
					break;
				}
			}
		}
		while(d.compareTo(BigInteger.ONE) == 0 && !isFactorFound);
		//}
		
		Edward P = new Edward(x, y, BigInteger.ONE, d, n); //Construct curve.
		temp1 = x; //Used to find the initial point again.
		temp2 = y;
		

		
		//--------------------------------1. Stage--------------------------------//
		  		
		//Optimizations: Another prime batch strategy, a better single scalar multiplication scheme, store NAF representation for primes. 
		BigInteger z = BigInteger.ONE;
		int index2stage=0; //Used to index correctly into the prime array in the 2. stage.
		BigInteger gcdx;
		if(!isFactorFound){
			int prime = 0;
			byte[] NAF; 
			int length;
			
			prime = primes.get(0);
			for(int i = 1; prime < B1+1; i++){
				NAF = computeNAF(prime);
				if(NAF[NAF.length-1] == 1){
					length = NAF.length-1;
				}
				else{
					length = NAF.length-2;
				} 
				for(int j = 1 ; Math.pow(prime,j) < B1+1; j++){
					for(int s = length; s > 0; s--){
						P.dup();
						if(NAF[s-1] == 1){
							P.add(x, y, z);
						}
						else if(NAF[s-1] == -1){
							P.sub(x, y, z);
						}
					}
					x = P.x;
					y = P.y;
					z = P.z;
				}
				prime = primes.get(i);
				index2stage++;
			}
			
			//Check for a factor in n.
			gcdx = n.gcd(P.x);
			if(gcdx.compareTo(BigInteger.ONE) != 0 && gcdx.compareTo(n) != 0){
				if(!isFactorFound){
					factorIsFoundFirstStage(temp1, temp2, d, n.gcd(P.x));
					//DEBUG System.out.println("Found7");
				}
			}		
		}
	
		
		//--------------------------------2. Stage--------------------------------//
		
		//Optimizations:  
		//Calculation method is bullshit!
		//Can't this be done with a more efficient datastructure?
		//Can we speed up by normalizing something?
		if(!isFactorFound){
			//DEBUG System.out.println("2. stage"); 
			Edward[] Earray;
			Earray = new Edward[maxPrimeGap];
			
			//Calculate [2*s]P  for all s such that 2*s<=maxPrimeGap and store in array 'Earray'.
			Edward secStageDup = new Edward(P.x, P.y, P.z, d, n);
			Edward secStage = new Edward(P.x, P.y, P.z, d, n);
			secStage.dup();
			secStageDup.dup();
			Earray[0] = new Edward(secStageDup.x, secStageDup.y, secStageDup.z, d, n);

			for(int i = 1; i < Earray.length; i++){
				secStage.add(secStageDup.x, secStageDup.y, secStageDup.z);
				Earray[i] = new Edward(secStage.x, secStage.y, secStage.z, d, n) ;
			}
			
			//Update modular operation counters.
			M += secStage.M; S += secStage.S; D += secStage.D; A += secStage.A;
			M += secStageDup.M; S += secStageDup.S; D += secStageDup.D; A += secStageDup.A;
			
			//We need to calculate [q_0]P
			//Normalizing. primes.get(index2stage) is with high probability much larger than I/M hence may optimize by normalizing the point P=[x,y,z].
			//If z != 0 and gcd(z,n) = 1	
			BigInteger gcdz = z.gcd(n);
			if(gcdz.compareTo(BigInteger.ONE) == 0){
				if(z.compareTo(BigInteger.ZERO) != 0){
					z = z.modInverse(n);
					x = x.multiply(z);
					y = y.multiply(z);
					I += 1; M += 2;
					for(int i = 0; i < primes.get(index2stage); i++){
						P.addMixed(x, y);
					}
				}
			}
			else if(gcdz.compareTo(n) != 0){
				if(!isFactorFound){
					factorIsFoundSecondStage(temp1, temp2, d, n.gcd(P.x));
					//DEBUG System.out.println("Found8");
				}
			}
			
			//Check for at factor in n.
			gcdx = n.gcd(P.x);
			if(gcdx.compareTo(BigInteger.ONE) != 0 && gcdx.compareTo(n) != 0){
				if(!isFactorFound){
					factorIsFoundSecondStage(temp1, temp2, d, n.gcd(P.x));
					//DEBUG System.out.println("Found9");
				}
			}		
			
			//V_0=p_1-q_0, V_i=p_{i+1}-p_i
			//Next we calculate [q_0+V_0+V_1+..+V_i]P=[q_0+V_0+..V_{i-1}]P+[p_{i+1}-p_i]P for increasing i. 
			//[q_0+V_0+..V_{i-1}]P is known and [p_{i+1}-p_i]P=[2*s]P for some s.
			int difference = 0;
			for(int i = index2stage+1; i < primes.size(); i++){
				difference = (primes.get(i) - primes.get(i-1))/2; //Difference between to odd primes is always even.
				secStage = Earray[difference-1]; //Earray[i]= [2*(i+1)]P.  Earray[difference-1] = [2*difference]P = [(primes.get(i) - primes.get(i-1))]P
				P.add(secStage.x, secStage.y, secStage.z);

				//Check for at factor in n.
				gcdx = n.gcd(P.x);
				if(gcdx.compareTo(BigInteger.ONE) != 0 && gcdx.compareTo(n) != 0){
					if(!isFactorFound){
						factorIsFoundSecondStage(temp1, temp2, d, n.gcd(P.x));
						//DEBUG System.out.println("Found10");
					}
				}		
			}
		}
				
		//Update modular operation counters.
		M += P.M; S += P.S; D += P.D; A += P.A;
		//Flag that this thread is done. 
		//numberOfThreadsRunning--;
		decrementNumberOfThreadsRunning();
	}
	
	

	
	
	///////////////////////////////////////////
	//                                       //
	//           Auxiliary methods           //
	//                                       //
	///////////////////////////////////////////
	
	//Following synchronized methods secures that threads do not corrupt each other when 
	//setting global fields. 
	private synchronized void decrementNumberOfThreadsRunning(){
		numberOfThreadsRunning--;
	}
	
	private synchronized static void incrementNumberOfThreadsRunning(){
		numberOfThreadsRunning++;
	}

	private synchronized static void incrementNumberOfCurvesUsed(){
		numberOfCurvesUsed++;
	}
	
	private synchronized void factorIsFoundFirstStage(BigInteger x, BigInteger y, BigInteger d, BigInteger factor){
		isFactorFound = true;
		factorFound = factor;
		boundSet = true;
		luckyX = x;
		luckyY = y;
		luckyD = d;
	}
	
	private synchronized void factorIsFoundSecondStage(BigInteger x, BigInteger y, BigInteger d, BigInteger factor){
		isFactorFound = true;
		factorFound = factor;
		luckyX = x;
		luckyY = y;
		luckyD = d;
	}

	/*
	 * Compute bounds B1 and B2.
	 * If Either B1 or B2 are set manually these are checked; should have B1 < B2. 
	 * If B1 or B2 are not set manually they are computed using a heuristic scheme from
	 * http://www.loria.fr/~zimmerma/records/ecm/params.html and with B2~10*B1. If 10*B1 exceeds 
	 * 1'400'000'000 we set B2 = 1'400'000'000. 
	 */
	static void Bounds(){
		if(boundSet){
			if(B1-1 > B2){
				System.out.println("Bound B1 is greater or equal to B2. Exiting..");
				System.exit(1);
			}
		}
		else if(FS != 0){
			if(1 <= FS && FS < 5){
				B1 = 1000; B2 = 10000;
			}
			else if(FS < 10){
				B1 = 2500; B2 = 25000;
			}
			else if(FS < 15){
				B1 = 5000; B2 = 50000;
			}
			else if(FS < 20){
				B1 = 7500; B2 = 75000;
			}
			else if(FS < 25){
				B1 = 11000; B2 = 110000;
			}
			else if(FS < 30){
				B1 = 50000; B2 = 500000;
			}
			else if(FS < 35){
				B1 = 250000; B2 = 2500000;
			}
			else if(FS < 40){
				B1 = 1000000; B2 = 10000000;
			}
			else if(FS < 45){
				B1 = 3000000; B2 = 30000000;
			}
			else if(FS < 50){
				B1 = 11000000; B2 = 110000000;
			}
			else if(FS < 55){
				B1 = 43000000; B2 = 430000000;
			}
			else if(FS < 60){
				B1 = 110000000; B2 = 1100000000;
			}
			else if(FS < 65){
				B1 = 260000000; B2 = 1400000000; //Maximal B2 is 1'400'000'00 due to the implementation of the sieve being used. 
			}
			else if(FS < 70){
				B1 = 850000000; B2 = 1400000000;  
			}
			//else if(FS < 75){
			//	B1 = new Long("2900000000"); B2 = new Long("2600000000");
			//}
			else{
				System.out.println("EECM_Torben does not support that high a factor size.");
				System.out.println("Try inputting bounds manually.");
				System.out.println("Exiting..");
				System.exit(1);
			}
		}
	}
	
	/*
	 * Implementation of Erasthothenes sieve. 
	 * This method updates the array 'primes' with primes from 3 to B2. 
	 */
	static void ErasthosthenesSieve(){
		int sievingArraySize = (int)((B2+1)/2);
		boolean[] sievingArray = new boolean[sievingArraySize]; 

		int prime ;
		int loop = 1;
		primes.add(2);
		while(loop < (int)(Math.sqrt(B2)/2+1)){
			if(sievingArray[loop] == false){ //changed 2*loop+2 to prime+1
				prime = 2*loop+1;
				for(int i = loop*(prime+1); i < sievingArraySize; i += prime){ //(2*loop+1)^2=2(2*loop^2+2*loop)+1=2(loop*(2*loop+2))+1. Could use i = loop + prime which is easier to see; 3(2*loop+1)= 6*loop+3=6*loop+2+1=2(3*loop+1)+1=2(2*loop+1+loop)2*(prime+loop)+1
					sievingArray[i] = true;
				}
				primes.add(prime);
			}	
			loop++;	
		}
		
		while(loop < sievingArraySize){
			if(sievingArray[loop] == false){
				primes.add(2*loop+1);
			}
			loop++;
		}
	}
	
	//The NAF representation is stored as little-Endian.
	static byte[] computeNAF(int prime){
		//31-Integer.numberOfLeadingZeroes(prime) is a hack to calculate log_2(n).
		byte [] NAF = new byte[(31-Integer.numberOfLeadingZeros(prime))+2];
		byte temp = 0;
		while(prime > 0){
			if(prime % 2 == 1){
				NAF[temp] = (byte)(2 - (prime % 4));
				prime = prime - NAF[temp];
			}
			else{
				NAF[temp] = 0;
			}
			prime = prime/2;
			temp++;			
		}
		return NAF;
	}
	
	//Compute the maximal prime gap for primes below B2
	//See http://users.cybercity.dk/~dsl522332/math/primegaps/maximal.htm
	static void computeMaxPrimeGap(){
		if(B2 < 15684){
			maxPrimeGap = 22;
		}
		else if(B2 < 19610){
			maxPrimeGap = 26;
		}
		else if(B2 < 31398){
			maxPrimeGap = 36;
		}
		else if(B2 < 155922){
			maxPrimeGap = 43;
		}
		else if(B2 < 360654){
			maxPrimeGap = 48;
		}
		else if(B2 < 370262){
			maxPrimeGap = 56;
		}
		else if(B2 < 492113){
			maxPrimeGap = 57;
		}
		else if(B2 < 1349534){
			maxPrimeGap = 59;
		}
		else if(B2 < 1357202){
			maxPrimeGap = 66;
		}
		else if(B2 < 2010734){
			maxPrimeGap = 74;
		}
		else if(B2 < 4652354){
			maxPrimeGap = 77;
		}
		else if(B2 < 17051708){
			maxPrimeGap = 90;
		}
		else if(B2 < 20831324){
			maxPrimeGap = 105;
		}
		else if(B2 < 47326694){
			maxPrimeGap = 110;
		}
		else if(B2 < 122164748){
			maxPrimeGap = 111;
		}
		else if(B2 < 189695660){
			maxPrimeGap = 117;
		}
		else if(B2 < 191912784){
			maxPrimeGap = 124;
		}
		else if(B2 < 387096134){
			maxPrimeGap = 150;
		}
		else if(B2 < 436273010){
			maxPrimeGap = 141;
		}
		else if(B2 < 1294268492){
			maxPrimeGap = 144;
		}
		else{
			maxPrimeGap = 160;
		}
	//	else if(B2 < 1453168142){ //Maximal B2 is 1'400'000'00 due to the implementation of the sieve being used. 
	//		maxPrimeGap = 146;
	//	}
	//	else if(B2 < Long.parseLong("2300942550")){
	//		maxPrimeGap = 160;
	//	}
	//	else if(B2 < Long.parseLong("3842610774")){
	//		maxPrimeGap = 168;
	//	}
	//	else{ //B2 < Long.parseLong("4302407359")
	//		maxPrimeGap = 177;
	//	}
	}
	
	static boolean trialDivide(){
		BigInteger gcd;
		for(int i = 0; i < smallPrimes.length; i++){
			gcd = n.gcd(new BigInteger(Integer.toString(smallPrimes[i])));
			if(gcd.compareTo(BigInteger.ONE) != 0){
				if(gcd.compareTo(n) != 0){
					isFactorFound = true;
					trialDivide = true;
					factorFound = new BigInteger(Integer.toString(smallPrimes[i]));
					return true;
				}
			}
		}
		return false;
	}
	
	/*
	 * Checks whether n is a perfect power. If 'yes' base contain b such that b^a=n.
	 * Is to slow for large n!
	 */
	/*
	static boolean checkPerfectPower(){
		BigInteger aSquared;
		do
		{
			BigInteger result;

			int power = Math.max((int) (log()/log(base) - 2),1);
			int comparison;
			
			do
			{
				power++;
				result = base.pow(power);
				comparison = n.compareTo(result);
			}
			while( comparison > 0 && power < Integer.MAX_VALUE );
			
			if( comparison == 0 )
			{
				factorFound = base;
				isFactorFound = true;
				return true;
			}
			base = base.add(BigInteger.ONE);
			aSquared = base.pow(2);
		}
		while (aSquared.compareTo(n) <= 0);		
		return false;
	}
	*/

	
	/*
	 * Compute log_2(n)
	 */
	/*
	static double log()
	{
		if ( logSave != -1 )
			return logSave;

		BigInteger b;

		int temp = n.bitLength() - 1000;
		if (temp > 0) 
		{
			b=n.shiftRight(temp); 
			logSave = (Math.log(b.doubleValue()) + temp)*Math.log(2);
		}
		else 
			logSave = (Math.log(n.doubleValue()))*Math.log(2);

		return logSave;
	}
	*/
	/*
	 * Compute log_2(n) with a parameter.
	 */
	/*
	static double log(BigInteger x)
	{
		BigInteger b;
		
	    int temp = x.bitLength() - 1000;
	    if (temp > 0) 
	    {
	    	b=x.shiftRight(temp); 
	        return (Math.log(b.doubleValue()) + temp)*Math.log(2);
	    }
	    else 
	    	return (Math.log(x.doubleValue())*Math.log(2));
	}
	*/
}