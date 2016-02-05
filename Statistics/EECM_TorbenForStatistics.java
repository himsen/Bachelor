/*
 * ECM using Edwards curves.
 * Implements the ECM algorithm using Edwards curves. Part of the authors bachelor thesis. 
 * Version 1.0.
 * Torben Hansen
 * Aarhus University, Institute of Mathematics.
 */

import java.io.*;
import java.math.*;
import java.util.*;

public class EECM_TorbenForStatistics{

	
	//Input 15 if running stats for the 15bit_sampleset.txt. 
	public static void main(String[] args){
		FileReader reader;
		BufferedReader bufReader;
		FileWriter writer, writerTime, writerMult, writerTotalOp, writerCurves;
		BufferedWriter bufWriter, bufWriterTime, bufWriterMult, bufWriterTotalOp, bufWriterCurves;
		try {
			//Setup writers
			writer = new FileWriter(args[0] + "bit_stats.txt");
			bufWriter = new BufferedWriter(writer);
			writerTime = new FileWriter(args[0] + "bit_time.txt");
			bufWriterTime = new BufferedWriter(writerTime);
			writerMult = new FileWriter(args[0] + "bit_mult.txt");
			bufWriterMult = new BufferedWriter(writerMult);
			writerTotalOp = new FileWriter(args[0] + "bit_totalOp.txt");
			bufWriterTotalOp = new BufferedWriter(writerTotalOp);
			writerCurves = new FileWriter(args[0] + "bit_curves.txt");
			bufWriterCurves = new BufferedWriter(writerCurves);
			
			
			//Stats
			double time = 0;
			BigInteger mult = BigInteger.ZERO;
			BigInteger totalOp = BigInteger.ZERO;
			//BigInteger curves = BigInteger.ZERO;
			double curves = 0;
			
			//Reader for reading file containg n.
			reader = new FileReader(args[0] + "bit_sampleset.txt");
			bufReader = new BufferedReader(reader);
			
			//Test
			bufWriter.write("Testing information: " + args[0] + " bit factor. Using manual bounds: No. Using options: -FS" + (Integer.parseInt(args[0])/3) + ", -NC2000");
			String[] input = new String[3];
			input[0] = "-NC2000";
			input[1] = "-FS" + (Integer.parseInt(args[0]))/3;
			int numberOfTimesRunning = 500;
			for(int j = 0; j < numberOfTimesRunning; j++){
				EECM_Torben ecm = new EECM_Torben();
				input[2] = bufReader.readLine();
				ecm.main(input);
				
				time += ecm.time;
				bufWriterTime.write(Double.toString(ecm.time));
				bufWriterTime.newLine();
				
				mult = mult.add(new BigInteger(Long.toString(ecm.M + ecm.S)));
				bufWriterMult.write(Long.toString(ecm.M + ecm.S));
				bufWriterMult.newLine();
				
				totalOp = totalOp.add(new BigInteger(Long.toString(ecm.M + ecm.S + ecm.A + ecm.I + ecm.D)));
				bufWriterTotalOp.write(Long.toString(ecm.M + ecm.S + ecm.A + ecm.I + ecm.D));
				bufWriterTotalOp.newLine();
				
				//curves = curves.add(new BigInteger(Integer.toString(ecm.numberOfCurvesUsed)));
				curves += ecm.numberOfCurvesUsed;
				bufWriterCurves.write(Integer.toString(ecm.numberOfCurvesUsed));
				bufWriterCurves.newLine();
				ecm.primes.clear();
				
			}
			bufWriter.newLine();
			bufWriter.write("Average time: " + time/numberOfTimesRunning);
			bufWriter.newLine();
			bufWriter.write("Average multiplications: " + mult.divide(new BigInteger(Integer.toString(numberOfTimesRunning))));
			bufWriter.newLine();
			bufWriter.write("Average total modular operations: " + totalOp.divide(new BigInteger(Integer.toString(numberOfTimesRunning))));
			bufWriter.newLine();
			//bufWriter.write("Avarages number of curves: " + curves.divide(new BigInteger(Integer.toString(numberOfTimesRunning))));
			bufWriter.write("Avarages number of curves: " + curves/numberOfTimesRunning);
			bufWriter.close();
			writer.close();
			bufReader.close();
			reader.close();
			bufWriterTime.close();
			writerTime.close();
			bufWriterMult.close();
			writerMult.close();
			bufWriterTotalOp.close();
			writerTotalOp.close();
			bufWriterCurves.close();
			writerCurves.close();
		} 
		catch (IOException e) {
			e.printStackTrace();
		} 		

	}
}