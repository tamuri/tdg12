package tdg.rg;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Random;

public class Fractions {

	Random randomGen = new Random();

	double pd;
	double pn;
	double pa;
	double[] pdB = new double[1000];
	double[] pnB = new double[1000];
	double[] paB = new double[1000];
	double[] pdK = new double[3598];
	double[] pnK = new double[3598];
	double[] paK = new double[3598];
	boolean printStuff = true;
	boolean doBootStrap = true;

	String prefix = "";


	char[] aaCode = {	'F', 'F', 'L', 'L', 'S', 'S', 'S', 'S', 'Y', 'Y', 
					 		'C', 'C', 'W', 'W', 'L', 'L', 'L', 'L', 
							'P', 'P', 'P', 'P', 'H', 'H', 'Q', 'Q', 'R', 'R', 
							'R', 'R', 'I', 'I', 'M', 'M', 'T', 'T', 'T', 'T', 
							'N', 'N', 'K', 'K', 'S', 'S', 'V', 'V', 
							'V', 'V', 'A', 'A', 'A', 'A', 'D', 'D', 'E', 'E', 
							'G', 'G', 'G', 'G'};


	public static void main(String[] args) {
		Fractions fractions = new Fractions();
		fractions.run();
	}

	Fractions() {
		readFiles();
		if (doBootStrap) bootStrap();
		if (printStuff) printOutput();
	}

	Fractions(String pref) {
		this.prefix = pref;
		readFiles();
		if (doBootStrap) bootStrap();
		if (printStuff) printOutput();
	}

	void run() {
	}

	void printOutput() {
		Arrays.sort(paB);
		Arrays.sort(pnB);
		Arrays.sort(pdB);
		System.out.println(pd + "\t" + pdB[25] + "\t" + pdB[975]);
		System.out.println(pn + "\t" + pnB[25] + "\t" + pnB[975]);
		System.out.println(pa + "\t" + paB[25] + "\t" + paB[975]);
	}

	void bootStrap() {
		// Construct various bootstraps
		for (int iBoot = 0; iBoot < 1000; iBoot++) {
			for (int iLocation = 0; iLocation < 3598; iLocation++) {
				int iPoint = randomGen.nextInt(3598);
				paB[iBoot] += paK[iPoint];
				pnB[iBoot] += pnK[iPoint];
				pdB[iBoot] += pdK[iPoint];
			}
		}
		// Normalise all distributions
		double summ = pa + pn + pd;
		pa /= summ;
		pn /= summ;
		pd /= summ;
		for (int iBoot = 0; iBoot < 1000; iBoot++) {
			summ = paB[iBoot] + pnB[iBoot] + pdB[iBoot];
			paB[iBoot] /= summ;
			pnB[iBoot] /= summ;
			pdB[iBoot] /= summ;
		}
	}

	void readFiles() {
		boolean ok = true;
		try {
			BufferedReader buffQ0 = new BufferedReader(new FileReader(prefix + "Q0.txt"));
			String[] stringQ0 = buffQ0.readLine().split("\\s");
			double[] q0 = new double[3600];
			for (int iPair = 0; iPair < 3600; iPair++) {
				q0[iPair] = Double.parseDouble(stringQ0[iPair]);
			}
			buffQ0.close();

			BufferedReader buffS = new BufferedReader(new FileReader(prefix + "S.txt"));
			BufferedReader buffQ = new BufferedReader(new FileReader(prefix + "QS.txt"));
			BufferedReader buffPi = new BufferedReader(new FileReader(prefix + "PiS.txt"));
			boolean eof = false;
			int iLine = 0;
			while (!eof) {
				String lineS = buffS.readLine();
				String lineQ = buffQ.readLine();
				String linePi = buffPi.readLine();
				if (lineS == null) {
					eof = true;
				} else {
					String[] stringS = lineS.split("\\s");
					String[] stringQ = lineQ.split("\\s");
					String[] stringPi = linePi.split("\\s");
					for (int iPair = 0; iPair < 3600; iPair++) {
						int iCodon = (iPair - iPair%60)/60;
						int jCodon = iPair%60;
						char iRes = aaCode[iCodon];
						char jRes = aaCode[jCodon];
						double piValue = Double.parseDouble(stringPi[iCodon]);
						double qValue = Double.parseDouble(stringQ[iPair]);
						if ( (iRes != jRes) ) {
//						if ( (iCodon != jCodon) ) {
							double deltaS = Double.parseDouble(stringS[iPair]);
							double q = q0[iPair];
							q = qValue;
							if (deltaS > 2.0) {
								pa += piValue*q;
								paK[iLine] += piValue*q;
							} else if (deltaS > -10.0) {
								pn += piValue*q;
								pnK[iLine] += piValue*q;
							} else {
								pd += piValue*q;
								pdK[iLine] += piValue*q;
							}
						}
					}
					iLine++;
				}
			}
			buffS.close();
			buffQ.close();
			buffPi.close();
		} catch (IOException ioe) {
			System.out.println("Error -- " + ioe.toString());
			System.exit(1);
		}
	}


}
