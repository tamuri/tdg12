package tdg.rg;

import java.io.*;
import java.util.*;
import pal.math.*;
import tdg.utils.GeneticCode;

public class DistS implements MultivariateFunction {

    int sites = 759;
    //char[] aaCode = GeneticCode.VERTEBRATE_MITOCHONDRIAL_CODE.toCharArray();
    char[] aaCode = GeneticCode.STANDARD_CODE.toCharArray();


	Random randomGen = new Random();



	double[] histo = new double[600]; // histogram for real data
	double[][] kHisto = new double[sites][600]; // histogram for each location
	double[][] bStrapHisto = new double[1000][600]; // histogram for each of 1000 bootstraps
	double[] avg = new double[600]; // average of the bootstrapped histograms
	double[] avg2 = new double[600]; // average of the bootstrapped histograms squared
	double[][] results = new double[5][600];
	int iMax = 0;  // max range of histogram was 0
	int iMin = 600; // min range of histogram was 600
	double binSize = 0.25;
	int oSet = 300;
	boolean printStuff = true;
	boolean doBootStrap = true;
	int dataSet = -1;
	int leftBound = 230;
	int rightBound = 300;

	String functionType = "Reverse Gamma";

	String prefix = "";

	/*char[] aaCode = { 'F', 'F', 'L', 'L', 'S', 'S', 'S', 'S', 'Y', 'Y',
							'*', '*', 'C', 'C', 'W', 'W', 'L', 'L', 'L', 'L',
							'P', 'P', 'P', 'P', 'H', 'H', 'Q', 'Q', 'R', 'R',
							'R', 'R', 'I', 'I', 'M', 'M', 'T', 'T', 'T', 'T',
							'N', 'N', 'K', 'K', 'S', 'S', '*', '*', 'V', 'V',
							'V', 'V', 'A', 'A', 'A', 'A', 'D', 'D', 'E', 'E',
							'G', 'G', 'G', 'G'};*/

        // TODO don't forget about sites below!!
        //char[] aaCode = GeneticCode.STANDARD_CODE.toCharArray();

	public static void main(String[] args) {
		DistS distS = new DistS();
		distS.run();
	}

	int toBin(double x) {
		return Math.round(Math.round(x/binSize)+oSet);
	}

	double toX(int iBin) {
		return (iBin-oSet)*binSize;
	}

	DistS() {
		readFiles();
		if (doBootStrap) bootStrap();
		if (printStuff) printOutput();
	}

	DistS(String pref) {
		this.prefix = pref;
		readFiles();
		if (doBootStrap) bootStrap();
		if (printStuff) printOutput();
	}

	void run() {
	}

	void setDataSet(int dataSet) {
		this.dataSet = dataSet;
	}

	void setBounds(double[] bounds) {
		leftBound = toBin(bounds[0]);
		rightBound = toBin(bounds[1]);
	}

	double reverseGamma(double[] params, double x) {
			double negX = -1.0*x;
			return gamma(params, negX);
	}

	double gamma(double[] params, double x) {
			return params[0] * Math.pow(x, (params[1]-1.0)) * Math.exp(-x*params[2]);
	}

	double exponential(double[] params, double x) {
			return params[0] * Math.exp(-x*params[1]);
	}


	public double evaluate(double[] params) {
		double error2 = 0.0;
		for (int iS = leftBound; iS < rightBound; iS++) {
			double x = toX(iS);
			double funct = 0.0;
			if (functionType.equals("Reverse Gamma")) {
				funct = reverseGamma(params, x);
			} else if (functionType.equals("Gamma")) {
				funct = gamma(params, x);
			} else if (functionType.equals("Exponential")) {
				funct = exponential(params, x);
			}
			double delta = funct - histo[iS];
			if (dataSet >= 0) {
				delta = funct - bStrapHisto[dataSet][iS];
			}
			error2 += delta*delta/(avg2[iS]*avg2[iS]);
		}
		return error2;
	}

	public double getLowerBound(int n) {
		double[] lBounds = {0, 0, 0};
		return lBounds[n];
	}

	public double getUpperBound(int n) {
		double[] uBounds = {1.0, 10.0, 10.0};
		return uBounds[n];
	}

	public int getNumArguments() {
		if (functionType.equals("Exponential")) {
			return 2;
		} else {
			return 3;
		}
	}

	public OrthogonalHints getOrthogonalHints() {
		return null;
	}


	void printOutput() {
		for (int iS = iMin; iS < iMax; iS++) {
			System.out.println(iS + "\t" + toX(iS) + "\t" + histo[iS] + "\t" + avg[iS] + "\t" + avg2[iS]);
		}
	}

	void bootStrap() {
		// Construct various bootstraps
		for (int iBoot = 0; iBoot < 1000; iBoot++) {
			for (int iLocation = 0; iLocation < sites; iLocation++) {
				int iPoint = randomGen.nextInt(sites);
				for (int iS = 0; iS < 600; iS++) {
					bStrapHisto[iBoot][iS] += kHisto[iPoint][iS];
				}
			}
		}
		// Normalise all distributions
		double summ = 0.0;
		for (int iBoot = 0; iBoot < 1000; iBoot++) {
			summ = 0.0;
			for (int iS = 0; iS < 600; iS++) {
				summ += binSize * bStrapHisto[iBoot][iS];
			}
			for (int iS = 0; iS < 600; iS++) {
				bStrapHisto[iBoot][iS]/= summ;
				avg[iS] += bStrapHisto[iBoot][iS];
				avg2[iS] += bStrapHisto[iBoot][iS] * bStrapHisto[iBoot][iS];
			}
		}
		summ = 0.0;
		for (int iS = 0; iS < 600; iS++) {
			summ += binSize * histo[iS];
			avg[iS] /= 1000.0;
			avg2[iS] /= 1000.0;
			avg2[iS] = Math.sqrt(avg2[iS]-avg[iS]*avg[iS]);
		}
		for (int iS = 0; iS < 600; iS++) {
			histo[iS]/= summ;
			results[0][iS] = toX(iS);
			results[1][iS] = histo[iS];
			results[2][iS] = avg[iS];
			results[3][iS] = avg2[iS];
		}
	}


	void fillFit(double[] params) {
		for (int iS = leftBound; iS < rightBound; iS++) {
			double x = toX(iS);
			double funct = 0.0;
			if (functionType.equals("Reverse Gamma")) {
				funct = reverseGamma(params, x);
			} else if (functionType.equals("Gamma")) {
				funct = gamma(params, x);
			} else if (functionType.equals("Exponential")) {
				funct = exponential(params, x);
			}
			results[4][iS] = funct;
		}
	}

	void readFiles() {
		boolean ok = true;
		try {
			BufferedReader buffQ0 = new BufferedReader(new FileReader(prefix + "Q0.txt"));
			String[] stringQ0 = buffQ0.readLine().split("\\s");
			double[] q0 = new double[4096];
			for (int iPair = 0; iPair < 4096; iPair++) {
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
					for (int iPair = 0; iPair < 4096; iPair++) {
						int iCodon = (iPair - iPair%64)/64;
						char iRes = aaCode[iCodon];
						char jRes = aaCode[iPair%64];
						double piValue = Double.parseDouble(stringPi[iCodon]);
						double qValue = Double.parseDouble(stringQ[iPair]);
						if ( (iRes != jRes) ) {
							double deltaS = Double.parseDouble(stringS[iPair]);
                            deltaS = Math.min(10, deltaS);
                            deltaS = Math.max(-10, deltaS);
							int iBin = toBin(deltaS);
							histo[iBin] += piValue*q0[iPair];
							kHisto[iLine][iBin] += piValue*q0[iPair];
							iMax = Math.max(iMax, iBin);
							iMin = Math.min(iMin, iBin);
						}
					}
					iLine++;
				}
			}
			iMin = Math.max(0, iMin-1);
			iMax = Math.min(iMax+2, 600);
			buffS.close();
			buffQ.close();
			buffPi.close();
		} catch (IOException ioe) {
			System.out.println("Error -- " + ioe.toString());
			System.exit(1);
		}
	}


}