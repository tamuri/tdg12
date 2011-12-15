package tdg.results;

import pal.math.ConjugateDirectionSearch;

import java.util.Arrays;

public class CurveFit {


	public static void main(String[] args) {
		CurveFit curveFit = new CurveFit();
		double[][] realResults = curveFit.run("");

		if (false) System.exit(1);

		for (int iX = 0; iX < 600; iX++) {
			System.out.println(iX + "\t" + realResults[0][iX] + "\t" + realResults[1][iX] + "\t" 
				+ realResults[2][iX] + "\t" + realResults[3][iX] + "\t" + realResults[4][iX]);
		}

		if (true) System.exit(1);

		double[][] sim1Results = curveFit.run("simulation/1.swMutSel/");
		double[][] sim2Results = curveFit.run("simulation/2.swMutSel/");
		double[][] sim3Results = curveFit.run("simulation/3.swMutSel/");
		System.out.println("\n\nResults");
		for (int iX = 0; iX < 600; iX++) {
			System.out.println(iX + "\t" + realResults[0][iX] + "\t" + realResults[1][iX] + "\t" 
//				+ realResults[2][iX] + "\t" + realResults[3][iX] + "\t" + realResults[4][iX]);
				+ realResults[2][iX] + "\t" + realResults[3][iX] + "\t" + realResults[4][iX] + "\t"
				+ sim1Results[1][iX] + "\t" + sim2Results[1][iX] + "\t" + sim3Results[1][iX]);
		}
	}

	double[][] run(String prefix) {

		System.out.println(prefix);
		DistS distS = new DistS(prefix);
		double[] bounds = {-7.0, -2.0};
		distS.setBounds(bounds);
		distS.functionType = "Reverse Gamma";

		ConjugateDirectionSearch cds = new ConjugateDirectionSearch();
		double[] xVector = {0.01, 1.0, 1.0};
		double tolfx = 0.00001;
		double tolx = 0.00001;
		cds.optimize(distS, xVector, tolfx, tolx);
		double realParam0 = xVector[0];
		double realParam1 = xVector[1];
		double realParam2 = xVector[2];

		distS.fillFit(xVector);
		double[][] results = distS.results;

		double[] simParam0 = new double[1000];
		double[] simParam1 = new double[1000];
		double[] simParam2 = new double[1000];
		for (int iDataBase = 0; iDataBase < 1000; iDataBase++) {
			distS.setDataSet(iDataBase);
			cds.optimize(distS, xVector, tolfx, tolx);
			simParam0[iDataBase] = xVector[0];
			simParam1[iDataBase] = xVector[1];
			simParam2[iDataBase] = xVector[2];
		}
		Arrays.sort(simParam0);
		Arrays.sort(simParam1);
		Arrays.sort(simParam2);
		System.out.println("Reverse Gamma");
		System.out.println(0 + "\t" + realParam0 + "\t" + simParam0[25] + "\t" + simParam0[974]);
		System.out.println(1 + "\t" + realParam1 + "\t" + simParam1[25] + "\t" + simParam1[974]);
		System.out.println(2 + "\t" + realParam2 + "\t" + simParam2[25] + "\t" + simParam2[974]);


		bounds[0] = 0.0;
		bounds[1] = 5.0;
		distS.setBounds(bounds);
		distS.functionType = "Exponential";

		cds = new ConjugateDirectionSearch();
		xVector[0] = 0.01;
		xVector[1] = 1.0;
		tolfx = 0.00001;
		tolx = 0.00001;
		cds.optimize(distS, xVector, tolfx, tolx);
		realParam0 = xVector[0];
		realParam1 = xVector[1];

		distS.fillFit(xVector);

		simParam0 = new double[1000];
		simParam1 = new double[1000];
		for (int iDataBase = 0; iDataBase < 1000; iDataBase++) {
			distS.setDataSet(iDataBase);
			cds.optimize(distS, xVector, tolfx, tolx);
			simParam0[iDataBase] = xVector[0];
			simParam1[iDataBase] = xVector[1];
		}
		Arrays.sort(simParam0);
		Arrays.sort(simParam1);
		System.out.println("Exponential");
		System.out.println(0 + "\t" + realParam0 + "\t" + simParam0[25] + "\t" + simParam0[974]);
		System.out.println(1 + "\t" + realParam1 + "\t" + simParam1[25] + "\t" + simParam1[974]);

		return distS.results;

	}
}