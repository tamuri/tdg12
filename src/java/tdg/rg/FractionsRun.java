package tdg.rg;

import java.io.*;
import java.util.*;
import pal.math.*;
import tdg.rg.Fractions;

public class FractionsRun {


	public static void main(String[] args) {
		FractionsRun fractionsRun = new FractionsRun();
		fractionsRun.run("");
		fractionsRun.run("simulation/1.swMutSel/");
		fractionsRun.run("simulation/2.swMutSel/");
		fractionsRun.run("simulation/3.swMutSel/");
	}

	void run(String prefix) {
		System.out.println(prefix);
		Fractions fractions = new Fractions(prefix);
	}
}

