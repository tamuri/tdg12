package tdg.results;

import java.io.*;

public class Convert {

	char[] aaCode = {	'F', 'F', 'L', 'L', 'S', 'S', 'S', 'S', 'Y', 'Y', 
					 		'*', '*', 'C', 'C', 'W', 'W', 'L', 'L', 'L', 'L', 
							'P', 'P', 'P', 'P', 'H', 'H', 'Q', 'Q', 'R', 'R', 
							'R', 'R', 'I', 'I', 'M', 'M', 'T', 'T', 'T', 'T', 
							'N', 'N', 'K', 'K', 'S', 'S', '*', '*', 'V', 'V', 
							'V', 'V', 'A', 'A', 'A', 'A', 'D', 'D', 'E', 'E', 
							'G', 'G', 'G', 'G'};


	public static void main(String[] args) {
		Convert convert = new Convert();
	}

	Convert() {
		readFiles();
	}


	void readFiles() {
		boolean ok = true;
		try {
			BufferedReader buffS = new BufferedReader(new FileReader("S.txt"));
			BufferedReader buffQ = new BufferedReader(new FileReader("QS.txt"));
			BufferedReader buffPi = new BufferedReader(new FileReader("PiS.txt"));
			PrintWriter newS = new PrintWriter(new BufferedWriter(new FileWriter("S.new")));
			PrintWriter newQ = new PrintWriter(new BufferedWriter(new FileWriter("QS.new")));
			PrintWriter newPi = new PrintWriter(new BufferedWriter(new FileWriter("PiS.new")));

			boolean eof = false;
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
						if ( (iRes != '*') && (jRes != '*') ) {
							newS.print(stringS[iPair] + " ");
							newQ.print(stringQ[iPair] + " ");
						}
					}
					for (int iRes = 0; iRes < 64; iRes++) {
						char jRes = aaCode[iRes];
						if ( (jRes != '*') ) {
							newPi.print(stringPi[iRes] + " ");
						}
					}
					newS.println();
					newQ.println();
					newPi.println();
				}
			}
			buffS.close();
			buffQ.close();
			buffPi.close();
			newS.close();
			newPi.close();
			newQ.close();
		} catch (IOException ioe) {
			System.out.println("Error -- " + ioe.toString());
			System.exit(1);
		}
	}
	

}
