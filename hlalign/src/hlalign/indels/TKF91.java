package hlalign.indels;

import java.util.ArrayList;

import hlalign.base.IndelModel;
import hlalign.base.Utils;

public class TKF91 extends IndelModel {
	
	public TKF91(double par[]){
		/** lambda, mu */
		params = par.clone(); 
		
		/** Start, Match, Delete, Insert, End */
		int[][] emit = new int[][] { {0, 0},
							 {1, 1},
							 {1, 0},
							 {0, 1},
							 {0, 0} };
		
		two = new TKF91HMM(emit);
	
		
		/* state definitions

		0  Start
		1  x y z
		2  x - z
		3  x y -
		4  x - - 
		5  - y - (after z)
		6  - y - (after gap in z)
		7  - - z
		8  End
		*/
		emit = new int[][] { {0, 0},
							 {1, 1},
							 {0, 1},
							 {1, 0},
							 {0, 0},
							 {1, 0},
							 {1, 0},
							 {0, 1},
							 {0, 0} };
		
		three = new TKF91HMM(emit);
		
		
		/* state definitions

		0  Start
		1  w x y z
		2  w x - z
		3  w x y -
		4  w x - -
		5  w - - - 
		6  - - y - (after z)
		7  - - y - (after gap in z)
		8  - - - z
		9  - x y z
	   10  - x - z
	   11  - x y -
	   12  - x - -
	   13  End  */
		
		emit = new int[][] { {0, 0, 0, 0},
							 {1, 1, 1, 1},
							 {1, 1, 0, 1},
							 {1, 1, 1, 0},
							 {1, 1, 0, 0},
							 {1, 0, 0, 0},
							 {0, 0, 1, 0},
							 {0, 0, 1, 0},
							 {0, 0, 0, 1},
							 {0, 1, 1, 1},
							 {0, 1, 0, 1},
							 {0, 1, 1, 0},
							 {0, 1, 0, 0},
							 {0, 0, 0, 0}
							 };
		
		four = new TKF91HMM(emit);
	}
	
	public class TKF91HMM extends HMM {
			
		public TKF91HMM(int[][] em) {
			super(em);
		}

		public void calcTrans(double t){
			
			double l  = params[0];
			double mu = params[1];
			
			double a1 = Math.exp(-mu*t);
			double b1 = l*(1-Math.exp((l-mu)*t))/(mu-l*Math.exp((l-mu)*t));
			double c1 = Math.max(0,1 - mu*(1-Math.exp((l-mu)*t))/((1-Math.exp(-mu*t))*(mu-l*Math.exp((l-mu)*t))));

			trans = new double[states][states];

			for(int i = 0; i < states; i++)
				trans[i][0] = 0;
			
			trans[0][1] = trans[1][1] = trans[3][1] = (1-b1)*l/mu*a1;
			trans[0][2] = trans[1][2] = trans[3][2] = (1-b1)*l/mu*(1-a1);
			trans[0][3] = trans[1][3] = trans[3][3] = b1;
			trans[0][4] = trans[1][4] = trans[3][4] = (1-b1)*(1-l/mu);

			trans[2][1] = (1-c1)*l/mu*a1;
			trans[2][2] = (1-c1)*l/mu*(1-a1);
			trans[2][3] = c1;
			trans[2][4] = (1-c1)*(1-l/mu);

			trans[4][1] = 0;
			trans[4][2] = 0;
			trans[4][3] = 0;
			trans[4][4] = 0;

			for(int i = 0; i < states; i++)
				for(int j = 0; j < states; j++)
					trans[i][j] = Math.log(trans[i][j]);
		}
		
		@Override
		public void calcTrans(double t1, double t2){
			
			double l  = params[0];
			double mu = params[1];
			double p = l / mu;
			
			double a1 = Math.exp(-mu*t1);
			double b1 = l*(1-Math.exp((l-mu)*t1))/(mu-l*Math.exp((l-mu)*t1));
			double c1 = Math.max(0,1 - mu*(1-Math.exp((l-mu)*t1))/((1-Math.exp(-mu*t1))*(mu-l*Math.exp((l-mu)*t1))));

			double a2 = Math.exp(-mu*t2);
			double b2 = l*(1-Math.exp((l-mu)*t2))/(mu-l*Math.exp((l-mu)*t2));
			double c2 = Math.max(0,1 - mu*(1-Math.exp((l-mu)*t2))/((1-Math.exp(-mu*t2))*(mu-l*Math.exp((l-mu)*t2))));

			
			trans = new double[states][states];

			for(int i = 1; i < states-1; i++){
				trans[i][0] = 0;
				trans[i][1] = a1*a2 * p;
				trans[i][2] = (1-a1)*a2 * p;
				trans[i][3] = a1*(1-a2) * p;
				trans[i][4] = (1-a1)*(1-a2) * p;
			}
			
			for(int i = 1; i < 5; i++)
				trans[1][i] *= (1-b1)*(1-b2);
			trans[1][5] = b1;
			trans[1][6] = 0;
			trans[1][7] = (1-b1)*b2;
			trans[1][8] = (1-b1)*(1-b2)*(1-p);
			
			for(int i = 0; i < states; i++)
				trans[0][i] = trans[1][i];
			
			for(int i = 1; i < 5; i++)
				trans[2][i] *= (1-c1)*(1-b2);
			trans[2][5] = c1;
			trans[2][6] = 0;
			trans[2][7] = (1-c1)*b2;
			trans[2][8] = (1-c1)*(1-b2)*(1-p);
			
			for(int i = 1; i < 5; i++)
				trans[3][i] *= (1-b1)*(1-c2);
			trans[3][5] = 0;
			trans[3][6] = b1;
			trans[3][7] = (1-b1)*c2;
			trans[3][8] = (1-b1)*(1-c2)*(1-p);
					
			for(int i = 1; i < 5; i++)
				trans[4][i] *= (1-c1)*(1-c2);
			trans[4][5] = 0;
			trans[4][6] = c1;
			trans[4][7] = (1-c1)*c2;
			trans[4][8] = (1-c1)*(1-c2)*(1-p);
								  
			for(int i = 1; i < 5; i++)
				trans[5][i] *= (1-b1)*(1-b2);
			trans[5][5] = b1;
			trans[5][6] = 0;
			trans[5][7] = (1-b1)*b2;
			trans[5][8] = (1-b1)*(1-b2)*(1-p);
			
			for(int i = 1; i < 5; i++)
				trans[6][i] *= (1-b1)*(1-c2);
			trans[6][5] = 0;
			trans[6][6] = b1;
			trans[6][7] = (1-b1)*c2;
			trans[6][8] = (1-b1)*(1-c2)*(1-p);

			for(int i = 1; i < 5; i++)
				trans[7][i] *= (1-b2);
			trans[7][5] = 0;
			trans[7][6] = 0;
			trans[7][7] = b2;
			trans[7][8] = (1-b2)*(1-p);
			
			double[] sums = new double[states];
			
			
			// force sums to 1 - tested with actual expressions above
			// and differences seemed to be numeric only
			for(int i = 0; i < states-1; i++){
				for(int j = 0; j < states-1; j++){
					sums[i] += trans[i][j];
				}
				trans[i][states-1] = 1 - sums[i];
			}
			
			for(int i = 0; i < states; i++)
				for(int j = 0; j < states; j++)
					trans[i][j] = Math.log(trans[i][j]);
		}

		
		@Override
		public void calcTrans(double t1, double t2, double t3){
			
			double l  = params[0];
			double mu = params[1];
			
			double a1 = Math.exp(-mu*t1);
			double b1 = l*(1-Math.exp((l-mu)*t1))/(mu-l*Math.exp((l-mu)*t1));
			double c1 = Math.max(0,1 - mu*(1-Math.exp((l-mu)*t1))/((1-Math.exp(-mu*t1))*(mu-l*Math.exp((l-mu)*t1))));

			double a2 = Math.exp(-mu*t2);
			double b2 = l*(1-Math.exp((l-mu)*t2))/(mu-l*Math.exp((l-mu)*t2));
			double c2 = Math.max(0,1 - mu*(1-Math.exp((l-mu)*t2))/((1-Math.exp(-mu*t2))*(mu-l*Math.exp((l-mu)*t2))));

			double a3 = Math.exp(-mu*t3);
			double b3 = l*(1-Math.exp((l-mu)*t3))/(mu-l*Math.exp((l-mu)*t3));
			double c3 = Math.max(0,1 - mu*(1-Math.exp((l-mu)*t3))/((1-Math.exp(-mu*t3))*(mu-l*Math.exp((l-mu)*t3))));
		
			trans = new double[states][states];
			
			// Form the matrix as the product of 3 matrices to take advantage of repetition
			double[][] T1 = new double[states][states];
			double[][] T2 = new double[states][states];
			double[][] T3 = new double[states][states];
			
			for(int i = 0; i < states; i++)
				for(int j = 0; j < states; j++)
					T1[i][j] = 1;
			
			for(int i = 0; i < states; i++)
				for(int j = 0; j < states; j++)
					T2[i][j] = 1;
			
			for(int i = 0; i < states; i++)
				for(int j = 0; j < states; j++)
					T3[i][j] = 1;
			
			// states divided into new w: l/mu*(1-b1)
			//            and inserted x: b1
			for(int i = 0; i < states - 1; i++){
				for(int j = 0; j < 6; j++)
					T1[i][j] = l/mu*(1-b1);
				for(int j = 9; j < 13; j++)
					if(i != 5)
						T1[i][j] = b1;
			}
	
			for(int i = 1; i < 6; i++)
				T1[5][i] = l / mu;
			for(int i = 6; i < 9; i++)
				T1[5][i] = 0;

			// elements common to columns (same in each row)
			// eg: every transition to full match state involves a1*a2*a3
			double[] row = new double[] { 0, a1*a2*a3, a1*a2*(1-a3), a1*(1-a2)*a3, a1*(1-a2)*(1-a3), (1-a1),
	                  1, 1, 1, a2*a3, a2*(1-a3), (1-a2)*a3, (1-a2)*(1-a3), (1-b1)*(1 - l/mu) };
			
			for(int i = 0; i < states-1; i++)
				T2[i] = Utils.copyOf(row);
			
			T2[5][13] = (1-l/mu);
			
			// elements common to rows (same in each column)
			double[] col = new double[] { (1-b2)*(1-b3), (1-b2)*(1-b3), (1-b2)*(1-c3), (1-c2)*(1-b3), (1-c2)*(1-c3), 
                    (1-c1), (1-b2)*(1-b3), (1-c2)*(1-b3), (1-b2), (1-b2)*(1-b3), (1-b2)*(1-c3), (1-c2)*(1-b3), 
                    (1-c2)*(1-c3), 0 };
			for(int i = 1; i < 6; i++)
				for(int j = 0; j < states; j++)
					T3[j][i] = col[j];
			for(int j = 0; j < states; j++)
				T3[j][13] = col[j];
			
			col = new double[] { b3, b3, c3, 0, 0, 0, b3, 0, 0, b3, c3, 0, 0, 0 };
			for(int i = 0; i < states; i++)
				T3[i][6] = col[i];
			
			col = new double[] { 0, 0, 0, b3, c3, 0, 0, b3, 0, 0, 0, b3, c3, 0 };
			for(int i = 0; i < states; i++)
				T3[i][7] = col[i];
			
			col = new double[] { (1-b3)*b2, (1-b3)*b2, (1-c3)*b2, (1-b3)*c2, (1-c3)*c2,
				      0, (1-b3)*b2, (1-b3)*c2, b2, (1-b3)*b2, 
				      (1-c3)*b2, (1-b3)*c2, (1-c3)*c2, 0 };
			for(int i = 0; i < states; i++)
				T3[i][8] = col[i];
			
			col = new double[] { (1-b2)*(1-b3), (1-b2)*(1-b3), (1-b2)*(1-c3), (1-c2)*(1-b3), (1-c2)*(1-c3), 
	                 c1, (1-b2)*(1-b3), (1-c2)*(1-b3), (1-b2), (1-b2)*(1-b3), (1-b2)*(1-c3), (1-c2)*(1-b3), 
	                 (1-c2)*(1-c3), 0 };
			for(int i = 9; i < 13; i++)
				for(int j = 0; j < states; j++)
					T3[j][i] = col[j];
			
			for(int i = 0; i < states; i++)
				for(int j = 0; j < states; j++)
					trans[i][j] = T1[i][j] * T2[i][j] * T3[i][j];
			
			trans[states-1][0] = 0;
			
			double[] sums = new double[states];

			// force sums to 1 - tested with actual expressions above
			// and differences seemed to be numeric only
			for(int i = 0; i < states-1; i++){
				for(int j = 0; j < states-1; j++)
					sums[i] += trans[i][j];
				trans[i][states-1] = 1 - sums[i];
			}
			
			for(int i = 0; i < states; i++)
				for(int j = 0; j < states; j++)
					trans[i][j] = Math.log(trans[i][j]); 			
		}
		// /HMM
	}
	
	/**
	 * 
	 * TKF91 states
	 * 0	Start
	 * 1	x y
	 * 2	x -
	 * 3	- y
	 * 4	End
	 * 
	 * @param triAlign
	 * @return
	 */
	
	public int[][] translateAlign(int[] triAlign){
		ArrayList<Integer> alignLeft = new ArrayList<Integer>(0);
		ArrayList<Integer> alignRight = new ArrayList<Integer>(0);
		
		for(int i = 0; i < triAlign.length; i++){
			switch(triAlign[i]){
			case 0:
				alignLeft.add(0);
				alignRight.add(0);
				break;
			case 1:
				alignLeft.add(1);
				alignRight.add(1);
				break;
			case 2:
				alignLeft.add(2);
				alignRight.add(1);
				break;
			case 3:
				alignLeft.add(1);
				alignRight.add(2);
				break;
			case 4:
				alignLeft.add(2);
				alignRight.add(2);
				break;
			case 5:
				alignLeft.add(3);
				break;
			case 6:
				alignLeft.add(3);
				break;
			case 7:
				alignRight.add(3);
				break;
			case 8:
				alignLeft.add(4);
				alignRight.add(4);
				break;			
			}
		}
		int[][] aligns = new int[2][];
		aligns[0] = new int[alignLeft.size()];
		for(int i = 0; i < aligns[0].length; i++)
			aligns[0][i] = alignLeft.get(i);
		aligns[1] = new int[alignRight.size()];
		for(int i = 0; i < aligns[1].length; i++)
			aligns[1][i] = alignRight.get(i);
		
		return aligns;
	}
	
	public int[] translate4Align(int[][] align){
		ArrayList<Integer> states = new ArrayList<Integer>(0);
		states.add(0);
		int bin;
		int[] temp = new int[four.emit.length];
		int[] map = new int[(int)Math.pow(2, four.emit[0].length)];
		for(int i = 0; i < four.emit.length; i++)
			for(int j = 0; j < four.emit[0].length; j++)
				temp[i] += Math.pow(2, j) * four.emit[i][j];
		
		System.out.println("map length: " + map.length);
		for(int i = 0; i < temp.length; i++)
			map[temp[i]] = i;
			
			
		for(int i = 0; i < align[0].length; i++){
			bin = 0;
			for(int j = 0; j < align.length; j++)
				bin += Math.pow(2, j) * align[j][i];
			if(bin != 0){
				if(bin == 2){
					if(four.trans[states.get(states.size())][6] > -Utils.log0)
						states.add(6);
					else
						states.add(7);
				}
				states.add(map[bin]);
			}
			
		}
		// states.add(four.states-1);
		int[] fourAlign = new int[states.size()];
		for(int i = 0; i < states.size(); i++)
			fourAlign[i] = states.get(i);
			
		return fourAlign;
	}
	// /TKF91
}
