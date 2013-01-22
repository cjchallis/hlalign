package halign.indels;

import hlalign.base.IndelModel;

public class TKF91 extends IndelModel {
	
	public TKF91(double l, double mu){
		/** Start, Match, Delete, Insert, End */
		emit = new int[][] { {0, 0},
							 {1, 1},
							 {1, 0},
							 {0, 1},
							 {0, 0} };
		
		states = emit.length;
		
		emitVal = new int[states];
		for(int i = 0; i < states; i++)
			emitVal[i] = 2*emit[i][0] + emit[i][0];
		
		/** lambda, mu */
		params = new double[] {l, mu}; 
		
		/** only 1 type in TKF91 */
		type = new int[states];
		
	}
	
	@Override
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
		trans[4][4] = 1;

		for(int i = 0; i < states; i++)
			for(int j = 0; j < states; j++)
				trans[i][j] = Math.log(trans[i][j]);
	}
}
