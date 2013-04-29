package hlalign.indels;
import hlalign.base.IndelModel;
import hlalign.base.Utils;

public class TKF92 extends IndelModel{
	TKF91 tkf91;
	
	public TKF92(double[] par){
		params = par.clone();
		double l = params[0];
		double mu = params[1];
		// r = params[2];
		
		tkf91 = new TKF91(new double[]{l, mu});
		two = new TKF92HMM(tkf91.two.emit);
		three = new TKF92HMM(tkf91.three.emit);
		four = new TKF92HMM(tkf91.four.emit);
	}
	
	public class TKF92HMM extends HMM{
		public TKF92HMM(int[][] em){
			super(em);
		}
		
		public void calcTrans(double t){
			tkf91.two.calcTrans(t);
			trans = tkf91.two.trans;
			double lr = Math.log(params[2]);
			double lmr = Math.log(1-params[2]);
			// Multiply all transition probabilities (except from start or from end) by (1-r)
			for(int i = 1; i < trans.length - 1; i++)
				for(int j = 1; j < trans.length; j++)
					trans[i][j] += lmr;
			// Add r to self-transition probabilities
			for(int i = 1; i < trans.length - 1; i++)
				trans[i][i] = Utils.logAdd(trans[i][i], lr);
			
			System.out.println("Internal: ");
			for(int i = 0; i < trans.length; i++){
				for(int j = 0; j < trans.length; j++){
					System.out.print(Math.exp(trans[i][j]) + "  ");
				}
				System.out.println();
			}
		}
		
		public void calcTrans(double t1, double t2){
			tkf91.three.calcTrans(t1, t2);
			trans = tkf91.three.trans;
			double lr = Math.log(params[2]);
			double lmr = Math.log(1-params[2]);
			for(int i = 1; i < trans.length - 1; i++)
				for(int j = 1; j < trans.length - 1; j++)
					trans[i][j] += lmr;
			for(int i = 1; i < trans.length - 1; i++)
				trans[i][i] = Utils.logAdd(trans[i][i], lr);
		}
	}
}
