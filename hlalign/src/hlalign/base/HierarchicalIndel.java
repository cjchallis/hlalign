package hlalign.base;

import hlalign.indels.*;

public class HierarchicalIndel extends IndelModel {
	/** Indel model governing super links */
	IndelModel upper;
	double[] upparams;
	
	/** Indel model governing normal links */
	IndelModel[] lower;
	double[][] lowparams;
	
	HierHMM two;
	
	int types;
	double[] pi;
	
	public HierarchicalIndel(IndelModel up, String low, int num){
		upper = up;
		lower = new IndelModel[num];
		switch(low){
		case "TKF91":
			for(int i = 0; i < num; i++)
				lower[i] = new TKF91(lowparams[i]);
			break;
		case "TKF92":
			for(int i = 0; i < num; i++)
				lower[i] = new TKF92(lowparams[i]);
			break;
		default:
			System.out.println("No such indel model: " + low);
			break;			
		}
		types = num;
	}
	
	public class HierHMM extends HMM{
		public HierHMM(int[][] em){
			super(em);
		}
		
		public void calcTrans(double t){
			trans = new double[5*types+3][5*types+3];
			for(int i = 0; i < types; i++){
				lower[i].two.calcTrans(t);
				double[][] subTrans = lower[i].two.trans;
				int shift = 5*i;
				for(int j = 1; j < 4; j++)
					for(int k = 1; k < 4; k++)
						trans[shift + j][shift + k] = subTrans[j][k];
				
			}
		}
	}
}
