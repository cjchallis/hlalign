package hlalign.base;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import hlalign.indels.*;

// import statalign.model.ext.plugins.StructAlign.MultiNormCholesky;

public class MCMC {
	public int B;
	public int N;
	
	String[] names;
	double[][][] coords;
	double[] sigma2;
	double[][] fullCovar;
	Tree tree;
	Alignment[] align;
	// int[][] alignArray;
	Structure structure;
	
	public MCMC(int b, int n, double[][][] c, String[] pnames){
		B = b;
		N = n;
		names = pnames;
		sigma2 = new double[b+n];
		//tree = new Tree[b+n];
		align = new Alignment[b+n];
		coords = c;
		
		for(int i = 0; i < coords.length; i++){
			RealMatrix temp = new Array2DRowRealMatrix(coords[i]);
			RealVector mean = Utils.meanVector(temp);
			for(int j = 0; j < coords[i].length; j++)
				coords[i][j]= temp.getRowVector(j).subtract(mean).toArray();
		}
	}
	
	public void run(){
		initialize();
		System.out.println("LL: " + structure.logLikelihood(tree));
		
		
		double[] pX = structure.marginal(coords[0]);
		double[] pY = structure.marginal(coords[1]);
		double[][] pXY = structure.joint(coords[0], coords[1], 1);
		
		for(int i = 0; i < 6; i++)
			System.out.println(pX[i]);
		for(int i = 0; i < 6; i++)
			System.out.println(pY[i]);
		for(int i = 0; i < 6; i++){
			for(int j = 0; j < 6; j++)
				System.out.print(pXY[i][j]);
			System.out.println();
		}
				
		TKF91 tkf91 = new TKF91(.03, .033);
		tkf91.two.calcTrans(1);
		
		PairDP dp = new PairDP(pX, pY, pXY, tkf91.two);
		
		dp.forward();
		System.out.println("ML: " + dp.getMarginalLikelihood());
		int[] align = dp.backward();
		System.out.println("Single align:");
		for(int i = 0; i < align.length; i++)
			System.out.print(align[i]);
		
		int[][] alignInds = structure.makeAlignInds(tree.alignArray);
		structure.marginTree(coords, alignInds, 5, true, tree);
		structure.marginTree(coords, alignInds, 5, false, tree);
		structure.jointTree(coords, alignInds, 5, false, tree);
		
		double[] sums = new double[tkf91.three.states];
		tkf91.three.calcTrans(3, 1);
		
		for(int i = 0; i < tkf91.three.states; i++){
			for(int j = 0; j < tkf91.three.states; j++){
				sums[i] += Math.exp(tkf91.three.trans[i][j]);
			}
			System.out.println(sums[i]);
		}
		
		System.out.println("Four:");
		
		tkf91.four.calcTrans(1,2,3);
		// TODO Figure out where these lines should actually go
		tkf91.four.silent[12] = true;
		tkf91.four.emitVal = new int[tkf91.four.states];
		for(int i = 1; i < tkf91.four.states - 3; i++)
			tkf91.four.emitVal[i] = 1;
		int[][] fourAlign = tree.vertex[5].get4Align();
		System.out.println("fourAlign: " + fourAlign[0].length);
		int[] stateAlign = tkf91.translate4Align(fourAlign);
		System.out.println("State align: " + stateAlign.length);
		SingleDP trainWreck = new SingleDP(tkf91.four, stateAlign, false);
		trainWreck.forward();
		int[] new4Align = trainWreck.backward();
		
		for(int i = 0; i < stateAlign.length; i++)
			System.out.print(stateAlign[i] + " ");
		System.out.println();
		for(int i = 0; i < new4Align.length; i++)
			System.out.print(new4Align[i] + " ");
		

	}
	
	public void initialize(){
		structure = new Structure(coords, .1, 100, .1);
		tree = new Tree(names, coords, this);
	}
}
