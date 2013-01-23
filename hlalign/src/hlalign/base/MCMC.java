package hlalign.base;

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
	}
	
	public void run(){
		initialize();
		System.out.println("LL: " + structure.logLikelihood(tree));
		
		OU ou = new OU(1, 10, 1);
		double[] pX = ou.marginal(coords[0]);
		double[] pY = ou.marginal(coords[1]);
		double[][] pXY = ou.joint(coords[0], coords[1], 1);
		
		TKF91 tkf91 = new TKF91(.03, .033);
		tkf91.calcTrans(1);
		
		PairDP dp = new PairDP(pX, pY, pXY, tkf91);
		
		dp.forward();
		System.out.println("ML: " + dp.getMarginalLikelihood());
		int[] align = dp.backward();
		for(int i = 0; i < align.length; i++)
			System.out.print(align[i]);
		
	}
	
	public void initialize(){
		tree = new Tree(names, coords);
		structure = new Structure(coords);
	}
}
