package hlalign.base;

// import statalign.model.ext.plugins.StructAlign.MultiNormCholesky;

public class MCMC {
	public int B;
	public int N;
	
	String[] names;
	double[][][] coords;
	double[] sigma2;
	double[][] indelPars;
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
		indelPars = new double[b+n][];
		coords = c;
	}
	
	public void run(){
		initialize();
		System.out.println("LL: " + structure.logLikelihood(tree));
	}
	
	public void initialize(){
		tree = new Tree(names, coords);
		structure = new Structure(coords);
	}
}
