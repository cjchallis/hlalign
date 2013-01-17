package hlalign.base;

public class MCMC {
	public int B;
	public int N;
	
	String[] names;
	double[][][] coords;
	double[] sigma2;
	double[][] indelPars;
	Tree[] tree;
	Alignment[] align;
	
	public MCMC(int b, int n, double[][][] c, String[] pnames){
		B = b;
		N = n;
		names = pnames;
		sigma2 = new double[b+n];
		tree = new Tree[b+n];
		align = new Alignment[b+n];
		indelPars = new double[b+n][];
		coords = c;
	}
	
	public void run(){
		initialize();
	}
	
	public void initialize(){
		sigma2[0] = 1;
		initTree();
		initAlign();
		
	}
	
	public Tree initTree(){
		return new Tree(names);
	}
	
	public Alignment initAlign(){
		return new Alignment();
	}
	
}
