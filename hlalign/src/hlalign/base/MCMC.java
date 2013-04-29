package hlalign.base;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import hlalign.indels.*;
import hlalign.heterogeneous.*;

// import statalign.model.ext.plugins.StructAlign.MultiNormCholesky;

public class MCMC {
	public int B;
	public int N;
	
	String[] names;
	double[][][] coords;
	char[][] seqs;
	double[] sigma2;
	double[][] fullCovar;
	Tree tree;
	Alignment[] align;
	// int[][] alignArray;
	Structure structure;
	Residue[][] resAlign;
	
	public MCMC(int b, int n, double[][][] c, char[][] s, String[] pnames){
		B = b;
		N = n;
		names = pnames;
		sigma2 = new double[b+n];
		//tree = new Tree[b+n];
		align = new Alignment[b+n];
		coords = c;
		seqs = s;
		
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
				
		TKF91 tkf91 = new TKF91(new double[]{.03, .033});
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


		System.out.println("TKF92 matrix sums:");
		TKF92 tkf92 = new TKF92(new double[] {.03, .033, .7});
		tkf92.two.calcTrans(1);
		
		sums = new double[tkf92.two.states];
		
		for(int i = 0; i < tkf92.two.states; i++){
			for(int j = 0; j < tkf92.two.states; j++){
				sums[i] += Math.exp(tkf92.two.trans[i][j]);
			}
			System.out.println(sums[i]);
		}
		
		
		for(int i = 0; i < tkf92.two.states; i++){
			for(int j = 0; j < tkf92.two.states; j++){
				System.out.print(Math.exp(tkf92.two.trans[i][j]) + "  ");
			}
			System.out.println();
		}
		
		convertAlign(tree.alignArray);


	}
	
	public void initialize(){
		structure = new Structure(coords, .1, 100, .1);
		tree = new Tree(names, coords, this);
	}
	
	public void convertAlign(int[][] align){
		int n = align.length;
		int[] lengths = new int[n];
		for(int i = 0; i < n; i++)
			for(int j = 0; j < align[0].length; j++)
				lengths[i] += align[i][j];
		resAlign = new Residue[n][];
		for(int i = 0; i < n; i++)
			resAlign[i] = new Residue[lengths[i]];
		
		for(int i = 0; i < n; i++)
			for(int j = 0; j < lengths[i]; j++)
				resAlign[i][j] = new Residue(j);
		
		int[] idx = new int[n];
		int[] parents = new int[n];
		int[] lefts = new int[n];
		int[] rights = new int[n];
		for(int i = 0; i < n; i++){
			if(tree.vertex[i].parent != null)
				parents[i] = tree.vertex[i].parent.index;
			else
				parents[i] = -1;
			if(tree.vertex[i].left != null)
				lefts[i] = tree.vertex[i].left.index;
			else
				lefts[i] = -1;
			if(tree.vertex[i].right != null)
				rights[i] = tree.vertex[i].right.index;
			else
				rights[i] = -1;
			
		}
			
		
		for(int j = 0; j < align[0].length; j++){
			for(int i = 0; i < align.length; i++){
				if(align[i][j] == 1){
					if(parents[i] > -1)
						if(align[parents[i]][j] == 1)
							resAlign[i][idx[i]].parent = resAlign[parents[i]][idx[parents[i]]];
					if(lefts[i] > -1){
						if(align[lefts[i]][j] == 1)
							resAlign[i][idx[i]].left = resAlign[lefts[i]][idx[lefts[i]]];
						if(align[rights[i]][j] == 1)
							resAlign[i][idx[i]].right = resAlign[rights[i]][idx[rights[i]]];
					}
				}
			}
		}		
	}
	
}
