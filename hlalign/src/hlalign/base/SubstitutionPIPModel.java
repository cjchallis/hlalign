package hlalign.base;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

public class SubstitutionPIPModel {
	public double[][] subQ; // substitution-only matrix, read from data
	public double[][] Q;	// substitution-deletion matrix, calculate from subQ and mu
	
	public double[][] u;	// eigenvector matrix of Q
	public double[] d;		// eigenvalues of Q
	public double[][] v;	// inverse eigenvector matrix of Q
	public EigenDecomposition eigen;
	
	public double[] pi;
	public double mu = .1;
	public double l = 15;
	
	public char[] alphabet;
	public String alphaString;
	
	int n;
	
	public SubstitutionPIPModel(double[][] matrix, double[] e){
		subQ = matrix;
		pi = e;
		Q = new double[subQ.length + 1][subQ.length + 1];
		d = new double[Q.length];
		for(int i = 0; i < subQ.length; i++)
			for(int j = 0; j < subQ.length; j++)
				Q[i][j] = subQ[i][j];
		n = Q.length;
		updateDecomposition();
	}

	public void updateDecomposition(){
		for(int i = 0; i < subQ.length; i++){
			Q[i][i] = subQ[i][i] - mu;
			Q[i][Q.length-1] = mu;
		}
		RealMatrix rmQ = new Array2DRowRealMatrix(Q);
		eigen = new EigenDecomposition(rmQ);
		u = eigen.getV().getData();
		RealMatrix dd = eigen.getD();
		for(int i = 0; i < Q.length; i++)
			d[i] = dd.getEntry(i,i);
		v = new LUDecomposition(eigen.getV()).getSolver().getInverse().getData();
	}
	
	public double[][] calcSubMatrix(double t){
		double exp_dtk;
		double [][] subMatrix = new double[n][n];
		for(int k = 0; k < n; k++){
			exp_dtk = Math.exp(d[k]*t);
			for(int i = 0; i < n; i++){
				for(int j = 0; j < n; j++){
					subMatrix[i][j] += u[i][k] * exp_dtk * v[k][j];
				}
			}
		}
		for(int i = 0; i < n; i++)
			for(int j = 0; j < n; j++)
				subMatrix[i][j] = Math.log(subMatrix[i][j]);
		return subMatrix;
	}
}

