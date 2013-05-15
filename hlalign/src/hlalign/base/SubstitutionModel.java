package hlalign.base;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.jblas.MatrixFunctions;

public class SubstitutionModel {
	public double[][] theta; // substitution-only matrix, read from data
	public double[][] Q;	// substitution-deletion matrix, calculate from theta and mu
	
	public double[][] u;	// eigenvector matrix of Q
	public double [] d;		// eigenvalues of Q
	public double [][] v;	// inverse eigenvector matrix of Q
	public EigenDecomposition eigen;
	
	public double[] pi;

	int n = d.length;

	public char[] alphabet;
	
	public SubstitutionModel(double[][] matrix, double mu){
		theta = matrix;
		Q = new double[theta.length + 1][theta.length + 1];
		d = new double[Q.length];
		for(int i = 0; i < theta.length; i++)
			for(int j = 0; j < theta.length; j++)
				Q[i][j] = theta[i][j];
		updateDecomposition(mu);
	}

	public void updateDecomposition(double mu){
		for(int i = 0; i < theta.length; i++){
			Q[i][i] = theta[i][i] - mu;
			Q[i][Q.length-1] = mu;
		}
		RealMatrix rmQ = new Array2DRowRealMatrix(Q);
		eigen = new EigenDecomposition(rmQ, 0);
		u = eigen.getV().getData();
		RealMatrix dd = eigen.getD();
		for(int i = 0; i < Q.length; i++)
			d[i] = dd.getEntry(i,i);
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
