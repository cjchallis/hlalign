package hlalign.base;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.distribution.NormalDistribution;

public class Structure {
	
	public double sigma2;
	public double epsilon;
	public double tau;
	public double curLogLike;
	public int[][] align;
	public double [][][] coords;
	public double[][] fullCovar;
	
	public Structure(double [][][] c, double s2, double t, double e){
		coords = c;
		sigma2 = s2;
		tau = t;
		epsilon = e;
	}
	
	public double logLikelihood(Tree tree){
		align = tree.alignArray;
		
		fullCovar = calcFullCovar(tree);
	
		
		double logli = calcAllColumnContrib();
		curLogLike = logli;
		
		// testing
		//System.out.println("Total log likelihood " + curLogLike);
		
		return curLogLike;
	}
	
	public double[][] calcFullCovar(Tree tree){
		double[][] distMat = new double[tree.names.length][tree.names.length];
		calcDistanceMatrix(tree.root, distMat);
		//System.out.print("Distance: " + distMat[0][1]);
		
		//System.out.println("Current tree:");
		//printTree(tree.root, "o");
		
		for(int i = 0; i < tree.names.length; i++)
			for(int j = i; j < tree.names.length; j++)
				distMat[j][i] = distMat[i][j] = tau * Math.exp(-distMat[i][j]);
		for(int i = 0; i < tree.names.length; i++)
			distMat[i][i] += epsilon;
		return distMat;
	}
	
	/**
	 * recursive algorithm to traverse tree and calculate distance matrix between leaves 
	 */		
	public int[] calcDistanceMatrix(Vertex v, double[][] distMat){
		int[] subTree = new int[distMat.length + 1];
		
		// either both left and right are null or neither is
		if(v.left != null){
			int[] subLeft  = calcDistanceMatrix(v.left, distMat);
			int[] subRight = calcDistanceMatrix(v.right, distMat);
			int i = 0;
			while(subLeft[i] > -1){
				subTree[i] = subLeft[i];
				i++;
			}
			for(int j = 0; i+j < subTree.length; j++)
				subTree[i+j] = subRight[j];
		}
		else{
			subTree[0] = v.index;
			for(int j = 1; j < subTree.length; j++)
				subTree[j] = -1;
		}

		addEdgeLength(distMat, subTree, v.edgeLength * sigma2 / (2*tau));	

		/*System.out.println();
		System.out.println("Distmat:");
		for(int i = 0; i < distMat.length; i++)
			for(int j = 0; j < distMat[0].length; j++)
				System.out.println(distMat[i][j]);*/
		return subTree;
	}
		
	/** adds the length of the current edge to the distance between all leaves
	 * of a subtree to all other leaves
	 * 'rows' contains the indices of vertices in the subtree */
	public void addEdgeLength(double[][] distMat, int[] subTree, double edgeLength){
		
		int i = 0;
		while(subTree[i] > -1){
			for(int j = 0; j < distMat.length; j++){  
				distMat[subTree[i]][j] += edgeLength;
				distMat[j][subTree[i]] += edgeLength;
			}
			i++;		
		}
			
		// edge length should not be added to distance between vertices in the subtree
		// subtract the value from these entries of the distance matrix
		i = 0;
		while(subTree[i] > -1){
			int j = 0;
			while(subTree[j] > -1){
				distMat[subTree[i]][subTree[j]] -= edgeLength;
				distMat[subTree[j]][subTree[i]] -= edgeLength;
				j++;
			}
			i++;
		}
	}

	
	
	public double calcAllColumnContrib(){
		double logli = 0;
		int l = coords.length;
		int[] inds = new int[l];		// current char indices
		int[] col = new int[l];  
		for(int i = 0; i < align[0].length; i++) {
			for(int j = 0; j < l; j++)
				col[j] = align[j][i] == 0 ? -1 : inds[j]++;
			double ll = columnContrib(col); 
			logli += ll;
			//System.out.println("Column: " + Arrays.toString(col) + "  ll: " + ll);
		}
		return logli;
	}
	
	/**
	 * Calculates the structural likelihood contribution of a single alignment column
	 * @param col the column, id of the residue for each sequence (or -1 if gapped in column)
	 * @return the likelihood contribution
	 */
	public double columnContrib(int[] col) {
		// count the number of ungapped positions in the column
		int numMatch = 0;
		for(int i = 0; i < col.length; i++){
			if(col[i] != -1)
				numMatch++;
		}
		if(numMatch == 0)  
			return 1;
		// collect indices of ungapped positions
		int[] notgap = new int[numMatch];
		int j = 0;
		for(int i = 0; i < col.length; i++)
			if(col[i] != -1)
				notgap[j++] = i;
		
		// extract covariance corresponding to ungapped positions
		double[][] subCovar = getSubMatrix(fullCovar, notgap, notgap);
		// create normal distribution with mean 0 and covariance subCovar
		MultiNormCholesky multiNorm = new MultiNormCholesky(new double[numMatch], subCovar);
		
		double logli = 0;
		double[] vals = new double[numMatch];
		// loop over all 3 coordinates
		
		/*System.out.println("Calculating log likelihood: ");
		System.out.println("Mean: " + Arrays.toString(multiNorm.getMeans()));
		System.out.println("Variance: " + Arrays.toString(subCovar[0]));*/
		for(j = 0; j < 3; j++){
			for(int i = 0; i < numMatch; i++)
				vals[i] = coords[notgap[i]][col[notgap[i]]][j];
			//System.out.println("Values: " + Arrays.toString(vals));
			logli += multiNorm.logDensity(vals);
			//System.out.println("LL: " + multiNorm.logDensity(vals));
		}
		return logli;
	}

	/**
	 * extracts the specified rows and columns of a 2d array
	 * @param matrix, 2d array from which to extract; rows, rows to extract; cols, columns to extract
	 * @return submatrix
	 */
	public double[][] getSubMatrix(double[][] matrix, int[] rows, int[] cols) {
		double[][] submat = new double[rows.length][cols.length];
		for(int i = 0; i < rows.length; i++)
			for(int j = 0; j < cols.length; j++)
				submat[i][j] = matrix[rows[i]][cols[j]];
		return submat;
	}

	/**
	 * calculates log probabilities from stationary OU process of the structure X
	 * @param X, a protein structure
	 * @return array of log probabilities, for use in PairDP 
	 */
	public double[] marginal(double[][] X){
		double[] pX = new double[X.length + 2];
		pX[0] = pX[1] = Utils.log0;
		
		NormalDistribution stationary = new NormalDistribution(0, tau);

		for(int i = 2; i < pX.length; i++)
			for(int j = 0; j < 3; j++)
				pX[i] += Math.log(stationary.density(X[i-2][j]));
		
		return pX;		    
	}

	/**
	 * calculates log probabilities from joint distribution of OU process
	 * @param 	X, Y: protein structures, 
	 * @param 	t: branch length between @X and @Y
	 * @return matrix of log probabilities - every pairwise match between @X and @Y, for use in PairDP
	 */
	public double[][] joint(double[][] X, double[][] Y, double t){
		
		int lx = X.length + 2;
		int ly = Y.length + 2;
		
		double[][] pXY = new double[lx][ly];
		
		double var = tau + epsilon, covar = tau * Math.exp(-t * sigma2 / (2 * tau));
		
		double[][] covMat = new double[][] { {var, covar}, {covar, var} };

		MultiNormCholesky joint = new MultiNormCholesky(new double[2], covMat);
		
		for(int i = 0; i < 2; i++)
			for(int j = 0; j < ly; j++)
				pXY[i][j] = Utils.log0;

		for(int i = 2; i < lx; i++)
			for(int j = 0; j < 2; j++)
				pXY[i][j] = Utils.log0;
			
		for(int i = 2; i < lx; i++)
			for(int j = 2; j < ly; j++)
				for(int k = 0; k < 3; k++)
					pXY[i][j] += joint.logDensity(new double[] {X[i-2][k], Y[j-2][k]});
		
		return pXY;
	}
	
	/**
	 * calculates log probabilities of a set of leaves in a subtree or complement of the tree
	 * needed for realigning an internal node and its parent
	 * @param coords   	Array containing all protein coordinates
	 * @param align		int array with 1's and 0's representing alignment
	 * @param vertex	int with the index of a vertex
	 * @param subTree	if true, use leaves in subtree of @vertex, false, use all other leaves
	 * @param tree  	tree
	 * @return array of log probabilities, for use in PairDP
	 */
	public double[] marginTree(double[][][] coords, int[][] alignInds, int vertex, boolean subTree, Tree tree){
		ArrayList<Integer> complement = collectLeaves(tree.vertex[vertex]);
		if(subTree){
			ArrayList<Integer> leaves = new ArrayList<Integer>(0);
			for(int i = 0; i < coords.length; i++)
				leaves.add(i);
			for(int i = 0; i < complement.size(); i++)
				leaves.remove(complement.get(i));
			complement = leaves;
		} else { vertex = tree.vertex[vertex].parent.index; }
		
		fullCovar = calcFullCovar(tree);
		
		int l = coords.length;
		int[] col = new int[l];
		
		int lx = 2;
		for(int i = 0; i < alignInds[0].length; i++)
			lx += alignInds[vertex][i] == -1 ? 0 : 1;
			
		double[] probs = new double[lx];
		probs[0] = probs[1] = Utils.log0;	// dynamic programming requires 2 leading log0's
		int k = 2;
		
		for(int i = 0; i < alignInds[0].length; i++){
			if(alignInds[vertex][i] != -1){
				for(int j = 0; j < l; j++){
					if(complement.contains(j))
						col[j] = -1;
					else
						col[j] = alignInds[j][i];
				}
				/*System.out.println("Column for i = " + i);
				for(int j = 0; j < col.length; j++)
					System.out.println(col[j]);*/
				probs[k] = columnContrib(col);
				k++;
				
			}
		}
		
		return probs;
	}
	
	/**
	 * calculates joint probabilities for realignment of an internal node and its parent
	 * 
	 * @param coords   	Array containing all protein coordinates
	 * @param align		int array with 1's and 0's representing alignment
	 * @param vertex	int with the index of a vertex
	 * @param subTree	if true, use leaves in subtree of @vertex, false, use all other leaves
	 * @param tree  	tree
	 * @return array of log probabilities, for use in PairDP
	 */
	public double[][] jointTree(double[][][] coords, int[][] alignInds, int vertex, boolean children, Tree tree){
		int v1;
		int v2;
		
		ArrayList<Integer> subTree1 = new ArrayList<Integer>(0);
		subTree1.add(-1);
		ArrayList<Integer> subTree2;
		
		if(!children){
			v1 = tree.vertex[vertex].parent.index;
			v2 = vertex;
			subTree2 = collectLeaves(tree.vertex[v2]);
		} else {
			v1 = tree.vertex[vertex].left.index;
			v2 = tree.vertex[vertex].right.index;
			subTree1 = collectLeaves(tree.vertex[v1]);
			subTree2 = collectLeaves(tree.vertex[v2]);
		}
		
		int lx = 2, ly = 2;		// 
		for(int i = 0; i < alignInds[v1].length; i++)
			lx += alignInds[v1][i] == -1 ? 0 : 1;
		for(int i = 0; i < alignInds[v2].length; i++)
			ly += alignInds[v2][i] == -1 ? 0 : 1;
		
		double[][] pXY = new double[lx][ly];
		
		System.out.println("lx: " + lx);
		System.out.println("ly: " + ly);
			
		for(int i = 0; i < 2; i++)
			for(int j = 0; j < ly; j++)
				pXY[i][j] = Utils.log0;

		for(int i = 2; i < lx; i++)
			for(int j = 0; j < 2; j++)
				pXY[i][j] = Utils.log0;
		
		fullCovar = calcFullCovar(tree);
		
		int l = coords.length;
		int[] col = new int[l];
		
		int x = 2, y = 2;

		for(int i = 0; i < alignInds[v1].length; i++){
			if(alignInds[v1][i] != -1){
				y = 2;
				for(int j = 0; j < alignInds[v2].length; j++){
					if(alignInds[v2][j] != -1){
						for(int k = 0; k < l; k++){
							if(subTree2.contains(k))
								col[k] = alignInds[k][j];
							else if(!children || subTree1.contains(k))
								col[k] = alignInds[k][i];
							else
								col[k] = -1;
						}
						//System.out.println("x: " + x + " y: " + y);
						pXY[x][y] = columnContrib(col);
						y++;
						
						//System.out.println("Column for i = " + i + " and j = " + j + ":");
						//for(int k = 0; k < col.length; k++)
						//	System.out.println(col[k]);
					}
				}
				x++;
			}
		}
		return pXY;
	}
	
	public int[][] makeAlignInds(int[][] align){	
		int[] inds = new int[align.length];
		int[][] alignInds = new int[align.length][align[0].length];
		
		for(int i = 0; i < align.length; i++)
			for(int j = 0; j < align[0].length; j++)
				alignInds[i][j] = align[i][j] == 0 ? -1 : inds[i]++;
		return alignInds;
	}
	
	public ArrayList<Integer> collectLeaves(Vertex v){
		ArrayList<Integer> inds = new ArrayList<Integer>(0);
		moveDown(v, inds);
		return inds;
	}
	
	public void moveDown(Vertex v, List<Integer> inds){
		if(v.left != null){
			moveDown(v.left, inds);
			moveDown(v.right, inds);
		}
		else
			inds.add(v.index);
	}
	
	
	
}
