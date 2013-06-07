package hlalign.base;

import java.util.ArrayList;

/**
 */

public class Vertex {

    static final int ROUNDING = 100; // tells the number of digits in rounding when printing the tree
    static final double SELECTING = 0.5; /*this is the probability for selecting a 
    										homologous column not to be changed during topology changing
	 */

    static final double EMPTY_WINDOW = 0.01; /* this is the probability that an empty window will be realigned*/

    Tree owner;

    Vertex old;

    /** The name of the sequence associated to the vertex.
     * Used only if the vertex is a leaf of the tree. */
    public String name;
    
    /** Is this vertex a leaf? */
    public boolean leaf;
    
    /** Index of this vertex in the Tree.vertex array */
    public int index;

    /** This reference points to the parent of the vertex */
    public Vertex parent;
    
    /** This reference points to the left child of the vertex. If the
     * vertex is a leaf, it is set to null. */
    public Vertex left;
    
    /** This reference points to the right child of the vertex. If the
     * vertex is a leaf, it is set to null. */
    public Vertex right;

    /** The length of the edge that connects this vertex with its parent. */
    public double edgeLength;

    int length;									// sequence length
    String seq;									// original sequence of this Vertex (given for leaves only)
    int[] intSeq;								// sequence converted to integers
    
    /** alignments between Left and Right children */
    public int[] alignL;
    public int[] alignR;
    
    public int descentOrder;
    
    public final int GAP;
    
    public ArrayList<Integer> descChars;
    
    public double[][] subMatrix;				// (log) substitution probabilities from parent to this vertex
    
    // the following are all column-specific and relate to the Poisson Indel Process
    public double survival;						// survival probability
    public double commonAncestor;				// if a common ancestor of the column, Utils.log0 if not
    public double[] fels;						// felsenstein vector, \tilde{f}_v(\sigma)
    public double felsSum;						// weighted sum by Pi over felsenstein vector of this vertex, \tilde{f}_v
    public double firstVertex;					// probability of alignment column given that this was the first vertex
    												// after insertion of the original character, f_v
    												// felsenstein sum adjusted by survival probability and zeroed 
    												// if this vertex is not ancestral to all characters in the column
    public double priorFirst;					// prior probability that this was the first vertex
    
    public int leafCount;

    Vertex(int i, String nm, Tree tree) {
    	// edgeLength = 1.0;
    	left = null;
    	right = null;
    	name = nm;
    	index = i;
    	owner = tree;
    	fels = new double[owner.owner.subModel.Q.length];
    	GAP = fels.length - 1;
    }
    
    /**
     * Calculate the felsenstein vector for this vertex (for a particular column)
     * from its children
     * Uses an alignment matrix containing indices of characters, not just +/-
     * For PIP, gap character is last index (20 for proteins, 0-19 amino acids)
     */
    public void felsenstein(int col){
    	if(leaf){
    		for(int i = 0; i < fels.length; i++)			
    			fels[i] = Utils.log0;
    		fels[owner.align.matrix[index][col]] = 0;  // set fels of leaf character to log(1)
    	} else{
    		for(int i = 0; i < fels.length; i++){
    			double l = Utils.log0;
    			double r = Utils.log0;
    			//System.out.println("Vertex " + index + "col " + col);
    			for(int j = 0; j < fels.length; j++){
    				//System.out.println(left.subMatrix[i][j] + " " + left.fels[j]);
    				l = Utils.logAdd(l, left.subMatrix[i][j] + left.fels[j]);
    				r = Utils.logAdd(r, right.subMatrix[i][j] + right.fels[j]);
    			}
    			fels[i] = l + r;
    		}  	
    	}
    }
    
    /**
     * Calculate the pi-weighted felsenstein sum, equivalent to substitution likelihood
     * if the alignment column originated at this vertex
     * This quantity is \tilde{f}_v from the PIP paper
     */
    public void calcFelsSum(){
    	felsSum = Utils.log0;
    	for(int i = 0; i < fels.length-1; i++)
    		felsSum = Utils.logAdd(felsSum, owner.owner.subModel.pi[i] + fels[i]);
    }
    
    /**
     * Calculate probability of column history given this as the first vertex
     * This quantity is f_v from the PIP paper
     * @param empty Calculate for empty column?
     */
    public void calcFirstVertex(boolean empty){
    	if(this == owner.root)
    		firstVertex = felsSum;
    	else if(!empty)
    		firstVertex = commonAncestor + survival + felsSum;
    	else
    		firstVertex = Math.log(1 + Math.exp(survival)*(1-Math.exp(felsSum)));			// calc probability of empty column
    																						// given this as the first non-empty vertex  		
    }
    
    /**
     * Calculate survival probability for this node.  This quantity is \beta(v) from PIP paper 
     */
    public void calcSurvival(){
    	if(this == owner.root)
    		survival = 0;
    	else
    		survival = -Math.log(edgeLength) - Math.log(owner.owner.subModel.mu) + Math.log(1 - Math.exp(-edgeLength * owner.owner.subModel.mu));
    }
    
    /**
     * Aggregate the descendant characters from each vertex into an ArrayList
     * Later used for checking which vertices are ancestral to all characters in a column
     * (step 4 of PIP appendix)
     */
    public void findDescChars(int col){
    	descChars = new ArrayList<Integer>(0);
    	if(leaf)
    		descChars.add(owner.align.matrix[index][col]);
    	else {
    		left.findDescChars(col);
    		right.findDescChars(col);
    		for(int i = 0; i < left.descChars.size(); i++)
    			descChars.add(left.descChars.get(i));
    		for(int i = 0; i < right.descChars.size(); i++)
    			descChars.add(right.descChars.get(i));
    	}
    }
    
    /**
     * 
     */
    public void checkAncestral(int nchar){
    	int count = 0;
    	for(int i = 0; i < descChars.size(); i++)
			count += descChars.get(i) < GAP ? 1 : 0;
    	if(this == owner.root)	// initial call always to root, nchar is set here to check against other vertices
    		nchar = count;
    	if(nchar == 0)
    		firstVertex = Math.log(1 + Math.exp(survival) * (Math.exp(felsSum) - 1));
    	else{
    		if(count < nchar)
    			firstVertex = Utils.log0;
    		else
    			firstVertex = felsSum + survival;
    	}
		if(left != null){
			left.checkAncestral(nchar);
			right.checkAncestral(nchar);
		}
    }
    
    /**
     * 
     */
    public void calcPrior(){
    	if(this != owner.root)
    		priorFirst = Math.log( edgeLength / (owner.totalBranchLength + 1.0 / owner.owner.subModel.mu));
    	else
    		priorFirst = Math.log( 1.0 / owner.owner.subModel.mu / (owner.totalBranchLength + 1.0 / owner.owner.subModel.mu));
    }
    
    public void edgeChangeUpdate(){
    	subMatrix = owner.owner.subModel.calcSubMatrix(edgeLength);
    	owner.calcBranchLength();
    }
}


