package hlalign.base;


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

    /**
     * The name of the sequence associated to the vertex.
     * Used only if the vertex is a leaf of the tree.
     */
    public String name;
    
    /** Index of this vertex in the Tree.vertex array */
    public int index;

    /** This reference points to the parent of the vertex */
    public Vertex parent;
    /**
     * This reference points to the left child of the vertex. If the
     * vertex is a leaf, it is set to null.
     */
    public Vertex left;
    /**
     * This reference points to the right child of the vertex. If the
     * vertex is a leaf, it is set to null.
     */
    public Vertex right;

    /** The length of the edge that connects this vertex with its parent. */
    public double edgeLength;          			// length of edge to parent vertex

    int length;									// sequence length
    String seq;									// original sequence of this Vertex (given for leaves only)

    /** alignments between Left and Right children */
    public int[] alignL;
    public int[] alignR;
    
    public int descentOrder;
    
    // the following are all column-specific and relate to the Poission Indel Process
    public double survival;						// survival probability
    public double commonAncestor;				// if a common ancestor of the column, Utils.log0 if not
    public double[] fels;						// felsenstein vector, \tilde{f}_v(\sigma)
    public double felsSum;						// weighted sum by Pi over felsenstein vector of this vertex, \tilde{f}_v
    public double firstVertex;					// probability of alignment column given that this was the first vertex
    												// after insertion of the original character, f_v
    public double[][] subMatrix;				// (log) substitution probabilities from parent to this vertex

    public int leafCount;

    Vertex(int i, String nm) {
    	edgeLength = 1.0;
    	left = null;
    	right = null;
    	name = nm;
    	index = i;
    }
    
    /**
     * Calculate the felsenstein vector for this vertex (for a particular column)
     * from its children
     */
    public void felsenstein(){   	
    	for(int i = 0; i < fels.length; i++){
    	   	double l = Utils.log0;
        	double r = Utils.log0;
    		for(int j = 0; j < fels.length; j++){
    			l = Utils.logAdd(l, left.subMatrix[i][j] + left.fels[j]);
    			r = Utils.logAdd(r, right.subMatrix[i][j] + right.fels[j]);
    		}
    		fels[i] = l + r;
    	}  	
    }
    
    /**
     * Calculate the pi-weighted felsenstein sum, equivalent to substitution likelihood
     * if the alignment column originated at this vertex
     * This quantity is \tilde{f}_v from the PIP paper
     */
    public void calcFelsSum(){
    	felsSum = Utils.log0;
    	for(int i = 0; i < fels.length; i++)
    		felsSum = Utils.logAdd(felsSum, owner.subs.pi[i] + fels[i]);
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
    
    public void calcSurvival(){
    	survival = -Math.log(edgeLength) - Math.log(owner.mu) + Math.log(1 - Math.exp(edgeLength * owner.mu));
    }
}