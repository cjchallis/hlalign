package hlalign.base;


/**
 * This is a vertex of the tree.
 * The 'hardcore' functions are implemented in this class, developers are suggested
 * not change functions in it. The implemented functions are quite unreadable, since we
 * opted for efficiency and not readability. You should be able to develop novel
 * functionality of the software package (postprocessing, substitution models, etc.)
 * without touching this class.
 * @author miklos, novak
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

    /** The alignment between this vertex and its left/right descendant, given as the 
     * sequence of states through a pair HMM, with this vertex as the ancestor
     */
    public int[] alignL;
    public int[] alignR;
    
    public int descentOrder;
    
    int length;					// sequence length
//    AlignColumn first;			// first alignment column of Vertex
//    AlignColumn last;			// last, virtual alignment column of Vertex (always present)
    String seq;					// original sequence of this Vertex (given for leaves only)

    int winLength;                    // length of window
 //   AlignColumn winFirst;        // first alignment column of window
 //   AlignColumn winLast;        // first alignment column past window end
    boolean selected;                // shows if vertex is part of the selected subtree

    /** The length of the edge that connects this vertex with its parent. */
    public double edgeLength;                            // length of edge to parent vertex
    double[][] charTransMatrix;            // precalculated character transition likelihoods (subst. model)
    double[][] charPropTransMatrix;        // precalculated character transition likelihoods for proposals (subst. model)
    double[][] hmm2TransMatrix;            // precalculated state transition likelihoods for 2-seq HMM (indel model)
    double[][] hmm2PropTransMatrix;        // precalculated state transition likelihoods for 2-seq HMM used for proposals (indel model)
    double[][] hmm3TransMatrix;            // precalculated state transition likelihoods for 3-seq HMM (indel model)
    double[][] hmm3RedTransMatrix;    // precalculated st. trans. likelihoods for 3-seq HMM, silent st. removed (indel model)

    /**
     * The log-sum of the Felsenstein's likelihoods of characters that are inserted into the
     * sequence of this vertex.
     */
    public double orphanLogLike;        // log-sum of the likelihood of each orphan column in subtree (incl. this vertex)
    /**
     * The log-sum of the cumulative insertion-deletion loglikelihoods up to this vertex (ie. summed over the
     * subtree below this vertex.).
     */
    public double indelLogLike;

    public int leafCount;

    Vertex(int i, String nm) {
    	edgeLength = 1.0;
    	left = null;
    	right = null;
    	name = nm;
    	index = i;
    }
}