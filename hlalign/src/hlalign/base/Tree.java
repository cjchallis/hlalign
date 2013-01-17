package hlalign.base;


/**
 * This is the current tree in the MCMC run.
 * It extends Stoppable, so the MCMC run can be terminated in graphical interface mode
 * while it calls a function in the Tree class.
 * @author miklos, novak
 */
public class Tree {

    /**
     * When two fragments of two sequences on two neighbor vertices of the tree
     * are realigned, the proposed new alignment is drawn from a forward-backward
     * sampling of this HMM. It is a simplified version of the TKF92 model.
     */
    // public Hmm2 hmm2;
    /**
     * When two fragments of two sequences on two sibling vertices are realigned
     * drawing a novel ancestral substring for their parent vertex, the proposed new
     * alignment is drawn from a forward-backward sampling of this HMM. It is
     * a pair-HMM having states emitting into the non-observable ancestral sequence.
     */
    // public HmmSilent hmm3;
    
	/** The array of vertices of the tree. */
    public Vertex vertex[];

    /** The root of the tree */
    public Vertex root;

    /** The name of the tree */
    public String title;

    /** The name of the sequences */
    public String[] names;
    
    public Tree(String[] pnames){
    	names = pnames;
    	vertex = new Vertex[2*names.length - 1];
    	for(int i = 0; i < vertex.length; i++)
    		if(i < names.length)
    			vertex[i] = new Vertex(i, names[i]);
    		else
    			vertex[i] = new Vertex(i, null);
   
    	title = "Tree";
    	System.out.println(vertex.length);
    	
    	
    	/* hard-coded tree topology for 4 proteins */
    	root = vertex[vertex.length - 1];
    	root.left = vertex[vertex.length - 2];
    	root.right = vertex[vertex.length - 3];
    	root.name = "root";
    	root.edgeLength = 0;
    	vertex[vertex.length - 2].left = vertex[0];
    	vertex[vertex.length - 2].right = vertex[1];
    	vertex[vertex.length - 3].left = vertex[2];
    	vertex[vertex.length - 3].right = vertex[3];
    	
    	vertex[2].parent = vertex[3].parent = vertex[4];
    	vertex[0].parent = vertex[1].parent = vertex[5];
    	vertex[4].parent = vertex[5].parent = vertex[6];
    	
    	for(int i = 0; i < vertex.length; i++)
    		System.out.println(vertex[i].index + " " + vertex[i].name + " " + vertex[i].edgeLength + " " 
    				+ vertex[i].parent + " " + vertex[i].left + " " + vertex[i].right);
    }
 
}