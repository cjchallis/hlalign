package hlalign.base;

import java.util.*;

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
    
    /** alpha-carbon coordinates */
    public double [][][] coords;
    
    public int[][] fullAlign;
    
    // match, delete, insert - 11, 10, 01 in binary
    public int[] emit = {3, 2, 1};
    
    /** indices of the leaf vertices */
    public int[] leaves = {3, 4, 5, 6};
    
    Alignment align;
    
    int[][] alignArray;
    
    public Tree(String[] pnames, double [][][] c){
    	coords = c;
    	names = pnames;
    	vertex = new Vertex[2*names.length - 1];
    	for(int i = 0; i < vertex.length; i++)
    		vertex[i] = new Vertex(i, null);
   
    	title = "Tree";
    	
    	// currently hard-coded topology
    	initTopology();
    	root.descentOrder = 0;
    	assignDescentOrder(root, 1);
    	// create coord array where coord[i] are the coordinates
    	// for the vertex with index i
//    	initCoords();
    	// currently hard-coded
    	initAlign();
 
    }
    
    public void initTopology(){
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
    		vertex[i].index = i;

    }
    
    public void initAlign(){
       	// Initialize alignment for testing
    	root.alignL = new int[]{0, 0, 2};
    	root.alignR = new int[]{0, 2, 1};
    	vertex[5].alignL = new int[]{2, 0, 1, 1};
    	vertex[5].alignR = new int[]{2, 1, 1, 1};
    	vertex[4].alignL = new int[]{1, 2, 0};
    	vertex[4].alignR = new int[]{2, 0, 0};
    	
    	printTree(root);
    	
    	createMultAlign();
    	printAlign();
    }
    
    /** For testing, prints array of 1's and 0's representing multiple alignment */
    public void printAlign(){
    	System.out.println("Print align");
    	AlignColumn ac = align.first.next;
    	System.out.println(align.first.size());
    	System.out.println(align.first.next.size());
    	System.out.println(ac.size());
    	for(int i = 0; i < 7; i++){
    		ac = align.first.next;
    		while(ac.next != null){
    			System.out.print(ac.get(i));
    			ac = ac.next;
    		}
    		System.out.println();
    	}
    	System.out.println(alignArray.length);
    	for(int i = 0; i < alignArray.length; i++){
    		for(int j = 0; j < alignArray[0].length; j++)
    			System.out.print(alignArray[i][j]);
    		System.out.println();
    	}
    	
    }
    
    /** Indices match order in which vertices are added to multiple alignment 
     * when composed from pairwise alignments
     * @param v
     * @param i
     * @return
     */
    public int assignDescentOrder(Vertex v, int i){
    	v.left.descentOrder = i;
    	i++;
    	v.right.descentOrder = i;
    	i++;
    	if(v.left.left != null){
    		i = assignDescentOrder(v.left, i);
    		i = assignDescentOrder(v.right, i);
    	}
    	return i;
    }
    
    /** Match coordinate indices to vertex indices */
    public void initCoords(){
    	// hard-coded position of input proteins
    	// input should be two pairs of closer-related proteins
    	double [][][] c = new double [7][][];
    	c[3] = coords[0];
    	c[4] = coords[1];
    	c[5] = coords[2];
    	c[6] = coords[3];
    	coords = c;
    	
    	String[] n = new String[7];
    	n[0] = "root";
    	n[3] = names[0];
    	n[4] = names[1];
    	n[5] = names[2];
    	n[6] = names[3];
    	names = n;
    }
    
    /** Recursive function to print across the tree.  Currently prints only vertex indices */
    public void printTree(Vertex v){
    	System.out.println(v.index);
    	if(v.left != null){
    		printTree(v.left);
    		printTree(v.right);
    	}
    }
    
    /** Composes pairwise alignments stored in vertices to multiple alignment
     * in form of chain of AlignColumns.
     */
    public void createMultAlign(){
    	align = new Alignment();
    	align.first = new AlignColumn();
    	align.first.next = new AlignColumn();
    	align.first.next.prev = align.first;
    	combinePairwise(root);
    	alignArray = alignToArray(align.first);
    }   
    
    public int[][] alignToArray(AlignColumn first){
    	System.out.println("Begin alignToArray");
    	AlignColumn ac = first.next;
    	int k = 0;
    	while(ac.size() > 0){
    		k++;
    		ac = ac.next;
    	}
    	ac = first.next;
    	System.out.print("ac size is " + ac.size());
    	int[][] alignArray = new int[ac.size()][k];
    	for(int j = 0; j < k; j++){
    		for(int i = 0; i < ac.size(); i++){
    			System.out.println("i " + i + " j " + j + " do " + vertex[i].descentOrder);
    			alignArray[i][j] = ac.get(vertex[i].descentOrder); 
    		}
    		ac = ac.next;
    	}
    	return alignArray;
    }
    
    /** Recursive function which combines the two pairwise alignments stored at each
     * internal node and appends to the multiple alignment.
     * @param v
     */
    private void combinePairwise(Vertex v){
    	int i = 0, j = 0;
    	AlignColumn ac = align.first.next;
    	while(i < v.alignL.length || j < v.alignR.length){
    		if(v != root){		// existing AlignColumns with 0 at the current vertex should have 0 at both children as well
    			while(ac.get(v.descentOrder) == 0){
    				ac.add(0);
    				ac.add(0);
    				ac = ac.next;
    			}
    		}
    		if(emit[v.alignL[i]] == 1 || emit[v.alignR[j]] == 1){	// if either alignment has an insertion, insert a new AlignColumn
    			AlignColumn insert = new AlignColumn();
    			insert.col = new ArrayList<Integer>(Collections.nCopies(ac.size(), 0));
    			insert.prev = ac.prev;
    			insert.next = ac;
    			ac.prev.next = insert;
    			ac.prev = insert;
    			
    			ac = insert;
    			if(v == root)
    				ac.add(0);
    			
    			if(emit[v.alignL[i]] == 1){		// check if insertion was left or right child
    				ac.add(1);
    				ac.add(0);
    				i++;
    			} else {
    				ac.add(0);
    				ac.add(1);
    				j++; 
    			}
    		} else{				// the next state is not an insertion in either alignment, so the residue is present in v
    			if(v == root)
    				ac.add(1);
    			ac.add( emit[v.alignL[i]] == 3 ? 1 : 0 );
    			ac.add( emit[v.alignR[j]] == 3 ? 1 : 0 );
    			i++;
    			j++;
    		}
    		
    		if(ac.next == null){	// move to the next AlignColumn and continue
    			ac.next = new AlignColumn();
    			ac.next.prev = ac;
    		}
    		ac = ac.next;
    	}
    	
    	while(ac.size() > 0){		// when the end of both alignments is reached, any additional AlignColumns 
    		ac.add(0);				// need to be filled out with zeroes
    		ac.add(0);
    		ac = ac.next;
    	}
    	
    	if(v.left.left != null){
    		combinePairwise(v.left);
    		combinePairwise(v.right);
    	}
    		
    }    
       
}