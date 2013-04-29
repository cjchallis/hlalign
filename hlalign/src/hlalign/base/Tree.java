package hlalign.base;

import java.util.*;
import hlalign.indels.*;


public class Tree {

	/** The MCMC object containing the tree */
	public MCMC owner;
	
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
    
    /** The indel model*/
    IndelModel indels;
    
    /** The indel model for initializing alignments */
    IndelModel initIndel;
    
    /** The substitution model */
    SubstitutionModel subs;
    
    /** Maps Vertex.descentOrder to Vertex.index  */
    ArrayList<Integer> indexMap;
    
    public double mu;
    
    Alignment align;
    
    public int[][] alignArray;
    
    public Tree(String[] pnames, double [][][] c, MCMC mcmc){
    	owner = mcmc;
    	coords = c;
    	names = pnames;
    	vertex = new Vertex[2*names.length - 1];
    	for(int i = 0; i < vertex.length; i++){
    		vertex[i] = new Vertex(i, null);
    		vertex[i].owner = this;
    	}
   
    	title = "Tree";
    	
    	// currently hard-coded topology
    	initTopology();
    	
    	initIndel();
    	
    
    	// create coord array where coord[i] are the coordinates
    	// for the vertex with index i
// 		initCoords();
    	// currently hard-coded
    	initAlign();
 
    }
    
    public void initIndel(){
    	indels = new TKF91(new double[] {.03, .033});
    	indels.two.calcTrans(1);
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
    	// root.alignL = new int[]{0, 0, 2};
    	// root.alignR = new int[]{0, 2, 1};
    	// vertex[5].alignL = new int[]{2, 0, 1, 1};
    	// vertex[5].alignR = new int[]{2, 1, 1, 1};
    	// vertex[4].alignL = new int[]{1, 2, 0};
    	// vertex[4].alignR = new int[]{2, 0, 0};
    	
    	// printTree(root);
    

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
    

    public int[][] alignToArrayKeepOrder(AlignColumn first){
    	AlignColumn ac = first.next;
    	int k = 0;
    	while(ac.size() > 0){
    		k++;
    		ac = ac.next;
    	}
    	ac = first.next;
    	int[][] alignArray = new int[ac.size()][k];
  
    	for(int j = 0; j < k; j++){
    		for(int i = 0; i < ac.size(); i++)
    			alignArray[i][j] = ac.get(i); 
    		ac = ac.next;
    	}
    	return alignArray;
    }
 
}