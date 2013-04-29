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
    
    	parentAlign(root);
    	
    	align = createMultAlign(root);
    	alignArray = alignToArray(align.first);
    	printAlign();
    }
    
 
   
    
    public void parentAlign(Vertex v){
    	if(v.left.left != null){
    		parentAlign(v.left);
    		parentAlign(v.right);
    	}
    	
    	int[][] leftInds = createPartialAlignInds(v.left);
    	
    	System.out.println("leftInds: " + v.left.index);
    	for(int i = 0; i < leftInds.length; i++){
    		for(int j = 0; j < leftInds[0].length; j++)
    			System.out.print(leftInds[i][j] + " ");
    		System.out.println();
    	}
    	
    	int[][] rightInds = createPartialAlignInds(v.right);
    	
    	System.out.println("rightInds: " + v.right.index);
    	for(int i = 0; i < rightInds.length; i++){
    		for(int j = 0; j < rightInds[0].length; j++)
    			System.out.print(rightInds[i][j] + " ");
    		System.out.println();
    	}
    	
    	int[][] alignInds = new int[leftInds.length][];
    	for(int i = 0; i < leftInds.length; i++)
    		alignInds[i] = Utils.copyOf(leftInds[i]);
    	
    	
    	
    	for(int i = 0; i < indexMap.size(); i++){
    		int j = indexMap.get(i);
    		System.out.println("iMap: " + j);
    		alignInds[j] = rightInds[j];
    	}
    	
    	System.out.println("alignInds: ");
    	for(int i = 0; i < alignInds.length; i++){
    		for(int j = 0; j < alignInds[i].length; j++)
    			System.out.print(alignInds[i][j] + " ");
    		System.out.println();
    	}
    	
    	Structure struc = owner.structure;
    	
    	double[] marginLeft = struc.marginTree(coords, leftInds, v.left.index, true, this);
    	double[] marginRight = struc.marginTree(coords, rightInds, v.right.index, true, this);
    	double[][] joint = struc.jointTree(coords, alignInds, v.index, true, this);
    	
    	TKF91 initIndel = new TKF91(new double[]{.03, .033});
    	initIndel.three.calcTrans(v.left.edgeLength, v.right.edgeLength);
    	
    	System.out.println("marginLeft: " + marginLeft.length);
    	System.out.println("marginRight: " + marginRight.length);
    	System.out.println("joint: " + joint.length + " x " + joint[0].length);
    	
    	PairDP dp = new PairDP(marginLeft, marginRight, joint, initIndel.three);
    	dp.forward2();
    	int[] triAlign = dp.backward();
    	
    	int[][] aligns = initIndel.translateAlign(triAlign);
    	v.alignL = aligns[0];
    	v.alignR = aligns[1];
    	
    }
    
    public int[][] createPartialAlignInds(Vertex v){
    	int[][] al3;
    	if(v.left != null){
    		Alignment al = createMultAlign(v);
    		int[][] al2 = alignToArrayKeepOrder(al.first);
    		al3 = new int[vertex.length][al2[0].length];
    		for(int i = 0; i < indexMap.size(); i++)
    			al3[indexMap.get(i)] = al2[i];
    	} else {
    		indexMap = new ArrayList<Integer>(0);
    		indexMap.add(v.index);
    		al3 = new int[vertex.length][coords[v.index].length];
    		for(int i = 0; i < al3[v.index].length; i++)
    			al3[v.index][i] = 1;
    	}
    	Structure struc = owner.structure;
    	System.out.println("Length al3: " + al3.length);
    	return struc.makeAlignInds(al3);
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
    
    /**
     * 
     */
    
    public void assignDescentOrder(Vertex v){
    	indexMap = new ArrayList<Integer>(0);
    	v.descentOrder = 0;
    	indexMap.add(v.index);
    	assignDescentOrder(v, 1);
    }
    
    /** Indices match order in which vertices are added to multiple alignment 
     * when composed from pairwise alignments
     * @param v
     * @param i
     * @return
     */
    public int assignDescentOrder(Vertex v, int i){
    	v.left.descentOrder = i;
    	indexMap.add(v.left.index);
    	i++;
    	v.right.descentOrder = i;
    	indexMap.add(v.right.index);
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
    public Alignment createMultAlign(Vertex v){
    	Alignment al = new Alignment();
    	al.first = new AlignColumn();
    	al.first.next = new AlignColumn();
    	al.first.next.prev = al.first;
    	assignDescentOrder(v);
    	combinePairwise(v, al);
    	return al;
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
    	
    	
    	for(int i = 0; i < vertex.length; i++)
    		System.out.println("Vertex " + vertex[i].index + " " + vertex[i].descentOrder);   		
    	
    	for(int i = 0; i < indexMap.size(); i++)
    		System.out.println("iMap: " + indexMap.get(i));
    	
    	for(int j = 0; j < k; j++){
    		for(int i = 0; i < ac.size(); i++){
    			System.out.println("i " + i + " j " + j + " do " + vertex[i].descentOrder);
    			alignArray[i][j] = ac.get(vertex[i].descentOrder); 
    		}
    		ac = ac.next;
    	}
    	return alignArray;
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
    
    /** Recursive function which combines the two pairwise alignments stored at each
     * internal node and appends to the multiple alignment.
     * @param v
     */
    private void combinePairwise(Vertex v, Alignment al){
    	int i = 0, j = 0;
    	AlignColumn ac = al.first.next;
    	while(i < v.alignL.length || j < v.alignR.length){
    		if(ac.size() > 0){		// existing AlignColumns with 0 at the current vertex should have 0 at both children as well
    			while(ac.get(v.descentOrder) == 0){
    				ac.add(0);
    				ac.add(0);
    				ac = ac.next;
    			}
    		}
    		if(indels.two.emitVal[v.alignL[i]] == 1 || indels.two.emitVal[v.alignR[j]] == 1){	// if either alignment has an insertion, insert a new AlignColumn
    			AlignColumn insert = new AlignColumn();
    			insert.col = new ArrayList<Integer>(Collections.nCopies(ac.size(), 0));
    			insert.prev = ac.prev;
    			insert.next = ac;
    			ac.prev.next = insert;
    			ac.prev = insert;
    			
    			ac = insert;
    			if(ac.size() == 0)
    				ac.add(0);
    			
    			if(indels.two.emitVal[v.alignL[i]] == 1){		// check if insertion was left or right child
    				ac.add(1);
    				ac.add(0);
    				i++;
    			} else {
    				ac.add(0);
    				ac.add(1);
    				j++; 
    			}
    		} else{				// the next state is not an insertion in either alignment, so the residue is present in v
    			if(ac.size() == 0)
    				ac.add(1);
    			ac.add( indels.two.emitVal[v.alignL[i]] == 3 ? 1 : 0 );
    			ac.add( indels.two.emitVal[v.alignR[j]] == 3 ? 1 : 0 );
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
    		combinePairwise(v.left, al);
    		combinePairwise(v.right, al);
    	}
    		
    }    
       
}