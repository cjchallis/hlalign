package hlalign.mcmc;

import hlalign.base.MCMC;
import hlalign.base.Tree;
import hlalign.base.Utils;
import hlalign.base.Vertex;

public class StochasticNNI extends McmcMove{
	
	// remember the index of vertex selected and
	// *which* of the three topologies was proposed
	public int index, which;
	Tree tree;
	
	// tuning parameter for branch length multipliers
	double lambda;
	// 5 branch length multipliers
	double[] m;

	public StochasticNNI(MCMC mcmc, String n){
		owner = mcmc;
		name = n;
		tree = owner.tree;
		lambda = 2*Math.log(1.2);
	}
	
	@Override
	public void copyState(Object externalState) {
		// draw uniformly from branches not attached to leaves
		// also exclude right branch of the root, because selecting either root branch
		// will result in same move
		int nLeaves = tree.names.length;
		// should be impossible to select leaf or root
		index = Utils.generator.nextInt(tree.vertex.length - nLeaves - 1) + nLeaves;
		// now make sure it is not right child of root
		while(tree.root.right == tree.vertex[index])
			index = Utils.generator.nextInt(tree.vertex.length - nLeaves - 1) + nLeaves;
		// decide whether to swap left (0) or right (1) with brother, or stay the same (2)
		which = Utils.generator.nextInt(3);
		
		oldll = owner.getLogLike();
		
		m = new double[5];
		for(int i = 0; i < m.length; i++)
			m[i] = Math.exp(lambda*(Utils.generator.nextDouble() - 0.5));
	}

	@Override
	public double proposal(Object externalState) {
		/* Tree under consideration for nearest neighbor
		 * grandpa, parent, brother, vertex, right, left
		 * 
		 *		g
		 *	 	|
		 * 		p -- b
		 * 		|
		 * 		v -- r
		 * 		|
		 * 		l
		 * 
		 *	Note that if the root is in the tree it is on either p--v or g--p
		 *
		 *	We will always swap brother with either right or left, this grandpa
		 *	need not be identified and we only need to worry about the root on p--v
		 */
		if(Utils.DEBUG){
			System.out.println("vert is " + index);
			System.out.println("Before NNI:");
			tree.printTree(tree.root);
		}
		
		Vertex vert = tree.vertex[index];
		Vertex left = vert.left;
		Vertex right = vert.right;
		Vertex parent, brother;
		if(vert.parent == tree.root){
			parent = tree.root.right;
			brother = parent.right;
		}
		else{
			parent = vert.parent;
			if(vert == parent.left)
				brother = parent.right;
			else
				brother = parent.left;
		}
		
		left.edgeLength *= m[0];
		right.edgeLength *= m[1];
		vert.edgeLength *= m[2];
		if(vert.parent == tree.root){
			tree.root.right.edgeLength *= m[2];
			brother.edgeLength *= m[3];
			parent.left.edgeLength *= m[4];
		}
		else {
			brother.edgeLength *= m[3];
			parent.edgeLength *= m[4];
			if(tree.root.left == parent)
				tree.root.right.edgeLength *= m[4];
			else if(tree.root.right == parent)
				tree.root.left.edgeLength *= m[4];
		}
				
		// if which == 2, don't change the topology at all
		if(which < 2){
			Vertex swap = (which == 0 ? left : right);
			// parent always becomes the parent of swap
			// and vert always becomes parent of brother
			swap.parent = parent;
			brother.parent = vert;
			if(brother == parent.right)
				parent.right = swap;
			else
				parent.left = swap;		
			if(which == 0)
				vert.left = brother;
			else
				vert.right = brother;
		}
		if(Utils.DEBUG){
			System.out.println("After NNI:");
			tree.printTree(tree.root);
		}
		// proposal ratio product of m
		double logProp = 0;
		for(int i = 0; i < m.length; i++)
			logProp += Math.log(m[i]);
		return logProp;
	}

	@Override
	public double logPriorDensity(Object externalState) {
		// uniform prior over trees
		return 0;
	}

	@Override
	public void updateLikelihood(Object externalState) {
		owner.setLogLike(owner.tree.calcML());
	}
	
	@Override
	public boolean isParamChangeAccepted(double logProposalRatio){
		boolean accept = getOwner().isParamChangeAccepted(logProposalRatio,this);
		if(accept & which == 2)
			acceptanceCount--;
		return accept; 
	}

	@Override
	public void restoreState(Object externalState) {
		Vertex vert = tree.vertex[index];
		Vertex left = vert.left;
		Vertex right = vert.right;
		Vertex parent, brother;
		if(vert.parent == tree.root){
			parent = tree.root.right;
			brother = parent.right;
		}
		else{
			parent = vert.parent;
			if(vert == parent.left)
				brother = parent.right;
			else
				brother = parent.left;
		}
		
		Vertex swap = (which == 0 ? left : right);
		
		if(which < 2){
			// parent always becomes the parent of swap
			// and vert always becomes parent of brother
			swap.parent = parent;
			brother.parent = vert;
			if(brother == parent.right)
				parent.right = swap;
			else
				parent.left = swap;		
			if(which == 0)
				vert.left = brother;
			else
				vert.right = brother;
		}
		
		vert.left.edgeLength /= m[0];
		vert.right.edgeLength /= m[1];
		vert.edgeLength /= m[2];
		
		if(vert.parent == tree.root){
			tree.root.right.edgeLength /= m[2];
			tree.root.right.right.edgeLength /= m[3];
			tree.root.right.left.edgeLength /= m[4];
		}
		else {
			if(parent.left == vert)
				parent.right.edgeLength /= m[3];
			else
				parent.left.edgeLength /= m[3];
			parent.edgeLength /= m[4];
			if(tree.root.left == parent)
				tree.root.right.edgeLength /= m[4];
			else if(tree.root.right == parent)
				tree.root.left.edgeLength /= m[4];
		}
		
		owner.setLogLike(oldll);
		if(Utils.DEBUG){
			System.out.println("Restored tree: ");
			tree.printTree(tree.root);
		}
	}
	
}


