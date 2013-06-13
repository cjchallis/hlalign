package hlalign.mcmc;

import hlalign.base.MCMC;
import hlalign.base.Tree;
import hlalign.base.Utils;
import hlalign.base.Vertex;

public class NearestNeighborInterchange extends McmcMove{
	
	// remember the index of vertex selected and which of left and right was swapped
	public int index, which;
	Tree tree;

	public NearestNeighborInterchange(MCMC mcmc, String n){
		owner = mcmc;
		name = n;
		tree = owner.tree;
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
		// decide whether to swap left or right with brother
		which = Utils.generator.nextInt(2);
		
		oldll = owner.getLogLike();
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
		
		// proposal ratio is always 1
		return 0;
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
		
		owner.setLogLike(oldll);
	}
	
}
