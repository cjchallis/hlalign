package hlalign.heterogeneous;

/** For representing multiple alignments in a column-free way */
public class Residue {
	/** surviving descendant of this residue on left branch */
	public Residue left;
	
	/** surviving descendant of this residue on right branch */
	public Residue right;
	
	/** surviving ancestor of this residue */
	public Residue parent;
	
	/** position of this residue within its protein */
	public int index;
	
	/** first residue inserted by this residue on left branch */
	public Residue insertLeft;
	
	/** first residue inserted by this residue on right branch */
	public Residue insertRight;
	
	/** residue which inserted this residue, if first insertion */
	public Residue insertParent;
	
	public Residue(int idx){
		left = right = parent = insertLeft = insertRight = insertParent = null;
		index = idx;
	}
}
