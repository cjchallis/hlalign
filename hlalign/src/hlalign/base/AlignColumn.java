package hlalign.base;

import java.util.ArrayList;

/** Contains an ArrayList<Integer> representing a column of the multiple alignment.
 * 1 for presence of residue, 0 for gap.  Index in the array is the vertex index.
 * @author Challis
 *
 */
public class AlignColumn {

	ArrayList<Integer> col;
	AlignColumn next;
	AlignColumn prev;

	public int size(){
		return col.size();
	}

	public void add(int k){
		col.add(k);
	}

	public int get(int k){
		return col.get(k);
	}
	public AlignColumn(){
		col = new ArrayList<Integer>();
		next = null;
		prev = null;
	}
}
    
