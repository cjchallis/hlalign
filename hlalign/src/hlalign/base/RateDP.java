package hlalign.base;

public class RateDP {
	HierarchicalIndel indel;
	
	/** insertion eligibility for super links
	 * 0: not eligible
	 * 1: eligible with surviving sibling (use \beta)
	 * 2: eligible with deceased sibling (use \gamma)
	 */
	int[] superEl;
	
	/** insertion eligibility for normal links, see above */
	int[] normalEl;
	
	/** order by which insertions appear in alignment:
	 * - left child before right (arbitrary choice)
	 * - child before parent (enforced by meaning of model)
	 * - lower number = higher priority
	 */
	int[] priority;
	
	/** number of steps from root, used to find root of a given alignment column */
	int[] seniority;
	
	/** full alignment of all nodes */
	int[][] align;
	
	
	public void forward(){
		// Initialize insertion eligibility arrays
		for(int i = 0; i < superEl.length; i++)
			superEl[i] = normalEl[i] = 1;
			
		
		// begin loop over align columns
		for(int i = 0; i < align[0].length; i++){
			//begin loop over states (s)
			for(int s = 0; s < indel.types; s++){
				// new align column
				int sroot = findRoot(i, true);
				// find super link root
				// if(super link root == old super link root)
					// calculate super link transition for staying in same super link
					// find normal link root
					// calculate transition to normal link
					// update eligibilities
					// calculate survival and death probabilities
					// update eligibilities // this doesn't need to be done for each state, will be the same for all
					// tr[s][s] = result
				// endif
				// else tr[s][s] = 0
		
				// get probability to Q state
				// new loop over states (t) (to calculate probabilities for new super link)
					// multiply probability to Q by probability from Q to super link root of type t
					// update eligibilities
					// calculate super- survival and death probabilities
					// update eligibilities
		
					// calculate transition to normal root
					// update eligibilities
					// calculate survival and death probabilities
					// update eligibilities
				// end for
				// tr[s][t] += result // sum of probabilities, use log_add for logs
			//end for
			}
			// for s
				// for t
					// DP[i][t] += tr[s][t]
		}
		//end for
	}
	
	public int findRoot(int column, boolean sup){
		int sen = align.length;
		int root = align.length;
		for(int i = 0; i < align.length; i++){
			if(align[i][column] > (sup ? -1 : 0)){
				if(seniority[i] < sen){
					sen = seniority[i];
					root = i;
				}
			}
		}
		return root;
	}
}
