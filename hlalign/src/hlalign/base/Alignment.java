package hlalign.base;

public class Alignment {
	AlignColumn first;
	
	int[][] matrix;
	
	int[] empty;
	
	public void print(){
		for(int i = 0; i < matrix.length; i++){
			for(int j = 0; j < matrix[0].length; j++)
				System.out.printf("%3d", matrix[i][j]);
			System.out.println();
		}
	}
}
