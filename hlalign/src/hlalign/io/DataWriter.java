package hlalign.io;

import java.io.FileWriter;
import java.io.IOException;

import hlalign.base.MCMC;
import hlalign.mcmc.McmcMove;

public class DataWriter {
	
	public MCMC owner;
	public FileWriter outputFile;
	public FileWriter ll;
	
	public DataWriter(MCMC mcmc){
		owner = mcmc;
	}
	
	public void writeTraces(){
		try{
			outputFile = new FileWriter("pipTraces.txt");
			int maxLength = 0;
			for (McmcMove mcmcMove : owner.scalarMoves) {
				outputFile.write(mcmcMove.name+"\t");
				maxLength = mcmcMove.sample.size() > maxLength ? mcmcMove.sample.size() : maxLength; 
			}	
			outputFile.write("\n");
			System.out.println(maxLength);
			for(int i = 0; i < maxLength; i++){
				for (McmcMove mcmcMove : owner.scalarMoves) {
					if(mcmcMove.sample.size() > i-2)
						outputFile.write(String.format("%.3f", mcmcMove.sample.get(i))+"\t");
					else
						outputFile.write("\t");
				}
				outputFile.write("\n");
			}
			outputFile.close();
			
			ll = new FileWriter("logLike.txt");
			for(int i = 0; i < owner.logLikeTrace.size(); i++)
				ll.write(String.format("%.3f", owner.logLikeTrace.get(i)) + "\n");
			ll.close();
		} catch (IOException e) {System.out.println("Exception writing samples");}
	}
}