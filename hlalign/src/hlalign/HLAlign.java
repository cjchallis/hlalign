package hlalign;
import hlalign.io.DataReader;
import hlalign.base.MCMC;

public class HLAlign {
	
	public static void main(String[] args){
		System.out.println("Hierarchical Links Alignment v0.0");
		int n = args.length;
		String[] coordFiles = new String[n/2];
		String[] seqFiles = new String[n/2];
		
		for(int i = 0; i < n/2; i++){
			coordFiles[i] = args[i];
			seqFiles[i] = args[i + n/2];
		}
		
		DataReader datareader = new DataReader();
		datareader.readAllCoords(coordFiles);
		datareader.readAllSeqs(seqFiles);
		try{datareader.readSub();}
		catch(Exception e){System.out.println("Error reading substitutions.");}
		datareader.charToInt();
		
		
		MCMC mcmc = new MCMC(1000, 1000, datareader, coordFiles);
		mcmc.run();
		System.out.println("Done");
	}
}
