package hlalign.io;

import java.io.*;
// import java.lang.*;
import java.util.*;

public class DataReader {
	
	public double[][][] coords;
	
	public char[][] seqs;
	
	public double[][] subQ;
	public double[] e;
	
	public DataReader(){}
	
	public void readAllCoords(String[] files){
		coords = new double[files.length][][];
		for(int i = 0; i < files.length; i++)
			try{
				coords[i] = readCoords(files[i]);
			} catch(Exception e){System.out.println("File '"+files[i]+ "' not found.");}
	}
	
	public double[][] readCoords(String file) 
			throws java.io.FileNotFoundException{
		Scanner input = new Scanner (new File(file));
		double[] temp = new double[1000];
		int i = 0;
		while(input.hasNextDouble()){
			temp[i] = input.nextDouble();
			i++;
		}
		input.close();
		double[][] x = new double [i / 3][3];
		for(int j = 0; j < i; j++)
			x[j / 3][j % 3] = temp[j];
		return x;
	}
	
	public void readAllSeqs(String[] files){
		seqs = new char[files.length][];
		for(int i = 0; i < files.length; i++)
			try{
				seqs[i] = readSeq(files[i]);
			} catch(Exception e){System.out.println("File '"+files[i]+ "' not found.");}
	}
	
	public char[] readSeq(String file)
			throws java.io.FileNotFoundException{
		String line;
		char[] seq = new char[0];;
		InputStream fis = new FileInputStream(file);
		BufferedReader br = new BufferedReader(new InputStreamReader(fis));
		try{
			while((line = br.readLine()) != null){
				seq = new char[line.length()];
				for(int i = 0; i < line.length(); i++)
					seq[i] = line.charAt(i);
			}
			br.close();
		}
		catch(Exception e){System.out.println("Exception during sequence reading");}	
		br = null;
		fis = null;
		return seq;
	}
	
	public void readSub()
			throws java.io.FileNotFoundException{
		BufferedReader bf;
		try{
			bf = new BufferedReader(new FileReader("data/alphabetA.dat"));
			String a = bf.readLine();
			char[] alphabet = new char[a.length()];
			for(int i = 0; i < alphabet.length; i++){
				alphabet[i] = a.charAt(i);
			}
			bf.close();
			bf = new BufferedReader(new FileReader("data/a.dat"));
			int size = alphabet.length;
			subQ = new double[size][size];
			e = new double[size];
			
			String[] temp;
			for(int i = 0; i < size; i++){
				temp = (bf.readLine()).split(" ");
				for(int j = 0; j < size; j++){
					subQ[i][j] = Double.parseDouble(temp[j]);
				}
			}
			bf.close();
			
			String t;
			bf = new BufferedReader(new FileReader("data/aPi.dat"));
			for(int i = 0; i < size; i++){
				t = bf.readLine();
				e[i] = Double.parseDouble(t);
			}
			bf.close();
		}
		catch(Exception e){System.out.println("Exception in reading substitution matrix"); 
		System.out.println(e.getClass());
		System.out.println(e.getMessage());}
	}
}
