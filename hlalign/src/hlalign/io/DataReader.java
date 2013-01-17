package hlalign.io;

import java.io.*;
// import java.lang.*;
import java.util.*;

public class DataReader {
	
	public double[][][] coords;
	
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
}
