package jackSolver;

public class Atom{
	/*
	 * JackRHFSolver
	 * 
	 * by Andrew Long 
	 * Advisor: Jason A.C. Clyburne 
	 * with additional assistance from Cory Pye.
	 * 
	 * based on the Fortran IV program by attila Szabo and Neil S. Ostlund
	 * in the book Modern Quantum Chemisty: introduction to advanced electronic 
	 * structure theory/ Attila Szabo, Neil S. Ostlund
	 *  
	 *  If you're just learning this stuff it is well worth your time to get the 
	 *  book. It's excellent.  The steps in the comments in the main part of this 
	 *  program match up with the steps on page 146 of the book.
	 *  
	 *  
	 */
	
	public double x,y,z;
	int element;
	
	public Atom(double xin,double yin,double zin,int elementin){
		x=xin;
		y=yin;
		z=zin;
		element=elementin;
	}
	
	public String toString(){
		return " atom x: "+x+" y: "+y+" z: "+z+" Za: "+element;
	}
	
	public static double distance(Atom a,Atom b){
		
		return Math.sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z));
		
	}
}