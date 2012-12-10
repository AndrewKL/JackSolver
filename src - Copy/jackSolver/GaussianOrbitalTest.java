package jackSolver;

import static org.junit.Assert.*;

import org.junit.Test;

public class GaussianOrbitalTest {

	@Test
	public void testToString() {
		System.out.println("===================================================");
		System.out.println("========ToString===================================");
		
		Atom atom1 = new Atom(0,0,1.463200,1);
		Atom atom2 = new Atom(0,0,0,2);
		Atom[] atoms = {atom2,atom1};
		System.out.println(" atom array:"+atoms);
		GaussianOrbital[] orbitals = JackRHFSolver.generateOrbitalArray(atoms,3);
		System.out.println(orbitals[0]);
		assertTrue(orbitals[0].toString().contains("GTO x: 0.0 y: 0.0 z: 0.0"));
	}

	@Test
	public void testCalculateR() {
		System.out.println("===================================================");
		System.out.println("========CalculateR=================================");
		Atom atom1 = new Atom(0,0,1.463200,1);
		Atom atom2 = new Atom(0,0,0,2);
		Atom[] atoms = {atom2,atom1};
		System.out.println(" atom array:"+atoms);
		GaussianOrbital[] orbitals = JackRHFSolver.generateOrbitalArray(atoms,3);
		System.out.println(" orbitals: "+orbitals);
		System.out.println(" atom1: "+atoms[0]+" atoms2: "+atoms[1]);
		double R = GaussianOrbital.calcR(orbitals[0], orbitals[1]);
		System.out.println( "R: "+R);
		//System.out.println( "end test");
		
		assertTrue(R==1.4632);
	}

	@Test
	public void testCalcRabsqrd() {
		System.out.println("===================================================");
		System.out.println("========CalculateRaRbdifsqrd=======================");
		
		Atom atom1 = new Atom(2,2,2,1);
		Atom atom2 = new Atom(0,0,0,2);
		Atom[] atoms = {atom2,atom1};
		System.out.println(" atom array:"+atoms);
		GaussianOrbital[] orbitals = JackRHFSolver.generateOrbitalArray(atoms,3);
		System.out.println(" orbitals: "+orbitals);
		System.out.println(" atom1: "+atoms[0]+" atoms2: "+atoms[1]);
		double R = GaussianOrbital.calcRabsqrd(orbitals[0], orbitals[1]);
		System.out.println( "R: "+R);
		assertTrue(R==12);
	}

	@Test
	public void testCalculateRpRcsqrd() {
		Atom atom1 = new Atom(1.463200,1.463200,1.463200,1);
		Atom atom2 = new Atom(0-1.463200,0-1.463200,-1.463200,1);
		Atom[] atoms = {atom2,atom1};
		System.out.println(" atom array:"+atoms);
		GaussianOrbital[] orbitals = JackRHFSolver.generateOrbitalArray(atoms,3);
		System.out.println(" orbitals: "+orbitals);
		System.out.println(" atom1: "+atoms[0]+" atoms2: "+atoms[1]);
		double Rpc = GaussianOrbital.calcRpcsqrd(orbitals[0], 0, orbitals[1], 0,2,2, 2);
		System.out.println( "Rpc: "+Rpc);
		assertTrue(Rpc==12);
	}

	@Test
	public void testCalculateRpqsqrd() {
		Atom atom1 = new Atom(0,0,0,1);
		Atom atom2 = new Atom(2,2,2,1);
		Atom[] atoms = {atom2,atom1};
		System.out.println(" atom array:"+atoms);
		GaussianOrbital[] orbitals = JackRHFSolver.generateOrbitalArray(atoms,3);
		System.out.println(" orbitals: "+orbitals);
		System.out.println(" atom1: "+atoms[0]+" atoms2: "+atoms[1]);
		double Rpq = GaussianOrbital.calcRpqsqrd(orbitals, 0, 0, 1, 1, 0, 0, 0, 0);
		System.out.println( "Rpc: "+Rpq);
		assertTrue(Rpq==12);
	}
	
	@Test
	public void testCalcRp() {
		System.out.println("===================================================");
		System.out.println("========CalculateR=================================");
		Atom atom1 = new Atom(1.463200,1.463200,1.463200,1);
		Atom atom2 = new Atom(0-1.463200,0-1.463200,-1.463200,1);
		Atom[] atoms = {atom2,atom1};
		System.out.println(" atom array:"+atoms);
		GaussianOrbital[] orbitals = JackRHFSolver.generateOrbitalArray(atoms,3);
		System.out.println(" orbitals: "+orbitals);
		System.out.println(" atom1: "+atoms[0]+" atoms2: "+atoms[1]);
		Vec3D Rp = GaussianOrbital.calcRp(orbitals[0], 0, orbitals[1], 0);
		System.out.println( "R: "+Rp);
		System.out.println( "end test");
		
		assertTrue(Rp.x==0&&Rp.y==0&&Rp.z==0);
	}
	
	

}
