/**
 * 
 */
package jackSolver;

import static org.junit.Assert.*;

import org.junit.Test;

import cern.colt.*;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import cern.colt.matrix.linalg.SingularValueDecomposition;

/**
 * @author Andrew
 *
 */
public class JackRHFSolverTest {
	
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


	/**
	 * Test method for {@link jackSolver.JHFSolver#symetricOrthogonalization(cern.colt.matrix.DoubleMatrix2D)}.
	 */
	/*@Test
	public void testSymetricOrthogonalization() {
		System.out.println("===================================================");
		System.out.println("========testSymetricOrthogonalizaton===============");
		
		double[][] Smat = {{1.00,0.4507704116},{0.4507704116,1.0}};
		DoubleMatrix2D SMatrix = new DenseDoubleMatrix2D(Smat);
		
		DoubleMatrix2D XTransformMatrix = JHFSolver.symetricOrthogonalization(SMatrix);
		System.out.println("===XTransformMatrix===");
		System.out.println(XTransformMatrix);
		
		fail("Not yet implemented");
	}*/
	
	/**
	 * Test method for {@link jackSolver.JackRHFSolver#adjoint(cern.colt.matrix.DoubleMatrix2D)}.
	 */
	/*@Test
	public void testDiagonalizeMatrix() {
		System.out.println("===================================================");
		System.out.println("========testDiagonalizeMatrix======================");
		// assume that A = P^-1*D*P
		
		//double[][] Smat = {{1.00,0.4507704116},{0.4507704116,1.0}};
		double[][] Smat = {{-1,-1,1},{0,-2,1},{0,0,-1}};
		DoubleMatrix2D SMatrix = new DenseDoubleMatrix2D(Smat);
		System.out.println(SMatrix);
		DoubleMatrix2D Sdiagonalized;
		DoubleMatrix2D UnitaryMatrixforS;
		
		
		EigenvalueDecomposition evd = new EigenvalueDecomposition(SMatrix);
		Sdiagonalized = evd.getD();
		//System.out.println(evd);
		System.out.println(Sdiagonalized);
		
		
		UnitaryMatrixforS = evd.getV();
		System.out.println(" "+UnitaryMatrixforS);
		
		System.out.println();
		System.out.println("===test 2x2 case===");
		double[][] Smat2 = {{1.00,0.4507704116},{0.4507704116,1.0}};
		
		//double [][] diagonalized = diagonalize2x2Matrix(Smat);
		SMatrix = new DenseDoubleMatrix2D(Smat2);
		
		
		DoubleMatrix2D diagonalizedMatrix = JHFSolver.DiagonalizeMatrix(SMatrix);
		System.out.println(diagonalizedMatrix);
		
		
		fail("Not yet implemented");
	}*/
	
	
	
	@Test
	public void testBookDiagonalized() {
		System.out.println("===================================================");
		System.out.println("========testdiagonalize2x2Matrix===================");
		System.out.println((Math.PI/4));
		System.out.println("===================================================");
		
		double[][] Smat = {{1.00,0.4507704116},{0.4507704116,1.0}};
		
		//double [][] diagonalized = diagonalize2x2Matrix(Smat);
		DoubleMatrix2D SMatrix = new DenseDoubleMatrix2D(Smat);
		
		
		DoubleMatrix2D DMatrix = JackRHFSolver.diagonalize2x2MatrixgetEigenValues(SMatrix);
		DoubleMatrix2D PMatrix = JackRHFSolver.diagonalize2x2MatrixgetEigenVectors(SMatrix);
		
		System.out.println(SMatrix);
		System.out.println(DMatrix);
		System.out.println(PMatrix);
		
		//assertTrue(diagonalizedMatrix.hashCode()==689021860);
		System.out.println("===================================================");
		double[][] mat1 = {{-4,-6},{3,5}};
		
		//double [][] diagonalized = diagonalize2x2Matrix(Smat);
		SMatrix = new DenseDoubleMatrix2D(mat1);
		
		System.out.println(SMatrix);
		DMatrix = JackRHFSolver.diagonalize2x2MatrixgetEigenValues(SMatrix);
		System.out.println(DMatrix);
		PMatrix = JackRHFSolver.diagonalize2x2MatrixgetEigenVectors(SMatrix);
		System.out.println(PMatrix);
		
		
		System.out.println("===================================================");

		double[][] mat2 = {{1,3},{3,4}};
		
		
		SMatrix = new DenseDoubleMatrix2D(mat2);
		
		
		System.out.println(SMatrix);
		DMatrix = JackRHFSolver.diagonalize2x2MatrixgetEigenValues(SMatrix);
		System.out.println(DMatrix);
		PMatrix = JackRHFSolver.diagonalize2x2MatrixgetEigenVectors(SMatrix);
		System.out.println(PMatrix);
		
		System.out.println("===================================================");

		double[][] mat3 = {{-2.439732411,-0.5158386047},{-0.5158386047,-1.538667186}};
		
		
		SMatrix = new DenseDoubleMatrix2D(mat3);
		
		
		System.out.println(SMatrix);
		DMatrix = JackRHFSolver.diagonalize2x2MatrixgetEigenValues(SMatrix);
		System.out.println(DMatrix);
		PMatrix = JackRHFSolver.diagonalize2x2MatrixgetEigenVectors(SMatrix);
		System.out.println(PMatrix);
		
		
	}
	
	
	
	@Test
	public void testDiagonalizefprime() {
		System.out.println("===================================================");
		System.out.println("========testdiagonalize2x2Matrix===================");
		double[][] fPrimemat = {{-2.43973,-0.515838},{-0.515838,-1.538663}};
		
		
		
		DoubleMatrix2D fPrimeMatrix = new DenseDoubleMatrix2D(fPrimemat);
		
		
		DoubleMatrix2D diagonalizedMatrix = JackRHFSolver.diagonalize2x2MatrixgetEigenValues(fPrimeMatrix);
		System.out.println(diagonalizedMatrix);
		System.out.println(diagonalizedMatrix.hashCode());
		
	}
	
	
	
	@Test
	public void testCOLTDiagonalized() {
		System.out.println("===================================================");
		System.out.println("========testcoltdiagonalize========================");
		double[][] Smat = {{-2.43973,-0.515838},{-0.515838,-1.538663}};
		DoubleMatrix2D SMatrix = new DenseDoubleMatrix2D(Smat);
		
		System.out.println("===SMatrix=====");
		System.out.println(SMatrix);
		
		
		
		EigenvalueDecomposition evd = new EigenvalueDecomposition(SMatrix);
		
		System.out.println("===diagonalized S=====");
		System.out.println(evd);
		System.out.println("======================");
		//if(DEBUG)System.out.println(evd.getD);

		
		
		DoubleMatrix2D diagonalizedMatrix = evd.getD();
		System.out.println(diagonalizedMatrix);

	}

	
	/**
	 * Test method for {@link jackSolver.JackRHFSolver#EstimateDenstiyMatrix(jackSolver.GaussianOrbital[], int)}.
	 */
	@Test
	public void testEstimateDenstiyMatrix() {
		System.out.println("===================================================");
		System.out.println("========testEstimateDenstiyMatrix==================");
		Atom atom1 = new Atom(0,0,1.463200,1);
		Atom atom2 = new Atom(0,0,0,2);
		Atom[] atoms = {atom2,atom1};
		GaussianOrbital[] orbitals = JackRHFSolver.generateOrbitalArray(atoms,3);
		
		
		
		DoubleMatrix2D densityMatrix = JackRHFSolver.EstimateDenstiyMatrix(orbitals, 2);
		boolean iseveryelementzero = true;
		for(int i=0;i<densityMatrix.rows();i++){
			for(int j=0;j<densityMatrix.columns();j++){
				if(densityMatrix.get(i,j)!=0){
					iseveryelementzero=false;
				}
			}
		}
		
		
		assertTrue(iseveryelementzero);
		
	}
	
	/**
	 * Test method for {@link jackSolver.JackRHFSolver#inversesqrt(cern.colt.matrix.DoubleMatrix2D)}.
	 */
	@Test
	public void testInversesqrt() {
		System.out.println("===================================================");
		System.out.println("========testInversesqrt============================");
		double[][] mat = {{1,0,0},{0,4,0},{0,0,0.25}};
		DoubleMatrix2D matrix = new DenseDoubleMatrix2D(mat);
		System.out.println(matrix);
		JackRHFSolver.inversesqrt(matrix);
		System.out.println(matrix);
		
		double[][] answerkey={{1,0,0},{0,0.5,0},{0,0,2}};
		DoubleMatrix2D answerkeymatrix = new DenseDoubleMatrix2D(answerkey);
		System.out.println(answerkeymatrix.hashCode());
		
		System.out.println(answerkeymatrix);
		assertTrue(matrix.equals(answerkeymatrix));
	}
	
	/**
	 * Test method for {@link jackSolver.JackRHFSolver#cannonicalOrthogonalization(DoubleMatrix2D SMatrix)}.
	 */
	@Test
	public void testcannonicalOrthogonalization(){
		System.out.println("===================================================");
		System.out.println("========testcannonicalOrthogonalization============");
		
		double[][] Smat = {{1.00,0.4507704116},{0.4507704116,1.0}};
		DoubleMatrix2D SMatrix = new DenseDoubleMatrix2D(Smat);
		
		System.out.println("===SMatrix===");
		System.out.println(SMatrix);
		
		DoubleMatrix2D XMatrix = JackRHFSolver.cannonicalOrthogonalization(SMatrix);
		System.out.println("===XMatrix===");
		System.out.println(XMatrix);
		double value = XMatrix.hashCode();
		System.out.println(value);
		

		
	}
	
	/**
	 * Test method for {@link jackSolver.JackRHFSolver#calculate4centerintegrals(jackSolver.GaussianOrbital[])}.
	 */
	@Test
	public void testCalculate4centerintegrals() {
		System.out.println("===================================================");
		System.out.println("======== testCalculate4centerintegrals=============");
		//boolean passtest=true;
		Atom atom1 = new Atom(0,0,1.463200,1);
		Atom atom2 = new Atom(0,0,0,2);
		Atom[] atoms = {atom2,atom1};
		int STOnG =1;  
		GaussianOrbital[] orbitals = JackRHFSolver.generateOrbitalArray(atoms,STOnG);
		double[][][][] fourcenterintegrals = JackRHFSolver.calculate4centerintegrals(orbitals);
		
		for(int i=0;i<fourcenterintegrals.length;i++){
			for(int j=0;j<fourcenterintegrals[i].length;j++){
				for(int k=0;k<fourcenterintegrals[i][j].length;k++){
					for(int l=0;l<fourcenterintegrals[i][j][k].length;l++){
						
						System.out.println("4CI ( "+(1+i)+" "+(1+j)+" "+(1+k)+" "+(1+l)+" ) value: "+fourcenterintegrals[i][j][k][l]);
					}
				}
			}
			
		}
		
		
		
	}
	
	/**
	 * Test method for {@link jackSolver.JackRHFSolver#calculate4centerintegral(jackSolver.GaussianOrbital[])}.
	 */
	@Test
	public void testCalculate4centerintegral() {
		System.out.println("===================================================");
		System.out.println("======== testCalculate4centerintegral=============");
		
		boolean passtest=true;
		
		Atom atom1 = new Atom(0,0,1.463200,1);
		Atom atom2 = new Atom(0,0,0,2);
		Atom[] atoms = {atom2,atom1};
		int STOnG =1;  
		GaussianOrbital[] orbitals = JackRHFSolver.generateOrbitalArray(atoms,STOnG);
		
		//=======================================================
		double fourci = JackRHFSolver.calculateFourCenterIntegral(orbitals[0].alpha[0], orbitals[0].alpha[0], orbitals[0].alpha[0], orbitals[0].alpha[0], 0, 0, 0);
		double d = orbitals[0].D[0]*orbitals[0].D[0]*orbitals[0].D[0]*orbitals[0].D[0];
		double total = fourci*d;
		System.out.println(" alpha: "+orbitals[0].alpha[0]+" beta: "+orbitals[0].alpha[0]+" gamma: "+orbitals[0].alpha[0]+" delta: "+orbitals[0].alpha[0]);
		System.out.println(" fourci: "+fourci);
		System.out.println(" d: "+d);
		System.out.println(" total: "+total);
		
		if(Math.abs(total-1.229037413723358)>=0.0001)passtest=false;
		
		System.out.println("=========================");
		
		fourci = JackRHFSolver.calculateFourCenterIntegral(orbitals[1].alpha[0], orbitals[1].alpha[0], orbitals[1].alpha[0], orbitals[1].alpha[0], 0, 0, 0);
		d = orbitals[1].D[0]*orbitals[1].D[0]*orbitals[1].D[0]*orbitals[1].D[0];
		total = fourci*d;
		System.out.println("  "+orbitals[1].alpha[0]+ orbitals[1].alpha[0]+ orbitals[1].alpha[0]+ orbitals[1].alpha[0]+ 0+ 0+ 0);
		System.out.println(" alpha: "+orbitals[1].alpha[0]+" beta: "+orbitals[1].alpha[0]+" gamma: "+orbitals[1].alpha[0]+" delta: "+orbitals[1].alpha[0]);
		System.out.println(" fourci: "+fourci);
		System.out.println(" d: "+d);
		System.out.println(" total: "+total);
		
		if(Math.abs(total-0.7283184673916195)>=0.0001)passtest=false;
		
		assertTrue(passtest);
	}
	
	@Test
	public void testCalculate4centerintegralcontractedgaussianloop() {
		System.out.println("===================================================");
		System.out.println("======== testCalculate4centerintegral=============");
		
		boolean passtest =true;
		
		Atom atom1 = new Atom(0,0,1.463200,1);
		Atom atom2 = new Atom(0,0,0,2);
		Atom[] atoms = {atom2,atom1};
		int STOnG =3;  
		GaussianOrbital[] orbitals = JackRHFSolver.generateOrbitalArray(atoms,STOnG);
		double fourci = JackRHFSolver.calculateFourCenterIntegrals(orbitals,1,1,1,1);
		
		System.out.println("==================");
		System.out.println(" Total fourci: "+fourci);
		
		if(Math.abs(fourci-0.7746083600328786)>=0.0001)passtest=false;
		
		
		assertTrue(passtest);
	}
	
	
	
	/**
	 * Test method for {@link jackSolver.JackRHFSolver#generateOrbitalArray(Atom[] atoms)}.
	 */
	@Test
	public void testgenerateOrbitalArray() {
		
		System.out.println("===================================================");
		System.out.println("========generateOrbitalArray=======================");
		Atom atom1 = new Atom(0,0,1.463200,1);
		Atom atom2 = new Atom(0,0,0,2);
		Atom[] atoms = {atom2,atom1};
		GaussianOrbital[] orbitals = JackRHFSolver.generateOrbitalArray(atoms,3);
		boolean passtest =true;
		for(int i=0;i<orbitals.length;i++){
			System.out.println("orbital i: "+i);
			if(orbitals[i]==null)passtest=false;
			else{for(int j=0;j<orbitals[i].D.length;j++){
				System.out.print(" d: "+orbitals[i].D[j]);
			}}
		}
		
		assertTrue(passtest);
	}


	
	/**
	 * Test method for {@link jackSolver.JackRHFSolver#calcNewDensityMatrix(PDensityMatrix,cCoefficientsMatrix,numberofelectrons)}.
	 */
	@Test
	public void testcalcNewDensityMatrix() {
		System.out.println();
		System.out.println("===================================================");
		System.out.println("========testcalcNewDensityMatrix====================");
		double[][] pmat = {{0,0},{0,0}};
		
		DoubleMatrix2D PMatrix = new DenseDoubleMatrix2D(pmat);
		System.out.println(PMatrix);
		
		/*double[][] cmat = {{0.929145, -0.625857},{0.139834, 1.11151}};
		DoubleMatrix2D CMatrix = new DenseDoubleMatrix2D(cmat);
		System.out.println(CMatrix);*/


 
		
		
		
		double[][] Coefmat = {{.9291467304,-0.6258569539},{0.1398330503,0.1111511265}};
		
		DoubleMatrix2D CMatrix = new DenseDoubleMatrix2D(Coefmat);
		System.out.println(CMatrix);
		
		JackRHFSolver.calcNewDensityMatrix(PMatrix,CMatrix,2);
		System.out.println(PMatrix);
		
		// 1.726627   0.25985084
		// 0.2598508  0.0391065639
		
		
	}

	

	/**
	 * Test method for {@link jackSolver.JackRHFSolver#calcGMatrix(DoubleMatrix2D gMatrix,
			DoubleMatrix2D pDensityMatrix, GaussianOrbital[] orbitals)}.
	 */
	@Test
	public void testcalcGMatrix() {
		System.out.println("===================================================");
		System.out.println("========testcalcGMatrix============================");
		double[][] pmat = {{1.726627, 0.25985084},{0.2598508, 0.0391065639}};
		
		
		DoubleMatrix2D PMatrix = new DenseDoubleMatrix2D(pmat);
		System.out.println(PMatrix);
		
		double[][] gmat = {{0,0},{0,0}};
		DoubleMatrix2D GMatrix = new DenseDoubleMatrix2D(gmat);
		System.out.println(GMatrix);
		
		Atom atom1 = new Atom(0,0,1.463200,1);
		Atom atom2 = new Atom(0,0,0,2);
		Atom[] atoms = {atom2,atom1};
		GaussianOrbital[] orbitals = JackRHFSolver.generateOrbitalArray(atoms,3);
		
		
		JackRHFSolver.calcGMatrix(GMatrix,PMatrix, orbitals, null);
		System.out.println(GMatrix);
		
	}

	/**
	 * Test method for {@link jackSolver.JHFSolver#calcConvergence(DoubleMatrix2D pDensityMatrix,
			DoubleMatrix2D oldpDensityMatrix)}.
	 */
	@Test
	public void testcalcConvergence() {
		System.out.println("===================================================");
		System.out.println("========testtestcalcConvergence====================");
		
		double[][] pmat = {{1.726627, 0.25985084},{0.2598508, 0.0391065639}};
		
		
		DoubleMatrix2D PMatrix = new DenseDoubleMatrix2D(pmat);
		System.out.println(PMatrix);
		
		double[][] oldpmat = {{0,0},{0,0}};
		DoubleMatrix2D oldpMatrix = new DenseDoubleMatrix2D(oldpmat);
		System.out.println(oldpMatrix);
		
		double delta = JackRHFSolver.calcConvergence(PMatrix,oldpMatrix);
		System.out.println(" delta: "+delta);
		
		assertTrue(Math.abs(delta-0.882866781821925)<=0.0001);
		
		
	}

	/**
	 * Test method for {@link jackSolver.JHFSolver#calcNuclearRepulsionEnergy(Atom[] atoms)}.
	 */
	@Test
	public void testcalcNuclearRepulsionEnergy() {
		System.out.println("===================================================");
		System.out.println("========testcalcNuclearRepulsionEnergy=============");
		Atom atom1 = new Atom(0,0,1.463200,1);
		Atom atom2 = new Atom(0,0,0,2);
		Atom[] atoms = {atom2,atom1};
		
		double nuclearRepulsionEnergy = JackRHFSolver.calcNuclearRepulsionEnergy(atoms);
		System.out.println(" nuclear repulsion energy: "+nuclearRepulsionEnergy);
		
		assertTrue(Math.abs(nuclearRepulsionEnergy-1.366867140513942)<=0.0001);
		
		
		System.out.println("===================================================");
		
		atom1 = new Atom(0,0,1.463200,1);
		atom2 = new Atom(0,0,0,2);
		nuclearRepulsionEnergy = JackRHFSolver.calcNuclearRepulsionEnergy(atoms);
		System.out.println(" nuclear repulsion energy: "+nuclearRepulsionEnergy);
		
	}

	/**
	 * Test method for {@link jackSolver.JHFSolver#HFSolver(Atom[] atoms, int numberofelectrons)}.
	 */
	@Test
	public void testHFSolver() {
		System.out.println("=================================================================");
		
		System.out.println("calculation RHF for H-He+ ");
		Atom atom1 = new Atom(0,0,1.463200,1);
		Atom atom2 = new Atom(0,0,0,2);
		//Atom atom3 = new Atom(0,2,2,1);
		//Atom atom4 = new Atom(0,2,3,1);
		Atom[] atoms = {atom2,atom1};
		
		int numberofelectrons =2;
		System.out.println("Calculating H-He energy");
		
		double energy = JackRHFSolver.HFSolver(atoms,numberofelectrons);
		System.out.println(" energy: "+energy );
		assertTrue(Math.abs(energy+2.86065870404886)<=0.000001);
	}
	/**
	 * Test method for {@link jackSolver.JHFSolver#HFSolver(Atom[] atoms, int numberofelectrons)}.
	 */
	@Test
	public void test2HFSolver() {
		System.out.println("=================================================================");
		System.out.println("calculation RHF for H-H ");
		Atom atom1 = new Atom(0,0,1.463200,1);
		Atom atom2 = new Atom(0,0,0,1);
		//Atom atom3 = new Atom(0,2,2,1);
		//Atom atom4 = new Atom(0,2,3,1);
		Atom[] atoms2 = {atom2,atom1};
		
		
		int numberofelectrons =2;
		System.out.println("Calculating H-H energy");
		
		double energy = JackRHFSolver.HFSolver(atoms2,numberofelectrons);
		System.out.println(" energy: "+energy );
		assertTrue(Math.abs(energy+1.1140149480073371)<=0.0001);
		
	}
	
	@Test
	public void test3HFSolver(){
		System.out.println("calculation RHF for H-He+ ");
		Atom atom1 = new Atom(0,0,1.463200,1);
		Atom atom2 = new Atom(0,0,0,2);
		Atom atom3 = new Atom(0,5,2,1);
		Atom atom4 = new Atom(0,5,3,1);
		Atom[] atoms = {atom2,atom1,atom3,atom4};
		
		int numberofelectrons =2;
		System.out.println("Calculating H-He energy");
		
		double energy = JackRHFSolver.HFSolver(atoms,numberofelectrons);
		System.out.println(" energy: "+energy );
		//assertTrue(Math.abs(energy-1.1140149480073371)<=0.0001);
	}
	
	
	
	/**
	 * Test method for {@link jackSolver.JHFSolver#HFSolver(Atom[] atoms, int numberofelectrons)}.
	 */
	@Test
	public void testH2GasTestHFSolver() {
		System.out.println("=================================================================");
		System.out.println("========H2 Gas Test==============================================");
		int x =2,y=2,z=2;
		
		Atom[] atoms=generateH2Gas(x,y,z);
 
		
		
		int numberofelectrons =2*x*y*z;
		System.out.println("Calculating gas energy");
		
		double energy = JackRHFSolver.HFSolver(atoms,numberofelectrons);
		System.out.println(" energy: "+energy );
		
	}
	
	
	
	private Atom[] generateH2Gas(int x, int y, int z) {
		System.out.print(" number of atoms: "+(2*x*y*z));
		Atom[] atoms = new Atom[(2*x*y*z)+1];
		int currentatom=0;
		double spacing = 20;
		for(int i=0;i<x;i++){
			System.out.println(" x loop");
			for(int j=0;j<y;j++){
				System.out.println(" y loop");
				for(int k=0;k<z;k++){
					System.out.println(" z loop z:"+z+" k"+k);
					System.out.println(" adding H2");
					atoms[currentatom]=new Atom(i*x*spacing,j*y*spacing,k*z*spacing,1);
					currentatom++;
					System.out.println(" adding H");
					atoms[currentatom]=new Atom(i*x*spacing+1.46,j*y*spacing,k*z*spacing,1);
					currentatom++;
				}
			}
		}
		atoms[currentatom]=new Atom(-10,-10,-10,1);
		return atoms;
	}



	@Test
	public void testCOLTcannonicalOrthogonalizationSolver() {
		System.out.println("===================================================");
		System.out.println("========testCOLTcannonicalOrthogonalization============");
		
		double[][] Smat = {{1.00,0.4507704116},{0.4507704116,1.0}};
		DoubleMatrix2D SMatrix = new DenseDoubleMatrix2D(Smat);
		
		System.out.println("===SMatrix===");
		System.out.println(SMatrix);
		
		DoubleMatrix2D XMatrix = JackRHFSolver.cannonicalOrthogonalization(SMatrix);
		System.out.println("===XMatrix old method===");
		System.out.println(XMatrix);
		
		
		System.out.println("===XMatrix new method===");
		//aka diagonalize
		XMatrix = cannonicalOrthogonalization(SMatrix);
		
		System.out.println(XMatrix);
		
		
		
	
	}
	
	public static DoubleMatrix2D cannonicalOrthogonalization(DoubleMatrix2D SMatrix) {
		//based on the book method for 2x2 matrices
		
		//X=Us^-1/2
		
		SingularValueDecomposition svd = new SingularValueDecomposition(SMatrix);
		//System.out.println(svd);
		
		DoubleMatrix2D SdiagMatrix = svd.getS();//diagonalize2x2MatrixgetEigenValues(SMatrix);
		DoubleMatrix2D UMatrix = svd.getU();//diagonalize2x2MatrixgetEigenVectors(SMatrix);
		//System.out.println("===SdiagMatrix===");
		//System.out.println(SdiagMatrix);
		
		Algebra alg = new Algebra();
		
				
		JackRHFSolver.inversesqrt(SdiagMatrix);
		DoubleMatrix2D XMatrix = alg.mult(UMatrix, SdiagMatrix);
		
		return XMatrix;
	}
	
	@Test
	public void testadjoint() {
		System.out.println("===================================================");
		System.out.println("========testadjoint================================");
		
		double[][] Smat = {{1,2},{3,4}};
		System.out.println("===SMatrix===");
		DoubleMatrix2D SMatrix = new DenseDoubleMatrix2D(Smat);
		System.out.println(SMatrix);
		System.out.println("===SadjointMatrix===");
		DoubleMatrix2D SadjointMatrix = JackRHFSolver.adjoint(SMatrix);
		System.out.print(SadjointMatrix);
		System.out.println("===================================================");
		double[][] mat2 = {{1,2,3},{4,5,6},{7,8,9}};
		System.out.println("===SMatrix===");
		SMatrix = new DenseDoubleMatrix2D(mat2);
		System.out.println(SMatrix);
		System.out.println("===SadjointMatrix===");
		SadjointMatrix = JackRHFSolver.adjoint(SMatrix);
		System.out.print(SadjointMatrix);
		
		
		
	}
	
	@Test
	public void testCOLTv2cannonicalOrthogonalizationSolver() {
		System.out.println("===================================================");
		System.out.println("========testCOLTcannonicalOrthogonalization============");
		
		double[][] Smat = {{1.00,0.4507704116},{0.4507704116,1.0}};
		DoubleMatrix2D SMatrix = new DenseDoubleMatrix2D(Smat);
		
		System.out.println("===SMatrix===");
		System.out.println(SMatrix);
		System.out.println("===XMatrix new method===");
		DoubleMatrix2D XMatrix = cannonicalOrthogonalizationEVD(SMatrix);
		
		System.out.println(XMatrix);
		
		
		System.out.println("===XMatrix old method===");
		//aka diagonalize
		XMatrix = cannonicalOrthogonalization(SMatrix);
		
		System.out.println(XMatrix);
		
		
		
	
	}
	
	public static DoubleMatrix2D cannonicalOrthogonalizationEVD(DoubleMatrix2D SMatrix) {
		//based on the book method for 2x2 matrices
		
		//X=Us^-1/2
		
		
		DoubleMatrix2D SadjointMatrix = JackRHFSolver.adjoint(SMatrix);
		EigenvalueDecomposition evd = new EigenvalueDecomposition(SadjointMatrix);
		System.out.println(evd);
		
		
		DoubleMatrix2D SdiagMatrix = evd.getD();//diagonalize2x2MatrixgetEigenValues(SMatrix);
		DoubleMatrix2D UMatrix = evd.getV();//diagonalize2x2MatrixgetEigenVectors(SMatrix);
		//System.out.println("===SdiagMatrix===");
		//System.out.println(SdiagMatrix);
		
		Algebra alg = new Algebra();
		
				
		JackRHFSolver.inversesqrt(SdiagMatrix);
		DoubleMatrix2D XMatrix = alg.mult(UMatrix, SdiagMatrix);
		
		return XMatrix;
	}

	

}
