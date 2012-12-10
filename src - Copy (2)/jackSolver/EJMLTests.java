package jackSolver;

import static org.junit.Assert.*;

import org.ejml.alg.dense.decomposition.eig.SwitchingEigenDecomposition;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.EigenOps;
import org.ejml.ops.MatrixFeatures;
import org.ejml.simple.SimpleEVD;
import org.ejml.simple.SimpleMatrix;
import org.junit.Test;

import cern.colt.matrix.DoubleMatrix2D;
//import cern.colt.matrix.impl.DenseDoubleMatrix2D;
//import cern.colt.matrix.linalg.Algebra;
//import cern.colt.matrix.linalg.EigenvalueDecomposition;

public class EJMLTests {

	@Test
	public void test() {
		fail("Not yet implemented");
	}
	
	@Test
	public void testcannonicalOrthogonalization(){
		System.out.println("===================================================");
		System.out.println("========testcannonicalOrthogonalization============");
		
		double[][] Smat = {{1.00,0.4507704116},{0.4507704116,1.0}};
		DenseMatrix64F SMatrix = new DenseMatrix64F(Smat);
		MatrixFeatures.isSymmetric(SMatrix,0);
		
		//DoubleMatrix2D SMatrix = new DenseDoubleMatrix2D(Smat);
		
		System.out.println("===SMatrix===");
		System.out.println(SMatrix);
		
		DoubleMatrix2D XMatrix = cannonicalOrthogonalization(SMatrix);
		System.out.println("===XMatrix===");
		System.out.println(XMatrix);
		double value = XMatrix.hashCode();
		System.out.println(value);
		

		
	}
	
	public static DoubleMatrix2D cannonicalOrthogonalization(DenseMatrix64F SMatrix) {
		//based on the book method for 2x2 matrices
		
		//X=Us^-1/2
		
		SwitchingEigenDecomposition evd = new SwitchingEigenDecomposition(SMatrix.numCols, true, 0);
		evd.decompose(SMatrix);
		System.out.println("===EVD===");
		System.out.println(evd);

		DenseMatrix64F SdiagMatrix = EigenOps.createMatrixD(evd);//SMatrix.eig().; //diagonalize2x2MatrixgetEigenValues(SMatrix);
		
		System.out.println("===SdiagMatrix===");
		System.out.println(SdiagMatrix);
		
		DenseMatrix64F UMatrix = EigenOps.createMatrixV(evd);//diagonalize2x2MatrixgetEigenVectors(SMatrix);
		//Algebra alg = new Algebra();
		System.out.println("===UMatrix===");
		System.out.println(UMatrix);
		
		//inversesqrt(SdiagMatrix);
		//DoubleMatrix2D XMatrix = alg.mult(UMatrix, SdiagMatrix);
		
		return null;
	}
	
	void decompositionExample( DenseMatrix64F A ) {
        //EigenvalueDecomposition svd = EigenvalueDecomposition.;

        //if( !svd.decompose(A) )
        //    throw new RuntimeException("Decomposition failed");
        
        //DenseMatrix64F D = svd.getD();

        //DenseMatrix64F U = svd.getU(null,false);
        //DenseMatrix64F W = svd.getW(null);
        //DenseMatrix64F V = svd.getV(null,false);

    }

}
