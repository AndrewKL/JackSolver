package jackSolver;

public class Vec3D{
	public double x;
	public double y;
	public double z;
	
	public Vec3D(double xin,double yin,double zin){
		x=xin;
		y=yin;
		z=zin;
	}
	
	public String toString(){
		return " x: "+x+" y: "+y+" z: "+z;
	}
}
