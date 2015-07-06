/*Taken from IvusSnakes ImageJ plugin http://ivussnakes.sourceforge.net/
	Licensed with the GPL 3.0 or above
*/

package org.doube.bonej.pqct.selectroi.liveWireEngine;

class Vector2d {

    private double x, y;
    public Vector2d(){
		x = 0.0;
		y = 0.0;
    }
    public Vector2d(double x, double y){
		this.x = x;
		this.y = y;
    }
    public Vector2d(Vector2d other){
		x = other.getX();
		y = other.getY();
    }
    public void setX(double x){
		this.x = x;
    }
    public void setY(double y){
		this.y = y;
    }
    public double getX(){
		return x;
    }
    public double getY(){
		return y;
    }

    //returns the modulus of the vector
    public double mod(){
		return Math.sqrt( x*x + y*y );
    }
    //returns this - other
    public Vector2d sub(Vector2d other){
		Vector2d ans = new Vector2d();
		ans.setX( getX() - other.getX());
		ans.setY( getY() - other.getY());
		return ans;
    }
    //returns a unit vector (versor) from this
    public Vector2d getUnit(){
		Vector2d ans= new Vector2d();
		if((Math.abs(x)<0.1)&&(Math.abs(y)<0.1)){
			//if the vector is null
			//we'll return the direction (1,0) for it
			//this is because of a problem that happens when 
			//we need to calculate the direction cost
			ans.setX(1.0);
			ans.setY(0.0);
			
		}
		else{
			ans.setX( getX()/mod());
			ans.setY( getY()/mod());	
		}
		return ans;
    }       
    public double dotProduct(Vector2d other){
	//find angle between vectors
		return (getX()*other.getX()+getY()*other.getY());
    }
    //returns the vector that is normal to this
    public Vector2d getNormal(){
		Vector2d ans = new Vector2d(this.getY(),-this.getX());
		return ans;
    }
    public static void main(String args[]){
		Vector2d a = new Vector2d(3,4);
		System.out.println("X " + a.getX());
		System.out.println("Y " + a.getY());
		System.out.println("Mod " + a.mod());
		System.out.println("Norm x " + a.getUnit().getX() + 
				   "Norm y " + a.getUnit().getY());
		Vector2d b = new Vector2d(3,7);
		System.out.println("dotProduct " + a.dotProduct(b));
    }
}