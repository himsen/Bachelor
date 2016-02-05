/*
 * Edwards curve class. Initializes the Edwards curve with curve parameter d and point [x,y,z].
 * Implements different additions and subtraction on Edwards curves.
 * Version 1.0.
 * Torben Hansen
 * Aarhus University, Institute of Mathematics. 
 */

import java.math.*;

public class Edward {
	
	/*
	 * TODO:
	 * - Check the correctness of all methods.
	 * - tripling og quadoubling
	 */
	
	/*
	 * Note: x.pow(2).mod(n) is faster than x.modPow(new BigInteger(2),n). See test class: 'WhatPowModMethodIsBest'
	 * If z=0 we are actually not allowed to add, but the probability that it happen is very low. 
	 * If dx_1x_2y_1y_2-1=0 or dx_2x_2y_1y_2+1=0 then additions i actually not allowed, but it happens with negligible probability.
	 */
	
	//Global fields
	BigInteger x,y,z; //Point P=[x,y,z].
	BigInteger d; //Curve parameter; x^2+y^2=1+dx^2y^2.
	BigInteger n; //Number we wish to find a factor in.
	BigInteger r1,r2; //Temporary registers. 
	BigInteger bigZero = BigInteger.ZERO; //A global ZERO BigInteger.
	int M,S,D,A; //Cost counting. (M multiplications, S squarings, D multiplications by d, A additions/subtractions)
	
	//Constructor
	public Edward(BigInteger x, BigInteger y, BigInteger z, BigInteger d, BigInteger n){
		this.x = x;
		this.y = y;
		this.z = z;
		this.d = d;
		this.n = n;
	}
	
	
	//NOTE hvis output er 0 svarer punktet ikke til et rationelt. 
	/*
	 * Addition.
	 * Input: Point Q=[u,v,w].
	 * Operation: Calculate the Edwards sum [x,y,z]+[u,v,w]. 
	 * Notes: Use r1 and r2 as temporary registers. The formula used is unified and even complete when d is not a square. 
	 * On platforms with S/M < 0.75 one should use another addition. See Bernstein andLanges paper: 'Faster addition and doubling on elliptic curves'.
	 * Update: this.x, this.y, this.z are updated with the sum. 
	 */
	public void add(BigInteger u, BigInteger v, BigInteger w){	
		//Addition.
		z = (z.multiply(w)).mod(n);
		r1 = (x.add(y)).mod(n);
		r2 = (u.add(v)).mod(n);
		x = (x.multiply(u)).mod(n);
		y = (y.multiply(v)).mod(n);
		r1 = (r1.multiply(r2)).mod(n);
		r1 = (r1.subtract(x)).mod(n);
		r1 = (r1.subtract(y)).mod(n);
		r1 = (r1.multiply(z)).mod(n);
		r2 = (x.multiply(y)).mod(n);
		r2 = (r2.multiply(d)).mod(n);
		y = (y.subtract(x)).mod(n);
		y = (y.multiply(z)).mod(n);
		z = (z.pow(2)).mod(n);
		x = (z.subtract(r2)).mod(n);
		z = (z.add(r2)).mod(n);
		y = (y.multiply(z)).mod(n);
		z = (z.multiply(x)).mod(n);
		x = (x.multiply(r1)).mod(n);
		//Count.
		M += 10; S += 1; D += 1; A += 7;		
		//Reset temporary registers.
		r1 = bigZero;
		r2 = bigZero;
	}

	/*
	 * Mixed addition.
	 * Input: Point [u,v,1]
	 * This is the same operation as in the method 'void add(BigInteger u, BigInteger v, BigInteger w)'. 
	 * Only difference is the knowledge that the input point is a rational point i.e. z-coordinate is 1. 
	 * In this case we save the multiplication A=z_1*z_2
	 * Update: this.x, this.y and thix.z is updated with [x,y,z]+[u,v,1].
	 */
	public void addMixed(BigInteger u, BigInteger v){
		//Mixed addition.
		r1 = (x.add(y)).mod(n);
		r2 = (u.add(v)).mod(n);
		x = (x.multiply(u)).mod(n);
		y = (y.multiply(v)).mod(n);
		r1 = (r1.multiply(r2)).mod(n);
		r1 = (r1.subtract(x)).mod(n);
		r1 = (r1.subtract(y)).mod(n);
		r1 = (r1.multiply(z)).mod(n);
		r2 = (x.multiply(y)).mod(n);
		r2 = (r2.multiply(d)).mod(n);
		y = (y.subtract(x)).mod(n);
		y = (y.multiply(z)).mod(n);
		z = (z.pow(2)).mod(n);
		x = (z.subtract(r2)).mod(n);
		z = (z.add(r2)).mod(n);
		y = (y.multiply(z)).mod(n);
		z = (z.multiply(x)).mod(n);
		x = (x.multiply(r1)).mod(n);
		//Count.
		M += 9; S += 1; D += 1; A += 7;		
		//Reset temporary registers.
		r1 = bigZero;
		r2 = bigZero;
	}
	
	/*
	 * Subtraction.
	 * Input:point [u,v,w]
	 * Operation: Compute the difference [x,y,z]-[u,v,w].
	 * Note: Any points P=[q1,q2,q3] has -P=[-q1,q2,q3].
	 * Update: this.x, this.y and this.z is updated with [x,y,z]-[u,v,w]=[x,y,z]+[-u,v,w].
	 */
	public void sub(BigInteger u,BigInteger v, BigInteger w){
		this.add((u.negate()).mod(n), v, w);
	}
	
	/*
	 * Mixed subtraction.
	 * Input: Point [u,v,1]
	 * Operation: Compute the difference [x,y,z]-[u,v,1].
	 * Note: Any points P=[q1,q2,q3] has -P=[-q1,q2,q3].
	 * Update: this.x, this.y and this.z is updated with [x,y,z]-[u,v,1]=[x,y,z]+[-u,v,1].
	 */
	public void subMixed(BigInteger u, BigInteger v){
		this.addMixed((u.negate()).mod(n),v);
	}

	public void addBothMixed(){
		
	}
	
	/*
	 * Doubling.
	 * Operations: Calculates a doubling; 2[x,y,z].
	 * Notes: This formula is unified and accept inverse and neutral element. r1 and r2 are temporary registers.
	 * Updates: this.x, this.y, this.z. are updated with 2[x,y,z].
	 */
	public void dup(){
		//Doubling.
		r1 = (x.add(y)).mod(n);
		x = (x.pow(2)).mod(n); 
		y = (y.pow(2)).mod(n);
		z = (z.pow(2)).mod(n);
		r1 = (r1.pow(2)).mod(n);
		z = (z.add(z)).mod(n);
		r2 = (x.add(y)).mod(n);
		y = (x.subtract(y)).mod(n);
		r1 = (r1.subtract(r2)).mod(n);
		z = (r2.subtract(z)).mod(n);
		x = (z.multiply(r1)).mod(n);
			r1 = z; //Special allocation used if z==0;
		z = (z.multiply(r2)).mod(n);
		y = (y.multiply(r2)).mod(n);
		//Count.
		M += 3; S += 4; A += 6;
		//Reset registers.
		r1 = bigZero;
		r2 = bigZero;
	}

	public void tripling(){
		
	}
}
