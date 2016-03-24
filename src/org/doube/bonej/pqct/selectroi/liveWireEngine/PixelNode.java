
package org.doube.bonej.pqct.selectroi.liveWireEngine;

/*
	Store pixel nodes in a Priority Queue
	so that Dijkstra can run on O(n log n)
	The interface Comparable is required so that the Java class PriorityQueue
	could be used
	 The code is licensed under GPL 3.0 or newer
*/
class PixelNode implements Comparable<PixelNode> {
    private int[] myIndex;
    private double myDistance;
    private int[] whereFrom;
    /**Constructor
    	@param index the index of the node
    	@param distance the cost to the node
    	@param whereFrom from which node we got to this node from
    */
    public PixelNode(int[] index, double distance, int[] whereFrom){
		myIndex = index;
		myDistance = distance;
		this.whereFrom = whereFrom;
    }
    public double getDistance(){
		return myDistance;
    }
    public int[] getIndex(){
		return myIndex;
    }
    public int[] getWhereFrom(){
		return whereFrom;
    }

    public int compareTo(PixelNode other){
		if( myDistance < other.getDistance()){
			return -1;
		}else{
			if( myDistance > other.getDistance()){
				return +1;
			}else{
				return 0;
			}
		}
	}

}
