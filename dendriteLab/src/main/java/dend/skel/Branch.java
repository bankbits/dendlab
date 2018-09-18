package dend.skel;

import dend.Edge;
import dend.Vertex;

public class Branch implements Comparable<Branch> {
	
	
	private Edge edge = null;
	// Vertex on graph
	private Vertex vertex = null;
	// path distances
	private double pathDist = 0.0;
	// euclidean distance
	private double dist = 0.0;
	// Strahler order
	private int order = 0;
	// branch length
	private double length = 0;
	
	
	/***
	 * Branch point constructors
	 * Essentially a vertex object but also contains an "id" field
	 * that allows us to identify Vertex objects by name
	 */

	public Branch(Edge e) {
		this.edge = e;
		this.vertex = e.getV1();
		this.length = e.getLength();
	}
	
	public Branch(Edge e, double pathDist, double dist, int order) {
		this.edge = e;
		this.vertex = e.getV1();
		this.length = e.getLength();
		this.pathDist = pathDist;
		this.dist = dist;
		this.order = order;
	}

	public Vertex getVertex() {
		return vertex;
	}

	public double getPathDist() {
		return pathDist;
	}

	public double getDist() {
		return dist;
	}

	public int getOrder() {
		return order;
	}

	public double getLength() {
		return length;
	}
	
	public Edge getEdge() {
		return edge;
	}
	
	public void setVertex(Vertex vertex) {
		this.vertex = vertex;
	}

	public void setPathDist(double pathDist) {
		this.pathDist = pathDist;
	}

	public void setDist(double dist) {
		this.dist = dist;
	}

	public void setOrder(int order) {
		this.order = order;
	}
	
	public void setLength(int length) {
		this.length = length;
	}
	
	public void setEdge(Edge edge) {
		this.edge = edge;
	}
	
	public String toString() {
		return "Order: " + this.order + " Branch Length: " + this.length + 
				" Path Dist: " + this.pathDist + " Euc Dist: " + this.dist;
	}
	
	public int compareTo(Branch other) {
		Integer thisOrder = new Integer(this.order);
		Integer otherOrder = new Integer(other.order);
		return thisOrder.compareTo(otherOrder);
	}
	
	
	
	

}
