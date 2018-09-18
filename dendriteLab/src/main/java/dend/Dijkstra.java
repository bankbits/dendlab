package dend;

import java.util.ArrayList;

import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;



public class Dijkstra {

	/* List of all vertices in our graph */
	private final List<Vertex> nodes;
	/* List of all edges */
    private final List<Edge> edges;
    
    /* Vertices marked as visited */
    private Set<Vertex> settledNodes;
    /* Vertices marked as unvisited */
    private Set<Vertex> unSettledNodes;
    /* List of previous vertices */
    private Map<Vertex, Vertex> predecessors;
    /* Map of distance from each Vertex to a specified "source" Vertex */
    private Map<Vertex, Double> distance;

    /**
     * Dijkstra constructor that simply requires an inputed Graph to execute on
     * @param graph the inputted graph (taken from a skeleton result object)
     */
    public Dijkstra(Graph graph) {
        // create a copy of the array so that we can operate on this array
        this.nodes = new ArrayList<Vertex>(graph.getVertices());
        this.edges = new ArrayList<Edge>(graph.getEdges());
    }
    
    /**
     * Runs the main algorithm to find shortest paths between all nodes and the source Vertex
     * @param source the Vertex to use as the source
     */
    public void execute(Vertex source) {
        settledNodes = new HashSet<Vertex>();
        unSettledNodes = new HashSet<Vertex>();
        distance = new HashMap<Vertex, Double>();
        predecessors = new HashMap<Vertex, Vertex>();
        
        distance.put(source, 0.0);
        unSettledNodes.add(source);
        
        while (unSettledNodes.size() > 0) {
            Vertex node = getMinimum(unSettledNodes);
            settledNodes.add(node);
            unSettledNodes.remove(node);
            findMinimalDistances(node);
        }
    }
    /**
     * Getter method to return the map of distances
     * @return distances
     */
    public Map<Vertex, Double> getDistances() {
    	return this.distance;
    }
    
    /**
     * Given a Vertex node, this method calculates the shortest distances between
     * the node Vertex and its neighbors
     * 
     * @param node
     */
    private void findMinimalDistances(Vertex node) {
        List<Vertex> adjacentNodes = getNeighbors(node);
        for (Vertex target : adjacentNodes) {
            if (getShortestDistance(target) > getShortestDistance(node)
                    + getDistance(node, target)) {
                distance.put(target, getShortestDistance(node)
                        + getDistance(node, target));
                predecessors.put(target, node);
                unSettledNodes.add(target);
            }
        }
    }
    
    /**
     * Get the distance between two vertices on opposite ends of an edge
     * 
     * @param node
     * @param target
     * @return edge length
     */
    private double getDistance(Vertex node, Vertex target) {
        for(Edge e : node.getBranches()) {
        	if(e.getV2().equals(target) || e.getV1().equals(target)) {
        		return e.getLength();
        	}
        }
        throw new RuntimeException("Should not happen");
    }

    /**
     * Gets a list of Vertices immediately connected to the node Vertex by
     * one edge
     * 
     * @param node
     * @return neighbors - List of neighbor vertices
     */
    private ArrayList<Vertex> getNeighbors(Vertex node) {
    	
    	ArrayList<Vertex> neighbors = new ArrayList<Vertex>();
    	ArrayList<Edge> edges = node.getBranches();
    	for(Edge e : edges) {
    		if(e.getV1().equals(node))
    			neighbors.add(e.getV2());
    		else if(e.getV2().equals(node))
    			neighbors.add(e.getV1());
    	}
    	return neighbors;
    }
    
    /**
     * Returns the Vertex with the shortest distance from a set of vertices
     * @param vertexes
     * @return
     */
    private Vertex getMinimum(Set<Vertex> vertices) {
        Vertex minimum = null;
        for (Vertex vertex : vertices) {
            if (minimum == null) {
                minimum = vertex;
            } else {
                if (getShortestDistance(vertex) < getShortestDistance(minimum)) {
                    minimum = vertex;
                }
            }
        }
        return minimum;
    }

    /**
     * Checks if a node is settled
     * @param vertex
     * @return true if settled, false otherwise
     */
    private boolean isSettled(Vertex vertex) {
        return settledNodes.contains(vertex);
    }
    
    /**
     * Returns the calculated shortest distance for a destination vertex
     * from the map of distances
     * 
     * @param destination
     * @return shortest distance
     */
    public double getShortestDistance(Vertex destination) {
        Double d = distance.get(destination);
        if (d == null) {
            return Integer.MAX_VALUE;
        } else {
            return d;
        }
    }
    
    /**
     * This method returns the path from the source node to the selected target
     * node and returns NULL if no path exists
     * 
     * @param target
     * @return path
     */
    public LinkedList<Vertex> getPath(Vertex target) {
        LinkedList<Vertex> path = new LinkedList<Vertex>();
        Vertex step = target;
        // check if a path exists
        if (predecessors.get(step) == null) {
            return null;
        }
        path.add(step);
        while (predecessors.get(step) != null) {
            step = predecessors.get(step);
            path.add(step);
        }
        // Put it into the correct order
        Collections.reverse(path);
        return path;
    }
    
}
