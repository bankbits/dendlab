/*
 * #%L
 * hIPNAT plugins for Fiji distribution of ImageJ
 * %%
 * Copyright (C) 2017 Tiago Ferreira
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package dend.skel;

import java.awt.Checkbox;



import java.awt.Choice;
import java.awt.Font;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.ColorModel;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Vector;

import javax.swing.JMenuItem;
import javax.swing.JPopupMenu;

import dend.AnalyzeSkeleton2_;
import dend.ColorMaps;
import dend.Dijkstra;
import dend.Edge;
import dend.Graph;
import dend.IPNAT;
import dend.Point;
import dend.SkeletonResult;
import dend.Utils;
import dend.Vertex;
import dend.processing.Binary;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.plugin.ZProjector;
import ij.process.ImageProcessor;
import sc.fiji.analyzeSkeleton.AnalyzeSkeleton_;
//import sc.fiji.analyzeSkeleton.AnalyzeSkeleton_;
//import sc.fiji.analyzeSkeleton.AnalyzeSkeleton_;
//import sc.fiji.analyzeSkeleton.AnalyzeSkeleton_;
//import sc.fiji.analyzeSkeleton.Edge;
//import sc.fiji.analyzeSkeleton.Graph;
//import sc.fiji.analyzeSkeleton.Point;
//import sc.fiji.analyzeSkeleton.Vertex;
//import sc.fiji.skeletonize3D.Skeletonize3D_;
//import sholl.gui.EnhancedGenericDialog;
import sc.fiji.skeletonize3D.Skeletonize3D_;
import sholl.gui.EnhancedGenericDialog;



/**
 * This class implements the ImageJ {@code Strahler Analysis} plugin. For more
 * information, visit the hIPNAT repository
 * {@literal https://github.com/tferr/hIPNAT} and the the plugin's documentation
 * page: {@literal http://imagej.net/Strahler_Analysis}
 *
 *
 * @author Tiago Ferreira
 */
public class Strahler implements PlugIn, DialogListener {

	protected static final String URL = "http://imagej.net/Strahler_Analysis";

	/* Default value for max. number of pruning cycles */
	int maxOrder = 30;

	/* Default option for loop detection */
	private int pruneChoice = AnalyzeSkeleton_.SHORTEST_BRANCH;
	
	private int treeChoice = 1;

	/* Default option for 'root-protection' ROI */
	private boolean protectRoot = true;
	
	/* Default option for showing the generated paths image */
	private boolean showPaths = false;

	/* Default option for 'iteration-stack' output */
	private boolean outIS = false;

	private /* Default option for verbose mode */
	boolean verbose = true;

	/* Default option for tabular option */
	private boolean tabular = false;

	/* Remove isolated pixels from thinned images? */
	private boolean erodeIsolatedPixels = true;

	/* Title of main results window */
	private static final String STRAHLER_TABLE = "Strahler_Table";

	/* Title of detailed results window */
	private static final String VERBOSE_TABLE = "Strahler_Iteration_Log";
	
	/* Title of detailed euclidean branch distances to root table */
	private static final String BRANCH_DIST = "Euclidean Branch Distances";
	
	/* Title of accurate branch distances to root table */
	private static final String BRANCH_PATHS = "Accurate Branch Distances"; 
	
	/* Title of comparison results table */
	private static final String COMPARE_TABLE = "Compare Exact and Estimated Distances";
	
	ImageStack shortpathImage = null;

	public static final String[] treeTypes = { "Root endpoint", "Root junction" };
	
	/*
	 * Grayscale image for intensity-based pruning of skel. loops. While it is
	 * unlikely that the iterative pruning of terminal branches will cause new
	 * loops on pre-existing skeletons, offering the option to resolve loops
	 * with intensity based methods remains useful specially when analyzing
	 * non-thinned grayscale images.
	 */
	ImagePlus grayscaleImp = null;
	int grayscaleImpChoice;

	ImagePlus srcImp; // Image to be analyzed (we'll be working on a copy)
	boolean validRootRoi; // Flag assessing validity of 'root-protective' ROI
	String title; // Title of active image
	Roi rootRoi; // Reference to the "root-protecting" ROI
	ImageProcessor ip;

	/**
	 * Calls {@link fiji.Debug#run(String, String) fiji.Debug.run()} so that the
	 * plugin can be debugged from an IDE
	 *
	 * @param args
	 *            the arguments as specified in {@code plugins.config}
	 */
	public static void main(final String[] args) {

		// Debug.run("Strahler Analysis...", null);
		new ImageJ(); // start ImageJ	
		
		ImagePlus imp = IJ.openImage(Strahler.class.getResource("/skeleton2.tif").getFile());
		//final ImagePlus imp = new LSystemsTree().createTreeStack("DebugImp");
		imp.setRoi((int)(470), (int)(530), 25, 25);
		
		// throw an exception if the user did not define an ROI around the root
		if(imp.getRoi() == null) {
			throw new RuntimeException("No ROI found");
		}
		imp.show();
		IJ.runPlugIn(imp, "dend.skel.Strahler", null);
		WindowManager.addWindow(imp.getWindow());
	}

	/**
	 * This method is called when the plugin is loaded.
	 *
	 * @param arg
	 *            the arguments as specified in {@code plugins.config}
	 *
	 */
	@Override
	public void run(final String arg) {

		// Retrieve analysis image and its ROI
		srcImp = WindowManager.getCurrentImage();
		if (!validRequirements(srcImp))
			return;

		title = srcImp.getTitle();
		rootRoi = srcImp.getRoi();
		validRootRoi = (rootRoi != null && rootRoi.getType() == Roi.RECTANGLE);

		// TODO: 3D Roots are special. We need to:
		// 1) Check if ROI is associated with all slices or just one
		// 2) Ignore counts above/below the ROI, as needed
		// 3) Extend ip operations to stack
		if (srcImp.getNSlices() > 1) {
			final String warning = "3D images are currently supported with the following limitations:\n"
					+ "    - 'Root-protecting' ROIs are not yet supported\n"
					+ "    - Lengths are estimated from Z-projections\n \n"
					+ "These issues will be addressed in future releases.";
			if (IJ.macroRunning())
				IJ.log(warning);
			else
				IJ.showMessage("Warning", warning);
			validRootRoi = false;
		}

		// Retrieve grayscale image for intensity-based pruning of skel. loops
		if (!getSettings())
			return;
		
		// Work on a skeletonized copy since we'll be modifying the image
		if (rootRoi != null)
			srcImp.killRoi();
		final ImagePlus imp = srcImp.duplicate();
		if (rootRoi != null)
			srcImp.setRoi(rootRoi);
		ip = imp.getProcessor();
		skeletonizeWithoutHermits(imp);

		// Initialize ResultsTable: main and detailed info
		final ResultsTable rt = Utils.getTable(STRAHLER_TABLE);
		final ResultsTable logrt = Utils.getTable(VERBOSE_TABLE);
		ResultsTable tt = Utils.getTable(COMPARE_TABLE);
		//ResultsTable bt = Utils.getTable(BRANCH_DIST);
		//ResultsTable dt = Utils.getTable(BRANCH_PATHS);
		//ResultsTable ct = Utils.getTable(COMPARE_TABLE);

		// Analyze root
		ImagePlus rootImp;
		ImageProcessor rootIp = null;
		dend.SkeletonResult rootResult = null;
		ArrayList<Point> rootEndpointsList = null;
		int nRootEndpoints = 0, nRootJunctions = 0;
		Vertex rootVertex = null;
		ArrayList<Point> juncts = new ArrayList<Point>();
		List<Branch> branchesList = new ArrayList<Branch>();

		AnalyzeSkeleton2_ root = null;
		
		
		// The code for finding distances to the root hinges on
		// having a valid root ROI
		// if none is found, the distances calculated will be inaccurate
		if (validRootRoi && verbose) {

			// Duplicate entire canvas. Ignore tree(s) outside ROI
			rootImp = imp.duplicate();
			rootIp = rootImp.getProcessor();
			rootIp.setValue(0.0);
			rootIp.fillOutside(rootRoi);

			// Get root properties
			root = new AnalyzeSkeleton2_();
			root.setup("", rootImp);
			
		
			rootResult = root.run(pruneChoice, false, false, grayscaleImp, true, false);
			rootImp.flush();

			// We assume ROI contains only end-point branches, slab voxels and
			// no junction points. We'll thus remove end-points at ROI
			// boundaries
		
			// Our array juncts will hold all Points inside the root ROI
			juncts.addAll(rootResult.getListOfJunctionVoxels());
			juncts.addAll(rootResult.getListOfEndPoints());
			
			nRootJunctions = sum(rootResult.getJunctions());
			
			rootEndpointsList = rootResult.getListOfEndPoints();
			final ListIterator<Point> it = rootEndpointsList.listIterator();
			final Rectangle r = rootRoi.getBounds();
			while (it.hasNext()) {
				final Point p = it.next();
				if (p.x == r.x || p.y == r.y || p.x == (int) (r.x + r.getWidth() - 1)
						|| p.y == (int) (r.y + r.getHeight() - 1))
					it.remove();
			}
			rootResult.setListOfEndPoints(rootEndpointsList);
			nRootEndpoints = rootEndpointsList.size();
		}
		

		

		// Initialize display images. Use Z-projections to populate
		// iteration stack when dealing with 3D skeletons
		final int nSlices = imp.getNSlices();
		ZProjector zp = null;

		ImageStack imgStack = srcImp.getImageStack();
		
		
		final ImageStack iterationStack = new ImageStack(imp.getWidth(), imp.getHeight());
		if (nSlices > 1) {
			zp = new ZProjector(imp);
			zp.setMethod(ZProjector.MAX_METHOD);
			zp.setStartSlice(1);
			zp.setStopSlice(nSlices);
		}

		// Initialize AnalyzeSkeleton_
		final AnalyzeSkeleton2_ as = new AnalyzeSkeleton2_();
		as.setup("", imp);

		// Perform the iterative pruning
		int order = 1, nEndpoints = 0, nJunctions = 0, nJunctions2 = 0;
		ArrayList<Point> endpointsList = null, junctionsList = null;
		
		int maxBranches = 0;
		
		// original skeleton result BEFORE pruning
		// perform dijkstra on this skeleton before branches are cut off
		// and original graph connects branches to the root
		dend.SkeletonResult ogSr = as.run(0, false, false, grayscaleImp, true, false);
		Graph graphOG = null;
		
		// Array to hold list of Vertices that are end points of branches
		ArrayList<Point> ends = ogSr.getListOfEndPoints();
		ArrayList<Vertex> endVertices = new ArrayList<Vertex>();
		
		
		// Find the grpah with the most number of branches which should
		// represent the bulk of the neuron skeleton structure
		for(Graph g : ogSr.getGraph()) {
			if(g.getVertices().size() > maxBranches) {
				maxBranches = g.getVertices().size();
				graphOG = g;
			}
		}
		
		// this loop adds vertices to our list that lie on the endpoints of the graph
		// aka, the outermost branches that we want to track
		if(graphOG != null) {
			for(Vertex v : graphOG.getVertices()) {
				for(Point p : v.getPoints()) {
					if(ends.contains(p)) {
						endVertices.add(v);
						break;
					}
				}
			}
		}
		
		// Finds the vertex on the graph associated to the correct root points
		// taken from junction and endpoint voxels in the root ROI
		Vertex[] verts = new Vertex[graphOG.getVertices().size()];
		verts = graphOG.getVertices().toArray(verts);
		
		rootVertex = null;
		ArrayList<Edge> rootBranches = new ArrayList<Edge>();
		for(Point j : juncts) {
			rootVertex = as.findPointVertex(verts, j);
			if(rootVertex != null) {
				rootBranches = rootVertex.getBranches();
				break;
			}
		}
		Edge[] rootEdges = new Edge[rootBranches.size()];
		rootEdges = rootBranches.toArray(rootEdges);
		
		
		
		Point rootPoint = null;
		
		// IF a root vertex / junction does not exist, we have to create one ourselves
		// REEEEEE
		//
		//
		//
		//
		if(rootVertex == null) {
			
			as.run(0, false, false, grayscaleImp, true, false);
			ArrayList<Point> junctions = as.listOfJunctionVoxels;
			
			Rectangle rect = rootRoi.getBounds();
			int x = rect.x + rect.width/2;
			int y = rect.y + rect.height/2;
			
			Point center = new Point(x, y, 0);
			
			double minDist = Double.MAX_VALUE;
			Point minPoint = null;
			Vertex minVertex = null;
			
			for(Point p : junctions) {
				double dist = calculateDistance(p, center, srcImp);
				if(dist < minDist) {
					minDist = dist;
					minPoint = p;
					
					
					minVertex = as.findPointVertex(verts, p);
					if(rootVertex != null) {
						rootBranches = rootVertex.getBranches();
						break;
					}
				}
			}
			
			rootVertex = minVertex;
			
			rootPoint = minPoint;
			
			for(Edge e : rootBranches) {
				rootBranches.addAll(e.getV1().getBranches());
				rootBranches.addAll(e.getV2().getBranches());
			}
			
			rootEdges = new Edge[rootBranches.size()];
			
			rootEdges = rootBranches.toArray(rootEdges);

		} else {
			rootPoint = rootVertex.getPoints().get(0);
		}
		
		
		
		// Copy input image 
		ColorModel cm = ColorModel.getRGBdefault();
		imgStack = new ImageStack(srcImp.getWidth(), srcImp.getHeight(), cm);
		for(int i=1; i<=srcImp.getStack().getSize(); i++)
			imgStack.addSlice(srcImp.getImageStack().getSliceLabel(i), srcImp.getImageStack().getProcessor(i).duplicate());

		// Finally run our path finding algorithm on the graph to find
		// paths and distances between every vertex in the graph and
		// the root vertex
		Map<Vertex, Double> branchDist = testDijkstra(imgStack, graphOG, rootVertex, showPaths);
		
		/**
		 * Create list of branches with properties and relevant data
		 */
		ArrayList<Edge> edges = graphOG.getEdges();
		for (Edge e : edges) {
			Branch branch = new Branch(e);
			
			// Find the path distance for this edge
			for(Map.Entry<Vertex, Double> entry : branchDist.entrySet()) {
				if(entry.getKey() == e.getV1()) {
					branch.setPathDist(entry.getValue());
					break;
				} else if (entry.getKey() == e.getV2())  {
					branch.setPathDist(entry.getValue());
					branch.setVertex(e.getV2());
				} else {
					// sets path disatnce to -1 if it couldn't be matched to an entry
					branch.setPathDist(-1);
				}
			}
			
			// Find the Euclidean distance for this edge
			Point p = branch.getVertex().getPoints().get(0);
			double dist = calculateDistance(p, rootPoint, srcImp);
			branch.setDist(dist);
			
			// add new Branch to our list of branches
			branchesList.add(branch);
		}
		
		
		
		// This do-while loop performs the function of pruning end point branches
		// at each order, re-skeletonizing the image and classifying each branch
		// by Strahler order
		do {

			IJ.showStatus("Retrieving measurements for order " + order + "...");
			IJ.showProgress(order, getMaxOrder());

			// (Re)skeletonize image
			if (order > 1)
				skeletonizeWithoutHermits(imp);

			// Get properties of loop-resolved tree(s)
			
			
			SkeletonResult sr = as.run(pruneChoice, false, false, grayscaleImp, true, false);
			
			ArrayList<Point> totalPoints = new ArrayList<Point>();
			totalPoints.addAll(sr.getListOfEndPoints());
			totalPoints.addAll(sr.getListOfSlabVoxels());
			totalPoints.addAll(sr.getListOfJunctionVoxels());
			totalPoints.addAll(sr.getListOfStartingSlabVoxels());
			
			for(Branch b : branchesList) {
				ArrayList<Point> branchPoints = new ArrayList<Point>();
				//branchPoints.addAll(b.getEdge().getV1().getPoints());
				//branchPoints.addAll(b.getEdge().getV2().getPoints());
				branchPoints.addAll(b.getEdge().getSlabs());
				
				for(Point p : branchPoints) {
					if(totalPoints.contains(p))
						b.setOrder(6 - order);
				}
			}
			
			nEndpoints = sum(sr.getEndPoints());
			nJunctions = sum(sr.getJunctions());

			if (order == 1) {
				// Remember initial properties
				endpointsList = sr.getListOfEndPoints();
				junctionsList = sr.getListOfJunctionVoxels();

				// Do not include root in 1st order calculations
				nEndpoints -= nRootEndpoints;
				nJunctions -= nRootJunctions;
			}

			// Is it worth proceeding?
			if (nEndpoints == 0 || nJunctions2 == nJunctions) {
				//errorMsg = "Error! Iteration " + order + " aborted: ";
				// errorMsg += (nEndpoints == 0) ? "No end-poins found" : "Unsolved loop(s) detected";
				break;
			}

			// Add current tree(s) to debug animation
			ImageProcessor ipd;
			if (nSlices > 1 && zp != null) {
				zp.doProjection();
				ipd = zp.getProjection().getProcessor();
			} else {
				ipd = ip.duplicate();
			}
			
			iterationStack.addSlice("Order " + IJ.pad(order, 2), ipd);

			// Report properties of pruned structures
			if (verbose) {
				logrt.incrementCounter();
				logrt.addValue("Image", title);
				logrt.addValue("Structure", "Skel. at iteration " + Integer.toString(order));
				logrt.addValue("# Trees", sr.getNumOfTrees());
				logrt.addValue("# Branches", sum(sr.getBranches()));
				logrt.addValue("# End-points", nEndpoints);
				logrt.addValue("# Junctions", nJunctions);
				logrt.addValue("# Triple points", sum(sr.getTriples()));
				logrt.addValue("# Quadruple points", sum(sr.getQuadruples()));
				logrt.addValue("Average branch length", average(sr.getAverageBranchLength()));
			}

			// Remember main results
			nJunctions2 = nJunctions;
			
			// Eliminate end-points
			as.run(pruneChoice, true, false, grayscaleImp, true, false, rootBranches);
			
			// Eliminate end-points
			//as.run(pruneChoice, true, false, grayscaleImp, true, false);

		} while (order++ <= getMaxOrder() && nJunctions > 0);
		
		// Set counter to the de facto order
		order -= 1;
		
		// set branches where the order wasn't found equal to the last order (eg 5)
		
		/*for(Branch b : branchesList) {
			if(b.getOrder() == 1)
				b.setOrder(2);
			if(b.getPathDist() == 0.0 || b.getDist() == 0.0)
				b.setOrder(1);
			if(b.getOrder() == 0 && order <= 5)
				b.setOrder(5);
			else if(b.getOrder() == 0 && order >= 6)
				b.setOrder(6);
			
			if(b.getOrder() == -1)
				b.setOrder(7);
			if(b.getOrder() == -2)
				b.setOrder(8);
		}*/
		
		// sort the branches by order, with the root order at the top
		Collections.sort(branchesList);

		// Append root properties to log table
		if (validRootRoi && verbose) {

			// Check if ROI contains unexpected structures
			
			logrt.incrementCounter();
			logrt.addValue("Image", title);
			logrt.addValue("Structure", "Root");
			logrt.addValue("# Trees", rootResult == null ? 0 : rootResult.getNumOfTrees());
			logrt.addValue("# Branches", rootResult == null ? 0 : sum(rootResult.getBranches()));
			logrt.addValue("# End-points", nRootEndpoints);
			logrt.addValue("# Junctions", nRootJunctions);
			logrt.addValue("# Triple points", rootResult == null ? 0 : sum(rootResult.getTriples()));
			logrt.addValue("# Quadruple points", rootResult == null ? 0 : sum(rootResult.getQuadruples()));
			logrt.addValue("Average branch length",
					rootResult == null ? Double.NaN : average(rootResult.getAverageBranchLength()));

		}
		
		// Safety check
		if (iterationStack.getSize() < 1) {
			error("Enable \"detailed\" mode and check " + VERBOSE_TABLE + " for details.");
			return;
		}
		
		
		
		// Create iteration stack
		final Calibration cal = srcImp.getCalibration();
		final ImagePlus imp2 = new ImagePlus("StrahlerIteration_" + title, iterationStack);
		
		System.out.println(imp2.isThreshold());
		
		imp2.setCalibration(cal);
		if (outIS) {
			if (validRootRoi) {
				iterationStack.addSlice("Root", rootIp);
				paintPoints(iterationStack, rootEndpointsList, 255, "Root end-points");
				imp2.setRoi(rootRoi);
			}
			paintPoints(iterationStack, endpointsList, 255, "End-points");
			paintPoints(iterationStack, junctionsList, 255, "Junction-points");
		}

		// Generate Strahler mask
		zp = new ZProjector(imp2);
		zp.setMethod(ZProjector.SUM_METHOD);
		zp.setStartSlice(1);
		zp.setStopSlice(order);
		zp.doProjection();
		final ImageProcessor ip3 = zp.getProjection().getProcessor().convertToShortProcessor(false);
		clearPoints(ip3, junctionsList); // disconnect branches
		ip3.multiply(1 / 255.0); // map intensities to Strahler orders
		final ImagePlus imp3 = new ImagePlus("StrahlerMask_" + title, ip3);
		imp3.setCalibration(cal);
		
		
		Graph[] graphs = ogSr.getGraph();

		// Analyze segmented order
		// Segment branches by order
		ImagePlus maskImp = imp3.duplicate(); // Calibration is
													// retained
		IJ.setThreshold(maskImp, 1, 1);
		IJ.run(maskImp, "Convert to Mask", "");

		// Analyze segmented order
		AnalyzeSkeleton2_ maskAs = new AnalyzeSkeleton2_();
		maskAs.setup("", maskImp);
		SkeletonResult maskSr = maskAs.run(pruneChoice, false, false, grayscaleImp, true, false);
		maskImp.flush();
		
		// trees if the user requested it
		int nBranches = (erodeIsolatedPixels) ? sum(maskSr.getBranches()) : maskSr.getNumOfTrees();
		maxBranches = nBranches;
		
		// array to hold estimated distances of branches
		double[][] distances = new double[maxBranches][order];
		graphs = maskSr.getGraph();
		
		// Measure segmented orders
		double prevNbranches = Double.NaN;
		for (int i = 0; i < order; i++) {

			// Segment branches by order
			maskImp = imp3.duplicate(); // Calibration is retained
			IJ.setThreshold(maskImp, i, i);
			IJ.run(maskImp, "Convert to Mask", "");

			// Analyze segmented order
			maskAs = new AnalyzeSkeleton2_();
			maskAs.setup("", maskImp);
			maskSr = maskAs.run(pruneChoice, false, false, grayscaleImp, true, false);
			
			if(i == order) {
				Graph[] maskGraphs = maskSr.getGraph();
			}
			
			maskImp.flush();
			nBranches = (erodeIsolatedPixels) ? sum(maskSr.getBranches()) : maskSr.getNumOfTrees();

			// Log measurements
			rt.incrementCounter();
			rt.addValue("Image", title);
			rt.addValue("Strahler Order", order - i + 1);
			rt.addValue("# Branches", nBranches);
			rt.addValue("Ramification ratios", prevNbranches / nBranches);
			rt.addValue("Average branch length", average(maskSr.getAverageBranchLength()));
			rt.addValue("Unit", cal.getUnit());
			String noteMsg = "";
			if (i == 1) {
				noteMsg = (erodeIsolatedPixels) ? "Ignoring" : "Including";
				noteMsg += " single-point arbors...";
			}
			//rt.addValue("Notes", noteMsg);

			// Remember results for previous order
			prevNbranches = nBranches;
		}
		
		
		// Display estimated and exact distances in outputted results table
		//displayBranchDistances(order, maxBranches, distances, bt);
		//displayExactDistances(branchDist, endVertices, dt);
		//displayComparisonTable(graphOG, branchDist, ct, rootPoint);
			
		// Display outputs
		if (!tabular) {
			if (outIS)
				imp2.show();
			ip3.setMinAndMax(1, order + 1);
			ColorMaps.applyViridisColorMap(imp3,  1, false);
			if (validRootRoi)
				imp3.setRoi(rootRoi);
			imp3.show();
			addCalibrationBar(imp3, 0, "White");
		}
		if (verbose)
			logrt.show(VERBOSE_TABLE);
		rt.show(STRAHLER_TABLE);
		
		displayAllResults(branchesList, tt);

		IJ.showProgress(0, 0);
		IJ.showTime(imp, imp.getStartTime(), "Strahler Analysis concluded... ");
		imp.flush();
	
	}
	
	/* -----------------------------------------------------------------------*/
	/**
	 * Set pixel in 3D image.
	 * 
	 * @param image 3D image
	 * @param p point coordinates
	 * @param value pixel value
	 */
	private void setPixel(ImageStack image, Point p, byte value)
	{
		if(p.x >= 0 && p.x < image.getWidth() && p.y >= 0 && p.y < image.getHeight() && p.z >= 0 && p.z < image.getBitDepth())
			((byte[]) image.getPixels(p.z + 1))[p.x + p.y * image.getWidth()] = value;
	} // end setPixel 
	
	
	/**
	 * Test our implementation of Dijkstra's algorithm using our 'Dijkstra' object
	 * Given a graph, a source and destination vertex, this function finds the shortest
	 * distance and the path between the two vertices
	 * 
	 * 
	 * @param graph
	 * @param source
	 * @param dest
	 */
	public Map<Vertex, Double> testDijkstra(ImageStack img, Graph graph, Vertex source, boolean showPaths) {
		
		Dijkstra dij = new Dijkstra(graph);
		dij.execute(source);
		
		// value to color the graph
		byte SHORTEST_PATH = 96;
		
		// we store the branch distsances in a Map that maps each vertex
		// to its respective distance to the root
		Map<Vertex, Double> distances = dij.getDistances();

		if(showPaths) {
		// Colors in the path from each vertex to the source
			for(Vertex v : graph.getVertices()) {
				
				LinkedList<Vertex> path = dij.getPath(v);
				
				if(path != null) {
					for(Vertex z : path) {
						for(Point p : z.getPoints()) {
							setPixel(img, p, SHORTEST_PATH);
						}
						for(Edge e : z.getBranches()) {
							for(Point s : e.getSlabs()) {
								setPixel(img, s, SHORTEST_PATH);
							}
						}
					}
				}
			}
		}
		
		return distances;
	}
	
	/**
	 * Displays all results for Branches in the skeleton
	 * Info includes: Strahler order, branch length, path distance, distance
	 * 
	 * @param branches List of Branch objects
	 * @param tt	   Results table to write to
	 */
	public void displayAllResults(List<Branch> branches, ResultsTable tt) {
		
		int counter = 0;
		for(Branch b : branches) {
			counter++;
			tt.incrementCounter();
			tt.addValue("Branch #", counter);
			tt.addValue("Strahler Order", b.getOrder());
			tt.addValue("Branch Length", b.getLength());
			tt.addValue("Path Distance", b.getPathDist());
			tt.addValue("Distance", b.getDist());
		}
		tt.show(COMPARE_TABLE);
		
	}
	
	
	
	/**
	 * Display the Euclidean distance from each branch endpoint to the root
	 * of the skeleton in a ResultsTable called "Branch Distances"
	 * 
	 * @param order			the Strahler order of the branches
	 * @param maxBranches	the maximum number of branches in one Strahler order
	 * @param distances 	the distances of each branch to the root
	 * @param bt			the results table to populate
	 */
	public void displayBranchDistances(int order, int maxBranches, double[][] distances, ResultsTable bt) {

		// Set up the heading of the results table
		bt.incrementCounter();
		bt.addValue("Branch #", "Average Distance");
		
		// calculate and display the average distance in the results table
		for(int i = 0; i < order; i++) {
			double average = 0;
			int counter = 0;
			for(int j = 0; j < maxBranches; j++) {
				if(distances[j][i] > 0.0) {
					counter++;
					average += distances[j][i];
				}
			}
			average /= counter;
			bt.addValue("Strahler Order " + (order - i), average);
		}
		bt.incrementCounter();
		
		// Loop through and display branch lengths in the results table
		for(int i = 0; i < maxBranches; i++) {
			bt.incrementCounter();
			bt.addValue("Branch #", i);
			for(int j = 0; j < order; j++) {
				
					bt.addValue("Strahler Order " + (order - j), distances[i][j]);
			}
		}

		bt.show(BRANCH_DIST);
	}
	
	/**
	 * Displays the exact distances from each vertex in the graph to the root
	 * in a results table
	 * 
	 * @param distances	Map of vertex to distances
	 * @param endVerts	list of endpoint vertices
	 * @param dt		our exact distances ResultsTable
	 */
	public void displayExactDistances(Map<Vertex, Double> distances, ArrayList<Vertex> endVerts, ResultsTable dt) {
		dt.incrementCounter();
		
		int counter = 0;
		
		for(Map.Entry<Vertex, Double> entry : distances.entrySet()) {
			counter++;
			dt.incrementCounter();
			
			Vertex key = entry.getKey();
			Double val = entry.getValue();
			
			dt.addValue("Vertex #", counter);
			dt.addValue("Distance", val);
			
			if(endVerts.contains(key))
				dt.addValue("EndPoint Vertex", "Yes");
			else
				dt.addValue("EndPoint Vertex", "No");
		}
		dt.show(BRANCH_PATHS);
	}
	
	
	/**
	 * Outputs the comparison values for the exact and estimated branch distances
	 * 
	 * @param g
	 * @param branchDist
	 * @param ct
	 * @param rootPoint
	 */
	public void displayComparisonTable(Graph g, Map<Vertex, Double> branchDist, ResultsTable ct, Point rootPoint) {

		ct.incrementCounter();
		int count = 1;
		
		for(Vertex v : g.getVertices()) {
			if(v.getPoints().size() >= 1) {
				Point branchPoint = v.getPoints().get(0);
				double dist = calculateDistance(branchPoint, rootPoint, srcImp);
				if(branchDist.get(v) != null) {
					count++;
					ct.incrementCounter();
					ct.addValue("Branch #", count);
					ct.addValue("Exact Distance", branchDist.get(v));
					ct.addValue("Estimated Distance", dist);
				}
			}
		}
		ct.show(COMPARE_TABLE);
		
	}	
	// -----------------------------------------------------------------------
	/**
	 * Calculate Euclidean distance between two points in 3D.
	 * 
	 * @param point1 first point coordinates
	 * @param point2 second point coordinates
	 * @return distance (in the corresponding units)
	 */
	private double calculateDistance(Point point1, Point point2, ImagePlus image) 
	{
		return Math.sqrt(  Math.pow( (point1.x - point2.x) * image.getCalibration().pixelWidth, 2) 
				          + Math.pow( (point1.y - point2.y) * image.getCalibration().pixelHeight, 2)
				          + Math.pow( (point1.z - point2.z) * image.getCalibration().pixelDepth, 2));
	}

	/**
	 * Checks if image to be analyzed fulfills analysis requirements and warns
	 * the user if required dependencies are present (i.e,, if all the required
	 * update sites have been subscribed).
	 *
	 * @param imp
	 *            the image to be analyzed
	 * @return {@code true}, if assessment was successful. If {@code false} a
	 *         macro friendly {@link Utils#error} is displayed.
	 */
	boolean validRequirements(final ImagePlus imp) {
		boolean validImp = imp != null && imp.getBitDepth() == 8;
		final boolean validSetup = Utils.validSkelDependencies();
		if (!validImp) {
			final String msg = (imp == null) ? "An 8-bit image is required but none was found."
					: imp.getTitle() + " is not an 8-bit image.";
			if (IJ.macroRunning()) {
				Utils.error("Invalid image", msg, imp);
			} else {
				final GenericDialog gd = new GenericDialog("Invalid Image");
				gd.addMessage(msg);
				gd.enableYesNoCancel("OK", "Analyze Sample Image");
				gd.hideCancelButton();
				gd.showDialog();
				if (!gd.wasOKed() && !gd.wasCanceled()) {
					final LSystemsTree lst = new LSystemsTree();
					this.srcImp = lst.sampleTree();
					this.srcImp.setRoi(58, 130, 25, 35);
					srcImp.show();
					new ij.plugin.Zoom().run("in");
					validImp = true;
				}
			}
		}
		return validSetup && validImp;
	}

	/**
	 * Displays an error message that will not disrupt macro calls. This is
	 * useful for batch processing of images: Even if the analysis of a
	 * particular image fails, remaining images can still be analyzed by the
	 * same macro
	 *
	 * @param errorMsg
	 *            the error message
	 */
	private void error(final String errorMsg) {
		Utils.error("Strahler Analysis", errorMsg, srcImp);
	}

	/**
	 * Gets the analysis parameters from the user.
	 *
	 * @return {@code true} if the dialog input is valid and dialog was not
	 *         dismissed.
	 */
	private boolean getSettings() {

		final EnhancedGenericDialog gd = new EnhancedGenericDialog("Strahler Analysis :: " + IPNAT.getVersion());
		final Font headerFont = new Font("SansSerif", Font.BOLD, 12);
		gd.setSmartRecording(true);

		// Part 1. Main Options
		gd.setInsets(0, 0, 0);
		gd.addMessage("Tree Classification:", headerFont);
		gd.addCheckbox("Infer root end-points from rectangular ROI", protectRoot);
		gd.addCheckbox("Ignore single-point arbors (Isolated pixels)", erodeIsolatedPixels);
		
		

		// Part 2: Loop elimination
		gd.setInsets(25, 0, 0);
		gd.addMessage("Elimination of Skeleton Loops:", headerFont);
		gd.addChoice("Method:", AnalyzeSkeleton_.pruneCyclesModes, AnalyzeSkeleton_.pruneCyclesModes[pruneChoice]);

		// Part 3: Tree type
		gd.addChoice("Tree Type:", treeTypes, treeTypes[1]); // deafult sets root as junction
		
		// 8-bit grayscale is the only image type recognized by
		// AnalyzeSkeleton_,
		// so we'll provide the user with a pre-filtered list of valid choices
		final ArrayList<Integer> validIds = new ArrayList<>();
		final ArrayList<String> validTitles = new ArrayList<>();
		final int[] ids = WindowManager.getIDList();
		for (int i = 0; i < ids.length; ++i) {
			final ImagePlus imp = WindowManager.getImage(ids[i]);
			if (imp.getBitDepth() == 8) { // TODO: ignore composites?
				validIds.add(ids[i]);
				validTitles.add(imp.getTitle());
			}
		}
		gd.addChoice("8-bit grayscale image:", validTitles.toArray(new String[validTitles.size()]), title);

		// Part 4: Output
		gd.setInsets(25, 0, 0);
		gd.addMessage("Output Options:", headerFont);
		gd.addCheckbox("Display_iteration stack", outIS);
		gd.addCheckbox("Show detailed information", verbose);
		gd.addCheckbox("Tabular data only (no image output)", tabular);

		gd.addDialogListener(this);
		dialogItemChanged(gd, null); // update prompt

		// Add More>> dropdown menu
		final JPopupMenu popup = new JPopupMenu();
		JMenuItem mi;
		mi = new JMenuItem("Online documentation");
		mi.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(final ActionEvent e) {
				IJ.runPlugIn("ij.plugin.BrowserLauncher", URL);
			}
		});
		popup.add(mi);
		popup.addSeparator();
		mi = new JMenuItem("List hIPNAT commands...");
		mi.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(final ActionEvent e) {
				IJ.runPlugIn("ij.plugin.BrowserLauncher", IPNAT.DOC_URL + "#List_of_commands");
			}
		});
		popup.add(mi);
		mi = new JMenuItem("About hIPNAT plugins...");
		mi.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(final ActionEvent e) {
				IJ.runPlugIn("ipnat.Help", "");
			}
		});
		popup.add(mi);
		gd.assignPopupToHelpButton(popup);

		gd.showDialog();

		// Set grayscale image for intensity-based pruning of skel. loops.
		if (pruneChoice == AnalyzeSkeleton_.LOWEST_INTENSITY_VOXEL
				|| pruneChoice == AnalyzeSkeleton_.LOWEST_INTENSITY_BRANCH) {
			grayscaleImp = WindowManager.getImage(validIds.get(grayscaleImpChoice));
		} else {
			grayscaleImp = null;
		}

		return gd.wasOKed();

	}

	/* Retrieve dialog options using the DialogListener interface */
	@Override
	public boolean dialogItemChanged(final GenericDialog gd, final java.awt.AWTEvent e) {

		protectRoot = gd.getNextBoolean();
		erodeIsolatedPixels = gd.getNextBoolean();
		pruneChoice = gd.getNextChoiceIndex();
		treeChoice = gd.getNextChoiceIndex();
		grayscaleImpChoice = gd.getNextChoiceIndex();
		outIS = gd.getNextBoolean();
		verbose = gd.getNextBoolean();
		tabular = gd.getNextBoolean();

		// Enable/Disable key components of GenericDialog
		if (!IJ.macroRunning()) {
			final Choice cImgChoice = (Choice) gd.getChoices().elementAt(1);
			final Choice treeImgChoice = (Choice) gd.getChoices().elementAt(2);
			final Vector<?> checkboxes = gd.getCheckboxes();
			final Checkbox roiOption = (Checkbox) checkboxes.elementAt(0);
			final Checkbox stackOption = (Checkbox) checkboxes.elementAt(3);

			cImgChoice.setEnabled(pruneChoice == AnalyzeSkeleton_.LOWEST_INTENSITY_VOXEL
					|| pruneChoice == AnalyzeSkeleton_.LOWEST_INTENSITY_BRANCH
					|| pruneChoice == AnalyzeSkeleton_.SHORTEST_BRANCH
					|| pruneChoice == AnalyzeSkeleton_.NONE);
			treeImgChoice.setEnabled(treeChoice == 0 || treeChoice == 1);
			roiOption.setEnabled(validRootRoi);
			stackOption.setEnabled(!tabular);

		}

		return !gd.wasCanceled();

	}

	/**
	 * Returns the sum of the values in the input array, or zero if the array is
	 * empty or {@code null}.
	 *
	 * @param array
	 *            array of values to be summed
	 * @return the sum of elements in the array. Returns zero if array is
	 *         {@code null} or empty.
	 */
	int sum(final int[] array) {
		int sum = 0;
		if (array != null)
			for (final int i : array)
				sum += i;
		return sum;
	}

	/**
	 * Returns the sum of the values in the input array, or zero if the array is
	 * empty or {@code null}.
	 *
	 * @param array
	 *            array of values to be summed
	 * @return the sum of elements in the array. Returns zero if array is
	 *         {@code null} or empty.
	 */
	double sum(final double[] array) {
		double sum = 0; // TODO Use org.apache.commons.math3.stat.StatUtils?
		if (array != null && array.length > 0)
			for (final double i : array)
				sum += i;
		return sum;
	}

	/**
	 * Returns the arithmetic mean of the values in the input array, or
	 * {@code Double.NaN} if the array is empty or {@code null}.
	 *
	 * @param array
	 *            array of values to be averaged
	 * @return the arithmetic mean of the array. Returns {@code Double.NaN} if
	 *         array is {@code null} or empty.
	 */
	double average(final double[] array) {
		if (array != null && array.length > 0)
			return sum(array) / array.length; // TODO Use
												// org.apache.commons.math3.stat.StatUtils?
		return Double.NaN;
	}

	/* Paints point positions. */
	void paintPoints(final ImageStack stack, final ArrayList<Point> points, final int value, final String sliceLabel) {
		if (points != null) {
			final ImageProcessor ipp = stack.getProcessor(1).createProcessor(stack.getWidth(), stack.getHeight());
					//ip.createProcessor(stack.getWidth(), stack.getHeight());
			for (int j = 0; j < points.size(); j++) {
				final Point point = points.get(j);
				ipp.putPixel(point.x, point.y, value);
			}
			stack.addSlice(sliceLabel, ipp);
		}
	}

	/* Clears point positions */
	private void clearPoints(final ImageProcessor processor, final ArrayList<Point> points) {
		if (points != null) {
			for (int j = 0; j < points.size(); j++) {
				final Point point = points.get(j);
				processor.putPixel(point.x, point.y, 0);
			}
		}
	}

	/*
	 * Skeletonization method that erodes the thinned structure in order to
	 * eliminate isolated pixels. Thinning and pruning may give rise to single
	 * point arbors. These 'debris' trees have 1 end-point but no branches or
	 * junctions. If present they overestimate the total number of end-points
	 */
	private void skeletonizeWithoutHermits(final ImagePlus imp) {
		final Skeletonize3D_ thin = new Skeletonize3D_();
		thin.setup("", imp);
		thin.run(null);

		if (erodeIsolatedPixels)
			Binary.removeIsolatedPixels(imp);
	}

	/**
	 * Runs {@link ij.plugin.CalibrationBar} on the specified image using
	 * sensible settings.
	 *
	 * @param imp
	 *            the image to processed
	 * @param nLabels
	 *            the n. of labels in the calibration bar
	 * @param color
	 *            Labels' foreground color as per {@link ij.plugin.Colors}
	 **/
	private void addCalibrationBar(final ImagePlus imp, final int nLabels, final String color) {
		final ImageCanvas ic = imp.getCanvas();
		double zoom = (imp.getHeight() > 200) ? 1.0 : 0.8;
		final double mag = (ic != null) ? ic.getMagnification() : 1.0;
		if (zoom <= 1 && mag < 1)
			zoom = 1.0 / mag;
		IJ.run(imp, "Calibration Bar...",
				"fill=None label=" + color + " number=" + nLabels + " zoom=" + zoom + " overlay");
	}

	/**
	 * Returns the maximum Strahler order being considered by the plugin.
	 *
	 * @return The maximum number of pruning cycles of end-point branches that
	 *         the plugin should perform
	 * @see #setMaxOrder(int)
	 */
	public int getMaxOrder() {
		return maxOrder;
	}

	/**
	 * Sets the maximum Strahler order to be considered by the plugin.
	 *
	 * @param maxOrder
	 *            The maximum number of pruning cycles of end-point branches
	 *            that the plugin should perform
	 * @see #getMaxOrder()
	 */
	public void setMaxOrder(final int maxOrder) {
		this.maxOrder = maxOrder;
	}

}
