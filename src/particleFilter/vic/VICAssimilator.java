/**
 * Copyright 2018 University of Pittsburgh
 * 
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
 * in compliance with the License. You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software distributed under the 
 * License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 * express or implied. See the License for the specific language governing permissions and
 * limitations under the License.
 */

package particleFilter.vic;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.time.Duration;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;

import org.apache.commons.io.FileUtils;

import particleFilter.Assimilator;
import particleFilter.Model;
import particleFilter.ModelResult;
import particleFilter.Particle;
import probDist.KernelDensity;
import probDist.Normal;
import probDist.multiVar.MultiVarKernelDensity;
import probDist.multiVar.tools.Sample;
import utilities.stat.ContSeries;
import vic.Forcing;
import vic.Soil;
import vic.routing.MuskingumNetwork;
import vic.routing.Simulation;
import vic.routing.State;

/**
 * Citation: Hernández, F. and Liang, X.: Hybridizing Bayesian and variational data assimilation
 * for high-resolution hydrologic forecasting, Hydrol. Earth Syst. Sci., 22, 5759-5779,
 * https://doi.org/10.5194/hess-22-5759-2018, 2018.
 * 
 * @author Felipe Hernández (developer)
 * @author Xu Liang (advisor)
 */
public class VICAssimilator implements Model
{

	// --------------------------------------------------------------------------------------------
	// Constants
	// --------------------------------------------------------------------------------------------
	
	public final static String	DATE_TIME_FORMAT_FOLDER		= "yyyy-MM-dd HH.mm";
	public final static String	MEAN_Q_FILE					= "/Streamflow.txt";
	
	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------

	private String						parametersFolder;
	private ArrayList<Soil>				soils;
	private ArrayList<String>			globalFileParams;
	private Duration					modelTimeStep;
	private ArrayList<Forcing>			cellForcings;
	private MuskingumNetwork			network;
	private Hashtable<String, Double>	areas;
	private ArrayList<String>			outputs;
	private Hashtable<String, Double>	directFractions;
	
	private String						vicExec;
	private long						simMaxTime;
	private long						forecastSimMaxTime;
	private boolean						removeFiles;
	
	private ArrayList<State>			initialStates;
	private Hashtable<String, Double>	currentStreamflow;
	
	private String						modelsFolder;
	
	private LocalDateTime				current;
	
	private HashSet<String>				toDelete;
	
	private KernelDensity				endStreamflow;
	private KernelDensity				endEvaporation;
	private KernelDensity				endSoilMoisture;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	public VICAssimilator(String parametersFolder, ArrayList<Soil> soils, 
				ArrayList<String> globalFileParams, Duration modelTimeStep, 
				ArrayList<Forcing> cellForcings, MuskingumNetwork network, 
				Hashtable<String, Double> areas, ArrayList<String> outputs,
				Hashtable<String, Double> directFractions, String vicExec, long simMaxTime,
				long forecastSimMaxTime, boolean removeFiles)
	{
		this.parametersFolder	= parametersFolder;
		this.soils				= soils;
		this.globalFileParams	= globalFileParams;
		this.modelTimeStep		= modelTimeStep;
		this.cellForcings		= cellForcings;
		this.network			= network;
		this.areas				= areas;
		this.outputs			= outputs;
		this.directFractions	= directFractions;
		this.vicExec			= vicExec;
		this.simMaxTime			= simMaxTime;
		this.forecastSimMaxTime	= forecastSimMaxTime;
		this.removeFiles		= removeFiles;
		this.current			= null;
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------

	public ArrayList<State> assimilate(LocalDateTime start, LocalDateTime end,
			LocalDateTime forecastEnd, ArrayList<State> initialStates, Hashtable<LocalDateTime,
			Double> obsQ, double obsError, boolean absoluteError, int ensembleSize,
			boolean resample, boolean perturb, boolean fClassKernels, String outputFolder,
			String modelsFolder, String forecastFolder) throws IOException
	{
		// Prepare data
		this.modelsFolder					= modelsFolder;
		this.initialStates					= initialStates;
		toDelete							= new HashSet<>();
		DateTimeFormatter formatter			= DateTimeFormatter.ofPattern(DATE_TIME_FORMAT_FOLDER);
		
		// Create mean streamflow file
		String meanQFile					= outputFolder + MEAN_Q_FILE;
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(meanQFile, false)));
		out.println("Date time\tObserved\tMean streamflow\tSt. dev.");
		out.close();
		
		// Initialize states
		ArrayList<Particle> currentState	= new ArrayList<>();
		MultiVarKernelDensity dist			= new MultiVarKernelDensity();
		int index							= 1;
		for (State state : initialStates)
		{
			ArrayList<Double> values		= state.toArray();
			dist.addSample(new Sample(1.0, values));
			currentState.add(new Particle("Root " + index, values, 1.0));
			index++;
		}
		
		// Generate additional states
		if (fClassKernels)
			dist.computeGaussianBW(true);
		else
			dist.computeGaussianDiagBW(true);
		ArrayList<Sample> samples = dist.sampleMultipleOb(ensembleSize - initialStates.size());
		index								= 1;
		for (Sample sample : samples)
			currentState.add(new Particle("Generated " + (index++), sample.getValues(), 1.0));
		
		// Perform sequential assimilation
		current								= start;
		while (current.isBefore(end))
		{
			// Obtain observation
			double observed					= obsQ.get(current.plus(modelTimeStep));
			double stDev					= absoluteError ? obsError : obsError*observed;
			Normal obs						= new Normal(observed, stDev);
			
			// Create folder
			String dateTime					= formatter.format(current);
			String stepFolder				= modelsFolder + "/" + dateTime;
			if (!Files.exists(FileSystems.getDefault().getPath(stepFolder)))
				Files.createDirectory(FileSystems.getDefault().getPath(stepFolder));
			System.out.println(dateTime);
			
			// Run assimilator
			currentStreamflow				= new Hashtable<>();
			Assimilator assimilator			= new Assimilator(this, obs);
			currentState					= assimilator.assimilate(currentState, ensembleSize, 
															resample, perturb, fClassKernels);
			current							= current.plus(modelTimeStep);
			
			// Store streamflow
			dateTime						= formatter.format(current);
			ContSeries streamflow			= new ContSeries(true);
			for (Particle particle : currentState)
			{
				String id					= particle.getId();
				String[] tokens				= id.split(" ");
				String id2					= tokens[0] + " " + tokens[1];
				double q					= currentStreamflow.get(id2);
				double weight				= particle.getWeight();
				streamflow.addValue(q, weight);
			}
			out = new PrintWriter(new BufferedWriter(new FileWriter(meanQFile, true)));
			out.println(dateTime + "\t" + obsQ.get(current) + "\t" + streamflow.getMean() 
							+ "\t" + streamflow.getStDev());
			out.close();
		}
		
		// Create final state
		ArrayList<State> endState			= new ArrayList<>(currentState.size());
		for (Particle particle : currentState)
			endState.add(new State(particle.getState(), current, initialStates.get(0)));		
		return endState;
	}
	
	public ArrayList<State> forecast(ArrayList<State> initialState, LocalDateTime start,
			LocalDateTime end, String folder, boolean removeFolder) throws IOException
	{
		// Create forecast folders
		System.out.println("");
		Files.createDirectory(FileSystems.getDefault().getPath(folder));
		String file				= folder + "/Streamflow.txt";
		PrintWriter out			= new PrintWriter(new BufferedWriter(new FileWriter(file, false)));
		out.println("Model\tWeight\tDischarge (m3/s)");
		
		// Prepare distributions
		ArrayList<State> endState		= new ArrayList<>(initialState.size());
		endStreamflow					= new KernelDensity();
		endEvaporation					= new KernelDensity();
		endSoilMoisture					= new KernelDensity();
		
		// Perform forecast
		int index						= 1;
		for (State state : initialState)
		{
			// Prepare model
			String id					= "Particle " + index;
			String runFolder			= folder + "/" + id;
			try
			{
				Simulation simulation	= new Simulation(runFolder, state, start, end, 
										modelTimeStep, parametersFolder, soils, network, areas, 
										outputs, directFractions, cellForcings, globalFileParams);
				simulation.run(vicExec, forecastSimMaxTime);
				
				// Retrieve and store hydrograph
				ArrayList<Double> streamF	= simulation.getStreamflow(outputs.get(0));
				ArrayList<Double> evap		= simulation.getEvaporation();
				ArrayList<Double> soilM		= simulation.getSoilMoisture();
				String line				= id + "\t" + 1.0 + "\t";
				for (Double value : streamF)
					line				+= value + "\t";
				out.println(line);
				
				// Store final values
				int endIndex			= streamF.size() - 1;
				endState.add(simulation.getEndState());
				endStreamflow.addSample(	streamF.get(	endIndex));
				endEvaporation.addSample(	evap.get(		endIndex));
				endSoilMoisture.addSample(	soilM.get(		endIndex));
				
				System.out.println("Completed forecast for " + id);
			} catch (Exception e)
			{
				System.out.println("Forecast for " + id + " failed: " + e.getMessage());
				//e.printStackTrace();
			} finally
			{
				// Stop executable
				String processName		= vicExec.substring(vicExec.lastIndexOf('/') + 1, 
													vicExec.length());
				try
				{
					Runtime.getRuntime().exec("taskkill /F /IM \"" + processName + "\"");
				} catch (IOException e1)
				{
					e1.printStackTrace();
				}
				index++;
			}
		}
		out.close();
		
		// Remove file folder
		if (removeFolder)
		{
			try
			{
				FileUtils.deleteDirectory(new File(folder));
			} catch (Exception e)
			{
				e.printStackTrace();
			}
		}
		
		return endState;
	}
	
	public KernelDensity getFinalStreamflow()
	{
		return endStreamflow;
	}
	
	public KernelDensity getFinalEvaporation()
	{
		return endEvaporation;
	}
	
	public KernelDensity getFinalSoilMoisture()
	{
		return endSoilMoisture;
	}
	
	@Override
	public ModelResult runModel(int index, ArrayList<Double> sourceValues)
	{
		// Prepare data
		LocalDateTime start				= current;
		Duration timeStep				= modelTimeStep;
		DateTimeFormatter formatter		= DateTimeFormatter.ofPattern(DATE_TIME_FORMAT_FOLDER);
		String runFolder				= modelsFolder + "/" + formatter.format(start) 
												+ "/Particle " + index;
		State sourceState				= new State(sourceValues, start, initialStates.get(0));
		String error					= "";
		State targetState				= null;
		ArrayList<Double> targetStateArray = null;
		double streamflow				= Double.NaN;
		
		// Run model
		LocalDateTime end				= start.plus(timeStep);
		Simulation simulation			= new Simulation(runFolder, sourceState, start, end, 
														modelTimeStep, parametersFolder, soils,
														network, areas,	outputs, directFractions, 
														cellForcings, globalFileParams);
		try
		{
			simulation.run(vicExec, simMaxTime);
			streamflow					= simulation.getStreamflow(outputs.get(0)).get(0);
			targetState					= simulation.getEndState();
			targetStateArray			= targetState.toArray();
		
			// Store current streamflow
			String id					= Assimilator.PREFIX + " " + index;
			currentStreamflow.put(id, streamflow);
		} catch (Exception e)
		{
			error						= e.getMessage();
			//e.printStackTrace();
			//if (error != null) System.out.println("Could not create \"" + id + "\": " + error);
		}
		
		// Stop executable
		String processName				= vicExec.substring(vicExec.lastIndexOf('/') + 1, 
											vicExec.length());
		try
		{
			Runtime.getRuntime().exec("taskkill /F /IM \"" + processName + "\"");
		} catch (IOException e1)
		{
			e1.printStackTrace();
		}
		
		// Remove simulation files
		if (removeFiles)
		{
			String deleteFolder			= runFolder;
			try
			{
				FileUtils.deleteDirectory(new File(deleteFolder));
				for (String folder : toDelete)
				{
					deleteFolder		= folder;
					FileUtils.deleteDirectory(new File(deleteFolder));
				}
			} catch (Exception e)
			{
				toDelete.add(deleteFolder);
			}
		}
		
		// Return
		return new ModelResult(targetStateArray, streamflow, error);
	}
	
}
