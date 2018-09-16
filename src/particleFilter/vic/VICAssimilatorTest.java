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
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.time.Duration;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Scanner;
import java.util.TreeSet;

import org.apache.commons.io.FileUtils;

import probDist.KernelDensity;
import utilities.geom.Point2D;
import vic.Forcing;
import vic.Soil;
import vic.routing.MuskingumNetwork;
import vic.routing.State;

/**
 * Citation: Hernández, F., & Liang, X. (2018). "Hybridizing Bayesian and variational data
 * assimilation for high-resolution hydrologic forecasting". Hydrol. Earth Syst. Sci.
 * 
 * @author Felipe Hernández (developer)
 * @author Xu Liang (advisor)
 */
public class VICAssimilatorTest
{
	
	// --------------------------------------------------------------------------------------------
	// Constants
	// --------------------------------------------------------------------------------------------
	
	public final static String STATES_FOLDER			= "Base states";
	public final static String MODELS_PREP_FOLDER		= "Preparation";
	public final static String DA_FOLDER 				= "Data assimilation";

	public final static String STATE_FILE_DT_FORMAT		= "yyyyMMdd HH-mm";

	public final static String OUT_FILE_STATS			= "Stats.txt";
	public final static String OUT_FILE_Q_VALS			= "Q - values.txt";
	public final static String OUT_FILE_EV_VALS			= "Ev - values.txt";
	public final static String OUT_FILE_SM_VALS			= "SM - values.txt";

	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	private VICAssimilator			assimilator;
	
	private Duration				modelTimeStep;
	private boolean					removeDAFiles;
	private boolean					removeDAModelFiles;
	private boolean					removeForecastFiles;
	
	private String					modelsFolder;
	private String					outputFolder;
	
	private DateTimeFormatter		formatter;
	
	private State					modelState;
	private TreeSet<LocalDateTime>	stateList;					
	private String					statesFolder;
	private String					stateFileHeader;
	
	private HashSet<String>			toDelete;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	public VICAssimilatorTest(String parameterFolder, ArrayList<Soil> soils, 
			ArrayList<String> globalFileParams, ArrayList<Forcing> cellForcings,
			MuskingumNetwork network, Duration modelTimeStep, Hashtable<String, Double> areas,
			ArrayList<String> outputs, Hashtable<String, Double> directFractions, String vicExec,
			long simMaxTime, long forecastSimMaxTime, boolean removeDAFiles,
			boolean removeDAModelFiles, boolean removeForecastFiles,
			Hashtable<LocalDateTime, Double> obsQ)
	{
		this.modelTimeStep			= modelTimeStep;
		this.removeDAFiles			= removeDAFiles;
		this.removeDAModelFiles		= removeDAModelFiles;
		this.removeForecastFiles	= removeForecastFiles;
		assimilator					= new VICAssimilator(parameterFolder, soils, globalFileParams,
							modelTimeStep, cellForcings, network, areas, outputs, directFractions,
							vicExec, simMaxTime, forecastSimMaxTime, removeDAFiles);
		formatter					= DateTimeFormatter.ofPattern(STATE_FILE_DT_FORMAT);
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------
	
	public void testAssimilator(int runIndex, String outputFolder, String modelsFolder,
			LocalDateTime forecastStart, LocalDateTime forecastEnd, ArrayList<Duration> leadTimes,
			ArrayList<State> baseState, ArrayList<String> variables, Hashtable<LocalDateTime,
			Double> obsQ, double obsError, boolean absoluteError, int ensembleSize,
			boolean resample, boolean perturb, boolean fClassKernels) throws IOException
	{
		this.outputFolder					= outputFolder;
		this.modelsFolder					= modelsFolder;
		
		modelState							= baseState.get(0);
		toDelete							= new HashSet<>();
		
		// Initialize report files
		Collections.sort(leadTimes);
		for (Duration leadTime : leadTimes)
		{
			// Create folder
			String ltStr					= leadTime.toString();
			String ltFolder					= outputFolder + "/Lead time = " + ltStr;
			Files.createDirectory(FileSystems.getDefault().getPath(ltFolder));
			
			// Create streamflow files
			String mFile					= ltFolder + "/" + OUT_FILE_STATS;
			String qFile					= ltFolder + "/" + OUT_FILE_Q_VALS;
			String eFile					= ltFolder + "/" + OUT_FILE_EV_VALS;
			String sFile					= ltFolder + "/" + OUT_FILE_SM_VALS;
			PrintWriter out	= new PrintWriter(new BufferedWriter(new FileWriter(mFile, false)));
			out.println("Date-time\tQ_Mean\tQ_0.05q\tQ_0.25q\tQ_median\tQ_0.75q\tQ_0.95q"
					+ "\tEv_Mean\tEv_0.05q\tEv_0.25q\tEv_median\tEv_0.75q\tEv_0.95q"
					+ "\tSM_Mean\tSM_0.05q\tSM_0.25q\tSM_median\tSM_0.75q\tSM_0.95q");
			out.close();
			out				= new PrintWriter(new BufferedWriter(new FileWriter(qFile, false)));
			out.println("Date-time\tStreamflow values");
			out.close();
			out				= new PrintWriter(new BufferedWriter(new FileWriter(eFile, false)));
			out.println("Date-time\tEvaporation values");
			out.close();
			out				= new PrintWriter(new BufferedWriter(new FileWriter(sFile, false)));
			out.println("Date-time\tSoil moisture values");
			out.close();
		}
		
		// Initialize state folders
		stateList							= new TreeSet<>();
		statesFolder						= outputFolder + "/" + STATES_FOLDER;
		String prepFolder					= modelsFolder + "/" + MODELS_PREP_FOLDER;
		StringBuilder builder				= new StringBuilder("Id\tWeight");
		for (String var : variables)
			builder.append("\t" + var);
		stateFileHeader						= builder.toString();
		Files.createDirectory(FileSystems.getDefault().getPath(prepFolder));
		Files.createDirectory(FileSystems.getDefault().getPath(statesFolder));
		writeStateToFile(baseState);
		
		// Perform test loop
		Files.createDirectory(FileSystems.getDefault().getPath(outputFolder + "/" + DA_FOLDER));
		LocalDateTime current				= forecastStart;
		while (!current.isAfter(forecastEnd))
		{
			// Prepare assimilation
			LocalDateTime daStart			= current.minus(modelTimeStep);
			LocalDateTime daEnd				= current;
			ArrayList<State> initState		= getState(daStart);
			String strDT					= formatter.format(current);
			String folder					= outputFolder + "/" + DA_FOLDER + "/" + strDT;
			String modelsiFolder			= modelsFolder + "/" + strDT;
			Files.createDirectory(FileSystems.getDefault().getPath(folder));
			Files.createDirectory(FileSystems.getDefault().getPath(modelsiFolder));
			
			// Assimilate and store target state
			System.out.println("\nPerforming assimilation at " + daStart + "...");
			ArrayList<State> targetState	= assimilator.assimilate(daStart, daEnd, null,
					initState, obsQ, obsError, absoluteError, ensembleSize, resample, perturb,
					fClassKernels, folder, modelsiFolder, "");
			writeStateToFile(targetState);
			
			// Delete assimilation model folders
			if (removeDAModelFiles)
			{
				try
				{
					FileUtils.deleteDirectory(new File(modelsiFolder));
				} catch (Exception e)
				{
					toDelete.add(folder);
				}
			}
			
			// Perform forecast
			LocalDateTime dateTime			= current;
			String forecastFolder			= folder + "/Forecast";
			Files.createDirectory(FileSystems.getDefault().getPath(forecastFolder));
			ArrayList<State> state			= targetState;
			for (Duration leadTime : leadTimes)
			{
				LocalDateTime end			= current.plus(leadTime);
				strDT						= formatter.format(dateTime);
				String forecastiFolder		= forecastFolder + "/" + leadTime.toString();
				System.out.println("\nPerforming forecast at " + dateTime + "...");
				state						= assimilator.forecast(state, dateTime, end,
														forecastiFolder, removeForecastFiles);
				KernelDensity streamflow	= assimilator.getFinalStreamflow();
				KernelDensity evaporation	= assimilator.getFinalEvaporation();
				KernelDensity soilMoisture	= assimilator.getFinalSoilMoisture();
				dateTime					= end;
				storeResults(dateTime, leadTime, state, streamflow, evaporation, soilMoisture);
			}
			
			// Delete folders
			if (removeDAFiles)
			{
				try
				{
					FileUtils.deleteDirectory(new File(folder));
				} catch (Exception e)
				{
					toDelete.add(folder);
				}
				try
				{
					FileUtils.deleteDirectory(new File(modelsiFolder));
				} catch (Exception e)
				{
					toDelete.add(modelsiFolder);
				}
			}
			
			// Advance to next time step
			current							= current.plus(modelTimeStep);
		}
		
		// Delete remaining folders
		if (removeDAFiles)
		{
			for (String folder : toDelete)
				try
				{
					FileUtils.deleteDirectory(new File(folder));
				} catch (Exception e) {}
		}
		
		System.out.println("\nFinished!");
	}
	
	private ArrayList<State> getState(LocalDateTime dateTime) throws IOException
	{
		String strDT				= formatter.format(dateTime);
		try
		{
			ArrayList<State> state	= readStateFromFile(dateTime);
			if (!stateList.contains(dateTime))
				stateList.add(dateTime);
			System.out.println("\nLoaded initial state (" + strDT + ")...");
			return state;
		} catch (Exception e)
		{
			// Continue
		}
		
		System.out.println("\nPreparing initial state for assimilation (" + strDT + ")...");
		LocalDateTime baseDT		= stateList.lower(dateTime);
		ArrayList<State> baseState	= readStateFromFile(baseDT);
		String folder		= modelsFolder + "/" + MODELS_PREP_FOLDER + "/" + strDT;
		//Files.createDirectory(FileSystems.getDefault().getPath(folder));
		ArrayList<State> state		= assimilator.forecast(baseState, baseDT, dateTime, folder,
										removeDAModelFiles);
		writeStateToFile(state);
		stateList.add(dateTime);
		return state;
	}

	private ArrayList<State> readStateFromFile(LocalDateTime dateTime)
	{
		DateTimeFormatter formatter	= DateTimeFormatter.ofPattern(STATE_FILE_DT_FORMAT);
		String fileName				= statesFolder + "/" + formatter.format(dateTime) + ".txt";
		Scanner scanner				= null;
		try
		{
			// Read samples from file
			scanner					= new Scanner(new FileInputStream(new File(fileName)));
			scanner.nextLine();
			ArrayList<State> states = new ArrayList<>();
			while (scanner.hasNextLine())
			{
				String[] tokens		= scanner.nextLine().split("\t");
				int valueCount		= tokens.length - 2;
				ArrayList<Double> values = new ArrayList<>(valueCount);
				for (int v = 2; v < valueCount + 2; v++)
					values.add(Double.valueOf(tokens[v]));
				State state		= new State(values, dateTime, modelState);
				states.add(state);
			}
			return states;
		} catch (Exception e)
		{
			throw new RuntimeException(e.getMessage());
		} finally
		{
			scanner.close();
		}
	}

	private void writeStateToFile(ArrayList<State> states) throws IOException
	{
		LocalDateTime dateTime		= states.get(0).dateTime;
		String fileName				= statesFolder + "/" + formatter.format(dateTime) + ".txt";
		PrintWriter out	= new PrintWriter(new BufferedWriter(new FileWriter(fileName, false)));
		out.println(stateFileHeader);
		int index					= 1;
		for (State state : states)
		{
			StringBuilder builder	= new StringBuilder(index + "\t1.0");
			for (Double value : state.toArray())
			{
				builder.append("\t");
				builder.append(value);
			}
			out.println(builder.toString());
		}
		out.close();
		stateList.add(dateTime);
	}

	private void storeResults(LocalDateTime dateTime, Duration leadTime, ArrayList<State> state,
			KernelDensity streamflow, KernelDensity evaporation, KernelDensity soilMoisture)
					throws IOException
	{
		// Get strings
		String ltStr		= leadTime.toString();
		String ltFolder		= outputFolder + "/Lead time = " + ltStr;
		String dtStr		= formatter.format(dateTime);
		
		//String mFile		= ltFolder + "/" + OUT_FILE_STATS;
		String qFile		= ltFolder + "/" + OUT_FILE_Q_VALS;
		String eFile		= ltFolder + "/" + OUT_FILE_EV_VALS;
		String sFile		= ltFolder + "/" + OUT_FILE_SM_VALS;
		
		// Write streamflow and weight values
		PrintWriter outQ	= new PrintWriter(new BufferedWriter(new FileWriter(qFile, true)));
		String lineQ		= dtStr;
		ArrayList<Point2D> samples = streamflow.getSamplesWeights();
		for (Point2D sample : samples)
		{
			lineQ			+= "\t" + sample.x;
		}
		outQ.println(lineQ);
		outQ.close();
		
		// Write evaporation values
		PrintWriter outE	= new PrintWriter(new BufferedWriter(new FileWriter(eFile, true)));
		String lineE		= dtStr;
		samples				= evaporation.getSamplesWeights();
		for (Point2D sample : samples)
			lineE			+= "\t" + sample.x;
		outE.println(lineE);
		outE.close();
		
		// Write soil moisture values
		PrintWriter outS	= new PrintWriter(new BufferedWriter(new FileWriter(sFile, true)));
		String lineS		= dtStr;
		samples				= soilMoisture.getSamplesWeights();
		for (Point2D sample : samples)
			lineS			+= "\t" + sample.x;
		outS.println(lineS);
		outS.close();
	}

}
