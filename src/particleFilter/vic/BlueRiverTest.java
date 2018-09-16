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

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.time.Duration;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Scanner;

import vic.Forcing;
import vic.Soil;
import vic.routing.MuskingumElement;
import vic.routing.MuskingumNetwork;
import vic.routing.State;

/**
 * Citation: Hernández, F., & Liang, X. (2018). "Hybridizing Bayesian and variational data
 * assimilation for high-resolution hydrologic forecasting". Hydrol. Earth Syst. Sci.
 * 
 * @author Felipe Hernández (developer)
 * @author Xu Liang (advisor)
 */
public class BlueRiverTest
{

	public static void main(String[] args) throws IOException
	{
		// General
		int runIndex					= 1;
		String outputFolder				= "data/Test " + runIndex + "/Output";
		String modelsFolder				= "data/Test " + runIndex + "/Models";
		String inputDataFolder			= "data/Blue River";
		boolean removeDAFiles			= false;
		boolean removeDAModelFiles		= true;
		boolean removeForecastFiles		= false;
		
		// Filter parameters
		int ensembleSize				= 30;
		boolean resample				= true;
		boolean perturb					= true;
		boolean fClassKernels			= false;
		
		// Determine times
		LocalDateTime forecastStart		= LocalDateTime.of(1996, 12, 15, 0, 0);
		LocalDateTime forecastEnd		= LocalDateTime.of(1997,  6,  1, 0, 0);
		Duration modelTimeStep			= Duration.ofDays(1);
		ArrayList<Duration> leadTimes	= new ArrayList<>();
		leadTimes.add(Duration.ofDays(1));
		leadTimes.add(Duration.ofDays(3));
		leadTimes.add(Duration.ofDays(6));
		leadTimes.add(Duration.ofDays(12));
		LocalDateTime baseStateTime		= LocalDateTime.of(1996, 12, 1, 0, 0);
		String initStateFile			= inputDataFolder + "/States/19961201 00-00.txt";
		
		// Simulation options
		String vicExec					= "data/VIC/vicNl.exe";
		long simMaxTime					= 6*1000;
		long forecastSimMaxTime			= 6*1000;
		
		// Load observed flow
		LocalDateTime qObsStart			= LocalDateTime.of(1996, 10, 16, 0, 0);
		Hashtable<LocalDateTime, Double> obsQ = loadObsHydrograph(inputDataFolder + "/obsQ.txt",
																	qObsStart, modelTimeStep);
		double obsError					= 0.2;
		boolean absoluteError			= false;
		
		// Read global parameters file and forcings
		String parameterFolder			= inputDataFolder + "/Parameters";
		String soilFile					= parameterFolder + "/soil.dat";
		ArrayList<Soil> soils			= Soil.readFromFile(soilFile, 3);
		String paramFile				= inputDataFolder + "/vic_global_file_val";
		ArrayList<String> globalFileParams = new ArrayList<>();
		Scanner scanner					= new Scanner(new FileInputStream(new File(paramFile)));
		while (scanner.hasNextLine())
			globalFileParams.add(scanner.nextLine());
		scanner.close();
		ArrayList<Forcing> cellForcings	= Forcing.loadFromFiles(inputDataFolder + "/Forcing");
		
		// Read network file
		String routingFile				= parameterFolder + "/routing.txt";
		Hashtable<String, Double> areas = new Hashtable<>();
		HashSet<String> inNetwork		= new HashSet<>();
		ArrayList<String> outputs		= new ArrayList<>();
		Hashtable<String, Double> directFractions = new Hashtable<>();
		MuskingumNetwork network		= new MuskingumNetwork();
		scanner							= new Scanner(new FileInputStream(new File(routingFile)));
		while (scanner.hasNextLine())
		{
			String line					= scanner.nextLine().trim();
			String[] tokens				= line.split("\t");
			String id					= tokens[0];
			double k					= Double.valueOf(tokens[1]);
				double x				= Double.valueOf(tokens[2]);
			String downstream			= null;
			if (tokens.length > 3)
			{
				downstream				= tokens[3];
				if (tokens.length > 4)
				{
					double area			= Double.valueOf(tokens[4]);
					areas.put(id, area);
					if (tokens.length > 5)
					{
						directFractions.put(id, Double.valueOf(tokens[5]));
						if (tokens.length > 6)
							if (Boolean.valueOf(tokens[6]))
								outputs.add(id);
					}
				if (!inNetwork.contains(downstream))
					downstream			= null;
				}
			}
			network.addElement(id, new MuskingumElement(k, x, 0.0, 0.0), downstream);
			inNetwork.add(id);
		}
		scanner.close();
		
		// Setup initial ensemble
		String exStateFolder			= inputDataFolder + "/states/scenario 1/";
		State modelState = new State(exStateFolder + "01", exStateFolder + "01_rout.txt");
		ArrayList<String> variables		= modelState.getVariableNames();
		ArrayList<State> baseState		= createInitialState(initStateFile, baseStateTime,
											modelState, variables);
		
		VICAssimilatorTest test			= new VICAssimilatorTest(parameterFolder, soils,
				globalFileParams, cellForcings, network, modelTimeStep, areas, outputs,
				directFractions, vicExec, simMaxTime, forecastSimMaxTime, removeDAFiles,
				removeDAModelFiles, removeForecastFiles, obsQ);
		test.testAssimilator(runIndex, outputFolder, modelsFolder, forecastStart, forecastEnd,
								leadTimes, baseState, variables, obsQ, obsError, absoluteError,
								ensembleSize, resample, perturb, fClassKernels);
	}
	
	private static ArrayList<State> createInitialState(String initStateFile,
			LocalDateTime dateTime, State modelState, ArrayList<String> variables)
					throws FileNotFoundException
	{
		Scanner scanner				= new Scanner(new FileInputStream(new File(initStateFile)));
		String line					= scanner.nextLine();
		
		// Verify variables match
		String[] tokens				= line.split("\t");
		for (int i = 0; i < variables.size(); i++)
		{
			String varName			= variables.get(i);
			String onFile			= tokens[i + 2];
			if (varName.compareTo(onFile) != 0)
			{
				scanner.close();
				throw new IllegalArgumentException("Variable name mismatch: \"" 
													+ varName + "\", \"" + onFile + "\"");
			}
		}
		
		ArrayList<State> states		= new ArrayList<>();
		while (scanner.hasNextLine())
		{
			line					= scanner.nextLine();
			tokens					= line.split("\t");
			ArrayList<Double> values = new ArrayList<>();
			for (int t = 2; t < tokens.length; t++)
				values.add(Double.valueOf(tokens[t]));
			states.add(new State(values, dateTime, modelState));
		}
		scanner.close();
		
		return states;
	}

	private static Hashtable<LocalDateTime, Double> loadObsHydrograph(String fileRoute, 
			LocalDateTime start, Duration timeStep) throws FileNotFoundException
	{
		Hashtable<LocalDateTime, Double> qObs = new Hashtable<>();
		Scanner scanner					= new Scanner(new FileInputStream(new File(fileRoute)));
		LocalDateTime dateTime			= start;
		while (scanner.hasNextLine())
		{
			Double value				= Double.valueOf(scanner.nextLine());
			qObs.put(dateTime, value);
			dateTime					= dateTime.plus(timeStep);
		}
		scanner.close();
		return qObs;
	}

}
