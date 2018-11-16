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
 * Citation: Hernández, F. and Liang, X.: Hybridizing Bayesian and variational data assimilation
 * for high-resolution hydrologic forecasting, Hydrol. Earth Syst. Sci., 22, 5759-5779,
 * https://doi.org/10.5194/hess-22-5759-2018, 2018.
 * 
 * @author Felipe Hernández (developer)
 * @author Xu Liang (advisor)
 */
public class BlueRiver
{

	public static void main(String[] args) throws IOException
	{
		// General configuration
		int scenario					= 1;
		int runIndex					= 1;
		String outputFolder				= "data/Tests/Scenario " + scenario + "/0" + runIndex;
		String modelsFolder				= outputFolder + "/Models";
		String forecastFolder			= outputFolder;
		String inputDataFolder			= "data/Blue River/";
		String vicExec					= "data/VIC/vicNl.exe";
		
		// Filter parameters
		int ensembleSize				= 100;
		boolean resample				= true;
		boolean perturb					= true;
		boolean fClassKernels			= true;
		
		// Determine times
		Duration modelTimeStep			= Duration.ofDays(1);
		LocalDateTime start				= null;
		LocalDateTime end				= null;
		LocalDateTime forecastEnd		= null;
		if (scenario == 1)
		{
			start						= LocalDateTime.of(1996, 10, 15, 0, 0);
			end							= LocalDateTime.of(1996, 10, 29, 0, 0);
			forecastEnd					= LocalDateTime.of(1996, 11, 13, 0, 0);
		}
		else if (scenario == 2)
		{
			start						= LocalDateTime.of(1997,  1, 15, 0, 0);
			end							= LocalDateTime.of(1997,  1, 29, 0, 0);
			forecastEnd					= LocalDateTime.of(1997,  2, 13, 0, 0);
		}
		else if (scenario == 3)
		{
			start						= LocalDateTime.of(1997,  6,  1, 0, 0);
			end							= LocalDateTime.of(1997,  6, 15, 0, 0);
			forecastEnd					= LocalDateTime.of(1997,  6, 29, 0, 0);
		}
		
		// Load observed flow
		LocalDateTime qObsStart			= start.plus(modelTimeStep);
		Hashtable<LocalDateTime, Double> obsQ = loadObsHydrograph(inputDataFolder + "obsQ_" 
													+ scenario + ".txt", qObsStart, modelTimeStep);
		double obsError					= 0.2;
		boolean absoluteError			= false;
		
		// Simulation options
		long simMaxTime					= 8*1000;
		long forecastSimMaxTime			= 45*1000;
		boolean removeFiles				= true;
		
		// Read global parameters file and forcings
		String parameterFolder			= inputDataFolder + "Parameters";
		String soilFile					= parameterFolder + "/soil.dat";
		ArrayList<Soil> soils			= Soil.readFromFile(soilFile, 3);
		String paramFile				= inputDataFolder + "vic_global_file_val";
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
		
		// Initial states
		ArrayList<State> initialStates	= new ArrayList<>();
		String initStateFolder			= inputDataFolder + "/States/scenario " + scenario + "/";
		initialStates.add(new State(initStateFolder + "01", initStateFolder + "01_rout.txt"));
		initialStates.add(new State(initStateFolder + "02", initStateFolder + "02_rout.txt"));
		initialStates.add(new State(initStateFolder + "03", initStateFolder + "03_rout.txt"));
		initialStates.add(new State(initStateFolder + "04", initStateFolder + "04_rout.txt"));
		initialStates.add(new State(initStateFolder + "05", initStateFolder + "05_rout.txt"));
		initialStates.add(new State(initStateFolder + "06", initStateFolder + "06_rout.txt"));
		initialStates.add(new State(initStateFolder + "07", initStateFolder + "07_rout.txt"));
		initialStates.add(new State(initStateFolder + "08", initStateFolder + "08_rout.txt"));
		
		// Perform assimilation
		VICAssimilator assimilator		= new VICAssimilator(parameterFolder, soils, 
									globalFileParams, modelTimeStep, cellForcings, network, areas, 
									outputs, directFractions, vicExec, simMaxTime, 
									forecastSimMaxTime, removeFiles);
		assimilator.assimilate(start, end, forecastEnd, initialStates, obsQ, obsError, 
									absoluteError, ensembleSize, resample, perturb, 
									fClassKernels, outputFolder, modelsFolder, forecastFolder);
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
