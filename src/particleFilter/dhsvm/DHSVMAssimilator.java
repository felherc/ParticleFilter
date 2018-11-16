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

package particleFilter.dhsvm;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Duration;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Scanner;
import java.util.concurrent.ArrayBlockingQueue;

import org.apache.commons.io.FileUtils;

import dhsvm.MetStation;
import dhsvm.Soil;
import dhsvm.Vegetation;
import dhsvm.grid.Input;
import dhsvm.grid.State;
import dhsvm.stream.StreamNetwork;
import particleFilter.Assimilator;
import particleFilter.Model;
import particleFilter.ModelResult;
import particleFilter.Particle;
import probDist.KernelDensity;
import probDist.Normal;
import probDist.multiVar.NonParametric;
import probDist.multiVar.tools.ContMultiSample;
import probDist.multiVar.tools.Sample;
import utilities.Utilities;
import utilities.geom.PointSD;
import utilities.stat.ContSeries;
import utilities.thread.Executor;
import utilities.thread.ExecutorThread;

/**
 * Citation: Hernández, F. and Liang, X.: Hybridizing Bayesian and variational data assimilation
 * for high-resolution hydrologic forecasting, Hydrol. Earth Syst. Sci., 22, 5759-5779,
 * https://doi.org/10.5194/hess-22-5759-2018, 2018.
 * 
 * @author Felipe Hernández (developer)
 * @author Xu Liang (advisor)
 */
public class DHSVMAssimilator implements Model, Executor
{

	// --------------------------------------------------------------------------------------------
	// Constants
	// --------------------------------------------------------------------------------------------
	
	public final static String	MEAN_Q_FILE					= "/Streamflow.txt";
	public final static String	PARTICLE_ID_PREFIX			= "Particle";
	public final static String	DATE_TIME_FORMAT_FOLDER		= "yyyy-MM-dd HH.mm";
	public final static String	DATE_TIME_FORMAT_FLOW_1		= "MM.dd.yyyy-HH:mm:ss";
	public final static String	DATE_TIME_FORMAT_FLOW_2		= "MM/dd/yyyy-HH:mm:ss";
	public final static String	DATE_TIME_FORMAT_REPORT		= "MM/dd/yyyy HH:mm";
	
	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	private ModelConfigurator						configurator;
	private DHSVMModel								defaultModel;
	private Input									input;
	private StreamNetwork							network;
	private boolean									defaultParameters;
	private ArrayList<String>						metFiles;
	private Duration								modelTimeStep;
	
	private String									modelsFolder;
	private String									dhsvmExec;
	
	private Duration								daTimeStep;
	
	private Hashtable<LocalDateTime, Double>		obsQ;
	
	private ArrayList<Particle>						currentStates;
	private LocalDateTime							current;
	private Hashtable<String, Double>				currentStreamflow;
	
	private String 									meanQFile;
	private LocalDateTime							forecastEnd;
	private String									forecastFolder;
	private ArrayBlockingQueue<Particle>			forecastQueue;
	private Hashtable<String, String>				forecastQ;
	private ArrayList<ContMultiSample>				forecastEndStates;
	private ArrayList<LocalDateTime>				statesToSave;
	private Hashtable<LocalDateTime, KernelDensity>	forecastStreamflow;
	private Hashtable<LocalDateTime, KernelDensity>	forecastEvaporation;
	private Hashtable<LocalDateTime, KernelDensity>	forecastSoilMoistureL1;
	private Hashtable<LocalDateTime, KernelDensity>	forecastSoilMoistureL2;
	private Hashtable<LocalDateTime, KernelDensity>	forecastSoilMoistureL3;
	
	private boolean									removeForecastFiles;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	public DHSVMAssimilator(ModelConfigurator configurator, Input input, ArrayList<Soil> soils,
			ArrayList<Vegetation> vegetations, StreamNetwork network,
			ArrayList<MetStation> stations, boolean defaultParameters, String optionsFile,
			String areaFile, String constantsFile,  ArrayList<String> metFiles,
			Duration modelTimeStep, String dhsvmExec) throws FileNotFoundException
	{
		this.configurator		= configurator;
		this.input				= input;
		this.network			= network;
		this.defaultParameters	= defaultParameters;
		this.metFiles			= metFiles;
		this.modelTimeStep		= modelTimeStep;
		this.dhsvmExec			= dhsvmExec;
		
		// Copy meteorological files
		metFiles					= new ArrayList<>();
		for (MetStation station : stations)
		{
			String file				= station.dataFile;
			metFiles.add(file);
			station.dataFile		= "../../Met" + file.substring(file.lastIndexOf("/"), 
																		file.length());
		}
		
		// Read configuration sections' files
		String[] optionsSection		= readSectionFile(optionsFile);
		String[] areaSection		= readSectionFile(areaFile);
		String[] constantsSection	= readSectionFile(constantsFile);
		
		// Copy inputs and create model
		String inputFolderCopy		= "../../Input";
		String demFile				= inputFolderCopy + "/DEM.bin";
		String maskFile				= inputFolderCopy + "/Mask.bin";
		String flowDirFile			= inputFolderCopy + "/Flow_dir.bin";
		String soilTypeFile			= inputFolderCopy + "/Soil.bin";
		String soilDepthFile		= inputFolderCopy + "/Soil_depth.bin";
		String vegTypeFile			= inputFolderCopy + "/Vegetation.bin";
		String streamClassFile		= inputFolderCopy + "/stream_class.txt";
		String streamNetworkFile	= inputFolderCopy + "/stream_network.txt";
		String streamMapFile		= inputFolderCopy + "/stream_map.txt";
		String surfaceRoutingFile	= inputFolderCopy + "/surface_routing.txt";
		this.defaultModel			= new DHSVMModel(optionsSection, areaSection, 
				constantsSection, demFile, maskFile, flowDirFile, soilTypeFile, soilDepthFile, 
				vegTypeFile, streamClassFile, streamNetworkFile, streamMapFile, 
				surfaceRoutingFile, true, soils, vegetations, stations);
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------

	public ArrayList<Particle> assimilate(NonParametric initialDist, String outputFolder,
			String modelsFolder, LocalDateTime start, LocalDateTime end, Duration daTimeStep,
			Hashtable<LocalDateTime, Double> obsQ, boolean absoluteError, double obsError,
			int ensembleSize, boolean resample, boolean perturb, boolean fClassKernels)
					throws IOException
	{
		// Prepare data and prepare models folder
		this.modelsFolder					= modelsFolder;
		this.daTimeStep						= daTimeStep;
		this.obsQ							= obsQ;
		DateTimeFormatter formatter			= DateTimeFormatter.ofPattern(DATE_TIME_FORMAT_FOLDER);
		createModelsFolders(modelsFolder);
		
		// Initialize states
		this.currentStates					= new ArrayList<>();
		int index							= 1;
		for (ContMultiSample sample : initialDist.getSamples(true, 1.0, ensembleSize))
		{
			ArrayList<Double> values		= sample.getValues();
			currentStates.add(new Particle("Root " + index, values, 1.0));
			index++;
		}
		int missing							= currentStates.size() - ensembleSize;
		if (missing > 0)
		{
			double[][] samples				= initialDist.sampleMultiple(missing);
			index							= 1;
			for (int s = 0; s < samples.length; s++)
			{
				ArrayList<Double> values	= Utilities.toArrayList(samples[s]);
				currentStates.add(new Particle("Generated " + index, values, 1.0));
				index++;
			}
		}
		
		// Create mean streamflow file
		String meanQFile					= outputFolder + MEAN_Q_FILE;
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(meanQFile, false)));
		out.println("Date time\tObserved\tMean streamflow\tSt. dev.");
		out.close();

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
			currentStates					= assimilator.assimilate(currentStates, ensembleSize,
																resample, perturb, fClassKernels);
			current							= current.plus(daTimeStep);
			
			// Store streamflow
			ContSeries streamflow			= new ContSeries(true);
			for (Particle particle : currentStates)
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

		System.out.println();
		return currentStates;
	}
	
	@Override
	public ModelResult runModel(int index, ArrayList<Double> valueArray)
	{
		// Prepare data
		String particleID							= PARTICLE_ID_PREFIX + " " + index;
		State initState						= null;
		LocalDateTime start					= current;
		LocalDateTime end					= current.plus(daTimeStep);
				
		// Prepare files
		DateTimeFormatter formatter			= DateTimeFormatter.ofPattern(DATE_TIME_FORMAT_FOLDER);
		String runFolder					= modelsFolder + "/" + formatter.format(start) 
												+ "/Particle " + index;
		String stateFolder					= runFolder + "/state";
		String outputFolder					= runFolder + "/output";
		String configFile					= runFolder + "/Configuration.txt";
		
		// Configure model
		ModelConfiguration config			= configurator.configure(valueArray, start,
																		defaultParameters);
		initState							= config.initialState;
		ArrayList<Soil> soils				= config.soils;
		ArrayList<Vegetation> vegetations	= config.vegetations;
		StreamNetwork network				= config.network;
		String classFileName				= "stream_class.txt";
		String networkFileName				= "stream_network.txt";
		String mapFileName					= "stream_map.txt";
		String surfaceRoutingFile			= "surface_routing.txt";
		DHSVMModel model					= new DHSVMModel(defaultModel.optionsSection,
				defaultModel.areaSection, defaultModel.constantsSection, defaultModel.demFile,
				defaultModel.maskFile, defaultModel.flowDirFile, defaultModel.soilTypeFile,
				defaultModel.soilDepthFile, defaultModel.vegTypeFile, classFileName,
				networkFileName, mapFileName, surfaceRoutingFile, defaultModel.d8FlowDirection,
				soils, vegetations, defaultModel.stations);
		
		// Write files
		try
		{
			if (Files.exists(FileSystems.getDefault().getPath(runFolder	)))
				FileUtils.deleteDirectory(new File(runFolder));
			
			Files.createDirectory(FileSystems.getDefault().getPath(runFolder	));
			Files.createDirectory(FileSystems.getDefault().getPath(stateFolder	));
			Files.createDirectory(FileSystems.getDefault().getPath(outputFolder	));
			network.writeToFiles(runFolder, classFileName, networkFileName, mapFileName,
									surfaceRoutingFile);
			
			initState.writeToFiles(stateFolder);
			model.writeConfigFile(configFile, start, end, modelTimeStep, "state", "output", true, null);
		} catch (IOException e)
		{
			throw new RuntimeException("Could not create model files: " + e.getMessage());
		}
		
		String error						= "";
		State endState						= null;
		ArrayList<Double> endArray			= null;
		Hashtable<LocalDateTime, Double> modeledQ = null;
		try
		{
			// Run model
			ProcessBuilder pb				= new ProcessBuilder(dhsvmExec, configFile);
			pb.directory(new File(runFolder));
			Process process					= pb.start();
			BufferedReader br				= new BufferedReader(new InputStreamReader(
																	process.getInputStream()));
			while (br.readLine() != null)
			{
				synchronized (this)	{try {wait(1);} catch (InterruptedException e) {}}
			}
			
			// Retrieve target state and streamflow
			int rows						= input.elevation.length;
			int cols						= input.elevation[0].length;
			endState						= new State(end, outputFolder, rows, cols);
			endArray						= configurator.toArray(new ModelConfiguration(
														soils, vegetations, network, endState));
			modeledQ 						= getHydrograph(outputFolder);
			currentStreamflow.put(particleID, modeledQ.get(end));
		} catch (IOException e)
		{
			String line						= formatter.format(start) + " - " + particleID + ": ";
			error							= e.getMessage();
			line							+= error;
			System.out.println(line);
		}
		
		double streamflow					= modeledQ != null ? modeledQ.get(end) : Double.NaN;
		return new ModelResult(endArray, streamflow, error);
	}
	
	private Hashtable<LocalDateTime, Double> getHydrograph(String outputFolder) 
																	throws FileNotFoundException
	{
		// Obtain modeled hydrograph
		String flowFile					= outputFolder + "/Stream.Flow";
		Scanner scanner					= new Scanner(new FileInputStream(new File(flowFile)));
		Hashtable<LocalDateTime, Double> hydrograph = new Hashtable<>();
		DateTimeFormatter formatter		= DateTimeFormatter.ofPattern(DATE_TIME_FORMAT_FLOW_1);
		while (scanner.hasNextLine())
		{
			String line					= scanner.nextLine();
			String[] tokens				= line.split("\\s+");
			LocalDateTime dateTime		= LocalDateTime.parse(tokens[0], formatter);
			double value				= Double.valueOf(tokens[4])/3.6;	// Convert to l/s
			hydrograph.put(dateTime, value);
		}
		scanner.close();
		return hydrograph;
	}
	
	private Hashtable<LocalDateTime, Double> getEvaporation(String outputFolder) 
																	throws FileNotFoundException
	{
		String outputFile				= outputFolder + "/Aggregated.Values";
		Scanner scanner					= new Scanner(new FileInputStream(new File(outputFile)));
		Hashtable<LocalDateTime, Double> evaporation = new Hashtable<>();
		DateTimeFormatter formatter		= DateTimeFormatter.ofPattern(DATE_TIME_FORMAT_FLOW_2);
		while (scanner.hasNextLine())
		{
			String line					= scanner.nextLine();
			String[] tokens				= line.split("\\s+");
			LocalDateTime dateTime		= LocalDateTime.parse(tokens[0], formatter);
			double value				= Double.valueOf(tokens[8]);
			evaporation.put(dateTime, value);
		}
		scanner.close();
		return evaporation;
	}
	
	private Hashtable<LocalDateTime, double[]> getSoilMoisture(String outputFolder) 
																	throws FileNotFoundException
	{
		String outputFile				= outputFolder + "/Aggregated.Values";
		Scanner scanner					= new Scanner(new FileInputStream(new File(outputFile)));
		Hashtable<LocalDateTime, double[]> soilMoisture = new Hashtable<>();
		DateTimeFormatter formatter		= DateTimeFormatter.ofPattern(DATE_TIME_FORMAT_FLOW_2);
		while (scanner.hasNextLine())
		{
			String line					= scanner.nextLine();
			String[] tokens				= line.split("\\s+");
			LocalDateTime dateTime		= LocalDateTime.parse(tokens[0], formatter);
			double[] sm					= new double[3];
			sm[0]						= Double.valueOf(tokens[30]);
			sm[1]						= Double.valueOf(tokens[31]);
			sm[2]						= Double.valueOf(tokens[32]);
			soilMoisture.put(dateTime, sm);
		}
		scanner.close();
		return soilMoisture;
	}
	
	private String[] readSectionFile(String textFile) throws FileNotFoundException
	{
		ArrayList<String> lines		= new ArrayList<>();
		Scanner scanner				= new Scanner(new File(textFile));
		while (scanner.hasNextLine())
			lines.add(scanner.nextLine());
		scanner.close();
		
		String[] array				= new String[lines.size()];
		for (int l = 0; l < lines.size(); l++)
			array[l]				= lines.get(l);
		return array;
	}
	
	public void prepareForecast(ArrayList<Particle> currentStates, LocalDateTime dateTime,
					String outputFolder, Hashtable<LocalDateTime, Double> obsQ) throws IOException
	{
		// Create mean streamflow file
		this.obsQ			= obsQ;
		meanQFile			= outputFolder + MEAN_Q_FILE;
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(meanQFile, false)));
		out.println("Date time\tObserved streamflow (l/s)\tMean streamflow (l/s)\tSt. dev. (l/s)");
		out.close();
		
		// Prepare distribution
		this.currentStates	= currentStates;
		current				= dateTime;
	}
	
	public ArrayList<Particle> forecast(LocalDateTime forecastEnd, String folder, int threadCount, 
						LocalDateTime performanceStart, ArrayList<LocalDateTime> statesToSave,
						boolean removeFolder, long timeLimit) throws IOException
	{
		System.out.println("Starting forecast...");
		this.forecastEnd			= forecastEnd;
		this.forecastFolder			= folder;
		this.removeForecastFiles	= removeFolder;
		if (statesToSave == null)
			statesToSave			= new ArrayList<>();
		if (statesToSave.size() == 0)
			statesToSave.add(forecastEnd);
		this.statesToSave			= statesToSave;
		if (!folder.equals(modelsFolder))
			createModelsFolders(folder);
		
		// Prepare distributions of forecasted outputs
		LocalDateTime forecastStart	= current.plus(modelTimeStep);
		LocalDateTime dateTime		= forecastStart;
		currentStreamflow			= new Hashtable<>();
		forecastStreamflow			= new Hashtable<>();
		forecastEvaporation			= new Hashtable<>();
		forecastSoilMoistureL1		= new Hashtable<>();
		forecastSoilMoistureL2		= new Hashtable<>();
		forecastSoilMoistureL3		= new Hashtable<>();
		while (!dateTime.isAfter(forecastEnd))
		{
			forecastStreamflow.put(		dateTime, new KernelDensity());
			forecastEvaporation.put(	dateTime, new KernelDensity());
			forecastSoilMoistureL1.put(	dateTime, new KernelDensity());
			forecastSoilMoistureL2.put(	dateTime, new KernelDensity());
			forecastSoilMoistureL3.put(	dateTime, new KernelDensity());
			dateTime				= dateTime.plus(modelTimeStep);
		}
		
		// Verify folders
		Path path					= FileSystems.getDefault().getPath(folder + "/Forecasts");
		if (!Files.exists(path))
			Files.createDirectory(path);
		path						= FileSystems.getDefault().getPath(folder + "/Input");
		if (!Files.exists(path))
			createModelsFolders(folder);
		
		// Populate forecast queue
		int sampleCount				= currentStates.size();
		forecastEndStates			= new ArrayList<>(sampleCount);
		forecastQueue				= new ArrayBlockingQueue<>(sampleCount);
		for (Particle particle : currentStates)
			forecastQueue.offer(particle);
		
		// Launch threads
		this.forecastQ				= new Hashtable<>();
		for (int t = 0; t < threadCount; t++)
		{
			ExecutorThread thread	= new ExecutorThread(this);
			thread.start("Perform forecast");
		}
		
		// Wait for threads to finish
		long start					= System.currentTimeMillis();
		boolean timeUp				= false;		
		while (forecastQ.size() < sampleCount && !timeUp)
			synchronized (this)
			{
				try
				{
					long time		= System.currentTimeMillis() - start;
					timeUp			= time > timeLimit;
					if (timeUp)
						System.out.println("Forecast time up!");
					wait(1);
				} catch (InterruptedException e) {}
			}
		forecastQueue.clear();
		
		// Compute bandwidth of distribution
		for (LocalDateTime dateTime2 : forecastStreamflow.keySet())
		{
			forecastStreamflow.get(		dateTime2).computeGaussianBandwidth();
			forecastEvaporation.get(	dateTime2).computeGaussianBandwidth();
			forecastSoilMoistureL1.get(	dateTime2).computeGaussianBandwidth();
			forecastSoilMoistureL2.get(	dateTime2).computeGaussianBandwidth();
			forecastSoilMoistureL3.get(	dateTime2).computeGaussianBandwidth();
		}
		
		// Write mean streamflow
		DateTimeFormatter formatter	= DateTimeFormatter.ofPattern(DATE_TIME_FORMAT_REPORT);
		PrintWriter out 	= new PrintWriter(new BufferedWriter(new FileWriter(meanQFile, true)));
		dateTime					= forecastStart;
		while (!dateTime.isAfter(forecastEnd))
		{
			KernelDensity dist		= forecastStreamflow.get(dateTime);
			out.println(formatter.format(dateTime) + "\t" + obsQ.get(dateTime) + "\t" 
									+ dist.getMean() + "\t" + dist.getStDev());
			dateTime				= dateTime.plus(modelTimeStep);
		}
		out.close();
		
		// Compute forecast performance metrics
		if (performanceStart != null)
		{
			// Obtain data
			dateTime					= performanceStart;
			ArrayList<Double> obs	= new ArrayList<>();
			ArrayList<Double> mod	= new ArrayList<>();
			ContSeries density		= new ContSeries(false);
			ContSeries rarity		= new ContSeries(false);
			while (!dateTime.isAfter(forecastEnd))
			{
				KernelDensity dist	= forecastStreamflow.get(dateTime);
				Double observed		= obsQ.get(dateTime);
				if (observed != null)
				{
					// Deterministic evaluation
					obs.add(observed);
					mod.add(dist.getMean());
					
					// Probabilistic evaluation
					density.addValue(dist.getpdf(observed));
					rarity.addValue(2*Math.abs(dist.getCDF(observed) - 0.5));
				}
				dateTime				= dateTime.plus(modelTimeStep);
			}
			
			// Compute performance
			double[] obsArr			= Utilities.toArray(obs);
			double[] modArr			= Utilities.toArray(mod);
			double nse_l2			= Utilities.computeNashSutcliffe(obsArr, modArr);
			double nse_l1			= Utilities.computeNashSutcliffe(obsArr, modArr, 1.0);
			double mare				= Utilities.computeMeanAbsRelativeError(obsArr, modArr);
			double meanDensity		= density.getMean();
			double meanRarity		= rarity.getMean();
			
			// Store performance
			String file				= folder + "/Forecasts/Performance.txt";
			out					= new PrintWriter(new BufferedWriter(new FileWriter(file, false)));
			out.println("Values =\t"			+ obsArr.length	);
			out.println("NSE_l2 =\t"			+ nse_l2		);
			out.println("NSE_l1 =\t"			+ nse_l1		);
			out.println("MARE =\t"				+ mare			);
			out.println("Mean density =\t"		+ meanDensity	);
			out.println("Mean rarity =\t"		+ meanRarity	);
			out.close();
			
			System.out.println("\nValues =\t"		+ obsArr.length	);
			System.out.println("NSE_l2 =\t"			+ nse_l2		);
			System.out.println("NSE_l1 =\t"			+ nse_l1		);
			System.out.println("MARE =\t"			+ mare			);
			System.out.println("Mean density =\t"	+ meanDensity	);
			System.out.println("Mean rarity =\t"	+ meanRarity	);
		}
		
		// Write forecasted streamflow file
		out						= null;
		String file				= folder + "/Forecasts/Streamflow.txt";
		out						= new PrintWriter(new BufferedWriter(new FileWriter(file, false)));
		out.println("Model\tWeight\tDischarge (l/s)");
		ArrayList<PointSD> sorter = new ArrayList<>();
		for (Particle particle : currentStates)
			sorter.add(new PointSD(particle.getId(), particle.getWeight(), false));
		Collections.sort(sorter);
		for (int p = sorter.size() - 1; p >= 0; p--)
			out.println(forecastQ.get(sorter.get(p).getX()));
		out.close();
		out.close();
		
		// Remove file folders
		if (removeFolder)
		{
			try
			{
				FileUtils.deleteDirectory(new File(folder + "/Forecasts"	));
				FileUtils.deleteDirectory(new File(folder + "/Input"		));
				FileUtils.deleteDirectory(new File(folder + "/Met"			));
			} catch (Exception e)
			{
				System.out.println("Error removing folder: " + e.getMessage());
			}
		}
		
		System.out.println("");
		
		// Return final distribution
		if (forecastEndStates.size() > 0)
		{
			try
			{
				synchronized (forecastEndStates)
				{
					return configurator.createDistribution(forecastEndStates);
				}
			} catch (Exception e)
			{
				e.printStackTrace();
				return null;
			}
		}
		else
		{
			System.out.println("No states to create final distribution");
			return null;
		}
	}
	
	private void createModelsFolders(String modelsFolder) throws IOException
	{
		String inputFolder					= modelsFolder + "/Input";
		String metFolder					= modelsFolder + "/Met";
		if (!Files.exists(			FileSystems.getDefault().getPath(inputFolder	)))
			Files.createDirectory(	FileSystems.getDefault().getPath(inputFolder	));
		if (!Files.exists(			FileSystems.getDefault().getPath(metFolder		)))
			Files.createDirectory(	FileSystems.getDefault().getPath(metFolder		));
		input.writeToFiles(inputFolder);
		network.writeToFiles(inputFolder, "stream_class.txt", "stream_network.txt", 
								"stream_map.txt", "surface_routing.txt");
		for (String origFile : metFiles)
		{
			String targetFile				= metFolder + origFile.substring(
												origFile.lastIndexOf("/"), origFile.length());
			File file						= new File(targetFile);
			if (!Files.exists(file.toPath()))
				Files.copy(new File(origFile).toPath(), file.toPath());
		}
	}
	
	@Override
	public void execute(String processID)	// Perform forecast
	{
		while (forecastQueue.size() > 0)
		{
			String runFolder				= null;
			try
			{
				Particle particle			= forecastQueue.poll();
				LocalDateTime end			= forecastEnd;
				
				// Prepare model
				String particleID			= particle.getId();
				ModelConfiguration config	= configurator.configure(particle.getState(), current, defaultParameters);
				State initialState			= config.initialState;
				LocalDateTime start			= initialState.dateTime;
				runFolder					= forecastFolder + "/Forecasts/" + particleID;
				String stateFolder			= runFolder + "/state";
				String outputFolder			= runFolder + "/output";
				String configFile			= runFolder + "/Configuration.txt";
				Files.createDirectory(FileSystems.getDefault().getPath(runFolder	));
				Files.createDirectory(FileSystems.getDefault().getPath(stateFolder	));
				Files.createDirectory(FileSystems.getDefault().getPath(outputFolder	));
				initialState.writeToFiles(stateFolder);
				
				String classFileName				= "stream_class.txt";
				String networkFileName				= "stream_network.txt";
				String mapFileName					= "stream_map.txt";
				String surfaceRoutingFile			= "surface_routing.txt";
				
				ArrayList<Soil> soils				= config.soils;
				ArrayList<Vegetation> vegetations	= config.vegetations;
				StreamNetwork network				= config.network;
				DHSVMModel model					= new DHSVMModel(defaultModel.optionsSection,
				defaultModel.areaSection, defaultModel.constantsSection, defaultModel.demFile,
				defaultModel.maskFile, defaultModel.flowDirFile, defaultModel.soilTypeFile,
				defaultModel.soilDepthFile, defaultModel.vegTypeFile, classFileName,
				networkFileName, mapFileName, surfaceRoutingFile, defaultModel.d8FlowDirection,
				soils, vegetations, defaultModel.stations);
				
				network.writeToFiles(runFolder, classFileName, networkFileName, mapFileName,
										surfaceRoutingFile);
				
				model.writeConfigFile(configFile, start, end, modelTimeStep, "state", "output", 
										false, statesToSave);
				
				// Run forecast
				ProcessBuilder pb		= new ProcessBuilder(dhsvmExec, configFile);
				pb.directory(new File(runFolder));
				Process process			= pb.start();
				BufferedReader br		= new BufferedReader(new InputStreamReader(
																		process.getInputStream()));
				while (br.readLine() != null)
				{
					synchronized (this)	{try {wait(1);} catch (InterruptedException e) {}}
				}
				
				// Retrieve target outputs
				String line				= particleID + "\t" + particle.getWeight() + "\t";
				Hashtable<LocalDateTime, Double> q		= getHydrograph(	outputFolder);
				Hashtable<LocalDateTime, Double> ev		= getEvaporation(	outputFolder);
				Hashtable<LocalDateTime, double[]> sm	= getSoilMoisture(	outputFolder);
				currentStreamflow.put(particleID, q.get(end));
				LocalDateTime current 	= start.plus(modelTimeStep);
				while (!current.isAfter(end))
				{
					line				+= q.get(current) + "\t";
					current				= current.plus(modelTimeStep);
				}
				current					= start.plus(modelTimeStep);
				double weight			= particle.getWeight();
				while (!current.isAfter(forecastEnd))
				{
					Double qi			= q.get(current);
					Double evi			= ev.get(current);
					double[] smc		= sm.get(current);
					if (qi != null)		forecastStreamflow.get(current).addSample(qi, weight);
					if (evi != null)	forecastEvaporation.get(current).addSample(evi, weight);
					if (smc != null)
					{
						forecastSoilMoistureL1.get(current).addSample(smc[0], weight);
						forecastSoilMoistureL2.get(current).addSample(smc[1], weight);
						forecastSoilMoistureL3.get(current).addSample(smc[2], weight);
					}
					current				= current.plus(modelTimeStep);
				}
				
				try
				{
					State finalState		= new State(forecastEnd, runFolder + "/output",
															input.rows, input.cols);
					ModelConfiguration conf	= new ModelConfiguration(soils, vegetations, network,
																		finalState);
					ArrayList<Double> vals	= configurator.toArray(conf);
					forecastEndStates.add(new Sample(weight, vals));
					
					System.out.println("Completed forecast for " + particleID);
				} catch (Exception e)
				{
					System.out.println("Completed partial forecast for "
										+ particleID + ": " + e.getMessage());
				}
				
				forecastQ.put(particleID, line);
			} catch (Exception e)
			{
				e.printStackTrace();
			} finally
			{
				// Delete folder
				if (removeForecastFiles)
					try
					{
						FileUtils.deleteDirectory(new File(runFolder));
					} catch (Exception e)
					{
						// Do nothing
					}
			}
		}
	}

	public Hashtable<LocalDateTime, KernelDensity> getForecastStreamflow()
	{
		return forecastStreamflow;
	}
	
	public Hashtable<LocalDateTime, KernelDensity> getForecastEvaporation()
	{
		return forecastEvaporation;
	}
	
	public Hashtable<LocalDateTime, KernelDensity> getForecastSoilMoistureL1()
	{
		return forecastSoilMoistureL1;
	}
	
	public Hashtable<LocalDateTime, KernelDensity> getForecastSoilMoistureL2()
	{
		return forecastSoilMoistureL2;
	}
	
	public Hashtable<LocalDateTime, KernelDensity> getForecastSoilMoistureL3()
	{
		return forecastSoilMoistureL3;
	}
	
}
