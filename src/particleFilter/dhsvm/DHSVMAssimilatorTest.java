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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
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

import dhsvm.MetStation;
import dhsvm.Soil;
import dhsvm.Vegetation;
import dhsvm.grid.Input;
import dhsvm.stream.StreamNetwork;
import maestro_mo.ContVar;
import particleFilter.Particle;
import probDist.ContProbDist;
import probDist.KernelDensity;
import probDist.multiVar.MultiVarKernelDensity;
import probDist.multiVar.tools.Sample;
import utilities.Utilities;
import utilities.geom.Point2D;
import utilities.stat.ContStats;

/**
 * Citation: Hernández, F. and Liang, X.: Hybridizing Bayesian and variational data assimilation
 * for high-resolution hydrologic forecasting, Hydrol. Earth Syst. Sci., 22, 5759-5779,
 * https://doi.org/10.5194/hess-22-5759-2018, 2018.
 * 
 * @author Felipe Hernández (developer)
 * @author Xu Liang (advisor)
 */
public class DHSVMAssimilatorTest
{

	// --------------------------------------------------------------------------------------------
	// Constants
	// --------------------------------------------------------------------------------------------
	
	public final static String STATES_FOLDER		= "Base states";
	public final static String MODELS_PREP_FOLDER	= "Preparation";
	public final static String DA_FOLDER 			= "Data assimilation";
	
	public final static String STATE_FILE_DT_FORMAT	= "yyyyMMdd HH-mm";
	
	public final static String OUT_FILE_STATS		= "Stats.txt";
	public final static String OUT_FILE_Q_VALS		= "Q.txt";
	public final static String OUT_FILE_EV_VALS		= "Ev.txt";
	public final static String OUT_FILE_SM_L1_VALS	= "SM1.txt";
	public final static String OUT_FILE_SM_L2_VALS	= "SM2.txt";
	public final static String OUT_FILE_SM_L3_VALS	= "SM3.txt";
	public final static String OUT_FILE_WEIGHTS		= "W.txt";
	public final static String OUT_FILE_PERF		= "Performance.txt";
	
	public final static int MAX_STATE_FILES			= 50;
	
	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	private DHSVMAssimilator					assimilator;
	
	private Duration							modelTimeStep;
	private Duration							daTimeStep;
	
	private Hashtable<LocalDateTime, Double>	obsQ;
	
	private boolean								removeDAFiles;
	private boolean								removeDAModelFiles;
	private boolean								removeForecastFiles;
	private String								modelsFolder;
	private int									threadCount;
	private long								timeLimit;
	
	private String								outputFolder;
	private DateTimeFormatter					formatter;
	
	private TreeSet<LocalDateTime>				stateList;					
	private String								statesFolder;
	private String								stateFileHeader;
	
	private ArrayList<Duration>					leadTimes;
	private ArrayList<LocalDateTime>			progress;
	
	private HashSet<String>						toDelete;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	public DHSVMAssimilatorTest(ModelConfigurator configurator, boolean defaultParameters, 
			Duration modelTimeStep, Input input, int layers, ArrayList<Soil> soils,
			ArrayList<Vegetation> vegetations, StreamNetwork network,
			ArrayList<MetStation> stations, String optionsFile, String areaFile,
			String constantsFile, String dhsvmExec, boolean removeDAFiles,
			boolean removeDAModelFiles, boolean removeForecastFiles) throws IOException
	{
		this.modelTimeStep			= modelTimeStep;
		
		this.removeDAFiles			= removeDAFiles;
		this.removeDAModelFiles		= removeDAModelFiles;
		this.removeForecastFiles	= removeForecastFiles;
		
		ArrayList<String> metFiles	= new ArrayList<>(stations.size());
		for (MetStation station : stations)
			metFiles.add(station.dataFile);
		
		assimilator					= new DHSVMAssimilator(configurator, input, soils, vegetations,
				network, stations, defaultParameters, optionsFile, areaFile, constantsFile,
				metFiles, modelTimeStep, dhsvmExec);
		
		formatter					= DateTimeFormatter.ofPattern(STATE_FILE_DT_FORMAT);
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------
	
	public void testAssimilator(String outputFolder,
			String modelsFolder, LocalDateTime forecastStart, LocalDateTime forecastEnd,
			Duration daTimeStep, ArrayList<Duration> leadTimes, double leadTimeRange,
			ArrayList<Particle> baseState, LocalDateTime baseStateTime,
			ArrayList<ContVar> variables, int ensembleSize, int dimLimit,
			double corrThreshold, double scaling, int threadCount, long timeLimit,
			int maxDARetries, Hashtable<LocalDateTime, Double> obsQ, boolean absoluteError,
			boolean resample, boolean perturb, boolean fClassKernels, double obsError)
					throws IOException
	{
		this.outputFolder					= outputFolder;
		this.modelsFolder					= modelsFolder;
		this.daTimeStep						= daTimeStep;
		this.leadTimes						= leadTimes;
		this.threadCount					= threadCount;
		this.timeLimit						= timeLimit;
		this.obsQ							= obsQ;
		
		toDelete							= new HashSet<>();
		
		// Prepare list with current progress times
		int leadTimeCount					= leadTimes.size();
		progress							= new ArrayList<>(leadTimeCount);
		for (int t = 0; t < leadTimeCount; t++)
			progress.add(forecastStart);
		
		// Initialize report files
		Collections.sort(leadTimes);
		for (int t = 0; t < leadTimes.size(); t++)
		{
			Duration leadTime				= leadTimes.get(t);
			
			// Determine file names
			String ltStr					= leadTime.toString();
			String ltFolder					= outputFolder + "/Lead time = " + ltStr;
			String mFile					= ltFolder + "/" + OUT_FILE_STATS;
			String qFile					= ltFolder + "/" + OUT_FILE_Q_VALS;
			String eFile					= ltFolder + "/" + OUT_FILE_EV_VALS;
			String s1File					= ltFolder + "/" + OUT_FILE_SM_L1_VALS;
			String s2File					= ltFolder + "/" + OUT_FILE_SM_L2_VALS;
			String s3File					= ltFolder + "/" + OUT_FILE_SM_L3_VALS;
			String wFile					= ltFolder + "/" + OUT_FILE_WEIGHTS;
			
			// Initialize output folder and determine if progress has already been made
			Path path						= FileSystems.getDefault().getPath(ltFolder);
			if (Files.exists(path))
			{
				Scanner scanner				= new Scanner(new FileInputStream(new File(mFile)));
				String line					= null;
				while (scanner.hasNextLine())
					line					= scanner.nextLine();
				scanner.close();
				if (line != null)
				{
					String[] tokens			= line.split("\t");
					LocalDateTime dateTime	= LocalDateTime.parse(tokens[0], formatter);
					dateTime				= dateTime.plus(modelTimeStep);
					if (dateTime.isAfter(progress.get(t)))
							progress.set(t, dateTime);
				}
				continue;
			}
			Files.createDirectory(path);
			
			// Create output files
			PrintWriter out	= new PrintWriter(new BufferedWriter(new FileWriter(mFile, false)));
			out.println("Date-time\tQ_Mean\tQ_stDev\tEv_Mean\tEv_stDev\tSM1_Mean\tSM1_stDev\t"
							+ "SM2_Mean\tSM2_stDev\tSM3_Mean\tSM3_stDev");
			out.close();
			out				= new PrintWriter(new BufferedWriter(new FileWriter(qFile, false)));
			out.println("Date-time\tStreamflow values");
			out.close();
			out				= new PrintWriter(new BufferedWriter(new FileWriter(eFile, false)));
			out.println("Date-time\tEvaporation values");
			out.close();
			out				= new PrintWriter(new BufferedWriter(new FileWriter(s1File, false)));
			out.println("Date-time\tSoil moisture values (layer 1)");
			out.close();
			out				= new PrintWriter(new BufferedWriter(new FileWriter(s2File, false)));
			out.println("Date-time\tSoil moisture values (layer 2)");
			out.close();
			out				= new PrintWriter(new BufferedWriter(new FileWriter(s3File, false)));
			out.println("Date-time\tSoil moisture values (layer 3)");
			out.close();
			out				= new PrintWriter(new BufferedWriter(new FileWriter(wFile, false)));
			out.println("Date-time\tParticle weights");
			out.close();
		}
		
		// Initialize state folders
		stateList							= new TreeSet<>();
		statesFolder						= modelsFolder + "/" + STATES_FOLDER;
		String prepFolder					= modelsFolder + "/" + MODELS_PREP_FOLDER;
		StringBuilder builder				= new StringBuilder("Id\tWeight");
		for (ContVar var : variables)
		{
			builder.append("\t");
			builder.append(var.getName());
		}
		stateFileHeader						= builder.toString();
		Path path							= FileSystems.getDefault().getPath(prepFolder);
		if (!Files.exists(path))
			Files.createDirectory(path);
		else
		{
			// Find existing files
			File statesDir					= new File(statesFolder);
			for (File file : statesDir.listFiles())
			{
				String name					= file.getName();
				name						= name.substring(0, name.length() - 4);
				try
				{
					LocalDateTime dateTime	= LocalDateTime.parse(name, formatter);
					stateList.add(dateTime);
				} catch (Exception e)
				{
					e.printStackTrace();
				}
			}
		}
		path								= FileSystems.getDefault().getPath(statesFolder);
		if (!Files.exists(path))			Files.createDirectory(path);
		writeStateToFile(baseStateTime, baseState);
		
		// Perform test loop
		path					= FileSystems.getDefault().getPath(outputFolder + "/" + DA_FOLDER);
		if (!Files.exists(path))			Files.createDirectory(path);
		boolean done						= false;
		long modelTSMillis					= modelTimeStep.toMillis();
		int retries							= 0;
		LocalDateTime oldDAStart			= LocalDateTime.MIN;
		while (!done)
		{
			// Determine assimilation start time
			LocalDateTime daStart			= LocalDateTime.MAX;
			for (int t = 0; t < leadTimeCount; t++)
			{
				LocalDateTime current		= progress.get(t);
				if (current.isAfter(forecastEnd))
					continue;
				double millis				= leadTimes.get(t).toMillis()*(1.0 - leadTimeRange);
				int modelTSMult				= (int)Math.ceil((millis)/modelTSMillis);
				Duration lower				= modelTimeStep.multipliedBy(modelTSMult);
				LocalDateTime start	= current.minus(lower).minus(modelTimeStep).minus(daTimeStep);
				daStart						= start.isBefore(daStart) ? start : daStart;
			}
			
			// Verify if enough retries have been attempted
			if (!daStart.isAfter(oldDAStart))
			{
				if (retries >= maxDARetries)
				{
					storeNullResults(daStart);
					daStart					= daStart.plus(modelTimeStep);
					retries					= 0;
				}
				else
					retries++;
			}
			else				
				retries						= 0;
			oldDAStart						= daStart;
			LocalDateTime daEnd				= daStart.plus(daTimeStep);
			
			// Prepare assimilation
			ArrayList<Particle> initStates	= getState(daStart);
			
			// TODO Remove flag
			System.out.println("Samples in loaded state: " + initStates.size());
			
			String strDT					= formatter.format(daEnd.plus(modelTimeStep));			
			String folder					= outputFolder + "/" + DA_FOLDER + "/" + strDT;
			String modelsiFolder			= modelsFolder + "/" + strDT;
			path							= FileSystems.getDefault().getPath(folder);
			if (retries > 0)
			{
				folder						+= "_" + retries;
				path						= FileSystems.getDefault().getPath(folder);
			}
			Files.createDirectory(path);
			Files.createDirectory(FileSystems.getDefault().getPath(modelsiFolder));
			
			// Assimilate and store target state
			MultiVarKernelDensity dist		= new MultiVarKernelDensity(true);
			for (Particle particle : initStates)
				dist.addSample(new Sample(particle.getWeight(), particle.getState()));
			if (fClassKernels)
				dist.computeGaussianBW(scaling, dimLimit, corrThreshold);
			else
				dist.computeGaussianDiagBW(scaling);
			ArrayList<Particle> targetStates = assimilator.assimilate(dist, folder,
					modelsiFolder, daStart, daEnd, daTimeStep, obsQ, absoluteError, obsError,
					ensembleSize, resample, perturb, fClassKernels);
					
			writeStateToFile(daEnd, targetStates);
			
			// Determine forecast times
			ArrayList<ArrayList<LocalDateTime>> forTimes = new ArrayList<>(leadTimeCount);
			HashSet<LocalDateTime> allTimes	= new HashSet<>();
			for (int t = 0; t < leadTimeCount; t++)
			{
				LocalDateTime dateTime		= daEnd.plus(modelTimeStep);
				ArrayList<LocalDateTime> times = new ArrayList<>();
				forTimes.add(times);
				LocalDateTime prog			= progress.get(t);
				Duration leadTime			= leadTimes.get(t);
				double rangeMillis			= leadTime.toMillis()*leadTimeRange;
				int modelTSMult				= (int)Math.ceil((rangeMillis)/modelTSMillis);
				Duration range				= modelTimeStep.multipliedBy(modelTSMult);
				Duration lower				= leadTime.minus(	range);
				Duration upper				= leadTime.plus(	range);
				
				Duration remaining			= Duration.between(dateTime.plus(lower), forecastEnd);
				if (dateTime.plus(leadTime).isBefore(prog)
						&& upper.minus(lower).compareTo(remaining) < 0)
					continue;
				
				Duration currLeadTime		= Duration.between(daEnd, dateTime);
				while (currLeadTime.compareTo(upper) <= 0 && !dateTime.isAfter(forecastEnd))
				{
					if (currLeadTime.compareTo(lower) >= 0 && !dateTime.isBefore(prog))
					{
						times.add(dateTime);
						allTimes.add(dateTime);
					}
					dateTime				= dateTime.plus(modelTimeStep);
					currLeadTime			= Duration.between(daEnd, dateTime);
				}
			}
			
			// Determine forecast duration
			ArrayList<LocalDateTime> statesToSave = new ArrayList<>(allTimes);
			Collections.sort(statesToSave);
			LocalDateTime forecastiEnd		= statesToSave.get(statesToSave.size() - 1);
			Duration forecastDuration		= Duration.between(daEnd, forecastiEnd);
			System.out.println("Forecast duration: " + forecastDuration.toString());
			for (int l = 0; l < leadTimes.size(); l++)
			{
				ArrayList<LocalDateTime> times = forTimes.get(l);
				String line					= leadTimes.get(l).toString() + ": " + times.size() + " hours";
				if (times.size() > 0)
					line					+= " starting at " + times.get(0).toString();
				System.out.println(line);
			}
			
			// Perform forecast
			assimilator.prepareForecast(targetStates, daEnd, modelsiFolder, obsQ);
			statesToSave.clear();
			statesToSave.add(forecastiEnd);
			assimilator.forecast(forecastiEnd, modelsiFolder, threadCount, null, statesToSave,
									removeForecastFiles, timeLimit);
			
			// Obtain and store results
			Hashtable<LocalDateTime, KernelDensity> q	= assimilator.getForecastStreamflow();
			Hashtable<LocalDateTime, KernelDensity> ev	= assimilator.getForecastEvaporation();
			Hashtable<LocalDateTime, KernelDensity> sm1	= assimilator.getForecastSoilMoistureL1();
			Hashtable<LocalDateTime, KernelDensity> sm2	= assimilator.getForecastSoilMoistureL2();
			Hashtable<LocalDateTime, KernelDensity> sm3	= assimilator.getForecastSoilMoistureL3();
			storeResults(forTimes, q, ev, sm1, sm2, sm3);
			
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
			
			// Keep running?
			boolean keepRunning				= false;
			for (LocalDateTime prog : progress)
				if (!prog.isAfter(forecastEnd))
					
				{
					keepRunning				= true;
					break;
				}
			done							= !keepRunning;
		}

		// Compute performance
		computeForecastPerformance(leadTimes, forecastStart, forecastEnd);
		
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
	
	private ArrayList<Particle> getState(LocalDateTime dateTime) throws IOException
	{
		String strDT				= formatter.format(dateTime);
		
		// Try loading initial state
		DateTimeFormatter formatter	= DateTimeFormatter.ofPattern(STATE_FILE_DT_FORMAT);
		String fileName				= statesFolder + "/" + formatter.format(dateTime) + ".txt";
		if (Files.exists(FileSystems.getDefault().getPath(fileName)))
		{
			ArrayList<Particle> state = readStateFromFile(dateTime);
			System.out.println("\nLoaded initial state (" + strDT + ")...");
			return state;
		}
		
		System.out.println("\nPreparing initial state for assimilation (" + strDT + ")...");
		LocalDateTime baseDT		= stateList.lower(dateTime);
		ArrayList<Particle> baseState = readStateFromFile(baseDT);
		String folder				= modelsFolder + "/" + MODELS_PREP_FOLDER + "/" + strDT;
		Files.createDirectory(FileSystems.getDefault().getPath(folder));
		assimilator.prepareForecast(baseState, baseDT, folder, obsQ);
		ArrayList<Particle> state	= assimilator.forecast(dateTime, folder, threadCount, null,
													null, removeForecastFiles, timeLimit);
		try
		{
			FileUtils.deleteDirectory(new File(folder));
		} catch (Exception e) {}
		return state;
	}
	
	private void writeStateToFile(LocalDateTime dateTime, ArrayList<Particle> particles)
					throws IOException
	{
		String fileName					= statesFolder + "/" + formatter.format(dateTime) + ".txt";
		PrintWriter out	= new PrintWriter(new BufferedWriter(new FileWriter(fileName, false)));
		out.println(stateFileHeader);
		
		// TODO Remove flag
		System.out.println("Particles in state to save: " + particles.size());
		
		int index						= 1;
		for (Particle particle : particles)
		{
			StringBuilder builder		= new StringBuilder(index + "\t" + particle.getWeight());
			for (Double value : particle.getState())
			{
				builder.append("\t");
				builder.append(value);
			}
			out.println(builder.toString());
			index++;
		}
		out.close();
		stateList.add(dateTime);
		
		// Limit number of state files
		while (stateList.size() > MAX_STATE_FILES)
		{
			try
			{
				ArrayList<LocalDateTime> selector = new ArrayList<>(stateList);
				Collections.shuffle(selector);
				LocalDateTime selected	= selector.get(0);
				fileName				= statesFolder + "/" + formatter.format(selected) + ".txt";
				stateList.remove(selected);
				Files.delete(FileSystems.getDefault().getPath(fileName));
			} catch (Exception e)
			{
				e.printStackTrace();
			}
		}
	}
	
	private ArrayList<Particle> readStateFromFile(LocalDateTime dateTime)
													throws FileNotFoundException
	{
		DateTimeFormatter formatter	= DateTimeFormatter.ofPattern(STATE_FILE_DT_FORMAT);
		String fileName				= statesFolder + "/" + formatter.format(dateTime) + ".txt";
		Scanner scanner				= null;

		// Read samples from file
		scanner						= new Scanner(new FileInputStream(new File(fileName)));
		scanner.nextLine();
		ArrayList<Particle> particles = new ArrayList<>();
		int index					= 1;
		while (scanner.hasNextLine())
		{
			String[] tokens			= scanner.nextLine().split("\t");
			double weight			= Double.valueOf(tokens[1]);
			int valueCount			= tokens.length - 2;
			ArrayList<Double> values = new ArrayList<>(valueCount);
			for (int v = 2; v < valueCount + 2; v++)
				values.add(Double.valueOf(tokens[v]));
			particles.add(new Particle(DHSVMAssimilator.PARTICLE_ID_PREFIX + " " + index,
										values, weight));
			index++;
		}
		scanner.close();
		
		return particles;
	}

	private void storeResults(	ArrayList<ArrayList<LocalDateTime>> forcastTimes,
								Hashtable<LocalDateTime, KernelDensity> q,
								Hashtable<LocalDateTime, KernelDensity> ev,
								Hashtable<LocalDateTime, KernelDensity> sm1,
								Hashtable<LocalDateTime, KernelDensity> sm2,
								Hashtable<LocalDateTime, KernelDensity> sm3) throws IOException
	{
		for (int l = 0; l < leadTimes.size(); l++)
		{
			Duration leadTime				= leadTimes.get(l);
			
			// Get strings
			String ltStr					= leadTime.toString();
			String ltFolder					= outputFolder + "/Lead time = " + ltStr;

			String mFile					= ltFolder + "/" + OUT_FILE_STATS;
			String qFile					= ltFolder + "/" + OUT_FILE_Q_VALS;
			String eFile					= ltFolder + "/" + OUT_FILE_EV_VALS;
			String s1File					= ltFolder + "/" + OUT_FILE_SM_L1_VALS;
			String s2File					= ltFolder + "/" + OUT_FILE_SM_L2_VALS;
			String s3File					= ltFolder + "/" + OUT_FILE_SM_L3_VALS;
			String wFile					= ltFolder + "/" + OUT_FILE_WEIGHTS;
			
			ArrayList<LocalDateTime> times	= forcastTimes.get(l);
			if (times.size() == 0)
				continue;
			
			LocalDateTime next				= progress.get(l);
			for (LocalDateTime forecastTime : times)
			{
				String dtStr				= formatter.format(forecastTime);
				
				KernelDensity streamflow	= q.get(	forecastTime);
				KernelDensity evaporation	= ev.get(	forecastTime);
				KernelDensity soilMoisture1	= sm1.get(	forecastTime);
				KernelDensity soilMoisture2	= sm2.get(	forecastTime);
				KernelDensity soilMoisture3	= sm3.get(	forecastTime);
				
				try
				{
					if (streamflow.getSamples().size() == 0	|| evaporation.getSamples().size() == 0
															|| soilMoisture1.getSamples().size() == 0)
						break;
				} catch (Exception e)
				{
					break;
				}
				
				// Write statistics
				Path mPath					= Paths.get(mFile);
				String line					= System.lineSeparator() + dtStr + "\t"
								+ streamflow.getMean() + "\t" + streamflow.getStDev() + "\t"
								+ evaporation.getMean() + "\t" + evaporation.getStDev() + "\t"
								+ soilMoisture1.getMean() + "\t" + soilMoisture1.getStDev() + "\t"
								+ soilMoisture2.getMean() + "\t" + soilMoisture2.getStDev() + "\t"
								+ soilMoisture3.getMean() + "\t" + soilMoisture3.getStDev();
				BufferedWriter writer = Files.newBufferedWriter(mPath, StandardOpenOption.APPEND);
				writer.write(line);
				writer.close();
				
				// Write streamflow and weights
				PrintWriter outQ			= new PrintWriter(new BufferedWriter(
															new FileWriter(qFile, true)));
				PrintWriter outW			= new PrintWriter(new BufferedWriter(
															new FileWriter(wFile, true)));
				String lineQ				= dtStr;
				String lineW				= dtStr;
				ArrayList<Point2D> samples	= streamflow.getSamplesWeights();
				for (Point2D sample : samples)
					if (sample != null)
					{
						lineQ				+= "\t" + sample.x;
						lineW				+= "\t" + sample.y;
					}
				outQ.println(lineQ);
				outW.println(lineW);
				outQ.close();
				outW.close();
				
				// Write evaporation values
				PrintWriter outE			= new PrintWriter(new BufferedWriter(
															new FileWriter(eFile, true)));
				String lineE				= dtStr;
				samples						= evaporation.getSamplesWeights();
				for (Point2D sample : samples)
					if (sample != null)
						lineE				+= "\t" + sample.x;
				outE.println(lineE);
				outE.close();
				
				// Write soil moisture values (layer 1)
				PrintWriter outS1			= new PrintWriter(new BufferedWriter(
															new FileWriter(s1File, true)));
				String lineS1				= dtStr;
				samples						= soilMoisture1.getSamplesWeights();
				for (Point2D sample : samples)
					if (sample != null)
						lineS1				+= "\t" + sample.x;
				outS1.println(lineS1);
				outS1.close();
				
				// Write soil moisture values (layer 2)
				PrintWriter outS2			= new PrintWriter(new BufferedWriter(
															new FileWriter(s2File, true)));
				String lineS2				= dtStr;
				samples						= soilMoisture2.getSamplesWeights();
				for (Point2D sample : samples)
					if (sample != null)
						lineS2				+= "\t" + sample.x;
				outS2.println(lineS2);
				outS2.close();
				
				// Write soil moisture values (layer 3)
				PrintWriter outS3			= new PrintWriter(new BufferedWriter(
															new FileWriter(s3File, true)));
				String lineS3				= dtStr;
				samples						= soilMoisture3.getSamplesWeights();
				for (Point2D sample : samples)
					if (sample != null)
						lineS3				+= "\t" + sample.x;
				outS3.println(lineS3);
				outS3.close();
				
				next						= forecastTime.plus(modelTimeStep);
			}
			
			// Update current progress
			progress.set(l, next);
		}
	}

	private void storeNullResults(LocalDateTime daStart) throws IOException
	{
		for (int l = 0; l < leadTimes.size(); l++)
		{
			Duration leadTime		= leadTimes.get(l);
			LocalDateTime prog		= progress.get(l);
			if (prog.isAfter(daStart.plus(daTimeStep).plus(leadTime)))
				continue;
			
			String ltStr			= leadTime.toString();
			String ltFolder			= outputFolder + "/Lead time = " + ltStr;
			String mFile			= ltFolder + "/" + OUT_FILE_STATS;
			String qFile			= ltFolder + "/" + OUT_FILE_Q_VALS;
			String eFile			= ltFolder + "/" + OUT_FILE_EV_VALS;
			String s1File			= ltFolder + "/" + OUT_FILE_SM_L1_VALS;
			String s2File			= ltFolder + "/" + OUT_FILE_SM_L2_VALS;
			String s3File			= ltFolder + "/" + OUT_FILE_SM_L3_VALS;
			String wFile			= ltFolder + "/" + OUT_FILE_WEIGHTS;
			
			String dtStr			= formatter.format(prog);
			
			Path mPath				= Paths.get(mFile);
			BufferedWriter writer	= Files.newBufferedWriter(mPath, StandardOpenOption.APPEND);
			writer.write(System.lineSeparator() + dtStr);
			writer.close();
			
			PrintWriter outQ		= new PrintWriter(new BufferedWriter(
																new FileWriter(qFile, true)));
			PrintWriter outW		= new PrintWriter(new BufferedWriter(
																new FileWriter(wFile, true)));
			PrintWriter outE		= new PrintWriter(new BufferedWriter(
																new FileWriter(eFile, true)));
			PrintWriter outS1		= new PrintWriter(new BufferedWriter(
																new FileWriter(s1File, true)));
			PrintWriter outS2		= new PrintWriter(new BufferedWriter(
																new FileWriter(s2File, true)));
			PrintWriter outS3		= new PrintWriter(new BufferedWriter(
																new FileWriter(s3File, true)));
			
			outQ.println(	dtStr);
			outW.println(	dtStr);
			outE.println(	dtStr);
			outS1.println(	dtStr);
			outS2.println(	dtStr);
			outS3.println(	dtStr);
			
			outQ.close();
			outW.close();
			outE.close();
			outS1.close();
			outS2.close();
			outS3.close();
			
			progress.set(l, prog.plus(modelTimeStep));
		}
	}
	
	private void computeForecastPerformance(ArrayList<Duration> leadTimes,
			LocalDateTime forecastStart, LocalDateTime forecastEnd) throws IOException
	{
		System.out.println("\nForecast performance:");
		System.out.println("Lead time\tNSE_l2\tNSE_l1\tMARE\tCRPS\tDensity");
		for (Duration leadTime : leadTimes)
		{
			String ltStr					= leadTime.toString();
			String ltFolder					= outputFolder + "/Lead time = " + ltStr;
			String valuesFile				= ltFolder + "/" + OUT_FILE_Q_VALS;
			String weightsFile				= ltFolder + "/" + OUT_FILE_WEIGHTS;
			ArrayList<ContProbDist> series	= new ArrayList<>();
			Scanner scannerV				= new Scanner(new FileInputStream(new File(valuesFile)));
			Scanner scannerW				= new Scanner(new FileInputStream(new File(weightsFile)));
			scannerV.nextLine();	// Header
			scannerW.nextLine();	// Header
			
			// Load distribution time series from files
			ArrayList<Double> meanQ			= new ArrayList<>();
			ArrayList<Double> obs			= new ArrayList<>();
			while (scannerV.hasNextLine())
			{
				String[] tokensV			= scannerV.nextLine().split("\t");
				String[] tokensW			= scannerW.nextLine().split("\t");
				
				LocalDateTime dateTime		= LocalDateTime.parse(tokensV[0], formatter);
				Double obsi					= obsQ.get(dateTime);
				if (obsi != null 	&& dateTime.isAfter(	forecastStart	)
									&& !dateTime.isAfter(	forecastEnd		))
				{
					KernelDensity dist			= new KernelDensity();
					dist.setWeighted(true);
					for (int v = 1; v < Math.min(tokensV.length, tokensW.length); v++)
					{
						double value			= Double.valueOf(tokensV[v]);
						double weight			= Double.valueOf(tokensW[v]);
						dist.addSample(value, weight);
					}
					if (dist.getSamples().size() == 0)
					{
						scannerV.close();
						scannerW.close();
						throw new RuntimeException("No forecasted values at time " + tokensV[0]);
					}
					dist.computeGaussianBandwidth();
					series.add(dist);
					meanQ.add(dist.getMean());
					obs.add(obsi);
				}
			}
			scannerV.close();
			scannerW.close();
			
			// Compute density and CRPS
			ContStats density		= new ContStats(false);
			ContStats crps			= new ContStats(false);
			for (int t = 0; t < obs.size(); t++)
			{
				double obsi			= obs.get(t);
				ContProbDist dist	= series.get(t);
				density.addValue(dist.getpdf(obsi));
				KernelDensity dist2	= (KernelDensity)dist;
				crps.addValue(dist2.computeEnsembleCRPS(obsi));
			}

			// Compute performance
			double[] obsArr			= Utilities.toArray(obs);
			double[] modArr			= Utilities.toArray(meanQ);
			double nse_l2			= Utilities.computeNashSutcliffe(obsArr, modArr);
			double nse_l1			= Utilities.computeNashSutcliffe(obsArr, modArr, 1.0);
			double mare				= Utilities.computeMeanAbsRelativeError(obsArr, modArr);
			double meanDensity		= density.getMean();
			double meanCRPS			= crps.getMean();

			// Store performance file
			String performanceFile	= ltFolder + "/" + OUT_FILE_PERF;
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(performanceFile)));
			out.println("Values\t"			+ obsArr.length	);
			out.println("NSE_l2\t"			+ nse_l2		);
			out.println("NSE_l1\t"			+ nse_l1		);
			out.println("MARE\t"			+ mare			);
			out.println("Mean CRPS\t"		+ meanCRPS		);
			out.println("Mean density\t"	+ meanDensity	);
			out.close();
			
			// Print to console
			System.out.println(ltStr + "\t" + nse_l2 + "\t" + nse_l1 + "\t" + mare + "\t"
								+ meanCRPS + "\t" + meanDensity);
		}
	}
	
}
