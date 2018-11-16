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
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.time.Duration;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;

import dhsvm.MetStation;
import dhsvm.Soil;
import dhsvm.Vegetation;

/**
 * Contains the information of a DHSVM model except its temporal span and its initial state
 * 
 * Citation: Hernández, F. and Liang, X.: Hybridizing Bayesian and variational data assimilation
 * for high-resolution hydrologic forecasting, Hydrol. Earth Syst. Sci., 22, 5759-5779,
 * https://doi.org/10.5194/hess-22-5759-2018, 2018.
 * 
 * @author Felipe Hernández (developer)
 * @author Xu Liang (advisor)
 */
public class DHSVMModel
{

	// --------------------------------------------------------------------------------------------
	// Constants
	// --------------------------------------------------------------------------------------------
	
	public final static String CONF_FILE_TIME_FORMAT	= "MM/dd/yyyy-HH:mm:ss";
	
	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	// Configuration-file section lines
	public String[]					optionsSection;
	public String[]					areaSection;
	public String[]					constantsSection;
	
	// File routes of input files
	public String					demFile;
	public String					maskFile;
	public String					flowDirFile;
	public String					soilTypeFile;
	public String					soilDepthFile;
	public String					vegTypeFile;
	public String					streamClassFile;
	public String					streamNetworkFile;
	public String					streamMapFile;
	public String					surfaceRoutingFile;
	
	// Additional model elements
	public boolean					d8FlowDirection;
	public ArrayList<Soil>			soils;
	public ArrayList<Vegetation>	vegetations;
	public ArrayList<MetStation>	stations;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	public DHSVMModel()
	{
		// Do nothing
	}

	public DHSVMModel(String[] optionsSection, String[] areaSection, String[] constantsSection, 
			String demFile,	String maskFile, String flowDirFile, String soilTypeFile, 
			String soilDepthFile, String vegTypeFile, String streamClassFile, 
			String streamNetworkFile, String streamMapFile, String surfaceRoutingFile,
			boolean d8FlowDirection, ArrayList<Soil> soils, ArrayList<Vegetation> vegetations, 
			ArrayList<MetStation> stations)
	{
		this.optionsSection		= optionsSection;
		this.areaSection		= areaSection;
		this.constantsSection	= constantsSection;
		this.demFile			= demFile;
		this.maskFile			= maskFile;
		this.flowDirFile		= flowDirFile;
		this.soilTypeFile		= soilTypeFile;
		this.soilDepthFile		= soilDepthFile;
		this.vegTypeFile		= vegTypeFile;
		this.streamClassFile	= streamClassFile;
		this.streamNetworkFile	= streamNetworkFile;
		this.streamMapFile		= streamMapFile;
		this.surfaceRoutingFile	= surfaceRoutingFile;
		this.d8FlowDirection	= d8FlowDirection;
		this.soils				= soils;
		this.vegetations		= vegetations;
		this.stations			= stations;
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------
	
	public void writeConfigFile(String fileRoute, LocalDateTime start, LocalDateTime end, 
					Duration timeStep, String stateFolder, String outputFolder, 
					boolean storeFinalState, ArrayList<LocalDateTime> extraStatesToSave)
							throws IOException
	{
		String configFile						= fileRoute;
		try
		{
			Files.delete(FileSystems.getDefault().getPath(configFile));
		} catch (IOException e) {}
		
		DateTimeFormatter formatter		= DateTimeFormatter.ofPattern(CONF_FILE_TIME_FORMAT);
		
		PrintWriter out	= new PrintWriter(new BufferedWriter(new FileWriter(configFile, true)));
		
		// Options and area section
		for (int l = 0; l < optionsSection.length; l++)
			out.println(optionsSection[l]													);
		out.println(	""																	);
		for (int l = 0; l < areaSection.length; l++)
			out.println(areaSection[l]														);
		out.println(	""																	);
		
		// Time section
		out.println(	"[TIME]"															);
		out.println(	"Time Step            = " + timeStep.toHours()						);
		out.println(	"Min Time Step        = 6"											);
		out.println(	"Model Start          = " + start.format(formatter)					);
		out.println(	"Model End            = " + end.format(formatter)					);
		out.println(	""																	);
		
		// Constants section
		for (int l = 0; l < constantsSection.length; l++)
			out.println(constantsSection[l]													);
		out.println(	""																	);
		
		// Terrain section
		out.println(	"[TERRAIN]"															);
		out.println(	"DEM File              = " + demFile								);
		out.println(	"Basin Mask File       = " + maskFile								);
		out.println(	"Flow Direction Method = " + (d8FlowDirection ? "D8" : "D4")		);
		out.println(	"Flow Direction File   = " + maskFile								);
		out.println(	""																	);
		
		// Routing section
		out.println(	"[ROUTING]"															);
		out.println(	"Stream Class File     = " + streamClassFile						);
		out.println(	"Stream Network File   = " + streamNetworkFile						);
		out.println(	"Stream Map File       = " + streamMapFile							);
		out.println(	""																	);
		
		// Meteorology section
		out.println(	"[METEOROLOGY]"														);
		out.println(	"Number of stations    = " + stations.size()						);
		out.println(	""																	);
		for (int s = 0; s < stations.size(); s++)
		{
			MetStation station					= stations.get(s);
			out.println("Station Name       " + (s + 1) + " = " + station.name				);
			out.println("North Coordinate   " + (s + 1) + " = " + station.northCoordinate	);
			out.println("East Coordinate    " + (s + 1) + " = " + station.eastCoordinate	);
			out.println("Elevation          " + (s + 1) + " = " + station.elevation			);
			out.println("Station File       " + (s + 1) + " = " + station.dataFile			);
			out.println(""																	);
		}
		
		// Soils section
		out.println(	"[SOILS]"															);
		out.println(	"Soil Map File         = " + soilTypeFile							);
		out.println(	"Soil Depth File       = " + soilDepthFile							);
		out.println(	"Number of Soil Types  = " + soils.size()							);
		out.println(	""																	);
		for (Soil soil : soils)
		{
			String[] soilParams					= soil.getFileConf();
			for (int s = 0; s < soilParams.length; s++)
				out.println(soilParams[s]													);
			out.println(""																	);
		}
		
		// Vegetation section
		out.println(	"[VEGETATION]"														);
		out.println(	"Vegetation Map File				= " + vegTypeFile				);
		out.println(	"IMPERVIOUS SURFACE ROUTING FILE	= " + surfaceRoutingFile		);
		out.println(	"Number of Vegetation Types			= " + vegetations.size()		);
		out.println(	""																	);
		for (Vegetation vegetation : vegetations)
		{
			String[] vegParams			= vegetation.getFileConf();
			for(int s = 0; s < vegParams.length; s++)
				out.println(vegParams[s]													);
			out.println(""																	);
		}
		
		// Output section
		out.println(	"[OUTPUT]"															);
		out.println(	"Output Directory            = " + outputFolder	+ "/"				);
		out.println(	"Initial State Directory     = " + stateFolder	+ "/"				);
		out.println(	""																	);
		int statesToSave		= (extraStatesToSave != null ? extraStatesToSave.size() : 0)
								+ (storeFinalState ? 1 : 0);
		if (statesToSave > 0)
		{
			out.println("Number of Model States		= " + statesToSave						);
			out.println(""																	);
			for (int t = 0; t < statesToSave - (storeFinalState ? 1 : 0); t++)
			{
				String timeStr			= formatter.format(extraStatesToSave.get(t));
				out.println("State Date               " + (t + 1) + "	= " + timeStr		);
			}
			if (storeFinalState)
				out.println("State Date               "
											+ statesToSave + "	= " + formatter.format(end)	);

			out.println(""																	);
		}
		else
			out.println("Number of Model States		= 0"									);
		out.println(	"Number of Map Variables     = 0"									);
		out.println(	"Number of Image Variables   = 0"									);
		out.println(	"Number of Graphics          = 0"									);
		out.println(	""																	);
		out.println(	"[End]"																);
		out.println(	""																	);
		
		out.close();
	}
	
}
