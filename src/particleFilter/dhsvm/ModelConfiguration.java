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

import java.util.ArrayList;

import dhsvm.Soil;
import dhsvm.Vegetation;
import dhsvm.grid.State;
import dhsvm.stream.StreamNetwork;

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
public class ModelConfiguration
{

	public ArrayList<Soil>			soils;
	public ArrayList<Vegetation>	vegetations;
	public StreamNetwork			network;
	public State					initialState;
	
	public ModelConfiguration(ArrayList<Soil> soils, ArrayList<Vegetation> vegetations,
								StreamNetwork network, State initialState)
	{
		this.soils			= soils;
		this.vegetations	= vegetations;
		this.network		= network;
		this.initialState	= initialState;
	}
	
}
