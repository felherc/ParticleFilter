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

import java.time.LocalDateTime;
import java.util.ArrayList;

import particleFilter.Particle;
import probDist.multiVar.tools.ContMultiSample;

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
public interface ModelConfigurator
{

	ModelConfiguration configure(ArrayList<Double> valueArray, LocalDateTime dateTime,
									boolean defaultParameters);
	
	ArrayList<Double> toArray(ModelConfiguration configuration);
	
	ArrayList<Particle> createDistribution(ArrayList<ContMultiSample> samples);
	
}