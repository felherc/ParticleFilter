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

package particleFilter;

import java.util.ArrayList;

/**
 * Citation: Hernández, F., & Liang, X. (2018). "Hybridizing Bayesian and variational data
 * assimilation for high-resolution hydrologic forecasting". Hydrol. Earth Syst. Sci.
 * 
 * @author Felipe Hernández (developer)
 * @author Xu Liang (advisor)
 */
public interface Model
{

	// --------------------------------------------------------------------------------------------
	// Method
	// --------------------------------------------------------------------------------------------
	
	public ModelResult runModel(int index, ArrayList<Double> sourceState);
	
}
