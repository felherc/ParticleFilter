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
import java.util.Collections;
import java.util.Enumeration;
import java.util.Hashtable;

import probDist.ContProbDist;
import probDist.multiVar.MultiVarKernelDensity;
import probDist.multiVar.MultiVarNormal;
import probDist.multiVar.tools.Sample;
import utilities.stat.ContSeries;

/**
 * Citation: Hernández, F. and Liang, X.: Hybridizing Bayesian and variational data assimilation
 * for high-resolution hydrologic forecasting, Hydrol. Earth Syst. Sci., 22, 5759-5779,
 * https://doi.org/10.5194/hess-22-5759-2018, 2018.
 * 
 * @author Felipe Hernández (developer)
 * @author Xu Liang (advisor)
 */
public class Assimilator
{

	// --------------------------------------------------------------------------------------------
	// Constants
	// --------------------------------------------------------------------------------------------
	
	public static final String PREFIX = "Particle";

	public static final String RSPL = " - resample ";
	
	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	private Model model;
	
	private ContProbDist observation;

	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	public Assimilator(Model model, ContProbDist observation)
	{
		this.model			= model;
		this.observation	= observation;
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------
	
	public ArrayList<Particle> assimilate(ArrayList<Particle> sourceState, int ensembleSize, 
										boolean resample, boolean perturb, boolean fClassKernels)
	{
		// Run model and compute weights
		ArrayList<Particle> particles		= new ArrayList<>();
		ArrayList<Particle> particles2		= new ArrayList<>();
		ContSeries partWeights				= new ContSeries(true);
		double weightSum					= 0.0;
		for (int p = 0; p < sourceState.size(); p++)
		{
			Particle particle 				= sourceState.get(p);
			String id						= PREFIX + " " + (p + 1);
			ModelResult result 				= model.runModel(p + 1, particle.getState());
			String error					= result.error;
			double weight					= Double.NaN;
			if (error.equals(""))
			{
				weight						= observation.getpdf(result.output);
				weightSum					+= weight;
				System.out.println(id + ": weight = " + weight);
			}
			else
			{
				weight						= 0.0;
				System.out.println(id + ": " + error);
			}
			particles.add(new Particle(id, result.state, weight));
			partWeights.addValue(p, weight);	
		}
		
		// Adjust weights if all are zero
		if (weightSum == 0.0)
		{
			partWeights						= new ContSeries(false);
			for (Particle particle : particles)
			{
				particle.setWeight(1.0);
				partWeights.addValue(1.0);
			}
		}
		
		if (!resample)
		{
			ArrayList<Integer> indices		= new ArrayList<>();
			for (int i = 0; i < particles.size(); i++)
				indices.add(i);
			Collections.shuffle(indices);
			ArrayList<Integer> indices2		= new ArrayList<>();
			for (int i = 0; i < ensembleSize ; i++)
				indices2.add(indices.get(i));
			Collections.sort(indices2);
			for (int i = 0; i < indices2.size(); i++)
				particles2.add(particles.get(indices2.get(i)));
			return particles2;
		}
		
		// Resample
		Hashtable<Integer, Integer> samples	= new Hashtable<>();
		for (int s = 0; s < ensembleSize; s++)
		{
			Integer index					= (int) partWeights.sample();
			if (samples.containsKey(index))
			{
				Integer count				= samples.get(index);
				samples.put(index, count + 1);
			}
			else
				samples.put(index, 1);
		}
		
		if (!perturb)
		{
			Enumeration<Integer> keys		= samples.keys();
			while (keys.hasMoreElements())
			{
				Integer index				= keys.nextElement();
				int count					= samples.get(index);
				Particle original			= particles.get(index);
				for (int r = 0; r < count; r++)
				{
					String id				= original.getId();
					id						+= r > 0 ? RSPL + r : "";
					particles2.add(new Particle(id, original.getState(), 1.0));
				}
			}
			return particles2;
		}
		
		// Determine sampling kernel
		MultiVarKernelDensity dist			= new MultiVarKernelDensity();
		for (Particle particle : particles)
			if (particle.getWeight() > 0.0)
				dist.addSample(new Sample(particle.getWeight(), particle.getState()));
		MultiVarNormal kernel				= null;
		int dimensions						= dist.getSamples().get(0).getValues().size();
		double[] zeroVec					= new double[dimensions];
		for (int d = 0; d < dimensions; d++)
			zeroVec[d]						= 0.0;
		if (fClassKernels)
		{
			dist.computeGaussianBW(true);
			kernel							= new MultiVarNormal(zeroVec, dist.getBandwidth());
		}
		else
		{
			dist.computeGaussianDiagBW(true);
			kernel							= new MultiVarNormal(zeroVec, dist.getDiagBandwidth());
		}
		
		// Perturb
		Enumeration<Integer> keys			= samples.keys();
		while (keys.hasMoreElements())
		{
			Integer index					= keys.nextElement();
			int count						= samples.get(index);
			Particle original				= particles.get(index);
			ArrayList<Double> center		= original.getState();
			particles2.add(new Particle(original.getId(), center, 1.0));
			for (int r = 0; r < count - 1; r++)
			{
				String id					= original.getId() + RSPL + (r + 1);
				double[] randVec			= kernel.sample();
				ArrayList<Double> values	= new ArrayList<>();
				for (int d = 0; d < dimensions; d++)
					values.add(center.get(d) + randVec[d]);
				particles2.add(new Particle(id, values, 1.0));
			}
		}
		return particles2;
	}
	
}
