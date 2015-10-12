package DWR.DMS.PTM;

import java.util.HashMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.Random;
import java.util.Calendar;
import java.io.File;

import javax.xml.bind.JAXBElement.GlobalScope;

import ncsa.hdf.hdf5lib.exceptions.HDF5SymbolTableException;
import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5SimpleReader;
import ch.systemsx.cisd.hdf5.IHDF5SimpleWriter;

/**
 * @author Doug Jackson
 * doug.jackson@noaa.gov
 */
public class BehavedParticle extends Particle
{	
	// Static fields
	public static String outputFilename = MainPTM.getBehaviorOutputFilename();
	public static IHDF5SimpleWriter writer = initializeWriter();
	public static IHDF5SimpleReader reader;
	public static String behaviorParameterFile = MainPTM.getBehaviorInputFilename();
	public static HashMap<String, double[]> ECVecHash = new HashMap<String, double[]>();
	public static HashMap<Integer, Integer> nodeECHash = new HashMap<Integer, Integer>();
	public static boolean ECEnabled;
	public static boolean echoedSetpoints = false;
	public static int tideCountThr;
	public static boolean swimTime = true;
	public static String sunriseTime, sunsetTime;
	public static boolean immortal = false;
	
	// Static initializer
	static 
	{
		double[] ECUpNodes;
		
		// Open the HDF5 file that contains the parameter values
		try
		{
			reader = HDF5Factory.openForReading(behaviorParameterFile);
			System.out.println("Opened " + behaviorParameterFile);
			
			// Create a HashMap to translate a node ID to an index into the EC vectors from QualData
			ECUpNodes = reader.readDoubleArray("QualData/upNodes");
			for(int i=0; i<ECUpNodes.length; i++)
			{
				nodeECHash.put((int)ECUpNodes[i], i);
			}
			
			// Read in the sunrise and sunset times and set the appropriate hours in MainPTM
			sunriseTime = reader.readString("sunriseTime");
			sunsetTime = reader.readString("sunsetTime");
			MainPTM.sunriseHour = Integer.parseInt(sunriseTime.substring(0, 2));
			MainPTM.sunriseMin = Integer.parseInt(sunriseTime.substring(2, 4));
			MainPTM.sunsetHour = Integer.parseInt(sunsetTime.substring(0, 2));
			MainPTM.sunsetMin = Integer.parseInt(sunsetTime.substring(2, 4));
					
			// Read in the tideCountThr and set the class variable in SmartChannel
			tideCountThr = reader.readInt("tideCountThr");
			SmartChannel.setTideCountThr(tideCountThr);	
			
			// Read in the immortal flag
			immortal = reader.readBoolean("immortal");
			
		} catch (HDF5SymbolTableException e)
		{
			System.out.println("Could not find either QualData/upNodes, tideCountThr, dielSwimPeriod, or immortal in the HDF5 file " + behaviorParameterFile);
			System.exit(2);
		}catch (Exception e)
		{
			System.out.println("Failed to open " + behaviorParameterFile + ". Does it exist?");
			System.out.println("Aborting execution");
			System.exit(1);
		} 
				
	}

	// Behavior parameters
	public float lastDecisionAttemptTime;
	public HashMap<Integer, ArrayList<Integer>> possibleChoicesHash;
	public int swimCode;
	public float swimSpeed;
	public float meanSwimSpeed;
	public float stdSwimSpeed;
	public float epsSwimSpeed;
	public boolean variableSwimSpeed;
	public float filterK;
	public float holdThr;
	public float daytimeSwimProb;
	
	// Parameters for velocity memory
	public float timeSinceDecision;
	public double[] velIntMemory; 
	public int velDecisionPeriod;
	
	// Parameters for "downstream" detection
	public float currentDirection;
	public float constProbConfusion;
	public float slopeProbConfusion;
	public double maxProbConfusion;
	public double probConfusion;
	public float confusionFactor;
	public boolean randAssess;
	public float probAssess;
	public float initProbConfusion;
	public boolean enteredSmartChannel;
	
	// tide stage change threshold for behaviors 5 and 6
	public float stageThresholdInc;
	public float stageThresholdDec;
	public float sumStageChanges;
	public float previousStage;
	public boolean stageInitialized = false;
	
	public boolean tideIncreasing = false;
	public double upNodeEC, downNodeEC;
	public Waterbody previousWB;
	
	// makeNodeDecision parameters
	public int nodeDecisionIndex;
	public int[] nodeDecisions;
	public HashMap<Integer, double[][]> weightsTransformation = new HashMap<Integer, double[][]>();
	public double[][] channelPars;
	
	// Channel-specific parameters
	public HashMap<Integer, Double> channelLambda = new HashMap<Integer, Double>();
	public HashMap<Integer, Double> channelOmega = new HashMap<Integer, Double>();
    public HashMap<Integer, Double> channelSwimSpeed = new HashMap<Integer, Double>();
    public HashMap<Integer, Double> channelHoldThr = new HashMap<Integer, Double>();
    public HashMap<Integer, Double> channelConstProbConfusion = new HashMap<Integer, Double>();
    public HashMap<Integer, Double> channelDaytimeSwimProb = new HashMap<Integer, Double>();
	
	// mortality parameters
	Random generator = new Random();
	public double realizedSurvProb;
	
	// Checkpoint parameters
	public int ChippsPassCount;
	public int ExitPassCount;
	public int SWPpassCount;
	public int CVPpassCount;
	public int[] checkpoints;
	public int[] checkpointsPassCount;
	
	////////////////////////////////////////////////////////////////////
	// Instance methods
	////////////////////////////////////////////////////////////////////
	
	public BehavedParticle(ParticleFixedInfo pFI)
	{		
		super(pFI);
		
		// lastDecisionTime == -999.0f indicates that makeNodeDecision has never been run
		lastDecisionAttemptTime = -999.0f;
		possibleChoicesHash = new HashMap<Integer, ArrayList<Integer>>();
		
		// Initialize the stage information
		sumStageChanges = 0.0f;
		previousStage = 0.0f;
		
		// Initialize the EC memory
		upNodeEC = downNodeEC = 0.0;		
				
		// Read in the parameter values
		try
		{
			swimCode = reader.readInt("swimCode");
			stdSwimSpeed = reader.readFloat("stdSwimSpeed");
			filterK = reader.readFloat("filterK");
			holdThr = reader.readFloat("holdThr"); //NA
			daytimeSwimProb = reader.readFloat("daytimeSwimProb"); // NA
			variableSwimSpeed = reader.readBoolean("variableSwimSpeed");
			stageThresholdInc = reader.readFloat("stageThresholdInc");
			stageThresholdDec = reader.readFloat("stageThresholdDec");
			velDecisionPeriod = reader.readInt("velDecisionPeriod");
			velIntMemory = new double[velDecisionPeriod];
			constProbConfusion = reader.readFloat("constProbConfusion"); //NA
			slopeProbConfusion = reader.readFloat("slopeProbConfusion");
			randAssess = reader.readBoolean("randAssess");
			probAssess = reader.readFloat("probAssess");
			initProbConfusion = reader.readFloat("initProbConfusion");
			currentDirection = 1.0f;
			
			nodeDecisionIndex = 0;
			nodeDecisions = reader.readIntArray("nodeDecisions");
						
			for(int i=0; i<nodeDecisions.length; i++)
			{
				if(!echoedSetpoints) System.out.println("Node decision " + Integer.toString(i) + " = " + Integer.toString(nodeDecisions[i]));
				weightsTransformation.put(nodeDecisions[i], reader.readDoubleMatrix("weightsTransformation/weightsTransformation_" 
						+ Integer.toString(nodeDecisions[i])));
				if(!echoedSetpoints) writer.writeDoubleMatrix("weightsTransformation/weightsTransformation_" + Integer.toString(nodeDecisions[i]), 
					weightsTransformation.get(nodeDecisions[i]));				
			}
			
			// See if we need to update EC
			ECEnabled=false;
			if(swimCode==7)	ECEnabled=true;
			else
			{
				for(int i=0; i<nodeDecisions.length; i++)
				{
					if(nodeDecisions[i]==2) ECEnabled=true;
				}
			}
				
            // channelPars
			// read in HDF5 file, store channel specific values in hashmaps
            // 0: channel
            // 1: lambda
            // 2: omega
            // 3: swimSpeed
			// 4: HoldThr
			// 5: ConstProbConfusion
			// 6: daytimeSwimProb
			
			channelPars = reader.readDoubleMatrix("channelPars");
			for(int i=0; i<channelPars.length; i++)
			{
				channelLambda.put((int)Math.floor(channelPars[i][0]), channelPars[i][1]);
				channelOmega.put((int)Math.floor(channelPars[i][0]), channelPars[i][2]);
                channelSwimSpeed.put((int)Math.floor(channelPars[i][0]), channelPars[i][3]);
                channelHoldThr.put((int)Math.floor(channelPars[i][0]), channelPars[i][4]);
                channelConstProbConfusion.put((int)Math.floor(channelPars[i][0]), channelPars[i][5]);
                channelDaytimeSwimProb.put((int)Math.floor(channelPars[i][0]), channelPars[i][6]);
			}
			
			ChippsPassCount = 0;
			ExitPassCount = 0;
			CVPpassCount = 0;
			SWPpassCount = 0;
			
			// Read the checkpoints and sort them so we can use Arrays.binarySearch() to see if the
			// list contains a particular checkpoint
			checkpoints = reader.readIntArray("checkpoints");
			Arrays.sort(checkpoints);
			checkpointsPassCount = new int[checkpoints.length];	

			// Write the parameter values to the output file
			if(!echoedSetpoints)
			{
				System.out.println("swimCode = " + swimCode + 
						", stdSwimSpeed=" + stdSwimSpeed + ", variableSwimSpeed=" + variableSwimSpeed +
						", velDecisionPeriod=" + velDecisionPeriod + 
						", constProbConfusion = " + constProbConfusion +
						", slopeProbConfusion=" + slopeProbConfusion +
						", randAssess=" + randAssess +
						", probAssess=" + probAssess +
						", initProbConfusion=" + initProbConfusion +
						", tideCountThr=" + tideCountThr + ", holdThr=" + holdThr +
						", stageThresholdInc=" + stageThresholdInc + ", stageThresholdDec=" + stageThresholdDec +
						", daytimeSwimProb=" + daytimeSwimProb + ", sunriseTime=" + sunriseTime + ", sunsetTime=" + sunsetTime +  // daytimeSwimProb NA
						", checkpoints = " + Arrays.toString(checkpoints) + ", immortal=" + immortal);
				
				writer.writeInt("swimCode", swimCode);
				writer.writeFloat("stdSwimSpeed", stdSwimSpeed);
				writer.writeBoolean("variableSwimSpeed", variableSwimSpeed);
				writer.writeInt("velDecisionPeriod", velDecisionPeriod);
				writer.writeFloat("constProbConfusion", constProbConfusion);
				writer.writeFloat("slopeProbConfusion", slopeProbConfusion);
				writer.writeBoolean("randAssess", randAssess);
				writer.writeFloat("probAssess", probAssess);
				writer.writeFloat("initProbConfusion", initProbConfusion);
				writer.writeInt("tideCountThr", tideCountThr);
				writer.writeFloat("filterK", filterK);
				writer.writeFloat("holdThr", holdThr);
				writer.writeFloat("stageThresholdInc", stageThresholdInc);
				writer.writeFloat("stageThresholdDec", stageThresholdDec);
				writer.writeFloat("daytimeSwimProb", daytimeSwimProb); // NA
				writer.writeString("sunriseTime", sunriseTime, 4);
				writer.writeString("sunsetTime", sunsetTime, 4);
				writer.writeIntArray("nodeDecisions", nodeDecisions);
				writer.writeDoubleMatrix("channelPars", channelPars);
				writer.writeIntArray("checkpoints", checkpoints);
				writer.writeBoolean("immortal", immortal);
			}
			
			echoedSetpoints = true; 
			
		} catch (HDF5SymbolTableException e)
		{
			System.out.println("Could not find one of the parameters in the HDF5 file " + behaviorParameterFile);
			System.exit(2);
		}
		
		// Initialize the swimSpeed to 0.0 just to be safe
		swimSpeed = 0.0f;
		meanSwimSpeed = 0.0f;
		epsSwimSpeed = 0.0f;
		
		// Write the realized swimSpeed for each particle to the output file, but 
		// only do this if stdSwimSpeed>0 and variableSwimSpeed==false
		if(stdSwimSpeed>0 && variableSwimSpeed==false)
		{
			writer.writeFloat("swimSpeed/particleNum/" + Integer.toString(this.getId()), swimSpeed);
		}	
		
		// Initialize confusionFactor (-1 is confused)
		if(generator.nextDouble()<initProbConfusion)
		{
			confusionFactor = -1.0f;
		}
		else
		{
			confusionFactor = 1.0f;
		}
				
		// Set maximum prob confusion
		maxProbConfusion = 0.5;
		
		// Initializations
		probConfusion = initProbConfusion;
		realizedSurvProb = 1.0;
		enteredSmartChannel = false;
	}

	// Main swim behavior stuff here
	/**
	 * Create some swimming behaviors for movement along the channel axis.
	 * 0: never move (used to test particle environment queries) 
	 * 1: passive drift (basic PTM)
	 * 2: swim downstream all the time
	 * 3: swim downstream when the flow is towards the bay (positive), otherwise drift 
	 * 4: hold still when flow is negative, swim towards bay when positive
	 * 5: swim with the flow when tide falls, otherwise drift
	 * 6: swim with the flow when tide falls, otherwise hold still
	 * 7: swim towards higher salinity (EC) all the time
	 * 8: swim downstream during certain times (e.g., nighttime), otherwise hold (Note: with the addition of dielSwimPeriod, swimCode 8==2)
	 * 9: swim "downstream" at a constant swimSpeed; with advection
	 * 10: swim "downstream" at a constant swimSpeed, with "downstream" varying by channel; with advection
	 * 11: swim "downstream" when the velocity is above -holdThr, with "downstream" varying by channel; hold otherwise
	 */
	
	/**
	 * Externally induced Deterministic (the effect of flow in the channel)
	 */
	@Override
	protected float calcXVelocityExtDeterministic()
	{
		float channelDir;
		float flowVelocity = getFlowVelocity();
		
		float particleVelocity = 0.0f;	
		
		if(wb instanceof SmartChannel)
		{
			channelDir = ((SmartChannel)wb).getChannelDir();
		}
		else
		{
			channelDir = 1.0f;
		}

		// Nocturnal and diurnal behaviors may differ
		if(swimTime)
		{
			switch (swimCode)
			{
				case 0:
					particleVelocity = 0.0f;
					break;
	
				case 1: case 2: case 3: case 5: case 8: case 7: case 9: case 10:
					particleVelocity = flowVelocity;
					break;
	
				case 4:
					if (flowVelocity > 0)
						particleVelocity = flowVelocity; // get into flow when positive
					else
						particleVelocity = 0.0f; // hold position while flows are negative
					break;
					
				case 6:
					if (tideIncreasing)
						particleVelocity = 0.0f; // hold still if tide is rising
					else
						particleVelocity = flowVelocity;  // on the falling tide, go with the flow
					break;
				
				case 11:
					if(flowVelocity*channelDir*confusionFactor>-holdThr)
					{
						particleVelocity = flowVelocity;
					}
					else
					{
						particleVelocity = 0.0f;
					}
					
					break;
					
				default:
					throw new IllegalArgumentException("Unrecognized swimCode");
				
			}
		}
		else
		{
			switch(swimCode)
			{
				case 0: case 2: case 3: case 4: case 5: case 6: case 7: case 8: case 9: case 10: case 11:
					particleVelocity = 0.0f;
					break;
				case 1:
					particleVelocity = flowVelocity;
					break;
				
				default:
					throw new IllegalArgumentException("Unrecognized swimCode");
			}
		}
		
		return particleVelocity;
	}
	
	/**
	 * Internally induced Deterministic (the particle's swimming behavior)
	 */
	@Override
	protected float calcXVelocityIntDeterministic()
	{	
		float flowVelocity = getFlowVelocity();
		float swimVelocity = 0.0f;
		float channelDir;
		
		if(wb instanceof SmartChannel) {
			channelDir = ((SmartChannel)wb).getChannelDir();

		} else {
			channelDir = 1.0f;
		}
		
		// Nocturnal and diurnal behaviors may differ
		if(swimTime)
		{
			swimSpeed = meanSwimSpeed + epsSwimSpeed;
			
			switch (swimCode)
			{
				case 0: case 1:
					swimVelocity = 0.0f;
					break;
	
				case 2: case 8:
					swimVelocity = swimSpeed;
					break;
				
				// swim with flow when positive
				case 3: case 4: 
					if (flowVelocity > 0.0f)
						swimVelocity = swimSpeed;
					else
						swimVelocity = 0.0f;
					break;
	
				// swim with the flow if the tide is falling
				case 5: case 6:
					if (tideIncreasing)
						swimVelocity = 0.0f; 
					else
					{
						// swim in the direction of the flow
						if (flowVelocity >= 0.0f)
							swimVelocity = swimSpeed; // swim "downstream"
						else
							swimVelocity = (-1.0f * swimSpeed); // swim "upstream" 'cause that's where the flow is going.
					}
					break;
				
				// Swim towards the node with highest salinity
				// In the interest of performance, upNodeEC and downNodeEC are only updated after a new waterbody
				// is entered.
				case 7:
					if(wb instanceof Channel)
					{
						if(downNodeEC >= upNodeEC)
							swimVelocity = swimSpeed;
						else
							swimVelocity = -swimSpeed;
					}
					else
						swimVelocity = 0.0f;
					break;
				
				// Swim "downstream" all the time, with "downstream" determined based on the direction of flow
				// experienced by the fish over the past velDecisionPeriod hours
				case 9:
					swimVelocity = swimSpeed*currentDirection*confusionFactor;			
					break;
					
				// Swim "downstream" all the time, with "downstream" determined by the average direction of flow
				// in the current channel over the last tideCountThr cycles
				case 10:
					swimVelocity = swimSpeed*channelDir*confusionFactor;
					break;
					
				// Swim "downstream" if the "upstream" flow is less than holdThr, with "downstream" determined by the average direction of flow
				// in the current channel over the last tideCountThr cycles. Hold otherwise.
				case 11:
					if(flowVelocity*channelDir*confusionFactor>-holdThr)
					{
						swimVelocity = swimSpeed*channelDir*confusionFactor;
					}
					else
					{
						swimVelocity = 0.0f;
					}
					break;
					
				default:
					throw new IllegalArgumentException("Unrecognized swimCode");

			}
		}
		else
		{
			swimVelocity = 0.0f;				
		}
		
		return swimVelocity;
	}

	@Override
	/**
	 * Decide which WaterBody to enter into next
	 */
	protected void makeNodeDecision()
	{

		ArrayList<Integer> possibleChoices;
		
		int choiceIndex = 0;
		int numWaterBodies = nd.getNumberOfWaterbodies();
		double [] weightVector;
		double sumWeightVector;
		float tempChannelLength;
		int [] indexVector;
		boolean madeDecision = false;
		float minFlow = 0.0f;
		
		previousWB = wb;
		
		// Send message to observer about change
		if (observer != null)
		{
			observer.observeChange(ParticleObserver.NODE_CHANGE, this);
		}
		
		// Loop until a decision is made
		do 
		{
			// Clear possibleChoicesHash, reset nodeDecisionIndex, and possibly check to see
			// if the fish becomes confused unless this is a retry of an unsuccessful choice
			if(((float) Globals.currentModelTime + tmLeft) != lastDecisionAttemptTime)
			{				
				nodeDecisionIndex = 0;
				possibleChoicesHash.clear();
				
				// If this is a junction, check to see if the fish is confused in this stretch
				if(nd.getNumChannels()>2 && !randAssess)
				{
					checkConfusion();
				}
			}
			
			// Remember the last currentModelTime && tmLeft combination when we attempt to make a decision
			lastDecisionAttemptTime = ((float) Globals.currentModelTime + tmLeft);
			
			// If this is the first time we've been to this node during this decision, create a 
			// new list
			if(!possibleChoicesHash.containsKey(nd.getEnvIndex()))
			{
				//System.out.println("Entered possibleChoicesHash creation for node " + 
					//	Integer.toString(nd.getEnvIndex()));
				possibleChoices = new ArrayList<Integer>();
				possibleChoicesHash.put(nd.getEnvIndex(), possibleChoices);
			}
			
			// Retrieve possibleChoices for this node
			possibleChoices = possibleChoicesHash.get(nd.getEnvIndex());
			
			// If the possibleChoices list for this node is empty but there's another decision
			// type to try, add all of the waterbodies to the possibleChoices list.
			if(possibleChoices.isEmpty())
			{
				// Is there another decision type to try?
				if(nodeDecisionIndex<nodeDecisions.length)
				{
					for(int i=0; i<numWaterBodies; i++)
					{
						// Prevent movement into boundaries, but allow entry into 901 and 915 (where water is taken out)
						if(!(nd.getWaterbody(i) instanceof Boundary) || nd.getWaterbodyEnvIndex(i)==901 || nd.getWaterbodyEnvIndex(i)==915)
						{
							possibleChoices.add(i);
						}
					}
					possibleChoicesHash.put(nd.getEnvIndex(), possibleChoices);
				}
				// If there are no choices, either move out into the channel a small amount if it's a 
				// dead end, or wait
				else
				{
					
					if(numWaterBodies==1)
					{
						x = getPerturbedXLocation();
					}
					else
					{
						particleWait = true;
					}
					return;	
				}
			}
									
			weightVector = new double[possibleChoices.size()];
			sumWeightVector = 0.0;
			indexVector = new int[possibleChoices.size()];
			
			switch(nodeDecisions[nodeDecisionIndex])
			{
			
			// Outflow-based decision
			case 0:
				
				for(int i=0 ; i<possibleChoices.size(); i++)
				{
					weightVector[i] = nd.getFilterOp(possibleChoices.get(i))*nd.getOutflow(possibleChoices.get(i));
					sumWeightVector += weightVector[i];
					indexVector[i] = possibleChoices.get(i);
				}
				break;
				
			// Channel width-based decision
			case 1:
				
				for(int i=0; i<possibleChoices.size(); i++)
				{
					if(nd.getWaterbody(possibleChoices.get(i)) instanceof Channel)
					{
						// Get channel width at either the beginning of the channel or the end, depending
						// on which end of the channel we're at
						if(((Channel)nd.getWaterbody(possibleChoices.get(i))).getUpNodeId()==nd.getEnvIndex())
						{
							weightVector[i] = ((Channel)nd.getWaterbody(possibleChoices.get(i))).getWidth(0.0f);
						}
						else
						{
							tempChannelLength = ((Channel)nd.getWaterbody(possibleChoices.get(i))).getLength();
							weightVector[i] =
									((Channel)nd.getWaterbody(possibleChoices.get(i))).getWidth(tempChannelLength);
						}
					}
					else
					{
						weightVector[i] = 0.0;
					}
					sumWeightVector += weightVector[i];
					indexVector[i] = possibleChoices.get(i);
				}
				break; 
			
			// Salinity (EC)-based decision
			case 2:
				
				for(int i=0; i<possibleChoices.size(); i++)
				{
					if(nd.getWaterbody(possibleChoices.get(i)) instanceof Channel)
					{
						// Get salinity at the opposite end of the channel
						if(((Channel)nd.getWaterbody(possibleChoices.get(i))).getUpNodeId()==nd.getEnvIndex())
						{
							// false = use downNode
							weightVector[i] = lookupEC((Channel)nd.getWaterbody(possibleChoices.get(i)), false);
						}
						else
						{
							// true = use upNode
							weightVector[i] = lookupEC((Channel)nd.getWaterbody(possibleChoices.get(i)), true);
						}
					}
					else
					{
						weightVector[i] = 0.0;
					}
					sumWeightVector += weightVector[i];
					indexVector[i] = possibleChoices.get(i);
				}
				break; 
			
			// Total flow-based decision
			// Choose waterbody based on flow relative to the channel with the lowest (or most negative)
			// flow. Exclude the waterbody that the particle is coming from.
			case 3:
				
				// Find the lowest (or most negative) flow
				for(int i=0; i<possibleChoices.size(); i++)
				{
					if(i==0)
					{
						minFlow = nd.getFilterOp(possibleChoices.get(i))*nd.getSignedOutflow(possibleChoices.get(i));
					}
					else if((nd.getFilterOp(possibleChoices.get(i))*nd.getSignedOutflow(possibleChoices.get(i)))<minFlow)
					{
						minFlow =  nd.getFilterOp(possibleChoices.get(i))*nd.getSignedOutflow(possibleChoices.get(i));
					}
				}
				
				// weightVector is equal to flow-minFlow+1.0f (the +1.0f is to ensure that sumWeightVector>0 if there's
				// a valid choice). Set weightVector = 0.0f for the waterbody the fish is coming from. Also set weightVector = 0.0f
				// if the outflow is zero, which indicates that the outflow is blocked by a gate.
				for(int i=0; i<possibleChoices.size(); i++)
				{
					if(nd.getWaterbody(possibleChoices.get(i))==previousWB)
					{
						weightVector[i] = 0.0f;
					}
					else if(nd.getSignedOutflow(possibleChoices.get(i))==0.0f)
					{
						weightVector[i] = 0.0f;
					}
					else
					{
						weightVector[i] = nd.getFilterOp(possibleChoices.get(i))*nd.getSignedOutflow(possibleChoices.get(i)) -
								minFlow + 1.0f;
						sumWeightVector += weightVector[i];
					}
					indexVector[i] = possibleChoices.get(i);
				}
				
				break;
				
			default:
					throw new IllegalArgumentException("Unrecognized decision type in nodeDecisions");
			}
					
			// If there are viable options, make a choice. Otherwise, clear possibleChoices and
			// increment nodeDecisionIndex (to trigger either the next decision type or the dead end logic)
			if(sumWeightVector > 0)
			{
				// Apply the user-defined transformation
				weightVector = transformWeightVector(weightVector, sumWeightVector);
				
				// Make the choice and remove it from the list of possible future choices for this node
				choiceIndex = weightedChoice(weightVector);
				possibleChoices.remove(possibleChoices.indexOf(indexVector[choiceIndex]));
				possibleChoicesHash.put(nd.getEnvIndex(), possibleChoices);
				madeDecision = true;
				
				// If we don't have any remaining options with this decision type, increment nodeDecisionIndex in 
				// case this attempt fails
				if(possibleChoices.isEmpty())
				{
					nodeDecisionIndex ++;
				}
			}
			else
			{
				possibleChoices.clear();
				nodeDecisionIndex ++;
			}
			
		} while(!madeDecision);		

		// Get a pointer to the water body that the particle entered
		wb = nd.getWaterbody(indexVector[choiceIndex]);
		
		// Update parameters, etc., when entering a new channel
		enterChannel();
		
		// Determine if the particle entered from the upNode side or the downNode side of 
		// a channel and set currentDirection accordingly
		if (wb.getPTMType() == Waterbody.CHANNEL)
		{
			if (((Channel) wb).getDownNodeId() == nd.getEnvIndex())
				currentDirection = -1.0f;
			else
			{
				currentDirection = 1.0f;
			}
		}
		else
		{
			currentDirection = 1.0f;
		}
		
		// See if the particle has passed any checkpoints
		checkCheckpoints();
				
		// Remember the upNode and downNode ECs (if necessary)
		if(ECEnabled)
		{
			if(wb instanceof Channel)
			{
				upNodeEC = lookupEC((Channel)wb, true);
				downNodeEC = lookupEC((Channel)wb, false);
			}
			else
			{
				upNodeEC = downNodeEC = 0.0;
			}
		}
		
		// Send message to observer about change
		if (observer != null)
			observer.observeChange(ParticleObserver.WATERBODY_CHANGE, this);
		
		// Set x as beginning of Channel...
		x = getXLocationInChannel();
	}
	
	// Transform the values in weightVector
	protected double[] transformWeightVector(double[] weightVector, double sumWeightVector)
	{
		double[] newWeightVector = new double[weightVector.length];
		System.arraycopy(weightVector, 0, newWeightVector, 0, weightVector.length);
		
		double sumNewWeightVector = 0.0;
		
		for(int i=0; i<newWeightVector.length; i++)
		{
			// normalize newWeightVector
		    newWeightVector[i] /= sumWeightVector;
		
			// apply the transformation
		    newWeightVector[i] *= interpLinear(weightsTransformation.get(nodeDecisions[nodeDecisionIndex])[0],
		    		weightsTransformation.get(nodeDecisions[nodeDecisionIndex])[1], newWeightVector[i]);
			
			sumNewWeightVector += newWeightVector[i];
		}
		
		sumNewWeightVector = sumNewWeightVector + 0.0;
		
		// renormalize newWeightVector
		for(int i=0; i<newWeightVector.length; i++)
		{
			newWeightVector[i] /= sumNewWeightVector;
		}
		 
		return newWeightVector;
	}
	
	public double interpLinear(double[] x, double[] y, double xi) throws IllegalArgumentException
	{

		double[] dx = new double[x.length-1];
		double[] dy = new double[x.length-1];
		double[] slope = new double[x.length-1];
		double[] intercept = new double[x.length-1];
		int location;
		double yi;
		
		if (x.length != y.length)
		{
			throw new IllegalArgumentException("weightsTransformation vectors must be the same length.");
		}
		if (x.length == 1)
		{
			throw new IllegalArgumentException("weightsTransformation vectors must contain more than one value.");
		}

		// Calculate the line equation between each point
		for (int i=0; i<x.length-1; i++)
		{
			dx[i] = x[i+1]-x[i];
			if (dx[i]==0)
			{
				throw new IllegalArgumentException("weightsTransformation vectors must be montotonic. A duplicate " +
						"x-value was found");
			}
			if (dx[i] < 0)
			{
				throw new IllegalArgumentException("X must be sorted");
			}
			dy[i] = y[i + 1]-y[i];
			slope[i] = dy[i]/dx[i];
			intercept[i] = y[i]-x[i]*slope[i];
		}

		// Perform the interpolation here
		if ((xi>x[x.length-1]) || (xi<x[0]))
		{
			yi = 1.0;
		}
		else
		{
			// binarySearch(double[], double) returns index of the search key, if it is contained in the array; 
			// otherwise, (-(insertion point) - 1). The insertion point is defined as the point at which the key 
			// would be inserted into the array: the index of the first element greater than the key, or a.length 
			// if all elements in the array are less than the specified key. Note that this guarantees that the 
			// return value will be >= 0 if and only if the key is found.
			location = Arrays.binarySearch(x, xi);
			if (location < -1)
			{
				location = -location-2;
				yi = slope[location]*xi+intercept[location];
			}
			else
			{
				yi = y[location];
			}
		}

		return yi;
	}
	
	public float getFlowVelocity()
	{
		float flowVelocity;
		
		if(wb instanceof SmartChannel)
			flowVelocity = (((Channel) wb).getVelocity(x, y, z, channelVave, channelWidth,
					channelDepth));
		else
			flowVelocity = 0.0f;
		
		if(Float.isNaN(flowVelocity))
		{
			System.out.println("x=" + x + ", y=" + y + ", z=" + z + ", ID=" + wb.getEnvIndex());
			System.exit(2);
		}
		return flowVelocity;
	}
	
	public void checkCheckpoints()
	{
		int checkpointIndex; 
		
		// Check to see if the particle has reached Chipps Island yet
		if(previousWB!=null)
		{
			if(wb instanceof Channel && 
					((wb.getEnvIndex()==422 || wb.getEnvIndex()==417) && 
							(previousWB.getEnvIndex()==275 || previousWB.getEnvIndex()==281 || previousWB.getEnvIndex()==278)))
			{
				ChippsPassCount++;
				recordCheckpoint(this, "Chipps", ChippsPassCount);
				
				// Write the realized survival to the output file
				writer.writeDouble("realizedSurvProb/" + this.getId(), realizedSurvProb);
			}
		}		
		
		// Check to see if the particle has exited the system
		if(nd.getEnvIndex()==412)
		{
			ExitPassCount++;
			recordCheckpoint(this, "Exit", ExitPassCount);
		}
		
		// Check to see if the particle was exported via SWP
		if(wb instanceof Reservoir)
		{
			if(((Reservoir)wb).getName().equals("clifton_court"))
			{
				SWPpassCount++;
				recordCheckpoint(this, "SWP", SWPpassCount);
				isDead = true;
			    recordDeath(this);
			}
		}
				
		if(wb.getEnvIndex()==204)
		{
			CVPpassCount++;
			recordCheckpoint(this, "CVP", CVPpassCount);
			isDead = true;
		    recordDeath(this);
		}
		
		// Check to see if the particle passed one of the other checkpoints
		checkpointIndex = Arrays.binarySearch(checkpoints, nd.getEnvIndex());
		if(checkpointIndex>=0)
		{
			checkpointsPassCount[checkpointIndex]++;
			recordCheckpoint(this, Integer.toString(nd.getEnvIndex()), checkpointsPassCount[checkpointIndex]);
		}
		
	}
	
	public void updateTideIncreasing()
	{
		// Keep track of the stage changes since the last tideIncreasing assessment
		// was made
		if(!stageInitialized)
		{
			previousStage = channelStage;
			stageInitialized = true;
		}
		sumStageChanges += (channelStage-previousStage);
				
		// If the sums of the tide increases or decreases exceed stageThreshold, make a 
		// tideIncreasing assessment, i.e., decide if the tide is increasing or decreasing.
		// If the fish has already decided that the tide is increasing (decreasing), any
		// further increases (decreases) just extend the assumed peak (valley) of the tide cycle;
		// if the tide backs off from the peak by stageThreshold, then the fish assumes that the
		// tide direction has changed.
		if((sumStageChanges>stageThresholdInc) || (sumStageChanges>0 && tideIncreasing))
		{
			tideIncreasing = true;
			sumStageChanges = 0.0f;
		}
		else if((sumStageChanges<(-stageThresholdDec)) || (sumStageChanges<0 && !tideIncreasing))
		{
			tideIncreasing = false;
			sumStageChanges = 0.0f;
		}
		// otherwise don't change tideIncreasing	

		previousStage = channelStage;
		
	}
	
	// Check to see if the fish is confused in this stretch in terms of their assessment
	// of which way is "downstream"
	public void checkConfusion()
	{
		// Fish becomes confused with probability=probConfusion
		if(generator.nextDouble()<probConfusion)
		{
			confusionFactor = -1.0f;
		}
		else
		{
			confusionFactor = 1.0f;
		}	
	}
	
	// Update parameters, etc., when a fish enters a new channel
	public void enterChannel()
	{
		// Update probConfusion whenever a new SmartChannel is entered
		updateProbConfusion();
        
		if(wb instanceof SmartChannel)
		{
			meanSwimSpeed = channelSwimSpeed.get(wb.getEnvIndex()).floatValue();
			holdThr = channelHoldThr.get(wb.getEnvIndex()).floatValue();
			constProbConfusion = channelConstProbConfusion.get(wb.getEnvIndex()).floatValue();
			daytimeSwimProb = channelDaytimeSwimProb.get(wb.getEnvIndex()).floatValue();
		}
	}
    
	// Update the probability of confusion based on the signalToNoise of the current channel
	public void updateProbConfusion()
	{
		double signalToNoise, lnSignalToNoise, term;
		
		// Only change probConfusion if the fish is currently in a smart channel.
		if(wb instanceof SmartChannel)
		{
			signalToNoise = ((SmartChannel) wb).getSignalToNoise();
			lnSignalToNoise = Math.log(Math.max(1E-10, signalToNoise));
			term = Math.exp(constProbConfusion + slopeProbConfusion*lnSignalToNoise);
			probConfusion = maxProbConfusion*term/(1+term);
		}
	}
	
	@Override
	public void insert()
	{
		super.insert();
		
		// Record the insertion time
		recordInsertion(this);		
	}
	@Override 
	// Override updatePosition to implement BehavedParticle actions that occur every 15 minutes
	public void updatePosition(float delT)
	{		
		// Update the probability of confusion the first time that the fish enters a SmartChannel
		if(enteredSmartChannel==false)
		{
			if(wb instanceof SmartChannel)
			{
				updateProbConfusion();
				enteredSmartChannel=true;
			}
		}
		
		// Determine if "downstream" should be assessed, i.e., whether the fish becomes
		// confused/unconfused during this time step
		if(randAssess)
		{
			if(generator.nextDouble()<probAssess)
			{
				checkConfusion();
			}
		}
		
		// Draw a new epsSwimSpeed every delT seconds (which will be 15 minutes)
		if(variableSwimSpeed)
		{
			// Draw epsSwimSpeed from a normal distribution
			epsSwimSpeed = (float) (generator.nextGaussian()*stdSwimSpeed);
		}
		
		// Clear the memory of time spent and movement in the previous time step
		movementTimeDistance.clear();
				
		super.updatePosition(delT);		
	}
	
	@Override
	public void updateParticleParameters(float timeStep)
	{
		int lastHoursSinceDecision = (int) Math.floor(timeSinceDecision/3600.0f);
		int hoursSinceDecision, lastMemIndex, memIndex, minHoursSinceDecision, thisMemIndex;
		float thisTimeRemainder, thisTimeStep, lastTimeStep;
		double velIntSum;
		
		super.updateParticleParameters(timeStep);
		
		timeSinceDecision += timeStep;
		hoursSinceDecision = (int) Math.floor(timeSinceDecision/3600.0f);
		
		// Calculate the amount of this time step that belongs to the current hour and the 
		// amount that belongs to the previous hour
		thisTimeRemainder = timeSinceDecision % 3600.0f;
		thisTimeStep = Math.min(timeStep, thisTimeRemainder);
		lastTimeStep = timeStep - thisTimeStep;
		
		// Calculate the indices into this and the previous memory location
		lastMemIndex = lastHoursSinceDecision % velDecisionPeriod;
		memIndex = hoursSinceDecision % velDecisionPeriod;
		
		// Reset the integrator if this is the first time we've added to it during the current hour
		if(hoursSinceDecision != lastHoursSinceDecision)
		{
			velIntMemory[memIndex] = 0.0; 
		}
		
		// Add +1*timeStep to the integrators if the fish is traveling in the "downstream" direction, i.e., with
		// the flow; subtract if they're going the "wrong" direction
		if(getFlowVelocity()*currentDirection>0)
		{
			velIntMemory[lastMemIndex] += (double)(lastTimeStep);
			velIntMemory[memIndex] += (double)(thisTimeStep);
		}
		else
		{
			velIntMemory[lastMemIndex] -= (double)(lastTimeStep);
			velIntMemory[memIndex] -= (double)(thisTimeStep);
		}
			
				
		// Calculate the integral. All memory locations except the current one are simply added since
		// they represent a time period of one hour.
		velIntSum = velIntMemory[memIndex]*(thisTimeRemainder/3600.0f);
		minHoursSinceDecision = Math.max(0, hoursSinceDecision - velDecisionPeriod + 1);
		for(int h=minHoursSinceDecision; h<hoursSinceDecision; h++)
		{
			thisMemIndex = h % velDecisionPeriod;
			velIntSum += velIntMemory[thisMemIndex];
		}
		
		// Change the direction if the integrated velocity is negative (which indicates that
		// the fish has been going "upstream") and hoursSinceDecision>=velDecisionPeriod
		if(velIntSum<0 && hoursSinceDecision>=velDecisionPeriod)
		{
			currentDirection = -currentDirection;
			timeSinceDecision = 0.0f;
		}		
	}
	
	@Override
	// Particle mortality
	protected void checkHealth()
	{
		double survivalProb = 1.0;
		int channelNum;
		double lambda=0.0, omega=0.0, time, distance;
		double[] timeDistance;
		
		// Loop through all of the movementTimeDistance entries
		for(Map.Entry<Integer, double[]> entry : movementTimeDistance.entrySet())
		{
			
			channelNum = entry.getKey();
			if(channelLambda.containsKey(channelNum))
			{
				lambda = channelLambda.get(channelNum);
			}
			else
			{
				System.out.println("Could not find lambda for channel " + Integer.toString(channelNum));
			}
			
			if(channelOmega.containsKey(channelNum))
			{
				omega = channelOmega.get(channelNum);
			}
			else
			{
				System.out.println("Could not find omega for channel " + Integer.toString(channelNum));
			}
			
			// Retrieve the amount of time and distance traveled in this channel
			timeDistance = entry.getValue();
			time = timeDistance[0];
			distance = timeDistance[1];
			
			// From Anderson, J. J., Gurarie, E., & Zabel, R. W. (2005). Mean free-path length theory of 
			// predator–prey interactions: Application to juvenile salmon migration. Ecological Modelling, 
			// 186(2), 196–211. doi:10.1016/j.ecolmodel.2005.01.014
			// Units of lambda are feet; units of omega are feet/sec.
			survivalProb *= 
					Math.exp((-1.0/lambda)*Math.sqrt((Math.pow(distance, 2.0) + (Math.pow(omega, 2.0)*Math.pow(time, 2.0)))));
			
		}
		
		realizedSurvProb *= survivalProb;
		
		// Particle dies with P(1-survivalProb) 
		if(generator.nextDouble()>survivalProb && !immortal)
		{
			isDead = true;
//		    observer.observeChange(ParticleObserver.DEATH,this);
		    recordDeath(this);
		}
	}
	
	public void checkSwimTime()
	{
		// If it's nighttime, the fish will swim.
		// If it's daytime, fish swims with probability daytimeSwimProb. 
		if(!MainPTM.isDaytime || (MainPTM.isDaytime && generator.nextDouble()<daytimeSwimProb)) {
			swimTime = true;
		} else {
			swimTime = false;
		}
	}
	
	////////////////////////////////////////////////////////////////////
	// Class methods
	////////////////////////////////////////////////////////////////////
	public static IHDF5SimpleWriter initializeWriter()
	{	
		// Delete outputFilename if it already exists
		File testFile = new File(outputFilename);
		try
		{
			if(testFile.exists())
			{
				testFile.delete();
			}

		} catch (Exception e)
		{
			System.out.println("Cannot delete old output file " + outputFilename);
			System.out.println("Try deleting it manually and restarting. Aborting execution");
			System.exit(1);
		}
		IHDF5SimpleWriter w = HDF5Factory.open(outputFilename);
		System.out.println("Opened " + outputFilename + " for writing.");
		return w;
	}
	
	public static void recordDeath(BehavedParticle bP)
	{
		int julianMin = Globals.currentModelTime;
		String modelDate, modelTime;
		
		modelDate = Globals.getModelDate(julianMin);
		modelTime = Globals.getModelTime(julianMin);
		
		// Write to the HDF5 file
		writer.writeString("died/particleNum/" + Integer.toString(bP.getId()) + "/modelDate", modelDate, 9);
		writer.writeInt("died/particleNum/" + Integer.toString(bP.getId()) + "/modelTime", new Integer(modelTime).intValue());
		writer.writeInt("died/particleNum/" + Integer.toString(bP.getId()) + "/waterBody", bP.getCurrentWaterbody().getEnvIndex());
		
	}
	
	public static void recordCheckpoint(BehavedParticle bP, String checkpoint, int passCount)
	{
		int julianMin = Globals.currentModelTime;
		String modelDate, modelTime;
		
		modelDate = Globals.getModelDate(julianMin);
		modelTime = Globals.getModelTime(julianMin);
		
		// Write to the HDF5 file
		writer.writeString(checkpoint + "/particleNum/" + Integer.toString(bP.getId()) + "/modelDate_" + passCount, modelDate, 9);
		writer.writeInt(checkpoint + "/particleNum/" + Integer.toString(bP.getId()) + "/modelTime_" + passCount, new Integer(modelTime).intValue());
	}
	
	public static void recordInsertion(BehavedParticle bP)
	{
		int julianMin = Globals.currentModelTime;
		String modelDate, modelTime;
		
		modelDate = Globals.getModelDate(julianMin);
		modelTime = Globals.getModelTime(julianMin);
		
		// Write to the HDF5 file
		writer.writeString("inserted/particleNum/" + Integer.toString(bP.getId()) + "/modelDate", modelDate, 9);
		writer.writeInt("inserted/particleNum/" + Integer.toString(bP.getId()) + "/modelTime", new Integer(modelTime).intValue());	
		writer.writeInt("inserted/particleNum/" + Integer.toString(bP.getId()) + "/insertionNode", bP.nd.getEnvIndex());
		
	}
	
	public static double lookupEC(Channel c, Boolean useUpNode)
	{
		double[] ECVec;
		double minEC=0.001;
		double EC = minEC;
		int node, tryCount=0;
		boolean success=false;
		int julianMin = Globals.currentModelTime;
		String modelDate, modelTime, modelDateTime;
		modelDate = Globals.getModelDate(julianMin);
		modelTime = Globals.getModelTime(julianMin);
		modelDateTime = modelDate+modelTime;
		
		if(ECVecHash.containsKey(modelDateTime))
		{
			ECVec = ECVecHash.get(modelDateTime);
		}
		else
		{
			ECVec = reader.readDoubleArray("QualData/" + modelDateTime);
			ECVecHash.put(modelDateTime, ECVec);
		}
		
		// Try the node indicated by useUpNode. If that doesn't work, use the opposite node. If that
		// doesn't work, return the default, arbitrarily small, EC
		while(tryCount<2 && !success)
		{
			if(useUpNode) node = c.getUpNodeId();
			else node = c.getDownNodeId();
			
			// See if this node is stored in QualData
			if(nodeECHash.containsKey(node))
			{
				// put an arbitrary lower bound on EC so the fish choose an option with EC==0.0 if there are no other 
				// options
				EC =  Math.max(minEC, ECVec[nodeECHash.get(node)]);
				success=true;
			}		
			else
			{
				useUpNode = !useUpNode;
				tryCount++;
			}
		}

		return EC;		
	}
	
	public static void destructor()
	{
		reader.close();
		System.out.println("Closed " + behaviorParameterFile);
		writer.close();
		System.out.println("Closed " + outputFilename);
	}

}
