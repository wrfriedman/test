// Create an HDF5 file to hold the behavior parameters.
// Some parameters are hardcoded, others (e.g., Qual data) are read from files.
// Doug Jackson
// doug.jackson@noaa.gov

package createBehaviorHDF5;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5SimpleWriter;

public class CreateBehaviorHDF5
{
	public int numberOfChannels;
	public int numberOfReservoirs;
	public int swimCode;
	public float meanSwimSpeed;
	public float stdSwimSpeed;
	public float filterK;
	public boolean variableSwimSpeed;
	public float stageThresholdInc;
	public float stageThresholdDec;
	public float daytimeSwimProb;
	public String sunriseTime;
	public String sunsetTime;
	public int velDecisionPeriod;
	public int tideCountThr;
	public float constProbConfusion;
	public float slopeProbConfusion;
	public boolean randAssess;
	public float probAssess;
	public float initProbConfusion;
	public float holdThr;
	public int[] nodeDecisions;
	public double[][] wT_0;
	public double[][] wT_1;
	public double[][] wT_2;
	public double[][] wT_3;
	public double[] upNodes;	
	public String upNodesFile;
	public String channelParsFile;
	public String outputFilename;
	public String qualDatafileName;
	public enum releaseLocations{FREEPORT, SUTTER, STEAMBOAT, J1, J2, GEORGIANA, DCC, MOK, RIO, CVO, NA, FREEPORTTRACK};
	public releaseLocations releaseLocation;
	public int[] checkpoints;
	public boolean immortal;
	public String outputSpecFile;
	
	public static DataInputStream inputStream;
	public static IHDF5SimpleWriter writer;
	
	public static void main(String[] args)
	{		
		// Idiom to avoid using all static methods and variables
		CreateBehaviorHDF5 thisObj = new CreateBehaviorHDF5();
		
		// Read the qual filename if provided
		if(args.length>0)
		{
			thisObj.channelParsFile = args[0];
			thisObj.outputFilename = args[1];
			thisObj.qualDatafileName = args[2];
			thisObj.swimCode = Integer.parseInt(args[3]);
			thisObj.meanSwimSpeed = Float.parseFloat(args[4]);
			thisObj.stdSwimSpeed = Float.parseFloat(args[5]);
			thisObj.filterK = Float.parseFloat(args[6]);
			thisObj.variableSwimSpeed = Boolean.parseBoolean(args[7]);
			thisObj.stageThresholdInc = Float.parseFloat(args[8]);
			thisObj.stageThresholdDec = Float.parseFloat(args[9]);
			thisObj.daytimeSwimProb = Float.parseFloat(args[10]);
			thisObj.velDecisionPeriod = Integer.parseInt(args[11]);
			thisObj.tideCountThr = Integer.parseInt(args[12]);
			thisObj.holdThr = Float.parseFloat(args[13]);
			thisObj.constProbConfusion = Float.parseFloat(args[14]);
			thisObj.slopeProbConfusion = Float.parseFloat(args[15]);
			thisObj.randAssess = Boolean.parseBoolean(args[16]);
			thisObj.probAssess = Float.parseFloat(args[17]);
			thisObj.initProbConfusion = Float.parseFloat(args[18]);
			thisObj.releaseLocation = releaseLocations.valueOf(args[19]);
			thisObj.immortal = Boolean.parseBoolean(args[20]);
			thisObj.upNodesFile = args[21];
			thisObj.outputSpecFile = args[22];
		}
				
		thisObj.createFile();
	}
	
	//////////////////////////////////////////////////////////////////////
	// Constructor
	//////////////////////////////////////////////////////////////////////	
	public CreateBehaviorHDF5()
	{
		outputFilename = "/Users/doug.jackson/Documents/Shared_Windows/PTM/data/PTM_calibration/PTM_sweep_sC11_stdSwimSpeed_1_15apr15/MCMC/PTM_behavior.h5";
		
		channelParsFile = "/Users/doug.jackson/Documents/Shared_Windows/PTM/data/PTM_calibration/PTM_sweep_sC11_stdSwimSpeed_1_15apr15/MCMC/Model_C_16April2015_0510_chain_926423.csv";
		qualDatafileName = "C:/Users/doug.jackson/Documents/Sacramento_PTM_dry_run/createBehaviorHDF5/createBehaviorHDF5/qual_PTM_output_D000000041.txt";
		upNodesFile = "/Users/doug.jackson/Documents/Shared_Windows/models/git/createBehaviorHDF5/createBehaviorHDF5/channels_v8_1_2.csv";
		outputSpecFile = "";
		
		numberOfChannels = 521;
		numberOfReservoirs = 7; //6 
		swimCode = 11;
		meanSwimSpeed = 0.5f;
		stdSwimSpeed = 0.0f;
		filterK = 0.1f;
		variableSwimSpeed = true;
		stageThresholdInc = 0.5f;
		stageThresholdDec = 0.0f;
		daytimeSwimProb = 1.0f; // Probability of swimming during the day
		
		// From http://www.esrl.noaa.gov/gmd/grad/solcalc/ for Dec. 31
		sunriseTime = "0727";
		sunsetTime = "1700";
		
		// Number of hours to average the velocity over when trying to determine if
		// a particular fish has been fighting the flow and should turn around
		velDecisionPeriod = 12;
		
		// Number of tide cycles to average the velocity over when determining which
		// direction is "downstream" for a particular channel
		tideCountThr = 2;
		
		// The probability of confusion is a linear function of the log of the signal to noise ratio,
		// which is mean(vel)/std(vel) for a specific channel.
		// If randAssess==false, confusion is assessed at every junction of 3 or more waterbodies, and persists until 
		// the next junction. If randAssess==true, confusion assessment is a Bernoulli process that occurs based on
		// a coin flip with probability=probAssess every 15 minute time step. initProbConfusion determines the initial
		// probability that a fish will be confused.
		constProbConfusion = 3f;
		slopeProbConfusion = -0.25f;
		randAssess = true;
		probAssess = 0.01f;
		initProbConfusion = 0.5f;		
		
		// Upstream flow velocity at which to begin holding. E.g., 0.1f means to start
		// holding when the velocity is 0.1 m/sec upstream
		holdThr = 1.0f;
		
		nodeDecisions = new int[]{3};
		wT_0 = new double[][]{{0, 0.5, 1.0}, {1.0, 1.0, 1.0}};
		wT_1 = new double[][]{{0, 0.5, 1.0}, {1.0, 1.0, 1.0}};
		wT_2 = new double[][]{{0, 0.5, 1.0}, {1.0, 1.0, 1.0}};
		wT_3 = new double[][]{{0, 0.5, 1.0}, {1.0, 1.0, 1.0}};
		
		releaseLocation = releaseLocations.NA;
		
		immortal = false;
		
	}	

	//////////////////////////////////////////////////////////////////////
	// Instance methods
	//////////////////////////////////////////////////////////////////////
	public void createFile()
	{
		double[][] channelPars;
		HashMap<String, double[]> EC = new HashMap<String, double[]>(); 
		String line;
		int thisChannelNumber = 0;
		String modelDate, modelTime, thisEC, thisModelDateTime;
		BufferedReader bReader;
		
		Matcher matcher;
		Pattern channelNumberStartPattern = Pattern.compile("(?<=CHANNEL_)\\d+");
		Pattern modelDatePattern = Pattern.compile("[0-9]{2}[A-Z]{3}[0-9]{4}");
		Pattern modelTimePattern = Pattern.compile("(?<=\\s)\\d{4}");
		Pattern ECPattern = Pattern.compile("\\d*\\.\\d*");

		writer = initializeWriter();
		
		// Set checkpoints
		switch(releaseLocation)
		{
			case FREEPORT:
				checkpoints = new int[]{300, 270, 273, 230, 259, 309, 304, 306, 271};
				break;
			
			case SUTTER:
				checkpoints = new int[]{270, 277, 274, 285, 313, 315};
				break;
				
			case STEAMBOAT:
				checkpoints = new int[]{305, 273, 272, 275, 285, 313, 315};
				break;
			
			case J1:
				checkpoints = new int[]{271, 270, 304, 273, 305, 230, 307, 308, 259, 309};
				break;
				
			case J2:
				checkpoints = new int[]{315, 230, 259, 309, 308};
				break;
				
			case GEORGIANA:
				checkpoints = new int[]{247, 245, 258, 259, 246, 307, 308};
				break;
				
			case DCC:
				checkpoints = new int[]{307, 308, 230, 251, 236, 243, 246, 259};
				break;
				
			case MOK:
				checkpoints = new int[]{246, 247, 243, 224};
				break;
			
			case RIO:
				checkpoints = new int[]{315, 42, 317};			
				break;
			
			case CVO:
				// Read the checkpoints from the outputSpecFile
				ArrayList<Integer> tempCheckpoints = new ArrayList<Integer>();
				try
				{
					bReader = new BufferedReader(new FileReader(outputSpecFile));

					while ((line = bReader.readLine()) != null)
					{
			            System.out.println(line);
			            try
			            {
			            	tempCheckpoints.add(Integer.parseInt(line));
			            }
			            catch(NumberFormatException e)
			            {
			            	System.out.println("Checkpoint " + line + "is either already included or not valid.");
			            }
					}
					bReader.close();
				}
				catch(Exception e)
				{
					System.out.println("Error while reading file " + outputSpecFile + ":" + e.getMessage()); 
				}
				
				checkpoints = new int[tempCheckpoints.size()];
				for(int i=0; i<tempCheckpoints.size(); i++)
				{
					checkpoints[i] = tempCheckpoints.get(i);	
				}
				break;
				
			case NA:
				checkpoints = new int[]{};
				break;
			
			case FREEPORTTRACK:
				checkpoints = new int[]{300, 305, 307, 315};
				break;
			
			default:
				checkpoints = new int[]{};
			
		}
		
		// Read the channelMortalitySlopes
		channelPars = readChannelPars();
		
		// Read the node numbers
		upNodes = readUpNodes();
		
		// Read the Qual data file until the end if swimCode==7
		if(swimCode==7)
		{
			try
			{
				// Clear the previous values (just in case)
				thisEC = thisModelDateTime = "";
				bReader = new BufferedReader(new FileReader(qualDatafileName));
				while((line = bReader.readLine())!=null)
				{
					// Find the channelNumber
					matcher = channelNumberStartPattern.matcher(line);
					if(matcher.find())
					{
						thisChannelNumber = Integer.parseInt(matcher.group());
						System.out.println("Qual channel  = " + thisChannelNumber);
					}
					
					// Find modelDate
					matcher = modelDatePattern.matcher(line);
					if(matcher.find())
					{
						modelDate = matcher.group();
						
						// Find modelTime
						matcher = modelTimePattern.matcher(line);
						if(matcher.find())
						{
							modelTime = matcher.group();
							thisModelDateTime = modelDate + modelTime;
						}
						
						// Find EC
						matcher = ECPattern.matcher(line);
						if(matcher.find())
						{
							thisEC = matcher.group();
						}

						// Create a new EC array for this modelDateTime if it doesn't already exist
						if(!EC.containsKey(thisModelDateTime))
						{
							EC.put(thisModelDateTime, new double[numberOfChannels]);	
						}
						EC.get(thisModelDateTime)[thisChannelNumber-1] = Double.parseDouble(thisEC);
					}			

				}
				bReader.close();
			} catch (IOException e)
			{
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		//////////////////////////////////////////////////////////////////////
		// Write data to HDF5 file
		//////////////////////////////////////////////////////////////////////
		
		writer.writeString("channelParsFile", channelParsFile, 250);
		writer.writeString("qualDatafileName", qualDatafileName, 250);
		writer.writeInt("swimCode", swimCode);
		writer.writeFloat("meanSwimSpeed", meanSwimSpeed);
		writer.writeFloat("stdSwimSpeed", stdSwimSpeed);
		writer.writeFloat("filterK", filterK);
		writer.writeBoolean("variableSwimSpeed", variableSwimSpeed);
		writer.writeFloat("stageThresholdInc", stageThresholdInc);
		writer.writeFloat("stageThresholdDec", stageThresholdDec);
		writer.writeFloat("daytimeSwimProb", daytimeSwimProb);
		writer.writeString("sunriseTime", sunriseTime, 4);
		writer.writeString("sunsetTime", sunsetTime, 4);
		writer.writeInt("velDecisionPeriod", velDecisionPeriod);
		writer.writeInt("tideCountThr", tideCountThr);
		writer.writeFloat("holdThr", holdThr);
		writer.writeFloat("constProbConfusion", constProbConfusion);
		writer.writeFloat("slopeProbConfusion", slopeProbConfusion);
		writer.writeBoolean("randAssess", randAssess);
		writer.writeFloat("probAssess", probAssess);
		writer.writeFloat("initProbConfusion", initProbConfusion);
		writer.writeIntArray("nodeDecisions", nodeDecisions);
		writer.writeDoubleMatrix("weightsTransformation/weightsTransformation_0", wT_0);
		writer.writeDoubleMatrix("weightsTransformation/weightsTransformation_1", wT_1);
		writer.writeDoubleMatrix("weightsTransformation/weightsTransformation_2", wT_2);
		writer.writeDoubleMatrix("weightsTransformation/weightsTransformation_3", wT_3);
		writer.writeDoubleMatrix("channelPars", channelPars);
		writer.writeIntArray("checkpoints", checkpoints);
		writer.writeBoolean("immortal", immortal);
		writer.writeString("upNodesFile", upNodesFile, 250);
		
		// Write the Qual data to the HDF5 file
		writer.writeString("QualData/qualFile", qualDatafileName, qualDatafileName.length());
		writer.writeDoubleArray("QualData/upNodes", upNodes);
		for(Map.Entry<String, double[]> entry : EC.entrySet())
		{
			writer.writeDoubleArray("QualData/" + entry.getKey(), entry.getValue());
		}
		writer.close();	
			
		System.out.println("Done");
	}

	public double[] readUpNodes()
	{
		double[] uN = new double[numberOfChannels];
		int index;
		String line;
		String[] values;
		
		try
		{
			BufferedReader bReader = new BufferedReader(new FileReader(upNodesFile));
			index = 0;
			// Burn the header line
			line = bReader.readLine();
			while((line = bReader.readLine())!=null)
			{
				// Find the channelNumber
				values = line.split(",");
				uN[index] = Double.parseDouble(values[1]);
				index++;
			}
			bReader.close();
		} catch (FileNotFoundException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return uN;
	}
	public double[][] readChannelPars()
	{
		String line;
		String[] values;
		int index;
		double[][] cP = new double[numberOfChannels + numberOfReservoirs][7];
		
		try
		{
			BufferedReader bReader = new BufferedReader(new FileReader(channelParsFile));
			index = 0;
			// Burn the header line
			line = bReader.readLine();
			while((line = bReader.readLine())!=null)
			{
				// channelPars
				// 0: channel
				// 1: lambda
				// 2: omega
				// 3: meanSwimSpeed
				// 4: holdThr
				// 5: constProbConfusion
				// 6: daytimeSwimProb
				values = line.split(",");
				cP[index][0] = Integer.parseInt(values[0]);
				cP[index][1] = Double.parseDouble(values[1]);
				cP[index][2] = Double.parseDouble(values[2]);
				cP[index][3] = Double.parseDouble(values[3]);
				cP[index][4] = Double.parseDouble(values[4]);
				cP[index][5] = Double.parseDouble(values[5]);
				cP[index][6] = Double.parseDouble(values[6]);
				index++;
			}
			bReader.close();
		} catch (FileNotFoundException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return cP;
	}
	
	public IHDF5SimpleWriter initializeWriter()
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

}
