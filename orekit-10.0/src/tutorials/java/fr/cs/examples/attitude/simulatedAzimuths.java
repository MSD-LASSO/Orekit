package fr.cs.examples.attitude;
// This program looks at TLE data for 9 dif Satellites to get azimuths
//3 ground stations
//Event detectors to see when all three are in view
//azimuth, elevation, and range with respect to one station
//time, az el, range then abs lat , long, altitude("elevation")
//to matlab



//possibly use circular field of view detector? or unecessary?



//TO DO FOR ME:

// needa pick sun synchronous orbit. Choose parameters close to one form the TLE. e close to 0. inclination and semi major linked.
// then do max error.
// use paper to get inclination and semi major axis.
// chose whatever I want. 
// do monte carlo 1000 or 10000 times. get mean, std deviation of all the measurements. Error from the nominal...

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.SortedSet;
import java.util.TreeSet;

import org.hipparchus.geometry.euclidean.threed.RotationOrder;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.util.FastMath;
import org.orekit.attitudes.AttitudeProvider;
import org.orekit.attitudes.AttitudesSequence;
import org.orekit.attitudes.LofOffset;
import org.orekit.attitudes.NadirPointing;
import org.orekit.bodies.BodyShape;
import org.orekit.bodies.OneAxisEllipsoid;
import org.orekit.bodies.CelestialBodyFactory;
import org.orekit.bodies.GeodeticPoint;
import org.orekit.bodies.OneAxisEllipsoid;
import org.orekit.data.DataProvidersManager;
import org.orekit.data.DirectoryCrawler;
import org.orekit.errors.OrekitException;
import org.orekit.frames.FactoryManagedFrame;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.LOFType;
import org.orekit.frames.TopocentricFrame;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.orbits.KeplerianOrbit;
import org.orekit.orbits.Orbit;
import org.orekit.orbits.PositionAngle;
import org.orekit.propagation.BoundedPropagator;
import org.orekit.propagation.Propagator;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.analytical.EcksteinHechlerPropagator;
import org.orekit.propagation.analytical.KeplerianPropagator;
import org.orekit.propagation.analytical.tle.SGP4;
import org.orekit.propagation.analytical.tle.TLE;
import org.orekit.propagation.analytical.tle.TLEPropagator;
import org.orekit.propagation.events.BooleanDetector;
import org.orekit.propagation.events.EclipseDetector;
import org.orekit.propagation.events.ElevationDetector;
import org.orekit.propagation.events.EventDetector;
import org.orekit.propagation.events.EventsLogger;
import org.orekit.propagation.events.EventsLogger.LoggedEvent;
import org.orekit.propagation.events.handlers.ContinueOnEvent;
import org.orekit.propagation.events.handlers.EventHandler;
import org.orekit.propagation.events.handlers.RecordAndContinue;
import org.orekit.propagation.events.handlers.StopOnEvent;
import org.orekit.propagation.sampling.OrekitFixedStepHandler;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScale;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.AngularDerivativesFilter;
import org.orekit.utils.Constants;
import org.orekit.utils.IERSConventions;
import org.orekit.utils.PVCoordinates;
import org.orekit.utils.PVCoordinatesProvider;
import org.orekit.frames.TopocentricFrame;

public class simulatedAzimuths {
	
	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {

		// Need to import this orekit data to get the UTC time scale.
		//orekit data includes: ephemerides (orbit estimations of planets), drag coefficients, densities, time scales
		File home       = new File(System.getProperty("user.home"));
		File orekitData = new File(home, "orekit-data");
		if (!orekitData.exists()) {
			System.err.format(Locale.US, "Failed to find %s folder%n",
					orekitData.getAbsolutePath());
			System.err.format(Locale.US, "You need to download %s from %s, unzip it in %s and rename it 'orekit-data' for this tutorial to work%n",
					"orekit-data-master.zip", "https://gitlab.orekit.org/orekit/orekit-data/-/archive/master/orekit-data-master.zip",
					home.getAbsolutePath());
			System.exit(1);
		}
		DataProvidersManager manager = DataProvidersManager.getInstance();
		manager.addProvider(new DirectoryCrawler(orekitData));
		final SortedSet<String> output = new TreeSet<String>();


		// ORBITAL PARAMETERS OF THE LIBERTAD 1
		/*
   double mu =  3.986004415e+14;
   double a = 7083000;                     // semi major axis in meters
   double e = 0.0097707;                   // eccentricity
   double i = FastMath.toRadians(98.2);        // inclination
   double omega = FastMath.toRadians(64.69);  // perigee argument
   double raan = FastMath.toRadians(231.0598);   // right ascension of ascending node
   double lM =  FastMath.toRadians(0);                           // mean anomaly
   //*/
		///*//Anthony's parameters  HIGH ELEVATION ~90 deg
		double mu =  3.986004415e+14;
		double a = 6871000;                     // semi major axis in meters
		double e = 0;//0.99999999;                   // eccentricity
		double i = 1.682-0.3;//-1+1;//FastMath.toRadians(98.2+1);        // inclination
		double omega = 0.563;//-.5;//FastMath.toRadians(64.69+0.5);  // perigee argument
		double raan = 2.096;//FastMath.toRadians(231.0598);   // right ascension of ascending node
		double lM = 0.971205 ;//FastMath.toRadians(0);                           // mean anomaly
		//*/
		/*    LOW ELEVATION ~30 deg
   double mu =  3.986004415e+14;
   double a = 6871000;                     // semi major axis in meters
   double e = 0;//0.99999999;                   // eccentricity
   double i = 1.682-.85;//-1+1;//FastMath.toRadians(98.2+1);        // inclination
   double omega = 0.2;//-.5;//FastMath.toRadians(64.69+0.5);  // perigee argument
   double raan = 1.6;//FastMath.toRadians(231.0598);   // right ascension of ascending node
   double lM = 0.971205 ;//FastMath.toRadians(0);                           // mean anomaly
   //*/

		// DEFINING INERTIAL FRAME, UTC TIME SCALE AND START TIME OF SIMULATION (first time stamp)
		Frame inertialFrame = FramesFactory.getEME2000();
		TimeScale utc = TimeScalesFactory.getUTC();

		//set initial date as October 30th, 2019 at 0:00
		AbsoluteDate initialDate = new AbsoluteDate(2019, 10, 30, 0, 0, 00.000, utc);

		AbsoluteDate epoch = new AbsoluteDate(2019, 10,30,10, 0, 00.000, utc);

		// Setting end time
		AbsoluteDate extrapDate = initialDate;
		AbsoluteDate finalDate = new AbsoluteDate(initialDate, 3.154e+7, utc);  //end time, 60000 seconds after start.
		//1 year=354 days=3.154e+7
		//AbsoluteDate ogTime=new AbsoluteDate(2000,1,1,12,0,00.000,utc);
		//System.out.println(initialDate.durationFrom(ogTime));


		//DEFINING ORBIT of the Libertad 1 in the EME2000 frame (=inertialFrame)
		Orbit initialOrbit = new KeplerianOrbit(a, e, i, omega, raan, lM, PositionAngle.MEAN,
				inertialFrame, epoch, mu);

		//DEFINING KEPLERIAN PROPAGATOR THAT WILL ITERATE/PROPAGATE SOLUTION OF SATELLITE MOTION

		KeplerianPropagator kepler = new KeplerianPropagator(initialOrbit);
		kepler.setSlaveMode();  //slave mode is default

		// Defining a ground station positions

		//CASE 1: SINGLE ROOF: ON INSTITUTE BLOCK
		//INSTITUTE1: 43.085242, -77.679250
		//INSTITUTE2: 43.085395, -77.678947
		//INSTITUTE3: 43.085563, -77.679124

		//CASE 2: THREE RIT BUILDINGS: INSTITUTE, RIT INN, RIVERWOOD
		//INSTITUTE: 43.085346, -77.679105
		//RIT INN:   43.048300, -77.658663
		//ELLINGSON: 43.086285  -77.668015

		//CASE 3: FAR LOCATIONS: INSTITUTE, BROCKPORT, BRISTOL
		//INSTITUTE: 43.085346, -77.679105
		//BROCKPORT: 43.209037, -77.950921
		//BRISTOL:   42.700192, -77.408628


		// input lat lon in degrees in the following arrays
		//INSTITUTE HALL: 43.085346, -77.679105   Common Location
		//RIT INN &CONFERENCE CENTER~: 43.048300, -77.658663;    
		// ELLINGSON: 43.086285, -77.668015 , 154 [m];
		// THE HILL: 43.063532, -77.689936, 154 [m];
		// 
		// RIVERWOOD: 43.057017, -77.692965;
		// BROCKPORT: 43.209037, -77.950921;  
		// UR OBSERVATORY BRISTOL: 42.700192, -77.408628;
		// WEBSTER SHCROEDER: 43.204291, -77.469981

		/* CASE 1: ONE ROOF
   double stationLatitudes[]= {43.085242, 43.085395, 43.085563};
   double stationLongitudes[]=  {-77.679250, -77.678947, -77.679124};
   double stationAltitudes[]=  {0,  0,  0 };  
   // */
		/* CASE 2: AT RIT
   double stationLatitudes[]= {43.085346, 43.048300, 43.086285};
   double stationLongitudes[]=  {-77.679105, -77.658663, -77.668015};
   double stationAltitudes[]=  {0,  0,  0 };  
   // */
		/* CASE 3: FAR AWAY
   double stationLatitudes[]= {43.085346, 43.209037, 42.700192};
   double stationLongitudes[]=  {-77.679105,-77.950921,-77.408628};
   double stationAltitudes[]=  {0,  0,  0 };  
   // */
		///* CASE 4: Mess Bristol , Brockport,   Webster High School:43.204291, -77.469981
		double stationLatitudes[]= {43.209037, 42.700192,43.204291 };
		double stationLongitudes[]=  {-77.950921,-77.408628,-77.469981};
		double stationAltitudes[]=  {0,  0,  0 };  
		// */
		/* Testing...  Hamamatsu: 34.711606, 137.724915
   int mod= 5;
   double stationLatitudes[]= {43.085346+1.2*mod, 43.209037-2*mod, 42.700192};
   double stationLongitudes[]=  {-77.679105+1.2*mod,-77.950921-2*mod,-77.408628};
   double stationAltitudes[]=  {0,  0,  0 };  
   // */
		//Call on function to create ground station array from input arrays
		GeodeticPoint[] groundStations=createStations(stationLatitudes,stationLongitudes,stationAltitudes);
		// System.out.println(groundStations[1].getLatitude());


		//Create earth- model as a one axis ellipsoid- WSG84, IERS-2010.

		Frame earthFrame = FramesFactory.getITRF(IERSConventions.IERS_2010, true);
		///*
		BodyShape earth = new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
				Constants.WGS84_EARTH_FLATTENING,
				earthFrame);
		//*/
		//SPHERICAL EARTH MODEL//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//BodyShape earth= new OneAxisEllipsoid(6371000,0,earthFrame);

		//Call on function to create reference frames for the stations, based on earth position.
		TopocentricFrame[] stationFrames=createStationFrames(groundStations,earth);

		// BUNCH OF SATELLITES FROM FOR 1 YEAR TLE PROPAGATION


		//from org.orekit.propagation.analytical.tle import TLE
		//TLEPropagator class
		// FalconSAT3 Example
		//1 30776U 07006E   19322.84223699  .00002052  00000-0  57149-4 0  9995
		//2 30776  35.4334 192.8565 0001675  36.1627 323.9213 15.36801110704254
		
		
		
		
		//Order: Falconsat-3, TW-1A F0-99, F0-29, S0-50, TW-1C, Reaktor hello World   //TLE's from N2YO on 11/21/19
								//Falcon
		String[][] tleLines= { {"1 30776U 07006E   19322.84223699  .00002052  00000-0  57149-4 0  9995" ,
								"2 30776  35.4334 192.8565 0001675  36.1627 323.9213 15.36801110704254"
							   },//TW-1A
							   {"1 40928U 15051D   19325.45785450  .00002238  00000-0  55211-4 0  9996",
								"2 40928  97.1313 353.2644 0009532 221.6069 244.2001 15.41546646233161"	
							   },//F0-99
							   {"1 43937U 19003F   19325.39826597  .00001109  00000-0  48541-4 0  9995",
								"2 43937  97.2576  18.4692 0022804 240.9369 181.0320 15.23890991 46791"
							   },//F0-29
							   {"1 24278U 96046B   19324.83254322 -.00000043 +00000-0 -75556-5 0  9998",
							    "2 24278 098.5611 097.3754 0350126 197.7028 161.1655 13.53094848148532"
							   },//SO-50
							   {"1 27607U 02058C   19325.22363176 -.00000028  00000-0  16912-4 0  9991",
							   	"2 27607  64.5548 294.0712 0071991 226.1009 133.4136 14.75589945909981"
							   },//TW-1C
							   {"1 40926U 15051B   19325.09613970  .00012267  00000-0  15862-3 0  9991", 
							    "2 40926  97.1168  11.2112 0008414 142.6972 271.5568 15.59031121234349" 
							   },//REAKTOR Hello World
							   {"1 43743U 18096AA  19325.35647626  .00001127  00000-0  45252-4 0  9997", 
							   	"2 43743  97.4449  35.0810 0019078 142.6241 279.8987 15.26926366 54489"
							   }, //PRISM
							   {"1 33493U 09002B   19325.23103822  .00000786  00000-0  74081-4 0  9999", 
							   	"2 33493  98.1273 252.1078 0016114  33.0885 327.1340 14.96814275588412"   
							   }, //LILACSAT
							   {"1 40908U 15049K   19339.41053193  .00000498  00000-0  31094-4 0  9993", 
							   	"2 40908  97.4848 337.0453 0015544 341.9934 157.4020 15.13907814232487"   
							   }
							  
							 };
		//Double[][] satelliteMasses= {54,};			   
		//System.out.println(tleLines[0][1]);
		//Falconsat-3
		//String tleLine1= "1 30776U 07006E   19322.84223699  .00002052  00000-0  57149-4 0  9995" ;
		//String tleLine2="2 30776  35.4334 192.8565 0001675  36.1627 323.9213 15.36801110704254" ;
		
		ArrayList<ArrayList<Double>> azimuthArray=new ArrayList<ArrayList<Double>>(1);
		for (int k=0;k<tleLines.length;k++)
		{
			azimuthArray.add(new ArrayList<Double>(1));
			String tleLine1=tleLines[k][0];
			String tleLine2=tleLines[k][1];
					
			
			/*
			TW-1A
			F0-99
			VZULSAT 1 (CANNOT FIND ONLINE)
			F0-29
			S0-50
			TW-1C
			Reaktor Hello World
			PRISM
			LILACSAT-2
			*/
			
			TLE orekitTLE = new TLE(tleLine1, tleLine2);
	
			double mass=100; // 54 kg is FalconSat3. 100kg assumption?
			
			NadirPointing nadirPointing = new NadirPointing(inertialFrame, earth);
			SGP4 oreTLEPropagator=new SGP4(orekitTLE,nadirPointing,mass);	
	
	
			// EVENT DETECTION USING ELEVATION DETECTORS/////////////////////////////////////////////////////////
	
			double maxCheck  = 60.0;  //"maximum checking interval"
			double threshold =  0.001; //convergence threshold value
			double minElevation = 0;     //min elevation (trigger elevation)
	
			//Call on function to create a detector that will check when satellite is in view of all stations.
			BooleanDetector stationVisOverlapDetector=createBooleanDetector(stationFrames, maxCheck,
					threshold, minElevation);
			EventsLogger booleanLogger=new EventsLogger(); //creating logger to get data from detector
			oreTLEPropagator.addEventDetector(booleanLogger.monitorDetector(stationVisOverlapDetector));  //add event detector to propagator
	
			//Propagation
			////// Possibly use ElevationExtremumDetector.....
	
	
			SpacecraftState initialState= oreTLEPropagator.getInitialState();
			oreTLEPropagator.propagate(initialDate, finalDate);
			List<EventsLogger.LoggedEvent> stationOverlap=booleanLogger.getLoggedEvents(); //getting event instances.
			//System.out.println(stationOverlap.get(0).getState().getDate().toString());
			//System.out.println(initialDate.toString());
	
			//oreTLEPropagator.
			//    Index 0 corresponds to First Entry time.
			//    Index 1 corresponds to First Exit time.

			
			//System.out.println(maxElevationArray);
			//ArrayList<Double> maxElevationArray=new ArrayList<Double>(1);
			//System.out.println(stationOverlap.size());
			//System.out.println(k);
			//System.out.println(tleLine1);
			//System.out.println(stationOverlap.size());
			for (int b=0;b<stationOverlap.size()-1;b=b+2)
			{
				//System.out.println(b);
				//oreTLEPropagator.
				int entryIndex=b; //must be even. CHOOSING WHICH PASS TO LOOK AT
	
				AbsoluteDate propagateTime=stationOverlap.get(entryIndex).getState().getDate(); //setting initial time for propagation.
				AbsoluteDate endTime=stationOverlap.get(entryIndex+1).getState().getDate();
	
				//Index integer to get az, elevation, and range of satellite wrt the chosen ground station.
				int stationReference=0;
				String stationNo=Integer.toString(stationReference);
				// propagation loop to get lat, lon, azi, elevation during one overhead pass in view of all stations/////////////////////////////////////////////////////////////
	
	
				
	
				while (propagateTime.compareTo(endTime)<=0) {
					PVCoordinates pvInert   = oreTLEPropagator.getPVCoordinates(propagateTime);
					Vector3D positionVectorSatellite=pvInert.getPosition();   //3D vector of satellite in earth EME2000 frame.
	
					// System.out.println(positionVectorSatellite.getNorm());
					// Get the azimuth, elevation, and range from the .get functions. 
					// in reference to station Reference#. Position vector is in the inertial earth frame.
					double azimuth=stationFrames[stationReference].getAzimuth(positionVectorSatellite, inertialFrame,propagateTime );
					double elevation=stationFrames[stationReference].getElevation(positionVectorSatellite, inertialFrame, propagateTime);
					double range=stationFrames[stationReference].getRange(positionVectorSatellite, inertialFrame, propagateTime);
					String currentTimeStamp=propagateTime.getDate().toString();
					
					
					//System.out.format("Elevation: %6.15f   Azimuth: %8.15f    Range: %12.15f Time: %s %n" ,elevation*180/Math.PI, azimuth*180/Math.PI,range,propagateTime.getDate().toString());      
					GeodeticPoint satLatLonAlt = earth.transform(pvInert.getPosition(), FramesFactory.getEME2000(),propagateTime);
					//System.out.println(satLatLonAlt);
					double latitude=satLatLonAlt.getLatitude()*180/Math.PI;
					double longitude=satLatLonAlt.getLongitude()*180/Math.PI;
					double altitude=satLatLonAlt.getAltitude();
					// System.out.format("Latitude: %3.15f	Longitude: %3.15f	Altitude: %3.15f %n", //%n
					//		   latitude,longitude,altitude);//*180/Math.PI);
					propagateTime = propagateTime.shiftedBy(60); //getting info every x seconds.
					
					//if (elevation>(20*Math.PI/180))
					//{
						azimuthArray.get(k).add(azimuth);
			//		}
					
				}
				//System.out.println(currentMaxElevation);
				

				//System.out.println(maxElevationArray);
				/*
				if (currentMaxElevation>30)
				{
				System.out.println(currentMaxElevation);
				}
				*/
			}  //end of propagation while loop
			//*/
		}   //end of for loop for each satellite
		 PrintWriter unWriter= new PrintWriter("/Users/mtruong/Desktop/JAVA/Histogram/azimuthdebug.txt", "UTF-8");
		 unWriter.printf("%d \n", tleLines.length);
		// unWriter.printf("%s\t%s\t%s\t%s\n","Latitude","Longitude","Altitude","Station");
		for (int s=0;s<tleLines.length;s++)
		{
		//System.out.format("%s %n",maxElevationArray.get(s).toString().replaceAll("[,\\[\\]]","")); 
		
		unWriter.printf("%s %n", azimuthArray.get(s).toString().replaceAll("[,\\[\\]]",""));
		
		System.out.format("%s %n",azimuthArray.get(s).toString().replaceAll("[,\\[\\]]","")); 
		//System.out.println(maxElevationArray);
		}
		unWriter.close();
	}	
	
	
	

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// FUNCTIONS
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//function to create array of stations
	// static means you can call the function without it being attached to an object instance.
	public static GeodeticPoint[] createStations(double[] latArray, double[] lonArray, double[] altArray) {
		
		int size=latArray.length;
		
		GeodeticPoint[] stationArray= new GeodeticPoint[size];
		
		for (int i=0;i<size;i++) {
			
			stationArray[i]= new GeodeticPoint(FastMath.toRadians(latArray[i]), FastMath.toRadians(lonArray[i]),
					altArray[i]);
		}
		return stationArray;
	}
	
	//function to create array of station frames 
	public static TopocentricFrame[] createStationFrames(GeodeticPoint[] stations, BodyShape earth) {
		
		int size=stations.length;
		TopocentricFrame[] stationFrames=new TopocentricFrame[size];
		
		for (int i=0; i<size;i++) {
			String identifier="Station "+Integer.toString(i);
			stationFrames[i]=new TopocentricFrame(earth, stations[i], identifier);
		}
		return stationFrames;
	}
	
	public static BooleanDetector createBooleanDetector(TopocentricFrame[] stationFrames,double maxCheck,
			double threshold,double minElevation) {
		
		   /*
		   EventDetector sta1Visi =
		     new ElevationDetector(maxcheck, threshold, staF1).
		     withConstantElevation(elevation).    //setting the minimum elevation to check for
		     withHandler(new RecordAndContinue());
		  */
		int size=stationFrames.length;
		EventDetector[] stationVisibilityDetectors=new EventDetector[size];
		for (int i=0; i<size;i++) {
			
			stationVisibilityDetectors[i]=new ElevationDetector(maxCheck,threshold,stationFrames[i]).
					withConstantElevation(minElevation).
					withHandler(new RecordAndContinue());
		}
		
		
		BooleanDetector stationOverlapDetector=BooleanDetector.andCombine(stationVisibilityDetectors);
		return stationOverlapDetector;
	}
	
}

/*
   
   // debugging with Anthony code
   /*
   double stationLat=groundStations[stationReference].getLatitude();
   double stationLon=groundStations[stationReference].getLongitude();
   
   
   System.out.format("\tSLat%.5f\tSLon%.5f\t %n",stationLat,stationLon);
   System.out.println(groundStations[stationReference].getNorth());
   System.out.println(groundStations[stationReference].getEast());
   System.out.println(groundStations[stationReference].getZenith());
   System.out.println(groundStations[stationReference].getAltitude());

   KeplerianPropagator keplerTest = new KeplerianPropagator(initialOrbit);
   keplerTest.setSlaveMode();  //slave mode is default
   PVCoordinates pvInertTest   = keplerTest.propagate(initialDate).getPVCoordinates();	 
   Vector3D PVSat=pvInertTest.getPosition(); 
   double X=PVSat.getX();
   double Y=PVSat.getY();
   double Z=PVSat.getZ();
   //positionWriter.printf("%.15f\t%.15f\t%.15f\n",X,Y,Z);
   System.out.format("%.15f\t%.15f\t%.15f%n",X,Y,Z);
   System.out.println(PVSat.getNorm() );
   */
   
