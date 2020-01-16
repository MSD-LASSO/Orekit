package fr.cs.examples.attitude;

//3 ground stations
//Event detectors to see when all three are in view
//azimuth, elevation, and range with respect to one station
//time, az el, range then abs lat , long, altitude("elevation")
//to matlab



//possibly use circular field of view detector? or unecessary?


import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
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

public class StationDetectionModular {
	
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
   /* currently used parameters
   double mu =  3.986004415e+14;
   double a = 7083000;                     // semi major axis in meters
   double e = 0.0097707;                   // eccentricity
   double i = FastMath.toRadians(98.2+1);        // inclination
   double omega = FastMath.toRadians(64.69+0.5);  // perigee argument
   double raan = FastMath.toRadians(231.0598);   // right ascension of ascending node
   double lM =  FastMath.toRadians(0);                           // mean anomaly
   //*/
   /*//Anthony's parameters  HIGH ELEVATION ~90 deg
   double mu =  3.986004415e+14;
   double a = 6871000;                     // semi major axis in meters
   double e = 0;//0.99999999;                   // eccentricity
   double i = 1.682-0.3;//-1+1;//FastMath.toRadians(98.2+1);        // inclination
   double omega = 0.563;//-.5;//FastMath.toRadians(64.69+0.5);  // perigee argument
   double raan = 2.096;//FastMath.toRadians(231.0598);   // right ascension of ascending node
   double lM = 0.971205 ;//FastMath.toRadians(0);                           // mean anomaly
   //*/
   ///*    LOW ELEVATION ~30 deg
   double mu =  3.986004415e+14;
   double a = 6871000;                     // semi major axis in meters
   double e = 0;//0.99999999;                   // eccentricity
   double i = 1.682-.85;//-1+1;//FastMath.toRadians(98.2+1);        // inclination
   double omega = 0.2;//-.5;//FastMath.toRadians(64.69+0.5);  // perigee argument
   double raan = 1.6;//FastMath.toRadians(231.0598);   // right ascension of ascending node
   double lM = 0.971205 ;//FastMath.toRadians(0);                           // mean anomaly
   //*/
   //System.out.println(i*180/Math.PI);
   // DEFINING INERTIAL FRAME, UTC TIME SCALE AND START TIME OF SIMULATION (first time stamp)
   Frame inertialFrame = FramesFactory.getEME2000();
   TimeScale utc = TimeScalesFactory.getUTC();
   
   //set initial date is October 20th, 2019 at 23:30 
   //AbsoluteDate initialDate = new AbsoluteDate(2019, 10, 20, 23, 30, 00.000, utc);
   AbsoluteDate initialDate = new AbsoluteDate(2019, 10, 30, 0, 0, 00.000, utc);
   
   
   AbsoluteDate epoch = new AbsoluteDate(2019, 10,30,10, 0, 00.000, utc);
   
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
 /* CASE 4:  Brockport, Mees Bristol ,  Webster High School:43.204291, -77.469981. Relastic Elevations
   double stationLatitudes[]= {43.209037, 42.700192,43.204291 };
   double stationLongitudes[]=  {-77.950921,-77.408628,-77.469981};
   double stationAltitudes[]=  {165+10,  701,  147 };     //from google earth/web sources  
   // */
///* CASE 5: Pavillion, Mees Bristol, Webster
   double stationLatitudes[]= {42.871390, 42.700192,43.204291 };
   double stationLongitudes[]=  {-78.018577,-77.408628,-77.469981};
   double stationAltitudes[]=  {313+10,  701,  147 };     //from google earth/web sources  
   //*/
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
   
   // EVENT DETECTION USING ELEVATION DETECTORS
   
   double maxCheck  = 60.0;  //"maximum checking interval"
   double threshold =  0.001; //convergence threshold value
   double minElevation = 0;     //min elevation (trigger elevation)
   
   //Call on function to create a detector that will check when satellite is in view of all stations.
   BooleanDetector stationVisOverlapDetector=createBooleanDetector(stationFrames, maxCheck,
			threshold, minElevation);
   EventsLogger booleanLogger=new EventsLogger(); //creating logger to get data from detector
   kepler.addEventDetector(booleanLogger.monitorDetector(stationVisOverlapDetector));  //add event detector to propagator
   
   
   // Setting end time
   AbsoluteDate extrapDate = initialDate;
   AbsoluteDate finalDate = new AbsoluteDate(initialDate, 60000, utc);  //end time, 60000 seconds after start.
   
   
   //Propagation
 
   kepler.propagate(finalDate);  //making satellite travel to final date.
   List<EventsLogger.LoggedEvent> stationOverlap=booleanLogger.getLoggedEvents(); //getting event instances.

   // USING EPHEMERIS MODE TO GET ONE GROUND TRACK OVER THE STATIONS. Index 0 corresponds to First Entry time.
																	//    Index 1 corresponds to First Exit time.
  ///*    PR
   
   int entryIndex=0; //must be even. CHOOSING WHICH PASS TO LOOK AT
   kepler.resetInitialState(stationOverlap.get(entryIndex).getState());
   kepler.setEphemerisMode();  //ephemeris mode will give us only info from time 0 to time 1. 
   kepler.propagate(stationOverlap.get(entryIndex+1).getState().getDate());
   BoundedPropagator ephemeris = kepler.getGeneratedEphemeris();  //creating ephermis, which contains data of first ground track.
   AbsoluteDate ephemerisPropagateTime=stationOverlap.get(entryIndex).getState().getDate(); //setting initial time for propagation.
   AbsoluteDate ephemerisEndTime=stationOverlap.get(entryIndex+1).getState().getDate();
   
   //Index integer to get az, elevation, and range of satellite wrt the chosen ground station.
   int stationReference=0;
   String stationNo=Integer.toString(stationReference);
   //defining writer object to write to text file. /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   PrintWriter unWriter= new PrintWriter("/Users/mtruong/Desktop/JAVA/realElevationsCase4and5/updatedCase5WGS84LowEle.txt", "UTF-8");
 
   // Printing Orbit Info

   unWriter.printf("semimajor axis, a [meters]: %.15f%n eccentricity, e [-]: %.15f%n inclination, i [rad]: %15f%n"
   		+ "Perigee argument, omega [rad]: %15f%n Right Ascension of Ascending Node, raan [rad], %15f%n Initial Mean anomaly [rad]: %15f%n", 
   		a, e, i ,omega, raan, lM);
    
   ///* Printing station info
   unWriter.printf("%s\t%s\t%s\t%s\n","Latitude","Longitude","Altitude","Station");
   for (int b=0; b<groundStations.length;b++) {
	   System.out.println(b);
	unWriter.printf("%.15f\t%.15f\t%.15f\t%d \n",groundStations[b].getLatitude()*180/Math.PI,groundStations[b].getLongitude()*180/Math.PI,groundStations[b].getAltitude(),b);
//	unWriter.printf("%.3f %.3f %.3f \n",groundStations[b].getLatitude(),groundStations[b].getLongitude(),groundStations[b].getAltitude());
   }
   //header text 
   unWriter.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%n",
		   "Elevation wrt "+stationNo,"Azimuth wrt "+stationNo,"Range wrt "+stationNo,"Latitude","Longitude","Altitude","Time");

   
   //DEFINING WRITER OBJECT FOR GETTING POSITIONS
   //PrintWriter positionWriter= new PrintWriter("/Users/mtruong/Desktop/JAVA/testOrbitDatafakePlaceholder.txt","UTF-8");
   //positionWriter.printf("%s\t%s\t%s\t%s\n","X","Y","Z","Time");
  
   // propagation loop to get lat, lon, azi, elevation during one overhead pass in view of all stations
   while (ephemerisPropagateTime.compareTo(ephemerisEndTime)<=0) {
  	 PVCoordinates pvInert   = ephemeris.propagate(ephemerisPropagateTime).getPVCoordinates();
  	 
      Vector3D positionVectorSatellite=pvInert.getPosition();   //3D vector of satellite in earth EME2000 frame.
       
     //  System.out.println(positionVectorSatellite.getNorm());
       // Get the azimuth, elevation, and range from the .get functions. 
       // in reference to station Reference#. Position vector is in the inertial earth frame.
       double azimuth=stationFrames[stationReference].getAzimuth(positionVectorSatellite, inertialFrame,ephemerisPropagateTime );
       double elevation=stationFrames[stationReference].getElevation(positionVectorSatellite, inertialFrame, ephemerisPropagateTime);
       double range=stationFrames[stationReference].getRange(positionVectorSatellite, inertialFrame, ephemerisPropagateTime);
       
 
       String currentTimeStamp=ephemerisPropagateTime.getDate().toString();
       
   	//getAzimuth(Vector3D extPoint, Frame frame, AbsoluteDate date)
      System.out.format("Elevation: %6.15f   Azimuth: %8.15f    Range: %12.15f Time: %s %n" ,elevation*180/Math.PI, azimuth*180/Math.PI,range,ephemerisPropagateTime.getDate().toString());

  
       
       GeodeticPoint satLatLonAlt = earth.transform(pvInert.getPosition(), FramesFactory.getEME2000(),ephemerisPropagateTime);
       //System.out.println(satLatLonAlt);
       double latitude=satLatLonAlt.getLatitude()*180/Math.PI;
       double longitude=satLatLonAlt.getLongitude()*180/Math.PI;
       double altitude=satLatLonAlt.getAltitude();
       System.out.format("Latitude: %3.15f	Longitude: %3.15f	Altitude: %3.15f %n", //%n
    		   latitude,longitude,altitude);//*180/Math.PI);
       
      // System.out.println(stationFrames[stationReference].pointAtDistance(azimuth, elevation, range));
       //unWriter.printf("%3.3f	%3.3f%n", satLatLonAlt.getLatitude()*180/Math.PI,satLatLonAlt.getLongitude()*180/Math.PI);
       unWriter.printf("%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%s\t%n",
    		   elevation,azimuth,range,latitude,longitude,altitude,currentTimeStamp);
       
       double X=positionVectorSatellite.getX();
       double Y=positionVectorSatellite.getY();
       double Z=positionVectorSatellite.getZ();
    //   positionWriter.printf("%.15f\t%.15f\t%.15f\t%s\n",X,Y,Z,currentTimeStamp);
   //    System.out.format("%.15f\t%.15f\t%.15f%n",X,Y,Z);
       
       
       
       ephemerisPropagateTime = ephemerisPropagateTime.shiftedBy(15); //getting info every x seconds.
   }
   
    kepler.resetInitialState(stationOverlap.get(entryIndex+1).getState()); //getting very last state 
    Vector3D positionVectorSat=kepler.getInitialState().getPVCoordinates().getPosition(); 
    double az=stationFrames[stationReference].getAzimuth(positionVectorSat, inertialFrame,ephemerisPropagateTime );
    double el=stationFrames[stationReference].getElevation(positionVectorSat, inertialFrame, ephemerisPropagateTime);
    double ran=stationFrames[stationReference].getRange(positionVectorSat, inertialFrame, ephemerisPropagateTime);
    String currentTimeStamp=stationOverlap.get(entryIndex+1).getState().getDate().toString();
    System.out.format("Elevation: %6.15f   Azimuth: %8.15f    Range: %12.15f Time: %s %n" ,el*180/Math.PI, az*180/Math.PI,ran,currentTimeStamp);
    
	unWriter.close();
	//positionWriter.close();
	
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	
	
    // Constants
	int c = 299792458; /// m/s

	// Parameters
	double azEl_weight = 1.0; // Will be normalized later (i.e divided by the number of observations)
	double azEl_sigma = 1.0; //Estimated covariance of the range measurements, in meters

	
	AbsoluteDate odDate =new AbsoluteDate(2019, 10, 30,0,0,00.000,utc);// # Beginning of the orbit determination
		
	double collectionDuration = 1/24; //# days
	
	AbsoluteDate startCollectionDate =new  AbsoluteDate(2019, 11, 30,0,0,00.000,utc);//odDate + timedelta(days=-collectionDuration)

	//Orbit propagator parameters
	double prop_min_step = 0.001;// # s
	double prop_max_step = 300.0;// # s
	double prop_position_error = 10.0;// # m

	//# Estimator parameters
	double estimator_position_scale = 1.0;// # m
	double estimator_convergence_thres = 1e-3;
	double estimator_max_iterations = 25;
	double estimator_max_evaluations = 35;
	
	
	
	
	
	
	
	
	
	
	
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
   
