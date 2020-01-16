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
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.SortedSet;
import java.util.TreeSet;
import org.hipparchus.geometry.euclidean.threed.RotationOrder;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.QRDecomposer;
import org.hipparchus.optim.nonlinear.vector.leastsquares.GaussNewtonOptimizer;
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
import org.orekit.estimation.leastsquares.BatchLSEstimator;
import org.orekit.estimation.measurements.AngularAzEl;
import org.orekit.estimation.measurements.GroundStation;
import org.orekit.estimation.measurements.ObservableSatellite;
import org.orekit.frames.FactoryManagedFrame;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.LOFType;
import org.orekit.frames.TopocentricFrame;
import org.orekit.models.earth.ReferenceEllipsoid;
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
import org.orekit.propagation.conversion.DormandPrince853IntegratorBuilder;
import org.orekit.propagation.conversion.NumericalPropagatorBuilder;
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
import org.orekit.propagation.integration.AbstractIntegratedPropagator;
import org.orekit.propagation.sampling.OrekitFixedStepHandler;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScale;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.AngularDerivativesFilter;
import org.orekit.utils.Constants;
import org.orekit.utils.IERSConventions;
import org.orekit.utils.PVCoordinates;
import org.orekit.utils.PVCoordinatesProvider;
import org.orekit.utils.TimeStampedPVCoordinates;
import org.orekit.frames.TopocentricFrame;

public class simulatedOrbitDetermination {
	
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
		   									0,    //let flattening = 0 (sphere)//Constants.WGS84_EARTH_FLATTENING,
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
   
	// FalconSAT Example
	//1 30776U 07006E   19322.84223699  .00002052  00000-0  57149-4 0  9995
	//2 30776  35.4334 192.8565 0001675  36.1627 323.9213 15.36801110704254
   
	String tleLine1= "1 30776U 07006E   19322.84223699  .00002052  00000-0  57149-4 0  9995" ;
	String tleLine2="2 30776  35.4334 192.8565 0001675  36.1627 323.9213 15.36801110704254" ;
	TLE orekitTLE = new TLE(tleLine1, tleLine2);
	
	double mass=54; //kg
	NadirPointing nadirPointing = new NadirPointing(inertialFrame, earth);
	SGP4 oreTLEPropagator=new SGP4(orekitTLE,nadirPointing,mass);
   
   
   //Call on function to create a detector that will check when satellite is in view of all stations.
   BooleanDetector stationVisOverlapDetector=createBooleanDetector(stationFrames, maxCheck,
			threshold, minElevation);
   EventsLogger booleanLogger=new EventsLogger(); //creating logger to get data from detector
   oreTLEPropagator.addEventDetector(booleanLogger.monitorDetector(stationVisOverlapDetector));  //add event detector to propagator
   
   // Setting end time
   AbsoluteDate extrapDate = initialDate;
   AbsoluteDate finalDate = new AbsoluteDate(initialDate, 60000, utc);  //end time, 60000 seconds after start.
   
   
   //Propagation
   SpacecraftState initialState= oreTLEPropagator.getInitialState();
   oreTLEPropagator.propagate(initialDate,finalDate);  //making satellite travel to final date.
   List<EventsLogger.LoggedEvent> stationOverlap=booleanLogger.getLoggedEvents(); //getting event instances.

   // USING EPHEMERIS MODE TO GET ONE GROUND TRACK OVER THE STATIONS. Index 0 corresponds to First Entry time.
																	//    Index 1 corresponds to First Exit time.
  ///*    PR
   
   int entryIndex=4; //must be even. CHOOSING WHICH PASS TO LOOK AT
   /*
   oreTLEPropagator.resetInitialState(stationOverlap.get(entryIndex).getState());
   kepler.setEphemerisMode();  //ephemeris mode will give us only info from time 0 to time 1. 
   kepler.propagate(stationOverlap.get(entryIndex+1).getState().getDate());
   BoundedPropagator ephemeris = kepler.getGeneratedEphemeris();  //creating ephermis, which contains data of first ground track.
   AbsoluteDate ephemerisPropagateTime=stationOverlap.get(entryIndex).getState().getDate(); //setting initial time for propagation.
   AbsoluteDate ephemerisEndTime=stationOverlap.get(entryIndex+1).getState().getDate();
   */
   //Index integer to get az, elevation, and range of satellite wrt the chosen ground station.
   int stationReference=0;
   String stationNo=Integer.toString(stationReference);

   	AbsoluteDate propagateTime=stationOverlap.get(entryIndex).getState().getDate(); //setting initial time for propagation.
	AbsoluteDate endTime=stationOverlap.get(entryIndex+1).getState().getDate();

  
   // propagation loop to get lat, lon, azi, elevation during one overhead pass in view of all stations/////////////////////////////////////////////////////////////
   
   ArrayList<Double> azArray=new ArrayList<Double>(1);
   ArrayList<Double> elArray=new ArrayList<Double>(1);
   ArrayList<AbsoluteDate> timeStamps= new ArrayList<AbsoluteDate>(1);
   
   while (propagateTime.compareTo(endTime)<=0) {
  	 PVCoordinates pvInert   = oreTLEPropagator.getPVCoordinates(propagateTime);
  	 
      Vector3D positionVectorSatellite=pvInert.getPosition();   //3D vector of satellite in earth EME2000 frame.
       
     //  System.out.println(positionVectorSatellite.getNorm());
       // Get the azimuth, elevation, and range from the .get functions. 
       // in reference to station Reference#. Position vector is in the inertial earth frame.
       double azimuth=stationFrames[stationReference].getAzimuth(positionVectorSatellite, inertialFrame,propagateTime );
       double elevation=stationFrames[stationReference].getElevation(positionVectorSatellite, inertialFrame, propagateTime);
       double range=stationFrames[stationReference].getRange(positionVectorSatellite, inertialFrame, propagateTime);
       
 
       String currentTimeStamp=propagateTime.getDate().toString();
       
   	  //getAzimuth(Vector3D extPoint, Frame frame, AbsoluteDate date)
      System.out.format("Elevation: %6.15f   Azimuth: %8.15f    Range: %12.15f Time: %s %n" ,elevation*180/Math.PI, azimuth*180/Math.PI,range,propagateTime.getDate().toString());

      
       
       GeodeticPoint satLatLonAlt = earth.transform(pvInert.getPosition(), FramesFactory.getEME2000(),propagateTime);
       //System.out.println(satLatLonAlt);
       double latitude=satLatLonAlt.getLatitude()*180/Math.PI;
       double longitude=satLatLonAlt.getLongitude()*180/Math.PI;
       double altitude=satLatLonAlt.getAltitude();
       System.out.format("Latitude: %3.15f	Longitude: %3.15f	Altitude: %3.15f %n", //%n
    		   latitude,longitude,altitude);//*180/Math.PI);

       double X=positionVectorSatellite.getX();
       double Y=positionVectorSatellite.getY();
       double Z=positionVectorSatellite.getZ();
    //   positionWriter.printf("%.15f\t%.15f\t%.15f\t%s\n",X,Y,Z,currentTimeStamp);
   //    System.out.format("%.15f\t%.15f\t%.15f%n",X,Y,Z);
       
       // adding az el and timestamps of flyby to arrays
       azArray.add(azimuth);
       elArray.add(elevation);
       timeStamps.add(propagateTime);
       
       propagateTime = propagateTime.shiftedBy(30); //getting info every x seconds.
   }
   

    //adding final az el and timestamp of flyby
  
	
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Orbit Determination Attempt
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	

    // Constants
	int c = 299792458; /// m/s

	// Parameters
	double azEl_weight = 1; // Will be normalized later (i.e divided by the number of observations)
	double azEl_sigma =  1; //Estimated covariance of the range measurements, in meters

	
	//Orbit propagator parameters
	double prop_min_step = 0.001;// # s
	double prop_max_step = 300.0;// # s
	double prop_position_error = 10.0;// # m

	//# Estimator parameters
	double estimator_position_scale = 1.0;// # m
	double estimator_convergence_thres = 1e-1;
	int estimator_max_iterations = 25;
	int estimator_max_evaluations = 100;
	
	
	SpacecraftState tleInitialState = oreTLEPropagator.getInitialState();
	
	AbsoluteDate tleEpoch = tleInitialState.getDate();
	Orbit tleOrbit_TEME = tleInitialState.getOrbit();
	TimeStampedPVCoordinates tlePV_ECI = tleOrbit_TEME.getPVCoordinates(inertialFrame);
	
	
	//extra ellipsoid definition just to get GM- attraction coefficient
	ReferenceEllipsoid wgs84Ellipsoid = ReferenceEllipsoid.getWgs84(earthFrame);
	
	//from org.orekit.orbits import CartesianOrbit
	CartesianOrbit tleOrbit_ECI =new CartesianOrbit(tlePV_ECI, inertialFrame, wgs84Ellipsoid.getGM());
	
	//from org.orekit.propagation.conversion import DormandPrince853IntegratorBuilder
	DormandPrince853IntegratorBuilder integratorBuilder = new DormandPrince853IntegratorBuilder(prop_min_step, prop_max_step, prop_position_error);

	//from org.orekit.propagation.conversion import NumericalPropagatorBuilder
	//from org.orekit.orbits import PositionAngle
	NumericalPropagatorBuilder propagatorBuilder = new NumericalPropagatorBuilder(tleOrbit_ECI,integratorBuilder, PositionAngle.MEAN, estimator_position_scale);
	propagatorBuilder.setMass(mass);
	propagatorBuilder.setAttitudeProvider(nadirPointing);
		
	//from org.hipparchus.linear import QRDecomposer
	QRDecomposer matrixDecomposer =new QRDecomposer(1e-11);
	//from org.hipparchus.optim.nonlinear.vector.leastsquares import GaussNewtonOptimizer
	GaussNewtonOptimizer optimizer =new GaussNewtonOptimizer(matrixDecomposer,false);

	//from org.orekit.estimation.leastsquares import BatchLSEstimator
	BatchLSEstimator estimator =new BatchLSEstimator(optimizer, propagatorBuilder);
	estimator.setParametersConvergenceThreshold(estimator_convergence_thres);
	estimator.setMaxIterations(estimator_max_iterations);
	estimator.setMaxEvaluations(estimator_max_evaluations);
	
	//azArray,elArray,timestamps
	//estimator.addMeasurement(azArray);
	GroundStation station0=new GroundStation(stationFrames[0]);
	double[] weight=new double[] {azEl_weight, azEl_weight};
	double[] sigma= new double[] {azEl_sigma, azEl_sigma};
	
	
	//what is "propagator index" for the ObservedSatellite constructor??????????
	
	ObservableSatellite satellite=new ObservableSatellite(0);
	for (int s=0;s<azArray.size();s++)
	{
		
		double[] currentAzEl=new double[]{azArray.get(s)*180/Math.PI,elArray.get(s)*180/Math.PI}; 
		AngularAzEl measurement=new AngularAzEl(station0,timeStamps.get(s),currentAzEl,weight,sigma,satellite);
		estimator.addMeasurement(measurement);
		System.out.println(currentAzEl[1]);
	}
	
	AbstractIntegratedPropagator[] estimatedPropagatorArray = estimator.estimate();
	AbstractIntegratedPropagator estimatedPropagator = estimatedPropagatorArray[0];
	SpacecraftState estimatedInitialState = estimatedPropagator.getInitialState();
	AbsoluteDate actualOdDate = estimatedInitialState.getDate();
    Orbit estimatedOrbit_init = estimatedInitialState.getOrbit();
	
    System.out.println("done");
	//estimatedOrbit_init.

	
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
   
