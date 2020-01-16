package fr.cs.examples.attitude;

//needa pick sun synchronous orbit. Choose parameters close to one form the TLE. e close to 0. inclination and semi major linked.
//then do max error.
//use paper to get inclination and semi major axis.
//chose whatever I want. 
//do monte carlo 1000 or 10000 times. get mean, std deviation of all the measurements. Error from the nominal...



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
import java.util.Random;

public class simODtestV2 {
	
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
   

   // DEFINING INERTIAL FRAME, UTC TIME SCALE AND START TIME OF SIMULATION (first time stamp)
   Frame inertialFrame = FramesFactory.getEME2000();
   TimeScale utc = TimeScalesFactory.getUTC();
   
   //set initial date as October 30th, 2019 at 0:00
   AbsoluteDate initialDate = new AbsoluteDate(2019, 10, 30, 0, 0, 00.000, utc);
   
   
   AbsoluteDate epoch = new AbsoluteDate(2019, 10,30,10, 0, 00.000, utc);
   
   //KeplerianOrbit
   ///*//Anthony's parameters  HIGH ELEVATION ~90 deg
   double mu =  3.986004415e+14;
   double a = 6871000;                     // semi major axis in meters
   double e = 0.01;//0.99999999;                   // eccentricity
   double i = 1.54;//-1+1;//FastMath.toRadians(98.2+1);        // inclination
   double omega = 0.5;//-.5;//FastMath.toRadians(64.69+0.5);  // perigee argument
   double raan = 2.3;//FastMath.toRadians(231.0598);   // right ascension of ascending node
   double lM = 0.971205 ;//FastMath.toRadians(0);                           // mean anomaly
   //*/
   
   KeplerianOrbit keplerOrbit = new KeplerianOrbit(a, e, i, omega, raan, lM, PositionAngle.MEAN,
           inertialFrame, epoch, mu);
  // keplerOrbit.
   //DEFINING KEPLERIAN PROPAGATOR THAT WILL ITERATE/PROPAGATE SOLUTION OF SATELLITE MOTION
  
   KeplerianPropagator kepler = new KeplerianPropagator(keplerOrbit);
   kepler.setSlaveMode();  //slave mode is default
   
   
   
   
   
 
 ///* Ground Station Triangle
   double stationLatitudes[]= {43.209037, 42.700192,43.204291 };
   double stationLongitudes[]=  {-77.950921,-77.408628,-77.469981};
   double stationAltitudes[]=  {0,  0,  0 };  

   //Call on function to create ground station array from input arrays
   GeodeticPoint[] groundStations=createStations(stationLatitudes,stationLongitudes,stationAltitudes);
  // System.out.println(groundStations[1].getLatitude());
   

   //Create earth- model as a one axis ellipsoid- WSG84, IERS-2010.
   
   Frame earthFrame = FramesFactory.getITRF(IERSConventions.IERS_2010, true);
  ///*
   BodyShape earth = new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
		   									0,    //let flattening = 0 (sphere)//Constants.WGS84_EARTH_FLATTENING,
                                          earthFrame);

   //Call on function to create reference frames for the stations, based on earth position.
   TopocentricFrame[] stationFrames=createStationFrames(groundStations,earth);
   
   // EVENT DETECTION USING ELEVATION DETECTORS
   
   double maxCheck  = 60.0;  //"maximum checking interval"
   double threshold =  0.001; //convergence threshold value
   double minElevation = 0;     //min elevation (trigger elevation)
   
	// FalconSAT-3 Example TLE
	//1 30776U 07006E   19322.84223699  .00002052  00000-0  57149-4 0  9995
	//2 30776  35.4334 192.8565 0001675  36.1627 323.9213 15.36801110704254
   
	//String tleLine1= "1 30776U 07006E   19322.84223699  .00002052  00000-0  57149-4 0  9995" ;
	//String tleLine2="2 30776  35.4334 192.8565 0001675  36.1627 323.9213 15.36801110704254" ;
	//TLE orekitTLE = new TLE(tleLine1, tleLine2);
	
	//double mass=54; //kg
	//NadirPointing nadirPointing = new NadirPointing(inertialFrame, earth);
	//SGP4 oreTLEPropagator=new SGP4(orekitTLE,nadirPointing,mass);
   

   //Call on function to create a detector that will check when satellite is in view of all stations.
   BooleanDetector stationVisOverlapDetector=createBooleanDetector(stationFrames, maxCheck,
			threshold, minElevation);
   EventsLogger booleanLogger=new EventsLogger(); //creating logger to get data from detector
   kepler.addEventDetector(booleanLogger.monitorDetector(stationVisOverlapDetector));  //add event detector to propagator
   
   // Setting end time
   AbsoluteDate extrapDate = initialDate;
   AbsoluteDate finalDate = new AbsoluteDate(initialDate, 60000, utc);  //end time, 60000 seconds after start.
   
   
   //Propagation
   SpacecraftState initialState= kepler.getInitialState();
   kepler.propagate(initialDate,finalDate);  //making satellite travel to final date.
   List<EventsLogger.LoggedEvent> stationOverlap=booleanLogger.getLoggedEvents(); //getting event instances.

   ArrayList<Double> azArray=new ArrayList<Double>(1);
   ArrayList<Double> elArray=new ArrayList<Double>(1);
   ArrayList<AbsoluteDate> timeStamps= new ArrayList<AbsoluteDate>(1);
   
   
for (int b=4;b<=4;b=b+2)
{
  int entryIndex=b; //must be even. CHOOSING WHICH PASS TO LOOK AT
   kepler.resetInitialState(stationOverlap.get(entryIndex).getState());
   kepler.setEphemerisMode();  //ephemeris mode will give us only info from time 0 to time 1. 
   kepler.propagate(stationOverlap.get(entryIndex+1).getState().getDate());
   BoundedPropagator ephemeris = kepler.getGeneratedEphemeris();  //creating ephermis, which contains data of first ground track.
   AbsoluteDate ephemerisPropagateTime=stationOverlap.get(entryIndex).getState().getDate(); //setting initial time for propagation.
   AbsoluteDate ephemerisEndTime=stationOverlap.get(entryIndex+1).getState().getDate();
   
   //Index integer to get az, elevation, and range of satellite wrt the chosen ground station.
   int stationReference=0;
   String stationNo=Integer.toString(stationReference);

   // propagation loop to get lat, lon, azi, elevation during one overhead pass in view of all stations/////////////////////////////////////////////////////////////
   

   
   Random rand=new Random();
   int dataCounter=0;
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
       
       //printing for debugging
      //System.out.format("Elevation: %6.15f   Azimuth: %8.15f    Range: %12.15f Time: %s %n" ,elevation*180/Math.PI, azimuth*180/Math.PI,range,ephemerisPropagateTime.getDate().toString());
       
       GeodeticPoint satLatLonAlt = earth.transform(pvInert.getPosition(), FramesFactory.getEME2000(),ephemerisPropagateTime);
       //System.out.println(satLatLonAlt);
       double latitude=satLatLonAlt.getLatitude()*180/Math.PI;
       double longitude=satLatLonAlt.getLongitude()*180/Math.PI;
       double altitude=satLatLonAlt.getAltitude();
     //  System.out.format("Latitude: %3.15f	Longitude: %3.15f	Altitude: %3.15f %n", //%n
    //		   latitude,longitude,altitude);   //printing for debugging

       // adding az el and timestamps of flyby to arrays
       double rand_dub1=rand.nextGaussian()*1*Math.PI/180;
       double rand_dub2=rand.nextGaussian()*1*Math.PI/180;
       //System.out.println(rand_dub1*180/Math.PI);
       azArray.add(azimuth);//+rand_dub1);
       elArray.add(elevation);//+rand_dub2);
       timeStamps.add(ephemerisPropagateTime);
       
       ephemerisPropagateTime = ephemerisPropagateTime.shiftedBy(5); //getting info every x seconds.
       dataCounter++;
   }
}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Orbit Determination Attempt
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	

    // Constants
	int c = 299792458; /// m/s

	// Parameters
	double azEl_weight = .1; // Will be normalized later (i.e divided by the number of observations)
	double azEl_sigma =  0.03; //Estimated covariance of the range measurements, in meters

	
	//Orbit propagator parameters
	double prop_min_step = 0.001;// # s
	double prop_max_step = 300.0;// # s
	double prop_position_error = 10.0;// # m

	//# Estimator parameters
	double estimator_position_scale = 1.0;// # m
	double estimator_convergence_thres = 1e-2;
	int estimator_max_iterations = 10000;
	int estimator_max_evaluations = 10000;
	
	
	//SpacecraftState tleInitialState = oreTLEPropagator.getInitialState();
	
	//AbsoluteDate tleEpoch = tleInitialState.getDate();
	//Orbit tleOrbit_TEME = tleInitialState.getOrbit();
	//TimeStampedPVCoordinates tlePV_ECI = tleOrbit_TEME.getPVCoordinates(inertialFrame);
	
	
	//extra ellipsoid definition just to get GM- attraction coefficient
	ReferenceEllipsoid wgs84Ellipsoid = ReferenceEllipsoid.getWgs84(earthFrame);
	
	//from org.orekit.orbits import CartesianOrbit
	//CartesianOrbit tleOrbit_ECI =new CartesianOrbit(tlePV_ECI, inertialFrame, wgs84Ellipsoid.getGM());
	
	//from org.orekit.propagation.conversion import DormandPrince853IntegratorBuilder
	DormandPrince853IntegratorBuilder integratorBuilder = new DormandPrince853IntegratorBuilder(prop_min_step, prop_max_step, prop_position_error);

	//from org.orekit.propagation.conversion import NumericalPropagatorBuilder
	//from org.orekit.orbits import PositionAngle
	
	double mass=54; //kg
	NadirPointing nadirPointing = new NadirPointing(inertialFrame, earth);
	
	NumericalPropagatorBuilder propagatorBuilder = new NumericalPropagatorBuilder(keplerOrbit,integratorBuilder, PositionAngle.MEAN, estimator_position_scale);
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
	
	
	//ObservableSatellite Index 0 is OK when only one satellite is being observed.
	
	ObservableSatellite satellite=new ObservableSatellite(0);
	for (int s=0;s<azArray.size();s++)
	{
		
		//double[] currentAzEl=new double[]{azArray.get(s)*180/Math.PI,elArray.get(s)*180/Math.PI}; 
		//double[] currentAzEl=new double[]{azArray.get(s),elArray.get(s)}; 
		double[] currentAzEl=new double[]{azArray.get(s),elArray.get(s)}; 
		
		AngularAzEl measurement=new AngularAzEl(station0,timeStamps.get(s),currentAzEl,sigma,weight,satellite);
		estimator.addMeasurement(measurement);
		System.out.println(currentAzEl[1]);
	}
	
	AbstractIntegratedPropagator[] estimatedPropagatorArray = estimator.estimate();
	AbstractIntegratedPropagator estimatedPropagator = estimatedPropagatorArray[0];
	SpacecraftState estimatedInitialState = estimatedPropagator.getInitialState();
	AbsoluteDate actualOdDate = estimatedInitialState.getDate();
    Orbit estimatedOrbit_init = estimatedInitialState.getOrbit();
	
    System.out.println("done");

    System.out.format("eccentricity: %.6f \n",estimatedOrbit_init.getE());
    System.out.format("semimajor axis: %.6f \n",estimatedOrbit_init.getA());
    System.out.format("Inclination: %.6f \n",estimatedOrbit_init.getI());
    
    //System.out.println(estimatedOrbit_init.get;

    KeplerianOrbit temp=(KeplerianOrbit) estimatedOrbit_init; 
    System.out.format("Mean Anomaly: %.6f \n",temp.getMeanAnomaly());
    System.out.format("Right Ascension: %.6f \n",temp.getRightAscensionOfAscendingNode());
    System.out.format("Argument of Perigee: %.6f \n",temp.getPerigeeArgument());
    
    /*//Anthony's parameters  HIGH ELEVATION ~90 deg
  double mu =  3.986004415e+14;
   double a = 6871000;                     // semi major axis in meters
   double e = 0.01;//0.99999999;                   // eccentricity
   double i = 1.54;//-1+1;//FastMath.toRadians(98.2+1);        // inclination
   double omega = 0.5;//-.5;//FastMath.toRadians(64.69+0.5);  // perigee argument
   double raan = 2.3;//FastMath.toRadians(231.0598);   // right ascension of ascending node
   double lM = 0.971205 
    //*/
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


