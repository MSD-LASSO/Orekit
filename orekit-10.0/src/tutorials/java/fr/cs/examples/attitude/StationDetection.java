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

public class StationDetection {
	
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
     double mu =  3.986004415e+14;
     double a = 7083000;                     // semi major axis in meters
     double e = 0.0097707;                   // eccentricity
     double i = FastMath.toRadians(98.2);        // inclination
     double omega = FastMath.toRadians(64.69);  // perigee argument
     double raan = FastMath.toRadians(231.0598);   // right ascension of ascending node
     double lM =  FastMath.toRadians(0);                           // mean anomaly
     
     // DEFINING INERTIAL FRAME, UTC TIME SCALE AND START TIME OF SIMULATION (first time stamp)
     Frame inertialFrame = FramesFactory.getEME2000();
     TimeScale utc = TimeScalesFactory.getUTC();
     AbsoluteDate initialDate = new AbsoluteDate(2019, 10, 20, 23, 30, 00.000, utc);
     
     
     //DEFINING ORBIT of the Libertad 1 in the EME2000 frame (=inertialFrame)
     
     Orbit initialOrbit = new KeplerianOrbit(a, e, i, omega, raan, lM, PositionAngle.MEAN,
             inertialFrame, initialDate, mu);
     
     //DEFINING KEPLERIAN PROPAGATOR THAT WILL ITERATE/PROPAGATE SOLUTION OF SATELLITE MOTION
     
     KeplerianPropagator kepler = new KeplerianPropagator(initialOrbit);
     kepler.setSlaveMode();  //slave mode is default
     
     // Defining a ground station positions
     
     double longitudeStation1 = FastMath.toRadians(40.);
     double latitudeStation1  = FastMath.toRadians(43.);
     double altitudeStation1  = 0.;
     GeodeticPoint station1 = new GeodeticPoint(latitudeStation1, longitudeStation1, altitudeStation1);
     
     double longitudeStation2 = FastMath.toRadians(50.);
     double latitudeStation2  = FastMath.toRadians(53.);
     double altitudeStation2  = 0.;
     GeodeticPoint station2 = new GeodeticPoint(latitudeStation2, longitudeStation2, altitudeStation2);
     
     /*
     double longitudeStation3 = FastMath.toRadians(50.);
     double latitudeStation3  = FastMath.toRadians(43.);
     double altitudeStation3  = 0.;
     GeodeticPoint station3 = new GeodeticPoint(latitudeStation3, longitudeStation3, altitudeStation3);
     */
     //make earth- model as a one axis ellipsoid- WSG84, IERS-2010.
     
     Frame earthFrame = FramesFactory.getITRF(IERSConventions.IERS_2010, true);
     BodyShape earth = new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
                                            Constants.WGS84_EARTH_FLATTENING,
                                            earthFrame);
     // Defining frame for the station
     
     TopocentricFrame staF1 = new TopocentricFrame(earth, station1, "Station 1");
     TopocentricFrame staF2 = new TopocentricFrame(earth, station2, "Station 2");
     
     // Setting end time
     AbsoluteDate extrapDate = initialDate;
     AbsoluteDate finalDate = new AbsoluteDate(initialDate, 60000, utc);  //stop time
     
     ////// EVENT DETECTION

     double maxcheck  = 60.0;
     double threshold =  0.001;
     double elevation = 0;
     
     EventDetector sta1Visi =
       new ElevationDetector(maxcheck, threshold, staF1).
       withConstantElevation(elevation).    //setting the minimum elevation to check for
       withHandler(new RecordAndContinue());
    
     EventsLogger logger1=new EventsLogger();   //Event logger to log when event is detected.
     
     EventDetector sta2Visi =
    	       new ElevationDetector(maxcheck, threshold, staF2).
    	       withConstantElevation(elevation).    //setting the minimum elevation to check for
    	       withHandler(new RecordAndContinue());
    	    
     EventsLogger logger2=new EventsLogger();   //Event logger to log when event is detected.
     
     kepler.addEventDetector(logger1.monitorDetector(sta1Visi)); 
     kepler.addEventDetector(logger2.monitorDetector(sta2Visi)); 
     
     BooleanDetector bothStations=BooleanDetector.andCombine(sta1Visi,sta2Visi);
     EventsLogger booleanLogger=new EventsLogger();
     kepler.addEventDetector(booleanLogger.monitorDetector(bothStations));
     
     kepler.propagate(finalDate);

     List<EventsLogger.LoggedEvent> allEventsStation1=logger1.getLoggedEvents();
     List<EventsLogger.LoggedEvent> allEventsStation2=logger2.getLoggedEvents();
     List<EventsLogger.LoggedEvent> stationOverlap=booleanLogger.getLoggedEvents();
     
     kepler.resetInitialState(stationOverlap.get(0).getState());
     kepler.setEphemerisMode();
     kepler.propagate(stationOverlap.get(1).getState().getDate());
     BoundedPropagator ephemeris = kepler.getGeneratedEphemeris();
     AbsoluteDate ephemerisPropagateTime=stationOverlap.get(0).getState().getDate();
     AbsoluteDate ephemerisEndTime=stationOverlap.get(1).getState().getDate();
     
     while (ephemerisPropagateTime.compareTo(ephemerisEndTime)<=0) {
    	 PVCoordinates pvInert   = ephemeris.propagate(ephemerisPropagateTime).getPVCoordinates();
         Vector3D positionVectorSatellite=pvInert.getPosition();
         
         // Get the azimuth, elevation, and range from the .get functions. 
         // in the staF (station frame). Get the position vector in the inertial frame
         double azimuth=staF1.getAzimuth(positionVectorSatellite, inertialFrame,ephemerisPropagateTime );
         double elevationSat=staF1.getElevation(positionVectorSatellite, inertialFrame, ephemerisPropagateTime);
         double range=staF1.getRange(positionVectorSatellite, inertialFrame, ephemerisPropagateTime);
         
         
     	//getAzimuth(Vector3D extPoint, Frame frame, AbsoluteDate date)
         System.out.format("Elevation: %6.3f   Azimuth: %8.3f    Range: %12.3f Time: %s %n" ,elevationSat*180/Math.PI, azimuth*180/Math.PI,range,ephemerisPropagateTime.getDate());
         
         
         GeodeticPoint satLatLonAlt = earth.transform(pvInert.getPosition(), FramesFactory.getEME2000(), extrapDate);
         //System.out.println(satLatLonAlt);
         
         System.out.format("%3.3f	%3.3f%n",satLatLonAlt.getLatitude()*180/Math.PI,satLatLonAlt.getLongitude()*180/Math.PI);
         //unWriter.printf("%3.3f	%3.3f%n", satLatLonAlt.getLatitude()*180/Math.PI,satLatLonAlt.getLongitude()*180/Math.PI);
         
         ephemerisPropagateTime = ephemerisPropagateTime.shiftedBy(60);
     }
     //kepler.get
     
     /*
     for(int b=0;b<allEventsStation1.size();b++){  
    	 System.out.println(b);  
    	 String station1Date=allEventsStation1.get(b).getState().getDate().toString();
    	 String station2Date=allEventsStation2.get(b).getState().getDate().toString();
    	 String overlapDate=stationOverlap.get(b).getState().getDate().toString();
    	 System.out.format("Station 1: %s%nStation 2: %s%nOverlap:   %s%n",station1Date,station2Date,overlapDate);
    	 //on entry, the later date is chosen for both in view start time
    	 //on exit, the earlier date is chosen for both in view end time
    	 
     }
     //*/
     /*
     while (extrapDate.compareTo(finalDate) <= 0)  {

         //Position and velocity coordinates of our satellite
         PVCoordinates pvInert   = kepler.propagate(extrapDate).getPVCoordinates();
         Vector3D positionVectorSatellite=pvInert.getPosition();
         
         
         // POSITION AND VELOCITY OF KEPLER OBJECT IN THE STATION FRAME.
         PVCoordinates pvStation = inertialFrame.getTransformTo(staF1, extrapDate).transformPVCoordinates(pvInert);
         
         // GET POSITION VECTOR OF SATELLITE WRT STATION
         Vector3D positionVectorWRTStation=pvStation.getPosition();
         
         
         // Get the azimuth, elevation, and range from the .get functions. 
         // in the staF (station frame). Get the position vector in the inertial frame
         double azimuth=staF1.getAzimuth(positionVectorSatellite, inertialFrame, extrapDate);
         double elevationSat=staF1.getElevation(positionVectorSatellite, inertialFrame, extrapDate);
         double range=staF1.getRange(positionVectorSatellite, inertialFrame, extrapDate);
         
     	//getAzimuth(Vector3D extPoint, Frame frame, AbsoluteDate date)
         System.out.format("Elevation: %3.3f   Azimuth: %3.3f    Range: %3.3f Time: %s %n" ,elevationSat*180/Math.PI, azimuth*180/Math.PI,range,extrapDate.getDate());
         extrapDate = extrapDate.shiftedBy(60);
     }
	//*/
	
	
	}	
	
	
	
	//public  {
	//} function to create array of stations from lat and long array
	// function to create event detectors for the stations
	// function to boolean stuff
	//print results
	
}
