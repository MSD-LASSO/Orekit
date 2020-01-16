package fr.cs.examples.attitude;

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
import org.orekit.propagation.Propagator;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.analytical.EcksteinHechlerPropagator;
import org.orekit.propagation.analytical.KeplerianPropagator;
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


public class TESTCLASS {
	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
		// TODO Auto-generated method stub
		
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
         
		
         // DEFINING INERTIAL FRAME AND FIRST TIME STAMP 
         Frame inertialFrame = FramesFactory.getEME2000();
         TimeScale utc = TimeScalesFactory.getUTC();
         AbsoluteDate initialDate = new AbsoluteDate(2019, 10, 20, 23, 30, 00.000, utc);
         
         // ORBITAL PARAMETERS OF THE LIBERTAD 1
         double mu =  3.986004415e+14;
         double a = 7083000;                     // semi major axis in meters
         double e = 0.0097707;                   // eccentricity
         double i = FastMath.toRadians(98.2);        // inclination
         double omega = FastMath.toRadians(64.69);  // perigee argument
         double raan = FastMath.toRadians(231.0598);   // right ascension of ascending node
         double lM = 0;                           // mean anomaly
         
         // Defining orbit from the libertad paramters
         
         //orbit is defined in the EME2000 frame= inertialFrame.
         Orbit initialOrbit = new KeplerianOrbit(a, e, i, omega, raan, lM, PositionAngle.MEAN,
                 inertialFrame, initialDate, mu);
         
         //using keplerian algorithm/integrator to propogate
         KeplerianPropagator kepler = new KeplerianPropagator(initialOrbit);
         kepler.setSlaveMode();  //slave mode is default
         
         // Defining a ground station position
         
         double longitude = FastMath.toRadians(43.);
         double latitude  = FastMath.toRadians(77.);
         double altitude  = 0.;
         GeodeticPoint station = new GeodeticPoint(latitude, longitude, altitude);
         
         //make earth- model as a one axis ellipsoid- WSG84, IERS-2010.
         
         Frame earthFrame = FramesFactory.getITRF(IERSConventions.IERS_2010, true);
         BodyShape earth = new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
                                                Constants.WGS84_EARTH_FLATTENING,
                                                earthFrame);
         // Defining frame for the station
         
         TopocentricFrame staF = new TopocentricFrame(earth, station, "station");
         
         // "time stamp" time to propagate to
         AbsoluteDate extrapDate = initialDate;
         AbsoluteDate finalDate = new AbsoluteDate(initialDate, 60000, utc);  //stop time

         // Event Detector 
         
         double maxcheck  = 60.0;
         double threshold =  0.001;
         double elevation = 0;
         EventDetector sta1Visi =
           new ElevationDetector(maxcheck, threshold, staF).
           withConstantElevation(elevation).    //setting the minimum elevation to check for
           withHandler(new RecordAndContinue());
        
         EventsLogger logger=new EventsLogger();   //Event logger to log when event is detected.
         kepler.addEventDetector(logger.monitorDetector(sta1Visi)); 
         kepler.propagate(finalDate);
         
         List<EventsLogger.LoggedEvent> allEvents=logger.getLoggedEvents();
         
         System.out.println(allEvents.get(1).getState().getDate().toString());
        // System.out.println(allEvents.get(2).getState().getDate().toString());
        // System.out.println(allEvents.get(3).getState().getDate().toString());
        // System.out.println(allEvents.get(4).getState().getDate().toString());

         
         
         
         KeplerianPropagator kepler2 = new KeplerianPropagator(initialOrbit);
         kepler2.setSlaveMode();  //slave mode is default
         
         
         EventDetector sta1Visi2 =
                 new ElevationDetector(maxcheck, threshold, staF).
                 withConstantElevation(elevation).    //setting the minimum elevation to check for
                 withHandler(new StopOnEvent());
              
               EventsLogger logger2=new EventsLogger();
               kepler2.addEventDetector(logger2.monitorDetector(sta1Visi2)); 
         SpacecraftState asdf=kepler2.propagate(finalDate);
         System.out.println(asdf.getDate().toString());
         
         
         //System.out.println(allEvents.get(10).getState().getDate().toString());
         
         /*
         kepler.addEventDetector(sta1Visi);
         SpacecraftState asdf= kepler.propagate(finalDate);
         
         System.out.println(asdf.getDate().toString());
         System.out.println(finalDate.toString());
         double hi= 0.000;
         */
         
         //defining writer object to write to text file.
         //PrintWriter unWriter= new PrintWriter("/Users/mtruong/Desktop/JAVA/testOrbitData.txt", "UTF-8");
        
        /* 
         while (extrapDate.compareTo(finalDate) <= 0)  {

             //Position and velocit coordinates of our satellite
             PVCoordinates pvInert   = kepler.propagate(extrapDate).getPVCoordinates();
             Vector3D positionVectorSatellite=pvInert.getPosition();
             
             
             // POSITION AND VELOCITY OF KEPLER OBJECT IN THE STATION FRAME.
             PVCoordinates pvStation = inertialFrame.getTransformTo(staF, extrapDate).transformPVCoordinates(pvInert);
             
             // GET POSITION VECTOR OF SATELLITE WRT STATION
             Vector3D positionVectorWRTStation=pvStation.getPosition();
             
             
             // Get the azimuth, elevation, and range from the .get functions. 
             // in the staF (station frame). Get the position vector in the inertial frame
             double azimuth=staF.getAzimuth(positionVectorSatellite, inertialFrame, extrapDate);
             double elevation=staF.getElevation(positionVectorSatellite, inertialFrame, extrapDate);
             double range=staF.getRange(positionVectorSatellite, inertialFrame, extrapDate);
             
         	//getAzimuth(Vector3D extPoint, Frame frame, AbsoluteDate date)
             System.out.format("Elevation: %3.3f   Azimuth: %3.3f    Range: %3.3f %n" ,elevation*180/Math.PI, azimuth*180/Math.PI,range);
             
             
             
             //System.out.format("%3.3f	%3.3f%n",pvStation.getLatitude()*180/Math.PI,pvStation.getLongitude()*180/Math.PI);
             
             //GeodeticPoint satLatLonAlt = earth.transform(pvInert.getPosition(), FramesFactory.getEME2000(), extrapDate);
             //System.out.println(satLatLonAlt);
             
             //System.out.format("%3.3f	%3.3f%n",satLatLonAlt.getLatitude()*180/Math.PI,satLatLonAlt.getLongitude()*180/Math.PI);
             //unWriter.printf("%3.3f	%3.3f%n", satLatLonAlt.getLatitude()*180/Math.PI,satLatLonAlt.getLongitude()*180/Math.PI);
             extrapDate = extrapDate.shiftedBy(60);
             
         }
		 */
		 
        // unWriter.close();
	}
}

//3 ground stations
//Event detectors to see when all three are in view
//azimuth, elevation, and range with respect to one station
//time, az el, range then abs lat , long, altitude("elevation")
//to matlab



