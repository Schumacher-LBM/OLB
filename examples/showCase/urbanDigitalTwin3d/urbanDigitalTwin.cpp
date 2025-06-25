/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Dennis Teutscher
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */
/*
  This configuration is based on the preprint "A Digital Urban Twin Enabling Interactive Pollution Predictions and Enhanced Planning"
  (DOI: https://doi.org/10.48550/arXiv.2502.13746), which is to be published in *Building and Environment*.

  You can download the required files using the "download_data.sh" script located in the "example" folder.
  For quicker testing, use "map_small.osm", which provides faster results. The full-scale case from the paper
  is represented by "map.osm".

  To run the small case, define `small`; to run the full case, define `original`.
  This example requires the PROJ library.
  Ensure it is installed on your system, then enable it by adding 'PROJ' to the FEATURES variable in config.mk:
  FEATURES := PROJ
*/

#define small//small//original

#ifdef FEATURE_PROJ

#include "tinyxml2.h"
#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <proj.h>
#include <limits>
#include <cmath>
#include "olb.h"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = float;
using DESCRIPTOR = D3Q19<AVERAGE_VELOCITY,POROSITY,collision::LES::SMAGORINSKY>;
using MOMENTA = momenta::AdvectionDiffusionBulkTuple;
using ADDESCRIPTOR = D3Q7<VELOCITY,OMEGA,AVERAGE_DENSITY,POROSITY,LATTICE_TIME,SCALAR,SOURCE>;
using ADBulkDynamics = dynamics::Tuple<T, ADDESCRIPTOR, MOMENTA, equilibria::FirstOrder, collision::ParameterFromCell<descriptors::OMEGA, collision::BGK>,AdvectionDiffusionExternalVelocityCollision>;



using BulkDynamics = SmagorinskyBGKdynamics<T,DESCRIPTOR>;
using PorousSmagorinskyRLBthirdOrder = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Porous<momenta::BulkTuple>,
  equilibria::ThirdOrder,
  collision::ParameterFromCell<collision::LES::SMAGORINSKY,collision::SmagorinskyEffectiveOmega<collision::RLBThirdOrder>>
>;

// Parameters for the simulation setup
#ifdef small
  const T update=20.; //time after which the simulation is updated with new wind and pollution data (Value is set for the small geometry "map_small.osm". For map.osm a value of 400 can be used.)
#elif defined(original)
  const T update=400.;
#endif

const int nOTimeStamps = 24; // the number of time stamps from the meassuring data. In this case each time stamp represents 1 h. (With the data for this example it can be extended up to 160 h)
const T maxPhysT = update*nOTimeStamps; // max. simulation time in s, SI unit (update * timeStamps from the meassuring data)


T rhoBurg;
T rhoLeder;
T rhoExhaust;

T winddirection=0;
T winddirectionOld=0;
int step =1;
const T D = 0.0000233;
T Kin = 0.02*0.02/150* std::pow(0.97,3)/std::pow(1.-0.97,2);
T windSpeedCurrent=0;
T windSpeedOld=0;
int counter = 0; // Initialisierter Zähler
std::map<std::string, std::map<std::string, std::map<std::string, std::string>>> airQualityDataLederstrasse;
std::map<std::string, std::map<std::string, std::map<std::string, std::string>>> airQualityDataAlteburgstrasse;
std::map<std::string, std::map<std::string, std::pair<std::string, std::string>>> weatherData;

template <typename T, typename S>
class Wind3D : public AnalyticalF3D<T,S> {

protected:
  T um, delta, alpha;
  int zDir, flowDir;

public:
  Wind3D(T um_, T delta_, T alpha_, int zDir_, int flowDir_) : AnalyticalF3D<T, S>(3)
  {
    um = um_;
    delta = delta_;
    alpha = alpha_;
    zDir = zDir_;
    flowDir = flowDir_;
  };

  bool operator()(T output[], const S input[])
  {
    output[0] = 0.;
    output[1] = 0.;
    output[2] = 0.;
    output[flowDir] = um*util::pow((input[zDir]/delta),alpha);

    return true;
  };
};
void readWeatherData(
    const std::string& filename,
    const std::string& dateCol,
    const std::string& timeCol,
    const std::string& datetimeCol,
    bool splitDatetime,
    const std::string& wdirCol,
    const std::string& wspdCol,
    std::map<std::string, std::map<std::string, std::pair<std::string, std::string>>>& weatherData
) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line;
    std::getline(file, line); // Read header
    std::stringstream headerStream(line);
    std::string column;
    std::vector<std::string> headers;
    std::map<std::string, size_t> columnIndices;

    while (std::getline(headerStream, column, ';')) {
        headers.push_back(column);
        columnIndices[column] = headers.size() - 1;
    }

    if (splitDatetime) {
        if (columnIndices.count(datetimeCol) == 0) {
            std::cerr << "Error: datetime column not found." << std::endl;
            return;
        }
    } else {
        if (columnIndices.count(dateCol) == 0 || columnIndices.count(timeCol) == 0) {
            std::cerr << "Error: date or time column not found." << std::endl;
            return;
        }
    }

    if (columnIndices.count(wdirCol) == 0 || columnIndices.count(wspdCol) == 0) {
        std::cerr << "Error: wind direction or speed column not found." << std::endl;
        return;
    }

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> row;
        while (std::getline(ss, cell, ';')) {
            row.push_back(cell);
        }

        if (row.size() < headers.size()) continue;

        std::string date, time;
        if (splitDatetime) {
            std::string datetime = row[columnIndices[datetimeCol]];
            if (datetime.length() >= 16) {
                date = datetime.substr(0, 10);
                time = datetime.substr(11, 5);
            } else {
                std::cerr << "Skipping invalid datetime format: " << datetime << std::endl;
                continue;
            }
        } else {
            date = row[columnIndices[dateCol]];
            time = row[columnIndices[timeCol]];
        }

        std::string wdir = row[columnIndices[wdirCol]];
        std::string wspd = row[columnIndices[wspdCol]];

        weatherData[date][time] = {wdir, wspd};
    }

    file.close();
}

void readAirQualityData(
    const std::string& filename,
    const std::string& dateCol,
    const std::string& timeCol,
    const std::string& datetimeCol,
    bool splitDatetime, // true = use datetimeCol and split it; false = use dateCol + timeCol
    const std::string& pollutantCol,
    const std::string& valueCol,
    std::map<std::string, std::map<std::string, std::map<std::string, std::string>>>& airQualityData
) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line;
    std::getline(file, line); // Read header
    std::stringstream headerStream(line);
    std::string column;
    std::vector<std::string> headers;
    std::map<std::string, size_t> columnIndices;

    while (std::getline(headerStream, column, ';')) {
        headers.push_back(column);
        columnIndices[column] = headers.size() - 1;
    }

    if (splitDatetime) {
        if (columnIndices.count(datetimeCol) == 0) {
            std::cerr << "Error: datetime column '" << datetimeCol << "' not found." << std::endl;
            return;
        }
    } else {
        if (columnIndices.count(dateCol) == 0 || columnIndices.count(timeCol) == 0) {
            std::cerr << "Error: date or time column not found." << std::endl;
            return;
        }
    }

    if (columnIndices.count(pollutantCol) == 0 || columnIndices.count(valueCol) == 0) {
        std::cerr << "Error: pollutant or value column not found." << std::endl;
        return;
    }

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> row;
        while (std::getline(ss, cell, ';')) {
            row.push_back(cell);
        }

        if (row.size() < headers.size()) continue;

        std::string date, time;
        if (splitDatetime) {
            std::string datetime = row[columnIndices[datetimeCol]];
            if (datetime.length() >= 16) {
                date = datetime.substr(0, 10);
                time = datetime.substr(11, 5);
            } else {
                std::cerr << "Skipping invalid datetime format: " << datetime << std::endl;
                continue;
            }
        } else {
            date = row[columnIndices[dateCol]];
            time = row[columnIndices[timeCol]];
        }

        std::string pollutant = row[columnIndices[pollutantCol]];
        std::string value = row[columnIndices[valueCol]];

        airQualityData[date][time][pollutant] = value;
    }

    file.close();
}

// Function to search for a road by name
std::list<typename OSMParser<T>::RoadStatistics>::const_iterator findRoadByName(
    const std::list<typename OSMParser<T>::RoadStatistics>& roadStatisticsList,
    const std::string& roadName)
{
    return std::find_if(roadStatisticsList.begin(), roadStatisticsList.end(),
        [&roadName](const OSMParser<T>::RoadStatistics& road) {
            return road.name == roadName; // Compare the road name
        });
}

// Adjusting date format from yyyy-mm-dd to dd.mm.yyyy
std::string convertDateFormat(const std::string& date) {
    // Date format conversion from yyyy-mm-dd to dd.mm.yyyy
    return date.substr(8, 2) + "." + date.substr(5, 2) + "." + date.substr(0, 4);
}


// Function now returns a tuple containing the date and a vector of pollutants with their concentrations
std::tuple<std::string, std::vector<std::pair<std::string, std::string>>> getConcentrationData(
    const std::map<std::string, std::map<std::string, std::map<std::string, std::string>>>& airQualityData,
    const std::string& weatherDate, const std::string& weatherTime)
{
    OstreamManager clout("getConcentration");
    // Convert weatherDate from yyyy-mm-dd to dd.mm.yyyy format
    std::string convertedDate = convertDateFormat(weatherDate);
    clout << convertedDate << std::endl;

    // Quote the weather time
    std::string quotedWeatherTime = "'" + weatherTime + "'";
    clout << quotedWeatherTime << std::endl;

    // Iterate through air quality data
    for (const auto& stationData : airQualityData) {
        const std::string& date = stationData.first;
        const auto& stationTimeData = stationData.second;

        // Check if the date matches
        if (convertedDate == date) {
            for (const auto& timeEntry : stationTimeData) {
                const std::string& time = timeEntry.first;

                // Check if the time matches
                if (time == quotedWeatherTime) {
                    const auto& pollutantData = timeEntry.second;

                    // Collect all pollutants and their concentrations
                    std::vector<std::pair<std::string, std::string>> pollutants;
                    for (const auto& pollutant : pollutantData) {
                        const std::string& type = pollutant.first;
                        std::string concentration = pollutant.second;

                        // Replace "-" with an empty string
                        if (concentration == "-") {
                            concentration = "";
                        }

                        // Add pollutant type and concentration to the vector
                        pollutants.emplace_back(type, concentration);
                        clout << "Pollutant: " << type << ", Concentration: " << concentration << std::endl;
                    }

                    // Return the converted date and the vector of pollutants
                    return {convertedDate, pollutants};
                }
            }
        }
    }

    // If no matching concentration data is found, return an empty tuple with an empty vector
    return { "", {} };
}

std::tuple<std::string, std::string, std::string, std::string> getWeatherData(
    const std::map<std::string, std::map<std::string, std::pair<std::string, std::string>>>& weatherData,
    int counter)
{
    // Make sure the counter is valid
    if (counter < 0) {
        throw std::out_of_range("Counter cannot be negative");
    }

    // Iterate through weather data entries
    int currentIndex = 0;

    for (const auto& dateData : weatherData) {
        const std::string& date = dateData.first;
        const auto& timeData = dateData.second;

        // Iterate through the time data for the current date
        for (const auto& timeEntry : timeData) {
            const std::string& time = timeEntry.first;
            const auto& windData = timeEntry.second;  // Contains wind direction and speed

            // Check if the current index matches the counter
            if (currentIndex == counter) {
                // Replace ',' with '.' for the values to ensure correct format
                std::string windDirection = windData.first;
                std::string windSpeed = windData.second;

                for (char& c : windDirection) {
                    if (c == ',') {
                        c = '.';
                    }
                }
                for (char& c : windSpeed) {
                    if (c == ',') {
                        c = '.';
                    }
                }

                return std::make_tuple(date, time, windDirection, windSpeed);
            }
            currentIndex++;  // Increment the counter
        }
    }

    // Throw an exception if the counter exceeds available weather data
    throw std::out_of_range("Counter exceeds available weather data");
}

// Stores data from stl file in geometry in form of material numbers
void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter, SuperGeometry<T,3>& superGeometry,  OSMParser<T> osmParser,IndicatorF3D<T>& bounding)
{

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0,2 );
  superGeometry.rename( 2,1,{1,1,1} );
  {
    auto origin = bounding.getMin();
    origin[2] = -converter.getPhysDeltaX();
    auto extent = bounding.getMax() - bounding.getMin();
    extent[2] = converter.getPhysDeltaX();
    extent[0] = bounding.getMax()[0]+bounding.getMax()[0];
    extent[1] = bounding.getMax()[1]+bounding.getMax()[0];
    IndicatorCuboid3D<T> floor(extent, origin);
    superGeometry.rename(2,4,1,floor);
  }
  {
    auto origin = bounding.getMax()-converter.getPhysDeltaX();
    origin[0]= bounding.getMin()[0];
    origin[1]= bounding.getMin()[1];
    auto extent = bounding.getMax() - bounding.getMin();
    extent[2] = 0.1*extent[2];
    extent[0] = bounding.getMax()[0]+bounding.getMax()[0];
    extent[1] = bounding.getMax()[1]+bounding.getMax()[0];
    IndicatorCuboid3D<T> sky(extent, origin);
    superGeometry.rename(2,5,1,sky);
  }
  osmParser.createGeometry(superGeometry,osmParser.BUILDING,1,4,converter.getPhysDeltaX());
  osmParser.createGeometry(superGeometry,osmParser.TREES,1,9,converter.getPhysDeltaX());
  osmParser.createGeometry(superGeometry,osmParser.AREAVEGITATION,1,9,converter.getPhysDeltaX());
  osmParser.createGeometry(superGeometry,osmParser.ROADS,4,6,converter.getPhysDeltaX(),"Alteburgstraße");
  osmParser.createGeometry(superGeometry,osmParser.ROADS,4,7,converter.getPhysDeltaX(),"Lederstraße");
  osmParser.generateBuildingExhaust(superGeometry,1*converter.getPhysDeltaX(),1,8);

//small sockel for trees. Instabel if to close to source
  {
    auto origin = bounding.getMin();
    auto extent = bounding.getMax() - bounding.getMin();
    extent[2] = converter.getPhysDeltaX();
    extent[0] = bounding.getMax()[0]+0.1*bounding.getMax()[0];
    extent[1] = bounding.getMax()[1]+0.1*bounding.getMax()[0];
    IndicatorCuboid3D<T> floor(extent, origin);
    superGeometry.rename(9,4,floor);
  }


  superGeometry.checkForErrors();
  superGeometry.print();
  clout << "Geometry build!" <<std::endl;
}

// Set up the geometry of the simulation
void prepareLattice( SuperLattice<T,DESCRIPTOR>& sLattice,
                     UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperGeometry<T,3>& superGeometry )
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();
  T latticeDist = 1. - converter.getLatticeLength(3.3e-3);

  // Material=0 -->do nothing
  sLattice.defineDynamics<NoDynamics>(superGeometry,0);

  // Material=1 -->bulk dynamics
  auto bulkIndicator = superGeometry.getMaterialIndicator({1,2,5,9});
  sLattice.defineDynamics<PorousSmagorinskyRLBthirdOrder::template wrap_collision<collision::TrackAverageVelocity>>(bulkIndicator);

  // Material=2 -->bounce back
  sLattice.defineDynamics<BounceBack>(superGeometry, 4);
  sLattice.defineDynamics<BounceBack>(superGeometry, 6);
  sLattice.defineDynamics<BounceBack>(superGeometry, 7);

  boundary::set<boundary::InterpolatedVelocity>(sLattice, superGeometry.getMaterialIndicator(2));
  boundary::set<boundary::FullSlip>(sLattice, superGeometry, 5);

  // Initial conditions
  AnalyticalConst3D<T,T> rhoF( 1 );
  Vector<T,3> velocityV;
  AnalyticalConst3D<T,T> uF(velocityV);
  // Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU( bulkIndicator, rhoF, uF );
  sLattice.iniEquilibrium( bulkIndicator, rhoF, uF );

  sLattice.setParameter<OMEGA>(omega);
  sLattice.setParameter<collision::LES::SMAGORINSKY>(0.12);

  sLattice.setParameter<LATTICE_TIME>(1);

  AnalyticalConst3D<T,T> initialPorosityF(1);
  sLattice.defineField<POROSITY>(superGeometry.getMaterialIndicator({0,1,2,4,5,6,7,8}), initialPorosityF);


  T h = converter.getPhysDeltaX();
  T d = 0.9992;
  clout<< "Porosity d = "<<d <<std::endl;
  //Defining porosity for geometry (0 not permeable)
  AnalyticalConst3D<T,T> solidPorosityF(d);
  sLattice.defineField<POROSITY>(superGeometry.getMaterialIndicator(9),solidPorosityF);
  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

void prepareLatticeAD(SuperLattice<T,ADDESCRIPTOR>& sLatticeAD,
                    SuperGeometry<T,3>& sGeometry,
                    const UnitConverter<T,DESCRIPTOR>& converter)
{

  OstreamManager clout( std::cout,"prepareLatticeAD" );
  clout << "Prepare Lattice for ADE..." << std::endl;

  const T omegaA = converter.getLatticeRelaxationFrequencyFromDiffusivity<ADDESCRIPTOR>(D);


  // Material=0 -->do nothing
  sLatticeAD.defineDynamics<NoDynamics>(sGeometry, 0);

  auto bulkIndicatorAD = sGeometry.getMaterialIndicator({1,6,7,8,9});
  sLatticeAD.defineDynamics<ADBulkDynamics::template wrap_collision<collision::TrackAverageDensity>>(bulkIndicatorAD);

  //outlet
 // setZeroGradientBoundary<T,ADDESCRIPTOR>(sLatticeAD,sGeometry.getMaterialIndicator({3,4}));
  boundary::set<boundary::ZeroDistribution>(sLatticeAD,sGeometry.getMaterialIndicator({2,5}));

  //walls
  auto bounceIndicator = sGeometry.getMaterialIndicator({4});
  sLatticeAD.defineDynamics<BounceBack>(bounceIndicator);

  //emitters
  sLatticeAD.defineDynamics<EquilibriumBoundaryFirstOrder>(sGeometry, 6);
  sLatticeAD.defineDynamics<EquilibriumBoundaryFirstOrder>(sGeometry, 7);
  sLatticeAD.defineDynamics<EquilibriumBoundaryFirstOrder>(sGeometry, 8);

  //boundary::set<boundary::AdvectionDiffusionDirichlet>(sLatticeAD,sGeometry.getMaterialIndicator({6,7}));
  sLatticeAD.setParameter<LATTICE_TIME>(1);


  AnalyticalConst3D<T,T> rho0( converter.getLatticeDensity(0) );
  AnalyticalConst3D<T,T> rho1( converter.getLatticeDensity(0) );

  AnalyticalConst3D<T,T> uF( T( 0 ), T( 0 ), T(0) );
  AnalyticalConst3D<T,T> u0( T( 0 ), T( 0 ), T(0) );
  //sLatticeAD.defineField<VELOCITY>(sGeometry.getMaterialIndicator({1,2,4,5,6,7,8,9}),uF);

  AnalyticalConst3D<T,T> om( omegaA );

  //Initialize all values of distribution functions to their local equilibrium
  sLatticeAD.defineRhoU( sGeometry, 6, rho1,u0 );
  sLatticeAD.iniEquilibrium(sGeometry, 6, rho1, u0 );
  sLatticeAD.defineRhoU( sGeometry, 7, rho1,u0 );
  sLatticeAD.iniEquilibrium(sGeometry, 7, rho1, u0 );
  sLatticeAD.defineRhoU( sGeometry, 8, rho1,u0 );
  sLatticeAD.iniEquilibrium(sGeometry, 8, rho1, u0 );

  sLatticeAD.defineRho( sGeometry.getMaterialIndicator({1,2,4,5,9}), rho0 );
  sLatticeAD.iniEquilibrium( sGeometry.getMaterialIndicator({1,2,4,5,9}), rho0, u0 );

  sLatticeAD.defineField<SOURCE>(sGeometry.getMaterialIndicator({1,2,4,5,6,7,8,9}), rho0);



  sLatticeAD.setParameter<OMEGA>(omegaA);
  AnalyticalConst3D<T,T> omegaAF(omegaA);
  sLatticeAD.defineField<OMEGA>(sGeometry.getMaterialIndicator({1,2,4,5,6,7,8,9}), omegaAF);

  AnalyticalConst3D<T,T> initialPorosityF(1);
  sLatticeAD.defineField<POROSITY>(sGeometry.getMaterialIndicator({0,1,2,3,4,5,6,7,8}), initialPorosityF);

  T h = converter.getPhysDeltaX();
  T d = 0.9992;
  clout<< "Porosity d = "<<d <<std::endl;
  //Defining porosity for geometry (0 not permeable)
  AnalyticalConst3D<T,T> solidPorosityF(d);
  sLatticeAD.defineField<POROSITY>(sGeometry.getMaterialIndicator(9),solidPorosityF);

  sLatticeAD.initialize();


  {
    auto& communicator = sLatticeAD.getCommunicator(stage::PostCoupling());
    communicator.requestField<VELOCITY>();
    communicator.requestField<OMEGA>();
    communicator.requestOverlap(2);
    communicator.exchangeRequests();
  }
  clout << "Omega ADE: " << 1/omegaA << std::endl;
  clout << "Prepare Lattice for ADE... OK" << std::endl;
}
// Function to interpolate between the old and new values
void interpolateWind(T &winddirectionOld, T &windSpeedOld,
                     T winddirection, T windSpeedCurrent,
                     int counter, int totalSteps) {
    std::ostream &clout = std::cout; // Output stream for logging (can be customized)

    // Interpolation step size for wind direction and wind speed
    T windDirectionStep = (winddirection - winddirectionOld) / totalSteps;
    T windSpeedStep = (windSpeedCurrent - windSpeedOld) / totalSteps;

    // Incremental interpolation based on the current counter
    winddirectionOld += windDirectionStep * counter; // Move towards new wind direction
    windSpeedOld += windSpeedStep * counter; // Move towards new wind speed

    // To keep the wind direction within 0°-360°, we use fmod for modular arithmetic
    winddirectionOld = fmod(winddirectionOld, 360.0);
    if (winddirectionOld < 0) winddirectionOld += 360.0;

    // After the loop, set the final values to the new values
    if (counter == totalSteps) {
        winddirectionOld = winddirection; // Final wind direction
        windSpeedOld = windSpeedCurrent; // Final wind speed
    }
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setBoundaryValues( SuperLattice<T, DESCRIPTOR>& sLattice,SuperLattice<T, ADDESCRIPTOR>& sLatticeAD,
                        UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                        SuperGeometry<T,3>& superGeometry, OSMParser<T> osmParser, int updateCounter)
{
  OstreamManager clout( std::cout,"setBoundaryValues" );

  // No of time steps for smooth start-up
  #ifdef small
    int iTmaxStart = converter.getLatticeTime(15.);
  #elif defined(original)
    int iTmaxStart = converter.getLatticeTime(45.);
  #endif

  if (iT % converter.getLatticeTime(update) == 0) {  // Every 60th timestep
    try {
      // Access weather data based on the counter
      auto [weatherDate, weatherTime, windDirection, windSpeed] = getWeatherData(weatherData, counter);
      T windDirectionValue = std::stod(windDirection);  // Wind direction
      T windSpeedValue = std::stod(windSpeed);    // Wind speed

      T PM25_burg =0.;
      T NO2_burg =0.;
      auto [date, pollutants] = getConcentrationData(airQualityDataAlteburgstrasse, weatherDate, weatherTime);
      for (const auto& [type, concentration] : pollutants) {
        //clout << "Alte Burgstrasse: "<<"Pollutant: " << type << ", Concentration: " << concentration << std::endl;
        if(type=="\"Feinstaub (PM₂,₅)\"")
        {
          if(concentration =="")
          {
            PM25_burg=0.;
          }
          else
          {
            PM25_burg=stod(concentration);
          }
        }
        else
        {
          if(concentration =="")
          {
            NO2_burg=0.;
          }
          else
          {
            NO2_burg=stod(concentration);
          }
        }
      }
      auto [date2, pollutants2] = getConcentrationData(airQualityDataLederstrasse, weatherDate, weatherTime);
      T PM25_leder =0.;
      T NO2_leder =0.;
      for (const auto& [type, concentration] : pollutants2) {
        //clout << "Lederstrasse: "<<"Pollutant: " << type << ", Concentration: " << concentration << std::endl;
        if(type=="\"Feinstaub (PM₂,₅)\"")
        {
          if(concentration =="")
          {
            PM25_leder=0.;
          }
          else
          {
            PM25_leder=stod(concentration);
          }
        }
        else
        {
          if(concentration =="")
          {
            NO2_leder=0.;
          }
          else
          {
            NO2_leder=stod(concentration);
          }
        }
      }
      ++counter;  // Increment the counter to process the next data set

    std::map<std::string, std::pair<T, T>> nitrogenConcentrations = {
      {"Alteburgstraße", {NO2_burg,PM25_burg}},{"Lederstraße",{NO2_leder,PM25_leder}}};
    auto stats = osmParser.calculateStreetStatistics(nitrogenConcentrations,converter.getPhysDeltaX());

    clout <<updateCounter <<"---------------------------------------------------------------------"<<std::endl;
    clout <<" Date: " <<weatherDate <<" Time: " << weatherTime << " Wind speed: "<<windSpeedValue<<" Wind direction: "<<windDirection<<std::endl;
    clout << "---------------------------------------------------------------------"<<std::endl;
    clout << "Streetname " << "Alteburgstrasse"<<std::endl;
    clout <<"Pollutant: PM25 " << ", measured concentration: " << PM25_burg << std::endl;
    clout <<"Pollutant: NO2 " << ", measured concentration: " << NO2_burg << std::endl;
    auto roadStats= findRoadByName(stats,"Alteburgstraße");
    clout <<"Estimated cars: " << roadStats->numberCars << std::endl;
    clout <<"PM25 from cars: " << roadStats->PM25_cars << std::endl;
    clout <<"PM25 from buildings: " << roadStats->PM25_buildings << std::endl;
    clout << "---------------------------------------------------------------------"<<std::endl;
    clout << "Streetname " << "Lederstrasse"<<std::endl;
    clout <<"Pollutant: PM25 " << ", measured concentration: " << PM25_leder << std::endl;
    clout <<"Pollutant: NO2 " << ", measured concentration: " << NO2_leder << std::endl;
    auto roadStats_= findRoadByName(stats,"Lederstraße");
    clout <<"Estimated cars: " << roadStats_->numberCars << std::endl;
    clout <<"PM25 from cars: " << roadStats_->PM25_cars << std::endl;
    clout <<"PM25 from buildings: " << roadStats_->PM25_buildings << std::endl;
    clout << "---------------------------------------------------------------------"<<std::endl;
    sLatticeAD.setProcessingContext(ProcessingContext::Evaluation);
    AnalyticalConst3D<T,T> rho6(roadStats->PM25_cars  );
    AnalyticalConst3D<T,T> rhoF(1);
    AnalyticalConst3D<T,T> rho7( roadStats_->PM25_cars  );
    AnalyticalConst3D<T,T> rho8( (roadStats->PM25_buildings + roadStats_->PM25_buildings)/2. );

    rhoBurg = roadStats->PM25_cars;
    rhoLeder= roadStats_->PM25_cars;
    rhoExhaust=roadStats->PM25_buildings + roadStats_->PM25_buildings;
    AnalyticalConst3D<T,T> u0( T( 0 ), T( 0 ), T(0) );
    sLatticeAD.defineRhoU( superGeometry, 6, rho6,u0 );
    sLatticeAD.iniEquilibrium(superGeometry, 6, rho6, u0 );
    sLatticeAD.defineRhoU( superGeometry, 7, rho7,u0 );
    sLatticeAD.iniEquilibrium(superGeometry, 7, rho7, u0 );
    sLatticeAD.defineRhoU( superGeometry, 8, rho8,u0 );
    sLatticeAD.iniEquilibrium(superGeometry, 8, rho8, u0 );
    sLatticeAD.setProcessingContext(ProcessingContext::Simulation);
    winddirectionOld=winddirection;
    winddirection=windDirectionValue;
    windSpeedOld=windSpeedCurrent;
    windSpeedCurrent=windSpeedValue;
    updateCounter=0;
    } catch (const std::out_of_range& e) {
        clout << "No more weather data available for counter: " << counter << std::endl;
    }

  }
  // Time allocated for adjustment to a new wind direction.
  // If instabilities occur during the simulation, consider increasing this value.
  #ifdef small
    const T window= converter.getLatticeTime(5.);
  #elif defined(original)
    const T window= converter.getLatticeTime(15.);
  #endif
  if (iT % converter.getLatticeTime(update) == 0)
  {
    if(step==window)
      step=0;

  }
  if (iT >= converter.getLatticeTime(update) &&
      converter.getLatticeTime(iT)%30==0 && converter.getLatticeTime(step) < window){

    // Smooth start curve, polynomial
    SinusStartScale<T,int> StartScale( window, T( 1 ) );

    // Creates and sets the Poiseuille inflow profile using functors
    int iTvec[1] = {step};
    T frac[1] = {};
    StartScale( frac,iTvec );
    AnalyticalConst3D<T,T> uD(T(0),T(0),converter.getLatticeVelocity(0.2054)*frac[0]);

    sLatticeAD.defineU( superGeometry, 6,uD );
    sLatticeAD.defineU( superGeometry, 7,uD );
    sLatticeAD.defineU( superGeometry, 8,uD );
    sLatticeAD.defineField<VELOCITY>(superGeometry.getMaterialIndicator({6,7,8}), uD);

    interpolateWind(winddirectionOld,windSpeedOld,winddirection,windSpeedCurrent,step,window);
    T windRad = M_PI/180.*winddirectionOld;
    AnalyticalWindProfileF3D<T> uF(converter.getLatticeVelocity(windSpeedOld),0.1,10.,0.4,0,0.,2,windRad);
    step++;

    sLattice.defineU(superGeometry,2,uF);
    sLattice.template setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
    sLatticeAD.template setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
    sLatticeAD.template setProcessingContext<Array<momenta::FixedDensity::RHO>>(
        ProcessingContext::Simulation);
  }


  if (iT<=iTmaxStart && iT%30==0&&counter!=0){//iT <= iTmaxStart && iT%30==0. ) {
    // Smooth start curve, polynomial
    SinusStartScale<T,int> StartScale( iTmaxStart, T( 1 ) );

    // Creates and sets the Poiseuille inflow profile using functors
    int iTvec[1] = {updateCounter};
    T frac[1] = {};
    StartScale( frac,iTvec );
    T windRad = M_PI/180.*winddirection;
    AnalyticalWindProfileF3D<T> uF(converter.getLatticeVelocity(windSpeedCurrent)*frac[0],0.1,10.,0.4,0,0.,2,windRad);
    AnalyticalConst3D<T,T> uD(T(0),T(0),converter.getLatticeVelocity(0.2054)*frac[0]);
    sLatticeAD.defineU(superGeometry, 6, uD);
    sLatticeAD.defineU(superGeometry, 7, uD);
    sLatticeAD.defineU(superGeometry, 8, uD);
    sLatticeAD.defineField<VELOCITY>(superGeometry.getMaterialIndicator({6,7,8}), uD);
    sLattice.defineU(superGeometry,2,uF);
    sLattice.template setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
    sLatticeAD.template setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
    sLatticeAD.template setProcessingContext<Array<momenta::FixedDensity::RHO>>(
        ProcessingContext::Simulation);
  }
}



// Computes the pressure drop between the voxels before and after the cylinder
void getResults( SuperLattice<T, DESCRIPTOR>& sLattice,SuperLattice<T, ADDESCRIPTOR>& sLatticeAD,IndicatorF3D<T>& bounding,
                 UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                 SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer, int updateCounter)
{

  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter3D<T> vtmWriter( "city" );
  SuperVTMwriter3D<T> vtmWriterADE( "cityAD" );

  const int vtkIter  = converter.getLatticeTime( update-3);
  const int statIter = converter.getLatticeTime( 5);

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank(sLattice);
    SuperGeometryF3D<T> geom (superGeometry);
    vtmWriter.write(geom);
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );

    vtmWriter.createMasterFile();
    vtmWriterADE.createMasterFile();
  }

  // Writes output on the console
  if ( iT%statIter == 0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
    sLatticeAD.getStatistics().print( iT,converter.getPhysTime( iT ) );
  }
  #ifdef small
    int iTstartAvg = converter.getLatticeTime(15.);
  #elif defined(original)
    int iTstartAvg = converter.getLatticeTime(45.);
  #endif
  if (iT == iTstartAvg) {
    SuperLatticeDensity3D<T,ADDESCRIPTOR> density(sLatticeAD);
    SuperLatticePhysVelocity3D<T,DESCRIPTOR> latticeVelocity(sLattice,converter);
    sLattice.defineField<AVERAGE_VELOCITY>(superGeometry.getMaterialIndicator({1,2,4,5,6,7,8,9}), latticeVelocity);
    sLatticeAD.defineField<AVERAGE_DENSITY>(superGeometry.getMaterialIndicator({1,2,4,5,6,7,8,9}),density);
  }
  if (iT < iTstartAvg) {
    sLatticeAD.setParameter<LATTICE_TIME>(2); // dummy just to prevent nan
    sLattice.setParameter<LATTICE_TIME>(2);
  }
  else {
    sLatticeAD.setParameter<LATTICE_TIME>(iT - iTstartAvg + 1);
    sLattice.setParameter<LATTICE_TIME>(iT - iTstartAvg + 1);
  }
  // Writes the vtk files
  if ( iT%vtkIter == 0) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    sLattice.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtmWriter("city");
      SuperLatticePhysVelocity3D velocity(sLattice, converter);
      SuperLatticePhysExternalVectorField3D<T,DESCRIPTOR,AVERAGE_VELOCITY> velocityAv(sLattice,
            converter.getConversionFactorVelocity(),"aVel");
      SuperLatticePhysPressure3D pressure(sLattice, converter);
      vtmWriter.addFunctor(velocity);
      vtmWriter.addFunctor(velocityAv);
      vtmWriter.addFunctor(pressure);
      task(vtmWriter, iT);
    });
    sLatticeAD.setProcessingContext(ProcessingContext::Evaluation);
    sLatticeAD.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtmWriterADE("cityAD");
      SuperLatticeDensity3D<T, ADDESCRIPTOR> concentration( sLatticeAD );
      SuperLatticeExternalScalarField3D<T,ADDESCRIPTOR,descriptors::AVERAGE_DENSITY> avgD(sLatticeAD);
      avgD.getName() = "average(Density)";
      vtmWriterADE.addFunctor(avgD);
      vtmWriterADE.addFunctor( concentration );
      task(vtmWriterADE, iT);
    });
  }
  if ( iT ==converter.getLatticeTime(30)) {
    SuperLatticePhysVelocity3D<T,DESCRIPTOR> velocity(sLattice, converter);
    SuperEuklidNorm3D<T> normVel( velocity );
    auto measure = (bounding.getMin()[0]+2*converter.getPhysDeltaX()+ bounding.getMax()[0]) / 2.0;

    auto line3d = Line3D<T>().originAt(measure).parallelTo({0,0,1});
    LineLattice3D<T> lineLattice(superGeometry.getCuboidDecomposition(),line3d);
    BlockReduction3D1D<T> line(normVel,lineLattice,BlockDataSyncMode::ReduceAndBcast,BlockDataReductionMode::Discrete);
    int n = line.getN();
    static Gnuplot<T> gplot( "inflow" );
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    for(int i =0; i < n; i++)
    {
      T res[3]{};
      line(res,i);
      auto physR= line.getPhysR(i);
      gplot.setData(res[0],T(physR[2]),"olb");

    }
    gplot.writePNG();
  }
}
void saveStepToFile(int iT, const std::string& filename) {
    std::ofstream outFile(filename, std::ios::trunc); // Open in truncate mode (overwrite)
    if (!outFile) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }
    outFile << "Step: " << iT << std::endl;
    outFile.close();
}

int loadStepFromFile(const std::string& filename) {
    std::ifstream inFile(filename);
    if (!inFile) {
        std::cerr << "Error: Could not open file " << filename << " for reading." << std::endl;
        return -1; // Return an invalid value if file can't be read
    }
    std::string label;
    int step;
    if (inFile >> label >> step && label == "Step:") {
        inFile.close();
        return step;
    }
    inFile.close();
    return -1; // Return invalid value if reading fails
}
int main( int argc, char* argv[] )
{

  // === 1st Step: Initialization ===
  initialize(&argc, &argv);
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );
  std::string fName( "urbanDigitalTwin.xml" );
  XMLreader config( fName );

//Weather data:
  readAirQualityData(
    "station_2024-11-06-2024-11-13_alteBurgStrasse.csv",
    "Datum", "Uhrzeit", "",  // separate columns
    false,
    "Schadstoff", "Messwert",
    airQualityDataAlteburgstrasse);
  readAirQualityData(
    "station_2024-11-06-2024-11-13_lederStrasse.csv",
    "Datum", "Uhrzeit", "",  // separate columns
    false,
    "Schadstoff", "Messwert",
    airQualityDataLederstrasse);

  // Provide the path to the CSV file
  readWeatherData(
      "winddataReutlingen6_14.CSV",
      "", "", "time",        // separate date/time not used → combined "time"
      true,                  // splitDatetime = true
      "wdir", "wspd",
      weatherData
  );

  UnitConverter<T,DESCRIPTOR>* converter = createUnitConverter<T,DESCRIPTOR>( config );

  //Load specific parameters from XML
  T logT = config["Output"]["Log"]["SaveTime"].get<T>();
  T vtkSave = config["Output"]["VisualizationVTK"]["SaveTime"].get<T>();
  int timerPrintMode = config["Output"]["Timer"]["PrintMode"].get<int>();
  std::string filenameVtk = config["Output"]["VisualizationVTK"]["Filename"].get<std::string>();

  // Prints the converter log as console output
  converter->print();
  // Writes the converter log in a file
  converter->write("city");
  #ifdef small
    OSMParser<T> osmParser("map_small.osm");
  #elif defined(original)
    OSMParser<T> osmParser("map.osm");
  #endif
  osmParser.convertToUTMLocal();
  osmParser.convertToUTMLocalTrees(osmParser.trees,true);
  auto extent = osmParser.getDomainSize(30.);
  extent += 0.4*extent;
  olb::Vector<T,3> origin(-0.2*extent[0]*converter->getPhysDeltaX(),-0.2*extent[1]*converter->getPhysDeltaX(),-converter->getPhysDeltaX());

  clout << extent << std::endl;
  olb::IndicatorCuboid3D<T> bounding(extent, origin);

  #ifdef PARALLEL_MODE_MPI
    const int noOfCuboids = singleton::mpi().getSize();
  #else
    const int noOfCuboids = 1;
  #endif
  olb::CuboidDecomposition3D<T> cGeometry(bounding,converter->getPhysDeltaX(),noOfCuboids);
  olb::BlockLoadBalancer<T> loadBalancer(cGeometry);


  // Instantiation of a superGeometry
  SuperGeometry<T,3> superGeometry( cGeometry, loadBalancer);

  prepareGeometry( *converter, superGeometry,osmParser,bounding);

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );
  SuperLattice<T, ADDESCRIPTOR> sLatticeAD( superGeometry );

  //prepareLattice and set boundaryCondition
  prepareLattice( sLattice, *converter, superGeometry );
  prepareLatticeAD(sLatticeAD,superGeometry,*converter);


  const T omega = converter->getLatticeRelaxationFrequencyFromDiffusivity<ADDESCRIPTOR>(D);
  SuperLatticeCoupling coupling(LESADECoupling<T>{},names::NavierStokes{}, sLattice, names::Concentration0{}, sLatticeAD);
  coupling.setParameter<LESADECoupling<T>::SMAGORINSKY_PREFACTOR>(0.2);
  coupling.setParameter<LESADECoupling<T>::SCHMIDT>(0.04);
  coupling.setParameter<LESADECoupling<T>::OMEGA_NSE>(converter->getLatticeRelaxationFrequency());
  coupling.setParameter<LESADECoupling<T>::OMEGA_ADE>(omega);

  SuperLatticeCoupling coupling2(
    PorousADECorrection<T>{},
    names::Concentration0{}, sLatticeAD);
  coupling2.setParameter<PorousADECorrection<T>::DIFFUSION>(D/converter->getCharPhysVelocity()*converter->getPhysDeltaX());

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converter->getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();
  int updateCounter=0;

  for (int iT =0; iT < converter->getLatticeTime( maxPhysT ); ++iT) {
    //=== 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( sLattice,sLatticeAD, *converter, iT, superGeometry,osmParser,updateCounter);

    // //// === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    coupling.execute();
    sLatticeAD.collideAndStream();
    coupling2.execute();

    updateCounter=updateCounter+1;
    // //=== 7th Step: Computation and Output of the Results ===
    getResults( sLattice,sLatticeAD,bounding,*converter, iT, superGeometry, timer,0);
  }

  timer.stop();
  timer.printSummary();
}

#else
int main(int argc, char* argv[]) { }
#endif
