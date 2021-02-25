//#include <vtkXMLUnstructuredGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkImageData.h>
//#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLImageDataReader.h>
#include <vtkHexahedron.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <cmath>
#include <vtkCellData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetSurfaceFilter.h>
//#include <experimental/filesystem>
//using namespace std::experimental::filesystem;

#include "boost/filesystem.hpp"
#include <iostream>
#include <iomanip>
using namespace boost::filesystem;

static constexpr unsigned outputPrecision = 10;

// the data for one timestep
struct CellData
{
  std::vector<double> u;
  std::vector<double> v;
  std::vector<double> p;

  std::vector< std::array<double, 3 > > cellCenters;

  size_t size;

  template<typename T>
  void resize( const T s )
  {
    assert( s > -1 );
    size = size_t(s);
    u.resize( size );
    v.resize( size );
    p.resize( size );

    cellCenters.resize( size );
  }

};

struct ParsedData
{
  CellData stokes;
  CellData darcy;

//  CellData stokesIter;
//  CellData darcyIter;
};



struct ErrorsBase
{
  double l1;
  double l2;
  double linf;

  ErrorsBase() : l1(0.), l2(0.), linf( 0.) {}
};


struct Errors
{
  ErrorsBase abs;
  ErrorsBase rel;
};


struct AllErrors
{
  Errors stokes;
  Errors darcy;
//  Errors total;
};

std::string getMonolithicName(const std::string& eq_name,
                              const std::string& n,
                              const std::string& alpha,
                              const std::string& perm,
                              const std::string& dp )
{
  return eq_name + "-" + n + "-" + alpha + "-" + perm + "-" + dp;
}

std::string getIterativeName(const std::string& eq_name,
                             const std::string& n,
                             const std::string& alpha,
                             const std::string& perm,
                             const std::string& dp,
                             const std::string& caseName,
                             const std::string& preciceRelTol )
{
  std::string name = eq_name + "-"
                     + preciceRelTol + "-"
                     + n + "-"
                     + alpha + "-" +
                     perm + "-"
                     + dp;
  if (caseName == "darcy-first")
  {
    name += "-darcy-first";
  }
  return name;
}

// the data of one directory, i.e. for multiple timesteps
//struct SimulationData
//{
//  std::vector<double> times;      //< time points of the parsed values
//  std::vector<CellData> values;   //< the values for a single time step, each, corresponding to times
//};

//void parseData(std::string filename, MeshData &meshData, double &t, const bool ignoreObstacle)
//{
//  // Create a reader
//  vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();

//  reader->SetFileName(filename.c_str());
//  reader->Update();


//  if (reader->GetNumberOfPointArrays() == 0)
//  {
//    std::cout << "Error: File \"" << filename << "\" contains no data." << std::endl;
//  }

//  vtkSmartPointer<vtkDataSet> dataSet = reader->GetOutputAsDataSet(0);
//  vtkSmartPointer<vtkFieldData> fieldData = dataSet->GetFieldData();

//  /*
//  std::cout << "file \"" << filename << "\" contains " << reader->GetNumberOfPointArrays() << " point arrays, "
//    << reader->GetNumberOfCellArrays() << "," << reader->GetNumberOfColumnArrays() << std::endl;
//  std::cout << " field has " << fieldData->GetNumberOfArrays() << " arrays.";
//  */

//  // get stored time
//  if (fieldData->GetNumberOfArrays() > 0)
//  {
//    vtkSmartPointer<vtkDataArray> array = fieldData->GetArray(0);
//    double *value = array->GetTuple(0);
//    t = *value;
//  }
//  else
//  {
//    std::cout << "Warning: File \"" << filename << "\" has no time." << std::endl;
//  }

//  // get stored velocity values
//  vtkSmartPointer<vtkPointData> pointData = dataSet->GetPointData();
//  bool hasTemperature = false;
//  for (int arrayNo = 0; arrayNo < pointData->GetNumberOfArrays(); arrayNo++)
//  {
//    vtkSmartPointer<vtkDoubleArray> array = vtkDoubleArray::SafeDownCast(pointData->GetAbstractArray(arrayNo));
//    int nEntries = array->GetNumberOfTuples();
//    //int nComponents = array->GetNumberOfComponents();
//    //std::cout << "  " << array->GetName() << ", nEntries: " << nEntries << ", nComponents: " << nComponents << std::endl;

//    if (std::string(array->GetName()) == "velocity")
//    {

//      meshData.u.resize(nEntries);
//      meshData.v.resize(nEntries);

//      // loop over values
//      for (int entryNo = 0; entryNo < nEntries; entryNo++)
//      {
//        std::array<double,3> values;
//        array->GetTuple(entryNo, values.data());
//        meshData.u[entryNo] = values[0];
//        meshData.v[entryNo] = values[1];

//        if ( std::isnan( meshData.u[entryNo] ) ) {
//          meshData.u[entryNo] = 0.;
//        }

//        if ( std::isnan( meshData.v[entryNo] ) ) {
//          meshData.v[entryNo] = 0.;
//        }
//      }
//      // Set values to zero that could have different values
//      // depending on how the values are set for obstacle cells
//      // Only works for Karman vortex test case and only for given mesh
//      if (ignoreObstacle)
//      {
//        //Hard coded for now
//        meshData.u[918] = 0.;
//        meshData.u[919] = 0.;

//        meshData.u[1019] = 0.;
//        meshData.u[1020] = 0.;
//        meshData.u[1021] = 0.;

//        meshData.u[1121] = 0.;
//        meshData.u[1122] = 0.;

//        meshData.v[918] = 0.;
//        meshData.v[919] = 0.;

//        meshData.v[1019] = 0.;
//        meshData.v[1020] = 0.;
//        meshData.v[1021] = 0.;

//        meshData.v[1121] = 0.;
//        meshData.v[1122] = 0.;
//      }
//    }


//    meshData.temperature.resize(nEntries, 0.);
//    if (std::string(array->GetName()) == "T")
//    {
//      hasTemperature = true;
//      // loop over values
//      for (int entryNo = 0; entryNo < nEntries; entryNo++)
//      {
//        std::array<double,1> values;
//        array->GetTuple(entryNo, values.data());
//        meshData.temperature[entryNo] = values[0];


//        if ( std::isnan( meshData.temperature[entryNo] ) ) {
//          meshData.temperature[entryNo] = 0.;
//        }
//      }
//    }
//  }
//  if ( hasTemperature == false )
//  {
//    std::cout << "File does not have temperature field!" << std::endl;
//  }
//}

//void parseDirectory(std::string directory, SimulationData &simulationData, const bool ignoreObstacle )
//{
//  // get files in directory
//  std::vector<std::string> filenames;

//  directory_iterator end_iter; // default construction yields past-the-end
//  for (directory_iterator iter(directory); iter != end_iter; ++iter)
//  {
//    if (is_directory(iter->status()))
//      continue;

//    if (is_regular_file(iter->status()))
//    {
//      std::string filename(iter->path().string());

//      if (filename.find(".vti") != std::string::npos)
//      {
//        filenames.push_back(filename);
//      }
//    }
//  }

//  // sort files by filename
//  std::sort(filenames.begin(), filenames.end());

//  //std::cout << "Parse data in \"" << directory << "\"." << std::endl;

//  // loop over filenames and parse files
//  for (std::vector<std::string>::iterator iter = filenames.begin(); iter != filenames.end(); iter++)
//  {
//    std::string filename = *iter;
//    //std::cout << "  File " << filename << std::endl;
//    double t;
//    simulationData.values.emplace_back();
//    parseData(filename, simulationData.values.back(), t, ignoreObstacle);
//    simulationData.times.push_back(t);
//  }
//}


std::string findLastVtu( const std::string& problemName, const std::string& directory )
{
  // get files in directory
  std::vector<std::string> filenames;

  directory_iterator end_iter; // default construction yields past-the-end
  for (directory_iterator iter(directory); iter != end_iter; ++iter)
  {
    if (is_directory(iter->status()))
      continue;

    if (is_regular_file(iter->status()))
    {
      std::string filename(iter->path().string());

      if (filename.find(".vtu") != std::string::npos
         && filename.find(problemName) != std::string::npos )
      {
        filenames.push_back(filename);
      }
    }
  }

  // sort files by filename
  std::sort(filenames.begin(), filenames.end());


//  std::cout << "Found files " << std::endl;
//  for (const auto& name : filenames) {
//    std::cout << "  " << name << std::endl;
//  }

  return filenames.back();
}

//double compareData(SimulationData &testData, SimulationData &referenceData )
//{
//  const bool output = false;


//  // loop over entries in test data
//  double differenceNormData = 0;
//  double differenceNormTData = 0;
//  int nFiles = testData.times.size();

//  struct {
//    int uvPointIndex = 0;
//    int uvFileIndex = 0;

//    int tPointIndex = 0;
//    int tFileIndex = 0;

////    double uDiff = 0.;
////    double vDiff = 0.;
//    double uvDiff = 0.;
//    double tDiff = 0.;

//  } differenceIndicator;

//  for (int i = 0; i < nFiles; i++)
//  {
//    double t = testData.times[i];
//    MeshData &testMeshData = testData.values[i];

//    if (output)
//      std::cout << "test t = " << t << std::endl;

//    // get corresponding data sets of reference Data
//    int firstReferenceIndex = 0;
//    int secondReferenceIndex = 0;
//    for (int j = 0; j < (int)referenceData.times.size(); j++)
//    {
//      if (referenceData.times[j] > t)
//      {
//        secondReferenceIndex = j;
//        break;
//      }
//      firstReferenceIndex = j;
//    }

//    double alpha = 0.0;

//    if (secondReferenceIndex == 0)
//    {
//      secondReferenceIndex = firstReferenceIndex;
//    }
//    else
//    {
//      alpha = (t - referenceData.times[firstReferenceIndex]) / (referenceData.times[secondReferenceIndex] - referenceData.times[firstReferenceIndex]);
//    }

//    if (output)
//    {
//      std::cout << "   corresponding reference datasets: " << std::endl
//        << "     t=" << referenceData.times[firstReferenceIndex] << " at index " << firstReferenceIndex << std::endl
//        << "     t=" << referenceData.times[secondReferenceIndex] << " at index " << secondReferenceIndex << ", alpha: " << alpha << std::endl;
//    }

//    // loop over values
//    double differenceNormDataset = 0;
//    double differenceNormTDataset = 0;
//    int nPoints = testMeshData.u.size();
//    for (int j = 0; j < nPoints; j++)
//    {
//      double testU = testMeshData.u[j];
//      double testV = testMeshData.v[j];
//      double testT = testMeshData.temperature[j];

//      double referenceU = (1.-alpha) * referenceData.values[firstReferenceIndex].u[j] + alpha * referenceData.values[secondReferenceIndex].u[j];
//      double referenceV = (1.-alpha) * referenceData.values[firstReferenceIndex].v[j] + alpha * referenceData.values[secondReferenceIndex].v[j];

//      double referenceT = (1.-alpha) * referenceData.values[firstReferenceIndex].temperature[j] + alpha * referenceData.values[secondReferenceIndex].temperature[j];

//      double differenceNorm = sqrt( (referenceU - testU)*(referenceU - testU) + (referenceV - testV)*(referenceV - testV));
//      const double temperatureNorm = sqrt( (referenceT - testT)*(referenceT - testT) );


//      if ( differenceIndicator.uvDiff < std::fabs(referenceU - testU) + std::fabs(referenceV - testV) )
//      {
//        differenceIndicator.uvFileIndex = i;
//        differenceIndicator.uvPointIndex = j;
//        differenceIndicator.uvDiff = std::fabs(referenceU - testU) + std::fabs(referenceV - testV) ;
////        differenceIndicator.tDiff = std::max(differenceIndicator.tDiff, std::fabs(referenceT - testT) );
//      }

//      if ( differenceIndicator.tDiff < std::fabs(referenceT - testT) )
//      {
//        differenceIndicator.tFileIndex = i;
//        differenceIndicator.tPointIndex = j;
//        differenceIndicator.tDiff = std::fabs(referenceT - testT);
//      }

//      //std::cout << "      j=" << j << "/" << nPoints << ", error: " << differenceNorm << std::endl;
//      differenceNormDataset += differenceNorm;
//      differenceNormTDataset += temperatureNorm;
//    }
//    differenceNormDataset /= nPoints;
//    differenceNormData += differenceNormDataset;

//    differenceNormTDataset /= nPoints;
//    differenceNormTData += differenceNormTDataset;
//  }
//  differenceNormData /= nFiles;
//  differenceNormTData /= nFiles;
//  std::cout << nFiles << " files with " << testData.values[0].u.size() << " entries each." << std::endl
//    << "average 2-norm error per velocity vector: " << differenceNormData << std::endl;

//  std::cout << "average 2-norm error per temperature: " << differenceNormTData << std::endl;

//  std::cout << "Biggest difference uv: " << "File " << differenceIndicator.uvFileIndex
//            << " at point " << differenceIndicator.uvPointIndex
//            << ", |u-u_ref|+|v-v_ref|: " << differenceIndicator.uvDiff << std::endl;

//  std::cout << "Biggest difference T: " << "File " << differenceIndicator.tFileIndex
//            << " at point " << differenceIndicator.tPointIndex
//            << ", |T-T_ref|: " << differenceIndicator.tDiff << std::endl;

//  if (fabs(differenceNormData) > 1e-4)
//  {
//    std::cout << "Velocity error is higher than tolerance of 1e-4!" << std::endl;
//  }

//  if (fabs(differenceNormTData) > 1e-4)
//  {
//    std::cout << "Temperature error is higher than tolerance of 1e-4!" << std::endl;
//  }

//  return differenceNormData;
//}



void FindAllData(vtkPolyData* polydata)
{
  std::cout << "Normals: " << polydata->GetPointData()->GetNormals() << std::endl;

  vtkIdType numberOfPointArrays = polydata->GetPointData()->GetNumberOfArrays();
  std::cout << "Number of PointData arrays: " << numberOfPointArrays << std::endl;

  vtkIdType numberOfCellArrays = polydata->GetCellData()->GetNumberOfArrays();
  std::cout << "Number of CellData arrays: " << numberOfCellArrays << std::endl;

  std::cout << "Type table/key: " << std::endl;;
  //more values can be found in <VTK_DIR>/Common/vtkSetGet.h
  std::cout << VTK_UNSIGNED_CHAR << " unsigned char" << std::endl;
  std::cout << VTK_UNSIGNED_INT << " unsigned int" << std::endl;
  std::cout << VTK_FLOAT << " float" << std::endl;
  std::cout << VTK_DOUBLE << " double" << std::endl;

  for(vtkIdType i = 0; i < numberOfPointArrays; i++)
  {
    // The following two lines are equivalent
    //arrayNames.push_back(polydata->GetPointData()->GetArray(i)->GetName());
    //arrayNames.push_back(polydata->GetPointData()->GetArrayName(i));
    int dataTypeID = polydata->GetPointData()->GetArray(i)->GetDataType();
    std::cout << "Array " << i << ": " << polydata->GetPointData()->GetArrayName(i)
              << " (type: " << dataTypeID << ")" << std::endl;
  }

  for(vtkIdType i = 0; i < numberOfCellArrays; i++)
  {
    // The following two lines are equivalent
    //polydata->GetPointData()->GetArray(i)->GetName();
    //polydata->GetPointData()->GetArrayName(i);
    int dataTypeID = polydata->GetCellData()->GetArray(i)->GetDataType();
    std::cout << "Array " << i << ": " << polydata->GetCellData()->GetArrayName(i)
              << " (type: " << dataTypeID << ")" << std::endl;
  }
}

void printData( const vtkUnstructuredGrid* unstructuredGrid )
{
//  unstructuredGrid
//unstructuredGrid->Get
}

void checkCellCenters(const CellData& cellDataA,
                      const CellData& cellDataB )
{
//  if ( pda->Get )
//vtkPolyData* asdf;
//asdf->GetPo
  if (cellDataA.size != cellDataB.size )
  {
    throw std::runtime_error("Datasets have different sizes!");
  }
  constexpr double prec = 1e-8;

  double normLinf = std::numeric_limits<double>::min();
  double normLinfB = std::numeric_limits<double>::min();

  for (size_t i = 0; i< cellDataA.size; ++i)
  {
    const auto& pa = cellDataA.cellCenters[i];
    const auto& pb = cellDataB.cellCenters[i];

    double normLinfLocal = std::numeric_limits<double>::min();
    for (int d = 0; d < 3; ++d)
    {
      normLinfLocal = std::max( std::fabs(pa[d] - pb[d]), normLinfLocal );
      normLinfB = std::max( std::fabs(pb[d]), normLinfB);
    }
    normLinf = std::max( normLinfLocal, normLinf );
    if ( normLinfLocal > prec )
    {
      std::cout << "Point " << i << "deviates!" << std::endl
                <<"  Dataset A: " << pa[0] << ", " << pa[1] << ", " << pa[2] << std::endl
                <<"  Dataset B: " << pb[0] << ", " << pb[1] << ", " << pb[2] << std::endl;
    }
  }
  const double normLinfRel = normLinf / normLinfB;
//  std::cout << "Coordinate deviations: " << std::endl
//            << "  linf-norm: " << normLinf << " (abs), "
//            << normLinfRel << "(rel)" << std::endl;

  if ( normLinfRel > prec )
  {
    throw std::runtime_error("Coordinates deviate too much!");
  }

}

vtkPolyData* getPolyDataFromFile( const std::string& filename )
{
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();

  reader->SetFileName(filename.c_str());
  reader->Update();

  const auto grid = reader->GetOutput();

  vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter =
      vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
  surfaceFilter->SetInputData(grid);
  surfaceFilter->Update();

  return surfaceFilter->GetOutput();
}

void parseData(const std::string& filename,
               CellData& cellData)
{

  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();

  reader->SetFileName(filename.c_str());
  reader->Update();

  const auto grid = reader->GetOutput();

  vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter =
      vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
  surfaceFilter->SetInputData(grid);
  surfaceFilter->Update();

  surfaceFilter->GetOutput();
  vtkPolyData* polydata = surfaceFilter->GetOutput();

  std::cout << "Output has " << polydata->GetNumberOfPoints() << " points." << std::endl;


//  FindAllData( polydata );


  vtkIdType pressureId = -1;
  vtkIdType velocityId = -1;
  {
    const vtkIdType numberOfCellArrays = polydata->GetCellData()->GetNumberOfArrays();
    for(vtkIdType i = 0; i < numberOfCellArrays; i++)
    {
      // The following two lines are equivalent
      //polydata->GetPointData()->GetArray(i)->GetName();
      //polydata->GetPointData()->GetArrayName(i);
      const std::string arrayName = polydata->GetCellData()->GetArrayName(i);
      if ( arrayName == "p" )
        pressureId = i;
      else if ( arrayName == "velocity_liq (m/s)" )
        velocityId = i;
//      int dataTypeID = polydata->GetCellData()->GetArray(i)->GetDataType();
//      std::cout << "Array " << i << ": " << polydata->GetCellData()->GetArrayName(i)
//                << " (type: " << dataTypeID << ")" << std::endl;
    }
  }
  assert( pressureId != -1 );
  assert( velocityId != -1 );

  //Parse data
  {
    const auto pressureArray = polydata->GetCellData()->GetArray(pressureId);
    const auto velocityArray = polydata->GetCellData()->GetArray(velocityId);

    const auto pressureArraySize = pressureArray->GetNumberOfTuples();
    const auto velocityArraySize = velocityArray->GetNumberOfTuples();
    assert( velocityArraySize == pressureArraySize );
    cellData.resize( velocityArraySize );

    std::array<double, 3> velocityVal;
    std::array<double, 1> pressureVal;

    for (vtkIdType i = 0; i < pressureArraySize; ++i)
    {
      pressureArray->GetTuple( i, pressureVal.data() );
      velocityArray->GetTuple( i, velocityVal.data() );

      cellData.p[i] = pressureVal[0];

      cellData.u[i] = velocityVal[0];
      cellData.v[i] = velocityVal[1];


      std::array<double, 6> bounds;
      polydata->GetCellBounds( i, bounds.data() );
      for (int j = 0; j < 3; ++j) {
        cellData.cellCenters[i][j] = bounds[j*2] + bounds[j*2+1];
        cellData.cellCenters[i][j] *= 0.5;
      }


    }

  }
}


//Errors computeError(const CellData& dataSetA,
//                    const CellData& dataSetB)
//{
//  assert( dataSetA.size == dataSetB.size );
//  Errors errorsVelocity;
//  Errors errorsPressure;


//  ErrorsBase errV;
//  ErrorsBase errP;
//  for (size_t i = 0; i < dataSetA.size; ++i)
//  {
//    errorsPressure.abs.l1 += std::fabs(dataSetA.p[i] - dataSetB.p[i]);
//    errorsVelocity.abs.l1 += std::fabs(dataSetA.u[i] - dataSetB.u[i])
//                         + std::fabs(dataSetA.v[i] - dataSetB.v[i]) ;

//    errP.l1 += std::fabs( dataSetB.p[i] );
//    errV.l1 += std::fabs( dataSetB.u[i] ) + std::fabs( dataSetB.v[i] );

//    errorsPressure.abs.l2 += std::pow(dataSetA.p[i] - dataSetB.p[i], 2);
//    errorsVelocity.abs.l2 += std::pow(dataSetA.u[i] - dataSetB.u[i], 2)
//                         + std::pow(dataSetA.v[i] - dataSetB.v[i], 2) ;

//    errP.l2 += std::pow( dataSetB.p[i], 2 );
//    errV.l2 += std::pow( dataSetB.u[i], 2 ) + std::pow( dataSetB.v[i], 2 );

//    errorsPressure.abs.linf = std::max( std::fabs(dataSetA.p[i] - dataSetB.p[i] ),
//                                   errorsPressure.abs.linf );

//    const double u = std::fabs(dataSetA.u[i] - dataSetB.u[i]);
//    const double v = std::fabs(dataSetA.v[i] - dataSetB.v[i]);
//    const double uvMax = std::max( u, v );
//    errorsVelocity.abs.linf = std::max( uvMax, errorsVelocity.abs.linf );

//    errP.linf = std::max( std::fabs(dataSetB.p[i]), errP.linf );
//    errV.linf = std::max( std::fabs(dataSetB.u[i]), errV.linf );
//    errV.linf = std::max( std::fabs(dataSetB.v[i]), errV.linf );

//    constexpr double tol = 1e-10;
//    if ( uvMax > tol )
//    {
//      std::cout << "Large error " << uvMax << " at Point " << i << std::endl
//                << "  Coords: " << dataSetA.cellCenters[i][0] << ", "
//                << dataSetA.cellCenters[i][1] << ", "
//                << dataSetA.cellCenters[i][2] << ", "
//                << std::endl;
//    }
//  }
//  errorsPressure.abs.l2 = std::sqrt(errorsPressure.abs.l2);
//  errorsVelocity.abs.l2 = std::sqrt(errorsVelocity.abs.l2);

//  errP.l2 = std::sqrt( errP.l2 );
//  errV.l2 = std::sqrt( errV.l2 );

//  errorsPressure.rel.l1 = errorsPressure.abs.l1 / errP.l1;
//  errorsPressure.rel.l2 = errorsPressure.abs.l2 / errP.l2;
//  errorsPressure.rel.linf = errorsPressure.abs.linf / errP.linf;

//  errorsVelocity.rel.l1 = errorsVelocity.abs.l1 / errV.l1;
//  errorsVelocity.rel.l2 = errorsVelocity.abs.l2 / errV.l2;
//  errorsVelocity.rel.linf = errorsVelocity.abs.linf / errV.linf;


//  std::cout << "Errors: " << std::endl
//            << "  p: " << std::endl
//            << "    " << errorsPressure.abs.l1 << " (l1), " << errorsPressure.rel.l1 << ", " << errP.l1 << std::endl
//            << "    " << errorsPressure.abs.l2 << " (l2), " << errorsPressure.rel.l2 << ", " << errP.l2 << std::endl
//            << "    " << errorsPressure.abs.linf << " (linf), " << errorsPressure.rel.linf << ", " << errP.linf << std::endl;
//  std::cout << "Errors: " << std::endl
//            << "  vel: " << std::endl
//            << "    " << errorsVelocity.abs.l1 << " (l1), " << errorsVelocity.rel.l1 << ", " << errV.l1 << std::endl
//            << "    " << errorsVelocity.abs.l2 << " (l2), " << errorsVelocity.rel.l2 << ", " << errV.l2 << std::endl
//            << "    " << errorsVelocity.abs.linf << " (linf), " << errorsVelocity.rel.linf << ", " << errV.linf << std::endl;

////  return errors;
//}

Errors computePressureError(const CellData& dataSetA,
                    const CellData& dataSetB)
{
  assert( dataSetA.size == dataSetB.size );
  Errors errorsPressure;
  ErrorsBase errP;

  for (size_t i = 0; i < dataSetA.size; ++i)
  {
    errorsPressure.abs.l1 += std::fabs(dataSetA.p[i] - dataSetB.p[i]);

    errP.l1 += std::fabs( dataSetB.p[i] );

    errorsPressure.abs.l2 += std::pow(dataSetA.p[i] - dataSetB.p[i], 2);
    errP.l2 += std::pow( dataSetB.p[i], 2 );

    errorsPressure.abs.linf = std::max( std::fabs(dataSetA.p[i] - dataSetB.p[i] ),
                                       errorsPressure.abs.linf );

    errP.linf = std::max( std::fabs(dataSetB.p[i]), errP.linf );

    constexpr double tol = 1e-10;
    if ( std::fabs(dataSetA.p[i] - dataSetB.p[i] ) > tol )
    {
//      std::cout << "Large error " << std::fabs(dataSetA.p[i] - dataSetB.p[i] ) << " (p) at Point " << i << std::endl
//                << "  Coords: " << dataSetA.cellCenters[i][0] << ", "
//                << dataSetA.cellCenters[i][1] << ", "
//                << dataSetA.cellCenters[i][2] << ", "
//                << std::endl;
    }
  }
  errorsPressure.abs.l2 = std::sqrt(errorsPressure.abs.l2);

  errP.l2 = std::sqrt( errP.l2 );

  errorsPressure.rel.l1 = errorsPressure.abs.l1 / errP.l1;
  errorsPressure.rel.l2 = errorsPressure.abs.l2 / errP.l2;
  errorsPressure.rel.linf = errorsPressure.abs.linf / errP.linf;

  std::cout << "Errors: " << std::endl
            << "  p: " << std::endl
            << "    " << errorsPressure.abs.l1 << " (l1), " << errorsPressure.rel.l1 << ", " << errP.l1 << std::endl
            << "    " << errorsPressure.abs.l2 << " (l2), " << errorsPressure.rel.l2 << ", " << errP.l2 << std::endl
            << "    " << errorsPressure.abs.linf << " (linf), " << errorsPressure.rel.linf << ", " << errP.linf << std::endl;

  return errorsPressure;
}

Errors computeVelocityError(const CellData& dataSetA,
                            const CellData& dataSetB)
{
  assert( dataSetA.size == dataSetB.size );
  Errors errorsVelocity;
  ErrorsBase errV;
  for (size_t i = 0; i < dataSetA.size; ++i)
  {
    errorsVelocity.abs.l1 += std::fabs(dataSetA.u[i] - dataSetB.u[i])
                             + std::fabs(dataSetA.v[i] - dataSetB.v[i]) ;
    errV.l1 += std::fabs( dataSetB.u[i] ) + std::fabs( dataSetB.v[i] );

    errorsVelocity.abs.l2 += std::pow(dataSetA.u[i] - dataSetB.u[i], 2)
                             + std::pow(dataSetA.v[i] - dataSetB.v[i], 2) ;

    errV.l2 += std::pow( dataSetB.u[i], 2 ) + std::pow( dataSetB.v[i], 2 );

    const double u = std::fabs(dataSetA.u[i] - dataSetB.u[i]);
    const double v = std::fabs(dataSetA.v[i] - dataSetB.v[i]);
    const double uvMax = std::max( u, v );
    errorsVelocity.abs.linf = std::max( uvMax, errorsVelocity.abs.linf );

    errV.linf = std::max( std::fabs(dataSetB.u[i]), errV.linf );
    errV.linf = std::max( std::fabs(dataSetB.v[i]), errV.linf );

    constexpr double tol = 1e-10;
    if ( uvMax > tol )
    {
//      std::cout << "Large error " << uvMax << " at Point " << i << std::endl
//                << "  Coords: " << dataSetA.cellCenters[i][0] << ", "
//                << dataSetA.cellCenters[i][1] << ", "
//                << dataSetA.cellCenters[i][2] << ", "
//                << std::endl;
    }
  }
  errorsVelocity.abs.l2 = std::sqrt(errorsVelocity.abs.l2);

  errV.l2 = std::sqrt( errV.l2 );

  errorsVelocity.rel.l1 = errorsVelocity.abs.l1 / errV.l1;
  errorsVelocity.rel.l2 = errorsVelocity.abs.l2 / errV.l2;
  errorsVelocity.rel.linf = errorsVelocity.abs.linf / errV.linf;

  std::cout << "Errors: " << std::endl
            << "  vel: " << std::endl
            << "    " << errorsVelocity.abs.l1 << " (l1), " << errorsVelocity.rel.l1 << ", " << errV.l1 << std::endl
            << "    " << errorsVelocity.abs.l2 << " (l2), " << errorsVelocity.rel.l2 << ", " << errV.l2 << std::endl
            << "    " << errorsVelocity.abs.linf << " (linf), " << errorsVelocity.rel.linf << ", " << errV.linf << std::endl;

  return errorsVelocity;
}

std::tuple<AllErrors, AllErrors> computeErrors(const ParsedData& dataA,
                     const ParsedData& dataB)
{
  AllErrors pressureErrors;
  AllErrors velocityErrors;
  //Stokes error
  {
    checkCellCenters( dataA.stokes, dataB.stokes );
    pressureErrors.stokes = computePressureError( dataA.stokes, dataB.stokes );
    velocityErrors.stokes = computeVelocityError( dataA.stokes, dataB.stokes );
  }

  //Darcy error
  {
    checkCellCenters( dataA.darcy, dataB.darcy );
    pressureErrors.darcy = computePressureError( dataA.darcy, dataB.darcy );
    velocityErrors.darcy = computeVelocityError( dataA.darcy, dataB.darcy );
  }

  return std::make_tuple( pressureErrors, velocityErrors );
}

int main(int argc, char *argv[])
{
  // parse command line arguments
//  if (argc != 4)
//  {
//    std::cout << "usage: " << argv[0] << " <testcasename> <directory to test> <directory with reference data>" << std::endl;
//    exit(-1);
//  }normLinf

  const std::vector<std::string> couplingTypes = {"serial-implicit", "parallel-implicit"};
  const std::vector<std::string> meshSizes = {"20", "40", "80" };
  const std::vector<std::string> pressureDifferences = { "1e-9" };
  const std::vector<std::string> permeabilities={"1e-6"};
  const std::vector<std::string> alpha_BJS={"1.0"};
  const std::vector<bool> hasInertiaTerms={true, false};
  const std::vector<std::string> preciceRelTols={"1e-2", "1e-4", "1e-6", "1e-8"};
  const std::vector<std::string> preciceCases={"serial-implicit", "serial-implicit-darcy-first"};
//  const std::vector<std::string> interfaceCases={"darcy", "stokes"};


//  CellData stokesMonolithicData;
//  ParsedData monolithicData;
//  {
//    const std::string stokesMonolithicFilename = "monolithic-stokes.vtu";
//    const std::string darcyMonolithicFilename = "monolithic-darcy.vtu";
//    parseData( stokesMonolithicFilename, monolithicData.stokes );
//    parseData( darcyMonolithicFilename, monolithicData.darcy );
//  }


//  ParsedData iterativeData;
//  {
//    const std::string stokesIterativeFilename = "iterative-stokes.vtu";
//    const std::string darcyIterativeFilename = "iterative-darcy.vtu";
//    parseData( stokesIterativeFilename, iterativeData.stokes );
//    parseData( darcyIterativeFilename, iterativeData.darcy );
//  }

//  checkCellCenters( monolithicData.stokes, iterativeData.stokes );
//  checkCellCenters( monolithicData.darcy, iterativeData.darcy );

//  AllErrors pressureErrors;
//  AllErrors velocityErrors;
//  std::tie( pressureErrors, velocityErrors ) = computeErrors( monolithicData, iterativeData );

  //Read file and process files
  const std::string monolithicRootDir = "../../results/monolithic";
  const std::string iterativeRootDir = "../../results/iterative";

  {
    for (const auto& couplingType: couplingTypes)
    {
      ParsedData monolithicData;
      ParsedData iterativeData;
      for (const auto isNavierStokes: hasInertiaTerms)
      {
        const std::string equationName = (isNavierStokes) ? "navier-stokes" : "stokes";
        for (const auto& meshSize : meshSizes)
        {
          const double h = 1./ std::atof(meshSize.c_str());
          for (const auto& alpha : alpha_BJS)
          {
            for (const auto& permeability : permeabilities)
            {
              for (const auto& dp : pressureDifferences)
              {
                //Parse
                {
                  const std::string& monolithicDir = monolithicRootDir + "/" + getMonolithicName(equationName,
                                                                                                 meshSize,
                                                                                                 alpha,
                                                                                                 permeability,
                                                                                                 dp);
                  std::cout << "Parsing monolithic data from " << monolithicDir << std::endl;
                  const std::string& flowFilename = findLastVtu( "_" + equationName, monolithicDir );
                  std::cout << "Loading monolithic free-flow data from: " << std::endl
                            << "  " << flowFilename << std::endl;
                  parseData( flowFilename, monolithicData.stokes );
                  const std::string& darcyFilename = findLastVtu( "_darcy", monolithicDir );
                  std::cout << "Loading monolithic porous-media-flow data from: " << std::endl
                            << "  " << darcyFilename << std::endl;
                  parseData( darcyFilename, monolithicData.darcy );
                }


                for ( const auto& preciceCase: preciceCases )
                {
                  const std::string caseName = (preciceCase == "serial-implicit") ? "stokes-first" : "darcy-first";
                  std::ofstream velocityErrorFile;
                  std::ofstream pressureErrorFile;
                  {
                    const std::string monolithicString = getMonolithicName(equationName,
                                                                           meshSize,
                                                                           alpha,
                                                                           permeability,
                                                                           dp);
                    {
                      const std::string resultFilename = couplingType + "-errors-u-" + monolithicString
                                                         + "-" + caseName + ".csv";
                      velocityErrorFile.open( resultFilename,  std::ofstream::trunc);

                      velocityErrorFile << "precice_rel_tol" << ","
                                        << "err_u_stokes_l1_abs" << "," << "err_stokes_l1_rel" << ","
                                        << "err_u_stokes_l2_abs" << "," << "err_u_stokes_l2_rel" << ","
                                        << "err_u_stokes_linf_abs" << "," << "err_u_stokes_linf_rel" << ","
                                        << "err_u_darcy_l1_abs" << "," << "err_u_darcy_l1_rel" << ","
                                        << "err_u_darcy_l2_abs" << "," << "err_u_darcy_l2_rel" << ","
                                        << "err_u_darcy_linf_abs" << "," << "err_u_darcy_linf_rel" << std::endl;
                    }

                    {
                      const std::string resultFilename = couplingType + "-errors-p-" + monolithicString
                                                         + "-" + caseName + ".csv";
                      pressureErrorFile.open( resultFilename,  std::ofstream::trunc);

                      pressureErrorFile << "precice_rel_tol" << ","
                                        << "err_p_stokes_l1_abs" << "," << "err_stokes_l1_rel" << ","
                                        << "err_p_stokes_l2_abs" << "," << "err_p_stokes_l2_rel" << ","
                                        << "err_p_stokes_linf_abs" << "," << "err_p_stokes_linf_rel" << ","
                                        << "err_p_darcy_l1_abs" << "," << "err_p_darcy_l1_rel" << ","
                                        << "err_p_darcy_l2_abs" << "," << "err_p_darcy_l2_rel" << ","
                                        << "err_p_darcy_linf_abs" << "," << "err_p_darcy_linf_rel" << std::endl;
                    }
                  }
                  for (const auto& preciceRelTol: preciceRelTols )
                  {
                    {
                      const std::string& iterativeDir = iterativeRootDir + "/" + couplingType + "/" + getIterativeName(equationName,
                                                                                                                       meshSize,
                                                                                                                       alpha,
                                                                                                                       permeability,
                                                                                                                       dp,
                                                                                                                       caseName,
                                                                                                                       preciceRelTol );

                      std::cout << "Parsing iterative data from " << iterativeDir << std::endl;
                      const std::string& flowFilename = findLastVtu( equationName+"-iterative", iterativeDir );
                      std::cout << "Loading monolithic free-flow data from: " << std::endl
                                << "  " << flowFilename << std::endl;
                      parseData( flowFilename, iterativeData.stokes );
                      const std::string& darcyFilename = findLastVtu( "darcy-iterative", iterativeDir );
                      std::cout << "Loading monolithic porous-media-flow data from: " << std::endl
                                << "  " << darcyFilename << std::endl;
                      parseData( darcyFilename, iterativeData.darcy );

                      AllErrors pressureErrors;
                      AllErrors velocityErrors;
                      std::tie( pressureErrors, velocityErrors ) = computeErrors( monolithicData, iterativeData );

                      pressureErrorFile << std::setprecision(outputPrecision) << std::scientific;
                      pressureErrorFile << preciceRelTol << ","
                                        << pressureErrors.stokes.abs.l1 << "," << pressureErrors.stokes.rel.l1 << ","
                                        << pressureErrors.stokes.abs.l2 << "," << pressureErrors.stokes.rel.l2 << ","
                                        << pressureErrors.stokes.abs.linf << "," << pressureErrors.stokes.rel.linf << ","
                                        << pressureErrors.darcy.abs.l1 << "," << pressureErrors.darcy.rel.l1 << ","
                                        << pressureErrors.darcy.abs.l2 << "," << pressureErrors.darcy.rel.l2 << ","
                                        << pressureErrors.darcy.abs.linf << "," << pressureErrors.darcy.rel.linf << std::endl;

                      velocityErrorFile << std::setprecision(outputPrecision) << std::scientific;
                      velocityErrorFile << preciceRelTol << ","
                                        << velocityErrors.stokes.abs.l1 << "," << velocityErrors.stokes.rel.l1 << ","
                                        << velocityErrors.stokes.abs.l2 << "," << velocityErrors.stokes.rel.l2 << ","
                                        << velocityErrors.stokes.abs.linf << "," << velocityErrors.stokes.rel.linf << ","
                                        << velocityErrors.darcy.abs.l1 << "," << velocityErrors.darcy.rel.l1 << ","
                                        << velocityErrors.darcy.abs.l2 << "," << velocityErrors.darcy.rel.l2 << ","
                                        << velocityErrors.darcy.abs.linf << "," << velocityErrors.darcy.rel.linf << std::endl;
                    }
                  }
                  velocityErrorFile.close();
                  pressureErrorFile.close();
                }
              }

            }
          }
        }
      }
    }
  }

  //Write file
  {
  }

  return EXIT_SUCCESS;
}
