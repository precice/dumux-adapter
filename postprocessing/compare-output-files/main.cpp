//#include <vtkXMLUnstructuredGrid.h>
#include <vtkImageData.h>
#include <vtkStructuredGrid.h>
//#include <vtkXMLStructuredGridWriter.h>
#include <vtkCellData.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkDoubleArray.h>
#include <vtkHexahedron.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <cmath>

#include <iomanip>
#include <iostream>
// the data for one timestep
struct CellData {
    std::vector<double> u;
    std::vector<double> v;
    std::vector<double> p;

    std::vector<std::array<double, 3> > cellCenters;

    size_t size;

    template<typename T>
    void resize(const T s)
    {
        assert(s > -1);
        size = size_t(s);
        u.resize(size);
        v.resize(size);
        p.resize(size);

        cellCenters.resize(size);
    }
};

struct ErrorsBase {
    double l1;
    double l2;
    double linf;

    ErrorsBase() : l1(0.), l2(0.), linf(0.) {}
};

struct Errors {
    ErrorsBase abs;
    ErrorsBase rel;

    ErrorsBase dataset_norm;
};

void checkCellCenters(const CellData &cellDataA, const CellData &cellDataB)
{
    //  if ( pda->Get )
    //vtkPolyData* asdf;
    //asdf->GetPo
    if (cellDataA.size != cellDataB.size) {
        throw std::runtime_error("Datasets have different sizes!");
    }
    constexpr double prec = 1e-8;

    double normLinf = std::numeric_limits<double>::min();
    double normLinfB = std::numeric_limits<double>::min();

    for (size_t i = 0; i < cellDataA.size; ++i) {
        const auto &pa = cellDataA.cellCenters[i];
        const auto &pb = cellDataB.cellCenters[i];

        double normLinfLocal = std::numeric_limits<double>::min();
        for (int d = 0; d < 3; ++d) {
            normLinfLocal = std::max(std::fabs(pa[d] - pb[d]), normLinfLocal);
            normLinfB = std::max(std::fabs(pb[d]), normLinfB);
        }
        normLinf = std::max(normLinfLocal, normLinf);
        if (normLinfLocal > prec) {
            std::cout << "Point " << i << "deviates!" << std::endl
                      << "  Dataset A: " << pa[0] << ", " << pa[1] << ", "
                      << pa[2] << std::endl
                      << "  Dataset B: " << pb[0] << ", " << pb[1] << ", "
                      << pb[2] << std::endl;
        }
    }
    const double normLinfRel = normLinf / normLinfB;
    //  std::cout << "Coordinate deviations: " << std::endl
    //            << "  linf-norm: " << normLinf << " (abs), "
    //            << normLinfRel << "(rel)" << std::endl;

    if (normLinfRel > prec) {
        throw std::runtime_error("Coordinates deviate too much!");
    }
}

vtkPolyData *getPolyDataFromFile(const std::string &filename)
{
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();

    reader->SetFileName(filename.c_str());
    reader->Update();

    const auto grid = reader->GetOutput();

    vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter =
        vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surfaceFilter->SetInputData(grid);
    surfaceFilter->Update();

    return surfaceFilter->GetOutput();
}

void parseData(const std::string &filename, CellData &cellData)
{
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();

    reader->SetFileName(filename.c_str());
    reader->Update();

    const auto grid = reader->GetOutput();

    vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter =
        vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surfaceFilter->SetInputData(grid);
    surfaceFilter->Update();

    surfaceFilter->GetOutput();
    vtkPolyData *polydata = surfaceFilter->GetOutput();

    //std::cout << "Output has " << polydata->GetNumberOfPoints() << " points."
    //          << std::endl;

    vtkIdType pressureId = -1;
    vtkIdType velocityId = -1;
    {
        const vtkIdType numberOfCellArrays =
            polydata->GetCellData()->GetNumberOfArrays();
        for (vtkIdType i = 0; i < numberOfCellArrays; i++) {
            // The following two lines are equivalent
            //polydata->GetPointData()->GetArray(i)->GetName();
            //polydata->GetPointData()->GetArrayName(i);
            const std::string arrayName =
                polydata->GetCellData()->GetArrayName(i);
            if (arrayName == "p")
                pressureId = i;
            else if (arrayName == "velocity_liq (m/s)")
                velocityId = i;
            //      int dataTypeID = polydata->GetCellData()->GetArray(i)->GetDataType();
            //      std::cout << "Array " << i << ": " << polydata->GetCellData()->GetArrayName(i)
            //                << " (type: " << dataTypeID << ")" << std::endl;
        }
    }
    assert(pressureId != -1);
    assert(velocityId != -1);

    //Parse data
    {
        const auto pressureArray =
            polydata->GetCellData()->GetArray(pressureId);
        const auto velocityArray =
            polydata->GetCellData()->GetArray(velocityId);

        const auto pressureArraySize = pressureArray->GetNumberOfTuples();
        const auto velocityArraySize = velocityArray->GetNumberOfTuples();
        assert(velocityArraySize == pressureArraySize);
        cellData.resize(velocityArraySize);

        std::array<double, 3> velocityVal;
        std::array<double, 1> pressureVal;

        for (vtkIdType i = 0; i < pressureArraySize; ++i) {
            pressureArray->GetTuple(i, pressureVal.data());
            velocityArray->GetTuple(i, velocityVal.data());

            cellData.p[i] = pressureVal[0];

            cellData.u[i] = velocityVal[0];
            cellData.v[i] = velocityVal[1];

            std::array<double, 6> bounds;
            polydata->GetCellBounds(i, bounds.data());
            for (int j = 0; j < 3; ++j) {
                cellData.cellCenters[i][j] = bounds[j * 2] + bounds[j * 2 + 1];
                cellData.cellCenters[i][j] *= 0.5;
            }
        }
    }
}

Errors computePressureError(const CellData &dataSetA, const CellData &dataSetB)
{
    assert(dataSetA.size == dataSetB.size);
    Errors pressureErrors;
    ErrorsBase errP;

    for (size_t i = 0; i < dataSetA.size; ++i) {
        pressureErrors.abs.l1 += std::fabs(dataSetA.p[i] - dataSetB.p[i]);

        errP.l1 += std::fabs(dataSetB.p[i]);

        pressureErrors.abs.l2 += std::pow(dataSetA.p[i] - dataSetB.p[i], 2);
        errP.l2 += std::pow(dataSetB.p[i], 2);

        pressureErrors.abs.linf = std::max(
            std::fabs(dataSetA.p[i] - dataSetB.p[i]), pressureErrors.abs.linf);

        errP.linf = std::max(std::fabs(dataSetB.p[i]), errP.linf);

        constexpr double tol = 1e-10;
        if (std::fabs(dataSetA.p[i] - dataSetB.p[i]) > tol) {
            //      std::cout << "Large error " << std::fabs(dataSetA.p[i] - dataSetB.p[i] ) << " (p) at Point " << i << std::endl
            //                << "  Coords: " << dataSetA.cellCenters[i][0] << ", "
            //                << dataSetA.cellCenters[i][1] << ", "
            //                << dataSetA.cellCenters[i][2] << ", "
            //                << std::endl;
        }
    }
    pressureErrors.abs.l2 = std::sqrt(pressureErrors.abs.l2);

    errP.l2 = std::sqrt(errP.l2);

    pressureErrors.rel.l1 = pressureErrors.abs.l1 / errP.l1;
    pressureErrors.rel.l2 = pressureErrors.abs.l2 / errP.l2;
    pressureErrors.rel.linf = pressureErrors.abs.linf / errP.linf;

    pressureErrors.dataset_norm.l1 = errP.l1;
    pressureErrors.dataset_norm.l2 = errP.l2;
    pressureErrors.dataset_norm.linf = errP.linf;

    //    std::cout << "Errors: " << std::endl
    //              << "  p: " << std::endl
    //              << "     (abs)        (rel)       (norm dataset 2)" << std::endl
    //              << "    " << pressureErrors.abs.l1 << " (l1), "
    //              << pressureErrors.rel.l1 << ", " << errP.l1 << std::endl
    //              << "    " << pressureErrors.abs.l2 << " (l2), "
    //              << pressureErrors.rel.l2 << ", " << errP.l2 << std::endl
    //              << "    " << pressureErrors.abs.linf << " (linf), "
    //              << pressureErrors.rel.linf << ", " << errP.linf << std::endl;

    //    std::cout << "Errors (inf-norm): " << std::endl
    //              << std::setw(8) << "Pressure " << std::setw(14) << "(abs)" << std::setw(14) << "(rel)" << std::setw(20) << "(norm of dataset)" << std::endl
    //              << std::setw(8) << " " << std::setw(14) << pressureErrors.abs.linf <<  std::setw(14) << pressureErrors.rel.linf << std::setw(20) << errP.linf << std::endl;

    return std::move(pressureErrors);
}

Errors computeVelocityError(const CellData &dataSetA, const CellData &dataSetB)
{
    assert(dataSetA.size == dataSetB.size);
    Errors velocityErrors;
    ErrorsBase errV;
    for (size_t i = 0; i < dataSetA.size; ++i) {
        velocityErrors.abs.l1 += std::fabs(dataSetA.u[i] - dataSetB.u[i]) +
                                 std::fabs(dataSetA.v[i] - dataSetB.v[i]);
        errV.l1 += std::fabs(dataSetB.u[i]) + std::fabs(dataSetB.v[i]);

        velocityErrors.abs.l2 += std::pow(dataSetA.u[i] - dataSetB.u[i], 2) +
                                 std::pow(dataSetA.v[i] - dataSetB.v[i], 2);

        errV.l2 += std::pow(dataSetB.u[i], 2) + std::pow(dataSetB.v[i], 2);

        const double u = std::fabs(dataSetA.u[i] - dataSetB.u[i]);
        const double v = std::fabs(dataSetA.v[i] - dataSetB.v[i]);
        const double uvMax = std::max(u, v);
        velocityErrors.abs.linf = std::max(uvMax, velocityErrors.abs.linf);

        errV.linf = std::max(std::fabs(dataSetB.u[i]), errV.linf);
        errV.linf = std::max(std::fabs(dataSetB.v[i]), errV.linf);

        constexpr double tol = 1e-10;
        if (uvMax > tol) {
            //      std::cout << "Large error " << uvMax << " at Point " << i << std::endl
            //                << "  Coords: " << dataSetA.cellCenters[i][0] << ", "
            //                << dataSetA.cellCenters[i][1] << ", "
            //                << dataSetA.cellCenters[i][2] << ", "
            //                << std::endl;
        }
    }
    velocityErrors.abs.l2 = std::sqrt(velocityErrors.abs.l2);

    errV.l2 = std::sqrt(errV.l2);

    velocityErrors.rel.l1 = velocityErrors.abs.l1 / errV.l1;
    velocityErrors.rel.l2 = velocityErrors.abs.l2 / errV.l2;
    velocityErrors.rel.linf = velocityErrors.abs.linf / errV.linf;

    velocityErrors.dataset_norm.l1 = errV.l1;
    velocityErrors.dataset_norm.l2 = errV.l2;
    velocityErrors.dataset_norm.linf = errV.linf;

    //    std::cout << "Errors: " << std::endl
    //              << "  vel: " << std::endl
    //              << "     (abs)        (rel)       (norm dataset 2)" << std::endl
    //              << "    " << velocityErrors.abs.l1 << " (l1), "
    //              << velocityErrors.rel.l1 << ", " << errV.l1 << std::endl
    //              << "    " << velocityErrors.abs.l2 << " (l2), "
    //              << velocityErrors.rel.l2 << ", " << errV.l2 << std::endl
    //              << "    " << velocityErrors.abs.linf << " (linf), "
    //              << velocityErrors.rel.linf << ", " << errV.linf << std::endl;

    //    std::cout << "Errors (inf-norm): " << std::endl
    //              << std::setw(8) << "Velocity " << std::setw(14) << "(abs)" << std::setw(14) << "(rel)" << //std::setw(20) << "(norm of dataset)" << std::endl
    //              << std::setw(8) << " " << std::setw(14) << velocityErrors.abs.linf <<  std::setw(14) << //velocityErrors.rel.linf << std::setw(20) << errV.linf << std::endl;

    return std::move(velocityErrors);
}

std::tuple<Errors, Errors> computeErrors(const CellData &dataA,
                                         const CellData &dataB)
{
    checkCellCenters(dataA, dataB);
    std::cout.precision(6);
    std::cout << std::scientific;
    const Errors pressureErrors = computePressureError(dataA, dataB);
    const Errors velocityErrors = computeVelocityError(dataA, dataB);

    std::cout << "Errors (inf-norm): " << std::endl
              << std::setw(10) << " " << std::setw(14) << "(abs)"
              << std::setw(14) << "(rel)" << std::setw(20)
              << "(norm of dataset)" << std::endl
              << std::setw(10) << "Pressure " << std::setw(14)
              << pressureErrors.abs.linf << std::setw(14)
              << pressureErrors.rel.linf << std::setw(20)
              << pressureErrors.dataset_norm.linf << std::endl;

    std::cout << std::setw(10) << "Velocity " << std::setw(14)
              << velocityErrors.abs.linf << std::setw(14)
              << velocityErrors.rel.linf << std::setw(20)
              << velocityErrors.dataset_norm.linf << std::endl;

    return std::make_tuple(pressureErrors, velocityErrors);
}

int main(int argc, char *argv[])
{
    // parse command line arguments
    if (argc != 3) {
        std::cout << "usage: " << argv[0]
                  << " <1. VTK filename> <2. VTK filename>" << std::endl;
        exit(-1);
    }

    const std::string filename1(argv[1]);
    const std::string filename2(argv[2]);

    CellData parsedData1;
    CellData parsedData2;

    parseData(filename1, parsedData1);
    parseData(filename2, parsedData2);

    Errors pressureErrors;
    Errors velocityErrors;
    std::tie(pressureErrors, velocityErrors) =
        computeErrors(parsedData1, parsedData2);

    return EXIT_SUCCESS;
}