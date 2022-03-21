#include <vector>
#include <virtual-casing.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;


PYBIND11_MODULE(virtual_casing, m) {

    py::class_<sctl::Vector<double>>(m, "SCTLDoubleVector")
         .def(py::init<>())
         .def(py::init<const std::vector<double> & >());

    py::class_<sctl::Vector<int>>(m, "SCTLIntVector")
        .def(py::init<>())
        .def(py::init<const std::vector<int> & >());

    py::enum_<biest::SurfType>(m, "SurfType")
        .value("AxisymCircleWide", biest::SurfType::AxisymCircleWide)           // 125 x 50
        .value("AxisymCircleNarrow", biest::SurfType::AxisymCircleNarrow)       // 250 x 50
        .value("AxisymWide", biest::SurfType::AxisymWide)                       // 125 x 50
        .value("AxisymNarrow", biest::SurfType::AxisymNarrow)                   // 250 x 50
        .value("RotatingEllipseWide", biest::SurfType::RotatingEllipseWide)     // 125 x 50
        .value("RotatingEllipseNarrow", biest::SurfType::RotatingEllipseNarrow) // 250 x 50
        .value("Quas3", biest::SurfType::Quas3)                                 // 250 x 50
        .value("LHD", biest::SurfType::LHD)                                     // 250 x 50
        .value("W7X", biest::SurfType::W7X)                                     // 250 x 45
        .value("Stell", biest::SurfType::Stell)                                 // 250 x 50
        .value("W7X_", biest::SurfType::W7X_);

    py::class_<VirtualCasing<double>>(m, "VirtualCasing")
        .def(py::init<>())
        .def("set_surface", &VirtualCasing<double>::SetSurface)
        .def("set_accuracy", &VirtualCasing<double>::SetAccuracy)
        .def("compute_external_B", &VirtualCasing<double>::ComputeBext);

    py::class_<VirtualCasingTestData<double>>(m, "VirtualCasingTestData")
        .def(py::init<>())
        .def_static("surface_coordinates", &VirtualCasingTestData<double>::SurfaceCoordinates)
        .def_static("magnetic_field_data", &VirtualCasingTestData<double>::BFieldData);

}

