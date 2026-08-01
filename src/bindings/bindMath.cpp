#include "bindings.h"

#include <pybind11/operators.h>

#include "../engine/math/Vector3.h"


void bindMath(py::module_& m)
{
    py::class_<Vector3>(m, "Vector3")
        // VECTOR DEFINITION
        .def(py::init<>())
        .def(py::init<float, float, float>(), py::arg("x"), py::arg("y"), py::arg("z"))
        .def_readwrite("x", &Vector3::x)
        .def_readwrite("y", &Vector3::y)
        .def_readwrite("z", &Vector3::z)

        // VECTOR NORMS
        .def("squaredNorm", &Vector3::squaredNorm)
        .def("norm", &Vector3::norm)
        .def("normalize", &Vector3::normalize)
        .def("normalized", &Vector3::normalized)

        // VECTOR OPERATIONS
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(-py::self)
        .def(py::self * float())
        .def(py::self / float())
        .def(py::self *= float())
        .def(float() * py::self)


        // REPRESENTATION
        .def("__repr__", [](const Vector3& vector)
            {
                return "Vector3(" +
                    std::to_string(vector.x) + ", " +
                    std::to_string(vector.y) + ", " +
                    std::to_string(vector.z) + ")";
            }
        );


    // FREE FUNCTIONS
    m.def("dot", &dot, py::arg("a"), py::arg("b"));
    m.def("cross", &cross, py::arg("a"), py::arg("b"));
}