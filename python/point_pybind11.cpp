#include "pybind11/pybind11.h"
#include "cartesian_kernel.h"
#include "solve_lp.h"
#include "pybind11/numpy.h"
#include <eigen3/Eigen/Eigen>
#include "polytopes.h"
#include <iostream>

typedef float                    NT;
typedef Cartesian<NT>             Kernel;
typedef typename Kernel::Point    Point;
typedef HPolytope<Point> Polytope;

namespace py=pybind11;

PYBIND11_MODULE(volesti, m) {
    py::class_<Point>(m, "Point", py::buffer_protocol())
        .def(py::init([](py::array_t<NT> b) {
            py::buffer_info info = b.request();
            Point* p = new Point(info.shape[0]);
            p->set_dimension(info.shape[0]);
            NT* buffer_data = static_cast<NT*>(info.ptr);
            for (int i=0; i<p->dimension(); i++) {
                p->set_coord(i, buffer_data[i]);
            }
            return p;
        }))
        .def_buffer([](Point &p) -> py::buffer_info {
            return py::buffer_info(
                p.data(),
                sizeof(NT),
                py::format_descriptor<NT>::format(),
                1,
                { p.dimension() },
                { sizeof(NT) }
            );
        })
        .def("__add__", [](Point &p, Point &p2) {
            return p+p2;
        }, py::is_operator());

    py::class_<Polytope>(m, "Polytope")
        .def(py::init([](py::array_t<NT> A, py::array_t<NT> b) {
            py::buffer_info A_info = A.request();
            py::buffer_info b_info = b.request();

            Eigen::Map<Polytope::MT> A_map((NT*) A_info.ptr, A_info.shape[0], A_info.shape[1]);
            Eigen::Map<Polytope::VT> b_map((NT*) b_info.ptr, b_info.shape[0], 1);
            Polytope* p = new Polytope(A_map, b_map);
            return p;
        }))
        .def("print", &Polytope::print)
        .def("create_point_representation", [](Polytope &p, py::array_t<NT> internalPoint) {
            py::buffer_info ip_info = internalPoint.request();
            p.create_point_representation((NT*) ip_info.ptr);
        })
        .def("is_in", [](Polytope &p, py::array_t<NT> point) -> bool {
            py::buffer_info point_info = point.request();
            return p.contains_point((NT* ) point_info.ptr);
        });

}
