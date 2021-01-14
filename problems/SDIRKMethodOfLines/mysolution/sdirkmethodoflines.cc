/** @file
 * @brief NPDE SDIRKMethodOfLines
 * @author Erick Schulz
 * @date 12/04/2019
 * @copyright Developed at ETH Zurich
 */

#include "sdirkmethodoflines.h"

namespace SDIRKMethodOfLines {

/** @brief This class implements a Lehrfem++ matrix provider defining a
LinFEMassMatrixProvider::Eval function returning the local MASS matrix for
linear first-order lagrange FE bases over triangular mesh (only!). Integration
over the triangular cells is performed using the trapezoidal rule.*/
class LinFEMassMatrixProvider {
 public:
  /** @brief default constructor */
  explicit LinFEMassMatrixProvider() = default;
  /** @brief Default implement: all cells are active */
  virtual bool isActive(const lf::mesh::Entity & /*cell*/) { return true; }
  /** @brief Main method for computing the element vector
   * @param cell refers to current cell for which the element vector is desired
   * The implementation uses an analytic formula defined over triangular cells*/
  Eigen::Matrix<double, 3, 3> Eval(const lf::mesh::Entity &tria);
};  // class LinFEMassMatrixProvider
/** Implementing member function Eval of class LinFEMassMatrixProvider*/
/* SAM_LISTING_BEGIN_8 */
Eigen::Matrix<double, 3, 3> LinFEMassMatrixProvider::Eval(
    const lf::mesh::Entity &tria) {
  Eigen::Matrix<double, 3, 3> elMat;
  LF_VERIFY_MSG(tria.RefEl() == lf::base::RefEl::kTria(), "unsupp type" << tria);
  const double area = lf::geometry::Volume(*(tria.Geometry()));
  //auto endpoints = lf::geometry::Corners(*(tria.Geometry()));
  //double edge_length = (endpoints.col(1) - endpoints.col(0)).norm();

  elMat << 2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0;
  elMat *= area / 12.0;
  //====================
  // Your code goes here
  //====================

  return elMat;  // return the local mass element matrix
}  // LinFEMassMatrixProvider::Eval
/* SAM_LISTING_END_8 */

/** @brief This class implements a Lehrfem++ matrix provider defining a
LinearMassEdgeMatrixProvider<FUNCTOR>::Eval function returning the local EDGE
MASS matrix for linear first-order lagrange FE over triangular mesh (only!).
Integration over the triangular cells is performed using the trapezoidal rule.*/
template <typename FUNCTOR>
class LinearMassEdgeMatrixProvider {
 public:
  /** @brief Constructor storing the right hand side function
   *  @param predicate is lambda function returning booleans on edge entities
   *         cool_coeff is the convective cooling coefficient */
  explicit LinearMassEdgeMatrixProvider(FUNCTOR predicate, double cool_coeff)
      : predicate_(predicate), cool_coeff_(cool_coeff) {}
  /** @brief Default implement: all edges are active */
  virtual bool isActive(const lf::mesh::Entity & /*edge*/) { return true; }
  /** @brief Main method for computing the element vector
   * @param edge is current entity for which the element vector is desired
   * The implementation uses simple vertex based quadrature */
  Eigen::Matrix<double, 2, 2> Eval(const lf::mesh::Entity &edge);

 private:
  /** predicate_ provides booleans for boundary edges */
  FUNCTOR predicate_;
  double cool_coeff_;
};  // class LinearMassEdgeMatrixProvider
/* Implementing member function Eval of class LinearMassEdgeMatrixProvider */
/* SAM_LISTING_BEGIN_9 */
template <typename FUNCTOR>
Eigen::Matrix<double, 2, 2> LinearMassEdgeMatrixProvider<FUNCTOR>::Eval(
    const lf::mesh::Entity &edge) {
  Eigen::Matrix<double, 2, 2> elBdyEdgeMat;
  if(predicate_(edge)){
  auto endpoints = lf::geometry::Corners(*(edge.Geometry()));
  double edge_length = (endpoints.col(1) - endpoints.col(0)).norm();
  elBdyEdgeMat << 2, 1, 1, 2;
  double coef = edge_length * cool_coeff_ / 6.0;
  elBdyEdgeMat *= coef;
  //====================
  // Your code goes here
  //====================
  return elBdyEdgeMat;
}else{
  return Eigen::Matrix<double, 2, 2>::Zero();
}  // LinearMassEdgeMatrixProvider<FUNCTOR>::Eval
/* SAM_LISTING_END_9 */
    }

/* SAM_LISTING_BEGIN_1 */
std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>>
assembleGalerkinMatrices(const lf::assemble::DofHandler &dofh,
                         double cool_coeff) {
  std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>>
      sparse_pair;

  Eigen::SparseMatrix<double> A_sparse;
  Eigen::SparseMatrix<double> M_sparse;
  auto mesh_p = dofh.Mesh();
  auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  double Numdofs = dofh.NumDofs();
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)}; //edges
  auto edges_predicate = [&bd_flags](const lf::mesh::Entity &edge) -> bool{return bd_flags(edge);};

  lf::assemble::COOMatrix<double> A(Numdofs, Numdofs);
  lf::assemble::COOMatrix<double> M(Numdofs, Numdofs);
  // auto alpha = [](Eigen::Vector2d x) -> double{return 1.0};
  // auto betta = [](Eigen::Vector2d x) -> double{return 0.0};
  // auto gamma = [](Eigen::Vector2d x) -> double{return 1.0};

  auto one = lf::mesh::utils::MeshFunctionConstant(1.0);
  auto zero = lf::mesh::utils::MeshFunctionConstant(0.0);
  //auto betta = mesh_p::utils::MeshFunctionConstant(1.0);
  
  // auto one = [](Eigen::Vector2d x) -> double{return 1.0;};
  // auto zero = [](Eigen::Vector2d x) -> double{return 0.0;};

 
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<double, decltype(one), decltype(zero)> elmat_builder_A(fe_space, one, zero); //for a
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<double, decltype(zero), decltype(one)> elmat_builder_M(fe_space, zero, one); //for m
  //lf::uscalfe::MassEdgeMatrixProvider<double, decltype(betta), decltype(edges_predicate)> edge_mat_builder(fe_space, betta, edges_predicate);
  LinearMassEdgeMatrixProvider loc_edge_mass(edges_predicate, cool_coeff);
  LinFEMassMatrixProvider mass_builder_M;
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder_A, A);
  lf::assemble::AssembleMatrixLocally(1, dofh, dofh, loc_edge_mass, A);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, mass_builder_M, M);

  A_sparse = A.makeSparse();
  M_sparse = M.makeSparse();
  sparse_pair.first = A_sparse;
  sparse_pair.second = M_sparse;

  //====================
  // Your code goes here
  //====================
  return sparse_pair;
}
/* SAM_LISTING_END_1 */

/* Implementation of class SDIRK2Timestepper */
// Implementation of SDIRK2Timestepper constructor
SDIRK2Timestepper::SDIRK2Timestepper(const lf::assemble::DofHandler &dofh,
                                     double tau /*nb. steps*/,
                                     double cool_coeff /*cooling coeff*/)
    : tau_(tau) {
  
  //====================
  // Your code goes here
  //====================
}  // SDIRK2Timestepper constructor

/* Implementation of SDIRK2Timestepper member function */
/* SAM_LISTING_BEGIN_9 */
Eigen::VectorXd SDIRK2Timestepper::discreteEvolutionOperator(
    const Eigen::VectorXd &mu) const {
  Eigen::VectorXd discrete_evolution_operator;

  //====================
  // Your code goes here
  //====================
  return discrete_evolution_operator;
}  // SDIRK2Timestepper::discreteEvolutionOperator
/* SAM_LISTING_END_9 */

/** @Brief Implementing the temperature evolution solver. The solver obtains the
discrete evolution operator from the SDIRK2Timestepper class and repeatedly
iterates its applicaiton starting from the initial condition argument
* @param cool_coeff is the convective cooling coefficient */
/* SAM_LISTING_BEGIN_6 */
std::pair<Eigen::VectorXd, Eigen::VectorXd> solveTemperatureEvolution(
    const lf::assemble::DofHandler &dofh, unsigned int m, double cool_coeff,
    Eigen::VectorXd initial_temperature_vec) {
  std::pair<Eigen::VectorXd, Eigen::VectorXd> solution_pair;

  //====================
  // Your code goes here
  //====================
  return solution_pair;
}  // solveTemperatureEvolution

/* SAM_LISTING_END_6 */

/** @brief thermalEnergy computes the integral of the discrete argument function
passed as a vector of coefficients whose index is the global index of the
degree of freedom (linear lagrange basis function) that they multiply */
/* SAM_LISTING_BEGIN_7 */
double thermalEnergy(const lf::assemble::DofHandler &dofh,
                     const Eigen::VectorXd &temperature_vec) {
  double thermal_energy = 0.0;
  auto mesh_p = dofh.Mesh();
  double themal_energy_loc;
  for(const lf::mesh::Entity *tria : mesh_p->Entities(0)){
    themal_energy_loc = 0;
    double area = lf::geometry::Volume(*(tria->Geometry()));
    for(const lf::assemble::gdof_idx_t &g_idx: dofh.GlobalDofIndices(*tria)){
      themal_energy_loc += temperature_vec[g_idx];
    }
    thermal_energy += themal_energy_loc * area / 3.0; 
    //auto endpoints = lf::geometry::Corners(*(tria.Geometry()));
    // auto dof_idx = dofh.GlobalDofIndices(*tria);
    // double num_dofs = dofh.NumLocalDofs(*tria);
    // //assert(dofh.NumLocalDofs(*tria) == 3);
    // for (int i = 0; i < num_dofs; i++){
    //   thermal_energy += area / num_dofs * temperature_vec(dof_idx[i]);
    // }
  }
  //====================
  // Your code goes here
  //====================
  return thermal_energy;
}  // thermalEnergy

/* SAM_LISTING_END_7 */

}  // namespace SDIRKMethodOfLines
