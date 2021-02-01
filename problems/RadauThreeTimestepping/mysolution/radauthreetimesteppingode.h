/** @file radauthreetimesteppingode.h
 * @brief NPDE RadauThreeTimestepping
 * @author Daniil Emtsev
 * @date 01/02/2021
 * @copyright Developed at ETH Zurich
 */

#include <vector>

namespace RadauThreeTimestepping {

/**
 * @brief Solves the ODE given in exercise 6-1.c
 * @param m The number of timesteps
 * @returns A vector of solutions at every timestep
 */
std::vector<double> twoStageRadauTimesteppingLinScalODE(unsigned int m);

/**
 * @brief Print a table of convergence rates
 */
void testConvergenceTwoStageRadauLinScalODE();

}  // namespace RadauThreeTimestepping
