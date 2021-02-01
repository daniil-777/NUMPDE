/**
 * @file radauthreetimesteppingode.cc
 * @brief NPDE homework RadauThreeTimestepping
 * @author Daniil Emtsev
 * @date 01/02/2021
 * @copyright Developed at ETH Zurich
 */

#include "radauthreetimesteppingode.h"

#include <cmath>
#include <iostream>
#include <vector>

namespace RadauThreeTimestepping {

/* SAM_LISTING_BEGIN_1 */
std::vector<double> twoStageRadauTimesteppingLinScalODE(unsigned int m) {
  std::vector<double> sol_vec;
  double tau = 5.0 / m;
  double y_current = 1.0;
  double mult = 1 - tau * (1 + tau / 6) / ((1 + 5 * tau / 12) * (1 + tau / 4) + tau * tau / 16);
  sol_vec.push_back(y_current);
  for (int i = 0; i < m; i++){
    y_current *= mult;
    sol_vec.push_back(y_current);
  }
  //====================
  // Your code goes here
  //====================
  return sol_vec;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
void testConvergenceTwoStageRadauLinScalODE() {
  constexpr int nIter = 10;       // total number of iterations
  double max_norm_errors[nIter];  // errors vector for all approx. sols
  double rates[nIter - 1];        // The rates of convergence
  double avg_rate = 0.0;  // The average rate of convergence over all iterations

  //====================
  // Your code goes here
  //====================
  /* SAM_LISTING_END_2 */

  // Printing results
  std::cout << "\n" << std::endl;
  std::cout << "*********************************************************"
            << std::endl;
  std::cout << "         Convergence of two-stage Radau Method           "
            << std::endl;
  std::cout << "*********************************************************"
            << std::endl;
  std::cout << "--------------------- RESULTS ---------------------------"
            << std::endl;
  std::cout << "Iteration"
            << "\t| Nsteps"
            << "\t| error"
            << "\t\t| rates" << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;
  for (int k = 0; k < nIter; k++) {
    std::cout << k << "\t"
              << "\t|" << 10 * std::pow(2, k) << "\t\t|" << max_norm_errors[k];
    if (k > 0) {
      std::cout << "\t|" << rates[k - 1];
    }
    std::cout << "\n";
  }
  std::cout << "---------------------------------------------------------"
            << std::endl;
  std::cout << "Average rate of convergence: " << avg_rate << "\n" << std::endl;
}

}  // namespace RadauThreeTimestepping
