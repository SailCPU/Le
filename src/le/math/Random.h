//
// Created by Sail Yang on 2019-09-14.
//

#ifndef PROJECT_RANDOM_H
#define PROJECT_RANDOM_H

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>

namespace
{
    boost::mt19937 generator;
    boost::uniform_real<> distributor;

    typedef boost::normal_distribution<double> dist_type;

    boost::variate_generator<boost::mt19937&, dist_type> normal_distribution(
            generator,
            boost::normal_distribution<double>(0, 1));
}

namespace le {
    namespace math {
//! @addtogroup math

//! @{
/**
 * @brief Generation of random numbers.
 */
        class Random {
        public:
            Random() = delete;
            ~Random() = delete;

            /**
             * @brief A random double in the range [0, 1[ using a uniform distribution.
             *
             * @note Uses boost::random
             */
            static double ran()
            {
                return distributor(generator);
            }

            /**
             * @brief Seeds the random number generator.
             *
             * @note Uses boost::random
             */
            static void seed(unsigned seed)
            {
                generator.seed(
                        static_cast<boost::mt19937::result_type>(
                                seed));
            }

            /**
             * @brief Seeds the random number generator with current time of day
             *
             * @note Uses boost::random
             */
            static void seed()
            {
                seed((unsigned)rw::common::TimerUtil::currentTimeMs());
            }

            /**
             * @brief A random double in the range [from, to[ using a uniform distribution.
             *
             * @note Uses boost::random
             */
            static double ran(double from, double to)
            {
                if(from>to) {
                    RW_THROW("From must be smaller than to: " << from << ">" << to);
                } else if(from==to){
                    return from;
                }

                double res = from;
                do {
                    res = from + (to - from) * Random::ran();
                } while (res >= to);

                return res;
            }

            /**
             * @brief A random integer in the range [from, to[ using a uniform distribution.
             *
             * @note Uses boost::random
             */
            static int ranI(int from, int to)
            {
                return (int)floor(Random::ran(from, to));
            }

            /**
             * @brief Returns a random sample around \b mean with standard deviation \b sigma  using the normal distribution.
             *
             * @note Uses boost::random
             *
             * @param mean [in] Means value
             * @param sigma [in] Standard deviation
             * @return Random sample
             */
            static double ranNormalDist(double mean, double sigma)
            {
                return mean + sigma * normal_distribution();
            }
        };
//! @}
    } /* namespace math */
} /* namespace rw */

#endif //PROJECT_RANDOM_H
