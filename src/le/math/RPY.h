//
// Created by sail on 19-9-7.
//

#ifndef LE_RPY_H
#define LE_RPY_H

/**
 * @file RPY
 */

#include "Rotation3DVector.h"
#include "Rotation3D.h"
#include "Vector3D.h"

namespace le { namespace math {

        /** @addtogroup math */
        /*@{*/

        /**
         * @brief A class for representing Roll-Pitch-Yaw Euler angle rotations.
         */
        template<class T = double>
        class RPY: public Rotation3DVector<T> {
        public:

            /**
             * @brief Constructs rotation in which all elements are initialized to 0
             */
            RPY() : _rpy(0.0, 0.0, 0.0) {}

            /**
             * @brief Constructs an initialized roll-pitch-yaw euler angle vector
             * @param roll Rotation around z
             * @param pitch Rotation around y
             * @param yaw Rotation around x
             */
            RPY(T roll, T pitch, T yaw) : _rpy(roll, pitch, yaw) {}

            /**
             * @brief Constructs an RPY object initialized according to the specified Rotation3D
             *
             * \f$ \beta = arctan2(-r_{31},\sqrt{r_{11}^2+r_{21}^2}) \f$
             * \f$ \alpha = arctan2(r_{21}/cos(\beta), r_{11}/cos(\beta)) \f$
             * \f$ \beta = arctan2(r_{32}/cos(\beta), r_{33}/cos(\beta))) \f$
             *
             * @param R [in] A 3x3 rotation matrix \f$ \mathbf{R} \f$
             *
             * @param epsilon [in] Value specifying the value for which \f$
             * cos(\beta)\f$ is assumed 0 and the special case solution assuming
             * \f$\alpha=0, \beta=\pi/2 and \gamma = arctan2(r_{21}, r_{22})\f$ is
             * to be used.
             */
            explicit RPY(const Rotation3D<T>& R, T epsilon = 1e-5)
            {
                const T r11 = R(0, 0);
                const T r12 = R(0, 1);
                const T r21 = R(1, 0);
                const T r22 = R(1, 1);
                const T r31 = R(2, 0);
                const T r32 = R(2, 1);
                const T r33 = R(2, 2);

                double sum = r11*r11 + r21*r21;

                // We check for rounding errors in input so that we don't take the sqrt of
                // some negative number.
                if (sum < 0) sum = 0;
                // TODO: Check that if sum < 0 then sum ~= 0 also.

                const double cos_beta = sqrt(sum);
                const double sin_beta = -r31;

                // If beta == 90 deg or beta == -90 deg:
                if (fabs(cos_beta) < epsilon) {

                    // If beta == -90 deg:
                    if (sin_beta < 0) {
                        _rpy(0) = 0;
                        _rpy(1) = static_cast<T>(-M_PI / 2);
                        _rpy(2) = - atan2(r12, r22);
                    }

                        // If beta == 90 deg:
                    else {
                        _rpy(0) = 0;
                        _rpy(1) = static_cast<T>(M_PI / 2);
                        _rpy(2) = atan2(r12, r22);
                    }

                } else {
                    _rpy(1) = static_cast<T>(atan2(sin_beta, cos_beta));
                    //_rpy(0) = static_cast<T>(atan2(r21 / cos_beta, r11 / cos_beta));
                    _rpy(0) = static_cast<T>(atan2(r21, r11));
                    //_rpy(2) = static_cast<T>(atan2(r32 / cos_beta, r33 / cos_beta));
                    _rpy(2) = static_cast<T>(atan2(r32, r33));
                }
            }

            /**
             * @copydoc Rotation3DVector::toRotation3D
             */
            const Rotation3D<T> toRotation3D() const
            {
                const T a = _rpy(0);
                const T b = _rpy(1);
                const T c = _rpy(2);

                const T ca = cos(a);
                const T sa = sin(a);
                const T cb = cos(b);
                const T sb = sin(b);
                const T cc = cos(c);
                const T sc = sin(c);

                return Rotation3D<T>(
                        ca * cb,
                        ca * sb * sc-sa * cc,
                        ca * sb * cc+sa * sc,

                        sa * cb,
                        sa * sb * sc+ca * cc,
                        sa * sb * cc-ca * sc,

                        -sb,
                        cb * sc,
                        cb * cc);
            }

            /**
             * @brief Returns reference to the element
             *
             * @param index [in] index of element
             *
             * @return reference to the element
             */
            T& operator()(size_t index){
                return _rpy(index);
            }

            /**
             * @brief Returns a const reference the an element
             *
             * @param index [in] index of element
             *
             * @return const reference to the element
             */
            const T& operator()(size_t index) const {
                return _rpy(index);
            }


            /**
             * @brief Returns a const reference the an element
             *
             * @param i [in] index of element
             *
             * @return const reference to the element
             */
            const T& operator[](size_t i) const { return (*this)(i); }

            /**
              * @brief Returns reference to the element
              *
              * @param i [in] index of element
              *
              * @return reference to the element
              */
            T& operator[](size_t i) { return (*this)(i); }

            /**
             * @brief Comparison operator.
             *
             * The comparison operator makes a element wise comparison.
             * Returns true only if all elements are equal.
             *
             * @param rhs [in] RPY to compare with
             * @return True if equal.
             */
            bool operator==(const RPY<T> &rhs) const {
                return (_rpy(0) == rhs(0) && _rpy(1) == rhs(1) && _rpy(2) == rhs(2));
            }

            /**
             * @brief Comparison operator.
             *
             * The comparison operator makes a element wise comparison.
             * Returns true if any of the elements are different.
             *
             * @param rhs [in] RPY to compare with
             * @return True if not equal.
             */
            bool operator!=(const RPY<T> &rhs) const {
                return !(*this == rhs);
            }

            /**
             * @brief size of this RPY.
             * @return the value 3
             */
            size_t size() const { return 3; }

            /**
             * @brief Ouputs RPY to stream
             *
             * @param os [in/out] stream to use
             *
             * @param rpy [in] rpy rotation
             *
             * @return the resulting stream
             */
            friend std::ostream& operator<<(std::ostream& os, const RPY<T>& rpy){
                return os <<"RPY {"<<rpy(0)<<", "<<rpy(1)<<", "<<rpy(2)<<"}";
                //return os << rpy._rpy;
            }

        private:
            Vector3D<T> _rpy;
        };

        /**
        * @brief Casts RPY<T> to RPY<Q>
        *
        * @param rpy [in] RPY with type T
        *
        * @return RPY with type Q
        */
        template<class Q, class T>
        const RPY<Q> cast(const RPY<T>& rpy) {
            return RPY<Q>(
                    static_cast<Q>(rpy(0)),
                    static_cast<Q>(rpy(1)),
                    static_cast<Q>(rpy(2)));
        }

        template class le::math::RPY<double>;
        template class le::math::RPY<float>;

        /*@}*/
    }} // end namespaces
    
#endif //LE_RPY_H
