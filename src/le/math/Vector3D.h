//
// Created by sail on 19-9-7.
//

#ifndef LE_VECTOR3D_H
#define LE_VECTOR3D_H

#include "Eigen/Eigen"

namespace le {
    namespace math {

        /** @addtogroup math */
        /*@{*/

        /**
         * @brief A 3D vector @f$ \mathbf{v}\in \mathbb{R}^3 @f$
         *
         * @f$ \robabx{i}{j}{\mathbf{v}} = \left[
         *  \begin{array}{c}
         *  v_x \\
         *  v_y \\
         *  v_z
         *  \end{array}
         *  \right]
         *  @f$
         *
         *  Usage example:
         *
         *  \code
         *  const Vector3D<> v1(1.0, 2.0, 3.0);
         *  const Vector3D<> v2(6.0, 7.0, 8.0);
         *  const Vector3D<> v3 = cross(v1, v2);
         *  const double d = dot(v1, v2);
         *  const Vector3D<> v4 = v2 - v1;
         *  \endcode
         */
        template<class T = double>
        class Vector3D {
        public:

            //! Eigen type equivalent to Vector3D
            typedef Eigen::Matrix<T, 3, 1> EigenVector3D;


            //! Value type.
            typedef T value_type;

            /**
             * @brief Creates a 3D vector initialized with 0's
             */
            Vector3D() {
                _vec[0] = 0;
                _vec[1] = 0;
                _vec[2] = 0;
            }

            /**
             * @brief Creates a 3D vector
             * @param x [in] @f$ x @f$
             * @param y [in] @f$ y @f$
             * @param z [in] @f$ z @f$
             */
            Vector3D(T x, T y, T z) {
                _vec[0] = x;
                _vec[1] = y;
                _vec[2] = z;
            }

            /**
             * @brief Creates a 3D vector from vector_expression
             *
             * @param r [in] an ublas vector_expression
             */
            template<class R>
            explicit Vector3D(const Eigen::MatrixBase <R> &r) {
                _vec[0] = r.row(0)(0);
                _vec[1] = r.row(1)(0);
                _vec[2] = r.row(2)(0);
            }

            /**
             * @brief Returns Eigen vector with the values of *this
             */
            const EigenVector3D e() const {
                EigenVector3D v;
                v(0) = _vec[0];
                v(1) = _vec[1];
                v(2) = _vec[2];
                return v;
            }

            /**
               @brief The dimension of the vector (i.e. 3).

               This method is provided to help support generic algorithms using
               size() and operator[].
            */
            size_t size() const { return 3; }

            //----------------------------------------------------------------------
            // Various operators

            /**
             * @brief Returns reference to vector element
             * @param i [in] index in the vector \f$i\in \{0,1,2\} \f$
             * @return const reference to element
             */
            const T &operator()(size_t i) const { return _vec[i]; }

            /**
             * @brief Returns reference to vector element
             * @param i [in] index in the vector \f$i\in \{0,1,2\} \f$
             * @return reference to element
             */
            T &operator()(size_t i) { return _vec[i]; }

            /**
             * @brief Returns reference to vector element
             * @param i [in] index in the vector \f$i\in \{0,1,2\} \f$
             * @return const reference to element
             */
            const T &operator[](size_t i) const { return _vec[i]; }

            /**
             * @brief Returns reference to vector element
             * @param i [in] index in the vector \f$i\in \{0,1,2\} \f$
             * @return reference to element
             */
            T &operator[](size_t i) { return _vec[i]; }


            /**
               @brief Scalar division.
             */
            const Vector3D<T> operator/(T s) const {
                return Vector3D<T>(_vec[0] / s, _vec[1] / s, _vec[2] / s);
            }

            /**
               @brief Scalar multiplication.
             */
            const Vector3D<T> operator*(T s) const {
                return Vector3D<T>(_vec[0] * s, _vec[1] * s, _vec[2] * s);
            }

            /**
               @brief Scalar multiplication.
             */
            friend const Vector3D<T> operator*(T s, const Vector3D<T> &v) {
                return Vector3D<T>(v[0] * s, v[1] * s, v[2] * s);
            }

            /**
               @brief Vector subtraction.
             */
            const Vector3D<T> operator-(const Vector3D<T> &b) const {
                return Vector3D<T>(_vec[0] - b[0], _vec[1] - b[1], _vec[2] - b[2]);
            }

            /**
               @brief Vector addition.
             */
            const Vector3D<T> operator+(const Vector3D<T> &b) const {
                return Vector3D<T>(_vec[0] + b[0], _vec[1] + b[1], _vec[2] + b[2]);
            }

            /**
               @brief Scalar multiplication.
             */
            Vector3D<T> &operator*=(T s) {
                _vec[0] *= s;
                _vec[1] *= s;
                _vec[2] *= s;
                return *this;
            }

            /**
               @brief Scalar division.
             */
            Vector3D<T> &operator/=(T s) {
                _vec[0] /= s;
                _vec[1] /= s;
                _vec[2] /= s;
                return *this;
            }

            /**
               @brief Vector addition.
             */
            Vector3D<T> &operator+=(const Vector3D<T> &v) {
                _vec[0] += v._vec[0];
                _vec[1] += v._vec[1];
                _vec[2] += v._vec[2];
                return *this;
            }

            /**
               @brief Vector subtraction.
             */
            Vector3D<T> &operator-=(const Vector3D<T> &v) {
                _vec[0] -= v._vec[0];
                _vec[1] -= v._vec[1];
                _vec[2] -= v._vec[2];
                return *this;
            }

            /**
               @brief Unary minus.
             */
            const Vector3D<T> operator-() const {
                return Vector3D<T>(-_vec[0], -_vec[1], -_vec[2]);
            }

            /**
               @brief Streaming operator.
             */
            friend std::ostream &operator<<(std::ostream &out, const Vector3D<T> &v) {
                return out
                        << "Vector3D("
                        << v[0] << ", " << v[1] << ", " << v[2]
                        << ")";
            }

            //----------------------------------------------------------------------
            // Various friend functions

            /**
             * @brief Returns the Euclidean norm (2-norm) of the vector
             * @return the norm
             */
            T norm2() const {
                return sqrt(_vec[0] * _vec[0] + _vec[1] * _vec[1] + _vec[2] * _vec[2]);
            }

            /**
             * @brief Returns the Manhatten norm (1-norm) of the vector
             * @return the norm
             */
            T norm1() const {
                return fabs(_vec[0]) + fabs(_vec[1]) + fabs(_vec[2]);
                //return norm_1(_vec);
            }

            /**
             * @brief Returns the infinte norm (\f$\inf\f$-norm) of the vector
             * @return the norm
             */
            T normInf() const {
                T res = fabs(_vec[0]);
                const T f1 = fabs(_vec[1]);
                if (f1 > res)
                    res = f1;
                const T f2 = fabs(_vec[2]);
                if (f2 > res)
                    res = f2;
                return res;
                //return Math::max(fabs(_vec[0]), Math::max(fabs(_vec[1]), fabs(_vec[2])));
                //return norm_inf(_vec);
            }

            /**
               @brief Compare with \b b for equality.
               @param b [in] other vector.
               @return True if a equals b, false otherwise.
            */
            bool operator==(const Vector3D<T> &b) const {
                return _vec[0] == b[0] && _vec[1] == b[1] && _vec[2] == b[2];
            }

            /**
               @brief Compare with \b b for inequality.
               @param b [in] other vector.
               @return True if a and b are different, false otherwise.
            */
            bool operator!=(const Vector3D<T> &b) const {
                return !(*this == b);
            }

            /**
             * @brief Get zero vector.
             * @return vector.
             */
            static Vector3D<T> zero() { return Vector3D<T>(0, 0, 0); }

            /**
             * @brief Get x vector (1,0,0)
             * @return vector.
             */
            static Vector3D<T> x() { return Vector3D<T>(1.0, 0, 0); }

            /**
             * @brief Get y vector (0,1,0)
             * @return vector.
             */
            static Vector3D<T> y() { return Vector3D<T>(0, 1.0, 0); }

            /**
             * @brief Get z vector (0,0,1)
             * @return vector.
             */
            static Vector3D<T> z() { return Vector3D<T>(0, 0, 1.0); }

        private:
            T _vec[3];
        };

        /**
         * @brief Calculates the 3D vector cross product @f$ \mathbf{v1} \times \mathbf{v2} @f$
         * @param v1 [in] @f$ \mathbf{v1} @f$
         * @param v2 [in] @f$ \mathbf{v2} @f$
         *
         * @return the 3D vector cross product @f$ \mathbf{v1} \times \mathbf{v2} @f$
         *
         * The 3D vector cross product is defined as:
         * @f$
         * \mathbf{v1} \times \mathbf{v2} = \left[\begin{array}{c}
         *  v1_y * v2_z - v1_z * v2_y \\
         *  v1_z * v2_x - v1_x * v2_z \\
         *  v1_x * v2_y - v1_y * v2_x
         * \end{array}\right]
         * @f$
         *
         * @relates Vector3D
         */
        template<class T>
        const Vector3D<T> cross(const Vector3D<T> &v1, const Vector3D<T> &v2) {
            return Vector3D<T>(
                    v1[1] * v2[2] - v1[2] * v2[1],
                    v1[2] * v2[0] - v1[0] * v2[2],
                    v1[0] * v2[1] - v1[1] * v2[0]);
        }

        /**
         * @brief Calculates the 3D vector cross product @f$ \mathbf{v1} \times \mathbf{v2} @f$
         * @param v1 [in] @f$ \mathbf{v1} @f$
         * @param v2 [in] @f$ \mathbf{v2} @f$
         * @param dst [out] the 3D vector cross product @f$ \mathbf{v1} \times \mathbf{v2} @f$
         *
         * The 3D vector cross product is defined as:
         * @f$
         * \mathbf{v1} \times \mathbf{v2} = \left[\begin{array}{c}
         *  v1_y * v2_z - v1_z * v2_y \\
         *  v1_z * v2_x - v1_x * v2_z \\
         *  v1_x * v2_y - v1_y * v2_x
         * \end{array}\right]
         * @f$
         *
         * @relates Vector3D
         */
        template<class T>
        void cross(const Vector3D<T> &v1, const Vector3D<T> &v2, Vector3D<T> &dst) {
            dst[0] = v1[1] * v2[2] - v1[2] * v2[1];
            dst[1] = v1[2] * v2[0] - v1[0] * v2[2];
            dst[2] = v1[0] * v2[1] - v1[1] * v2[0];
        }


        /**
         * @brief Calculates the dot product @f$ \mathbf{v1} . \mathbf{v2} @f$
         * @param v1 [in] @f$ \mathbf{v1} @f$
         * @param v2 [in] @f$ \mathbf{v2} @f$
         *
         * @return the dot product @f$ \mathbf{v1} . \mathbf{v2} @f$
         *
         * @relates Vector3D
         */
        template<class T>
        T dot(const Vector3D<T> &v1, const Vector3D<T> &v2) {
            return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
            //return inner_prod(v1.m(), v2.m());
        }

        /**
         * @brief Returns the normalized vector \f$\mathbf{n}=\frac{\mathbf{v}}{\|\mathbf{v}\|} \f$.
         * In case \f$ \|mathbf{v}\| = 0\f$ the zero vector is returned.
         * @param v [in] \f$ \mathbf{v} \f$ which should be normalized
         * @return the normalized vector \f$ \mathbf{n} \f$
         *
         * @relates Vector3D
         */
        template<class T>
        const Vector3D<T> normalize(const Vector3D<T> &v) {
            T length = v.norm2();
            if (length != 0)
                return Vector3D<T>(v(0) / length, v(1) / length, v(2) / length);
            else
                return Vector3D<T>(0, 0, 0);
        }

        /**
         * @brief Calculates the angle from @f$ \mathbf{v1}@f$ to @f$ \mathbf{v2} @f$
         * around the axis defined by @f$ \mathbf{v1} \times \mathbf{v2} @f$ with n
         * determining the sign.
         * @param v1 [in] @f$ \mathbf{v1} @f$
         * @param v2 [in] @f$ \mathbf{v2} @f$
         * @param n [in] @f$ \mathbf{n} @f$
         *
         * @return the angle
         *
         * @relates Vector3D
         */
        template<class T>
        double angle(const Vector3D<T> &v1, const Vector3D<T> &v2, const Vector3D<T> &n) {
            const Vector3D<T> nv1 = normalize(v1);
            const Vector3D<T> nv2 = normalize(v2);
            const Vector3D<T> nn = normalize(n);
            return atan2(dot(nn, cross(nv1, nv2)), dot(nv1, nv2));
        }

        /**
         * @brief Calculates the angle from @f$ \mathbf{v1}@f$ to @f$ \mathbf{v2} @f$
         * around the axis defined by @f$ \mathbf{v1} \times \mathbf{v2} @f$
         * @param v1 [in] @f$ \mathbf{v1} @f$
         * @param v2 [in] @f$ \mathbf{v2} @f$
         *
         * @return the angle
         *
         * @relates Vector3D
         */
        template<class T>
        double angle(const Vector3D<T> &v1, const Vector3D<T> &v2) {
            Vector3D<T> n = cross(v1, v2);
            return angle(v1, v2, n);
        }

        /**
         * @brief Casts Vector3D<T> to Vector3D<Q>
         * @param v [in] Vector3D with type T
         * @return Vector3D with type Q
         *
         * @relates Vector3D
         */
        template<class Q, class T>
        const Vector3D<Q> cast(const Vector3D<T> &v) {
            return Vector3D<Q>(
                    static_cast<Q>(v(0)),
                    static_cast<Q>(v(1)),
                    static_cast<Q>(v(2)));
        }

        template class le::math::Vector3D<double>;

        template class le::math::Vector3D<float>;

        /**@}*/
    }
} // end namespaces

#endif //LE_VECTOR3D_H
