//
// Created by sail on 19-9-7.
//

#ifndef LE_EAA_H
#define LE_EAA_H

namespace le { namespace math {
        /** @addtogroup math */
        /*@{*/

        namespace
        {
            template<class T> Vector3D<T> angleAxis(const Rotation3D<T>& R)
            {
                typedef Vector3D<T> V;

                const T eps = static_cast<T>(1e-6);

                T cos_theta = (T)0.5 * (R(0, 0) + R(1, 1) + R(2, 2) - 1);

                // Numerical rounding errors force us to make this check:
                if (cos_theta > 1) cos_theta = 1;
                else if (cos_theta < -1) cos_theta = -1;

                // ... because otherwise this often yields NaN:
                const T angle = acos(cos_theta);
                const V axis(R(2, 1) - R(1, 2),
                             R(0, 2) - R(2, 0),
                             R(1, 0) - R(0, 1));

                if (fabs(angle) < eps) {   // Is angle close to 0 degree
                    return 0.5 * axis;
                }

                if (fabs(angle - Pi) < eps) { // Is the angle close to 180 degree
                    V values(
                            (T)0.5 * (R(0, 0) + (T)1.0),
                            (T)0.5 * (R(1, 1) + (T)1.0),
                            (T)0.5 * (R(2, 2) + (T)1.0)
                    );
                    values[0] = (values[0] < 0) ? 0 : sqrt(values[0]);
                    values[1] = (values[1] < 0) ? 0 : sqrt(values[1]);
                    values[2] = (values[2] < 0) ? 0 : sqrt(values[2]);
                    // Find index the axis element with largest value.
                    std::size_t k = 0;
                    if (std::fabs(axis[1]) > std::fabs(axis[0]))
                        k = 1;
                    if (std::fabs(axis[2]) > std::fabs(axis[1]) && std::fabs(axis[2]) > std::fabs(axis[0]))
                        k = 2;
                    const int sign_k = (axis[k] >= 0) ? 1 : -1;
                    values *= static_cast<T>(sign_k);
                    // Determine signs
                    if (k == 0) {
                        T v1v2 = R(0,1)+R(1,0);
                        T v1v3 = R(0,2)+R(2,0);
                        if (v1v2 < 0.)
                            values[1] = -values[1];
                        if (v1v3 < 0.)
                            values[2] = -values[2];
                    } else if (k == 1) {
                        T v1v2 = R(0,1)+R(1,0);
                        T v2v3 = R(1,2)+R(2,1);
                        if (v1v2 < 0.)
                            values[0] = -values[0];
                        if (v2v3 < 0.)
                            values[2] = -values[2];
                    } else if (k == 2) {
                        T v1v3 = R(0,2)+R(2,0);
                        T v2v3 = R(1,2)+R(2,1);
                        if (v1v3 < 0.)
                            values[0] = -values[0];
                        if (v2v3 < 0.)
                            values[1] = -values[1];
                    }
                    return values*(static_cast<T>(Pi)-axis.norm2()/2);
                }

                return normalize(axis)*angle;
            }
        }

        // Forward declare cross function
        template<class T> class EAA;
        template<class T> const Vector3D<T> cross(const Vector3D<T>& v, const EAA<T>& eaa);

        /**
         * @brief A class for representing an equivalent angle-axis rotation
         *
         * This class defines an equivalent-axis-angle orientation vector also known
         * as an @f$ \thetak @f$ vector or "axis+angle" vector
         *
         * The equivalent-axis-angle vector is the product of a unit vector @f$
         * \hat{\mathbf{k}} @f$ and an angle of rotation around that axis @f$ \theta
         * @f$
         *
         * @note given two EAA vectors @f$ \theta_1\mathbf{\hat{k}}_1 @f$ and @f$
         * \theta_2\mathbf{\hat{k}}_2 @f$ it is generally not possible to subtract
         * or add these vectors, except for the special case when @f$
         * \mathbf{\hat{k}}_1 == \mathbf{\hat{k}}_2 @f$ this is why this class does
         * not have any subtraction or addition operators
         */
        template<class T = double>
        class EAA : public Rotation3DVector<T>
        {
        public:
            /**
             * @brief Extracts Equivalent axis-angle vector from Rotation matrix
             *
             * @param R [in] A 3x3 rotation matrix @f$ \mathbf{R} @f$
             *
             * @f$
             * \theta = arccos(\frac{1}{2}(Trace(\mathbf{R})-1)=arccos(\frac{r_{11}+r_{22}+r_{33}-1}{2})
             * @f$
             *
             * @f$
             * \thetak=log(\mathbf{R})=\frac{\theta}{2 sin \theta}(\mathbf{R}-\mathbf{R}^T) =
             * \frac{\theta}{2 sin \theta}
             * \left[
             * \begin{array}{c}
             * r_{32}-r_{23}\\
             * r_{13}-r_{31}\\
             * r_{21}-r_{12}
             * \end{array}
             * \right]
             * @f$
             *
             * @f$
             * \thetak=
             * \left[
             * \begin{array}{c}
             * 0\\
             * 0\\
             * 0
             * \end{array}
             * \right]
             * @f$ if @f$ \theta = 0 @f$
             *
             * @f$
             * \thetak=\pi
             * \left[
             * \begin{array}{c}
             * \sqrt{(R(0,0)+1.0)/2.0}\\
             * \sqrt{(R(1,1)+1.0)/2.0}\\
             * \sqrt{(R(2,2)+1.0)/2.0}
             * \end{array}
             * \right]
             * @f$ if @f$ \theta = \pi @f$
             *
             */
            explicit EAA(const Rotation3D<T>& R):_eaa(angleAxis(R)){}

            /**
             * @brief Constructs an EAA vector initialized to \f$\{0,0,0\}\f$
             */
            EAA():
                    _eaa(0,0,0)
            {}

            /**
             * @brief Constructs an initialized EAA vector
             * @param axis [in] \f$ \mathbf{\hat{k}} \f$
             * @param angle [in] \f$ \theta \f$
             * @pre norm_2(axis) = 1
             */
            EAA(const Vector3D<T>& axis, T angle) :
                    _eaa(axis * angle)
            {}

            /**
             * @brief Constructs an initialized EAA vector
             * @f$ \thetak =
             * \left[\begin{array}{c}
             *    \theta k_x\\
             *    \theta k_y\\
             *    \theta k_z
             * \end{array}\right]
             * @f$
             * @param thetakx [in] @f$ \theta k_x @f$
             * @param thetaky [in] @f$ \theta k_y @f$
             * @param thetakz [in] @f$ \theta k_z @f$
             */
            EAA(T thetakx, T thetaky, T thetakz) :
                    _eaa(Vector3D<T>(thetakx, thetaky, thetakz)){}

            /**
             * @brief Constructs an EAA vector that will rotate v1 into
             * v2. Where v1 and v2 are normalized and described in the same reference frame.
             * @param v1 [in] normalized vector
             * @param v2 [in] normalized vector
             */
            EAA(const Vector3D<T>& v1, const Vector3D<T>& v2) :
                    _eaa(0,0,0)
            {
                /* HMM, this seems incorrect. The projection of v1 onto v2 in dot(v1,v2) will be close
                 * to zero when the vectors are perpendicular, and not parallel. (JAJ, 06-09-2010)
                 *
                const T epsilon = (T)0.00001;
                T dval = dot(v1,v2);
                if(fabs(dval)<epsilon){
                    // if the angle is 0 then do nothing, if its 180 degrees then the
                    // rotation axis must be choosen to be perpendicular to v1 or v2
                    if(dval<epsilon){
                        _eaa = Vector3D<T>(v1(2)*(T)Pi,v1(0)*(T)Pi,v1(1)*(T)Pi);
                    }
                } else {
                    T cosangle = acos( dval );
                    _eaa = normalize( cross(v1,v2) )*cosangle;
                }
                */

                const T epsilon = (T)1e-15;
                T dval = dot(v1,v2);
                if(fabs(dval-1)<epsilon) {
                    /*
                     * if the projection is close to 1 then the angle between the vectors are almost 0
                     * and we cannot reliably determine the perpendicular axis.
                     *
                     * A good approximation is therefore just to set the EAA equal to 0.
                     * */
                    _eaa = Vector3D<T>(0,0,0);
                } else if(fabs(dval+1)<epsilon){
                    // if the projection is close to -1 then the angle between the vectors are almost 180 and we choose
                    // a rotation axis perpendicular to the vector
                    int idx = 0;
                    if( fabs(v1(0))>fabs(v1(1)) )
                        idx = 1;
                    if( fabs(v1(idx))>fabs(v1(2)) )
                        idx = 2;
                    Vector3D<T> v3(0,0,0);
                    v3(idx) = 1;
                    _eaa = normalize( cross(v1,v3) ) * (T)Pi;
                } else {
                    T cosangle = acos( dval );
                    _eaa = normalize( cross(v1,v2) )*cosangle;
                }

            }

            /**
             * @brief Constructs an initialized EAA vector
             *
             * The angle of the EAA are \f$\|eaa\|\f$ and the axis is \f$\frac{eaa}{\|eaa\|}\f$
             * @param eaa [in] Values to initialize the EAA
             */
            explicit EAA(Vector3D<T> eaa) : _eaa(eaa) {}

            //! @brief destructor
            virtual ~EAA(){}

            /**
             * @copydoc Rotation3DVector::toRotation3D()
             *
             * @f$
             * \mathbf{R} = e^{[\mathbf{\hat{k}}],\theta}=\mathbf{I}^{3x3}+[\mathbf{\hat{k}}] sin\theta+[{\mathbf{\hat{k}}}]^2(1-cos\theta) =
             *  \left[
             *    \begin{array}{ccc}
             *      k_xk_xv\theta + c\theta & k_xk_yv\theta - k_zs\theta & k_xk_zv\theta + k_ys\theta \\
             *      k_xk_yv\theta + k_zs\theta & k_yk_yv\theta + c\theta & k_yk_zv\theta - k_xs\theta\\
             *      k_xk_zv\theta - k_ys\theta & k_yk_zv\theta + k_xs\theta & k_zk_zv\theta + c\theta
             *    \end{array}
             *  \right]
             * @f$
             *
             * where:
             * - @f$ c\theta = cos \theta @f$
             * - @f$ s\theta = sin \theta @f$
             * - @f$ v\theta = 1-cos \theta @f$
             */
            virtual const Rotation3D<T> toRotation3D() const
            {
                T theta = angle();
                T ca = cos(theta);
                T sa = sin(theta);
                T va = 1-ca;

                Vector3D<T> k = axis();
                T kx = k[0];
                T ky = k[1];
                T kz = k[2];

                return Rotation3D<T>(
                        kx * kx * va + ca , kx * ky * va - kz * sa , kx * kz * va + ky * sa,
                        kx * ky * va + kz * sa, ky * ky * va + ca, ky * kz * va - kx * sa,
                        kx * kz * va - ky * sa, ky * kz * va + kx * sa, kz * kz * va + ca);
            }

            /**
             * @brief Extracts the angle of rotation @f$ \theta @f$
             * @return @f$ \theta @f$
             */
            T angle() const
            {
                return norm_2(_eaa.m());
            }

            /**
             * @brief Extracts the axis of rotation vector @f$ \mathbf{\hat{\mathbf{k}}} @f$
             * @return @f$ \mathbf{\hat{\mathbf{k}}} @f$
             */
            const Vector3D<T> axis() const
            {
                T theta = angle();
                if (theta < 1e-6)
                    return Vector3D<T>(0, 0, 0);
                else
                    return _eaa / theta;
            }

            /**
             * @brief Returns element of EAA
             * @param i [in] index (@f$ 0 < i < 3 @f$)
             * @return the @f$ i @f$'th element
             */
            const T& operator[](size_t i) const
            {
                assert(i < 3);
                return _eaa[i];
            }

            /**
             * @brief Returns element of EAA
             * @param i [in] index (@f$ 0 < i < 3 @f$)
             * @return the @f$ i @f$'th element
             */
            T& operator[](size_t i)
            {
                assert(i < 3);
                return _eaa[i];
            }


            /**
             * @brief Returns element of EAA
             * @param i [in] index (@f$ 0 < i < 3 @f$)
             * @return the @f$ i @f$'th element
             */
            const T& operator()(size_t i) const
            {
                assert(i < 3);
                return _eaa[i];
            }

            /**
             * @brief Returns element of EAA
             * @param i [in] index (@f$ 0 < i < 3 @f$)
             * @return the @f$ i @f$'th element
             */
            T& operator()(size_t i)
            {
                assert(i < 3);
                return _eaa[i];
            }

            /**
             * @brief Comparison operator.
             *
             * The comparison operator makes a element wise comparison.
             * Returns true only if all elements are equal.
             *
             * @param rhs [in] EAA to compare with
             * @return True if equal.
             */
            bool operator==(const EAA<T> &rhs) const
            {
                return (_eaa(0) == rhs(0) && _eaa(1) == rhs(1) && _eaa(2) == rhs(2));
            }

            /**
             * @brief Comparison operator.
             *
             * The comparison operator makes a element wise comparison.
             * Returns true if any of the elements are different.
             *
             * @param rhs [in] EAA to compare with
             * @return True if not equal.
             */
            bool operator!=(const EAA<T> &rhs) const
            {
                return !(*this == rhs);
            }

            /**
             * @brief Get the size of the EAA.
             * @return the size (always 3).
             */
            size_t size() const { return 3; }

            /**
             * @brief Calculates \f$ \robabx{a}{c}{\thetak} =
             * \robabx{a}{b}{\mathbf{R}} \robabx{b}{c}{\mathbf{\thetak}} \f$
             *
             * @param aRb [in] \f$ \robabx{a}{b}{\mathbf{R}} \f$
             * @param bTKc [in] \f$ \robabx{b}{c}{\thetak} \f$
             * @return \f$ \robabx{a}{c}{\thetak} \f$
             */
            friend const EAA operator*(const Rotation3D<T>& aRb, const EAA& bTKc)
            {
                return EAA(aRb * bTKc._eaa);
                /* return Vector3D<T>(prod(aRb.m(), bTKc._eaa.m()))); */
            }

            /**
             * @brief Ouputs EAA to stream
             * @param os [in/out] stream to use
             * @param eaa [in] equivalent axis-angle
             * @return the resulting stream
             */
            friend std::ostream& operator<<(std::ostream& os, const EAA<T>& eaa){
                return os <<" EAA { "<<eaa(0)<<", "<<eaa(1)<<", "<<eaa(2)<<"}";
                //return os << eaa._eaa;
            }

            /**
             * @brief Calculates the cross product
             * @param v [in] a 3D vector
             * @param eaa [in] a 3D eaa vector
             * @return the resulting 3D vector
             */
            friend const Vector3D<T> cross<T>(const Vector3D<T>& v, const EAA<T>& eaa);

        private:
            Vector3D<T> _eaa;
        };

        /**
        * @brief Calculates the cross product
        * @param v [in] a 3D vector
        * @param eaa [in] a 3D eaa vector
        * @return the resulting 3D vector
        */
        template<class T>
        const Vector3D<T> cross(const Vector3D<T>& v, const EAA<T>& eaa)
        {
            return cross(v, eaa._eaa);
        }

        /**
        * @brief Casts EAA<T> to EAA<Q>
        * @param eaa [in] EAA with type T
        * @return EAA with type Q
        */
        template<class Q, class T>
        const EAA<Q> cast(const EAA<T>& eaa) {
            return EAA<Q>(
                    static_cast<Q>(eaa(0)),
                    static_cast<Q>(eaa(1)),
                    static_cast<Q>(eaa(2)));
        }

        template class le::math::EAA<double>;
        template class le::math::EAA<float>;

        /*@}*/

    }} // end namespaces

#endif //LE_EAA_H
