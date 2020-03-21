//
// Created by Sail Yang on 2019-09-14.
//

#ifndef PROJECT_ROTATION3DVECTOR_H
#define PROJECT_ROTATION3DVECTOR_H

namespace le { namespace math {

        /** @addtogroup math */
        /* @{*/

        /**
         * @brief An abstract base class for Rotation3D parameterisations
         *
         * Classes that represents a parametrisation of a 3D rotation may inherit
         * from this class
         */
        template<class T = double>
        class Rotation3DVector{
        public:
            /**
             * @brief Virtual destructor
             */
            virtual ~Rotation3DVector(){}

            /**
             * @brief Returns the corresponding @f$ 3\times 3 @f$ Rotation matrix
             * @return The rotation matrix
             */
            virtual const Rotation3D<T> toRotation3D() const = 0;

        protected:
            /**
             * @brief Copy Constructor
             *
             * We allow subclasses of this class to be copied.
             */
            Rotation3DVector(const Rotation3DVector&) {}

            /**
             * @brief Assignment operator is protected to force subclasses to
             * implement it by themself.
             */
            Rotation3DVector& operator=(const Rotation3DVector&) { return *this;}

        protected:
            /**
             * @brief Default Constructor
             */
            Rotation3DVector() {}
        };

        template class le::math::Rotation3DVector<double>;
        template class le::math::Rotation3DVector<float>;

        /**@}*/
    }} // end namespaces

#endif //PROJECT_ROTATION3DVECTOR_H
