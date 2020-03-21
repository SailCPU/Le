//
// Created by sail on 19-9-7.
//

#ifndef LE_Q_H
#define LE_Q_H

#endif //LE_Q_H

#include "Eigen/Eigen"

namespace le {
    namespace math {

        class Q {
        public:
            //! Eigen vector used as internal datastructure.
            typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Base;

            /**
             * @brief A configuration of vector of length \b dim.
             */
            explicit Q(size_t dim) : _vec(dim) {}

            /**
             * @brief Default constructor.
             *
             * The vector will be of dimension zero.
             */
            Q() : _vec((Base::Index) 0) {}

            /**
             * @brief Creates a Q of length \b n and initialized with values from \b values
             *
             * The method reads n values from \b values and do not check whether reading out of bounds.
             *
             * @param n [in] Length of q.
             * @param values [in] Values to initialize with
             */
            Q(size_t n, const double *values):_vec(n){
                for (size_t i = 0; i<n; i++)
                    _vec(i) = values[i];
            }

            /**
             * @brief Creates a Q of length \b n and initialize all values in Q to \b value
             *
             * @param n [in] Length of q.
             * @param value [in] Value to initialize
             */
            Q(size_t n, double value):_vec(n){
                for (size_t i = 0; i<n; i++)
                    _vec(i) = value;
            }

            /**
             * @brief Creates a Q of length \b n and initialize all values in Q to the values specified after \b n
             *
             * The number of arguments after \b n must match the number n.
             *
             * @param n [in] Length of q.
             * @param a0 [in] Value to initialize q(0)
             * @param a1 [in] Value to initialize q(1)
             * @param ... [in] Values to initialize [q(2);q(n-1)]
             *
             */

            Q(size_t n, double a0, double a1):_vec(n){
                if(n<2) assert("Vector size must be >= 2");
                _vec[0] = a0; _vec[1] = a1;
            }

            //! @copydoc Q(size_t,double,double)
            Q(size_t n, double a0, double a1, double a2):_vec(n){
                if(n<3) assert("Vector size must be >= 3");
                _vec[0] = a0; _vec[1] = a1; _vec[2] = a2;
            }

            //! @copydoc Q(size_t,double,double)
            Q(size_t n, double a0, double a1, double a2, double a3):_vec(n){
                if(n<4) assert("Vector size must be >= 4");
                _vec[0] = a0; _vec[1] = a1; _vec[2] = a2; _vec[3] = a3;
            }

            //! @copydoc Q(size_t,double,double)
            Q(size_t n, double a0, double a1, double a2, double a3, double a4):_vec(n){
                if(n<5) assert("Vector size must be >= 5");
                _vec[0] = a0; _vec[1] = a1; _vec[2] = a2; _vec[3] = a3; _vec[4] = a4;
            }

            //! @copydoc Q(size_t,double,double)
            Q(size_t n, double a0, double a1, double a2, double a3, double a4, double a5):_vec(n){
                if(n<6) assert("Vector size must be >= 6");
                _vec[0] = a0; _vec[1] = a1; _vec[2] = a2; _vec[3] = a3; _vec[4] = a4; _vec[5] = a5;
            }

            //! @copydoc Q(size_t,double,double)
            Q(size_t n, double a0, double a1, double a2, double a3, double a4, double a5, double a6):_vec(n){
                if(n<7) assert("Vector size must be >= 7");
                _vec[0] = a0; _vec[1] = a1; _vec[2] = a2; _vec[3] = a3; _vec[4] = a4; _vec[5] = a5; _vec[6] = a6;
            }

            //! @copydoc Q(size_t,double,double)
            Q(size_t n, double a0, double a1, double a2, double a3, double a4, double a5, double a6, double a7):_vec(n){
                if(n<8) assert("Vector size must be >= 8");
                _vec[0] = a0; _vec[1] = a1; _vec[2] = a2; _vec[3] = a3; _vec[4] = a4; _vec[5] = a5; _vec[6] = a6; _vec[7] = a7;
            }

            //! @copydoc Q(size_t,double,double)
            Q(size_t n, double a0, double a1, double a2, double a3, double a4, double a5, double a6, double a7,
              double a8):_vec(n){
                if(n<9) assert("Vector size must be >= 9");
                _vec[0] = a0; _vec[1] = a1; _vec[2] = a2; _vec[3] = a3; _vec[4] = a4; _vec[5] = a5; _vec[6] = a6; _vec[7] = a7; _vec[8] = a8;

            }

            //! @copydoc Q(size_t,double,double)
            Q(size_t n, double a0, double a1, double a2, double a3, double a4, double a5, double a6, double a7,double a8, double a9):_vec(n){
                if(n<10) assert("Vector size must be >= 10");
                _vec[0] = a0; _vec[1] = a1; _vec[2] = a2; _vec[3] = a3; _vec[4] = a4; _vec[5] = a5; _vec[6] = a6; _vec[7] = a7; _vec[8] = a8; _vec[9] = a9;

            }

            /**
             * @brief Construct a configuration vector from a std::vector
             * expression.
             *
             * @param r [in] An expression for a vector of doubles
             */
            Q(const std::vector<double> &r) :
                    _vec(r.size()) {
                for (size_t i = 0; i < r.size(); i++)
                    _vec(i) = r[i];
            }

            /**
             * @brief Construct from Eigen base.
             * @param q [in] Eigen base.
             */
            Q(const Base &q) :
                    _vec(q.rows()) {
                for (int i = 0; i < q.size(); i++)
                    _vec(i) = q(i, 0);
            }

            //! @brief Destructor.
            virtual ~Q(){}

            /**
             * @brief Returns Q of length \b n initialized with 0's
             */
            static Q zero(std::size_t n) {
                return Q(Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(n));
            }

            /**
             * @brief The dimension of the configuration vector.
             */
            size_t size() const {
                return _vec.rows();
            }


            /**
               @brief True if the configuration is of dimension zero.
             */
            bool empty() const { return size() == 0; }


            /**
             * @brief Accessor for the internal Eigen vector state.
             */
            const Base &e() const {
                return _vec;
            }

            /**
             * @brief Accessor for the internal Eigen vector state.
             */
            Base &e() {
                return _vec;
            }

            /**
             * @brief Extracts a sub part (range) of this Q.
             * @param start [in] Start index
             * @param cnt [in] the number of elements to include
             * @return
             */
            const Q getSubPart(size_t start, size_t cnt) const {
                assert(start + cnt <= size());

                Q res(cnt);
                for (size_t i = 0; i < cnt; i++) {
                    res(i) = (*this)[start + i];
                }
                return res;
            }

            /**
             * @brief Set subpart of vector.
             * @param index [in] the initial index.
             * @param part [in] the part to insert beginning from \b index.
             */
            void setSubPart(size_t index, const Q &part) {
                assert(index + part.size() <= size());
                for (size_t i = 0; i < part.size(); i++) {
                    (*this)[index + i] = part(i);
                }
            }

            //----------------------------------------------------------------------
            // Norm utility methods

            /**
             * @brief Returns the Euclidean norm (2-norm) of the configuration
             * @return the norm
             */
            double norm2() const {
                return _vec.norm();
                //return norm_2(m());
            }

            /**
             * @brief Returns the Manhatten norm (1-norm) of the configuration
             * @return the norm
             */
            double norm1() const {
                return _vec.lpNorm<1>();
                //return norm_1(m());
            }

            /**
             * @brief Returns the infinte norm (\f$\inf\f$-norm) of the configuration
             * @return the norm
             */
            double normInf() const {
                return _vec.lpNorm<Eigen::Infinity>();
                //return norm_inf(m());
            }

            //----------------------------------------------------------------------
            // Various operators

            /**
             * @brief Returns reference to vector element
             * @param i [in] index in the vector
             * @return const reference to element
             */
            const double &operator()(size_t i) const { return _vec(i); }

            /**
             * @brief Returns reference to vector element
             * @param i [in] index in the vector
             * @return reference to element
             */
            double &operator()(size_t i) { return _vec(i); }

            /**
             * @brief Returns reference to vector element
             * @param i [in] index in the vector
             * @return const reference to element
             */
            const double &operator[](size_t i) const { return _vec(i); }

            /**
             * @brief Returns reference to vector element
             * @param i [in] index in the vector
             * @return reference to element
             */
            double &operator[](size_t i) { return _vec(i); }

            /**
               @brief Scalar division.
             */
            const Q operator/(double s) const {
                return Q(_vec / s);
            }

            /**
             * @brief Scalar multiplication.
             */
            const Q operator*(double s) const {
                return Q(_vec * s);
            }

            /**
             * @brief Scalar multiplication.
             */
            friend const Q operator*(double s, const Q &v) {
                return Q(s * v.e());
            }

            /**
             * @brief Scalar division.
             */
            friend const Q operator/(double s, const Q &v) {
                Q res = v;
                for (size_t i = 0; i < v.size(); i++)
                    res(i) = s / v(i);
                return res;
            }

            /**
             * @brief Vector subtraction.
             */
            const Q operator-(const Q &b) const {
                return Q(_vec - b.e());
            }

            /**
             * @brief Vector addition.
             */
            const Q operator+(const Q &b) const {
                return Q(_vec + b.e());
            }

            /**
             * @brief Scalar multiplication.
             */
            Q &operator*=(double s) {
                _vec *= s;
                return *this;
            }

            /**
             * @brief Scalar division.
             */
            Q &operator/=(double s) {
                _vec /= s;
                return *this;
            }

            /**
             * @brief Vector addition.
             */
            Q &operator+=(const Q &v) {
                _vec += v.e();
                return *this;
            }

            /**
             * @brief Vector subtraction.
             */
            Q &operator-=(const Q &v) {
                _vec -= v.e();
                return *this;
            }

            /**
             * @brief Unary minus.
             */
            Q operator-() const {
                return Q(-_vec);
            }

            /**
             * @brief Compares whether this is less than \b q
             *
             * The less operator is defined such that the first index is the most significant. That is
             * if (*this)[0] < q[0] then true is returned. If (*this)[0] > q[0] false is returned and
             * only if (*this)[0] == q[0] is the next index considered.
             */
            bool operator<(const Q &q) const {
                assert(size() == q.size());
                for (size_t i = 0; i < size(); i++) {
                    if (_vec[i] < q[i])
                        return true;
                    else if (_vec[i] > q[i])
                        return false;
                }
                return false;
            }

            /**
             * @brief Convert to a standard vector.
             * @param v [out] the result.
             */
            void toStdVector(std::vector<double> &v) const {
                v.resize(size());
                for (size_t i = 0; i < size(); i++) {
                    v[i] = _vec[i];
                }
            }

            /**
             * @brief Convert to a standard vector.
             * @return the result.
             */
            std::vector<double> toStdVector() const {
                std::vector<double> v(size());
                toStdVector(v);
                return v;
            }

        private:
            void init(size_t n, const double *values) {
                for (size_t i = 0; i < n; i++)
                    _vec(i) = values[i];
            }

        private:
            Base _vec;
        };

        /**
          * @brief Compares \b q1 and \b q2 for equality.
          *
          * \b q1 and \b q2 are considered equal if and only if they have equal
          * length and if q1(i) == q2(i) for all i.
          *
          * @relates Q
          *
          * @param q1 [in]
          * @param q2 [in]
          * @return True if q1 equals q2, false otherwise.
          */
        bool operator==(const Q& q1, const Q& q2){
            if (q1.size() != q2.size())
                return false;

            for (size_t i = 0; i < q1.size(); i++)
                if (q1(i) != q2(i))
                    return false;
            return true;
        }

        /**
           @brief Inequality operator

           The inverse of operator==().
         */
        inline bool operator!=(const Q& q1, const Q& q2) { return !(q1 == q2); }

        /**
         * @brief Streaming operator.
         *
         * @relates Q
         */
        std::ostream& operator<<(std::ostream& out, const Q& v){
            if (v.size() == 0)
                return out << "Q[0]{}";
            else {
                out << "Q[" << (int)v.size() << "]{";
                for (size_t i = 0; i < v.size() - 1; i++)
                    out << v[i] << ", ";
                return out << v[v.size() - 1] << "}";
            }
        }

        /**
         * @brief Input streaming operator
         *
         * Parse input stream according to how operator<< streams out
         *
         * @relates Q
         * @param in [in] Input stream
         * @param q [in] Target of q read in
         * @return reference to \b in
         */
        std::istream& operator>>(std::istream& in, Q& q){
            char ch1, ch2;
            do {
                in.get(ch1);
            } while (ch1 == ' ' || ch1 == '\t'); //Ignore space and tab, but not line changes.


            int size = -1;

            if (ch1 == 'Q') {
                in.get(ch2);
                if (ch1 != 'Q' || ch2 != '[')
                    assert("Content of input stream does not match format of Q");
                in >> size;

                in.get(ch1);
                in.get(ch2);
                if (ch1 != ']' || ch2 != '{')
                    assert("Content of input stream does not match format of Q");
            } else if (ch1 != '{') {
                assert("Content of input stream does not match format of Q");
            }

            std::vector<double> res;
            while (ch1 != '}') {
                double d;
                in >> d;
                if (!in.eof()) {
                    res.push_back(d);
                }
                in.get(ch1);
            }

            if (ch1 != '}')
                assert("Content of input stream does not match format of Q");

            if (size > -1 && (int)res.size() != size) {
                assert("Length of Q does not match device");
            }

            q = Q(res.size(), &res[0]);
            return in;
        }

        /**
           @brief The dot product (inner product) of \b a and \b b.

           @relates Q
        */
        double dot(const Q& a, const Q& b){
            return a.e().dot(b.e());
        }

        /**
         * @brief concatenates q1 onto q2 such that the returned q has
         * the configurations of q1 in [0;q1.size()[ and has q2 in
         * [q1.size();q1.size()+q2.size()[
         * @param q1 [in] the first Q
         * @param q2 [in] the second Q
         * @return the concatenation of q1 and q2
         */
        Q concat(const Q& q1, const Q& q2){
            Q q(q1.size()+q2.size());
            for(size_t i=0;i<q1.size();i++)
                q(i) = q1(i);
            for(size_t i=0;i<q2.size();i++)
                q(q1.size()+i) = q2(i);
            return q;
        }

        /*@}*/
    }
}