//
// Created by Sail Yang on 2019-09-08.
//

#include "le/math/Q.h"
#include "le/math/Transform3D.h"
#include "le/math/RPY.h"

#include <iostream>
#include <math.h>

int main(){
    le::math::Q speed(6,1.0);
    le::math::Transform3D<> trans(le::math::RPY<>(0.1,0.1,0.5).toRotation3D());
    le::math::Vector3D<> eye(0.1,0.2,0.7);
    le::math::Rotation3D<> rot = le::math::RPY<>(0.1,0.1,0.5).toRotation3D();
    le::math::Transform3D<> invTrans = le::math::inverse(le::math::Transform3D<>(rot*eye,rot));
    std::cout << speed/M_PI*180.0 << "\n";
    std::cout << trans << "\n";
    std::cout << invTrans << "\n";
    return 0;
}

