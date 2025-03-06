#include "../include/GLProgram.h"
#include <string>
#include "../include/Algorithm.h"

#define WINDOW_WIDTH 1600
#define WINDOW_HEIGHT 1200

// declare static members for use in callback functions
int GLProgram::windowWidth = WINDOW_WIDTH;
int GLProgram::windowHeight = WINDOW_HEIGHT;
Camera GLProgram::camera;
bool GLProgram::mousePressed = false;
double GLProgram::prevMouseX, GLProgram::prevMouseY;
glm::mat4 GLProgram::modelMatrix = glm::mat4(1.0f);

typedef glm::vec<3, float> tvecn;
typedef glm::vec<4, float> tvecnp1;

int main()
{
    GLProgram program;
    tinynurbs::RationalCurve<float> contour;
    contour.control_points = {
        glm::fvec3(0, -0.5, 0.5), glm::fvec3(0, -0.5, 1), glm::fvec3(0, 1, 1),
        glm::fvec3(0, 1, 0.7), glm::fvec3(0, 0, 0.5),
        glm::fvec3(0, 1, 0.3), glm::fvec3(0, 1, 0), glm::fvec3(0, -0.5, 0),
        glm::fvec3(0, -0.5, 0.5)
    };
    contour.knots = {0, 0, 0, 0.15, 0.3, 0.4, 0.6, 0.7, 0.85, 1, 1, 1};
    contour.degree = 2;
    contour.weights = {1, 1, 1, 1, 1, 1, 1, 1, 1};

    tinynurbs::RationalCurve<float> trace;
    trace.control_points = {
        glm::fvec3(0, 0, 0), glm::fvec3(1, 0.5, 0.5), glm::fvec3(2, 1, 0),
        glm::fvec3(3, 0.5, -0.5), glm::fvec3(4, 0, 0), glm::fvec3(5, -1, 1)
    };
    // reverse(trace.control_points.begin(), trace.control_points.end());
    trace.knots = {0, 0, 0, 0, 0.33, 0.66, 1, 1, 1, 1};
    trace.degree = 3;
    trace.weights = {0.9, 0.8, 0.7, 0.6, 0.5, 0.5};
    // reverse(trace.weights.begin(), trace.weights.end());

    std::vector<tinynurbs::RationalCurve<float> > profile_curves = {
        tinynurbs::RationalCurve<float>(2, {0, 0, 0, 0.15, 0.3, 0.4, 0.6, 0.7, 0.85, 1, 1, 1,},
                                        {
                                            glm::fvec3(0, -0.5, 0.5), glm::fvec3(0, -0.5, 1), glm::fvec3(0, 1, 1),
                                            glm::fvec3(0, 1, 0.7),
                                            glm::fvec3(0, 0, 0.5), glm::fvec3(0, 1, 0.3), glm::fvec3(0, 1, 0),
                                            glm::fvec3(0, -0.5, 0),
                                            glm::fvec3(0, -0.5, 0.5)
                                        },
                                        {1, 1, 1, 1, 1, 1, 1, 1, 1}),
        tinynurbs::RationalCurve<float>(3, {0, 0, 0, 0, 0.2, 0.4, 0.6, 0.8, 1, 1, 1, 1},
                                        {
                                            glm::fvec3(0, 0, 0),
                                            glm::fvec3(1, 2, 0),
                                            glm::fvec3(2, -1, 0),
                                            glm::fvec3(3, 1, 0),
                                            glm::fvec3(4, -2, 0),
                                            glm::fvec3(5, 2, 0),
                                            glm::fvec3(6, -1, 0),
                                            glm::fvec3(7, 0, 0),
                                        },
                                        {1.0, 1.1, 1.0, 0.9, 0.9, 1.0, 1.0, 1.0}),
        tinynurbs::RationalCurve<float>(2, {0, 0, 0, 0.15, 0.3, 0.4, 0.6, 0.7, 0.85, 1, 1, 1,},
                                        {
                                            glm::fvec3(2.11214, 0.300734, 0.445626),
                                            glm::fvec3(2.27571, 0.737865, 0.266298),
                                            glm::fvec3(2.26004, 0.173603, -1.12344),
                                            glm::fvec3(2.1619, -0.0886756, -1.01584),
                                            glm::fvec3(2.10691, 0.112647, -0.0176185),
                                            glm::fvec3(2.03104, -0.438381, -0.872377),
                                            glm::fvec3(1.93289, -0.700659, -0.76478),
                                            glm::fvec3(1.94856, -0.136397, 0.624955),
                                            glm::fvec3(2.11214, 0.300734, 0.445626),
                                        },
                                        {1, 1, 1, 1, 1, 1, 1, 1, 1,}),
        tinynurbs::RationalCurve<float>(2, {0, 0, 0, 0.15, 0.3, 0.4, 0.6, 0.7, 0.85, 1, 1, 1,},
                                        {
                                            glm::fvec3(2.78925, 0.103418, 0.213457),
                                            glm::fvec3(2.98301, 0.159616, -0.244035),
                                            glm::fvec3(2.93749, -1.32605, -0.445812),
                                            glm::fvec3(2.82124, -1.35977, -0.171317),
                                            glm::fvec3(2.77408, -0.391805, 0.146198),
                                            glm::fvec3(2.66623, -1.40473, 0.194676),
                                            glm::fvec3(2.54997, -1.43845, 0.469171),
                                            glm::fvec3(2.59549, 0.04722, 0.670949),
                                            glm::fvec3(2.78925, 0.103418, 0.213457),
                                        },
                                        {1, 1, 1, 1, 1, 1, 1, 1, 1,}),
        tinynurbs::RationalCurve<float>(2, {0, 0, 0, 0.15, 0.3, 0.4, 0.6, 0.7, 0.85, 1, 1, 1,},
                                        {
                                            glm::fvec3(3.60077, -0.33267, 0.378142),
                                            glm::fvec3(4.00971, -0.335648, 0.0904634),
                                            glm::fvec3(3.91209, -1.82738, -0.0328582),
                                            glm::fvec3(3.66673, -1.82559, 0.139749),
                                            glm::fvec3(3.56823, -0.829914, 0.337035),
                                            glm::fvec3(3.33958, -1.82321, 0.369892),
                                            glm::fvec3(3.09421, -1.82142, 0.542499),
                                            glm::fvec3(3.19183, -0.329692, 0.66582),
                                            glm::fvec3(3.60077, -0.33267, 0.378142),
                                        },
                                        {1, 1, 1, 1, 1, 1, 1, 1, 1,}),
        tinynurbs::RationalCurve<float>(2, {0, 0, 0, 0.15, 0.3, 0.4, 0.6, 0.7, 0.85, 1, 1, 1,},
                                        {
                                            glm::fvec3(5, -1.5, 1.5), glm::fvec3(5.49053, -1.43151, 1.43151),
                                            glm::fvec3(5.62048, -2.84555, 0.948186),
                                            glm::fvec3(5.32616, -2.88664, 0.989278),
                                            glm::fvec3(5.04332, -1.97135, 1.33889),
                                            glm::fvec3(4.93373, -2.94144, 1.04407),
                                            glm::fvec3(4.63942, -2.98253, 1.08516),
                                            glm::fvec3(4.50947, -1.56849, 1.56849), glm::fvec3(5, -1.5, 1.5),
                                        },
                                        {1, 1, 1, 1, 1, 1, 1, 1, 1,}),
    };

    std::vector<tinynurbs::RationalCurve<float> > trace_curves = {
        tinynurbs::RationalCurve<float>(3, {0, 0, 0, 0, 0.33, 0.66, 1, 1, 1, 1},
                                        {
                                            glm::fvec3(0, 0, 0), glm::fvec3(1, 0.5, 0.5), glm::fvec3(2, 1, 0),
                                            glm::fvec3(3, 0.5, -0.5), glm::fvec3(4, 0, 0), glm::fvec3(5, -1, 1)
                                        },
                                        {1, 1, 1, 1, 1, 1}),
        tinynurbs::RationalCurve<float>(2, {0, 0, 0, 0.2, 0.4, 0.6, 0.8, 1, 1, 1},
                                        {
                                            glm::fvec3(2, 3, 4), glm::fvec3(0, 3, 4), glm::fvec3(0, 3, 0),
                                            glm::fvec3(3, 3, 0), glm::fvec3(3, 0, 0), glm::fvec3(3, 0, 2),
                                            glm::fvec3(5, 0, 2)
                                        },
                                        {1.0, 1.2, 1.2, 1.0, 1.0, 1.0, 1.0}),
        tinynurbs::RationalCurve<float>(3, {0, 0, 0, 0, 0.5, 0.5, 0.5, 0.75, 1, 1, 1, 1},
                                        {
                                            glm::fvec3(0, 0, 0), glm::fvec3(1, 0.5, 0.5), glm::fvec3(2, 0, 0.5),
                                            glm::fvec3(2, -1, 1), glm::fvec3(2, -1, 1), glm::fvec3(4, 0, 2),
                                            glm::fvec3(5, 1, 1.5), glm::fvec3(6, 0, 1),
                                        },
                                        {1.0, 1.2, 1.2, 1.0, 1.0, 1.0, 1.0, 1.1}),
        tinynurbs::RationalCurve<float>(2, {0, 0, 0, 0.5, 0.5, 0.75, 1, 1, 1},
                                        {
                                            glm::fvec3(0, 0, 0), glm::fvec3(2, 1, 1), glm::fvec3(4, 0, 2),
                                            glm::fvec3(4, 0, 2), glm::fvec3(6, -1, 1), glm::fvec3(8, 0, 0)
                                        },
                                        {0.9, 0.8, 0.7, 0.6, 0.5, 0.5})
    };

    std::vector<float> v_bar = {
        0, 0.116865, 0.206807, 0.311239, 0.412281, 0.487043, 0.570073, 0.630555, 0.704791,
        0.8212, 1,
    };
    std::vector<std::vector<glm::vec3> > frames(v_bar.size());
    std::fill(frames.begin(), frames.end(),
              std::vector<glm::vec3>({glm::fvec3(1, 0, 0), glm::fvec3(0, 1, 0), glm::fvec3(0, 0, 1)}));

    std::vector<tinynurbs::RationalSurface<float> > mySweepSurf(10);
    std::vector<tinynurbs::RationalCurve<float> > section(10);
    std::vector<std::vector<glm::vec3> > frames1{};
    std::vector<std::vector<glm::vec3> > frames2{};
    std::vector<std::vector<glm::vec3> > frames3{};
    std::vector<float> v_bar1{};
    std::vector<float> v_bar2{};
    std::vector<float> v_bar3{};

    // approximate method
    createSweepSurfaceWithInterpolation(profile_curves[1], trace_curves[3], 6, mySweepSurf[0],
                                        v_bar1, frames1, 2, 2);
    // offset method
    createSweepSurfaceWithOffset(profile_curves[1], trace_curves[3], 0.3f, 10, mySweepSurf[1], v_bar2, frames2, 1, 2);
    // skinning method
    createSweepSurfaceWithSkinning(profile_curves[1], trace_curves[3], 6, mySweepSurf[2],
                                   v_bar3, frames3, 2, 0);

    if (hasSelfIntersection(mySweepSurf[2]))
    {
        std::cout << "Self-intersection detected!" << std::endl;
    }
    else
    {
        std::cout << "No self-intersection detected!" << std::endl;
    }

    program.init({profile_curves[1]}, {trace_curves[3]}, {profile_curves}, {frames1}, {mySweepSurf[2]});
    program.setClearColor(0.05f, 0.18f, 0.25f, 1.0f);
    program.run({profile_curves[1]}, {trace_curves[3]}, {profile_curves}, {frames1}, {v_bar1}, {mySweepSurf[2]});
    program.cleanup();

    return 0;
}
