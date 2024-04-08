#include "pr.hpp"
#include <iostream>
#include <memory>

template <class T>
inline T ss_max(T x, T y)
{
    return (x < y) ? y : x;
}
template <class T>
inline T ss_min(T x, T y)
{
    return (y < x) ? y : x;
}
template <class T>
inline T ss_abs(T x)
{
    return x < T(0) ? -x : x;
}

float sign_of(float x)
{
    return 0.0f < x ? 1.0f : -1.0f;
}

// ax^2 + bx + c == 0
int solve_quadratic( float xs[2], float a, float b, float c)
{
    // bx + c = 0 case
    if( a == 0.0f )
    {
        if( b == 0.0f )
        {
            return 0;
        }
        xs[0] = -c / b;
        return 1;
    }

    float det = b * b - 4.0f * a * c;
    if( det < 0.0f )
    {
        return 0;
    }

    float k = (-b - sign_of(b) * std::sqrtf(det)) / 2.0f;
    float x0 = k / a;
    float x1 = c / k;
    xs[0] = ss_min( x0, x1 );
    xs[1] = ss_max( x0, x1 );
    return 2;
}
float quadratic(float x, float a, float b, float c) {
    return x * x * a + x * b + c;
};
float cubic(float x, float a, float b, float c, float d) {
    return x * x * x * a + x * x * b + x * c + d;
};

// High-Performance Polynomial Root Finding for Graphics

// f (x)= a x^3 + b x^2 + c x + d
// f'(x)= aD x^2 + bD x + cD
float search_range( float LSign, float lBound, float rBound, int iterations, float a, float b, float c, float d, float aD, float bD, float cD )
{
    float x = (lBound + rBound) * 0.5f;
    float y = cubic( x, a, b, c, d );

    for (int i = 0; i < iterations; i++)
    {
        if( LSign * y < 0.0f )
        {
            rBound = x;
        }
        else
        {
            lBound = x;
        }

        //pr::DrawLine({ lBound, 0.1f + i * 0.1f, 0 }, { rBound,  0.1f + i * 0.1f, 0 }, {255, 0, 0}, 4);
        //pr::DrawSphere({ x, y, 0 }, 0.05f, { 128, 128, 255 });

        float dx = -y / quadratic(x, aD, bD, cD);
        float xn = x + dx;

        //pr::DrawLine({ x, y, 0 }, { xn, 0, 0 }, { 255, 255, 0 }, 2);

        if( xn < lBound || rBound < xn )
        {
            //pr::DrawLine({ x, y, 0 }, { xn, 0, 0 }, { 64, 64, 0 }, 1);
            xn = ( lBound + rBound ) * 0.5f;
            //pr::DrawLine({ x, y, 0 }, { xn, 0, 0 }, { 255, 255, 0 }, 2);
        }
        //else
        //{
        //    pr::DrawLine({ x, y, 0 }, { xn, 0, 0 }, { 255, 255, 0 }, 2);
        //}

        y = cubic(xn, a, b, c, d);
        x = xn;
    }
    return x;
}

int solve_cubic( float xs[3], float xMinCubic, float xMaxCubic, int iterations, float a, float b, float c, float d )
{
    if( a == 0.f )
    {
        return solve_quadratic(xs, b, c, d);
    }

    float aD = 3.0f * a;
    float bD = 2.0f * b;
    float cD = c;

    float borders[2];

    if( solve_quadratic( borders, aD, bD, cD ) ) // aD != 0.0f
    {
        borders[0] = ss_max( borders[0], xMinCubic );
        borders[1] = ss_min( borders[1], xMaxCubic );

        float yL = cubic(xMinCubic, a, b, c, d);
        float y0 = cubic(borders[0], a, b, c, d);
        float y1 = cubic(borders[1], a, b, c, d);
        float yR = cubic(xMaxCubic, a, b, c, d);

        float sY;
        float lBound;
        float rBound;

        if( yL * y0 < 0.0f ) // Left
        {
            sY = sign_of(yL);
            lBound = xMinCubic;
            rBound = borders[0];
        }
        else if( y0 * y1 < 0.0f ) // middle
        {
            sY = sign_of(y0);
            lBound = borders[0];
            rBound = borders[1];
        }
        else if( y1 * yR < 0.0f ) // Right
        {
            sY = sign_of(y1);
            lBound = borders[1];
            rBound = xMaxCubic;
        }
        else
        {
            return 0; 
        }
        
        // if there is one of solution in this range, the others can be found by deflation even if they are outside of this range.

        float x = search_range(sY, lBound, rBound, iterations, a, b, c, d, aD, bD, cD);
        xs[0] = x;

        float aPrime = a;
        float bPrime = b + x * aPrime;
        float cPrime = c + x * bPrime;
        if( solve_quadratic(xs + 1, aPrime, bPrime, cPrime) ) // aPrime != 0.0f
        {
            return 3;
        }
        return 1;
    }
    else
    {
        float yL = cubic(xMinCubic, a, b, c, d);
        float yR = cubic(xMaxCubic, a, b, c, d);
        if (yL * yR < 0.0f)
        {
            float lBound = xMinCubic;
            float rBound = xMaxCubic;
            xs[0] = search_range(sign_of(yL), lBound, rBound, iterations, a, b, c, d, aD, bD, cD);
            return 1;
        }
    }
    return 0;
}

int main() {
    using namespace pr;

    Config config;
    config.ScreenWidth = 1920;
    config.ScreenHeight = 1080;
    config.SwapInterval = 1;
    Initialize(config);

    Camera3D camera;
    camera.origin = { 0, 0, 4 };
    camera.lookat = { 0, 0, 0 };

    double e = GetElapsedTime();

    while (pr::NextFrame() == false) {
        if (IsImGuiUsingMouse() == false) {
            UpdateCameraBlenderLike(&camera);
        }

        ClearBackground(0.1f, 0.1f, 0.1f, 1);

        BeginCamera(camera);

        PushGraphicState();

        enum
        {
            SOLVER_QUADRATIC,
            SOLVER_CUBIC,
        };
        static int solver = SOLVER_CUBIC;

        static float a = 1.0f;
        static float b = 2.0f;
        static float c = -2.0f;
        static float d = -1.0f;

        static float xMinCubic = -16.f;
        static float xMaxCubic = +16.f;

        DrawGrid(GridAxis::XY, 1.0f, 32, { 128, 128, 128 });
        DrawXYZAxis(1.0f);


        if (solver == SOLVER_QUADRATIC)
        {
            PrimBegin( pr::PrimitiveMode::LineStrip );
            int N = 1000;
            LinearTransform itoX(0, N, -5, 5);
            for (int i = 0; i < N; i++)
            {
                float x = itoX(i);
                float y = quadratic(x, a, b, c);
                PrimVertex({ x, y, 0 }, { 255,255,255 });
            }
            PrimEnd();

            float xs[2];
            if (solve_quadratic(xs, a, b, c))
            {
                DrawText({ xs[0], 0, 0 }, "x0");
                DrawText({ xs[1], 0, 0 }, "x1" );

                DrawSphere({ xs[0], quadratic(xs[0], a, b, c), 0 }, 0.01f, { 255, 0, 0 });
                DrawSphere({ xs[1], quadratic(xs[1], a, b, c), 0 }, 0.01f, { 255, 0, 0 });
            }
        }

        if (solver == SOLVER_CUBIC)
        {
            PrimBegin(pr::PrimitiveMode::LineStrip);
            int N = 1000;
            LinearTransform itoX(0, N, xMinCubic, xMaxCubic);
            for (int i = 0; i < N; i++)
            {
                float x = itoX(i);
                float y = cubic(x, a, b, c, d);
                PrimVertex({ x, y, 0 }, { 255,255,255 });
            }
            PrimEnd();

            DrawLine({ xMinCubic, -16, 0 }, { xMinCubic, 16, 0 }, { 128, 0, 0 }, 4);
            DrawLine({ xMaxCubic, -16, 0 }, { xMaxCubic, 16, 0 }, { 128, 0, 0 }, 4);


            float xs[3];
            int nSolutions = solve_cubic(xs, xMinCubic, xMaxCubic, 16, a, b, c, d);
            
            for (int i = 0; i < nSolutions; i++)
            {
                char buffer[256];
                sprintf(buffer, "[%d] y=%f", i, cubic(xs[i], a, b, c, d));
                DrawText({ xs[i], 0, 0 }, buffer);
                DrawSphere({ xs[i], cubic(xs[i], a, b, c, d), 0}, 0.01f, {255, 0, 0});
            }
        }

        PopGraphicState();
        EndCamera();

        BeginImGui();

        ImGui::SetNextWindowSize({ 500, 800 }, ImGuiCond_Once);
        ImGui::Begin("Panel");
        ImGui::Text("fps = %f", GetFrameRate());

        ImGui::SliderFloat("a", &a, -10, 10);
        ImGui::SliderFloat("b", &b, -10, 10);
        ImGui::SliderFloat("c", &c, -10, 10);
        ImGui::SliderFloat("d", &d, -10, 10);

        ImGui::SliderFloat("xMinCubic", &xMinCubic, -16, 16);
        ImGui::SliderFloat("xMaxCubic", &xMaxCubic, -16, 16);

        ImGui::RadioButton("quadratic", &solver, SOLVER_QUADRATIC);
        ImGui::RadioButton("cubic", &solver, SOLVER_CUBIC);

        if (ImGui::Button("randomize"))
        {
            a = glm::mix( -2.0f, 2.0f, (float)rand() / RAND_MAX );
            b = glm::mix( -2.0f, 2.0f, (float)rand() / RAND_MAX );
            c = glm::mix( -2.0f, 2.0f, (float)rand() / RAND_MAX );
            d = glm::mix( -2.0f, 2.0f, (float)rand() / RAND_MAX );
        }

        ImGui::End();

        EndImGui();
    }

    pr::CleanUp();

    return 0;
}
