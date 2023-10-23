#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <functional>
#include  <string>

using std::min;
using std::cout;
using std::cin;
using std::vector;
using std::string;

const double pi = 3.14159265358979311600;

int nPoints = 100;
int padding = 1;
double x0;
double x1;
double tStart = 0.0;
double tStop;
double a = 1;
double dx;
double C = 1;
//define a vector to store discretised data
vector<double> u(nPoints + 2*padding); //u[0] and u[nPoints+1] are padding
vector<double> f(nPoints + 2*padding - 1); //f[i] holds the flux from u[i] to u[i+1]

//u.resize(nPoints + 2*padding);
//f.resize(nPoints + 2*padding - 1);


void boundary_condition(vector<double> &values, int pad)
{
    for (int index = 0; index < pad; index++)
    {
        values[index] = values[pad];
        values[pad + nPoints + index] = values[pad + nPoints];
    }
}

double (*u_0)(double);

double (*flux)(double);

double (*nummerical_flux)(double, double, double, double);

double advection(double u)
{
    return a*u;
}

double Burgers(double u)
{
    return (0.5 * u * u);
}

double shock(double x)
{
    if (x <= 0.5)
    {
        return 2.0;
    }
    else
    {
        return 1.0;
    }
}

double rarefraction(double x)
{
    if (x <= 0.5)
    {
        return 1;
    }
    else
    {
        return 2;
    }
}

double Toro(double x)
{
    if (x <= 0.5)
    {
        return (-0.5);
    }
    else
    {
        if (x <= 1)
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
}

double cosine(double x)
{
    return cos(2*pi*x);
}


vector<double> PDE_Step(const vector<double> &current, const vector<double> &flux, double time_step, int pad)
{
    vector<double> uPlus1 = current;
    for (int index = pad; index < (nPoints + pad); index++)
    {
        uPlus1[index] = current[index] - ((time_step/dx)*(flux[index] - flux[index-1]));
    }
    boundary_condition(uPlus1, pad);
    return uPlus1;
}

double Lax_Friedrichs_flux(const double velocity, const double velocity_next, double d_x, double d_t)
{
    double result = (d_x / d_t)*(velocity - velocity_next);
    result += flux(velocity);
    result += flux(velocity_next);
    return 0.5*result;
}

double Richtmyer_flux(const double velocity, const double velocity_next, double d_x, double d_t)
{
    double velocity_half_step = 0.5*(velocity + velocity_next);
    velocity_half_step += (-0.5)*(d_t/d_x)*(flux(velocity_next)-flux(velocity));
    return flux(velocity_half_step);
}

double FORCE_flux(const double velocity, const double velocity_next, double d_x, double d_t)
{
    double result = Lax_Friedrichs_flux(velocity, velocity_next, d_x, d_t);
    result += Richtmyer_flux(velocity, velocity_next, d_x, d_t);
    return 0.5*result;
}

double Godunov_flux(const double velocity, const double velocity_next, double d_x, double d_t)
{
    double s = (velocity + velocity_next)/2;
    double velocity_half;
    if (velocity > velocity_next)
    {
        if (s > 0)
        {
            velocity_half = velocity;
        }
        else
        {
            if (s < 0)
            {
                velocity_half = velocity_next;
            }
            else
            {
                std::cout << "Zero shock velocity\n";
                exit(3);
            }
        }
    }
    else
    {
        if (velocity > 0)
        {
            velocity_half = velocity;
        }
        else
        {
            if ((velocity <= 0)&&(0 <= velocity_next))
            {
                velocity_half = 0;
            }
            else
            {
                velocity_half = velocity_next;
            }
        }
    }
    return flux(velocity_half);
}

void set_initial_value(vector<double> &vel, int pad)
{
    double x;
    double integral = 0;
    for (int i = 0; i < nPoints ; i ++)
    {
        x = x0 + ((i + 0.5) * dx);
        // integral = 0.25 * u_0(x);
        // x += (dx/2);
        // integral += 0.5 * u_0(x);
        // x += (dx/2);
        // integral += 0.25 * u_0(x);
        vel[i+pad] = u_0(x);
    }
    boundary_condition(u, pad);
}

void evaluate_flux(const vector<double> &vel, vector<double> &flu, double space_step, double time_step)
{
    for (int i = 0; i < (nPoints+1); i++)
    {
        flu[i] = nummerical_flux(vel[i], vel[i+1], space_step, time_step);
    }
}

void print_result(vector<double> data)
{
    for (int i = 0; i < nPoints; i++)
    {
        std::cout << data[i+padding] << " ";
    }
    std::cout << std::endl;
}

double maximal(const vector<double>& values)
{
    int number = values.size();
    double result = 0;
    for (int i = 0; i < number; i++)
    {
        if (abs(values[i]) > result)
        {
            result = abs(values[i]);
        }
    }
    return result;
}

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        "Remember to specify the equation, initial condition, and method\n";
        exit(1);
    }
    std::cout << argv[1][0] << " " << argv[2] << " " << argv[3] << std::endl;
    // the first argument is the equation
    if (argv[1][0] == 'a') // advection
    {
        flux = advection;
        std::cout << "a?\n";
        std::cin >> a;
    }
    else
    {
        if (argv[1][0] == 'B') // Burgers
        {
            flux = Burgers;
        }
        else
        {
            std::cout << "Invalid equation\n";
            exit(2);
        }
    }
    // the second argument is the initial condition
    if (argv[2][0] == 's') // shock
    {
        u_0 = shock;
        tStop = 0.2;
        x0 = 0;
        x1 = 1;
    }
    else
    {
        if (argv[2][0] == 'r') // rarefraction
        {
            u_0 = rarefraction;
            tStop = 0.2;
            x0 = 0;
            x1 = 1;
        }
        else
        {
            if (argv[2][0] == 'T') // Toro
            {
                u_0 = Toro;
                tStop = 1.5;
                C = 0.8;
                x0 = 0;
                x1 = 1.5;
            }
            else
            {
                if (argv[2][0] == 'c') // cosine
                {
                    u_0 = cosine;
                    tStop = 0.2;
                    x0 = 0;
                    x1 = 1;
                }
                else
                {
                    std::cout << "Invalid initial condition\n";
                    exit(2);
                }
            }
        }
    }
    // the third argument is the numerical method
    if ((argv[3][0] == 'L')&&(argv[3][4] = 'F')) // Lax_Friedrichs
    {
        nummerical_flux = Lax_Friedrichs_flux;
    }
    else
    {
        if (argv[3][0] == 'F') // FORCE
        {
            nummerical_flux = FORCE_flux;
        }
        else
        {
            if (argv[3][0] == 'G') //Godunov
            {
                nummerical_flux = Godunov_flux;
            }
            else
            
            {
                std::cout << "Invalid method\n";
                exit(2);
            }
        }
    }

    dx = (x1 - x0) / nPoints;

    set_initial_value(u, padding);
    print_result(u);

    double t = tStart;
    double dt;
    double a_max;
    while (t < tStop)
    {
        a_max = maximal(u);
        dt = C * dx / a_max;
        // avoid overshooting the final time
        if ((t + dt) > tStop)
        {
            dt = tStop - t;
        }
        
        evaluate_flux(u, f, dx, dt);

        u = PDE_Step(u, f, dt, padding);

        t = t + dt;
    }

    print_result(u);
    return 0;
}
