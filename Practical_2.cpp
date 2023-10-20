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

// setup output
std::ofstream output_file("./practical_1_output1.dat");

int nPoints;
int padding;
double x0;
double x1;
double tStart = 0.0;
double tStop;
double a;
double dx;
double dt;
double C;
//define a vector to store discretised data
vector<double> u;


void boundary_condition(vector<double> &values, int pad)
{
    for (int index = 0; index < pad; index++)
    {
        values[index] = values[index + nPoints];
        values[pad + nPoints + index] = values[pad + index];
    }
}

double (*u_0)(double);

double square_function(double x)
{
    if ((x >= 0.3) && (x <= 0.7))
    {
        return 1.0;
    }
    else
    {
        return 0.0;
    }
}

double Gauss_function(double x)
{
    return exp((-8)*x*x);
}

vector<double> (*PDE_Step)(const vector<double>&, double, int);

vector<double> forward_difference_Step(const vector<double> &current, double time_step, int pad)
{
    vector<double> uPlus1 = current;
    for (int index = pad; index < nPoints+pad; index++)
    {
        uPlus1[index] = current[index] - C*(current[index+1] - current[index]);
    }
    boundary_condition(uPlus1, pad);
    return uPlus1;
}

vector<double> backward_difference_Step(const vector<double> &current, double time_step, int pad)
{
    vector<double> uPlus1 = current;
    for (int index = pad; index < nPoints+pad; index++)
    {
        uPlus1[index] = current[index] - C*(current[index] - current[index-1]);
    }
    boundary_condition(uPlus1, pad);
    return uPlus1;
}

vector<double> upwind_Step(const vector<double> &current, double time_step, int pad)
{
    if (a >= 0)
    {
        return backward_difference_Step(current, time_step, pad);
    }
    else
    {
        return forward_difference_Step(current, time_step, pad);
    }
}

vector<double> Lax_Friedrichs_Step(const vector<double> &current, double time_step, int pad)
{
    vector<double> uPlus1 = current;
    cout << "Remember to implement Lax_Friedrichs_Step\n";
    return uPlus1;
}

vector<double> Lax_Wandroff_Step(const vector<double> &current, double time_step, int pad)
{
    vector<double> uPlus1 = current;
    cout << "Remember to implement Lax_Wandroff_Step\n";
    return uPlus1;
}

vector<double> Warming_Beam_Step(const vector<double> &current, double time_step, int pad)
{
    vector<double> uPlus1 = current;
    cout << "Remember to implement Warming_Beam_Step\n";
    return uPlus1;
}

void read()
{
    cout << "Time:\n";
    cin >> tStop;
    cout << "Number of Points:\n";
    cin >> nPoints;
    cout << "Equation:\n";
    cin >> a;
    cout << "Initial value:\n";
    string help;
    cin >> help;
    if (help == "square")
    {
        u_0 = &square_function;
        x0 = 0;
        x1 = 1;
    }
    else
    {
        if (help == "Gauss")
        {
            u_0 = &Gauss_function;
            x0 = -1;
            x1 = 1;
        }
        else
        {
            cout << "This initial value function has not been implemented (yet)\n";
            exit(1);
        }
    }
    cout << "Mathod:\n";
    cin >> help;
    if (help == "forward_difference")
    {
        PDE_Step = forward_difference_Step;
        padding = 1;
    }
    else
    {
        if (help == "backward_difference")
        {
            PDE_Step = backward_difference_Step;
            padding = 1;
        }
        else
        {
            if (help == "upwind")
            {
                PDE_Step = upwind_Step;
                padding = 1;
            }
            else
            {
                if (help == "Lax_Friedrichs")
                {
                    PDE_Step = Lax_Friedrichs_Step;
                    padding = 1;
                }
                else
                {
                    if (help == "Lax_Wandroff")
                    {
                        PDE_Step = Lax_Wandroff_Step;
                        padding = 1; 
                    }
                    else
                    {
                        if (help == "Warming_Beam")
                        {
                            PDE_Step = Warming_Beam_Step;
                            padding = 2;
                        }
                        else
                        {
                            cout << "This initial method has not been implemented (yet)\n";
                            exit(1);
                        }

                    }
        
                }
            }
        }
    }
    u.resize(nPoints + (2*padding));
    dx = (x1 - x0) / nPoints;
    dt = dx;
    C = (-1)*a*dt/dx;
}

void set_initial_value(int pad)
{
    double x;
    for (int i = 0; i < nPoints ; i ++)
    {
        x = x0 + (i * dx);
        u[i+pad] = u_0(x);
    }
    boundary_condition(u, pad);
}

int print_result(vector<double> data)
{
    for (int i = 0; i < nPoints; i++)
    {
        output_file << data[i+padding] << " ";
    }
    output_file << std::endl;
    return 0;
}


int main()
{
    read();
    set_initial_value(padding);

    double t = tStart;
    double local_dt;
    while (t < tStop)
    {
        local_dt = min(dt, tStop-tStart); // avoid overshooting the final time

        u = PDE_Step(u, local_dt, padding);

        t = t + local_dt;
    }

    print_result(u);
    output_file.close();
    return 0;
}
