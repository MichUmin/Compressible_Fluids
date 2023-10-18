#include <vector>
#include <cmath>
#include <fstream>


using std::min;
using std::vector;

// setup output
std::ofstream output_file("./practical_1_output1.dat");

int nPoints = 100; // Chosen such that distance
int padding = 1;
double x0 = 0.0; // between points is 0.1
double x1 = 1.0;
double tStart = 0.0;
double tStop = 1.0;
double a = -1.0;
double dx = (x1 - x0) / nPoints;
// See problem sheet for more options here:
double dt = dx;
//define a vector with enough space to store discretised data
// at each point
std::vector<double> u (nPoints + (2*padding));
std::vector<double> uPlus1 = u;

double u_0(double x)
{
    return sin(x);
}

void set_initial_value()
{
    
    double x;
    for (int i = 0; i < padding ; i ++)
    {
        x = x0 + (i + nPoints - padding)*dx;
        u[i] = u_0(x);
        x = x0 + i*dx;
        u[i+padding+nPoints] = u_0(x);
    }
    for (int i = 0; i < nPoints; i++)
    {
        x = x0 + i*dx;
        u[i+padding] = u_0(x);
    }
}

double update(int index)
{
    return (u[index] - a*(dt/dx)*(u[index+1] - u[index]));
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
    double t = tStart;
    double local_dt;
    set_initial_value();
    while (t < tStop)
    {
        local_dt = min(dt, tStop-tStart); // avoid overshooting the final time

        //Periodic boundary conditions
        for (int i = 0; i < padding; i++)
        {
            u[i] = u[nPoints+i];
            u[nPoints+padding+i] = u[padding+i];
        }
        //Update the data
        //Make sure the limits on this loop correspond only to the
        // true domain
        for (int i = padding; i < nPoints+padding; i++)
        {

            uPlus1[i] = update(i);
        }
        // Now replace u with the updated data for the next time step
        u = uPlus1;
        t = t + local_dt;
    }

    print_result(u);
    output_file.close();
    return 0;
}
