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

int nPoints = 10000;
int padding = 1;
double x0 = 0;
double x1 = 1;
double tStart = 0.0;
double tStop;
double a = 1;
double dx = (x1 - x0)/nPoints;
double C = 0.8;
double gas_coef = 1.8;
int num_variables = 3;
double v_l, v_r, rho_l, rho_r, p_l, p_r;
//define a vector to store discretised data
vector<vector<double>> u; //u[0] and u[nPoints+1] are padding
vector<vector<double>> f; //f[i] holds the flux from u[i] to u[i+1]





void boundary_condition(vector<vector<double>> &values, int pad)
{
    for (int i = 0; i < pad; i++)
    {
        for (int j = 0; j < num_variables; j++)
        {
            values[i][j] = values[pad][j];
            values[pad + nPoints + i][j] = values[pad + nPoints][j];
        }
    }
}

//vector<double> (*u_0)(vector<double>);

vector<double> (flux)(const vector<double> &);

vector<double> (*nummerical_flux)(const vector<double>&, const vector<double>&, double, double);


vector<vector<double>> PDE_Step(const vector<vector<double>> &current, const vector<vector<double>> &flux, double time_step, int pad)
{
    vector<vector<double>> uPlus1 = current;
    for (int index = pad; index < (nPoints + pad); index++)
    {
        for (int j = 0; j < num_variables; j++)
        {
            uPlus1[index][j] = current[index][j] - ((time_step/dx)*(flux[index][j] - flux[index-1][j]));
        }
    }
    boundary_condition(uPlus1, pad);
    return uPlus1;
}

vector<double> Lax_Friedrichs_flux(const vector<double> & velocity, const vector<double> & velocity_next, double d_x, double d_t)
{
    vector<double> result;
    result.resize(num_variables);
    vector<double> f_left = flux(velocity);
    vector<double> f_right = flux(velocity_next);
    for (int i = 0; i < num_variables; i++)
    {
        result[i] = (d_x / d_t)*(velocity[i] - velocity_next[i]);
        result[i] += f_left[i];
        result[i] += f_right[i];
        result[i] = result[i] / 2;
    }
    return result;
}

vector<double> Richtmyer_flux(const vector<double> & velocity, const vector<double> & velocity_next, double d_x, double d_t)
{
    vector<double> velocity_half_step;
    velocity_half_step.resize(num_variables);
    vector<double> f_left = flux(velocity);
    vector<double> f_right = flux(velocity_next);
    for (int i = 0; i < num_variables; i++)
    {
        velocity_half_step[i]= 0.5*(velocity[i] + velocity_next[i]);
        velocity_half_step[i] += (-0.5)*(d_t/d_x)*(f_right[i]-f_left[i]);
    }
    return flux(velocity_half_step);
}

vector<double> FORCE_flux(const vector<double> & velocity, const vector<double> & velocity_next, double d_x, double d_t)
{
    vector<double> result = Lax_Friedrichs_flux(velocity, velocity_next, d_x, d_t);
    vector<double> R = Richtmyer_flux(velocity, velocity_next, d_x, d_t);
    for (int i = 0; i < num_variables; i++)
    {
        result[i] = 0.5*(result[i] + R[i]);
    }
    return result;
}

void set_initial_value(vector<vector<double>> &vel, int pad)
{
    double rho_v_l = rho_l*v_l;
    double E_l = p_l/(gas_coef - 1) + (rho_l*v_l*v_l/2);
    double rho_v_r = rho_r*v_r;
    double E_r = p_r/(gas_coef - 1) + (rho_r*v_r*v_r/2);
    double x;
    double integral = 0;
    for (int i = 0; i <= nPoints/2 ; i ++)
    {
        vel[i+pad][0] = rho_l;
        vel[i+pad][1] = rho_v_l;
        vel[i+pad][2] = E_l;
        vel[nPoints-i+pad][0] = rho_r;
        vel[nPoints-i+pad][1] = rho_v_r;
        vel[nPoints-i+pad][2] = E_r;
    }
    if (nPoints % 2 == 1)
    {
        vel[(nPoints/2)+1+pad][0] = (rho_r + rho_l)/2.0;
        vel[(nPoints/2)+1+pad][1] = (rho_v_r + rho_v_l)/2.0;
        vel[(nPoints/2)+1+pad][2] = (E_r + E_l)/2.0;
    }
    boundary_condition(vel, pad);
}

void print_result(vector<vector<double>> data)
{
    for (int i = 0; i < nPoints; i++)
    {
        for (int j = 0; j < num_variables; j++)
        {
            std::cout << data[i+padding][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout <<std::endl << std::endl;
}

double pressure(const vector<double> & state)
{
    double rho_v = state[1];
    double rho = state[0];
    return ((gas_coef - 1)*(state[2] - (rho_v*rho_v)/(2*rho)));
}

double speed_of_sound(const vector<double> & state)
{
    double rho = state[0];
    double p = pressure(state);
    return (sqrt(gas_coef*p/rho));
}

double maximal_v(const vector<vector<double>>& values)
{
    int number = values.size();
    double result = 0;
    double v_c;
    for (int i = 0; i < number; i++)
    {
        v_c = values[i][1] / values[i][0]; // rho_v / rho
        v_c = abs(v_c) + speed_of_sound(values[i]);
        if (v_c > result)
        {
            result = v_c;
        }
    }
    return result;
}

vector<double> flux(const vector<double> & state)
{
    vector<double> result;
    result.resize(3);
    double rho = state[0];
    double rho_v = state[1];
    double E = state[2];
    double p = pressure(state);

    result[0] = rho_v;
    result[1] = p + rho_v*rho_v/rho;
    result[2]= (E + p)*(rho_v/rho);
    
    return result;
}

void evaluate_flux(const vector<vector<double>> &vel, vector<vector<double>> &flu, double space_step, double time_step)
{
    for (int i = 0; i < (nPoints + (2*padding) - 1); i++)
    {
        //std::cout << i << std::endl;
        flu[i] = nummerical_flux(vel[i], vel[i+1], space_step, time_step);
    }
}

int main(int argc, char *argv[])
{
    nummerical_flux = FORCE_flux;
    if (argc != 2)
    {
        std::cout << "Remember to specify the test case you are running\n";
        exit(1);
    }
    // the first argument is the test case
    if (argv[1][4] == '1') // Test1
    {
        rho_l = 1.0;
        v_l = 0.0;
        p_l = 1.0;
        rho_r = 0.125;
        v_r = 0.0;
        p_r = 0.1;
        tStop = 0.25;
    }
    else
    {
        if (argv[1][4] == '2') // Test2
        {
            rho_l = 1.0;
            v_l = -2.0;
            p_l = 0.4;
            rho_r = 1.0;
            v_r = 2.0;
            p_r = 0.4;
            tStop = 0.15;
        }
        else
        {
            if (argv[1][4] == '3') // Test3
            {
                rho_l = 1.0;
                v_l = 0.0;
                p_l = 1000.0;
                rho_r = 1.0;
                v_r = 0.0;
                p_r = 0.01;
                tStop = 0.012;
            }
            else
            {
                if (argv[1][4] == '4') // Test4
                {
                    rho_l = 1.0;
                    v_l = 0.0;
                    p_l = 0.01;
                    rho_r = 1.0;
                    v_r = 0.0;
                    p_r = 100.0;
                    tStop = 0.035;
                }
                else
                {
                    if (argv[1][4] == '5') // Test5
                    {
                        rho_l =  5.99924;
                        v_l =  19.5975;
                        p_l =  460.894;
                        rho_r = 5.99242;
                        v_r = -6.19633;
                        p_r =  46.095;
                        tStop = 0.035;
                    }
                    else
                    {
                        std::cout << "invalid test case\n";
                        exit(1);
                    }
                }
            }
        }
    }
    

    dx = (x1 - x0) / nPoints;
    u.resize(nPoints + 2*padding);
    f.resize(nPoints + 2*padding - 1);

    for (int i = 0; i < (nPoints+(2*padding)-1); i++)
    {
        u[i].resize(num_variables);
        f[i].resize(num_variables);
    }
    u[nPoints+(2*padding)-1].resize(3);

    set_initial_value(u, padding);
    print_result(u);
    std::cout << "\t\n";

    double t = tStart;
    double dt;
    double a_max;
    int step = 1;
    while (t < tStop)
    {
        a_max = maximal_v(u);
        dt = C * dx / a_max;
        // avoid overshooting the final time
        if ((t + dt) > tStop)
        {
            dt = tStop - t;
        }
        // if (step % 100 == 1)
        // {
        //     std::cout << "Running "<< step << "th step, t = " << t << " dt = " << dt << std::endl;
        // }
        
        evaluate_flux(u, f, dx, dt);

        //std::cout << "Flux evaluated\n";

        u = PDE_Step(u, f, dt, padding);

        //std::cout << "PDE Step done\n";

        t = t + dt;
        step += 1;
    }

    print_result(u);
    return 0;
}
