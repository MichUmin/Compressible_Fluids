#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <functional>
#include <string>
#include "tables.cpp"

using std::min;
using std::cout;
using std::cin;
using std::vector;
using std::string;

const double pi = 3.14159265358979311600;

int nPoints_x = 101;
int nPoints_y = 101;
int padding = 1;
double x_0 = 0.0;
double x_1 = 1.0;
double y_0 = 0.0;
double y_1 = 1.0;
double tStart = 0.0;
double tStop;
double dx = (x_1 - x_0)/nPoints_x;
double dy = (y_1 - y_0)/nPoints_y;
double C = 0.8;
double gas_coef = 1.4;
int num_dimensions = 2;
int num_variables = num_dimensions + 2;
int step = 0;
double v_l, v_r, rho_l, rho_r, p_l, p_r;
bool use_halfstep = true;
//define a vector to store discretised data
table<vector<double>> u(nPoints_x + 2*padding, nPoints_y + 2*padding); // computed values
table<vector<double>> f(nPoints_x + 2*padding, nPoints_y + 2*padding); // f[i, j] holds the flux from u[i, j] to u[i+1, j] or to u[i, j+1]
table<vector<double>> u_plus(nPoints_x + 2*padding, nPoints_y + 2*padding);
table<vector<double>> u_minus(nPoints_x + 2*padding, nPoints_y + 2*padding);



void (*boundary_condition)(table<vector<double>> &values, int pad);

void transmissive_boundary_condition(table<vector<double>> &values, int pad)
{
    for (int i = 0; i < pad; i++)
    {
        for (int j = pad; j < nPoints_x + pad; j++)
        {
            values(j, i) = values(j, pad);
            values(j, pad + nPoints_y + i) = values(j, pad + nPoints_y - 1);
        }
        for (int j = pad; j < nPoints_y + pad; j++)
        {
            values(i, j) = values(padding,j);
            values(pad + nPoints_x + i, j) = values(padding + nPoints_x - 1, j);
        }
    }
}

vector<double> (flux)(int ,const vector<double> &);

double vec_len(const vector<double> & v)
{
    int num = v.size();
    double result = 0;
    for (int i = 0; i < num; i++)
    {
        result += (v[i]*v[i]);
    }
    return (sqrt(result));
}

vector<double> (*nummerical_flux)(int direction, const vector<double>&, const vector<double>&, double, double);

double square(double num)
{
    return num*num;
}

table<vector<double>> (*PDE_Step)(table<vector<double>> & current, table<vector<vector<double>>> & flux, double time_step, int pad);
/*
table<vector<double>> directional_PDE_Step(int direction, const table<vector<double>> &current, const table<vector<double>> &flux, double space_step, double time_step)
{
    int next_x = 0;
    int next_y = 0;
    switch (direction)
    {
    case 1:
        next_x = 1;
        break;
    case 2:
        next_y = 1;
        break;
    default:
        std::cout << "Invalid direction\n";
        exit(3);
    }

    table<vector<double>> uPlus1(nPoints_x + 2*padding, nPoints_y + 2*padding);
    uPlus1 = current;
    for (int index_x = padding; index_x < (nPoints_x + padding); index_x++)
    {
        for (int index_y = padding; index_y < (nPoints_y + padding); index_y++)
        {
            for (int j = 0; j < num_variables; j++)
            {
                uPlus1(index_x, index_y)[j] = current(index_x, index_y)[j] - ((time_step / space_step)*(flux(index_x, index_y)[j] - flux(index_x - next_x, index_y - next_y)[j]));
            }
        }
    }
    boundary_condition(uPlus1, pad);
    return uPlus1;
}
*/
vector<double> Lax_Friedrichs_flux(int direction, const vector<double> & velocity, const vector<double> & velocity_next, double d_x, double d_t)
{
    vector<double> result;
    result.resize(num_variables);
    vector<double> f_left = flux(direction, velocity);
    vector<double> f_right = flux(direction, velocity_next);
    for (int i = 0; i < num_variables; i++)
    {
        result[i] = (d_x / d_t)*(velocity[i] - velocity_next[i]);
        result[i] += f_left[i];
        result[i] += f_right[i];
        result[i] = result[i] / 2.0;
    }
    return result;
}

vector<double> Richtmyer_flux(int direction, const vector<double> & velocity, const vector<double> & velocity_next, double d_x, double d_t)
{
    vector<double> velocity_half_step;
    velocity_half_step.resize(num_variables);
    vector<double> f_left = flux(direction, velocity);
    vector<double> f_right = flux(direction, velocity_next);
    for (int i = 0; i < num_variables; i++)
    {
        velocity_half_step[i]= 0.5*(velocity[i] + velocity_next[i]);
        velocity_half_step[i] += (-0.5)*(d_t/d_x)*(f_right[i]-f_left[i]);
    }
    //std::cout << velocity_half_step[0] << " " << velocity_half_step[1] << " " << velocity_half_step[2] << "\n";
    return flux(direction, velocity_half_step);
}

vector<double> FORCE_flux(int direction, const vector<double> & velocity, const vector<double> & velocity_next, double d_x, double d_t)
{
    vector<double> result = Lax_Friedrichs_flux(direction, velocity, velocity_next, d_x, d_t);
    vector<double> R = Richtmyer_flux(direction, velocity, velocity_next, d_x, d_t);
    for (int i = 0; i < num_variables; i++)
    {
        result[i] = 0.5*(result[i] + R[i]);
    }
    return result;
}

void (*set_initial_value)(table<vector<double>> &);

void set_X(table<vector<double>> &vel)
{
    double rho_v_l = rho_l*v_l;
    double E_l = p_l/(gas_coef - 1.0) + (rho_l*v_l*v_l/2.0);
    double rho_v_r = rho_r*v_r;
    double E_r = p_r/(gas_coef - 1.0) + (rho_r*v_r*v_r/2.0);
    double x;
    double integral = 0;
    for (int j = 0; j < nPoints_y + 2*padding; j++)
    {
        for (int i = 0; i <= nPoints_x/2; i++)
        {
            vel(i+padding, j)[0] = rho_l;
            vel(i+padding, j)[1] = rho_v_l;
            vel(i+padding, j)[2] = 0;
            vel(i+padding, j)[num_variables-1] = E_l;
            vel(nPoints_x-i+padding, j)[0] = rho_r;
            vel(nPoints_x-i+padding, j)[1] = rho_v_r;
            vel(nPoints_x-i+padding, j)[2] = 0.0;
            vel(nPoints_x-i+padding, j)[num_variables-1] = E_r;
        }
        if (nPoints_x % 2 == 1)
        {
            vel((nPoints_x/2)+1+padding, j)[0] = (rho_r + rho_l)/2.0;
            vel((nPoints_x/2)+1+padding, j)[1] = (rho_v_r + rho_v_l)/2.0;
            vel((nPoints_x/2)+1+padding, j)[2] = 0.0;
            vel((nPoints_x/2)+1+padding, j)[num_variables-1] = (E_r + E_l)/2.0;
        }
    }
    //std::cout << "main part done\n";
    boundary_condition(vel, padding);
}

void set_Y(table<vector<double>> &vel)
{
    double rho_v_l = rho_l*v_l;
    double E_l = p_l/(gas_coef - 1.0) + (rho_l*v_l*v_l/2.0);
    double rho_v_r = rho_r*v_r;
    double E_r = p_r/(gas_coef - 1.0) + (rho_r*v_r*v_r/2.0);

    for (int j = 0; j < nPoints_x + 2*padding; j++)
    {
        for (int i = 0; i <= nPoints_y/2; i++)
        {
            vel(j, i+padding)[0] = rho_l;
            vel(j, i+padding)[1] = 0.0;
            vel(j, i+padding)[2] = rho_v_l;
            vel(j, i+padding)[num_variables-1] = E_l;
            vel(j, nPoints_y-i+padding)[0] = rho_r;
            vel(j, nPoints_y-i+padding)[1] = 0.0;
            vel(j, nPoints_y-i+padding)[2] = rho_v_r;
            vel(j, nPoints_y-i+padding)[num_variables-1] = E_r;
        }
        if (nPoints_y % 2 == 1)
        {
            vel(j, (nPoints_y/2)+1+padding)[0] = (rho_r + rho_l)/2.0;
            vel(j, (nPoints_y/2)+1+padding)[1] = 0.0;
            vel(j, (nPoints_y/2)+1+padding)[2] = (rho_v_r + rho_v_l)/2.0;
            vel(j, (nPoints_y/2)+1+padding)[num_variables-1] = (E_r + E_l)/2.0;
        }
    }
    boundary_condition(vel, padding);
}

void set_diagonal(table<vector<double>> &vel)
{
    double rho_vx_l = rho_l * v_l * nPoints_x / sqrt(nPoints_x*nPoints_x + nPoints_y*nPoints_y);
    double rho_vy_l = rho_l * v_l * nPoints_y / sqrt(nPoints_x*nPoints_x + nPoints_y*nPoints_y);
    double E_l = p_l/(gas_coef - 1.0) + (rho_l*v_l*v_l/2.0);
    double rho_vx_r = rho_r * v_r * nPoints_x / sqrt(nPoints_x*nPoints_x + nPoints_y*nPoints_y);
    double rho_vy_r = rho_r * v_r * nPoints_y / sqrt(nPoints_x*nPoints_x + nPoints_y*nPoints_y);
    double E_r = p_r/(gas_coef - 1.0) + (rho_r*v_r*v_r/2.0);

    int k;
    for (int i = 0; i < nPoints_x; i++)
    {
        for (int j = 0; j <= nPoints_y; j++)
        {
            if ((i + j) < ((nPoints_x + nPoints_y)/2))
            {
                vel(i+padding, j+padding)[0] = rho_l;
                vel(i+padding, j+padding)[1] = rho_vx_l;
                vel(i+padding, j+padding)[2] = rho_vy_l;
                vel(i+padding, j+padding)[num_variables-1] = E_l;
            }
            else
            {
                vel(i+padding, j+padding)[0] = rho_r;
                vel(i+padding, j+padding)[1] = rho_vx_r;
                vel(i+padding, j+padding)[2] = rho_vy_r;
                vel(i+padding, j+padding)[num_variables-1] = E_r;
            }
        }
        if ((nPoints_x+nPoints_y) % 2 == 1)
        {
            k = (nPoints_x+nPoints_y) / 2 - i;
            vel(i+padding, k+padding)[0] = (rho_r + rho_l)/2.0;
            vel(i+padding, k+padding)[1] = (rho_vx_r + rho_vx_l)/2.0;
            vel(i+padding, k+padding)[2] = (rho_vy_r + rho_vy_l)/2.0;
            vel(i+padding, k+padding)[num_variables-1] = (E_r + E_l)/2.0;
        }
    }
    boundary_condition(vel, padding);
}

void set_cylindrical(table<vector<double>> &vel)
{
    double E_l = p_l/(gas_coef - 1.0) + (rho_l*v_l*v_l/2.0);
    double E_r = p_r/(gas_coef - 1.0) + (rho_r*v_r*v_r/2.0);

    for (int i = 0; i < nPoints_x; i++)
    {
        for (int j = 0; j <= nPoints_y; j++)
        {
            vel(i+padding, j+padding)[1] = 0;
            vel(i+padding, j+padding)[2] = 0;
            if (((2*i - nPoints_x)*(2*i - nPoints_x)*nPoints_y*nPoints_y + (2*j - nPoints_y)*(2*j - nPoints_y)*nPoints_x*nPoints_x) < 0.16*nPoints_x*nPoints_x*nPoints_y*nPoints_y)
            {
                vel(i+padding, j+padding)[0] = rho_l;
                vel(i+padding, j+padding)[num_variables-1] = E_l;
            }
            else
            {
                vel(i+padding, j+padding)[0] = rho_r;
                vel(i+padding, j+padding)[num_variables-1] = E_r;
            }
        }
    }
    boundary_condition(vel, padding);
}

void (*print_result)(table<vector<double>> & data);

void print_horizontal(table<vector<double>> & data)
{
    int j = (nPoints_y / 2) + padding;
    for (int i = padding; i < nPoints_x+padding; i++)
    {
        std::cout << data(i, j) << std::endl;
    }
    std::cout <<std::endl;
}

void print_vertical(table<vector<double>> & data)
{
    int i = (nPoints_x / 2) + padding;
    for (int j = padding; j < nPoints_y+padding; j++)
    {
        std::cout << data(i, j) << std::endl;
    }
    std::cout <<std::endl;
}

void print_diagonal(table<vector<double>> & data)
{
    for (int i = padding; i < nPoints_x+padding; i++)
    {
        std::cout << data(i, i) << std::endl;
    }
    std::cout <<std::endl;
}

double pressure(const vector<double> & state)
{
    double rho = state[0];
    double E = state[num_variables-1];
    vector<double> v;
    v.resize(num_dimensions);
    double p;
    if (rho == 0.0)
    {
        /*
        std::cout << "Zero density when calculating the pressure\n";
        exit(2);
        */
       p = ((gas_coef - 1)*E);
    }
    else
    {
        for (int j = 0; j < num_dimensions; j++)
        {
            v[j] = state[1+j]/rho;
        }
        p = ((gas_coef - 1)*(E - (square(vec_len(v))*rho*0.5)));
    }
    /*
    if (p < 0.0)
    {
        std::cout << "Negative pressure\n";
        exit(3);
    }
    */
    return p;
}

double speed_of_sound(const vector<double> & state)
{
    double rho = state[0];
    double p = pressure(state);
    double cs;
    if (p < 0.0)
    {
        /*
        std::cout << "Negative pressure when calculating the speed of sound\n";
        exit(3);
        */
        p = (-1.0)*p;
    }
    cs = sqrt(gas_coef*p/rho);
    /*
    if (cs > 100)
    {
        std::cout << "Extreme cs: p = " << p << ", rho = " << rho << std::endl;
    }
    */
    return (cs);
}

double maximal_v(table<vector<double>>& values)
{
    double result = 0;
    double rho, v_c, cs;
    vector<double> state, v;
    v.resize(num_dimensions);
    for (int i = padding; i < nPoints_x+padding; i++)
    {
        for (int j = padding; j < nPoints_y+padding; j++)
        {
            state = values(i,j);
            rho = state[0];
            if (rho != 0.0)
            {
                for (int k = 0; k < num_dimensions; k++)
                {
                    v[k] = state[1+k]/rho;
                }
                cs = speed_of_sound(values(i, j));
                v_c = vec_len(v) + cs;
                if (v_c > result)
                {
                    result = v_c;
                }
            }
            else
            {
                std::cout << "rho equal 0 at i = " << i << ", j = " << j << std::endl;
                exit(3);
            }
        }
    }
    return result;
}

vector<double> flux(int direction, const vector<double> & state)
{
    // direction: x - 1, y - 2, z - 3
    vector<double> result;
    result.resize(num_variables);
    double rho = state[0];
    double rho_v = state[direction];
    double E = state[num_variables-1];
    double p = pressure(state);
    double v;

    result[0] = rho_v;
    if (rho == 0.0)
    {
        /*
        std::cout << "Zero density\n";
        exit(2);
        */
       for (int i = 1; i <= num_dimensions; i++)
       {
            result[i] = 0;
       }
       result[direction] = p;
       result[num_variables-1] = 0;
    }
    else
    {
        v = rho_v / rho;
        for (int i=1; i <= num_dimensions; i++)
        {
            result[i] = state[i]*v;
        }
        result[direction] += p;
        result[num_variables-1]= (E + p)*v;
    }
    return result;
}

double (*limiter)(double r);

double minbee(double r)
{
    if (r <= 0.0)
    {
        return 0.0;
    }
    else
    {
        if ((r <= 1.0)&&(r > 0.0))
        {
            return r;
        }
        else
        {
            return (2.0 / (1.0 + r));
        }
    }
}

double zero_limiter(double r)
{
    return 0.0;
}

double superbee(double r)
{
    if (r <= 0.0)
    {
        return 0.0;
    }
    else
    {
        if (r <= 0.5)
        {
            return 2.0*r;
        }
        else
        {
            if (r <= 1)
            {
                return 1.0;
            }
            else
            {
                return (2.0 / (1.0 + r));
            }
        }
    }
}

double Van_Leer(double r)
{
    if (r <= 0.0)
    {
        return 0.0;
    }
    else
    {
        if (r <= 1.0)
        {
            return ((2.0 * r) / (1.0 + r));
        }
        else
        {
            return (2.0 / (1.0 + r));
        }
    }
}

void SLICK_Step(int direction, table<vector<double>> &current, table<vector<double>> &fluxes, double space_step, double time_step)
{
    int next_x = 0;
    int next_y = 0;
    switch (direction)
    {
    case 1:
        next_x = 1;
        break;
    case 2:
        next_y = 1;
        break;
    default:
        std::cout << "Invalid direction\n";
        exit(3);
    }
    u_minus = current;
    u_plus = current;
    double r;
    double lim;
    double grad;
    for (int index_x = next_x; index_x < (nPoints_x + 2*padding - next_x); index_x++)
    {
        for (int index_y = next_y; index_y < (nPoints_y + 2*padding - next_y); index_y++)
        {
            if (current(index_x, index_y)[num_variables-1] != current(index_x + next_x, index_y + next_y)[num_variables-1])
            {
                r = ((current(index_x, index_y)[num_variables-1] - current(index_x - next_x, index_y - next_y)[num_variables-1]) / (current(index_x + next_x, index_y + next_y)[num_variables-1] - current(index_x, index_y)[num_variables-1]));
                lim = limiter(r);
            }
            else
            {
                lim = 0.0;
            }
            for (int j = 0; j < num_variables; j++)
            {
                grad = (current(index_x+next_x, index_y+next_y)[j] - current(index_x-next_x, index_y-next_y)[j])/2.0;
                u_minus(index_x, index_y)[j] -= grad*lim/2.0;
                u_plus(index_x, index_y)[j] += grad*lim/2.0;
            }
        }
    }
    boundary_condition(u_minus, padding);
    boundary_condition(u_plus, padding);
    if (use_halfstep)
    {
        vector<double> flux_plus, flux_minus;
        double change;
        for (int index_x = padding; index_x < nPoints_x + padding; index_x++)
        {
            for (int index_y = padding; index_y < nPoints_y + padding; index_y++)
            {
                flux_minus = flux(direction, u_minus(index_x, index_y));
                flux_plus = flux(direction, u_plus(index_x, index_y));
                for (int j = 0; j < num_variables; j++)
                {
                    change = (flux_plus[j] - flux_minus[j]) * time_step / space_step / 2.0;
                    u_minus(index_x, index_y)[j] -= change;
                    u_plus(index_x, index_y)[j] -= change;
                }
            }
        }
        boundary_condition(u_minus, padding);
        boundary_condition(u_plus, padding);
    }

    for (int i = 0; i < (nPoints_x + (2*padding) - next_x); i++)
    {
        for (int j = 0; j < (nPoints_y + (2*padding) - next_y); j++)
        {
            fluxes(i, j) = nummerical_flux(direction, u_plus(i, j), u_minus(i + next_x, j + next_y), space_step, time_step);
        }
    }
    
    for (int index_x = padding; index_x < (nPoints_x + padding); index_x++)
    {
        for (int index_y = padding; index_y < (nPoints_y + padding); index_y++)
        {
            for (int j = 0; j < num_variables; j++)
            {
                current(index_x, index_y)[j] = current(index_x, index_y)[j] - ((time_step / space_step)*(fluxes(index_x, index_y)[j] - fluxes(index_x - next_x, index_y - next_y)[j]));
            }
        }
    }
    boundary_condition(current, padding);
}

int main(int argc, char *argv[])
{
    nummerical_flux = FORCE_flux;
    limiter = Van_Leer;
    boundary_condition = transmissive_boundary_condition;

    if (argc < 2)
    {
        std::cout << "Remember to specify the test case you are running\n";
        exit(1);
    }
    // the first argument is the test case
    if (argv[1][0] == 'T')
    {
        if (argv[1][4] == 'C') // ToroCylindrical
        {
            x_1 = 2.0;
            y_1 = 2.0;
            rho_l = 1.0;
            v_l = 0.0;
            p_l = 1.0;
            rho_r = 0.125;
            v_r = 0.0;
            p_r = 0.1;
            tStop = 0.25;
            set_initial_value = set_cylindrical;
            print_result = print_horizontal;
        }
        else
        {
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
                    std::cout << "invalid test case\n";
                    exit(1);
                }
            }
            if (argv[2][0] == 'X') // 
            {
                set_initial_value = set_X;
                print_result = print_horizontal;
            }
            else
            {
                if (argv[2][0] == 'Y')
                {
                    set_initial_value = set_Y;
                    print_result = print_vertical;
                }
                else
                {
                    if (argv[2][0] == 'D') // Diagonal
                    {
                        set_initial_value = set_diagonal;
                        print_result = print_diagonal;
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
    else
    {
        if (argv[1][0] == 'F') // Full2D
        {
            std::cout << "Remember to implement this test case\n";
            exit(0);

        }
        else
        {
            std::cout << "invalid test case\n";
            exit(1);
        }
    }

    dx = (x_1 - x_0) / nPoints_x;
    dy = (y_1 - y_0) / nPoints_y;

    for (int i = 0; i < (nPoints_x+(2*padding)); i++)
    {
        for (int j = 0; j < (nPoints_y+(2*padding)); j++)
        {
            u(i, j).resize(num_variables);
            f(i, j).resize(num_variables);
            u_minus(i, j).resize(num_variables);
            u_plus(i, j).resize(num_variables);
        }
    }
    //std::cout << "resized\n";

    set_initial_value(u);

    //std::cout << "initial falue set\n";
    print_result(u);
    std::cout << "\t\n";

    double t = tStart;
    double dt;
    double a_max;
    step = 0;
    //tStop = 0.001;

    while (t < tStop)
    {
        //std::cout << "step: " << step << std::endl; 
        a_max = maximal_v(u); 
        dt = C * min(dx, dy) / a_max;
        // avoid overshooting the final time
        if ((t + 2*dt) > tStop)
        {
            dt = (tStop - t) / 2.0;
        }
        //std::cout << "dt = " << dt << " t = " << t << " a_max = " << a_max << std::endl;      

        SLICK_Step(1, u, f, dx, dt);
        SLICK_Step(2, u, f, dy, dt);
        SLICK_Step(2, u, f, dy, dt);
        SLICK_Step(1, u, f, dx, dt);

        t += (2.0 * dt);
        step += 2;

    }

    print_result(u);
    //debug_result(u);
    return 0;
}
