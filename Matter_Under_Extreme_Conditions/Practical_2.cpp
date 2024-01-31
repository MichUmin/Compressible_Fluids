#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <functional>
#include  <string>
#include <cassert>

using std::min;
using std::cout;
using std::cin;
using std::vector;
using std::string;

const double pi = 3.14159265358979311600;

int nPoints = 1000;
const int padding = 2;
double x0 = 0;
double x1 = 1;
double tStart = 0.0;
double tStop;
double dx = (x1 - x0)/nPoints;
double CFL = 0.8;
double gas_coef = 1.4; 
const int num_dimensions = 1;
const int num_variables = num_dimensions + 2;
int step = 0;
double v_l, v_r, rho_l, rho_r, p_l, p_r, rho_v_l, rho_v_r, epsilon_l, epsilon_r, E_l, E_r;
const bool use_halfstep = true;
//const bool use_halfstep = false;

vector<vector<double>> u(nPoints + 2*padding); //u[0] u[1] u[nPoints+2] u[nPoints+3] are padding
//vector<vector<double>> u_copy(nPoints + 2*padding);
vector<vector<double>> f(nPoints + 2*padding - 1); //f[i] holds the flux from u[i] to u[i+1]
vector<vector<double>> u_minus(nPoints + 2*padding);
vector<vector<double>> u_plus(nPoints + 2*padding);

vector<double> (flux)(const vector<double> &);

vector<double> (*source)(const vector<double> &, double);

void set_size(vector<vector<double>> & my_vector)
{
    int length = my_vector.size();
    for (int it = 0; it != length; it++)
    {
        my_vector[it].resize(num_variables);
    }
}

double (pressure)(const vector<double> &);

double (speed_of_sound)(const vector<double> &);

void boundary_condition(vector<vector<double>> &values, int pad)
{
    for (int i = 0; i < pad; i++)
    {
        for (int j = 0; j < num_variables; j++)
        {
            values[i][j] = values[((2*pad)-i)-1][j];
            values[pad + nPoints + i][j] = values[pad + nPoints - 1][j];
        }
        values[i][1] = (-1.0)*values[((2*pad)-i)-1][1]; // reflective boundary condition in r direction at r=0
    }
}

vector<double> (*nummerical_flux)(const vector<double>&, const vector<double>&, double, double);

vector<double> Lax_Friedrichs_flux(const vector<double> & velocity, const vector<double> & velocity_next, double space_step, double time_step)
{
    vector<double> result;
    result.resize(num_variables);
    vector<double> f_left = flux(velocity);
    vector<double> f_right = flux(velocity_next);
    for (int i = 0; i < num_variables; i++)
    {
        result[i] = (space_step / time_step)*(velocity[i] - velocity_next[i]);
        result[i] += f_left[i];
        result[i] += f_right[i];
        result[i] = result[i] / 2.0;
    }
    return result;
}

vector<double> Richtmyer_flux(const vector<double> & velocity, const vector<double> & velocity_next, double space_step, double time_step)
{
    vector<double> velocity_half_step;
    velocity_half_step.resize(num_variables);
    vector<double> f_left = flux(velocity);
    vector<double> f_right = flux(velocity_next);
    for (int i = 0; i < num_variables; i++)
    {
        velocity_half_step[i]= 0.5*(velocity[i] + velocity_next[i]);
        velocity_half_step[i] += (-0.5)*(time_step/space_step)*(f_right[i]-f_left[i]);
    }
    //std::cout << velocity_half_step[0] << " " << velocity_half_step[1] << " " << velocity_half_step[2] << "\n";
    return flux(velocity_half_step);
}

vector<double> FORCE_flux(const vector<double> & velocity, const vector<double> & velocity_next, double space_step, double time_step)
{
    vector<double> result = Lax_Friedrichs_flux(velocity, velocity_next, space_step, time_step);
    vector<double> R = Richtmyer_flux(velocity, velocity_next, space_step, time_step);
    for (int i = 0; i < num_variables; i++)
    {
        result[i] = 0.5*(result[i] + R[i]);
    }
    return result;
}

void set_initial_value(vector<vector<double>> &vel, int pad)
{
    double rho_v_l = rho_l*v_l;
    double E_l = p_l/(gas_coef - 1.0) + (rho_l*v_l*v_l/2.0);
    double rho_v_r = rho_r*v_r;
    double E_r = p_r/(gas_coef - 1.0) + (rho_r*v_r*v_r/2.0);

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

void print_result(const vector<vector<double>> & data)
{
    for (int i = 0; i < nPoints; i++)
    {
        //std::cout << i << " ";
        for (int j = 0; j < num_variables; j++)
        {
            //std::cout << j << " ";
            std::cout << data[i+padding][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout <<std::endl << std::endl;
}

void debug_result(const vector<vector<double>> & data)
{
    for (auto i = data.begin(); i != data.end(); i++)
    {
        for (int j = 0; j < num_variables; j++)
        {
            std::cout << (*i)[j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout <<std::endl << std::endl;
}

double internal_energy(const vector<double> & state)
{
    double rho_v = state[1];
    double rho = state[0];
    double v;
    double E = state[2];
    double epsilon;
    if (rho == 0.0)
    {
        cout << "Zero density at step: " << step  <<"\n";
        exit(3);
    }
    else
    {
        v = rho_v / rho;
        epsilon = (E - (0.5*rho_v*v))/rho;
    }
    return epsilon;
}

double pressure(const vector<double> & state)
{
    double rho_v = state[1];
    double rho = state[0];
    double E = state[2];
    double p;
    if (rho <= 0.0001)
    {
        /*
        std::cout << "Zero density when calculating the pressure\n";
        std::cout << rho << " " << rho_v << "\n";
        exit(2);
        */
       p = ((gas_coef - 1)*E);
    }
    else
    {
        p = ((gas_coef - 1)*(E - (rho_v*rho_v)/(2*rho)));
    }
    
    if (p < 0.0)
    {    
        cout << "Negative pressure\n";  
        //debug_result(u);
        //debug_result(u_left);
        //debug_result(u_right);
        
        //cout << "p = " << p << " E = " << E << " rho_v = " << rho_v << " rho = " << rho << " step = " << step << "\n";
        //exit(3);
    }
    
    return p;
}

double speed_of_sound(const vector<double> & state)
{
    double rho = state[0];
    if (rho == 0)
    {
        std::cout << "Zero density when calculating the speed of sound\n";
        std::cout << rho << "\n";
        exit(2);
    }
    double p = pressure(state);
    if (p < 0.0)
    {
        /*
        debug_result(u);
        std::cout << "Negative pressure when calculating the speed of sound\n";
        std::cout << p << "\n";
        exit(3);
        */
       return 0.0;
    }
    else
    {
        return (sqrt(gas_coef*p/rho));
    }
}


double maximal_v(const vector<vector<double>> & values)
{
    double result = 0.0;
    double v_c;
    for (int index = padding; index < nPoints + padding; index++)
    {
        v_c = speed_of_sound(values[index]);
        v_c += abs(values[index][1] / values[index][0]); // rho_v / rho

        if (v_c > result)
        {
            result = v_c;
        }
    }
    return result;
}

vector<double> flux(const vector<double> & state)
{
    vector<double> result(num_variables);
    double rho = state[0];
    double rho_v = state[1];
    double E = state[2];
    double p = pressure(state);
    double v;

    result[0] = rho_v;
    if (rho == 0.0)
    {
        /*
        std::cout << "Zero density\n";
        std::cout << rho << " " << rho_v << " " << E << "\n";
        exit(2);
        */
       result[1] = p;
       result[2] = 0;
    }
    else
    {
        v = rho_v / rho;
        result[1] = p + rho_v*v;
        result[2]= (E + p)*v;
    }
    
    return result;
}

vector<double> source_cylindrical(const vector<double> & state, double r)
{
    vector<double> result(num_variables);
    double rho = state[0];
    if (rho == 0.0)
    {
        cout << "Zero density when evaluating the source term\n";
        exit(1);
    }
    double rho_vr = state[1];
    double vr = rho_vr/rho;
    double rho_vz, vz;
    if (num_variables > 3)
    {
        rho_vz = state[2];
        vz = rho_vz/rho;
        result[2] = (-1.0)*rho_vr * vz / r;
    }
    double E = state[num_variables-1];
    double p = pressure(state);

    result[0] = (-1.0)*rho_vr / r;
    result[1] = (-1.0)*rho_vr * vr / r;

    result[num_variables - 1] = (-1.0)*(E + p) * vr / r;

    return result;
}

vector<double> source_spherical(const vector<double> & state, double r)
{
    vector<double> result(num_variables);
    double rho = state[0];
    if (rho == 0.0)
    {
        cout << "Zero density when evaluating the source term\n";
        exit(1);
    }
    double rho_vr = state[1];
    double vr = rho_vr/rho;
    double E = state[num_variables-1];
    double p = pressure(state);

    result[0] = (-2.0)*rho_vr / r;
    result[1] = (-2.0)*rho_vr * vr / r;
    result[num_variables - 1] = (-2.0) * (E + p) * vr / r;

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

void (*slope_reconstruction)(vector<vector<double>> const &, vector<vector<double>> &, vector<vector<double>> &);

void (slope_reconstruction_energy)(vector<vector<double>> const & middle, vector<vector<double>> & left, vector<vector<double>> & right)
{
    double r;
    double lim;
    double grad;
    for (int index_x = 1; index_x < (nPoints + 2*padding - 1); index_x++)
    {
        if (middle[index_x][num_variables-1] != middle[index_x + 1][num_variables-1])
        {
            r = ((middle[index_x][num_variables-1] - middle[index_x - 1][num_variables-1]) / (middle[index_x + 1][num_variables-1] - middle[index_x][num_variables-1]));
            lim = limiter(r);
        }
        else
        {
            lim = 0.0;
        }
        assert((lim<=1.0)&&(lim>=0.0));
        for (int j = 0; j < num_variables; j++)
        {
            grad = (middle[index_x+1][j] - middle[index_x-1][j])*0.5;
            left[index_x][j] = middle[index_x][j] - (grad*lim*0.5);
            right[index_x][j] = middle[index_x][j] + (grad*lim*0.5);
        }
    }
    boundary_condition(left, padding);
    boundary_condition(right, padding);
}

void (slope_reconstruction_min)(vector<vector<double>> const & middle, vector<vector<double>> & left, vector<vector<double>> & right)
{
    double r;
    double lim;
    double grad;
    for (int index_x = 1; index_x < (nPoints + 2*padding - 1); index_x++)
    {
        lim = 1.0;
        for (int j = 0; j < num_variables; j++)
        {
            if (middle[index_x][j] != middle[index_x + 1][j])
            {
                r = ((middle[index_x][j] - middle[index_x - 1][j]) / (middle[index_x + 1][j] - middle[index_x][j]));
                lim = min(lim, limiter(r));
            }
            else
            {
                if (middle[index_x - 1][j] != middle[index_x][j])
                {
                    lim = 0.0;
                }
                else
                {   
                    lim = min(1.0, lim);
                }
            }
        }
        assert((lim<=1.0)&&(lim>=0.0));
        for (int j = 0; j < num_variables; j++)
        {
            grad = (middle[index_x+1][j] - middle[index_x-1][j])*0.5;
            left[index_x][j] = middle[index_x][j] - (grad*lim*0.5);
            right[index_x][j] = middle[index_x][j] + (grad*lim*0.5);
        }
    }
    boundary_condition(left, padding);
    boundary_condition(right, padding);
}

void (*Conservative_Step)(vector<vector<double>> &, vector<vector<double>> &, const double, const double);

void SLICK_Step(vector<vector<double>> &current, vector<vector<double>> &fluxes, const double space_step, const double time_step)
{
    slope_reconstruction(current, u_minus, u_plus);
    
    if (use_halfstep)
    {
        vector<double> flux_plus, flux_minus;
        double change;
        for (int index_x = padding; index_x < nPoints+ padding; index_x++)
        {
            flux_minus = flux(u_minus[index_x]);
            flux_plus = flux(u_plus[index_x]);
            for (int j = 0; j < num_variables; j++)
            {
                change = (flux_plus[j] - flux_minus[j]) * (time_step / space_step) * 0.5;
                u_minus[index_x][j] -= change;
                u_plus[index_x][j] -= change;
            }
        }
        boundary_condition(u_minus, padding);
        boundary_condition(u_plus, padding);
    }
    //cout << "halfstep done\n";
    for (int i = 0; i < (nPoints + (2*padding) - 1); i++)
    {
        fluxes[i] = FORCE_flux(u_plus[i], u_minus[i + 1], space_step, time_step);
    }
    //cout << "fluxes evaluated\n";
    for (int index_x = padding; index_x < (nPoints + padding); index_x++)
    {
        for (int j = 0; j < num_variables; j++)
        {
            current[index_x][j] = current[index_x][j] - ((time_step / space_step)*(fluxes[index_x][j] - fluxes[index_x - 1][j]));
        }
    }

    boundary_condition(current, padding);
}

void (*ODE_Step)(vector<vector<double>> &, const double, const double);

void Euler_Step(vector<vector<double>> & states, const double space_step, const double time_step)
{
    vector<double> source_contribution(num_variables);
    double r = 0.5 * space_step;
    for (int position = padding; position < nPoints+padding ; position++)
    {
        source_contribution = source(states[position], r);
        for (int index = 0; index < num_variables; index++)
        {
            states[position][index] += time_step*source_contribution[index];
        }
        r += space_step;
    }
    boundary_condition(states, padding);
}

void Runge_Kutta_Step(vector<vector<double>> & states, const double space_step, const double time_step)
{
    vector<double> local_state(num_variables);
    vector<double> help(num_variables);
    vector<double> K1(num_variables);
    vector<double> K2(num_variables);
    double r = 0.5*space_step;
    for (int position = padding; position < nPoints + padding ; position++)
    {
        K1 = source(states[position], r);
        for (int index = 0; index < num_variables; index++)
        {
            K1[index] *= time_step;
            help[index] = states[position][index] + K1[index];
        }
        K2 = source(help, r);
        for (int index = 0; index < num_variables; index++)
        {
            K2[index] *= time_step;
            states[position][index] += (0.5*(K1[index] + K2[index]));
        }
        r += space_step;
    }
    boundary_condition(states, padding);
}


int main(int argc, char *argv[])
{
    nummerical_flux = FORCE_flux;
    limiter = Van_Leer;
    //limiter = superbee;
    //limiter = zero_limiter;
    Conservative_Step = SLICK_Step;
    //ODE_Step = Euler_Step;
    ODE_Step = Runge_Kutta_Step;
    //slope_reconstruction = slope_reconstruction_min;
    slope_reconstruction = slope_reconstruction_energy;

    if (argc != 3)
    {
        std::cout << "Remember to specify the test case you are running\n";
        exit(1);
    }
    //CFL = atof(argv[3]);
    // the first argument is the symmetry
    if (argv[1][0] == 'C') //Cylindrical
    {
        source = source_cylindrical;
    }
    else
    {
        if (argv[1][0] == 'S')
        {
            source = source_spherical;
        }
        else
        {
            cout << "Invlaid symmetry\n";
            exit(1);
        }
    }

    if (argv[2][0] == 'N') // Nothing
    {
        tStop = 1.0;

        rho_l = 1.0;
        v_l = 0.0;
        p_l = 1.0;
        rho_r = rho_l;
        v_r = v_l;
        p_r = p_l;
    }
    else
    {
        if (argv[2][0] == 'E') // Explosion
        {
            tStop = 0.25;

            rho_l = 1.0;
            v_l = 0.0;
            p_l = 1.0;
            rho_r = 0.125;
            v_r = 0.0;
            p_r = 0.1;
        }
        else
        {
            if (argv[2][0] == 'I') // Implosion
            {
                tStop = 0.25;

                rho_l = 0.125;
                v_l = 0.0;
                p_l = 0.1;
                rho_r = 1.0;
                v_r = 0.0;
                p_r = 1.0;

            }
            else
            {
                std::cout << "invalid test case\n";
                exit(1);
            }
        }
    }
    
    dx = (x1 - x0) / nPoints;

    set_size(u);
    //set_size(u_copy);
    set_size(u_minus);
    set_size(u_plus);
    set_size(f);

    set_initial_value(u, padding);
    //set_initial_value(u_copy, padding);
    print_result(u);
    std::cout << "\t\n";

    double t = tStart;
    double dt;
    double a_max;
    step = 0;
 
    while (t < tStop)
    {
        a_max = maximal_v(u);
        dt = CFL * dx / a_max;
        // avoid overshooting the final time
        if ((t + dt) > tStop)
        {
            dt = tStop - t;
        }     

        ODE_Step(u, dx, 0.5*dt);
        Conservative_Step(u, f, dx, dt);
        ODE_Step(u, dx, 0.5*dt);

        /*
        ODE_Step(u_copy, dx, dt);
        Conservative_Step(u_copy, f, dx, dt);

        for (int position = padding; position < nPoints+padding; position++)
        {
            for (int index = 0; index < num_variables; index++)
            {
                u[position][index] = 0.5*(u[position][index] + u_copy[position][index]);
                u_copy[position][index] = u[position][index];
            }
        }

        boundary_condition(u, padding);
        boundary_condition(u_copy, padding);
        */
        

        //cout << "done " << step << "\n";
        //print_result(u);
        t = t + dt;
        step += 1;
    }

    print_result(u);
    //debug_result(u);
    return 0;
}
