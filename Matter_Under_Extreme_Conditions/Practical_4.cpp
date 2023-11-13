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
double x0 = 0;
double x1 = 1;
double tStart = 0.0;
double tStop;
double a = 1;
double dx = (x1 - x0)/nPoints;
double C = 0.8;
double gas_coef = 1.4;
double gamma = 7.15;
double p_inf = 300000000.0;
double Gamma = 0.25;
double epsiolon_delta = 0.0;
double rho_0 = 1840.0;
double A = 854500000000.0;
double B = 20500000000.0;
double R_1 = 4.6;
double R_2 = 1.35; 
int num_dimensions = 1;
int num_variables = num_dimensions + 2;
int step = 0;
double v_l, v_r, rho_l, rho_r, p_l, p_r;
bool use_halfstep = true;
//define a vector to store discretised data
vector<vector<double>> u; //u[0] and u[nPoints+1] are padding
vector<vector<double>> f; //f[i] holds the flux from u[i] to u[i+1]
vector<vector<double>> center_deltas;
vector<vector<double>> between_deltas;
vector<vector<double>> limiters;
vector<vector<double>> u_left;
vector<vector<double>> u_right;
vector<vector<double>> uPlusOne_left;
vector<vector<double>> uPlusOne_right;

vector<double> (flux)(const vector<double> &);

double (*pressure)(const vector<double> & state);

double (*speed_of_sound)(const vector<double> & state);


void boundary_condition(vector<vector<double>> &values, int pad)
{
    for (int i = 0; i < pad; i++)
    {
        for (int j = 0; j < num_variables; j++)
        {
            values[i][j] = values[pad][j];
            values[pad + nPoints + i][j] = values[pad + nPoints - 1][j];
        }
    }
}

vector<double> (*nummerical_flux)(const vector<double>&, const vector<double>&, double, double);

vector<vector<double>> PDE_Step(const vector<vector<double>> &current, const vector<vector<double>> &flux, double time_step, int pad)
{
    vector<vector<double>> uPlus1 = current;
    for (int index = pad; index < (nPoints + pad); index++)
    {
        for (int j = 0; j < num_variables; j++)
        {
            uPlus1[index][j] = current[index][j] - ((time_step / dx)*(flux[index][j] - flux[index-1][j]));
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
        result[i] = result[i] / 2.0;
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
    //std::cout << velocity_half_step[0] << " " << velocity_half_step[1] << " " << velocity_half_step[2] << "\n";
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
    double E_l = p_l/(gas_coef - 1.0) + (rho_l*v_l*v_l/2.0);
    double rho_v_r = rho_r*v_r;
    double E_r = p_r/(gas_coef - 1.0) + (rho_r*v_r*v_r/2.0);
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
        epsilon = e;
    }
    else
    {
        v = rho_v / rho;
        epsilon = E - (0.5*rho_v*v);
    }
    return epsilon;
}

double ideal_pressure(const vector<double> & state)
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
        /*
        debug_result(u);
        debug_result(u_left);
        debug_result(u_right);
        std::cout << "Negative pressure\n";
        std::cout << "p = " << p << " E = " << E << " rho_v = " << rho_v << " rho = " << rho << " step = " << step << "\n";
        */
        exit(3);
    }
    
    return p;
}

double ideal_speed_of_sound(const vector<double> & state)
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

double stiffened_pressure(const vector<double> & state)
{
    double rho_v = state[1];
    double rho = state[0];
    double E = state[2];
    double p;
    if (rho == 0.0)
    {
        /*
        std::cout << "Zero density when calculating the pressure\n";
        std::cout << rho << " " << rho_v << "\n";
        exit(2);
        */
       p = (((gamma - 1)*E) - (gamma * p_inf));
    }
    else
    {
        p = (((gamma - 1)*(E - (rho_v*rho_v)/(2*rho))) - (gamma * p_inf));
    }

    return p;
}

double stiffened_speed_of_sound(const vector<double> & state)
{
    double rho = state[0];
    double p = pressure(state);
    double cs_2;
    if (p < ((-1.0) * p_inf))
    {
        std::cout << "Too low pressure when calculating the speed of sound\n";
        std::cout << p << "\n";
        exit(3);
    }
    else
    {
        cs_2 = gamma * (p + p_inf) / rho;
        return (sqrt(cs_2));
    }
}

double p_ref(double rho)
{
    double exponent1, exponent2;
    exponent1 = (-1.0)*R_1*rho_0 / rho;
    exponent2 = (-1.0)*R_2*rho_0 / rho;
    return (A*exp(exponent1) + B*exp(exponent2));

}

double epsilon_ref(double rho)
{
    double exponent1, exponent2;
    double result;
    exponent1 = (-1.0)*R_1*rho_0 / rho;
    exponent2 = (-1.0)*R_2*rho_0 / rho;
    result = (A/R_1)*exp(exponent1) + (B/R_2)*exp(exponent2);
    result = (result/rho_0) - epsilon_delta;
    return result;

}

double JWL_pressure(const vector<double> & state)
{
    double rho_v = state[1];
    double rho = state[0];
    double E = state[2];
    double p;
    p = p_ref(rho) + Gamma*rho*(internal_energy(state) - epsilon_ref(rho));
    
    return p;
}

double JWL_speed_of_sound(const vector<double> & state)
{
    double rho = state[0];
    double p = pressure(state);
    double exponent1, exponent2;
    double epsilon = internal_energy(state);
    double result;
    exponent1 = (-1.0)*R_1*rho_0 / rho;
    exponent2 = (-1.0)*R_2*rho_0 / rho;
    result = ((-1.0)*exponent1*A/rho)*exp(exponent1) + ((-1.0)*exponent2*B/rho)*exp(exponent2);
    result += Gamma*(epsilon - epsilon_ref(rho));
    result += Gamma*(p/rho + A/R_1*exponent1*exp(exponent1) + B/R_2*exponent2*exp(exponent2));
    return sqrt(result);
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

void evaluate_flux(const vector<vector<double>> & values_left, const vector<vector<double>> & values_right, vector<vector<double>> &flu, double space_step, double time_step)
{
    for (int i = 0; i < (nPoints + (2*padding) - 1); i++)
    {
        flu[i] = nummerical_flux(values_right[i], values_left[i+1], space_step, time_step);
    }
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

void find_deltas(const vector<vector<double>> & values, vector<vector<double>> & between, vector<vector<double>> & center)
{
    for (int i = 1; i < nPoints + (2*padding)-1; i++)
    {
        for (int j = 0; j < num_variables; j++)
        {
            center[i][j] = 0.5*(values[i+1][j] - values[i-1][j]);
        } 
    }
    for (int j = 0; j < num_variables; j++)
    {
        center[0][j] = 0.0;
        center[nPoints + (2*padding) - 1][j] = 0.0;
    }
}

void (*find_limiters)(const vector<vector<double>> & between, vector<vector<double>> & limiters);

void find_limiters_energy(const vector<vector<double>> & values, vector<vector<double>> & limiters)
{
    double r;
    double val;
    for (int i = padding; i < nPoints+padding; i++)
    {
        if (values[i+1][num_variables-1] != values[i][num_variables-1])
            {
                //std::cout << between[i + padding][j] << std::endl;
                r = ((values[i][num_variables] - values[i-1][num_variables-1]) / (values[i+1][num_variables-1] - values[i][num_variables-1]));
                val = limiter(r);
            }
            else
            {
                val = 0.0;
            }
        for (int j = 0; j < num_variables; j++)
        {
            limiters[i][j] = val;
        }
    }
    for (int i = 0; i < padding; i++)
    {
        for (int j = 0; j < num_variables; j++)
        {
            limiters[i][j] = 0.0;
            limiters[i + nPoints + padding][j] = 0.0;
        }
    }
}

int main(int argc, char *argv[])
{
    nummerical_flux = FORCE_flux;
    //limiter = zero_limiter;
    limiter = Van_Leer;
    find_limiters = find_limiters_energy;

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
    u_left.resize(nPoints + 2*padding);
    u_right.resize(nPoints + 2*padding);
    uPlusOne_left.resize(nPoints + 2*padding);
    uPlusOne_right.resize(nPoints + 2*padding);
    f.resize(nPoints + 2*padding - 1);
    center_deltas.resize(nPoints + 2*padding);
    limiters.resize(nPoints + 2*padding);
    between_deltas.resize(nPoints + 2*padding - 1);

    for (int i = 0; i < (nPoints+(2*padding)-1); i++)
    {
        u[i].resize(num_variables);
        f[i].resize(num_variables);
        between_deltas[i].resize(num_variables);
        limiters[i].resize(num_variables);
        center_deltas[i].resize(num_variables);
        u_left[i].resize(num_variables);
        u_right[i].resize(num_variables);
        uPlusOne_left[i].resize(num_variables);
        uPlusOne_right[i].resize(num_variables);
    }
    int last_index = nPoints+(2*padding)-1;
    u[last_index].resize(num_variables);
    limiters[last_index].resize(num_variables);
    center_deltas[last_index].resize(num_variables);
    u_left[last_index].resize(num_variables);
    u_right[last_index].resize(num_variables);
    uPlusOne_left[last_index].resize(num_variables);
    uPlusOne_right[last_index].resize(num_variables);

    set_initial_value(u, padding);
    print_result(u);
    std::cout << "\t\n";

    double t = tStart;
    double dt;
    double a_max;
    step = 1;
    vector<double> flux_left, flux_right;
    flux_left.resize(num_variables);
    flux_right.resize(num_variables);
    while (t < tStop)
    {
        a_max = maximal_v(u);
        dt = C * dx / a_max;
        // avoid overshooting the final time
        if ((t + dt) > tStop)
        {
            dt = tStop - t;
        }     
        //std::cout << dt << std::endl; 

        find_deltas(u, between_deltas, center_deltas);
        //std::cout << "Deltas ";
        //find_limiters(between_deltas, limiters);
        find_limiters(u, limiters);
        //std::cout << "limiters\n";
        //debug_result(limiters);

        double  change;
        for (int i = 0; i < (nPoints + (2*padding)); i++)
        {
            for (int j = 0; j < num_variables; j++)
            {
                change = 0.5*(limiters[i][j] * center_deltas[i][j]);
                //std::cout << change << " ";
                u_left[i][j] = u[i][j] - change; 
                u_right[i][j] = u[i][j] + change;
            }
            //std::cout << std::endl;
        }
        //debug_result(u);
        //debug_result(u_left);
        if (use_halfstep)
        {
            for (int i = 0; i < nPoints; i++)
            {
                flux_left = flux(u_left[i + padding]);
                flux_right = flux(u_right[i + padding]);
                for (int j = 0; j < num_variables; j++)
                {
                    change = (flux_right[j] - flux_left[j]) * dt / dx / 2.0;
                    uPlusOne_left[i+padding][j] = u_left[i+padding][j] - change;
                    uPlusOne_right[i+padding][j] = u_right[i+padding][j] - change;
                }
            }
            for (int i = 0; i < padding; i++)
            {
                for (int j = 0; j < num_variables; j++)
                {
                    uPlusOne_left[i][j] = u_left[i][j];
                    uPlusOne_right[i][j] = u_right[i][j];
                    uPlusOne_right[nPoints+padding + i][j] = u_right[nPoints+padding+i][j];
                    uPlusOne_left[nPoints+padding + i][j] = u_left[nPoints+padding+i][j];
                }
            }
            u_left = uPlusOne_left;
            u_right = uPlusOne_right;  
        }
        //boundary_condition(u_left, padding);
        //boundary_condition(u_right, padding);
        //debug_result(u_left);
        //debug_result(u_right);

        evaluate_flux(u_left, u_right, f, dx, dt);

        //debug_result(f);

        //std::cout << "Flux evaluated\n";

        u = PDE_Step(u, f, dt, padding);

        //print_result(u);
        //std::cout << "PDE Step done\n";

        t = t + dt;
        step += 1;
    }

    print_result(u);
    //debug_result(u);
    return 0;
}
