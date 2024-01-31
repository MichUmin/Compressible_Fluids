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
int padding = 2;
double x0 = 0;
double x1 = 1;
double tStart = 0.0;
double tStop;
double a = 1;
double dx = (x1 - x0)/nPoints;
double C = 0.8;
double gas_coef = 1.4;
double p_inf = 300000000.0;
double Gamma = 0.25;
double epsilon_delta = 0.0;
double rho_0 = 1840.0;
double A = 854500000000.0;
double B = 20500000000.0;
double R_1 = 4.6;
double R_2 = 1.35; 
int num_dimensions = 1;
int num_variables = num_dimensions + 2;
int step = 0;
double v_l, v_r, rho_l, rho_r, p_l, p_r, rho_v_l, rho_v_r, epsilon_l, epsilon_r, E_l, E_r;
bool use_halfstep = true;

vector<vector<double>> u; //u[0] and u[nPoints+1] are padding
vector<vector<double>> f; //f[i] holds the flux from u[i] to u[i+1]
vector<vector<double>> u_minus;
vector<vector<double>> u_plus;

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
    /*
    double rho_v_l = rho_l*v_l;
    double E_l = p_l/(gas_coef - 1.0) + (rho_l*v_l*v_l/2.0);
    double rho_v_r = rho_r*v_r;
    double E_r = p_r/(gas_coef - 1.0) + (rho_r*v_r*v_r/2.0);
    */
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
       p = (((gas_coef - 1)*E) - (gas_coef * p_inf));
    }
    else
    {
        p = (((gas_coef - 1)*(E - (rho_v*rho_v)/(2*rho))) - (gas_coef * p_inf));
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
        cs_2 = gas_coef * (p + p_inf) / rho;
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

double p_ref_prime(double rho)
{
    double exponent1, exponent2;
    exponent1 = (-1.0)*R_1*rho_0 / rho;
    exponent2 = (-1.0)*R_2*rho_0 / rho;
    return ((A*R_1*rho_0 / (rho*rho))*exp(exponent1) + (B*R_2*rho_0 / (rho*rho))*exp(exponent2));
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

double epsilon_ref_prime(double rho)
{
    double exponent1, exponent2;
    double result;
    exponent1 = (-1.0)*R_1*rho_0 / rho;
    exponent2 = (-1.0)*R_2*rho_0 / rho;
    result = (A/(rho*rho))*exp(exponent1) + (B/(rho*rho))*exp(exponent2);
    return result;
}

double JWL_pressure(const vector<double> & state)
{
    double rho_v = state[1];
    double rho = state[0];
    double E = state[2];
    double eps = (E - (0.5*rho_v*rho_v/rho))/rho;
    double p;
    p = p_ref(rho) + Gamma*rho*(eps - epsilon_ref(rho));
    
    return p;
}

/*
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
*/

double JWL_speed_of_sound(const vector<double> & state)
{
    double rho = state[0];
    if (rho == 0.0)
    {
        cout << "Zero density when evaluating speed of sound\n";
        exit(2);
    }
    double p = pressure(state);
    double result;
    result = p*Gamma/rho + p_ref_prime(rho) - Gamma*rho*epsilon_ref_prime(rho) + (p - p_ref(rho))/rho;
    if (result <= 0.0)
    {
        cout << "Negative cs^2 at step: " << step << "\n";
        exit(2);
    }
    return sqrt(result);
}


double maximal_v(const vector<vector<double>>& values)
{
    int number = values.size();
    double result = 0.0;
    double v_c;
    for (int i = padding; i < nPoints + padding; i++)
    {
        v_c = speed_of_sound(values[i]);
        if (values[i][0] >= 0.0000001)
        {
            v_c += abs(values[i][1]) / values[i][0]; // rho_v / rho
        }
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
/*
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
*/

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

void SLICK_Step(vector<vector<double>> &current, vector<vector<double>> &fluxes, double space_step, double time_step)
{
    int next_x = 1;

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
                change = (flux_plus[j] - flux_minus[j]) * time_step / space_step / 2.0;
                u_minus[index_x][j] -= change;
                u_plus[index_x][j] -= change;
            }
        }
        boundary_condition(u_minus, padding);
        boundary_condition(u_plus, padding);
    }

    for (int i = 0; i < (nPoints + (2*padding) - 1); i++)
    {
        fluxes[i] = FORCE_flux(u_plus[i], u_minus[i + 1], space_step, time_step);
    }
    
    for (int index_x = padding; index_x < (nPoints + padding); index_x++)
    {
        for (int j = 0; j < num_variables; j++)
        {
            current[index_x][j] = current[index_x][j] - ((time_step / space_step)*(fluxes[index_x][j] - fluxes[index_x - 1][j]));
        }
    }
    boundary_condition(current, padding);
}

int main(int argc, char *argv[])
{
    nummerical_flux = FORCE_flux;
    limiter = zero_limiter;
    //limiter = Van_Leer;
    //find_limiters = find_limiters_energy;
    slope_reconstruction = slope_reconstruction_energy;

    if (argc < 2)
    {
        std::cout << "Remember to specify the test case you are running\n";
        exit(1);
    }
    // the first argument is the test case
    if (argv[1][0] == 'J') // JWL
    {
        tStop = 0.000012;
        pressure = JWL_pressure;
        speed_of_sound = JWL_speed_of_sound;

        rho_l = 1700.0;
        v_l = 0.0;
        p_l = 1000000000000.0;
        rho_v_l = rho_l * v_l;
        rho_r = 1000.0;
        v_r = 0.0;
        p_r = 50000000000.0;
        rho_v_r = rho_r * v_r;
        epsilon_l = ((p_l - p_ref(rho_l)) / (Gamma * rho_l))+ epsilon_ref(rho_l);
        E_l = rho_l * (epsilon_l + 0.5*v_l*v_l);
        epsilon_r = ((p_r - p_ref(rho_r)) / (Gamma * rho_r)) + epsilon_ref(rho_r);
        E_r = rho_r * (epsilon_r + 0.5*v_r*v_r);
    }
    else
    {
        if (argv[1][9] == '1') // Stiffened1
        {
            tStop = 0.00005;
            pressure = stiffened_pressure;
            speed_of_sound = stiffened_speed_of_sound;
            gas_coef = 7.15;

            rho_l = 1500.0;
            v_l = 0.0;
            p_l = 3000.0 * 101325.0;
            rho_r = 1000.0;
            v_r = 0.0;
            p_r = 101325.0;
            rho_v_l = rho_l * v_l;
            rho_v_r = rho_r * v_r;
            E_l = (p_l + gas_coef*p_inf)/(gas_coef - 1.0) + 0.5*rho_l*v_l*v_l;
            E_r = (p_r + gas_coef*p_inf)/(gas_coef - 1.0) + 0.5*rho_r*v_r*v_r;
        }
        else
        {
            if (argv[1][9] == '2') // Stiffened2
            {
                tStop = 0.00005;
                pressure = stiffened_pressure;
                speed_of_sound = stiffened_speed_of_sound;
                gas_coef = 7.15;

                rho_l = 1000.0;
                v_l = -100.0;
                p_l = 2.0 * 101325.0;
                rho_r = 1000.0;
                v_r = 100.0;
                p_r = 2.0*101325.0;
                rho_v_l = rho_l * v_l;
                rho_v_r = rho_r * v_r;
                E_l = (p_l + gas_coef*p_inf)/(gas_coef - 1.0) + 0.5*rho_l*v_l*v_l;
                E_r = (p_r + gas_coef*p_inf)/(gas_coef - 1.0) + 0.5*rho_r*v_r*v_r;
            }
            else
            {
                std::cout << "invalid test case\n";
                exit(1);
            }
        }
    }
    
    int nSteps = atoi(argv[2]);

    dx = (x1 - x0) / nPoints;
    u.resize(nPoints + 2*padding);
    u_minus.resize(nPoints + 2*padding);
    u_plus.resize(nPoints + 2*padding);
    f.resize(nPoints + 2*padding - 1);

    for (int i = 0; i < (nPoints+(2*padding)-1); i++)
    {
        u[i].resize(num_variables);
        f[i].resize(num_variables);
        u_minus[i].resize(num_variables);
        u_plus[i].resize(num_variables);
    }
    int last_index = nPoints+(2*padding)-1;
    u[last_index].resize(num_variables);
    u_minus[last_index].resize(num_variables);
    u_plus[last_index].resize(num_variables);

    set_initial_value(u, padding);
    print_result(u);
    std::cout << "\t\n";

    double t = tStart;
    double dt;
    double a_max;
    step = 0;
 
    while ((t < tStop) && (step < nSteps))
    {
        a_max = maximal_v(u);
        dt = C * dx / a_max;
        // avoid overshooting the final time
        if ((t + dt) > tStop)
        {
            dt = tStop - t;
        }     
        
        SLICK_Step(u, f, dx, dt);

        t = t + dt;
        if (step % 100 == 0)
        {
            std::cout << step << " " << t <<"\n";
        }
        step += 1;
    }

    print_result(u);
    //debug_result(u);
    return 0;
}
