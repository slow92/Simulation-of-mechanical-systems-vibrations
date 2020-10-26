#ifndef SIMULATIONS_V1_H
#define SIMULATIONS_V1_H

#include "iSimulations.h"
#include <vector>

typedef std::vector<double> double_vector;

class Simulations_v1
: public iSimulations
{
    public:
        static Simulations_v1& get_singleSimulation()
        {
            static Simulations_v1 singleSimulation;
            return singleSimulation;
        }

        virtual bool one_freedom_degree() override;
        virtual bool two_freedom_degree() override;
        virtual bool tree_freedom_degree()override;
        virtual bool save_data()override;

    protected:

    private:
        Simulations_v1();
        virtual ~Simulations_v1();

        const double Gravity = 9.81;
        double dt;

        bool reset_vector();
        bool init_vector(const int& freedom_degree_number, const long& iteration_number);
        double_vector v_u2, v_v2, v_a2;
        double_vector v_u3, v_v3, v_a3;
        double_vector v_u4, v_v4, v_a4;
};

#endif // SIMULATIONS_V1_H
